      module b2mod_blnc
      use b2mod_types
      use b2mod_diag
      use b2mod_time
      use b2mod_b2cmpa
      implicit none

c----------------------------------------------------------------------
c*** debug_flags usage:
c***  90  >0 : debugging printout on
c***  90, 91 : 'non-default surfaces' (90) to (91)
c***  92, 93 : strata (92) to (93); 0 means the sum over all strata
c***  94, 95 : more detailed output for strata (94) to (95)
c----------------------------------------------------------------------
c*** Output data (for each species):
c***
c***    header   array
c***                                        --- Global particle balance
c***   ion_core  flxic  ion flow across the core boundary
c***   ntr_core  flxnc  neutral flow across the core boundary
c***   ion_targ  flxit  ion-driven net flow through the targets
c***   ntr_targ  flxnt  neutral-driven net flow through the targets
c***   ion_wall  flxiw  ion-driven net flow through the side walls
c***   ntr_wall  flxnw  neutral-driven net flow through the side walls
c***   ntr_puff  flxnp  flow from the gas puff
c***   flux_tot  flxps  total flow coming to the volume
c***   flux_ttr  flxrl  the same in terms of the inverse time scales
c***   flx_d/dt  flxid  total flow equivalent to the time derivative
c***   flux_bln  flxbl  error in particle balance (should be 0)
c***   flux_blr  flxbr  the same in terms of the inverse time scales
c***
c***                                        --- Ion balance in B2
c***                                            (fluxes in terms of the
c***                                                inverse time scales)
c***   ions_tot  blnna   the total number of ions
c***   bln_d/dt  blnnd   the time derivative
c***   bln_core  blnnc   ion flow from the core
c***   bln_targ  blnnt   ion flow to the targets
c***   bln_wall  blnnw   ion flow to the side walls
c***   bln_srce  blnns   integral sources
c***   bln_blnc  blnne   error in particle balance
c***   bln_blnr  blnnr   the same divided by the sum of the abs. values
c***                     of the other terms (relative error)
c----------------------------------------------------------------------

c
c*** lv    :  version number to appear in the file
c*** np_tmp:  number of output data for each species
c
      integer lv, np_tmp
      parameter (lv=961229, np_tmp=20)
      integer ltrc, lun
      parameter (ltrc=ntrc_bln)
      integer, save :: nc

      real (kind=R8), save :: stim = 0.0_R8
      real (kind=R8), allocatable :: hlpit_lp(:),hlpit_rp(:),hlpit_pl(:)
      real (kind=R8), allocatable, save :: hlp_neu(:,:)
      real (kind=R8), allocatable :: hlpic_s(:),hlpiw_s(:),hlpit_s(:)
c
      real (kind=R8), allocatable ::
     *     flxic_tmp(:),flxnc_tmp(:),flxit_tmp(:),
     ,     flxnt_tmp(:),flxiw_tmp(:),flxnw_tmp(:),
     ,     flxnp_tmp(:),flxps_tmp(:),blnna_tmp(:),
     ,     blnnd_tmp(:),blnnc_tmp(:),blnnt_tmp(:),
     ,     blnnw_tmp(:),blnns_tmp(:),blnne_tmp(:),
     ,     blnnr_tmp(:),
c*** from 960727 on
     ,     flxid_tmp(:),flxbl_tmp(:),flxbr_tmp(:),
c*** from 960808 on
     ,     flxrl_tmp(:)

      character f*(lenfnm), hlp_h(np_tmp)*8, hlp_frm*80

      save f, hlp_frm, hlp_h, lun,
     , hlpit_lp,hlpit_rp,hlpit_pl,hlpic_s,hlpiw_s,hlpit_s,
     , flxic_tmp,flxnc_tmp,flxit_tmp,flxnt_tmp,
     , flxiw_tmp,flxnw_tmp,flxnp_tmp,flxps_tmp,
     , blnna_tmp,blnnd_tmp,blnnc_tmp,blnnt_tmp,
     , blnnw_tmp,blnns_tmp,blnne_tmp,blnnr_tmp,
     , flxid_tmp,flxbl_tmp,flxbr_tmp,
     , flxrl_tmp

c*** Column headers:
      data hlp_h /
     1          'ion_core','ntr_core','ion_targ','ntr_targ','ion_wall',
     2          'ntr_wall','ntr_puff','flux_tot','ions_tot','bln_d/dt',
     3          'bln_core','bln_targ','bln_wall','bln_srce','bln_blnc',
     4          'bln_blnr','flx_d/dt','flux_bln','flux_blr','flux_ttr'/

      contains

      subroutine b2blnc(nx,ny,ns,dtim)
      use b2mod_geo
      use b2mod_plasma
      use b2mod_sources
      use b2mod_indirect
      use b2mod_numerics_namelist
      use b2mod_neutrals_namelist
      use b2mod_boundary_namelist
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif
#include "DIMENSIONS.F"
c
c  version : 19.07.2001 18:46
c
c======================================================================
c*** Computes the global particle flows across the boundaries in
c*** B2-Eirene and the global particle balance in B2
c----------------------------------------------------------------------
      implicit none
c
      integer nx, ny, ns
      real (kind=R8) :: dtim
c
      integer i, j, k, l, m, n, ic, ii, is, ix, iy, js, ks, ls,
     ,        kks
      logical first, hlp_lw, hlp_lt
c<<<
      logical hlp_tt
      real (kind=R8) :: hlp_wldpp
      real (kind=R8) :: hlpit_l,hlpit_r
      real (kind=R8) :: hlp_pls(nbcd,ns)
      integer kst_hlp
      integer use_recyceir, use_eirene
      real (kind=R8) :: recyctmp
c>>>
c
      real (kind=R8) ::
     *     hlpic,hlpiw,hlpit,hlpia,hlpid,hlpnc,hlpnw,hlpnt,hlpis,hlpnp,
     ,     hlpdt,hlp,hlpids
c----------------------------------------------------------------------
c
      character hlp_prt*26
c
c----------------------------------------------------------------------
c*** Taken from the common blocks:
c***
c*** lcbs: boundary segments representing the core boundary in B2
c*** lcns: 'non-default surfaces' for the core boundary (Eirene)
c*** ltns: 'non-default surfaces' for the targets (Eirene)
c*** lsns: list of Eirene surfaces comprising the recycling sources
c*** msns: msns(1,k): index of first surface in lsns list for Eirene
c***                  recycling source k
c***       msns(2,k): number of surfaces for Eirene recycling source k
c*** mlcmp: composition of molecules.
c***        Molecule of kind j contains mlcmp(i,j) atoms of kind i
c*** gpfc: composition of gas puff sources.
c***       gpfc(i,j)=(fraction of atoms of kind i in stratum number j)
c*** sclstra: scaling factors for species and strata
c*** lmolscl: pointer to the species with which these molecules scale
c----------------------------------------------------------------------
*   ..initialisation
      data use_recyceir /1/, use_eirene /0/
c
      data first / .true. /
c----------------------------------------------------------------------
      save i, k, l, m, n, ic, is, ix, iy,
     1 js, ks, ls, kks, use_recyceir,
     2 first, hlp_lw, hlp_lt, hlp_tt,
     3 hlpit_l,hlpit_r,
     4 hlpic,hlpiw,hlpit,hlpia,hlpid,hlpnc,hlpnw,hlpnt,hlpis,hlpnp,
     5 hlpdt,hlp,hlpids
c======================================================================

      if(first) then
c<<<
        call ipgeti ('eirene_use_recyceir', use_recyceir)
        call ipgeti ('b2mndr_eirene', use_eirene)
        print *,'=== Regions of the grid:'
        l=ny+1
        n=mod(l-1,10)+1
        i=(nx+2)/10+1
        write(hlp_prt,'(a,i3,a)') '(1x,i4,1x,i1,',i,'(1x,10i1))'
        do k=0,ny-1
          m=l
          l=max(0,m-n+1)
          n=10
          do j=m-1,l-1,-1
            write(*,hlp_prt) j,(region(i,j,0),i=-1,nx)
          end do
          l=l-1
        end do
        write(*,'()')
c>>>
c
c*** Initialisation
c
        do j=1,nstrat
          hlp=0.0_R8
          do i=1,nnatmi
            hlp=hlp+gpfc(i,j)
          end do
          if(hlp.gt.0.0_R8) then
            do i=1,nnatmi
              gpfc(i,j)=gpfc(i,j)/hlp
            end do
          else
            do i=1,nnatmi
              gpfc(i,j)=0.0_R8
            end do
          end if
        end do
        do j=1,DEF_NSRFS*nstrat
          if(lsns(j).lt.0) lsns(j)=nlim-lsns(j)
        end do
        first=.false.
      end if

c*** Re-scale the neutral fluxes

csw 2012 correct nlimps --> nnstsi+nnlimi (known after 1st Eirene call)
c      do i=1,nlimps
      do i=1,nnstsi+nnlimi
        do j=1,nnatmi
          wldna(i,j,0)=0.0_R8
          wldra(i,j,0)=0.0_R8
          wldpa(i,j,0)=0.0_R8
          do k=1,nstrat
            wldna(i,j,0)=wldna(i,j,0)+wldna(i,j,k)*sclstra(j,k)
            wldra(i,j,0)=wldra(i,j,0)+wldra(i,j,k)*sclstra(j,k)
            wldpa(i,j,0)=wldpa(i,j,0)+wldpa(i,j,k)*sclstra(j,k)
          end do
        end do
        do j=1,nnmoli
          l=lmolscl(j)
          wldnm(i,j,0)=0.0_R8
          wldrm(i,j,0)=0.0_R8
          wldpm(i,j,0)=0.0_R8
          do k=1,nstrat
            if(l.gt.0 .and. l.le.nnatmi) then
              wldnm(i,j,0)=wldnm(i,j,0)+wldnm(i,j,k)*sclstra(l,k)
              wldrm(i,j,0)=wldrm(i,j,0)+wldrm(i,j,k)*sclstra(l,k)
              wldpm(i,j,0)=wldpm(i,j,0)+wldpm(i,j,k)*sclstra(l,k)
            else
              wldnm(i,j,0)=wldnm(i,j,0)+wldnm(i,j,k)
              wldrm(i,j,0)=wldrm(i,j,0)+wldrm(i,j,k)
              wldpm(i,j,0)=wldpm(i,j,0)+wldpm(i,j,k)
            end if
          end do
        end do
      end do

c<<<
      if(debug_flags(90).gt.0.and.use_eirene.ne.0) then
        print *
        print *,'%%% Eirene data in B2 %%%'
        do k=debug_flags(92),debug_flags(93)
          print *
          print *,'wldpp:',k
          write (hlp_prt,'(a,i3,a)') '(1p,i4,',nnplsi,'e9.2)'
          do j=nnlimi+debug_flags(90),nnlimi+debug_flags(91)
            write(*,hlp_prt) j,(wldpp(j,i,k),i=1,nnplsi)
          end do
          print *
          if (nnmoli.gt.0) then
            print *,'wldpeb, wldpa, wldpm, hlpit_pl:',k
          else
            print *,'wldpeb, wldpa, hlpit_pl:',k
          end if
          write (hlp_prt,'(a,i3,a)')
     .     '(1p,i4,',1+nnatmi+nnmoli+nspecies,'e9.2)'
          do j=nnlimi+debug_flags(90),nnlimi+debug_flags(91)
            ks=0
            do ls=1,nspecies
              hlpit_pl(ls)=0.0_R8
              do is=1,nfluids(ls)
                ks=ks+1
                hlpit_pl(ls)=hlpit_pl(ls)+wldpp(j,ks,k)
              end do
            end do
            write(*,hlp_prt) j,wldpeb(j,k),
     ,         (wldpa(j,i,k),i=1,nnatmi),(wldpm(j,i,k),i=1,nnmoli),
     ,                                       (hlpit_pl(i),i=1,nspecies)
           end do
         end do
       end if
c>>>

c*** The plasma flows across the boundaries

      ks=0
      do ls=1,nspecies
        hlpic =0.0_R8
        hlpit =0.0_R8
        hlpiw =0.0_R8
        hlpit_lp(ls)=0.0_R8
        hlpit_rp(ls)=0.0_R8
c<<<
        do ic=1,nbc
          hlp_pls(ic,ls)=0.0_R8
          if (bcchar(ic).eq.'S' .or. bcchar(ic).eq.'N' .or.
     .        bcchar(ic).eq.'E' .or. bcchar(ic).eq.'W') kst_hlp=ic
        end do
c>>>
        do is=1,nfluids(ls)
          ks=ks+1
          do 1100 ic = 1, nbc
            if (bcchar(ic).eq.'S') then
              do j = 1, ncbs
                if(ic.eq.lcbs(j)) then
c***     core boundary
                  do i = 1, bc_list_size(ic)
                    hlpic=hlpic+
     .               fna(topix(bc_list_x(i,ic),bc_list_y(i,ic)),
     .                   topiy(bc_list_x(i,ic),bc_list_y(i,ic)),
     .                   1,ions_list(ks))
c<<<
                    hlp_pls(ic,ls)=hlp_pls(ic,ls)+
     .               fna(topix(bc_list_x(i,ic),bc_list_y(i,ic)),
     .                   topiy(bc_list_x(i,ic),bc_list_y(i,ic)),
     .                   1,ions_list(ks))
c>>>
                  end do
                  go to 1100
                end if
              end do
c***     wall in PFR
              do i = 1, bc_list_size(ic)
                hlpiw=hlpiw+
     .           fna(topix(bc_list_x(i,ic),bc_list_y(i,ic)),
     .               topiy(bc_list_x(i,ic),bc_list_y(i,ic)),
     .               1,ions_list(ks))
c<<<
                hlp_pls(ic,ls)=hlp_pls(ic,ls)+
     .           fna(topix(bc_list_x(i,ic),bc_list_y(i,ic)),
     .               topiy(bc_list_x(i,ic),bc_list_y(i,ic)),
     .               1,ions_list(ks))
c>>>
              end do
            else if (bcchar(ic).eq.'N') then
c***     outer wall
              do i = 1, bc_list_size(ic)
                hlpiw=hlpiw
     .           -fna(bc_list_x(i,ic),bc_list_y(i,ic),1,ions_list(ks))
c<<<
                hlp_pls(ic,ls)=hlp_pls(ic,ls)-
     .           fna(bc_list_x(i,ic),bc_list_y(i,ic),1,ions_list(ks))
c>>>
              end do
            else if (bcchar(ic).eq.'W') then
c***     left target
              do i = 1, bc_list_size(ic)
                hlpit=hlpit+
     .           fna(rightix(bc_list_x(i,ic),bc_list_y(i,ic)),
     .               rightiy(bc_list_x(i,ic),bc_list_y(i,ic)),
     .               0,ions_list(ks))
                hlpit_lp(ls)=hlpit_lp(ls)+
     .           fna(rightix(bc_list_x(i,ic),bc_list_y(i,ic)),
     .               rightiy(bc_list_x(i,ic),bc_list_y(i,ic)),
     .               0,ions_list(ks))
c<<<
                hlp_pls(ic,ls)=hlp_pls(ic,ls)+
     .           fna(rightix(bc_list_x(i,ic),bc_list_y(i,ic)),
     .               rightiy(bc_list_x(i,ic),bc_list_y(i,ic)),
     .               0,ions_list(ks))
c>>>
              end do
            else if (bcchar(ic).eq.'E') then
c***     right target
              do i = 1, bc_list_size(ic)
                hlpit=hlpit-
     .           fna(bc_list_x(i,ic),bc_list_y(i,ic),0,ions_list(ks))
                hlpit_rp(ls)=hlpit_rp(ls)-
     .           fna(bc_list_x(i,ic),bc_list_y(i,ic),0,ions_list(ks))
c<<<
                hlp_pls(ic,ls)=hlp_pls(ic,ls)-
     .           fna(bc_list_x(i,ic),bc_list_y(i,ic),0,ions_list(ks))
c>>>
              end do
            end if
 1100     continue
        end do
        hlpic_s(ls)=hlpic
        hlpit_s(ls)=hlpit
        hlpiw_s(ls)=hlpiw
      end do
c<<<
      if(debug_flags(90).gt.0) then
c       write (hlp_prt,'(a,i2,a)') '(1x,a,1p,',2*nspecies,'e12.4)'
c       write (*,hlp_prt) '%%% hlpit_lp, hlpit_rp: ',
c     ,                           (hlpit_lp(i),hlpit_rp(i),i=1,nspecies)
        write (hlp_prt,'(a,i2,a)') '(1x,a,1p,',nspecies,'e12.4)'
        write (*,hlp_prt) '%%# hlpiw_s : ',(hlpiw_s(i),i=1,nspecies)
      end if
c>>>

      ks=0
      do ls=1,nspecies
        kks=eb2spcr(ls)+1
c<<<
        if(debug_flags(90).gt.0) then
          print *,'%%# Species # ',ls
          print *,'%%# eb2spcr(ls) = ',eb2spcr(ls)
          print *,'%%# nfluids(ls) = ',nfluids(ls)
        end if
c>>>
c
c*** Ion fluxes and average densities
c
        hlpid =0.0_R8
        hlpids=0.0_R8
        hlpia =0.0_R8
        hlpis =0.0_R8
c<<<
        hlpit_l =0.0_R8
        hlpit_r =0.0_R8
c>>>
        do is=1,nfluids(ls)
          ks=ks+1
c          print *,'%%% Fluid # ',is,ks     !###
c*** total number of particles and time derivative
          do iy=0,ny-1
            do ix=0,nx-1
              j=region(ix,iy,0)
              if(j.ne.0 .and.
     &           leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
                hlpia=hlpia+na(ix,iy,ions_list(ks))*vol(ix,iy)
                hlpis=hlpis+
     1               sna(ix,iy,0,ions_list(ks))+
     2               sna(ix,iy,1,ions_list(ks))*na(ix,iy,ions_list(ks))
                if (solveco(ions_list(ks),j)) then
                  hlpdt=dtim*dtco(ions_list(ks),j)*time_factor(ix,iy)
                else
                  hlpdt=0.0_R8
                endif
                if(hlpdt.gt.0.0_R8) then
                  hlpdt=vol(ix,iy)/hlpdt
                  hlpid=hlpid+na(ix,iy,ions_list(ks))*hlpdt
                  hlpids=hlpids+na0(ix,iy,ions_list(ks))*hlpdt
                end if
              end if
            end do
          end do
c<<<
c          hlpit_l=hlpit_l+wldpp(nnlimi+7,ks,0)
c          hlpit_r=hlpit_r+wldpp(nnlimi+8,ks,0)
c>>>
        end do

        hlpid=hlpid-hlpids
        hlpic=hlpic_s(ls)
        hlpit=hlpit_s(ls)
        hlpiw=hlpiw_s(ls)
c*** Correct the source for d/dt terms
        hlpis=hlpis+hlpid
c<<<
c        print '(1x,a,1p,4e12.4)',
c     ,      '%%% Ion fluxes: hlpit, hlpiw, hlpic, hlpid =',
c     ,                       hlpit, hlpiw, hlpic, hlpid
c        print *,'%%% hlpit_l     hlpit_r     hlpit_t'
c        print '(i3,a)',ls,' hlpit_lp    hlpit_rp    hlpit'
c        print '(1p,3e12.4)',hlpit_l,hlpit_r,hlpit_r+hlpit_l,
c     ,                                  hlpit_lp(ls),hlpit_rp(ls),hlpit
c>>>
        blnna_tmp(ls)=hlpia
        blnnd_tmp(ls)=hlpid/hlpia
        blnnc_tmp(ls)=hlpic/hlpia
        blnnt_tmp(ls)=hlpit/hlpia
        blnnw_tmp(ls)=hlpiw/hlpia
        blnns_tmp(ls)=hlpis/hlpia
        blnne_tmp(ls)=blnns_tmp(ls)+blnnc_tmp(ls)+blnnt_tmp(ls)+
     +                                      blnnw_tmp(ls)-blnnd_tmp(ls)
        blnnr_tmp(ls)=blnne_tmp(ls)/(abs(blnns_tmp(ls))+
     +        abs(blnnc_tmp(ls))+abs(blnnt_tmp(ls))+abs(blnnw_tmp(ls))+
     +                                              abs(blnnd_tmp(ls)))

c
c*** Neutral sources - assume the same order of species in B2 and Eirene
c
        hlpnc=0.0_R8
        hlpnp=0.0_R8
        hlpnw=0.0_R8
        hlpnt=0.0_R8
        hlp_neu=0.0_R8
c*** Flux across the core boundary
        do is=1,ncns
          js=lcns(is)+nnlimi
          do j=1,nnatmi
            if(latmscl(j).eq.ls) then
              hlpnc=hlpnc-wldna(js,j,0)+wldra(js,j,0)+wldpa(js,j,0)
              do ic=1,nnmoli
                if(mlcmp(j,ic).gt.0) hlpnc=hlpnc-mlcmp(j,ic)*
     .                (wldnm(js,ic,0)-wldrm(js,ic,0)-wldpm(js,ic,0))
              end do
            end if
          end do
        end do
c*** Flux across the targets
        do is=1,ntns
          js=ltns(is)+nnlimi
          do j=1,nnatmi
            if(latmscl(j).eq.ls) then
              hlpnt=hlpnt-wldna(js,j,0)+wldra(js,j,0)+wldpa(js,j,0)
              do ic=1,nnmoli
                if(mlcmp(j,ic).gt.0) hlpnt=hlpnt-mlcmp(j,ic)*
     .                (wldnm(js,ic,0)-wldrm(js,ic,0)-wldpm(js,ic,0))
              end do
            end if
          end do
        end do
c*** Flux across the side walls
        do js=1,nnlimi
          do j=1,nnatmi
            if(latmscl(j).eq.ls) then
              hlpnw=hlpnw-wldna(js,j,0)+wldra(js,j,0)+wldpa(js,j,0)
              do ic=1,nnmoli
                if(mlcmp(j,ic).gt.0) hlpnw=hlpnw-mlcmp(j,ic)*
     .           (wldnm(js,ic,0)-wldrm(js,ic,0)-wldpm(js,ic,0))
              end do
            end if
          end do
        end do
c*** Gas puffing source and recycling
        do ic=1,nstrat
          hlp_tt = debug_flags(90).gt.0 .and.
     .                ic.ge.debug_flags(94) .and. ic.le.debug_flags(95)
          if(use_recyceir.gt.0) then
            recyctmp=recyceir(ic)
          else
            recyctmp=recyc(kks,ic)
          endif
          if (crcstra(ic).eq.'C') then
            do j=1,nnatmi
              if(latmscl(j).eq.ls)
     >          hlpnp=hlpnp+userfluxparm(ic,1)*gpfc(j,ic)*recyctmp
            end do
          else
            hlp_lw = crcstra(ic).eq.'S' .or. crcstra(ic).eq.'N'
            hlp_lt = crcstra(ic).eq.'E' .or. crcstra(ic).eq.'W'
c<<<
            if(hlp_tt) print *,'%%% Recycling? - ',hlp_lt,hlp_lw
c>>>
            if(hlp_lw .or. hlp_lt) then
              l=msns(1,ic)
              m=msns(2,ic)
              hlp=0.0_R8
c<<<
              if(hlp_tt) write(*,'(1x,a,2i5,1x,a1,1x,3i5)')
     ,                           '%%# ls,ic,crcstra,l,m,kks : ',
     ,                                ls,ic,crcstra(ic),l,m,kks
c>>>
              do i=l,l+m-1
                k=lsns(i)
                if(k.lt.0) k=nnlimi-k
                if(k.gt.0) then
c<<<
                  if(hlp_tt) then
                    hlp_wldpp=0.0_R8
                    do ii=kks,kks+nfluids(ls)-1
                      hlp_wldpp=hlp_wldpp+wldpp(k,ii,ic)
                    end do
                    print '(1x,a,2i5,1p,40e12.4)',
     ,                 '%%# i,k,hlp,wldpp= ',
     ,                  i,k,hlp,hlp_wldpp
                  end if
c>>>
                  do j=1,nnatmi
                    if (latmscl(j).eq.ls) then
                      hlp=hlp+wldpa(k,j,ic)*sclstra(j,ic)
c<<<
                      if(hlp_tt)
     ,                 print '(1x,a,i5,a,i5,a,i5,a,1p,1e12.4)',
     ,                 '%%% wldpa(',k,',',j,',',ic,') =',
     ,                      wldpa(k,j,ic)*sclstra(j,ic)
c>>>
                    end if
                  end do
                  do j=1,nnmoli
                    do ii=1,nnatmi
                      if (latmscl(ii).eq.ls) then
                        hlp=hlp+mlcmp(ii,j)*wldpm(k,j,ic)*sclstra(ii,ic)
c<<<
                        if(hlp_tt) print '(1x,a,2i3,1p,40e12.4)',
     ,                       '%%# j,mlcmp,wldpm,sclstra= ',
     ,                       j,mlcmp(ii,j),wldpm(k,j,ic),sclstra(ii,ic)
                        if(hlp_tt .and. mlcmp(ii,j).gt.0)
     ,                   print '(1x,a,i5,a,i5,a,i5,a,1p,1e12.4)',
     ,                       '%%% wldpm(',k,',',j,',',ic,') =',
     ,                            wldpm(k,j,ic)*sclstra(ii,ic)
c>>>
                      end if
                    end do
                  end do
                end if
              end do
c<<<
              hlp_neu(ic,ls)=hlp*recyctmp
              if(hlp_tt) print '(1x,a,2l3,2i3,1p,4e12.4)',
     ,          '%%# hlp_lw,hlp_lt,kks,is,hlpiw,hlpit,hlp,recyc=',
     ,               hlp_lw,hlp_lt,kks,is,hlpiw,hlpit,hlp,recyctmp
c>>>
              if(hlp_lw) hlpiw=hlpiw+hlp*recyctmp
              if(hlp_lt) hlpit=hlpit+hlp*recyctmp
            end if
          end if
        end do
c<<<
        if(debug_flags(90).gt.0) then
          print '(1x,a,1p,5e12.4)',
     ,      '%%# Total fluxes: hlpit, hlpiw, hlpic, hlpid, hlpnp =',
     ,                                hlpit, hlpiw, hlpic, hlpid, hlpnp
          print *,'%%# Recycling sources. ls= ',ls
          print '(a4,8x,20(5x,i3,4x))',' %%#',(i,i=1,kst_hlp)
          write(hlp_prt,'(a,i3,a)') '(1x,a,i3,a,1p,',nstrat,'e12.5)'
          do j=1,nnatmi
            if (latmscl(j).eq.ls) then
              write(*,hlp_prt) '%%# sclrtio(',j,':)= ',
     ,                                       (sclrtio(j,i),i=1,nstrat)
              write(*,hlp_prt) '%%# sclstra(',j,':)= ',
     ,                                       (sclstra(j,i),i=1,nstrat)
            end if
          end do
          write(hlp_prt,'(a,i3,a)') '(1x,a,1p,',kst_hlp,'e12.5)'
          write(*,hlp_prt) '%%# plasma  ',(hlp_pls(i,ls),i=1,kst_hlp)
          write(hlp_prt,'(a,i3,a)') '(1x,a,1p,',nstrat,'e12.5)'
          write(*,hlp_prt) '%%# neutrals',(hlp_neu(i,ls),i=1,nstrat)
        end if
c>>>
        flxid_tmp(ls)=hlpid
        flxic_tmp(ls)=hlpic
        flxit_tmp(ls)=hlpit
        flxiw_tmp(ls)=hlpiw
        flxnc_tmp(ls)=hlpnc
        flxnt_tmp(ls)=hlpnt
        flxnw_tmp(ls)=hlpnw
        flxnp_tmp(ls)=hlpnp
c
c*** Sum of all the fluxes...
c
        flxps_tmp(ls)=flxic_tmp(ls)+flxnc_tmp(ls)+flxit_tmp(ls)+
     +          flxnt_tmp(ls)+flxiw_tmp(ls)+flxnw_tmp(ls)+flxnp_tmp(ls)
        flxrl_tmp(ls)=flxps_tmp(ls)/blnna_tmp(ls)
c
c*** ...minus the time derivative
c
        flxbl_tmp(ls)=flxps_tmp(ls)-flxid_tmp(ls)
        flxbr_tmp(ls)=flxbl_tmp(ls)/blnna_tmp(ls)
c
      end do
c
      if(ank_tracing.eq.0) return
      if(filedata(ltrc)) then
c
c*** Write the data in a separate file for plotting...
c
        if(nc.gt.np_tmp .or. nc.le.0) then
          write(*,*) 'b2blnc: unexpected error. nc =',nc
          filedata(ltrc)=.false.
        else
          write(lun,hlp_frm) tim,
     .      (flxic_tmp(i),flxnc_tmp(i),flxit_tmp(i),
     ,       flxnt_tmp(i),flxiw_tmp(i),flxnw_tmp(i),
     ,       flxnp_tmp(i),flxps_tmp(i),blnna_tmp(i),
     ,       blnnd_tmp(i),blnnc_tmp(i),blnnt_tmp(i),
     ,       blnnw_tmp(i),blnns_tmp(i),blnne_tmp(i),
     ,       blnnr_tmp(i),
     ,       flxid_tmp(i),flxbl_tmp(i),flxbr_tmp(i),
     ,       flxrl_tmp(i),i=1,nspecies)
        end if
      else
c
c*** ... or just print out the data
c
        write(*,*) 'particle_fluxes ',nspecies
        write(*,'(30(3x,a8,1x))') (hlp_h(i),i=1,np_tmp)
        do i=1,nspecies
          write(*,'(1p,30e12.3)')
     ,     flxic_tmp(i),flxnc_tmp(i),flxit_tmp(i),
     ,     flxnt_tmp(i),flxiw_tmp(i),flxnw_tmp(i),
     ,     flxnp_tmp(i),flxps_tmp(i),blnna_tmp(i),
     ,     blnnd_tmp(i),blnnc_tmp(i),blnnt_tmp(i),
     ,     blnnw_tmp(i),blnns_tmp(i),blnne_tmp(i),
     ,     blnnr_tmp(i),
     ,     flxid_tmp(i),flxbl_tmp(i),flxbr_tmp(i),
     ,     flxrl_tmp(i)
        end do
      end if
      return
      end subroutine b2blnc
c======================================================================
c
      subroutine b2blnci
      implicit none
      character, save :: tab
      integer :: i, j, nv
      character :: hlp_c
      logical :: ex
c----------------------------------------------------------------------
c*** Initialise the file for tracing data
c*** Only called with filedata(ntrc_bln)=.true.
c----------------------------------------------------------------------
c
      allocate(hlpit_lp(nspecies),hlpit_rp(nspecies),hlpit_pl(nspecies))
      allocate(hlp_neu(nstra,nspecies))
      allocate(hlpic_s(nspecies),hlpiw_s(nspecies),hlpit_s(nspecies))
      allocate(
     , flxic_tmp(nspecies),flxnc_tmp(nspecies),flxit_tmp(nspecies),
     , flxnt_tmp(nspecies),flxiw_tmp(nspecies),flxnw_tmp(nspecies),
     , flxnp_tmp(nspecies),flxps_tmp(nspecies),blnna_tmp(nspecies),
     , blnnd_tmp(nspecies),blnnc_tmp(nspecies),blnnt_tmp(nspecies),
     , blnnw_tmp(nspecies),blnns_tmp(nspecies),blnne_tmp(nspecies),
     , blnnr_tmp(nspecies),
c*** from 960727 on
     , flxid_tmp(nspecies),flxbl_tmp(nspecies),flxbr_tmp(nspecies),
c*** from 960808 on
     , flxrl_tmp(nspecies))
      call ipgetr ('b2mndr_stim', stim)
      if(.not.filedata(ntrc_bln)) return

c*** lun   :  logical unit number for the output file (=ntrc(ltrc))
c*** f     :  the associated file name
      tab=char(9)
      nc=0
      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)
      call find_file(f,ex)
      open(lun,file=f,err=90)
      rewind(lun)
      if(ex) read(lun,*,end=91) hlp_c
      goto 93
   91 ex = .false.
   93 continue
      rewind(lun)
      if(ex.and.stim.lt.0.0_R8) then
        read(lun,*,err=40) hlp_c,nv,j,nc
        read(lun,*) hlp_c
        go to 42
   40   read(lun,*,err=92) hlp_c,i,j
        lv_trc(ltrc)=lv
        if(i.lt.960727) then
          nc=16
        else
          go to 92
        end if
   42   if(nv.ne.lv) then
          write(*,'(a,i8,3a,i8)')
     ,     '!!! b2blnc: current version ',lv,
     ,     ' differs from that in the file ',trim(f),' : ',nv
          if(nv .lt. 960727) then
c*** It was a major change of the file format - impossible to continue
            write(*,*) 'No ',trim(f),' tracing file - ',
     ,                                'remove the existing file first!'
            filedata(ltrc)=.false.
            return
          end if
        end if
        if(nc.ne.np_tmp) then
          write(*,'(a,i4,3a,i4)')
     ,     '!!! b2blnc: current number of columns ',np_tmp,
     ,     ' differs from that in the file ',trim(f),' : ',nc
          if(nc.gt.np_tmp) then
            write(*,*) 'No ',trim(f),' tracing file - ',
     ,                                'remove the existing file first!'
            filedata(ltrc)=.false.
            return
          end if
          write(*,*) ' You should remove the file!'
        end if
        if(j.ne.nspecies) then
          write(*,'(a,i2,3a,i2)')
     ,     '!!! b2blnc: current number of species ',nspecies,
     ,      ' differs from that in the file ',trim(f),' : ',j
          write(*,*) 'No ',trim(f),' tracing file - ',
     ,                              'remove the existing file first!'
          filedata(ltrc)=.false.
          return
        end if
   50   read(lun,'(a1)',end=55) hlp_c
        go to 50
   55   backspace lun
      else
        nc=np_tmp
        write(lun,700) lv,nspecies,np_tmp,label
        write(hlp_frm,'(9h(5h#time,,i3,16h(1h'//tab//',a8,1h_,a2)))')
     ,                                                  np_tmp*nspecies
c       print *,'b2blnc: nspecies= ',nspecies                   !###
c       print *,'species_list :',(species_list(j),':',j=1,nspecies) !###
        write(lun,hlp_frm) ((hlp_h(i),species_list(j),i=1,np_tmp),
     ,                                                    j=1,nspecies)
      end if
      write(hlp_frm,'(10h(1p,e16.7,,i3,12h(1h'//tab//',e12.3)))')
     ,                                                  np_tmp*nspecies
      return
c======================================================================
   90 write(*,*) 'b2blnc: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      return
c
   92 write(*,*) 'b2blnc: wrong format of the file ',trim(f)
      filedata(ltrc)=.false.
      return
c======================================================================
  700 format('#   ',3i8/'#  ',a)
c======================================================================
      end subroutine b2blnci
c
      subroutine b2blncc(iline)
      implicit none
      integer, intent(out) :: iline
      real (kind=R8) :: tdummy
      character :: hlp_c
      logical :: ex
c----------------------------------------------------------------------
c*** - to be called by b2co for locating truncation for savefile
c----------------------------------------------------------------------
c
      if(.not.filedata(ntrc_bln)) return
      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)
      call find_file(f,ex)
      open(lun,file=f,err=190)
      rewind(lun)
      iline = 0
      read(lun,'(a1)',err=155,end=155) hlp_c
      do while (hlp_c.eq.'#')
        iline = iline+1
        read(lun,'(a1)',err=155,end=155) hlp_c
      end do
      backspace lun
      read(lun,*,err=155,end=155) tdummy
      do while (tdummy.le.tim)
        iline = iline+1
        read(lun,*,err=155,end=155) tdummy
      end do
  155 continue
      return
c======================================================================
  190 write(*,*) 'b2blnc: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      iline=-1
      return
c======================================================================
      end subroutine b2blncc
c
      subroutine b2blncs
c----------------------------------------------------------------------
c*** - to be called BEFORE each timestep
c----------------------------------------------------------------------
c
      return
      end subroutine b2blncs
c======================================================================
c
      subroutine b2blncf
c----------------------------------------------------------------------
c*** - finalization step
c----------------------------------------------------------------------
c
      return
      end subroutine b2blncf
c======================================================================
      end module b2mod_blnc
