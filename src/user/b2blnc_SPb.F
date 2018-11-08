      subroutine b2blnc_SPb(nx,ny,ns,dtim,iline)
      use b2mod_types
      use b2mod_geo
      use b2mod_time
      use b2mod_diag
      use b2mod_b2cmpa
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
c  version : 22.06.2018 by Ilya Senichenkov
c
c======================================================================
c*** Computes the global particle flows across the boundaries in
c*** B2-Eirene and the global particle balance in B2
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c*** debug_flags usage:
c***  90  >0 : debugging printout on
c***  90, 91 : 'non-default surfaces' (90) to (91)
c***  92, 93 : strata (92) to (93); 0 means the sum over all strata
c***  94, 95 : more detailed output for strata (94) to (95)
c----------------------------------------------------------------------


c----------------------------------------------------------------------
c*** Output data (for each species):
c***
c***    header   array
c***                                        --- Global particle balance
c***   ion_core        flxic       ion flow across the core boundary (based on fna array)
c***   ion_core_mdf    flxic_mdf   ion flow across the core boundary (based on fna_mdf array)
c***   ion_targ        flxit       ion flow to the targets (based on fna array)                                    ! IYS 12.01.2018
c***   ion_wall        flxiw       ion flow to the side walls (based on fna array)                                 ! IYS 12.01.2018
c***   ion_wall        flxi_out    ion flow to the side walls and divertors  (based on fna_mdf array)              ! IYS 12.01.2018
c***   ion_dndt        flxid       time derivative of the ion content in the whole domain
c***   src_ioniz       flxi_out    ion source due to ionization from neutral (ground) state

c***   ntr_core        flxnc       neutral flow across the core boundary
c***   ntr_targ_tot    flxnt       neutral net flow through the targets including recycling source
c***   ntr_targ_rec    flxnt_r     neutral source due to ion recycling on targets       ! IYS 12.01.2018
c***   ntr_targ_pmp    flxnt_pump  neutral-driven net flow through the targets without the source due to ion recycling            ! IYS 12.01.2018
c***   ntr_wall_tot    flxnw       neutral net flow through the side walls including recycling source
c***   ntr_wall_rec    flxnw_r     neutral source due to recycling on side walls       ! IYS 12.01.2018
c***   ntr_wall_pmp    flxnw_pump  neutral-driven net flow through the side walls without the source due to ion recycling         ! IYS 12.01.2018
c***   ntr_puff        flxnp       flow from the gas puff
c***   pump            Npump_tmp   neutral pumping
c***   flux_tot        flxps       total flow coming to the volume
c***   Ftt_inv_time    flxrl       the same in terms of the inverse time scales
c***   flux_balance    flxbl       error in particle balance (should be 0)
c***   Bal_inv_time    flxbr       the same in terms of the inverse time scales
c***


c***                                        --- Ion balance in B2
c***  { Excluded from the output variables list by IYS 12.01.2018
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
!      IYS 12.01.2018 }

c----------------------------------------------------------------------

      implicit none
c
      integer nx, ny, ns, iline
      real (kind=R8) :: dtim
c
      integer lun, ltrc
      parameter (ltrc=ntrc_bln_SPb)

      integer i, j, k, l, m, n, ic, ii, is, ix, iy, js, ks, ls, lv,
     ,        nc, nv, kks, np_tmp

      logical first, ex !, hlp_lw, hlp_lt
      real (kind=R8) :: tdummy
      integer, save :: mdf_fnb = 0, style_b2mndt = 0                ! IYS 12.01.2018

      real (kind=R8), allocatable, save :: hlpit_lp(:),hlpit_rp(:),
     &                                                 hlpit_pl(:)
      real (kind=R8) :: hlp_pls(nbcd,ns)
      integer use_recyceir
      real (kind=R8) :: recyctmp
      real (kind=R8), save :: stim = 0.0_R8

c>>>
c
c*** lv    :  version number to appear in the file
c*** np_tmp:  number of output data for each species
c
      parameter (lv=20180622, np_tmp=16)                            ! IYS 22.06.2018
c
      real (kind=R8) :: hlpdt, hlp 

      real (kind=R8), allocatable, save :: 
     *     flxic(:),flxnc(:),flxit(:),                                            ! IYS 12.01.2018
     ,     flxnt(:),flxiw(:),flxnw(:),
     ,     flxnp(:),flxps(:),flxnw_rec(:),flxnw_pump(:),flxnw_sput(:),            ! IYS 12.01.2018
     ,     flxnt_rec(:),flxnt_pump(:),flxnt_sput(:),                              ! IYS 12.01.2018
c*** from 960727 on
     ,     flxid(:),!flxbl(:),flxbr(:),
c*** from 960808 on
!     ,     flxrl(:),
c*** from 20180604 on
     ,     flxis(:), Nions(:)

c----------------------------------------------------------------------
c
      character f*(lenfnm), hlp_c, hlp_h(np_tmp)*12, hlp_frm*256, tab
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
c*** Column headers:
      data hlp_h /
     1    '    ion_core','    ntr_core','    ion_targ','ntr_targ_tot',
     2    'ntr_targ_rec','ntr_targ_spt','ntr_targ_pmp','    ion_wall',
     3    'ntr_wall_tot','ntr_wall_rec','ntr_wall_spt','ntr_wall_pmp',
     4    '   src_ioniz','    ntr_puff','    flux_tot','   ion_dn/dt'/
!     5    'Ftt_inv_time',,'flux_balance','Bal_inv_time'/

*   ..initialisation
      data use_recyceir /1/
c
      data first / .true. /
c----------------------------------------------------------------------
      save lun, i, j, k, l, m, n, ic, is, ix, iy,
     1 js, ks, ls, nc, nv, kks, use_recyceir,
     2 first, ex
c======================================================================

      if(first_time_step) return

      if(ank_tracing.eq.0) return

      if(first) then

        call ipgeti('b2tfnb_mdf',mdf_fnb)
        call ipgeti('b2mndt_style',style_b2mndt)

#ifdef B25_EIRENE
c<<<
        call ipgeti ('eirene_use_recyceir', use_recyceir)
        print *,'=== Regions of the grid:'
        print *
        l=ny+1
        n=mod(l-1,10)+1
        do k=0,ny-1
          m=l
          l=max(0,m-n+1)
          n=10
          do j=m-1,l-1,-1
            print '(1x,i4,1x,i1,40(1x,10i1))',j,
     ,                                         (region(i,j,0),i=-1,nx)
          end do
          l=l-1
        end do
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

csw 2012 correct nlimps --> nnstsi+nnlimi (known after 1st eirene call)
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
      debug_flags(90) = 1
      debug_flags(91) = 6
      debug_flags(92) = 1
      debug_flags(93) = nstrat
 
      if(debug_flags(90).gt.0) then
        print *
        print *,'%%% Eirene data in B2 %%%'
        do k=debug_flags(92),debug_flags(93)
          print *
          print *,'wldpp:',k
          do j=nnlimi+debug_flags(90),nnlimi+debug_flags(91)
            print '(1p,i4,20e9.2)',j,(wldpp(j,i,k),i=1,nnplsi)
          end do
          print *
          print *,'wldpeb, wldpa, wldpm, hlpit_pl:',k
          do j=nnlimi+debug_flags(90),nnlimi+debug_flags(91)
            ks=0
            do ls=1,nspecies
              hlpit_pl(ls)=0.0_R8
              do is=1,nfluids(ls)
                ks=ks+1
                hlpit_pl(ls)=hlpit_pl(ls)+wldpp(j,ks,k)
              end do
            end do
            print '(1p,i4,20e12.4)',j,wldpeb(j,k),
     ,         (wldpa(j,i,k),i=1,nnatmi),(wldpm(j,i,k),i=1,nnmoli),
     ,                                       (hlpit_pl(i),i=1,nspecies)
          end do
        end do
      end if

#else
        first=.false.
      end if
#endif
c>>>

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! {IYS 04.06.2018
! The logics is changed completely

! Set all output arrays to zero

      flxic        = 0.0_R8                 ! Ion flux through core boundary based on fna array (sum over all core boundary segments and over all ionization states of given isonuclear sequence)

      flxit        = 0.0_R8                 ! Ion flux through targets (sum over all targets and over all ionization states of given isonuclear sequence)

      flxiw        = 0.0_R8                 ! Ion flux through wall boundaries based on fna array 

      flxid        = 0.0_R8                 ! Time derivative of the ion content in the whole domain (sum over all cells and over all ionization states of given isonuclear sequence)

      flxis        = 0.0_R8                 ! Ion source due to ionization from the neutral state (sum over all cells for first ionized state of given isonuclear sequence)

      flxnc        = 0.0_R8                 ! Neutral flux through core boundary (sum over all core boundary segments for neutrals corresponding to given isonuclear sequence)
      flxnt        = 0.0_R8                 ! Neutral flux through targets (sum over all targets and over all kind of molecules containing given atom)
      flxnt_rec    = 0.0_R8                 ! Neutral flux emitted from the target due to plasma ion recycling (part of the source from b2stbr for fluid neutrals),  IYS 12.01.2018, 04.06.2018
      flxnt_pump   = 0.0_R8                 ! Pumped/absorbed neutral flux at targets (source from b2stbc for fluid neutrals),  IYS 12.01.2018, 04.06.2018   
      flxnt_sput   = 0.0_R8                 ! Sputtered neutral flux at targets (part of the source from b2stbr for fluid neutrals),  IYS 22.06.2018 
      flxnw        = 0.0_R8                 ! Neutral flux through wall boundary (sum over all side wall segments for neutrals corresponding to given isonuclear sequence)
      flxnw_rec    = 0.0_R8                 ! Neutral flux emitted from the side wall due to plasma ion recycling (part of the source from b2stbr for fluid neutrals),  IYS 12.01.2018, 04.06.2018
      flxnw_pump   = 0.0_R8                 ! Pumped/absorbed neutral flux at wall (source from b2stbc for fluid neutrals),  IYS 12.01.2018, 04.06.2018
      flxnw_sput   = 0.0_R8                 ! Sputtered neutral flux at wall (part of the source from b2stbr for fluid neutrals),  IYS 22.06.2018, 04.06.2018
      flxnp        = 0.0_R8                 ! Neutral puff (for fluid neutrals -- the flux prescribed by BCCON=5,8 on walls)

      flxps        = 0.0_R8                 ! Sum of ion and neutral fluxes into the domain with sign (positive: inside the domain; negative: outside the domain)

c**********************************
! The flows of B2 fluids
c**********************************

!     Loop over isonuclear sequence
      do ls = 1,nspecies

c<<<
        if(debug_flags(90).gt.0) then
          print *,'%%# Species # ',ls
          print *,'%%# eb2spcr(ls) = ',eb2spcr(ls)
          print *,'%%# nfluids(ls) = ',nfluids(ls)
        end if
c>>>

!       Loop over species within given isonuclear sequence
c***        EB2SPCR(ls) points to the first B2 fluid (i.e. the neutral fluid) corresponding to the isonuclear sequence labeled by ls
c***        NFLUIDS(ls) is the number of charges states (== B2 plasma fluids) corresponding to the isonuclear sequence labeled by ls 
        do is = EB2SPCR(ls), EB2SPCR(ls) + NFLUIDS(ls)

!         Loop over B2 boundaries
          do ic = 1,nbc

!           Loop over poloidal coordinate along this boundary
            do ix = 1, BC_LIST_SIZE(ic)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             Check whether it is a southern boundary
              if (BCCHAR(ic).EQ.'S') then

!             Check whether it is a closed surface, i.e. is the part of the core boundary
                if (on_closed_surface(BC_LIST_X(1,ic), 
     .                                        BC_LIST_Y(1,ic))) then
#ifndef B25_EIRENE
! Count fluid neutrals if the Eirene is not coupled
                  if (is_neutral(is)) then
                     flxnc(ls) = flxnc(ls) + 
     &     fna(topix(BC_LIST_X(ix,ic),BC_LIST_Y(ix,ic)),
     &               topiy(BC_LIST_X(ix,ic),BC_LIST_Y(ix,ic)),1,is)
                  else  ! not neutral
#endif
                    if (mdf_fnb .eq. 0 .and. style_b2mndt .eq. 0) then
                      flxic(ls) = flxic(ls) +
     .     fna(topix(bc_list_x(ix,ic),bc_list_y(ix,ic)),
     .                 topiy(bc_list_x(ix,ic),bc_list_y(ix,ic)),1,is)
                    else 
                      flxic(ls) = flxic(ls) +
     .     fna_mdf(topix(bc_list_x(ix,ic),bc_list_y(ix,ic)),
     .                  topiy(bc_list_x(ix,ic),bc_list_y(ix,ic)),1,is)
                    endif

                    if(debug_flags(90).gt.0) then
                       hlp_pls(ic,ls)=hlp_pls(ic,ls)+
     .      fna(topix(bc_list_x(ix,ic),bc_list_y(ix,ic)),
     .                  topiy(bc_list_x(ix,ic),bc_list_y(ix,ic)),1,is)
                    endif
#ifndef B25_EIRENE
                  endif    ! if (is_neutral(is)) ...
#endif
                else        ! .not. (on_closed_surface(:,:)), i.e. it is the Privatte Flux Region

#ifndef B25_EIRENE
                  if (is_neutral(is)) then
                    flxnw(ls) = flxnw(ls) +
     .      fna(topix(bc_list_x(ix,ic),bc_list_y(ix,ic)),
     .                   topiy(bc_list_x(ix,ic),bc_list_y(ix,ic)),1,is)

                    if (BCCON(is,ic) .eq. 5) flxnp(ls) = flxnp(ls) + 
     &   CONPAR(IS,IC,1) * gs(topix(BC_LIST_X(IX,IC),BC_LIST_Y(IX,IC)), 
     &              topiy(BC_LIST_X(IX,IC),BC_LIST_Y(IX,IC)),1)              

                  else
#endif

                    flxiw(ls) = flxiw(ls) +
     .      fna(topix(bc_list_x(ix,ic),bc_list_y(ix,ic)),
     .                   topiy(bc_list_x(ix,ic),bc_list_y(ix,ic)),1,is)
c<<<
                    if(debug_flags(90).gt.0) then
                       hlp_pls(ic,ls)=hlp_pls(ic,ls) +
     .      fna(topix(bc_list_x(ix,ic),bc_list_y(ix,ic)),
     .                   topiy(bc_list_x(ix,ic),bc_list_y(ix,ic)),1,is)
                    endif
#ifndef B25_EIRENE
                  endif
#endif
              
                endif      ! if (on_closed_surface()) ... 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             Northern boundary
!             Do the same as for PRF boundary, but add fluxes with opposite sign due to opposite surface normal orientation
              elseif (BCCHAR(ic).EQ.'N') then

#ifndef B25_EIRENE
                  if (is_neutral(is)) then
                    flxnw(ls) = flxnw(ls) -
     .      fna(bc_list_x(ix,ic),bc_list_y(ix,ic),1,is)

                    if (BCCON(is,ic) .eq. 5) flxnp(ls) = flxnp(ls) + 
     &      CONPAR(IS,IC,1) * gs(BC_LIST_X(IX,IC),BC_LIST_Y(IX,IC),1)               

                  else
#endif

                    flxiw(ls) = flxiw(ls) -
     .      fna(bc_list_x(ix,ic),bc_list_y(ix,ic),1,is)

                    if(debug_flags(90).gt.0) then
                      hlp_pls(ic,ls)=hlp_pls(ic,ls) -
     .      fna(bc_list_x(ix,ic),bc_list_y(ix,ic),1,is)
                    endif
#ifndef B25_EIRENE
                  endif
#endif
           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             Western boundary (left target(s))
              elseif (BCCHAR(ic).EQ.'W') then

#ifndef B25_EIRENE
                if (is_neutral(is)) then
                  flxnt(ls) = flxnt(ls) +
     .      fna(rightix(bc_list_x(ix,ic),bc_list_y(ix,ic)),
     .                 rightiy(bc_list_x(ix,ic),bc_list_y(ix,ic)),0,is)

                  if (BCCON(is,ic) .eq. 5) flxnp(ls) = flxnp(ls) + 
     & CONPAR(IS,IC,1) * gs(rightix(bc_list_x(ix,ic),bc_list_y(ix,ic)),
     &                    rightiy(bc_list_x(ix,ic),bc_list_y(ix,ic)),0)               

                else
#endif
                  flxit(ls) = flxit(ls) +
     .      fna(rightix(bc_list_x(ix,ic),bc_list_y(ix,ic)),
     .                 rightiy(bc_list_x(ix,ic),bc_list_y(ix,ic)),0,is)

                  if(debug_flags(90).gt.0) then
                  
                    hlpit_lp(ls) = hlpit_lp(ls) +
     .               fna(rightix(bc_list_x(ix,ic),bc_list_y(ix,ic)),
     .                 rightiy(bc_list_x(ix,ic),bc_list_y(ix,ic)),0,is)

                    hlp_pls(ic,ls) = hlp_pls(ic,ls) +
     .               fna(rightix(bc_list_x(ix,ic),bc_list_y(ix,ic)),
     .                 rightiy(bc_list_x(ix,ic),bc_list_y(ix,ic)),0,is)
                  endif
#ifndef B25_EIRENE
                endif
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             Eastern boundary (right target(s))
!             Do the same as for left target, but add fluxes with opposite sign due to opposite surface normal orientation
              elseif (BCCHAR(ic).EQ.'E') then

#ifndef B25_EIRENE
                if (is_neutral(is)) then
                  flxnt(ls) = flxnt(ls) -
     .      fna(bc_list_x(ix,ic),bc_list_y(ix,ic),0,is)

                  if (BCCON(is,ic) .eq. 5) flxnp(ls) = flxnp(ls) + 
     &  CONPAR(is,ic,1) * gs(bc_list_x(ix,ic),bc_list_y(ix,ic),0)                

                else
#endif
                  flxit(ls) = flxit(ls) -
     .      fna(bc_list_x(ix,ic),bc_list_y(ix,ic),0,is)

                  if(debug_flags(90).gt.0) then
                  
                    hlpit_lp(ls) = hlpit_lp(ls) -
     .      fna(bc_list_x(ix,ic),bc_list_y(ix,ic),0,is)

                    hlp_pls(ic,ls) = hlp_pls(ic,ls) -
     .      fna(bc_list_x(ix,ic),bc_list_y(ix,ic),0,is)
                  endif
#ifndef B25_EIRENE
                endif
#endif
              endif     ! if (BCCHAR(ic).EQ.'S')

!           End of the loop over poloidal coordinate along the boundary
            enddo

#ifndef B25_EIRENE
            if (is_neutral(is) ) then 

! Compute the puffing rate
              if (BCCON(is,ic).eq.8 .or. 
     &            BCCON(is,ic).eq.16 .or. BCCON(is,ic).eq.26) then
                flxnp(ls) = flxnp(ls) + CONPAR(is,ic,1)
              endif

! Compute neutral flux due ion recycling
! and neutral pumped flux
              if (BCCHAR(ic).EQ.'W' .or. BCCHAR(ic).EQ.'E') then

!FIXME
                  flxnt_rec(ls) = flxnt_rec(ls) + sna_stbr_neut(ls,ic)

                  flxnt_pump(ls) = flxnt_pump(ls) + 
     &                                    sna_stbc_neut(ls,ic)

                  flxnt_sput(ls) = flxnt_sput(ls) + 
     &                                    sna_stbr_sput(ls,ic)

              elseif ( BCCHAR(ic).EQ.'N' .or. (BCCHAR(ic).EQ.'S' .and.
     & .not. on_closed_surface(BC_LIST_X(1,ic),BC_LIST_Y(1,ic))) ) then

!FIXME
                    flxnw_rec(ls) = flxnw_rec(ls) + sna_stbr_neut(ls,ic)

                    flxnw_sput(ls) = flxnw_sput(ls) + 
     &                                    sna_stbr_sput(ls,ic)

!             Exclude surfaces with puffing/seeding (for fluid neutrals)
                    if (BCCON(is,ic).ne.5 .and. BCCON(is,ic).ne.8
     &  .and. BCCON(is,ic).ne.16 .and. BCCON(is,ic).ne.26) then

                      flxnw_pump(ls) = flxnw_pump(ls) + 
     &                                    sna_stbc_neut(ls,ic)
                    endif

              endif
            endif
#endif

!         End of the loop over B2 boundaries
          enddo

!!!!!
! Compute the sum of ionization sources for all ionized states of given isonuclear sequence and the time derivative term
!!!!!
          do iy=0,ny-1
            do ix=0,nx-1

              j=region(ix,iy,0)

              if (j.ne.0 .and. leftix(ix,iy).ne.-2 
     &                    .and. rightix(ix,iy).ne.nx+1) then

                if (solveco(is,j)) then
                  hlpdt = dtim * dtco(is,j) * time_factor(ix,iy)     ! hlpdt is temporary variable - BE CAREFUL!!!!
                else
                  hlpdt = 0.0_R8
                endif
                if (hlpdt .gt. 0.0_R8 .and. .not. is_neutral(is)) then
                  flxid(ls) = flxid(ls) + 
     &        (na(ix,iy,is) - na0(ix,iy,is)) * vol(ix,iy) / hlpdt     ! hlpdt is temporary variable - BE CAREFUL!!!!
                endif

                if (is .eq. EB2SPCR(ls)+1) then  ! means the first ionized state
#ifdef B25_EIRENE
!     If Eirene is coupled, just take the ionization source from Eirene for given isonuclear sequence
                  flxis(ls) = sna_eir(ls)
#else
!     For standalone run take the source of neutrals with the opposite sign
                  flxis(ls) = flxis(ls) - sna(ix,iy,0,is-1) -
     1             sna(ix,iy,1,is-1) * na(ix,iy,is-1) 

!     But be careful: the source due to time derivative should be counted
                  if (hlpdt .gt. 0.0_R8) then
                    flxis(ls) = flxis(ls) - 
     &       (na(ix,iy,is-1) - na0(ix,iy,is-1)) * vol(ix,iy) / hlpdt     ! hlpdt is temporary variable - BE CAREFUL!!!!
                  endif
#endif
                endif


              endif ! if (j.ne.0 .and. leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)

            enddo
          enddo

!       End of the loop over species within given isonuclear sequence (is)
        enddo

c**********************************
! The flows of B2 fluids are now computed
c**********************************

c****************************************************************
! Now compute the neutral fluxes from Eirene, if Eirene is coupled

#ifdef B25_EIRENE

c*** Flux across the core boundary

!       Loof over segments of Eirene core boundary
        do is = 1,ncns

c**       Here is is not a species index, but the index of eirene surface,
!!!    because
c**       ncns - Number of 'non-default' eirene surfaces which are part of the core boundary
c**       lcns() - index of these 'non-default' surfaces

          js = lcns(is) + nnlimi
c**       js - index of the eirene surface, js > nnlimi means that this is the 'non-default' one

          do j = 1,nnatmi
            if (latmscl(j) .eq. ls) then
              flxnc(ls) = flxnc(ls) 
     &               - wldna(js,j,0) + wldra(js,j,0) + wldpa(js,j,0)
              do ic=1,nnmoli
                if (mlcmp(j,ic) .gt. 0) flxnc(ls) = flxnc(ls) - 
     & mlcmp(j,ic) * (wldnm(js,ic,0) - wldrm(js,ic,0) - wldpm(js,ic,0))
              end do
            end if
          end do
        end do

c*** Flux across the targets

!       Loop over segments of Eirene target boundaries
        do is = 1,ntns

c**       ntns - Number of 'non-default' eirene surfaces which are part of the target boundaries
c**       ltns() - index of these 'non-default' surfaces

          js=ltns(is)+nnlimi
c**       js - index of the eirene surface, js > nnlimi means that this is the 'non-default' one

          do j = 1,nnatmi

            if (latmscl(j) .eq. ls) then

              flxnt(ls) = flxnt(ls) - 
     &             wldna(js,j,0) + wldra(js,j,0) + wldpa(js,j,0)
              flxnt_rec(ls) = flxnt_rec(ls) + wldpa(js,j,0)                        ! IYS 12.01.2018
              flxnt_pump(ls) = flxnt_pump(ls) + wlpump(j,js)                       ! IYS 22.06.2018
              flxnt_sput(ls) = flxnt_sput(ls) + wldspta(js,j,0)                    ! IYS 22.06.2018

              do ic = 1,nnmoli

                if (mlcmp(j,ic) .gt. 0) then
 
                  flxnt(ls) = flxnt(ls) - mlcmp(j,ic) *
     &               (wldnm(js,ic,0) - wldrm(js,ic,0) - wldpm(js,ic,0))
                  flxnt_rec(ls) = flxnt_rec(ls) + 
     &                            mlcmp(j,ic) * wldpm(js,ic,0)
                  flxnt_pump(ls) = flxnt_pump(ls) + 
     &                            mlcmp(j,ic) * wlpump(nnatmi+ic,js)  ! IYS 22.06.2018
                  flxnt_sput(ls) = flxnt_sput(ls) + 
     &                            mlcmp(j,ic) * wldsptm(js,ic,0)      ! IYS 22.06.2018

                endif

              end do ! end of the loop over molecules

            end if

          end do ! end of the loop over atoms

        end do ! end of the loop over target segments

c*** Flux across the side walls

!       Loop over segments of the side wall
        do js = 1,nnlimi

          do j = 1,nnatmi

            if (latmscl(j) .eq. ls) then

              flxnw(ls) = flxnw(ls) - 
     &                 wldna(js,j,0) + wldra(js,j,0) + wldpa(js,j,0)
              flxnw_rec(ls) = flxnw_rec(ls) + wldpa(js,j,0)           ! IYS 12.01.2018, wldpa(js,j,0) should be zero, since no ion reaches side wall
              flxnw_pump(ls) = flxnw_pump(ls) + wlpump(j,js)                       ! IYS 22.06.2018
              flxnw_sput(ls) = flxnw_sput(ls) + wldspta(js,j,0)                    ! IYS 22.06.2018

              do ic = 1,nnmoli

                if (mlcmp(j,ic) .gt. 0) then

                  flxnw(ls) = flxnw(ls) - mlcmp(j,ic) *
     .              (wldnm(js,ic,0) - wldrm(js,ic,0) - wldpm(js,ic,0))        ! wldpm(js,j,0) should be zero, since no ion reaches side wall
                  flxnw_rec(ls) = flxnw_rec(ls) + 
     &                         mlcmp(j,ic) * wldpm(js,ic,0)
                  flxnw_pump(ls) = flxnw_pump(ls) +
     &                 mlcmp(j,ic) * wlpump(nnatmi+ic,js)  ! IYS 22.06.2018
                  flxnw_sput(ls) = flxnw_sput(ls) +
     &                 mlcmp(j,ic) * wldsptm(js,ic,0)      ! IYS 22.06.2018

                endif

              end do ! end of the loop over molecules

            end if

          end do ! end of the loop over atoms

        end do ! end of the loop over wall segments

c*** Puffing/seeding rate and recycling source

!       Loop over strata

        do ic = 1,nstrat

         if(use_recyceir.gt.0) then
              recyctmp=recyceir(ic)
          else
              recyctmp=recyc(kks,ic)
          endif

          if (crcstra(ic) .eq. 'C') then  ! Puffing or seeding

            do j = 1,nnatmi
              if (latmscl(j) .eq. ls) flxnp(ls) = flxnp(ls) +
     &                 userfluxparm(ic,1) * gpfc(j,ic) * recyctmp
            end do

          elseif (crcstra(ic) .eq. 'S' .or. crcstra(ic) .eq. 'N') then ! Recycling on the B2 outer boundaries except targets

            hlp = 0.0_R8

            do i = msns(1,ic), msns(1,ic) + msns(2,ic) - 1
 
               k = lsns(i)
               if (k .lt. 0) k=nnlimi-k

               if (k .gt. 0) then

                 do j=1,nnatmi                   
                   if (latmscl(j).eq.ls) then
                      hlp = hlp + wldpa(k,j,ic) * sclstra(j,ic)  ! BE CAREFUL: hlp is a temporary variable
                   end if
                 end do

                 do j=1,nnmoli
                   do ii=1,nnatmi
                     if (latmscl(ii).eq.ls) then
                       hlp = hlp + 
     &                     mlcmp(ii,j) * wldpm(k,j,ic) * sclstra(ii,ic)  ! BE CAREFUL: hlp is a temporary variable
                     end if
                   enddo
                 enddo

                endif

              enddo
          
              flxnw_rec(ls) = flxnw_rec(ls) + hlp * recyctmp

            endif

!       End of the loop over strata
        enddo

c**********************************
! The neutral fluxes from Eirene are now computed
c**********************************
#endif

c**** Total flux coming to the volume =
c       Ion and neutral flux through the core +
c       Ion and neutral flux through targets +
c       Ion and neutral flux through side walls

c     For fluid neutrals "side walls" means northern boundary plus PRF southerm boundary.
c     Puffing (BCCONs=5,8,16,26) and leackage (BCCONs=10,15 etc) fluxes are counted automatically,
c     as well as recycling and sputtering?

c     For coupled runs "side walls" means the real vessel boundary.
c     Then the pumping is counted automatically, but puffing (which is located
c     somewhere in between the vessel walls and B2 boundary) should be added manually.
c     Also, since no ion reaches the vessel boundary (except targets), the neutral source
c     due to ion recycling on B2 northern boundary and B2 PRF southerm boundary should
c     be added to a global balance as a wall recycling source

! FIXME

        flxps(ls) =  flxic(ls) + flxnc(ls) + flxit(ls) + flxnt(ls) +
     &          flxiw(ls) + flxnw(ls) 
#ifdef B25_EIRENE
     &          + flxnp(ls) + flxnw_rec(ls)
#endif

!     End of the loop over isonuclear sequence (ls)
      enddo


c********* Debug output **************
      if(debug_flags(90).gt.0) then
c        print '(a,1p,10e12.4)','%%% hlpit_lp, hlpit_rp: ',
c     ,                           (hlpit_lp(i),hlpit_rp(i),i=1,nspecies)
        print '(a,1p,10e12.4)','%%# hlpiw_s : ',
     ,                                        (flxiw(i),i=1,nspecies)
      end if
c********* Debug output **************
c

c*** Write the data in a separate file for plotting...
      if(filedata(ltrc)) then
c
c
        if(nc.gt.np_tmp .or. nc.le.0) then
          write(*,*) 'b2blnc_SPb: unexpected error. nc =',nc
          filedata(ltrc)=.false.
        else
          write(lun,hlp_frm) tim,
     .      (flxic(i), flxnc(i), flxit(i),                                 ! IYS 12.01.2018
     ,       flxnt(i), flxnt_rec(i), flxnt_sput(i), flxnt_pump(i),         ! IYS 12.01.2018
     ,       flxiw(i), flxnw(i), flxnw_rec(i), flxnw_sput(i),              ! IYS 12.01.2018
     ,       flxnw_pump(i),                                                ! IYS 12.01.2018
     ,       flxis(i), flxnp(i), flxps(i),
     ,       flxid(i),
     ,       i=1,nspecies)
        end if

      else
c
c*** ... or just print out the data
c
        write(*,*) 'particle_fluxes ',nspecies
        write(*,'(30(3x,a8,1x))') (hlp_h(i),i=1,np_tmp)
        do i=1,nspecies
          write(*,'(1p,30e12.3)') 
     ,     flxic(i), flxnc(i), flxit(i),                                   ! IYS 12.01.2018
     ,       flxnt(i), flxnt_rec(i), flxnt_sput(i), flxnt_pump(i),         ! IYS 12.01.2018
     ,       flxiw(i), flxnw(i), flxnw_rec(i), flxnw_sput(i),              ! IYS 12.01.2018
     ,       flxnw_pump(i),                                                ! IYS 12.01.2018
     ,       flxis(i),flxnp(i),  flxps(i),
     ,       flxid(i)
        end do
      end if
      return

c======================================================================
c
      entry b2blnci_SPb(nx,ny,ns,dtim,iline)
c----------------------------------------------------------------------
c*** Initialise the file for tracing data
c*** Only called with filedata(ntrc_bln_SPb)=.true.
c----------------------------------------------------------------------
c
      if (.not. allocated(hlpit_lp))
     &   allocate(hlpit_lp(nspecies), hlpit_rp(nspecies), 
     &            hlpit_pl(nspecies))

      if (.not. allocated(flxic)) 
     &      allocate(
     &   flxic(nspecies), flxnc(nspecies), flxit(nspecies),
     &   flxnt(nspecies), flxiw(nspecies), flxnw(nspecies),
     &   flxnp(nspecies), flxps(nspecies), 
     & flxnw_rec(nspecies), flxnw_pump(nspecies), flxnw_sput(nspecies),                                      ! IYS 12.01.2018
     & flxnt_rec(nspecies), flxnt_pump(nspecies), flxnt_sput(nspecies),                                       ! IYS 12.01.2018
c*** from 960727 on 
     &   flxid(nspecies), !flxbl(nspecies), flxbr(nspecies),
c*** from 960808 on
!     &   flxrl(nspecies),
c*** from 20180604 on
     &   flxis(nspecies), Nions(nspecies))


      call ipgetr ('b2mndr_stim', stim)

      if(.not.filedata(ntrc_bln_SPb)) return

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
     ,     '!!! b2blnc_SPb: current version ',lv,
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
     ,     '!!! b2blnc_SPb: current number of columns ',np_tmp,
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
     ,     '!!! b2blnc_SPb: current number of species ',nspecies,
     ,      ' differs from that in the file ',trim(f),' : ',j
          write(*,*) 'No ',trim(f),' tracing file - ',
     ,                              'remove the existing file first!'
          filedata(ltrc)=.false.
          return
        end if
   50   read(lun,'(a1)',end=55) hlp_c
        go to 50
   55   backspace lun

      else ! file 'blnc.trc' does not exist

        nc = np_tmp

        write(lun,700) lv, nspecies, np_tmp, label
        write(hlp_frm,
     &      '(18h(1h#,7x,4htime,4x,i3,17h(1h'//tab//',a16,1h_,a2)))')
     &                                                  np_tmp*nspecies

c       print *,'b2blnc_SPb: nspecies= ',nspecies                   !###
c       print *,'species_list :',(species_list(j),':',j=1,nspecies) !###

        write(lun,hlp_frm) ((hlp_h(i),species_list(j),i=1,np_tmp),
     &                                                    j=1,nspecies)

      end if

      write(hlp_frm,'(10h(1p,e16.7,,i3,12h(1h'//tab//',e16.7)))')
     &                                                  np_tmp*nspecies
      return
c======================================================================
   90 write(*,*) 'b2blnc_SPb: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      return
c
   92 write(*,*) 'b2blnc_SPb: wrong format of the file ',trim(f)
      filedata(ltrc)=.false.
      return
c======================================================================
c
      entry b2blncc_SPb(nx,ny,ns,dtim,iline)
c----------------------------------------------------------------------
c*** - to be called by b2co for locating truncation for savefile
c----------------------------------------------------------------------
c
      if(.not.filedata(ntrc_bln_SPb)) return
      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)
      call find_file(f,ex)
      open(lun,file=f,err=190)
      rewind(lun)
      iline = 0
      read(lun,'(a1)',end=155) hlp_c
      do while (hlp_c.eq.'#')
        iline = iline+1
        read(lun,'(a1)',end=155) hlp_c
      end do
      backspace lun
      read(lun,*,end=155) tdummy
      do while (tdummy.le.tim)
        iline = iline+1
        read(lun,*,end=155) tdummy
      end do
  155 continue
      return
c======================================================================
  190 write(*,*) 'b2blnc_SPb: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      iline=-1
      return
c======================================================================
c
      entry b2blncs_SPb(nx,ny,ns,dtim,iline)
c----------------------------------------------------------------------
c*** - to be called BEFORE each timestep
c----------------------------------------------------------------------
c
      return
c======================================================================
c
      entry b2blncf_SPb(nx,ny,ns,dtim)
c----------------------------------------------------------------------
c*** - finalization step
c----------------------------------------------------------------------
c
      if (allocated(flxic)) 
     &  deallocate(
     &    flxic, flxnc, flxit,
     &    flxnt, flxiw, flxnw,
     &    flxnp, flxps, !blnna_tmp,
     &    flxnw_rec, flxnw_pump, flxnw_sput,                                       ! IYS 12.01.2018
     &    flxnt_rec, flxnt_pump, flxnt_sput,                                       ! IYS 12.01.2018
c*** from 960727 on 
     &    flxid, !flxbl, flxbr,
c*** from 960808 on
!     &    flxrl,
c*** from 20180604 on
     &    flxis, Nions)

      if (allocated(hlpit_pl)) 
     &  deallocate(hlpit_pl, hlpit_lp, hlpit_rp)


      return

c======================================================================
  700 format('#   ',3i8/'#  ',a)
  710 format(a1,3x,3i8/3x,a)
c======================================================================
      end subroutine b2blnc_SPb
