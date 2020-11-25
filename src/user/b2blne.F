      subroutine b2blne(nx,ny,ns,BoRiS,iline)
      use b2mod_types
      use b2mod_geo
      use b2mod_time
      use b2mod_diag
      use b2mod_rates
      use b2mod_b2cmpa
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_constants
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_numerics_namelist
      use b2mod_external
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif

c
c  version : 29.07.2000 23:52
c
c=======================================================================
c*** Computes the global energy flows across the boundaries and the
c*** sources in B2-Eirene, and the global energy balance in B2
c-----------------------------------------------------------------------
      implicit none
      integer nx,ny,ns,iline
      real (kind=R8) :: BoRiS
c
c*** lv    :  version number to appear in the file
c
      integer lv
      parameter (lv=20000704)
c
c-----------------------------------------------------------------------
c*** Power balance:
c***
c***      header  var       Legend
c***                                                            (global)
c***     pwr_totl pwtot   sum of the inputs and losses
c***                      (should be zero in steady state)
c***     pwr_tscl pwttr   the same divided by the total energy content
c***     tot_core pwcrt   total power input from the core
c***     tot_targ pwtgt   total power flow to the targets
c***     tot_wall pwwlt   total power flow to the side walls
c***     tot_rad  pwrdt   total radiated power
c***     heat_cor pwcrh   heat from the core
c***     ptnt_cor pwcrp   potential energy from the core
c***     ntrl_cor pwcrn   neutral power loss to the core
c***     plsm_trg pwtgp   plasma power to the targets
c***     ntrl_trg pwtgn   neutral power to the targets
c***     plsm_wll pwwlp   plasma to the side walls
c***     ntrl_wll pwwln   neutrals to the side walls
c***     brms_rad radbr   Bremsstrahlung losses
c***     ntrl_rad radnu   power radiated by neutrals
c***     tot_enrg enrtt   total plasma energy content
c***     tot_enre enrte   electron energy content
c***     tot_enri enrti   ion energy content
c***                                                          (species)
c***     imp_rad_ radim   power radiated by impurity ions
c***                                                          (targets)
c***     pwr_totl pwttt   total power flow to the target
c***     pwr_plsm pwtpl   plasma to the target
c***     pwr_neut pwtnu   neutrals to the target
c***     pwr_ionz pwtio   ionisation energy to the target
c***     pwr_diss pwtds   dissociation energy to the target
c***   suffices: il, ol, iu, ou, wl for inner and outer targets and wall
c-----------------------------------------------------------------------
c*** Global values
      integer nglb_tmp, ng_tmp
      parameter (nglb_tmp=19)
      real (kind=R8) ::
     *     equg_tmp(nglb_tmp)
      real (kind=R8) ::
     *     pwtot_tmp, pwttr_tmp, pwcrt_tmp, pwtgt_tmp, pwwlt_tmp,
     ,     pwrdt_tmp, pwcrh_tmp, pwcrp_tmp, pwcrn_tmp, pwtgp_tmp,
     ,     pwtgn_tmp, pwwlp_tmp, pwwln_tmp, radbr_tmp, radnu_tmp,
     ,     rdtim_tmp, enrtt_tmp, enrti_tmp, enrte_tmp
      equivalence
     ,     (equg_tmp(1),  pwtot_tmp), (equg_tmp(2),  pwttr_tmp),
     ,     (equg_tmp(3),  pwcrt_tmp), (equg_tmp(4),  pwtgt_tmp),
     ,     (equg_tmp(5),  pwwlt_tmp), (equg_tmp(6),  pwrdt_tmp),
     ,     (equg_tmp(7),  pwcrh_tmp), (equg_tmp(8),  pwcrp_tmp),
     ,     (equg_tmp(9),  pwcrn_tmp), (equg_tmp(10), pwtgp_tmp),
     ,     (equg_tmp(11), pwtgn_tmp), (equg_tmp(12), pwwlp_tmp),
     ,     (equg_tmp(13), pwwln_tmp), (equg_tmp(14), radbr_tmp),
     ,     (equg_tmp(15), radnu_tmp), (equg_tmp(16), rdtim_tmp),
     ,     (equg_tmp(17), enrtt_tmp), (equg_tmp(18), enrti_tmp),
     ,     (equg_tmp(19), enrte_tmp)
      character hlp_hg(nglb_tmp)*8
c-----------------------------------------------------------------------
c*** Values by species
      integer nspc_tmp, nimp_tmp, ni_tmp, ns_tmp
      integer, allocatable :: limp_tmp(:)
      parameter (nspc_tmp=1)
      real (kind=R8), allocatable :: radim_tmp(:)
      character*2, allocatable :: hlp_spc(:)
      character*8 hlp_hs(nspc_tmp)
      real (kind=R8), save :: stim = 0.0_R8
      real (kind=R8) :: tdummy
c-----------------------------------------------------------------------
c*** Values by surface. Here ntrg_tmp is the number of values,
c*** and ntrgs_tmp is the maximum number of different surfaces.
c*** At present, there can be up to 5 surfaces: 4 targets and side wall
c*** (to be generalised later)
c*** ntrgf_tmp is the actual number of surfaces for the geometry used
c*** nz_tmp is the number of surfaces read and to be written
c*** nt_tmp is the number of values read and to be written

      integer ntrg_tmp, ntrgs_tmp, nt_tmp, nz_tmp, ntrgf_tmp
      parameter (ntrg_tmp=5, ntrgs_tmp=5)
      real (kind=R8) ::
     *     equt_tmp(ntrgs_tmp,ntrg_tmp)
      real (kind=R8) ::
     *     pwttt_tmp(ntrgs_tmp), pwtpl_tmp(ntrgs_tmp),
     ,     pwtnu_tmp(ntrgs_tmp), pwtio_tmp(ntrgs_tmp),
     ,     pwtds_tmp(ntrgs_tmp)
      equivalence
     ,     (equt_tmp(1,1), pwttt_tmp), (equt_tmp(1,2), pwtpl_tmp),
     ,     (equt_tmp(1,3), pwtnu_tmp), (equt_tmp(1,4), pwtio_tmp),
     ,     (equt_tmp(1,5), pwtds_tmp)
      character hlp_trg(ntrgs_tmp)*2, hlp_ht(ntrg_tmp)*8,
     ,          hlp_trgs(ntrgs_tmp,2)*2
      integer, allocatable :: ltgs(:,:)
      integer indt_hlp(nbcd,ntrgs_tmp), nsrf_hlp(ntrgs_tmp)
      logical ltarg_hlp(ntrgs_tmp)
c-----------------------------------------------------------------------

      integer lun, ltrc
      parameter (ltrc=ntrc_ble)
      integer i, j, k, l, n, ic, is, ix, iy,
     .        ixr, iyr, ixf, iyf, nv
      real (kind=R8) :: hlp_h, hlp_p
      logical lcross_hlp, lminus_hlp, lbound_hlp, lcore_hlp
      logical ex, warn_tmp
      character f*(lenfnm), hlp_c, hlp_frm*80, tab, hlp_frf*160
#ifdef B25_EIRENE
      integer, save :: m, js
      integer, save :: ntgs(ntrgs_tmp)
      data ntgs / ntrgs_tmp*1 /
      real (kind=R8), save :: hlp
#endif
c-----------------------------------------------------------------------
c*** Local lists:
c***  ltgs      Eirene surfaces corresponding to the target surfaces
c***  indt_hlp  B2 boundaries corresponding to the target surfaces
c***
c*** Taken from the common blocks:
c***
c*** lcbs: boundary segments representing the core boundary in B2
c*** lcns: 'non-default surfaces' for the core boundary (Eirene)
c*** ltns: 'non-default surfaces' for the targets (Eirene)
c*** lsns: list of Eirene surfaces comprising the recycling sources
c*** msns: msns(1,k): index of first surface in lsns list for Eirene
c***                  recycling source k
c***       msns(2,k): number of surfaces for Eirene recycling source k
c-----------------------------------------------------------------------
c*** Identification suffices:
      data hlp_trgs / 'il','ol','wl','  ','  ',
     ,                'il','iu','ou','ol','wl'/
c*** Column headers:
      data hlp_hg /
     1          'pwr_totl','pwr_tscl','tot_core','tot_targ','tot_wall',
     2          'tot_rad ','heat_cor','ptnt_cor','ntrl_cor','plsm_trg',
     3          'ntrl_trg','plsm_wll','ntrl_wll','brms_rad','ntrl_rad',
     4          'imp_rad ','tot_enrg','tot_enre','tot_enri'/
      data hlp_hs / 'imp_rad_' /
      data hlp_ht /
     1          'pwr_totl','pwr_plsm','pwr_neut','pwr_ionz','pwr_diss'/
c-----------------------------------------------------------------------
c
      save
     +  ng_tmp,
     +  equg_tmp,
     +  pwtot_tmp, pwttr_tmp, pwcrt_tmp, pwtgt_tmp, pwwlt_tmp,
     ,     pwrdt_tmp, pwcrh_tmp, pwcrp_tmp, pwcrn_tmp, pwtgp_tmp,
     ,     pwtgn_tmp, pwwlp_tmp, pwwln_tmp, radbr_tmp, radnu_tmp,
     ,     rdtim_tmp, enrtt_tmp, enrti_tmp, enrte_tmp,
     +  hlp_hg,
     +  nimp_tmp, ni_tmp, ns_tmp, limp_tmp, ntrgf_tmp,
     +  radim_tmp,
     +  hlp_spc, hlp_hs,
     +  nt_tmp, nz_tmp,
     +  equt_tmp,
     +  pwttt_tmp, pwtpl_tmp,
     ,     pwtnu_tmp, pwtio_tmp,
     ,     pwtds_tmp,
     +  hlp_trg, hlp_ht,
     +  ltgs,
     +  indt_hlp, nsrf_hlp,
     +  ltarg_hlp,
     +  lun,
     +  i, j, k, l, ic, is, ix, iy, nv,
     +  hlp_h, hlp_p,
     +  lcross_hlp, lminus_hlp, lbound_hlp, lcore_hlp, ex, warn_tmp,
     +  f, hlp_c, hlp_frm, tab, hlp_frf
c=======================================================================

c*** Initialisation
      do i=1,nglb_tmp !{
        equg_tmp(i)=0.0_R8
      end do !}
      do i=1,nspecies !{
        radim_tmp(i)=0.0_R8
      end do !}
      do i=1,ntrgs_tmp !{
        do j=1,ntrg_tmp !{
          equt_tmp(i,j)=0.0_R8
        end do !}
      end do !}

#ifdef B25_EIRENE
c*** Correct the neutral fluxes

      do j=1,nnlimi+nnstsi !{
        wldpeb(j,0)=0.0_R8
        wldnek(j,0)=0.0_R8
        wldnep(j,0)=0.0_R8
        do k=1,nstrat !{
          wldpeb(j,0)=wldpeb(j,0)+wldpeb(j,k)*sclstra(0,k)
          wldnek(j,0)=wldnek(j,0)+wldnek(j,k)*sclstra(0,k)
          wldnep(j,0)=wldnep(j,0)+wldnep(j,k)*sclstra(0,k)
        end do !}
      end do !}
#endif

c*** Plasma flows across the boundaries

      do ic = 1, nbc !{

c*** First, check the boundary kind.
c*** Does it belong to the core boundary?

        lcore_hlp=.false.
        do j = 1, ncbs !{
          if(ic.eq.lcbs(j)) then !{
            lcore_hlp=.true.
            go to 10
          end if !}
        end do !}
c
   10   do 100 j = 1,ntrgf_tmp !{
c*** Belongs to the target j?
          ltarg_hlp(j)=.false.
          do k = 1,nsrf_hlp(j) !{
            if(ic.eq.indt_hlp(k,j)) then !{
              ltarg_hlp(j)=.true.
              go to 100
            end if !}
          end do !}
  100   continue !}
c
c*** Across or along the field?
        lcross_hlp = bcchar(ic).eq.'S' .or. bcchar(ic).eq.'N'
        lminus_hlp = bcchar(ic).eq.'E' .or. bcchar(ic).eq.'N'
        lbound_hlp = bcchar(ic).eq.'S' .or. bcchar(ic).eq.'N' .or.
     .               bcchar(ic).eq.'E' .or. bcchar(ic).eq.'W'
c
        do i = 1, bc_list_size(ic)
          ix = bc_list_x(i,ic)
          iy = bc_list_y(i,ic)
          if (bcchar(ic).eq.'W') then
            ixr=rightix(bc_list_x(i,ic),bc_list_y(i,ic))
            iyr=rightiy(bc_list_x(i,ic),bc_list_y(i,ic))
          else if (bcchar(ic).eq.'E') then
            ixr=leftix(bc_list_x(i,ic),bc_list_y(i,ic))
            iyr=leftiy(bc_list_x(i,ic),bc_list_y(i,ic))
          else if (bcchar(ic).eq.'S') then
            ixr=topix(bc_list_x(i,ic),bc_list_y(i,ic))
            iyr=topiy(bc_list_x(i,ic),bc_list_y(i,ic))
          else if (bcchar(ic).eq.'N') then
            ixr=bottomix(bc_list_x(i,ic),bc_list_y(i,ic))
            iyr=bottomiy(bc_list_x(i,ic),bc_list_y(i,ic))
          endif
          if (lminus_hlp) then
            ixf=ix
            iyf=iy
          else
            ixf=ixr
            iyf= iyr
          endif
          if (lcross_hlp) then !{
            hlp_h=fht(ixf,iyf,1)-fhj(ixf,iyf,1)+fhi_ext(ixf,iyf,1)
            hlp_p=0.0_R8
            do is=0,ns-1 !{
              hlp_h=hlp_h-fhp(ixf,iyf,1,is)
              hlp_p=hlp_p+fhp(ixf,iyf,1,is)
            end do !}
            do is=0,ns_ext-1 !{
              hlp_p=hlp_p+fa_ext(ixf,iyf,1,is)*ev*
     .         (pt_ext(ix,iy,is)*hy1(ixr,iyr)+
     .          pt_ext(ixr,iyr,is)*hy1(ix,iy))/
     .         (hy1(ixr,iyr)+hy1(ix,iy))
              hlp_h=hlp_h+fa_ext(ixf,iyf,1,is)*(1.0_R8-BoRiS)*
     .          0.5_R8*am_ext(is)*mp*
     .         (ua_ext(ix,iy,is)**2*hy1(ixr,iyr)+
     .          ua_ext(ixr,iyr,is)**2*hy1(ix,iy))/
     .         (hy1(ixr,iyr)+hy1(ix,iy))
              hlp_h=hlp_h+fa_ext(ixf,iyf,1,is)*
     .         (ta_ext(ix,iy,is)*hy1(ixr,iyr)+
     .          ta_ext(ixr,iyr,is)*hy1(ix,iy))/
     .         (hy1(ix,iy)+hy1(ixr,iyr))
c
            end do !}
          else if (lbound_hlp) then !}{
            hlp_h=fht(ixf,iyf,0)-fhj(ixf,iyf,0)+fhi_ext(ixf,iyf,0)
            hlp_p=0.0_R8
            do is=0,ns-1 !{
              hlp_h=hlp_h-fhp(ixf,iyf,0,is)
              hlp_p=hlp_p+fhp(ixf,iyf,0,is)
            end do !}
            do is=0,ns_ext-1 !{
              hlp_p=hlp_p+fa_ext(ixf,iyf,0,is)*ev*
     .         (pt_ext(ix,iy,is)*hx(ixr,iyr)+
     .          pt_ext(ixr,iyr,is)*hx(ix,iy))/
     .         (hx(ix,iy)+hx(ixr,iyr))
              hlp_h=hlp_h+fa_ext(ixf,iyf,0,is)*(1.0_R8-BoRiS)*
     .          0.5_R8*am_ext(is)*mp*
     .         (ua_ext(ix,iy,is)**2*hx(ixr,iyr)+
     .          ua_ext(ixr,iyr,is)**2*hx(ix,iy))/
     .         (hx(ix,iy)+hx(ixr,iyr))
              hlp_h=hlp_h+fa_ext(ixf,iyf,0,is)*
     .         (ta_ext(ix,iy,is)*hx(ixr,iyr)+
     .          ta_ext(ixr,iyr,is)*hx(ix,iy))/
     .         (hx(ix,iy)+hx(ixr,iyr))
c
            end do !}
          end if !}
          if (lminus_hlp) then !{
            hlp_h=-hlp_h
            hlp_p=-hlp_p
          end if !}
c<<<
c          print *,'b2blne: ic,bcchar,lminus,hlp_h=',ic,' ',bcchar(ic),
c     ,                                             ' ',lminus_hlp,hlp_h
c>>>
c
c*** Now add the data to the proper accumulators
          if (lcore_hlp) then !{
            pwcrh_tmp=pwcrh_tmp+hlp_h
            pwcrp_tmp=pwcrp_tmp+hlp_p
          else !}{
            if (lcross_hlp) then !{
              pwwlp_tmp=pwwlp_tmp+hlp_h+hlp_p
            else if (lbound_hlp) then !}{
              pwtgp_tmp=pwtgp_tmp+hlp_h+hlp_p
            end if !}
          end if !}
          do j = 1,ntrgf_tmp !{
            if (ltarg_hlp(j)) then !{
              pwtpl_tmp(j)=pwtpl_tmp(j)+hlp_h
              pwtio_tmp(j)=pwtio_tmp(j)+hlp_p
            end if !}
          end do !}
        end do !}
      end do !}
c
c*** Correction for reflected neutrals
c
c<<<
c      print *,'b2blne: ntrgf_tmp = ',ntrgf_tmp
c      print *,'indt =',(indt_hlp(1,i),i=1,ntrgf_tmp)
c      print '(a,99a4)',' crcstra : ',crcstra
c      print '(a,99a4)',' bcchar  : ',bcchar
c      print '(a,99i4)',' rcpos   : ',rcpos
c      print '(a,99i4)',' bcpos   : ',bcpos
c>>>
#ifdef B25_EIRENE
      do ic=1,nstrat !{
        if(crcstra(ic).eq.'E' .or. crcstra(ic).eq.'W' .or.
     .     crcstra(ic).eq.'N' .or. crcstra(ic).eq.'S') then !{
          l=msns(1,ic)
          m=msns(2,ic)
          hlp=0.0_R8
          do i=l,l+m-1 !{
            k=lsns(i)
            if(k.lt.0) k=nnlimi-k
            if(k.gt.0) hlp=hlp+wldpeb(k,ic)
          end do !}
          k=0
          if(crcstra(ic).eq.'E' .or. crcstra(ic).eq.'W') then !{
            do i=1,nbc !{
              if(bcchar(i).eq.crcstra(ic) .and.
     &           bc_list_x(1,i).eq.rc_list_x(1,ic)) k=i
            end do !}
          end if !}
c<<<
c          print *,'b2blne: ic,k= ',ic,k
c>>>
          if(k.gt.0) then !{
            j=0
            do i=1,ntrgf_tmp-1 !{
              if(indt_hlp(1,i).eq.k) j=i
            end do !}
            if(j.gt.0) then !{
              k=j
              pwtgp_tmp=pwtgp_tmp+hlp
            else !}{
              k=0
              write(*,*) 'b2blne: failed to relate a stratum',ic,
     ,                                              ' to a B2 surface.'
            end if !}
          else !}{
            k=ntrgf_tmp
            pwwlp_tmp=pwwlp_tmp+hlp
          end if !}
          if(k.gt.0) pwtpl_tmp(k)=pwtpl_tmp(k)+hlp
        end if !}
      end do !}
#endif
c
c*** Radiation
c
      radnu_tmp=-rad_neutr
      radbr_tmp=-b2brem
      do is = 1, ni_tmp !{
        radim_tmp(is)=-rad_imp(is)
      end do !}
      if(ni_tmp.le.0) radim_tmp(1)=0.0_R8

#ifdef B25_EIRENE
c
c*** Neutral fluxes
c
csw 18nov2013 check for nnatmi (known after 1st EIRENE call)
      if(nnatmi > 0) then

c*** flux across the core boundary
       do is=1,ncns !{
        js=lcns(is)+nnlimi
        pwcrn_tmp=pwcrn_tmp+wldnek(js,0)+wldnep(js,0)
       end do !}
c*** flux across the targets
       do is=1,ntns !{
        js=ltns(is)+nnlimi
        pwtgn_tmp=pwtgn_tmp+wldnek(js,0)+wldnep(js,0)
       end do !}
c*** flux across the side walls
       do js=1,nnlimi !{
        pwwln_tmp=pwwln_tmp+wldnek(js,0)+wldnep(js,0)
       end do !}
c*** flux across the specified zones
       do j=1,ntrgf_tmp !{
        do l=1,ntgs(j) !{
          js=ltgs(l,j)
          if(js.lt.0) js=nnlimi-js
          if(js.gt.0) then !{
            pwtnu_tmp(j)=pwtnu_tmp(j)+wldnek(js,0)
            pwtds_tmp(j)=pwtds_tmp(j)+wldnep(js,0)
          else !}{
            do js=1,nnlimi !{
              pwtnu_tmp(j)=pwtnu_tmp(j)+wldnek(js,0)
              pwtds_tmp(j)=pwtds_tmp(j)+wldnep(js,0)
            end do !}
          end if !}
        end do !}
       end do !}
      endif !nnatmi
#endif
c
c*** Energy content
c
      do ix=0,nx-1 !{
        do iy=0,ny-1 !{
          if(region(ix,iy,0).ne.0) then !{
            enrte_tmp=enrte_tmp+ne(ix,iy)*te(ix,iy)*vol(ix,iy)
            enrti_tmp=enrti_tmp+ni(ix,iy,1)*ti(ix,iy)*vol(ix,iy)
            do is =0,ns_ext-1 !{
              enrti_tmp = enrti_tmp +
     .         na_ext(ix,iy,is)*ta_ext(ix,iy,is)*vol(ix,iy)
            end do !}
          end if !}
        end do !}
      end do !}
      enrte_tmp=1.5_R8*enrte_tmp
c***            kinetic energy of the plasma flow is still neglected!!!
      enrti_tmp=1.5_R8*enrti_tmp
c
c*** Totals
c
      pwtgp_tmp=-pwtgp_tmp
      pwwlp_tmp=-pwwlp_tmp
c
      pwcrt_tmp=pwcrh_tmp+pwcrp_tmp-pwcrn_tmp
      pwwlt_tmp=pwwlp_tmp+pwwln_tmp
      pwtgt_tmp=pwtgp_tmp+pwtgn_tmp
      pwrdt_tmp=radbr_tmp+radnu_tmp
      do k=1,nimp_tmp !{
        rdtim_tmp=rdtim_tmp+radim_tmp(k)
        pwrdt_tmp=pwrdt_tmp+radim_tmp(k)
      end do !}
      pwtot_tmp=pwcrt_tmp-(pwwlt_tmp+pwtgt_tmp+pwrdt_tmp)
      enrtt_tmp=enrti_tmp+enrte_tmp
      pwttr_tmp=pwtot_tmp/enrtt_tmp
c
      do k=1,ntrgf_tmp !{
        pwtpl_tmp(k)=-pwtpl_tmp(k)
        pwtio_tmp(k)=-pwtio_tmp(k)
        pwttt_tmp(k)=pwtpl_tmp(k)+pwtnu_tmp(k)+pwtio_tmp(k)+pwtds_tmp(k)
      end do !}
c
      n=max(ni_tmp,1)
      if(ank_tracing.eq.0) return
      if(filedata(ltrc)) then !{
c
c*** Write the data in a separate file for plotting...
c
        write(lun,hlp_frm) tim,equg_tmp,
     ,                     (radim_tmp(i),i=1,n),
     ,                    ((equt_tmp(i,j),j=1,nt_tmp),i=1,nz_tmp)
      else !}{
c
c*** ... or just print out the data
c
        write(6,'(8x,20(1x,a8,1x))') '   sum  ','targ_tot',
     ,         'wall_tot',' rad_tot',' plsm s ',' ioniz s',' plsm tg',
     ,         ' plsm wl',' neut tg',' neut wl','Bremsstr',' rad neu',
     ,         'rad imp:'
        write(6,'(a8,1p,20e10.2)') 'power >>', pwtot_tmp, pwtgt_tmp,
     ,        pwwlt_tmp, pwrdt_tmp, pwcrh_tmp, pwcrp_tmp, pwtgp_tmp,
     ,        pwwlp_tmp, pwtgn_tmp, pwwln_tmp, radbr_tmp, radnu_tmp,
     ,                                     (radim_tmp(i),i=1,n)
      end if !}
      return
c=======================================================================
c
      entry b2blnei(nx,ny,ns,BoRiS,iline)
c-----------------------------------------------------------------------
c*** Initialise the file for tracing data
c*** Only called with filedata(ntrc_ble)=.true.
c-----------------------------------------------------------------------
c
      allocate(ltgs(nlim,ntrgs_tmp))
      allocate(limp_tmp(nspecies))
      allocate(radim_tmp(nspecies))
      allocate(hlp_spc(nspecies))
      call ipgetr ('b2mndr_stim', stim)
c
      tab=char(9)
      nimp_tmp=0
      do is = 1, nspecies !{
        if(nfluids(is).gt.1) then !{
          nimp_tmp=nimp_tmp+1
          limp_tmp(nimp_tmp)=is
          hlp_spc(nimp_tmp)=species_list(is)
        end if !}
      end do !}
      if(nimp_tmp.le.0) hlp_spc(1)='-'

c*** Determine the target surfaces.

c*** First, targets in B2

      warn_tmp=.true.
      k=0
      do j=1,nbc !{
        if(bcchar(j).eq.'W' .or. bcchar(j).eq.'E') then !{
          if(minval(bc_list_y(1:bc_list_size(j),j)).le.0 .and.
     &       maxval(bc_list_y(1:bc_list_size(j),j)).ge.ny-1) then !{
            k=k+1
            if(k.le.ntrgs_tmp) then !{
              nsrf_hlp(k)=1
              indt_hlp(1,k)=j
              indt_hlp(2,k)=k
            else !}{
              if(warn_tmp) write(0,*) '!!! Too many targets (',k,' ) ',
     ,                'found. Check the value of NTRGS_TMP in b2blne.F'
              warn_tmp=.false.
            end if !}
          end if !}
        end if !}
      end do !}

c*** Sort the targets
      do l=1,k !{
        do i=l+1,k !{
          if(bc_list_x(1,indt_hlp(1,l)).gt.bc_list_x(1,indt_hlp(1,i)))
     >     then !{
            j=indt_hlp(1,l)
            indt_hlp(1,l)=indt_hlp(1,i)
            indt_hlp(1,i)=j
          end if !}
        end do !}

c*** Eirene surfaces corresponding to the target surfaces
        ltgs(1,l)=-abs(ltns(l))
      end do !}

c*** Assign the identification suffices

      if(k.le.2) then !{
        l=1
      else if(k.le.4) then !}{
        l=2
      end if !}
      if(k.le.4) then !{
        do i=1,k !{
          hlp_trg(i)=hlp_trgs(i,l)
        end do !}
      else !}{
        write(*,*) '@ b2blne: too many targets. Check the hlp_trgs ',
     ,                                                          'array'
        do i=1,k !{
          write(hlp_trg(i),'(i2.2)') i
        end do !}
      end if !}

c<<<
c      print *,'b2blne: targets ',(indt_hlp(1,i),':',hlp_trg(i),':',
c     ,                                             ltgs(1,i),' ',i=1,k)
c>>>

c*** Walls in B2

      if(k.lt.ntrgs_tmp) then !{
        k=k+1
        l=0
        do j=1,nbc !{
          ex=.false.
          do i=1,k-1 !{
            ex= ex .or. indt_hlp(1,i).eq.j
          end do !}
          do i=1,ncbs !{
            ex= ex .or. lcbs(i).eq.j
          end do !}
          if(.not.ex) then !{
            l=l+1
            indt_hlp(l,k)=j
          end if !}
        end do !}
        nsrf_hlp(k)=l
        hlp_trg(k)='wl'
        ltgs(1,k)=0
      else !}{
        write(*,*) '@ b2blne: no room for the wall data. ',
     ,                                                'Check ntrgs_hlp'
      end if !}
      ntrgf_tmp=k

c
c*** lun   :  logical unit number for the output file (=ntrc(ltrc))
c*** f     :  the associated file name
      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)
c
      call find_file(f,ex)
      open(lun,file=f,err=90)
      rewind(lun)
      if(ex) read(lun,*,end=91) hlp_c
      goto 93
   91 ex = .false.
   93 continue
      rewind(lun)
      if(ex.and.stim.lt.0.0_R8) then !{
        read(lun,*,err=92) hlp_c,nv,ng_tmp,nt_tmp,nz_tmp,ni_tmp,ns_tmp
        lv_trc(ltrc)=nv
        warn_tmp=.false.
        if(nv.ne.lv) then !{
          warn_tmp=.true.
          write(*,'(a,i8,3a,i8)')
     ,     '!!! b2blne: current version ',lv,
     ,     ' differs from that in the file ',trim(f),' : ',nv
          if(nv.lt.960825 .or. (nv.le.961229 .and. ni_tmp.eq.0)) then !{
c*** It was a major change of the file format - impossible to continue
            write(*,*) 'No ',trim(f),' tracing file - ',
     ,                                'remove the existing file first!'
            filedata(ltrc)=.false.
            return
          end if !}
        end if !}
        if(ni_tmp.ne.nimp_tmp) then !{
          warn_tmp=.true.
          write(*,'(a,i2,3a,i2)')
     ,     '!!! b2blne: current number of impurity species ',nimp_tmp,
     ,     ' differs from that in the file ',trim(f),' : ',ni_tmp
          filedata(ltrc)=.false.
        end if !}
        if(nz_tmp.ne.ntrgf_tmp) then !{
          warn_tmp=.true.
          write(*,'(a,i2,3a,i2)')
     ,     '!!! b2blne: current number of targets ',ntrgf_tmp,
     ,     ' differs from that in the file ',trim(f),' : ',nz_tmp
          filedata(ltrc)=.false.
        end if !}
        if(ng_tmp.ne.nglb_tmp) then !{
          warn_tmp=.true.
          write(*,'(2a,i4,3a,i4)')
     ,     '!!! b2blne: current number of columns',
     ,     ' for general data ',nglb_tmp,
     ,     ' differs from that in the file ',trim(f),' : ',ng_tmp
          if(ng_tmp.gt.nglb_tmp) then !{
            filedata(ltrc)=.false.
          end if !}
        end if !}
        if(ns_tmp.ne.nspc_tmp) then !{
          warn_tmp=.true.
          write(*,'(2a,i4,3a,i4)')
     ,     '!!! b2blne: current number of columns',
     ,     ' for impurity data ',nspc_tmp,
     ,     ' differs from that in the file ',trim(f),' : ',ns_tmp
          if(ns_tmp.gt.nspc_tmp) then !{
            filedata(ltrc)=.false.
          end if !}
        end if !}
        if(nt_tmp.ne.ntrg_tmp) then !{
          warn_tmp=.true.
          write(*,'(2a,i4,3a,i4)')
     ,     '!!! b2blne: current number of columns',
     ,     ' for impurity data ',ntrg_tmp,
     ,     ' differs from that in the file ',trim(f),' : ',nt_tmp
          if(nt_tmp.gt.ntrg_tmp) then !{
            filedata(ltrc)=.false.
          end if !}
        end if !}
        if(.not.filedata(ltrc)) then !{
          write(*,*) 'No ',trim(f),' tracing file - ',
     ,                                'remove the existing file first!'
          return
        end if !}
        if(warn_tmp) write(*,*) ' You should remove the file!'
c*** Read to the end of the file
   50   read(lun,'(a1)',end=55) hlp_c
        go to 50
   55   backspace lun
      else !}{
c*** New file
        ng_tmp=nglb_tmp
        nt_tmp=ntrg_tmp
        nz_tmp=ntrgf_tmp
        ni_tmp=nimp_tmp
        n=max(ni_tmp,1)
        ns_tmp=nspc_tmp
        write(lun,700) lv,ng_tmp,nt_tmp,nz_tmp,ni_tmp,ns_tmp,label
        hlp_frf='(9h(5h#time,,i3,9h(1h'//tab//',a8),,' //
     /                       'i3,12h(1h'//tab//',a8,a2),,' //
     /                       'i3,16h(1h'//tab//',a8,1h_,a2)))'
        write(hlp_frm,hlp_frf) ng_tmp,ns_tmp*n,nt_tmp*nz_tmp
        write(lun,hlp_frm) (hlp_hg(i),i=1,ng_tmp),
     ,                  ((hlp_hs(i),hlp_spc(j),i=1,ns_tmp),j=1,n),
     ,                  ((hlp_ht(i),hlp_trg(j),i=1,nt_tmp),j=1,nz_tmp)
      end if !}
      n=max(ni_tmp,1)
      write(hlp_frm,'(10h(1p,e16.7,,i3,12h(1h'//tab//',e12.3)))')
     -                                    ng_tmp+ns_tmp*n+nt_tmp*nz_tmp
      return
c=======================================================================
   90 write(*,*) 'b2blne: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      return
c
   92 write(*,*) 'b2blne: wrong format of the file ',trim(f)
      filedata(ltrc)=.false.
      return
c======================================================================
c
      entry b2blnec(nx,ny,ns,BoRiS,iline)
c----------------------------------------------------------------------
c*** - to be called by b2co for locating truncation for savefile
c----------------------------------------------------------------------
c
      if(.not.filedata(ntrc_ble)) return
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
  190 write(*,*) 'b2blne: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      iline=-1
      return
c=======================================================================
      entry b2blnes(nx,ny,ns,BoRiS,iline)
c-----------------------------------------------------------------------
c*** Save the values of the plasma fluxes across the boundaries
c-----------------------------------------------------------------------
c
      return
c=======================================================================
      entry b2blnef(nx,ny,ns,BoRiS,iline)
c-----------------------------------------------------------------------
c*** Finalization step
c-----------------------------------------------------------------------
c
      return
c=======================================================================
  700 format('#   ',6i8/'#  ',a)
c=======================================================================
      end subroutine b2blne
