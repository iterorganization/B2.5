*-----------------------------------------------------------------------
*.specification
      subroutine b2stbc_fb_post (nx, ny, ns, ismain,
     1  na, ne, ni, te, ti, fna, fhe, fhi, fch,
     &  sna0, smo0, she0, shi0, sch0, sne0,
     &  main_call, diagno)
      use b2mod_types
      use b2mod_geo
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_feedback
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, diagno
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  te(-1:nx, -1:ny), ti(-1:nx, -1:ny)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     *  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1)
      logical main_call                                                  !xpb
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3), sne0(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBC_FB provides linearised expressions for the sources of
*     density, parallel momentum, electron heat, atom heat, and charge,
*     due to feedback boundary conditions after the normal boundary
*     conditions have been applied.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, ivac, ireg
      real (kind=R8) ::
     4  cum_pressure, cum_density, cum_area,
     5  ave_pressure, ave_density, ave_temp,
     6  vacuum_corr
      save ncall
*   ..procedures
      external xertst
      external b2xvsg
*   ..initialisation
      data ncall/0/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbc_fb')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call read_b2mod_feedback_control(99)
       call read_b2mod_feedback(99,'    ',feedback_namelist_used)
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ni, ne, te, ti
       call b2xvsg (n2, ni(-1,-1,0), 1, 'ni0', '.gt.')
       call b2xvsg (n2, ni(-1,-1,1), 1, 'ni1', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
      endif

! process feedbacks in b2.feedback_control.parameters
      if(vacuum_communication.gt.0) then
! this allows for a communication of particle fluxes across vacuum regions
        do ivac=1, vacuum_communication
          do is = 0, ns-1
! this option only applies to neutrals
            if(is_neutral(is)) then
! calculate an average pressure and an average temperature
              cum_pressure=0.0_R8
              cum_density=0.0_R8
              cum_area=0.0_R8
              do ireg=1, vacuum_communication_nreg(ivac)
                iy=vacuum_communication_iy(ireg, ivac)
                do ix=vacuum_communication_ix1(ireg, ivac),
     1                vacuum_communication_ix2(ireg, ivac)
                  cum_pressure = cum_pressure +
     1             ti(ix,iy)*na(ix,iy,is)*gs(ix,iy,2)
                  cum_density = cum_density +
     1             na(ix,iy,is)*gs(ix,iy,2)
                  cum_area = cum_area + gs(ix,iy,2)
                enddo
              enddo
              ave_pressure = cum_pressure / cum_area
              ave_density = cum_density / cum_area
              ave_temp = ave_pressure / ave_density
              write(*,'(a,a,1x,i2,1x,i3,1x,1p,3g12.3)')
     1         'Vacuum-Communication: ivac, is, ave_pressure, ',
     2         'ave_density, ave_temp ',
     1         ivac, is, ave_pressure, ave_density, ave_temp/ev
! adjust the density boundary condition based on the difference between the average pressure and the local pressure
              do ireg=1, vacuum_communication_nreg(ivac)
                iy=vacuum_communication_iy(ireg, ivac)
                do ix=vacuum_communication_ix1(ireg, ivac),
     1           vacuum_communication_ix2(ireg, ivac)
                  if(vacuum_communication_method(ivac).eq.1) then
! try to set a flux
                    vacuum_corr = (vacuum_communication_beta(is,ivac) *
     1               ave_pressure - ti(ix,iy)*na(ix,iy,is)) / ave_temp *
     2               vacuum_communication_alpha(is,ivac)
                    sna0(bottomix(ix,iy),bottomiy(ix,iy),0,is) =
     .               sna0(bottomix(ix,iy),bottomiy(ix,iy),0,is) +
     1               vacuum_corr * gs(ix,iy,2)
                  else if(vacuum_communication_method(ivac).eq.2) then
! try to set a density based on pressure equality
                    vacuum_corr = vacuum_communication_beta(is,ivac) *
     1               ave_pressure / ti(bottomix(ix,iy),bottomiy(ix,iy))
                    sna0(bottomix(ix,iy),bottomiy(ix,iy),0,is) =
     .               sna0(bottomix(ix,iy),bottomiy(ix,iy),0,is) +
     1               vacuum_communication_alpha(is,ivac) *
     2               vacuum_corr * gs(ix,iy,2)
                    sna0(bottomix(ix,iy),bottomiy(ix,iy),1,is) =
     .               sna0(bottomix(ix,iy),bottomiy(ix,iy),1,is) -
     1               vacuum_communication_alpha(is,ivac) * gs(ix,iy,2)
                  else if(vacuum_communication_method(ivac).eq.3) then
! try to set a density based on density equality
                    vacuum_corr = vacuum_communication_beta(is,ivac) *
     1               ave_density
                    sna0(bottomix(ix,iy),bottomiy(ix,iy),0,is) =
     .               sna0(bottomix(ix,iy),bottomiy(ix,iy),0,is) +
     1               vacuum_communication_alpha(is,ivac) *
     2               vacuum_corr * gs(ix,iy,2)
                    sna0(bottomix(ix,iy),bottomiy(ix,iy),1,is) =
     .               sna0(bottomix(ix,iy),bottomiy(ix,iy),1,is) -
     1               vacuum_communication_alpha(is,ivac) * gs(ix,iy,2)
                  else
                    write(*,*) 'Vacuum communication method ',
     1               vacuum_communication_method(ivac),
     2               ' not implemented yet'
                  endif
                enddo
              enddo
            endif
          enddo
        enddo
      endif

* ..return
      if(ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      if(ncall.gt.1 .or. .not.main_call) return
      call b2stbc_fb_final
      return
*-----------------------------------------------------------------------
*.end b2stbc_fb_post

      end subroutine b2stbc_fb_post

* ..finalise feedback
      subroutine b2stbc_fb_final
      use b2mod_types
      use b2mod_feedback
      use b2mod_usrspb
      implicit none
      feedback_namelist_used = feedback_namelist_used .or.
     .        max(saved_cbshe_core, saved_cbshi_core,
     .            saved_cbsna_core, saved_cbsch_core,
     .            saved_cbsna_pfr1, saved_cbsna_pfr2,
     .            saved_cbsna_sol, maxval(saved_na_feedback_actuator))
     .            .gt.0.0_R8
      if(feedback_namelist_used) then
        call write_b2mod_feedback(99,'    ')
      endif

      return
      end subroutine b2stbc_fb_final

*-----------------------------------------------------------------------
*.specification
      subroutine b2stbc_fb_pre (nx, ny, ns, ismain, BoRiS,
     1  vol, gs, na, ne, ni, te, ti, kinrgy,
     &  fna, fna_mdf, fna_fcor, fna_nodrift,
     &  fna_he, fnaPSch, fna_eir, fna_32, fna_52,
     &  fhe, fhi, fhm, fch,
     &  sna0, smo0, she0, shi0, sch0, sne0,
     &  main_call, diagno)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_feedback
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_diag
      use b2mod_time
      use b2mod_neutrals_namelist
      use b2mod_numerics_namelist
      use b2mod_boundary_namelist
      use b2mod_user_namelist
      use b2mod_usrspb
     ,   , only : Q_pneut, calc_pneut, calc_pneut_init
      use b2mod_geo
     ,   , only : hy,hz,hx
#ifdef B25_EIRENE
      use eirmod_wneutrals
      use eirmod_braeir
#endif
      use b2mod_subsys
      use b2mod_external
      use b2mod_plasma
     , ,only: na0, ne0, fna0, kinrgy0, ni0,
     ,        nap, nep, fnap, kinrgyp, fht
      use b2mod_sources
     ,   , only : rqrad, rqbrm
      use b2mod_dimensions
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, diagno
      real (kind=R8) :: BoRiS,
     *  vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  te(-1:nx, -1:ny), ti(-1:nx, -1:ny), kinrgy(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1),
     &  fhm(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_mdf(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_fcor(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_nodrift(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_he(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_32(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_52(-1:nx,-1:ny,0:1,0:ns-1),
     &  fnaPSch(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_eir(-1:nx,-1:ny,0:1,0:ns-1),
     &  fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1)
      logical main_call                                                  !xpb
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3), sne0(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBC_FB_PRE provides a change to the plasma or boundary conditions
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, ispecies, is0, 
     & mdf_fnb, style_b2mndt, 
#ifdef B25_EIRENE
     &        i, j, k, l, npump_list,
     &        pump_list(DEF_NLIM+DEF_NSTS,0:DEF_NATM-1),
#endif
     &        ib, iss, ixl1, ixr1, ixl2, ixr2, istrai, icount
      integer, save :: is_start, is_end
      integer, save :: ncall_old, ncall_gap, n_flux_zero
      integer, save :: fhiycore_kinetic_energy, ncallfeedback
      integer, save :: coreregno, coreregn2, solregno,
     &                 pfrregno1, pfrregno2,
     &                 isfeedback, nesepm_istra, use_eirene,
     &                 nxtl, nxtr, puff_strata(0:DEF_NATM-1)
      real (kind=R8) ::
     &  na_feedback_current, na_feedback_rescale, vol_feedback,
     &  na_feedback_prev(0:DEF_NATM-1), na_feedback_save(0:DEF_NATM-1),
     &  dt_prev, sum_e, sum_s,
     &  total_flux, flux_old, target_flux,
     &  fnaycore, fnaycore_alpha, fnaycurr,
     1  fheycore, fheycore_alpha, fheycurr,
     1  fhiycore, fhiycore_alpha, fhiycurr,
     2  fchycore, fchycore_alpha, fchycurr,
     3  nesepm, nesepm_sol, nepedm_sol, nesepm_pfr,
     3  nesepm_alpha, nesepm_beta, nesepm_gamma, nesepm_curr,
     3  nesepm_minpuff, nesepm_maxpuff, nesepm_overshoot,
     3  nesepm_error, nesepm_error_old, nesepm_integral_error,
     4  volrec_sol, volrec_alpha, volrec_beta, volrec_overshoot,
     4  volrec_curr,
     5  ndes, ndes_sol, ntot, den(0:1), rcbsna(0:1),
     6  cbsnafac, cbsnainc, eps,
#ifdef B25_EIRENE
     &  flxic, flxnc, flxex, flxsp, flxnp, sna_pump, absr, 
     &  absr_save(0:DEF_NATM-1),
#endif
     7  private_flux_puff, private_flux_area,
     8  neutral_sources_rescale,
     9  save_core_flux_while_overshot, save_gas_puff_while_overshot,
     9  save_pfr1_puff_while_overshot, save_pfr2_puff_while_overshot
      real (kind=R8), allocatable :: charge_frac(:)
      integer, save :: jsep, jxa, jxi, iyped
      integer, save :: b2sral_style = 2
      character*72 message
      logical, save :: consistent
      logical consistent_core, consistent_sol, consistent_pfr1,
     .        consistent_pfr2, zeroslope
      save ncall, dt_prev, na_feedback_prev, na_feedback_save,
     1  fnaycore, fnaycore_alpha, fheycore, fheycore_alpha,
     2  fhiycore, fhiycore_alpha, fchycore, fchycore_alpha,
     3  nesepm, nesepm_sol, nepedm_sol, nesepm_pfr,
     3  nesepm_alpha, nesepm_beta, nesepm_gamma, nesepm_error_old,
     3  nesepm_minpuff, nesepm_maxpuff, nesepm_overshoot,
     3  nesepm_integral_error,
     4  volrec_sol, volrec_alpha, volrec_beta, volrec_overshoot,
     5  ndes, ndes_sol, den, rcbsna,
     6  cbsnafac, 
#ifdef B25_EIRENE
     &  absr_save, pump_list,
#endif
     7  private_flux_puff, private_flux_area, flux_old,
     8  neutral_sources_rescale,
     9  save_core_flux_while_overshot, save_gas_puff_while_overshot,
     9  save_pfr1_puff_while_overshot, save_pfr2_puff_while_overshot
      save zeroslope
*   ..procedures
      intrinsic abs, mod, nint
      external xertst, ipgetr, ipgeti
      external b2xvsg
*   ..initialisation
      parameter (npump_list=DEF_NATM*(DEF_NLIM+DEF_NSTS))
      data ncall/0/, ncall_old/0/, ncall_gap/5/, ncallfeedback/0/,
     &  n_flux_zero/0/, na_feedback_save/DEF_NATM*0.0_R8/,
     &  na_feedback_prev/DEF_NATM*0.0_R8/, use_eirene/0/,
     1  fnaycore/-1.0e30_R8/, fnaycore_alpha/1.0e-1_R8/,
     1  fheycore/0.0_R8/, fheycore_alpha/1.0e-1_R8/,
     2  fhiycore/0.0_R8/, fhiycore_alpha/1.0e-1_R8/,
     2  fhiycore_kinetic_energy/0/,
     2  fchycore/-1.0e30_R8/, fchycore_alpha/1.0e-1_R8/,
     3  nesepm/0.0_R8/, nesepm_sol/0.0_R8/, nepedm_sol/0.0_R8/,
     3  nesepm_pfr/0.0_R8/, nesepm_alpha/1.0e-3_R8/,
     3  nesepm_beta/0.0_R8/, nesepm_gamma/0.99_R8/,
     3  nesepm_minpuff/0.0_R8/, nesepm_maxpuff/0.0_R8/,
     3  nesepm_overshoot/0.0_R8/, nesepm_error_old/0.0_R8/,
     3  nesepm_integral_error/0.0_R8/,
     4  volrec_sol/0.0_R8/, volrec_alpha/1.0e-3_R8/,
     4  volrec_beta/1.0_R8/, volrec_overshoot/0.0_R8/,
     5  ndes/0.0_R8/, ndes_sol/0.0_R8/,
     6  cbsnafac/0.01_R8/,
#ifdef B25_EIRENE
     &  absr_save/DEF_NATM*0.0_R8/, pump_list/NPUMP_LIST*0/,
#endif
     7  private_flux_puff/0.0_R8/,
     8  neutral_sources_rescale/1.0_R8/,
     9  save_core_flux_while_overshot/0.0_R8/,
     9  save_gas_puff_while_overshot/0.0_R8/,
     9  save_pfr1_puff_while_overshot/0.0_R8/,
     9  save_pfr2_puff_while_overshot/0.0_R8/,
     i  coreregn2/4/, solregno/3/,
     i  pfrregno1/0/, pfrregno2/2/, isfeedback/0/, nesepm_istra/-1/,
     l  zeroslope/.false./, puff_strata/DEF_NATM*0/, eps/1.0E-30_R8/
#ifdef COMPLICATED
#ifndef B25_EIRENE
      integer i
#endif
      real (kind=R8) ::
     .  teave1, teave2, neave1, neave2, tenext,
     .  tiave1, tiave2, niave1, niave2, tinext
#endif
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbc_fb_pre')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call get_jsep(nx,ny,jxi,jxa,jsep)
       call read_b2mod_feedback_control(99)
       call read_b2mod_feedback(99,'    ',feedback_namelist_used)
       dt_prev = dteff
*   ..test for input consistency, i.e. are the cbs.. arrays used in the b.c settings
       call ipgeti ('b2sral_style', b2sral_style)
       call get_boundary_namelist
*   ..fill gas puff strata data (used for the gas puff actuator only)
       do ispecies = 0, nspecies-1
         if (na_feedback_actuator(ispecies).eq.1) then
           if (na_feedback_ib(ispecies).lt.0) then
             istrai = -na_feedback_ib(ispecies)
             if ((crcstra(istrai).eq.'C') .and.
     &           (species_start(istrai).le.eb2spcr(ispecies+1)) .and.
     &           (species_end(istrai).ge.eb2spcr(ispecies+1)))
     .         puff_strata(ispecies) = istrai
           else
             do istrai = 1, nstrat
               if ((crcstra(istrai).eq.'C') .and.
     &             (species_start(istrai).le.eb2spcr(ispecies+1)) .and.
     &             (species_end(istrai).ge.eb2spcr(ispecies+1)))
     .           puff_strata(ispecies) = istrai
             end do
           end if
         end if
       end do 
       consistent = .true.
       do ispecies = 0, nspecies-1
         if (na_feedback_ib(ispecies).le.0) cycle
         is0 = eb2spcr(ispecies+1)
         is_start = eb2spcr(ispecies+1)+1
         is_end = eb2spcr(ispecies+1)+nfluids(ispecies+1)
         consistent = consistent .and.
     &    ((na_feedback_choice(ispecies) .eq. 0 .or.
     &      na_feedback_option(ispecies) .eq. 0 .or.
     &      na_feedback_actuator(ispecies) .eq. 0) .or.
     &      na_feedback_actuator(ispecies) .eq. 2 .or.
     &      na_feedback_actuator(ispecies) .eq. 4 .or.
     &     (na_feedback_actuator(ispecies) .eq. 1 .and.
     &      bccon(is0,na_feedback_ib(ispecies)) .eq. 11) .or.
     &     (na_feedback_actuator(ispecies) .eq. 1 .and.
     &      bccon(is0,na_feedback_ib(ispecies)) .eq. 13) .or.
     &     (na_feedback_actuator(ispecies) .eq. 3 .and.
     &      maxval(bccon(is_start:is_end,
     .                   na_feedback_ib(ispecies))).eq.11 .and.
     &      minval(bccon(is_start:is_end,
     .                   na_feedback_ib(ispecies))).eq.11))
       end do
       consistent = consistent .or. lbndusr .or. b2sral_style.eq.1 .or.
     &  boundary_namelist.lt.1
       do ispecies = 0, nspecies-1
         if (na_feedback_ib(ispecies).ge.0) cycle
         consistent = consistent .and.
     &    ((na_feedback_choice(ispecies) .eq. 0 .or.
     &      na_feedback_option(ispecies) .eq. 0 .or.
     &      na_feedback_actuator(ispecies) .eq. 0) .or.
     &     (.not.(na_feedback_actuator(ispecies).eq.3) .and.
     &      (.not.(na_feedback_actuator(ispecies).eq.1) .or.
     &       puff_strata(ispecies).ne.0)))
       end do
       do ispecies = 0, nspecies-1
         if (na_feedback_ib(ispecies).ne.0) cycle
         consistent = consistent .and.
     &    ((na_feedback_choice(ispecies) .eq. 0 .or.
     &      na_feedback_option(ispecies) .eq. 0 .or.
     &      na_feedback_actuator(ispecies) .eq. 0) .or.
     &     (na_feedback_actuator(ispecies).eq.2 .or.
     &      na_feedback_actuator(ispecies).eq.4))
       end do
       do ispecies = 0, nspecies-1
         if (na_feedback_actuator(ispecies).ne.2) cycle
         do is = is_start, is_end
           do ib = 1, nbc
             consistent = consistent .and.
     &        ((na_feedback_choice(ispecies) .eq. 0 .or.
     &          na_feedback_option(ispecies) .eq. 0 .or.
     &          na_feedback_actuator(ispecies) .eq. 0) .or.
     &         (bccon(is,ib).eq.0 .or. bccon(is,ib).eq.2 .or.
     &          bccon(is,ib).eq.3 .or. bccon(is,ib).eq.4 .or.
     &          bccon(is,ib).eq.9 .or. bccon(is,ib).eq.10.or.
     &          bccon(is,ib).eq.14.or. bccon(is,ib).eq.15.or.
     &         (bccon(is,ib).eq.5 .and.conpar(is,ib,1).eq.0.0_R8).or.
     &         (bccon(is,ib).eq.8 .and.conpar(is,ib,1).eq.0.0_R8)))
           end do
         end do
       end do
       do ispecies = 0, nspecies-1
         consistent = consistent .and.
     &    ((na_feedback_choice(ispecies) .eq. 0 .or.
     &      na_feedback_option(ispecies) .eq. 0 .or.
     &      na_feedback_actuator(ispecies) .eq. 0) .or.
     &      (.not.(na_feedback_actuator(ispecies).eq.4) .or.
     &       (na_feedback_actuator(ispecies).eq.4 .and.
     &        na_feedback_ib(ispecies).le.0) .or.
     &       (na_feedback_actuator(ispecies).eq.4 .and.
     &        na_feedback_ib(ispecies).gt.0 .and.
     &        na_feedback_ib(ispecies).le.n_spcsrf .and.
     &        i_spcsrf(na_feedback_ib(ispecies)) .le.
     &        j_spcsrf(na_feedback_ib(ispecies)))))
       end do
       call xertst (consistent,
     .  'feedback scheme not consistent with B.C.''s!')
*   ..test if feedback option is consistent with the actuator
       do ispecies = 0, nspecies-1
         consistent = consistent .and.
     &    ((na_feedback_choice(ispecies) .eq. 0 .or.
     &      na_feedback_option(ispecies) .eq. 0 .or.
     &      na_feedback_actuator(ispecies) .eq. 0) .or.
     &      ((na_feedback_actuator(ispecies).ne.1 .and.
     &       (na_feedback_option(ispecies).ne.4 .and.
     &        na_feedback_option(ispecies).ne.5 .and.
     &        na_feedback_option(ispecies).ne.8)) .or.
     &       na_feedback_actuator(ispecies).eq.1))
       end do
       call xertst (consistent,
     .  'feedback option not consistent with feedback actuator')
*   ..test if feedback maximum is greater than minimum
       do ispecies = 0, nspecies-1
         consistent = consistent .and.
     &    ((na_feedback_choice(ispecies) .eq. 0 .or.
     &      na_feedback_actuator(ispecies) .eq. 0 .or.
     &      na_feedback_actuator(ispecies) .eq. 0) .or.
     &      (na_feedback_actuator_max(ispecies).le.0.0_R8 .or.
     &       (na_feedback_actuator_max(ispecies) .gt.
     &        na_feedback_actuator_min(ispecies))))
       end do
       call xertst (consistent,
     .  'feedback actuator maximum <= minimum')
*   ..test if neutral pressure can be computed if it is required
       do ispecies = 0, nspecies-1
         if (na_feedback_choice(ispecies).ne.10) cycle
         call xertst (lpfrs_pmp.ne.0,
     &     'Neutral pressure cannot be computed and used for feedback')
       end do
*  ..now read data for other feedback schemes given directly in b2mn.dat
       call ipgetr ('b2stbc_fnaycore', fnaycore)
       call ipgetr ('b2stbc_fnaycore_alpha', fnaycore_alpha)
       call xertst (0.0_R8.le.fnaycore_alpha,
     1   'faulty internal parameter fnaycore_alpha')
       call ipgetr ('b2stbc_fheycore', fheycore)
       call xertst (0.0_R8.le.fheycore,
     1   'faulty internal parameter fheycore')
       call ipgetr ('b2stbc_fheycore_alpha', fheycore_alpha)
       call xertst (0.0_R8.le.fheycore_alpha,
     1   'faulty internal parameter fheycore_alpha')
       call ipgetr ('b2stbc_fhiycore', fhiycore)
       call xertst (0.0_R8.le.fhiycore,
     1   'faulty internal parameter fhiycore')
       call ipgetr ('b2stbc_fhiycore_alpha', fhiycore_alpha)
       call xertst (0.0_R8.le.fhiycore_alpha,
     1   'faulty internal parameter fhiycore_alpha')
       call ipgeti ('b2tsbc_fhiycore_kinetic_energy',
     1  fhiycore_kinetic_energy)
       call ipgetr ('b2stbc_fchycore', fchycore)
       call ipgetr ('b2stbc_fchycore_alpha', fchycore_alpha)
       call xertst (0.0_R8.le.fchycore_alpha,
     1   'faulty internal parameter fchycore_alpha')
       call ipgetr ('b2stbc_nesepm', nesepm)
       call xertst (0.0_R8.le.nesepm,
     1   'faulty internal parameter nesepm')
       call ipgetr ('b2stbc_nesepm_sol', nesepm_sol)
       call xertst (0.0_R8.le.nesepm_sol,
     1   'faulty internal parameter nesepm_sol')
       call ipgetr ('b2stbc_nepedm_sol', nepedm_sol)
       call xertst (0.0_R8.le.nepedm_sol,
     1   'faulty internal parameter nepedm_sol')
       call ipgetr ('b2stbc_nesepm_pfr', nesepm_pfr)
       call xertst (0.0_R8.le.nesepm_pfr,
     1   'faulty internal parameter nesepm_pfr')
       call ipgetr ('b2stbc_volrec_sol', volrec_sol)
       call xertst (0.0_R8.le.volrec_sol,
     1   'faulty internal parameter volrec_sol')
       call ipgetr ('b2stbc_nesepm_alpha', nesepm_alpha)
       call xertst (0.0_R8.le.nesepm_alpha,
     1   'faulty internal parameter nesepm_alpha')
       call ipgetr ('b2stbc_nesepm_beta', nesepm_beta)
       call ipgetr ('b2stbc_nesepm_gamma', nesepm_gamma)
       call xertst (0.0_R8.le.nesepm_gamma,
     1   'faulty internal parameter nesepm_gamma')
       call ipgetr ('b2stbc_nesepm_minpuff', nesepm_minpuff)
       call xertst (0.0_R8.le.nesepm_minpuff,
     1   'faulty internal parameter nesepm_minpuff')
       call ipgetr ('b2stbc_nesepm_maxpuff', nesepm_maxpuff)
       call xertst (0.0_R8.le.nesepm_maxpuff,
     1   'faulty internal parameter nesepm_maxpuff')
       call ipgetr ('b2stbc_nesepm_overshoot', nesepm_overshoot)
       call xertst (1.0_R8.le.nesepm_overshoot.or.
     1  0.0_R8.eq.nesepm_overshoot,
     1  'faulty internal parameter nesepm_overshoot')
       call ipgetr ('b2stbc_volrec_alpha', volrec_alpha)
       call xertst (0.0_R8.le.volrec_alpha,
     1   'faulty internal parameter volrec_alpha')
       call ipgetr ('b2stbc_volrec_beta', volrec_beta)
       call ipgetr ('b2stbc_volrec_overshoot', volrec_overshoot)
       call xertst (1.0_R8.le.volrec_overshoot.or.
     1  0.0_R8.eq.volrec_overshoot,
     1  'faulty internal parameter volrec_overshoot')
       if (nnreg(0).le.2) then
         coreregno = 0
       else
         coreregno = 1
       end if
       call ipgeti ('b2stbc_coreregno', coreregno)
       call ipgeti ('b2stbc_coreregn2', coreregn2)
       call ipgeti ('b2stbc_solregno', solregno)
       call ipgetr ('b2stbc_ndes_sol', ndes_sol)
       call xertst (0.0_R8.le.ndes_sol,
     .  'faulty internal parameter ndes_sol')
       if(fnaycore.ne.-1.0e30_R8.or.fheycore.gt.0.0_R8.or.
     &    fhiycore.gt.0.0_R8.or.fchycore.ne.-1.0e30_R8.or.
     2    nesepm.gt.0.0_R8) then
         call xertst (cbirso.le.coreregno.and.
     1          coreregno.lt.cbirso+cbnrso-1,
     2    'faulty internal parameter coreregno')
         if (nncut.ge.2)
     >    call xertst (cbirso.le.coreregn2.and.
     1    coreregn2.lt.cbirso+cbnrso-1,
     2    'faulty internal parameter coreregn2')
       endif
       if(nesepm_sol.gt.0.0_R8.or.volrec_sol.gt.0.0_R8.or.
     1    nepedm_sol.gt.0.0_R8.or.ndes_sol.gt.0.0_R8) then
         call xertst (cbirno.le.solregno.and.
     1    solregno.le.cbirno+cbnrno-1,
     2    'faulty internal parameter solregno')
         call xertst(
     .    (ndes_sol.gt.0.0_R8 .and. volrec_sol.eq.0.0_R8 .and.
     .     nesepm_sol.eq.0.0_R8 .and. nepedm_sol.eq.0.0_R8) .or.
     .    (ndes_sol.eq.0.0_R8 .and. volrec_sol.gt.0.0_R8 .and.
     .     nesepm_sol.eq.0.0_R8 .and. nepedm_sol.eq.0.0_R8) .or.
     .    (ndes_sol.eq.0.0_R8 .and. volrec_sol.eq.0.0_R8 .and.
     .     nesepm_sol.gt.0.0_R8 .and. nepedm_sol.eq.0.0_R8) .or.
     .    (ndes_sol.eq.0.0_R8 .and. volrec_sol.eq.0.0_R8 .and.
     .     nesepm_sol.eq.0.0_R8 .and. nepedm_sol.gt.0.0_R8),
     .    'Simultaneous SOL feedback schemes not allowed!')
       endif
       if(nesepm.gt.0.0_R8 .or. nesepm_sol.gt.0.0_R8 .or.
     &    nesepm_pfr.gt.0.0_R8 .or. nepedm_sol.gt.0.0_R8) then
          call get_jsep(nx,ny,jxi,jxa,jsep)
          iyped = jsep/2
       endif
       if(nepedm_sol.gt.0.0_R8) then
         call ipgeti ('b2stbc_iyped', iyped)
         call xertst (iyped.gt.-2 .and. iyped.le.ny,
     .    'faulty input iyped')
         call xertst (on_closed_surface(jxa,iyped),
     .    'location of top of pedestal outside core region!')
       endif
       call ipgeti ('b2stbc_ncallfeedback', ncallfeedback)
       call ipgetr ('b2stbc_ndes', ndes)
       call ipgeti('b2tfnb_mdf', mdf_fnb)
       call ipgeti('b2mndt_style', style_b2mndt)
       call xertst (0.0_R8.le.ndes, 'faulty internal parameter ndes')
       call ipgetr ('b2stbc_cbsnafac', cbsnafac)
       call xertst (0.0_R8.lt.cbsnafac,
     1   'faulty internal parameter cbsnafac')
       call ipgetr ('b2stbc_private_flux_puff', private_flux_puff)
       call xertst (0.0_R8.le.private_flux_puff,
     1  'faulty internal parameter private_flux_puff')
       call ipgeti ('b2stbc_pfrregno1', pfrregno1)
       call ipgeti ('b2stbc_pfrregno2', pfrregno2)
       call ipgeti ('b2stbc_isfeedback', isfeedback)
       call xertst (0.le.isfeedback.and.isfeedback.lt.ns,
     1  'faulty internal parameter isfeedback')
       if (ndes.gt.0.0_R8.or.private_flux_puff.gt.0.0_R8.or.
     &     nesepm_pfr.gt.0.0_R8) then
         call xertst (cbirso.le.pfrregno1.and.
     1    pfrregno1.le.cbirso+cbnrso-1,
     2    'faulty internal parameter pfrregno1')
         call xertst (cbirso.le.pfrregno2.and.
     1    pfrregno2.le.cbirso+cbnrso-1,
     2    'faulty internal parameter pfrregno2')
         if (pfrregno1.gt.cbirso) then
           ixl1 = max(0, nint(cbrbrk(pfrregno1-1)*nx))
         else
           ixl1 = 0
         endif
         ixr1 = min(nx, nint(cbrbrk(pfrregno1)*nx))-1
         if (pfrregno2.gt.cbirso) then
           ixl2 = max(0, nint(cbrbrk(pfrregno2-1)*nx))
         else
           ixl2 = 0
         endif
         ixr2 = min(nx, nint(cbrbrk(pfrregno2)*nx))-1
         call xertst (pfrregno1.eq.pfrregno2 .or.
     &    rightix(ixr1,-1).eq.ixl2 .or. rightix(ixr2,-1).eq.ixl1,
     .    'pfrregno1 and pfrregno2 are not connected!')
         call xertst(
     .    (ndes.gt.0.0_R8 .and. private_flux_puff.eq.0.0_R8 .and.
     .     nesepm_pfr.eq.0.0_R8) .or.
     .    (ndes.eq.0.0_R8 .and. private_flux_puff.gt.0.0_R8 .and.
     .     nesepm_pfr.eq.0.0_R8) .or.
     .    (ndes.eq.0.0_R8 .and. private_flux_puff.eq.0.0_R8 .and.
     .     nesepm_pfr.gt.0.0_R8),
     .    'Simultaneous PFR feedback schemes not allowed!')
       endif
       call ipgetr ('b2mndr_rescale_neutrals_sources',
     *       neutral_sources_rescale)
       call xertst (0.0_R8.le.neutral_sources_rescale,
     1   'faulty input neutral_sources_rescale')
       call ipgeti ('b2mndr_eirene', use_eirene)
       call ipgeti ('eirene_nesepm_istra', nesepm_istra)
       if(nesepm_istra.gt.0) then
         call xertst (nesepm_istra.le.nstrat,
     1     'faulty internal parameter nesepm_istra')
         call xertst (crcstra(nesepm_istra).eq.'C',
     1     'Stratum nesepm_istra is not declared as a puff stratum!')
       else
         do istrai = 1, nstrat
           if(nesepm_istra.le.0.and.CRCSTRA(istrai).EQ.'C') then
             nesepm_istra=istrai
           endif
         enddo
         if (ndes.gt.0.0_R8.or.nesepm_pfr.gt.0.0_R8.or.
     &    private_flux_puff.gt.0.0_R8.or.nesepm_sol.gt.0.0_R8.or.
     *    volrec_sol.gt.0.0_R8.or.ndes_sol.gt.0.0_R8.or.
     &    nepedm_sol.gt.0.0_R8) then
           icount = 0
           do istrai = 1, nstrat
             if (crcstra(istrai).eq.'C' .and.
     .           species_start(istrai).le.isfeedback .and.
     .           species_end(istrai).ge.isfeedback) icount = icount+1
           end do
           if (icount.eq.0) then
             message =
     .        'No puff stratum found for feedback scheme!'
           else if (icount.gt.1) then
             message =
     .        'Multiple definitions of puff strata, use nesepm_istra'
           end if
           call xertst (icount.eq.1, trim(message))
         endif
       endif
       if (ndes.gt.0.0_R8.or.nesepm_pfr.gt.0.0_R8.or.
     &    private_flux_puff.gt.0.0_R8.or.nesepm_sol.gt.0.0_R8.or.
     *    volrec_sol.gt.0.0_R8.or.ndes_sol.gt.0.0_R8.or.
     &    nepedm_sol.gt.0.0_R8) then
         call xertst (nesepm_istra.gt.0,
     &     'No "C" stratum for feedback scheme specified !')
       endif
*   ..test for input consistency, i.e. are the cbs.. arrays used in the b.c settings
       consistent_core = fnaycore.eq.-1.0e30_R8 .and. fheycore.eq.0.0_R8
     &             .and. fhiycore.eq.0.0_R8 .and. nesepm.eq.0.0_R8
     &             .and. fchycore.eq.-1.0e30_R8
       consistent_pfr1 = ndes.eq.0.0_R8 .and. nesepm_pfr.eq.0.0_R8
     &             .and. private_flux_puff.eq.0.0_R8
       consistent_pfr2 = consistent_pfr1 .or. pfrregno1.eq.pfrregno2
       consistent_sol = nesepm_sol.eq.0.0_R8 .and. nepedm_sol.eq.0.0_R8
     &             .and. ndes_sol.eq.0.0_R8 .and. volrec_sol.eq.0.0_R8
       do ib = 1, nbc
         if (bcchar(ib).eq.'N') then
           if(nesepm_sol.gt.0.0_R8 .or. nepedm_sol.gt.0.0_R8 .or.
     .        volrec_sol.gt.0.0_R8) then
             consistent_sol = consistent_sol .or.
     .        bccon(isfeedback,ib).eq.11
           end if
           if(ndes_sol.gt.0.0_R8) then
             consistent_sol = consistent_sol .or.
     .        bccon(isfeedback,ib).eq.12
           end if
         else if (bcchar(ib).eq.'S') then
           if (core_boundary(ib)) then
             if (fnaycore.ne.-1.0e30_R8 .or. nesepm.gt.0.0_R8) then
               consistent_core = consistent_core .or.
     .          bccon(isfeedback,ib).eq.11
             end if
             if(fheycore.gt.0.0_R8) then
               consistent_core = consistent_core .or. bcene(ib).eq.10
             end if
             if(fhiycore.gt.0.0_R8) then
               consistent_core = consistent_core .or. bceni(ib).eq.10
             end if
             if(fchycore.ne.-1.0e30_R8) then
               consistent_core = consistent_core .or. bcpot(ib).eq.10
             end if
           else
             if(ndes.gt.0.0_R8) then
               if (mod(region(bc_list_x(1,ib),
     >                        bc_list_y(1,ib),0),4).eq.3) then
                 consistent_pfr1 = consistent_pfr1 .or.
     .            bccon(isfeedback,ib).eq.12
               else
                 consistent_pfr2 = consistent_pfr2 .or.
     .            bccon(isfeedback,ib).eq.12
               end if
             end if
             if(nesepm_pfr.gt.0.0_R8 .or. private_flux_puff.gt.0.0_R8)
     >        then
               if (mod(region(bc_list_x(1,ib),
     >                        bc_list_y(1,ib),0),4).eq.3) then
                 consistent_pfr1 = consistent_pfr1 .or.
     .            bccon(isfeedback,ib).eq.11
               else
                 consistent_pfr2 = consistent_pfr2 .or.
     .            bccon(isfeedback,ib).eq.11
               end if
             end if
           end if
         end if
       end do
       consistent = consistent_core .and. consistent_sol .and.
     &              consistent_pfr1 .and. consistent_pfr2
       consistent = consistent .or. lbndusr .or. b2sral_style.eq.1 .or.
     &  boundary_namelist.lt.1
       call xertst (consistent,
     .  'feedback scheme not consistent with B.C.''s!')
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
*    ..test sign of ni, ne, te, ti
       call b2xvsg (n2, ni(-1,-1,0), 1, 'ni0', '.gt.')
       call b2xvsg (n2, ni(-1,-1,1), 1, 'ni1', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
      endif

      if(maxval(na_feedback_target(0:nspecies-1)).gt.0.0_R8 .or.
     &   maxval(na_feedback_choice(0:nspecies-1)).gt.0) then
! we have a request for feedback control of na
! the idea here is to break the actions into 3 phases
! 1/ calculate something based on the present state of the plasma depending on a
!    species block based on "na_feedback_choice"
! 2/ compare this with the reference value and calculate a change
! 3/ apply this in some way to the species block based on "na_feedback_actuator"

        do ispecies = 0, nspecies-1

          if(na_feedback_choice(ispecies) .gt. 0 .and.
     &       na_feedback_option(ispecies) .gt. 0 .and.
     &       na_feedback_actuator(ispecies) .gt. 0) then

            if (main_call) then
              na_feedback_prev(ispecies) = na_feedback_save(ispecies)
            end if
! phase 1:
            select case(na_feedback_choice(ispecies))

            case (0) ! no action

            case (1) ! local density of isonuclear sequence
                     ! (including contributions from Eirene species)

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'b2stbc_fb_pre_choice: using ',
     &           'local_density for sequence ',
     &            ispecies+1
              endif

              na_feedback_current = 0.0_R8
              vol_feedback = 0.0_R8
              is_start=-1
              is_end=-1
              do is = 0, ns-1
                if(b2espcr(is).eq.ispecies+1) then
                  if(is_start.eq.-1) is_start=is
                  do iy = na_feedback_iy1(ispecies),
     &                    na_feedback_iy2(ispecies)
                    do ix = na_feedback_ix1(ispecies),
     &                      na_feedback_ix2(ispecies)
                      if(region(ix,iy,0).ne.0 .and.
     &                   leftix(ix,iy).ne.-2 .and.
     &                   rightix(ix,iy).ne.nx+1 .and.
     &                   bottomiy(ix,iy).ne.-2 .and.
     &                   topiy(ix,iy).ne.ny+1) then
                        na_feedback_current =
     &                   na_feedback_current +
     &                   na(ix,iy,is) * vol(ix, iy)
                        if(is.eq.is_start) then
                          vol_feedback = vol_feedback +
     &                                   vol(ix, iy)
                        endif
                      endif
                    enddo
                  enddo
                  is_end=is
                endif
              enddo
#ifdef B25_EIRENE
              do j=1,nnatmi !{
                if (latmscl(j).eq.ispecies+1) then !{
                  do iy = na_feedback_iy1(ispecies),
     &                    na_feedback_iy2(ispecies) !{
                    do ix = na_feedback_ix1(ispecies),
     &                      na_feedback_ix2(ispecies) !{
                      if(region(ix,iy,0).ne.0 .and.
     &                   leftix(ix,iy).ne.-2 .and.
     &                   rightix(ix,iy).ne.nx+1 .and.
     &                   bottomiy(ix,iy).ne.-2 .and.
     &                   topiy(ix,iy).ne.ny+1) then !{
                        na_feedback_current = na_feedback_current +
     &                   dab2(ix+1,iy+1,j,1) * vol(ix,iy)
                        do k=1,nnmoli !{
                          if (mlcmp(j,k).eq.0) cycle
                          na_feedback_current =
     &                     na_feedback_current + mlcmp(j,k) *
     &                     dmb2(ix+1,iy+1,k,1) * vol(ix,iy)
                        end do !}
                        do k=1,nnioni !{
                          if (micmp(j,k).eq.0) cycle
                          na_feedback_current =
     &                     na_feedback_current + micmp(j,k) *
     &                      dib2(ix+1,iy+1,k,1) * vol(ix,iy)
                        end do !}
                      end if !}
                    end do !}
                  end do !}
                end if !}
              end do !}
#endif
              na_feedback_current = na_feedback_current /
     &                 vol_feedback

            case (2) ! local electron density

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'b2stbc_fb_pre_choice: using ',
     &            'local_electron_density for sequence ',
     &             ispecies+1
              endif

              na_feedback_current = 0.0_R8
              vol_feedback = 0.0_R8
              is_start=-1
              is_end=-1
              do is = 0, ns-1
                 if(b2espcr(is).eq.ispecies+1) then
                   if(is_start.eq.-1) is_start=is
                   is_end=is
                 endif
              enddo
              do iy = na_feedback_iy1(ispecies),
     &                na_feedback_iy2(ispecies)
                do ix = na_feedback_ix1(ispecies),
     &                  na_feedback_ix2(ispecies)
                  if(region(ix,iy,0).ne.0 .and.
     &               leftix(ix,iy).ne.-2 .and.
     &               rightix(ix,iy).ne.nx+1 .and.
     &               bottomiy(ix,iy).ne.-2 .and.
     &               topiy(ix,iy).ne.ny+1) then
                    na_feedback_current = na_feedback_current +
     &               ne(ix, iy) * vol(ix, iy)
                    vol_feedback = vol_feedback + vol(ix, iy)
                 endif
                enddo
              enddo
              na_feedback_current = na_feedback_current /
     &                 vol_feedback

            case (3) ! nesepm electron density

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'b2stbc_fb_pre_choice: using ',
     &           'nesepm_electron_density for sequence ',
     &            ispecies+1
              endif

              na_feedback_current = 0.5_R8 *
     *          (ne(jxa,jsep)+
     *           ne(topix(jxa,jsep),topiy(jxa,jsep)))
              is_start=-1
              is_end=-1
              do is = 0, ns-1
                if(b2espcr(is).eq.ispecies+1) then
                  if(is_start.eq.-1) is_start=is
                  is_end=is
                endif
              enddo

            case (4) ! total particle content for that isonuclear sequence
                     ! (including contributions from Eirene species)

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'b2stbc_fb_pre_choice: using ',
     &            'total particle content for sequence ',
     &             ispecies+1
              endif
              na_feedback_current = 0.0_R8
              is_start=-1
              is_end=-1
              do is = 0, ns-1
                if(b2espcr(is).eq.ispecies+1) then
                  if(is_start.eq.-1) is_start=is
                  do iy = na_feedback_iy1(ispecies),
     &                    na_feedback_iy2(ispecies)
                    do ix = na_feedback_ix1(ispecies),
     &                      na_feedback_ix2(ispecies)
                      if(region(ix,iy,0).ne.0 .and.
     &                   leftix(ix,iy).ne.-2 .and.
     &                   rightix(ix,iy).ne.nx+1 .and.
     &                   bottomiy(ix,iy).ne.-2 .and.
     &                   topiy(ix,iy).ne.ny+1) then
                        na_feedback_current =
     &                    na_feedback_current +
     &                    na(ix,iy,is) * vol(ix, iy)
                      endif
                    enddo
                  enddo
                  is_end=is
                endif
              enddo
#ifdef B25_EIRENE
              do j=1,nnatmi !{
                if (latmscl(j).eq.ispecies+1) then !{
                  do iy = na_feedback_iy1(ispecies),
     &                    na_feedback_iy2(ispecies) !{
                    do ix = na_feedback_ix1(ispecies),
     &                      na_feedback_ix2(ispecies) !{
                      if(region(ix,iy,0).ne.0 .and.
     &                   leftix(ix,iy).ne.-2 .and.
     &                   rightix(ix,iy).ne.nx+1 .and.
     &                   bottomiy(ix,iy).ne.-2 .and.
     &                   topiy(ix,iy).ne.ny+1) then !{
                        na_feedback_current = na_feedback_current +
     &                   dab2(ix+1,iy+1,j,1) * vol(ix,iy)
                        do k=1,nnmoli !{
                          if (mlcmp(j,k).eq.0) cycle
                          na_feedback_current =
     &                     na_feedback_current + mlcmp(j,k) *
     &                     dmb2(ix+1,iy+1,k,1) * vol(ix,iy)
                        end do !}
                        do k=1,nnioni !{
                          if (micmp(j,k).eq.0) cycle
                          na_feedback_current =
     &                     na_feedback_current + micmp(j,k) *
     &                     dib2(ix+1,iy+1,k,1) * vol(ix,iy)
                        end do !}
                      end if !}
                    end do !}
                  end do !}
                end if !}
              end do !}
#endif

            case (5) ! total ion content for that isonuclear sequence
                     ! (not including neutrals)

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'b2stbc_fb_pre_choice: using ',
     &            'total ion content for sequence ',
     &             ispecies+1
              endif
              na_feedback_current = 0.0_R8
              is_start=-1
              is_end=-1
              do is = 0, ns-1
                if(is_neutral(is)) cycle
                if(b2espcr(is).eq.ispecies+1) then
                  if(is_start.eq.-1) is_start=is
                  do iy = na_feedback_iy1(ispecies),
     &                    na_feedback_iy2(ispecies)
                    do ix = na_feedback_ix1(ispecies),
     &                      na_feedback_ix2(ispecies)
                      if(region(ix,iy,0).ne.0 .and.
     &                   leftix(ix,iy).ne.-2 .and.
     &                   rightix(ix,iy).ne.nx+1 .and.
     &                   bottomiy(ix,iy).ne.-2 .and.
     &                   topiy(ix,iy).ne.ny+1) then
                        na_feedback_current =
     &                    na_feedback_current +
     &                    na(ix,iy,is) * vol(ix, iy)
                      endif
                    enddo
                  enddo
                  is_end=is
                endif
              enddo

            case (6) ! neutral particle flux through the core boundary

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'b2stbc_fb_pre_choice: using ',
     &            'core boundary flux for sequence ',
     &             ispecies+1
              endif
              na_feedback_current = na_feedback_const(ispecies)

#ifdef B25_EIRENE
              do k=1,ncns !{
                l=nlim+lcns(k)
                do j=1,nnatmi !{
                  if (latmscl(j).eq.ispecies+1) then !{
                    na_feedback_current=na_feedback_current
     .                +wldna(l,j,0)-wldra(l,j,0)-wldpa(l,j,0)
                    do i=1,nnmoli !{
                      na_feedback_current=na_feedback_current
     .             +mlcmp(j,i)*(wldnm(l,i,0)-wldrm(l,i,0)-wldpm(l,i,0))
                    end do !}
                  end if !}
                end do !}
              end do !}
#endif
              is_start = eb2spcr(ispecies+1)
              is_end = eb2spcr(ispecies+1)+nfluids(ispecies+1)
              allocate(charge_frac(1:is_end-is_start+1))
              charge_frac = 0.0_R8
              sum_s = 0.0_R8
              do is = is_start, is_end
                iss = is-is_start+1
                if (is_neutral(is).and.use_eirene.ne.0) cycle
                do iy = na_feedback_iy1(ispecies),
     .                  na_feedback_iy2(ispecies)
                  do ix = na_feedback_ix1(ispecies),
     .                    na_feedback_ix2(ispecies)
                    if (on_closed_surface(ix,iy)) then
                      sum_s = sum_s + na(ix,iy,is) * vol(ix,iy)
                      charge_frac(iss) = charge_frac(iss) +
     .                  na(ix,iy,is) * vol(ix,iy)
                    end if
                  end do
                end do
              end do
              charge_frac = charge_frac / sum_s

            case (7) ! relative average concentration of this species at the separatrix

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'b2stbc_fb_pre_choice: using ',
     &            'average separatrix concentration for sequence ',
     &             ispecies+1
              endif
              na_feedback_current = 0.0_R8
              sum_e = 0.0_R8
              is_start = -1
              is_end = -1
              do is = 0, ns-1
                if(b2espcr(is).eq.ispecies+1) then
                  if(is_start.eq.-1) is_start=is
                  if(is_neutral(is).and.use_eirene.ne.0) cycle
                  do ix = 0, nx-1
                    if(on_closed_surface(ix,jsep)) then
                      na_feedback_current =
     &                  na_feedback_current + (na(ix,jsep,is) +
     &                  na(topix(ix,jsep),topiy(ix,jsep),is)) *
     &                  gs(topix(ix,jsep),topiy(ix,jsep),1)
                      if (is.eq.is_start+1) then
                        sum_e = sum_e + (ne(ix,jsep) +
     &                    ne(topix(ix,jsep),topiy(ix,jsep))) *
     &                    gs(topix(ix,jsep),topiy(ix,jsep),1)
                      end if
                    endif
                  enddo
                  is_end = is
                endif
              enddo
              na_feedback_current = na_feedback_current/sum_e

            case (8) ! relative average concentration of this species on desired domain

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'b2stbc_fb_pre_choice: using average ',
     &            'concentration for sequence ', ispecies+1
              endif
              na_feedback_current = 0.0_R8
              sum_e = 0.0_R8
              is_start=-1
              is_end=-1
              do is = 0, ns-1
                if(b2espcr(is).eq.ispecies+1) then
                  if(is_start.eq.-1) is_start=is
                  if(is_neutral(is).and.use_eirene.ne.0) cycle
                  do iy = na_feedback_iy1(ispecies),
     &                    na_feedback_iy2(ispecies)
                    do ix = na_feedback_ix1(ispecies),
     &                      na_feedback_ix2(ispecies)
                      na_feedback_current =
     &                  na_feedback_current +
     &                  na(ix,iy,is) * vol(ix,iy)
                      if (is.eq.is_start+1) then
                        sum_e = sum_e +
     &                    ne(ix,iy) * vol(ix,iy)
                      end if
                    enddo
                  enddo
                  is_end=is
                endif
              enddo
              na_feedback_current = na_feedback_current/sum_e

            case (9) ! inner midplane electron density

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'b2stbc_fb_pre_choice: using ',
     &            'inner_midplance_electron_density'//
     &            'for sequence ',
     &             ispecies+1
              endif

              na_feedback_current = 0.5_R8 *
     *                 (ne(jxi,jsep)+
     *                 ne(topix(jxi,jsep),topiy(jxi,jsep)))
              is_start=-1
              is_end=-1
              do is = 0, ns-1
                if(b2espcr(is).eq.ispecies+1) then
                   if(is_start.eq.-1) is_start=is
                   is_end=is
                endif
              enddo
              
            case (10)  ! neutral pressure at PFR southern boundary (by A. Kukushkin)
            
              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'b2stbc_fb_pre_choice: using ',
     &            'neutral pressure for sequence ',
     &             ispecies+1
                call get_nxt(nx, nxtl, nxtr)
                call calc_pneut_init(nx, ny, nxtl, nxtr)
              endif
              
              call calc_pneut(ns, ny)
              
              if (Q_pneut.ne.Q_pneut) then
                  na_feedback_current = na_feedback_target(ispecies)
              else
                  na_feedback_current = Q_pneut
              endif
              
              is_start=-1
              is_end=-1
              do is = 0, ns-1
                if(b2espcr(is).eq.ispecies+1) then
                  if(is_start.eq.-1) is_start=is
                  is_end=is
                endif
              enddo
               
            case (11) ! total particle content for the isonuclear sequence
                      ! (including contributions from Eirene species)
                      ! in the edge (i.e. Krasheninnikov&Kukushkin style)

              if(ncall.eq.0) then !{
                write(*,'(a,a,i3)')
     &            'b2stbc_fb_pre_choice: using ',
     &            'total particle content in the edge for sequence ',
     &             ispecies+1
              endif !}

              na_feedback_current = 0.0_R8
              is_start=-1
              is_end=-1
              do is = 0, ns-1 !{
                if(b2espcr(is).eq.ispecies+1) then !{
                  if(is_start.eq.-1) is_start=is
                    do iy = 0,ny-1 !{
                      do ix = 0,nx-1 !{
                        if(.not.ON_CLOSED_SURFACE(ix,iy))!{
     &                    na_feedback_current =
     &                     na_feedback_current +
     &                     na(ix,iy,is) * vol(ix, iy) !}
                      enddo !}
                    enddo !}
                  is_end=is
                endif !}
              enddo !}
#ifdef B25_EIRENE
              do j=1,nnatmi !{
                if (latmscl(j).eq.ispecies+1) then !{
                  do iy = 0,ny-1 !{
                    do ix = 0,nx-1 !{
                      if(.not.ON_CLOSED_SURFACE(ix,iy)) then !{
                        na_feedback_current = na_feedback_current +
     &                   dab2(ix+1,iy+1,j,1) * vol(ix,iy)
                        do k=1,nnmoli !{
                          if (mlcmp(j,k).eq.0) cycle
                            na_feedback_current =
     &                       na_feedback_current + mlcmp(j,k) *
     &                       dmb2(ix+1,iy+1,k,1) * vol(ix,iy)
                        end do !}
                        do k=1,nnioni !{
                          if (micmp(j,k).eq.0) cycle
                            na_feedback_current =
     &                       na_feedback_current + micmp(j,k) *
     &                       dib2(ix+1,iy+1,k,1) * vol(ix,iy)
                        end do !}
                      end if !}
                    end do !}
                  end do !}
                end if !}
              end do !}
#endif

            case (12)  ! total radiation on desired domain
            
              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'b2stbc_fb_pre_choice: using ',
     &            'total radiation in area for sequence ',
     &             ispecies+1
              endif

              na_feedback_current = 0.0_R8

              do iy = na_feedback_iy1(ispecies),
     &                    na_feedback_iy2(ispecies)
                do ix = na_feedback_ix1(ispecies),
     &                      na_feedback_ix2(ispecies)
                  if (leftix(ix,iy).ne.-2 .and.
     &                rightix(ix,iy).ne.nx+1 .and.
     &                bottomiy(ix,iy).ne.-2 .and.
     &                topiy(ix,iy).ne.ny+1) then
                    na_feedback_current = na_feedback_current +
     .               sum(rqbrm(ix,iy,0:ns-1)) + sum(rqrad(ix,iy,0:ns-1))
#ifdef B25_EIRENE
                    na_feedback_current = na_feedback_current -
     .                         sum(eneutrad(ix+1,iy+1,1:nnatmi,0)) -
     .                          sum(emolrad(ix+1,iy+1,1:nnmoli,0)) -
     .                           sum(eionrad(ix+1,iy+1,1:nnioni,0))
#endif
                  endif
                enddo
              enddo

              is_start=-1
              is_end=-1
              do is = 0, ns-1
                if(b2espcr(is).eq.ispecies+1) then
                  if(is_start.eq.-1) is_start=is
                  is_end=is
                endif
              enddo
               
            case (13)  ! total radiation outside the separatrix
            
              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'b2stbc_fb_pre_choice: using ',
     &            'total radiation outside the separatrix ',
     &             ispecies+1
              endif

              na_feedback_current = 0.0_R8

              do iy = 0, ny-1
                do ix = 0, nx-1
                  if (.not. ON_CLOSED_SURFACE(ix,iy)) then
                    na_feedback_current = na_feedback_current +
     .               sum(rqbrm(ix,iy,0:ns-1)) + sum(rqrad(ix,iy,0:ns-1))
#ifdef B25_EIRENE
                    na_feedback_current = na_feedback_current -
     .                         sum(eneutrad(ix+1,iy+1,1:nnatmi,0)) -
     .                          sum(emolrad(ix+1,iy+1,1:nnmoli,0)) -
     .                           sum(eionrad(ix+1,iy+1,1:nnioni,0))
#endif
                  endif  
                enddo
              enddo

              is_start=-1
              is_end=-1
              do is = 0, ns-1
                if(b2espcr(is).eq.ispecies+1) then
                  if(is_start.eq.-1) is_start=is
                  is_end=is
                endif
              enddo
              
            case default
                  
              write(*,*) 'na feedback_choice not coded ',
     &                 na_feedback_choice(ispecies)

            end select
! phase 2:

c           Safeguard for na_feedback_options with multiplicative rescaling:
            if (na_feedback_option(ispecies).ne.4 .and.
     &          na_feedback_option(ispecies).ne.5) then
              if (na_feedback_current.lt.eps) then
                write(*,'(a,a,i3,a,i3,a,es9.1,a,es9.1)') 
     .            'b2stbc_fb_pre: Warning!   ',
     .            'na_feedback_current for sequence ',ispecies+1,
     .            ' using multiplicative rescale na_feedback_option',
     .            na_feedback_option(ispecies),' is below ',eps,
     .            ' and will be set to ',eps
                na_feedback_current = eps
              endif
            endif

            select case(na_feedback_option(ispecies))

            case (0) ! no action

            case (1) ! rescale slowed by na_feedback_alpha

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'b2stbc_fb_pre_option: using ',
     &            'slowed_rescale for sequence ',
     &             ispecies+1
              endif

              na_feedback_rescale =
     &           (1.0_R8 + na_feedback_alpha(ispecies) *
     &            na_feedback_target(ispecies) /
     &            na_feedback_current) /
     &           (1.0_R8 + na_feedback_alpha(ispecies))

            case (2) ! pure rescale

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'b2stbc_fb_pre_option: using ',
     &           'rescale for sequence ',
     &            ispecies+1
              endif

              na_feedback_rescale = na_feedback_target(ispecies) /
     &           na_feedback_current

            case (3) ! rescale slowed by tanh_log
                     ! 2**(tanh(log(x)/b)*log(a)/log(2))

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'b2stbc_fb_pre_option: using ',
     &           'tanh_log_rescale for sequence ',
     &            ispecies+1
              endif

              na_feedback_rescale = 2.0_R8 **
     &           (tanh(log(na_feedback_target(ispecies) /
     &            na_feedback_current)/
     &            na_feedback_beta(ispecies))*
     &            log(na_feedback_alpha(ispecies))/log(2.0_R8))

            case (4) ! rescale done according to SOLPS4 formula:
                     ! Feedback target changes as:
                     ! N = A + B*(time-T)
                     ! The actuator (should be gas puff only) changes as
                     ! S --> S + D, where
                     ! D = A((N_tar - N_cur)/dt + (N_prev - N_cur))/dt_prev

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'b2stbc_fb_pre_option: using ',
     &           'SOLPS4 linear increase for sequence ',
     &            ispecies+1
                if (na_feedback_const(ispecies).le.0.0_R8)
     &            na_feedback_const(ispecies) =
     &                        na_feedback_target(ispecies)
                na_feedback_prev(ispecies) = na_feedback_current
              endif

              na_feedback_target(ispecies) =
     .          na_feedback_const(ispecies) +
     .          na_feedback_beta(ispecies)*
     .         (tim-na_feedback_time(ispecies))
              na_feedback_rescale =
     .          na_feedback_alpha(ispecies) *
     .        ((na_feedback_target(ispecies) -
     .          na_feedback_current)/dteff +
     .         (na_feedback_prev(ispecies) -
     .          na_feedback_current)/dt_prev)
              na_feedback_save(ispecies) = na_feedback_current

            case (5) ! rescale done according to SOLPS4 formula:
                     ! Feedback target changes as:
                     ! Nt = A*exp((time-T)*B)
                     ! The actuator (should be gas puff only) changes as
                     ! S --> S + D, where
                     ! D = A((N_tar - N_cur)/dt + (N_prev - N_cur))/dt_prev

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'b2stbc_fb_pre_option: using ',
     &           'SOLPS4 exponential decay for sequence ',
     &            ispecies+1
                if (na_feedback_const(ispecies).le.0.0_R8)
     &            na_feedback_const(ispecies) =
     &                        na_feedback_target(ispecies)
                na_feedback_prev(ispecies) = na_feedback_current
              endif

              na_feedback_target(ispecies) =
     .           na_feedback_const(ispecies) *
     .           exp(na_feedback_beta(ispecies)*
     .          (tim-na_feedback_time(ispecies)))
              na_feedback_rescale =
     .           na_feedback_alpha(ispecies) *
     .         ((na_feedback_target(ispecies) -
     .           na_feedback_current)/dteff +
     .          (na_feedback_prev(ispecies) -
     .           na_feedback_current)/dt_prev)
              na_feedback_save(ispecies) = na_feedback_current

            case (6) ! rescale slowed by na_feedback_alpha (SOLPS4 style)

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'b2stbc_fb_pre_option: using ',
     &            'linear rescale for sequence ',
     &             ispecies+1
              endif

              na_feedback_rescale =
     &           na_feedback_alpha(ispecies) * 
     .           na_feedback_target(ispecies) / na_feedback_current

            case (7) ! combination of pure rescaling and rescaling slowed by tanh(log(x))
                     ! linear rescaling if |na_feedback_rescale - 1.| >= na_feedback_gamma
                     ! otherwise tanh slowing is applied {2**(tanh(log(x)/b)*log(a)/log(2))}
                      
              if(ncall.eq.0) then
                write(*,'(a,a,a,i3)')
     &             'b2stbc_fb_pre_option: using ',
     &             'combined pure rescale and tanh_log_rescale ',
     &             'for sequence ',
     &              ispecies+1
                if((na_feedback_gamma(ispecies).lt.0.0_R8) .or.
     &             (na_feedback_gamma(ispecies).gt.1.0_R8))
     .            na_feedback_gamma = 0.01_R8
              endif

              na_feedback_rescale =
     &          (na_feedback_target(ispecies) /
     &          na_feedback_current)
              if(abs(na_feedback_rescale-1.0_R8).lt.
     &           na_feedback_gamma(ispecies))
     .          na_feedback_rescale = 2.0_R8 **
     .           (tanh(log(na_feedback_rescale)/
     .           na_feedback_beta(ispecies))*
     .           log(na_feedback_alpha(ispecies))/log(2.0_R8))

            case (8)   ! Same as (7) but scaled with the
                       ! (pumping - core source - external source)
                       ! intended to be used with na_feedback_actuator=1
                       ! currently only applicable to EIRENE coupled runs 
            
              if(ncall.eq.0) then
                write(*,'(a,a,a,a,i3)')
     &             'b2stbc_fb_pre_option: using ',
     &             'combined pure rescale and tanh_log_rescale, ',
     &             'scaled with pumping and seeding balance ',
     &             'for sequence ',
     &              ispecies+1
                if((na_feedback_gamma(ispecies).lt.0.0_R8) .or.
     &             (na_feedback_gamma(ispecies).gt.1.0_R8))
     .            na_feedback_gamma = 0.01_R8
                if (use_eirene.eq.0) then
                  write(message,'(a,i3,a)')
     &              'Feedback option ',na_feedback_choice(ispecies),
     &              ' cannot be applied to standalone B25 runs'
                  call xerrab( trim(message) )
                endif
              endif

#ifdef B25_EIRENE
              na_feedback_rescale =
     &          (na_feedback_target(ispecies) /
     &          na_feedback_current)
              if(abs(na_feedback_rescale-1.0_R8).lt.
     &           na_feedback_gamma(ispecies))
     .          na_feedback_rescale = 2.0_R8 **
     .           (tanh(log(na_feedback_rescale)/
     .           na_feedback_beta(ispecies))*
     .           log(na_feedback_alpha(ispecies))/log(2.0_R8))

              if (mod(ncall,eirene_mod) .eq. 0) then
                call calc_pump_and_flux(nx, ny, ispecies,
     &                mdf_fnb, style_b2mndt, sna_pump,
     &                flxnp, flxic, flxnc, flxex, flxsp)
                if (flxsp.gt.0.0_R8) then
                  write(message,'(a,i3,a,i3,a)')
     &              'Sputtering detected for sequence ',ispecies+1,
     &              ', feedback option ',na_feedback_option(ispecies),
     &              'cannot be applied.'
                  call xertst(flxsp.gt.0.0_R8, trim(message) )
                else
                  na_feedback_rescale = (sna_pump - flxnc - flxic - 
     .                                   flxex) * na_feedback_rescale
                endif
              else
                na_feedback_rescale =
     .            saved_na_feedback_actuator(ispecies)
              endif
#endif

            end select
               
            write(*,'(a,1x,1p,i3,7g15.6)')
     &       'na_feedback_rescale ', ispecies+1,
     &        na_feedback_current,
     &        na_feedback_target(ispecies),
     &        saved_na_feedback_actuator(ispecies),
     &        na_feedback_rescale,
     &        na_feedback_prev(ispecies),
     &        dteff, dt_prev

! phase 3:

            select case(na_feedback_actuator(ispecies))

            case (0) ! no action

            case (1) ! gas puff via boundary condition

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'b2stbc_fb_pre_actuator: using ',
     &           'gas_puff for sequence ',
     &           ispecies+1
              endif

              if(ncall.eq.0 .and. 
     &           saved_na_feedback_actuator(ispecies).eq.0.0_R8 .and.
     &           .not.(feedback_namelist_used .and. 
     &                 (na_feedback_option(ispecies).eq.4 .or.
     &                  na_feedback_option(ispecies).eq.5 .or.
     &                  na_feedback_option(ispecies).eq.8))) then
                if(na_feedback_ib(ispecies).gt.0) then
                  write(message,'(a,i3,a,i2,a)')
     .              'Initial gas puff strength.le.0, set by CONPAR(',
     .              is_start,',',na_feedback_ib(ispecies),',1) !'
                  call xertst(
     .              conpar(is_start, na_feedback_ib(ispecies), 1)
     .                .gt.0.0_R8, trim(message) )
                  saved_na_feedback_actuator(ispecies) =
     .             conpar(is_start, na_feedback_ib(ispecies), 1)
                else  ! (na_feedback_ib(ispecies).lt.0 at this point) 
                  write(message,'(2a,i3,a)')
     .              'Initial gas puff strength.le.0, set by ',
     .               'USERFLUXPARM(',-na_feedback_ib(ispecies),',1) !'
                  call xertst(
     .              userfluxparm(-na_feedback_ib(ispecies), 1)
     .                .gt.0.0_R8, trim(message) )
                  saved_na_feedback_actuator(ispecies) =
     .             userfluxparm(-na_feedback_ib(ispecies),1)
                end if
              end if

              if ((na_feedback_option(ispecies).eq.4) .or.
     &            (na_feedback_option(ispecies).eq.5)) then
                saved_na_feedback_actuator(ispecies) =
     .               saved_na_feedback_actuator(ispecies) +
     .               na_feedback_rescale
              elseif (na_feedback_option(ispecies).eq.8) then
                saved_na_feedback_actuator(ispecies) = 
     .                                   na_feedback_rescale
              else
                saved_na_feedback_actuator(ispecies) =
     .               saved_na_feedback_actuator(ispecies) *
     .               na_feedback_rescale
              end if

              if(na_feedback_actuator_min(ispecies) .ge. 0.0_R8)
     &          saved_na_feedback_actuator(ispecies) =
     &            max(na_feedback_actuator_min(ispecies),
     &            saved_na_feedback_actuator(ispecies))
              if(na_feedback_actuator_max(ispecies) .gt. 0.0_R8)
     &          saved_na_feedback_actuator(ispecies) =
     &            min(na_feedback_actuator_max(ispecies),
     &            saved_na_feedback_actuator(ispecies))

              if(na_feedback_overshoot(ispecies) .gt. 1.0_R8 .and.
     &           (na_feedback_current .gt.
     &            na_feedback_target(ispecies) *
     &            na_feedback_overshoot(ispecies))) 
     &          saved_na_feedback_actuator(ispecies) = 0.0_R8

              if (use_eirene.ne.0) then
                userfluxparm(puff_strata(ispecies),1) =
     &            saved_na_feedback_actuator(ispecies)
              else
                conpar(is_start,na_feedback_ib(ispecies),1) =
     &            saved_na_feedback_actuator(ispecies)
              endif

            case (2) ! rescale na

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'b2stbc_fb_pre_actuator: using ',
     &           'rescale for sequence ',
     &           ispecies+1
              endif

              do is = is_start, is_end
                na(:,:,is) = na(:,:,is) * na_feedback_rescale
                kinrgy(:,:,is) =
     &           kinrgy(:,:,is) * na_feedback_rescale
                ne(:,:)=ne(:,:)+(na_feedback_rescale-1.0_R8)*
     &           na(:,:,is)*rza(:,:,is)
                fna(:,:,0:1,is) =
     &           fna(:,:,0:1,is)*na_feedback_rescale
                fna_mdf(:,:,0:1,is) =
     &           fna_mdf(:,:,0:1,is)*na_feedback_rescale
                fna_fcor(:,:,0:1,is) =
     &           fna_fcor(:,:,0:1,is)*na_feedback_rescale
                fna_nodrift(:,:,0:1,is) =
     &           fna_nodrift(:,:,0:1,is)*na_feedback_rescale
                fna_he(:,:,0:1,is) =
     &           fna_he(:,:,0:1,is)*na_feedback_rescale
                fnaPSch(:,:,0:1,is) =
     &           fnaPSch(:,:,0:1,is)*na_feedback_rescale
                fna_eir(:,:,0:1,is) =
     &           fna_eir(:,:,0:1,is)*na_feedback_rescale
                fna_32(:,:,0:1,is) =
     &           fna_32(:,:,0:1,is)*na_feedback_rescale
                fna_52(:,:,0:1,is) =
     &           fna_52(:,:,0:1,is)*na_feedback_rescale
              end do
*   ..initialise previous state and fluxes
              call b2scopy (n2*ns, na, 1, na0, 1)
              call b2scopy (n2, ne, 1, ne0, 1)
              call b2scopy (n2*2*ns, fna, 1, fna0, 1)
              call b2scopy (n2*ns, kinrgy, 1, kinrgy0, 1)
              call b2scopy (n2*ns, na, 1, nap, 1)
              call b2scopy (n2, ne, 1, nep, 1)
              call b2scopy (n2*2*ns, fna, 1, fnap, 1)
              call b2scopy (n2*ns, kinrgy, 1, kinrgyp, 1)
*   ..compute initial ni
              call b2xpni (nx, ny, ns, na, ni)
*   ..compute ni0
              call b2xpni (nx, ny, ns, na0, ni0)

            case (3) ! core boundary flux condition

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'b2stbc_fb_pre_actuator: using ',
     &           'core boundary flux condition for sequence ',
     &            ispecies+1
              endif

              is_start = eb2spcr(ispecies+1)+1
              is_end = eb2spcr(ispecies+1)+nfluids(ispecies+1)

              if(feedback_namelist_used .and.
     &             ncall .eq. 0 .and.
     &             saved_na_feedback_actuator(ispecies) .eq. 0.0_R8)
     &             saved_na_feedback_actuator(ispecies) =
     &             conpar(is_end, na_feedback_ib(ispecies), 1)

              saved_na_feedback_actuator(ispecies) =
     &             (1.0_R8 - na_feedback_alpha(ispecies)) *
     &             saved_na_feedback_actuator(ispecies) +
     &             na_feedback_rescale

              if(na_feedback_actuator_min(ispecies) .ge. 0.0_R8)
     &             saved_na_feedback_actuator(ispecies) =
     &             max(na_feedback_actuator_min(ispecies),
     &             saved_na_feedback_actuator(ispecies))
              if(na_feedback_actuator_max(ispecies) .gt. 0.0_R8)
     &             saved_na_feedback_actuator(ispecies) =
     &             min(na_feedback_actuator_max(ispecies),
     &             saved_na_feedback_actuator(ispecies))

              if(na_feedback_overshoot(ispecies) .gt. 1.0_R8 .and.
     &             na_feedback_current .gt.
     &             na_feedback_target(ispecies) *
     &             na_feedback_overshoot(ispecies)) then
                 total_flux = 0.0_R8
              else
                 total_flux = saved_na_feedback_actuator(ispecies)
              endif

              do is = is_start, is_end
                iss = is-is_start+1
                conpar(is, na_feedback_ib(ispecies), 1) =
     &                total_flux * charge_frac(iss)
              end do
                
            case (4)  ! rescale absorption coefficient of surfaces that
                      ! belong to a special surface group = na_feedback_ib(ispecies)
                      ! in this case saved_na_feedback_actuator is a new
                      ! absorption value (1 - recyct)
            
              if(ncall.eq.0) then
                  write(*,'(a,a,i3)')
     &                    'b2stbc_fb_pre_actuator: using ',
     &                    'specie absorption for sequence',
     &                    ispecies+1
                if (use_eirene.eq.0) then
                  write(message,'(a,i3,a)')
     &              'Feedback actuator ',na_feedback_actuator(ispecies),
     &              ' cannot be applied to standalone B25 runs'
                  call xerrab( trim(message) )
                endif
#ifdef B25_EIRENE
                if (saved_na_feedback_actuator(ispecies).le.0.0_R8 .or.
     &             saved_na_feedback_actuator(ispecies).gt.1.0_R8) then
                  write(*,'(a,a,a,i3,a)')
     .             'Warning!  ',
     .             'Saved_na_feedback_actuator has to be .ge.0  ',
     .             'and .le.1 for na_feedback_actuator = ',
     .             na_feedback_actuator(ispecies),
     .             ', default value will be used instead' 
                  if(sps_absr(na_feedback_ib(ispecies)).gt.0.0_R8 .and.
     &              sps_absr(na_feedback_ib(ispecies)).le.1.0_R8) then
                    saved_na_feedback_actuator(ispecies) =
     .                             sps_absr(na_feedback_ib(ispecies))
                  else
                    saved_na_feedback_actuator(ispecies) = 0.5_R8
                  endif
                endif
                absr_save(ispecies) =
     .                    saved_na_feedback_actuator(ispecies)
#endif
              endif

#ifdef B25_EIRENE
              if (ncall.eq.1) then
                call get_pump_list(ispecies+1, na_feedback_ib(ispecies),
     .                      pump_list(1:size(pump_list,DIM=1),ispecies))
              endif
              call calc_pump_and_flux(nx, ny, ispecies,
     &              mdf_fnb, style_b2mndt, sna_pump,
     &              flxnp, flxic, flxnc, flxex, flxsp)
              if (flxsp.ne.0.0_R8) then
                write(message,'(a,i3,a,i3,a)')
     &            'Sputtering detected for sequence: ',ispecies+1,
     &            ', feedback actuator ',na_feedback_actuator(ispecies),
     &            ' cannot be applied.'
                call xertst(flxsp.gt.0.0_R8, trim(message) )
              else
                if (mod(ncall,eirene_mod).eq.0 .and. ncall .gt. 0) then
                  if(abs(na_feedback_rescale*(sna_pump - flxnc)).gt.eps)
     &            then
                    na_feedback_rescale = (flxnp + flxic + flxex) / 
     .                         (sna_pump - flxnc) / na_feedback_rescale
                  else
                    na_feedback_rescale = (flxnp + flxic + flxex) / eps
                  endif
                  absr = saved_na_feedback_actuator(ispecies) *
     .                   na_feedback_rescale
                  if (absr.le.0.0_R8) then
                    if (saved_na_feedback_actuator(ispecies).gt.0.0_R8)
     .                absr_save(ispecies) = 
     .                          saved_na_feedback_actuator(ispecies)
                    if (na_feedback_rescale.gt.1.0_R8) then
                      saved_na_feedback_actuator(ispecies) = 
     .                                           absr_save(ispecies)
                    else
                      saved_na_feedback_actuator(ispecies) = 0.0_R8
                    end if
                  else if (absr.ge.1.0_R8) then
                    saved_na_feedback_actuator(ispecies) = 1.0_R8
                  else
                    saved_na_feedback_actuator(ispecies) = absr
                  end if
                  call eirene_recyct_update(ispecies+1,
     &                  saved_na_feedback_actuator(ispecies),
     &                  pump_list(1:size(pump_list,DIM=1),ispecies))
                endif
              endif
#endif
            case default
    
              write(*,*) 'na feedback_actuator not coded ',
     &                 na_feedback_actuator(ispecies)
     
            end select

            if (na_feedback_actuator(ispecies).eq.1 .or.
     &          na_feedback_actuator(ispecies).eq.3 .or.
     &              na_feedback_actuator(ispecies).eq.4)
     >        write(*,'(a,1p,i3,g15.6)') 'na_feedback_actuator ',
     &              ispecies+1, saved_na_feedback_actuator(ispecies)

            if (na_feedback_choice(ispecies).eq.6)
     >          deallocate(charge_frac)

          endif  ! na_feedback_choice > 0 .and. na_feedback_option > 0 .and. na_feedback_actuator > 0

        enddo ! ispecies

        feedback_namelist_used = feedback_namelist_used .or.
     &   na_feedback_actuator(ispecies).eq.1 .or.
     &   na_feedback_actuator(ispecies).eq.3
        
        dt_prev = dteff

      endif ! maxval(na_feedback_target) > 0.0 .and. maxval(na_feedback_choice) > 0

* DPC special treatment of core boundary condition
      if(fnaycore.ne.-1.0e30_R8) then
        fnaycurr=0.0_R8
        do ib = 1, nbc
          if(core_boundary(ib)) then
            do ix = 1, bc_list_size(ib)
              fnaycurr=fnaycurr+
     1         fna(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     2             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,
     3             isfeedback)
            end do
          endif
        enddo
        if(ncall.eq.0) then
          if(saved_cbsna_core.ne.0.0_R8) then
            cbsna(0,isfeedback,coreregno)=saved_cbsna_core
            if(nncut.ge.2)
     1       cbsna(0,isfeedback,coreregn2)=saved_cbsna_core
          endif
        endif
        cbsna(0,isfeedback,coreregno)=cbsna(0,isfeedback,coreregno)*
     *   (1.0_R8+fnaycore_alpha*
     1   (fnaycore-fnaycurr)/(abs(fnaycore)+abs(fnaycurr)+1.0e16_R8))
        if(nncut.ge.2)
     1   cbsna(0,isfeedback,coreregn2)=cbsna(0,isfeedback,coreregn2)*
     *    (1.0_R8+fnaycore_alpha*
     1    (fnaycore-fnaycurr)/(abs(fnaycore)+abs(fnaycurr)+1.0e16_R8))
        saved_cbsna_core=cbsna(0,isfeedback,coreregno)
        if(diagno.ge.1) then
          if(nncut.lt.2)
     >     write(*,*) 'fnaycurr= ',fnaycurr,
     1     ' ; target fnaycore= ',fnaycore,
     *     ' ; rescaling cbsna to ',cbsna(0,isfeedback,coreregno)
          if(nncut.ge.2)
     >     write(*,*) 'fnaycurr= ',fnaycurr,
     1     ' ; target fnaycore= ',fnaycore,
     *     ' ; rescaling cbsna to ',cbsna(0,isfeedback,coreregno),
     *     cbsna(0,isfeedback,coreregn2)
        endif
      endif

      if(fheycore.gt.0.0_R8) then
#ifdef COMPLICATED
        teave1=0.0_R8
        neave1=0.0_R8
        teave2=0.0_R8
        neave2=0.0_R8
#endif
        fheycurr=0.0_R8
        do ib = 1, nbc
          if(core_boundary(ib)) then
            do ix = 1, bc_list_size(ib)
#ifdef COMPLICATED
              teave1=teave1+te(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .                      ne(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .                     vol(bc_list_x(ix,ib),bc_list_y(ix,ib))
              neave1=neave1+ne(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .                     vol(bc_list_x(ix,ib),bc_list_y(ix,ib))
              teave2=teave2+te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                         topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .                     *ne(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                         topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .                    *vol(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                         topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
              neave2=neave2+ne(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                         topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .                    *vol(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                         topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
#endif
              fheycurr=fheycurr+
     .                 fhe(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                     topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)
            enddo
          endif
        enddo
#ifdef COMPLICATED
        teave1=teave1/neave1
        teave2=teave2/neave2
        tenext=teave1+abs(teave1-teave2)*
     1   (fheycore-max(0.0_R8,fheycurr))/fheycore*0.5_R8
        if(ncall.lt.ncallfeedback.and.main_call.and.diagno.ge.1) then
          write(*,*) ncall,': CORETEFEEDBACK (-)',fheycore,fheycurr,
     1     teave1/ev,teave2/ev,tenext/ev
        else
          if(diagno.ge.1)
     >     write(*,*) ncall,': CORETEFEEDBACK (+)',fheycore,fheycurr,
     1     teave1/ev,teave2/ev,tenext/ev
          call xertst (0.lt.tenext,
     1     'faulty feedback for core electron temperature')
          do i=0,7
            cbshe(i,coreregno)=0.0_R8
            if(nncut.ge.2) cbshe(i,coreregn2)=0.0_R8
          enddo
          cbshe(0,coreregno)=tenext*1.0e30_R8
          if(nncut.ge.2) cbshe(0,coreregn2)=tenext*1.0e30_R8
          cbshe(1,coreregno)=-1.0e30_R8
          if(nncut.ge.2) cbshe(1,coreregn2)=-1.0e30_R8
        endif
#else
        if(ncall.eq.0) then
          if(saved_cbshe_core.ne.0.0_R8) then
            cbshe(0,coreregno)=saved_cbshe_core
            if(nncut.ge.2) cbshe(0,coreregn2)=saved_cbshe_core
          endif
        endif
        cbshe(0,coreregno)=cbshe(0,coreregno)*
     *       (1.0_R8+fheycore_alpha*fheycore/fheycurr)/
     *       (1.0_R8+fheycore_alpha)
        if(nncut.ge.2) cbshe(0,coreregn2)=cbshe(0,coreregn2)*
     *       (1.0_R8+fheycore_alpha*fheycore/fheycurr)/
     *       (1.0_R8+fheycore_alpha)
        saved_cbshe_core=cbshe(0,coreregno)
        if(diagno.ge.1)
     >   write(*,*) 'Rescaling cbshe to ',cbshe(0,coreregno)
#endif
      endif

      if(fhiycore.gt.0.0_R8) then
#ifdef COMPLICATED
        tiave1=0.0_R8
        niave1=0.0_R8
        tiave2=0.0_R8
        niave2=0.0_R8
#endif
        fhiycurr=0.0_R8
        do ib=1,nbc
          if(core_boundary(ib)) then
            do ix=1,bc_list_size(ib)
#ifdef COMPLICATED
              tiave1=tiave1+ti(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .                      ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)*
     .                     vol(bc_list_x(ix,ib),bc_list_y(ix,ib))
              niave1=niave1+ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)*
     .                     vol(bc_list_x(ix,ib),bc_list_y(ix,ib))
              tiave2=tiave2+
     .               ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                  topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .              *ni(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                  topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),0)
     .             *vol(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                  topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
              niave2=niave2+
     .               ni(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                  topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),0)*
     .              vol(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                  topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
#endif
              fhiycurr=fhiycurr+
     .                 fhi(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                     topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)
              if (fhiycore_kinetic_energy.eq.1) then
                do is = 0, ns-1
                  fhiycurr=fhiycurr+(1.0_R8-BoRiS)*
     .             fhm(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,is)
                enddo
              endif
            enddo
          endif
        enddo
#ifdef COMPLICATED
        tiave1=tiave1/niave1
        tiave2=tiave2/niave2
        tinext=tiave1+abs(tiave1-tiave2)*
     1         (fhiycore-max(0.0_R8,fhiycurr))/fhiycore*0.5_R8
        if(ncall.lt.ncallfeedback.and.main_call.and.diagno.ge.1) then
          write(*,*) ncall,': CORETIFEEDBACK (-)',fhiycore,fhiycurr,
     1     tiave1/ev,tiave2/ev,tinext/ev
        else
          if(diagno.ge.1)
     >     write(*,*) ncall,': CORETIFEEDBACK (+)',fhiycore,fhiycurr,
     1     tiave1/ev,tiave2/ev,tinext/ev
          call xertst (0.0_R8.lt.tinext,
     1           'faulty feedback for core ion temperature')
          do is=0,ns-1
            do i=0,7
              cbshi(i,is,coreregno)=0.0_R8
              if(nncut.ge.2) cbshi(i,is,coreregn2)=0.0_R8
            enddo
            cbshi(0,is,coreregno)=tinext*1.0e30_R8
            if(nncut.ge.2) cbshi(0,is,coreregn2)=tinext*1.0e30_R8
            cbshi(1,is,coreregno)=-1.0e30_R8
            if(nncut.ge.2) cbshi(1,is,coreregn2)=-1.0e30_R8
          enddo
        endif
#else
        if(ncall.eq.0) then
          if(saved_cbshi_core.ne.0.0_R8) then
            cbshi(0,0:ns-1,coreregno)=saved_cbshi_core
            if(nncut.ge.2) cbshi(0,0:ns-1,coreregn2)=saved_cbshi_core
          endif
        endif
        cbshi(0,0:ns-1,coreregno)=
     *       cbshi(0,0:ns-1,coreregno)*
     *       (1.0_R8+fhiycore_alpha*fhiycore/fhiycurr)/
     *       (1.0_R8+fhiycore_alpha)
        if(nncut.ge.2) cbshi(0,0:ns-1,coreregn2)=
     *       cbshi(0,0:ns-1,coreregn2)*
     *       (1.0_R8+fhiycore_alpha*fhiycore/fhiycurr)/
     *       (1.0_R8+fhiycore_alpha)
        saved_cbshi_core=cbshi(0,0,coreregno)
        if(diagno.ge.1)
     >   write(*,*) 'Rescaling cbshi to ',cbshi(0,0:ns-1,coreregno)
#endif
      endif

cxpb current feedback not yet tested!
      if(fchycore.ne.-1.0e30_R8) then
        fchycurr=0.0_R8
        do ib=1,nbc
          if(core_boundary(ib)) then
            do ix = 1, bc_list_size(ib)
              fchycurr=fchycurr+
     .                 fch(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                     topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)
            enddo
          endif
        enddo
        if(ncall.eq.0) then
          if(saved_cbsch_core.ne.0.0_R8) then
            cbsch(0,coreregno)=saved_cbsch_core
            if(nncut.ge.2) cbsch(0,coreregn2)=saved_cbsch_core
          endif
        endif
        cbsch(0,coreregno)=
     *       cbsch(0,coreregno)*
     *       (1.0_R8+fchycore_alpha*fchycore/fchycurr)/
     *       (1.0_R8+fchycore_alpha)
        if(nncut.ge.2) cbsch(0,coreregn2)=
     *       cbsch(0,coreregn2)*
     *       (1.0_R8+fchycore_alpha*fchycore/fchycurr)/
     *       (1.0_R8+fchycore_alpha)
        saved_cbsch_core=cbsch(0,coreregno)
        if(diagno.ge.1)
     >   write(*,*) 'Rescaling cbsch to ',cbsch(0,coreregno)
      endif

      if(ndes.gt.0.0_R8) then
        if (ncall.eq.0) then
          if (saved_cbsna_pfr1.ne.0.0_R8)
     >      cbsna(0,isfeedback,pfrregno1)=saved_cbsna_pfr1
     .       *neutral_sources_rescale
          if (saved_cbsna_pfr2.ne.0.0_R8)
     >      cbsna(0,isfeedback,pfrregno2)=saved_cbsna_pfr2
     .       *neutral_sources_rescale
        endif
        ntot = 0.0_R8
        do iy=-1,ny
          do ix=-1,nx
            do is = 0, ns-1
              if (nint(zn(is)).eq.nint(zn(isfeedback)).and.
     &            abs(am(is)-am(isfeedback)).lt.amtol)
     >       ntot = ntot + na(ix,iy,is)*vol(ix,iy)
            enddo
          enddo
        enddo
        cbsnainc=(ndes-ntot)*cbsnafac
        if(ncall.gt.2) then
          if((ntot-den(1))*(den(1)-den(0)).lt.0.0_R8) then
            zeroslope=.true.
            rcbsna(0)=cbsna(0,isfeedback,pfrregno1)
            if (pfrregno1.ne.pfrregno2)
     >       rcbsna(1)=cbsna(0,isfeedback,pfrregno2)
            if (pfrregno1.ne.pfrregno2) then
              if (diagno.ge.1) write(*,*)
     1         'ZEROSLOPE:ncall,cbsna(0,',isfeedback,',ireg)',ncall,
     1          cbsna(0,isfeedback,pfrregno1),
     1          cbsna(0,isfeedback,pfrregno2)
            else
              if (diagno.ge.1) write(*,*)
     1         'ZEROSLOPE:ncall,cbsna(0,',isfeedback,',ireg)',ncall,
     1          cbsna(0,isfeedback,pfrregno1)
            endif
          endif
          if(zeroslope.and.(ntot-ndes)*(den(1)-ndes).lt.0) then
            cbsna(0,isfeedback,pfrregno1)=rcbsna(0)
            if (pfrregno1.ne.pfrregno2)
     >       cbsna(0,isfeedback,pfrregno2)=rcbsna(1)
            cbsnainc=0
            zeroslope=.false.
            if (pfrregno1.ne.pfrregno2) then
              if(diagno.ge.1) write(*,*)
     >         'ZEROCROSS:ncall,cbsna(0,',isfeedback,',ireg)',ncall,
     1         cbsna(0,isfeedback,pfrregno1),
     1         cbsna(0,isfeedback,pfrregno2)
            else
              if (diagno.ge.1) write(*,*)
     >         'ZEROCROSS:ncall,cbsna(0,',isfeedback,',ireg)',ncall,
     1         cbsna(0,isfeedback,pfrregno1)
            endif
          endif
        endif
        if(ncall.lt.ncallfeedback.and.diagno.ge.1.and.main_call) then
          write(*,*) ncall,': PFLUXPUFFFEEDBACK (-)',ntot,ndes,cbsnainc
        else
          if(diagno.ge.1)
     >     write(*,*)
     1     ncall,': PFLUXPUFFFEEDBACK (+)',ntot,ndes,cbsnainc
          cbsna(0,isfeedback,pfrregno1)=
     .     max(0.0_R8,cbsna(0,isfeedback,pfrregno1)+
     &     neutral_sources_rescale * cbsnainc)
          if (pfrregno1.ne.pfrregno2)
     >     cbsna(0,isfeedback,pfrregno2)=
     .     max(0.0_R8,cbsna(0,isfeedback,pfrregno2)+
     &     neutral_sources_rescale * cbsnainc)
        endif
        if (pfrregno1.ne.pfrregno2) then
          if (diagno.ge.1)
     >     write(*,*) 'cbsna(0,',isfeedback,',',pfrregno1,'),',
     .                'cbsna(0,',isfeedback,',',pfrregno2,') ',
     1     cbsna(0,isfeedback,pfrregno1),cbsna(0,isfeedback,pfrregno2)
        else
          if (diagno.ge.1)
     >     write(*,*) 'cbsna(0,',isfeedback,',',pfrregno1,'),',
     1     cbsna(0,isfeedback,pfrregno1)
        endif
        den(0)=den(1)
        den(1)=ntot
      endif

      if(nesepm.gt.0.0_R8) then
        if(ncall.eq.0) then
          if(saved_cbsna_core.ne.0.0_R8) then
            cbsna(0,isfeedback,coreregno)=saved_cbsna_core
            if(nncut.ge.2)
     >      cbsna(0,isfeedback,coreregn2)=saved_cbsna_core
          endif
        endif
        nesepm_curr = 0.5_R8 *
     *   (ne(jxa,jsep)+ne(topix(jxa,jsep),topiy(jxa,jsep)))
        if(save_core_flux_while_overshot.ne.0.0_R8) then
          cbsna(0,isfeedback,coreregno)=save_core_flux_while_overshot
          if(nncut.ge.2)
     >    cbsna(0,isfeedback,coreregn2)=save_core_flux_while_overshot
          save_core_flux_while_overshot=0.0_R8
        endif
        cbsna(0,isfeedback,coreregno)=cbsna(0,isfeedback,coreregno)*
     *   (1.0_R8+nesepm_alpha*nesepm/nesepm_curr)/
     *   (1.0_R8+nesepm_alpha)
        if(nncut.ge.2)
     *   cbsna(0,isfeedback,coreregn2)=cbsna(0,isfeedback,coreregn2)*
     *   (1.0_R8+nesepm_alpha*nesepm/nesepm_curr)/
     *   (1.0_R8+nesepm_alpha)
        saved_cbsna_core=cbsna(0,isfeedback,coreregno)
        if(nesepm_overshoot.ge.1.0_R8.and.
     1   nesepm_curr.gt.nesepm*nesepm_overshoot) then
          save_core_flux_while_overshot=cbsna(0,ismain,coreregno)
          cbsna(0,isfeedback,coreregno)=0.0_R8
          if (nncut.ge.2) cbsna(0,isfeedback,coreregn2)=0.0_R8
          if (diagno.ge.1) write(*,*) 'Overshot density! ',
     1     nesepm_curr,nesepm*nesepm_overshoot
        endif
        if(diagno.ge.1) then
          if(nncut.lt.2)
     >     write(*,'(a,1p,g14.7,a,g14.7)')
     1     'nesepm_curr= ',nesepm_curr,
     *     ' ; rescaling cbsna to ',cbsna(0,isfeedback,coreregno)
          if(nncut.ge.2)
     >     write(*,'(a,1p,g14.7,a,2g14.7)')
     1     'nesepm_curr= ',nesepm_curr,
     *     ' ; rescaling cbsna to ',cbsna(0,isfeedback,coreregno),
     *     cbsna(0,isfeedback,coreregn2)
        endif
      endif

      if (nesepm_sol.gt.0.0_R8.or.ndes_sol.gt.0.0_R8) then
c nesepm contains the desired separatrix density
        if (ncall.eq.0) then
c do this only once at the beginning
          if (saved_cbsna_sol.ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=saved_cbsna_sol
     .       *neutral_sources_rescale
c cbsna(0,isfeedback,solregno) contains the boundary condition we are
c changing under feedback control
          elseif (cbsna(0,isfeedback,solregno).ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=cbsna(0,isfeedback,solregno)
     .       *neutral_sources_rescale
          elseif (userfluxparm(nesepm_istra,1).ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=userfluxparm(nesepm_istra,1)
     .       *neutral_sources_rescale
          else
            cbsna(0,isfeedback,solregno)=
     1       1.0e18_R8*neutral_sources_rescale
            write(*,*) 'Warning!  ',
     1       'The gas puff was started with an initial rate of 1e18!'
          endif
c make sure that we save the feedback value before the code stops
        endif
        ntot = 0.0_R8
        nesepm_curr = 0.0_R8
        if (nesepm_sol.gt.0.0_R8) then
          nesepm_curr = 0.5_R8 *
     *     (ne(jxa,jsep)+ne(topix(jxa,jsep),topiy(jxa,jsep)))
        else if (ndes_sol.gt.0.0_R8) then
          do iy = -1, ny
            do ix = -1, nx
              do is = 0, ns-1
                if (nint(zn(is)).eq.nint(zn(isfeedback)).and.
     &              abs(am(is)-am(isfeedback)).lt.amtol)
     >            ntot = ntot + na(ix,iy,is)*vol(ix,iy)
              enddo
            enddo
          enddo
        endif
c the current value of the separatrix density
        if(save_gas_puff_while_overshot.ne.0.0_R8) then
          cbsna(0,isfeedback,solregno)=save_gas_puff_while_overshot
          save_gas_puff_while_overshot=0.0_R8
        endif
        if (nesepm_sol.gt.0.0_R8) then
          if(nesepm_beta.ne.0.0_R8) then
c new attempt
            nesepm_error = (nesepm_sol - nesepm_curr) / nesepm_sol
            if(nesepm_error * nesepm_error_old .le. 0.0_R8
     &          .and. ncall .gt. ncall_old+ncall_gap) then
              write(*,*) 'sign change in nesepm_error; ncall, flux = ',
     &          ncall, cbsna(0,isfeedback,solregno)
              if(n_flux_zero .eq. 1) then
                cbsna(0,isfeedback,solregno) = 0.5_R8 *
     &            (cbsna(0,isfeedback,solregno) + flux_old)
                write(*,*) 'Reset gas puff to ',
     &            cbsna(0,isfeedback,solregno)
                n_flux_zero = 0
                write(*,*) 'Approximate oscillation period is ',
     &            2*(ncall-ncall_old)
              else
                n_flux_zero = 1
                flux_old = cbsna(0,isfeedback,solregno)
                write(*,*) 'Saved old flux ',
     &            cbsna(0,isfeedback,solregno)
              endif
              ncall_old = ncall
            endif
            nesepm_error_old = nesepm_error
            nesepm_integral_error =
     1       nesepm_gamma * nesepm_integral_error +
     1       (1-nesepm_gamma)* nesepm_error
            write(*,'(a,1p,2e15.7)')
     .       'nesepm_error, nesepm_integral_error= ',
     1        nesepm_error, nesepm_integral_error
            cbsna(0,isfeedback,solregno) =
     .       cbsna(0,isfeedback,solregno)*
     *        (1.0_R8 + nesepm_alpha*nesepm_error +
     +         nesepm_beta*nesepm_integral_error)
          else
c old approach
            cbsna(0,isfeedback,solregno)=
     .       cbsna(0,isfeedback,solregno)*
     *        (1.0_R8+nesepm_alpha*nesepm_sol/nesepm_curr)/
     *        (1.0_R8+nesepm_alpha)
          endif
        else if (ndes_sol.gt.0.0_R8) then
          cbsna(0,isfeedback,solregno)=
     .     cbsna(0,isfeedback,solregno)*
     *      (1.0_R8+nesepm_alpha*ndes_sol/ntot)/
     *      (1.0_R8+nesepm_alpha)
        endif
        if(nesepm_minpuff.gt.0.0_R8) then
          cbsna(0,isfeedback,solregno)=
     1     max(cbsna(0,isfeedback,solregno),nesepm_minpuff)
        endif
        if(nesepm_maxpuff.gt.0.0_R8) then
          cbsna(0,isfeedback,solregno)=
     1     min(cbsna(0,isfeedback,solregno),nesepm_maxpuff)
        endif
c the key part: update the gas puff based on the ratio of desired
c value to the current value; nesepm_alpha is used to stabilise the
c system and defaults to 1e-3
        saved_cbsna_sol=cbsna(0,isfeedback,solregno)/
     .   neutral_sources_rescale
c we save the value
        if (nesepm_overshoot.ge.1.0_R8 .and.
     1    ((nesepm_sol.gt.0.0_R8 .and.
     2      nesepm_curr.gt.nesepm_sol*nesepm_overshoot) .or.
     3     (ndes_sol.gt.0.0_R8 .and.
     4      ntot.gt.ndes_sol*nesepm_overshoot))) then
c if we have overshot the density by more than
c nesepm_sol*nesepm_overshoot then we force the gas puff to 0
          save_gas_puff_while_overshot=cbsna(0,isfeedback,solregno)
          cbsna(0,isfeedback,solregno)=0.0_R8
          if (diagno.ge.1 .and. nesepm_sol.gt.0.0_R8)
     >     write(*,*) 'Overshot density! ',
     1     nesepm_curr,nesepm_sol*nesepm_overshoot
          if (diagno.ge.1 .and. ndes_sol.gt.0.0_R8)
     >     write(*,*) 'Overshot density! ',
     1     ntot,ndes_sol*nesepm_overshoot
        endif
        if (diagno.ge.1 .and. nesepm_sol.gt.0.0_R8)
     >   write(*,'(a,1p,g14.7,a,g14.7,a,g14.7)')
     1    ' nesepm_curr= ',nesepm_curr,
     *    ' target = ',nesepm_sol,
     *    ' ; rescaling cbsna to ',cbsna(0,isfeedback,solregno)
        if (diagno.ge.1 .and. ndes_sol.gt.0.0_R8)
     >   write(*,'(a,1p,g14.7,a,g14.7,a,g14.7)')
     1    ' ntot= ',ntot,
     *    ' target = ',ndes_sol,
     *    ' ; rescaling cbsna to ',cbsna(0,isfeedback,solregno)
      endif

      if (nesepm_pfr.gt.0.0_R8) then
        if (ncall.eq.0) then
          if (saved_cbsna_pfr1.ne.0.0_R8)
     >      cbsna(0,isfeedback,pfrregno1)=saved_cbsna_pfr1
     .       *neutral_sources_rescale
          if (saved_cbsna_pfr2.ne.0.0_R8)
     >      cbsna(0,isfeedback,pfrregno2)=saved_cbsna_pfr2
     .       *neutral_sources_rescale
        endif
      endif
      if(private_flux_puff.gt.0.0_R8.or.nesepm_pfr.gt.0.0_R8) then
        if(ncall.eq.0) then
          private_flux_area=0.0_R8
          if (saved_cbsna_pfr1.ne.0.0_R8) then
            cbsna(0,isfeedback,pfrregno1)=saved_cbsna_pfr1
     .       *neutral_sources_rescale
          elseif (cbsna(0,isfeedback,pfrregno1).ne.0.0_R8) then
            cbsna(0,isfeedback,pfrregno1)=cbsna(0,isfeedback,pfrregno1)
     .       *neutral_sources_rescale
          elseif (userfluxparm(nesepm_istra,1).ne.0.0_R8) then
            if (pfrregno1.eq.pfrregno2) then
              cbsna(0,isfeedback,pfrregno1)=userfluxparm(nesepm_istra,1)
     .         *neutral_sources_rescale
            else
              cbsna(0,isfeedback,pfrregno1)=userfluxparm(nesepm_istra,1)
     .         *neutral_sources_rescale/2.0_R8
            endif
          else
            cbsna(0,isfeedback,pfrregno1)=
     1       1.0e18_R8*neutral_sources_rescale
            write(*,*) 'Warning!  ',
     1       'The gas puff was started with an initial rate of 1e18!'
          endif
          if (pfrregno1.ne.pfrregno2) then
            if (saved_cbsna_pfr2.ne.0.0_R8) then
              cbsna(0,isfeedback,pfrregno2)=saved_cbsna_pfr2
     .         *neutral_sources_rescale
            elseif (cbsna(0,isfeedback,pfrregno2).ne.0.0_R8) then
              cbsna(0,isfeedback,pfrregno2)=
     .         cbsna(0,isfeedback,pfrregno2)*neutral_sources_rescale
            elseif (userfluxparm(nesepm_istra,1).ne.0.0_R8) then
              cbsna(0,isfeedback,pfrregno2)=userfluxparm(nesepm_istra,1)
     .         *neutral_sources_rescale/2.0_R8
            else
              cbsna(0,isfeedback,pfrregno2)=
     1         1.0e18_R8*neutral_sources_rescale
            endif
          endif
          do ix=-1,nx
            if((mod(region(ix,-1,0),4).eq.3 .or.
     &         (mod(region(ix,-1,0),4).eq.0 .and. region(ix,-1,0).ne.0))
     &       .and.((ix.ge.ixl1 .and. ix.le.ixr1) .or.
     &             (ix.ge.ixl2 .and. ix.le.ixr2)))
     >       private_flux_area=private_flux_area+
     2        gs(topix(ix,-1),topiy(ix,-1),1)
          enddo
          if(diagno.ge.1) then
            write(*,*) 'private_flux_area = ',private_flux_area
            if (pfrregno1.ne.pfrregno2) then
              write(*,'(a,1p,2g14.7)')
     .        'Original cbsna(0,isfeedback,pfrregno1&2) = ',
     1         cbsna(0,isfeedback,pfrregno1),
     2         cbsna(0,isfeedback,pfrregno2)
            else
              write(*,'(a,1p,g14.7)')
     .        'Original cbsna(0,isfeedback,pfrregno1) = ',
     1         cbsna(0,isfeedback,pfrregno1)
            endif
          endif
        endif
        if (nesepm_pfr.gt.0.0_R8) then
          nesepm_curr = 0.5_R8 *
     *     (ne(jxa,jsep)+ne(topix(jxa,jsep),topiy(jxa,jsep)))
          if(save_pfr1_puff_while_overshot.ne.0.0_R8) then
            cbsna(0,isfeedback,pfrregno1)=save_pfr1_puff_while_overshot
            save_pfr1_puff_while_overshot=0.0_R8
          endif
          if(save_pfr2_puff_while_overshot.ne.0.0_R8) then
            cbsna(0,isfeedback,pfrregno2)=save_pfr2_puff_while_overshot
            save_pfr2_puff_while_overshot=0.0_R8
          endif
          cbsna(0,isfeedback,pfrregno1)=
     .     cbsna(0,isfeedback,pfrregno1)*
     *     (1.0_R8+nesepm_alpha*nesepm_pfr/nesepm_curr)/
     *     (1.0_R8+nesepm_alpha)
          saved_cbsna_pfr1=cbsna(0,isfeedback,pfrregno1)/
     .     neutral_sources_rescale
          if (pfrregno1.ne.pfrregno2) then
            cbsna(0,isfeedback,pfrregno2)=
     .       cbsna(0,isfeedback,pfrregno2)*
     *       (1.0_R8+nesepm_alpha*nesepm_pfr/nesepm_curr)/
     *       (1.0_R8+nesepm_alpha)
            saved_cbsna_pfr2=cbsna(0,isfeedback,pfrregno2)/
     .       neutral_sources_rescale
          endif
          if(nesepm_overshoot.ge.1.0_R8.and.
     1     nesepm_curr.gt.nesepm_pfr*nesepm_overshoot) then
            save_pfr1_puff_while_overshot=cbsna(0,isfeedback,pfrregno1)
            save_pfr2_puff_while_overshot=cbsna(0,isfeedback,pfrregno2)
            cbsna(0,isfeedback,pfrregno1)=0.0_R8
            cbsna(0,isfeedback,pfrregno2)=0.0_R8
            if (diagno.ge.1) write(*,*) 'Overshot density! ',
     1       nesepm_curr,nesepm_pfr*nesepm_overshoot
          endif
          if (diagno.ge.1)
     >     write(*,'(a,1p,g14.7,a,g14.7,a,2g14.7)')
     1     'nesepm_curr= ',nesepm_curr,
     *      ' target = ',nesepm_pfr,
     *      ' ; rescaling cbsna to ',cbsna(0,isfeedback,pfrregno1),
     *      cbsna(0,isfeedback,pfrregno2)
        elseif (private_flux_puff.gt.0.0_R8) then
          target_flux=0.0_R8
          do ix=-1,nx
            if (leftix(ix,-1).eq.-2.and.
     &            ((ix.ge.ixl1 .and. ix.le.ixr1) .or.
     &             (ix.ge.ixl2 .and. ix.le.ixr2))) then
              do iy=-1,ny
                target_flux=target_flux
     1           -rza(rightix(ix,iy),rightiy(ix,iy),ismain)
     2           *fna(rightix(ix,iy),rightiy(ix,iy),0,ismain)
              enddo
            elseif (rightix(ix,-1).eq.nx+1.and.
     &            ((ix.ge.ixl1 .and. ix.le.ixr1) .or.
     &             (ix.ge.ixl2 .and. ix.le.ixr2))) then
              do iy=-1,ny
                target_flux=target_flux
     1           +rza(ix,iy,ismain)*fna(ix,iy,0,ismain)
              enddo
            endif
          enddo
          if(diagno.ge.1) write(*,*) 'target_flux = ',target_flux
          cbsna(0,isfeedback,pfrregno1)= neutral_sources_rescale*
     .     target_flux*private_flux_puff/private_flux_area
          if (pfrregno1.ne.pfrregno2)
     >     cbsna(0,isfeedback,pfrregno2)= neutral_sources_rescale*
     .     target_flux*private_flux_puff/private_flux_area
        endif
        if (pfrregno1.ne.pfrregno2) then
          if (diagno.ge.1)
     >     write(*,'(a,1p,2g14.7)')
     .     'cbsna(0,isfeedback,pfrregno1&2) = ',
     1      cbsna(0,isfeedback,pfrregno1),cbsna(0,isfeedback,pfrregno2)
        else
          if (diagno.ge.1)
     >     write(*,'(a,1p,g14.7)') 'cbsna(0,isfeedback,pfrregno1) = ',
     1     cbsna(0,isfeedback,pfrregno1)
        endif
      endif

      if (volrec_sol.gt.0.0_R8) then
        if (ncall.eq.0) then
          if (saved_cbsna_sol.ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=saved_cbsna_sol
     .       *neutral_sources_rescale
          elseif (cbsna(0,isfeedback,solregno).ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=cbsna(0,isfeedback,solregno)
     .       *neutral_sources_rescale
          elseif (userfluxparm(nesepm_istra,1).ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=userfluxparm(nesepm_istra,1)
     .       *neutral_sources_rescale
          else
            cbsna(0,isfeedback,solregno)=
     1       1.0e18_R8*neutral_sources_rescale
            write(*,*) 'Warning!  ',
     1       'The gas puff was started with an initial rate of 1e18!'
          endif
        endif
        volrec_curr = cum_volrec
        if(volrec_curr.le.0.0_R8) then
          write(*,*) 'volrec_curr <= 0 ', volrec_curr
        else
          if(save_gas_puff_while_overshot.ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=save_gas_puff_while_overshot
            save_gas_puff_while_overshot=0.0_R8
          endif
          cbsna(0,isfeedback,solregno)=
     .     cbsna(0,isfeedback,solregno)*
     *     (1.0_R8+volrec_alpha*(volrec_sol/volrec_curr)**volrec_beta)/
     *     (1.0_R8+volrec_alpha)
          saved_cbsna_sol=cbsna(0,isfeedback,solregno)/
     .     neutral_sources_rescale
          if(volrec_overshoot.ge.1.0_R8.and.
     1     volrec_curr.gt.volrec_sol*volrec_overshoot) then
            save_gas_puff_while_overshot=cbsna(0,isfeedback,solregno)
            cbsna(0,isfeedback,solregno)=0.0_R8
            if (diagno.ge.1) write(*,*) 'Overshot density! ',
     1       volrec_curr,volrec_sol*volrec_overshoot
          endif
        endif
        if (diagno.ge.1)
     >   write(*,*) 'volrec_curr= ',volrec_curr,
     *    ' target = ',volrec_sol,
     *    ' ; rescaling cbsna to ',cbsna(0,isfeedback,solregno)
      endif

      if (nepedm_sol.gt.0.0_R8) then
c nepedm contains the desired density at the top of the pedestal
        if (ncall.eq.0) then
c do this only once at the beginning
          if (saved_cbsna_sol.ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=saved_cbsna_sol
     .       *neutral_sources_rescale
c cbsna(0,isfeedback,solregno) contains the boundary condition we are
c changing under feedback control
          elseif (cbsna(0,isfeedback,solregno).ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=cbsna(0,isfeedback,solregno)
     .       *neutral_sources_rescale
          elseif (userfluxparm(nesepm_istra,1).ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=userfluxparm(nesepm_istra,1)
     .       *neutral_sources_rescale
          else
            cbsna(0,isfeedback,solregno)=
     1       1.0e18_R8*neutral_sources_rescale
            write(*,*) 'Warning!  ',
     1       'The gas puff was started with an initial rate of 1e18!'
          endif
c make sure that we save the feedback value before the code stops
        endif
        ntot = 0.0_R8
c the current value of the density at the top of the pedestal
        nesepm_curr = ne(jxa,iyped)
        if(save_gas_puff_while_overshot.ne.0.0_R8) then
          cbsna(0,isfeedback,solregno)=save_gas_puff_while_overshot
          save_gas_puff_while_overshot=0.0_R8
        endif
        if(nesepm_beta.ne.0.0_R8) then
c new attempt
          nesepm_error = (nesepm_sol - nesepm_curr) / nepedm_sol
          if(nesepm_error * nesepm_error_old .le. 0.0_R8
     &       .and. ncall .gt. ncall_old+ncall_gap) then
            write(*,*) 'sign change in nesepm_error; ncall, flux = ',
     &       ncall, cbsna(0,isfeedback,solregno)
            if(n_flux_zero .eq. 1) then
              cbsna(0,isfeedback,solregno) = 0.5_R8 *
     &         (cbsna(0,isfeedback,solregno) + flux_old)
              write(*,*) 'Reset gas puff to ',
     &         cbsna(0,isfeedback,solregno)
              n_flux_zero = 0
              write(*,*) 'Approximate oscillation period is ',
     &         2*(ncall-ncall_old)
            else
              n_flux_zero = 1
              flux_old = cbsna(0,isfeedback,solregno)
              write(*,*) 'Saved old flux ',
     &         cbsna(0,isfeedback,solregno)
            endif
            ncall_old = ncall
          endif
          nesepm_error_old = nesepm_error
          nesepm_integral_error =
     1     nesepm_gamma * nesepm_integral_error +
     1      (1-nesepm_gamma)* nesepm_error
          write(*,'(a,1p,2e14.7)')
     .     'nesepm_error, nesepm_integral_error= ',
     1      nesepm_error, nesepm_integral_error
          cbsna(0,isfeedback,solregno) =
     .     cbsna(0,isfeedback,solregno)*
     *      (1.0_R8 + nesepm_alpha*nesepm_error +
     +       nesepm_beta*nesepm_integral_error)
        else
c old approach
          cbsna(0,isfeedback,solregno)=cbsna(0,isfeedback,solregno)*
     *     (1.0_R8+nesepm_alpha*nepedm_sol/nesepm_curr)/
     *     (1.0_R8+nesepm_alpha)
        endif
        if(nesepm_minpuff.gt.0.0_R8) cbsna(0,isfeedback,solregno)=
     1   max(cbsna(0,isfeedback,solregno),nesepm_minpuff)
        if(nesepm_maxpuff.gt.0.0_R8) cbsna(0,isfeedback,solregno)=
     1   min(cbsna(0,isfeedback,solregno),nesepm_maxpuff)
c the key part: update the gas puff based on the ratio of desired
c value to the current value; nesepm_alpha is used to stabilise the
c system and defaults to 1e-3
        saved_cbsna_sol=cbsna(0,isfeedback,solregno)/
     .   neutral_sources_rescale
c we save the value
        if (nesepm_overshoot.ge.1.0_R8 .and.
     2      nesepm_curr.gt.nepedm_sol*nesepm_overshoot) then
c if we have overshot the density by more than
c nepedm_sol*nesepm_overshoot then we force the gas puff to 0
          save_gas_puff_while_overshot=cbsna(0,isfeedback,solregno)
          cbsna(0,isfeedback,solregno)=0.0_R8
          if (diagno.ge.1)
     >     write(*,*) 'Overshot density! ',
     1     nesepm_curr,nepedm_sol*nesepm_overshoot
        endif
        if (diagno.ge.1)
     >   write(*,'(a,1p,g20.10,a,g20.10,a,g20.10)')
     1    ' nesepm_curr= ',nesepm_curr,
     *    ' target = ',nepedm_sol,
     *    ' ; rescaling cbsna to ',cbsna(0,isfeedback,solregno)
      endif
      
* ..return
      if(ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      return


*-----------------------------------------------------------------------
*.end b2stbc_fb_pre

      end subroutine b2stbc_fb_pre

#ifdef B25_EIRENE
*.flux balance computation routine 
c    Relays on the EIRENE data, cannot be used in standalone B25 runs
      subroutine calc_pump_and_flux(nx, ny, ispecies,
     &     mdf_fnb, style_b2mndt, sna_pump,
     &     flxnp, flxic, flxnc, flxex, flxsp)
      use b2mod_types
      use eirmod_wneutrals
      use b2mod_eirene
     ,   , only : eirpump, eirspta, eirsptm
      use eirmod_braeir
      use eirmod_comusr
     ,   , only : lkindp
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_diag
      use b2mod_b2cmpa
     ,   , only : is_neutral
      use b2mod_plasma
     ,   , only : fna_mdf, fna
      use b2mod_geo
     ,   , only : gs
      use b2mod_sources
     ,   , only : ext_sna
       
      implicit none
       
      integer ::  j, k, is, ispecies, ic, idir, ia, im, ip,
     &    ix, ixf, iy, iyf, mdf_fnb, style_b2mndt,
     &    js, nx, ny
      real (kind=R8) :: flxnp, flxic, flxnc, flxex, flxsp,
     &    sna_pump, bsign
        
      sna_pump = 0.0_R8  ! pumping
      flxic = 0.0_R8     ! ion through the core
      flxnc = 0.0_R8     ! neutral through the core
      flxnp = 0.0_R8     ! puffing
      flxex = 0.0_R8     ! external source
      flxsp = 0.0_R8     ! sputtering flux
!     flx.. - positive if source particles into the system
!     sna_pump - positive if pump particles out of the system

c Gather eirene pumping data (at this point includes neutral to the core)

      do ia = 1, nnatmi
        if (latmscl(ia).ne.ispecies+1) cycle
        sna_pump = sna_pump + eirpump(ia)
        flxsp = flxsp + eirspta(ia)
        do im = 1, nnmoli
          if (mlcmp(ia,im).gt.0) then
            sna_pump = sna_pump + mlcmp(ia,im) *eirpump(nnatmi+im)
            flxsp = flxsp + mlcmp(ia,im) *eirsptm(im)
          endif
        enddo
      end do
     
      do ip = 1, nnplsi
        if (latmscl(lkindp(ip)).ne.ispecies+1) cycle
        sna_pump = sna_pump +
     &            eirpump(nnatmi+nnmoli+nnioni+ip)
      enddo
     
c  Calculate ion flux through the core
      do is = EB2SPCR(ispecies+1), EB2SPCR(ispecies+1) + ! Loop over the isonuclear sequence 
     &       NFLUIDS(ispecies+1)

        do ic = 1,nbc ! Loop over B2 boundaries

          if (BCCHAR(IC).eq.'S') then
            idir  = 1
            bsign = 1.0_R8
          else if (BCCHAR(IC).eq.'N') then
            idir  = 1
            bsign =-1.0_R8
          else if (BCCHAR(IC).eq.'W') then
            idir  = 0
            bsign = 1.0_R8
          else if (BCCHAR(IC).eq.'E') then
            idir  = 0
            bsign =-1.0_R8
          else
            cycle
          end if

          do ix = 1, BC_LIST_SIZE(ic) ! Loop over poloidal coordinate along this boundary

            if (BCCHAR(ic).EQ.'S') then
              ixf = topix(BC_LIST_X(ix,ic),BC_LIST_Y(ix,ic))
              iyf = topiy(BC_LIST_X(ix,ic),BC_LIST_Y(ix,ic))
            else if (BCCHAR(ic).EQ.'N') then
              ixf = BC_LIST_X(ix,ic)
              iyf = BC_LIST_Y(ix,ic)
            else if (BCCHAR(ic).EQ.'W') then
              ixf = rightix(BC_LIST_X(ix,ic),BC_LIST_Y(ix,ic))
              iyf = rightiy(BC_LIST_X(ix,ic),BC_LIST_Y(ix,ic))
            else if (BCCHAR(ic).eq.'E') then
              ixf = BC_LIST_X(ix,ic)
              iyf = BC_LIST_Y(ix,ic)
            end if

            if (core_boundary(ic)) then ! Check whether this is part of the core boundary
              if (.not.is_neutral(is)) then
                if (mdf_fnb .eq. 0 .and. style_b2mndt .le. 0) then
                  flxic = flxic + bsign*fna(ixf,iyf,idir,is)
                else
                  flxic = flxic + bsign*fna_mdf(ixf,iyf,idir,is)
                end if
              end if
            end if  

          end do ! End Loop over poloidal coordinate along the boundary
        end do ! End loop over B2 boundaries 

c  Calculate the external source
        do iy=0,ny-1
          do ix=0,nx-1
            j=region(ix,iy,0)
            if ((j.ne.0) .and. (leftix(ix,iy).ne.-2) .and.
     &           (rightix(ix,iy).ne.nx+1) .and. 
     &            (is.eq.EB2SPCR(ispecies+1)))
     .             flxex = flxex + ext_sna(ix,iy,is)
          end do
        end do
      end do ! End of the loop over species within given isonuclear sequence (is)

c  Calculate neutral flux to the core   
      do is = 1,ncns  ! Loop over segments of Eirene core boundary
        js = lcns(is) + nnlimi
        do j = 1,nnatmi
          if (latmscl(j) .eq. ispecies+1) then
            flxnc = flxnc - wldna(js,j,0) + wldra(js,j,0) 
     &             + wldpa(js,j,0)
            do ic=1,nnmoli
              if (mlcmp(j,ic) .gt. 0) flxnc = 
     &                flxnc - mlcmp(j,ic) * (wldnm(js,ic,0)
     &               - wldrm(js,ic,0) - wldpm(js,ic,0))
            end do
          end if
        end do
      end do ! End loop over segments of Eirene core boundary

c  Calculate seeding
      do ic = 1,nstrat ! Loop over strata to gather Eirene seeding
        if (crcstra(ic) .eq. 'C') then  ! Puffing or seeding
          do ia = 1,nnatmi
            if (latmscl(ia).ne.ispecies+1) cycle
            flxnp = flxnp +
     &                 userfluxparm(ic,1) * gpfc(ia,ic)
          end do
        end if
      end do ! End loop over strata

c  substract neutral flux to the core from pumping 
      sna_pump = sna_pump + flxnc
*-----------------------------------------------------------------------
*.end of flux balance calculation routine
      
      end subroutine calc_pump_and_flux
      
*.eirene wall recycling update routine
      subroutine eirene_recyct_update(ispecies, absr, spump_list)
        use eirmod_clgin
        use b2mod_neutrals_namelist
        use b2mod_dimensions
        implicit none
        integer, intent(in) :: ispecies,   ! Index of B2.5 species (or isonuclear sequence)
     .                         spump_list(DEF_NLIM+DEF_NSTS) ! list of surfaces to update
        real (kind=R8), intent(in) :: absr ! new value of the absorbtion coefficient
        integer ja, l, jm, k
        
      do l = 1, DEF_NLIM+DEF_NSTS
        k = spump_list(l)
        if (k.eq.0) exit
        if (k.lt.0) k = nlim - k
        do ja = 1,nnatmi
          if (latmscl(ja) .eq. ispecies) then
            recyct(ja,k) = 1.0_R8 - absr
            do jm = 1,nnmoli                  
              if (mlcmp(ja,jm) .gt. 0) then
                recyct(jm+nnatmi,k) = 1.0_R8 - absr
              endif
            enddo
          endif
        enddo
      enddo

*-----------------------------------------------------------------------
*.end of eirene wall recycling update routine

      end subroutine eirene_recyct_update

*.eirene pump element list construction routine
      subroutine get_pump_list(ispecies, spcsrf_id, spump_list)
        use eirmod_clgin
        use b2mod_neutrals_namelist
        use b2mod_dimensions
        implicit none
        integer, intent(in) :: ispecies, spcsrf_id
        integer, intent(inout) :: spump_list(DEF_NLIM+DEF_NSTS)
        integer ja, l, jm, k, i
        
      i = 1
      if (spcsrf_id.gt.0) then !take surfaces list form spcsrf group directly
        do l=i_spcsrf(spcsrf_id), j_spcsrf(spcsrf_id)
          spump_list(i) = l_spcsrf(l)
          i = i + 1
        enddo
      else  !take all the pumping surfaces for given ispecies
        do l = 1, nnlimi + ncns - 1
          do ja = 1,nnatmi
            if (latmscl(ja).eq.ispecies) then
              if (recyct(ja,l).lt.1.0_R8) then
                spump_list(i) = l
                i = i + 1
              endif
            endif
          enddo
        enddo
      endif

      if (spump_list(1).eq.0) 
     &  write(*,'(2a,i3)') "get_pump_list: Warning! ",
     &      "No pumping surfaces found for sequence: ",ispecies
          
*-----------------------------------------------------------------------
*.end of eirene pump element list construction routine

      end subroutine get_pump_list
#endif
