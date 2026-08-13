*-----------------------------------------------------------------------
*.specification

      subroutine b2stfm (nCv, nFc, ns,
     &  switch, geo, mpg, pl, dv, rt, rtw, sr_mol, sr)
      use b2mod_types
      use b2mod_diag
      use b2mod_tallies
      use b2mod_constants
      use b2mod_eirene_globals
      use b2mod_b2cmpa
      use b2mod_b2cmrc
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_subsys
      use b2mod_dimensions
      use b2mod_fluid_molecules
      use b2mod_neutrals_namelist
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, ns, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2RatesWork), intent (in) :: rtw
      type (B2SourceMol), intent (inout) :: sr_mol
      type (B2Source), intent (inout) :: sr

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STFM provides linearised expressions for the sources of
*     density, parallel momentum, electron heat and atom heat, due to
*     interactions with fluid molecules, as well as the sources for those
*     molecules themselves.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, is, ism, is_a, is_i
      character chns*3!,
      real (kind=R8) :: pit(nCv), pit2(nCv)
      real (kind=R8) :: !! rate, number of reaction events per cell per second
     *  rate_d(nCv),  !! dissociation
     *  rate_nd(nCv), !! non-dissociative ionization
     *  rate_di1(nCv),!! dissociative ionization (from neutral molecule)
     *  rate_ic(nCv), !! ion conversion
     *  rate_ec_m(nCv), !! elastic collisions (mom)
     *  rate_ec_E(nCv), !! elastic collisions (energy)
     *  rate_dr(nCv), !! dissociative recombination
     *  rate_de(nCv), !! dissociative excitation
     *  rate_di2(nCv) !! dissociative ionization (from molecular ion)
cc sources for plasma and fluid atoms
      real (kind=R8) ::
     *  sna0_d(nCv,0:1,0:ns-1),  !! dissociation
     *  sna0_nd(nCv,0:1,0:ns-1), !! non-dissociative ionization
     *  sna0_di1(nCv,0:1,0:ns-1),!! dissociative ionization (from neutral molecule)
     *  sna0_ic(nCv,0:1,0:ns-1), !! ion conversion
     *  sna0_dr(nCv,0:1,0:ns-1), !! dissociative recombination
     *  sna0_de(nCv,0:1,0:ns-1), !! dissociative excitation
     *  sna0_di2(nCv,0:1,0:ns-1) !! dissociative ionization (from molecular ion)
      real (kind=R8) ::
     *  smo0_d(nCv,0:1,0:ns-1),  smo0_nd(nCv,0:1,0:ns-1), 
     *  smo0_di1(nCv,0:1,0:ns-1),smo0_ic(nCv,0:1,0:ns-1), 
     *  smo0_dr(nCv,0:1,0:ns-1), smo0_de(nCv,0:1,0:ns-1), 
     *  smo0_di2(nCv,0:1,0:ns-1), smo0_ec(nCv,0:1,0:ns-1) 
      real (kind=R8) :: !! here still per species, summation for B2.5 energy equations later on
     *  shi0_d(nCv,0:1,0:ns-1),  shi0_nd(nCv,0:1,0:ns-1), 
     *  shi0_di1(nCv,0:1,0:ns-1),shi0_ic(nCv,0:1,0:ns-1), 
     *  shi0_dr(nCv,0:1,0:ns-1), shi0_de(nCv,0:1,0:ns-1), 
     *  shi0_di2(nCv,0:1,0:ns-1), shi0_ec(nCv,0:1,0:ns-1) 
      real (kind=R8) :: 
     *  she0_d(nCv,0:1),  she0_nd(nCv,0:1), she0_di1(nCv,0:1), 
     *  she0_ic(nCv,0:1), she0_dr(nCv,0:1), she0_de(nCv,0:1), 
     *  she0_di2(nCv,0:1)
      real (kind=R8) ::  
     *  sna0_tot(nCv,0:1,0:ns-1), smo0_tot(nCv,0:1,0:ns-1),
     *  shi0_tot(nCv,0:1,0:ns-1), she0_tot(nCv,0:1)
cc molecule sources
      real (kind=R8) :: 
     *  sn_mol_d(nCv), sn_mol_nd(nCv), sn_mol_di1(nCv),
     *  sn_mol_ic(nCv), sn_mol_ec(nCv), sn_mol_tot(nCv)
      real (kind=R8) :: 
     *  sm_mol_d(nCv), sm_mol_nd(nCv), sm_mol_di1(nCv),
     *  sm_mol_ic(nCv), sm_mol_ec(nCv), sm_mol_tot(nCv),
     *  sm_mol_rad_d(nCv,0:1), sm_mol_rad_nd(nCv,0:1), 
     *  sm_mol_rad_di1(nCv,0:1), sm_mol_rad_ic(nCv,0:1),
     *  sm_mol_rad_ec(nCv,0:1), sm_mol_rad_tot(nCv,0:1),
     *  sm_mol_dia_d(nCv,0:1), sm_mol_dia_nd(nCv,0:1), 
     *  sm_mol_dia_di1(nCv,0:1), sm_mol_dia_ic(nCv,0:1),
     *  sm_mol_dia_ec(nCv,0:1), sm_mol_dia_tot(nCv,0:1)
      real (kind=R8) :: 
     *  sh_mol_d(nCv), sh_mol_nd(nCv), sh_mol_di1(nCv),
     *  sh_mol_ic(nCv), sh_mol_ec(nCv), sh_mol_tot(nCv)
      real (kind=R8) :: 
     *  se_mol_d(nCv), se_mol_nd(nCv), se_mol_di1(nCv),
     *  se_mol_ic(nCv), se_mol_ec(nCv), se_mol_tot(nCv),
     *  se_ec_fr(nCv), se_ec_kin(nCv), tkin, tfr
cc molecule sources
      real (kind=R8) :: F_p, F_d, F_r, Vcom_p, Vcom_d, Vcom_r, 
     *                  wrk(nCv), vsq(nCv), vsqi(nCv)
      real (kind=R8) :: t0, freac, num_a, den_a, a0, mr, s0, t0_i, t0_a,
     *                  t0_mol, tau_inv, beta_num, beta_denom, m_d2,
     *                  total_mol_sink, mfp(nCv), vth, t0_rad, t0_dia
      real (kind=R8), allocatable :: alpha(:,:), nmi_new(:,:)
      real (kind=R8), parameter :: eps = 1e-10_R8, nm_fac = 1000.0_R8 
      real (kind=R8) :: fac_nlim(nCv,0:ns_mol-1), nm
*   ..procedures
      external xertst, sfill
      external b2xvsg
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stfm')

*   ..initialise rates to 0. TODO add nms dimension
      call sfill (nCv, 0.0_R8, rate_d, 1)
      call sfill (nCv, 0.0_R8, rate_nd, 1)
      call sfill (nCv, 0.0_R8, rate_di1, 1)
      call sfill (nCv, 0.0_R8, rate_ic, 1)
      call sfill (nCv, 0.0_R8, rate_ec_m, 1)
      call sfill (nCv, 0.0_R8, rate_ec_E, 1)
      call sfill (nCv, 0.0_R8, rate_dr, 1)
      call sfill (nCv, 0.0_R8, rate_de, 1)
      call sfill (nCv, 0.0_R8, rate_di2, 1)

*   ..initialise sources to 0
      call sfill (nCv*2*ns, 0.0_R8, sna0_d, 1)
      call sfill (nCv*2*ns, 0.0_R8, sna0_nd, 1)
      call sfill (nCv*2*ns, 0.0_R8, sna0_di1, 1)
      call sfill (nCv*2*ns, 0.0_R8, sna0_ic, 1)
      call sfill (nCv*2*ns, 0.0_R8, sna0_dr, 1)
      call sfill (nCv*2*ns, 0.0_R8, sna0_de, 1)
      call sfill (nCv*2*ns, 0.0_R8, sna0_di2, 1)
      call sfill (nCv*2*ns, 0.0_R8, sna0_tot, 1)

      call sfill (nCv*2*ns, 0.0_R8, smo0_d, 1)
      call sfill (nCv*2*ns, 0.0_R8, smo0_nd, 1)
      call sfill (nCv*2*ns, 0.0_R8, smo0_di1, 1)
      call sfill (nCv*2*ns, 0.0_R8, smo0_ic, 1)
      call sfill (nCv*2*ns, 0.0_R8, smo0_dr, 1)
      call sfill (nCv*2*ns, 0.0_R8, smo0_de, 1)
      call sfill (nCv*2*ns, 0.0_R8, smo0_di2, 1)
      call sfill (nCv*2*ns, 0.0_R8, smo0_ec, 1)
      call sfill (nCv*2*ns, 0.0_R8, smo0_tot, 1)

      call sfill (nCv*2*ns, 0.0_R8, shi0_d, 1)
      call sfill (nCv*2*ns, 0.0_R8, shi0_nd, 1)
      call sfill (nCv*2*ns, 0.0_R8, shi0_di1, 1)
      call sfill (nCv*2*ns, 0.0_R8, shi0_ic, 1)
      call sfill (nCv*2*ns, 0.0_R8, shi0_dr, 1)
      call sfill (nCv*2*ns, 0.0_R8, shi0_de, 1)
      call sfill (nCv*2*ns, 0.0_R8, shi0_di2, 1)
      call sfill (nCv*2*ns, 0.0_R8, shi0_ec, 1)
      call sfill (nCv*2*ns, 0.0_R8, shi0_tot, 1)

      call sfill (nCv*2, 0.0_R8, she0_d, 1)
      call sfill (nCv*2, 0.0_R8, she0_nd, 1)
      call sfill (nCv*2, 0.0_R8, she0_di1, 1)
      call sfill (nCv*2, 0.0_R8, she0_ic, 1)
      call sfill (nCv*2, 0.0_R8, she0_dr, 1)
      call sfill (nCv*2, 0.0_R8, she0_de, 1)
      call sfill (nCv*2, 0.0_R8, she0_di2, 1)
      call sfill (nCv*2, 0.0_R8, she0_tot, 1)

      call sfill (nCv, 0.0_R8, sn_mol_d, 1)
      call sfill (nCv, 0.0_R8, sn_mol_nd, 1)
      call sfill (nCv, 0.0_R8, sn_mol_di1, 1)
      call sfill (nCv, 0.0_R8, sn_mol_ic, 1)
      call sfill (nCv, 0.0_R8, sn_mol_ec, 1)

      call sfill (nCv, 0.0_R8, sm_mol_d, 1)
      call sfill (nCv, 0.0_R8, sm_mol_nd, 1)
      call sfill (nCv, 0.0_R8, sm_mol_di1, 1)
      call sfill (nCv, 0.0_R8, sm_mol_ic, 1)
      call sfill (nCv, 0.0_R8, sm_mol_ec, 1)

      call sfill (nCv*2, 0.0_R8, sm_mol_rad_d, 1)
      call sfill (nCv*2, 0.0_R8, sm_mol_rad_nd, 1)
      call sfill (nCv*2, 0.0_R8, sm_mol_rad_di1, 1)
      call sfill (nCv*2, 0.0_R8, sm_mol_rad_ic, 1)
      call sfill (nCv*2, 0.0_R8, sm_mol_rad_ec, 1)

      call sfill (nCv*2, 0.0_R8, sm_mol_dia_d, 1)
      call sfill (nCv*2, 0.0_R8, sm_mol_dia_nd, 1)
      call sfill (nCv*2, 0.0_R8, sm_mol_dia_di1, 1)
      call sfill (nCv*2, 0.0_R8, sm_mol_dia_ic, 1)
      call sfill (nCv*2, 0.0_R8, sm_mol_dia_ec, 1)

      call sfill (nCv, 0.0_R8, sh_mol_d, 1)
      call sfill (nCv, 0.0_R8, sh_mol_nd, 1)
      call sfill (nCv, 0.0_R8, sh_mol_di1, 1)
      call sfill (nCv, 0.0_R8, sh_mol_ic, 1)
      call sfill (nCv, 0.0_R8, sh_mol_ec, 1)

      call sfill (nCv, 0.0_R8, sn_mol_tot, 1)
      call sfill (nCv, 0.0_R8, sm_mol_tot, 1)
      call sfill (nCv, 0.0_R8, sh_mol_tot, 1)
      call sfill (nCv*2, 0.0_R8, sm_mol_rad_tot, 1)
      call sfill (nCv*2, 0.0_R8, sm_mol_dia_tot, 1)


cc temporary, total energy source, for benchmarking with matlab
      call sfill (nCv, 0.0_R8, se_mol_d, 1)
      call sfill (nCv, 0.0_R8, se_mol_nd, 1)
      call sfill (nCv, 0.0_R8, se_mol_di1, 1)
      call sfill (nCv, 0.0_R8, se_mol_ic, 1)
      call sfill (nCv, 0.0_R8, se_mol_ec, 1)
      call sfill (nCv, 0.0_R8, se_mol_tot, 1)

cc set diamagnetic velocity
      pit = dabs(geo%cvBb(:,0))/geo%cvBb(:,3)  !! pitch, b_theta / B
      pit2 = dabs(geo%cvBb(:,2))/geo%cvBb(:,3) !! b_phi / B
      if (switch%b2tfmol_3D_mom.eq.0) then
        u_dia_mol(:,0) = (u_pol_mol(:,0)-u_mol(:,0)*pit)/pit2
      endif

      m_D2 = am_d2*mp
      mr = (am(1)*am_d2)/(am(1)+am_d2)*mp
      if (switch%b2stfm_ekin_perp.eq.0) then
        vsq = u_mol(:,0)**2
        vsqi = u_moli(:,0)**2
      else
        vsq = u_mol(:,0)**2+u_dia_mol(:,0)**2+u_rad_mol(:,0)**2
        vsqi = u_moli(:,0)**2+u_dia_mol(:,0)**2+u_rad_mol(:,0)**2
      endif


* ..compute source terms
cc This part here is kept as an example of how to parallelize the sources. But current implementation is not yet parallel.
#if ( !defined(NO_OPENMP_B2STFM) && defined(_OPENMP) )
!svdk Note that we parallellized over iCv here instead of is, since
!svdk both is and is+1 indices are used and hence race conditions may
!svdk exist when parallellizing over is, which are non-trivial to 
!svdk eliminate
!$omp parallel default(none)
!$omp& private(iCv, is)
!$omp& shared(ns, mpg, switch, pl, geo)
!$omp& shared(zn, am)
!$omp& shared(sna0_d, sna0_nd, sna0_di1, sna0_ic)
!$omp& shared(sna0_dr, sna0_de, sna0_di2)
!$omp& shared(is_neutral)
#endif
      do is = 0, ns-2
*     ..particle source and heat source
#if ( !defined(NO_OPENMP_B2STEL) && defined(_OPENMP) )
!$omp do 
#endif
        do iCv = 1, mpg%nCi ! internal cells only
        enddo
#if ( !defined(NO_OPENMP_B2STEL) && defined(_OPENMP) )
!$omp end do 
!$omp do
#endif
*     ..parallel momentum source
        do iCv = 1, mpg%nCi ! internal cells only
        enddo
#if ( !defined(NO_OPENMP_B2STEL) && defined(_OPENMP) )
!$omp end do 
#endif
      enddo
#if ( !defined(NO_OPENMP_B2STEL) && defined(_OPENMP) )
!$omp end parallel
#endif

      do is = 1, ns-1
*     ..particle source and heat source
        do iCv = 1, mpg%nCi ! internal cells only
        enddo
*     ..parallel momentum source
        do iCv = 1, mpg%nCi ! internal cells only
        enddo
      enddo

cc TODO !! species are now hardcoded for testing !! create general scheme to see which particles go where

cc compute the equilibrium factor for the molecular ion density. nmi = alpha*nm
cc TODO compare this nmi with the one directly from EIRENE. Should match if all rates consistent. 
      if (.not.allocated(alpha)) then
        allocate(alpha(nCv,0:ns_mol-1))
        allocate(nmi_new(nCv,0:ns_mol-1))
      endif
      do ism = 0, ns_mol-1
        do iCv = 1, mpg%nCi
          num_a = (dv%ne(iCv)*rtw%r_nd(iCv,ism)*switch%b2stfm_freac_nd + 
     .          pl%na(iCv,1)*rtw%r_ic(iCv,ism)*switch%b2stfm_freac_ic)
          den_a = (dv%ne(iCv)*rtw%r_dr(iCv,ism)*switch%b2stfm_freac_dr +
     .             dv%ne(iCv)*rtw%r_de(iCv,ism)*switch%b2stfm_freac_de +
     .          dv%ne(iCv)*rtw%r_di2(iCv,ism)*switch%b2stfm_freac_di2)
          if (num_a.gt.0.0_R8.and.den_a.eq.0.0_R8) then
            call xerrab('molecular ions produced but not removed!')
          endif
          a0 = num_a/den_a
          alpha(iCv,ism) =  a0
          nmi_new(iCv,ism) = n_mol(iCv,ism)*a0
        enddo
      enddo


      if (switch%b2stfm_recompute_molion.ne.0) then
        n_moli(:,0) = nmi_new(:,0)
        u_moli(:,0) = u_mol(:,0)
        t_moli(:,0) = t_mol(:,0)
      endif

cc compute a factor to reduce specific source terms when the density is close to the minimum numerical denisity
      do ism = 0, ns_mol-1
        do iCv = 1, nCv
          nm = n_mol(iCv,ism)
          if (nm.gt.nm_min_ext.and.nm.lt.nm_min_ext*nm_fac) then 
            fac_nlim(iCv,ism) = (log(nm)-log(nm_min))/
     .                          (log(nm_fac*nm_min_ext)-log(nm_min_ext))
          else
            fac_nlim(iCv,ism) = 1.0_R8
          endif
        enddo
      enddo

* ..Particle sources

cc non-dissociative ionization
      freac = switch%b2stfm_freac_nd
      if (freac.ge.0.0_R8) then
        if (freac.ne.1.0_R8) then
          write(*,*) 'freac_nd /= 1 only for testing !!'
        endif
        do ism = 0, ns_mol-1
          do iCv = 1, mpg%nCi
            t0 =  n_mol(iCv,0)*dv%ne(iCv)*rtw%r_nd(iCv,ism)*
     .            freac*geo%cvVol(iCv)
            rate_nd(iCv) = t0
            !! creates a molecular ion -> no direct contribution to atoms or ions
            sn_mol_nd(iCv) = -t0
          enddo
        enddo
      endif

cc dissociation, particle source
      freac = switch%b2stfm_freac_d
      if (freac.ge.0.0_R8) then
        if (freac.ne.1.0_R8) then
          write(*,*) 'freac_d /= 1 only for testing !!'
        endif
        do ism = 0, ns_mol-1
          is_a = mol_is_rp_atom(ism) !! B2.5 species index of atom to which reaction product belongs
          is_i = mol_is_rp_ion(ism)  !! B2.5 species index of atom to which reaction product belongs
          do iCv = 1, mpg%nCi
            t0 =  n_mol(iCv,0)*dv%ne(iCv)*rtw%r_d(iCv,ism)*
     .            freac*geo%cvVol(iCv)
            rate_d(iCv) = t0
            sna0_d(iCv,0,is_a) = 2.0_R8*t0
            sn_mol_d(iCv) = -t0
          enddo
        enddo
      endif

cc dissociative ionization, particle source
      freac = switch%b2stfm_freac_di1
      if (freac.ge.0.0_R8) then
        if (freac.ne.1.0_R8) then
          write(*,*) 'freac_di1 /= 1 only for testing !!'
        endif
        do ism = 0, ns_mol-1
          is_a = mol_is_rp_atom(ism) 
          is_i = mol_is_rp_ion(ism)  
          do iCv = 1, mpg%nCi 
            t0 = n_mol(iCv,0)*dv%ne(iCv)*rtw%r_di1(iCv,ism)*freac*
     .           geo%cvVol(iCv)
            rate_di1(iCv) = t0
            sna0_di1(iCv,0,is_a) = t0
            sna0_di1(iCv,0,is_i) = t0
            sn_mol_di1(iCv) = -t0
          enddo
        enddo
      endif

cc ion conversion, particle source
      freac = switch%b2stfm_freac_ic
      if (freac.ge.0.0_R8) then
        if (freac.ne.1.0_R8) then
          write(*,*) 'freac_ic /= 1 only for testing !!'
        endif
        do ism = 0, ns_mol-1
          is_a = mol_is_rp_atom(ism) 
          is_i = mol_is_rp_ion(ism)  
          do iCv = 1, mpg%nCi 
            t0 = n_mol(iCv,0)*pl%na(iCv,1)*rtw%r_ic(iCv,ism)*freac*
     .           geo%cvVol(iCv)
            rate_ic(iCv) = t0
            sna0_ic(iCv,0,is_a) = t0
            sna0_ic(iCv,0,is_i) =-t0
            sn_mol_ic(iCv) = -t0
          enddo
        enddo
      endif

cc elastic collisions, no particle source, but set rate
      freac = switch%b2stfm_freac_ec
      if (freac.ge.0.0_R8) then
        if (freac.ne.1.0_R8) then
          write(*,*) 'freac_ec /= 1 only for testing !!'
        endif
        do ism = 0, ns_mol-1
          do iCv = 1, mpg%nCi 
            t0 = n_mol(iCv,0)*pl%na(iCv,1)*rtw%r_ec_m(iCv,ism)*freac*
     .           geo%cvVol(iCv)
            rate_ec_m(iCv) = t0
            t0 = n_mol(iCv,0)*pl%na(iCv,1)*rtw%r_ec_E(iCv,ism)*freac*
     .           geo%cvVol(iCv)
            rate_ec_E(iCv) = t0
          enddo
        enddo
      endif


cc dissociative recombination, particle source
      freac = switch%b2stfm_freac_dr
      if (freac.ge.0.0_R8) then
        if (freac.ne.1.0_R8) then
          write(*,*) 'freac_dr /= 1 only for testing !!'
        endif
        do ism = 0, ns_mol-1
          is_a = mol_is_rp_atom(ism) 
          do iCv = 1, mpg%nCi 
            t0 = n_moli(iCv,0)*dv%ne(iCv)*rtw%r_dr(iCv,ism)*freac*
     .           geo%cvVol(iCv)
            rate_dr(iCv) = t0
            sna0_dr(iCv,0,is_a) = 2.0_R8*t0
          enddo
        enddo
      endif

cc dissociative excitation, particle source
      freac = switch%b2stfm_freac_de
      if (freac.ge.0.0_R8) then
        if (freac.ne.1.0_R8) then
          write(*,*) 'freac_de /= 1 only for testing !!'
        endif
        do ism = 0, ns_mol-1
          is_a = mol_is_rp_atom(ism) 
          is_i = mol_is_rp_ion(ism)  
          do iCv = 1, mpg%nCi 
            t0 = n_moli(iCv,0)*dv%ne(iCv)*rtw%r_de(iCv,ism)*freac*
     .           geo%cvVol(iCv)
            rate_de(iCv) = t0
            sna0_de(iCv,0,is_a) = t0
            sna0_de(iCv,0,is_i) = t0
          enddo
        enddo
      endif

cc dissociative ionization (2), particle source
      freac = switch%b2stfm_freac_di2
      if (freac.ge.0.0_R8) then
        if (freac.ne.1.0_R8) then
          write(*,*) 'freac_di2 /= 1 only for testing !!'
        endif
        do ism = 0, ns_mol-1
          is_i = mol_is_rp_ion(ism)  
          do iCv = 1, mpg%nCi 
            t0 = n_moli(iCv,0)*dv%ne(iCv)*rtw%r_di2(iCv,ism)*freac*
     .           geo%cvVol(iCv)
            rate_di2(iCv) = t0
            sna0_di2(iCv,0,is_i) = 2.0_R8*t0
          enddo
        enddo
      endif



* ..Parallel momentum sources

cc non-dissociative ionization
      !! creates a molecular ion, no direct momentum transfer to ions or atoms
      do ism = 0, ns_mol-1
        do iCv = 1, mpg%nCi
          t0 =  am_d2*mp*rate_nd(iCv)*u_mol(iCv,0)*geo%cvHz(iCv)
          sm_mol_nd(iCv) = -t0
          !! perpendicular directions
          t0_rad =  am_d2*mp*rate_nd(iCv)
          t0_dia =  am_d2*mp*rate_nd(iCv)
          sm_mol_rad_nd(iCv,1) = -t0_rad
          sm_mol_dia_nd(iCv,1) = -t0_dia
        enddo
      enddo

cc dissociation, parallel momentum source
      do ism = 0, ns_mol-1
        is_a = mol_is_rp_atom(ism) 
        do iCv = 1, mpg%nCi
          t0 =  am_d2*mp*rate_d(iCv)*u_mol(iCv,0)*geo%cvHz(iCv)
          smo0_d(iCv,0,is_a) = t0
          sm_mol_d(iCv) = -t0
          !! perpendicular directions
          t0_rad =  am_d2*mp*rate_d(iCv)
          t0_dia =  am_d2*mp*rate_d(iCv)
          sm_mol_rad_d(iCv,1) = -t0_rad
          sm_mol_dia_d(iCv,1) = -t0_dia
        enddo
      enddo

cc dissociative ionization (1), parallel momentum source
      do ism = 0, ns_mol-1
        is_a = mol_is_rp_atom(ism) 
        is_i = mol_is_rp_ion(ism)  
        do iCv = 1, mpg%nCi
          t0 =  am_d2*mp*rate_di1(iCv)*u_mol(iCv,0)*geo%cvHz(iCv)
          smo0_di1(iCv,0,is_a) = 0.5_R8*t0 !! momentum equally split over ion and atom
          smo0_di1(iCv,0,is_i) = 0.5_R8*t0
          sm_mol_di1(iCv) = -t0
          !! perpendicular direction
          t0_rad =  am_d2*mp*rate_di1(iCv)
          t0_dia =  am_d2*mp*rate_di1(iCv)
          sm_mol_rad_di1(iCv,1) = -t0_rad
          sm_mol_dia_di1(iCv,1) = -t0_dia
        enddo
      enddo

cc ion conversion, parallel momentum source
      !! the ion becomes atom
      do ism = 0, ns_mol-1
        is_a = mol_is_rp_atom(ism) 
        is_i = mol_is_rp_ion(ism)  
        do iCv = 1, mpg%nCi
          t0 =   am(1)*mp*rate_ic(iCv)*pl%ua(iCv,1)*geo%cvHz(iCv) !! TODO generalize indices
          smo0_ic(iCv,0,is_a) =  t0
          smo0_ic(iCv,0,is_i) = -t0
          t0 =  am_d2*mp*rate_ic(iCv)*u_mol(iCv,0)*geo%cvHz(iCv)
          sm_mol_ic(iCv) = -t0
          !! perpendicular
          t0_rad =  am_d2*mp*rate_ic(iCv)
          t0_dia =  am_d2*mp*rate_ic(iCv)
          sm_mol_rad_ic(iCv,1) = -t0_rad
          sm_mol_dia_ic(iCv,1) = -t0_dia
        enddo
      enddo

cc elastic collisions, parallel momentum source
      do ism = 0, ns_mol-1
        is_i = mol_is_rp_ion(ism)  
        do iCv = 1, mpg%nCi
          t0 =  mr*rate_ec_m(iCv)*(u_mol(iCv,0)-pl%ua(iCv,1))*
     .          geo%cvHz(iCv) 
          smo0_ec(iCv,0,is_i) = t0
          sm_mol_ec(iCv) = -t0
          !! perpendicular
          t0_rad =  mr*rate_ec_m(iCv)
          t0_dia =  mr*rate_ec_m(iCv)
          sm_mol_rad_ec(iCv,1) = -t0_rad
          sm_mol_dia_ec(iCv,1) = -t0_dia
        enddo
      enddo

cc dissociative recombination, parallel momentum source
      do ism = 0, ns_mol-1
        is_a = mol_is_rp_atom(ism) 
        do iCv = 1, mpg%nCi
          t0 =  am_d2*mp*rate_dr(iCv)*u_moli(iCv,0)*geo%cvHz(iCv) 
          smo0_dr(iCv,0,is_a) =  t0
        enddo
      enddo


cc dissociative excitation, parallel momentum source
      do ism = 0, ns_mol-1
        is_a = mol_is_rp_atom(ism) 
        is_i = mol_is_rp_ion(ism)  
        do iCv = 1, mpg%nCi
          t0 =  am_d2*mp*rate_de(iCv)*u_moli(iCv,0)*geo%cvHz(iCv) 
          smo0_de(iCv,0,is_a) =  0.5_R8*t0
          smo0_de(iCv,0,is_i) =  0.5_R8*t0
        enddo
      enddo

cc dissociative ionization (2), parallel momentum source
      do ism = 0, ns_mol-1
        is_i = mol_is_rp_ion(ism)  
        do iCv = 1, mpg%nCi
          t0 =  am_d2*mp*rate_di2(iCv)*u_moli(iCv,0)*geo%cvHz(iCv) 
          smo0_di2(iCv,0,is_i) =  t0
        enddo
      enddo

* ..Ion / atom energy sources
!! TODO conversion total to internal !!!!

cc non-dissociative ionization
      if (switch%b2stfm_freac_nd.ge.0.0_R8) then
      !! creates a molecular ion, no direct energy transfer to ions or atoms
        do ism = 0, ns_mol-1
          do iCv = 1, mpg%nCi
            !! for the molecule itself 
            t0_mol = -rate_nd(iCv)*(1.5_R8*t_mol(iCv,ism) +         
     .          0.5_R8*vsq(iCv)*am_d2*mp)
            se_mol_nd(iCv) = t0_mol
            !! conversion total -> internal for molecule
            t0_mol = t0_mol - u_mol(iCv,0)*sm_mol_nd(iCv)/geo%cvHz(iCv)
     .         + 0.5_R8*am_d2*mp*sn_mol_nd(iCv)*u_mol(iCv,0)**2  
            sh_mol_nd(iCv) = t0_mol
          enddo
        enddo
      else
        se_mol_nd = 0.0_R8
        sh_mol_nd = 0.0_R8
      endif

cc dissociation, ion energy source
      if (switch%b2stfm_freac_d.ge.0.0_R8) then
        do ism = 0, ns_mol-1
          is_a = mol_is_rp_atom(ism) 
          is_i = mol_is_rp_ion(ism)  
          do iCv = 1, mpg%nCi
            !! total energy
            t0 =  rate_d(iCv)*(1.5_R8*t_mol(iCv,ism) +         
     .          0.5_R8*vsq(iCv)*am_d2*mp +
     .          switch%b2stfm_erp_d*ev)
            !! conversion total -> internal for atom
            t0 = t0 - pl%ua(iCv,0)*smo0_d(iCv,0,0)/geo%cvHz(iCv) +
     .         0.5_R8*am(0)*mp*sna0_d(iCv,0,0)*pl%ua(iCv,0)**2       
            shi0_d(iCv,0,is_a) = t0
          
            !! for the molecule itself !! TODO generalize reaction product energies to guarantee energy conversvation between ions, electrons and molecules
            t0_mol = -rate_d(iCv)*(1.5_R8*t_mol(iCv,ism) +         
     .          0.5_R8*vsq(iCv)*am_d2*mp)
            t0_mol = t0_mol
            se_mol_d(iCv) = t0_mol
            !! conversion total -> internal for molecule
            t0_mol = t0_mol - u_mol(iCv,0)*sm_mol_d(iCv)/geo%cvHz(iCv)+
     .         0.5_R8*am_d2*mp*sn_mol_d(iCv)*u_mol(iCv,0)**2  
            sh_mol_d(iCv) = t0_mol*fac_nlim(iCv,ism)
          enddo
        enddo
      else
        shi0_d = 0.0_R8
        se_mol_d = 0.0_R8
        sh_mol_d = 0.0_R8
      endif

cc dissociative ionization (1), ion energy source
      if (switch%b2stfm_freac_di1.ge.0.0_R8) then
        do ism = 0, ns_mol-1
          is_a = mol_is_rp_atom(ism) 
          is_i = mol_is_rp_ion(ism)  
          do iCv = 1, mpg%nCi
            !! total energy
            t0 =  rate_di1(iCv)*(1.5_R8*t_mol(iCv,ism) + 
     .          0.5_R8*vsq(iCv)*am_d2*mp + 
     .          switch%b2stfm_erp_di1*ev)
            !! split
            t0_i = 0.5_R8*t0 
            t0_a = 0.5_R8*t0
            !! conversion total -> internal
            t0_i = t0_i  - pl%ua(iCv,1)*smo0_di1(iCv,0,1)/geo%cvHz(iCv)
     .         + 0.5_R8*am(1)*mp*sna0_di1(iCv,0,1)*pl%ua(iCv,1)**2  
            t0_a = t0_a  - pl%ua(iCv,0)*smo0_di1(iCv,0,0)/geo%cvHz(iCv)
     .         + 0.5_R8*am(0)*mp*sna0_di1(iCv,0,0)*pl%ua(iCv,0)**2  
            shi0_di1(iCv,0,is_a) = t0_a 
            shi0_di1(iCv,0,is_i) = t0_i 

            !! for the molecule itself
            t0_mol =  -rate_di1(iCv)*(1.5_R8*t_mol(iCv,ism) + 
     .          0.5_R8*vsq(iCv)*am_d2*mp)
            t0_mol = t0_mol
            se_mol_di1(iCv) = t0_mol
            !! conversion total -> internal for molecule
            t0_mol = t0_mol -u_mol(iCv,0)*sm_mol_di1(iCv)/geo%cvHz(iCv)
     .         + 0.5_R8*am_d2*mp*sn_mol_di1(iCv)*u_mol(iCv,0)**2  
            sh_mol_di1(iCv) = t0_mol*fac_nlim(iCv,ism)
          enddo
        enddo
      else
        shi0_di1 = 0.0_R8
        se_mol_di1 = 0.0_R8
        sh_mol_di1 = 0.0_R8
      endif

cc ion conversion, ion energy source
      !! the ion becomes atom
      if (switch%b2stfm_freac_ic.ge.0.0_R8) then
        do ism = 0, ns_mol-1
          is_a = mol_is_rp_atom(ism) 
          is_i = mol_is_rp_ion(ism)  
          do iCv = 1, mpg%nCi
            !! total energy
            t0 =  rate_ic(iCv)*(1.5_R8*pl%ti(iCv) + 
     .          0.5_R8*pl%ua(iCv,1)**2*am(1)*mp)
            !! the ion becomes atom
            t0_i = -t0 
            t0_a =  t0
            !! conversion total -> internal
            t0_i = t0_i  - pl%ua(iCv,1)*smo0_ic(iCv,0,1)/geo%cvHz(iCv)
     .         + 0.5_R8*am(1)*mp*sna0_ic(iCv,0,1)*pl%ua(iCv,1)**2  
            t0_a = t0_a  - pl%ua(iCv,0)*smo0_ic(iCv,0,0)/geo%cvHz(iCv)
     .         + 0.5_R8*am(0)*mp*sna0_ic(iCv,0,0)*pl%ua(iCv,0)**2 
            shi0_ic(iCv,0,is_a) =  t0_a
            shi0_ic(iCv,0,is_i) =  t0_i

            !! for the molecule
            t0_mol = -rate_ic(iCv)*(1.5_R8*t_mol(iCv,0) + 
     .          0.5_R8*vsq(iCv)*am_d2*mp)
            t0_mol = t0_mol
            se_mol_ic(iCv) = t0_mol
            !! conversion total -> internal for molecule
            t0_mol = t0_mol - u_mol(iCv,0)*sm_mol_ic(iCv)/geo%cvHz(iCv)
     .         + 0.5_R8*am_d2*mp*sn_mol_ic(iCv)*u_mol(iCv,0)**2 
            sh_mol_ic(iCv) = t0_mol*fac_nlim(iCv,ism)
          enddo
        enddo
      else
        shi0_ic = 0.0_R8
        se_mol_ic = 0.0_R8
        sh_mol_ic = 0.0_R8
      endif

cc elastic collisions, ion energy source
      if (switch%b2stfm_freac_ec.ge.0.0_R8) then
        do ism = 0, ns_mol-1
          is_i = mol_is_rp_ion(ism)  
          do iCv = 1, mpg%nCi
            if (switch%b2stfm_debug_icv.eq.iCv) then !! option to simulate energy transfer due to elastic collisions in detail for one cell. 
             call integrate_elastic_collisions_cell(mp*am_d2, mp*am(1), 
     *        t_mol(iCv,ism)/ev, pl%ti(iCv)/ev, pl%ua(iCv,1), 0.0_R8,
     *        0.0_R8,u_mol(iCv,ism), u_rad_mol(iCv,ism),
     *        u_dia_mol(iCv,ism))
            endif
            !! total energy, thermal part
            t0 =  rate_ec_E(iCv)*1.5_R8*(t_mol(iCv,ism)-pl%ti(iCv))  
            !! total energy, contributions from mechanical work
            F_p = mr*rate_ec_m(iCv)*(u_mol(iCv,ism)-pl%ua(iCv,1)) !! force on ions, parallel
              !! center-of-mass parallel velocity
            Vcom_p = (am_d2*u_mol(iCv,ism) + am(1)*pl%ua(iCv,1))/
     *               (am_d2+am(1)) 
            if (switch%b2stfm_ekin_perp.eq.0) then
              F_d = 0.0_R8
              F_r = 0.0_R8
              Vcom_d = 0.0_R8
              Vcom_r = 0.0_R8
            else !! take also into account perpendicular mechanical work
              F_d = mr*rate_ec_m(iCv)*u_dia_mol(iCv,ism) !! force on ions, diamagnetic 
              F_r = mr*rate_ec_m(iCv)*u_rad_mol(iCv,ism) !! force on ions, radial
              Vcom_d = (am_d2*u_dia_mol(iCv,ism))/(am_d2+am(1))  !! center-of-mass diamagnetic velocity, neglecting ion contribution
              Vcom_r = (am_d2*u_rad_mol(iCv,ism))/(am_d2+am(1))  !! center-of-mass radial velocity, neglecting ion contribution
            endif
            wrk(iCv) = F_p*Vcom_p + F_d*Vcom_d + F_r*Vcom_r
            !! add mechanical work to thermal part
            t0 = t0 + wrk(iCv)
            t0_mol = -t0 !! save for step below
            !! conversion total -> internal (no particle source in this case)
            t0 = t0 - pl%ua(iCv,1)*smo0_ec(iCv,0,1)/geo%cvHz(iCv) 
            shi0_ec(iCv,0,is_i) = t0


            !! for the molecule
            se_mol_ec(iCv) = t0_mol
            !! conversion total -> internal for molecule
            t0_mol = t0_mol - u_mol(iCv,0)*sm_mol_ec(iCv)/geo%cvHz(iCv) 
            sh_mol_ec(iCv) = t0_mol
          enddo
        enddo
      else
        shi0_ec = 0.0_R8
        se_mol_ec = 0.0_R8
        sh_mol_ec = 0.0_R8
      endif

cc dissociative recombination, ion energy source
      if (switch%b2stfm_freac_dr.ge.0.0_R8) then
        do ism = 0, ns_mol-1
          is_a = mol_is_rp_atom(ism) 
          is_i = mol_is_rp_ion(ism)  
          do iCv = 1, mpg%nCi
            !! total energy
            t0 =  rate_dr(iCv)*(1.5_R8*t_moli(iCv,ism) + 
     .          0.5_R8*vsqi(iCv)*am_d2*mp + 
     .          switch%b2stfm_erp_dr*ev)
            !! special treatment, see H8 AMJUEL
            s0 = -1.681368547011e+01_R8+
     .          3.964355004318e-01_R8*log(pl.te(iCv)/ev)
            s0 = exp(s0)/1.0e6_R8*geo%cvVol(iCv)*
     .         dv%ne(iCv)*n_moli(iCv,ism)*ev
            t0 = t0 + s0
            !! conversion total -> internal
            t0 = t0 - pl%ua(iCv,0)*smo0_dr(iCv,0,0)/geo%cvHz(iCv) +
     .         0.5_R8*am(0)*mp*sna0_dr(iCv,0,0)*pl%ua(iCv,0)**2            
            shi0_dr(iCv,0,is_a) =  t0
          enddo
        enddo
      else
        shi0_dr = 0.0_R8
      endif


cc dissociative excitation, ion energy source
      if (switch%b2stfm_freac_dr.ge.0.0_R8) then
        do ism = 0, ns_mol-1
          is_a = mol_is_rp_atom(ism) 
          is_i = mol_is_rp_ion(ism)  
          do iCv = 1, mpg%nCi
            !! total energy
            t0 =  rate_de(iCv)*(1.5_R8*t_moli(iCv,ism) + 
     .          0.5_R8*vsqi(iCv)*am_d2*mp + 
     .          switch%b2stfm_erp_de*ev)
            !! split
            t0_i = 0.5_R8*t0 
            t0_a = 0.5_R8*t0
            !! conversion total -> internal
            t0_i = t0_i  - pl%ua(iCv,1)*smo0_de(iCv,0,1)/geo%cvHz(iCv) +
     .         0.5_R8*am(1)*mp*sna0_de(iCv,0,1)*pl%ua(iCv,1)**2  
            t0_a = t0_a  - pl%ua(iCv,0)*smo0_de(iCv,0,0)/geo%cvHz(iCv) +
     .         0.5_R8*am(0)*mp*sna0_de(iCv,0,0)*pl%ua(iCv,0)**2  
            shi0_de(iCv,0,is_a) = t0_a 
            shi0_de(iCv,0,is_i) = t0_i 
          enddo
        enddo
      else
        shi0_de = 0.0_R8
      endif

cc dissociative ionization (2), ion energy source
      if (switch%b2stfm_freac_di2.ge.0.0_R8) then
        do ism = 0, ns_mol-1
          is_i = mol_is_rp_ion(ism)  
          do iCv = 1, mpg%nCi
            !! total energy
            t0 =  rate_di2(iCv)*(1.5_R8*t_moli(iCv,ism) + 
     .          0.5_R8*vsq(iCv)*am_d2*mp + 
     .          switch%b2stfm_erp_di2*ev)
            !! conversion total -> internal
            t0 = t0 - pl%ua(iCv,0)*smo0_di2(iCv,0,1)/geo%cvHz(iCv) +
     .         0.5_R8*am(1)*mp*sna0_di2(iCv,0,1)*pl%ua(iCv,1)**2            
            shi0_di2(iCv,0,is_i) =  t0
          enddo
        enddo
      else
        shi0_di2 = 0.0_R8
      endif


* ..Electron energy sources

cc dissociation , electron energy source
      do ism = 0, ns_mol-1
        do iCv = 1, mpg%nCi
          t0 =  rate_d(iCv)*10.5_R8*ev
          she0_d(iCv,0) = -t0
        enddo
      enddo

cc non-dissociative ionization , electron energy source
      do ism = 0, ns_mol-1
        do iCv = 1, mpg%nCi
          t0 =  rate_nd(iCv)*15.4_R8*ev
          she0_nd(iCv,0) = -t0
        enddo
      enddo

cc dissociative ionization (1), electron energy source
      do ism = 0, ns_mol-1
        do iCv = 1, mpg%nCi
          t0 =  rate_di1(iCv)*25.0_R8*ev
          she0_di1(iCv,0) = -t0
        enddo
      enddo

cc ion conversion
  !! no electron energy sink ? -> seems not

cc elastic collisions
  !! no electron energy sink

cc dissociative recombination, electron energy source
      do ism = 0, ns_mol-1
        do iCv = 1, mpg%nCi
          !! special treatment, see H8 AMJUEL
          s0 = -1.681368547011e+01_R8+
     .          3.964355004318e-01_R8*log(pl.te(iCv)/ev)
          s0 = exp(s0)/1.0e6_R8*geo%cvVol(iCv)*
     .         dv%ne(iCv)*n_moli(iCv,ism)*ev
          t0 = s0
          she0_dr(iCv,0) = -t0
        enddo
      enddo


cc dissociative excitation, electron energy source
      do ism = 0, ns_mol-1
        do iCv = 1, mpg%nCi
          t0 =  rate_de(iCv)*10.4_R8*ev
          she0_de(iCv,0) = -t0
        enddo
      enddo

cc dissociative ionization (2), electron energy source
      do ism = 0, ns_mol-1
        do iCv = 1, mpg%nCi
          t0 =  rate_di2(iCv)*15.5_R8*ev
          she0_di2(iCv,0) = -t0
        enddo
      enddo

* ..Sum over all reactions

!! for the plasma and fluid atoms
      do is = 0, ns-1
        !! particles
        sna0_tot(:,0,is) = sna0_nd(:,0,is) + sna0_d(:,0,is) + 
     .                     sna0_di1(:,0,is) + sna0_ic(:,0,is) +
     .                     sna0_dr(:,0,is) + sna0_de(:,0,is) +
     .                     sna0_di2(:,0,is) 
        !! parallel momentum
        smo0_tot(:,0,is) = smo0_nd(:,0,is) + smo0_d(:,0,is) + 
     .                     smo0_di1(:,0,is) + smo0_ic(:,0,is) +
     .                     smo0_dr(:,0,is) + smo0_de(:,0,is) +
     .                     smo0_di2(:,0,is) + smo0_ec(:,0,is)
        !! ion / atom energy
        shi0_tot(:,0,is) = shi0_nd(:,0,is) + shi0_d(:,0,is) + 
     .                     shi0_di1(:,0,is) + shi0_ic(:,0,is) +
     .                     shi0_dr(:,0,is) + shi0_de(:,0,is) +
     .                     shi0_di2(:,0,is) + shi0_ec(:,0,is)
      enddo 
      she0_tot(:,0) = she0_nd(:,0) + she0_d(:,0) + 
     .                she0_di1(:,0) + she0_ic(:,0) +
     .                she0_dr(:,0) + she0_de(:,0) +
     .                she0_di2(:,0)

!! for the molecules themselves
      sn_mol_tot = sn_mol_d + sn_mol_nd + sn_mol_di1 + sn_mol_ic
      sm_mol_tot = sm_mol_d + sm_mol_nd + sm_mol_di1 + sm_mol_ic +
     .             sm_mol_ec
      sh_mol_tot = sh_mol_d + sh_mol_nd + sh_mol_di1 + sh_mol_ic +
     .             sh_mol_ec
      se_mol_tot = se_mol_d + se_mol_nd + se_mol_di1 + se_mol_ic +
     .             se_mol_ec
      sm_mol_rad_tot = sm_mol_rad_d + sm_mol_rad_nd + sm_mol_rad_di1 + 
     .             sm_mol_rad_ic + sm_mol_rad_ec
      sm_mol_dia_tot = sm_mol_dia_d + sm_mol_dia_nd + sm_mol_dia_di1 + 
     .             sm_mol_dia_ic + sm_mol_dia_ec

      total_mol_sink = sum(sn_mol_tot)
      write(*,*) 'total volumetric sink of molecules'
      write(*,*) total_mol_sink

!! Fill in general source structure
      sr_mol%sn_mol(:,0) = sr_mol%sn_mol(:,0) + 0.5_R8*sn_mol_tot
      sr_mol%sm_mol(:,0) = sr_mol%sm_mol(:,0) + 0.5_R8*sm_mol_tot
      sr_mol%sh_mol(:,0) = sr_mol%sh_mol(:,0) + 0.5_R8*sh_mol_tot
      sr_mol%sm_mol_rad(:,0) = sr_mol%sm_mol_rad(:,0) +  !! here everytinng in 1?
     .                         sm_mol_rad_tot(:,0)
      sr_mol%sm_mol_dia(:,0) = sr_mol%sm_mol_dia(:,0) + 
     .                         sm_mol_dia_tot(:,0)

      sr_mol%sn_mol(:,1) = sr_mol%sn_mol(:,1) +
     .                     0.5_R8*sn_mol_tot/n_mol(:,0)
      sr_mol%sm_mol(:,1) = sr_mol%sm_mol(:,1) + 
     .                     0.5_R8*sm_mol_tot/(u_mol(:,0)+eps)
      sr_mol%sh_mol(:,1) = sr_mol%sh_mol(:,1) + 
     .                     0.5_R8*sh_mol_tot/t_mol(:,0)
      sr_mol%sm_mol_rad(:,1) = sr_mol%sm_mol_rad(:,1) + 
     .                     sm_mol_rad_tot(:,1)
      sr_mol%sm_mol_dia(:,1) = sr_mol%sm_mol_dia(:,1) + 
     .                     sm_mol_dia_tot(:,1)


!! Calculate time constant tau for use in the transport model
      do iCv = 1, mpg%nCi
        tau_inv = 0.0_R8
        t0 =  pl%na(iCv,1)*rtw%r_ec_m(iCv,0)*switch%b2stfm_freac_ec 
        tau_inv = tau_inv + t0
        t0 =  dv%ne(iCv)*rtw%r_nd(iCv,0)*switch%b2stfm_freac_nd
        tau_inv = tau_inv + t0
        t0 =  dv%ne(iCv)*rtw%r_d(iCv,0)*switch%b2stfm_freac_d
        tau_inv = tau_inv + t0
        t0 =  dv%ne(iCv)*rtw%r_di1(iCv,0)*switch%b2stfm_freac_di1
        tau_inv = tau_inv + t0
        t0 =  pl%na(iCv,1)*rtw%r_ic(iCv,0)*switch%b2stfm_freac_ic
        tau_inv = tau_inv + t0
        tau_mol(iCv) = 1.0_R8/tau_inv
      enddo

!! Calculate the pressure-diffusion coefficient
      do iCv = 1, mpg%nCv ! -> these also in guard cells, for computation in guard face conductivities
        vth = sqrt(t_mol(iCv,0)/m_D2)
        t0 = 0.0_R8
        t0 =  t0 + mr*pl%na(iCv,1)*rtw%r_ec_m(iCv,0)*   !! mr instead of m_D2 here!
     .                       switch%b2stfm_freac_ec
        mfp(iCV) = vth/
     .      (pl%na(iCv,1)*rtw%r_ec_m(iCv,0)*  
     .                       switch%b2stfm_freac_ec)     !! mean free path
        t0 =  t0 + m_D2*dv%ne(iCv)*rtw%r_nd(iCv,0)*
     .                       switch%b2stfm_freac_nd
        t0 =  t0 + m_D2*dv%ne(iCv)*rtw%r_d(iCv,0)*
     .                       switch%b2stfm_freac_d
        t0 =  t0 + m_D2*dv%ne(iCv)*rtw%r_di1(iCv,0)*
     .                       switch%b2stfm_freac_di1
        t0 =  t0 + m_D2*pl%na(iCv,1)*rtw%r_ic(iCv,0)*
     .                       switch%b2stfm_freac_ic
        diff_coeff(iCv) = 1.0_R8/t0
        coeff_vi(iCv) = mr*pl%na(iCv,1)*rtw%r_ec_m(iCv,0)* 
     .               switch%b2stfm_freac_ec/t0
        if (switch%b2stfm_diff_coeff_ene_style.ne.0) then
          t0 = t0 - mr*pl%na(iCv,1)*rtw%r_ec_m(iCv,0)*   
     .                       switch%b2stfm_freac_ec
          t0 = t0 + mr*pl%na(iCv,1)*rtw%r_ec_E(iCv,0)*   
     .                       switch%b2stfm_freac_ec
          diff_coeff_ene(iCv) = 1.0_R8/t0
        else
          diff_coeff_ene(iCv) = diff_coeff(iCv)
        endif
      enddo

!! Calculate the equilibrium density factor needed in the pressure-diffusion equation
      do iCv = 1, mpg%nCv  ! -> these also in guard cells, for computation in guard face conductivities
            beta_num = 0.0_R8        !! numerator
            beta_denom = 0.0_R8      !! denominator
            t0 =  dv%ne(iCv)*rtw%r_nd(iCv,0)*switch%b2stfm_freac_nd
            beta_denom = beta_denom + t0
            t0 =  dv%ne(iCv)*rtw%r_d(iCv,0)*switch%b2stfm_freac_d
            beta_denom = beta_denom + t0
            t0 =  dv%ne(iCv)*rtw%r_di1(iCv,0)*switch%b2stfm_freac_di1
            beta_denom = beta_denom + t0
            t0 = pl%na(iCv,1)*rtw%r_ec_m(iCv,0)*switch%b2stfm_freac_ec
            beta_denom = beta_denom + t0
            beta_num = beta_num + t0
            t0 =  pl%na(iCv,1)*rtw%r_ic(iCv,0)*switch%b2stfm_freac_ic
            beta_denom = beta_denom + t0
            beta_num = beta_num + t0
            beta_eq(iCv) = beta_num/beta_denom
      enddo

* ..Add to main source contributions
      do is = 0, ns-1
        call b2saxpy (nCv, 1.0_R8, sna0_tot(:,0,is), 1,
     .                             sr%sna(:,0,is), 1)
        call b2saxpy (nCv, 1.0_R8, smo0_tot(:,0,is), 1,
     .                             sr%smq(:,0,is), 1)
        if (switch%tn_style.eq.0) then
          call b2saxpy (nCv, 1.0_R8, shi0_tot(:,0,is), 1,
     .                             sr%shi(:,0), 1)
        elseif (switch%tn_style.eq.1) then 
          if (.not.is_neutral(is)) then 
            call b2saxpy (nCv, 1.0_R8, shi0_tot(:,0,is), 1,
     .                             sr%shi(:,0), 1)
          endif
        else
          if (.not.is_neutral(is)) then 
            call b2saxpy (nCv, 1.0_R8, shi0_tot(:,0,is), 1,
     .                             sr%shi(:,0), 1)
          else
            call b2saxpy (nCv, 1.0_R8, shi0_tot(:,0,is), 1,
     .                             sr%shn(:,0), 1)
          endif
        endif
      enddo
      call b2saxpy (nCv, 1.0_R8,she0_tot(:,0) , 1, sr%she(:,0), 1)



* ..Write output

      if (switch%b2stfm_iout.eq.1) then    
        call my_out_us(70,nCv,0,n_mol(:,0),'b2stfm_n_mol')     
        call my_out_us(70,nCv,0,u_mol(:,0),'b2stfm_u_mol')                        
        call my_out_us(70,nCv,0,t_mol(:,0),'b2stfm_t_mol')                        
        call my_out_us(70,nCv,0,n_moli(:,0),'b2stfm_n_moli')                        
        call my_out_us(70,nCv,0,u_moli(:,0),'b2stfm_u_moli')                        
        call my_out_us(70,nCv,0,t_moli(:,0),'b2stfm_t_moli')  
        call my_out_us(70,nCv,0,nmi_new(:,0),'b2stfm_nmi_new')  
        call my_out_us(70,nCv,0,mfp,'b2stfm_mfp')            
        call my_out_us(70,nCv,0,fac_nlim,'b2stfm_fac_nlim')            

                                                                                                        
        do is = 0, ns-1
          write (chns,'(i3.3)') is
cc Particles
          call my_out_us(70,nCv,0,sna0_nd(:,0,is),
     .                            'b2stfm_sna0_nd_'//chns)
          call my_out_us(70,nCv,0,sna0_d(:,0,is),
     .                            'b2stfm_sna0_d_'//chns)
          call my_out_us(70,nCv,0,sna0_di1(:,0,is),
     .                            'b2stfm_sna0_di1_'//chns)
          call my_out_us(70,nCv,0,sna0_ic(:,0,is),
     .                            'b2stfm_sna0_ic_'//chns)
          call my_out_us(70,nCv,0,sna0_dr(:,0,is),
     .                            'b2stfm_sna0_dr_'//chns)
          call my_out_us(70,nCv,0,sna0_de(:,0,is),
     .                            'b2stfm_sna0_de_'//chns)
          call my_out_us(70,nCv,0,sna0_di2(:,0,is),
     .                            'b2stfm_sna0_di2_'//chns)
          call my_out_us(70,nCv,0,sna0_tot(:,0,is),
     .                            'b2stfm_sna0_tot_'//chns)
cc Parallel momentum
          call my_out_us(70,nCv,0,smo0_nd(:,0,is),
     .                            'b2stfm_smo0_nd_'//chns)
          call my_out_us(70,nCv,0,smo0_d(:,0,is),
     .                            'b2stfm_smo0_d_'//chns)
          call my_out_us(70,nCv,0,smo0_di1(:,0,is),
     .                            'b2stfm_smo0_di1_'//chns)
          call my_out_us(70,nCv,0,smo0_ic(:,0,is),
     .                            'b2stfm_smo0_ic_'//chns)
          call my_out_us(70,nCv,0,smo0_ec(:,0,is),
     .                            'b2stfm_smo0_ec_'//chns)
          call my_out_us(70,nCv,0,smo0_dr(:,0,is),
     .                            'b2stfm_smo0_dr_'//chns)
          call my_out_us(70,nCv,0,smo0_de(:,0,is),
     .                            'b2stfm_smo0_de_'//chns)
          call my_out_us(70,nCv,0,smo0_di2(:,0,is),
     .                            'b2stfm_smo0_di2_'//chns)
          call my_out_us(70,nCv,0,smo0_tot(:,0,is),
     .                            'b2stfm_smo0_tot_'//chns)
cc Ion/atom energy
          call my_out_us(70,nCv,0,shi0_nd(:,0,is),
     .                            'b2stfm_shi0_nd_'//chns)
          call my_out_us(70,nCv,0,shi0_d(:,0,is),
     .                            'b2stfm_shi0_d_'//chns)
          call my_out_us(70,nCv,0,shi0_di1(:,0,is),
     .                            'b2stfm_shi0_di1_'//chns)
          call my_out_us(70,nCv,0,shi0_ic(:,0,is),
     .                            'b2stfm_shi0_ic_'//chns)
          call my_out_us(70,nCv,0,shi0_ec(:,0,is),
     .                            'b2stfm_shi0_ec_'//chns)
          call my_out_us(70,nCv,0,shi0_dr(:,0,is),
     .                            'b2stfm_shi0_dr_'//chns)
          call my_out_us(70,nCv,0,shi0_de(:,0,is),
     .                            'b2stfm_shi0_de_'//chns)
          call my_out_us(70,nCv,0,shi0_di2(:,0,is),
     .                            'b2stfm_shi0_di2_'//chns)
          call my_out_us(70,nCv,0,shi0_tot(:,0,is),
     .                            'b2stfm_shi0_tot_'//chns)
        enddo
cc Electron energy
        call my_out_us(70,nCv,0,she0_nd(:,0),
     .                            'b2stfm_she0_nd')
        call my_out_us(70,nCv,0,she0_d(:,0),
     .                            'b2stfm_she0_d')
        call my_out_us(70,nCv,0,she0_di1(:,0),
     .                            'b2stfm_she0_di1')
        call my_out_us(70,nCv,0,she0_ic(:,0),
     .                            'b2stfm_she0_ic')
        call my_out_us(70,nCv,0,she0_dr(:,0),
     .                            'b2stfm_she0_dr')
        call my_out_us(70,nCv,0,she0_de(:,0),
     .                            'b2stfm_she0_de')
        call my_out_us(70,nCv,0,she0_di2(:,0),
     .                            'b2stfm_she0_di2')
        call my_out_us(70,nCv,0,she0_tot(:,0),
     .                            'b2stfm_she0_tot')

cc Molecules 
        call my_out_us(70,nCv,0,sn_mol_tot(:),
     .                            'b2stfm_sn_mol')
        call my_out_us(70,nCv,0,sm_mol_tot(:),
     .                            'b2stfm_sm_mol')
        call my_out_us(70,nCv,0,sh_mol_tot(:),
     .                            'b2stfm_sh_mol')
        call my_out_us(70,nCv,0,se_mol_tot(:),
     .                            'b2stfm_se_mol')
        call my_out_us(70,nCv,0,sm_mol_rad_tot(:,0),
     .                            'b2stfm_sm_rad0_mol')
        call my_out_us(70,nCv,0,sm_mol_rad_tot(:,1),
     .                            'b2stfm_sm_rad1_mol')
        call my_out_us(70,nCv,0,sm_mol_dia_tot(:,0),
     .                            'b2stfm_sm_dia0_mol')
        call my_out_us(70,nCv,0,sm_mol_dia_tot(:,1),
     .                            'b2stfm_sm_dia1_mol')
c for debugging, maybe remove later 
        call my_out_us(70,nCv,0,sh_mol_d(:),
     .                            'b2stfm_sh_mol_d')
        call my_out_us(70,nCv,0,sh_mol_nd(:),
     .                            'b2stfm_sh_mol_nd')
        call my_out_us(70,nCv,0,sh_mol_di1(:),
     .                            'b2stfm_sh_mol_di1')
        call my_out_us(70,nCv,0,sh_mol_ic(:),
     .                            'b2stfm_sh_mol_ic')
        call my_out_us(70,nCv,0,sh_mol_ec(:),
     .                            'b2stfm_sh_mol_ec')
        call my_out_us(70,nCv,0,se_mol_d(:),
     .                            'b2stfm_se_mol_d')
        call my_out_us(70,nCv,0,se_mol_nd(:),
     .                            'b2stfm_se_mol_nd')
        call my_out_us(70,nCv,0,se_mol_di1(:),
     .                            'b2stfm_se_mol_di1')
        call my_out_us(70,nCv,0,se_mol_ic(:),
     .                            'b2stfm_se_mol_ic')
        call my_out_us(70,nCv,0,se_mol_ec(:),
     .                            'b2stfm_se_mol_ec')
        call my_out_us(70,nCv,0,wrk(:),
     .                            'b2stfm_wrk')

c rates
        call my_out_us(70,nCv,0,rate_d,
     .                            'b2stfm_rate_d')
        call my_out_us(70,nCv,0,rate_nd,
     .                            'b2stfm_rate_nd')
        call my_out_us(70,nCv,0,rate_di1,
     .                            'b2stfm_rate_di1')
        call my_out_us(70,nCv,0,rate_ic,
     .                            'b2stfm_rate_ic')
        call my_out_us(70,nCv,0,rate_ec_m,
     .                            'b2stfm_rate_ec_m')
        call my_out_us(70,nCv,0,rate_dr,
     .                            'b2stfm_rate_dr')
        call my_out_us(70,nCv,0,rate_de,
     .                            'b2stfm_rate_de')
        call my_out_us(70,nCv,0,rate_di2,
     .                            'b2stfm_rate_di2')

c write out molecule sinks per reactions
      total_mol_sink = sum(sn_mol_tot)
      write(*,*) 'total volumetric sink of molecules'
      write(*,*) total_mol_sink
      t0 = sum(sn_mol_d)
      write(*,*) 'molecules sink due to d', t0
      t0 = sum(sn_mol_nd)
      write(*,*) 'molecules sink due to nd', t0
      t0 = sum(sn_mol_di1)
      write(*,*) 'molecules sink due to di1', t0
      t0 = sum(sn_mol_ic)
      write(*,*) 'molecules sink due to ic', t0

      endif
  

* ..return
      call subend ()
      return

      end subroutine b2stfm
