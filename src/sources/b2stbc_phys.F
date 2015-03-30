*-----------------------------------------------------------------------
*.specification
      subroutine b2stbc_phys (nx, ny, ns, ismain, ismain0,
     1  BoRiS, facdrift, fac_ExB,
     1  vol, crx, hx, hy, hz, qz, qc, gs, pbs, pbshz, bb, na, ua,     !srv 24.08.08
     1  uadia, vedia, vadia, veecrb, vaecrb,
     &  te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc, fna, fna_mdf, !srv 28.10.02
     1  fnaPSch, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,      !srv 13.10.06
     &  fhe, fhi, fhe_mdf, fhi_mdf, fna_fcor, fhePSch, fhiPSch,
     &  fch, fne, fni,                                                !srv 13.10.06
     &  fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,!srv 21.10.02
     2  cdna, chce, chve, chci, chvi, calf,                           !srv 25.09.09
     &  sna0, smo0, she0, shi0, sch0, sne0,
     &  wrong_flow, main_call, diagno,
     &  b2stbc_fchy_dia, b2stbc_neoclassical, praverage,
     1  b2stbc_fchy_dia_coreonly, secmodel, area_fix,
     2  sheath_drift_fix)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_anomalous_transport
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_time
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_external
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif
#ifdef PAR_OPT
      use b2mod_par_opt
#endif
      use b2mod_astra_to_b2                                          !srv 28.07.11
      implicit none
#     include "DIMENSIONS.F"
#ifndef B25_EIRENE
#     include "COUPLE/EIRDIAG.F"
#endif
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, ismain0, diagno
      real (kind=R8) ::
     *  BoRiS, facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny),
     &  vol(-1:nx,-1:ny), crx(-1:nx,-1:ny,0:3),
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), hz(-1:nx,-1:ny),           !srv 24.08.08
     &  qz(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1),
     &  pbshz(-1:nx,-1:ny,0:1),
     1  bb(-1:nx,-1:ny,0:3), gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     1  uadia(-1:nx,-1:ny,0:1,0:ns-1),
     &  vedia(-1:nx,-1:ny,0:1), vadia(-1:nx,-1:ny,0:1,0:ns-1),
     2  veecrb(-1:nx,-1:ny,0:1), vaecrb(-1:nx,-1:ny,0:1,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  kinrgy(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1), fna_32(-1:nx,-1:ny,0:1,0:ns-1),
     *  fna_mdf(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_52(-1:nx,-1:ny,0:1,0:ns-1), fni_32(-1:nx,-1:ny,0:1),
     &  fni_52(-1:nx,-1:ny,0:1), fne_32(-1:nx,-1:ny,0:1),
     &  fne_52(-1:nx,-1:ny,0:1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fhe_mdf(-1:nx,-1:ny,0:1),
     &  fhi_mdf(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     &  fch_p(-1:nx,-1:ny,0:1), fchvispar(-1:nx,-1:ny,0:1),           !srv 11.10.02 21.10.02 {
     &  fchvisper(-1:nx,-1:ny,0:1), fchin(-1:nx,-1:ny,0:1),
     &  fchvisq(-1:nx,-1:ny,0:1), fchdia(-1:nx,-1:ny,0:1),
     &  fchinert(-1:nx,-1:ny,0:1),                                    !srv 21.10.02 }
     &  fne(-1:nx,-1:ny,0:1), fni(-1:nx,-1:ny,0:1),
     &  fna_fcor(-1:nx,-1:ny,0:1,0:ns-1),                            !srv 25.02.00
     &  floe_noc(-1:nx,-1:ny,0:1), floi_noc(-1:nx,-1:ny,0:1),        !srv 28.10.02
     &  fhePSch(-1:nx,-1:ny,0:1), fhiPSch(-1:nx,-1:ny,0:1),          !srv 13.10.06
     &  fnaPSch(-1:nx,-1:ny,0:1,0:ns-1)                              !srv 13.10.06
      real (kind=R8) ::
     1  chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),
     2  chci(-1:nx,-1:ny,0:1), chvi(-1:nx,-1:ny,0:1),
     3  calf(-1:nx,-1:ny,0:1), cdna(-1:nx,-1:ny,0:1,0:ns-1)          !srv 25.09.09
      real (kind=R8) ::
     &  b2stbc_fchy_dia, b2stbc_neoclassical, praverage(2)
      logical main_call                                                  !xpb
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3), sne0(-1:nx,-1:ny,0:1)
      integer b2stbc_fchy_dia_coreonly, secmodel, area_fix,
     &  sheath_drift_fix
      logical wrong_flow                                                 !xpb
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBC_PHYS provides linearised expressions for the sources of
*     density, parallel momentum, electron heat, atom heat, and charge,
*     due to the boundary conditions.
*     Here the boundary conditions are expressed by use of the old B2
*     namelist format.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, ib, ireg, coreregno, coreregn2,
     1  bc_type13_ref, bc_type13_ix, bc_type13_iy, bc_type13_i,
     2  bc_type13_ixx,
     1  bc_type16_ref, bc_type16_ix, bc_type16_iy, bc_type16_i,
     2  bc_type16_ixx,
     1  bc_type20_ref, bc_type20_ix, bc_type20_iy, bc_type20_i,
     2  bc_type20_ixx,
     1  bc_type21_ref, bc_type21_ix, bc_type21_iy, bc_type21_i,
     2  bc_type21_ixx,
     1  bc_type16_kinetic_energy
      integer       :: num_inner, k, isb, isa, ism, ibw                 !srv 26.11.02
      integer, save :: istyle_average_pr = 0, istyle_fchi = 0           !srv 08.08.02 14.02.05
      integer, save :: iform_corr_dia = 0                               !srv 25.07.05
      integer, save :: mdf_fnb = 0, mdf_fhe = 0, mdf_fhi = 0, iout = 0, !srv 13.10.06
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0, iout_b2npmo =0
      integer, save :: istyle_cur_contr_on_S_and_N = 0                  !srv 01.02.09
      integer, save :: bcene_16_style =0, bceni_16_style = 0            !srv 16.02.12
      integer, save :: bcene_15_style =0                                !srv 19.04.13
      integer, save :: pot_eq = 1                                       !srv 27.02.13
      integer, save :: use_astra
      integer, save :: drift_style
      logical once_bccon, once_bcene, once_bceni                        !srv 16.02.12
      real (kind=R8) , save :: cbc = 0.0e0_R8, bc_ref = 0.01            !srv 14.10.99
      real (kind=R8) , save :: bc_ref_te = 0.01, bc_ref_ti=0.01         !IYS 07.03.2013 11.04.13
      real (kind=R8) , save :: senepar, senipar, sconpar                !srv 16.02.12
      real (kind=R8) ::
     *  cutll, cutlo, vte, s0, s1, s0hz, s1hz, fche, fchi, t0, t1, t2,      !srv 15.11.07
     1  t3, totflux, weight,
     1  naave, volave, teave, neave, tiave, niave, poave,
     2  vwest, veast, totfluxn, totfluxi
     &  ,pzaverage, pzsum, volsum                                     !srv 24.03.99 18.01.02
     &  ,ua2sum, ua2average(0:ns-1)                                   !srv 07.04.00 15.02.02
     &  ,uasum, uaaverage(0:ns-1)                                     !srv 15.02.02
     &  ,nesum, neaverage, tisum, tiaverage, tesum, teaverage         !srv 15.12.00 25.01.02
     &  ,nasum, naaverage(0:ns-1), bbsum, bbaverage                   !srv 22.01.02 18.10.04
      real (kind=R8) :: bzhzsum, bzhzaverage                          !lk 20.11.07
      real (kind=R8) ::
     *  pr(-1:nx,-1:ny), wdia(-1:nx,-1:ny),
     *  uu(-1:nx,-1:ny,0:ns+ns_ext-1),
     *  sx(-1:nx,-1:ny), ud(-1:nx,-1:ny,0:ns-1),
     *  sxhz(-1:nx,-1:ny), pz(-1:nx,-1:ny),         !srv 24.08.08
     &  sna0_no_mdf(-1:nx,-1:ny,0:1,0:ns-1)                          !srv 23.09.08
      real (kind=R8), save ::                                        !srv 16.02.12
     *  NUEOUT, UBIG, VBIG, PBIG, UP00, EPSLON, NE00, S2, US, FS, VS,
     *  NAS, TES, TIS, NIVS, FFAC, DFS
      real (kind=R8) , save :: phm0 = 1.0e0_R8, phm1 = 1.0e0_R8,
     &                         phm2 = 1.0e0_R8,
     &                         corr_flux   = 0.0e0_R8                 !srv 25.07.05
      character chns*3, chk*1                                        !srv 11.07.05
      logical done, decay_length_ok
      real (kind=R8) ::
     *  ro(-1:nx,-1:ny), cs, csb, uuout, upout, zbar, zarg, uutest,  !srv 01.02.09
     *  fchy_dia, fchy_inertia, seec                                 !srv 07.04.00
      integer isfeedback, fix_fch_in_fhe_sheath
      integer, save :: fix_sx_in_fhe_sheath = 1                       !srv 08.05.09
      integer, save :: use_eirene
      real (kind=R8) :: nesepm, nesepm_overshoot,
     1 bc_type13_norm, bc_type13_fac, na0eps
      real (kind=R8) :: fch_anomalous
      character*26 boundary_location(nbcd)
      character*24 hlp_frm
      save ncall, coreregno, coreregn2,  bc_type13_ref,
     1  bc_type16_ref, bc_type20_ref, bc_type13_norm, bc_type13_fac,
     2  na0eps, fch_anomalous, bc_type21_ref,
     *  isfeedback, nesepm, nesepm_overshoot, fix_fch_in_fhe_sheath
*   ..procedures
      intrinsic min, max, abs, sqrt, exp, log, dim, mod
      real (kind=R8) ::
     *  machsfr, expu2, hy1, sy, bcint, syhz                   !srv 24.08.08
      external subini, subend, xertst, ipgeti, machsfr, bcint, eseec0
      external ipgetr, b2xvsg, b2xvff, b2xvps
      expu2(t0) = dim(exp(min(20.0_R8,max(t0,cutll))),cutlo)
*     (expu2 is the exponential function with underflow set to 0 and
*     large values limited at a somewhat arbitrary level.)
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
      sy(ix,iy)=gs(topix(ix,iy),topiy(ix,iy),1)
      syhz(ix,iy)=gs(topix(ix,iy),topiy(ix,iy),1)*                    !srv 24.08.08
     &            hz(ix,iy)
!     &            hz(topix(ix,iy),topiy(ix,iy))
*   ..initialisation
      data ncall/0/, coreregno/1/, coreregn2/4/, bc_type13_ref/1/,
     1 bc_type16_ref/1/, bc_type20_ref/1/, na0eps/1.0e-30_R8/,
     * isfeedback/0/, nesepm/0.0_R8/, nesepm_overshoot/0.0_R8/,
     5 fix_fch_in_fhe_sheath/2/, bc_type16_kinetic_energy/0/,
     & bc_type13_norm/1.0e15_R8/, bc_type13_fac/1.0_R8/,
     & fch_anomalous/1.0_R8/, bc_type21_ref/1/,
     & use_astra/0/, use_eirene/0/, drift_style/0/
!! FIXME: Maybe change default value of drift_style ???
#ifdef WANT_THIS
      real (kind=R8) ::
     .  Tguess, fy, ehy, fylimt, prefactor,
     .  fyt, gam0, gam1, tlimit, cone0, coni0, flo, a,
     .  cum_flo_vadia, cum_flo_vaecrb, flo_vedia, flo_veecrb
      real (kind=R8) :: hybr, harm, pit, upwind
      integer ibnd, no_current, conduction_only, no_hybr_e, no_hybr_i
      data no_current/0/, conduction_only/0/, no_hybr_e/0/, no_hybr_i/0/
      save no_current, conduction_only, no_hybr_e, no_hybr_i
      hybr(flo,t1) = max(t1,abs(flo/2.),sqrt(t1**2+flo**2/4.))
      harm(t1,t2)= 2.0_R8*t1*t2/(t1+t2)
      pit(ix,iy)=abs(bb(ix,iy,0)/bb(ix,iy,3))
      upwind(flo,t1,t2) = max(flo,0.0_R8)*t1+min(flo,0.0_R8)*t2
#endif
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbc_phys')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2stbc_phm0', phm0)                                  !srv 17.01.07
       call ipgetr ('b2stbc_phm1', phm1)                                  !srv 17.01.07
       call ipgetr ('b2stbc_phm2', phm2)                                  !srv 17.01.07
        call ipgeti ('b2stbc_coreregno', coreregno)
        call xertst (cbirso.le.coreregno.and.
     1   coreregno.le.cbirso+cbnrso-1,
     2   'faulty internal parameter coreregno')
        if (nncut.ge.2) then
          call ipgeti ('b2stbc_coreregn2', coreregn2)
          call xertst (cbirso.le.coreregn2.and.
     1     coreregn2.le.cbirso+cbnrso-1,
     2     'faulty internal parameter coreregn2')
        endif
#ifdef WANT_THIS
        call ipgeti ('b2tfhe_conduction_only', conduction_only)
        call ipgeti ('b2tfhe_no_hybr', no_hybr_e)
        call ipgeti ('b2tfhi_no_hybr', no_hybr_i)
#endif
        call ipgeti ('b2mndr_eirene', use_eirene)
        call ipgeti ('b2mndr_astra', use_astra)
        call ipgeti ('b2tfnb_drift_style', drift_style)
        call ipgeti ('b2stbc_type13_ref', bc_type13_ref)
        call ipgeti ('b2stbc_type16_ref', bc_type16_ref)
        call ipgeti ('b2stbc_type20_ref', bc_type20_ref)
        call ipgeti ('b2stbc_type21_ref', bc_type21_ref)
        call ipgetr ('b2stbc_type13_norm', bc_type13_norm)
        call ipgetr ('b2stbc_type13_fac', bc_type13_fac)
        call ipgetr ('b2mndr_na0eps', na0eps)
        call xertst (0.lt.na0eps, 'faulty parameter na0eps')
        call ipgeti ('b2stbc_fix_fch_in_fhe_sheath',
     1   fix_fch_in_fhe_sheath)
        call ipgeti ('b2stbc_fix_sx_in_fhe_sheath',
     1   fix_sx_in_fhe_sheath)                                       !srv 08.05.09
        call ipgeti ('b2stbc_type16_kinetic_energy',
     1   bc_type16_kinetic_energy)
        call ipgetr ('b2stbc_cbc', cbc)                              !srv 14.10.99
        call ipgeti ('b2stbc_average_pr', istyle_average_pr)         !srv 08.08.02
        call ipgeti ('b2stbc_istyle_fchi', istyle_fchi)              !srv 14.02.05
        call ipgetr ('b2stbc_corr_flux_in_core', corr_flux)          !srv 25.07.05
        call ipgeti ('b2stbc_iform_corr_dia', iform_corr_dia)        !srv 25.07.05
        call ipgeti ('b2stbc_iout', iout)                            !srv 11.07.05
        call ipgeti ('b2tfnb_mdf', mdf_fnb)                          !srv 13.10.06
        call ipgeti ('b2tfhe_mdf', mdf_fhe)                          !srv 13.10.06
        call ipgeti ('b2tfhi_mdf', mdf_fhi)                          !srv 13.10.06
        call ipgeti ('b2stbc_istyle_cur_contr_on_S_and_N',
     &                       istyle_cur_contr_on_S_and_N)            !srv 01.02.09
        call ipgeti ('b2npc9_iout', iout_b2npc9)                     !srv 11.09.09 {
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)
        call ipgeti ('b2npp7_iout', iout_b2npp7)
        call ipgeti ('b2nppo_iout', iout_b2nppo)                     !srv 11.09.09 }
        call ipgetr ('b2stbc_bc_ref', bc_ref)
        call ipgetr ('b2stbc_bc_ref_te', bc_ref_te)                  !srv 11.04.13
        call ipgetr ('b2stbc_bc_ref_ti', bc_ref_ti)                  !srv 11.04.13
        call ipgeti ('b2stbc_bcene_15_style', bcene_15_style)
        call ipgeti ('b2stbc_bcene_16_style', bcene_16_style)
        call ipgeti ('b2stbc_bceni_16_style', bceni_16_style)
        call ipgeti ('b2tfhe_no_current', no_current)
        call ipgeti ('b2news_poteq', pot_eq)                         !srv 27.02.13
        call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &   'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
        call ipgetr ('b2tfhe_anomalous', fch_anomalous)
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
      call xertst (0.le.ismain0.and.ismain0.lt.ns,
     . 'invalid main neutral species index ismain0')
      call xertst (is_neutral(ismain0).or.ismain.eq.ismain0,
     . 'invalid main neutral species ismain0; must be neutral')
*   ..test facdrift (not on first call)                                  !xpb
      if (ncall.gt.0)                                                    !xpb
     > call xertst (0.0_R8.le.minval(facdrift) .and.
     &       1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol, hx, hy
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
*    ..test state
       call b2xvps (nx, ny, ns,
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
*    ..test edge values of fne, fni, chce, chve, chci, chvi, calf
       call b2xvff (nx, ny, fne, 'fne')
       call b2xvff (nx, ny, fni, 'fni')
       call b2xvff (nx, ny, chce, 'chce')
       call b2xvff (nx, ny, chve, 'chve')
       call b2xvff (nx, ny, chci, 'chci')
       call b2xvff (nx, ny, chvi, 'chvi')
       call b2xvff (nx, ny, calf, 'calf')
*    ..test sign of chce, chci
       call b2xvsg (n2*2, chce, 1, 'chce', '.ge.')
       call b2xvsg (n2*2, chci, 1, 'chci', '.ge.')
      endif
*   ..set cutll, cutlo
      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))

*  -- INITIALISE CERTAIN CONSTANTS --
      EPSLON=1.e-10_R8
*  -- SET NE00, UP00, UBIG, VBIG --
      NE00=1.0e20_R8
      UP00=SQRT(100*EV/MP)
      UBIG=UP00 / EPSLON
      VBIG=1.0e-2_R8*UBIG
      PBIG=UP00

      sna0_no_mdf = 0.0_R8

      call b2xppz (nx, ny, ns, ne, na, te, ti, pz)               !srv 18.01.02
*   ..find the average pressure for inner flux surface            !srv 24.03.99
      if (istyle_average_pr .eq.0) then                            !srv 08.08.02
       iy = -1                                                     !srv 24.03.99
       pzsum =0.0e0_R8                                             !srv 24.03.99 18.01.02
       volsum=0.0e0_R8                                             !srv 24.03.99
       do ix = -1,nx                                               !srv 26.11.02
        if (mod(region(ix,iy,0),4).eq.1) then                      !srv 26.11.02
         pzsum=pzsum + pz(ix,iy)*vol(ix,iy)                        !srv 24.03.99 18.01.02
         volsum=volsum + vol(ix,iy)                                !srv 24.03.99
        endif                                                      !srv 26.11.02
       enddo                                                       !srv 24.03.99
       pzaverage=pzsum/volsum                                      !srv 24.03.99 18.01.02
      else                                                         !srv 08.08.02
       iy = 0                                                      !srv 24.03.99  07.06.02 {
       pzsum =0.0e0_R8                                             !srv 24.03.99 18.01.02
       volsum=0.0e0_R8                                             !srv 24.03.99
       num_inner=0                                                 !srv 26.11.02
       do ix = -1,nx                                               !srv 26.11.02
        if (mod(region(ix,iy,0),4).eq.1) then                      !srv 26.11.02
         pzsum=pzsum
     *   + (pz(ix,iy)+pz(bottomix(ix,iy),bottomiy(ix,iy)))         !srv 24.03.99 18.01.02
         volsum=volsum + vol(ix,iy)                                !srv 24.03.99
         num_inner=num_inner+1                                     !srv 26.11.02
        endif                                                      !srv 26.11.02
       enddo                                                       !srv 24.03.99
       pzaverage=pzsum*0.5e0_R8/num_inner                         !srv 24.03.99 18.01.02 07.06.02 } 26.11.02
      endif                                                        !srv 08.08.02
*   ..find the average ne, te, ti, bb for inner flux surface      !srv 15.12.00 18.10.04
      bbsum = 0.0e0_R8                                            !srv 18.10.04
      nesum = 0.0e0_R8                                            !srv 15.12.00 25.01.02
      tesum = 0.0e0_R8                                            !srv 15.12.00
      tisum = 0.0e0_R8                                            !srv 15.12.00
      bzhzsum  = 0.0e0_R8                                         !lk 20.11.07
      bbaverage  = 0.0e0_R8                                       !lk 20.11.07
      do ix = -1,nx                                               !srv 26.11.02
       if (mod(region(ix,iy,0),4).eq.1) then                      !srv 26.11.02
        bbsum=bbsum + bb(ix,iy,3)*vol(ix,iy)                      !srv 18.10.04
        nesum=nesum + ne(ix,iy)*vol(ix,iy)                        !srv 15.12.00 25.01.02
        tisum=tisum + ti(ix,iy)*vol(ix,iy)                        !srv 15.12.00
        tesum=tesum + te(ix,iy)*vol(ix,iy)                        !srv 15.12.00
        bzhzsum = bzhzsum+abs(bb(ix,iy,2))*vol(ix,iy)**2/gs(ix,iy,2)   !lk 20.11.07 !srv 07.05.13
       endif                                                      !srv 26.11.02
      enddo                                                       !srv 15.12.00
      bbaverage=bbsum/volsum                                      !srv 18.10.04
      neaverage=nesum/volsum                                      !srv 15.12.00 25.01.02
      tiaverage=tisum/volsum                                      !srv 15.12.00
      teaverage=tesum/volsum                                      !srv 15.12.00
      bzhzaverage=bzhzsum/volsum                                  !lk 20.11.07
*   ..find the average ua**2, na for inner flux surface           !srv 22.01.02 {
      do is = 0, ns-1
       uasum  = 0.0e0_R8                                           !srv 15.02.02
       ua2sum = 0.0e0_R8                                           !srv 07.04.00 15.02.02
       nasum  = 0.0e0_R8
       do ix = -1,nx                                               !srv 26.11.02
        if (mod(region(ix,iy,0),4).eq.1) then                      !srv 26.11.02
         ua2sum = ua2sum + ua(ix,iy,is)**2*vol(ix,iy)              !srv 07.04.00 15.02.02
         uasum  = uasum  + ua(ix,iy,is)*vol(ix,iy)                 !srv 15.02.02
         nasum  = nasum  + na(ix,iy,is)*vol(ix,iy)
        endif                                                      !srv 26.11.02
       enddo
       ua2average(is) = ua2sum/volsum                              !srv 07.04.00 15.02.02
       uaaverage(is) = uasum/volsum                                !srv 07.04.00 15.02.02 16.04.02
       naaverage(is) = nasum/volsum
      enddo                                                        !srv 22.01.02 }

      call b2xpro (nx, ny, ns, mp, am, na, ro)                       !srv 01.02.09
      call b2xppr (nx, ny, ns, ne, na, te, ti, pr)
      wdia(:,:)=(ne(:,:)*te(:,:)+ni(:,:,1)*ti(:,:))*bb(:,:,2)
      do is = 0, ns_ext-1
        if (.not.is_ext_neutral(is)) wdia(:,:)=wdia(:,:)+
     .   na_ext(:,:,is)*ta_ext(:,:,is)*bb(:,:,2)
      end do
      do iy = -1, ny
        do ix = -1, nx
          if (rightix(ix,iy).ne.nx+1) then
            sx(ix,iy)=gs(rightix(ix,iy),rightiy(ix,iy),0)
            if (area_fix.ge.3)
     >       sx(ix,iy) = sx(ix,iy)*qc(rightix(ix,iy),rightiy(ix,iy))
            sxhz(ix,iy)=gs(rightix(ix,iy),rightiy(ix,iy),0)*         !srv 24.08.08
     &                  hz(ix,iy)
!     &                  hz(rightix(ix,iy),rightiy(ix,iy))
            if (area_fix.ge.3)
     >       sxhz(ix,iy) = sxhz(ix,iy)*qc(rightix(ix,iy),rightiy(ix,iy))
          else
            sx(ix,iy)=0.0_R8
            sxhz(ix,iy)=0.0_R8                                       !srv 24.08.08
          endif
          if (leftix(ix,iy).ne.-2) then
            do is = 0, ns-1
              uu(ix,iy,is)=fna(ix,iy,0,is)/gs(ix,iy,0)/
     .         ((na(ix,iy,is)*hx(leftix(ix,iy),leftiy(ix,iy))+
     .           na(leftix(ix,iy),leftiy(ix,iy),is)*hx(ix,iy))/
     .          (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy))))
            end do
            do is = 0, ns_ext-1
              uu(ix,iy,is+ns)=fa_ext(ix,iy,0,is)/gs(ix,iy,0)/
     .         ((na_ext(ix,iy,is)*hx(leftix(ix,iy),leftiy(ix,iy))+
     .           na_ext(leftix(ix,iy),leftiy(ix,iy),is)*hx(ix,iy))/
     .          (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy))))
            end do
            if (area_fix.ge.1) uu(ix,iy,:) = uu(ix,iy,:)/qc(ix,iy)
          else
            uu(ix,iy,:)=0.0_R8
          endif
        enddo
      enddo
cxpb  Following suggestion from Russians, the drift velocity
cxpb  should only include ExB and not diamagnetic contribution
      if (sheath_drift_fix.eq.0) then
        ud(:,:,:) = uadia(:,:,0,:)
      else
        ud(:,:,:) = vaecrb(:,:,0,:)
      endif
      if (main_call .and. diagno.ge.1) then                            !xpb
        write(*,'(a,1p,(10000e12.4))') 'dna0',
     1   (dna0(nx/2,ny,is),is=0,ns-1)
        write(*,'(a,1p,(10000e12.4))') 'dpa0',
     1   (dpa0(nx/2,ny,is),is=0,ns-1)
        write(*,'(a,1p,(10000e12.4))') 'vsa0',
     1   (vsa0(nx/2,ny,is),is=0,ns-1)
        write(*,'(a,1p,(10000e12.4))') 'hce0',
     1   hce0(nx/2,ny)
        write(*,'(a,1p,(10000e12.4))') 'hci0',
     1   hci0(nx/2,ny)
        write(*,'(a,1p,(10000e12.4))') 'sig0',
     1   sig0(nx/2,ny)
      endif                                                            !xpb

      if(ncall.eq.0) write(*,*) 'Doing B2-style bc'
      once_bccon = .false.                                             !srv 16.02.12 {
      once_bcene = .false.
      once_bceni = .false.                                             !srv 16.02.12 }

* -- LOOP OVER ALL BOUNDARIES
      DO IB=1, NBC
        if (bcpos(ib).eq.-2) then
          write(boundary_location(ib),'(a23,i2)')
     .     ' on boundary cell list ',IB
        else
          write(boundary_location(ib),'(a4,i4,a6,i4,a4,i4)')
     .     ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
        endif

* -- Moving core tally according to B.C.
* -- Must be done only once!

        IF (ncall.eq.0 .and. bc_type16_ref.gt.1) THEN
          IF (BCCHAR(IB).EQ.'S' .AND.
     &     (BCENE(IB).EQ.16.or.BCENE(IB).EQ.17.or.BCENE(IB).EQ.20.or.
     &      BCENI(IB).EQ.16.or.BCENI(IB).EQ.17.or.BCENI(IB).EQ.20)) THEN
            DO IX = 1, bc_list_size(ib)
              bc_type16_ix=bc_list_x(ix,ib)
              bc_type16_iy=bc_list_y(ix,ib)
              do bc_type16_i=1,bc_type16_ref      ! how many steps away from the boundary?
                bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
                bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
                bc_type16_ix=bc_type16_ixx
                call xertst (-2.lt.bc_type16_ix.and.
     &                       bc_type16_ix.lt.nx+1,
     1           'stepped off x-grid in type 16 bc')
                call xertst (-2.lt.bc_type16_iy.and.
     &                       bc_type16_iy.lt.ny+1,
     1           'stepped off y-grid in type 16 bc')
              enddo
              region(bc_type16_ix,bc_type16_iy,2) =
     .         region(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2)
              region(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2) = 0
            ENDDO
            if (bc_type16_ref.gt.2) then
              write(*,'(a,i2,a)')
     .         'Moved core radial flux tally position by ',
     .          bc_type16_ref-1,' radial steps.'
            else if (bc_type16_ref.eq.2) then
              write(*,'(a)')
     .         'Moved core radial flux tally position by '//
     .         '1 radial step.'
            endif
          ENDIF
        ENDIF

        DO IS=0, NS-1

* -- BOUNDARY CONDITIONS FOR CONTINUITY EQUATIONS --
          if(BCCON(IS,IB).eq.0) goto 10000
          GOTO(100,200,300,400,500,600,700,800,900,1000,
     1     1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,
     2     2100,2200,2300,2400),
     .     BCCON(IS,IB)
          write(*,*) BCCON(IS,IB)
          call xerrab ('BCCON out of bounds')

* -- BCCON=1 -- PRESCRIBE THE VALUE OF THE DENSITY

  100     CONTINUE
          if(ncall.eq.0) then
            IF(BCCHAR(IB).EQ.'S') THEN
            if(is.eq.ismain.and.use_astra.ne.0) then
              CONPAR(IS,IB,1) = neb_astra                               !srv 28.07.11
            endif
            ENDIF
            if (CONPAR(IS,IB,1).eq.0.0_R8) then
              CONPAR(IS,IB,1) = na0eps
              write(*,*) 'Zero density B.C. no longer allowed!'
              write(*,*) 'Changed to na0eps.'
            endif
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'constant density ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1       'BCCON = 1, CONPAR(,,1) !> 0')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*CONPAR(IS,IB,1)*
     .         SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .         -VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*CONPAR(IS,IB,1)*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=-VBIG*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            end do
          ELSEIF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         UBIG*CONPAR(IS,IB,1)*
     .         SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .         -UBIG*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .          UBIG*CONPAR(IS,IB,1)*
     .         SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .         -UBIG*
     .         SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            end do
          ELSEIF(BCCHAR(IB).EQ.'X') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)= UBIG*
     .         CONPAR(IS,IB,1)*VOL(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=-UBIG*
     .         VOL(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'Y') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)= UBIG*
     .         CONPAR(IS,IB,1)*VOL(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=-UBIG*
     .         VOL(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=1, INVALID OPTION.')
          ENDIF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                     !srv 23.09.08 }
          GOTO 10000

* -- BCCON=2 -- PRESCRIBE THE GRADIENT OF THE DENSITY

  200     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'density gradient ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .          VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         (-CONPAR(IS,IB,1)*0.5_R8*
     .         (HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     .          HY1(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))) +
     .         na(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .         -VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)= VBIG*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .         (CONPAR(IS,IB,1)*0.5_R8*
     .         (HY1(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     .          HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))) +
     .         na(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=-VBIG*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            end do
          ELSEIF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         UBIG*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     .         (-CONPAR(IS,IB,1)*0.5_R8*
     .         (HX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))+
     .          HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .             RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))) +
     .         na(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .         -UBIG*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)= UBIG*
     .         SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))*
     .         (CONPAR(IS,IB,1)*0.5_R8*
     .         (HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .             LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))+
     .          HX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))) +
     .         na(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=-UBIG*
     .         SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=2, INVALID OPTION.')
          END IF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                     !srv 23.09.08 }
          GOTO 10000

* -- BCCON=3 -- SHEATH CONDITIONS

  300     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,a,a,a,i3)')
     .       'sheath density bc on ',BCCHAR(IB),
     1       boundary_location(ib),
     2       ' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)= UBIG*
     .         na(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .         SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .         -UBIG*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)= UBIG*
     .         na(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .         SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=-UBIG*
     .         SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=3, NO SHEATH IN RADIAL DIR.')
          END IF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                     !srv 23.09.08 }
          GOTO 10000

* -- BCCON=4 -- PRESCRIBE THE VALUE OF THE DENSITY, WEAKLY
*                 A MIXED BOUNDARY CONDITION

  400     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'weak density ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1       'BCCON = 4, CONPAR(,,1) !> 0')
            call xertst(CONPAR(IS,IB,2).gt.0.0_R8,
     1       'BCCON = 4, CONPAR(,,2) !> 0')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .          CONPAR(IS,IB,2)
     .         *CONPAR(IS,IB,1)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     1         -CONPAR(IS,IB,2)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         CONPAR(IS,IB,2)*
     .         CONPAR(IS,IB,1)*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .         -CONPAR(IS,IB,2)*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            end do
          ELSEIF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .          CONPAR(IS,IB,2)
     .         *CONPAR(IS,IB,1)*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     1         -CONPAR(IS,IB,2)*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .          CONPAR(IS,IB,2)*
     .         CONPAR(IS,IB,1)*
     .         SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .         -CONPAR(IS,IB,2)*
     .         SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=4, NOT YET AVAILABLE')
          END IF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                     !srv 23.09.08 }
          GOTO 10000

* -- BCCON=5 -- PRESCRIBE THE PARTICLE FLUX PER UNIT AREA

  500     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'particle flux per area ',CONPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     1         CONPAR(IS,IB,1)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            end do
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         CONPAR(IS,IB,1)*
     1         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .           BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            end do
          ELSEIF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     1         CONPAR(IS,IB,1)*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
            end do
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         CONPAR(IS,IB,1)*
     1         SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .           LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=5, NOT YET AVAILABLE')
          END IF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                     !srv 23.09.08 }

          GOTO 10000

* -- BCCON=6 -- PRESCRIBE THE TOTAL PARTICLE FLUX FOR A
*                 CONSTANT DENSITY

cghgh
  600     CONTINUE
          if(ncall.eq.0) then
           write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .      'total particle flux ',CONPAR(IS,IB,1),
     1      ' with const density on ',BCCHAR(IB),
     2      boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'N') THEN
            s1=0.0_R8
            s2=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              s1=s1+SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                 BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/
     .         (0.5_R8*
     .         (HY1(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     .          HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
              s2=s2+SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                 BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .          na(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .             BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS)/
     .         (0.5_R8*
     .         (HY1(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     .          HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
            ENDDO
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)= VBIG*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .         (CONPAR(IS,IB,1)/dna0(bc_list_x(ix,ib),
     .                               bc_list_y(ix,ib),is) + s2) / s1
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=-VBIG*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            end do
          ELSEIF(BCCHAR(IB).EQ.'S') THEN
            s1=0.0_R8
            s2=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              s1=s1+SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .         (0.5_R8*(HY1(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                      TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     .                  HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
              s2=s2+SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         na(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS)/
     .         (0.5_R8*(HY1(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                      TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     .                  HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
            ENDDO
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .          VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         (CONPAR(IS,IB,1)/dna0(bc_list_x(ix,ib),
     .                               bc_list_y(ix,ib),is) + s2) / s1
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .          - VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            end do
          ELSE
            CALL XERRAB ('B2STBC -- BCCON=6, NOT AVAILABLE FOR W,E')
          ENDIF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                     !srv 23.09.08 }
          GOTO 10000
cghgh

* -- BCCON=7 -- PRESCRIBE THE DENSITY AS A FUNCTION OF OTHER
*                 PLASMA PARAMETERS

  700     CONTINUE
          CALL XERRAB ('B2STBC -- BCCON=7, NOT YET AVAILABLE')
          GOTO 10000

* -- BCCON=8 -- PRESCRIBE THE TOTAL PARTICLE FLUX WITH
*                 CONSTANT FLUX DENSITY

  800     CONTINUE
          if(ncall.eq.0) then
            if(is.eq.ismain.and.BCCHAR(IB).eq.'S'.and.use_astra.ne.0)
     >       CONPAR(IS,IB,1) = fneb_astra                               !srv 04.10.11
            write(*,'(a20,1p,g14.7,a4,a1,a,a13,i3)')
     .       'total particle flux ', CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            US=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              US=US+SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            ENDDO
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)+
     1         CONPAR(IS,IB,1)/US*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &            sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+
     &            sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &            sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &            sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) !srv 23.09.08
            if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) -
     &           fnaPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                   topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &            sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) +
     &            vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &            hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &            vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                   topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
            endif                                                    !srv 12.09.08 }
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            US=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              US=US+SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                 BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            ENDDO
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)+
     1         CONPAR(IS,IB,1)/US*
     2         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &          sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+
     &          sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &          sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &          sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) !srv 23.09.08
            if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) +
     &           fnaPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) -
     &           vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &           vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
            endif                                                    !srv 12.09.08 }
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'W') THEN
            US=0.0_R8
            DO IY = 1, BC_LIST_SIZE(IB)
              US=US+SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            ENDDO
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)+
     1         CONPAR(IS,IB,1)/US*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
              sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     &         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     &         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) !srv 23.09.08
            if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     &           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) -
     &           fnaPSch(
     &            rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &            rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     &           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) +
     &           vol(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &           hx(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     &           vaecrb(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                  rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)
            endif                                                    !srv 12.09.08 }
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            US=0.0_R8
            DO IY = 1, BC_LIST_SIZE(IB)
              US=US+SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            ENDDO
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)+
     1         CONPAR(IS,IB,1)/US*
     2         SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
              sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     .         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) +
     &           fnaPSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) -
     &           vol(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &           hx(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     &           vaecrb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)
            endif                                                    !srv 12.09.08 }
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'X') THEN
            US=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              US = US + SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            ENDDO
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)+
     1         CONPAR(IS,IB,1)/US*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &          sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+
     &          sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &          sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &          sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) +
     &          (fnaPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &           fnaPSch(
     &            topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &            topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is))/2.0_R8
                if (drift_style.eq.0) then
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) -
     &              vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &              hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &              vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
                else if (drift_style.eq.1) then
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) -
     &              vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &              hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*0.5_R8*
     &             (vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &              vaecrb(
     &               topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &               topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is))
                endif
              endif                                                    !srv 12.09.08 }
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'Y') THEN
            US=0.0_R8
            DO IY = 1, BC_LIST_SIZE(IB)
              US = US + SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            ENDDO
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)+
     1         CONPAR(IS,IB,1)/US*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
              sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     .         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) +
     &          (fnaPSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)+
     &           fnaPSch(
     &            rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &            rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is))/
     &           2.0_R8
                if (drift_style.eq.0) then
                  sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .             sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) -
     &             vol(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &             hx(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     &             vaecrb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)
                else if (drift_style.eq.1) then
                  sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .             sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) -
     &             vol(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &             hx(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*0.5_R8*
     &            (vaecrb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)+
     &             vaecrb(
     &              rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &              rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is))
                endif                                                   !srv 12.09.08 }
              endif
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCCON=8, NOT AVAILABLE')
          ENDIF
          GOTO 10000

* -- BCCON=9 -- PRESCRIBE THE DECAY LENGTH FOR THE DENSITY

  900     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'density decay length ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).ne.0.0_R8,
     1       'B2STBC -- BCCON=9, CONPAR(,,1) = 0')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=0.0_R8
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     1         -(dna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)+
     2           dpa0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)*
     3           (rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .             te(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .             ti(bc_list_x(ix,ib),bc_list_y(ix,ib))))/
     2         CONPAR(IS,IB,1)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     .         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 08.09.08 {
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     .           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) -
     &           fnaPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                   topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     .           vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .           vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
              endif                                                    !srv 08.09.08 }
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=0.0_R8
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     1         -(dna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)+
     2           dpa0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)*
     3           (rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .             te(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .             ti(bc_list_x(ix,ib),bc_list_y(ix,ib))))/
     2         CONPAR(IS,IB,1)*
     2         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     .         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 08.09.08 {
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) +
     &           fnaPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)-
     &           vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &           vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
              endif                                                    !srv 08.09.08 }
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCCON=9, NOT AVAILABLE FOR W,E')
          ENDIF
          GOTO 10000

* -- BCCON=10 -- leakage option for density

 1000     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'leakage ',CONPAR(IS,IB,1),' on ',BCCHAR(IB),
     1       boundary_location(ib),
     2       ' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              s0 = vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .             hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))                 !srv 26.02.09
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=0.0_R8
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=s0*
     2         sqrt((rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .                te(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .                ti(bc_list_x(ix,ib),bc_list_y(ix,ib)))/
     3         (am(is)*mp))*
     2         CONPAR(IS,IB,1)
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,is) =
     &         sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,is)+
     &         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,is) !srv 23.09.08
*    ..drifts to the wall                                                !srv 02.03.05 {
            if(mdf_fnb.eq.0) then                                    !srv 09.11.06 08.09.08 {
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     *           min(0.0e0_R8,s0*(
     .            vadia(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)+
     &           vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)))
            else
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     .           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) -
     &           fnaPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                   topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &            s0*min(0.0e0_R8,
     &           vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is))   !srv 02.03.05 }
            endif                                                    !srv 09.11.06 08.09.08 }
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              s0=vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))                  !srv 26.02.09
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=0.0_R8
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=s0*
     2         sqrt((rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .                te(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .                ti(bc_list_x(ix,ib),bc_list_y(ix,ib)))/
     3         (am(is)*mp))*
     2         CONPAR(IS,IB,1)
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,is) =
     &         sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,is)+
     &         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,is) !srv 23.09.08
*    ..drifts to the wall                                                !srv 02.03.05 {
            if(mdf_fnb.eq.0) then                                    !srv 09.11.06 08.09.08 {
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)-
     *           max(0.0e0_R8,
     .               s0*(vadia(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     *                  vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)))
            else
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     .           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) +
     &           fnaPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)-
     &           max(0.0e0_R8,
     .               s0*vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is))    !lk 12.12.07 }
            endif                                                    !srv 09.11.06 08.09.08 }
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              s0 = vol(bc_list_x(iy,ib),bc_list_y(iy,ib))/
     .              hx(bc_list_x(iy,ib),bc_list_y(iy,ib))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=0.0_R8
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=s0*
     2         sqrt((rza(bc_list_x(iy,ib),bc_list_y(iy,ib),is)*
     .                te(bc_list_x(iy,ib),bc_list_y(iy,ib))+
     .                ti(bc_list_x(iy,ib),bc_list_y(iy,ib)))/
     3         (am(is)*mp))*
     2         CONPAR(IS,IB,1)
              sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,is) =
     &         sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,is)+
     &         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,is) !srv 23.09.08
*    ..drifts to the wall                                                !srv 02.03.05 {
            if(mdf_fnb.eq.0) then                                    !srv 09.11.06 08.09.08 {
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is)-
     *           max(0.0e0_R8,s0*(vadia(
     .            rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)+
     *     vaecrb(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)))
            else
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) +
     &           fnaPSch(
     &            rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is)-
     &           max(0.0e0_R8,s0*vaecrb(
     .               rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is))    !lk 12.12.07 }
            endif                                                    !srv 09.11.06 08.09.08 }
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              s0 = vol(bc_list_x(iy,ib),bc_list_y(iy,ib))/
     .              hx(bc_list_x(iy,ib),bc_list_y(iy,ib))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=0.0_R8
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=s0*
     2         sqrt((rza(bc_list_x(iy,ib),bc_list_y(iy,ib),is)*
     .                te(bc_list_x(iy,ib),bc_list_y(iy,ib))+
     .                ti(bc_list_x(iy,ib),bc_list_y(iy,ib)))/
     3         (am(is)*mp))*
     2         CONPAR(IS,IB,1)
              sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,is) =
     &         sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,is)+
     &         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,is) !srv 23.09.08
*    ..drifts to the wall                                                !srv 02.03.05 {
            if(mdf_fnb.eq.0) then                                    !srv 09.11.06 08.09.08 {
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is)-
     *           max(0.0e0_R8,
     .               s0*(vadia(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)+
     *                  vaecrb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)))
            else
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) +
     &           fnaPSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is)-
     &           max(0.0e0_R8,
     .               s0*vaecrb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is))    !lk 12.12.07 }
            endif                                                    !srv 09.11.06 08.09.08 }
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCCON=10, NOT AVAILABLE')
          ENDIF
          GOTO 10000

* -- BCCON=11 -- particle flux feedback boundary condition

 1100     CONTINUE
          if(ncall.eq.0) then
           write(*,'(a,a,a,a,i3)')
     .      'particle flux feedback boundary condition on ',
     1      BCCHAR(IB),boundary_location(ib),' for species ',IS
           call ipgeti ('b2stbc_isfeedback', isfeedback)
           call xertst (0.le.isfeedback.and.isfeedback.lt.ns,
     1      'faulty internal parameter isfeedback')
           call xertst(is.eq.isfeedback,
     1       'BCCON(,IS,1) = 11 for IS <> ISFEEDBACK!')
           call xertst(lfeedback, 'LFEEDBACK must be turned on!')
          endif
          IF (BCCHAR(IB).EQ.'S') THEN
*   ..initialise ireg
            ireg = cbirso
*   ..update ireg
 1101       if (cbrbrk(ireg)*nx.lt.
     .          maxval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8) then
              call xertst (ireg.lt.cbirso+cbnrso-1,
     &         'error in computing region on south boundary')
              ireg = ireg+1
              goto 1101
            endif
            if (ireg.gt.cbirso) then
              call xertst (cbrbrk(ireg-1)*nx.lt.
     &         minval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8,
     &         'feedback bc covers more than one region!')
            endif
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         cbsna(0,IS,ireg)*
     .         SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            ENDDO
          ELSEIF (BCCHAR(IB).EQ.'N') THEN
*   ..initialise ireg
            ireg = cbirno
*   ..update ireg
 1102       if (cbrbrk(ireg)*nx.lt.
     .          maxval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8) then
              call xertst (ireg.lt.cbirno+cbnrno-1,
     &         'error in computing region on north boundary')
              ireg = ireg+1
              goto 1102
            endif
            if (ireg.gt.cbirno) then
              call xertst (cbrbrk(ireg-1)*nx.lt.
     &         minval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8,
     &         'feedback bc covers more than one region!')
            endif
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         cbsna(0,IS,ireg)*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            ENDDO
          ELSE
            CALL XERRAB
     .       ('B2STBC -- BCCON=11, NOT AVAILABLE FOR W,E')
          ENDIF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                     !srv 23.09.08 }
          GOTO 10000

* -- BCCON=12 -- particle density feedback boundary condition

 1200     CONTINUE
          if(ncall.eq.0) then
           write(*,'(a,a,a,a,i3)')
     .      'particle density feedback boundary condition on ',
     1      BCCHAR(IB),boundary_location(ib),' for species ',IS
           call ipgeti ('b2stbc_isfeedback', isfeedback)
           call xertst (0.le.isfeedback.and.isfeedback.lt.ns,
     1      'faulty internal parameter isfeedback')
           call xertst(is.eq.isfeedback,
     1       'BCCON(,IS,1) = 12 for IS <> ISFEEDBACK!')
           call xertst(lfeedback, 'LFEEDBACK must be turned on!')
           call ipgetr ('b2stbc_nesepm', nesepm)
           call ipgetr ('b2stbc_nesepm_overshoot', nesepm_overshoot)
           call xertst (0.0_R8.eq.nesepm.or.
     1      0.0_R8.eq.nesepm_overshoot,
     1      'BCCON=12 coding incompatible with nesepm_overshoot!')
          endif
          IF (BCCHAR(IB).EQ.'S') THEN
*   ..initialise ireg
            ireg = cbirso
*   ..update ireg
 1201       if (cbrbrk(ireg)*nx.lt.
     .          maxval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8) then
              call xertst (ireg.lt.cbirso+cbnrso-1,
     &         'error in computing region on south boundary')
              ireg = ireg+1
              goto 1201
            endif
            if (ireg.gt.cbirso) then
              call xertst (cbrbrk(ireg-1)*nx.lt.
     .          minval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8,
     &         'feedback bc covers more than one region!')
            endif
            call xertst (cbsna(0,is,ireg).gt.0.0_R8,
     &       'BCCON=12, trying to set a negative density!')
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*cbsna(0,IS,ireg)*
     .         SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            ENDDO
          ELSEIF (BCCHAR(IB).EQ.'N') THEN
*   ..initialise ireg
            ireg = cbirno
*   ..update ireg
 1202       if (cbrbrk(ireg)*nx.lt.
     .          maxval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8) then
              call xertst (ireg.lt.cbirno+cbnrno-1,
     &         'error in computing region on north boundary')
              ireg = ireg+1
              goto 1202
            endif
            if (ireg.gt.cbirno) then
              call xertst (cbrbrk(ireg-1)*nx.lt.
     .          minval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8,
     &         'feedback bc covers more than one region!')
            endif
            call xertst (cbsna(0,is,ireg).gt.0.0_R8,
     &       'BCCON=12, trying to set a negative density!')
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*cbsna(0,IS,ireg)*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=-VBIG*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            ENDDO
          ELSE
            CALL XERRAB
     .       ('B2STBC -- BCCON=12, NOT AVAILABLE FOR W,E')
          ENDIF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                     !srv 23.09.08 }
          GOTO 10000

* -- BCCON=13 -- Feedback constant density for a given particle flux

 1300     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,a,a,a,i3,a,1p,2g14.7)')
     .       'particle density to achieve specified total flux on ',
     1       BCCHAR(IB),boundary_location(ib),
     2       ' for species ',IS, ' with parameters ',
     3       CONPAR(IS,IB,1),CONPAR(IS,IB,2)
          endif
c dpc
          IF(BCCHAR(IB).EQ.'S') THEN
            FS=0.0_R8
            VS=0.0_R8
            NAS=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
               bc_type13_ix=bc_list_x(ix,ib)
               bc_type13_iy=bc_list_y(ix,ib)
               do bc_type13_i=1,bc_type13_ref ! how many steps away from the boundary?
                  bc_type13_ixx=topix(bc_type13_ix,bc_type13_iy)
                  bc_type13_iy=topiy(bc_type13_ix,bc_type13_iy)
                  bc_type13_ix=bc_type13_ixx
                  call xertst (-2.lt.bc_type13_ix.and.
     &                         bc_type13_ix.lt.nx+1,
     1                 'stepped off x-grid in type 13 bc')
                  call xertst (-2.lt.bc_type13_iy.and.
     &                         bc_type13_iy.lt.ny+1,
     1                 'stepped off y-grid in type 13 bc')
               enddo
               FS=FS+FNA(bc_type13_ix,bc_type13_iy,1,is)
               VS=VS+VOL(bc_list_x(ix,ib),bc_list_y(ix,ib))
               NAS=NAS+NA(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .          vol(bc_list_x(ix,ib),bc_list_y(ix,ib))
            ENDDO
            NAS=NAS/VS
            DFS=CONPAR(IS,IB,1)+CONPAR(IS,IB,3)
            FFAC=(1.0_R8+CONPAR(IS,IB,2)*
     1       (DFS-FS)/
     2       (abs(DFS)+abs(FS)+
     3        abs(bc_type13_norm)+abs(NAS*bc_type13_fac)))
            NAS=NAS*FFAC
            write(*,'(a,i3,1p,3(1x,e15.6))')
     1       'b2stbc_conpar=13: is,curr_flux,desired_flux,density ',
     2       IS,FS,DFS,NAS
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*NAS*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            ENDDO
          ELSE IF(BCCHAR(IB).EQ.'W') THEN
            FS=0.0_R8
            VS=0.0_R8
            NAS=0.0_R8
            DO IY = 1, BC_LIST_SIZE(IB)
               bc_type13_ix=bc_list_x(ib,iy)
               bc_type13_iy=bc_list_y(ib,iy)
               do bc_type13_i=1,bc_type13_ref ! how many steps away from the boundary?
                  bc_type13_ixx=rightix(bc_type13_ix,bc_type13_iy)
                  bc_type13_iy=rightiy(bc_type13_ix,bc_type13_iy)
                  bc_type13_ix=bc_type13_ixx
                  call xertst (-2.lt.bc_type13_ix.and.
     &                         bc_type13_ix.lt.nx+1,
     1                 'stepped off x-grid in type 13 bc')
                  call xertst (-2.lt.bc_type13_iy.and.
     &                         bc_type13_iy.lt.ny+1,
     1                 'stepped off y-grid in type 13 bc')
               enddo
               FS=FS+FNA(bc_type13_ix,bc_type13_iy,0,is)
               VS=VS+VOL(bc_list_x(ib,iy),bc_list_y(ib,iy))
               NAS=NAS+NA(bc_list_x(ib,iy),bc_list_y(ib,iy),is)*
     .          vol(bc_list_x(ib,iy),bc_list_y(ib,iy))
            ENDDO
            NAS=NAS/VS
            DFS=CONPAR(IS,IB,1)+CONPAR(IS,IB,3)
            FFAC=(1.0_R8+CONPAR(IS,IB,2)*
     1       (DFS-FS)/
     2       (abs(DFS)+abs(FS)+
     3        abs(bc_type13_norm)+abs(NAS*bc_type13_fac)))
            NAS=NAS*FFAC
            write(*,'(a,i3,1p,3(1x,e15.6))')
     1       'b2stbc_conpar=13: is,curr_flux,desired_flux,density ',
     2       IS,FS,DFS,NAS
            DO IY = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IB,IY),BC_LIST_Y(IB,IY),0,IS)=
     .         VBIG*NAS*SX(BC_LIST_X(IB,IY),BC_LIST_Y(IB,IY))
              SNA0(BC_LIST_X(IB,IY),BC_LIST_Y(IB,IY),1,IS)=
     .        -VBIG*SX(BC_LIST_X(IB,IY),BC_LIST_Y(IB,IY))
            ENDDO
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=13, ONLY CODED FOR SOUTH AND WEST')
          ENDIF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                     !srv 23.09.08 }
          GOTO 10000

* -- BCCON=14 -- Sound speed flux

 1400     CONTINUE
          if (ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       ' sound speed flux ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'W') THEN
              DO IY = 1, BC_LIST_SIZE(IB)
                S1 = ABS(
     1           PBS(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0))
                CS=SQRT(PR(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     .                  RO(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
!                CS=SQRT(BCINT(PR,NX,NY,IY,IB)/BCINT(RO,NX,NY,IY,IB))
                SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=0.0_R8
                SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .           SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)+
     .           S1*CONPAR(IS,IB,1)*CS
                sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,IS)=
     &           sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,IS) +
     &           SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,IS)
              if(mdf_fnb.ne.0) then                                   !srv 13.10.06 {
                  sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     &             sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)-
     &          fnaPSch(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                  rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)
              endif                                                   !srv 13.10.06 }
              ENDDO
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              DO IY = 1, BC_LIST_SIZE(IB)
                S1 = ABS(PBS(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))
                CS=SQRT(PR(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     .                  RO(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
!                CS=SQRT(BCINT(PR,NX,NY,IY,IB)/BCINT(RO,NX,NY,IY,IB))
                SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=0.0_R8
                SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .           SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)+
     .           S1*CONPAR(IS,IB,1)*CS
                sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,IS)=
     &           sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,IS) +
     .           SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,IS)
              if(mdf_fnb.ne.0) then                                   !srv 13.10.06 {
                  sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     &             sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)+
     &             fnaPSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)
              endif                                                   !srv 13.10.06 }
              ENDDO
            ELSE
              CALL XERRAB ('B2STBC -- BCCON=14, NOT IN RADIAL DIR')
          ENDIF
          GOTO 10000

* -- BCCON=15 -- prescribe the decay length for density for drifts switched on

 1500     CONTINUE                                                  !srv 25.09.09 16.02.12{
          if(ncall.eq.0) then
            write(*,'(2a,1p,g14.7,a,a,a,a,i3)')
     .       'prescribe the decay length for the density',
     &       ' for drifs switched on ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN                                !srv 23.07.10 {
            DO IX = 1, BC_LIST_SIZE(IB)
              csb=sqrt((rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &                   te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &                   ti(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/
     &                 (am(is)*mp))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)= VBIG*
     &         CDNA(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &              TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,IS)*
     .           NA(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=-VBIG*
     &         (CDNA(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &               TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,IS)-
     &          CONPAR(IS,IB,1)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &          csb)
            END DO
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              csb=sqrt((rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &                   te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &                   ti(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/
     &                 (am(is)*mp))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)= VBIG*
     &         CDNA(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*
     .         NA(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)= -VBIG*
     &         (CDNA(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)-
     &          CONPAR(IS,IB,1)*
     .          SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &             BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*csb)
            END DO
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=15, NO IN POLOIDAL DIR.')
          END IF                                                    !srv 23.07.10 }
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                     !srv 23.09.08 }
          GOTO 10000                                                !srv 25.09.09 }

* -- BCCON=16 -- Feedback boundary condition: given total particle flux with constant density  !srv 16.02.12 {

 1600     CONTINUE
          if(ncall.eq.0) then
            call xertst(is.eq.ismain,
     1       'BCCON = 16 is applied for main ion only ')
            write(*,'(a,a,i3)') 
     .       'given total particle flux with constant',
     &       ' density for main ion is =', ismain
!           write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
            if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     1           boundary_location(ibw),' for species ',IS
                call xertst(
     .         mod(region(BC_LIST_X(1,IBW),BC_LIST_Y(1,IBW),0),4).eq.1,
     .           'BCCON = 16 only applicable to core boundaries !')
                call xertst(CONPAR(IS,IBW,1).gt.0.0_R8,
     1           'BCCON = 16, CONPAR(,,1) !> 0')
              endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.19.or.
     1            BCCON(IS,IBW).EQ.21.or.BCCON(IS,IBW).EQ.22.or.
     1            BCCON(IS,IBW).EQ.23)) then
                call xerrab
     .           ('BCCON types 16, 19, 21, 22 and 23 cannot be mixed!')
            endif
          ENDDO
            if (use_astra.ne.0) then
#ifdef ASTRA
              SCONPAR = fneb_astra
#else
              call xerrab
     .         ('Compile with -DASTRA option to couple to ASTRA!')
#endif
            else
              SCONPAR = 0.0_R8
              DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                SCONPAR = SCONPAR + CONPAR(IS,IBW,1)
              endif
            ENDDO
            endif
            write(*,*) ' Summed value is ', SCONPAR
          endif
          NAS = 0.0_R8
          VS = 0.0_R8
          DO IBW=IB,NBC
          if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN
               DO IX = 1, BC_LIST_SIZE(IBW)
               NAS = NAS+NA(BC_LIST_X(IBW,IX),BC_LIST_Y(IBW,IX),is)*
     .                     VOL(BC_LIST_X(IBW,IX),BC_LIST_Y(IBW,IX))
                 VS = VS + VOL(BC_LIST_X(IBW,IX),BC_LIST_Y(IBW,IX))
               END DO
             endif
          END DO
          NAS=NAS/VS
          if(.not.once_bccon) then
            IF(BCCHAR(IB).EQ.'S') THEN
            totflux = 0.0_R8
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                if(mdf_fnb.ne.0) then
                    DO IX = 1, BC_LIST_SIZE(IBW)
                    totflux = totflux +
     &         fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,is)
                    ENDDO
                else
                    DO IX = 1, BC_LIST_SIZE(IBW)
                    totflux = totflux +
     &            fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,is)
                    ENDDO
                endif
              ENDIF
            ENDDO
              NAS = NAS *
     &         (1.0_R8-bc_ref*(totflux-SCONPAR)/
     &                 max(abs(totflux),SCONPAR))
!              T0 = (totflux-CONPAR(IS,IB,1))/
!     &           max(abs(totflux),CONPAR(IS,IB,1))
!              NAS = NAS *
!     &         (1.0_R8-max(-bc_ref,min(bc_ref,T0)))
              write(*,"(' BCCON=16: NAS=',es14.7,' totflux=',es14.7,
     *         ' SCONPAR=',es14.7,' IS=',i3)") NAS, totflux, SCONPAR, IS
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     *                VBIG*NAS*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     *               -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  END DO
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)=
     .               SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  END DO                                              !srv 23.09.08 }
              END IF
            END DO
            once_bccon = .true.
            ELSE                                                       !srv 16.02.12 }
              CALL XERRAB
     1         ('B2STBC -- BCCON=16, NOT AVAILABLE FOR N,W,E')
            ENDIF
        endif                                                          !srv 16.02.12 }
          GOTO 10000

* -- BCCON=17 -- FILLER

 1700     CONTINUE
          CALL XERRAB ('B2STBC -- BCCON=17, NOT YET AVAILABLE')
          GOTO 10000

#ifdef ASTRA
* -- BCCON=18 -- PRESCRIBE THE MODIFIED TOTAL PARTICLE FLUX

 1800     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,a,a,a,i3)') 'modified total particle flux on ',
     1         BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(is.eq.ismain,
     1       'BCCON = 18 is applied for main ion only ')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
          t0 = 0.0_R8
            US=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              US=US+SY(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX))
              DO ISB = 0, NS-1
                if (is_neutral(isb)) then
                  t0 = t0 +
     .             fna(topix(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX)),
     .                 topiy(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX)),1,isb)
                endif
              END DO
            END DO
            CONPAR(IS,IB,1) = (fSeNBI_astra - t0)                     !srv 29.11.11

!            if(ncall.eq.0) then
               write (*,'(a,a,i3,1p,3(a,g11.4))')
     &        ' BCCON=18', ' IS =', IS,
     &        ' CONPAR=',CONPAR(IS,IB,1), ' fSeNBI_astra=',fSeNBI_astra,
     &        ' integrated neutral flux =',t0
!            endif

            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX),0,IS)=
     &         SNA0(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX),0,IS)+         !srv 01.02.09
     1          CONPAR(IS,IB,1)/US*SY(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX))
              SNA0(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX),1,IS)=0.0_R8
              sna0_no_mdf(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX),0,is) =
     &         sna0_no_mdf(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX),0,is)+
     &         sna0(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX),0,is) !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX),1,is) =
     &         sna0_no_mdf(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX),1,is)+
     &         sna0(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX),1,is) !srv 23.09.08
            if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                sna0(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX),0,is) =
     &           sna0(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX),0,is) -
     &           fnaPSch(topix(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX)),
     &                   topiy(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX)),1,is)
               sna0(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX),1,is) =
     &          sna0(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX),1,is) +
     &           vol(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX))/
     &           hy1(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX))*
     &           vaecrb(topix(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX)),
     &                  topiy(BC_LIST_X(IB,IX),BC_LIST_Y(IB,IX)),1,is)
            endif                                                    !srv 12.09.08 }
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCCON=18, ONLY FOR SOUTH')
          ENDIF
          GOTO 10000
#else
* -- BCCON=18 -- FILLER

 1800     CONTINUE
          CALL XERRAB ('B2STBC -- BCCON=18, NOT YET AVAILABLE')
          GOTO 10000
#endif

* -- BCCON=19 -- Feedback boundary condition: given sum of integrated neutrals and
*                main ion particle flux with constant ion density

 1900     CONTINUE
          if(ncall.eq.0) then
            call xertst(is.eq.ismain,
     1       'BCCON = 19 is applied for main ion only ')
          write(*,*) 'given sum of integrated neutrals and main ',
     &                 'ion particle flux with constant ion density '
!           write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
            if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     1           boundary_location(ibw),' for species ',IS
                call xertst(
     .         mod(region(BC_LIST_X(1,IBW),BC_LIST_Y(1,IBW),0),4).eq.1,
     .            'BCCON = 19 only applicable to core boundaries !')
                call xertst(CONPAR(IS,IBW,1).gt.0.0_R8,
     1           'BCCON = 19, CONPAR(,,1) !> 0')
            endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.21.or.BCCON(IS,IBW).EQ.22.or.
     1            BCCON(IS,IBW).EQ.23)) then
                call xerrab
     .           ('BCCON types 16, 19, 21, 22 and 23 cannot be mixed!')
              endif
          ENDDO
            if (use_astra.ne.0) then
#ifdef ASTRA
              SCONPAR = fSeNBI_astra
#else
              call xerrab
     .         ('Compile with -DASTRA option to couple to ASTRA!')
#endif
            else
              SCONPAR = 0.0_R8
              DO IBW=IB,NBC
                if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  SCONPAR = SCONPAR + CONPAR(IS,IBW,1)
              endif
            ENDDO
            endif
          write(*,*) ' Summed value is ', SCONPAR
            NAS = 0.0_R8
            VS = 0.0_R8
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
              NAS = NAS+NA(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)*
     .                   VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                VS = VS + VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
              endif
            ENDDO
            NAS = NAS/VS
          endif
        if(.not.once_bccon) then
            IF(BCCHAR(IB).EQ.'S') THEN
            totfluxn = 0.0_R8
            totfluxi = 0.0_R8
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
#ifdef B25_EIRENE
                    isa = b2espcr(ismain0)
                    if (isa.ne.-1) then
                      totfluxn = totfluxn + rfluxa(
     &             topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+1,
     &             topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+1,isa,1)*
     &                    sy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                      do ism = 1,nnmoli
                        if (mlcmp(isa,ism).ne.0) then
                    totfluxn = totfluxn + mlcmp(isa,ism)*rfluxm(
     &             topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+1,
     &             topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+1,ism,1)*
     &                sy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                        endif
                      enddo
                    endif
#endif
                   isb = ismain0
                   if (is_neutral(isb))
     >              totfluxn = totfluxn + fna(
     &             topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &             topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,isb)
                  isb = ismain
                  if(mdf_fnb.ne.0) then
                    totfluxi = totfluxi + fna_mdf(
     &                 topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,isb)
                  else
                    totfluxi = totfluxi + fna(
     &                 topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,isb)
                  endif
                ENDDO
              ENDIF
            ENDDO
            NAS = NAS *
     &       (1.0_R8-bc_ref*(totfluxn+totfluxi-SCONPAR)/
     &           max(abs(totfluxn+totfluxi),SCONPAR))
            write(*,"(' BCCON=19: NAS=',es14.7,' totfluxn=',es14.7,
     &         ' totfluxi=',es14.7,' SCONPAR=',es14.7,
     &         ' totfluxn+totfluxi=',es14.7)")
     *         NAS, totfluxn, totfluxi, SCONPAR, totfluxn+totfluxi
            DO IBW=IB,NBC
              IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &             VBIG*NAS*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &            -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                ENDDO
                DO IX = 1, BC_LIST_SIZE(IBW)
                  SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                        BC_LIST_Y(IX,IBW),:,IS)=
     &             SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                ENDDO                                                !srv 23.09.08 }
              ENDIF
            ENDDO
            once_bccon = .true.
            ELSE                                                       !srv 16.02.12 }
              CALL XERRAB
     1         ('B2STBC -- BCCON=19, NOT AVAILABLE FOR N,W,E')
            ENDIF
        endif                                                          !srv 16.02.12 }
          GOTO 10000

* -- BCCON=20 -- constant density feedback condition

 2000     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'constant density (neo) ',
     1       CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1       'BCCON = 1, CONPAR(,,1) !> 0')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            naave = 0.0_R8
            volave = 0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type20_ix=bc_list_x(ix,ib)
              bc_type20_iy=bc_list_y(ix,ib)
              do bc_type20_i=1,bc_type20_ref      ! how many steps away from the boundary?
                bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
                bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
                bc_type20_ix=bc_type20_ixx
                call xertst (-2.lt.bc_type20_ix.and.
     &                       bc_type20_ix.lt.nx+1,
     1           'stepped off x-grid in type 20 bc')
                call xertst (-2.lt.bc_type20_iy.and.
     &                       bc_type20_iy.lt.ny+1,
     1           'stepped off y-grid in type 20 bc')
              enddo
              naave=naave+na(bc_type20_ix,bc_type20_iy,is)*
     1         vol(bc_type20_ix,bc_type20_iy)
              volave=volave+vol(bc_type20_ix,bc_type20_iy)
            end do
            naave=naave/volave
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type20_ix=bc_list_x(ix,ib)
              bc_type20_iy=bc_list_y(ix,ib)
              do bc_type20_i=1,bc_type20_ref      ! how many steps away from the boundary?
                bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
                bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
                bc_type20_ix=bc_type20_ixx
                call xertst (-2.lt.bc_type20_ix.and.
     &                       bc_type20_ix.lt.nx+1,
     1           'stepped off x-grid in type 20 bc')
                call xertst (-2.lt.bc_type20_iy.and.
     &                       bc_type20_iy.lt.ny+1,
     1           'stepped off y-grid in type 20 bc')
              enddo
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*CONPAR(IS,IB,1)*
     .         SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         na(bc_type20_ix,bc_type20_iy,is)/naave
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=20, NOT AVAILABLE FOR N,W,E')
          ENDIF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                     !srv 23.09.08 }
          GOTO 10000

!! FIXME: Deal with BC number clash here in a better way
cxpb Careful here because of BC number clash
* From 5.0: starting at version 01.001.030
* -- BCCON=21 -- constant density scaled by density on the ring bc_type21_ref away
* to be moved up to next available number

* From 5.2, after version number 01.001.078:
* -- BCCON=21 -- PRESCRIBE THE VALUE OF THE DENSITY  !lk and IYS 01.03.2013{
* AND ADD A DENSITY PERTURBATION TO GET A SOLUTION
* WHICH IS AS CLOSE AS POSSIBLE TO NEOCLASSICAL THEORY

* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH
* CORRESPONDING CONDITION ON IONS TEMPERATURE (BCENI=24)

 2100     CONTINUE
          if(ncall.eq.0) then
            if (boundary_version.ge.'01.001.030' .and.
     &          boundary_version.lt.'01.001.078') then
              write(*,*) 'BCCON=21 type B.C. moved to 24 ',
     &         ' because of clash with 5.2 numbering. '
              write(*,*) 'Please modify b2.boundary.parameters file'
              call xerrab('B2STBC -- BCCON=21, DISPLACED')
            endif
            write(*,'(a,1p,g14.7)')
     .       'constant density',CONPAR(IS,IB,1)
            call xertst(is.eq.ismain,
     1       'BCCON = 21 is applied for main ion only ')
            DO IBW=IB,NBC
            if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                call xertst(
     .         mod(region(BC_LIST_X(1,IBW),BC_LIST_Y(1,IBW),0),4).eq.1,
     .           'BCCON = 21 only applicable to core boundaries !')
                call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1           'BCCON = 21, CONPAR(,,1) !> 0')
              endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.16.or.
     1           BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.22.or.
     1           BCCON(IS,IBW).EQ.23)) then
                call xerrab
     .           ('BCCON types 16, 19, 21, 22 and 23 cannot be mixed!')
              endif
              if (bceni(ibw).ne.23.and.bceni(ibw).ne.24) then
                write(0,*)
     &    'Should use BCCON(IS,IB).EQ.21 together with BCENI=23 or 24'
                write(*,*)
     &    'Should use BCCON(IS,IB).EQ.21 together with BCENI=23 or 24'
                write(*,*)
     &           'Please ensure that this boundary condition '//
     &           'is consistent with '//
     &           'ion energy boundary condition'  !srv 21.09.12
              endif
            END DO
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
          if(is.eq.ismain.and.use_astra.ne.0) then
#ifdef ASTRA
            CONPAR(IS,IB,1) = neb_astra                               !srv 28.07.11
#else
              call xerrab
     .         ('Compile with -DASTRA option to couple to ASTRA!')
#endif
            else
              SCONPAR = 0.0_R8
              DO IBW=IB,NBC
                if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  SCONPAR = SCONPAR + CONPAR(IS,IBW,1)
              endif
            ENDDO
            endif
          write(*,*) ' Summed value is ', SCONPAR
          ENDIF
          if (.not.once_bccon) then
            IF(BCCHAR(IB).EQ.'S' .and.
     $       (mod(region(BC_LIST_X(1,IB),BC_LIST_Y(1,IB),0),4).eq.1))
     >       THEN
!            DO IX = 1, BC_LIST_SIZE(IB)
!              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
!     .         VBIG*CONPAR(IS,IB,1)*
!     .         SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
!              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
!     .        -VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
!            end do
            t0=0.0e0_R8   ! Volume between two surfaces                  !lk 12.12.12 {
            t1=0.0e0_R8   ! Sum of volume-averaged electrons and ions temperatures
              do ix = -1,nx
                if (mod(region(ix,0,0),4).eq.1) then
                t0=t0 + vol(ix,0)
                t1=t1 + (ti(ix,-1)+te(ix,-1))*vol(ix,0)
                endif
            enddo
            t1=t1/t0        ! Sum of volume-averaged electrons and ions temperatures
              DO IBW=IB,NBC
                do IX = 1, BC_LIST_SIZE(IBW)
                  if(mod(region(BC_LIST_X(IX,IBW),
     &                          BC_LIST_Y(IX,IBW),0),4).eq.1) then
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*CONPAR(IS,IBW,1)*t1
     &                 /(ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &                 te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &               SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                endif
              enddo                                    !lk 12.12.12 }
                DO IX = 1, BC_LIST_SIZE(IBW)                             !srv 23.09.08
                  SNA0_NO_MDF(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)=
     &             SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                ENDDO
              ENDDO
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCCON=21, VALID FOR CORE ONLY')
            ENDIF
            once_bccon = .true.
          endif
          GOTO 10000 !lk and IYS 01.03.2013}

* -- BCCON=22 -- Feedback boundary condition: given total particle flux with constant
*    average density. The poloidal density variation is added to make the solution
*    as close as possible to the neoclassical value. It is recommended to use this
*    boundary condition with corresponding boundary condition on ions temperature
*    (BCENI=24) !srv 16.02.12, IYS 01.03.2013 {

 2200     CONTINUE
          if(ncall.eq.0) then
            call xertst(is.eq.ismain,
     1       'BCCON = 22 is applied for main ion only ')
            write(*,'(a,a,i3)') 
     &       'given total particle flux with constant',
     &       ' density for main ion is =', ismain
!           write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
            if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                call xertst(
     .         mod(region(BC_LIST_X(1,IBW),BC_LIST_Y(1,IBW),0),4).eq.1,
     .           'BCCON = 22 only applicable to core boundaries !')
                call xertst(CONPAR(IS,IBW,1).gt.0.0_R8,
     1          'BCCON = 22, CONPAR(,,1) !> 0')
            endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.21.or.
     1            BCCON(IS,IBW).EQ.23)) then
                call xerrab
     .           ('BCCON types 16, 19, 21, 22 and 23 cannot be mixed!')
              endif
              if (bceni(ibw).ne.23.and.bceni(ibw).ne.24) then
                write(0,*)
     &    'Should use BCCON(IS,IB).EQ.22 together with BCENI=23 or 24'
                write(*,*)
     &    'Should use BCCON(IS,IB).EQ.22 together with BCENI=23 or 24'
                write(*,*)
     &           'Please ensure that this boundary condition '//
     &           'is consistent with '//
     &           'ion energy boundary condition'  !srv 21.09.12
              endif
           ENDDO
            if (use_astra.ne.0) then
#ifdef ASTRA
              SCONPAR = fneb_astra
#else
              call xerrab
     .         ('Compile with -DASTRA option to couple to ASTRA!')
#endif
            else
              SCONPAR = 0.0_R8
              DO IBW=IB,NBC
                if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  SCONPAR = SCONPAR + CONPAR(IS,IBW,1)
              endif
            ENDDO
            endif
            write(*,*) ' Summed value is ', SCONPAR
          endif
          if(.not.once_bccon) then
            IF(BCCHAR(IB).EQ.'S') THEN
            totflux = 0.0_R8
            t0=0.0e0_R8   ! Volume between two surfaces                 !lk 12.12.12 {
            t1=0.0e0_R8   ! Sum of volume-averaged electrons and ions temperatures
              NAS=0.0e0_R8
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                if(mdf_fnb.ne.0) then
                    DO IX = 1, BC_LIST_SIZE(IBW)
                    totflux = totflux +
     &         fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,is)
                    t0=t0 + vol(BC_LIST_X(IX,IBW),0)
                    t1=t1 + (ti(BC_LIST_X(IX,IBW),-1)+
     &                         te(BC_LIST_X(IX,IBW),-1))*
     .                        vol(BC_LIST_X(IX,IBW),0)
                      NAS=NAS+na(BC_LIST_X(IX,IBW),-1,is)*
     .                       vol(BC_LIST_X(IX,IBW),0)
                    END DO
                else
                    DO IX = 1, BC_LIST_SIZE(IBW)
                    totflux = totflux +
     &         fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &             topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,is)
                    t0=t0 + vol(BC_LIST_X(IX,IBW),0)
                    t1=t1 + (ti(BC_LIST_X(IX,IBW),-1)+
     &                         te(BC_LIST_X(IX,IBW),-1))*
     &                        vol(BC_LIST_X(IX,IBW),0)
                      NAS=NAS+na(BC_LIST_X(IX,IBW),-1,is)*
     &                       vol(BC_LIST_X(IX,IBW),0)
                    ENDDO
                endif
              ENDIF
            ENDDO
            t1=t1/t0        ! Sum of volume-averaged electrons ans ions temperatures
              NAS=NAS/t0      ! Volume-averaged density

              NAS = NAS*
     &         (1.0_R8-bc_ref*(totflux-SCONPAR)/
     &           max(abs(totflux),SCONPAR))
!              T0 = (totflux-CONPAR(IS,IB,1))/
!     &           max(abs(totflux),CONPAR(IS,IB,1))
!              NAS = NA((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB),is) *
!     &         (1.0_R8-max(-bc_ref,min(bc_ref,T0)))
              write(*,"(' BCCON=22: NAS=',es14.7,' totflux=',es14.7,
     *       '  SCONPAR=',es14.7,' IS=',i3)") NAS, totflux, SCONPAR, IS
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*NAS*t1
     &                   /(ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &                 te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
     &               *SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))  ! lk
                   SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                 ENDDO
                 DO IX = 1, BC_LIST_SIZE(IBW)                           !srv 23.09.08 {
                   SNA0_NO_MDF(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
     &             =SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                 ENDDO                                                  !srv 23.09.08 }
             ENDIF
           ENDDO
             once_bccon = .true.
           ELSE                                                       !srv 16.02.12 }
             CALL XERRAB
     1        ('B2STBC -- BCCON=22 IS AVAILABLE FOR CORE ONLY')
           ENDIF
         endif                                                          !srv 16.02.12 }
         GOTO 10000                            !IYS 01.03.2013 }

* -- BCCON=23 -- Feedback boundary condition: given total particle flux with constant
*    average density. The poloidal density variation is added to make the solution
*    as close as possible to the neoclassical value. It is recommended to use this boundary
*    condition with corresponding boundary condition on ions temperature (BCENI=24) !srv 16.02.12, IYS 04.03.2013 {

 2300    CONTINUE
         if(ncall.eq.0) then
           call xertst(is.eq.ismain,
     1      'BCCON = 23 is applied for main ion only ')
           write(*,'(a,a,i3)') 
     .      'given total particle flux with constant',
     &      ' density for main ion is =', ismain
!          write all boundaries with this type of boundary condition
           DO IBW=IB,NBC
             if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >        THEN
               write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .          boundary_location(ibw),' for species ',IS
               call xertst(
     .          mod(region(BC_LIST_X(1,IBW),BC_LIST_Y(1,IBW),0),4).eq.1,
     .          'BCCON = 23 only applicable to core boundaries !')
               call xertst(CONPAR(IS,IBW,1).gt.0.0_R8,
     1          'BCCON = 23, CONPAR(,,1) !> 0')
             endif
             if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.16.or.
     1          BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.21.or.
     1          BCCON(IS,IBW).EQ.22)) then
               call xerrab
     .          ('BCCON types 16, 19, 21, 22 and 23 cannot be mixed!')
             endif
             if (bceni(ibw).ne.23..and.bceni(ibw).ne.24) then
               write(0,*)
     &     'Should use BCCON(IS,IB).EQ.23 together with BCENI=23 or 24'
               write(*,*)
     &     'Should use BCCON(IS,IB).EQ.23 together with BCENI=23 or 24'
               write(*,*)
     &          'Please ensure that this boundary condition '//
     &          'is consistent with '//
     &          'ion energy boundary condition'  !srv 21.09.12
             endif
           ENDDO
         endif

         if(.not.once_bccon) then
           if (use_astra.ne.0) then
#ifdef ASTRA
             SCONPAR = fSeNBI_astra
#else
             call xerrab
     .        ('Compile with -DASTRA option to couple to ASTRA!')
#endif
           else
             SCONPAR = 0.0_R8
             DO IBW=IB,NBC
               if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >          THEN
                 SCONPAR = SCONPAR + CONPAR(IS,IBW,1)
               endif
             ENDDO
           endif
           write(*,*) ' Summed value is ', SCONPAR
           IF(BCCHAR(IB).EQ.'S') THEN
             totflux = 0.0_R8
             t0=0.0e0_R8   ! Volume between two surfaces                  !lk 12.12.12 {
             t1=0.0e0_R8   ! Sum of volume-averaged electrons and ions temperatures
             NAS=0.0e0_R8
             DO IBW=IB,NBC
               IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >          THEN
                 DO IX = 1, BC_LIST_SIZE(IBW)
#ifdef B25_EIRENE
                  isa = b2espcr(ismain0)
                  if (isa.ne.-1) then
                    totfluxn = totfluxn + rfluxa(
     &             topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+1,
     &             topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+1,isa,1)
     &               *sy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    do ism = 1,nnmoli
                      if (mlcmp(isa,ism).ne.0) then
                  totfluxn = totfluxn + mlcmp(isa,ism)*rfluxm(
     &             topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+1,
     &             topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+1,ism,1)*
     &                sy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                      endif
                    enddo
                  endif
#endif
                  isb = ismain0
                  if (is_neutral(isb))
     >             totfluxn = totfluxn +
     &           fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,isb)
                  isb = ismain
                  if(mdf_fnb.ne.0) then
                    totfluxi = totfluxi +
     &       fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,isb)
                  else
                    totfluxi = totfluxi +
     &           fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,isb)
                  endif
                  t0=t0+vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                      topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t1=t1 + (ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     .                     te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     .              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  NAS=NAS+na(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)*
     .              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                ENDDO
              ENDIF
            ENDDO
            t1=t1/t0        ! Sum of volume-averaged electrons and ions temperatures
            NAS=NAS/t0      ! Volume-averaged density
            NAS = NAS*
     &         (1.0_R8-bc_ref*(totfluxn+totfluxi-SCONPAR)/
     &           (abs(totfluxn)+abs(totfluxi)+abs(SCONPAR)))
            write(*,"(' BCCON=23: NAS=',es14.7,' totfluxn=',es14.7,
     &        ' totfluxi=',es14.7,' SCONPAR=',es14.7,
     &        ' totfluxn+totfluxi=',es14.7)")
     *      NAS, totfluxn, totfluxi, SCONPAR, totfluxn+totfluxi
            DO IBW=IB,NBC
              IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*NAS*t1
     &               /(ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &               te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &               SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))  ! lk
                  SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &             -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                ENDDO
                DO IX =1, BC_LIST_SIZE(IBW)                          !srv 23.09.08 {
                  SNA0_NO_MDF(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)=
     &             SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                ENDDO                                                 !srv 23.09.08 }
            ENDIF
          ENDDO
            once_bccon = .true.
          ELSE                                                       !srv 16.02.12 }
            CALL XERRAB
     1       ('B2STBC -- BCCON=23 IS AVAILABLE FOR CORE ONLY')
          ENDIF
      endif                                                          !srv 16.02.12 }
        GOTO 10000                            !IYS 04.03.2013 }

* -- BCCON=24 -- constant density scaled by density on the ring bc_type21_ref away

 2400     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'constant density (modified by nearby ring) ',
     1       CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1       'BCCON = 24, CONPAR(,,1) !> 0')
            call xertst(CONPAR(is,ib,2).gt.0.0_R8,
     &           'faulty parameter conpar2')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            naave = 0.0_R8
            volave = 0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type21_ix=bc_list_x(ix,ib)
              bc_type21_iy=bc_list_y(ix,ib)
              do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
                bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
                bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
                bc_type21_ix=bc_type21_ixx
                call xertst (-2.lt.bc_type21_ix.and.
     &                       bc_type21_ix.lt.nx+1,
     1           'stepped off x-grid in type 24 bc')
                call xertst (-2.lt.bc_type21_iy.and.
     &                       bc_type21_iy.lt.ny+1,
     1           'stepped off y-grid in type 24 bc')
              enddo
              naave=naave+na(bc_type21_ix,bc_type21_iy,is)*
     1         vol(bc_type21_ix,bc_type21_iy)
              volave=volave+vol(bc_type21_ix,bc_type21_iy)
            end do
            naave=naave/volave
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type21_ix=bc_list_x(ix,ib)
              bc_type21_iy=bc_list_y(ix,ib)
              do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
                bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
                bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
                bc_type21_ix=bc_type21_ixx
                call xertst (-2.lt.bc_type21_ix.and.
     &                       bc_type21_ix.lt.nx+1,
     1           'stepped off x-grid in type 24 bc')
                call xertst (-2.lt.bc_type21_iy.and.
     &                       bc_type21_iy.lt.ny+1,
     1           'stepped off y-grid in type 24 bc')
              enddo
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*CONPAR(IS,IB,1)*
     .         SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         ((1.0_R8-conpar(is,ib,2))+
     .           conpar(is,ib,2)*na(bc_type21_ix,bc_type21_iy,is)/naave)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            end do
          ELSE IF(BCCHAR(IB).EQ.'N') THEN
            naave = 0.0_R8
            volave = 0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type21_ix=bc_list_x(ix,ib)
              bc_type21_iy=bc_list_y(ix,ib)
              do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
                bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
                bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
                bc_type21_ix=bc_type21_ixx
                call xertst (-2.lt.bc_type21_ix.and.
     &                       bc_type21_ix.lt.nx+1,
     1           'stepped off x-grid in type 24 bc')
                call xertst (-2.lt.bc_type21_iy.and.
     &                       bc_type21_iy.lt.ny+1,
     1           'stepped off y-grid in type 24 bc')
              enddo
              naave=naave+na(bc_type21_ix,bc_type21_iy,is)*
     1         vol(bc_type21_ix,bc_type21_iy)
              volave=volave+vol(bc_type21_ix,bc_type21_iy)
            end do
            naave=naave/volave
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type21_ix=bc_list_x(ix,ib)
              bc_type21_iy=bc_list_y(ix,ib)
              do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
                bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
                bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
                bc_type21_ix=bc_type21_ixx
                call xertst (-2.lt.bc_type21_ix.and.
     &                       bc_type21_ix.lt.nx+1,
     1           'stepped off x-grid in type 24 bc')
                call xertst (-2.lt.bc_type21_iy.and.
     &                       bc_type21_iy.lt.ny+1,
     1           'stepped off y-grid in type 24 bc')
              enddo
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*CONPAR(IS,IB,1)*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .         ((1.0_R8-conpar(is,ib,2))+
     .           conpar(is,ib,2)*na(bc_type21_ix,bc_type21_iy,is)/naave)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=-VBIG*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=24, NOT AVAILABLE FOR W,E')
          ENDIF
          GOTO 10000

* -- BOUNDARY CONDITIONS FOR PARALLEL MOMENTUM EQUATIONS --

10000     CONTINUE
          if(BCMOM(IS,IB).eq.0) goto 20000
        GOTO(10100,10200,10300,10400,10500,
     1         10600,10700,10800,10900,11000,
     1     11100,11200,11300,11400,11500,11600),BCMOM(IS,IB)      !srv 01.02.09 19.04.13
          write(*,*) BCMOM(IS,IB)
          call xerrab ('BCMOM out of bounds')

* -- BCMOM=1 -- PRESCRIBE THE VALUE OF THE PARALLEL VELOCITY

10100     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'specified parallel velocity ',MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*NE00*AM(IS)*
     .         MOMPAR(IS,IB,1)*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -VBIG*NE00*AM(IS)*
     .         SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*NE00*AM(IS)*
     .         MOMPAR(IS,IB,1)*
     .         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -VBIG*NE00*AM(IS)*
     .         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            end do
          ELSEIF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         UBIG*NE00*AM(IS)*
     .         MOMPAR(IS,IB,1)*SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .        -UBIG*NE00*AM(IS)*
     .         SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         UBIG*NE00*AM(IS)*
     .         MOMPAR(IS,IB,1)*
     .         SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .        -UBIG*NE00*AM(IS)*
     .         SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            end do
          ELSEIF(BCCHAR(IB).EQ.'X') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         UBIG*NE00*AM(IS)*
     .         MOMPAR(IS,IB,1)*VOL(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -UBIG*NE00*AM(IS)*
     .         VOL(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'Y') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         UBIG*NE00*AM(IS)*
     .         MOMPAR(IS,IB,1)*VOL(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .        -UBIG*NE00*AM(IS)*
     .         VOL(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            end do
          ELSE
            CALL XERRAB ('B2STBC -- BCMOM=1, INVALID OPTION.')
          ENDIF
          GOTO 20000

* -- BCMOM=2 -- PRESCRIBE THE GRADIENT OF THE PARALLEL VELOCITY

10200     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'parallel velocity gradient ',MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)= VBIG*
     1         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     2         AM(IS)*na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)*
     .         (MOMPAR(IS,IB,1)*
     1         (HY1(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     1          HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/2 +
     2         UA(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=-VBIG*
     .         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))) *
     .         AM(IS)*na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            end do
          ELSEIF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*AM(IS) *
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS) *
     .         (-MOMPAR(IS,IB,1)*
     1         (HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     2          HY1(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))/2 +
     .         UA(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        - VBIG*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)) *
     .         AM(IS) *
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            end do
          ELSEIF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         VBIG*SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*AM(IS) *
     .         na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),IS) *
     .         (-MOMPAR(IS,IB,1)*
     1         (HX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))+
     2          HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .             RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))/2 +
     .         UA(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .         - VBIG*SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)) *
     .         AM(IS) *
     .         na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),IS)
            end do
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)= VBIG*
     .         SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))*AM(IS) *
     .         na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),IS) *
     .         (MOMPAR(IS,IB,1)*
     1         (HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .             LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))+
     2          HX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/2 +
     .         UA(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=-VBIG*
     .         SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))) *
     .         AM(IS) *
     .         na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),IS)
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCMOM=2, INVALID OPTION.')
          END IF
          GOTO 20000

* -- BCMOM=3 -- SHEATH CONDITIONS, MACH NUMBER AS INPUT

10300     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'sheath velocity ',MOMPAR(IS,IB,2),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              IF(MOMPAR(IS,IB,2).LT.0.5_R8) THEN
                CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))  !srv 01.02.09
                UUOUT=MOMPAR(IS,IB,1)*CS*
     .           abs(wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .               wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3))
                if(.not.is_neutral(is)) then                           !xpb
                  UPOUT=max(0.0_R8, UUOUT/
     .             abs(wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .                 wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)) +    !xpb
     .             bcint(ud(-1,-1,is),nx,ny,iy,ib)*
     .             abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/
     .                 wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)))     !xpb
                else                                                   !xpb
                  UPOUT = UUOUT/
     .             abs(wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .                 wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3))      !xpb
                endif                                                  !xpb
cxpb  Added multiplicative factor to UD according to the Chodura criterion
cxpb  and max(,,) statement for stability
              ELSE
                ZBAR=bcint(na(-1,-1,is),nx,ny,iy,ib)/
     .           bcint(NE,nx,ny,IY,IB)
                ZARG=(GAMMAI*bcint(ti,nx,ny,iy,ib)
     .           +bcint(rz2(-1,-1,is),nx,ny,iy,ib)*ZBAR
     .           *bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                CS=SQRT(ZARG)
                UUTEST=
     .           ua(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .           (1.0_R8+
     .            HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .    HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .       RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .    ua(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .       RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .            HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .    HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .       RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
                UUOUT=MAX(MOMPAR(IS,IB,1)*CS,-UUTEST)
                if(.not.is_neutral(is)) then                           !xpb
                  UPOUT=max(0.0_R8, UUOUT+                             !xpb
     .             bcint(ud(-1,-1,is),nx,ny,iy,ib)*
     .             abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/
     .                 wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)))     !xpb
                else                                                   !xpb
                  UPOUT = UUOUT                                        !xpb
                endif                                                  !xpb
cxpb  Added multiplicative factor to UD according to the Chodura criterion
cxpb  and max(,,) statement for stability
              ENDIF
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         (UBIG/EPSLON**2)*AM(IS)*mp*
     1         bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .         SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     2         (-UPOUT)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .         -(UBIG/EPSLON**2)*AM(IS)*mp*
     .         bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .         SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              IF(MOMPAR(IS,IB,2).LT.0.5_R8) THEN
                CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB)) !srv 01.02.09
                UUOUT=MOMPAR(IS,IB,1)*CS*
     .           abs(wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .               wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3))
                if(.not.is_neutral(is)) then                           !xpb
                  UPOUT=max(0.0_R8, UUOUT/                             !xpb
     .              abs(wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .                  wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3))     !xpb
     .             -bcint(ud(-1,-1,is),nx,ny,iy,ib)*
     .              abs(wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/
     .                  wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)))    !xpb
                else                                                   !xpb
                  UPOUT = UUOUT/
     .             abs(wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .                 wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3))      !xpb
                endif                                                  !xpb
cxpb  Added multiplicative factor to UD according to the Chodura criterion
cxpb  and max(,,) statement for stability
              ELSE
                ZBAR=bcint(na(-1,-1,is),nx,ny,iy,ib)/
     .           bcint(NE,nx,ny,IY,IB)
                ZARG=(GAMMAI*bcint(ti,nx,ny,iy,ib)
     .           +bcint(rz2(-1,-1,is),nx,ny,iy,ib)*ZBAR
     .           *bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                CS=SQRT(ZARG)
                UUTEST=
     .           ua(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .           (1.0_R8+
     .            HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .       HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .       ua(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .            HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .       HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
                UUOUT=MAX(MOMPAR(IS,IB,1)*CS,UUTEST)
                if(.not.is_neutral(is)) then                           !xpb
                  UPOUT=max(0.0_R8, UUOUT -                            !xpb
     .             bcint(ud(-1,-1,is),nx,ny,iy,ib)*
     .             abs(wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/
     .                 wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)))     !xpb
                else                                                   !xpb
                  UPOUT = UUOUT                                        !xpb
                endif                                                  !xpb
cxpb  Added multiplicative factor to UD according to the Chodura criterion
cxpb  and max(,,) statement for stability
              ENDIF
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         (UBIG/EPSLON**2)*AM(IS)*mp*
     1         bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     2         SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))*UPOUT
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .        -(UBIG/EPSLON**2)*AM(IS)*mp*
     .         bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .         SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCMOM=3, NO SHEATH IN RADIAL DIR.')
          END IF
          GOTO 20000

* -- BCMOM=4 -- PRESCRIBE THE VALUE OF THE VELOCITY, WEAKLY
*                 A MIXED BOUNDARY CONDITION

10400     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'weak velocity ',MOMPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(MOMPAR(IS,IB,2).gt.0.0_R8,
     1       'BCMOM = 4, MOMPAR(,,2) !> 0')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         MOMPAR(IS,IB,2)
     .         *MOMPAR(IS,IB,1)*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -MOMPAR(IS,IB,2)
     .         *SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         MOMPAR(IS,IB,2)*
     .         MOMPAR(IS,IB,1)*
     .         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -MOMPAR(IS,IB,2)*
     1         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            end do
          ELSEIF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         MOMPAR(IS,IB,2)
     .         *MOMPAR(IS,IB,1)*SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     1         -MOMPAR(IS,IB,2)*SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         MOMPAR(IS,IB,2)*
     .         MOMPAR(IS,IB,1)*
     .         SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .        -MOMPAR(IS,IB,2)*
     .         SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCMOM=4, OPTION NOT AVAILABLE')
          ENDIF
          GOTO 20000

* -- BCMOM=5 -- PRESCRIBE THE PARALLEL MOMENTUM FLUX PER UNIT AREA

10500     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'parallel momentum flux per unit area ',MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         MOMPAR(IS,IB,1)*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            end do
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         MOMPAR(IS,IB,1)*
     1         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            end do
          ELSEIF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     1         MOMPAR(IS,IB,1)*SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
            end do
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         MOMPAR(IS,IB,1)*
     1         SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCMOM=5, OPTION NOT AVAILABLE')
          ENDIF
          GOTO 20000

* -- BCMOM=6 -- PRESCRIBE THE TOTAL PARALLEL MOMENTUM FLUX FOR A
*                 CONSTANT PARALLEL VELOCITY

10600     CONTINUE
          CALL XERRAB ('B2STBC -- BCMOM=6, NOT YET AVAILABLE')
          GOTO 20000

* -- BCMOM=7 -- PRESCRIBE THE PARALLEL MOMENTUM AS A FUNCTION OF OTHER
*                 PLASMA PARAMETERS

10700     CONTINUE
          CALL XERRAB ('B2STBC -- BCMOM=7, NOT YET AVAILABLE')
          GOTO 20000

* -- BCMOM=8 -- SPECIAL : LIMITED SHEAR

10800     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,a,a,a,i3)')
     .       'limited shear velocity on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
#ifdef WANT_THIS
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              A=MIN(1.0_R8,
     .         SQRT(TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .              TI(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                 TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*AM(IS) *
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)*A*
     .         UA(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        - VBIG*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*AM(IS) *
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            end do
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              A=MIN(1.0_R8,
     .         SQRT(TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .          TI(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .             BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)= VBIG*
     .         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .         AM(IS) * na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)*A*
     .         UA(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=-VBIG*
     .         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .         AM(IS) * na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCMOM=8, NOT AVAILABLE FOR W,E')
          END IF
#else
          CALL XERRAB ('B2STBC -- BCMOM=8, NOT IMPLEMENTED.')
#endif
          GOTO 20000

* -- BCMOM=9 -- PRESCRIBE THE TOTAL PARALLEL MOMENTUM FLUX WITH
*                 CONSTANT FLUX DENSITY

10900     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'total parallel momentum flux ',MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            US=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              US=US+SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            ENDDO
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS) =
     1         MOMPAR(IS,IB,1)/US*
     1         SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            US=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              US=US+SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                   BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            ENDDO
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS) =
     1         MOMPAR(IS,IB,1)/US*
     2         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'W') THEN
            US=0.0_R8
            DO IY = 1, BC_LIST_SIZE(IB)
              US=US+SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            ENDDO
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         MOMPAR(IS,IB,1)/US*
     1         SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            US=0.0_R8
            DO IY = 1, BC_LIST_SIZE(IB)
              US=US+SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            ENDDO
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         MOMPAR(IS,IB,1)/US*
     1         SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'X') THEN
            US=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              US=US+SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            ENDDO
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         MOMPAR(IS,IB,1)/US*
     1         SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'Y') THEN
            US=0.0_R8
            DO IY = 1, BC_LIST_SIZE(IB)
              US=US+SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            ENDDO
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         MOMPAR(IS,IB,1)/US*
     1         SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCMOM=9, NOT AVAILABLE')
          ENDIF
          GOTO 20000

* -- BCMOM=10 -- PRESCRIBE THE DECAY LENGTH FOR THE PARALLEL MOMENTUM

11000     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'parallel momentum decay length ',MOMPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(MOMPAR(IS,IB,1).ne.0.0_R8,
     1       'BCMOM = 10, MOMPAR(,,1) = 0')
            if (bccon(is,ib).ne.9) then
              write(0,*)
     &         'Should use BCCON(IS,IB).EQ.9 together with BCMOM=10'
              write(*,*)
     &         'Should use BCCON(IS,IB).EQ.9 together with BCMOM=10'
              write(*,*) 'Please ensure that this boundary condition '//
     &                   'is consistent with density boundary condition'  !srv 21.09.12
            endif
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)=
     .         AM(IS)*mp*ua(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
     .         *SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -vsa0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)/
     1         MOMPAR(IS,IB,1)*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     .         AM(IS)*mp*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         AM(IS)*mp*ua(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
     .         *SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -vsa0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)/
     1         MOMPAR(IS,IB,1)*
     2         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))) +
     .         AM(IS)*mp*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCMOM=10, NOT AVAILABLE FOR W,E')
          ENDIF
          GOTO 20000

* -- BCMOM=11 -- ROZHANSKY VISCOSITY CONDITION FOR THE PARALLEL MOMENTUM

11100     CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'ROZ parallel momentum bc on ',BCCHAR(IB),
     1       boundary_location(ib),
     2       ' for species ',IS
            write(*,*) ' BCMOM=11 no longer supported !'
            write(0,*) ' BCMOM=11 no longer supported !'
          endif
#ifdef WANT_THIS
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              IF (BC_LIST_X(IX,IB).GE.0) THEN
                smo0(bc_list_x(ix,ib),bc_list_y(ix,ib),0,is)=
     1           vsa0(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .           na(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*am(is)*
     2           sqrt(1-pit(bc_list_x(ix,ib),bc_list_y(ix,ib))**2)*
     .           ua(bc_list_x(ix,ib),bc_list_y(ix,ib),is)
                smo0(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is)=
     1           -dna0(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*am(is)*
     2           (na(bc_list_x(ix,ib),bc_list_y(ix,ib),is)-
     2            na(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),is))/
     3           hy1(bc_list_x(ix,ib),bc_list_y(ix,ib))-
     3           vsa0(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .           na(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     4           am(is)*
     .           sqrt(1-pit(bc_list_x(ix,ib),bc_list_y(ix,ib))**2)/
     .           hy1(bc_list_x(ix,ib),bc_list_y(ix,ib))
              end if
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCMOM=11, NOT AVAILABLE FOR N,W,E')
          ENDIF
#else
          CALL XERRAB ('B2STBC -- BCMOM=11, NOT IMPLEMENTED.')
#endif
          GOTO 20000

* -- BCMOM=12 -- CONDITION from b2stbc_spb FOR THE PARALLEL MOMENTUM   !srv 01.02.09 {

11200     CONTINUE
          if(ncall.eq.0) then !srv added accumulation in order to account twice b.c. for corner cells
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a,a,i3)')
     .       'from b2stbc_spb parallel momentum bc with mompar1 =',
     1       MOMPAR(IS,IB,1),' and mompar2 =', MOMPAR(IS,IB,2),' on ',
     2       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              s0hz = hz(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &              vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &              hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))             !srv 09.01.01
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),2,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3,IS)=
     &         SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3,IS)+
     &             MOMPAR(IS,IB,2)*s0hz +
!     &         MOMPAR(IS,IB,2)*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &         MOMPAR(IS,IB,1)*
     &         max(0.0e0_R8,-fna_fcor(
     &          TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &          TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is))/
     &          na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
            END DO
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              s0hz = hz(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &              vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &              hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))             !srv 09.01.01
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),2,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3,IS)=
     &         SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3,IS)+
     &         MOMPAR(IS,IB,2)*s0hz+
!     &           MOMPAR(IS,IB,2)*
!     &           SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
!     &                BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     &         MOMPAR(IS,IB,1)*
     &         max(0.0e0_R8,
     &          fna_fcor(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is))/
     &                na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCMOM=12, NOT AVAILABLE FOR W,E')
          ENDIF
          GOTO 20000

* -- BCMOM=13 -- CONDITION from b2stbc_spb FOR THE PARALLEL MOMENTUM

11300     CONTINUE
          if(ncall.eq.0) then   !srv added accumulation in order to account twice b.c. for corner cells
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a,a,i3)')
     .       'from b2stbc_spb parallel momentum bc with mompar1 =',
     1       MOMPAR(IS,IB,1),' and mompar2 =', MOMPAR(IS,IB,2),' on ',
     2       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
!              CS=SQRT(BCINT(PR,NX,NY,IY,IB)/BCINT(RO,NX,NY,IY,IB))
              cs=sqrt(pr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &                ro(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              s1hz = pbshz(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                     rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0) !srv 11.01.13
              vwest=-cs*s1hz -                          !srv 17.02.04 19.02.08 11.01.13
     &         hz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     &         vol(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &             rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     &          hx(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &             rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))*  !*signmf*bb(bc_list_x(iy,ib),bc_listy(iy,ib),2)/bb(bc_list_x(iy,ib),bc_list_y(iy,ib),3)   !srv 05.02.05
     &      vaecrb(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &             rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)*
     &           cbc*          !srv 17.02.04 19.02.08
     &          (abs(bb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))/
     &               bb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))        !srv 11.01.13
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS) =
     &         SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS) +  !srv 07.05.09 {
     &         MOMPAR(IS,IB,1)*(UBIG/EPSLON**2)*AM(IS)*mp*
     &         na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*vwest
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     &         SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)-
     &         (UBIG/EPSLON**2)*AM(IS)*mp*
     &          na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*abs(s1hz) !srv 07.05.09 }
!             SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),2,IS)=
!     &        SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),2,IS)+      !srv 07.05.09 { equivalence form
!     &        MOMPAR(IS,IB,1)*(UBIG/EPSLON**2)*vwest
!             SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3,IS)=
!      &       SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3,IS)-     !srv 07.05.09 }
!     &        (UBIG/EPSLON**2)*S1HZ
            END DO
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
!             CS=SQRT(BCINT(PR,NX,NY,IY,IB)/BCINT(RO,NX,NY,IY,IB))
              cs=sqrt(pr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &                ro(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              s1hz = pbshz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)      !srv 11.01.13
             veast = cs*s1hz-                                      !srv 17.02.04 19.02.08 25.04.08 { 11.01.13
     &          hz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     &          vol(leftix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),iy)/
     &          hx(leftix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),iy)        !*signmf*bb(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/bb(bc_list_x(iy,ib),bc_list_y(iy,ib),3) !srv 05.02.05
     &          *vaecrb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)*cbc*    !srv 17.02.04 19.02.08
     &          (abs(bb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))/
     &               bb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))         !srv 11.01.13
               SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS) =
     &          SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS) +   !srv 07.05.09 {
     &          MOMPAR(IS,IB,1)*(UBIG/EPSLON**2)*AM(IS)*mp*
     &          na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*veast
               SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     &          SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)-
     &          (UBIG/EPSLON**2)*AM(IS)*mp*
     &          na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*abs(s1hz)  !srv 07.05.09 }
!              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),2,IS)=
!     &         SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),2,IS)+      !srv 07.05.09 { equivalence form
!     &         MOMPAR(IS,IB,1)*(UBIG/EPSLON**2)*veast
!              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3,IS)=
!               SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3,IS) -    !srv 07.05.09 }
!     &         (UBIG/EPSLON**2)*S1HZ
            END DO
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCMOM=13, NOT AVAILABLE FOR S,N')
          ENDIF
          GOTO 20000                                                 !srv 01.02.09 }

* -- BCMOM=14 -- CONDITION from b2stbc_spb FOR THE PARALLEL MOMENTUM  !srv 01.02.09 {

11400     CONTINUE
          if(ncall.eq.0) then !srv added accumulation in order to account twice b.c. for corner cells
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a,a,i3)')
     .       'from b2stbc_spb momentum bc with mompar1 =',
     1       MOMPAR(IS,IB,1),' and mompar2 =', MOMPAR(IS,IB,2),' on ',
     2       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              s0hz = hz(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &              vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &              hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))              !srv 09.01.01
              if (.not.is_neutral(is)) then                          !srv 15.02.02 {
              select case (iform_corr_dia)                          !srv 25.07.05
              case (0)
                  smo0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &             smo0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+  !srv 09.01.01 01.07.05
     &             s0hz*(am(is)*mp)*naaverage(is)*                   !srv 15.02.02
     &            (-2.0e0_R8*
     &             facdrift(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*phm2*  !srv 07.06.00 17.05.99 25.04.00 15.02.02
     &             tiaverage*bb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),2)/
     &             (qe*rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &             hx(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*     !srv 15.02.02
     &            (OnedBsqc(rightix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                      rightiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))-
     &              OnedBsqc(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))        !srv 05.03.99 15.02.02
     &             )*uaaverage(is) + corr_flux*s0hz                  !srv 15.02.02 08.09.03 25.07.05
                case (1)
                 smo0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &            smo0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+  !srv 09.01.01 01.07.05
     &            s0hz*(am(is)*mp)*naaverage(is)*                   !srv 15.02.02
     &           (-2.0e0_R8*facdrift(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &            phm2* !srv 07.06.00 17.05.99 25.04.00 15.02.02
     &            tiaverage*bb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),2)/
     &            (qe*rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &            hx(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*               !srv 15.02.02
     &            (OnedBsqc(rightix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                      rightiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))-
     &             OnedBsqc(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))) !srv 05.03.99 15.02.02
     &            )*ua(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)+corr_flux
     &             *s0hz               !srv 08.09.03 25.07.05
          end select                                           !srv 25.07.05
        endif                                                 !srv 15.02.02
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),2,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3,IS)=
     &         SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3,IS)+
     &         MOMPAR(IS,IB,2)*s0hz +
!     &         MOMPAR(IS,IB,2)*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &         MOMPAR(IS,IB,1)*
     &         max(0.0e0_R8,
     &       -fna_fcor(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is))/
     &                    na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
            END DO
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCMOM=14, NOT AVAILABLE FOR N,W,E')
          ENDIF
          GOTO 20000

* -- BCMOM=15 -- PRESCRIBE THE VALUE OF THE PARALLEL VELOCITY AS A FUNCTION OF !srv 02.04.09 {
*                POLOIDAL COORDINATE

11500        CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,a,1p,g14.7,a,a,a,a,i3)')
     .       'specified parallel velocity as a function of',
     &       ' poloidal coordinate ',MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*NE00*AM(IS)*MOMPAR(IS,IB,1)*
     .         (bbaverage/bb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3))*
     &         SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))                  !srv 24.08.08
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     &        -VBIG*NE00*AM(IS)*
     .         SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))                  !srv 24.08.08
            END DO
          ELSE
            CALL XERRAB ('B2STBC -- BCMOM=15, ONLY FOR S')
          ENDIF
          GOTO 20000                                                !srv 02.04.09 }

* -- BCMOM=16 -- PRESCRIBE THE AVERAGE VALUE OF THE PARALLEL VELOCITY !lk 12.12.12 {

11600        CONTINUE
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'specified average parallel velocity ',MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*NE00*AM(IS)*
     .         MOMPAR(IS,IB,1)*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))  !srv 24.08.08
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -VBIG*NE00*AM(IS)*
     .         SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))                  !srv 24.08.08
            END DO
            t0=0.0e0_R8
            t1=0.0e0_R8
            do ix = -1,nx
              if (mod(region(ix,0,0),4).eq.1) then
              t0=t0 + vol(ix,0)
              t1=t1 + UA(IX,0,IS)*vol(ix,0)
              endif
          enddo
            t1=t1/t0
            do ix = 1, BC_LIST_SIZE(IB)
            if(mod(region(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0),4).eq.1)
     >       then
                SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .           SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)+
     .             VBIG*NE00*AM(IS)*
     .         (UA(IX,TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS)-t1)*
     .           SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            endif
          END DO                        !lk 12.12.12 }
          ELSE
            CALL XERRAB ('B2STBC -- BCMOM=16, IS VALID FOR CORE ONLY')
          ENDIF
          GOTO 20000

* -- END LOOP OVER SPECIES
20000     continue
        enddo

* -- BOUNDARY CONDITIONS FOR ELECTRON ENERGY EQUATION --

        if(BCENE(IB).eq.0) goto 30000
        GOTO(20100,20200,20300,20400,20500,
     1   20600,20700,20800,20900,21000,
     2   21100,21200,21300,21400,21500,
     3   21600,21700,21800,21900,22000,
     4   22100,22200), BCENE(IB)
        write(*,*) BCENE(IB)
        call xerrab ('BCENE out of bounds')

* -- BCENE=1 -- PRESCRIBE THE VALUE OF THE ELECTRON TEMPERATURE

20100   CONTINUE
        if(ncall.eq.0) then
          IF(BCCHAR(IB).EQ.'S'.and.use_astra.ne.0) THEN
          ENEPAR(IB,1) = teb_astra*1.0e3_R8                        !srv 14.06.12
          ENDIF
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified electron temperature ',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENEPAR(IB,1).gt.0.0_R8,
     1     'BCENE = 1, ENEPAR(,1) !> 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*ENEPAR(IB,1)*EV
     .       *SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        VBIG*NE00*ENEPAR(IB,1)*EV*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       UBIG*NE00*ENEPAR(IB,1)*EV
     .       *SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .       -UBIG*NE00*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .        UBIG*NE00*ENEPAR(IB,1)*EV*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-UBIG*NE00*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'X') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= UBIG*NE00*
     .       ENEPAR(IB,1)*EV*VOL(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-UBIG*NE00*
     .       VOL(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'Y') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= UBIG*NE00*
     .       ENEPAR(IB,1)*EV*VOL(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-UBIG*NE00*
     .       VOL(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCENE=1, INVALID OPTION.')
        ENDIF

        GOTO 30000

* -- BCENE=2 -- PRESCRIBE THE GRADIENT OF THE ELECTRON TEMPERATURE

20200   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'electron temperature gradient ',ENEPAR(IB,1),
     1     ' on',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*EV*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .       (-ENEPAR(IB,1)*0.5_R8*
     .       (HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     .        HY1(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))) +
     .       TE(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       - VBIG*NE00*EV*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*EV*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .       (ENEPAR(IB,1)*0.5_R8*
     .       (HY1(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     .        HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))) +
     .       TE(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*EV*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .        VBIG*NE00*EV*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     .       (-ENEPAR(IB,1)*0.5_R8*
     .       (HX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))+
     .        HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .           RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))) +
     .       TE(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .       -VBIG*NE00*EV*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= VBIG*NE00*EV*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))*
     .       (ENEPAR(IB,1)*0.5_R8*
     .       (HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .           LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))+
     .        HX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))) +
     .       TE(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-VBIG*NE00*EV*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          end do
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCENE=2, INVALID OPTION.')
        END IF
        GOTO 30000

* -- BCENE=3 -- SHEATH CONDITIONS, ELECTRON ENERGY TRANSMISSION

20300   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,g14.7,a,g14.7,a,a,a)')
     .     'sheath electron with MOMPAR(ismain,IB,2) =',
     1     MOMPAR(ismain,IB,2),' and MOMPAR(ismain,IB,1) =',
     2     MOMPAR(ismain,IB,1),' and POTPAR(IB,2) =',POTPAR(IB,2),
     3     ' and ENEPAR(IB,1) =',ENEPAR(IB,1),
     4     ' on ',BCCHAR(IB),boundary_location(ib)
          if(fix_fch_in_fhe_sheath.eq.0 .and.
     &       fix_sx_in_fhe_sheath.eq.0) then !srv 08.05.09
            write(*,*) 'WARNING: Wrong (old) implementation of ',
     2       'the FCH contribution to FHE used'
          else if (fix_fch_in_fhe_sheath.eq.1 .and.
     &       fix_sx_in_fhe_sheath.eq.0) then
            write(*,*) 'WARNING: Wrong (second) implementation of ',
     2       'the FCH contribution to FHE used'
          endif
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
              CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))  !srv 01.02.09
              UUOUT=abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .                  wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3))*
     .         MOMPAR(ismain,IB,1)*CS
              if(fix_fch_in_fhe_sheath.eq.0) then
                NUEOUT=MIN(0.0_R8,-UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/QE)
              elseif(fix_fch_in_fhe_sheath.eq.1 .and.
     .               fix_sx_in_fhe_sheath.eq.0) then
                NUEOUT=MIN(0.0_R8,-UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)/QE)
              elseif(fix_fch_in_fhe_sheath.eq.2 .or.
     .               fix_sx_in_fhe_sheath.eq.1) then                     !srv 08.05.09 {
                NUEOUT=MIN(0.0_R8,-UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)/QE/
     .            sx(bc_list_x(iy,ib),bc_list_y(iy,ib)))
              endif
            ELSE
              NUEOUT=0.0_R8
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  ZBAR=bcint(na(-1,-1,IS),nx,ny,IY,IB)/
     .             bcint(NE,nx,ny,IY,IB)
                  ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .             bcint(rz2(-1,-1,IS),nx,ny,IY,IB)*ZBAR*
     .             bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                  CS=SQRT(ZARG)
                  UUTEST=
     .             UU(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .             (1.0_R8+
     .              HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .      HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .         RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .      UU(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .         RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .              HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .      HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .         RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
                  UUOUT=MAX( MOMPAR(IS,IB,1)*CS*
     .             abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .                 wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .             -UUTEST)
                  NUEOUT=NUEOUT - UUOUT*bcint(rza(-1,-1,IS),nx,ny,IY,IB)
     .             *bcint(na(-1,-1,IS),nx,ny,IY,IB)
                endif
              end do
              DO IS=0,NS_EXT-1
                if (.not.is_ext_neutral(is)) then
                  NUEOUT=NUEOUT +
     .             UU(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS+NS)
     .             *bcint(za_ext(-1,-1,IS),nx,ny,IY,IB)
     .             *bcint(na_ext(-1,-1,IS),nx,ny,IY,IB)
                endif
              end do
              if(fix_fch_in_fhe_sheath.eq.0) then
                NUEOUT=MIN(0.0_R8, NUEOUT
     .           - fch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/QE )
              elseif(fix_fch_in_fhe_sheath.eq.1 .and.
     .               fix_sx_in_fhe_sheath.eq.0) then
                NUEOUT=MIN(0.0_R8, NUEOUT -
     &           fch(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)/QE )
              elseif(fix_fch_in_fhe_sheath.eq.2 .or.
     .               fix_sx_in_fhe_sheath.eq.1) then
                NUEOUT=MIN(0.0_R8, NUEOUT -
     &           fch(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)/QE/
     &            sx(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)) )
              endif
            ENDIF
            if(enepar(ib,1).gt.0.0_R8) then                             !xpb
              she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .        she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) +
     .         (ENEPAR(IB,1)+
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT*
     .          SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            else           ! xpb new pressure linearisation as per B2.5-style
              she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)=
     .        she0(bc_list_x(iy,ib),bc_list_y(iy,ib),3) +               !xpb
     .         (-ENEPAR(IB,1)+                                          !xpb
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT*
     .          SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     .          ne(bc_list_x(iy,ib),bc_list_y(iy,ib))
            endif                                                       !xpb
            if(mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &         she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &         fhePSch(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                 rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)
            endif                                                !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
              CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))  !srv 01.02.09
              UUOUT=MOMPAR(ismain,IB,1)*CS*
     .         abs(wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .             wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),3))
              if(fix_fch_in_fhe_sheath.eq.0) then
                NUEOUT=MAX(0.0_R8,UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)/QE)
              elseif(fix_fch_in_fhe_sheath.eq.1 .and.
     .               fix_sx_in_fhe_sheath.eq.0) then
                NUEOUT=MAX(0.0_R8,UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/QE)
              elseif(fix_fch_in_fhe_sheath.eq.2 .or.
     .               fix_sx_in_fhe_sheath.eq.1) then
                NUEOUT=MAX(0.0_R8,UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/QE/
     .            sx(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
              endif
            ELSE
              NUEOUT=0.0_R8
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  ZBAR=bcint(na(-1,-1,IS),nx,ny,IY,IB)/
     .             bcint(NE,nx,ny,IY,IB)
                  ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .             bcint(rz2(-1,-1,IS),nx,ny,IY,IB)*ZBAR*
     .             bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                  CS=SQRT(ZARG)
                  UUTEST=
     .             UU(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .             (1.0_R8+
     .              HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .       HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .       UU(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .              HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .       HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
                  UUOUT=MAX( MOMPAR(IS,IB,1)*CS*
     .             abs(wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .                 wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .             +UUTEST)
                  NUEOUT=NUEOUT + UUOUT*bcint(rza(-1,-1,IS),nx,ny,IY,IB)
     .             *bcint(na(-1,-1,IS),nx,ny,IY,IB)
                endif
              end do
              DO IS=0,NS_EXT-1
                if (.not.is_ext_neutral(is)) then
                  NUEOUT=NUEOUT +
     .              UU(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),IS+NS)
     .             *bcint(za_ext(-1,-1,IS),nx,ny,IY,IB)
     .             *bcint(na_ext(-1,-1,IS),nx,ny,IY,IB)
                endif
              end do
              if(fix_fch_in_fhe_sheath.eq.0) then
                NUEOUT=MAX(0.0_R8, NUEOUT -
     .           fch(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)/QE )
              elseif(fix_fch_in_fhe_sheath.eq.1 .and.
     .               fix_sx_in_fhe_sheath.eq.0) then
                NUEOUT=MAX(0.0_R8,
     .           NUEOUT - fch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/QE )
              elseif(fix_fch_in_fhe_sheath.eq.2 .or.
     .               fix_sx_in_fhe_sheath.eq.1) then
                NUEOUT=MAX(0.0_R8,
     .           NUEOUT - fch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/QE/
     .           SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))) )
              endif
            ENDIF
            if(enepar(ib,1).gt.0.0_R8) then                             !xpb
              she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .        she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) -
     .         (ENEPAR(IB,1) +
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT*
     .          SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .             LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            else           ! xpb new pressure linearisation as per B2.5-style
              she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)=
     .        she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3) -              !xpb
     .         (-ENEPAR(IB,1) +                                        !xpb
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT*
     .          SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .             LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
     .         /ne(bc_list_x(iy,ib),bc_list_y(iy,ib))
            endif                                                      !xpb
            if(mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &         she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) +
     &         fhePSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)
            endif                                                !srv 12.09.08 }
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=3, NO SHEATH IN RADIAL DIR.')
        END IF
        GOTO 30000

* -- BCENE=4 -- PRESCRIBE THE VALUE OF THE ELECTRON TEMPERATURE,
*                 WEAKLY A MIXED BOUNDARY CONDITION

20400   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'weak electron temperature ',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENEPAR(IB,1).gt.0.0_R8,
     1     'BCENE = 4, ENEPAR(,1) !> 0')
          call xertst(ENEPAR(IB,2).gt.0.0_R8,
     1     'BCENE = 4, ENEPAR(,2) !> 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       ENEPAR(IB,2)*ENEPAR(IB,1)*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -ENEPAR(IB,2)/EV*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       ENEPAR(IB,2)*ENEPAR(IB,1)*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-ENEPAR(IB,2)/EV*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       ENEPAR(IB,2)*ENEPAR(IB,1)*
     1       SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .       -ENEPAR(IB,2)/EV*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
           DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       ENEPAR(IB,2)*ENEPAR(IB,1)*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-ENEPAR(IB,2)/EV*
     2       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=4, NOT YET AVAILABLE')
        END IF
        GOTO 30000

* -- BCENE=5 -- PRESCRIBE THE ELECTRON ENERGY FLUX PER UNIT AREA

20500   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'electron energy flux per area ',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENEPAR(IB,1)*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENEPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       ENEPAR(IB,1)*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=ENEPAR(IB,1)*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=5, NOT YET AVAILABLE')
        END IF
        GOTO 30000

* -- BCENE=6 -- PRESCRIBE THE TOTAL ELECTRON ENERGY FLUX FOR A
*                 CONSTANT ELECTRON TEMPERATURE

cghgh
20600   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'total electron energy flux ',ENEPAR(IB,1),
     1     ' with const Te on ',BCCHAR(IB),boundary_location(ib)
        endif
#ifdef WANT_THIS
        IF(BCCHAR(IB).EQ.'N') THEN
          Tguess=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            Tguess=Tguess+TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
          Tguess=Tguess/BC_LIST_SIZE(IB)
          if(diagno.ge.2) write(*,*) 'Initial Tguess:',Tguess/ev
          T0=Tguess*0.95_R8
          T1=Tguess*1.01_R8
          Tlimit=Tguess*0.0001_R8
          fylimt=ENEPAR(IB,1)*0.0001_R8
          ibnd=-1
20610     ibnd=ibnd+1
          if (ibnd.eq.0) then
            Tguess=T0
          elseif (ibnd.eq.1) then
            Tguess=T1
          else
            Tguess=T1-GAM1*(T1-T0)/(EPSLON + GAM1 - GAM0)
          endif
          fyt=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
cxpb code segment analog to b2tfhe
            if(ncall.eq.0) then
             vedia(bc_list_x(ix,ib),bc_list_y(ix,ib),1) =
     .        ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .           bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))*
     &        bb(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .           bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2)       !xpb
     &        /(qe*hx(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))   !xpb
     &        *(2.0_R8/
     &          bb(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)**3) !xpb
     &        *(wbbr(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)-
     &          wbbl(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3))
     &        *facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))                                   !xpb
             veecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1) =
     .         1.0_R8/hx(bc_list_x(ix,ib),bc_list_y(ix,ib))*          !xpb
     1          bb(bc_list_x(ix,ib),bc_list_y(ix,ib),2)/
     .          bb(bc_list_x(ix,ib),bc_list_y(ix,ib),3)**2*               !xpb
     2         (po(rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             rightiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     3          po(leftix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             leftiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))*0.5_R8*
     &         fac_ExB(bc_list_x(ix,ib),bc_list_y(ix,ib))                                    !xpb
            endif
!xpb            if(no_current.eq.1) then
!xpb              ehy = 0.0_R8
!xpb            else
            ehy = -(po(bc_list_x(ix,ib),bc_list_y(ix,ib))-
     &     po(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .        bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))+
     &    (ne(bc_list_x(ix,ib),bc_list_y(ix,ib))*Tguess-
     &     ne(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .        bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))*
     &     te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .        bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))*
     &    (2.0_R8/qe)/(ne(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                    bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &                +ne(bc_list_x(ix,ib),bc_list_y(ix,ib)))
!xpb            endif
            if (conduction_only.eq.1) then
              flo = 0.0_R8
            else
              prefactor=
     &         (vol(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &          vol(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &         /(hy1(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &           hy1(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &         *(ne(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &          +ne(bc_list_x(ix,ib),bc_list_y(ix,ib)))/2.0_R8
              flo_vedia=prefactor*
     .         vedia(bc_list_x(ix,ib),bc_list_y(ix,ib),1)
              flo_veecrb=prefactor*
     .         veecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1)
              flo = (1.5_R8+BoRiS)*
     .         fne_32(bc_list_x(ix,ib),bc_list_y(ix,ib),1)+
     1         2.5_R8*fne_52(bc_list_x(ix,ib),bc_list_y(ix,ib),1)+
     &       chve(bc_list_x(ix,ib),bc_list_y(ix,ib),1)*
     .       (ne(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     &       ne(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .          bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))/2.0_R8+
     &         2.5_R8*flo_vedia+1.5_R8*flo_veecrb
            endif
            if (no_hybr_e.eq.0) then
              cone0 =
     .         hybr(flo,chce(bc_list_x(ix,ib),bc_list_y(ix,ib),1))
              fy=flo*(te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &         Tguess)/2.0_R8-cone0*(Tguess-
     &         te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .            bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))+
     &         calf(bc_list_x(ix,ib),bc_list_y(ix,ib),1)*ehy*
     &         harm(te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     &         Tguess)
              if (fch_anomalous.ne.0.0_R8) fy = fy + fch_anomalous*
     &         calf_an(bc_list_x(ix,ib),bc_list_y(ix,ib),1)*ehy*
     &         harm(te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     &         Tguess)
            else
              cone0 = chce(bc_list_x(ix,ib),bc_list_y(ix,ib),1)
              t2 = dim ( cone0,
     .         abs(flo) *
     .         max(0.0_R8,
     .          min(Tguess,
     .              te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))))/
     .         max(cutlo, Tguess+
     .          te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))) )
              fy =
     .         upwind(flo,
     .                te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     .                Tguess)
     .         - t2 * (Tguess-te(bottomix(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)),
     .                           bottomiy(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)))) +
     &         calf(bc_list_x(ix,ib),bc_list_y(ix,ib),1)*
     .         harm(te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     &                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     &                            Tguess)*ehy
              if (fch_anomalous.ne.0.0_R8) fy = fy + fch_anomalous*
     &         calf_an(bc_list_x(ix,ib),bc_list_y(ix,ib),1)*
     .         harm(te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     &                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     &                            Tguess)*ehy
            endif
            fyt=fyt + fy
          ENDDO
          if(diagno.ge.2)
     >     write(*,'(a,i4,f15.6,1pe14.6)') 'ibnd,Tg,fyt(e):',
     .     ibnd,Tguess/ev,fyt
          if (ibnd.eq.0) then
            GAM0 = fyt - ENEPAR(IB,1)
            goto 20610
          elseif (ibnd.eq.1) then
            GAM1 = fyt - ENEPAR(IB,1)
            goto 20610
          else
            T0 = T1
            GAM0 = GAM1 - ENEPAR(IB,1)
            T1 = Tguess
            GAM1 = fyt
            if(diagno.ge.2) write(*,*) 'GAM0, GAM1 = ',GAM0,GAM1
            if(diagno.ge.2) write(*,*) 'fyt,t1,t0 = ',fyt,t1/ev,t0/ev
            if (ibnd .gt. 100)
     >       call xerrab ('B2STBC -- BCENE=6, NOT CONVERGING')
            if (abs(t1-t0).gt.Tlimit .and.
     &          abs(fyt-ENEPAR(IB,1)).gt.fylimt) goto 20610
          endif
          if(diagno.ge.1)
     >     write(*,'(a,f15.6,i4,1pe14.6)') ' tebound:',T1/ev,ibnd,fyt
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*T1*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'S') THEN
          Tguess=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            Tguess=Tguess+TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
          Tguess=Tguess/BC_LIST_SIZE(IB)
          if(diagno.ge.2) write(*,*) 'Initial Tguess:',Tguess/ev
          T0=Tguess*0.95_R8
          T1=Tguess*1.01_R8
          Tlimit=Tguess*0.0001_R8
          fylimt=ENEPAR(IB,1)*0.0001_R8
          ibnd=-1
20620     ibnd=ibnd+1
          if (ibnd.eq.0) then
            Tguess=T0
          elseif (ibnd.eq.1) then
            Tguess=T1
          else
            Tguess=T1-GAM1*(T1-T0)/(EPSLON + GAM1 - GAM0)
          endif
          fyt=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
cxpb code segment analog to b2tfhe
            if(ncall.eq.0) then
             vedia(bc_list_x(ix,ib),bc_list_y(ix,ib),1) =                 !xpb
     &        ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .           topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))*               !xpb
     &        bb(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .           topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2)/             !xpb
     &        (qe*hx(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))           !xpb
     &        *(2.0_R8/bb(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .          topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)**3) !xpb
     &        *(wbbr(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)-         !xpb
     &          wbbl(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3))         !xpb
     &        *facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))                                    !xpb
             veecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1) =
     .         1.0_R8/hx(bc_list_x(ix,ib),bc_list_y(ix,ib))*           !xpb
     1          bb(bc_list_x(ix,ib),bc_list_y(ix,ib),2)/
     .          bb(bc_list_x(ix,ib),bc_list_y(ix,ib),3)**2*                !xpb
     2         (po(rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             rightiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     3          po(leftix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             leftiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))*0.5_R8*
     &         fac_ExB(bc_list_x(ix,ib),bc_list_y(ix,ib))                                     !xpb
            endif
!xpb            if(no_current.eq.1) then
!xpb               ehy = 0.0_R8
!xpb            else
            ehy = -(po(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     &                po(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &         (ne(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))*
     &          te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     &          ne(bc_list_x(ix,ib),bc_list_y(ix,ib))*Tguess)*
     &         (2.0_R8/qe)/(ne(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                         topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &                     +ne(bc_list_x(ix,ib),bc_list_y(ix,ib)))
!xpb            endif
            if (conduction_only.eq.1) then
              flo = 0.0_R8
            else
              prefactor=
     &         (vol(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &          vol(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &         /(hy1(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &           hy1(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &         *(ne(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))                    !xpb
     &          +ne(bc_list_x(ix,ib),bc_list_y(ix,ib)))/2.0_R8                                       !xpb
              flo_vedia=prefactor*
     .         vedia(bc_list_x(ix,ib),bc_list_y(ix,ib),1)
              flo_veecrb=prefactor*
     .         veecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1)
              flo = (1.5_R8+BoRiS)*
     .         fne_32(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)+
     1         2.5_R8*fne_52(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)+
     &         chve(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)*
     &        (ne(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     &         ne(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .            topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))/2.0_R8+
     &         2.5_R8*flo_vedia+1.5_R8*flo_veecrb
            endif
            if (no_hybr_e.eq.0) then
              cone0 = hybr(flo,
     &         chce(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1))
              fy=flo*(te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &         Tguess)/2.0_R8-cone0*
     &         (te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-Tguess)+
     &         calf(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)*ehy*
     &         harm(te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),Tguess)
              if (fch_anomalous.ne.0.0_R8) fy = fy + fch_anomalous*
     &         calf_an(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)*ehy*
     &         harm(te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),Tguess)
            else
              cone0 = chce(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                     topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)
              t2 = dim ( cone0,
     .         abs(flo) *
     .         max(0.0_R8,
     .          min(Tguess,
     .              te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))))/
     .         max(cutlo, Tguess+
     .          te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))) )
              fy =
     .         upwind(flo, Tguess,
     .                te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .         - t2 * (te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                    topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     .                 Tguess) +
     &         calf(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)*
     &         harm(te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     &              Tguess)*ehy
              if (fch_anomalous.ne.0.0_R8) fy = fy + fch_anomalous*
     &         calf_an(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     &                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)*
     &         harm(te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     &              Tguess)*ehy
            endif
            fyt=fyt + fy
          ENDDO
          if(diagno.ge.2)
     >     write(*,'(a,i4,f15.6,1pe14.6)') 'ibnd,Tg,fyt(e):',
     .     ibnd,Tguess/ev,fyt
          if (ibnd.eq.0) then
            GAM0 = fyt - ENEPAR(IB,1)
            goto 20620
          elseif (ibnd.eq.1) then
            GAM1 = fyt - ENEPAR(IB,1)
            goto 20620
          else
            T0 = T1
            GAM0 = GAM1 - ENEPAR(IB,1)
            T1 = Tguess
            GAM1 = fyt
            if (ibnd .gt. 100)
     >       call xerrab ('B2STBC -- BCENE=6, NOT CONVERGING')
            if (abs(t1-t0).gt.Tlimit .and.
     &          abs(fyt-ENEPAR(IB,1)).gt.fylimt) goto 20620
          endif
          if(diagno.ge.1)
     >     write(*,'(a,f15.6,i4,1pe14.6)') ' tebound:',T1/ev,ibnd,fyt
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=6, NOT AVAILABLE FOR W,E')
        ENDIF
#else
        CALL XERRAB ('B2STBC -- BCENE=6, NOT IMPLEMENTED.')
#endif
        GOTO 30000
cghgh

* -- BCENE=7 -- PRESCRIBE THE ELECTRON TEMPERATURE AS A FUNCTION OF
*                 OTHER PLASMA PARAMETERS

20700   CONTINUE
        CALL XERRAB ('B2STBC -- BCENE=7, NOT YET AVAILABLE')
        GOTO 30000

* -- BCENE=8 -- PRESCRIBE THE TOTAL ELECTRON HEAT FLUX WITH CONSTANT
*                 FLUX DENSITY

20800   CONTINUE
        if(ncall.eq.0) then
          if (use_astra.ne.0) ENEPAR(IB,1) = fheb_astra           !srv 28.07.11
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'total electron energy flux ',ENEPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          US=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            US=US+SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       ENEPAR(IB,1)/US*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
            if(mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     .         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &         fhePSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     .         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     &         1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &         veecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)*
     &                   ne(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          endif                                                !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          US=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            US=US+SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .               BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENEPAR(IB,1)/US*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
            if(mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &         fhePSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) -
     &         1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &         veecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)*
     &         ne(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          endif                                                !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          US=0.0_R8
          DO IY = 1, BC_LIST_SIZE(IB)
            US=US+SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       ENEPAR(IB,1)/US*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          US=0.0_R8
          DO IY = 1, BC_LIST_SIZE(IB)
            US=US+SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          end do
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=ENEPAR(IB,1)/US*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=8, NOT AVAILABLE')
        ENDIF
        GOTO 30000

* -- BCENE=9 -- PRESCRIBE THE DECAY LENGTH FOR THE ELECTRON TEMPERATURE

20900   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'electron temperature decay length ',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENEPAR(IB,1).ne.0.0_R8,
     1     'BCENE = 9, ENEPAR(,1) = 0')
!          call xertst(all(BCCON(0:NS-1,IB).EQ.9),
!     1     'Should use BCCON(IS,IB).EQ.9 together with BCENE=9 ')       !srv 16.07.12
          ism = 0
          do is = 0, ns-1
            if (bccon(is,ib).ne.9.and..not.is_neutral(is)) then
              ism = ism + 1
              write(0,*)
     &         'Should use BCCON(IS,IB).EQ.9 together with BCENE=9'
              write(*,*)
     &         'Should use BCCON(IS,IB).EQ.9 together with BCENE=9'
            endif
            if (ism.gt.0) then
              write(*,*)
     &         'Please ensure that this boundary condition '//
     &         'is consistent with density boundary conditions'  !srv 21.09.12
            end if
          enddo
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -hce0(bc_list_x(ix,ib),bc_list_y(ix,ib))/ENEPAR(IB,1)*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            DO IS=0,NS-1
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .         rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*2.5_R8*
     .         TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .        she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     .         rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*2.5_R8*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            ENDDO
            if(mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &         fhePSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     &         1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &         veecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)*
     &                   ne(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            endif                                                !srv 12.09.08 }
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -hce0(bc_list_x(ix,ib),bc_list_y(ix,ib))/ENEPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            DO IS=0,NS-1
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .         rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*2.5_R8*
     .         TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .        she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     .         rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*2.5_R8*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            ENDDO
            if(mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     .         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &         fhePSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     .         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) -
     &         1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &         veecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)*
     &             ne(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            endif                                                !srv 12.09.08 }
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=9, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 30000

* -- BCENE=10 -- Feedback option for core

21000   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,a,a)')
     .     'electron heat flux feedback bc on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(lfeedback, 'LFEEDBACK must be turned on!')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if (region(bc_list_x(ix,ib),bc_list_y(ix,ib),0).eq.1) then
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .         cbshe(0,coreregno)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            elseif (region(bc_list_x(ix,ib),bc_list_y(ix,ib),0).eq.5)
     .       then
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .         cbshe(0,coreregn2)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            endif
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=10, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 30000

* -- BCENE=11 -- FROM B2STBC_SPB                                       !srv 01.02.09 {

21100   CONTINUE
        if(ncall.eq.0) then   !srv added accumulation in order to account twice b.c. for corner cells
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     ' FROM B2STBC_SPB with enepar1 =',ENEPAR(IB,1),
     1     ' and enepar2 =',ENEPAR(IB,2),
     2     ' on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
!            vte = sqrt(bcint(te,nx,ny,ix,ib)/me)
            vte = sqrt(te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/me)
            s0 = vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            if(mdf_fhe.ne.0) then                                    !srv 21.03.09 {
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &          fhePSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)
            endif                                                    !srv 21.03.09 }
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &       she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)+ENEPAR(IB,1)*
     &       max(0.0e0_R8,-floe_noc(
     &         topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &         topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1))       !srv 28.10.02
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)=
     &       she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)+ENEPAR(IB,2)*
     &       vte*s0
!     &       vte*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
!            vte = sqrt(bcint(te,nx,ny,ix,ib)/me)
            vte = sqrt(te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/me)
            s0 = vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            if(mdf_fhe.ne.0) then                                    !srv 21.03.09 {
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &         fhePSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)
          endif                                                    !srv 21.03.09 }
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &       she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)+ENEPAR(IB,1)*
     &       max(0.0e0_R8,floe_noc(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1))                  !srv 28.10.02
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)=
     &       she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)+ENEPAR(IB,2)*
     &       vte*s0
!     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
!     2          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=11, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 30000                                                    !srv 01.02.09 }

* -- BCENE=12 -- SHEATH CONDITIONS, EL ENERGY TRANSMISSION COEFFICIENT

21200   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,a,a)')
     .     'sheath (12) el on ',BCCHAR(IB),
     1     boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(bc_list_x(ix,ib),bc_list_y(ix,ib),3) =
     .      she0(bc_list_x(ix,ib),bc_list_y(ix,ib),3)-
     1       ENEPAR(IB,1)*max(0.0_R8,
     1        -fne(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1))/
     2       ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(bc_list_x(ix,ib),bc_list_y(ix,ib),3) =
     .      she0(bc_list_x(ix,ib),bc_list_y(ix,ib),3)-
     1             ENEPAR(IB,1)*max(0.0_R8,
     .       fne(bc_list_x(ix,ib),bc_list_y(ix,ib),1))/
     2       ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(bc_list_x(iy,ib),bc_list_y(iy,ib),3) =
     .      she0(bc_list_x(iy,ib),bc_list_y(iy,ib),3)-
     1       ENEPAR(IB,1)*max(0.0_R8,
     1        -fne(rightix(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .             rightiy(bc_list_x(iy,ib),bc_list_y(iy,ib)),0))/
     2       ne(bc_list_x(iy,ib),bc_list_y(iy,ib))
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(bc_list_x(iy,ib),bc_list_y(iy,ib),3) =
     .      she0(bc_list_x(iy,ib),bc_list_y(iy,ib),3)-
     1             ENEPAR(IB,1)*max(0.0_R8,
     .       fne(bc_list_x(iy,ib),bc_list_y(iy,ib),0))/
     2       ne(bc_list_x(iy,ib),bc_list_y(iy,ib))
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=12, INVALID OPTION')
        END IF
        GOTO 30000

* -- BCENE=13 -- PRESCRIBE THE ELECTRON ENERGY FLUX PER UNIT AREA
*             -- PROPORTIONAL TO TEMPERATURE

21300   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     1     'electron energy flux per area proportional to temperature ',
     2     ENEPAR(IB,1),' on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       ENEPAR(IB,1)*TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       ENEPAR(IB,1)*TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       ENEPAR(IB,1)*TE(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     1       SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       ENEPAR(IB,1)*TE(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=13, NOT YET AVAILABLE')
        END IF
        GOTO 30000

* -- BCENE=14 -- leakage option for electron energy

21400   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'electron energy leakage',ENEPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            SHE0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            CS=SQRT(bcint(PR,nx,ny,IX,IB)/bcint(RO,nx,ny,IX,IB))   !srv 01.02.09
            SHE0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     1       vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     2       hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*cs*ENEPAR(IB,1)
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            SHE0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            CS=SQRT(bcint(PR,nx,ny,IX,IB)/bcint(RO,nx,ny,IX,IB))   !srv 01.02.09
            SHE0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     1       vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     2       hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*cs*ENEPAR(IB,1)
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=14, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 30000

* -- BCENE=15 -- FROM B2STBC_SPB                                     !srv 01.02.09 {

21500   CONTINUE
        if(bcene_15_style.eq.0) then
          if(ncall.eq.0) then   !srv added accumulation in order to account twice b.c. for corner cells
            write(*,'(a,1p,g14.7,a,a,a)')
     .       'bcene_15_style = 0 with GAMMAE =',GAMMAE,
     1       ' on ',BCCHAR(IB),boundary_location(ib)
          endif
          IF(BCCHAR(IB).EQ.'W') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
!            CS=SQRT(bcint(PR,nx,ny,IX,IB)/bcint(RO,nx,ny,IX,IB))
              CS=SQRT(PR(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                RO(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              t0 = max(0.0e0_R8,                                        !srv 05.05.09
     &         (abs(pbs(rightix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  rightiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),0))*
     &                       ne(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*cs +   !srv 25.12.06 11.01.13
     &          fch_p(rightix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                rightiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),0)/qe))                  !srv 25.12.06 11.03.09 }
              if(pot_eq.eq.1) then                                      !srv 27.02.13 {
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-             !srv 11.03.09 {
     &                (1.0e0_R8+GAMMAE)/(1.0e0_R8-GAMMAE)*t0                !srv 19.04.13
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &           qe*po(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*t0            !srv 07.11.12
              else                            !srv 27.02.13
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-             !srv 11.03.09 {
     &                ((1.0e0_R8+GAMMAE)/(1.0e0_R8-GAMMAE)+3.1_R8)*t0   !srv 19.04.13
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)
              endif                                                     !srv 27.02.13 }
              if(mdf_fhe.ne.0) then                                     !srv 01.06.06 {
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &           fhePSch(rightix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                   rightiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),0)
            endif                                                       !srv 01.06.06 }
            END DO
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
!              CS=SQRT(bcint(PR,nx,ny,IX,IB)/bcint(RO,nx,ny,IX,IB))
              CS=SQRT(PR(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                RO(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              t0 = max(0.0e0_R8,                                        !srv 05.05.09
     &         (abs(pbs(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0))*
     &               ne(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*cs -         !srv 25.12.06 11.01.13
     &            fch_p(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)/qe))       !srv 25.12.06 11.03.09 }
              if(pot_eq.eq.1) then                                        !srv 27.02.13 {
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-             !srv 11.03.09 {
     &                (1.0e0_R8+GAMMAE)/(1.0e0_R8-GAMMAE)*t0            !srv 19.04.13
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &           qe*po(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*t0             !srv 07.11.12
              else
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-             !srv 11.03.09 {
     &               ((1.0e0_R8+GAMMAE)/(1.0e0_R8-GAMMAE)+3.1_R8)*t0        !srv 19.04.13
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)
              endif                                                     !srv 27.02.13 }
              if(mdf_fhe.ne.0) then                                     !srv 01.06.06 {
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &           fhePSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)
            endif                                                     !srv 01.06.06 }
            end do
          ELSE
            CALL XERRAB 
     .       ('B2STBC -- BCENE=15, not available for S and N')
          END IF

        elseif(bcene_15_style.eq.1) then
          if(ncall.eq.0) then   !srv added accumulation in order to account twice b.c. for corner cells
            write(*,'(a,1p,g14.7,a,a,a)')
     &       'bcene_15_style = 1 with GAMMAE =',GAMMAE,
     1       ' on ',BCCHAR(IB),boundary_location(ib)
          endif
          IF(BCCHAR(IB).EQ.'W') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
!              CS=SQRT(bcint(PR,nx,ny,IX,IB)/bcint(RO,nx,ny,IX,IB))
              s1=abs(pbs(rightix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                   rightiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),0)) !srv 18.04.13 {
              vte=sqrt(te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/me)
              t0 =s1*ne(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &          sqrt(1.0_R8/(2.0_R8*pi))*vte*
     &          expu2(max(-50.0_R8,
     &         -qe*(po(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))-potpar(ib,2))/
     &              te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
              if(pot_eq.eq.1) then                                      !srv 27.02.13 {
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &            (1.0e0_R8+GAMMAE)*t0
                 she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &            (1.0e0_R8-GAMMAE)*
     &            qe*po(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*t0
              else                        !srv 27.02.13
                cs=sqrt(pr(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                  ro(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &               s1*((1.0e0_R8+GAMMAE)+(1.0e0_R8-GAMMAE)*3.1_R8)*
     &           ne(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*cs
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)              !srv 18.04.13 }
              endif                                                     !srv 27.02.13 }
              if(mdf_fhe.ne.0) then                                     !srv 01.06.06 {
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &           fhePSch(rightix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                   rightiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),0)
            endif                                                     !srv 01.06.06 }
            END DO
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
!              CS=SQRT(bcint(PR,nx,ny,IX,IB)/bcint(RO,nx,ny,IX,IB))
              s1=abs(pbs(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0))                               !srv 18.04.13 {
              vte=sqrt(te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/me)
              t0 =s1*ne(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &          sqrt(1.0_R8/(2.0_R8*pi))*vte*
     &          expu2(max(-50.0_R8,
     &         -qe*(po(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))-potpar(ib,2))/
     &              te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
              if(pot_eq.eq.1) then                                      !srv 27.02.13 {
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &               (1.0e0_R8+GAMMAE)*t0
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &           (1.0e0_R8-GAMMAE)*qe*
     &           po(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*t0
              else                               !srv 27.02.13
                cs=sqrt(pr(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                  ro(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &               s1*((1.0e0_R8+GAMMAE)+(1.0e0_R8-GAMMAE)*3.1_R8)*
     &           ne(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*cs
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)            !srv 18.04.13 }
              endif                                                   !srv 27.02.13 }
              if(mdf_fhe.ne.0) then                                   !srv 01.06.06 {
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &           fhePSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)
            endif                                                     !srv 01.06.06 }
            end do
          ELSE
            CALL XERRAB 
     .       ('B2STBC -- BCENE=15, not available for S and N')
          END IF
        else
          CALL XERRAB
     .     ('B2STBC -- BCENE=15, incorrect bcene_15_style value')
        endif
        GOTO 30000                                                   !srv 01.02.09 }

* -- BCENE=16 -- Feedback boundary condition with constant temperature

21600   CONTINUE

        if(bcene_16_style.eq.0) then

          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .       'electron temperature feedback by flux',
     1       ENEPAR(IB,1),' with ',ENEPAR(IB,2),' on ',
     1       BCCHAR(IB),boundary_location(ib)
            call xertst(enepar(ib,1).ne.0, 'faulty parameter enepar1')
            call xertst(enepar(ib,2).gt.0, 'faulty parameter enepar2')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            totflux = 0.0_R8
            t0 = 0.0_R8
            weight = 0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type16_ix=bc_list_x(ix,ib)
              bc_type16_iy=bc_list_y(ix,ib)
              do bc_type16_i=1,bc_type16_ref     ! how many steps away from the boundary?
                bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
                bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
                bc_type16_ix=bc_type16_ixx
                call xertst (-2.lt.bc_type16_ix.and.
     &                             bc_type16_ix.lt.nx+1,
     1           'stepped off x-grid in type 16 bc')
                call xertst (-2.lt.bc_type16_iy.and.
     &                             bc_type16_iy.lt.ny+1,
     1           'stepped off y-grid in type 16 bc')
              enddo
              totflux=totflux+
     .         fhe(bc_type16_ix,bc_type16_iy,1)
              t0 = t0 +
     .          te(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .         vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .          ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
              weight = weight + vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .                          ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
            END DO
            t0 = t0 / weight
            if(diagno.ge.2) write(*,*) 'Totflux :',totflux
            T1 = (1.0_R8 -
     .       enepar(ib,2)*(totflux-enepar(ib,1))/enepar(ib,1)) * t0
            if(diagno.ge.2)
     1       write(*,*) 'Old/New electron temp :',t0/qe, t1/qe
            DO IX = 1, BC_LIST_SIZE(IB)
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .         VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .        -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            END DO
          ELSE
            CALL XERRAB ('B2STBC -- BCENE=16, NOT AVAILABLE FOR N,W,E')
          ENDIF
          GOTO 30000

        elseif(bcene_16_style.eq.1)then

          if(ncall.eq.0) then
          write(*,*)
     &       'Given integrated electron heat flux with constant Te '
!          write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCENE(IBW).EQ.BCENE(IB))THEN
            write(*,'(a,a,a)') ' on ',BCCHAR(IBW),boundary_location(ibw)
                call xertst(
     .         mod(region(BC_LIST_X(1,IBW),BC_LIST_Y(1,IBW),0),4).eq.1,
     .           'BCENE = 16 only applicable to core boundaries !')
                call xertst(enepar(ibw,1).ne.0.0_R8,
     .           'faulty parameter enepar1')
            endif
          END DO
            if (use_astra.ne.0) then
#ifdef ASTRA
              SENEPAR = fheb_astra
#else
              call xerrab(
     .         'Compile with -DASTRA option to couple to ASTRA!')
#endif
            else
              SENEPAR = 0.0_R8
              DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCENE(IBW).EQ.BCENE(IB))THEN
                  SENEPAR = SENEPAR + ENEPAR(IBW,1)
              endif
            END DO
            endif
          write(*,*) ' Summed value is ', SENEPAR
          endif
        if(.not.once_bcene) then
            TES = 0.0_R8
            volsum = 0.0_R8
            DO IBW=IB,NBC
            if(BCCHAR(IBW).EQ.'S'.and.BCENE(IBW).EQ.BCENE(IB))THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                TES = TES +
     &                  TE(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                 VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  VOLSUM = VOLSUM +
     &                 VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                END DO
              endif
            END DO
            TES=TES/VOLSUM
          IF(BCCHAR(IB).EQ.'S') THEN
            totflux = 0.0_R8
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCENE(IBW).EQ.BCENE(IB)) THEN
                if(mdf_fhe.ne.0) then
                    DO IX = 1, BC_LIST_SIZE(IBW)
                    totflux = totflux +
     &           fhe_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
     &                +fac_ExB(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*      !srv 11.04.13 {
     &                 ((po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))-
     &            po(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))/
     &                  (hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &            hx(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))+
     &          (po(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &              rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))-
     &                   po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))/
     &          (hx(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &              rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))+
     &                   hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))*
     &                   ne(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                   te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                   SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                   bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),2)/
     &                   bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),3)**2   !srv 11.04.13 }
                    ENDDO
                else
                    DO IX = 1, BC_LIST_SIZE(IBW)
                totflux = totflux +
     .             fhe(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
                    ENDDO
                endif
              ENDIF
            ENDDO
              TES = TES *
     &         (1.0_R8-bc_ref_te*(totflux-SENEPAR)/                     !srv 11.04.13
     &           max(abs(totflux),SENEPAR))
              write(*,"(' BCENE=16: TES=',es14.7,' totflux=',es14.7,
     *         ' SENEPAR=',es14.7)") TES/ev, totflux, SENEPAR
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCENE(IBW).EQ.BCENE(IB)) THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    she0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)= 
     .               VBIG*NE00*TES*
     .               SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    she0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)=
     .              -VBIG*NE00*
     .               SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  ENDDO
              ENDIF
            END DO
            once_bcene = .true.
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCENE=16, NOT AVAILABLE FOR N,W,E')
            ENDIF
          endif                                                         !srv 16.02.12 }
            GOTO 30000
        else
          CALL XERRAB
     1       (' Incorrect style of boundary condition')
        endif

* -- BCENE=17 -- Feedback boundary condition with constant shared temperature

21700   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'shared temperature feedback by flux ',
     1     ENEPAR(IB,1),' with ',ENEPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(enepar(ib,1).ne.0, 'faulty parameter enepar1')
          call xertst(enepar(ib,2).gt.0, 'faulty parameter enepar2')
          call xertst(bcene(ib).eq.bceni(ib),
     1     'both bcene and bceni should be set to 17')
          call xertst(enipar(ib,1).ne.0, 'faulty parameter enipar1')
          call xertst(enipar(ib,2).gt.0, 'faulty parameter enipar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          totflux = 0.0_R8
          t0 = 0.0_R8
          weight = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type16_ix=bc_list_x(ix,ib)
            bc_type16_iy=bc_list_y(ix,ib)
            do bc_type16_i=1,bc_type16_ref      ! how many steps away from the boundary?
              bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
              bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
              bc_type16_ix=bc_type16_ixx
              call xertst (-2.lt.bc_type16_ix.and.bc_type16_ix.lt.nx+1,
     1         'stepped off x-grid in type 16 bc')
              call xertst (-2.lt.bc_type16_iy.and.bc_type16_iy.lt.ny+1,
     1         'stepped off y-grid in type 16 bc')
            enddo
            totflux=totflux+fhe(bc_type16_ix,bc_type16_iy,1)+
     .       fhi(bc_type16_ix,bc_type16_iy,1)+
     .       fhi_ext(bc_type16_ix,bc_type16_iy,1)
            if (bc_type16_kinetic_energy.eq.1) then
              do is = 0, ns-1
                totflux=totflux+(1.0_R8-BoRiS)*
     .           fna(bc_type16_ix,bc_type16_iy,1,is)*0.5_R8*
     .           (kinrgy(bc_type16_ix,bc_type16_iy,is)+
     .            kinrgy(bottomix(bc_type16_ix,bc_type16_iy),
     .                   bottomiy(bc_type16_ix,bc_type16_iy),is))
              enddo
              do is = 0, ns_ext-1
                totflux=totflux+(1.0_R8-BoRiS)*
     .           fa_ext(bc_type16_ix,bc_type16_iy,1,is)*0.5_R8*
     .           0.5_R8*am_ext(is)*mp*
     .           (ua_ext(bc_type16_ix,bc_type16_iy,is)**2+
     .            ua_ext(bottomix(bc_type16_ix,bc_type16_iy),
     .                   bottomiy(bc_type16_ix,bc_type16_iy),is)**2)
              enddo
            endif
            t0 = t0 +
     .        te(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .        ne(bc_list_x(ix,ib),bc_list_y(ix,ib)) +
     .        ti(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .        ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)
            weight = weight +
     1       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .        ne(bc_list_x(ix,ib),bc_list_y(ix,ib)) +
     2       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .        ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)
          ENDDO
          t0 = t0 / weight
          if(diagno.ge.2) write(*,*) 'Totflux :',totflux
          T1 = (1.0_R8 -
     .     enepar(ib,2)*
     1       (totflux*enepar(ib,1)/(enepar(ib,1)+enipar(ib,1))-
     2         enepar(ib,1))/enepar(ib,1) -
     .     enipar(ib,2)*
     1       (totflux*enipar(ib,1)/(enepar(ib,1)+enipar(ib,1))-
     2         enipar(ib,1))/enipar(ib,1)
     .     ) * t0
          if(diagno.ge.2)
     1     write(*,*) 'Old/New shared temp :',t0/qe, t1/qe
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=17, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 30000

* -- BCENE=18 -- FRACTIONAL DROP CONDITION (FOR EDGE2D COMPARISON)

21800   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'electron temperature fractional drop ',
     1     ENEPAR(IB,1),' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENEPAR(IB,1).lt.1.0_R8,
     1     'BCENE = 18, ENEPAR(,1) => 1')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENEPAR(IB,1)) *
     .       TE(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENEPAR(IB,1)) *
     .       TE(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*T1*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENEPAR(IB,1)) *
     .       TE(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .        VBIG*NE00*T1*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .       -VBIG*NE00*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENEPAR(IB,1)) *
     .       TE(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= VBIG*NE00*T1*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-VBIG*NE00*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=18, INVALID OPTION')
        ENDIF
        GOTO 30000

* -- BCENE=19 -- PRESCRIBE THE DECAY LENGTH FOR THE ELECTRON TEMPERATURE

21900   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'electron temperature decay length (new) ',
     .     ENEPAR(IB,1), ' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(ENEPAR(IB,1).ne.0.0_R8,
     1     'BCENE = 19, ENEPAR(,1) = 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -hce0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/ENEPAR(IB,1)*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            t0=0.0_R8
            DO IS=0,NS-1
               t0=t0-bcint(rza(-1,-1,is),nx,ny,IX,IB)*(
     .              1.5_R8*fna_32(topix(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),
     .                            topiy(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),1,is)+
     .              2.5_R8*fna_52(topix(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),
     .                            topiy(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),1,is))
            ENDDO
            DO IS=0,NS_EXT-1
               t0=t0-bcint(za_ext(-1,-1,is),nx,ny,IX,IB)*
     .              2.5_R8*fa_ext(topix(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),
     .                            topiy(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),1,is)
            ENDDO
            if(t0.gt.0.0_R8) then
               she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .              t0 * TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            else
               she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) + t0
            endif
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
             she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
             she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -hce0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/ENEPAR(IB,1)*
     1            SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .               BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
             t0=0.0_R8
             DO IS=0,NS-1
                t0=t0+bcint(rza(-1,-1,is),nx,ny,IX,IB)*(
     .               1.5_R8*fna_32(bc_list_x(ix,ib),
     .                             bc_list_y(ix,ib),1,is)+
     .               2.5_R8*fna_52(bc_list_x(ix,ib),
     .                             bc_list_y(ix,ib),1,is))
             ENDDO
             DO IS=0,NS_EXT-1
                t0=t0+bcint(za_ext(-1,-1,is),nx,ny,IX,IB)*
     .               2.5_R8*fa_ext(bc_list_x(ix,ib),
     .                             bc_list_y(ix,ib),1,is)
             ENDDO
             if(t0.gt.0.0_R8) then
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .          she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .              t0 * TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
             else
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .          she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) + t0
             endif
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=19, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 30000

* -- BCENE=20 -- Feedback boundary condition with constant temperature

22000   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'electron temperature feedback by flux ',
     1     ENEPAR(IB,1),' with ',ENEPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(enepar(ib,1).ne.0, 'faulty parameter enepar1')
          call xertst(enepar(ib,2).gt.0, 'faulty parameter enepar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          totflux = 0.0_R8
          t0 = 0.0_R8
          weight = 0.0_R8
          teave = 0.0_R8
          neave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type16_ix=bc_list_x(ix,ib)
            bc_type16_iy=bc_list_y(ix,ib)
            do bc_type16_i=1,bc_type16_ref      ! how many steps away from the boundary?
              bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
              bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
              bc_type16_ix=bc_type16_ixx
              call xertst (-2.lt.bc_type16_ix.and.bc_type16_ix.lt.nx+1,
     1         'stepped off x-grid in type 16 bc')
              call xertst (-2.lt.bc_type16_iy.and.bc_type16_iy.lt.ny+1,
     1         'stepped off y-grid in type 16 bc')
            enddo
            totflux=totflux+
     .       fhe(bc_type16_ix,bc_type16_iy,1)
            t0 = t0 +
     .       te(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
            weight = weight + vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .                        ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
            bc_type20_ix=bc_list_x(ix,ib)
            bc_type20_iy=bc_list_y(ix,ib)
            do bc_type20_i=1,bc_type20_ref      ! how many steps away from the boundary?
              bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
              bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
              bc_type20_ix=bc_type20_ixx
              call xertst (-2.lt.bc_type20_ix.and.bc_type20_ix.lt.nx+1,
     1         'stepped off x-grid in type 20 bc')
              call xertst (-2.lt.bc_type20_iy.and.bc_type20_iy.lt.ny+1,
     1         'stepped off y-grid in type 20 bc')
            enddo
            teave=teave+te(bc_type20_ix,bc_type20_iy)*
     1        ne(bc_type20_ix,bc_type20_iy)*
     2       vol(bc_type20_ix,bc_type20_iy)
            neave=neave+ne(bc_type20_ix,bc_type20_iy)*
     1       vol(bc_type20_ix,bc_type20_iy)
          ENDDO
          t0 = t0 / weight
          teave=teave/neave
          if(diagno.ge.2) write(*,*) 'Totflux :',totflux
          T1 = (1.0_R8 -
     .     enepar(ib,2)*(totflux-enepar(ib,1))/enepar(ib,1)) * t0
          if(diagno.ge.2)
     1     write(*,*) 'Old/New electron temp :',t0/qe, t1/qe
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type20_ix=bc_list_x(ix,ib)
            bc_type20_iy=bc_list_y(ix,ib)
            do bc_type20_i=1,bc_type20_ref      ! how many steps away from the boundary?
              bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
              bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
              bc_type20_ix=bc_type20_ixx
              call xertst (-2.lt.bc_type20_ix.and.bc_type20_ix.lt.nx+1,
     1         'stepped off x-grid in type 20 bc')
              call xertst (-2.lt.bc_type20_iy.and.bc_type20_iy.lt.ny+1,
     1         'stepped off y-grid in type 20 bc')
            enddo
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*T1*
     1       te(bc_type20_ix,bc_type20_iy)/teave
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=20, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 30000

* -- BCENE=21 -- constant temperature scaled by temperature on the ring bc_type21_ref away

22100   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'fixed electron temperature (modified by nearby ring) ',
     1     ENEPAR(IB,1),' with ',ENEPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(enepar(ib,1).ne.0.0_R8,
     1     'faulty parameter enepar1')
          call xertst(enepar(ib,2).gt.0.0_R8,
     1     'faulty parameter enepar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          teave = 0.0_R8
          neave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            teave=teave+te(bc_type21_ix,bc_type21_iy)*
     1       ne(bc_type21_ix,bc_type21_iy)*
     2       vol(bc_type21_ix,bc_type21_iy)
            neave=neave+ne(bc_type21_ix,bc_type21_iy)*
     1       vol(bc_type21_ix,bc_type21_iy)
          ENDDO
          teave=teave/neave
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*ENEPAR(IB,1)*EV*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .       ((1.0_R8-enepar(ib,2))+
     .           enepar(ib,2)*te(bc_type21_ix,bc_type21_iy)/teave)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE IF(BCCHAR(IB).EQ.'N') THEN
          teave = 0.0_R8
          neave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            teave=teave+te(bc_type21_ix,bc_type21_iy)*
     1       ne(bc_type21_ix,bc_type21_iy)*
     2       vol(bc_type21_ix,bc_type21_iy)
            neave=neave+ne(bc_type21_ix,bc_type21_iy)*
     1       vol(bc_type21_ix,bc_type21_iy)
          ENDDO
          teave=teave/neave
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*ENEPAR(IB,1)*EV*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .       ((1.0_R8-enepar(ib,2))+
     .           enepar(ib,2)*te(bc_type21_ix,bc_type21_iy)/teave)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=21, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 30000

* -- BCENE=22 -- THE MODIFIED DECAY LENGTH FOR THE ELECTRON TEMPERATURE   !srv 23.07.10 {

22200   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     ' BCENE=22 with ENEPAR(IB,1)=',ENEPAR(IB,1),
     1     ' on',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            vte = sqrt(te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/me)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=VBIG*NE00*EV*
     &           CHCE(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)*
     &             TE(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=- VBIG*NE00*EV*
     &           (CHCE(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)-
     &            ENEPAR(IB,1)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &            VTE*NE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            vte = sqrt(te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/me)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=VBIG*NE00*EV*
     &           CHCE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)*
     &             TE(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=- VBIG*NE00*EV*
     &       (CHCE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-ENEPAR(IB,1)*
     &          SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &             BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     &               VTE*NE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=22, NOT AVAILABLE FOR W,E')
        END IF                                                          !srv 23.07.10 }
        GOTO 30000

* -- BOUNDARY CONDITIONS FOR ION ENERGY EQUATIONS --

30000   continue
        if(BCENI(IB).eq.0) goto 40000
        GOTO(30100,30200,30300,30400,30500,
     1     30600,30700,30800,30900,31000,
     2     31100,31200,31300,31400,31500,
     3     31600,31700,31800,31900,32000,
     4     32100,32200,32300,32400), BCENI(IB)                         !srv 01.02.09 19.04.13
        write(*,*) BCENI(IB)
        call xerrab ('BCENI out of bounds')

* -- BCENI=1 -- PRESCRIBE THE VALUE OF THE ION TEMPERATURE

30100   CONTINUE
        if(ncall.eq.0) then
          IF(BCCHAR(IB).EQ.'S'.and.use_astra.ne.0) THEN
            ENIPAR(IB,1) = tib_astra*1.0e3_R8                           !srv 14.06.12
          ENDIF
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified ion temperature ',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(ENIPAR(IB,1).gt.0.0_R8,
     1     'BCENI = 1, ENIPAR(,1) !> 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*ENIPAR(IB,1)*EV*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*ENIPAR(IB,1)*EV*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       UBIG*NE00*ENIPAR(IB,1)*EV*
     .       SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .      -UBIG*NE00*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       UBIG*NE00*ENIPAR(IB,1)*EV*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-UBIG*NE00*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'X') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= UBIG*NE00*
     .       ENIPAR(IB,1)*EV*VOL(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-UBIG*NE00*
     .       VOL(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'Y') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= UBIG*NE00*
     .       ENIPAR(IB,1)*EV*VOL(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-UBIG*NE00*
     .       VOL(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCENI=1, INVALID OPTION.')
        END IF
        GOTO 40000

* -- BCENI=2 -- PRESCRIBE THE GRADIENT OF THE ION TEMPERATURE

30200   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'ion temperature gradient ',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        VBIG*NE00*EV*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .       (-ENIPAR(IB,1)*0.5_R8*
     .         (HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     .          HY1(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))) +
     .       TI(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -VBIG*NE00*EV*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*EV*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .       (ENIPAR(IB,1)*0.5_R8*
     .         (HY1(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     .          HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))) +
     .       TI(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*EV*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .        UBIG*NE00*EV*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     .       (-ENIPAR(IB,1)*0.5_R8*
     .         (HX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))+
     .          HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .             RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))) +
     .       TI(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .       -UBIG*NE00*EV*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= UBIG*NE00*EV*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))*
     .       (ENIPAR(IB,1)*0.5_R8*
     .         (HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .             LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))+
     .          HX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))) +
     .       TI(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-UBIG*NE00*EV*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          end do
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCENI=2, INVALID OPTION.')
        END IF
        GOTO 40000

* -- BCENI=3 -- SHEATH CONDITIONS, ION ENERGY TRANSMISSION COEFFICIENT

30300   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,g14.7,a,g14.7,a,a,a)')
     .     'sheath ion with MOMPAR(ismain,IB,2) =',
     1     MOMPAR(ismain,IB,2),' and MOMPAR(ismain,IB,1) =',
     2     MOMPAR(ismain,IB,1),
     .     ' and enipar(ib,1) =',enipar(ib,1),               !srv 24.02.10 {
     1     ' and enipar(ib,2) =',enipar(ib,2),
     2     ' on ',BCCHAR(IB),boundary_location(ib)           !srv 24.02.10 }
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
              CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))  !srv 01.02.09
              UUOUT=MOMPAR(ismain,IB,1)*CS*
     .         abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .             wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3))
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  if(enipar(ib,1).gt.0.0_R8) then                        !xpb
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) - UUOUT *
     .               ENIPAR(IB,1)*bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .               SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
                  else     ! xpb new pressure linearisation as per B2.5-style
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)=
     .              SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3) - UUOUT *  !xpb
     .               (-ENIPAR(IB,1))*bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .                SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     .                ni(bc_list_x(iy,ib),bc_list_y(iy,ib),0)
                  endif                                                  !xpb
                endif
              end do
            ELSE
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  ZBAR=bcint(na(-1,-1,is),nx,ny,IY,IB)/
     .             bcint(NE,nx,ny,IY,IB)
                  ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .             bcint(rz2(-1,-1,is),nx,ny,IY,IB)*ZBAR*
     .             bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                  CS=SQRT(ZARG)
                  UUTEST=
     .             UU(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .             (1.0_R8+
     .              HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .      HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .         RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .      UU(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .         RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .              HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .      HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .         RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
                  UUOUT=MAX(MOMPAR(IS,IB,1)*CS*
     .             abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .                 wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .             -UUTEST)
                  if(enipar(ib,1).gt.0.0_R8) then                        !xpb
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) - UUOUT *
     .               ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .               SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
                  else     ! xpb new pressure linearisation as per B2.5-style
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)=
     .              SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3) - UUOUT *  !xpb
     .               (-ENIPAR(IB,1))*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .                SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     .                ni(bc_list_x(iy,ib),bc_list_y(iy,ib),0)
                  endif                                                  !xpb
                endif
              end do
            ENDIF
            if(mdf_fhi.ne.0) then                              !srv 12.09.08 {
              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &         shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &         fhiPSch(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                 rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)
            endif                                              !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
              CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))  !srv 01.02.09
              UUOUT=MOMPAR(ismain,IB,1)*CS*
     .         abs(wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .             wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),3))
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  if(enipar(ib,1).gt.0.0_R8) then                        !xpb
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) - UUOUT *
     .               ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .               SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                  LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
                  else     ! xpb new pressure linearisation as per B2.5-style
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)=
     .              SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3) - UUOUT *  !xpb
     .               (-ENIPAR(IB,1))*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .                SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .                ni(bc_list_x(iy,ib),bc_list_y(iy,ib),0)
                  endif                                                  !xpb
                endif
              end do
            ELSE
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  ZBAR=bcint(na(-1,-1,IS),nx,ny,IY,IB)/
     .             bcint(NE,nx,ny,IY,IB)
                  ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .             bcint(rz2(-1,-1,IS),nx,ny,IY,IB)*ZBAR*
     .             bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                  CS=SQRT(ZARG)
                  UUTEST=
     .             UU(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     1             (1.0_R8+
     2              HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .       HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .       UU(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .              HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .       HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
                  UUOUT=MAX(MOMPAR(IS,IB,1)*CS*
     .             abs(wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .                 wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .             UUTEST)
                  if(enipar(ib,1).gt.0.0_R8) then                        !xpb
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) - UUOUT *
     .               ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .               SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
                  else     ! xpb new pressure linearisation as per B2.5-style
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)=
     .              SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3) - UUOUT *  !xpb
     .               (-ENIPAR(IB,1))*
     .                 bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .                 SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                    LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .                ni(bc_list_x(iy,ib),bc_list_y(iy,ib),0)
                  endif                                                  !xpb
                endif
              end do
            ENDIF
            if(mdf_fhi.ne.0) then                              !srv 12.09.08 {
              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &         shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) +
     &         fhiPSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)
            endif                                              !srv 12.09.08 }
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=3, NO SHEATH IN RADIAL DIR.')
        END IF
        GOTO 40000

* -- BCENI=4 -- PRESCRIBE THE VALUE OF THE ION TEMPERATURE, WEAKLY
*                 A MIXED BOUNDARY CONDITION

30400   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'weak ion temperature ',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(ENIPAR(IB,1).gt.0.0_R8,
     1     'BCENI = 4, ENIPAR(,1) !> 0')
          call xertst(ENIPAR(IB,2).gt.0.0_R8,                !srv 02.02.11
     1     'BCENI = 4, ENIPAR(,2) !> 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= ENIPAR(IB,2)
     .       *ENIPAR(IB,1)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-ENIPAR(IB,2)
     .       /EV*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= ENIPAR(IB,2)*
     .       ENIPAR(IB,1)*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-ENIPAR(IB,2)/EV*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= ENIPAR(IB,2)*
     1       ENIPAR(IB,1)*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-ENIPAR(IB,2)/EV*
     .       SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= ENIPAR(IB,2)*
     .       ENIPAR(IB,1)*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-ENIPAR(IB,2)/EV*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=4, NOT YET AVAILABLE')
        END IF
        GOTO 40000

* -- BCENI=5 -- PRESCRIBE THE ION ENERGY FLUX PER UNIT AREA

30500   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'ion energy flux per area ',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENIPAR(IB,1)*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENIPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=ENIPAR(IB,1)*
     .        SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=ENIPAR(IB,1)*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=5, NOT YET AVAILABLE')
        END IF
        GOTO 40000

* -- BCENI=6 -- PRESCRIBE THE TOTAL ION ENERGY FLUX FOR A
*                 CONSTANT ION TEMPERATURE

cghgh
30600   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'total ion energy flux ',ENIPAR(IB,1),
     1     ' with const Ti on ',BCCHAR(IB),boundary_location(ib)
        endif
#ifdef WANT_THIS
        IF(BCCHAR(IB).EQ.'N') THEN
          Tguess=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            Tguess=Tguess+TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
          Tguess=Tguess/BC_LIST_SIZE(IB)
          if(diagno.ge.2) write(*,*) 'Initial Tguess:',Tguess/ev
          T0=Tguess*0.95_R8
          T1=Tguess*1.01_R8
          Tlimit=Tguess*0.0001_R8
          fylimt=ENIPAR(IB,1)*0.0001_R8
          ibnd=-1
30610     ibnd=ibnd+1
          if (ibnd.eq.0) then
            Tguess=T0
          elseif (ibnd.eq.1) then
            Tguess=T1
          else
            Tguess=T1-GAM1*(T1-T0)/(EPSLON + GAM1 - GAM0)
          endif
          fyt=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
cxpb code segment analog to b2tfhi
            if(ncall.eq.0) then
             do is = 0, ns-1
              if(.not.is_neutral(is)) then
               vadia(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is) =
     .          ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))*
     .          bb(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2)
     &          /(qe*
     &            rza(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),is)*
     *            hx(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     &          *(2.0_R8/
     &          bb(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)**3)
     &          *(wbbr(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)
     &           -wbbl(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3))
     &          *facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))
               vaecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is) =
     .          1.0_R8/hx(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     1          bb(bc_list_x(ix,ib),bc_list_y(ix,ib),2)/
     .          bb(bc_list_x(ix,ib),bc_list_y(ix,ib),3)**2*0.5_R8*
     2         (po(rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             rightiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     2          po( leftix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              leftiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))*
     &         fac_ExB(bc_list_x(ix,ib),bc_list_y(ix,ib))
              else
               vadia(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is)=0.0_R8
               vaecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is)=0.0_R8
              endif
             enddo
            endif
            if (conduction_only.eq.1) then
              flo = 0.0_R8
            else
              cum_flo_vadia=0.0_R8
              cum_flo_vaecrb=0.0_R8
              do is=0,ns-1
                prefactor=
     &           (vol(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &            vol(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &           /(hy1(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &             hy1(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &           *(na(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),is)
     &            +na(bc_list_x(ix,ib),bc_list_y(ix,ib),is))/2.0_R8
                cum_flo_vadia=cum_flo_vadia+
     1           prefactor*vadia(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is)
                cum_flo_vaecrb=cum_flo_vaecrb+prefactor*
     1           vaecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is)
              enddo
              flo = (1.5_R8+BoRiS)*
     .         fni_32(bc_list_x(ix,ib),bc_list_y(ix,ib),1)+
     1         2.5_R8*fni_52(bc_list_x(ix,ib),bc_list_y(ix,ib),1)+
     &         chvi(bc_list_x(ix,ib),bc_list_y(ix,ib),1)*
     &         (ni(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),0)+
     &          ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0))/2.0_R8+
     1         2.5_R8*cum_flo_vadia+1.5_R8*cum_flo_vaecrb
            endif
            if (no_hybr_i.eq.0) then
              coni0 = hybr(flo,
     .         chci(bc_list_x(ix,ib),bc_list_y(ix,ib),1))
              fy = flo*
     .         (ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &         +Tguess)/2.0_R8-coni0*
     1         (Tguess-ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     2                    bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
            else
              coni0 = chci(bc_list_x(ix,ib),bc_list_y(ix,ib),1)
              t2 = dim ( coni0,
     .         abs(flo) *
     .         max(0.0_R8,
     .             min(Tguess,ti(bottomix(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)),
     .                           bottomiy(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)))))/
     .         max(cutlo,
     .                 Tguess+ti(bottomix(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)),
     .                           bottomiy(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)))) )
              fy =
     .         upwind(flo,
     .                ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     .                Tguess)
     .         - t2 * (Tguess-ti(bottomix(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)),
     .                           bottomiy(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib))))
            endif
            fyt=fyt + fy
          ENDDO
          if(diagno.ge.2)
     >     write(*,'(a,i4,1pe15.6,1pe14.6)') 'ibnd,Tg,fyt(i):',
     .     ibnd,Tguess/ev,fyt
          if (ibnd.eq.0) then
            GAM0 = fyt - ENIPAR(IB,1)
            goto 30610
          elseif (ibnd.eq.1) then
            GAM1 = fyt - ENIPAR(IB,1)
            goto 30610
          else
            T0 = T1
            GAM0 = GAM1 - ENIPAR(IB,1)
            T1 = Tguess
            GAM1 = fyt
            if (ibnd .gt. 100)
     >       call xerrab ('B2STBC -- BCENI=6, NOT CONVERGING')
            if (abs(t1-t0).gt.Tlimit .and.
     &          abs(fyt-ENIPAR(IB,1)).gt.fylimt) goto 30610
          endif
          if(diagno.ge.1)
     >     write(*,'(a,f15.6,i4,1pe14.6)') ' tibound:',T1/ev,ibnd,fyt
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*T1*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'S') THEN
          Tguess=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            Tguess=Tguess+TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
          Tguess=Tguess/BC_LIST_SIZE(IB)
          if(diagno.ge.2) write(*,*) 'Initial Tguess:',Tguess/ev
          T0=Tguess*0.95_R8
          T1=Tguess*1.01_R8
          Tlimit=Tguess*0.0001_R8
          fylimt=ENIPAR(IB,1)*0.0001_R8
          ibnd=-1
30620     ibnd=ibnd+1
          if (ibnd.eq.0) then
            Tguess=T0
          elseif (ibnd.eq.1) then
            Tguess=T1
          else
            Tguess=T1-GAM1*(T1-T0)/(EPSLON + GAM1 - GAM0)
          endif
          fyt=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
cxpb code segment analog to b2tfhi
            if(ncall.eq.0) then
             do is = 0, ns-1
              if(.not.is_neutral(is)) then
               vadia(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,is) =
     &          ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))*
     &          bb(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2)/
     &          (qe*rza(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                  topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),is)*
     *           hx(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     &          *(2.0_R8/
     &          bb(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)**3)
     &          *(wbbr(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)-
     &            wbbl(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3))
     &          *facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))
               vaecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is) =
     .          1.0_R8/hx(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     1          bb(bc_list_x(ix,ib),bc_list_y(ix,ib),2)/
     .          bb(bc_list_x(ix,ib),bc_list_y(ix,ib),3)**2*0.5_R8*
     2         (po(rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             rightiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     2          po( leftix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              leftiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))*
     &         fac_ExB(bc_list_x(ix,ib),bc_list_y(ix,ib))
              else
               vadia(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,is)=
     .                                                        0.0_R8
               vaecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is)=0.0_R8
              endif
             enddo
            endif
            if (conduction_only.eq.1) then
              flo = 0.0_R8
            else
             cum_flo_vadia = 0.0_R8
             cum_flo_vaecrb = 0.0_R8
             do is=0,ns-1
               prefactor=
     &          (vol(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &           vol(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &          /(hy1(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &            hy1(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &          *(na(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),is)
     &           +na(bc_list_x(ix,ib),bc_list_y(ix,ib),is))/2.0_R8
               cum_flo_vadia=cum_flo_vadia+prefactor*
     1          vadia(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     2                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,is)
               cum_flo_vaecrb=cum_flo_vaecrb+
     1          prefactor*vaecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is)
             enddo
             flo = (1.5_R8+BoRiS)*
     1         fni_32(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)+
     1         2.5_R8*fni_52(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)+
     &         chvi(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)*
     &         (ni(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),0)+
     &          ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0))/2.0_R8+
     1         2.5_R8*cum_flo_vadia+1.5_R8*cum_flo_vaecrb
            endif
            if (no_hybr_i.eq.0) then
              coni0 = hybr(flo,
     &         chci(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1))
              fy = flo*(ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                     topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &                 +Tguess)/2.0_R8-coni0*
     1         (ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-Tguess)
            else
              coni0 = chci(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                     topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)
              t2 = dim ( coni0,
     .         abs(flo) *
     .         max(0.0_R8,
     .             min(Tguess,ti(topix(bc_list_x(ix,ib),
     .                                 bc_list_y(ix,ib)),
     .                           topiy(bc_list_x(ix,ib),
     .                                 bc_list_y(ix,ib)))))/
     .         max(cutlo,
     .                 Tguess+ti(topix(bc_list_x(ix,ib),
     .                                 bc_list_y(ix,ib)),
     .                           topiy(bc_list_x(ix,ib),
     .                                 bc_list_y(ix,ib)))) )
              fy =
     .         upwind(flo, Tguess,
     .                ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .         - t2 * (ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                    topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     .                 Tguess)
            endif
            fyt=fyt + fy
          ENDDO
          if(diagno.ge.2)
     >     write(*,'(a,i4,f15.6,1pe14.6)') 'ibnd,Tg,fyt(i):',
     .     ibnd,Tguess/ev,fyt
          if (ibnd.eq.0) then
            GAM0 = fyt - ENIPAR(IB,1)
            goto 30620
          elseif (ibnd.eq.1) then
            GAM1 = fyt - ENIPAR(IB,1)
            goto 30620
          else
            T0 = T1
            GAM0 = GAM1 - ENIPAR(IB,1)
            T1 = Tguess
            GAM1 = fyt
            if (ibnd .gt. 100)
     >       call xerrab ('B2STBC -- BCENI=6, NOT CONVERGING')
            if (abs(t1-t0).gt.Tlimit .and.
     &          abs(fyt-ENIPAR(ib,1)).gt.fylimt) goto 30620
          endif
          if(diagno.ge.1)
     >     write(*,'(a,f15.6,i4,1pe14.6)') ' tibound:',T1/ev,ibnd,fyt

          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .          VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .         -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=6, NOT AVAILABLE FOR W,E')
        ENDIF
#else
        CALL XERRAB ('B2STBC -- BCENI=6, NOT IMPLEMENTED.')
#endif
        GOTO 40000
cghgh

* -- BCENI=7 -- PRESCRIBE THE ION TEMPERATURE AS A FUNCTION OF OTHER
*                 PLASMA PARAMETERS

30700   CONTINUE
        CALL XERRAB ('B2STBC -- BCENI=7, NOT YET AVAILABLE')
        GOTO 40000

* -- BCENI=8 -- PRESCRIBE THE TOTAL ION HEAT FLUX WITH CONSTANT
*                 FLUX DENSITY

30800   CONTINUE
        if(ncall.eq.0) then
          if (use_astra.ne.0) ENIPAR(IB,1) = fhib_astra           !srv 28.07.11
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'total ion energy flux',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          US=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            US=US+SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       ENIPAR(IB,1)/US*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          if(mdf_fhi.ne.0) then                                !srv 12.09.08 {
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &         fhiPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)
              do is=0,ns-1
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     &           1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                  hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &           vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)*
     &           na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
              end do
          endif                                                !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          US=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            US=US+SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .               BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENIPAR(IB,1)/US*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          if(mdf_fhi.ne.0) then                                !srv 12.09.08 {
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     .         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &         fhiPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)
              do is=0,ns-1
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     .           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) -
     &           1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                  hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &           vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)*
     &               na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
              end do
          endif                                                !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          US=0.0_R8
          DO IY = 1, BC_LIST_SIZE(IB)
            US=US+SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .         ENIPAR(IB,1)/US*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          US=0.0_R8
          DO IY = 1, BC_LIST_SIZE(IB)
            US=US+SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          end do
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=ENIPAR(IB,1)/US*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=8, NOT AVAILABLE')
        ENDIF
        GOTO 40000

* -- BCENI=9 -- PRESCRIBE THE DECAY LENGTH FOR THE ION TEMPERATURE

30900   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'ion temperature decay length ',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(ENIPAR(IB,1).ne.0.0_R8,
     1     'BCENI = 9, ENIPAR(,1) = 0')
          decay_length_ok = .true.
          do is = 0, ns-1
            decay_length_ok = decay_length_ok .and.
     1       (bccon(is,ib).eq.9 .or.
     2       (is_neutral(is) .and. use_eirene.eq.1))
          end do
          if (.not.decay_length_ok) then
!          call xertst(all(BCCON(0:NS-1,IB).EQ.9),
!     1     'Should use BCCON(IS,IB).EQ.9 together with BCENI=9 ')       !srv 16.07.12
            write(0,*)
     &       'Should use BCENI.EQ.9 together with BCCON=9'
            write(*,*)
     &       'Should use BCENI.EQ.9 together with BCCON=9'
            write(*,*) 'Please ensure that this boundary condition '//
     &                 'is consistent with density boundary condition'  !srv 21.09.12
          endif
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .         -hci0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/ENIPAR(IB,1)*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            DO IS=0,NS-1
cxpb          if (.not.is_neutral(is)) then   ! dpc --- need to check
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .           2.5_R8*TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .           SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS) !srv 23.09.08
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     .           2.5_R8*
     .           SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*   !srv 23.09.08
     &                    na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
cxpb          endif
            ENDDO
            if(mdf_fhi.ne.0) then                                !srv 12.09.08 {
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &        shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &        fhiPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)
              do is=0,ns-1
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     &           1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                  hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &           vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)*
     &                     na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
              end do
            endif                                                !srv 12.09.08 }
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .         -hci0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/ENIPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            DO IS=0,NS-1
cxpb          if (.not.is_neutral(is)) then   ! dpc --- need to check
cxpb            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
cxpb     .       shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
cxpb     .           rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*2.5_R8*
cxpb     .            TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
cxpb     .           SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)
cxpb            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
cxpb     .       shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
cxpb     .           rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*2.5_R8*
cxpb     .       SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*
cxpb     .           na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
cxpb          endif
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .           2.5_R8*TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS) !srv 23.09.08
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     .           2.5_R8*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*     !srv 23.09.08
     &                  na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            ENDDO
            if(mdf_fhi.ne.0) then                                !srv 12.09.08 {
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     .         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &         fhiPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)
              do is=0,ns-1
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) -
     &           1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                  hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &           vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)*
     &               na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
             end do
            endif                                                !srv 12.09.08 }
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=9, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 40000

* -- BCENI=10 -- Feedback option for core

31000   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'ion heat flux feedback bc on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(lfeedback, 'LFEEDBACK must be turned on!')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if (region(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0).eq.1) then
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     1         cbshi(0,ismain,coreregno)*
     .         SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            elseif (region(bc_list_x(ix,ib),bc_list_y(ix,ib),0).eq.5)
     .       then
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     1         cbshi(0,ismain,coreregn2)*
     .         SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            endif
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=10, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 40000

* -- BCENI=11 -- SHEATH CONDITIONS, ION ENERGY TRANSMISSION COEFFICIENT

31100   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,a,a)')
     .     'sheath (11) ion on ',BCCHAR(IB),
     1     boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            DO IS=0,NS-1
              if (.not.is_neutral(is)) then
                ZBAR=bcint(na(-1,-1,IS),nx,ny,IY,IB)/
     .           bcint(NE,nx,ny,IY,IB)
                ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .           bcint(rz2(-1,-1,IS),nx,ny,IY,IB)*ZBAR*
     .           bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                CS=SQRT(ZARG)
                UUTEST=
     .           UU(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .           (1.0_R8+
     .            HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .    HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .       RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .    UU(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .       RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .            HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .    HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .       RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
                UUOUT=MAX(MOMPAR(IS,IB,1)*CS*
     .           abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .               wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .           -UUTEST)
                shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .           shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) - UUOUT *
     .           ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .           SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              endif
            end do
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            DO IS=0,NS-1
              if (.not.is_neutral(is)) then
                ZBAR=bcint(na(-1,-1,IS),nx,ny,IY,IB)/
     .           bcint(NE,nx,ny,IY,IB)
                ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .           bcint(rz2(-1,-1,IS),nx,ny,IY,IB)*ZBAR*
     .           bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                CS=SQRT(ZARG)
                UUTEST=UU(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                    LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     1        (1.0_R8+HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .         HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .            LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .         UU(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .            LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .                HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .         HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .            LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
                UUOUT=MAX(MOMPAR(IS,IB,1)*CS*
     .           abs(wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .               wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .           UUTEST)
                shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .           shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) - UUOUT *
     .           ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .           SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              endif
            end do
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=11, NO SHEATH IN RADIAL DIR.')
        END IF
        GOTO 40000

* -- BCENI=12 -- SHEATH CONDITIONS, ION ENERGY TRANSMISSION COEFFICIENT

31200   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,a,a)')
     .     'sheath (12) ion on ',BCCHAR(IB),
     1     boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            DO IS=0,NS-1
              shi0(bc_list_x(ix,ib),bc_list_y(ix,ib),3) =
     .        shi0(bc_list_x(ix,ib),bc_list_y(ix,ib),3)-
     1         ENIPAR(IB,1)*max(0.0_R8,
     1          -fna(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,is))/
     2         ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)
            end do
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            DO IS=0,NS-1
              shi0(bc_list_x(ix,ib),bc_list_y(ix,ib),3) =
     .        shi0(bc_list_x(ix,ib),bc_list_y(ix,ib),3)-
     1         ENIPAR(IB,1)*max(0.0_R8,
     .         fna(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is))/
     2         ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)
            end do
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            DO IS=0,NS-1
              shi0(bc_list_x(iy,ib),bc_list_y(iy,ib),3) =
     .        shi0(bc_list_x(iy,ib),bc_list_y(iy,ib),3)-
     1         ENIPAR(IB,1)*max(0.0_R8,
     1          -fna(rightix(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .               rightiy(bc_list_x(iy,ib),bc_list_y(iy,ib)),0,is))/
     2         ni(bc_list_x(iy,ib),bc_list_y(iy,ib),0)
            end do
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            DO IS=0,NS-1
              shi0(bc_list_x(iy,ib),bc_list_y(iy,ib),3) =
     .        shi0(bc_list_x(iy,ib),bc_list_y(iy,ib),3)-
     1         ENIPAR(IB,1)*max(0.0_R8,
     .         fna(bc_list_x(iy,ib),bc_list_y(iy,ib),0,is))/
     2         ni(bc_list_x(iy,ib),bc_list_y(iy,ib),0)
            end do
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=12, INVALID OPTION')
        END IF
        GOTO 40000

* -- BCENI=13 -- PRESCRIBE THE ION ENERGY FLUX PER UNIT AREA
*             -- PROPORTIONAL TO TEMPERATURE

31300   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     1     'ion energy flux per area proportional to temperature ',
     2     ENIPAR(IB,1),' on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENIPAR(IB,1)*
     .       TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENIPAR(IB,1)*
     .       TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=ENIPAR(IB,1)*
     .       TI(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     1       SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=ENIPAR(IB,1)*
     .       TI(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=13, INVALID OPTION')
        END IF
        GOTO 40000

* -- BCENI=14 -- leakage option for ion energy

31400   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'ion energy leakage ',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            CS=SQRT(bcint(PR,nx,ny,IX,IB)/bcint(RO,nx,ny,IX,IB))  !srv 01.02.09
            SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     1       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))/
     .       hy1(bc_list_x(ix,ib),bc_list_y(ix,ib))*cs*ENIPAR(IB,1)
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            CS=SQRT(bcint(PR,nx,ny,IX,IB)/bcint(RO,nx,ny,IX,IB))  !srv 01.02.09
            SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     1       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))/
     .       hy1(bc_list_x(ix,ib),bc_list_y(ix,ib))*cs*ENIPAR(IB,1)
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=14, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 40000

* -- BCENI=15 -- from b2stbc_spb                              !srv 01.02.09 {

31500   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     &     'from b2stbc_spb for BCENI with enipar1 =',
     &     ENIPAR(IB,1),' and enipar2 =', ENIPAR(IB,2),' on ',
     &     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
!            SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=0.0_R8
!            CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))
            CS=SQRT(PR(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &              RO(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            t0 = 0.0_R8
            do is = 0, ns-1
              if(.not.is_neutral(is)) t0 = t0 +
     &         na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)       !srv 11.03.09
            enddo
            SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &       SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) +                !srv 11.03.09
     &       ENIPAR(IB,1)*t0*cs*
     &       abs(pbs(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &               rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)) !srv 26.11.02 28.01.03
            if(mdf_fhi.ne.0) then                                        !srv 13.10.06 {
              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &         shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &         fhiPSch(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                 rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)
            endif                                                         !srv 13.10.06 }
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
!            SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=0.0_R8
!            CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))
            CS=SQRT(PR(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     .              RO(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            t0 = 0.0_R8
            do is = 0, ns-1
              if(.not.is_neutral(is)) t0 = t0 +
     .         na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)   !srv 11.03.09
            enddo
            SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     .       SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) +            !srv 11.03.09
     &       ENIPAR(IB,1)*t0*cs*
     &       abs(pbs(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))         !srv 26.11.02 28.01.03
            if(mdf_fhi.ne.0) then                                    !srv 13.10.06 {
              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     .         shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) +
     &         fhiPSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)
            endif                                                    !srv 13.10.06 }
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=15, NOT AVAILABLE FOR S,N')
        ENDIF
        GOTO 40000                                                   !srv 01.02.09 }

* -- BCENI=16 -- Feedback boundary condition with constant temperature

31600   CONTINUE
        if(bceni_16_style.eq.0) then
         if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .       'ion temperature feedback by flux',
     1       ENIPAR(IB,1),' with ',ENIPAR(IB,2),' on ',
     1       BCCHAR(IB),boundary_location(ib)
           call xertst(enipar(ib,1).ne.0, 'faulty parameter enipar1')
           call xertst(enipar(ib,2).gt.0, 'faulty parameter enipar2')
         endif
         IF(BCCHAR(IB).EQ.'S') THEN
           totflux = 0.0_R8
           t0 = 0.0_R8
           weight = 0.0_R8
           DO IX = 1, BC_LIST_SIZE(IB)
             bc_type16_ix=bc_list_x(ix,ib)
             bc_type16_iy=bc_list_y(ix,ib)
             do bc_type16_i=1,bc_type16_ref      ! how many steps away from the boundary?
               bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
               bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
               bc_type16_ix=bc_type16_ixx
               call xertst (-2.lt.bc_type16_ix.and.bc_type16_ix.lt.nx+1,
     1          'stepped off x-grid in type 16 bc')
               call xertst (-2.lt.bc_type16_iy.and.bc_type16_iy.lt.ny+1,
     1          'stepped off y-grid in type 16 bc')
             enddo
             totflux = totflux +
     .        fhi(bc_type16_ix,bc_type16_iy,1)
             if (bc_type16_kinetic_energy.eq.1) then
               do is = 0, ns-1
                 totflux=totflux+(1.0_R8-BoRiS)*
     .            fna(bc_type16_ix,bc_type16_iy,1,is)*0.5_R8*
     .            (kinrgy(bc_type16_ix,bc_type16_iy,is)+
     .             kinrgy(bottomix(bc_type16_ix,bc_type16_iy),
     .                    bottomiy(bc_type16_ix,bc_type16_iy),is))
               enddo
             endif
             t0 = t0 +
     .         ti(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .        vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         ni(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)
             weight = weight +
     .        vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         ni(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)
           ENDDO
           t0 = t0 / weight
           if(diagno.ge.2) write(*,*) 'Totflux :',totflux
           T1 = (1.0_R8 -
     .      enipar(ib,2)*(totflux-enipar(ib,1))/enipar(ib,1)) * t0
           if(diagno.ge.2)
     1      write(*,*) 'Old/New ion temp :',t0/qe, t1/qe
           DO IX = 1, BC_LIST_SIZE(IB)
             shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
             shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
           ENDDO
         ELSE
           CALL XERRAB ('B2STBC -- BCENI=16, NOT AVAILABLE FOR N,W,E')
         ENDIF
         GOTO 40000

        elseif(bceni_16_style.eq.1)then                                 !srv 16.02.12   {

          if(ncall.eq.0) then
          write(*,*) 'Given integrated ion heat flux',
     &                 ' with constant Ti '
!          write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
            if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
            write(*,'(a,a,a)') ' on ',BCCHAR(IBW),boundary_location(ibw)
                call xertst(
     .         mod(region(BC_LIST_X(1,IBW),BC_LIST_Y(1,IBW),0),4).eq.1,
     .         'BCENI = 16 only applicable to core boundaries !')
                call xertst(enipar(ibw,1).ne.0,
     &           'faulty parameter enipar1')
              endif
            if(BCCHAR(IBW).EQ.'S'.and.
     &         (BCENI(IBW).EQ.23.or.BCENI(IBW).EQ.24)) then
                call xerrab
     .           ('BCENI types 16, 23 and 24 cannot be mixed!')
            endif
          ENDDO
            if (use_astra.ne.0) then
#ifdef ASTRA
              SENIPAR = fhib_astra
#else
              call xerrab
     &         ('Compile with -DASTRA option to couple to ASTRA!')
#endif
            else
              SENIPAR = 0.0_R8
              DO IBW=IB,NBC
                if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
                  SENIPAR = SENIPAR + ENIPAR(IBW,1)
                endif
              ENDDO
            endif
            write(*,*) ' Summed value is ', SENIPAR
          endif
          if(.not.once_bceni) then
            TIS = 0.0_R8
            NIVS = 0.0_R8
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  TIS = TIS + TI(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     .                        NI(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     .                       VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  NIVS = NIVS +
     .                   NI(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     .                  VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                END DO
              endif
            END DO
            TIS = TIS/NIVS
            IF(BCCHAR(IB).EQ.'S') THEN
            totflux = 0.0_R8
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCENI(IBW).EQ.BCENI(IB)) THEN
                if(mdf_fhi.ne.0) then
                    DO IX = 1, BC_LIST_SIZE(IBW)
                totflux = totflux +
     &          fhi_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
     &       +fac_ExB(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*        !srv 11.04.13 {
     &      ((po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))-
     & po(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &    leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))/
     &       (hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     & hx(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &    leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))+
     &       (po(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &           rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))-
     &        po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))/
     &       (hx(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &           rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))+
     &        hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))*
     &        ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)*
     &        ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &        SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &        bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),2)/
     &        bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),3)**2                !srv 11.04.13 }
                    ENDDO
                else
                    DO IX = 1, BC_LIST_SIZE(IBW)
                    totflux = totflux +
     &             fhi(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
                    ENDDO
                endif
              ENDIF
            ENDDO
              TIS = TIS *
     &         (1.0_R8-bc_ref_ti*(totflux-SENIPAR)/                      !srv 11.04.13
     &           max(abs(totflux),SENIPAR))
              write(*,"(' BCENI=16: TIS=',es14.7,' totflux=',es14.7,
     *         ' SENIPAR=',es14.7)") TIS/ev, totflux, SENIPAR
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCENE(IBW).EQ.BCENE(IB)) THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)=
     .               VBIG*NE00*TIS*
     .               SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)=
     .              -VBIG*NE00*
     .               SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  ENDDO
              ENDIF
            ENDDO
              once_bceni = .true.
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCENI=16, NOT AVAILABLE FOR N,W,E')
            ENDIF
        endif                                                          !srv 16.02.12 }
          GOTO 40000
      else
        CALL XERRAB (' Incorrect style of boundary condition')
        endif


* -- BCENI=17 -- Feedback boundary condition with constant temperature

31700   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'ion temperature feedback by flux ',
     1     ENIPAR(IB,1),' with ',ENIPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(enipar(ib,1).ne.0, 'faulty parameter enipar1')
          call xertst(enipar(ib,2).gt.0, 'faulty parameter enipar2')
          call xertst(bcene(ib).eq.bceni(ib),
     1     'both bcene and bceni should be set to 17')
          call xertst(enepar(ib,1).ne.0, 'faulty parameter enepar1')
          call xertst(enepar(ib,2).gt.0, 'faulty parameter enepar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
c already handled when the electron heat bc was handled
c see label 2980
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=16, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 40000

* -- BCENI=18 -- FRACTIONAL DROP CONDITION (FOR EDGE2D COMPARISON)

31800   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'Ion temperature fractional drop',
     1     ENIPAR(IB,1),' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENIPAR(IB,1).lt.1.0_R8,
     1     'BCENI = 18, ENIPAR(,1) => 1')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENIPAR(IB,1)) *
     .       TI(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .           VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .          -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENIPAR(IB,1)) *
     .       TI(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*T1*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENIPAR(IB,1)) *
     .       TI(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .           VBIG*NE00*T1*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .          -VBIG*NE00*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENIPAR(IB,1)) *
     .       TI(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= VBIG*NE00*T1*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-VBIG*NE00*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=18, INVALID OPTION')
        ENDIF
        GOTO 40000

* -- BCENI=19 -- PRESCRIBE THE DECAY LENGTH FOR THE ION TEMPERATURE

31900   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'ion temperature decay length (new)',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(ENIPAR(IB,1).ne.0.0_R8,
     1     'BCENI = 19, ENIPAR(,1) = 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .          -hci0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/ENIPAR(IB,1)*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            t0=0.0_R8
            DO IS=0,NS-1
               t0=t0-
     .              1.5_R8*fna_32(topix(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),
     .                            topiy(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),1,is)-
     .              2.5_R8*fna_52(topix(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),
     .                            topiy(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),1,is)
            ENDDO
            if(t0.gt.0.0_R8) then
               shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .          shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .              t0 * TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            else
               shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .          shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) + t0
            endif
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
             shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
             shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .          -hci0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/ENIPAR(IB,1)*
     1            SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .               BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
             t0=0.0_R8
             DO IS=0,NS-1
                t0=t0+
     .           1.5_R8*fna_32(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is)+
     .           2.5_R8*fna_52(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is)
             ENDDO
            if(t0.gt.0.0_R8) then
               shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .          shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .              t0 * TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
             else
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .          shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) + t0
             endif
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=19, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 40000

* -- BCENI=20 -- Feedback boundary condition with constant temperature scaled by temperature on the ring bc_type20_ref away

32000   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'ion temperature feedback by flux (neo) ',
     1     ENIPAR(IB,1),' with ',ENIPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(enipar(ib,1).ne.0, 'faulty parameter enipar1')
          call xertst(enipar(ib,2).gt.0, 'faulty parameter enipar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          totflux = 0.0_R8
          t0 = 0.0_R8
          weight = 0.0_R8
          tiave = 0.0_R8
          niave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type16_ix=bc_list_x(ix,ib)
            bc_type16_iy=bc_list_y(ix,ib)
            do bc_type16_i=1,bc_type16_ref      ! how many steps away from the boundary?
              bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
              bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
              bc_type16_ix=bc_type16_ixx
              call xertst (-2.lt.bc_type16_ix.and.bc_type16_ix.lt.nx+1,
     1         'stepped off x-grid in type 16 bc')
              call xertst (-2.lt.bc_type16_iy.and.bc_type16_iy.lt.ny+1,
     1         'stepped off y-grid in type 16 bc')
            enddo
            totflux = totflux +
     .       fhi(bc_type16_ix,bc_type16_iy,1)
            t0 = t0 +
     .       ti(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)
            weight = weight + vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)
            bc_type20_ix=bc_list_x(ix,ib)
            bc_type20_iy=bc_list_y(ix,ib)
            do bc_type20_i=1,bc_type20_ref      ! how many steps away from the boundary?
              bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
              bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
              bc_type20_ix=bc_type20_ixx
              call xertst (-2.lt.bc_type20_ix.and.bc_type20_ix.lt.nx+1,
     1         'stepped off x-grid in type 20 bc')
              call xertst (-2.lt.bc_type20_iy.and.bc_type20_iy.lt.ny+1,
     1         'stepped off y-grid in type 20 bc')
            enddo
            tiave=tiave+ti(bc_type20_ix,bc_type20_iy)*
     1       ni(bc_type20_ix,bc_type20_iy,0)*
     2       vol(bc_type20_ix,bc_type20_iy)
            niave=niave+ni(bc_type20_ix,bc_type20_iy,0)*
     1       vol(bc_type20_ix,bc_type20_iy)
          ENDDO
          t0 = t0 / weight
          tiave=tiave/niave
          if(diagno.ge.2) write(*,*) 'Totflux :',totflux
          T1 = (1.0_R8 -
     .     enipar(ib,2)*(totflux-enipar(ib,1))/enipar(ib,1)) * t0
          if(diagno.ge.2)
     1     write(*,*) 'Old/New ion temp :',t0/qe, t1/qe
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type20_ix=bc_list_x(ix,ib)
            bc_type20_iy=bc_list_y(ix,ib)
            do bc_type20_i=1,bc_type20_ref      ! how many steps away from the boundary?
              bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
              bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
              bc_type20_ix=bc_type20_ixx
              call xertst (-2.lt.bc_type20_ix.and.bc_type20_ix.lt.nx+1,
     1         'stepped off x-grid in type 20 bc')
              call xertst (-2.lt.bc_type20_iy.and.bc_type20_iy.lt.ny+1,
     1         'stepped off y-grid in type 20 bc')
            enddo
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*T1*
     1       ti(bc_type20_ix,bc_type20_iy)/tiave
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=20, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 40000

!! FIXME: Clash here with 5.0/5.1 BCENI=21 type

* -- BCENI=21 -- from b2stbc_spb                              !srv 01.02.09 {

32100   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     &     'from b2stbc_spb for BCENI with enipar1 =',
     &     ENIPAR(IB,1),' and enipar2 =', ENIPAR(IB,2),' on ',
     &     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            s0 = vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            if(mdf_fhi.ne.0) then                                    !srv 21.03.09 {
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &         fhiPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)
          endif                                                    !srv 21.03.09 }
            do is = 0, ns-1
            if(.not.is_neutral(is)) then                            !srv 19.03.09 {
                t0 = na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                csb = 
     .           sqrt((rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &                  te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &                  ti(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/
     &               (am(is)*mp))
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)+ENIPAR(IB,2)*
     &                            t0*csb*s0
!     &                            t0*csb*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              endif
            enddo
            SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &       SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) -
     &       max(0.0e0_R8,-floi_noc(
     &        topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &        topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)) !srv 28.10.02
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            s0 = vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            if(mdf_fhi.ne.0) then                                    !srv 21.03.09 {
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &         fhiPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)
          endif                                                    !srv 21.03.09 }
            do is = 0, ns-1
            if(.not.is_neutral(is)) then                            !srv 19.03.09 }
                t0 = na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                csb = sqrt((rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &                       te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &                       ti(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/
     &                     (am(is)*mp))
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)+ENIPAR(IB,2)*
     &          t0*csb*s0
!     &          t0*csb*SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
!     &                    BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
             endif                                                      !srv 19.03.09 {
            enddo
            SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &       SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &       max(0.0e0_R8,
     &        floi_noc(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1))            !srv 28.10.02
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=21, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 40000                                                   !srv 01.02.09 }

* -- BCENI=22 -- THE MODIFIED DECAY LENGTH FOR THE ION TEMPERATURE   !srv 23.07.10 {

32200   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     &     ' BCENI=22 with ENIPAR(IB,1) =',ENIPAR(IB,1),
     &     ' on ', BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
          t0 = 0.0_R8
            do is = 0, ns-1
            if(.not.is_neutral(is)) then                            !srv 19.03.09 {
                csb = sqrt((rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &                       te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &                       ti(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/
     &                     (am(is)*mp))
              t0 = t0 + csb*na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
              endif                                                   !srv 19.03.09 }
            enddo
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=VBIG*NE00*EV*
     &       CHCI(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &            TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)*
     &         TI(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &            TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=- VBIG*NE00*EV*
     &       (CHCI(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &             TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)-
     &       ENIPAR(IB,1)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*T0)
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
          t0 = 0.0_R8
            do is = 0, ns-1
            if(.not.is_neutral(is)) then                            !srv 19.03.09 {
                csb = sqrt((rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &                       te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &                       ti(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/
     &                     (am(is)*mp))
              t0 = t0 + csb*na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
              endif                                                   !srv 19.03.09 }
            enddo
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=VBIG*NE00*EV*
     &       CHCI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)*
     &       TI(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=- VBIG*NE00*EV*
     &       (CHCI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-ENIPAR(IB,1)*   !srv 02.02.11
     &        SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &           BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*T0)
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=22, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 40000                                                   !srv 23.07.10 }

* -- BCENI=23 -- PRESCRIBE THE POLOIDALLY AVERAGED VALUE OF THE ION TEMPERATURE
* AND INTRODUCE A POLOIDAL VARIATION AS CLOSE AS POSSIBLE TO NEOCLASSICAL SOLUTION.
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH CORRESPONDING
* CONDITION ON IONS DENSITY

32300   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified ion temperature',ENIPAR(IB,1)
!          write all boundaries with this type of boundary condition
          DO IBW=IB,NBC
            IF (BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
              write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                  boundary_location(ibw)
              call xertst(
     .         mod(region(BC_LIST_X(1,IBW),BC_LIST_Y(1,IBW),0),4).eq.1,
     .         'BCENI = 23 only applicable to core boundaries !')
              call xertst(use_astra.ne.0.or.
     .         ENIPAR(IB,1).eq.ENIPAR(IBW,1),
     1         'BCENI = 23, inconsistent ENIPAR(,1) values!')
              if (.not.(bccon(ismain,ib).eq.21 .or. 
     .                  bccon(ismain,ib).eq.22 .or.
     .                  bccon(ismain,ib).eq.23)) then
                write(0,'(a,i3,a,i2,a)') 'Should use BCCON(',ismain,
     &             ',',IB,').EQ.21,22,23 together with BCENI=23'
                write(*,'(a,i3,a,i2,a)') 'Should use BCCON(',ismain,
     &             ',',IB,').EQ.21,22,23 together with BCENI=23'
                write(*,*)
     &           'Please ensure that this boundary condition '//
     &           'is consistent with '//
     &           'ion density boundary condition'  !srv 21.09.12
              end if
              if (ns.gt.1) then
               write(0,*) 
     &       'Use of BCENI=23 only meant for runs without impurities!'
               write(6,*) 
     &       'Use of BCENI=23 only meant for runs without impurities!'
              endif
            endif
            if(BCCHAR(IBW).EQ.'S'.and.
     &         (BCENI(IBW).EQ.16.or.BCENI(IBW).EQ.24)) then
              call xerrab
     .           ('BCENI types 16, 23 and 24 cannot be mixed!')
            endif
          ENDDO
          call xertst(use_astra.ne.0.or.ENIPAR(IBW,1).gt.0.0_R8,
     1     'BCENI = 23, ENIPAR(,1) !> 0')
        endif

        if(.not.once_bceni) then
          if (use_astra.ne.0) THEN
#ifdef ASTRA
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
                ENIPAR(IBW,1) = tib_astra*1.0e3_R8                          !srv 14.06.12
              endif
            ENDDO
#else
            call xerrab
     .       ('Compile with -DASTRA option to couple with ASTRA!')
#endif
          ENDIF

          IF(BCCHAR(IB).EQ.'S' .and.
     $     (mod(region(BC_LIST_X(1,IB),BC_LIST_Y(1,IB),0),4).eq.1)) 
     >     THEN
!          DO IX = 1, BC_LIST_SIZE(IB)
!            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*ENIPAR(IB,1)
!     .       *EV*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
!            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
!     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
!          end do
            t0=0.0e0_R8   ! Volume                                              !lk 12.12.12 {
            t1=0.0e0_R8   ! Poloidally volume-averaged ions temperature in the cell next to boundary
            t2=0.0e0_R8   ! Poloidally volume-averaged electrons density in the cell next to boundary
            t3=0.0e0_R8   ! Poloidally volume-averaged electrons density in the boundary cell
            do ix = -1,nx
              if (mod(region(ix,0,0),4).eq.1) then
                t0=t0 + vol(ix,0)
                t1=t1 + ti(ix,0)*vol(ix,0)
                t2=t2 + ne(ix,0)*vol(ix,0)
                t3=t3 + ne(ix,-1)*vol(ix,0)
              endif
            enddo
            t1=t1/t0
            t2=t2/t0
            t3=t3/t0
            DO IBW=IB,NBC
              IF(BCCHAR(IBW).EQ.'S'.AND.BCENI(IBW).EQ.BCENI(IB)) THEN
                do IX = 1, BC_LIST_SIZE(IBW)
                  if(mod(region(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),0),4).eq.1) then
                    shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)=
     &               shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)
     &                +VBIG*NE00*(ENIPAR(IBW,1)*EV
     &                +(ti(ix,0)-t1)*t1*t2/t3/ENIPAR(IBW,1)/EV)
     &                *SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)=
     &               -VBIG*NE00*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  endif
                enddo  
              ENDIF                                                 !lk 12.12.12 }
            ENDDO
            once_bceni = .true.
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCENI=23 IS AVAILABLE FOR CORE BOUNDARY ONLY.')
          END IF
        endif
        GOTO 40000

* -- BCENI=24 -- Feedback boundary condition with prescribed total ions flux,
*                constant poloidally averaged ions temperature and a
*        POLOIDAL VARIATION AS CLOSE AS POSSIBLE TO NEOCLASSICAL SOLUTION.
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH CORRESPONDING
* CONDITION ON IONS DENSITY

32400   CONTINUE
!IYS 04.03.2013   {
        if(ncall.eq.0) then
          write(*,*) 'Given integrated ion heat flux',
     &               ' with constant Ti '
!          write all boundaries with this type of boundary condition
          DO IBW=IB,NBC
            if (BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
              write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                  boundary_location(ibw)
              call xertst(
     .         mod(region(BC_LIST_X(1,IBW),BC_LIST_Y(1,IBW),0),4).eq.1,
     .         'BCENI = 24 only applicable to core boundaries !')
              call xertst(enipar(ibw,1).ne.0,'faulty parameter enipar1')
              if (.not.(bccon(ismain,ib).eq.21 .or. 
     .                  bccon(ismain,ib).eq.22 .or.
     .                  bccon(ismain,ib).eq.23)) then
                  write(0,'(a,i3,a,i2,a)') 'Should use BCCON(',ismain,
     &             ',',IB,').EQ.21,22,23 together with BCENI=24'
                  write(*,'(a,i3,a,i2,a)') 'Should use BCCON(',ismain,
     &             ',',IB,').EQ.21,22,23 together with BCENI=24'
                  write(*,*)
     &           'Please ensure that this boundary condition '//
     &           'is consistent with '//
     &           'ion density boundary condition'  !srv 21.09.12
              end if
              if (ns.gt.1) then
               write(0,*) 
     &       'Use of BCENI=24 only meant for runs without impurities!'
               write(6,*) 
     &       'Use of BCENI=24 only meant for runs without impurities!'
              endif
            endif
            if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.16
     &                           .and.BCENI(IBW).EQ.23) then
              call xerrab
     .           ('BCENI types 16, 23 and 24 cannot be mixed!')
            endif
          ENDDO
        endif

        if(.not.once_bceni) then
          if (use_astra.ne.0) then
#ifdef ASTRA
            SENIPAR = fhib_astra
#else
            call xerrab
     .       ('Compile with -DASTRA option to couple with ASTRA!')
#endif
          else
            SENIPAR = 0.0_R8
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
                SENIPAR = SENIPAR + ENIPAR(IBW,1)
              endif
            ENDDO
          endif
          write(*,*) ' Summed value is ', SENIPAR

          TIS = 0.0_R8
          NIVS = 0.0_R8
          DO IBW=IB,NBC
            if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
              DO IX = 1, BC_LIST_SIZE(IBW)
! TIS Poloidally volume-averaged ions temperature in the boundary cell
                TIS = TIS + TI(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     .                      NI(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     .                     VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                NIVS = NIVS + NI(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     .                       VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
              END DO
            endif
          END DO
          TIS = TIS/NIVS

          IF(BCCHAR(IB).EQ.'S' .and.
     $     (mod(region(BC_LIST_X(1,IB),BC_LIST_Y(1,IB),0),4).eq.1)) 
     >     THEN
            totflux = 0.0_R8
            t0=0.0e0_R8   ! Volume                                              !lk 12.12.12 {
            t1=0.0e0_R8   ! Poloidally volume-averaged ions temperature in the cell next to boundary
            t2=0.0e0_R8   ! Poloidally volume-averaged electrons density in the cell next to boundary
            t3=0.0e0_R8   ! Poloidally volume-averaged electrons density in the boundary cell
            DO IBW=IB,NBC
              IF(BCCHAR(IBW).EQ.'S'.AND.BCENI(IBW).EQ.BCENI(IB)) THEN
                if(mdf_fhi.ne.0) then
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    totflux = totflux +
     &        fhi_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
     &       +fac_ExB(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*            !srv 11.04.13 {
     &      ((po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))-
     &        po(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &           leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))/
     &       (hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &        hx(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &           leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))+
     &       (po(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &           rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))-
     &        po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))/
     &       (hx(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &           rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))+
     &        hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))*
     &        ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)*
     &        ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &        SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &        bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),2)/
     &        bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),3)**2             !srv 11.04.13 }
                    t0=t0 + vol(ix,0)
                    t1=t1 + ti(ix,0)*vol(ix,0)
                    TIS=TIS + ti(ix,-1)*vol(ix,0)
                    t2=t2 + ne(ix,0)*vol(ix,0)
                    t3=t3 + ne(ix,-1)*vol(ix,0)
                  ENDDO
                else
                  DO IX = 1, BC_LIST_SIZE(IB)
                    totflux = totflux +
     &           fhi(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
                    t0=t0 + vol(ix,0)
                    t1=t1 + ti(ix,0)*vol(ix,0)
                    TIS=TIS + ti(ix,-1)*vol(ix,0)
                    t2=t2 + ne(ix,0)*vol(ix,0)
                    t3=t3 + ne(ix,-1)*vol(ix,0)
                  ENDDO
                endif
              ENDIF
            ENDDO
!             TIS = TIS *
!     &         (1.0_R8-bc_ref*(totflux-SENIPAR)/
!     &           max(abs(totflux),SENIPAR))
            t1=t1/t0
            t2=t2/t0
            t3=t3/t0
            TIS=TIS/t0*
     &          (1.0_R8-bc_ref_ti*(totflux-SENIPAR)/
     &            max(abs(totflux),SENIPAR))
            write(*,"(' BCENI=24: TIS=',es14.7,' totflux=',es14.7,
     *        ' SENIPAR=',es14.7)") TIS/ev, totflux, SENIPAR
            DO IBW=IB,NBC
              IF(BCCHAR(IBW).EQ.'S'.AND.BCENE(IBW).EQ.BCENE(IB)) THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)=
     .             VBIG*NE00*
     $             (TIS+(ti(ix,0)-t1)*t1*t2/t3/TIS)*
     .              SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)=
     .             -VBIG*NE00*
     .              SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                ENDDO
              ENDIF
            ENDDO
            once_bceni = .true.
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCENI=24 IS AVAILABLE FOR CORE ONLY')
          ENDIF
        endif                                                          !srv 16.02.12 }
        GOTO 40000

* -- BOUNDARY CONDITIONS FOR POTENTIAL EQUATION --

40000   CONTINUE
        if(BCPOT(IB).eq.0) goto 50000
        GOTO(40100,40200,40300,40400,40500,
     1   40600,40700,40800,40900,41000,
     2   41100,41200,41300,41400,41500,
     3   41600,41700,41800,41900,42000,
     4   42100),BCPOT(IB)
        write(*,*) BCPOT(IB)
        call xerrab ('BCPOT out of bounds')

* -- BCPOT=1 -- PRESCRIBE THE VALUE OF THE POTENTIAL

40100   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified potential ',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       PBIG*POTPAR(IB,1)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -PBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       PBIG*POTPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-PBIG*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       UBIG*POTPAR(IB,1)*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .      -UBIG*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       UBIG*POTPAR(IB,1)*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-UBIG*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'X') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= UBIG*
     .       POTPAR(IB,1)*VOL(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-PBIG*
     .       VOL(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'Y') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= UBIG*
     .       POTPAR(IB,1)*VOL(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-PBIG*
     .       VOL(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCPOT=1, INVALID OPTION.')
        ENDIF
        GOTO 50000

* -- BCPOT=2 -- PRESCRIBE THE GRADIENT OF THE POTENTIAL

40200   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'potential gradient ',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       PBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .       (-POTPAR(IB,1)*0.5_R8*
     .       (HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     .        HY1(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))) +
     .       PO(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -PBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= PBIG*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .       (POTPAR(IB,1)*0.5_R8*
     .       (HY1(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     .        HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))) +
     .       PO(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-PBIG*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       UBIG*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     .       (-POTPAR(IB,1)*0.5_R8*
     .       (HX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))+
     .        HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .           RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))) +
     .        PO(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .           RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .      -UBIG*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= UBIG*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))*
     .       (POTPAR(IB,1)*0.5_R8*
     .       (HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .           LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))+
     .        HX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))) +
     .       PO(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-UBIG*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          end do
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCPOT=2, INVALID OPTION.')
        END IF
        GOTO 50000

* -- BCPOT=3 -- SHEATH CONDITIONS, ELECTRON ENERGY TRANSMISSION

40300   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,a,a,1p,g14.7,a)')
     .     'sheath potential bc on ',BCCHAR(IB),
     1     boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            s1 = abs(pbs(rightix(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .                   rightiy(bc_list_x(iy,ib),bc_list_y(iy,ib)),0))
            vte = sqrt(bcint(te,nx,ny,iy,ib)/me)
            t0 = 0.0_R8
            do is = 0, ns-1
              t0 = t0-bcint(rza(-1,-1,is),nx,ny,iy,ib)*
     1         fna(rightix(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .             rightiy(bc_list_x(iy,ib),bc_list_y(iy,ib)),0,is)
            enddo
            do is = 0, ns_ext-1
              t0 = t0-bcint(za_ext(-1,-1,is),nx,ny,iy,ib)*
     1         fa_ext(rightix(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .                rightiy(bc_list_x(iy,ib),bc_list_y(iy,ib)),0,is)
            enddo
            fchi = qe*t0
            if(t0.lt.0.0_R8 .and.
     .       bottomiy(bc_list_x(iy,ib),bc_list_y(iy,ib)).ne.-2 .and.
     &       topiy(bc_list_x(iy,ib),bc_list_y(iy,ib)).ne.ny+1) then
              wrong_flow = .true.                                     !xpb
              if(diagno.ge.0) then
                write(*,'(a,2i4,1x,1p,g14.7)') 'b2stbc_west:ix,iy,fchi',
     1           bc_list_x(iy,ib),bc_list_y(iy,ib),fchi
                write(hlp_frm,'(a,i3,a)') '(a,1p,',ns,'(g14.7,1x))'
                write(*,hlp_frm) 'rza : ',
     .           (rza(rightix(bc_list_x(iy,ib),
     .                        bc_list_y(iy,ib)),
     .                rightiy(bc_list_x(iy,ib),
     .                        bc_list_y(iy,ib)),is),is=0,ns-1)
                write(*,hlp_frm) 'fna (A) : ',
     &           (fna(rightix(bc_list_x(iy,ib),
     .                        bc_list_y(iy,ib)),
     .                rightiy(bc_list_x(iy,ib),
     .                        bc_list_y(iy,ib)),0,is)*qe,is=0,ns-1)
                write(*,'(a,1p,4g14.7)') 'Bx,Bz,Ud,-Ua =',
     &           wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0),
     .           wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),2),
     &           ud(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain),
     .          -ua(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain)
               endif
c              fchi=0.0_R8
            endif
            fche = sqrt(1/(2*pi))*qe*bcint(ne,nx,ny,iy,ib)*vte*s1*
     &       expu2(-qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
     .       bcint(te,nx,ny,iy,ib))
            if (secmodel.eq.1) then
              if (ne(bc_list_x(iy,ib),bc_list_y(iy,ib))/
     .           (rza(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain)*
     .            na(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain))
     .           .gt.1.1_R8) then
                write(*,*) 'Large ne contribution from impurities!'
                write(*,*) 'Kinetic s.e.c. model may not be valid!'
              endif
              call eseec0 (bcint(te,nx,ny,iy,ib),
     &         bcint(ti,nx,ny,iy,ib),bcint(ne,nx,ny,iy,ib),
     &         bcint(po,nx,ny,iy,ib),
     .         wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0),
     &         wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3),0,
     .         bc_list_x(iy,ib),bc_list_y(iy,ib),
     &         nx,ny,ns,ismain,seec,diagno)                           !cfs
            else
              seec = gammae
            endif
            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/bcint(te,nx,ny,iy,ib)
*    ..compute charge source
            sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),0) =
     .       sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),0)+
     &       (1.0_R8-seec)*fche-fchi+
     .       t0*po(bc_list_x(iy,ib),bc_list_y(iy,ib))
            sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),1) =
     .      sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),1)-t0
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            s1 = abs(pbs(bc_list_x(iy,ib),bc_list_y(iy,ib),0))
            vte = sqrt(bcint(te,nx,ny,iy,ib)/me)
            t0 = 0.0_R8
            do is = 0, ns-1
              t0 = t0+
     &         bcint(rza(-1,-1,is),nx,ny,iy,ib)*
     .         fna(bc_list_x(iy,ib),bc_list_y(iy,ib),0,is)
            enddo
            do is = 0, ns_ext-1
              t0 = t0+bcint(za_ext(-1,-1,is),nx,ny,iy,ib)*
     &         fa_ext(bc_list_x(iy,ib),bc_list_y(iy,ib),0,is)
            enddo
            fchi = qe*t0
            if(t0.lt.0.0_R8 .and.
     .       bottomiy(bc_list_x(iy,ib),bc_list_y(iy,ib)).ne.-2 .and.
     &       topiy(bc_list_x(iy,ib),bc_list_y(iy,ib)).ne.ny+1) then
              wrong_flow = .true.                                     !xpb
              if(diagno.ge.0) then
                write(*,'(a,2i4,1x,1p,g14.7)') 'b2stbc_east:ix,iy,fchi',
     1           bc_list_x(iy,ib),bc_list_y(iy,ib),fchi
                write(hlp_frm,'(a,i3,a)') '(a,1p,',ns,'(g14.7,1x))'
                write(*,hlp_frm) 'rza : ',
     &           (rza(bc_list_x(iy,ib),bc_list_y(iy,ib),is),is=0,ns-1)
                write(*,hlp_frm) 'fna (A) : ',
     .           (fna(bc_list_x(iy,ib),bc_list_y(iy,ib),0,is)*qe,
     .            is=0,ns-1)
                write(*,'(a,1p,4g14.7)') 'Bx,Bz,Ud,Ua =',
     &           wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0),
     .           wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),2),
     &           ud(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain),
     .           ua(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain)
               endif
c              fchi=0.0_R8
            endif
            fche = sqrt(1/(2*pi))*qe*bcint(ne,nx,ny,iy,ib)*vte*s1*
     &       expu2(-qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
     .       bcint(te,nx,ny,iy,ib))
            if (secmodel.eq.1) then
              if (ne(bc_list_x(iy,ib),bc_list_y(iy,ib))/
     .           (rza(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain)*
     .            na(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain))
     .           .gt.1.1_R8) then
                write(*,*) 'Large ne contribution from impurities!'
                write(*,*) 'Kinetic s.e.c. model may not be valid!'
              endif
              call eseec0 (bcint(te,nx,ny,iy,ib),
     &         bcint(ti,nx,ny,iy,ib),bcint(ne,nx,ny,iy,ib),
     &         bcint(po,nx,ny,iy,ib),
     .         wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0),
     &         wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),3),0,
     .         bc_list_x(iy,ib),bc_list_y(iy,ib),
     &         nx,ny,ns,ismain,seec,diagno)                           !cfs
            else
              seec = gammae
            endif
            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/bcint(te,nx,ny,iy,ib)
*    ..compute charge source
            sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),0) =
     .      sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),0)+
     &       (1.0_R8-seec)*fche-fchi+
     .       t0*po(bc_list_x(iy,ib),bc_list_y(iy,ib))
            sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),1) =
     .      sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),1)-t0
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=3, NO SHEATH IN RADIAL DIR.')
        END IF
        GOTO 50000

* -- BCPOT=4 -- PRESCRIBE THE VALUE OF THE POTENTIAL,
*                 WEAKLY A MIXED BOUNDARY CONDITION

40400   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'weak potential ',POTPAR(IB,1),' on ',BCCHAR(IB),
     1     boundary_location(ib)
          call xertst(POTPAR(IB,2).gt.0.0_R8,
     1     'BCPOT = 4, POTPAR(,2) !> 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= POTPAR(IB,2)*
     .       POTPAR(IB,1)*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-POTPAR(IB,2)*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= POTPAR(IB,2)*
     .       POTPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-POTPAR(IB,2)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= POTPAR(IB,2)*
     .       POTPAR(IB,1)*
     1       SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-POTPAR(IB,2)*
     .       SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= POTPAR(IB,2)*
     1       POTPAR(IB,1)*SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                       LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-POTPAR(IB,2)*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=4, NOT YET AVAILABLE')
        END IF
        GOTO 50000

* -- BCPOT=5 -- PRESCRIBE THE CURRENT FLUX DENSITY PER UNIT AREA

40500   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'current flux density per area ',POTPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=POTPAR(IB,1)*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=POTPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=POTPAR(IB,1)*
     .       SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=POTPAR(IB,1)*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=5, INVALID OPTION')
        END IF
        GOTO 50000

* -- BCPOT=6 -- PRESCRIBE THE TOTAL CURRENT FLUX DENSITY FOR A
*                 CONSTANT POTENTIAL

40600   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'total current flux density',POTPAR(IB,1),
     1     ' with const potential on ',BCCHAR(IB),boundary_location(ib)
        endif
        CALL XERRAB ('B2STBC -- BCPOT=6, NOT YET AVAILABLE')
        GOTO 50000

* -- BCPOT=7 -- PRESCRIBE THE POTENTIAL AS A FUNCTION OF
*                 OTHER PLASMA PARAMETERS

40700   CONTINUE
        CALL XERRAB ('B2STBC -- BCPOT=7, NOT YET AVAILABLE')
        GOTO 50000

* -- BCPOT=8 -- PRESCRIBE THE TOTAL ELECTRON CURRENT WITH CONSTANT
*                 FLUX DENSITY

40800   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'total current ',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          US=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            US=US+SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          end do
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=POTPAR(IB,1)/US*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          US=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            US=US+SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .               BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          end do
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=POTPAR(IB,1)/US*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          US=0.0_R8
          DO IY = 1, BC_LIST_SIZE(IB)
            US=US+SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          end do
          DO IY = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=POTPAR(IB,1)/US*
     .       SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          US=0.0_R8
          DO IY = 1, BC_LIST_SIZE(IB)
            US=US+SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          end do
          DO IY = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=POTPAR(IB,1)/US*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=8, INVALID OPTION')
        ENDIF
        GOTO 50000

* -- BCPOT=9 -- PRESCRIBE THE DECAY LENGTH FOR
*                 THE POTENTIAL

40900   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'potential decay length ',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(POTPAR(IB,1).ne.0.0_R8,
     1     'BCPOT = 9, POTPAR(,1) = 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -sig0(bc_list_x(ix,ib),bc_list_y(ix,ib))/POTPAR(IB,1)*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/qe/
     .       ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -sig0(bc_list_x(ix,ib),bc_list_y(ix,ib))/POTPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/qe/
     2       ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=9, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 50000

* -- BCPOT=10 -- Feedback option for core
*                NOT YET TESTED!!!!!!!!!

41000   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,a,a)') 'current feedback bc on ',
     1     BCCHAR(IB),boundary_location(ib)
          write(*,*) 'NOT YET TESTED !!!!!!!!!!'
          call xertst(lfeedback, 'LFEEDBACK must be turned on!')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if (region(bc_list_x(ix,ib),bc_list_y(ix,ib),0).eq.1) then
              sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        cbsch(0,coreregno)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            elseif (region(bc_list_x(ix,ib),bc_list_y(ix,ib),0).eq.5)
     .       then
              sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        cbsch(0,coreregn2)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            endif
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=10, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 50000

* -- BCPOT=11 -- SHEATH CONDITIONS, ELECTRON ENERGY TRANSMISSION FROM b2stbc_spb !srv 01.02.09 {

41100   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,a,1p,g14.7,a,g14.7,a,a,a)')
     &     'sheath potential bc from b2stbc_spb',
     &     ' with POTPAR(IB,2) =',POTPAR(IB,2),' gammae =', gammae,
     &     'on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            s1 = abs(pbs(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                   rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0))
!            vte = sqrt(bcint(te,nx,ny,iy,ib)/me)
            vte = sqrt(te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/me)
            wrong_flow=.false.
*    ..set boundary condition in explicit form (the same form as fna)
            CS=SQRT(PR(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &              RO(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            t0 = 0.0e0_R8
            do is = 0, ns-1
              if(istyle_fchi.eq.1 .and. BCCON(IS,IB).eq.14) then        !srv 14.02.05 {
                t0 = t0-rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*
     *           CONPAR(IS,IB,1)*cs*s1*
     *           na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
              else
                t0 = t0-rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*
     &           fna(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &               rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)   !srv 26.11.02
              endif                                                     !srv 14.02.05 }
            enddo
            fchi = qe*t0
!            t0 = 0
!            do is = 0, ns-1
!              t0 = t0-bcint(rza(-1,-1,is),nx,ny,iy,ib)*
!     1         fna(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
!     2             rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)
!            enddo
!            fchi = qe*t0
!            if(t0.lt.0 .and. bottomiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)).ne.-2 .and.
!     &       topiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)).ne.ny+1) then
!              wrong_flow = .true.                                     !xpb
!              if(diagno.ge.0) then
!                write(*,*) 'b2stbc_west:ix,iy,fchi',
!     1           BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),fchi
!                write(*,*)
!     .         rza(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
!     .             rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),ismain),
!     &         fna(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
!                   rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,ismain)
!                write(*,*) 'Bx,Bz,Uadia,-Ua =',
!     &           wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0),
!     &           wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),2),
!     &           uadia(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,ismain),
!     &          -ua(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),ismain)
!              endif
c              fchi=0.0_R8
!            endif
!            fche = sqrt(1/(2*pi))*qe*bcint(ne,nx,ny,iy,ib)*vte*s1*
!     &       expu2(-qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
!     .       bcint(te,nx,ny,iy,ib))
            fche = sqrt(1/(2*pi))*qe*
     &       ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*vte*s1*
     &       expu2(-qe*(po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))-
     &                  potpar(ib,2))/
     &       te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            if (secmodel.eq.1) then
              if (ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     .          (rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),ismain)*
     .            na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),ismain))
     >         .gt.1.1_R8) then
                write(*,*) 'Large ne contribution from impurities!'
                write(*,*) 'Kinetic s.e.c. model may not be valid!'
              endif
              call eseec0 (bcint(te,nx,ny,iy,ib),
     &         bcint(ti,nx,ny,iy,ib),bcint(ne,nx,ny,iy,ib),
     &         bcint(po,nx,ny,iy,ib),
     &         wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0),
     &         wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3),0,
     &         BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),
     &         nx,ny,ns,ismain,seec,diagno)                           !cfs
            else
              seec = gammae
            endif
!            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/bcint(te,nx,ny,iy,ib)
            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/
     &       te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
*    ..compute charge source
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &       sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)+
     &       (1.0_R8-seec)*fche-fchi+
     &       t0*po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &       sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-t0
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            s1 = abs(pbs(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))
!            vte = sqrt(bcint(te,nx,ny,iy,ib)/me)
!            t0 = 0
!            do is = 0, ns-1
!              t0 = t0+
!     &         bcint(rza(-1,-1,is),nx,ny,iy,ib)*
!     &         fna(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)
!            enddo
!            fchi = qe*t0
!            if(t0.lt.0 .and. bottomiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)).ne.-2 .and.
!     &       topiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)).ne.ny+1) then
!              wrong_flow = .true.                                     !xpb
!              if(diagno.ge.0) then
!                write(*,*) 'b2stbc_east:ix,iy,fchi',
!     1           BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),fchi
!                write(*,*)
!     &           rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),ismain),
!     &           fna(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,ismain)
!                write(*,*) 'Bx,Bz,Uadia,Ua =',
!     &           wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0),
!     &           wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),2),
!     &           uadia(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,ismain),
!     &           ua(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),ismain)
!              endif
c              fchi=0.0_R8
!            endif
!            fche = sqrt(1/(2*pi))*qe*bcint(ne,nx,ny,iy,ib)*vte*s1*
!     &       expu2(-qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
!     .       bcint(te,nx,ny,iy,ib))
            vte = sqrt(te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/me)
            wrong_flow=.false.
            CS=SQRT(PR(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &              RO(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            t0 = 0.0e0_R8
            do is = 0, ns-1
              if(istyle_fchi.eq.1 .and. BCCON(IS,IB).eq.14) then        !srv 14.02.05 {
*    ..set boundary condition in explicit form (the same form as fna)
                t0 = t0-rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*
     &           CONPAR(IS,IB,1)*cs*s1*
     &           na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
              else                                                      !srv 14.02.05 {
                t0 = t0+rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*
     &                  fna(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)     !srv 26.11.02
              endif                                                     !srv 14.02.05 }
            enddo
            fchi = qe*t0
            fche = sqrt(1/(2*pi))*qe*
     &       ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*vte*s1*
     &       expu2(-qe*(po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))-
     &                  potpar(ib,2))/
     &       te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            if (secmodel.eq.1) then
              if (ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &          (rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),ismain)*
     .            na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),ismain))
     >        .gt.1.1_R8) then
                write(*,*) 'Large ne contribution from impurities!'
                write(*,*) 'Kinetic s.e.c. model may not be valid!'
              endif
              call eseec0 (bcint(te,nx,ny,iy,ib),
     &         bcint(ti,nx,ny,iy,ib),bcint(ne,nx,ny,iy,ib),
     &         bcint(po,nx,ny,iy,ib),
     &         wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0),
     &         wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3),0,
     &         BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),
     &         nx,ny,ns,ismain,seec,diagno)                           !cfs
            else
              seec = gammae
            endif
!            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/bcint(te,nx,ny,iy,ib)
            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/
     &       te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
*    ..compute charge source
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &       sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)+
     &       (1.0_R8-seec)*fche-fchi+
     &       t0*po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &       sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-t0
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=11, NO SHEATH IN RADIAL DIR.')
        END IF
        GOTO 50000                                                    !srv 01.02.09 }

* -- BCPOT=12 -- Nothing implemented

41200   CONTINUE
        CALL XERRAB ('B2STBC -- BCPOT=12, NOT AVAILABLE FOR S,N,W,E')
        GOTO 50000

* -- BCPOT=13 -- Nothing implemented

41300   CONTINUE
        CALL XERRAB ('B2STBC -- BCPOT=13, NOT AVAILABLE FOR S,N,W,E')
        GOTO 50000

* -- BCPOT=14 -- Nothing implemented

41400   CONTINUE
        CALL XERRAB ('B2STBC -- BCPOT=14, NOT AVAILABLE FOR S,N,W,E')
        GOTO 50000

* -- BCPOT=15 -- Nothing implemented

41500   CONTINUE
        CALL XERRAB ('B2STBC -- BCPOT=15, NOT AVAILABLE FOR S,N,W,E')
        GOTO 50000

* -- BCPOT=16 -- Nothing implemented

41600   CONTINUE
        CALL XERRAB ('B2STBC -- BCPOT=16, NOT AVAILABLE FOR S,N,W,E')
        GOTO 50000

* -- BCPOT=17 -- Nothing implemented

41700   CONTINUE
        CALL XERRAB ('B2STBC -- BCPOT=17, NOT AVAILABLE FOR S,N,W,E')
        GOTO 50000

* -- BCPOT=18 -- Nothing implemented

41800   CONTINUE
        CALL XERRAB ('B2STBC -- BCPOT=18, NOT AVAILABLE FOR S,N,W,E')
        GOTO 50000

* -- BCPOT=19 -- Nothing implemented

41900   CONTINUE
        CALL XERRAB ('B2STBC -- BCPOT=19, NOT AVAILABLE FOR S,N,W,E')
        GOTO 50000

* -- BCPOT=20 -- Nothing implemented

42000   CONTINUE
        CALL XERRAB ('B2STBC -- BCPOT=20, NOT AVAILABLE FOR S,N,W,E')
        GOTO 50000

* -- BCPOT=21 -- constant potential scaled by potential on the ring bc_type21_ref away

42100   CONTINUE
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'fixed potential (modified by nearby ring) ',
     1     POTPAR(IB,1),' with ',POTPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(potpar(ib,2).gt.0.0_R8,
     1     'faulty parameter potpar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          poave = 0.0_R8
          neave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            poave=poave+po(bc_type21_ix,bc_type21_iy)*
     1       ne(bc_type21_ix,bc_type21_iy)*
     2       vol(bc_type21_ix,bc_type21_iy)
            neave=neave+ne(bc_type21_ix,bc_type21_iy)*
     1       vol(bc_type21_ix,bc_type21_iy)
          ENDDO
          poave=poave/neave
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*POTPAR(IB,1)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .       ((1.0_R8-potpar(ib,2))+
     .           potpar(ib,2)*po(bc_type21_ix,bc_type21_iy)/poave)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE IF(BCCHAR(IB).EQ.'N') THEN
          poave = 0.0_R8
          neave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            poave=poave+po(bc_type21_ix,bc_type21_iy)*
     1       ne(bc_type21_ix,bc_type21_iy)*
     2       vol(bc_type21_ix,bc_type21_iy)
            neave=neave+ne(bc_type21_ix,bc_type21_iy)*
     1       vol(bc_type21_ix,bc_type21_iy)
          ENDDO
          poave=poave/neave
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*POTPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .       ((1.0_R8-potpar(ib,2))+
     .           potpar(ib,2)*po(bc_type21_ix,bc_type21_iy)/poave)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=21, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 50000

50000   continue

        if(istyle_cur_contr_on_S_and_N.eq.0) then                       !srv 01.02.09
!xpb
!xpb  Adding fchy_dia contributions on south and north boundaries
!xpb  following "sv 01.07.99" coding (see B2.5-style implementation)
!xpb  Special attention is paid to the case where there is a discontinuity in B
!xpb
         IF (BCPOT(IB).EQ.1) THEN
          if (ncall.eq.0) write(*,'(a,i2,a)')
     &          'Disabled special treatment of fchy_dia on boundary ',
     &          ib,' because bcpot there is 1'
         ELSE IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if(maxval(facdrift) .ne. 0.0_R8) then
              if(nnreg(0).gt.1 .and.
     .         mod(region(bc_list_x(ix,ib),bc_list_y(ix,ib),0),4).eq.1)
     .         then
                if(periodic_bc.ne.1 .or.
     .             rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)).eq.
     .             bc_list_x(ix,ib)+1 .or. nncut.eq.0) then
                  fchy_dia = -praverage(1)*
     .             facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))
     &             *b2stbc_fchy_dia                                   !sv 01.07.99
     &             *(OnedBsqc(rightix(topix(bc_list_x(ix,ib),
     .                                      bc_list_y(ix,ib)),
     1                                topiy(bc_list_x(ix,ib),
     .                                      bc_list_y(ix,ib))),
     &                        rightiy(topix(bc_list_x(ix,ib),
     .                                      bc_list_y(ix,ib)),
     1                                topiy(bc_list_x(ix,ib),
     .                                      bc_list_y(ix,ib))))
     &              -OnedBsqc(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                        topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))      !sv
                else
      if(ncall.eq.0)
     > write(*,*) 'b2stbc: special treatment for cell ix,iy = ',
     . bc_list_x(ix,ib),bc_list_y(ix,ib)
                  fchy_dia = -praverage(1)
     .             *facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))
     &             *b2stbc_fchy_dia  !sv 01.07.99
     &             *2.0_R8*(1.0_R8/
     1              wbbv(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)**2
     &             -OnedBsqc(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                       topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))   !xpb
                endif
              else
                if(periodic_bc.ne.1 .or. nncut.eq.0 .or.
     &     rightix(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))).eq.
     .             bc_list_x(ix,ib)+1)
     >           then
                  fchy_dia = -bcint(wdia,nx,ny,ix,ib)
     .             *(vol(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     .               vol(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .             /(gs(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                  topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2)+
     .               gs(bc_list_x(ix,ib),bc_list_y(ix,ib),2))
     .      *(OnedBsqc(rightix(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     .                 rightiy(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .             -OnedBsqc(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                       topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .             *facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .              b2stbc_fchy_dia
                else
      if (ncall.eq.0)
     > write(*,*) 'b2stbc: special treatment for cell ix,iy = ',
     . bc_list_x(ix,ib),bc_list_y(ix,ib)
                 fchy_dia = -bcint(wdia,nx,ny,ix,ib)
     .             *(vol(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     .               vol(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .             /(gs(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                  topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2)+
     .               gs(bc_list_x(ix,ib),bc_list_y(ix,ib),2))*2.0_R8
     .     *(1.0_R8/wbbv(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)**2
     .     -OnedBsqc(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .             *facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .              b2stbc_fchy_dia
                endif
              endif
            else
              if(b2stbc_neoclassical .ne. 0.0_R8 .and. nnreg(0).gt.1
     &            .and. mod(region(bc_list_x(ix,ib),
     .                             bc_list_y(ix,ib),0),4).eq.1) then
                if(periodic_bc.ne.1 .or. nncut.eq.0 .or.
     &     rightix(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))).eq.
     .             bc_list_x(ix,ib)+1)
     >           then
                  fchy_dia = -praverage(1)
     .    *(OnedBsqc(rightix(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                       topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     .               rightiy(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                       topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .              -OnedBsqc(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                        topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .              *b2stbc_neoclassical
                else
                  fchy_dia = -praverage(1)*2.0_R8*
     .     *(1.0_R8/wbbv(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)**2
     .              -OnedBsqc(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                        topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .              *b2stbc_neoclassical
                endif
              else
                fchy_dia = 0.0_R8
              endif
            endif
            if(mod(region(bc_list_x(ix,ib),bc_list_y(ix,ib),0),4).ne.1
     .                      .and.
     1       b2stbc_fchy_dia_coreonly.eq.1) fchy_dia = 0.0_R8
            sch0(bc_list_x(ix,ib),bc_list_y(ix,ib),0) =
     .      sch0(bc_list_x(ix,ib),bc_list_y(ix,ib),0) + fchy_dia
          END DO
         ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if(maxval(facdrift) .ne. 0.0_R8) then
              if(region(bc_list_x(ix,ib),bc_list_y(ix,ib),0).eq.5) then
                fchy_dia = -praverage(2)
     .           *facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))
     &           *b2stbc_fchy_dia     !sv 01.07.99
     &           *(OnedBsqc(rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                      rightiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     1             OnedBsqc(bc_list_x(ix,ib),bc_list_y(ix,ib)))   !sv
              else
                fchy_dia = -bcint(wdia,nx,ny,ix,ib)
     .           *(vol(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .             vol(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .           /(gs(bc_list_x(ix,ib),bc_list_y(ix,ib),2)+
     .             gs(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2))
     .           *(OnedBsqc(rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                      rightiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .            -OnedBsqc(bc_list_x(ix,ib),bc_list_y(ix,ib)))*
     .             facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))
     &            *b2stbc_fchy_dia
              endif
            else
              if(b2stbc_neoclassical .ne. 0.0_R8 .and.
     &           region(bc_list_x(ix,ib),bc_list_y(ix,ib),0).eq.5) then
                fchy_dia = -praverage(2)
     .           *(OnedBsqc(rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                      rightiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .            -OnedBsqc(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .           *b2stbc_neoclassical
              else
                fchy_dia = 0.0_R8
              endif
            endif
            if(mod(region(bc_list_x(ix,ib),bc_list_y(ix,ib),0),4).ne.1
     .                      .and.
     1       b2stbc_fchy_dia_coreonly.eq.1) fchy_dia = 0.0_R8
            sch0(bc_list_x(ix,ib),bc_list_y(ix,ib),0) =
     .      sch0(bc_list_x(ix,ib),bc_list_y(ix,ib),0) - fchy_dia
          END DO
         ENDIF
!xpb
       else                                                         !srv 01.02.09 {
         IF (BCPOT(IB).EQ.1) THEN
          if (ncall.eq.0) write(*,'(a,i2,a)')
     &          'Disabled special treatment of fchy_dia on boundary ',
     &          ib,' because bcpot there is 1'
         ELSE IF(BCCHAR(IB).EQ.'S') THEN
           DO IX = 1, BC_LIST_SIZE(IB)
             if(mod(region(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0),4).eq.1)
     >        then
*        for inner surfaces
               if(maxval(facdrift) .ne. 0.0_R8) then
                fchy_dia =                                             !srv 01.07.99
     &          -pzaverage
     &          *facdrift(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*phm0      ! *signmf        !srv 24.03.99 01.07.99 21.04.00 18.01.02 2
     &          *bzhzaverage                                           !lk 20.11.07 }
     &          *(OnedBsqc(
     &            rightix(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))),
     &            rightiy(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))-
     &            OnedBsqc(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                     topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))   !srv 01.07.99
                sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) + fchy_dia
                fchy_inertia = 0.0_R8                                   !srv 25.06.03
              do is = 0, ns-1                                           !srv 25.06.03
                if(.not.is_neutral(is))then                             !srv 25.06.03
                    fchy_inertia = fchy_inertia                         !srv 07.04.00 25.06.03
     &               -0.5e0_R8*am(is)*mp*
     &                facdrift(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*phm1  ! *signmf                 !srv 07.04.00 21.04.00 21.10.02
     &               *naaverage(is)                                     !lk 20.11.07                       !srv 07.04.00
     &               *ua2average(is)                                    !srv 07.04.00 15.02.02
     &               *bzhzaverage                                       !lk 20.11.07
     &               *(OnedBsqc(
     &            rightix(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))),
     &            rightiy(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))-
     &           OnedBsqc(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))    !srv 07.04.00
                endif
              enddo
                sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)+fchy_inertia
               endif
             else
*        for PR region
               sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &          sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &           fchdia(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)+
     &           fchin(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)+
     &           fchvispar(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                     topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)+
     &           fchvisq(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                   topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)+
     &           fchvisper(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                     topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)+
     &           fchinert(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)
           endif
           END DO
         ELSEIF(BCCHAR(IB).EQ.'N') THEN
           DO IX = 1, BC_LIST_SIZE(IB)
             sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &        sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &             fchdia(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &             fchin(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &             fchvispar(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &             fchvisq(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &             fchvisper(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &             fchinert(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)
           END DO
         ENDIF
        endif                                                        !srv 01.02.09 }

      enddo

* ---------------------------------------------------------------------
*  INTERIOR BOUNDARY AT DEAD REGIONS
*  ---------------------------------------------------------------------
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            DO IS=0, NS-1
              SNA0(IX,IY,0,IS)= UBIG*NIISO(IS)
              SNA0(IX,IY,1,IS)=-UBIG
              smo0(ix,iy,0,is)= ubig*ne00*am(is)*0.0_R8*
     .         sx(min(nx-1,ix),iy)
              smo0(ix,iy,1,is)=-ubig*ne00*am(is)*
     .         sx(min(nx-1,ix),iy)
            end do
            she0(IX,IY,0)= UBIG*NE(IX,IY)*TEISO*EV
            she0(IX,IY,1)=-UBIG*NE(IX,IY)
            shi0(IX,IY,0)= UBIG*NE(IX,IY)*TIISO*EV
            shi0(IX,IY,1)=-UBIG*NE(IX,IY)
            sch0(IX,IY,0)= UBIG*PHIISO
            sch0(IX,IY,1)=-UBIG
          end do
        end do
      end do

      if (iout.eq.1 .or. iout_b2npmo.eq.1) then                       !srv 11.09.09 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
!        do k=0,3
!         write (chk,'(i1)') k
!         call my_out(70,nx,ny,smo0(-1,-1,k,is),
!     &               'b2stbc_phys_smo'//chk//chns)
!        enddo
       call my_out(70,nx,ny,smo0(:,:,0,is)+
     &        smo0(:,:,1,is)*ua(:,:,is)+
     &        smo0(:,:,2,is)*na(:,:,is)*mp*am(is)+
     &        smo0(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is),
     &        'b2stbc_phys_smo'//chns)
       enddo
      endif
      if (iout.eq.1 .or. (iout_b2npc9.eq.1 .or. iout_b2npco.eq.1)) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
!        do k=0,1
!         write (chk,'(i1)') k
!         call my_out(70,nx,ny,sna0(-1,-1,k,is),
!     &               'b2stbc_phys_sna'//chk//chns)
!        enddo
       call my_out(70,nx,ny,sna0(:,:,0,is)+
     &             sna0(:,:,1,is)*na(:,:,is),'b2stbc_phys_sna'//chns)
       enddo
      endif
      if (iout.eq.1 .or. (iout_b2npht.eq.1 .or. iout_b2nph9.eq.1)) then
!       do k=0,3
!        write (chk,'(i1)') k
!        call my_out(70,nx,ny,she0(-1,-1,k),'b2stbc_phys_she'//chk)
!        call my_out(70,nx,ny,shi0(-1,-1,k),'b2stbc_phys_shi'//chk)
!       enddo
       call my_out(70,nx,ny,she0(:,:,0)+she0(:,:,1)*te+
     &   she0(:,:,2)*ne + she0(:,:,3)*te*ne,'b2stbc_phys_she')
       call my_out(70,nx,ny,shi0(:,:,0)+shi0(:,:,1)*ti+
     &   shi0(:,:,2)*ni(:,:,0)+shi0(:,:,3)*ti*ni(:,:,0),
     &   'b2stbc_phys_shi')
      endif
      if (iout.eq.1 .or. (iout_b2npp7.eq.1 .or. iout_b2nppo.eq.1)) then
       do k=0,3
        write (chk,'(i1)') k
        call my_out(70,nx,ny,sch0(-1,-1,k),'b2stbc_phys_sch'//chk)
       enddo
!       call my_out(70,nx,ny,sch0(:,:,0)+sch0(:,:,1)*po+
!     &   sch0(:,:,2)*ne+sch0(:,:,3)*po*ne(:,:),'b2stbc_phys_sch')
      endif                                                           !srv 11.09.09 }

* ..return
      if(ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbc_phys

      end subroutine b2stbc_phys
