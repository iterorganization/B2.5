*-----------------------------------------------------------------------
*.specification
      subroutine b2stbc_phys (nx, ny, ns, ismain, ismain0,
     1  BoRiS, facdrift, fac_ExB,
     1  na, ua, uadia, vedia, vadia, veecrb, vaecrb,
     &  te, ti, po, ne, ni, floe_noc, floi_noc, fna, fna_mdf,         !srv 28.10.02
     1  fnaPSch, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,      !srv 13.10.06
     &  fhe, fhi, fhm, fhe_mdf, fhi_mdf, fna_fcor, fhePSch, fhiPSch,
     &  fch, fne, fni,                                                !srv 13.10.06
     &  fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,!srv 21.10.02
     &  fchanml,
     2  cdna, chce, chve, chci, chvi, calf,                           !srv 25.09.09
     &  sna0, smo0, she0, shi0, sch0, sne0,
     &  wrong_flow, main_call, diagno,
     &  b2stbc_fchy_dia, b2stbc_neoclassical, praverage,
     1  b2stbc_fchy_dia_coreonly, secmodel, area_fix,
     2  sheath_drift_fix)
      use b2mod_types
      use b2mod_geo
      use b2mod_math
      use b2mod_rates
      use b2mod_indirect
      use b2mod_feedback
      use b2mod_anomalous_transport
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_time
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_external
#ifdef B25_EIRENE
      use eirmod_wneutrals
      use eirmod_braeir
#endif
#ifdef PAR_OPT
      use b2mod_par_opt
#endif
      use b2mod_astra_to_b2                                          !srv 28.07.11
      use b2mod_subsys
      implicit none
#     include "DIMENSIONS.F"
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, ismain0, diagno
      real (kind=R8) ::
     *  BoRiS, facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     1  uadia(-1:nx,-1:ny,0:1,0:ns-1),
     &  vedia(-1:nx,-1:ny,0:1), vadia(-1:nx,-1:ny,0:1,0:ns-1),
     2  veecrb(-1:nx,-1:ny,0:1), vaecrb(-1:nx,-1:ny,0:1,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1), fna_32(-1:nx,-1:ny,0:1,0:ns-1),
     *  fna_mdf(-1:nx,-1:ny,0:1,0:ns-1), fhm(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_52(-1:nx,-1:ny,0:1,0:ns-1), fni_32(-1:nx,-1:ny,0:1),
     &  fni_52(-1:nx,-1:ny,0:1), fne_32(-1:nx,-1:ny,0:1),
     &  fne_52(-1:nx,-1:ny,0:1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fhe_mdf(-1:nx,-1:ny,0:1),
     &  fhi_mdf(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     &  fch_p(-1:nx,-1:ny,0:1), fchvispar(-1:nx,-1:ny,0:1),          !srv 11.10.02 21.10.02 {
     &  fchvisper(-1:nx,-1:ny,0:1), fchin(-1:nx,-1:ny,0:1),
     &  fchvisq(-1:nx,-1:ny,0:1), fchdia(-1:nx,-1:ny,0:1),
     &  fchinert(-1:nx,-1:ny,0:1), fchanml(-1:nx,-1:ny,0:1),         !srv 21.10.02 }
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
      logical main_call                                              !xpb
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3), sne0(-1:nx,-1:ny,0:1)
      integer b2stbc_fchy_dia_coreonly, secmodel, area_fix,
     &  sheath_drift_fix
      logical wrong_flow                                             !xpb
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
     1  bc_type13_ref, bc_type13_ix, bc_type13_iy, bc_type13_ixx,
     1  bc_type16_ref, bc_type16_ix, bc_type16_iy, bc_type16_i,
     2  bc_type16_ixx,
     1  bc_type20_ref, bc_type20_ix, bc_type20_iy, bc_type20_ixx,
     1  bc_type21_ref, bc_type21_ix, bc_type21_iy, bc_type21_ixx,
     1  bc_type16_kinetic_energy
      integer       :: num_inner, k                                  !srv 26.11.02
      integer, save :: istyle_average_pr = 0, istyle_fchi = 0        !srv 08.08.02 14.02.05
      integer, save :: iform_corr_dia = 0                            !srv 25.07.05
      integer, save :: mdf_fnb = 0, mdf_fhe = 0, mdf_fhi = 0,
     &                 iout = 0,                                     !srv 13.10.06
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2npc7 = 0, iout_b2npc11 = 0,            !srv 02.04.18
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0,
     &                 iout_b2npmo = 0, iout_b2wdat = 0
      integer, save :: istyle_cur_contr_on_S_and_N = 2               !srv 01.02.09 03.03.15
      integer, save :: bcene_16_style = 0, bceni_16_style = 0        !srv 16.02.12
      integer, save :: bcene_15_style = 0                            !srv 19.04.13
      integer, save :: pot_eq = 1                                    !srv 27.02.13
      integer, save :: use_astra
      integer, save :: drift_style
      logical once_bccon(0:ns-1), once_bcene, once_bceni, once_bcpot !srv 16.02.12 12.02.15 29.09.15
      real (kind=R8) , save :: cbc = 1.0e0_R8, bc_ref = 0.01         !srv 14.10.99
      real (kind=R8) , save :: bc_ref_te = 0.01, bc_ref_ti=0.01      !IYS 07.03.2013 11.04.13
      real (kind=R8) , save :: senepar, senipar,                     !srv 16.02.12 29.09.15
     &                         sconpar(0:DEF_NSD-1), scurpar
      real (kind=R8), save :: prev_cur_delta, cur_delta, po_step,    !Solovyev 01.05.14
     * po_prev, po_curr, coeff_16                                    !Solovyev 01.05.14, 24.05.14
      integer, save :: bcpot_16_step = 50                            !Solovyev 23.05.14 !srv 29.09.15 }
      integer, save :: eirene_sheath_pot = 1
      real (kind=R8) ::
     *  vte, s0, s1, s0hz, s1hz, fche, fchi, t0, t1, t2,             !srv 15.11.07
     1  t3, totflux, weight,
     1  naave, volave, teave, neave, tiave, niave, poave,
     2  vwest, veast, totfluxn, totfluxi
     &  ,pzaverage, pzsum, volsum                                    !srv 24.03.99 18.01.02
     &  ,ua2sum, ua2average(0:ns-1)                                  !srv 07.04.00 15.02.02
     &  ,uasum, uaaverage(0:ns-1)                                    !srv 15.02.02
     &  ,nesum, tisum, tiaverage, tesum                              !srv 15.12.00 25.01.02
     &  ,nasum, naaverage(0:ns-1), bbsum, bbaverage                  !srv 22.01.02 18.10.04
      real (kind=R8) :: bzhzsum, bzhzaverage                         !lk 20.11.07
      real (kind=R8) ::
     *  pr(-1:nx,-1:ny), wdia(-1:nx,-1:ny), wrk0(-1:nx,-1:ny),
     *  uu(-1:nx,-1:ny,0:ns+ns_ext-1),
     *  sx(-1:nx,-1:ny), ud(-1:nx,-1:ny,0:ns-1),
     *  sxhz(-1:nx,-1:ny), pz(-1:nx,-1:ny), rz(-1:nx,-1:ny),         !srv 24.08.08
     &  sna0_no_mdf(-1:nx,-1:ny,0:1,0:ns-1)                          !srv 23.09.08
      real (kind=R8), save ::                                        !srv 16.02.12
     *  NUEOUT, UBIG, VBIG, PBIG, UP00, EPSLON, NE00, S2, US, FS, VS,
     *  NAS, TES, TIS, NIVS, FFAC, DFS, POS                          !srv 29.09.15
      real (kind=R8), save :: ycur = 1.0e0_R8
      real (kind=R8), save :: vis_per = 1.0_R8, vis_q = 1.0_R8       !srv 02.10.17
      real (kind=R8), save :: phm0 = 1.0e0_R8, phm1 = 1.0e0_R8,
     &                        phm2 = 1.0e0_R8,
     &                        corr_flux   = 0.0e0_R8                 !srv 25.07.05
      real (kind=R8), save :: alfTeEh = 0.0_R8
      character chns*3, chk*1                                        !srv 11.07.05
      logical decay_length_ok, fedback_species, ishigh
      logical, save :: BCCON14(NBCD)                                 !srv 20.09.17
      logical, save :: consistent
      real (kind=R8) ::
     *  ro(-1:nx,-1:ny), cs, csb, uuout, upout, zbar, zarg, uutest,  !srv 01.02.09
     *  fchy_dia, fchy_inertia, seec                                 !srv 07.04.00
      integer isfeedback, fix_fch_in_fhe_sheath
      integer, save :: fix_sx_in_fhe_sheath = 1                      !srv 08.05.09
      integer, save :: use_eirene, no_current, ionising_core
      real (kind=R8) :: nesepm, nesepm_overshoot, ndes_sol,
     1 bc_type13_norm, bc_type13_fac, na_min, na_new
      character*26 boundary_location(nbcd)
      character*24 hlp_frm
      save ncall, coreregno, coreregn2,  bc_type13_ref,
     1  bc_type16_ref, bc_type20_ref, bc_type13_norm, bc_type13_fac,
     2  na_min, na_new, bc_type21_ref, ndes_sol,
     *  isfeedback, nesepm, nesepm_overshoot, fix_fch_in_fhe_sheath
*   ..procedures
      intrinsic abs, sqrt, mod
      real (kind=R8) ::
     *  bcint
      external xertst, ipgeti, bcint, eseec0
      external ipgetr, b2xvsg, b2xvff, b2xvps
*   ..initialisation
      data ncall/0/, coreregno/1/, coreregn2/4/, bc_type13_ref/1/,
     1 bc_type16_ref/1/, bc_type20_ref/1/, ionising_core/0/,
     2 na_min/1.0e4_R8/, na_new/1.0e14_R8/,
     * isfeedback/0/, nesepm/0.0_R8/, nesepm_overshoot/0.0_R8/,
     5 fix_fch_in_fhe_sheath/2/, bc_type16_kinetic_energy/0/,
     & bc_type13_norm/1.0e15_R8/, bc_type13_fac/1.0_R8/,
     & bc_type21_ref/1/, no_current/0/,
     & use_astra/0/, use_eirene/0/, ndes_sol/0.0_R8/
      real (kind=R8) ::
     .  Tguess, fy, ehy, fylimt, prefactor, w1,
     .  fyt, gam0, gam1, tlimit, cone0, coni0, flo, a,
     .  cum_flo_vadia, cum_flo_vaecrb, flo_vedia, flo_veecrb
      integer ibnd, conduction_only, no_hybr_e, no_hybr_i
      data conduction_only/0/, no_hybr_e/0/, no_hybr_i/0/
      save conduction_only, no_hybr_e, no_hybr_i
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbc_phys')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2stbc_phm0', phm0)                             !srv 17.01.07
       call ipgetr ('b2stbc_phm1', phm1)                             !srv 17.01.07
       call ipgetr ('b2stbc_phm2', phm2)                             !srv 17.01.07
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
       call ipgeti ('b2tfhe_conduction_only', conduction_only)
       call ipgeti ('b2tfhe_no_hybr', no_hybr_e)
       call ipgeti ('b2tfhi_no_hybr', no_hybr_i)
       call ipgeti ('b2mndr_eirene', use_eirene)
       call ipgeti ('b2mndr_astra', use_astra)
       if (redef_gmtry.eq.0) then
         drift_style = 1
       else
         drift_style = 2
       end if
       call ipgeti ('b2tfnb_drift_style', drift_style)
       call ipgeti ('b2stbc_type13_ref', bc_type13_ref)
       call ipgeti ('b2stbc_type16_ref', bc_type16_ref)
       call ipgeti ('b2stbc_type20_ref', bc_type20_ref)
       call ipgeti ('b2stbc_type21_ref', bc_type21_ref)
       call ipgetr ('b2stbc_type13_norm', bc_type13_norm)
       call ipgetr ('b2stbc_type13_fac', bc_type13_fac)
       call ipgetr ('b2mndr_na_min', na_min)
       call ipgetr ('b2mndr_na_new', na_new)
       call xertst (0.0_R8.lt.na_min, 'faulty parameter na_min')
       call xertst (na_new.gt.0.0_R8, 'faulty parameter na_new')
       call xertst (na_new.gt.na_min,
     &  'na_new must be larger than na_min!')
       call ipgeti ('b2stbc_fix_fch_in_fhe_sheath',
     1   fix_fch_in_fhe_sheath)
       call ipgeti ('b2stbc_fix_sx_in_fhe_sheath',
     1   fix_sx_in_fhe_sheath)                                       !srv 08.05.09
       call ipgeti ('b2stbc_type16_kinetic_energy',
     1   bc_type16_kinetic_energy)
       call ipgetr ('b2stbc_cbc', cbc)                               !srv 14.10.99
       call ipgeti ('b2stbc_average_pr', istyle_average_pr)          !srv 08.08.02
       call ipgeti ('b2stbc_istyle_fchi', istyle_fchi)               !srv 14.02.05
       call ipgetr ('b2stbc_corr_flux_in_core', corr_flux)           !srv 25.07.05
       call ipgeti ('b2stbc_iform_corr_dia', iform_corr_dia)         !srv 25.07.05
       call ipgeti ('b2stbc_iout', iout)                             !srv 11.07.05
       call ipgeti ('b2tfnb_mdf', mdf_fnb)                           !srv 13.10.06
       call ipgeti ('b2tfhe_mdf', mdf_fhe)                           !srv 13.10.06
       call ipgeti ('b2tfhi_mdf', mdf_fhi)                           !srv 13.10.06
       call ipgetr ('b2tfhe_alfTeEh', alfTeEh)                       !srv 28.04.08
       call ipgeti ('b2stbc_istyle_cur_contr_on_S_and_N',
     &                      istyle_cur_contr_on_S_and_N)             !srv 01.02.09
       call ipgeti ('b2wdat_iout', iout_b2wdat)                      !iys 17.04.17
       call ipgeti ('b2npc9_iout', iout_b2npc9)                      !srv 11.09.09 {
       call ipgeti ('b2npco_iout', iout_b2npco)
       call ipgeti ('b2npc7_iout', iout_b2npc7)                      !srv 02.04.18
       call ipgeti ('b2npc11_iout', iout_b2npc11)                    !srv 02.04.18
       call ipgeti ('b2npmo_iout', iout_b2npmo)
       call ipgeti ('b2nph9_iout', iout_b2nph9)
       call ipgeti ('b2npht_iout', iout_b2npht)
       call ipgeti ('b2npp7_iout', iout_b2npp7)
       call ipgeti ('b2nppo_iout', iout_b2nppo)                      !srv 11.09.09 }
       call ipgetr ('b2stbc_bc_ref', bc_ref)
       call ipgetr ('b2stbc_bc_ref_te', bc_ref_te)                   !srv 11.04.13
       call ipgetr ('b2stbc_bc_ref_ti', bc_ref_ti)                   !srv 11.04.13
       call ipgeti ('b2stbc_bcene_15_style', bcene_15_style)
       call ipgeti ('b2stbc_bcene_16_style', bcene_16_style)
       call ipgeti ('b2stbc_bceni_16_style', bceni_16_style)
       call ipgeti ('b2tfhe_no_current', no_current)
       call ipgeti ('b2news_poteq', pot_eq)                          !srv 27.02.13
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &   'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
       call ipgeti ('b2stbc_bcpot_16_step', bcpot_16_step)           !Solovyev 23.05.14 !srv 29.09.15
       call ipgetr ('b2tfhi_ycur', ycur)
       call ipgetr ('b2tfnb_vis_per', vis_per)                       !srv 02.10.17
       call ipgetr ('b2tfnb_vis_q', vis_q)                           !srv 02.10.17
       call ipgeti ('eirene_sheath_pot',eirene_sheath_pot)
       call ipgeti ('eirene_ionising_core', ionising_core)
       if (ionising_core.eq.0)     ! for backward compatibility
     >   call ipgeti ('eirene_ionizing_core', ionising_core)
       call xertst (ionising_core.ge.-nbc.and.ionising_core.le.1,
     1   'faulty input ionising_core, abs(ionising_core) > NBC')
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
      if (.not.b2mod_math_initialised) call b2mod_math_init

*  -- INITIALISE CERTAIN CONSTANTS --
      EPSLON=1.e-10_R8
*  -- SET NE00, UP00, UBIG, VBIG --
      NE00=1.0e20_R8
      UP00=SQRT(100*EV/MP)
      UBIG=UP00 / EPSLON
      VBIG=1.0e-2_R8*UBIG
      PBIG=UP00

      sna0_no_mdf = 0.0_R8

      call b2xppz (nx, ny, ns, ne, na, te, ti, pz)                !srv 18.01.02
*   ..find the average pressure for inner flux surface            !srv 24.03.99
      if (istyle_average_pr .eq.0) then                           !srv 08.08.02
       iy = -1                                                    !srv 24.03.99
       pzsum =0.0e0_R8                                            !srv 24.03.99 18.01.02
       volsum=0.0e0_R8                                            !srv 24.03.99
       do ix = -1, nx                                             !srv 26.11.02
        if (on_closed_surface(ix,iy)) then                        !srv 26.11.02
         pzsum=pzsum + pz(ix,iy)*vol(ix,iy)                       !srv 24.03.99 18.01.02
         volsum=volsum + vol(ix,iy)                               !srv 24.03.99
        endif                                                     !srv 26.11.02
       enddo                                                      !srv 24.03.99
       if (volsum.gt.0.0_R8) pzaverage=pzsum/volsum               !srv 24.03.99 18.01.02
      else                                                        !srv 08.08.02
       iy = 0                                                     !srv 24.03.99  07.06.02 {
       pzsum =0.0e0_R8                                            !srv 24.03.99 18.01.02
       volsum=0.0e0_R8                                            !srv 24.03.99
       num_inner=0                                                !srv 26.11.02
       do ix = -1, nx                                             !srv 26.11.02
        if (on_closed_surface(ix,iy)) then                        !srv 26.11.02
         pzsum=pzsum
     *   + (pz(ix,iy)+pz(bottomix(ix,iy),bottomiy(ix,iy)))        !srv 24.03.99 18.01.02
         volsum=volsum + vol(ix,iy)                               !srv 24.03.99
         num_inner=num_inner+1                                    !srv 26.11.02
        endif                                                     !srv 26.11.02
       enddo                                                      !srv 24.03.99
       if (num_inner.gt.0) pzaverage=pzsum*0.5e0_R8/num_inner     !srv 24.03.99 18.01.02 07.06.02 } 26.11.02
      endif                                                       !srv 08.08.02
*   ..find the average ne, te, ti, bb for inner flux surface      !srv 15.12.00 18.10.04
      bbsum = 0.0e0_R8                                            !srv 18.10.04
      nesum = 0.0e0_R8                                            !srv 15.12.00 25.01.02
      tesum = 0.0e0_R8                                            !srv 15.12.00
      tisum = 0.0e0_R8                                            !srv 15.12.00
      bzhzsum  = 0.0e0_R8                                         !lk 20.11.07
      bbaverage  = 0.0e0_R8                                       !lk 20.11.07
      do ix = -1, nx                                              !srv 26.11.02
       if (on_closed_surface(ix,iy)) then                         !srv 26.11.02
        bbsum=bbsum + bb(ix,iy,3)*vol(ix,iy)                      !srv 18.10.04
        nesum=nesum + ne(ix,iy)*vol(ix,iy)                        !srv 15.12.00 25.01.02
        tisum=tisum + ti(ix,iy)*vol(ix,iy)                        !srv 15.12.00
        tesum=tesum + te(ix,iy)*vol(ix,iy)                        !srv 15.12.00
        bzhzsum = bzhzsum+bb(ix,iy,2)*vol(ix,iy)**2/gs(ix,iy,2)   !lk 20.11.07 !srv 07.05.13
       endif                                                      !srv 26.11.02
      enddo                                                       !srv 15.12.00
      if (volsum.gt.0.0_R8) then
       bbaverage=bbsum/volsum                                     !srv 18.10.04
       tiaverage=tisum/volsum                                     !srv 15.12.00
       bzhzaverage=bzhzsum/volsum                                 !lk 20.11.07
      end if
*   ..find the average ua**2, na for inner flux surface           !srv 22.01.02 {
      do is = 0, ns-1
       uasum  = 0.0e0_R8                                           !srv 15.02.02
       ua2sum = 0.0e0_R8                                           !srv 07.04.00 15.02.02
       nasum  = 0.0e0_R8
       do ix = -1, nx                                              !srv 26.11.02
        if (on_closed_surface(ix,iy)) then                         !srv 26.11.02
         ua2sum = ua2sum + ua(ix,iy,is)**2*vol(ix,iy)              !srv 07.04.00 15.02.02
         uasum  = uasum  + ua(ix,iy,is)*vol(ix,iy)                 !srv 15.02.02
         nasum  = nasum  + na(ix,iy,is)*vol(ix,iy)
        endif                                                      !srv 26.11.02
       enddo
       if (volsum.gt.0.0_R8) then
        ua2average(is) = ua2sum/volsum                             !srv 07.04.00 15.02.02
        uaaverage(is) = uasum/volsum                               !srv 07.04.00 15.02.02 16.04.02
        naaverage(is) = nasum/volsum
       end if
      enddo                                                        !srv 22.01.02 }

      call b2xpro (nx, ny, ns, mp, am, na, ro)                     !srv 01.02.09
      call b2xprz (nx, ny, ns, mp, am, na, rz)
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
            sxhz(ix,iy)=gs(rightix(ix,iy),rightiy(ix,iy),0)*       !srv 24.08.08
     &                  hz(ix,iy)
!     &                  hz(rightix(ix,iy),rightiy(ix,iy))
            if (area_fix.ge.3)
     >       sxhz(ix,iy) = sxhz(ix,iy)*qc(rightix(ix,iy),rightiy(ix,iy))
          else
            sx(ix,iy)=0.0_R8
            sxhz(ix,iy)=0.0_R8                                     !srv 24.08.08
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
cwdk  With redefinition of drifts as being the poloidal projections,
cwdk  now need to recompute full diamagnetic component for ud
      if (isymm.eq.2.or.isymm.eq.4) then
        ud(:,:,:) = 0.0_R8
      else if (sheath_drift_fix.eq.0) then
        do is = 0,ns-1
          ud(:,:,is) = uadia(:,:,0,is)*bb(:,:,3)/bb(:,:,2)*signmf
        enddo
      else
        do is=0,ns-1
          ud(:,:,is) = vaecrb(:,:,0,is)*bb(:,:,3)/bb(:,:,2)*signmf
        enddo
      endif
      if (main_call .and. diagno.ge.1) then                           !xpb
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
      endif                                                           !xpb

      if(ncall.eq.0) write(*,*) 'Doing B2-style bc'

      once_bccon = .false.
      once_bcene = .false.
      once_bceni = .false.
      once_bcpot = .false.

* -- LOOP OVER ALL BOUNDARIES
      DO IB = 1, NBC

        if (bcpos(ib).eq.-2) then
          write(boundary_location(ib),'(a23,i2)')
     .     ' on boundary cell list ',IB
        else
          write(boundary_location(ib),'(a4,i4,a6,i4,a4,i4)')
     .     ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
        endif

* -- Moving core tally according to B.C.
* -- Must be done only once!

        IF (ncall.eq.0) THEN
          IF (BCCHAR(IB).EQ.'S' .AND. bc_type16_ref.gt.1 .AND.
     &       (BCENE(IB).EQ.16.or.BCENE(IB).EQ.17.or.BCENE(IB).EQ.20.or.
     &        BCENI(IB).EQ.16.or.BCENI(IB).EQ.17.or.BCENI(IB).EQ.20))
     &     THEN
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

* -- Consistency checks
          consistent = .true.
          bccon14(ib) = .true.
          IF (BCENE(IB).EQ.15 .OR. BCENI(IB).EQ.15 .OR. BCPOT(IB).EQ.11)
     .     THEN
            consistent = BCENE(IB).eq.15 .and. BCENI(IB).eq.15 .and.
     &                   BCPOT(IB).eq.11
            DO IS = 0, NS-1
              IF (IS_NEUTRAL(IS)) CYCLE
              consistent = consistent .and. BCMOM(IS,IB).EQ.13
              consistent = consistent .and. BCCON(IS,IB).EQ.14
              bccon14(ib) = bccon14(ib) .and. BCCON(IS,IB).EQ.14
            END DO
          END IF
          IF (.not.consistent) then
            write(*,*) 'DRIFT SHEATH BOUNDARY CONDITIONS AT BOUNDARY ',
     &       IB,' MAY NOT BE CONSISTENT!'
            WRITE(*,*) 'IT IS RECOMMENDED TO IMPOSE TOGETHER : '
            WRITE(*,*) 'BCENE/I=15, BCPOT=11, ',
     &       'with BCMOM=13 and BCCON=14 for all ion species'
            write(0,*) 'DRIFT SHEATH BOUNDARY CONDITIONS AT BOUNDARY ',
     &       IB,' MAY NOT BE CONSISTENT!'
            WRITE(0,*) 'IT IS RECOMMENDED TO IMPOSE TOGETHER : '
            WRITE(0,*) 'BCENE/I=15, BCPOT=11, ',
     &       'with BCMOM=13 and BCCON=14 for all ion species'
          END IF
        ENDIF

        call b2stbc_phys_species(IB)
        call b2stbc_phys_ene(IB)
        call b2stbc_phys_eni(IB)
        call b2stbc_phys_pot(IB)
      END DO

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

      if (iout_b2npmo.eq.1 .or. iout_b2wdat.eq.4) then               !srv 11.09.09 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        wrk0=smo0(:,:,0,is)+
     &       smo0(:,:,1,is)*ua(:,:,is)+
     &       smo0(:,:,2,is)*na(:,:,is)*mp*am(is)+
     &       smo0(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is)
        call my_out(70,nx,ny,wrk0,'b2stbc_phys_smo'//chns)
       enddo
      endif
      if (iout.eq.1 .or. (iout_b2npc9.eq.1 .or. iout_b2npco.eq.1
     &              .or.  iout_b2npc7.eq.1 .or. iout_b2npc11.eq.1 )  !srv 02.04.18
     &        .or. iout_b2wdat.eq.4 ) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
!        do k=0,1
!         write (chk,'(i1)') k
!         call my_out(70,nx,ny,sna0(-1,-1,k,is),
!     &               'b2stbc_phys_sna'//chk//chns)
!        enddo
        wrk0=sna0(:,:,0,is)+sna0(:,:,1,is)*na(:,:,is)
        call my_out(70,nx,ny,wrk0,'b2stbc_phys_sna'//chns)
       enddo
      endif
      if ((iout_b2npht.eq.1 .or. iout_b2nph9.eq.1)
     &        .or. iout_b2wdat.eq.4 ) then
       wrk0=she0(:,:,0)+she0(:,:,1)*te+
     &      she0(:,:,2)*ne + she0(:,:,3)*te*ne
       call my_out(70,nx,ny,wrk0,'b2stbc_phys_she')
       wrk0=shi0(:,:,0)+shi0(:,:,1)*ti+
     &      shi0(:,:,2)*ni(:,:,0)+shi0(:,:,3)*ti*ni(:,:,0)
       call my_out(70,nx,ny,wrk0,'b2stbc_phys_shi')
      endif
      if ((iout_b2npp7.eq.1 .or. iout_b2nppo.eq.1)
     &        .or. iout_b2wdat.eq.4 ) then
       wrk0=sch0(:,:,0)+sch0(:,:,1)*po+
     &      sch0(:,:,2)*ne+sch0(:,:,3)*po*ne(:,:)
       call my_out(70,nx,ny,wrk0,'b2stbc_phys_sch')
      endif

      if (iout.eq.1) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        do k=0,3
          write (chk,'(i1)') k
          call my_out(70,nx,ny,smo0(-1,-1,k,is),
     &               'b2stbc_phys_smo'//chk//chns)
        enddo
        call my_out(70,nx,ny,ua(:,:,is),'b2stbc_phys_ua'//chns)
        call my_out(70,nx,ny,na(:,:,is),'b2stbc_phys_na'//chns)
        do k=0,1
         write (chk,'(i1)') k
         call my_out(70,nx,ny,sna0(-1,-1,k,is),
     &               'b2stbc_phys_sna'//chk//chns)
        enddo
       enddo
       do k=0,3
         write (chk,'(i1)') k
         call my_out(70,nx,ny,sch0(-1,-1,k),'b2stbc_phys_sch'//chk)
         call my_out(70,nx,ny,she0(-1,-1,k),'b2stbc_phys_she'//chk)
         call my_out(70,nx,ny,shi0(-1,-1,k),'b2stbc_phys_shi'//chk)
       enddo
       call my_out(70,nx,ny,te,'b2stbc_phys_te')
       call my_out(70,nx,ny,ti,'b2stbc_phys_ti')
       call my_out(70,nx,ny,po,'b2stbc_phys_po')
      endif                                                           !srv 11.09.09 }

* ..return
      if(ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbc_phys

      contains

      subroutine b2stbc_phys_species(IB)
      integer ib
      integer ix, iy, is, isa, isb, ibw
#ifdef B25_EIRENE
      integer k, ism
#endif
      integer bc_type13_i, bc_type20_i, bc_type21_i

      if (bcpos(ib).eq.-2) then
        write(boundary_location(ib),'(a23,i2)')
     .   ' on boundary cell list ',IB
      else
        write(boundary_location(ib),'(a4,i4,a6,i4,a4,i4)')
     .   ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
      endif

      DO IS=0, NS-1
        if (is.eq.ns-1) then
          ishigh = .true.
        else
          ishigh = .not.lnext(is,is+1)
        end if

* -- BOUNDARY CONDITIONS FOR CONTINUITY EQUATIONS --
        select case (BCCON(IS,IB))
        case (0)
          continue

* -- BCCON=1 -- PRESCRIBE THE VALUE OF THE DENSITY

        case (1)
          if(ncall.eq.0) then
            IF(core_boundary(ib)) THEN
             if(is.eq.ismain.and.use_astra.ne.0) then
              CONPAR(IS,IB,1) = neb_astra                               !srv 28.07.11
             endif
             if (ishigh) call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 1')
            ENDIF
            if (CONPAR(IS,IB,1).eq.0.0_R8) then
              CONPAR(IS,IB,1) = na_min
              write(*,*) 'Zero density B.C. no longer allowed!'
              write(*,*) 'Changed to na_min.'
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

* -- BCCON=2 -- PRESCRIBE THE GRADIENT OF THE DENSITY

        case (2)
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'density gradient ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            if (core_boundary(ib).and.ishigh)
     >       call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 2')
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

* -- BCCON=3 -- SHEATH CONDITIONS

        case (3)
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

* -- BCCON=4 -- PRESCRIBE THE VALUE OF THE DENSITY, WEAKLY
*                 A MIXED BOUNDARY CONDITION

        case (4)
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'weak density ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1       'BCCON = 4, CONPAR(,,1) !> 0')
            call xertst(CONPAR(IS,IB,2).gt.0.0_R8,
     1       'BCCON = 4, CONPAR(,,2) !> 0')
            if (core_boundary(ib).and.ishigh)
     .       call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 4')
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

* -- BCCON=5 -- PRESCRIBE THE PARTICLE FLUX PER UNIT AREA

        case (5)
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

* -- BCCON=6 -- PRESCRIBE THE TOTAL PARTICLE FLUX FOR A
*                 CONSTANT DENSITY

        case (6)
          if(ncall.eq.0) then
           write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .      'total particle flux ',CONPAR(IS,IB,1),
     1      ' with const density on ',BCCHAR(IB),
     2      boundary_location(ib),' for species ',IS
           if (core_boundary(ib).and.ishigh)
     >      call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 6')
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

* -- BCCON=7 -- PRESCRIBE THE DENSITY AS A FUNCTION OF OTHER
*                 PLASMA PARAMETERS

        case (7)
          CALL XERRAB ('B2STBC -- BCCON=7, NOT YET AVAILABLE')

* -- BCCON=8 -- PRESCRIBE THE TOTAL PARTICLE FLUX WITH
*                 CONSTANT FLUX DENSITY

        case (8)
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
     &            sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)         !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &            sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &            sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)         !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) -
     &           fnaPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                   topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
                t0 =
     &            vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &            hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &            vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                   topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) +
     &           max(0.0_R8,t0)*na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) +
     &           min(0.0_R8,t0)
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
     &          sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)           !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &          sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &          sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)           !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) +
     &           fnaPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
                t0 =
     &           vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &           vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) +
     &           max(0.0_R8,t0)*na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) +
     &           min(0.0_R8,t0)
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
     &         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)            !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     &         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is)            !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     &           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) -
     &           fnaPSch(
     &            rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &            rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)
                t0 =
     &           vol(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &           hx(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     &           vaecrb(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                  rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     &           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) +
     &           max(0.0_R8,t0)*na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     &           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) +
     &           min(0.0_R8,t0)
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
     .         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)            !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is)            !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) +
     &           fnaPSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)
                t0 =
     &           vol(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &           hx(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     &           vaecrb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     &           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) +
     &           max(0.0_R8,t0)*na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     &           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) +
     &           min(0.0_R8,t0)
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
     &          sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)           !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &          sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &          sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)           !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) +
     &          (fnaPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &           fnaPSch(
     &            topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &            topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is))/2.0_R8
                if (drift_style.eq.0) then
                  t0 =
     &              vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &              hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &              vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
                else
                  t0 =
     &              vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &              hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*0.5_R8*
     &             (vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &              vaecrb(
     &               topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &               topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is))
                endif
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) +
     &           max(0.0_R8,t0)*na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) +
     &           min(0.0_R8,t0)
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
     .         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)            !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is)            !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) +
     &          (fnaPSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)+
     &           fnaPSch(
     &            rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &            rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is))/
     &           2.0_R8
                if (drift_style.eq.0) then
                  t0 =
     &             vol(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &             hx(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     &             vaecrb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)
                else
                  t0 =
     &             vol(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &             hx(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*0.5_R8*
     &            (vaecrb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)+
     &             vaecrb(
     &              rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &              rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is))
                endif                                                  !srv 12.09.08 }
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     &           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) +
     &           max(0.0_R8,t0)*na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     &           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) +
     &           min(0.0_R8,t0)
              endif
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCCON=8, NOT AVAILABLE')
          ENDIF

* -- BCCON=9 -- PRESCRIBE THE DECAY LENGTH FOR THE DENSITY

        case (9)
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'density decay length ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).ne.0.0_R8,
     1       'B2STBC -- BCCON=9, CONPAR(,,1) = 0')
            if (core_boundary(ib).and.ishigh)
     >       call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 9')
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
     .         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)            !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)            !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 08.09.08 {
                t0 =
     .           -fnaPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
                if (t0.ge.0.0_R8) then
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     .             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) + t0
                else
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) + t0/
     &               na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                end if
                t1 =
     .           vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
                if (t1.le.0.0_R8) then
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     .             vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .             hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*t1
                else
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     .             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+
     .             vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .             hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*t1*
     .              na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                end if
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
     .         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)            !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)            !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 08.09.08 {
                t0 = fnaPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
                if (t0.ge.0.0_R8) then
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) + t0
                else
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) + t0/
     &               na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                end if
                t1 = -vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
                if (t1.le.0.0_R8) then
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &             vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &             hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*t1
                else
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+
     &             vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &             hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*t1*
     &              na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                end if
              endif                                                    !srv 08.09.08 }
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCCON=9, NOT AVAILABLE FOR W,E')
          ENDIF

* -- BCCON=10 -- leakage option for density, recommended for cases with drifts

        case (10)
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'leakage ',CONPAR(IS,IB,1),' on ',BCCHAR(IB),
     1       boundary_location(ib),
     2       ' for species ',IS
            if (core_boundary(ib).and.ishigh)
     >       call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 10')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              s0 = vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .             hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))              !srv 26.02.09
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=0.0_R8
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=s0*
     2         sqrt((rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .                te(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .                ti(bc_list_x(ix,ib),bc_list_y(ix,ib)))/
     3         (am(is)*mp))*
     2         CONPAR(IS,IB,1)
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,is) =
     &         sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,is)+
     &         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,is)            !srv 23.09.08
*    ..drifts to the wall                                              !srv 02.03.05 {
              if(mdf_fnb.eq.0) then                                    !srv 09.11.06 08.09.08 {
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     *           min(0.0e0_R8,s0*(
     .            vadia(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)+
     &           vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)))
              else
                t0 =
     .          -fnaPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                   topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
                if (t0.ge.0.0_R8) then
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     .             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) + t0
                else
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) + t0/
     .             na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                end if
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &            s0*min(0.0e0_R8,
     &           vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)) !srv 02.03.05 }
              endif                                                     !srv 09.11.06 08.09.08 }
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              s0=vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))                !srv 26.02.09
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=0.0_R8
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=s0*
     2         sqrt((rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .                te(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .                ti(bc_list_x(ix,ib),bc_list_y(ix,ib)))/
     3         (am(is)*mp))*
     2         CONPAR(IS,IB,1)
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,is) =
     &         sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,is)+
     &         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,is)            !srv 23.09.08
*    ..drifts to the wall                                              !srv 02.03.05 {
              if(mdf_fnb.eq.0) then                                    !srv 09.11.06 08.09.08 {
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)-
     *           max(0.0e0_R8,
     .               s0*(vadia(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     *                  vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)))
              else
                t0 = fnaPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
                if (t0.ge.0.0_R8) then
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     .             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) + t0
                else
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) + t0/
     .             na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                end if
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)-
     &           max(0.0e0_R8,
     .               s0*vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)) !lk 12.12.07 }
              endif                                                     !srv 09.11.06 08.09.08 }
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
     &         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,is)            !srv 23.09.08
*    ..drifts to the wall                                              !srv 02.03.05 {
              if(mdf_fnb.eq.0) then                                    !srv 09.11.06 08.09.08 {
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is)-
     *           max(0.0e0_R8,s0*(vadia(
     .            rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)+
     *     vaecrb(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)))
              else
                t0 = -fnaPSch(
     &            rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)
                if (t0.ge.0.0_R8) then
                  sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     .             sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) + t0
                else
                  sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .             sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) + t0/
     .             na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
                end if
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is)-
     &           max(0.0e0_R8,s0*vaecrb(
     .               rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)) !lk 12.12.07 }
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
     &         sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,is)            !srv 23.09.08
*    ..drifts to the wall                                              !srv 02.03.05 {
              if(mdf_fnb.eq.0) then                                    !srv 09.11.06 08.09.08 {
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is)-
     *           max(0.0e0_R8,
     .               s0*(vadia(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)+
     *                  vaecrb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)))
              else
                t0 = fnaPSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)
                if (t0.ge.0.0_R8) then
                  sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     .             sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) + t0
                else
                  sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .             sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) + t0/
     .             na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
                endif
                sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     .           sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is)-
     &           max(0.0e0_R8,
     .               s0*vaecrb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)) !lk 12.12.07 }
              endif                                                     !srv 09.11.06 08.09.08 }
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCCON=10, NOT AVAILABLE')
          ENDIF

* -- BCCON=11 -- particle flux feedback boundary condition

        case (11)
          if(ncall.eq.0) then
           write(*,'(a,a,a,a,i3)')
     .      'particle flux feedback boundary condition on ',
     1      BCCHAR(IB),boundary_location(ib),' for species ',IS
           call ipgeti ('b2stbc_isfeedback', isfeedback)
           call xertst (0.le.isfeedback.and.isfeedback.lt.ns,
     1      'faulty internal parameter isfeedback')
           fedback_species = is.eq.isfeedback
           do isa = 0, nspecies-1
             if (na_feedback_ib(isa).eq.ib) then
               fedback_species = fedback_species .or.
     .                          (is.ge.eb2spcr(isa+1) .and.
     .                           is.le.eb2spcr(isa+1)+nfluids(isa+1))
             end if
           end do
           call xertst(fedback_species,
     1       'BCCON(,IS,1) = 11 for non-feedback species!')
           call xertst(lfeedback, 'LFEEDBACK must be turned on!')
          endif
          fedback_species = .false.
          do isa = 0, nspecies-1
            if (na_feedback_target(isa).gt.0.0_R8 .and.
     &          na_feedback_ib(isa).eq.ib .and.
     &          is.ge.eb2spcr(isa+1)+1 .and.
     .          is.le.eb2spcr(isa+1)+nfluids(isa+1)) then
              t0 = conpar(is,ib,1)
              fedback_species = .true.
            end if
          end do
          IF (BCCHAR(IB).EQ.'S') THEN
            if (.not.fedback_species) then
*   ..initialise ireg
              ireg = cbirso
*   ..update ireg
 1101         if (cbrbrk(ireg)*nx.lt.
     .          maxval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8) then
                call xertst (ireg.lt.cbirso+cbnrso-1,
     &           'error in computing region on south boundary')
                ireg = ireg+1
                goto 1101
              endif
              if (ireg.gt.cbirso) then
                call xertst (cbrbrk(ireg-1)*nx.lt.
     &           minval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8,
     &           'feedback bc covers more than one region!')
              endif
              t0 = cbsna(0,IS,ireg)
            end if
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         t0*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            ENDDO
          ELSEIF (BCCHAR(IB).EQ.'N') THEN
            if (.not.fedback_species) then
*   ..initialise ireg
              ireg = cbirno
*   ..update ireg
 1102         if (cbrbrk(ireg)*nx.lt.
     .          maxval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8) then
                call xertst (ireg.lt.cbirno+cbnrno-1,
     &           'error in computing region on north boundary')
                ireg = ireg+1
                goto 1102
              endif
              if (ireg.gt.cbirno) then
                call xertst (cbrbrk(ireg-1)*nx.lt.
     &           minval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8,
     &           'feedback bc covers more than one region!')
              endif
              t0 = cbsna(0,IS,ireg)
            end if
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         t0*SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .               BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            ENDDO
          ELSE
            CALL XERRAB
     .       ('B2STBC -- BCCON=11, NOT AVAILABLE FOR W,E')
          ENDIF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                       !srv 23.09.08 }

* -- BCCON=12 -- particle density feedback boundary condition

        case (12)
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
           call ipgetr ('b2stbc_ndes_sol', ndes_sol)
           call ipgetr ('b2stbc_nesepm_overshoot', nesepm_overshoot)
           call xertst ((0.0_R8.eq.nesepm.and.0.0_R8.eq.ndes_sol).or.
     1      0.0_R8.eq.nesepm_overshoot,
     1      'BCCON=12 coding incompatible with nesepm_overshoot!')
           if (core_boundary(ib).and.ishigh)
     >      call xertst (ionising_core.eq.0,
     .      'ionising_core switch usage inconsistent with BCCON = 12')
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
            t0=cbsna(0,IS,ireg)
            call xertst (t0.gt.0.0_R8,
     &       'BCCON=12, trying to set a negative density!')
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*t0*
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
            t0=cbsna(0,IS,ireg)
            call xertst (t0.gt.0.0_R8,
     &       'BCCON=12, trying to set a negative density!')
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*t0*
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
          ENDDO                                                       !srv 23.09.08 }

* -- BCCON=13 -- Feedback constant density for a given particle flux

        case (13)
          if(ncall.eq.0) then
            write(*,'(a,a,a,a,i3,a,1p,2g14.7)')
     .       'particle density to achieve specified total flux on ',
     1       BCCHAR(IB),boundary_location(ib),
     2       ' for species ',IS, ' with parameters ',
     3       CONPAR(IS,IB,1),CONPAR(IS,IB,2)
            if (core_boundary(ib).and.ishigh)
     .       call xertst (ionising_core.le.0,
     .       'ionising_core switch usage inconsistent with BCCON = 13')
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
               bc_type13_ix=bc_list_x(iy,ib)
               bc_type13_iy=bc_list_y(iy,ib)
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
               VS=VS+VOL(bc_list_x(iy,ib),bc_list_y(iy,ib))
               NAS=NAS+NA(bc_list_x(iy,ib),bc_list_y(iy,ib),is)*
     .          vol(bc_list_x(iy,ib),bc_list_y(iy,ib))
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
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         VBIG*NAS*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .        -VBIG*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            ENDDO
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=13, ONLY CODED FOR SOUTH AND WEST')
          ENDIF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                       !srv 23.09.08 }

* -- BCCON=14 -- Sound speed flux

        case (14)
          if (ncall.eq.0) then
            call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1       'BCCON = 14, CONPAR(,,1) <= 0 not allowed!')
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'sound speed flux ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              S1 = ABS(
     1           PBS(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0))
              CS=SQRT(PZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     .                RZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
!             CS=SQRT(BCINT(PZ,NX,NY,IY,IB)/BCINT(RZ,NX,NY,IY,IB))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=0.0_R8
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .           SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)-
     .           S1*CONPAR(IS,IB,1)*CS
              sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,IS)=
     &           sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,IS) +
     &           SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,IS)
              if(mdf_fnb.ne.0) then                                   !srv 13.10.06 {
                t0 = -
     &          fnaPSch(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                  rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)
                if (t0.ge.0.0_R8) then
                  sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     &             sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)+t0
                else
                  sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     &             sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is)+t0/
     &               na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
                endif
              endif                                                   !srv 13.10.06 }
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              S1 = ABS(PBS(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))
              CS=SQRT(PZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     .                RZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
!              CS=SQRT(BCINT(PZ,NX,NY,IY,IB)/BCINT(RZ,NX,NY,IY,IB))
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=0.0_R8
              SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .           SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)-
     .           S1*CONPAR(IS,IB,1)*CS
              sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,IS)=
     &           sna0_no_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,IS) +
     .           SNA0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),:,IS)
              if(mdf_fnb.ne.0) then                                   !srv 13.10.06 {
                t0 = fnaPSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)
                if(t0.ge.0.0_R8) then
                  sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is) =
     &             sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)+t0
                else
                  sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is) =
     &             sna0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,is)+t0/
     &               na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
                endif
              endif                                                   !srv 13.10.06 }
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCCON=14, NOT IN RADIAL DIR')
          ENDIF

* -- BCCON=15 -- prescribe a radial leakage velocity, as a multiplier to the local thermal velocity

        case (15)                                                   !srv 25.09.09 16.02.12
          if(ncall.eq.0) then
            write(*,'(2a,1p,g14.7,a,a,a,a,i3)')
     .       'prescribe the decay length for the density',
     &       ' for drifs switched on ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            if (core_boundary(ib).and.ishigh)
     .       call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 15')
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

* -- BCCON=16 -- Feedback boundary condition: given total particle flux with constant density  !srv 16.02.12 {

        case (16)
          if(ncall.eq.0) then
#ifdef ASTRA
            ! IYS 08.04.2014 temporary solution while we do not have impurities in ASTRA
            call xertst(is.eq.ismain,
     1         'BCCON = 16 is applied for main ion only for ASTRA') !srv 29.09.15
#endif
            write(*,'(a,a,i3)')
     .       'given total particle flux with constant',
     &       ' density for ion is =', is
!           write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     1           boundary_location(ibw),' for species ',IS
                call xertst(core_boundary(ibw),
     .           'BCCON = 16 only applicable to core boundaries !')
              endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.19.or.
     1            BCCON(IS,IBW).EQ.21.or.BCCON(IS,IBW).EQ.22.or.
     1            BCCON(IS,IBW).EQ.23.or.BCCON(IS,IBW).EQ.25.or.
     1            BCCON(IS,IBW).EQ.26.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab
     .('BCCON types 16, 19, 21, 22, 23, 25, 26 and 27 cannot be mixed!')
              endif
            ENDDO
            if (core_boundary(ib).and.ishigh)
     .       call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 16')
          endif
          if(.not.once_bccon(is)) then
            IF(BCCHAR(IB).EQ.'S') THEN
              if (use_astra.ne.0) then
#ifdef ASTRA
                SCONPAR(IS) = fneb_astra
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SCONPAR(IS) = 0.0_R8
                DO IBW=IB,NBC
                  if(BCCHAR(IBW).EQ.'S'.and.
     >               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    SCONPAR(IS) = SCONPAR(IS) + CONPAR(IS,IBW,1)
                  endif
                ENDDO
              endif
              write(*,*) 'BCCON=16: Summed value is ', SCONPAR(IS)
              NAS = 0.0_R8
              VS = 0.0_R8
              DO IBW=IB,NBC
                if(BCCHAR(IBW).EQ.'S'.and.
     &             BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    NAS = NAS+NA(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)
     .                      *VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    VS = VS +VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  END DO
                endif
              END DO
              NAS=NAS/VS
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
              IF (max(abs(totflux),abs(SCONPAR(IS))).GT.0.0_R8)
     >         NAS = NAS *
     &         (1.0_R8-bc_ref*(totflux-SCONPAR(IS))/
     &                 max(abs(totflux),abs(SCONPAR(IS))))
!              T0 = (totflux-CONPAR(IS,IB,1))/
!     &           max(abs(totflux),CONPAR(IS,IB,1))
!              NAS = NAS *
!     &         (1.0_R8-max(-bc_ref,min(bc_ref,T0)))
              write(*,'(a,es14.7,a,es14.7,a,es14.7,a,i3)')
     &         'NAS=',NAS,' totflux=',totflux,
     &         ' SCONPAR=',SCONPAR(IS),' IS=',IS
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
              once_bccon(is) = .true.
            ELSE                                                      !srv 16.02.12 }
              CALL XERRAB
     1         ('B2STBC -- BCCON=16, NOT AVAILABLE FOR N,W,E')
            ENDIF
          endif                                                       !srv 16.02.12 }

* -- BCCON=17 -- FILLER

        case (17)
          CALL XERRAB ('B2STBC -- BCCON=17, NOT YET AVAILABLE')

#ifdef ASTRA
* -- BCCON=18 -- PRESCRIBE THE MODIFIED TOTAL PARTICLE FLUX

        case (18)
          if(ncall.eq.0) then
            write(*,'(a,a,a,a,i3)') 'modified total particle flux on ',
     1         BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(is.eq.ismain,
     1       'BCCON = 18 is applied for main ion only ')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            t0 = 0.0_R8
            US=0.0_R8
            isb = ismain0
            DO IX = 1, BC_LIST_SIZE(IB)
              US=US+SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              if (is_neutral(isb)) t0 = t0 +
     .             fna(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                 topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,isb)
            END DO
            CONPAR(IS,IB,1) = (fSeNBI_astra - t0)                     !srv 29.11.11

!            if(ncall.eq.0) then
               write (*,'(a,a,i3,1p,3(a,g11.4))')
     &        ' BCCON=18', ' IS =', IS,
     &        ' CONPAR=',CONPAR(IS,IB,1), ' fSeNBI_astra=',fSeNBI_astra,
     &        ' integrated neutral flux =',t0
!            endif

            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     &         SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)+          !srv 01.02.09
     1          CONPAR(IS,IB,1)/US*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &         sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+
     &         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)           !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &         sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)           !srv 23.09.08
              if(mdf_fnb.ne.0) then                                   !srv 12.09.08 {
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) -
     &           fnaPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                   topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
                sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &           sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) +
     &            vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &            hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &           vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
              endif                                                   !srv 12.09.08 }
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCCON=18, ONLY FOR SOUTH')
          ENDIF
#else
* -- BCCON=18 -- FILLER

        case (18)
          CALL XERRAB ('B2STBC -- BCCON=18, NOT YET AVAILABLE')
#endif

* -- BCCON=19 -- Feedback boundary condition: given sum of integrated neutrals and
*                main ion particle fluxes with constant ion density.

        case (19)
          if(ncall.eq.0) then
            if (is.lt.ns-1) call xertst(.not.lnext(is,is+1),
     1       'BCCON = 19 can only be applied to ions in their '//
     2       'highest ionization stage')
            write(*,*) 'given sum of particle flux for '//
     .                 'isonuclear sequence with constant highest '//
     .                 'ionisation stage density '
!           write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     1           boundary_location(ibw),' for species ',IS
                call xertst(core_boundary(ibw),
     .            'BCCON = 19 only applicable to core boundaries !')
              endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.21.or.BCCON(IS,IBW).EQ.22.or.
     1            BCCON(IS,IBW).EQ.23.or.BCCON(IS,IBW).EQ.25.or.
     1            BCCON(IS,IBW).EQ.26.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab
     .('BCCON types 16, 19, 21, 22, 23, 25, 26 and 27 cannot be mixed!')
              endif
            ENDDO
            if (ishigh) call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 19')
          endif
          if(.not.once_bccon(is)) then
            IF(BCCHAR(IB).EQ.'S') THEN
              if (use_astra.ne.0) then
#ifdef ASTRA
                SCONPAR(IS) = fSeNBI_astra
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SCONPAR(IS) = 0.0_R8
                DO IBW=IB,NBC
                  if(BCCHAR(IBW).EQ.'S'.and.
     >               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    SCONPAR(IS) = SCONPAR(IS) + CONPAR(IS,IBW,1)
                  endif
                ENDDO
              endif
              write(*,*) 'BCCON=19: Summed value is ', SCONPAR(IS)
              NAS = 0.0_R8
              VS = 0.0_R8
              DO IBW=IB,NBC
                if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    NAS = NAS+NA(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)
     .                      *VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    VS = VS +VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  END DO
                endif
              ENDDO
              NAS = NAS/VS
              totfluxn = 0.0_R8
              totfluxi = 0.0_R8
#ifdef B25_EIRENE
cxpb  Add contribution from Eirene core boundary fluxes
              do isa = 1, nnatmi
                if (latmscl(isa).eq.b2espcr(is)) then
                  do k = 1, ncns
                    ix = nlim + lcns(k)
                    totfluxn = totfluxn -
     &               wldna(ix,isa,0) + wldra(ix,isa,0) + wldpa(ix,isa,0)
                    do ism = 1, nnmoli
                      if (mlcmp(isa,ism).gt.0) then
                        totfluxn = totfluxn - mlcmp(isa,ism) *
     &                   (wldnm(ix,ism,0) - wldrm(ix,ism,0)
     &                                    - wldpm(ix,ism,0))
                      end if
                    end do
                  end do
                end if
              end do
#endif
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
! totfluxn here means the sum of fluxes over all ions (including fluid
! neutrals) of this isotope except for the highest ionization stage.
! For coupled runs, both fna and fna_mdf are rescaled for neutrals
! Note that for neutrals fna = fna_mdf
                    do isb = 0, ns-1
                      if(b2espcr(is).eq.b2espcr(isb).and.is.ne.isb) then
                        if(.not.is_neutral(isb).or.use_eirene.eq.0) then
                          if(mdf_fnb.ne.0) then
                            totfluxn = totfluxn +
     &       fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,isb)
                          else
                            totfluxn = totfluxn +
     &           fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,isb)
                          endif
                        endif
                      endif
                    enddo
! totfluxi here means the flux of the ion in the highest ionization stage.
                    if(mdf_fnb.ne.0) then
                      totfluxi = totfluxi +
     &       fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,is)
                    else
                      totfluxi = totfluxi +
     &           fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,is)
                    endif
                  ENDDO
                ENDIF
              ENDDO
              IF(max(abs(totfluxn+totfluxi),abs(SCONPAR(IS))).GT.0.0_R8)
     >         NAS = NAS *
     &         (1.0_R8-bc_ref*(totfluxn+totfluxi-SCONPAR(IS))/
     &           max(abs(totfluxn+totfluxi),abs(SCONPAR(IS))))
              write(*,'(a,es14.7,a,es14.7,a,es14.7,a,es14.7,a,es14.7)')
     *         'NAS=',NAS,' totfluxn=',totfluxn,
     &         ' totfluxi=',totfluxi,' SCONPAR=',SCONPAR(IS),
     &         ' totfluxn+totfluxi=',totfluxn+totfluxi
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*NAS*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  ENDDO
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)=
     &               SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  ENDDO                                              !srv 23.09.08 }
                ENDIF
              ENDDO
              once_bccon(is) = .true.
            ELSE                                                     !srv 16.02.12 }
              CALL XERRAB
     1         ('B2STBC -- BCCON=19, NOT AVAILABLE FOR N,W,E')
            ENDIF
          endif                                                      !srv 16.02.12 }

* -- BCCON=20 -- constant density feedback condition

        case (20)
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'constant density (neo) ',
     1       CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1       'BCCON = 1, CONPAR(,,1) !> 0')
            if (core_boundary(ib).and.ishigh)
     >       call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 20')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            naave = 0.0_R8
            volave = 0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type20_ix=bc_list_x(ix,ib)
              bc_type20_iy=bc_list_y(ix,ib)
              do bc_type20_i=1,bc_type20_ref  ! how many steps away from the boundary?
                bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
                bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
                bc_type20_ix=bc_type20_ixx
                call xertst (-2.lt.bc_type20_ix.and.
     &                             bc_type20_ix.lt.nx+1,
     1           'stepped off x-grid in type 20 bc')
                call xertst (-2.lt.bc_type20_iy.and.
     &                             bc_type20_iy.lt.ny+1,
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
              do bc_type20_i=1,bc_type20_ref  ! how many steps away from the boundary?
                bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
                bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
                bc_type20_ix=bc_type20_ixx
                call xertst (-2.lt.bc_type20_ix.and.
     &                             bc_type20_ix.lt.nx+1,
     1           'stepped off x-grid in type 20 bc')
                call xertst (-2.lt.bc_type20_iy.and.
     &                             bc_type20_iy.lt.ny+1,
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

cxpb Careful here because of BC number clash
* From 5.0: starting at version 01.001.030
* -- BCCON=21 -- constant density scaled by density on the ring bc_type21_ref away
*                has been moved up to next available number (BCCON=24)

*  From 5.2, after version number 01.001.078:
*            -- BCCON=21 -- PRESCRIBE THE VALUE OF THE DENSITY !lk and IYS 01.03.2013
*                           AND ADD A DENSITY PERTURBATION TO GET A SOLUTION
*                           WHICH IS AS CLOSE AS POSSIBLE TO NEOCLASSICAL THEORY
*                           DERIVED FOR PURE H/D/T PLASMA ONLY.
*
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH
* CORRESPONDING CONDITION ON ION TEMPERATURE (BCENI=23,24)

        case (21)
          if(ncall.eq.0) then
            if (boundary_version.ge.'01.001.030' .and.
     &          boundary_version.lt.'02.000.000') then
              write(*,*) 'BCCON=21 type B.C. moved to 24 ',
     &         ' because of clash with 5.2 numbering. '
              write(*,*) 'Please modify b2.boundary.parameters file'
              call xerrab('B2STBC -- BCCON=21, DISPLACED')
            endif
            write(*,'(a,1p,g14.7)')
     .       'constant density',CONPAR(IS,IB,1)
            call xertst(is.eq.ismain,
     1       'BCCON = 21 is applied for main ion only ')
            if (.not.(ns.eq.1.or.(ns.eq.2.and.is_neutral(0)))) then
              write(0,*)
     &         'Use of BCCON=21 only meant for single-fluid runs!'
              write(0,*) 'Use BCCON=25 instead!'
              write(6,*)
     &         'Use of BCCON=21 only meant for single-fluid runs!'
              write(6,*) 'Use BCCON=25 instead!'
            endif
            DO IBW=IB,NBC
              if(core_boundary(ibw).and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1           'BCCON = 21, CONPAR(,,1) !> 0')
                call xertst(CONPAR(IS,IB,1).eq.CONPAR(IS,IBW,1),  ! IYS 22.01.2015
     1           'BCCON = 21, Average density on all parts ' //
     1             'of the surface should be the same')
                if (bceni(ibw).ne.23.and.bceni(ibw).ne.24) then
                  write(0,*)
     &    'Should use BCCON(IS,IB).EQ.21 together with BCENI=23 or 24'
                  write(*,*)
     &    'Should use BCCON(IS,IB).EQ.21 together with BCENI=23 or 24'
                  write(*,'(a)')
     &             'Please ensure that this boundary condition '//
     &             'is consistent with '//
     &             'ion energy boundary condition'  !srv 21.09.12
                endif
              endif
              if(core_boundary(ibw).and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.22.or.
     1            BCCON(IS,IBW).EQ.23.or.BCCON(IS,IBW).EQ.25.or.
     1            BCCON(IS,IBW).EQ.26.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab
     .('BCCON types 16, 19, 21, 22, 23, 25, 26 and 27 cannot be mixed!')
              endif
            END DO
            if (ishigh) call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 21')
          endif
          IF(core_boundary(ib)) THEN
            if (.not.once_bccon(is)) then
              if(is.eq.ismain.and.use_astra.ne.0) then
#ifdef ASTRA
                DO IBW=IB,NBC
                  if(core_boundary(ibw).and.
     &               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    CONPAR(IS,IBW,1) = neb_astra                       !srv 28.07.11
                  endif
                ENDDO
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              endif
              t0=0.0e0_R8   ! Volume between two surfaces               !lk 12.12.12 {
              t1=0.0e0_R8   ! Sum of volume-averaged electron and ion temperatures
              DO IBW=IB,NBC
                IF (core_boundary(ibw).and.
     &              BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    t0=t0 +
     .                vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    t1=t1 +
     .                (ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &                 te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &                vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  ENDDO
                ENDIF
              ENDDO
              if (t0.gt.0.0_R8) t1=t1/t0        ! Sum of volume-averaged electron and ion temperatures
              DO IBW=IB,NBC
                IF (core_boundary(ibw).and.
     &              BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                  do IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &                 VBIG*CONPAR(IS,IBW,1)*t1
     &                 /(ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &                   te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &                   SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &                -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  enddo                                                !lk 12.12.12 }
                  DO IX = 1, BC_LIST_SIZE(IBW)                         !srv 23.09.08
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)=
     &               SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  ENDDO
                ENDIF
              ENDDO
              once_bccon(is) = .true.
            endif
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=21, VALID FOR CORE ONLY')
          ENDIF

* -- BCCON=22 -- Feedback boundary condition: given total particle flux with
*    constant average density. The poloidal density variation is added to make
*    the solution as close as possible to the neoclassical value derived for
*    pure H/D/T plasma.
*    It is recommended to use this boundary condition
*    together with corresponding boundary condition on ion temperature
*    (BCENI=23,24) !srv 16.02.12, IYS 01.03.2013 {

        case (22)
          if(ncall.eq.0) then
            call xertst(is.eq.ismain,
     1       'BCCON = 22 is applied for main ion only ')
            write(*,'(a,a,i3)')
     &       'given total particle flux with constant',
     &       ' density for ion is =', is
!           write all boundaries with this type of boundary condition
            if (.not.(ns.eq.1.or.(ns.eq.2.and.is_neutral(0)))) then
              write(0,*)
     &         'Use of BCCON=22 only meant for single-fluid runs!'
              write(0,*) 'Use BCCON=26 instead!'
              write(6,*)
     &         'Use of BCCON=22 only meant for single-fluid runs!'
              write(6,*) 'Use BCCON=26 instead!'
            endif
            DO IBW=IB,NBC
              if(core_boundary(ibw).and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                if (bceni(ibw).ne.23.and.bceni(ibw).ne.24) then
                  write(0,*)
     &    'Should use BCCON(IS,IB).EQ.22 together with BCENI=23 or 24'
                  write(*,*)
     &    'Should use BCCON(IS,IB).EQ.22 together with BCENI=23 or 24'
                  write(*,'(a)')
     &             'Please ensure that this boundary condition '//
     &             'is consistent with '//
     &             'ion energy boundary condition'  !srv 21.09.12
                endif
              endif
              if(core_boundary(ibw).and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.21.or.
     1            BCCON(IS,IBW).EQ.23.or.BCCON(IS,IBW).EQ.25.or.
     1            BCCON(IS,IBW).EQ.26.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab
     .('BCCON types 16, 19, 21, 22, 23, 25, 26 and 27 cannot be mixed!')
              endif
            ENDDO
            if (ishigh) call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 22')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            if(.not.once_bccon(is)) then
              if (use_astra.ne.0) then
#ifdef ASTRA
                SCONPAR(IS) = fneb_astra
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SCONPAR(IS) = 0.0_R8
                DO IBW=IB,NBC
                  if(core_boundary(ibw).and.
     >               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    SCONPAR(IS) = SCONPAR(IS) + CONPAR(IS,IBW,1)
                  endif
                ENDDO
              endif
              write(*,*) 'BCCON=22: Summed value is ', SCONPAR(IS)
              totflux = 0.0_R8
              t0=0.0e0_R8   ! Volume between two surfaces             !lk 12.12.12 {
              t1=0.0e0_R8   ! Sum of volume-averaged electron and ion temperatures
              NAS=0.0e0_R8
              DO IBW=IB,NBC
                IF(core_boundary(ibw).AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    t0=t0 +
     &                vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    t1=t1 +
     &                (ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &                 te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     .                vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    NAS=NAS+na(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    if(mdf_fnb.ne.0) then
                      totflux = totflux +
     &         fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,is)
                    else
                      totflux = totflux +
     &         fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &             topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,is)
                    endif
                  END DO
                ENDIF
              ENDDO
              t1=t1/t0        ! Sum of volume-averaged electron and ion temperatures
              NAS=NAS/t0      ! Volume-averaged density
              IF (max(abs(totflux),abs(SCONPAR(IS))).GT.0.0_R8)
     >         NAS = NAS*
     &         (1.0_R8-bc_ref*(totflux-SCONPAR(IS))/
     &           max(abs(totflux),abs(SCONPAR(IS))))
!              T0 = (totflux-CONPAR(IS,IB,1))/
!     &           max(abs(totflux),CONPAR(IS,IB,1))
!              NAS = NA((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB),is) *
!     &         (1.0_R8-max(-bc_ref,min(bc_ref,T0)))
              write(*,'(a,es14.7,a,es14.7,a,es14.7,a,i3)')
     .         'NAS=',NAS,' totflux=',totflux,
     &         ' SCONPAR=',SCONPAR(IS),' IS=',IS
              DO IBW=IB,NBC
                IF(core_boundary(ibw).AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*NAS*t1
     &                   /(ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &                     te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
     &                    *SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))   ! lk
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  ENDDO
                  DO IX = 1, BC_LIST_SIZE(IBW)                       !srv 23.09.08 {
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)
     &                =SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  ENDDO                                              !srv 23.09.08 }
                ENDIF
              ENDDO
              once_bccon(is) = .true.
            endif                                                    !srv 16.02.12 }
          ELSE                                                       !srv 16.02.12 }
            CALL XERRAB
     1       ('B2STBC -- BCCON=22 IS AVAILABLE FOR CORE ONLY')
          ENDIF

* -- BCCON=23 -- Feedback boundary condition: given sum of integrated neutrals and
*    main ion particle fluxes with constant average density.
*    The poloidal density variation is added to make the solution
*    as close as possible to the neoclassical value derived for pure H/D/T plasma.
*    It is recommended to use this boundary condition with corresponding
*    boundary condition on ion temperature (BCENI=23,24) !srv 16.02.12, IYS 04.03.2013
*    Use this boundary condition for pure H/D/T plasma only,
*    If a minor impurity is present, please use the BCCON=27 instead. ! IYS 21.01.2016

       case (23)
          if(ncall.eq.0) then
            call xertst(is.eq.ismain,
     1       'BCCON = 23 is applied for main ion only ')
            write(*,'(a,a,i3)')
     .       'given total particle flux with constant',
     &       ' density for main ion is =', ismain
!          write all boundaries with this type of boundary condition
            if (.not.(ns.eq.1.or.(ns.eq.2.and.is_neutral(0)))) then
              write(0,*)
     &         'Use of BCCON=23 only meant for single-fluid runs!'
              write(0,*) 'Use BCCON=27 instead!'
              write(6,*)
     &         'Use of BCCON=23 only meant for single-fluid runs!'
              write(0,*) 'Use BCCON=27 instead!'
            endif
            DO IBW=IB,NBC
              if(core_boundary(ibw).and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                call xertst(CONPAR(IS,IBW,1).ge.0.0_R8,
     1           'BCCON = 23, CONPAR(,,1) < 0')
                if (bceni(ibw).ne.23.and.bceni(ibw).ne.24) then
                  write(0,*)
     &     'Should use BCCON(IS,IB).EQ.23 together with BCENI=23 or 24'
                  write(*,*)
     &     'Should use BCCON(IS,IB).EQ.23 together with BCENI=23 or 24'
                  write(*,'(a)')
     &             'Please ensure that this boundary condition '//
     &             'is consistent with '//
     &             'ion energy boundary condition'                   !srv 21.09.12
                endif
              endif
              if(core_boundary(ibw).and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.21.or.
     1            BCCON(IS,IBW).EQ.22.or.BCCON(IS,IBW).EQ.25.or.
     1            BCCON(IS,IBW).EQ.26.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab
     .('BCCON types 16, 19, 21, 22, 23, 25, 26 and 27 cannot be mixed!')
              endif
            ENDDO
            if (ishigh) call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 23')
          endif

          IF(core_boundary(ib)) THEN
            if(.not.once_bccon(is)) then
              if (use_astra.ne.0) then
#ifdef ASTRA
                SCONPAR(IS) = fSeNBI_astra
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SCONPAR(IS) = 0.0_R8
                DO IBW=IB,NBC
                  if(core_boundary(ibw).and.
     &               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    SCONPAR(IS) = SCONPAR(IS) + CONPAR(IS,IBW,1)
                  endif
                ENDDO
              endif
              write(*,*) 'BCCON=23: Summed value is ', SCONPAR(IS)
              totfluxn = 0.0_R8
              totfluxi = 0.0_R8
              t0=0.0e0_R8   ! Volume between two surfaces               !lk 12.12.12 {
              t1=0.0e0_R8   ! Sum of volume-averaged electron and ion temperatures
              NAS=0.0e0_R8
#ifdef B25_EIRENE
cxpb  Add contribution from Eirene core boundary fluxes
              do isa = 1, nnatmi
                if (latmscl(isa).eq.b2espcr(ismain0)) then
                  do k = 1, ncns
                    ix = nlim + lcns(k)
                    totfluxn = totfluxn -
     &               wldna(ix,isa,0) + wldra(ix,isa,0) + wldpa(ix,isa,0)
                    do ism = 1, nnmoli
                      if (mlcmp(isa,ism).gt.0) then
                        totfluxn = totfluxn - mlcmp(isa,ism) *
     &                   (wldnm(ix,ism,0) - wldrm(ix,ism,0)
     &                                    - wldpm(ix,ism,0))
                      end if
                    end do
                  end do
                end if
              end do
#endif
              DO IBW=IB,NBC
                IF(core_boundary(ibw).AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    isb = ismain0
                    if (is_neutral(isb)) totfluxn = totfluxn +
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
                    t0=t0+vol(topix(BC_LIST_X(IX,IBW),
     .                              BC_LIST_Y(IX,IBW)),
     &                        topiy(BC_LIST_X(IX,IBW),
     .                              BC_LIST_Y(IX,IBW)))
                    t1=t1 + (ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     .                       te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     .                vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    NAS=NAS+na(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  ENDDO
                ENDIF
              ENDDO
              t1=t1/t0        ! Sum of volume-averaged electron and ion temperatures
              NAS=NAS/t0      ! Volume-averaged density
              IF(abs(totfluxn)+abs(totfluxi)+abs(SCONPAR(IS)).GT.0.0_R8)
     >         NAS = NAS*
     &         (1.0_R8-bc_ref*(totfluxn+totfluxi-SCONPAR(IS))/
     &           (abs(totfluxn)+abs(totfluxi)+abs(SCONPAR(IS))))
              write(*,'(5(a,es14.7))')
     *          'NAS=',NAS,' totfluxn=',totfluxn,
     &          ' totfluxi=',totfluxi,' SCONPAR=', SCONPAR(IS),
     &          ' totfluxn+totfluxi=',totfluxn+totfluxi
              DO IBW=IB,NBC
                IF(core_boundary(ibw).AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*NAS*t1
     &               /(ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &                 te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &                 SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))         ! lk
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  ENDDO
                  DO IX =1, BC_LIST_SIZE(IBW)                          !srv 23.09.08 {
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)=
     &               SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  ENDDO                                                !srv 23.09.08 }
                ENDIF
              ENDDO
              once_bccon(is) = .true.
            endif                                                      !srv 16.02.12 }
          ELSE                                                         !srv 16.02.12 }
            CALL XERRAB
     1        ('B2STBC -- BCCON=23 IS AVAILABLE FOR CORE ONLY')
          ENDIF

* -- BCCON=24 -- constant density scaled by density on the ring bc_type21_ref away

        case (24)
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'constant density (modified by nearby ring) ',
     1       CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1       'BCCON = 24, CONPAR(,,1) !> 0')
            call xertst(CONPAR(is,ib,2).gt.0.0_R8,
     &           'faulty parameter conpar2')
            if (ishigh.and.core_boundary(ib))
     >       call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 24')
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
     &                             bc_type21_ix.lt.nx+1,
     1           'stepped off x-grid in type 24 bc')
                call xertst (-2.lt.bc_type21_iy.and.
     &                             bc_type21_iy.lt.ny+1,
     1           'stepped off y-grid in type 24 bc')
              enddo
              naave=naave+na(bc_type21_ix,bc_type21_iy,is)*
     1                   vol(bc_type21_ix,bc_type21_iy)
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
     &                             bc_type21_ix.lt.nx+1,
     1           'stepped off x-grid in type 24 bc')
                call xertst (-2.lt.bc_type21_iy.and.
     &                             bc_type21_iy.lt.ny+1,
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
     &                             bc_type21_ix.lt.nx+1,
     1           'stepped off x-grid in type 24 bc')
                call xertst (-2.lt.bc_type21_iy.and.
     &                             bc_type21_iy.lt.ny+1,
     1           'stepped off y-grid in type 24 bc')
              enddo
              naave=naave+na(bc_type21_ix,bc_type21_iy,is)*
     1                   vol(bc_type21_ix,bc_type21_iy)
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
     &                             bc_type21_ix.lt.nx+1,
     1           'stepped off x-grid in type 24 bc')
                call xertst (-2.lt.bc_type21_iy.and.
     &                             bc_type21_iy.lt.ny+1,
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

* -- BCCON=25  Feedback boundary condition: PRESCRIBE THE AVERAGE VALUE OF THE DENSITY  !lk and IYS 01.03.2013
*              AND ADD A DENSITY PERTURBATION FROM NEIGHBOURING RADIAL CELL.
*              THIS BOUNDARY CONDITION IS SUITABLE FOR ANY SPECIES OF
*              MULTI-SPECIES PLASMA (I.E. IN THE CASE WHEN THE CONDITION USED
*              IN BCCON=21 FAILS)
*              IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH
*              CORRESPONDING CONDITION ON ION TEMPERATURE (BCENI=26,27)

        case (25)
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7)')
     .       'constant average density',CONPAR(IS,IB,1)
            DO IBW=IB,NBC
              if(core_boundary(ibw).and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                call xertst(CONPAR(IS,IBW,1).gt.0.0_R8,
     1           'BCCON = 25, CONPAR(,,1) !> 0')
                call xertst(CONPAR(IS,IB,1).eq.CONPAR(IS,IBW,1),  ! IYS 22.01.2015
     1           'BCCON = 25, Average density on all parts ' //
     1             'of the surface should be the same')
                if (bceni(ibw).ne.26.and.bceni(ibw).ne.27) then
                  write(0,*)
     &    'Should use BCCON(IS,IB).EQ.25 together with BCENI=26 or 27'
                  write(*,*)
     &    'Should use BCCON(IS,IB).EQ.25 together with BCENI=26 or 27'
                  write(*,'(a)')
     &             'Please ensure that this boundary condition '//
     &             'is consistent with '//
     &             'ion energy boundary condition'  !srv 21.09.12
                endif
              endif
              if(core_boundary(ibw).and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.21.or.
     1            BCCON(IS,IBW).EQ.22.or.BCCON(IS,IBW).EQ.23.or.
     1            BCCON(IS,IBW).EQ.26.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab
     .('BCCON types 16, 19, 21, 22, 23, 25, 26 and 27 cannot be mixed!')
              endif
            END DO
            if (ishigh) call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 25')
          endif

          IF(core_boundary(ib)) THEN
            if (.not.once_bccon(is)) then
              if(is.eq.ismain.and.use_astra.ne.0) then
#ifdef ASTRA
                DO IBW=IB,NBC
                  if(core_boundary(ibw).and.
     &               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    CONPAR(IS,IBW,1) = neb_astra                       !srv 28.07.11
                  endif
                ENDDO
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              endif
              t0=0.0e0_R8   ! Volume between two surfaces              !lk 12.12.12 {
              t1=0.0e0_R8   ! Sum of volume-averaged density
              DO IBW=IB,NBC
                if(core_boundary(ibw).and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  do IX = 1, BC_LIST_SIZE(IBW)
                    t0=t0 +
     .                vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    t1=t1 +
     .                 na(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),is)
     .               *vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  enddo
                endif
              enddo
              t1=t1/t0
              DO IBW=IB,NBC
                if(core_boundary(ibw)) then
                  do IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*(CONPAR(IS,IBW,1) +
     &                 na(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),IS)
     .                -t1)*
     &               SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  enddo
                endif                                                  !lk 12.12.12 }
                DO IX = 1, BC_LIST_SIZE(IBW)                           !srv 23.09.08
                  SNA0_NO_MDF(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)=
     &             SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                ENDDO
              ENDDO
              once_bccon(is) = .true.
            endif
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=25, VALID FOR CORE ONLY')
          ENDIF

* -- BCCON=26  Feedback boundary condition: PRESCRIBE THE TOTAL ION FLUX AND
*              FIND THE AVERAGE DENSITY.  !lk and IYS 01.03.2013
*              A DENSITY PERTURBATION IS TAKEN FROM NEIGHBOURING RADIAL CELL.
*              THIS BOUNDARY CONDITION IS SUITABLE FOR ANY SPECIES OF
*              MULTI-SPECIES PLASMA (I.E. IN THE CASE WHEN THE CONDITION USED
*              IN BCCON=22 FAILS)
*              IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH
*              A CORRESPONDING CONDITION ON THE ION TEMPERATURE (BCENI=26,27)

        case (26)
          if(ncall.eq.0) then
            write(*,'(a,a,i3)')                                     !IYS 21.01.2016
     .       'given total particle flux with constant',             !IYS 21.01.2016
     &       ' average density for ion is =', is                    !IYS 21.01.2016
!           write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
              if(core_boundary(ibw).and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                if (bceni(ibw).ne.26.and.bceni(ibw).ne.27) then
                  write(0,*)
     &    'Should use BCCON(IS,IB).EQ.26 together with BCENI=26 or 27'
                  write(*,*)
     &    'Should use BCCON(IS,IB).EQ.26 together with BCENI=26 or 27'
                  write(*,'(a)')
     &           'Please ensure that this boundary condition '//
     &           'is consistent with '//
     &           'ion energy boundary condition'  !srv 21.09.12
                endif
              endif
              if(core_boundary(ibw).and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.21.or.
     1            BCCON(IS,IBW).EQ.22.or.BCCON(IS,IBW).EQ.23.or.
     1            BCCON(IS,IBW).EQ.25.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab
     .('BCCON types 16, 19, 21, 22, 23, 25, 26 and 27 cannot be mixed!')
              endif
            ENDDO
            if (ishigh) call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 26')
          endif  ! ncall == 0
          IF(core_boundary(ib)) THEN
            if(.not.once_bccon(is)) then
              if (use_astra.ne.0) then
#ifdef ASTRA
                SCONPAR(IS) = fneb_astra
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SCONPAR(IS) = 0.0_R8
                DO IBW=IB,NBC
                  if(core_boundary(ibw).and.
     >               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    SCONPAR(IS) = SCONPAR(IS) + CONPAR(IS,IBW,1)
                  endif
                ENDDO
              endif
              write(*,*) 'BCCON=26: Summed value is ', SCONPAR(IS)
              totflux = 0.0_R8
              t0=0.0e0_R8   ! Volume between two surfaces               !lk 12.12.12 {
              t1=0.0e0_R8   ! Volume-averaged density on previous surface
              NAS=0.0e0_R8
              DO IBW=IB,NBC
                IF(core_boundary(ibw).AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    t0=t0 +
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    t1=t1 +
     .                na(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),IS)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    NAS=NAS+na(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    if(mdf_fnb.ne.0) then
                      totflux = totflux +
     &         fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,is)
                    else
                      totflux = totflux +
     &         fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &             topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,is)
                    endif
                  ENDDO
                ENDIF
              ENDDO
              t1=t1/t0        ! Sum of volume-averaged electron and ion temperatures
              NAS=NAS/t0      ! Volume-averaged density
              if (max(abs(totflux),abs(SCONPAR(IS))).GT.0.0_R8)
     >         NAS = NAS*
     &           (1.0_R8-bc_ref*(totflux-SCONPAR(IS))/
     &             max(abs(totflux),abs(SCONPAR(IS))))
              write(*,'(a,es14.7,a,es14.7,a,es14.7,a,i3)')
     .         'NAS=',NAS,' totflux=',totflux,
     .         ' SCONPAR=',SCONPAR(IS),' IS=',IS
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*(NAS +
     &                 na(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),IS)
     &                 -t1)
     &               *SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))          ! lk
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  ENDDO
                  DO IX = 1, BC_LIST_SIZE(IBW)                         !srv 23.09.08 {
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)
     &               =SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  ENDDO                                                !srv 23.09.08 }
                ENDIF
              ENDDO
              once_bccon(is) = .true.
            endif                                                      !srv 16.02.12 }
          ELSE                                                         !srv 16.02.12 }
            CALL XERRAB
     1        ('B2STBC -- BCCON=26 IS AVAILABLE FOR CORE ONLY')
          ENDIF

* -- BCCON=27  Feedback boundary condition: prescribe the particle flux sum
*              for all neutrals and ions belonging to a given isonuclear sequence
*              and find the average density of the highest ionization stage.
*              A DENSITY PERTURBATION IS TAKEN FROM THE NEIGHBOURING RADIAL CELL.
*              IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH
*              A CORRESPONDING CONDITION ON THE ION TEMPERATURE (BCENI=26,27)

        case (27)
          if(ncall.eq.0) then
            if (is.lt.ns-1) call xertst(.not.lnext(is,is+1),
     1       'BCCON = 27 is applied to highest ionisation stages only')
            write(*,'(a,a,i3)')
     &       'given total particle flux with constant',
     &       ' density for highest ionisation stage ion is =', is
            DO IBW=IB,NBC
              if(core_boundary(ibw).and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                call xertst(CONPAR(IS,IBW,1).ge.0.0_R8,
     1           'BCCON = 27, CONPAR(,,1) < 0')
                if (bceni(ibw).ne.26.and.bceni(ibw).ne.27) then
                  write(0,*)
     &    'Should use BCCON(IS,IB).EQ.27 together with BCENI=26 or 27'
                  write(*,*)
     &    'Should use BCCON(IS,IB).EQ.27 together with BCENI=26 or 27'
                  write(*,'(a)')
     &             'Please ensure that this boundary condition '//
     &             'is consistent with '//
     &             'ion energy boundary condition'                    !srv 21.09.12
                endif
              endif
              if(core_boundary(ibw).and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.21.or.
     1            BCCON(IS,IBW).EQ.22.or.BCCON(IS,IBW).EQ.23.or.
     1            BCCON(IS,IBW).EQ.25.or.BCCON(IS,IBW).EQ.26)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab
     .('BCCON types 16, 19, 21, 22, 23, 25, 26 and 27 cannot be mixed!')
              endif
            ENDDO
            call xertst (ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 27')
          endif
          IF(core_boundary(ib)) THEN
            if(.not.once_bccon(is)) then
              if (use_astra.ne.0) then
#ifdef ASTRA
                SCONPAR(IS) = fSeNBI_astra
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SCONPAR(IS) = 0.0_R8
                DO IBW=IB,NBC
                  if(core_boundary(ibw).and.
     >               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    SCONPAR(IS) = SCONPAR(IS) + CONPAR(IS,IBW,1)
                  endif
                ENDDO
              endif
              write(*,*) 'BCCON=27: Summed value is ', SCONPAR(IS)
              totfluxn = 0.0_R8
              totfluxi = 0.0_R8
              t0=0.0e0_R8   ! Volume between two surfaces              !lk 12.12.12 {
              t1=0.0e0_R8   ! Volume-averaged on neighbouring surface
              NAS=0.0e0_R8
#ifdef B25_EIRENE
cxpb  Add contribution from Eirene core boundary fluxes
              do isa = 1, nnatmi
                if (latmscl(isa).eq.b2espcr(is)) then
                  do k = 1, ncns
                    ix = nlim + lcns(k)
                    totfluxn = totfluxn -
     &                 (wldna(ix,isa,0) - wldra(ix,isa,0)
     &                                  - wldpa(ix,isa,0))
                    do ism = 1, nnmoli
                      if (mlcmp(isa,ism).gt.0) then
                        totfluxn = totfluxn - mlcmp(isa,ism) *
     &                   (wldnm(ix,ism,0) - wldrm(ix,ism,0)
     &                                    - wldpm(ix,ism,0))
                      end if
                    end do
                  end do
                end if
              end do
#endif
              DO IBW=IB,NBC
                IF(core_boundary(ibw).AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
! totfluxn here means the sum of fluxes over all ions (including fluid
! neutrals) of this isotope except for the highest ionization stage.
! For coupled runs, both fna and fna_mdf are rescaled for neutrals
! Note that for neutrals fna = fna_mdf
                    do isb = 0, ns-1
                      if(b2espcr(is).eq.b2espcr(isb).and.is.ne.isb) then
                        if(.not.is_neutral(isb).or.use_eirene.eq.0) then
                          if(mdf_fnb.ne.0) then
                            totfluxn = totfluxn +
     &       fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,isb)
                          else
                            totfluxn = totfluxn +
     &           fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,isb)
                          endif
                        endif
                      endif
                    enddo
! totfluxi here means the flux of the ion in the highest ionization stage.
                    if(mdf_fnb.ne.0) then
                      totfluxi = totfluxi +
     &       fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,is)
                    else
                      totfluxi = totfluxi +
     &           fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,is)
                    endif
                    t0=t0 +
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    t1=t1 +
     .                na(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   IS)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    NAS=NAS+na(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  ENDDO
                ENDIF
              ENDDO
              t1=t1/t0        !
              NAS = NAS/t0    ! Volume-averaged density
              IF(abs(totfluxn)+abs(totfluxi)+abs(SCONPAR(IS)).GT.0.0_R8)
     >         NAS = NAS*
     &         (1.0_R8-bc_ref*(totfluxn+totfluxi-SCONPAR(IS))/
     &           (abs(totfluxn)+abs(totfluxi)+abs(SCONPAR(IS))))
              write(*,'(5(a,es14.7))')
     *         'NAS=',NAS,' totfluxn=',totfluxn,
     &         ' totfluxi=',totfluxi,' SCONPAR=',SCONPAR(IS),
     &         ' totfluxn+totfluxi=',totfluxn+totfluxi
              DO IBW=IB,NBC
                IF(core_boundary(ibw).AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*(NAS +
     &                 na(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),IS)
     &                -t1)
     &               *SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))          ! lk
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  ENDDO
                  DO IX =1, BC_LIST_SIZE(IBW)                          !srv 23.09.08 {
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)=
     &               SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  ENDDO                                                !srv 23.09.08 }
                ENDIF
              ENDDO
              once_bccon(is) = .true.
            endif                                                      !srv 16.02.12 }
          ELSE                                                         !srv 16.02.12 }
            CALL XERRAB
     1       ('B2STBC -- BCCON=27 IS AVAILABLE FOR CORE ONLY')
          ENDIF

        case default
          write(*,*) BCCON(IS,IB)
          call xerrab ('BCCON out of bounds')
        end select
          
* -- BOUNDARY CONDITIONS FOR PARALLEL MOMENTUM EQUATIONS --

        select case (BCMOM(IS,IB))
        case (0)
          continue

* -- BCMOM=1 -- PRESCRIBE THE VALUE OF THE PARALLEL VELOCITY

        case (1)
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

* -- BCMOM=2 -- PRESCRIBE THE GRADIENT OF THE PARALLEL VELOCITY

        case (2)
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

* -- BCMOM=3 -- SHEATH CONDITIONS, MACH NUMBER AS INPUT

        case (3)
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'sheath velocity bc ',MOMPAR(IS,IB,2),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              IF(MOMPAR(IS,IB,2).LT.0.5_R8) THEN
                CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))
                UUOUT=MOMPAR(IS,IB,1)*CS*
     .           abs(wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .               wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3))
                if(.not.is_neutral(is)) then                           !xpb
                  UPOUT=max(0.0_R8, UUOUT/
     .             abs(wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .                 wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)) +    !xpb
     .             bcint(ud(-1,-1,is),nx,ny,iy,ib)*
     .                (wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/
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
                UUTEST=abs(
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
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))
                UUOUT=MAX(MOMPAR(IS,IB,1)*CS,UUTEST)
                if(.not.is_neutral(is)) then                           !xpb
                  UPOUT=max(0.0_R8, UUOUT+                             !xpb
     .             bcint(ud(-1,-1,is),nx,ny,iy,ib)*
     .                (wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/
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
     2         (-UPOUT)*
     .         SIGN(1.0_R8,BB(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))    !csc
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .         -(UBIG/EPSLON**2)*AM(IS)*mp*
     .         bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .         SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              IF(MOMPAR(IS,IB,2).LT.0.5_R8) THEN
                CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))
                UUOUT=MOMPAR(IS,IB,1)*CS*
     .           abs(wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .               wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3))
                if(.not.is_neutral(is)) then                           !xpb
                  UPOUT=max(0.0_R8, UUOUT/                             !xpb
     .              abs(wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .                  wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3))     !xpb
     .             -bcint(ud(-1,-1,is),nx,ny,iy,ib)*
     .                 (wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/
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
                UUTEST=abs(
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
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))
                UUOUT=MAX(MOMPAR(IS,IB,1)*CS,UUTEST)
                if(.not.is_neutral(is)) then                           !xpb
                  UPOUT=max(0.0_R8, UUOUT -                            !xpb
     .             bcint(ud(-1,-1,is),nx,ny,iy,ib)*
     .                (wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/
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
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))*UPOUT*
     .         SIGN(1.0_R8,BB(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))    !csc
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

* -- BCMOM=4 -- PRESCRIBE THE VALUE OF THE VELOCITY, WEAKLY
*                 A MIXED BOUNDARY CONDITION

        case (4)
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

* -- BCMOM=5 -- PRESCRIBE THE PARALLEL MOMENTUM FLUX PER UNIT AREA

        case (5)
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

* -- BCMOM=6 -- PRESCRIBE THE TOTAL PARALLEL MOMENTUM FLUX FOR A
*                 CONSTANT PARALLEL VELOCITY

        case (6)
          CALL XERRAB ('B2STBC -- BCMOM=6, NOT YET AVAILABLE')

* -- BCMOM=7 -- PRESCRIBE THE PARALLEL MOMENTUM AS A FUNCTION OF OTHER
*                 PLASMA PARAMETERS

        case (7)
          CALL XERRAB ('B2STBC -- BCMOM=7, NOT YET AVAILABLE')

* -- BCMOM=8 -- SPECIAL : LIMITED SHEAR

        case (8)
          if(ncall.eq.0) then
            write(*,'(a,a,a,a,i3)')
     .       'limited shear velocity on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
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

* -- BCMOM=9 -- PRESCRIBE THE TOTAL PARALLEL MOMENTUM FLUX WITH
*                 CONSTANT FLUX DENSITY

        case (9)
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

* -- BCMOM=10 -- PRESCRIBE THE DECAY LENGTH FOR THE PARALLEL MOMENTUM

        case (10)
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
              write(*,'(a)')
     &                   'Please ensure that this boundary condition '//
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

* -- BCMOM=11 -- ROZHANSKY VISCOSITY CONDITION FOR THE PARALLEL MOMENTUM

        case (11)
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'ROZ parallel momentum bc on ',BCCHAR(IB),
     1       boundary_location(ib),
     2       ' for species ',IS
            write(*,*) ' BCMOM=11 no longer supported !'
            write(0,*) ' BCMOM=11 no longer supported !'
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              IF (BC_LIST_X(IX,IB).GE.0) THEN
                smo0(bc_list_x(ix,ib),bc_list_y(ix,ib),0,is)=
     1           vsa0(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .           na(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*am(is)*
     2           sqrt(1.0_R8-pit(bc_list_x(ix,ib),bc_list_y(ix,ib))**2)*
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
     .           sqrt(1.0_R8-pit(bc_list_x(ix,ib),bc_list_y(ix,ib))**2)/
     .           hy1(bc_list_x(ix,ib),bc_list_y(ix,ib))
              end if
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCMOM=11, NOT AVAILABLE FOR N,W,E')
          ENDIF

* -- BCMOM=12 -- CONDITION from b2stbc_spb FOR THE PARALLEL MOMENTUM   !srv 01.02.09

        case (12)
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

* -- BCMOM=13 -- CONDITION from b2stbc_spb FOR THE PARALLEL MOMENTUM

        case (13)
          if(ncall.eq.0) then   !srv added accumulation in order to account twice b.c. for corner cells
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a,a,i3)')
     .       'from b2stbc_spb parallel momentum bc with mompar1 =',
     1       MOMPAR(IS,IB,1),' and mompar2 =', MOMPAR(IS,IB,2),' on ',
     2       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
!              CS=SQRT(BCINT(PZ,NX,NY,IY,IB)/BCINT(RZ,NX,NY,IY,IB))
              cs=sqrt(pz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &                rz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
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
!             CS=SQRT(BCINT(PZ,NX,NY,IY,IB)/BCINT(RZ,NX,NY,IY,IB))
              cs=sqrt(pz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &                rz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              s1hz = pbshz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)      !srv 11.01.13
             veast = cs*s1hz-                                      !srv 17.02.04 19.02.08 25.04.08 { 11.01.13
     &          hz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     &          vol(leftix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &              leftiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     &          hx(leftix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &             leftiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))        !*signmf*bb(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/bb(bc_list_x(iy,ib),bc_list_y(iy,ib),3) !srv 05.02.05
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

* -- BCMOM=14 -- CONDITION from b2stbc_spb FOR THE PARALLEL MOMENTUM  !srv 01.02.09 {

        case (14)
          if(ncall.eq.0) then !srv added accumulation in order to account twice b.c. for corner cells
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a,a,i3)')
     .       'from b2stbc_spb momentum bc with mompar1 =',
     1       MOMPAR(IS,IB,1),' and mompar2 =', MOMPAR(IS,IB,2),' on ',
     2       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(core_boundary(ib),
     &       'BCMOM = 14 intended for core boundaries only!')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              s0hz = hz(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &              vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &              hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))              !srv 09.01.01
              if (.not.is_neutral(is)) then                             !srv 15.02.02 {
                select case (iform_corr_dia)                            !srv 25.07.05
                case (0)
                  smo0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &             smo0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+        !srv 09.01.01 01.07.05
     &             s0hz*(am(is)*mp)*naaverage(is)*                      !srv 15.02.02
     &            (-2.0e0_R8*
     &             facdrift(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*phm2*    !srv 07.06.00 17.05.99 25.04.00 15.02.02
     &             tiaverage*bb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),2)/
     &             (qe*rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &             hx(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*              !srv 15.02.02
     &            (OnedBsqc(rightix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                      rightiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))-
     &              OnedBsqc(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))        !srv 05.03.99 15.02.02
     &             )*uaaverage(is) + corr_flux*s0hz                     !srv 15.02.02 08.09.03 25.07.05
                case (1)
                  smo0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &             smo0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+        !srv 09.01.01 01.07.05
     &             s0hz*(am(is)*mp)*naaverage(is)*                      !srv 15.02.02
     &           (-2.0e0_R8*facdrift(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &             phm2*                                                !srv 07.06.00 17.05.99 25.04.00 15.02.02
     &             tiaverage*bb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),2)/
     &            (qe*rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &             hx(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*              !srv 15.02.02
     &            (OnedBsqc(rightix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                      rightiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))-
     &             OnedBsqc(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))         !srv 05.03.99 15.02.02
     &            )*ua(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)+corr_flux
     &             *s0hz                                                !srv 08.09.03 25.07.05
                end select                                              !srv 25.07.05
              endif                                                     !srv 15.02.02
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

* -- BCMOM=15 -- PRESCRIBE THE VALUE OF THE PARALLEL VELOCITY, SCALED WITH
*                B_average/B_local

        case (15)
          if(ncall.eq.0) then
            write(*,'(a,a,1p,g14.7,a,a,a,a,i3)')
     .       'specified parallel velocity as a function of',
     &       ' poloidal coordinate ',MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(core_boundary(ib),
     &       'BCMOM = 15 intended for core boundaries only!')
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

* -- BCMOM=16 -- PRESCRIBE THE AVERAGE VALUE OF THE PARALLEL VELOCITY   !lk 12.12.12 {

        case (16)
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
            do ibw = 1, nbc
              if (core_boundary(ibw)) then
                do ix = 1, bc_list_size(ibw)
                  t0=t0+vol(topix(bc_list_x(ix,ibw),bc_list_y(ix,ibw)),
     .                      topiy(bc_list_x(ix,ibw),bc_list_y(ix,ibw)))
                  t1=t1+
     .              ua(topix(bc_list_x(ix,ibw),bc_list_y(ix,ibw)),
     .                 topiy(bc_list_x(ix,ibw),bc_list_y(ix,ibw)),is)*
     .             vol(topix(bc_list_x(ix,ibw),bc_list_y(ix,ibw)),
     .                 topiy(bc_list_x(ix,ibw),bc_list_y(ix,ibw)))
                end do
              endif
            enddo
            if (t0.gt.0.0_R8) t1=t1/t0
            if (core_boundary(ib)) then
              do ix = 1, BC_LIST_SIZE(IB)
                SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .           SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)+
     .             VBIG*NE00*AM(IS)*
     .         (UA(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .             TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS)-t1)*
     .              SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              end do                        !lk 12.12.12 }
            endif
          ELSE
            CALL XERRAB ('B2STBC -- BCMOM=16, IS VALID FOR CORE ONLY')
          ENDIF

* -- BCMOM=17 -- leakage option for parallel momentum

        case (17)
          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'parallel momentum leakage ',MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     1         SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     2         sqrt((rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .                te(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .                ti(bc_list_x(ix,ib),bc_list_y(ix,ib)))/
     3         (am(is)*mp))*MOMPAR(IS,IB,1)*(am(is)*mp)*
     4         NA(bc_list_x(ix,ib),bc_list_y(ix,ib),is)
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     1         SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     2         sqrt((rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .                te(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .                ti(bc_list_x(ix,ib),bc_list_y(ix,ib)))/
     3         (am(is)*mp))*MOMPAR(IS,IB,1)*(am(is)*mp)*
     4         NA(bc_list_x(ix,ib),bc_list_y(ix,ib),is)
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=0.0_R8
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     1         SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     2         sqrt((rza(bc_list_x(iy,ib),bc_list_y(iy,ib),is)*
     .                te(bc_list_x(iy,ib),bc_list_y(iy,ib))+
     .                ti(bc_list_x(iy,ib),bc_list_y(iy,ib)))/
     3         (am(is)*mp))*MOMPAR(IS,IB,1)*(am(is)*mp)*
     4         NA(bc_list_x(iy,ib),bc_list_y(iy,ib),is)
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=0.0_R8
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     1         SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     2         sqrt((rza(bc_list_x(iy,ib),bc_list_y(iy,ib),is)*
     .                te(bc_list_x(iy,ib),bc_list_y(iy,ib))+
     .                ti(bc_list_x(iy,ib),bc_list_y(iy,ib)))/
     3         (am(is)*mp))*MOMPAR(IS,IB,1)*(am(is)*mp)*
     4         NA(bc_list_x(iy,ib),bc_list_y(iy,ib),is)
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCMOM=17, NOT AVAILABLE')
          ENDIF

        case default
          write(*,*) BCMOM(IS,IB)
          call xerrab ('BCMOM out of bounds')
        end select
          
* -- END LOOP OVER SPECIES
      enddo

      return
      end subroutine b2stbc_phys_species

      subroutine b2stbc_phys_ene(IB)
      integer ib
      integer ix, iy, is, ibw, ism
      integer bc_type16_i, bc_type20_i, bc_type21_i

* -- BOUNDARY CONDITIONS FOR ELECTRON ENERGY EQUATION --
      select case (BCENE(IB))
      case (0)
        continue
  
* -- BCENE=1 -- PRESCRIBE THE VALUE OF THE ELECTRON TEMPERATURE

      case (1)
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

* -- BCENE=2 -- PRESCRIBE THE GRADIENT OF THE ELECTRON TEMPERATURE

      case (2)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'electron temperature gradient ',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
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

* -- BCENE=3 -- SHEATH CONDITIONS, ELECTRON ENERGY TRANSMISSION

      case (3)
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
              CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))  !srv 01.02.09
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
                  UUTEST=abs(
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
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))
                  UUOUT=MAX( MOMPAR(IS,IB,1)*CS*
     .             abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .                 wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .             UUTEST)
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
              CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))  !srv 01.02.09
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
                  UUTEST=abs(
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
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))
                  UUOUT=MAX( MOMPAR(IS,IB,1)*CS*
     .             abs(wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .                 wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .             UUTEST)
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

* -- BCENE=4 -- PRESCRIBE THE VALUE OF THE ELECTRON TEMPERATURE,
*                 WEAKLY A MIXED BOUNDARY CONDITION

      case (4)
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

* -- BCENE=5 -- PRESCRIBE THE ELECTRON ENERGY FLUX PER UNIT AREA

      case (5)
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

* -- BCENE=6 -- PRESCRIBE THE TOTAL ELECTRON ENERGY FLUX FOR A
*                 CONSTANT ELECTRON TEMPERATURE

      case (6)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'total electron energy flux ',ENEPAR(IB,1),
     1     ' with const Te on ',BCCHAR(IB),boundary_location(ib)
        endif
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
     .          bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))/2.0_R8
              if(mdf_fhe.eq.0) then
                flo = flo+2.5_R8*flo_vedia+1.5_R8*flo_veecrb
              else
                flo = flo+1.5_R8*flo_veecrb
              endif
            endif
            if (no_hybr_e.eq.0) then
              cone0 =
     .         hybr(flo,chce(bc_list_x(ix,ib),bc_list_y(ix,ib),1))
              fy=flo*(te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &         Tguess)/2.0_R8-cone0*(Tguess-
     &         te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .            bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))+
     &         alfTeEh*calf(bc_list_x(ix,ib),bc_list_y(ix,ib),1)*ehy*
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
     &         alfTeEh*calf(bc_list_x(ix,ib),bc_list_y(ix,ib),1)*
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
     &          +ne(bc_list_x(ix,ib),bc_list_y(ix,ib)))/2.0_R8                   !xpb
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
     .            topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))/2.0_R8
              if(mdf_fhe.eq.0) then
                flo = flo+2.5_R8*flo_vedia+1.5_R8*flo_veecrb
              else
                flo = flo+1.5_R8*flo_veecrb
              endif
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
     &         alfTeEh*calf(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)*ehy*
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
     .                 Tguess) + alfTeEh *
     &         calf(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)*
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

* -- BCENE=7 -- PRESCRIBE THE ELECTRON TEMPERATURE AS A FUNCTION OF
*                 OTHER PLASMA PARAMETERS

      case (7)
        CALL XERRAB ('B2STBC -- BCENE=7, NOT YET AVAILABLE')

* -- BCENE=8 -- PRESCRIBE THE TOTAL ELECTRON HEAT FLUX WITH CONSTANT
*                 FLUX DENSITY

      case (8)
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

* -- BCENE=9 -- PRESCRIBE THE DECAY LENGTH FOR THE ELECTRON TEMPERATURE

      case (9)
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
            if (bccon(is,ib).eq.1) then
              write(0,*)
     &         'Should use BCENE=19 in conjunction with BCCON=1'
            end if
            if (ism.gt.0) then
              write(*,'(a)')
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

* -- BCENE=10 -- Feedback option for core

      case (10)
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

* -- BCENE=11 -- FROM B2STBC_SPB                                       !srv 01.02.09

      case (11)
        if(ncall.eq.0) then   !srv added accumulation in order to account twice b.c. for corner cells
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     ' from b2stbc_spb with ENEPAR(,1) =',ENEPAR(IB,1),
     1     ' and ENEPAR(,2) =',ENEPAR(IB,2),
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
     &         topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1))          !srv 28.10.02
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
          endif                                                      !srv 21.03.09 }
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &       she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)+ENEPAR(IB,1)*
     &       max(0.0e0_R8,floe_noc(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)) !srv 28.10.02
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)=
     &       she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)+ENEPAR(IB,2)*
     &       vte*s0
!     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
!     2          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=11, NOT AVAILABLE FOR W,E')
        ENDIF

* -- BCENE=12 -- SHEATH CONDITIONS, EL ENERGY TRANSMISSION COEFFICIENT

      case (12)
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

* -- BCENE=13 -- PRESCRIBE THE ELECTRON ENERGY FLUX PER UNIT AREA
*             -- PROPORTIONAL TO TEMPERATURE

      case (13)
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

* -- BCENE=14 -- leakage option for electron energy

      case (14)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'electron energy leakage ',ENEPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            SHE0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            CS=SQRT(bcint(PZ,nx,ny,IX,IB)/bcint(RZ,nx,ny,IX,IB))   !srv 01.02.09
            SHE0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)=
     1       vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     2       hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*cs*ENEPAR(IB,1)
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            SHE0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            CS=SQRT(bcint(PZ,nx,ny,IX,IB)/bcint(RZ,nx,ny,IX,IB))   !srv 01.02.09
            SHE0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)=
     1       vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     2       hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*cs*ENEPAR(IB,1)
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=14, NOT AVAILABLE FOR W,E')
        ENDIF

* -- BCENE=15 -- FROM B2STBC_SPB                                        !srv 01.02.09

      case (15)
        if(bcene_15_style.eq.0) then
          if(ncall.eq.0) then   !srv added accumulation in order to account twice b.c. for corner cells
            write(*,'(a,1p,g14.7,a,a,a)')
     .       'bcene_15_style = 0 with GAMMAE =',GAMMAE,
     1       ' on ',BCCHAR(IB),boundary_location(ib)
          endif
          IF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)                                 !srv 20.09.17 {
!            CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))
              CS=SQRT(PZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &                RZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              if(bccon14(ib)) then                                      !srv 20.09.17 {
*    .. apply bcene=15 according to the form of bccon=14
                t0 = 0.0_R8
                do is = 0, ns-1
                  t0 = t0 + rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*
     *            CONPAR(IS,IB,1)*
     *            na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
                enddo
              else
               t0 = ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              endif                                                     !srv 20.09.17
              t0 = max(0.0e0_R8,                                        !srv 05.05.09
     &          (abs(pbs(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &           rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0))*t0*cs + !srv 25.12.06 11.01.13 20.09.17
     &           fch_p(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &              rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)/qe))  !srv 25.12.06 11.03.09 }
              if(pot_eq.eq.1) then                                      !srv 27.02.13 {
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-             !srv 11.03.09 {
     &                (1.0e0_R8+GAMMAE)/(1.0e0_R8-GAMMAE)*t0            !srv 19.04.13
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &           qe*po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*t0            !srv 07.11.12
              else                            !srv 27.02.13
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-             !srv 11.03.09 {
     &                ((1.0e0_R8+GAMMAE)/(1.0e0_R8-GAMMAE)+3.1_R8)*t0   !srv 19.04.13
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)
              endif                                                     !srv 27.02.13 }
              if(mdf_fhe.ne.0) then                                     !srv 01.06.06 {
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &           fhePSch(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                   rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)
            endif                                                       !srv 01.06.06 }
            END DO
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)                                 !srv 20.09.17
!              CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))
              CS=SQRT(PZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &                RZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              if(bccon14(ib)) then                                      !srv 20.09.17 {
*    .. apply bcene=15 according to the form of bccon=14
                t0 = 0.0_R8
                do is = 0, ns-1
                  t0 = t0+rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*
     &            CONPAR(IS,IB,1)*
     &            na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
                enddo
              else
                t0 = ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              endif                                                     !srv 20.09.17 }
              t0 = max(0.0e0_R8,                                        !srv 05.05.09
     &          (abs(pbs(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))*t0*cs -  !srv 25.12.06 11.01.13 20.09.17
     &            fch_p(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/qe))       !srv 25.12.06 11.03.09 }
              if(pot_eq.eq.1) then                                      !srv 27.02.13 {
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-             !srv 11.03.09 {
     &                (1.0e0_R8+GAMMAE)/(1.0e0_R8-GAMMAE)*t0            !srv 19.04.13
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &           qe*po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*t0            !srv 07.11.12
              else
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-             !srv 11.03.09 {
     &               ((1.0e0_R8+GAMMAE)/(1.0e0_R8-GAMMAE)+3.1_R8)*t0    !srv 19.04.13
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)
              endif                                                     !srv 27.02.13 }
              if(mdf_fhe.ne.0) then                                     !srv 01.06.06 {
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) +
     &           fhePSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)
            endif                                                       !srv 01.06.06 }
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
            DO IY = 1, BC_LIST_SIZE(IB)                                 !srv 20.09.17
!              CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))
              s1=abs(pbs(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                   rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)) !srv 18.04.13 {
              vte=sqrt(te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/me)
              if(bccon14(ib)) then                                      !srv 20.09.17 {
*    .. apply bcene=15 according to the form of bccon=14
                t0 = 0.0_R8
                do is = 0, ns-1
                  t0 = t0 + rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*
     *            CONPAR(IS,IB,1)*
     *            na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
                enddo
              else
                t0 = ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              endif                                                     !srv 20.09.17
              t0 =s1*t0*sqrt(1.0_R8/(2.0_R8*pi))*vte*                   !srv 20.09.17
     &          expu2(max(-50.0_R8,
     &         -qe*(po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))-potpar(ib,2))/
     &              te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
              if(pot_eq.eq.1) then                                      !srv 27.02.13 {
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-
     &            (1.0e0_R8+GAMMAE)*t0
                 she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &            (1.0e0_R8-GAMMAE)*
     &            qe*po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*t0
              else                                                      !srv 27.02.13
                cs=sqrt(pz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &                  rz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-
     &               s1*((1.0e0_R8+GAMMAE)+(1.0e0_R8-GAMMAE)*3.1_R8)*
     &           ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*cs
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)              !srv 18.04.13 }
              endif                                                     !srv 27.02.13 }
              if(mdf_fhe.ne.0) then                                     !srv 01.06.06 {
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &           fhePSch(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                   rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)
            endif                                                       !srv 01.06.06 }
            END DO
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)                                 !srv 20.09.17
!              CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))
              s1=abs(pbs(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))          !srv 18.04.13 {
              vte=sqrt(te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/me)
              if(bccon14(ib)) then                                      !srv 20.09.17 {
*    .. apply bcene=15 according to the form of bccon=14
                t0 = 0.0_R8
                do is = 0, ns-1
                  t0 = t0+rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*
     &            CONPAR(IS,IB,1)*
     &            na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
                enddo
              else
                t0 = ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              endif                                                     !srv 20.09.17 }
              t0 =s1*t0*sqrt(1.0_R8/(2.0_R8*pi))*vte*                   !srv 20.09.17
     &          expu2(max(-50.0_R8,
     &         -qe*(po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))-potpar(ib,2))/
     &              te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
              if(pot_eq.eq.1) then                                      !srv 27.02.13 {
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-
     &               (1.0e0_R8+GAMMAE)*t0
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &           (1.0e0_R8-GAMMAE)*qe*
     &           po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*t0
              else                               !srv 27.02.13
                cs=sqrt(pz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &                  rz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-
     &               s1*((1.0e0_R8+GAMMAE)+(1.0e0_R8-GAMMAE)*3.1_R8)*
     &           ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*cs
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)              !srv 18.04.13 }
              endif                                                     !srv 27.02.13 }
              if(mdf_fhe.ne.0) then                                     !srv 01.06.06 {
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) +
     &           fhePSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)           !srv 20.09.17 }
            endif                                                       !srv 01.06.06 }
            end do
          ELSE
            CALL XERRAB
     .       ('B2STBC -- BCENE=15, not available for S and N')
          END IF
        else
          CALL XERRAB
     .     ('B2STBC -- BCENE=15, incorrect bcene_15_style value')
        endif

* -- BCENE=16 -- Feedback boundary condition with constant temperature

      case (16)

        if(bcene_16_style.eq.0) then

          if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .       'electron temperature feedback by flux ',
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

        elseif(bcene_16_style.eq.1)then

          if(ncall.eq.0) then
            write(*,*)
     &       'Given integrated electron heat flux with constant Te '
!          write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
              if(core_boundary(ibw).and.BCENE(IBW).EQ.BCENE(IB)) THEN
                write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                    boundary_location(ibw)
                call xertst(enepar(ibw,1).ne.0.0_R8,
     .           'faulty parameter enepar1')
              endif
            END DO
          endif
          IF(core_boundary(ib)) THEN
            if(.not.once_bcene) then
              if (use_astra.ne.0) then
#ifdef ASTRA
                SENEPAR = fheb_astra
#else
                call xerrab(
     .           'Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SENEPAR = 0.0_R8
                DO IBW=IB,NBC
                  if(core_boundary(ibw).and.BCENE(IBW).EQ.BCENE(IB))THEN
                    SENEPAR = SENEPAR + ENEPAR(IBW,1)
                  endif
                END DO
              endif
              write(*,*) 'BCENE=16: Summed value is ', SENEPAR
              TES = 0.0_R8
              volsum = 0.0_R8
              DO IBW=IB,NBC
                if(core_boundary(ibw).and.BCENE(IBW).EQ.BCENE(IB))THEN
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
              totflux = 0.0_R8
              DO IBW=IB,NBC
                IF(core_boundary(ibw).AND.BCENE(IBW).EQ.BCENE(IB)) THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    if(mdf_fhe.ne.0) then
                      totflux = totflux +
     &           fhe_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
     &                +fac_ExB(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*   !srv 11.04.13 {
     &                    ((po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))-
     &               po(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))/
     &                     (hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &               hx(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))+
     &             (po(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))-
     &                      po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))/
     &             (hx(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))+
     &                      hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))*
     &                      ne(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                      te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                      SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                      bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),2)/
     &                      bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),3)**2 !srv 11.04.13 }
                    else
                      totflux = totflux +
     .             fhe(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
                    endif
                  ENDDO
                ENDIF
              ENDDO
              TES = TES *
     &         (1.0_R8-bc_ref_te*(totflux-SENEPAR)/                     !srv 11.04.13
     &           max(abs(totflux),SENEPAR))
              write(*,'(3(a,es14.7))')
     &         'TES=',TES/ev,' totflux=',totflux,
     &         ' SENEPAR=',SENEPAR
              DO IBW=IB,NBC
                IF(core_boundary(ibw).AND.BCENE(IBW).EQ.BCENE(IB)) THEN
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
            endif                                                      !srv 16.02.12 }
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCENE=16, ONLY AVAILABLE FOR CORE!')
          ENDIF
        else
          CALL XERRAB
     1       (' Incorrect style of boundary condition')
        endif

* -- BCENE=17 -- Feedback boundary condition with constant shared temperature

      case (17)
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
            do bc_type16_i=1,bc_type16_ref    ! how many steps away from the boundary?
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
     .           fhm(bc_type16_ix,bc_type16_iy,1,is)
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

* -- BCENE=18 -- FRACTIONAL DROP CONDITION (FOR EDGE2D COMPARISON)

      case (18)
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

* -- BCENE=19 -- PRESCRIBE THE DECAY LENGTH FOR THE ELECTRON TEMPERATURE

      case (19)
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

* -- BCENE=20 -- Feedback boundary condition with constant temperature

      case (20)
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
            do bc_type16_i=1,bc_type16_ref    ! how many steps away from the boundary?
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
            do bc_type20_i=1,bc_type20_ref    ! how many steps away from the boundary?
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
            do bc_type20_i=1,bc_type20_ref    ! how many steps away from the boundary?
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

* -- BCENE=21 -- constant temperature scaled by temperature on the ring bc_type21_ref away

      case (21)
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
            do bc_type21_i=1,bc_type21_ref    ! how many steps away from the boundary?
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
            do bc_type21_i=1,bc_type21_ref    ! how many steps away from the boundary?
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
            do bc_type21_i=1,bc_type21_ref    ! how many steps away from the boundary?
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
            do bc_type21_i=1,bc_type21_ref    ! how many steps away from the boundary?
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

* -- BCENE=22 -- RADIAL LEAKAGE CONDITION FOR THE ELECTRON TEMPERATURE   !srv 23.07.10
*                SPECIFIES A TEMPERATURE GRADIENT SUCH THAT THE DIFFUSIVE
*                FLUX IS SET TO BE EQUIVALENT TO A LEAKAGE VELOCITY, IN
*                UNITS OF THE ELECTRON THERMAL VELOCITY

      case (22)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENE=22 with ENEPAR(IB,1)= ',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
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
        END IF                                                         !srv 23.07.10 }

      case default
        write(*,*) BCENE(IB)
        call xerrab ('BCENE out of bounds')
      end select

      return
      end subroutine b2stbc_phys_ene

      subroutine b2stbc_phys_eni(IB)
      integer ib
      integer ix, iy, is, ibw
      integer bc_type16_i, bc_type20_i, bc_type21_i

* -- BOUNDARY CONDITIONS FOR ION ENERGY EQUATIONS --
      select case (BCENI(IB))
      case (0)
        continue

* -- BCENI=1 -- PRESCRIBE THE VALUE OF THE ION TEMPERATURE

      case (1)
        if(ncall.eq.0) then
          IF(BCCHAR(IB).EQ.'S'.and.use_astra.ne.0) THEN
            ENIPAR(IB,1) = tib_astra*1.0e3_R8                          !srv 14.06.12
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

* -- BCENI=2 -- PRESCRIBE THE GRADIENT OF THE ION TEMPERATURE

      case (2)
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

* -- BCENI=3 -- SHEATH CONDITIONS, ION ENERGY TRANSMISSION COEFFICIENT

      case (3)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,g14.7,a,g14.7,a,a,a)')
     .     'sheath ion with MOMPAR(ismain,IB,2) =',
     1     MOMPAR(ismain,IB,2),' and MOMPAR(ismain,IB,1) =',
     2     MOMPAR(ismain,IB,1),
     .     ' and ENIPAR(IB,1) =',enipar(ib,1),               !srv 24.02.10 {
     1     ' and ENIPAR(IB,2) =',enipar(ib,2),
     2     ' on ',BCCHAR(IB),boundary_location(ib)           !srv 24.02.10 }
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
              CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))  !srv 01.02.09
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
                  UUTEST=abs(
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
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))
                  UUOUT=MAX(MOMPAR(IS,IB,1)*CS*
     .             abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .                 wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .             UUTEST)
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
            if(mdf_fhi.ne.0) then                                        !srv 12.09.08 {
              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &         shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &         fhiPSch(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                 rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)
            endif                                                        !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
              CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))       !srv 01.02.09
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
                  UUTEST=abs(
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
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))
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
            if(mdf_fhi.ne.0) then                                        !srv 12.09.08 {
              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &         shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) +
     &         fhiPSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)
            endif                                                        !srv 12.09.08 }
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=3, NO SHEATH IN RADIAL DIR.')
        END IF

* -- BCENI=4 -- PRESCRIBE THE VALUE OF THE ION TEMPERATURE, WEAKLY
*                 A MIXED BOUNDARY CONDITION

      case (4)
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

* -- BCENI=5 -- PRESCRIBE THE ION ENERGY FLUX PER UNIT AREA

      case (5)
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

* -- BCENI=6 -- PRESCRIBE THE TOTAL ION ENERGY FLUX FOR A
*                 CONSTANT ION TEMPERATURE

      case (6)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'total ion energy flux ',ENIPAR(IB,1),
     1     ' with const Ti on ',BCCHAR(IB),boundary_location(ib)
        endif
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
     &          ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0))/2.0_R8
              if (mdf_fhi.eq.0) then
                flo = flo + 2.5_R8*cum_flo_vadia+1.5_R8*cum_flo_vaecrb
              else
                flo = flo + 1.5_R8*cum_flo_vaecrb
              endif
            endif
            w1 = 0.0_R8
            if (BoRiS.eq.1.0_R8) then
              do is = 0, ns-1
                t2 = dim ( coni0,
     .           abs(fhm(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is))*
     .           max(0.0_R8,
     .            min(ti(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))))/
     .           max(cutlo,ti(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .                ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))))
                w1 = w1 +
     .            upwind(fhm(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is),
     .             ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     .             ti(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .             - t2 *
     .            (ti(bc_list_x(ix,ib),bc_list_y(ix,ib)) -
     .             ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
              end do
            end if
            if (no_hybr_i.eq.0) then
              coni0 = hybr(flo,
     .         chci(bc_list_x(ix,ib),bc_list_y(ix,ib),1))
              fy = w1 + flo*
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
              fy = w1 +
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
     &          ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0))/2.0_R8
              if (mdf_fhi.eq.0) then
                flo = flo + 2.5_R8*cum_flo_vadia+1.5_R8*cum_flo_vaecrb
              else
                flo = flo + 1.5_R8*cum_flo_vaecrb
              endif
            endif
            w1 = 0.0_R8
            if (BoRiS.eq.1.0_R8) then
              do is = 0, ns-1
                t2 = dim ( coni0,
     .           abs(fhm(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,is))
     .           *max(0.0_R8,
     .            min(ti(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))))/
     .           max(cutlo,ti(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .                ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))))
                w1 = w1 +
     .            upwind(fhm(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                       topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                       1,is),
     .             ti(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .             - t2 *
     .            (ti(bc_list_x(ix,ib),bc_list_y(ix,ib)) -
     .             ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
              end do
            end if
            if (no_hybr_i.eq.0) then
              coni0 = hybr(flo,
     &         chci(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1))
              fy = w1 +
     .          flo*(ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
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
              fy = w1 +
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

* -- BCENI=7 -- PRESCRIBE THE ION TEMPERATURE AS A FUNCTION OF OTHER
*                 PLASMA PARAMETERS

      case (7)
        CALL XERRAB ('B2STBC -- BCENI=7, NOT YET AVAILABLE')

* -- BCENI=8 -- PRESCRIBE THE TOTAL ION HEAT FLUX WITH CONSTANT
*                 FLUX DENSITY

      case (8)
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
            if(mdf_fhi.ne.0) then                                 !srv 12.09.08 {
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &         fhiPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) + 1.5e0_R8/qe*
     &         (fchinert(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                   topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)+
     &          fchvispar(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)+
     &          fchanml(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)+
     &          vis_per*                                          !srv 02.10.17
     &          fchvisper(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)+
     &          vis_q*                                            !srv 02.10.17
     &          fchvisq(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1))*
     &          ycur
              do is=0,ns-1
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     &           1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                  hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &           vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)*
     &           na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
              end do
            endif                                                 !srv 12.09.08 }
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
            if(mdf_fhi.ne.0) then                                 !srv 12.09.08 {
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     .         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &         fhiPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) - 1.5e0_R8/qe*
     &         (fchinert(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)+
     &          fchvispar(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)+
     &          fchanml(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)+
     &          vis_per*fchvisper(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)+   !srv 02.10.17
     &          vis_q*fchvisq(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1))*ycur  !srv 02.10.17
              do is=0,ns-1
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     .           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) -
     &           1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                  hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &           vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)*
     &               na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
              end do
            endif                                                 !srv 12.09.08 }
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

* -- BCENI=9 -- PRESCRIBE THE DECAY LENGTH FOR THE ION TEMPERATURE

      case (9)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'ion temperature decay length ',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(ENIPAR(IB,1).ne.0.0_R8,
     1     'BCENI = 9, ENIPAR(,1) = 0')
          decay_length_ok = .true.
          do is = 0, ns-1
            decay_length_ok = decay_length_ok .and.
     1       (bccon(is,ib).eq.9 .or. bccon(is,ib).eq.0 .or.
     2       (is_neutral(is) .and. use_eirene.eq.1))
            if (bccon(is,ib).eq.1) then
              write(0,*)
     &         'Should use BCENI=19 in conjunction with BCCON=1'
            end if
          end do
          if (.not.decay_length_ok) then
!          call xertst(all(BCCON(0:NS-1,IB).EQ.9),
!     1     'Should use BCCON(IS,IB).EQ.9 together with BCENI=9 ')       !srv 16.07.12
            write(0,*)
     &       'Should use BCENI.EQ.9 together with BCCON=9'
            write(*,*)
     &       'Should use BCENI.EQ.9 together with BCCON=9'
            write(*,'(a)')
     &                 'Please ensure that this boundary condition '//
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
     .           SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)    !srv 23.09.08
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     .           2.5_R8*
     .           SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*   !srv 23.09.08
     &                    na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
cxpb          endif
            ENDDO
            if(mdf_fhi.ne.0) then                                       !srv 12.09.08 {
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
            endif                                                       !srv 12.09.08 }
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
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)      !srv 23.09.08
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     .           2.5_R8*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*     !srv 23.09.08
     &                  na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            ENDDO
            if(mdf_fhi.ne.0) then                                       !srv 12.09.08 {
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
            endif                                                       !srv 12.09.08 }
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=9, NOT AVAILABLE FOR W,E')
        ENDIF

* -- BCENI=10 -- Feedback option for core

      case (10)
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

* -- BCENI=11 -- SHEATH CONDITIONS, ION ENERGY TRANSMISSION COEFFICIENT

      case (11)
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

* -- BCENI=12 -- SHEATH CONDITIONS, ION ENERGY TRANSMISSION COEFFICIENT

      case (12)
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

* -- BCENI=13 -- PRESCRIBE THE ION ENERGY FLUX PER UNIT AREA
*             -- PROPORTIONAL TO TEMPERATURE

      case (13)
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

* -- BCENI=14 -- leakage option for ion energy

      case (14)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'ion energy leakage ',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            CS=SQRT(bcint(PZ,nx,ny,IX,IB)/bcint(RZ,nx,ny,IX,IB))  !srv 01.02.09
            SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)=
     1       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))/
     .       hy1(bc_list_x(ix,ib),bc_list_y(ix,ib))*cs*ENIPAR(IB,1)
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            CS=SQRT(bcint(PZ,nx,ny,IX,IB)/bcint(RZ,nx,ny,IX,IB))  !srv 01.02.09
            SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)=
     1       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))/
     .       hy1(bc_list_x(ix,ib),bc_list_y(ix,ib))*cs*ENIPAR(IB,1)
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=14, NOT AVAILABLE FOR W,E')
        ENDIF

* -- BCENI=15 -- from b2stbc_spb                              !srv 01.02.09

      case (15)
        if(ncall.eq.0) then
          call xertst(ENIPAR(IB,1).gt.0.0_R8,
     1       'BCENI = 15, ENIPAR(I,,1) <= 0 not allowed!')              !srv 20.09.17
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     &     'from b2stbc_spb for ion temperature with ENIPAR(IB,1) =',
     &     ENIPAR(IB,1),' on ', BCCHAR(IB),boundary_location(ib)        !srv 20.09.17
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
!            SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=0.0_R8
!            CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))
            CS=SQRT(PZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &              RZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            if(bccon14(ib)) then                                        !srv 20.09.17 {
*    .. apply bceni=15 according to the form of bccon=14
              t0 = 0.0_R8
              do is = 0, ns-1
                if(.not.is_neutral(is)) t0 = t0 + CONPAR(IS,IB,1)*      !srv 20.09.17
     &           na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
              enddo
            else                                                        !srv 20.09.17 }
              t0 = 0.0_R8
              do is = 0, ns-1
                if(.not.is_neutral(is)) t0 = t0 +
     &            na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)              !srv 11.03.09
              enddo
            endif                                                       !srv 20.09.17
            SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &       SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) -                !srv 11.03.09 20.09.17
     &       ENIPAR(IB,1)*t0*cs*
     &       abs(pbs(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &               rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0))     !srv 26.11.02 28.01.03
            if(mdf_fhi.ne.0) then                                       !srv 13.10.06 {
              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &         shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &         fhiPSch(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                 rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)
            endif                                                       !srv 13.10.06 }
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
!            SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=0.0_R8
!            CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))
            CS=SQRT(PZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     .              RZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            if(bccon14(ib)) then                                        !srv 20.09.17 {
*    .. apply bceni=15 according to the form of bccon=14
              t0 = 0.0_R8
              do is = 0, ns-1
                if(.not.is_neutral(is)) t0 = t0 +CONPAR(IS,IB,1)*       !srv 20.09.17
     .            na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
              enddo
            else                                                        !srv 20.09.17 }
              t0 = 0.0_R8
              do is = 0, ns-1
                if(.not.is_neutral(is)) t0 = t0 +
     .            na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)              !srv 11.03.09
              enddo
            endif                                                       !srv 20.09.17
            SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     .       SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) -                !srv 11.03.09 20.09.17
     &       ENIPAR(IB,1)*t0*cs*
     &       abs(pbs(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))              !srv 26.11.02 28.01.03
            if(mdf_fhi.ne.0) then                                       !srv 13.10.06 {
              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     .         shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) +
     &         fhiPSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)
            endif                                                        !srv 13.10.06 }
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=15, NOT AVAILABLE FOR S,N')
        ENDIF

* -- BCENI=16 -- Feedback boundary condition with constant temperature

      case (16)
        if(bceni_16_style.eq.0) then
         if(ncall.eq.0) then
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .       'ion temperature feedback by flux ',
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
     .            fhm(bc_type16_ix,bc_type16_iy,1,is)
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

        elseif(bceni_16_style.eq.1)then                                 !srv 16.02.12 {

         if(ncall.eq.0) then
           write(*,*) 'Given integrated ion heat flux',
     &                 ' with constant Ti '
!          write all boundaries with this type of boundary condition
           DO IBW=IB,NBC
             if(core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB))THEN
               write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                   boundary_location(ibw)
               call xertst(enipar(ibw,1).ne.0,
     &           'faulty parameter enipar1')
             endif
             if(core_boundary(ibw).and.
     &         (BCENI(IBW).EQ.23.or.BCENI(IBW).EQ.24.or.
     &          BCENI(IBW).EQ.26.or.BCENI(IBW).EQ.27)) then  ! IYS added BCENI= 26 and 27
               call xerrab
     .          ('BCENI types 16, 23, 24, 26 and 27 cannot be mixed!')
             endif
           ENDDO
         endif
         IF(core_boundary(ib)) THEN
           if(.not.once_bceni) then
             if (use_astra.ne.0) then
#ifdef ASTRA
               SENIPAR = fhib_astra
#else
               call xerrab
     &          ('Compile with -DASTRA option to couple to ASTRA!')
#endif
             else
               SENIPAR = 0.0_R8
               DO IBW=IB,NBC
                 if(core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB))THEN
                   SENIPAR = SENIPAR + ENIPAR(IBW,1)
                 endif
               ENDDO
             endif
             write(*,*) 'BCENI=16: Summed value is ', SENIPAR
             TIS = 0.0_R8
             NIVS = 0.0_R8
             DO IBW=IB,NBC
               if(core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB))THEN
                 DO IX = 1, BC_LIST_SIZE(IBW)
                   TIS = TIS +
     .                   TI(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     .                   NI(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     .            VOL(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                   NIVS = NIVS +
     .                   NI(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     .            VOL(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                 END DO
               endif
             END DO
             TIS = TIS/NIVS
             totflux = 0.0_R8
             DO IBW=IB,NBC
               IF(core_boundary(ibw).AND.BCENI(IBW).EQ.BCENI(IB)) THEN
                 DO IX = 1, BC_LIST_SIZE(IBW)
                   if(mdf_fhi.ne.0) then
                      totflux = totflux +
     &          fhi_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
     &               +fac_ExB(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*    !srv 11.04.13 {
     &                   ((po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))-
     &              po(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))/
     &                    (hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &              hx(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))+
     &            (po(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))-
     &                     po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))/
     &            (hx(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))+
     &                     hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))*
     &                     ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)*
     &                     ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                     SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                     bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),2)/
     &                     bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),3)**2 !srv 11.04.13 }
                   else
                     totflux = totflux +
     &             fhi(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
                   endif
                 ENDDO
               ENDIF
             ENDDO
             TIS = TIS *
     &         (1.0_R8-bc_ref_ti*(totflux-SENIPAR)/                     !srv 11.04.13
     &           max(abs(totflux),SENIPAR))
             write(*,'(3(a,es14.7))')
     &         'TIS=',TIS/ev,' totflux=',totflux,
     &         ' SENIPAR=',SENIPAR
             DO IBW=IB,NBC
               IF(core_boundary(ibw).AND.BCENE(IBW).EQ.BCENE(IB)) THEN
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
           endif                                                       !srv 16.02.12 }
         ELSE
           CALL XERRAB
     1         ('B2STBC -- BCENI=16, ONLY AVAILABLE FOR CORE!')
         ENDIF
        else
         CALL XERRAB (' Incorrect style of boundary condition')
        endif

* -- BCENI=17 -- Feedback boundary condition with constant temperature

      case (17)
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
          CALL XERRAB ('B2STBC -- BCENI=17, NOT AVAILABLE FOR N,W,E')
        ENDIF

* -- BCENI=18 -- FRACTIONAL DROP CONDITION (FOR EDGE2D COMPARISON)

      case (18)
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

* -- BCENI=19 -- PRESCRIBE THE DECAY LENGTH FOR THE ION TEMPERATURE

      case (19)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'ion temperature decay length (new) ',ENIPAR(IB,1),' on ',
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

* -- BCENI=20 -- Feedback boundary condition with constant temperature scaled by temperature on the ring bc_type20_ref away

      case (20)
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

!! FIXME: Clash here with 5.0/5.1 BCENI=21 type
cxpb Careful here because of BC number clash
* From 5.0: starting at version 01.001.030
* -- BCCON=21 -- constant ion temperature scaled by ion temperature on the ring bc_type21_ref away
* has been moved up to next available number (BCENI=25)

* From 5.2, after version number 01.001.078:
* -- BCENI=21 -- from b2stbc_spb                                     !srv 01.02.09

      case (21)
        if(ncall.eq.0) then
          if (boundary_version.ge.'01.001.030' .and.
     &        boundary_version.lt.'02.000.000') then
            write(*,*) 'BCENI=21 type B.C. moved to 25 ',
     &       ' because of clash with 5.2 numbering. '
            write(*,*) 'Please modify b2.boundary.parameters file'
            call xerrab('B2STBC -- BCENI=25, DISPLACED')
          endif
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     &     'from b2stbc_spb for BCENI with ENIPAR(IB,1) =',
     &     ENIPAR(IB,1),' and ENIPAR(IB,2) =', ENIPAR(IB,2),' on ',
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
          endif                                                      !srv 21.03.09 }
            do is = 0, ns-1
              if(.not.is_neutral(is)) then                           !srv 19.03.09 {
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
     &        topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1))           !srv 28.10.02
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
              if(.not.is_neutral(is)) then                           !srv 19.03.09 }
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
              endif                                                  !srv 19.03.09 {
            enddo
            SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &       SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &       max(0.0e0_R8,
     &        floi_noc(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1))         !srv 28.10.02
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=21, NOT AVAILABLE FOR W,E')
        ENDIF

* -- BCENI=22 -- RADIAL LEAKAGE CONDITION FOR THE ION TEMPERATURE   !srv 23.07.10
*                SPECIFIES A TEMPERATURE GRADIENT SUCH THAT THE DIFFUSIVE
*                FLUX IS SET TO BE EQUIVALENT TO A LEAKAGE VELOCITY, IN
*                UNITS OF THE COLLECTIVE ION THERMAL VELOCITY

      case (22)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     &     'BCENI=22 with ENIPAR(IB,1)= ',ENIPAR(IB,1),
     &     ' on ', BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            t0 = 0.0_R8
            do is = 0, ns-1
              if(.not.is_neutral(is)) then                           !srv 19.03.09 {
                csb = sqrt((rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &                       te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &                       ti(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/
     &                     (am(is)*mp))
                t0 = t0 + csb*na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
              endif                                                  !srv 19.03.09 }
            enddo
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*EV*
     &       CHCI(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &            TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)*
     &         TI(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &            TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*EV*
     &       (CHCI(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &             TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)-
     &       ENIPAR(IB,1)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*T0)
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            t0 = 0.0_R8
            do is = 0, ns-1
              if(.not.is_neutral(is)) then                           !srv 19.03.09 {
                csb = sqrt((rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &                       te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &                       ti(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/
     &                     (am(is)*mp))
                t0 = t0 + csb*na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
              endif                                                  !srv 19.03.09 }
            enddo
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*EV*
     &       CHCI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)*
     &       TI(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*EV*
     &       (CHCI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-ENIPAR(IB,1)* !srv 02.02.11
     &        SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &           BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*T0)
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=22, NOT AVAILABLE FOR W,E')
        ENDIF

* -- BCENI=23 -- PRESCRIBE THE POLOIDALLY AVERAGED VALUE OF THE ION TEMPERATURE
* AND INTRODUCE A POLOIDAL VARIATION AS CLOSE AS POSSIBLE TO NEOCLASSICAL SOLUTION.
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH CORRESPONDING
* CONDITION ON ION DENSITY BCCON = 21,22,23

      case (23)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7)')
     .     'specified ion temperature ',ENIPAR(IB,1)
!          write all boundaries with this type of boundary condition
          DO IBW=IB,NBC
            IF (core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB)) THEN
              write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                  boundary_location(ibw)
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
                write(*,'(a)')
     &           'Please ensure that this boundary condition '//
     &           'is consistent with '//
     &           'ion density boundary condition'  !srv 21.09.12
              end if
              if (.not.(ns.eq.1.or.(ns.eq.2.and.is_neutral(0)))) then
                write(0,*)
     &           'Use of BCENI=23 only meant for single-fluid runs!'
                write(0,*) 'Use BCENI=26 instead!'
                write(6,*)
     &           'Use of BCENI=23 only meant for single-fluid runs!'
                write(6,*) 'Use BCENI=26 instead!'
              endif
            endif
            if(core_boundary(ibw).and.
     &         (BCENI(IBW).EQ.16.or.BCENI(IBW).EQ.24.or.
     &          BCENI(IBW).EQ.26.or.BCENI(IBW).EQ.27)) then  ! IYS added BCENI= 26 and 27
                call xerrab
     .           ('BCENI types 16, 23, 24, 26 and 27 cannot be mixed!')
            endif
          ENDDO
          call xertst(use_astra.ne.0.or.ENIPAR(IB,1).gt.0.0_R8,
     1     'BCENI = 23, ENIPAR(,1) !> 0')
        endif

        if(.not.once_bceni) then
          if (use_astra.ne.0) THEN
#ifdef ASTRA
            DO IBW=IB,NBC
              if(core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB)) THEN
                ENIPAR(IBW,1) = tib_astra*1.0e3_R8                     !srv 14.06.12
              endif
            ENDDO
#else
            call xerrab
     .       ('Compile with -DASTRA option to couple with ASTRA!')
#endif
          ENDIF

          IF(core_boundary(ib)) THEN
!          DO IX = 1, BC_LIST_SIZE(IB)
!            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*ENIPAR(IB,1)
!     .       *EV*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
!            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
!     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
!          end do
            t0=0.0e0_R8   ! Volume                                     !lk 12.12.12 {
            t1=0.0e0_R8   ! Poloidally volume-averaged ion temperature in the cell next to boundary
            t2=0.0e0_R8   ! Poloidally volume-averaged electron density in the cell next to boundary
            t3=0.0e0_R8   ! Poloidally volume-averaged electron density in the boundary cell
            DO IBW=IB,NBC
              IF(core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB)) THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  t0=t0 +
     &              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t1=t1 +
     &               ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t2=t2 +
     &               ne(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t3=t3 + ne(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &             vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                END DO
              ENDIF
            ENDDO
            t1=t1/t0
            t2=t2/t0
            t3=t3/t0
            DO IBW=IB,NBC
              IF(core_boundary(ibw).AND.BCENI(IBW).EQ.BCENI(IB)) THEN
                do IX = 1, BC_LIST_SIZE(IBW)
                  shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)=
     &             shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)
     &              +VBIG*NE00*(ENIPAR(IBW,1)*EV
     &              +(ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
     &                -t1)*t1*t2/t3/ENIPAR(IBW,1)/EV)
     &              *SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)=
     &             -VBIG*NE00*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                enddo
              ENDIF                                                 !lk 12.12.12 }
            ENDDO
            once_bceni = .true.
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCENI=23 IS AVAILABLE FOR CORE BOUNDARY ONLY.')
          END IF
        endif

* -- BCENI=24 -- Feedback boundary condition with prescribed total ion heat flux,
*                constant poloidally averaged ion temperature and a
*                poloidal variation as close as possible to neoclassical solution.
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH CORRESPONDING
* CONDITION ON ION DENSITY BCCON = 21,22,23

      case (24)
!IYS 04.03.2013 {
        if(ncall.eq.0) then
          write(*,*) 'Given integrated ion heat flux',
     &               ' with constant Ti '
!          write all boundaries with this type of boundary condition
          DO IBW=IB,NBC
            if (core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB)) THEN
              write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                  boundary_location(ibw)
              call xertst(enipar(ibw,1).ne.0,'faulty parameter enipar1')
              if (.not.(bccon(ismain,ib).eq.21 .or.
     .                  bccon(ismain,ib).eq.22 .or.
     .                  bccon(ismain,ib).eq.23)) then
                  write(0,'(a,i3,a,i2,a)') 'Should use BCCON(',ismain,
     &             ',',IB,').EQ.21,22,23 together with BCENI=24'
                  write(*,'(a,i3,a,i2,a)') 'Should use BCCON(',ismain,
     &             ',',IB,').EQ.21,22,23 together with BCENI=24'
                  write(*,'(a)')
     &             'Please ensure that this boundary condition '//
     &             'is consistent with '//
     &             'ion density boundary condition'  !srv 21.09.12
              end if
              if (.not.(ns.eq.1.or.(ns.eq.2.and.is_neutral(0)))) then
               write(0,*)
     &          'Use of BCENI=24 only meant for single-fluid runs!'
               write(0,*) 'Use BCENI=27 instead!'
               write(6,*)
     &          'Use of BCENI=24 only meant for single-fluid runs!'
               write(6,*) 'Use BCENI=27 instead!'
              endif
            endif
            if(core_boundary(ibw).and.
     &         (BCENI(IBW).EQ.16.or.BCENI(IBW).EQ.23.or.
     &          BCENI(IBW).EQ.26.or.BCENI(IBW).EQ.27)) then  ! IYS added BCENI= 26 and 27
                call xerrab
     .           ('BCENI types 16, 23, 24, 26 and 27 cannot be mixed!')
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
              if(core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB)) THEN
                SENIPAR = SENIPAR + ENIPAR(IBW,1)
              endif
            ENDDO
          endif
          write(*,*) 'BCENI=24: Summed value is ', SENIPAR

          IF(core_boundary(ib)) THEN
            TIS = 0.0_R8  ! Poloidally density-averaged ion temperature in the boundary cell     !IYS 21.01.2016
            NIVS = 0.0_R8
            totflux = 0.0_R8
            t0=0.0e0_R8   ! Volume                                     !lk 12.12.12 {
            t1=0.0e0_R8   ! Poloidally volume-averaged ion temperature in the cell next to boundary
            t2=0.0e0_R8   ! Poloidally volume-averaged electron density in the cell next to boundary
            t3=0.0e0_R8   ! Poloidally volume-averaged electron density in the boundary cell
            DO IBW=IB,NBC
              IF(core_boundary(ibw).AND.BCENI(IBW).EQ.BCENI(IB)) THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  t0=t0 +
     &               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t1=t1 +
     &                ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  TIS=TIS + ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                      ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     &               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  NIVS=NIVS +
     &                     ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     &              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t2=t2 +
     &                ne(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t3=t3 + ne(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
     &              *vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  if(mdf_fhi.ne.0) then
                    totflux = totflux +
     &        fhi_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
     &             +fac_ExB(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*      !srv 11.04.13 {
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
     &                   ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)*
     &                   ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                   SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                   bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),2)/
     &                   bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),3)**2  !srv 11.04.13 }
                  else
                    totflux = totflux +
     &               fhi(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
                  endif
                ENDDO
              ENDIF
            ENDDO
            TIS=TIS/NIVS
            t1=t1/t0
            t2=t2/t0
            t3=t3/t0
            TIS=TIS*
     &          (1.0_R8-bc_ref_ti*(totflux-SENIPAR)/
     &            max(abs(totflux),SENIPAR))
            write(*,'(3(a,es14.7))')
     &        'TIS=',TIS/ev,' totflux=',totflux,
     &        ' SENIPAR=',SENIPAR
            DO IBW=IB,NBC
              IF(core_boundary(ibw).AND.BCENE(IBW).EQ.BCENE(IB)) THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)=
     .             VBIG * NE00 * (TIS +
     .             (ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))-t1)*
     .              t1*t2/t3/TIS)*
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

* -- BCENI=25 -- constant ion temperature scaled by temperature on the ring bc_type21_ref away
* was BCENI=21 in 5.0/5.1

      case (25)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'fixed ion temperature (modified by nearby ring) ',
     1     ENIPAR(IB,1),' with ',ENIPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(enipar(ib,1).ne.0.0_R8,
     1     'faulty parameter enipar1')
          call xertst(enipar(ib,2).gt.0.0_R8,
     1     'faulty parameter enipar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          tiave = 0.0_R8
          niave = 0.0_R8
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
            tiave=tiave+ti(bc_type21_ix,bc_type21_iy)*
     1       ni(bc_type21_ix,bc_type21_iy,0)*
     2       vol(bc_type21_ix,bc_type21_iy)
            niave=niave+ni(bc_type21_ix,bc_type21_iy,0)*
     1       vol(bc_type21_ix,bc_type21_iy)
          ENDDO
          tiave=tiave/niave
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
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*ENIPAR(IB,1)*EV*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .       ((1.0_R8-enipar(ib,2))+
     .           enipar(ib,2)*ti(bc_type21_ix,bc_type21_iy)/tiave)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE IF(BCCHAR(IB).EQ.'N') THEN
          tiave = 0.0_R8
          niave = 0.0_R8
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
            tiave=tiave+ti(bc_type21_ix,bc_type21_iy)*
     1       ni(bc_type21_ix,bc_type21_iy,0)*
     2       vol(bc_type21_ix,bc_type21_iy)
            niave=niave+ni(bc_type21_ix,bc_type21_iy,0)*
     1       vol(bc_type21_ix,bc_type21_iy)
          ENDDO
          tiave=tiave/niave
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
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*ENIPAR(IB,1)*EV*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .       ((1.0_R8-enipar(ib,2))+
     .           enipar(ib,2)*ti(bc_type21_ix,bc_type21_iy)/tiave)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=25, NOT AVAILABLE FOR W,E')
        ENDIF

* -- BCENI=26 -- PRESCRIBE THE POLOIDALLY AVERAGED VALUE OF THE ION TEMPERATURE
* AND INTRODUCE A POLOIDAL VARIATION IN A SIMPLIFIED MANNER.
* THIS BOUNDARY CONDITION IS SUITABLE FOR ANY PLASMA COMPOSITION (i.e. when BCENI=23 fails).
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH CORRESPONDING
* CONDITION ON ION DENSITY BCCON = 25,26,27

      case (26)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7)')
     .     'specified ion temperature ',ENIPAR(IB,1)
!          write all boundaries with this type of boundary condition
          DO IBW=IB,NBC
            IF (core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB)) THEN
              write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                  boundary_location(ibw)
              call xertst(use_astra.ne.0.or.ENIPAR(IBW,1).gt.0.0_R8,
     1         'BCENI = 26, ENIPAR(,1) !> 0')
              call xertst(use_astra.ne.0.or.
     .         ENIPAR(IB,1).eq.ENIPAR(IBW,1),
     1         'BCENI = 26, inconsistent ENIPAR(,1) values!')
            endif
            if(core_boundary(ibw).and.
     &         (BCENI(IBW).EQ.16.or.BCENI(IBW).EQ.23.or.
     &          BCENI(IBW).EQ.24.or.BCENI(IBW).EQ.27)) then  ! IYS added BCENI= 26 and 27
                call xerrab
     .           ('BCENI types 16, 23, 24, 26 and 27 cannot be mixed!')
            endif
          ENDDO
        endif

        if(.not.once_bceni) then
          if (use_astra.ne.0) THEN
#ifdef ASTRA
            DO IBW=IB,NBC
              if(core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB)) THEN
                ENIPAR(IBW,1) = tib_astra*1.0e3_R8                     !srv 14.06.12
              endif
            ENDDO
#else
            call xerrab
     .       ('Compile with -DASTRA option to couple with ASTRA!')
#endif
          ENDIF
          IF(core_boundary(ib)) THEN
            t0=0.0e0_R8   ! Volume                                     !lk 12.12.12 {
            t1=0.0e0_R8   ! Poloidally volume-averaged ion temperature in the cell next to boundary
            DO IBW=IB,NBC
              IF(core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB)) THEN
                DO IX= 1, BC_LIST_SIZE(IBW)
                  t0 = t0 +
     .              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t1 = t1 +
     .               ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     .              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                END DO
              END IF
            END DO
            t1=t1/t0
            DO IBW=IB,NBC
              IF(core_boundary(ibw).AND.BCENI(IBW).EQ.BCENI(IB)) THEN
                do IX = 1, BC_LIST_SIZE(IBW)
                  shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)=
     &             shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)
     &              +VBIG*NE00*(ENIPAR(IBW,1)*EV
     &              +(ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
     &               -t1))
     &              *SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)=
     &             -VBIG*NE00*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                enddo
              ENDIF                                                 !lk 12.12.12 }
            ENDDO
            once_bceni = .true.
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCENI=26 IS AVAILABLE FOR CORE BOUNDARY ONLY.')
          END IF
        endif

* -- BCENI=27 -- Feedback boundary condition with prescribed total ion heat flux,
*                constant poloidally averaged ion temperature and a
*                poloidal variation in a simplified manner.
* THIS BOUNDARY CONDITION IS SUITABLE FOR ANY PLASMA COMPOSITION (i.e. when BCENI=24 fails).
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH CORRESPONDING
* CONDITION ON ION DENSITY BCCON = 25,26,27

      case (27)
        if(ncall.eq.0) then
          write(*,*) 'Given integrated ion heat flux',
     &               ' with constant Ti '
!          write all boundaries with this type of boundary condition
          DO IBW=IB,NBC
            if (core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB)) THEN
              write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                  boundary_location(ibw)
              call xertst(enipar(ibw,1).ne.0,'faulty parameter enipar1')
            endif
            if(core_boundary(ibw).and.
     &         (BCENI(IBW).EQ.16.or.BCENI(IBW).EQ.23.or.
     &          BCENI(IBW).EQ.24.or.BCENI(IBW).EQ.26)) then  ! IYS added BCENI= 26 and 27
                call xerrab
     .           ('BCENI types 16, 23, 24, 26 and 27 cannot be mixed!')
            endif
          ENDDO
        endif

        IF(core_boundary(ib)) THEN
          if(.not.once_bceni) then
            if (use_astra.ne.0) then
#ifdef ASTRA
              SENIPAR = fhib_astra
#else
              call xerrab
     .         ('Compile with -DASTRA option to couple with ASTRA!')
#endif
            else
              SENIPAR = 0.0_R8
              DO IBW=IB,NBC
                if(core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB)) THEN
                  SENIPAR = SENIPAR + ENIPAR(IBW,1)
                endif
              ENDDO
            endif
            write(*,*) 'BCENI=27: Summed value is ', SENIPAR
            TIS = 0.0_R8 ! Poloidally density-averaged ion temperature in the boundary cell
            NIVS = 0.0_R8
            t0=0.0e0_R8   ! Volume                                      !lk 12.12.12 {
            t1=0.0e0_R8   ! Poloidally volume-averaged ion temperature in the cell next to boundary
            totflux = 0.0_R8   ! IYS 22.03.2017
            DO IBW=IB,NBC
              if(core_boundary(ibw).and.BCENI(IBW).EQ.BCENI(IB)) THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  TIS=TIS + ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     .                      ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  NIVS=NIVS +
     .                      ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t0=t0 +
     &               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t1=t1 +
     &                ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  if(mdf_fhi.ne.0) then
                    totflux = totflux +
     &        fhi_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
     &             +fac_ExB(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*      !srv 11.04.13 {
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
     &                   ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)*
     &                   ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                   SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                   bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),2)/
     &                   bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),3)**2  !srv 11.04.13 }
                  else
                    totflux = totflux +
     &           fhi(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1)
                  endif
                END DO
              endif
            END DO
            TIS=TIS/NIVS
            t1=t1/t0
            TIS=TIS*
     &          (1.0_R8-bc_ref_ti*(totflux-SENIPAR)/
     &            max(abs(totflux),SENIPAR))
            write(*,'(3(a,es14.7))')
     &        'TIS=',TIS/ev,' totflux=',totflux,
     &        ' SENIPAR=',SENIPAR
            DO IBW=IB,NBC
              IF(core_boundary(ibw).AND.BCENE(IBW).EQ.BCENE(IB)) THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)=
     .             VBIG*NE00*
     $             (TIS+(ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                      topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
     &                  -t1))*
     .              SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)=
     .             -VBIG*NE00*
     .              SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                ENDDO
              ENDIF
            ENDDO
            once_bceni = .true.
          endif                                                        !srv 16.02.12 }
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCENI=27 IS AVAILABLE FOR CORE ONLY')
        ENDIF

      case default
        write(*,*) BCENI(IB)
        call xerrab ('BCENI out of bounds')
      end select

      return
      end subroutine b2stbc_phys_eni

      subroutine b2stbc_phys_pot(IB)
      integer ib
      integer ix, iy, is, ibw
      integer bc_type21_i

* -- BOUNDARY CONDITIONS FOR POTENTIAL EQUATION --
      select case (BCPOT(IB))
      case (0)
        continue

* -- BCPOT=1 -- PRESCRIBE THE VALUE OF THE POTENTIAL

      case (1)
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

* -- BCPOT=2 -- PRESCRIBE THE GRADIENT OF THE POTENTIAL

      case (2)
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

* -- BCPOT=3 -- SHEATH CONDITIONS, ELECTRON ENERGY TRANSMISSION

      case (3)
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
#ifdef B25_EIRENE
*     ..set sheath potential drop for Eirene
            if (eirene_sheath_pot.eq.1) then
              ! B2 indexing for delta_sheathxb
              delta_sheathxb(BC_LIST_X(IY,IB)+1,BC_LIST_Y(IY,IB)+1) =
     .          qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
     .          bcint(te,nx,ny,iy,ib)
            endif
#endif
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
#ifdef B25_EIRENE
*     ..set sheath potential drop for Eirene
            if (eirene_sheath_pot.eq.1) then
              ! B2 indexing for delta_sheathxb
              delta_sheathxb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)+1) =
     .          qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
     .          bcint(te,nx,ny,iy,ib)
            endif
#endif
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=3, NO SHEATH IN RADIAL DIR.')
        END IF

* -- BCPOT=4 -- PRESCRIBE THE VALUE OF THE POTENTIAL,
*                 WEAKLY A MIXED BOUNDARY CONDITION

      case (4)
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

* -- BCPOT=5 -- PRESCRIBE THE CURRENT FLUX DENSITY PER UNIT AREA

      case (5)
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

* -- BCPOT=6 -- PRESCRIBE THE TOTAL CURRENT FLUX DENSITY FOR A
*                 CONSTANT POTENTIAL

      case (6)
        if(ncall.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'total current flux density',POTPAR(IB,1),
     1     ' with const potential on ',BCCHAR(IB),boundary_location(ib)
        endif
        CALL XERRAB ('B2STBC -- BCPOT=6, NOT YET AVAILABLE')

* -- BCPOT=7 -- PRESCRIBE THE POTENTIAL AS A FUNCTION OF
*                 OTHER PLASMA PARAMETERS

      case (7)
        CALL XERRAB ('B2STBC -- BCPOT=7, NOT YET AVAILABLE')

* -- BCPOT=8 -- PRESCRIBE THE TOTAL ELECTRON CURRENT WITH CONSTANT
*                 FLUX DENSITY

      case (8)
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

* -- BCPOT=9 -- PRESCRIBE THE DECAY LENGTH FOR
*                 THE POTENTIAL

      case (9)
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

* -- BCPOT=10 -- Feedback option for core
*                NOT YET TESTED!!!!!!!!!

      case (10)
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

* -- BCPOT=11 -- SHEATH CONDITIONS, CURRENT CONTINUITY; FROM b2stbc_spb !srv 01.02.09

      case (11)
        if(ncall.eq.0) then
          write(*,'(a,a,1p,g14.7,a,g14.7,a,a,a)')
     &     'sheath potential bc from b2stbc_spb',
     &     ' with POTPAR(IB,2) =',POTPAR(IB,2),' gammae =', gammae,
     &     ' on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            s1 = abs(pbs(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                   rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0))
!            vte = sqrt(bcint(te,nx,ny,iy,ib)/me)
            vte = sqrt(te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/me)
            wrong_flow=.false.
*    ..set boundary condition in explicit form (the same form as fna)
            CS=SQRT(PZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &              RZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            t0 = 0.0e0_R8
            do is = 0, ns-1
              if(istyle_fchi.eq.1 .and. BCCON(IS,IB).eq.14) then        !srv 14.02.05 {
                t0 = t0+rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*
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
#ifdef B25_EIRENE
*     ..set sheath potential drop for Eirene
            if (eirene_sheath_pot.eq.1) then
              ! B2 indexing for delta_sheathxb
              delta_sheathxb(BC_LIST_X(IY,IB)+1,BC_LIST_Y(IY,IB)+1) =
     &          qe*(po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))-potpar(ib,2))/
     &              te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            endif
#endif
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
            CS=SQRT(PZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &              RZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            t0 = 0.0e0_R8
            do is = 0, ns-1
              if(istyle_fchi.eq.1 .and. BCCON(IS,IB).eq.14) then        !srv 14.02.05 {
*    ..set boundary condition in explicit form (the same form as fna)
                t0 = t0+rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*
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
#ifdef B25_EIRENE
*     ..set sheath potential drop for Eirene
            if (eirene_sheath_pot.eq.1) then
              ! B2 indexing for delta_sheathxb
              delta_sheathxb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)+1) =
     &          qe*(po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))-potpar(ib,2))/
     &              te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            endif
#endif
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=11, NO SHEATH IN RADIAL DIR.')
        END IF

* -- BCPOT=12 -- Test condition for south core region ! Solovyev 07.04.14   !srv 03.03.15
c code here is a part of code from 50000 for istyle_cur_contr_on_S_and_N.eq.1

      case (12)
        if(ncall.eq.0) then
          write(*,'(a,a,a,a,a)') 'Adding contributions ',
     &     'from diamagnetic and inertial currents ',
     &     'on ',BCCHAR(IB),boundary_location(ib)
          if(istyle_cur_contr_on_S_and_N.eq.0 .or.
     &       istyle_cur_contr_on_S_and_N.eq.1) then
            write(*,*) ' BCPOT=12 cannot be used when ',
     &                 'istyle_cur_contr_on_S_and_N = 0 or 1'
            call xerrab
     .       ('faulty parameter istyle_cur_contr_on_S_and_N !')
          else
            call xertst(core_boundary(ib),
     .       'BCPOT=12 intended for core boundaries only!')
          endif
        endif
        if (BCCHAR(ib).eq.'S') then
          do ix=1, BC_LIST_SIZE(IB)
*        for inner surfaces
            if(maxval(facdrift) .ne. 0.0_R8) then
              fchy_dia =                                               !srv 01.07.99
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
     &         sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) + fchy_dia
              fchy_inertia = 0.0_R8                                     !srv 25.06.03
              do is = 0, ns-1                                           !srv 25.06.03
                if(.not.is_neutral(is))then                             !srv 25.06.03
                  fchy_inertia = fchy_inertia                           !srv 07.04.00 25.06.03
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
     &         sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)+fchy_inertia
            endif
          enddo
        else
          call xerrab ('B2STBC -- BCPOT=12, NOT AVAILABLE FOR N,W,E')
        endif

* -- BCPOT=13 -- Test condition for south private & north regions ! Solovyev 07.04.14
c code here is a part of code from 50000 for istyle_cur_contr_on_S_and_N.eq.1

      case (13)
        if(ncall.eq.0) then
          write(*,'(a,a,a,a,a)') 'Adding contributions ',
     &     'from all currents ',
     &     'on ',BCCHAR(IB),boundary_location(ib)
          if(istyle_cur_contr_on_S_and_N.eq.0 .or.
     &       istyle_cur_contr_on_S_and_N.eq.1) then
            write(*,*) ' BCPOT=13 cannot be used when ',
     &                 'istyle_cur_contr_on_S_and_N = 0 or 1'
            call xerrab
     .       ('faulty parameter istyle_cur_contr_on_S_and_N !')
          else
            call xertst(.not.core_boundary(ib),
     .       'BCPOT=13 not intended for core boundaries!')
          endif
        endif
        if (BCCHAR(ib).eq.'S') then
          do ix=1, BC_LIST_SIZE(IB)
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
          enddo
        elseif (BCCHAR(ib).eq.'N') then
          do ix=1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &        sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &             fchdia(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &             fchin(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &             fchvispar(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &             fchvisq(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &             fchvisper(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &             fchinert(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)
          enddo
        else
          call xerrab ('B2STBC -- BCPOT=13, NOT AVAILABLE FOR W,E')
        endif

* -- BCPOT=14 -- Nothing implemented

      case (14)
        CALL XERRAB ('B2STBC -- BCPOT=14, NOT AVAILABLE FOR S,N,W,E')

* -- BCPOT=15 -- Nothing implemented

      case (15)
        CALL XERRAB ('B2STBC -- BCPOT=15, NOT AVAILABLE FOR S,N,W,E')

* -- BCPOT=16 -- Feedback boundary condition with constant potential    !Solovyev 22.03.14 !srv 29.09.15

      case (16)
        if (ncall.eq.0) then
          write(*,'(a,a,a)')
     .     'electric potential feedback by current on ',
     1     BCCHAR(IB),boundary_location(ib)
! add smth about #ifdef ASTRA here
          scurpar = 0.0_R8
c Init variables for BCPOT=16                                  { Solovyev 01.05.14
          cur_delta = 0.0_R8              ! square of delta between desired and calculated current
          prev_cur_delta = 0.0_R8         ! previous value of upper value
          po_step = 0.0_R8                ! step of approximation (set after 1st iteration)
          coeff_16 = 1.0_R8
c************************************************************* } Solovyev 01.05.14
          do ibw=ib,nbc
            if (BCCHAR(ibw).eq.'S'.and.BCPOT(ibw).eq.BCPOT(ib)) then
              scurpar = scurpar + potpar(ibw,1)
            endif
          enddo
          write (*,*) 'SCURPAR = ', scurpar
          POS = 0.0_R8                                                  !srv 29.09.15
          VS  = 0.0_R8
          DO IBW=IB,NBC
            if(BCCHAR(IBW).EQ.'S'.and.BCPOT(IBW).EQ.BCPOT(IB)) THEN
              DO IX=1, BC_LIST_SIZE(IBW)
                POS = POS + PO(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     *                     VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                VS = VS + VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
              ENDDO
            endif
          ENDDO
          POS=POS/VS                                                    !srv 29.09.15 }
          po_prev = POS
          write (*,*) 'POS at 0 = ', POS
        endif
        if (.not.once_bcpot) then
          if (BCCHAR(ib).eq.'S') then
            totflux = 0.0_R8
            do ibw=ib,nbc
              if (BCCHAR(ibw).eq.'S'.and.BCPOT(ibw).eq.BCPOT(ib)) then
                do ix=1, bc_list_size(ibw)
                  totflux = totflux +
     &              fch(topix(bc_list_x(ix,ibw),bc_list_y(ix,ibw)),
     &                  topiy(bc_list_x(ix,ibw),bc_list_y(ix,ibw)),1)
                enddo
              endif
            enddo
c The implementation is based on one step of Newton method
            cur_delta = totflux-scurpar                                 ! { Solovyev 03.05.14
            if (ncall.eq.0) then
              write (72, *) '# Newton'
              prev_cur_delta = cur_delta
            else
              if (ncall.eq.1) then                                      ! Solovyev 23.05.14
                po_curr = 0.0_R8                                        !srv 29.09.15
                VS  = 0.0_R8
                DO IBW=IB,NBC
                  if(BCCHAR(IBW).EQ.'S'.and.BCENE(IBW).EQ.BCENE(IB))THEN
                    DO IX=1, BC_LIST_SIZE(IBW)
                      po_curr=po_curr+
     &                   PO(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
     &                 *VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                      VS = VS + VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    ENDDO
                  endif
                ENDDO
                po_curr=po_curr/VS                                      !srv 29.09.15 }
              endif                                                     !Solovyev 23.05.14
              if (mod(ncall,bcpot_16_step).eq.0) then                   !Solovyev 23.05.14
                if (abs(cur_delta-prev_cur_delta).gt.1.0e-20_R8) then   !{ Solovyev 24.05.14
                  coeff_16 = (po_prev - po_curr)/
     &              (cur_delta - prev_cur_delta)
                endif
                po_step = coeff_16*cur_delta                            !} Solovyev 24.05.14
                POS = po_curr + po_step
                prev_cur_delta = cur_delta
                po_prev = po_curr
                po_curr = POS
                write (72, *) ncall, totflux, POS   ! Solovyev 29.04.14 - test output to b2fdebug
              endif                                                     ! Solovyev 23.05.14
            endif                                                       ! } Solovyev 03.05.14
c end of implementation
            write (*, *) 'ncall, totflux, POS', ncall,totflux, POS      !Solovyev 29.04.14 - test output to b2fdebug
            do ibw=ib,nbc
              if (BCCHAR(ibw).eq.'S'.and.BCPOT(ibw).eq.BCPOT(ib)) then
                do ix=1, BC_LIST_SIZE(IBW)
                  sch0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0) =
     &                 PBIG*POS*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)) !Solovyev 27.04.14
                  sch0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1) =
     &                -PBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))     !Solovyev 27.04.14
                enddo
              endif
            enddo
            once_bcpot = .true.
          else
            call xerrab('B2STBC -- BCPOT=16, NOT AVAILABLE FOR N,W,E')
          endif
        endif

* -- BCPOT=17 -- Nothing implemented

      case (17)
        CALL XERRAB ('B2STBC -- BCPOT=17, NOT AVAILABLE FOR S,N,W,E')

* -- BCPOT=18 -- Nothing implemented

      case (18)
        CALL XERRAB ('B2STBC -- BCPOT=18, NOT AVAILABLE FOR S,N,W,E')

* -- BCPOT=19 -- Nothing implemented

      case (19)
        CALL XERRAB ('B2STBC -- BCPOT=19, NOT AVAILABLE FOR S,N,W,E')

* -- BCPOT=20 -- Nothing implemented

      case (20)
        CALL XERRAB ('B2STBC -- BCPOT=20, NOT AVAILABLE FOR S,N,W,E')

* -- BCPOT=21 -- constant potential scaled by potential on the ring bc_type21_ref away

      case (21)
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

      case default
        write(*,*) BCPOT(IB)
        call xerrab ('BCPOT out of bounds')
      end select

      if(istyle_cur_contr_on_S_and_N.eq.0) then                       !srv 01.02.09
!xpb
!xpb  Adding fchy_dia contributions on south and north boundaries
!xpb  following "sv 01.07.99" coding (see B2.5-style implementation)
!xpb  Special attention is paid to the case where there is a discontinuity in B
!xpb
        IF (BCPOT(IB).NE.0.AND.BCPOT(IB).NE.5.AND.BCPOT(IB).NE.8) THEN
          if (ncall.eq.0) write(*,'(a,i2,a,i2)')
     &          'Disabled special treatment of fchy_dia on boundary ',
     &          ib,' because potential is imposed by BCPOT =',BCPOT(IB)
        ELSE IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if(maxval(facdrift) .ne. 0.0_R8) then
              if(core_boundary(ib)) then
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
                  if(ncall.eq.0) write(*,*)
     >             'b2stbc: special treatment for cell ix,iy = ',
     .              bc_list_x(ix,ib),bc_list_y(ix,ib)
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
                  if (ncall.eq.0) write(*,*)
     .             'b2stbc: special treatment for cell ix,iy = ',
     .              bc_list_x(ix,ib),bc_list_y(ix,ib)
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
              if(b2stbc_neoclassical .ne. 0.0_R8 .and.
     &           core_boundary(ib)) then
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
            if(.not.core_boundary(ib) .and.
     &         b2stbc_fchy_dia_coreonly.eq.1) fchy_dia = 0.0_R8
            sch0(bc_list_x(ix,ib),bc_list_y(ix,ib),0) =
     .      sch0(bc_list_x(ix,ib),bc_list_y(ix,ib),0) + fchy_dia
          END DO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if(maxval(facdrift) .ne. 0.0_R8) then
              if(core_boundary(ib)) then
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
     &           core_boundary(ib)) then
                fchy_dia = -praverage(2)
     .           *(OnedBsqc(rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                      rightiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .            -OnedBsqc(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .           *b2stbc_neoclassical
              else
                fchy_dia = 0.0_R8
              endif
            endif
            if(.not.core_boundary(ib).and.
     &         b2stbc_fchy_dia_coreonly.eq.1) fchy_dia = 0.0_R8
            sch0(bc_list_x(ix,ib),bc_list_y(ix,ib),0) =
     .      sch0(bc_list_x(ix,ib),bc_list_y(ix,ib),0) - fchy_dia
          END DO
        ENDIF
!xpb
      elseif(istyle_cur_contr_on_S_and_N.eq.1) then                   !srv 01.02.09 { 27.03.14

        IF (BCPOT(IB).NE.0.AND.BCPOT(IB).NE.5.AND.BCPOT(IB).NE.8) THEN
          if (ncall.eq.0) write(*,'(a,i2,a,i2)')
     &          'Disabled special treatment of fchy_dia on boundary ',
     &          ib,' because potential there is imposed by BCPOT=',
     &          BCPOT(IB)
        ELSE IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if(core_boundary(ib)) then
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
                do is = 0, ns-1                                         !srv 25.06.03
                  if(.not.is_neutral(is))then                           !srv 25.06.03
                    fchy_inertia = fchy_inertia                         !srv 07.04.00 25.06.03
     &               -0.5e0_R8*am(is)*mp*
     &                facdrift(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*phm1  ! *signmf     !srv 07.04.00 21.04.00 21.10.02
     &               *naaverage(is)                                     !lk 20.11.07  !srv 07.04.00
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
      else if (pot_eq.eq.1.and.
     &        (BCPOT(IB).eq.0.or.BCPOT(IB).eq.5.or.BCPOT(IB).eq.8)) then !srv 27.03.14
        if(ncall.eq.0) then
          IF (BCCHAR(IB).EQ.'S') THEN
            IF (core_boundary(ib)) THEN
              write(0,'(a/a)')
     &        'Core boundary potential B.C. does not '//
     &        'ensure drift current continuity',
     &        'Use BCPOT=12 to impose it.'                           !srv 03.03.15
            ELSE
              write(0,'(a/a)')
     &        'South boundary potential B.C. does not '//
     &        'ensure drift current continuity',
     &        'Use BCPOT=13 to impose it.'                           !srv 03.03.15
            ENDIF
          ELSE IF (BCCHAR(IB).EQ.'N') THEN
            write(0,'(a/a)')
     &        'North boundary potential B.C. does not '//
     &        'ensure drift current continuity',
     &        'Use BCPOT=13 to impose it.'                           !srv 03.03.15
          ENDIF
        endif
      endif                                                          !srv 01.02.09 }

      return
      end subroutine b2stbc_phys_pot

      end subroutine b2stbc_phys
