*-----------------------------------------------------------------------
*.specification

      subroutine b2news_m (nx, ny, ns, nscx, iscx,
     &  ismain, ismain0, dtim, rxf, ierr, ramp_call) !srv 23.11.10
      use b2mod_types
      use b2mod_geo
      use b2mod_ma28
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_transport_nspecies
      use b2mod_anomalous_transport
      use b2mod_work
      use b2mod_indirect
      use b2mod_constants
      use b2mod_numerics_namelist
      use b2mod_b2cmpa
      use b2mod_external
      use b2mod_time
#ifdef WITH_FTIMINGS
      use timings
#endif
      use b2mod_ma28_for_7diag
      use b2mod_ma28_for_9diag
      use b2mod_ma28_for_11diag
      use b2mod_ysmp_sdrv
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      integer nscx, iscx(0: nscxmax-1), ismain, ismain0
      real (kind=R8) ::
     *   dtim, rxf
*   ..output arguments                                                    !xpb
      integer ierr
      logical wrong_flow, ramp_call
*   ..common blocks

*-----------------------------------------------------------------------
*
*  1. purpose
*
*     B2NEWS_M performs one inner iteration for the B2 equations.
*     This routine defines the core structure of the solution
*     algorithm.
*
*
*  3. description (see also routine b2cdca)
*
*     For a detailed description of the equations and the method of
*     solution, see the references cited in the main program unit.
*
*     This routine performs one basic iteration for the B2 equations.
*     An outline of the procedure follows.
*     ..test input data
*     ..compute auxiliary quantities
*     ..compute transport coefficients
*     ..iterate on particle balance                   !srv 09.03.21
*     ..iterate on parallel momentum balance          !srv 09.03.21
*     ..iterate on charge balance
*     ..iterate on heat balance
*     Convergence monitoring must be done outside this routine.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer potit, area_fix, nsmin, nsmax, no_current
      integer, save :: style = 3, pot_eq = 1                             !srv 29.05.02 13.07.05
      integer, save :: mdf_fnb = 0, mdf_fhe = 0, mdf_fhi = 0             !srv 13.10.06
      integer, save :: guard_flows = 2
      integer, save :: no_solve = 0                                      !srv 02.07.08
!      integer, save :: re_compute_pb = 1                                !srv 15.02.10
      integer, save :: iout = 0,
     &                 iout_b2npc9 = 0, iout_b2npco = 0,                 !srv 11.09.09
     &                 iout_b2npc7 = 0, iout_b2npc11 = 0,                !srv 02.04.18
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npmo = 0, iout_b2npp7 = 0,
     &                 iout_b2nppo = 0, iout_b2wdat = 0
      integer, save :: ion_vlct_restrict = 0                             !iyv 06.02.14
      integer, save :: ue_style = 2                                      !srv 17.11.16
      integer, save :: get_residuals = 0                                 !srv 22.05.18
      integer, save :: visper_all_ions = 1
      real (kind=R8) ::
     *   potok
      integer ncall, n2, is, ier0, i, ix, iy, ireg
      character charns*3                                                 !srv 30.06.03
      real (kind=R8), save :: BoRiS = 0.0_R8
      real (kind=R8), save :: na_min = 1.0e4_R8, na_new = 1.0e14_R8
      real (kind=R8) ::
     *  xfm0, xfm1, xfm2, xfm3, rxf0, rxf1, rxf2, rxf3,
     1  t0, vis_per                                                      !srv 09.07.99 06.07.00
      real (kind=R8) ::                                                  !xpb
     &  lnlam(-1:nx,-1:ny), zeff(-1:nx,-1:ny),                           !srv 05.07.17
     &  pz(-1:nx,-1:ny), rz(-1:nx,-1:ny), cssb(-1:nx,-1:ny,0:1),         !iyv 06.02.14
     &  zero(-1:nx,-1:ny),                                               !lk 31.01.08 {
     &  fna_53(-1:nx,-1:ny,0:1,0:ns-1),
     &  fne_53(-1:nx,-1:ny,0:1)
      logical solving(4), solvedum(0:NNREGMAX), solvedm2(0:NNREGMAX)
      logical first_pass_5, first_pass_7, first_pass_9, first_pass_11
      logical inverted_ua, new_matrix
      save ncall, xfm0, xfm1, xfm2, xfm3, potit, potok                   !srv 09.07.99
     &     , vis_per, solving, area_fix, no_current                      !srv 06.07.00
     &     , nsmin, nsmax
*   ..procedures
      intrinsic mod
      real (kind=R8) ::
     *  smin, smax, b2sasum
      external xertst, ipgetr, smin, smax,
     &  b2scopy, b2saxpy, b2sasum, ipgeti
      external b2xvsg, b2xvps, b2xpne, b2xpni, b2xpfe, b2xpfi, b2xppb,
     &  b2xpnm, b2tral, b2treq, b2tfnb,
     &  b2npmo, b2npco, b2nppo, b2npht, b2nxdp, b2tfcc, b2npp7, b2xppz  !srv 06.07.00 18.01.02
*   ..initialization
      data ncall/0/, xfm0/1.0e0_R8/, xfm1/1.0e0_R8/, xfm2/1.0e0_R8/,
     1   xfm3/1.0e0_R8/, potit/50/, potok/1e-2_R8/, area_fix/3/,        !srv 09.07.99 09.03.09
     2   vis_per/1.0e0_R8/, no_current/0/                               !srv 06.07.00
      integer, save :: no_neutr_scl = 0, no_b2sral_call = 1
      integer, save :: do_2nd_b2npco_call = 0
      integer, save :: re_eval_prtls_fluxes = 0                         !srv 28.06.13
      integer, save :: nhist = 0, neutr_avg = 0, nrelax = 0, nstg2 = 0
      integer, save :: use_eirene = 0
      integer, save :: ramp_slow = 0
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2news_m')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       nsmin=0
       nsmax=ns
       call ipgetr ('b2news_BoRiS', BoRiS)
       call ipgeti ('b2ux7p_style', style)                             !srv 29.05.02
       call ipgetr ('b2news_xfm0', xfm0)
       call ipgetr ('b2news_xfm1', xfm1)
       call ipgetr ('b2news_xfm2', xfm2)
       call ipgetr ('b2news_xfm3', xfm3)
       call xertst (0.0_R8.le.xfm0.and.xfm0.le.1.0_R8,
     &   'faulty internal parameter xfm0')
       call xertst (0.0_R8.le.xfm1.and.xfm1.le.1.0_R8,
     &   'faulty internal parameter xfm1')
       call xertst (0.0_R8.le.xfm2.and.xfm2.le.1.0_R8,
     &   'faulty internal parameter xfm2')
       call xertst (0.0_R8.le.xfm3.and.xfm3.le.1.0_R8,
     &   'faulty internal parameter xfm3')
       call ipgetr ('b2mndr_na_min', na_min)
       call ipgetr ('b2mndr_na_new', na_new)
       call xertst (na_min.gt.0.0_R8,
     &  'faulty internal parameter na_min')
       call xertst (na_new.gt.0.0_R8,
     &  'faulty internal parameter na_new')
       call xertst (na_new.gt.na_min,
     &  'na_new must be larger than na_min!')
       call ipgeti ('b2news_ramp_slow', ramp_slow)
       call ipgeti ('b2news_potit', potit)
       call xertst (0.lt.potit,
     &   'faulty internal parameter potit')
       call ipgetr ('b2news_potok', potok)
       call xertst (0.0_R8.lt.potok,
     &   'faulty internal parameter potok')
       call ipgetr ('b2tfhe_vis_per', vis_per)                           !srv 06.07.00
       call ipgeti ('b2news_poteq', pot_eq)                              !srv 13.07.05
       call ipgeti ('b2tfhe_no_current', no_current)
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
       call ipgeti ('b2tfnb_mdf', mdf_fnb)                               !srv 13.10.06 {
       call ipgeti ('b2tfhe_mdf', mdf_fhe)
       call ipgeti ('b2tfhi_mdf', mdf_fhi)
       call ipgeti ('b2news_iout', iout)
       call ipgeti ('b2npc9_iout', iout_b2npc9)                          !srv 16.07.10 {
       call ipgeti ('b2npco_iout', iout_b2npco)
       call ipgeti ('b2npc7_iout', iout_b2npc7)                          !srv 02.04.18
       call ipgeti ('b2npc11_iout', iout_b2npc11)                        !srv 02.04.18
       call ipgeti ('b2npmo_iout', iout_b2npmo)
       call ipgeti ('b2npp7_iout', iout_b2npp7)
       call ipgeti ('b2nppo_iout', iout_b2nppo)
       call ipgeti ('b2nph9_iout', iout_b2nph9)
       call ipgeti ('b2npht_iout', iout_b2npht)                          !srv 16.07.10 }
       call ipgeti ('b2wdat_iout', iout_b2wdat)                          !iys 17.04.17
!       call ipgeti ('b2news_re_compute_pb', re_compute_pb)              !srv 15.02.10
       call ipgeti ('b2news_guard_flows', guard_flows)                   !xpb
       call ipgeti ('b2news_area_fix', area_fix)                         !srv 09.03.09
       call ipgeti ('b2npmo_ion_vlct_restrict', ion_vlct_restrict)       !iyv 06.02.14
       call ipgeti ('b2news_no_solve', no_solve)                         !srv 02.07.08 {
       call xertst (-15.le.no_solve.and.no_solve.le.1,
     &   'faulty input no_solve')
*    ..the table below shows which combination of equations will be solved
*      depending on the value of switch 'b2news_no_solve'
*----------------------------------------------------
*        | value |  co  |  mo  |  he  |  hi  |  po  |
*----------------------------------------------------
*        |   0   |  +   |  +   |  +   |  +   |  +   |
*----------------------------------------------------
*        |   1   |  -   |  -   |  -   |  -   |  -   |
*----------------------------------------------------
*        |  -1   |  -   |  +   |  -   |  -   |  -   |
*----------------------------------------------------
*        |  -2   |  +   |  -   |  -   |  -   |  -   |
*----------------------------------------------------
*        |  -3   |  +   |  +   |  -   |  -   |  -   |
*----------------------------------------------------
*        |  -4   |  -   |  -   |  -   |  -   |  +   |
*----------------------------------------------------
*        |  -5   |  -   |  +   |  -   |  -   |  +   |
*----------------------------------------------------
*        |  -6   |  +   |  -   |  -   |  -   |  +   |
*----------------------------------------------------
*        |  -7   |  +   |  +   |  -   |  -   |  +   |
*----------------------------------------------------
*        |  -8   |  -   |  -   |  +   |  +   |  -   |
*----------------------------------------------------
*        |  -9   |  -   |  +   |  +   |  +   |  -   |
*----------------------------------------------------
*        | -10   |  +   |  -   |  +   |  +   |  -   |
*----------------------------------------------------
*        | -11   |  +   |  +   |  +   |  +   |  -   |
*----------------------------------------------------
*        | -12   |  -   |  -   |  +   |  +   |  +   |
*----------------------------------------------------
*        | -13   |  -   |  +   |  +   |  +   |  +   |
*----------------------------------------------------
*        | -14   |  +   |  -   |  +   |  +   |  +   |
*----------------------------------------------------
* An alternative to the 'no_solve' switch for controlling of the
* solving of equations is to use the SOLVEMO, SOLVECO, SOLVEEE, SOLVEEI,
* SOLVEPO arrays from b2.numerics.parameters.
* In time-dependent mode, the arrays SOLVEET and SOLVEMT are not used, as
* the total momentum and energy equations are not solved.
* Please consult the code documentation.

       solving(1:4)=no_solve.eq.0
       do i=1,4
         solving(i)=solving(i).or.
     &    (no_solve.lt.0.and.mod(abs(no_solve),2**i).ge.2**(i-1))
       enddo
       if(.not.solving(1))
     >   write(*,*) 'Not solving parallel momentum equations!'
       if(.not.solving(2))
     >   write(*,*) 'Not solving continuity equations!'
       if(.not.solving(3))
     >   write(*,*) 'Not solving potential equation!'
       if(.not.solving(4))
     >   write(*,*) 'Not solving heat equations!'
       call set_facdrift_fac_ExB_init(nx, ny, ramp_slow)
       call ipgeti ('b2news_no_neutr_scl', no_neutr_scl)
       call ipgeti ('b2news_no_b2sral_call', no_b2sral_call)
       call ipgeti ('b2news_do_2nd_b2npco_call', do_2nd_b2npco_call)
       call ipgeti ('b2news_re_eval_prtls_fluxes', re_eval_prtls_fluxes) !srv 28.06.13

csw 26feb2014 some extra variables
       call ipgeti ('b2mndr_eirene', use_eirene)
       call ipgeti ('eirene_neutr_avg', neutr_avg)
       call ipgeti ('eirene_underrelax', nrelax)
       call ipgeti ('b2stbr_eir_src_nhist', nhist)
       call ipgeti ('b2mndt_nstg2', nstg2)
       if(nstg2 > 1 .and. no_neutr_scl .ne. 0) then
         write(*,*) 'WARNING, b2mndt_nstg2>1 .and. no_neutr_scl.ne.0'
       endif
       call xertst (.not.(
     .    (neutr_avg.gt.1 .and. (nhist.gt.1 .or. nrelax.gt.1)) .or.
     .    (nhist.gt.1 .and. (neutr_avg.gt.1 .or. nrelax.gt.1)) .or.
     .    (nrelax.gt.1 .and. (neutr_avg.gt.1 .or. nhist.gt.1))),
     .   'Cannot have simultaneous multiple Eirene averaging schemes!')
       call ipgeti ('b2news_nsmin', nsmin)
       if(nsmin.gt.ns) then
          nsmin=ns
       elseif(nsmin.lt.0) then
          nsmin=0
       endif
       call ipgeti ('b2news_nsmax', nsmax)
       if(nsmax.gt.ns) then
          nsmax=ns
       elseif(nsmax.lt.0) then
          nsmax=0
       endif
       write(*,*) 'nsmin, nsmax = ', nsmin, nsmax
       call ipgeti ('b2sigp_style', ue_style)
       call ipgeti ('get_residuals', get_residuals)                      !srv 22.05.18
       call ipgeti ('b2tvspr_visper_all_ions', visper_all_ions)
      endif                                                              !srv 02.07.08
*   ..set n2, zero
      n2 = (nx+2)*(ny+2)
      zero(-1:nx,-1:ny) = 0.0_R8
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
*   ..test facdrift
      call xertst (0.0_R8.le.minval(facdrift) .and.
     &     1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
      call xertst (0.0_R8.le.minval(fac_ExB) .and.
     &     1.0_R8.ge.maxval(fac_ExB), 'faulty argument fac_ExB')
      call xertst (0.0_R8.le.minval(fac_vis) .and.
     &     1.0_R8.ge.maxval(fac_vis), 'faulty argument fac_vis')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of bb
       call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*    ..test vol, hx, hy, qz, qc
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
       call xertst (-1.0_R8.lt.smin(n2,qz(-1,-1,0),1).and.
     &              smax(n2,qz(-1,-1,0),1).lt.1.0_R8,
     .              'faulty argument range qz(,,0)')
       call xertst (0.0_R8.lt.smin(n2,qz(-1,-1,1),1).and.
     &              smax(n2,qz(-1,-1,1),1).le.1.0_R8,
     .              'faulty argument range qz(,,1)')
       call xertst (0.0_R8.lt.smin(n2,qc,1).and.
     &              smax(n2,qc,1).le.1.0_R8,
     .              'faulty argument range qc')
*    ..test present state
       call b2xvps (nx, ny, ns,
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
       call b2xvsg (n2*ns, kinrgy, 1, 'kinrgy', '.ge.')
      endif
*   ..set rxf0/1/2/3
*     Note: The potential equation needs an accurate solution on each
*     iteration; hence the different setting of rxf2.
      rxf0 = rxf*xfm0
      rxf1 = rxf*xfm1
      rxf2 = 1.0e0_R8 ! (rxf*xfm2)**0.2e0_R8                             !srv 15.01.03
      if (rxf .eq. 0.0e0_R8) rxf2 = 0.0e0_R8
      rxf3 = rxf*xfm3
cxpb  Force a restart of the matrix builds on the first instance of every call
      first_pass_5 = .true.
      first_pass_7 = .true.
      first_pass_9 = .true.
      first_pass_11 = .true.

csw 01oct2012 initialize potential if necessary
      if (pot_eq.eq.0) then                                              !srv 13.07.05 {
        po = 0.0_R8
      else if (pot_eq.eq.2) then
        po(:,:)=3.1_R8*te(:,:)/qe
      endif

      if (ramp_call .or. ncall.eq.0) call set_facdrift_fac_ExB_inc
*   ..test facdrift, fac_ExB, fac_vis
      call xertst (0.0_R8.le.minval(facdrift) .and.
     &     1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
      call xertst (0.0_R8.le.minval(fac_ExB) .and.
     &     1.0_R8.ge.maxval(fac_ExB), 'faulty argument fac_ExB')
      call xertst (0.0_R8.le.minval(fac_vis) .and.
     &     1.0_R8.ge.maxval(fac_vis), 'faulty argument fac_vis')
* ..main computation
*   ..compute ni, ne, ne2
      call b2xpni (nx, ny, ns, na, ni)
      call b2xpne (nx, ny, ns, rza, na, ne_ext, ne)
      call b2xpne (nx, ny, ns, rz2, na, ne2_ext, ne2)
!*   ..compute ceqp                                                      !srv 03.07.03 {
!      call b2xpnm (nx, ny, ns, rz2, am, na, wrk0)
!      call b2treq (nx, ny, vol, te, ti, ne, ni, wrk0, lnlam, ceqp)      !srv 03.07.03 }
*   ..compute transport coefficients
      call b2tral (nx, ny, ns, nscx, nscxmax, iscx, ismain,              !srv 09.07.04 23.11.10
     1 bb, conn, vol, gs, hx, hy, hz, qz, qc, pbs, crx, cry, bzb,        !srv 09.01.01 13.10.06 01.07.09
     2 lnlam, fch, na, ua, te, ti, po, ne, ni, ne2, chvemx, chvimx,      !srv 03.06.03
     3 cdna, cdpa, cddi, cdde, cvla, cvsa, chce, chve, chci,             !srv 04.10.99 23.12.03 13.10.06 !som 03.11.21
     5 chvi, csig, csigin, calf, cthe, cthi,
     4 cdnahz, cdpahz, cvlahz, cvmahz, cvsahz, cvsa_cl, cvsahz_cl,
     * cvsahz_drho, cvsa_drho,                                           !som 18.08.21
     * chci_al_ast, chci_a, cvsahz_hAdp_albe, cvsahz_hBdp_al,            !som 08.10.21
     * cvsa_hAdp_albe, cvsa_hBdp_al,                                     !som 20.10.21
     * cvsahz_uAdp_albe, cvsahz_uBdp_al,                                 !som 22.12.23
     * cvsa_uAdp_albe, cvsa_uBdp_al,                                     !som 26.04.24
     5 fllim0fhi, fllimvisc, csig_cl, calf_cl, csig_stoch, chce_stoch)   !srv 09.01.01 15.12.05 24.08.21
*                                                                        !srv 09.03.21 {
*   ..recompute currents
      call b2tanml (nx, ny, ns, ismain,                                  !srv 20.09.06 }
     &              csig_an, po, ne, na, fchanml_a, fchanml)             !srv 08.09.21 !som 02.11.21
      call b2tvspa (nx, ny, ns, ismain, hy1, bb,                         !srv 02.11.21
     &              na, ua, ti, cvsa_cl, cvsa_drho, fac_vis,
     &              cvsa_uAdp_albe, cvsa_uBdp_al,                        !som 26.04.24
     &              fchvispar_a, fchvispar)                              !som 18.08.21
      call b2tiner (nx, ny, ns, ismain, vol, gs, bb, na, ua,             !srv 08.09.21
     &              facdrift, fchinert_a, fchinert)
      call b2tvspr (nx, ny, ns, ismain, fac_vis, vol, hy1,               !srv 10.10.17 02.11.21
     &              po, ti, na, avisper_a, bvisper_a, cvisper_a,
     &              fchvisper_a, fchvisper)
      call b2tvsq (nx, ny, ns, ismain, fac_vis,                          !srv 10.10.17 08.09.21 } 02.11.21
     &             ne2, lnlam, te, ti, ni(-1,-1,1), chci, na,            !som 21.10.21
     &             chci_al_ast, chci_a, cvsa_hAdp_albe, cvsa_hBdp_al,    !som 20.10.21
     &             fchvisq_a, fchvisq)
      call b2tinnt (nx, ny, ns, nscx, iscx,
     &              facdrift, fac_ExB,
     &              vol, gs, hx, hy1, hz, qz, qc, bb,
     &              te, po, ti, na, fna, csigin, fchin_a, fchin)         !IYS 17.10.2023

*     ..iterate on particle balance

      if(ion_vlct_restrict .ne. 0 ) then                                 !iyv 06.02.14
        call b2xprz (nx, ny, ns, mp, am, na, rz)
        call b2xppz (nx, ny, ns, ne, na, te, ti, pz)
        wrk0(:,:)=sqrt(pz(:,:)/rz(:,:))*abs(bb(:,:,0)/bb(:,:,3))
        call intfaceh(nx,ny,leftix,leftiy,vol,wrk0,cssb(-1,-1,0))
        call intfacev(nx,ny,bottomix,bottomiy,vol,wrk0,cssb(-1,-1,1))
      else
        cssb = 0.0_R8
      endif                                                              !iyv 06.02.14

#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2NEWS_UNDERSCORE_LOOP2) )
#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
!!FR Possible loop (2) for use with OpenMP
!$OMP PARALLEL DO DEFAULT(none)
!$OMP& PRIVATE(ier0, solvedum, solvedm2, ix, iy, ireg)
!$OMP& PRIVATE(new_matrix)
!$OMP& FIRSTPRIVATE(pb, flob, conb)
!$OMP& LASTPRIVATE(pb, flob, conb)
!$OMP& FIRSTPRIVATE(first_pass_5, first_pass_7)
!$OMP& FIRSTPRIVATE(first_pass_9, first_pass_11)
!$OMP& SHARED(last_solve_5, last_solve_7)
!$OMP& SHARED(last_solve_9, last_solve_11)
!$OMP& SHARED(nsmin, nsmax, na_min, ismain, mdf_fnb, visper_all_ions)
!$OMP& SHARED(ns, nx, ny, region, leftix, rightix, bottomiy)
!$OMP& SHARED(topiy, na, rza, te, ti, am, ua, ne)
!$OMP& SHARED(zn, zamax)
!$OMP& SHARED(is_neutral, uadia, ncall, resmo, corua)
!$OMP& SHARED(guard_flows, facdrift, fac_ExB, fac_vis, vol, pbs, pbshz)
!$OMP& SHARED(cssb, vaecrb, vadia)
!$OMP& SHARED(wadia, bb, hx, hy, hz, qz, gs, qc, po)
!$OMP& SHARED(cvla, cdna, cdpa, cddi, cdde, cvlahz, cdnahz, cdpahz)
!$OMP& SHARED(fchvispar, fchvisper, fchinert, fchanml, fchvisq, fchin)
!$OMP& SHARED(fchvispar_a,fchvisper_a, fchinert_a, fchanml_a, fchvisq_a)
!$OMP& SHARED(fchin_a, fna, fna_nodrift, fna_fcor, fna_he, fnaPSch)
!$OMP& SHARED(fna_mdf, fna_32, fna_52, fna_32nd, fna_52nd, fna_53)
!$OMP& SHARED(kinrgy, fllim0fna, solving, nnreg, solveco, solvemo)
!$OMP& SHARED(rxf0, dtim, sna, pcca, resco, corpa, b2stbr_sna)
!$OMP& SHARED(b2stbm_sna,ext_sna)
!$OMP& COPYIN(iwrk_5p, aaf_5p, wrk0_5p, wrk1_5p, wrk2_5p)
!$OMP& COPYIN(iwrk_7p, aaf_7p, wrk0_7p, wrk1_7p, wrk2_7p)
!$OMP& COPYIN(iwrk_9p, aaf_9p, wrk0_9p, wrk1_9p, wrk2_9p)
!$OMP& COPYIN(iwrk_11p, aaf_11p, wrk0_11p, wrk1_11p, wrk2_11p)
!$OMP& COPYIN(nc0_5p, nc0_7p, nc0_9p, nc0_11p)
!$OMP& COPYIN(nltrsol_5p, solvdim1_5p, solvdim_5p, n_nonzeros_5p)
!$OMP& COPYIN(nltrsol_7p, solvdim1_7p, solvdim_7p, n_nonzeros_7p)
!$OMP& COPYIN(nltrsol_9p, solvdim1_9p, solvdim_9p, n_nonzeros_9p)
!$OMP& COPYIN(nltrsol_11p, solvdim1_11p, solvdim_11p, n_nonzeros_11p)
!$OMP& COPYIN(r1_5p, s1_5p, a1h_5p, w_5p, res_5p, a1h_save_5p)
!$OMP& COPYIN(irn_5p, icn_5p, irnh_5p, icnh_5p)
!$OMP& COPYIN(ish_5p, jsh_5p, ijsh_5p, iw28_5p)
!$OMP& COPYIN(icnhkeep_5p, ikeep_5p, irnh_save_5p, icnh_save_5p)
!$OMP& COPYIN(r1_7p, s1_7p, a1h_7p, w_7p, res_7p, a1h_save_7p)
!$OMP& COPYIN(irn_7p, icn_7p, irnh_7p, icnh_7p)
!$OMP& COPYIN(ish_7p, jsh_7p, ijsh_7p, iw28_7p)
!$OMP& COPYIN(icnhkeep_7p, ikeep_7p, irnh_save_7p, icnh_save_7p)
!$OMP& COPYIN(r1_9p, s1_9p, a1h_9p, w_9p, res_9p, a1h_save_9p)
!$OMP& COPYIN(irn_9p, icn_9p, irnh_9p, icnh_9p)
!$OMP& COPYIN(ish_9p, jsh_9p, ijsh_9p, iw28_9p)
!$OMP& COPYIN(icnhkeep_9p, ikeep_9p, irnh_save_9p, icnh_save_9p)
!$OMP& COPYIN(r1_11p, s1_11p, a1h_11p, w_11p, res_11p, a1h_save_11p)
!$OMP& COPYIN(irn_11p, icn_11p, irnh_11p, icnh_11p)
!$OMP& COPYIN(ish_11p, jsh_11p, ijsh_11p, iw28_11p)
!$OMP& COPYIN(icnhkeep_11p, ikeep_11p, irnh_save_11p, icnh_save_11p)
!$OMP& COPYIN(r1_ysmp, s1_ysmp, a, rsp, lfirst_sdrv)
!$OMP& COPYIN(ia, ja, p, ip, isp, indexes, iasave, jasave)
!$OMP& COPYIN(indexes_9p, iasave_9p, jasave_9p)
!$OMP& COPYIN(indexes_11p, iasave_11p, jasave_11p)
!$OMP& COPYIN(lfirst_ilu, lfirst_ilu7, lfirst_ilu9, lfirst_ilu11)
#endif
!      do is = 0, ns-1
      do is = nsmin, nsmax-1
*    ..DPC
        do iy=-1,ny
          do ix=-1,nx
            na(ix,iy,is)=max(na(ix,iy,is),na_min) ! abs(na(ix,iy,is))    !srv 24.04.09
          enddo
        enddo
*    ..compute partial pressure
        call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)
*    ..compute particle flux
        call b2tfnb (nx, ny, is, ismain,
     &   guard_flows, facdrift, fac_ExB, am(is)*mp, vol, cssb,           !srv 07.06.00 09.01.01 08.09.01 19.05.03 !iyv 06.02.14
     1   rza(-1,-1,is), pbs, pbshz, na(-1,-1,is), pb, ua(-1,-1,is),      !srv 09.01.01 20.09.06
     1   uadia(-1,-1,0,is), vaecrb(-1,-1,0,is), vadia(-1,-1,0,is),       !srv 17.05.99
     1   wadia(-1,-1,0,is), ne, te, ti, bb, hx, hy, hz, gs, qc,          !srv 19.05.03 08.10.04 04.02.05
     &   po, cvla(-1,-1,0,is), cdna(-1,-1,0,is), cdpa(-1,-1,0,is),
     &   cddi(-1,-1,0,is),                                               !srv 13.10.06
     &   cvlahz(-1,-1,0,is),                                             !srv 09.01.01 15.12.05
     &   cdnahz(-1,-1,0,is), cdpahz(-1,-1,0,is),
     &   fchvispar_a(-1,-1,0,is), fchvisper_a(-1,-1,0,is),               !srv 08.09.21
     &   fchinert_a(-1,-1,0,is), fchanml_a(-1,-1,0,is),
     &   fchvisq_a(-1,-1,0,is), fchin_a(-1,-1,0,is),                     !srv 20.09.06, IYS 17.10.2023
     &   fna(-1,-1,0,is), fna_nodrift(-1,-1,0,is),
     &   fna_fcor(-1,-1,0,is),
     &   fna_he(-1,-1,0,is), fnaPSch(-1,-1,0,is),
     &   fna_mdf(-1,-1,0,is),                                            !srv 13.10.06
     &   flob, conb, fna_32(-1,-1,0,is), fna_52(-1,-1,0,is),
     &   fna_32nd(-1,-1,0,is), fna_52nd(-1,-1,0,is),
     &   fna_53(-1,-1,0,is),
     &   kinrgy(-1,-1,is), fllim0fna(-1,-1,0,is),
     &   .true.)                                                         !djm Jan2017
!        call b2tfco (nx, ny, is, na(-1,-1,is), te, ti, pb,
!     &      cdpa(-1,-1,0,is), flob, conb)
*    ..perform iteration

csw 16oct2012 init ier0
        ier0=0
!       if (solving(2).and.all(solveco(is,0:nnreg(0)))) then             !srv 26.09.12 22.05.18
          solvedum(0:nnreg(0))=solveco(is,0:nnreg(0))
          solvedm2(0:nnreg(0))=solvemo(is,0:nnreg(0))
          new_matrix = .false.
          if(mdf_fnb.eq.0) then                                          !srv 13.10.06
            if((visper_all_ions.ne.0 .and. .not.is_neutral(is)) .or.
     &        (visper_all_ions.eq.0 .and. is.eq.ismain)) then
              new_matrix = first_pass_7
              do ireg = 0, nnreg(0)
                new_matrix = new_matrix.or.
     &           solveco(is,ireg).neqv.last_solve_7(is,ireg)
              end do
              if(new_matrix.or.in_parallel()) call restart_ma28_for_7
              call b2npc7 (nx, ny, fac_vis, guard_flows, rza(-1,-1,is),
     &         nnreg(0), solving(2), solvedum, solving(1),               !srv 22.05.18
     &         solvedm2, ncall, is, rxf0, dtim, vol, pbs, flob, conb,    !srv 22.01.02 01.07.05 24.05.21
     &         sna(-1,-1,0,is), pcca(-1,-1,0,is),
     &         ti, pb, na(-1,-1,is), ua(-1,-1,is), fna(-1,-1,0,is),
     &         resco(-1,-1,is), corpa(-1,-1,is), ier0)
              call xertst (ier0.eq.0, 'error return from b2npc7')
              last_solve_7(is,0:nnreg(0))=solveco(is,0:nnreg(0))
              first_pass_7 = .false.
            else                                                         !srv 02.04.18
              new_matrix = first_pass_5
              do ireg = 0, nnreg(0)
                new_matrix = new_matrix.or.
     &           solveco(is,ireg).neqv.last_solve_5(is,ireg)
              end do
              if(new_matrix.or.in_parallel()) call restart_ma28
              call b2npco (nx, ny, guard_flows, nnreg(0), solving(2),    !srv 22.05.18
     &         solvedum, solving(1), solvedm2, ncall, is, rxf0, dtim,    !srv 24.05.21
     &         vol, pbs, flob, conb,                                     !srv 22.01.02 01.07.05
     &         sna(-1,-1,0,is), pcca(-1,-1,0,is), pb, na(-1,-1,is),
     &         ua(-1,-1,is), fna(-1,-1,0,is),
     &         resco(-1,-1,is), corpa(-1,-1,is), ier0)
              call xertst (ier0.eq.0, 'error return from b2npco')
              last_solve_5(is,0:nnreg(0))=solveco(is,0:nnreg(0))
              first_pass_5 = .false.
            endif                                                        !srv 02.04.18
          else                                                           !srv 13.10.06 {
            if(.not.is_neutral(is)) then
              if(visper_all_ions.ne.0 .or. is.eq.ismain) then
                new_matrix = first_pass_11
                do ireg = 0, nnreg(0)
                  new_matrix = new_matrix.or.
     &             solveco(is,ireg).neqv.last_solve_11(is,ireg)
                end do
                if(new_matrix.or.in_parallel()) call restart_ma28_for_11
                call b2npc11 (nx, ny, fac_vis, guard_flows,
     &           rza(-1,-1,is), nnreg(0), solving(2), solvedum,          !srv 22.05.18
     &           solving(1), solvedm2, ncall, is, rxf0, dtim,            !srv 24.05.21
     &           vol, pbs, flob, conb, sna(-1,-1,0,is),
     &           pcca(-1,-1,0,is), po, ti, pb, na(-1,-1,is),
     &           ua(-1,-1,is), fna_mdf(-1,-1,0,is),
     &           resco(-1,-1,is), corpa(-1,-1,is), ier0)
                call xertst (ier0.eq.0, 'error return from b2npc11')
                last_solve_11(is,0:nnreg(0))=solveco(is,0:nnreg(0))
                first_pass_11 = .false.
              else                                                       !srv 02.04.18
                new_matrix = first_pass_9
                do ireg = 0, nnreg(0)
                  new_matrix = new_matrix.or.
     &             solveco(is,ireg).neqv.last_solve_9(is,ireg)
                end do
                if(new_matrix.or.in_parallel()) call restart_ma28_for_9
                call b2npc9 (nx, ny, guard_flows, nnreg(0), solving(2),  !srv 22.05.18
     &           solvedum, solving(1), solvedm2, ncall, is, rxf0, dtim,  !srv 24.05.21
     &           vol, pbs, flob, conb,
     &           sna(-1,-1,0,is), pcca(-1,-1,0,is), pb, na(-1,-1,is),
     &           ua(-1,-1,is), fna_mdf(-1,-1,0,is),
     &           resco(-1,-1,is), corpa(-1,-1,is), ier0)
                call xertst (ier0.eq.0, 'error return from b2npc9')
                last_solve_9(is,0:nnreg(0))=solveco(is,0:nnreg(0))
                first_pass_9 = .false.
              endif
            else
              new_matrix = first_pass_5
              do ireg = 0, nnreg(0)
                new_matrix = new_matrix.or.
     &           solveco(is,ireg).neqv.last_solve_5(is,ireg)
              end do
              if(new_matrix.or.in_parallel()) call restart_ma28
              call b2npco (nx, ny, guard_flows, nnreg(0), solving(2),    !srv 22.05.18
     &         solvedum, solving(1), solvedm2, ncall, is, rxf0, dtim,    !srv 24.05.21
     &         vol, pbs, flob, conb,
     &         sna(-1,-1,0,is), pcca(-1,-1,0,is), pb, na(-1,-1,is),
     &         ua(-1,-1,is), fna_mdf(-1,-1,0,is),
     &         resco(-1,-1,is), corpa(-1,-1,is), ier0)
              call xertst (ier0.eq.0, 'error return from b2npco')
              last_solve_5(is,0:nnreg(0))=solveco(is,0:nnreg(0))
              first_pass_5 = .false.
            endif
          endif
!       endif                                                            !srv 13.10.06 }  22.05.18

        do iy=-1,ny
          do ix=-1,nx
            na(ix,iy,is)= max(na(ix,iy,is),na_min) !abs(na(ix,iy,is))    !srv 24.04.09
          enddo
        enddo
!        if(re_compute_pb.eq.1) then                                     !srv 15.02.10
!*    ..re-compute partial pressure                                      !srv 01.10.99
!        call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)   !srv 01.10.99
!        endif                                                           !srv 15.02.10
*    ..re-compute particle flux
        call b2tfnb (nx, ny, is, ismain,
     &   guard_flows, facdrift, fac_ExB, am(is)*mp, vol, cssb,           !srv 07.06.00 09.01.01 08.09.01 19.05.03 !iyv 06.02.14
     1   rza(-1,-1,is), pbs, pbshz, na(-1,-1,is), pb, ua(-1,-1,is),      !srv 09.01.01 20.09.06
     1   uadia(-1,-1,0,is), vaecrb(-1,-1,0,is), vadia(-1,-1,0,is),       !srv 17.05.99
     1   wadia(-1,-1,0,is), ne, te, ti, bb, hx, hy, hz, gs, qc,          !srv 19.05.03 08.10.04 04.02.05
     &   po, cvla(-1,-1,0,is), cdna(-1,-1,0,is), cdpa(-1,-1,0,is),
     &   cddi(-1,-1,0,is),                                               !srv 13.10.06
     &   cvlahz(-1,-1,0,is),                                             !srv 09.01.01 15.12.05
     &   cdnahz(-1,-1,0,is), cdpahz(-1,-1,0,is),
     &   fchvispar_a(-1,-1,0,is), fchvisper_a(-1,-1,0,is),               !srv 08.09.21
     &   fchinert_a(-1,-1,0,is), fchanml_a(-1,-1,0,is),
     &   fchvisq_a(-1,-1,0,is), fchin_a(-1,-1,0,is),                     !srv 20.09.06, IYS 17.10.2023
     1   fna(-1,-1,0,is), fna_nodrift(-1,-1,0,is),
     &   fna_fcor(-1,-1,0,is),
     &   fna_he(-1,-1,0,is), fnaPSch(-1,-1,0,is),
     &   fna_mdf(-1,-1,0,is),                                            !srv 13.10.06
     &   flob, conb, fna_32(-1,-1,0,is), fna_52(-1,-1,0,is),
     &   fna_32nd(-1,-1,0,is), fna_52nd(-1,-1,0,is),
     &   fna_53(-1,-1,0,is),
     &   kinrgy(-1,-1,is), fllim0fna(-1,-1,0,is),
     &   .false.)                                                        !djm Jan2017
      enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2NEWS_UNDERSCORE_LOOP2) )
!$OMP END PARALLEL DO
!! FR End of loop (2)
#ifdef WITH_FTIMINGS
      call timer%stop("parallel")
#endif
#endif
*      ..compute electron drift velocities
      call b2tfed (nx, ny, facdrift, fac_ExB, te, po, vedia, veecrb)

*   ..re-compute ne, ni, lnlam
      call b2xpne (nx, ny, ns, rza, na, ne_ext, ne)
      call nancheck(nx, ny, ne, 'b2news(ne): after call to b2xpne')
      call b2xpni (nx, ny, ns, na, ni)
      call b2xpne (nx, ny, ns, rz2, na, ne2_ext, ne2)                    !srv 11.07.99
      call b2xppz (nx, ny, ns, ne, na, te, ti, pz)                       !srv 18.01.02

*   ..compute ceqp                                                       !srv 03.07.03 {
      call b2xpnm (nx, ny, ns, rz2, am, na, wrk0)
      call b2treq (nx, ny, vol, te, ti, ne, ni, wrk0, lnlam, ceqp)       !srv 03.07.03 }
      if(iout .ne. 0) then                                               !srv 12.07.05 {
        do is = 0,ns-1
          write(charns,'(i3.3)') is
          call my_out(70,nx,ny,na(-1,-1,is),'b2news_m_na'//charns)
        enddo
        call my_out(70,nx,ny,ni(-1,-1,0),'b2news_m_ni0')
        call my_out(70,nx,ny,ni(-1,-1,1),'b2news_m_ni1')
        call my_out(70,nx,ny,ne,'b2news_m_ne')
        call my_out(70,nx,ny,pz,'b2news_m_pz')
        call my_out(70,nx,ny,ceqp,'b2news_m_ceqp')
      endif                                                              !srv 12.07.05 }

      call b2tanml (nx, ny, ns, ismain,                                  !srv 20.09.06 }!som 02.11.21
     &              csig_an, po, ne, na, fchanml_a, fchanml)             !srv 08.09.21
      call b2tvspa (nx, ny, ns, ismain, hy1, bb,                         !srv 19.02.21 !srv 02.11.21
     &              na, ua, ti, cvsa_cl, cvsa_drho, fac_vis,
     &              cvsa_uAdp_albe, cvsa_uBdp_al,                        !som 26.04.24
     &              fchvispar_a, fchvispar)                              !som 18.08.21
      call b2tiner (nx, ny, ns, ismain, vol, gs, bb,                     !srv 08.09.21
     &              na, ua, facdrift, fchinert_a, fchinert)
      call b2tvspr (nx, ny, ns, ismain, fac_vis,                         !srv 10.10.17 02.11.21
     &              vol, hy1, po, ti, na,
     &              avisper_a, bvisper_a, cvisper_a,
     &              fchvisper_a, fchvisper)
      call b2tvsq (nx, ny, ns, ismain, fac_vis,                          !srv 10.10.17 08.09.21 } 02.11.21
     &             ne2, lnlam, te, ti, ni(-1,-1,1), chci, na,            !som 21.10.21
     &             chci_al_ast, chci_a, cvsa_hAdp_albe, cvsa_hBdp_al,    !som 20.10.21
     &             fchvisq_a, fchvisq)
      call b2tinnt (nx, ny, ns, nscx, iscx,
     &              facdrift, fac_ExB,
     &              vol, gs, hx, hy1, hz, qz, qc, bb,
     &              te, po, ti, na, fna, csigin, fchin_a, fchin)         !IYS 17.10.2023

*   ..iterate on momentum balance
!      if (ncall.eq.0) then                                              !srv 01.10.99
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2NEWS_UNDERSCORE_LOOP1) )
#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
!! FR possible loop (1) for use with OpenMP
!$OMP PARALLEL DO DEFAULT(none) FIRSTPRIVATE(flob, conb, pb)
!$OMP& SHARED(nx, ny, rza, na, te, ti, ns, nsmin, nsmax, ismain)
!$OMP& SHARED(guard_flows, facdrift, fac_ExB, am, cssb)
!$OMP& SHARED(vol, pbs, pbshz, ua, uadia, vaecrb, vadia, wadia)
!$OMP& SHARED(ne, bb, hx, hy, hz, qz, gs, qc, po, cvla, cdna, cdpa)
!$OMP& SHARED(cddi, cdde, cvlahz, cdnahz, cdpahz, fchvispar)
!$OMP& SHARED(fchvisper, fchinert, fchanml, fchvisq, fchin)
!$OMP& SHARED(fchvispar_a,fchvisper_a, fchinert_a, fchanml_a, fchvisq_a)
!$OMP& SHARED(fchin_a, fna, fna_nodrift, fna_fcor, fna_he, fnaPSch)
!$OMP& SHARED(fna_mdf, fna_32, fna_52, fna_32nd, fna_52nd, fna_53)
!$OMP& SHARED(kinrgy, fllim0fna)
#endif
!        do is=0,ns-1                                                    !srv 01.10.99
        do is = nsmin, nsmax-1
*    ..compute partial pressure
          call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)  !srv 01.10.99
*    ..compute particle flux
          call b2tfnb (nx, ny, is, ismain,
     &   guard_flows, facdrift, fac_ExB, am(is)*mp, vol, cssb,           !srv 07.06.00 09.01.01 08.09.01 19.05.03 !iyv 06.02.14
     1   rza(-1,-1,is), pbs, pbshz, na(-1,-1,is), pb, ua(-1,-1,is),      !srv 09.01.01 20.09.06
     1   uadia(-1,-1,0,is), vaecrb(-1,-1,0,is), vadia(-1,-1,0,is),
     1   wadia(-1,-1,0,is), ne, te, ti, bb, hx, hy, hz, gs, qc,          !srv 19.05.03 08.10.04 04.02.05
     &   po, cvla(-1,-1,0,is), cdna(-1,-1,0,is), cdpa(-1,-1,0,is),
     &   cddi(-1,-1,0,is),                                               !srv 13.10.06
     &   cvlahz(-1,-1,0,is),                                             !srv 09.01.01 15.12.05
     &   cdnahz(-1,-1,0,is), cdpahz(-1,-1,0,is),
     &   fchvispar_a(-1,-1,0,is), fchvisper_a(-1,-1,0,is),               !srv 08.09.21
     &   fchinert_a(-1,-1,0,is), fchanml_a(-1,-1,0,is),
     &   fchvisq_a(-1,-1,0,is), fchin_a(-1,-1,0,is),                     !srv 20.09.06, IYS 17.10.2023
     &   fna(-1,-1,0,is), fna_nodrift(-1,-1,0,is),
     &   fna_fcor(-1,-1,0,is),
     &   fna_he(-1,-1,0,is), fnaPSch(-1,-1,0,is),
     &   fna_mdf(-1,-1,0,is),                                            !srv 13.10.06
     &   flob, conb, fna_32(-1,-1,0,is), fna_52(-1,-1,0,is),
     &   fna_32nd(-1,-1,0,is), fna_52nd(-1,-1,0,is),
     &   fna_53(-1,-1,0,is),
     &   kinrgy(-1,-1,is), fllim0fna(-1,-1,0,is),
     &   .false.)                                                         !djm Jan2017
        enddo                                                             !srv 01.10.99
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2NEWS_UNDERSCORE_LOOP1) )
!$OMP END PARALLEL DO
!! FR end of loop (1)
#ifdef WITH_FTIMINGS
        call timer%stop("parallel")
#endif
#endif
*      ..compute electron drift velocities
        call b2tfed (nx, ny, facdrift, fac_ExB, te, po, vedia, veecrb)
!      end if                                                            !srv 01.10.99

csw 10dec2012: added extra call if wanted/not-wanted (SOLPS4 has an
c   extra call at this place within internal iteration)
      if(no_b2sral_call .eq. 0) then
        write(*,*) ' re-compute source '
*     ..re-compute source linearisation                                  !xpb
        call b2sral (nx, ny, ns,
     &    nscx, nscxmax, nsmin, nsmax, iscx, ismain, ismain0,
     &    dtim, BoRiS, facdrift, fac_ExB, fac_vis,                       !xpb
     &    vol, hx, hy, hy1, hz, qz, qc, gs, pbs, bb, lnlam,              !srv 09.01.01
     &    na, ua, uadia, vedia, vadia, wadia, veecrb, vaecrb, ve, wedia,
     &    te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,                !srv 28.10.02
     &    fna, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &    fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,      !srv 13.10.06 25.12.06 23.10.09
     &    fhe, fhi, fhm, fht, fnaPSch, fhePSch, fhiPSch, fch,
     &    fchanml_a, fchinert_a, fchvispar_a, fchvisper_a, fchvisq_a,    !srv 08.09.21
     &    fchin_a, fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq,  !srv 21.10.02, IYS 17.10.2023
     &    fchinert, fchanml, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &    cdna, cdpa, cvsa_cl, cvsa_drho,                                !som 18.08.21
     &    chci_al_ast, chci_a, cvsa_hAdp_albe, cvsa_hBdp_al,             !som 20.10.21
     &    cvsa_uAdp_albe, cvsa_uBdp_al,                                  !som 26.04.24
     &    cvla, chce, chve, chci, chvi, calf, csigin,                    !xpb
     &    rlsa, rlra, rlqa, rlcx, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &    rza, rz2, rpt, rpi, sna, smo, smq, she, shi, sch, sne,
     &    wrong_flow, .false.)
        if(wrong_flow) then
          inverted_ua = .true.
          write(*,*) 'b2news_m: wrong_flow returned from b2sral'
          do iy = -1,ny
            if(topiy(-1,iy).eq.ny+1) cycle
            if(bottomiy(-1,iy).eq.-2) cycle
            if(region(-1,iy,0).eq.0) cycle
            if(.not.on_closed_surface(-1,iy)) then
              write(*,*) iy,'west =',
     >           fna(rightix(-1,iy),rightiy(-1,iy),0,ismain)
              inverted_ua = inverted_ua .and.
     .         ua(-1,iy,ismain)*
     .         pbs(rightix(-1,iy),rightiy(-1,iy),0).gt.0.0_R8
            end if
          end do
          do iy = -1,ny
            if(topiy(nx,iy).eq.ny+1) cycle
            if(bottomiy(nx,iy).eq.-2) cycle
            if(region(nx,iy,0).eq.0) cycle
            if(.not.on_closed_surface(nx,iy)) then
              write(*,*) iy,'east =',fna(nx,iy,0,ismain)
              inverted_ua = inverted_ua .and.
     .         ua(nx,iy,ismain)*pbs(nx,iy,0).lt.0.0_R8
            end if
          end do
          if (inverted_ua) then
            write(*,*) 'All target velocities have the wrong sign'
            write(*,*)
     .       'Please consider using the b2mndr_inverse_ua switch'
            call xerrab ('Inverted velocity field!')
          end if
        end if

*     ..smooth out sources
        call b2srsm (nx, ny, ns, dtim,
     &    vol, na, ua, te, ti, ne, ni,
     &    sna, smo, she, shi, sne, .false.)
*     ..include contributions due to timestep                            !xpb
        call b2srdt (nx, ny, ns, guard_flows, dtim, BoRiS, hz,           !xpb !srv 09.01.01 01.07.08
     &    vol, na0, ua0, te0, ti0, ne0, ni0, kinrgy0,
     &    na, ua, te, ti, ne, ni, kinrgy,                                !srv 16.07.10
     &    sna, smo, she, shi, sch, sne, .false.)                         !xpb

*     ..stabilise the source coefficients                                !xpb
        call b2srst (nx, ny, ns, guard_flows, na, ua, te, ti, po, ne,    !xpb !srv 01.07.08
     &    ni, sna, smo, she, shi, sch, sne)                              !xpb
      endif

      if (ue_style.eq.2 .or. iout.ne.0)
*   ..compute ue, electron parallel velocity
     > call b2xpve (nx, ny, ns, qe, rza, rz2, alfx_c, sigx_c, po,        !srv 05.07.17
     &              ne, te, na, ua, ue_ext, ne_ext, ue)
      if (ue_style.ne.2 .or. iout.ne.0) then
        do iy=-1,ny
          do ix=-1,nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
              wrk0(ix,iy)=-(fch_p(ix,iy,0)/pbs(ix,iy,0)+
     &          fch_p(rightix(ix,iy),rightiy(ix,iy),0)/
     &          pbs(rightix(ix,iy),rightiy(ix,iy),0))/qe/2.0_R8
            elseif(leftix(ix,iy).eq.-2) then
              wrk0(ix,iy)=-fch_p(rightix(ix,iy),rightiy(ix,iy),0)/
     &         pbs(rightix(ix,iy),rightiy(ix,iy),0)/qe
            elseif(rightix(ix,iy).eq.nx+1) then
              wrk0(ix,iy)=-fch_p(ix,iy,0)/pbs(ix,iy,0)/qe
            endif
!           wrk0(ix,iy)=wrk0(ix,iy)*prl_cur                              !srv 03.04.15
          enddo
        enddo
        do is=0,ns-1
          wrk0(:,:)=wrk0(:,:)+rza(:,:,is)*ua(:,:,is)*na(:,:,is)
        enddo
        wrk0(:,:)=wrk0(:,:)+ue_ext(:,:)*ne_ext(:,:)
        wrk0(:,:)=wrk0(:,:)/ne(:,:)                                      !srv 05.03.11
      endif
      if(iout.ne.0) then                                                 !srv 05.07.17
       call b2xzef (nx, ny, ns, rz2, na, ne, zeff)
       call b2xpic (nx, ny, ns, qe, rza, bb, na, ua, zeff, fch_pi_c)
       call b2xpif (nx, ny, hx, vol, fch_pi_c, fch_pi_f)
       call my_out(70,nx,ny,fch_pi_c,'b2news_m_fch_pi_c')
       call my_out(70,nx,ny,fch_pi_f,'b2news_m_fch_pi_f')
       call my_out(70,nx,ny,ue,'b2news_m_ve')
       call my_out(70,nx,ny,wrk0,'b2news_m_ue')
      endif                                                              !srv 05.07.17
      if (ue_style.ne.2) ue=wrk0

csw 11oct2012 added check on Coriolis force (taken from b2news)
      ier0=0
      call b2npmo (nx, ny, ns, ismain, ncall, solving(1), rxf1, dtim,    !srv 13.01.17 24.05.21
     &   vol, hx, hz, bb, pbshz, lnlam,                                  !srv 23.07.99 09.01.01 01.07.09
     &   cvsahz, cvsa, cvsahz_cl, cvsahz_drho, kt_neo,                   !srv 19.07.02
     &   alpha1, chci_al_ast, chci_a, cvsahz_hAdp_albe, cvsahz_hBdp_al,  !som 08.10.21
     *   cvsahz_uAdp_albe, cvsahz_uBdp_al,                               !som 22.12.23
     &   smo, smq, cthe, cthi, chci, smof,                               !srv 23.07.99 01.10.99 04.10.99 09.01.01
     &   resmo, resmt, corua, corut, pcca, pccm, ier0)                   !srv 23.07.99 10.11.02 23.12.03

      call xertst (ier0.eq.0, 'error return from b2npmo')

      if (pot_eq.eq.0) then                                              !srv 13.07.05 {
        po   =0.0e0_R8
      else if (pot_eq.eq.2) then                                         !srv 11.01.12
        po(:,:)=3.1_R8*te(:,:)/qe
      else                                                               !srv 13.07.05 }
*   ..correct po
        if (solving(3)) then
          call b2nxdp (nx, ny, nnreg(0), solvepo, qe, te, nep, ne, wrk0)
          call b2saxpy (n2, 1.0_R8, wrk0, 1, po, 1)
          if(iout .ne. 0) then                                           !srv 12.07.05 {
            call my_out(70,nx,ny,po,'b2nxdp_po')
          endif                                                          !srv 12.07.05 }
        endif
*   ..compute currents
        call b2tfch_ (nx, ny, ns, nscx, ismain, iscx,
     &    facdrift, fac_ExB, fac_vis,                                    !srv 09.01.20
     &    vol, gs, hx, hy1, hz, qz, qc, bb, te, po, ne,
     1    ne2, lnlam, pz, ti, ni(-1,-1,1), na, fna, uadia, ua,
     1    csig, csigin, csig_stoch, chci, calf, cvsa_cl, cvsa_drho,      !som 18.08.21 !srv 24.08.21
     &    chci_al_ast, chci_a, cvsa_hAdp_albe, cvsa_hBdp_al,             !som 20.10.21
     &    cvsa_uAdp_albe, cvsa_uBdp_al,                                  !som 26.04.24
     &    conc, avisper_a, bvisper_a, cvisper_a,
     &    fchanml_a, fchinert_a, fchvispar_a, fchvisper_a, fchvisq_a,    !srv 08.09.21
     &    fchin_a, fch, fch_p, fchdia, fchin,                            !IYS 17.10.2023
     &    fchvispar, fchvisper, fchvisq, fchinert, fchanml, fchstoch)
        if (solving(3)) then
*   ..solve the potential equation
          if (vis_per.eq.0.0_R8) then                                    !srv 22.12.19 06.07.00
            new_matrix = .false.
            do ireg = 0, nnreg(0)
              new_matrix = new_matrix.or.
     &         solvepo(ireg).neqv.last_solve_5(ns,ireg)
            end do
            if (new_matrix.or.in_parallel()) call restart_ma28
            call b2nppo (nx, ny, nnreg(0), solving(3), solvepo, ncall,   !srv 22.05.18
     &        rxf2, ne, te, conc, sch, po, fch, respo, corpo, ier0)
            last_solve_5(ns,0:nnreg(0))=solvepo(0:nnreg(0))
          else
            call b2npp7 (nx, ny, ns, ismain, nnreg(0), solving(3),       !srv 27.01.00 22.05.18
     &       solvepo, ncall, rxf2, ne, te, conc, avisper_a, bvisper_a,
     &       cvisper_a, sch, po, fch, respo, corpo, ier0)
          endif                                                          !srv 06.07.00
          t0=0.0_R8
crfs      do iy = 0, ny-1
          do iy = -1, ny
crfs        t0 = t0+b2sasum(nx,respo(0,iy),1)
            t0 = t0+b2sasum(nx+2,respo(-1,iy),1)
          enddo
          write(*,*) 'respo before solving',t0
          call xertst (ier0.eq.0, 'error return from b2nppo')
*   ..recompute currents
          call b2tfch_ (nx, ny, ns, nscx, ismain, iscx,
     &      facdrift, fac_ExB, fac_vis,                                  !srv 09.01.20
     &      vol, gs, hx, hy1, hz, qz, qc, bb, te, po, ne,
     1      ne2, lnlam, pz, ti, ni(-1,-1,1), na, fna, uadia, ua,
     1      csig, csigin, csig_stoch, chci, calf, cvsa_cl, cvsa_drho,    !som 18.08.21 !srv 24.08.21
     &      chci_al_ast, chci_a, cvsa_hAdp_albe, cvsa_hBdp_al,           !som 20.10.21
     &      cvsa_uAdp_albe, cvsa_uBdp_al,                                !som 26.04.24
     &      conc, avisper_a, bvisper_a, cvisper_a,
     &      fchanml_a, fchinert_a, fchvispar_a, fchvisper_a, fchvisq_a,  !srv 08.09.21
     &      fchin_a, fch, fch_p, fchdia, fchin,                          !IYS 17.10.2023
     &      fchvispar, fchvisper, fchvisq, fchinert, fchanml, fchstoch)
*   ..compute new residual following the tradition
          call b2ursd (nx, ny, ne, po, sch, fch, wrk0)
          t0=0.0_R8
crfs      do iy = 0, ny-1
          do iy = -1, ny
crfs        t0 = t0+b2sasum(nx,wrk0(0,iy),1)
            t0 = t0+b2sasum(nx+2,wrk0(-1,iy),1)
          enddo
          write(*,*) 'respo after solving',t0
        else
          if(get_residuals.eq.1)
     &       call b2ursd (nx, ny, ne, po, sch, fch, respo)
        endif
      endif           ! if pot_eq                                        !srv 13.07.05

      if(iout .ne. 0) then                                               !srv 12.07.05 {
        write (charns,'(i3.3)') i
        call my_out(70,nx,ny,respo,'b2news_m_respo'//charns)
        call my_out(70,nx,ny,po,'b2news_m_po'//charns)
      endif                                                              !srv 12.07.05 }

      if(ion_vlct_restrict .ne. 0 ) then                                 !iyv 06.02.14
        call b2xprz (nx, ny, ns, mp, am, na, rz)
        call b2xppz (nx, ny, ns, ne, na, te, ti, pz)
        wrk0(:,:)=sqrt(pz(:,:)/rz(:,:))*abs(bb(:,:,0)/bb(:,:,3))
        call intfaceh(nx,ny,leftix,leftiy,vol,wrk0,cssb(-1,-1,0))
        call intfacev(nx,ny,bottomix,bottomiy,vol,wrk0,cssb(-1,-1,1))
      else
        cssb = 0.0_R8
      endif                                                              !iyv 06.02.14

#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2NEWS_UNDERSCORE_LOOP3) )
#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
!! FR possible loop (3) for use with OpenMP
!$OMP PARALLEL DO DEFAULT(none) LASTPRIVATE(flob, conb, pb)
!$OMP& SHARED(nx, ny, rza, na, te, ti, ns, nsmin, nsmax, ismain)
!$OMP& SHARED(zamin, zamax)
!$OMP& SHARED(guard_flows, facdrift, fac_ExB, am, cssb)
!$OMP& SHARED(vol, pbs, pbshz, ua, uadia, vaecrb, vadia, wadia)
!$OMP& SHARED(ne, bb, hx, hy, hz, qz, gs, qc, po, cvla, cdna, cdpa)
!$OMP& SHARED(cddi, cdde, cvlahz, cdnahz, cdpahz, fchvispar)
!$OMP& SHARED(fchvisper, fchinert, fchanml, fchvisq, fchin)
!$OMP& SHARED(fchvispar_a,fchvisper_a, fchinert_a, fchanml_a, fchvisq_a)
!$OMP& SHARED(fchin_a, fna, fna_nodrift, fna_fcor, fna_he, fnaPSch)
!$OMP& SHARED(fna_mdf, fna_32, fna_52, fna_32nd, fna_52nd, fna_53)
!$OMP& SHARED(kinrgy, fllim0fna)
#endif
      do is = 0, ns-1
*    ..compute partial pressure
        call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)
*    ..re-compute particle flux
#ifdef DBG
        write(*,*) 'is, za(is) ',is, (zamin(is)+zamax(is))/2.0_R8
#endif
        call b2tfnb (nx, ny, is, ismain,
     &     guard_flows, facdrift, fac_ExB, am(is)*mp, vol, cssb,         !srv 07.06.00 09.01.01 08.09.01 19.05.03 !iyv 06.02.14
     1     rza(-1,-1,is), pbs, pbshz, na(-1,-1,is), pb, ua(-1,-1,is),    !srv 09.01.01 20.09.06
     1     uadia(-1,-1,0,is), vaecrb(-1,-1,0,is),
     1     vadia(-1,-1,0,is),                                            !srv 17.05.99
     1     wadia(-1,-1,0,is), ne, te, ti, bb, hx, hy, hz, gs, qc,        !srv 19.05.03 08.10.04 04.02.05
     &     po, cvla(-1,-1,0,is), cdna(-1,-1,0,is),
     &     cdpa(-1,-1,0,is), cddi(-1,-1,0,is),                           !srv 13.10.06
     &     cvlahz(-1,-1,0,is),                                           !srv 09.01.01 15.12.05
     &     cdnahz(-1,-1,0,is), cdpahz(-1,-1,0,is),
     &     fchvispar_a(-1,-1,0,is), fchvisper_a(-1,-1,0,is),             !srv 08.09.21
     &     fchinert_a(-1,-1,0,is), fchanml_a(-1,-1,0,is),
     &     fchvisq_a(-1,-1,0,is), fchin_a(-1,-1,0,is),                   !srv 20.09.06, IYS 17.10.2023
     &     fna(-1,-1,0,is), fna_nodrift(-1,-1,0,is),
     &     fna_fcor(-1,-1,0,is),
     &     fna_he(-1,-1,0,is), fnaPSch(-1,-1,0,is),
     &     fna_mdf(-1,-1,0,is),                                          !srv 13.10.06
     &     flob, conb, fna_32(-1,-1,0,is), fna_52(-1,-1,0,is),
     &     fna_32nd(-1,-1,0,is), fna_52nd(-1,-1,0,is),
     &     fna_53(-1,-1,0,is),
     &     kinrgy(-1,-1,is), fllim0fna(-1,-1,0,is),
     &     .false.)                                                      !djm Jan2017
      enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2NEWS_UNDERSCORE_LOOP3) )
!$OMP END PARALLEL DO
!!FR End of loop (3)
#ifdef WITH_FTIMINGS
      call timer%stop("parallel")
#endif
#endif

*   ..compute fne, fni
c      with diamagnetic drift part                                       !srv 13.05.99
      call b2xpfe (nx, ny, ns, ns_ext, qe,
     . rza, za_ext, fna_he, fa_ext, fch_p, fne)                          !srv 13.05.99 21.01.00 11.10.02
      call b2xpfe (nx, ny, ns, ns_ext, qe,
     . rza, za_ext, fna_53, fa_ext, fch_p, fne_53)                       !ajs 06.04.20
      call b2xpfe (nx, ny, ns, ns_ext, qe,
     1 rza, za_ext, fna_32nd, fa_ext, fch_32, fne_32)
      call b2xpfe (nx, ny, ns, 0, qe,
     1 rza, zero, fna_52nd, zero, fch_52, fne_52)
      call b2xpfi (nx, ny, ns, fna_he, fni)                              !srv 13.05.99 21.01.00
      call b2xpfi (nx, ny, ns, fna_32nd, fni_32)                         !xpb
      call b2xpfi (nx, ny, ns, fna_52nd, fni_52)                         !xpb
*   ..compute electron drift velocities
      call b2tfed (nx, ny, facdrift, fac_ExB, te, po, vedia, veecrb)
*   ..compute heat fluxes
      call b2tfhe_ (nx, ny, ns, ismain, facdrift,
     & vol, hx, hy1,                                                     !srv 29.01.20
     1 te, po, ne, na, vedia, fne, chce, chve, calf,                     !srv 11.07.99 18.01.02 01.07.09
     2 chvemx, cdde, fch_p, fchstoch,                                    !som 03.11.21
     & floe, cone, fhe, fhe_mdf, fhePSch, .true.)                        !djm Jan2017
csw 30apr2014 added fna_he, 17sep2014 added ua,cvsa
      call b2tfhi_ (nx, ny, ns, facdrift, vol, hx, hy1, bb,              !som 13.07.21
     1  ti, ni, na, fni, fna_he, ua,
     &  chci, cddi, chvi, chvimx, cvsa, vadia,                           !srv 17.05.99 27.12.01 25.01.02 13.10.06
     &  fchvispar, fchvisper, fchvisq, fchinert, fchanml,
     &  fhi, fhi_mdf, fhiPSch, floi, coni,                               !srv 13.10.06
     &  .true.)                                                          !djm Jan2017
       floi_noc = floi                                                   !srv 24.10.02
       floe_noc = floe                                                   !srv 24.10.02
      call b2xpen (nx, ny, ns,
     1 fna, fna_53, fne_53, fch, fhe, fhi, kinrgy, ua, cvsa,
     2 rpt, te, ti, po, BoRiS, fhm, fhp, fhj, fht)
      if (iout_b2wdat.eq.4) then
       call my_out(70,nx,ny,fht(-1,-1,0),'b2news_m_fhtx')
       call my_out(70,nx,ny,fht(-1,-1,1),'b2news_m_fhty')
      endif
*   ..add contribution from cdpa to floe, cone, floi, coni
      call b2tfcc (nx, ny, ns, vol, na, te, ti, cdpa,
     &  floe, cone, floi, coni)
      if(ue_style.ne.2) then
*   ..compute ue, electron parallel velocity                             !srv 05.03.11 {
       do iy=-1,ny
         do ix=-1,nx
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
             ue(ix,iy)=-(fch_p(ix,iy,0)/pbs(ix,iy,0)+
     &        fch_p(rightix(ix,iy),rightiy(ix,iy),0)/
     &          pbs(rightix(ix,iy),rightiy(ix,iy),0))/qe/2.0_R8
           elseif(leftix(ix,iy).eq.-2) then
             ue(ix,iy)=-fch_p(rightix(ix,iy),rightiy(ix,iy),0)/
     &         pbs(rightix(ix,iy),rightiy(ix,iy),0)/qe
           elseif(rightix(ix,iy).eq.nx+1) then
             ue(ix,iy)=-fch_p(ix,iy,0)/pbs(ix,iy,0)/qe
           endif
         enddo
       enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2NEWS_UNDERSCORE_LOOP4) )
#ifdef WITH_FTIMINGS
       call timer%start("parallel")
#endif
!!FR loop (4) over ns suitable for OpenMP parallelisation
!$OMP PARALLEL DO DEFAULT(none) SHARED(rza,ua,na,ns) REDUCTION(+:ue)
#endif
       do is=0,ns-1
        ue(:,:)=ue(:,:)+rza(:,:,is)*ua(:,:,is)*na(:,:,is)
       enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2NEWS_UNDERSCORE_LOOP4) )
!$OMP END PARALLEL DO
!!FR end of loop (4)
#ifdef WITH_FTIMINGS
       call timer%stop("parallel")
#endif
#endif
       ue(:,:)=ue(:,:)+ue_ext(:,:)*ne_ext(:,:)
       ue(:,:)=ue(:,:)/ne(:,:)                                           !srv 05.03.11 }
      else
       call b2xpve (nx, ny, ns, qe, rza, rz2, alfx_c, sigx_c, po,        !srv 17.11.16 {
     &               ne, te, na, ua, ue_ext, ne_ext, ue)
      endif

*   ..iterate on heat balance
      if(mdf_fhe.eq.0 .and. mdf_fhi.eq.0) then                           !srv 13.10.06
       call b2npht (nx, ny, ns, nnreg(0), ncall, ismain, solving,        !srv 05.07.17
     &  solveee, solveei, solveet, solvepo, solvemo,                     !srv 22.05.18
     &  rxf3, dtim, BoRiS, fac_ExB, hx, hy, hy1, gs, qc, vol, bb,        !srv 07.06.00 02.06.99 24.05.21
     &  pbs, lnlam, floe, cone, floi, coni, ceqp, pccm, she, shi,        !srv 07.06.00 01.07.09
     &  fne, fch, fchstoch, cvsa, cvsa_cl, cvsa_drho,                    !som 16.08.21 !srv 24.08.21
     &  cvsa_uAdp_albe, cvsa_uBdp_al,                                    !som 08.05.24
     &  fchdia, uadia, wadia, vadia, vaecrb,                             !srv 27.12.01 28.06.03 30.06.03 07.07.09 03.11.17
     &  na, ua, ue, te, ti, po, ne, ni, fhe, fhi,                        !srv 31.05.99
     &  shef, shif, reshe, reshi, wrk0, corte, corti, ier0)              !srv 31.05.99 10.11.02
      endif                                                              !srv 13.10.06 }
      if(mdf_fhe.ne.0 .and. mdf_fhi.ne.0) then                           !srv 13.10.06
        call b2nph9 (nx, ny, ns, nnreg(0), ncall, ismain, solving,       !srv 05.07.17 22.05.18
     &   solveee, solveei, solveet, solvepo, solvemo,                    !srv 22.05.18
     &   rxf3, dtim, BoRiS, fac_ExB, hx, hy, hy1, gs, qc, vol, bb,       !srv 07.06.00 02.06.99 24.05.21
     &   pbs, lnlam, floe, cone, floi, coni, ceqp, pccm, she, shi,       !srv 07.06.00 01.07.09
     &   fne, fch, fchstoch, cvsa,                                       !srv 24.08.21
     &   cvsa_cl, cvsa_drho,                                             !som 16.08.21
     &   cvsa_uAdp_albe, cvsa_uBdp_al,                                   !som 08.05.24
     &   fchdia, uadia, wadia, vadia, vaecrb, na, ua, ue,                !srv 27.12.01 28.06.03 30.06.03 07.07.09 03.11.17
     &   te, ti, po, ne, ni, fhe_mdf, fhi_mdf, shef, shif,               !srv 31.05.99
     &   reshe, reshi, wrk0, corte, corti, ier0)                         !srv 31.05.99 10.11.02
      endif                                                              !srv 13.10.06 }
      if(mdf_fhe.eq.0 .and. mdf_fhi.ne.0) then                           !srv 13.10.06
       call b2npht (nx, ny, ns, nnreg(0), ncall, ismain, solving,        !srv 05.07.17
     &  solveee, solveei, solveet, solvepo, solvemo,                     !srv 22.05.18
     &  rxf3, dtim, BoRiS, fac_ExB, hx, hy, hy1, gs, qc, vol, bb,        !srv 07.06.00 02.06.99 24.05.21
     &  pbs, lnlam, floe, cone, floi, coni, ceqp, pccm, she, shi,        !srv 07.06.00 01.07.09
     &  fne, fch, fchstoch, cvsa, cvsa_cl, cvsa_drho,                    !som 16.08.21 !srv 24.08.21
     &  cvsa_uAdp_albe, cvsa_uBdp_al,                                    !som 08.05.24
     &  fchdia, uadia, wadia, vadia, vaecrb,                             !srv 27.12.01 28.06.03 30.06.03 07.07.09 03.11.17
     &  na, ua, ue, te, ti, po, ne, ni, fhe, fhi_mdf,                    !srv 31.05.99
     &  shef, shif, reshe, reshi, wrk0, corte, corti, ier0)              !srv 31.05.99 10.11.02
      endif                                                              !srv 13.10.06 }
      if(mdf_fhe.ne.0 .and. mdf_fhi.eq.0) then                           !srv 13.10.06
       call b2npht (nx, ny, ns, nnreg(0), ncall, ismain, solving,        !srv 05.07.17
     &  solveee, solveei, solveet, solvepo, solvemo,                     !srv 22.05.18
     &  rxf3, dtim, BoRiS, fac_ExB, hx, hy, hy1, gs, qc, vol, bb,        !srv 07.06.00 02.06.99 24.05.21
     &  pbs, lnlam, floe, cone, floi, coni, ceqp, pccm, she, shi,        !srv 07.06.00 01.07.09
     &  fne, fch, fchstoch, cvsa, cvsa_cl, cvsa_drho,                    !som 16.08.21 !srv 24.08.21
     &  cvsa_uAdp_albe, cvsa_uBdp_al,                                    !som 08.05.24
     &  fchdia, uadia, wadia, vadia, vaecrb,                             !srv 27.12.01 28.06.03 30.06.03 07.07.09 03.11.17
     &  na, ua, ue, te, ti, po, ne, ni, fhe_mdf, fhi,                    !srv 31.05.99
     &  shef, shif, reshe, reshi, wrk0, corte, corti, ier0)              !srv 31.05.99 10.11.02
      endif                                                              !srv 13.10.06 }

csw 16aug2011 init potential
      if(pot_eq.eq.0) po=0.0_R8                                          !srv 13.07.05
      if(pot_eq.eq.2) po(:,:)=3.1_R8*te(:,:)/qe
      call xertst (ier0.eq.0, 'error return from b2npht')

csw 12oct2012, extra call to relax particle balance if wanted (as in SOLPS4)
      if(do_2nd_b2npco_call.ne.0) then

c   recompute particle fluxes
       do is=0,ns-1                                                      !srv 01.10.99
*    ..DPC
        do iy=-1,ny
          do ix=-1,nx
            na(ix,iy,is)=max(na(ix,iy,is),na_min)
          enddo
        enddo
*    ..compute partial pressure
        call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)
*    ..compute particle flux
        call b2tfnb (nx, ny, is, ismain,
     &   guard_flows, facdrift, fac_ExB, am(is)*mp, vol, cssb,           !srv 07.06.00 09.01.01 08.09.01 19.05.03 06.02.14
     1   rza(-1,-1,is), pbs, pbshz, na(-1,-1,is), pb, ua(-1,-1,is),      !srv 09.01.01 20.09.06
     1   uadia(-1,-1,0,is), vaecrb(-1,-1,0,is), vadia(-1,-1,0,is),       !srv 17.05.99
     1   wadia(-1,-1,0,is), ne, te, ti, bb, hx, hy, hz, gs, qc,          !srv 19.05.03 08.10.04  04.02.05
     &   po, cvla(-1,-1,0,is), cdna(-1,-1,0,is),
     &   cdpa(-1,-1,0,is), cddi(-1,-1,0,is),                             !srv 13.10.06
     &   cvlahz(-1,-1,0,is),                                             !srv 09.01.01 15.12.05
     &   cdnahz(-1,-1,0,is), cdpahz(-1,-1,0,is),
     &   fchvispar_a(-1,-1,0,is), fchvisper_a(-1,-1,0,is),               !srv 08.09.21
     &   fchinert_a(-1,-1,0,is), fchanml_a(-1,-1,0,is),
     &   fchvisq_a(-1,-1,0,is), fchin_a(-1,-1,0,is),                     !srv 20.09.06, IYS 17.10.2023
     &   fna(-1,-1,0,is), fna_nodrift(-1,-1,0,is),
     &   fna_fcor(-1,-1,0,is),
     &   fna_he(-1,-1,0,is), fnaPSch(-1,-1,0,is),
     &   fna_mdf(-1,-1,0,is),                                            !srv 13.10.06
     &   flob, conb, fna_32(-1,-1,0,is), fna_52(-1,-1,0,is),
     &   fna_32nd(-1,-1,0,is), fna_52nd(-1,-1,0,is),
     &   fna_53(-1,-1,0,is),
     &   kinrgy(-1,-1,is), fllim0fna(-1,-1,0,is),
     &   .true.)                                                         !djm Jan2017

*    ..perform iteration
        new_matrix = .false.
        if(mdf_fnb.eq.0) then                                            !srv 13.10.06
!         if (solving(2).and.all(solveco(is,0:nnreg(0)))) then           !srv  22.05.18
            solvedum(0:nnreg(0))=solveco(is,0:nnreg(0))
            solvedm2(0:nnreg(0))=solvemo(is,0:nnreg(0))
            if((visper_all_ions.ne.0 .and. .not.is_neutral(is)) .or.
     &        (visper_all_ions.eq.0 .and. is.eq.ismain)) then
              do ireg = 0, nnreg(0)
                new_matrix = new_matrix.or.
     &           solveco(is,ireg).neqv.last_solve_7(is,ireg)
              end do
              if(new_matrix.or.in_parallel()) call restart_ma28_for_7
              call b2npc7 (nx, ny, fac_vis, guard_flows, rza(-1,-1,is),
     &         nnreg(0), solving(2), solvedum, solving(1),               !srv 22.05.18
     &         solvedm2, ncall, is, rxf0, dtim, vol, pbs, flob, conb,    !srv 22.01.02 01.07.05 24.05.21
     &         sna(-1,-1,0,is), pcca(-1,-1,0,is),
     &         ti, pb, na(-1,-1,is), ua(-1,-1,is), fna(-1,-1,0,is),
     &         resco(-1,-1,is), corpa(-1,-1,is), ier0)
              call xertst (ier0.eq.0, 'error return from b2npc7')
              last_solve_7(is,0:nnreg(0))=solveco(is,0:nnreg(0))
            else                                                         !srv 02.04.18
              do ireg = 0, nnreg(0)
                new_matrix = new_matrix.or.
     &           solveco(is,ireg).neqv.last_solve_5(is,ireg)
              end do
              if(new_matrix.or.in_parallel()) call restart_ma28
              call b2npco (nx, ny, guard_flows, nnreg(0), solving(2),    !srv 22.05.18
     &         solvedum, solving(1), solvedm2, ncall, is, rxf0, dtim,    !srv 22.01.02 01.07.05
     &         vol, pbs, flob, conb, sna(-1,-1,0,is),
     &         pcca(-1,-1,0,is),                                         !srv 24.05.21
     &         pb, na(-1,-1,is), ua(-1,-1,is), fna(-1,-1,0,is),
     &         resco(-1,-1,is), corpa(-1,-1,is), ier0)
              call xertst (ier0.eq.0, 'error return from b2npco')
              last_solve_5(is,0:nnreg(0))=solveco(is,0:nnreg(0))
            endif
!         endif                                                          !srv 22.05.18
        else                                                             !srv 13.10.06 {
!         if (solving(2).and.all(solveco(is,0:nnreg(0)))) then           !srv 22.05.18
            solvedum(0:nnreg(0))=solveco(is,0:nnreg(0))
            solvedm2(0:nnreg(0))=solvemo(is,0:nnreg(0))
            if(.not.is_neutral(is)) then
              if(visper_all_ions.ne.0 .or. is.eq.ismain) then
                do ireg = 0, nnreg(0)
                  new_matrix = new_matrix.or.
     &             solveco(is,ireg).neqv.last_solve_11(is,ireg)
                end do
                if(new_matrix.or.in_parallel()) call restart_ma28_for_11
                call b2npc11 (nx, ny, fac_vis, guard_flows,
     &           rza(-1,-1,is), nnreg(0), solving(2), solvedum,          !srv 22.05.18
     &           solving(1), solvedm2, ncall, is, rxf0, dtim,            !srv 24.05.21
     &           vol, pbs, flob, conb, sna(-1,-1,0,is),
     &           pcca(-1,-1,0,is), po, ti, pb, na(-1,-1,is),
     &           ua(-1,-1,is), fna_mdf(-1,-1,0,is),
     &           resco(-1,-1,is), corpa(-1,-1,is), ier0)
                call xertst (ier0.eq.0, 'error return from b2npc11')
                last_solve_11(is,0:nnreg(0))=solveco(is,0:nnreg(0))
              else                                                       !srv 02.04.18
                do ireg = 0, nnreg(0)
                  new_matrix = new_matrix.or.
     &             solveco(is,ireg).neqv.last_solve_9(is,ireg)
                end do
                if(new_matrix.or.in_parallel()) call restart_ma28_for_9
                call b2npc9 (nx, ny, guard_flows, nnreg(0), solving(2),  !srv 22.05.18
     &           solvedum, solving(1), solvedm2, ncall, is, rxf0, dtim,  !srv 24.05.21
     &           vol, pbs, flob, conb,
     &           sna(-1,-1,0,is), pcca(-1,-1,0,is), pb, na(-1,-1,is),
     &           ua(-1,-1,is), fna_mdf(-1,-1,0,is),
     &           resco(-1,-1,is), corpa(-1,-1,is), ier0)
                call xertst (ier0.eq.0, 'error return from b2npc9')
                last_solve_9(is,0:nnreg(0))=solveco(is,0:nnreg(0))
              endif
            else
              do ireg = 0, nnreg(0)
                new_matrix = new_matrix.or.
     &           solveco(is,ireg).neqv.last_solve_5(is,ireg)
              end do
              if(new_matrix.or.in_parallel()) call restart_ma28
              call b2npco (nx, ny, guard_flows, nnreg(0), solving(2),    !srv 22.05.18
     &         solvedum, solving(1), solvedm2, ncall, is, rxf0, dtim,    !srv 24.05.21
     &         vol, pbs, flob, conb,
     &         sna(-1,-1,0,is), pcca(-1,-1,0,is), pb, na(-1,-1,is),
     &         ua(-1,-1,is), fna_mdf(-1,-1,0,is),
     &         resco(-1,-1,is), corpa(-1,-1,is), ier0)
              last_solve_5(is,0:nnreg(0))=solveco(is,0:nnreg(0))
              call xertst (ier0.eq.0, 'error return from b2npco')
            endif
!         endif                                                          !srv 22.05.18
        endif                                                            !srv 13.10.06 }
       enddo
*      ..compute electron drift velocities
       call b2tfed (nx, ny, facdrift, fac_ExB, te, po, vedia, veecrb)
      endif ! do_2nd_b2npco

csw 06dec2012 finally re-evaluate particle fluxes
c   (necessary to work out new rescaling factors sclrtio/sclalph in neutr_scl)
      if(re_eval_prtls_fluxes.eq.1 .or. use_eirene.ne.0) then            !srv 28.06.13 {
       do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            na(ix,iy,is)= max(na(ix,iy,is),na_min)   !abs(na(ix,iy,is))  !srv 24.04.09
          enddo
        enddo
*    ..compute partial pressure
        call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)    !srv 01.10.99
*    ..compute particle flux
        call b2tfnb (nx, ny, is, ismain,
     &   guard_flows, facdrift, fac_ExB, am(is)*mp, vol, cssb,           !srv 07.06.00 09.01.01 08.09.01 19.05.03 !iyv 06.02.14
     1   rza(-1,-1,is), pbs, pbshz, na(-1,-1,is), pb, ua(-1,-1,is),      !srv 09.01.01 20.09.06
     1   uadia(-1,-1,0,is), vaecrb(-1,-1,0,is), vadia(-1,-1,0,is),
     1   wadia(-1,-1,0,is), ne, te, ti, bb, hx, hy, hz, gs, qc,          !srv 19.05.03 08.10.04 04.02.05
     &   po, cvla(-1,-1,0,is), cdna(-1,-1,0,is),
     &   cdpa(-1,-1,0,is), cddi(-1,-1,0,is),                             !srv 13.10.06
     &   cvlahz(-1,-1,0,is),                                             !srv 09.01.01 15.12.05
     &   cdnahz(-1,-1,0,is), cdpahz(-1,-1,0,is),
     &   fchvispar_a(-1,-1,0,is), fchvisper_a(-1,-1,0,is),               !srv 08.09.21
     &   fchinert_a(-1,-1,0,is), fchanml_a(-1,-1,0,is),
     &   fchvisq_a(-1,-1,0,is), fchin_a(-1,-1,0,is),                     !srv 20.09.06, IYS 17.10.2023
     &   fna(-1,-1,0,is), fna_nodrift(-1,-1,0,is),
     &   fna_fcor(-1,-1,0,is),
     &   fna_he(-1,-1,0,is), fnaPSch(-1,-1,0,is),
     &   fna_mdf(-1,-1,0,is),                                            !srv 13.10.06
     &   flob, conb, fna_32(-1,-1,0,is), fna_52(-1,-1,0,is),
     &   fna_32nd(-1,-1,0,is), fna_52nd(-1,-1,0,is),
     &   fna_53(-1,-1,0,is),
     &   kinrgy(-1,-1,is), fllim0fna(-1,-1,0,is),
     &   .false.)                                                        !djm Jan2017
       enddo                                                             !srv 01.10.99
*      ..compute electron drift velocities
       call b2tfed (nx, ny, facdrift, fac_ExB, te, po, vedia, veecrb)
      endif                                                              !srv 28.06.13 }

csw 26oct2011 recalculate rescaling factors sclrtio/sclalph (SOLPS4)
      if(no_neutr_scl .eq. 0) then
        write(*,*) 'Calling b2stbr_neutr_scl'
        call b2stbr_neutr_scl(nx, ny, ns)                                !srv 25.10.08
        write(*,*) 'Called b2stbr_neutr_scl'
      endif


* ..return
      ierr = 0

      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2news_m

      end subroutine b2news_m
