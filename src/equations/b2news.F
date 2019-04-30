*-----------------------------------------------------------------------
*.specification

      subroutine b2news (nx, ny, ns, nxtl, nxtr,
     &  nscx, iscx, ismain, ismain0,
     &  dtim, rxf, ierr, wrong_flow, ramp_call)
      use b2mod_types
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_anomalous_transport
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_numerics_namelist
      use b2mod_work
      use b2mod_b2cmpa
      use b2mod_external
      use b2mod_neoclassical
#ifdef WITH_FTIMINGS
      use timings
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nxtl, nxtr
      integer nscx, iscx(0:nscxmax-1), ismain, ismain0
      real (kind=R8) ::
     *   dtim, rxf
*   ..output arguments                                                   !xpb
      integer ierr
      logical wrong_flow, ramp_call                                      !xpb
*   ..common blocks

*-----------------------------------------------------------------------
*
*  1. purpose
*
*     B2NEWS performs one inner iteration for the b2 equations.
*     This routine defines the core structure of the solution
*     algorithm.
*
*
*  3. description (see also routine b2cdca)
*
*     For a detailed description of the equations and the method of
*     solution, see the references cited in the main program unit.
*
*     This routine performs one basic iteration for the b2 equations.
*     An outline of the procedure follows.
*     ..test input data
*     ..compute auxiliary quantities
*     ..compute transport coefficients
*     ..iterate on parallel momentum balance
*     ..iterate on particle balance
*     ..iterate on charge balance
*     ..iterate on heat balance
*     Convergence monitoring must be done outside this routine.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer potit, ipass, recalculate_contributions, pot_eq,           !xpb,dpc
     *  potitmin, nsmin, nsmax, no_current
      integer, save :: style = 3                                         !srv 29.05.02
      integer, save :: mdf_fnb = 0, mdf_fhe = 0, mdf_fhi = 0             !srv 13.10.06
      integer, save :: use_b2npc9 = 1, use_b2nph9 = 1                    !srv 06.11.08 11.08.08
      integer, save :: iout = 0
      integer, save :: ion_vlct_restrict = 0                             !iyv 06.02.14
      integer, save :: ue_style = 2                                      !srv 17.11.16 13.01.17
      real (kind=R8) ::
     *  potok, mach10, BoRiS, vis_per, na_min, na_new
      integer ncall, n2, is, ier0, i, ix, iy, ncallout, ixmax, iymax,    !xpb
     1  guard_flows, ramp_slow, no_solve, ireg
      real (kind=R8) ::
     *  xfm0, xfm1, xfm2, xfm3, rxf0, rxf1, rxf2, rxf3, t0, tmax
      real (kind=R8) ::                                                  !xpb
     1  sna0(-1:nx,-1:ny,0:1,0:ns-1),                                    !xpb
     1  sne0(-1:nx,-1:ny,0:1),                                           !xpb
     2  shi0(-1:nx,-1:ny,0:3), she0(-1:nx,-1:ny,0:3),                    !xpb
     3  smo0(-1:nx,-1:ny,0:3,0:ns-1),                                    !xpb
     4  smq0(-1:nx,-1:ny,0:3,0:ns-1),                                    !xpb
     5  sch0(-1:nx,-1:ny,0:3),                                           !xpb
     6  ue(-1:nx,-1:ny), zeff(-1:nx,-1:ny), pz(-1:nx,-1:ny),             !srv 05.07.17
     7  lnlam(-1:nx,-1:ny), zero(-1:nx,-1:ny),
     &  rz(-1:nx,-1:ny), cssb(-1:nx,-1:ny,0:1),                          !iyv 06.02.14
     &  fchstoch(-1:nx,-1:ny)                                            !srv 17.12.13
      character charnc*6, hlp_frm*4
      logical solving(4), solvedum(0:NNREGMAX), solvedm2(0:NNREGMAX)
      logical new_matrix
      save ncall, xfm0, xfm1, xfm2, xfm3, potit, potok, pot_eq,
     1  recalculate_contributions, ramp_slow, guard_flows,
     3  potitmin, BoRiS, vis_per, no_solve, solving, na_min, na_new,
     *  nsmin, nsmax, no_current
*   ..procedures
      intrinsic abs, min, mod
      real (kind=R8) :: smin, smax
      integer len_of_digits
      external len_of_digits
      external subini, subend, xertst, ipgetr, smin, smax,
     &  b2scopy, b2saxpy, ipgeti
      external b2xvsg, b2xvps, b2xpne, b2xpni, b2xpfe, b2xpfi, b2xppb,
     &  b2xpnm, b2tral, b2trql, b2treq, b2tfnb, b2tfhe, b2tfhi,
     &  b2npmo, b2npco, b2nppo, b2npht, b2nxdp, b2tfcc, b2srsm,
     &  b2tfch, b2xpen, b2xppz
*   ..initialisation
      data ncall/0/, xfm0/1.0_R8/, xfm1/1.0_R8/, xfm2/1.0_R8/,
     1  xfm3/1.0_R8/, potit/50/, potok/1.0e-2_R8/, pot_eq/1/,
     2  ncallout/-1/, guard_flows/2/, no_solve/0/,
     5  recalculate_contributions/1/, ramp_slow/0/, no_current/0/,       !dpc
     6  potitmin/0/, BoRiS/0.0_R8/, vis_per/1.0_R8/, na_min/1.0e4_R8/,
     7  na_new/1.0e14_R8/
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2news')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       nsmin=0
       nsmax=ns
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
       call ipgetr ('b2news_BoRiS', BoRiS)
       call ipgeti ('b2ux7p_style', style)                               !srv 29.05.02
       call ipgetr ('b2tfhe_vis_per', vis_per)                           !srv 06.07.00
       call ipgetr ('b2mndr_na_min', na_min)
       call ipgetr ('b2mndr_na_new', na_new)
       call xertst (na_min.gt.0.0_R8,
     &  'faulty internal parameter na_min')
       call xertst (na_new.gt.0.0_R8,
     &  'faulty internal parameter na_new')
       call xertst (na_new.gt.na_min,
     &  'na_new must be larger than na_min!')
       call ipgeti ('b2news_ramp_slow', ramp_slow)
       call ipgeti ('b2news_guard_flows', guard_flows)                   !xpb
       call ipgeti ('b2news_ncallout', ncallout)
       call ipgeti ('b2news_potit', potit)
       call ipgeti ('b2news_potitmin', potitmin)
       call xertst (0.lt.potit,
     &   'faulty internal parameter potit')
       call xertst (potitmin.lt.potit,
     &   'faulty internal parameter potitmin')
       call ipgetr ('b2news_potok', potok)
       call xertst (0.0_R8.lt.potok,
     &   'faulty internal parameter potok')
       call ipgeti ('b2news_poteq', pot_eq)
       call ipgeti ('b2tfhe_no_current', no_current)
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
       call ipgeti ('b2npmo_ion_vlct_restrict', ion_vlct_restrict)       !iyv 06.02.14
       call ipgeti ('b2news_no_solve', no_solve)
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
* SOLVEPO, SOLVEMT, SOLVEET arrays from b2.numerics.parameters.
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
       call ipgeti ('b2news_recalculate_contributions',                  !dpc
     1  recalculate_contributions)                                       !dpc
       call b2news_facdrift_fac_ExB_init(nx, ny, ramp_slow)
       call ipgeti ('b2tfnb_mdf', mdf_fnb)                               !srv 13.10.06 {
       call ipgeti ('b2tfhe_mdf', mdf_fhe)
       call ipgeti ('b2tfhi_mdf', mdf_fhi)                               !srv 13.10.06 }
       call xertst (mdf_fhe.eq.0 .eqv. mdf_fhi.eq.0,
     &   ' Both mdf_fhe.eq.0 and mdf_fhe.eq.0 should be true,'//
     &   ' or both should be false')
       call ipgeti ('b2news_use_b2npc9', use_b2npc9)                     !srv 06.11.08
       call ipgeti ('b2news_use_b2nph9', use_b2nph9)                     !srv 06.11.08
       call ipgeti ('b2news_iout', iout)
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
       call ipgeti ('b2sigp_style', ue_style)                            !srv 13.01.17
      endif
*   ..set n2, zero
      n2 = (nx+2)*(ny+2)
      zero(-1:nx,-1:ny) = 0.0_R8
*   ..test nx, ny, ns, ismain
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
*    ..test vol, hx, hy, qz
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
      rxf2 =(rxf*xfm2)**0.2_R8
      rxf3 = rxf*xfm3

      if (ramp_call .or. ncall.eq.0) then
        call b2news_facdrift_fac_ExB_inc(nx,ny,ramp_slow)
      endif
*   ..test facdrift, fac_ExB, fac_vis
      call xertst (0.0_R8.le.minval(facdrift) .and.
     &     1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
      call xertst (0.0_R8.le.minval(fac_ExB) .and.
     &     1.0_R8.ge.maxval(fac_ExB), 'faulty argument fac_ExB')
      call xertst (0.0_R8.le.minval(fac_vis) .and.
     &     1.0_R8.ge.maxval(fac_vis), 'faulty argument fac_vis')
* ..main computation
      ipass = 0                                                          !xpb
*  ..loop for check of flow behaviour                                    !xpb
      do while(ipass.eq.0 .or. wrong_flow)                               !xpb
*   ..compute ni, ne, ne2
        call b2xpni (nx, ny, ns, na, ni)
        call b2xpne (nx, ny, ns, rza, na, ne_ext, ne)
        call b2xpne (nx, ny, ns, rz2, na, ne2_ext, ne2)
*   ..compute ceqp
        call b2xpnm (nx, ny, ns, rz2, am, na, wrk0)
        call b2treq (nx, ny, vol, te, ti, ne, ni, wrk0, lnlam, ceqp)     !xpb !srv 01.07.09 20.09.11
*   ..compute flux limit coefficients
        call b2trql (nx, ny, ns, bb, gs, vol, hx, hy, qz, qc, pbs,
     &   na, ua, te, ti, po, ne, ni, ne2, chvemx, chvimx)                !srv 16.06.08
*   ..compute transport coefficients
        call b2tral (nx, ny, ns, nscx, nscxmax, iscx, ismain,            !srv 23.11.10
     1   bb, conn, vol, gs, hx, hy, hz, qz, qc, pbs, crx, cry,
     .   bzb, lnlam,                                                     !srv 13.10.06 01.07.09
     1   fch_p, na, ua, te, ti, po, ne, ni, ne2, chvemx, chvimx,         !srv 03.06.03
     2   cdna, cdpa, cddi, cvla, cvsa, chce, chve, chci,                 !srv 13.10.06
     3   chvi, csig, csigin, calf, cthe, cthi,
     4   cdnahz, cdpahz, cvlahz, cvmahz, cvsahz, cvsa_cl, cvsahz_cl,     !srv 09.01.01 15.12.05
     5   fllim0fhi, fllimvisc, csig_cl, calf_cl)
        if (ncall.eq.ncallout) then
          write(*,*) 'ncall,ipass',ncall,ipass
          call hbout('b2ne_npmo   ',nx,ny,ns,te,'te',ti,'ti',
     &     na,'na',ua,'ua',resmo,'resmo',corua,'corua')
        endif
*   ..re-compute source terms if necessary                               !xpb
        if((wrong_flow .and. ipass.ne.0) .or. ncall.eq.0) then           !xpb
*      ..compute electron drift velocities
          call b2tfed (nx, ny, facdrift, fac_ExB, vol, hx, hy, qz, bb,   !srv 05.06.18
     &     te, po, vedia, veecrb)
*     ..re-compute source linearisation                                  !xpb
          call b2sral (nx, ny, ns, nxtl, nxtr,
     &     nscx, nscxmax, nsmin, nsmax, iscx, ismain, ismain0,
     &     dtim, BoRiS, facdrift, fac_ExB, fac_vis,                      !xpb
     &     vol, crx, hx, hy, hz, qz, qc, gs, pbs, pbshz, bb, lnlam, na,  !srv 09.01.01 02.04.18
     &     ua, uadia, vedia, vadia, wadia, veecrb, vaecrb, ve, wedia,
     &     te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,               !srv 28.10.02
     &     fna, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &     fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,     !srv 13.10.06 25.12.06 23.10.09
     &     fhe, fhi, fnaPSch, fhePSch, fhiPSch, fch, fchdia, fchin,
     &     fch_p, fchvispar, fchvisper, fchvisq, fchinert,               !srv 21.10.02
     &     fchanml, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &     cdna, cdpa, cvsa_cl, cvla, chce, chve, chci, chvi, calf,      !xpb
     &     rlsa, rlra, rlqa, rlcx, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &     rza, rz2, rpt, rpi, sna, smo, smq, she, shi, sch, sne,
     &     wrong_flow, .false.)
          if(wrong_flow) write(*,*)
     &         'b2news [1]: wrong_flow returned from b2sral'
*     ..smooth out sources
          call b2srsm (nx, ny, ns, dtim,
     &     vol, na, ua, te, ti, ne, ni,
     &     sna, smo, she, shi, sne, .false.)
*     ..include contributions due to timestep                            !xpb
          call b2srdt (nx, ny, ns, guard_flows, dtim, BoRiS, hz,         !xpb !srv 09.01.01 01.07.08
     &     vol, na0, ua0, te0, ti0, ne0, ni0, kinrgy0,
     &     na, ua, te, ti, ne, ni, kinrgy,                               !srv 16.07.10
     &     sna, smo, she, shi, sch, sne, .false.)                        !xpb
*     ..stabilise the source coefficients                                !xpb
          call b2srst (nx, ny, ns, guard_flows, na, ua, te, ti, po, ne,  !xpb !srv 01.07.08
     &     ni, sna, smo, she, shi, sch, sne)                             !xpb
        endif                                                            !xpb
        call b2tanml (nx, ny, ns, vol, hx, hy, csig_an, po, fchanml)     !srv 20.09.06 }
        call b2tvspa (nx, ny, ns, vol, hx, hy, qz, bb, ua, cvsa_cl,      !srv 20.09.06 {
     &                       fac_vis, fchvispar)
        call b2tiner (nx, ny, ns, vol,
     &                       gs, hx, hy, bb, na, ua, facdrift, fchinert)
        call b2tvspr (nx, ny, ns, fac_vis, vol, hx, hy, qz, bb,         !srv 10.10.17
     &          po, ti, na, avisper, bvisper, cvisper, fchvisper)
        call b2tvsq (nx, ny, ns, fac_vis, vol, hx, hy, qz, bb, ne2,     !srv 10.10.17
     &               lnlam, te, ti, ni, chci, fchvisq)

        if(ion_vlct_restrict .ne. 0 ) then                               !iyv 06.02.14
          call b2xprz (nx, ny, ns, mp, am, na, rz)
          call b2xppz (nx, ny, ns, ne, na, te, ti, pz)
          wrk0(:,:)=sqrt(pz(:,:)/rz(:,:))*abs(bb(:,:,0)/bb(:,:,3))
          call intfaceh(nx,ny,leftix,leftiy,vol,wrk0,cssb(-1,-1,0))
          call intfacev(nx,ny,bottomix,bottomiy,vol,wrk0,cssb(-1,-1,1))
        else
          cssb = 0.0_R8
        endif                                                            !iyv 06.02.14
#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
!! FR possible loop (1) for use with OpenMP
!$OMP PARALLEL DO DEFAULT(none) FIRSTPRIVATE(flob, conb, pb)
!$OMP& SHARED(nx, ny, rza, na, te, ti, ns, nsmin, nsmax, ismain)
!$OMP& SHARED(guard_flows, facdrift, fac_ExB, am, cssb)
!$OMP& SHARED(vol, pbs, pbshz, ua, uadia, vaecrb, vadia, wadia)
!$OMP& SHARED(ne, bb, hx, hy, hz, qz, gs, qc, po, cvla, cdna, cdpa)
!$OMP& SHARED(cddi, cvlahz, cvmahz, cdnahz, cdpahz, fchvispar)
!$OMP& SHARED(fchvisper, fchinert, fchanml, fchvisq)
!$OMP& SHARED(fna, fna_nodrift, fna_fcor, fna_he, fnaPSch, fna_mdf)
!$OMP& SHARED(fna_32, fna_52, fna_32nd, fna_52nd)
!$OMP& SHARED(kinrgy, fllim0fna) 
!        do is = 0, ns-1
        do is = nsmin, nsmax-1
*    ..compute partial pressure
         call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)
*    ..compute particle flux
         call b2tfnb (nx, ny, is, ismain,                                !srv 08.09.01
     1    guard_flows, facdrift, fac_ExB, am(is)*mp, vol, cssb,          !iyv 06.02.14
     1    rza(-1,-1,is), pbs, pbshz, na(-1,-1,is), pb, ua(-1,-1,is),
     1    uadia(-1,-1,0,is), vaecrb(-1,-1,0,is), vadia(-1,-1,0,is),      !xpb
     &    wadia(-1,-1,0,is), ne, te, ti, bb, hx, hy, hz, qz, gs, qc,     !srv 09.01.01
     &    po, cvla(-1,-1,0,is), cdna(-1,-1,0,is), cdpa(-1,-1,0,is),
     &    cddi(-1,-1,0,is),                                              !srv 13.10.06
     &    cvlahz(-1,-1,0,is), cvmahz(-1,-1,0,is),  cdnahz(-1,-1,0,is),   !srv 09.01.01 15.12.05
     &    cdpahz(-1,-1,0,is),
     &    fchvispar, fchvisper, fchinert, fchanml, fchvisq,              !srv 20.09.06
     &    fna(-1,-1,0,is), fna_nodrift(-1,-1,0,is),
     &    fna_fcor(-1,-1,0,is), fna_he(-1,-1,0,is),                      !srv 01.10.99 21.01.00
     &    fnaPSch(-1,-1,0,is), fna_mdf(-1,-1,0,is),                      !srv 13.10.06
     1    flob, conb, fna_32(-1,-1,0,is), fna_52(-1,-1,0,is),
     1    fna_32nd(-1,-1,0,is), fna_52nd(-1,-1,0,is),
     2    kinrgy(-1,-1,is), fllim0fna(-1,-1,0,is), .true.)               !djm Jan2017
        enddo
!$OMP END PARALLEL DO
!! FR end of loop (1)
#ifdef WITH_FTIMINGS
        call timer%stop("parallel")
#endif
*      ..compute electron drift velocities
        call b2tfed (nx, ny, facdrift, fac_ExB, vol, hx, hy, qz, bb,
     &    te, po, vedia, veecrb)
*     ..re-compute source linearisation                                  !xpb
        call b2sral (nx, ny, ns, nxtl, nxtr,
     &    nscx, nscxmax, nsmin, nsmax, iscx, ismain, ismain0,
     &    dtim, BoRiS, facdrift, fac_ExB, fac_vis,                       !xpb
     &    vol, crx, hx, hy, hz, qz, qc, gs, pbs, pbshz, bb, lnlam,       !srv 09.01.01 02.04.18
     &    na, ua, uadia, vedia, vadia, wadia, veecrb, vaecrb, ve, wedia,
     &    te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,                !srv 28.10.02
     &    fna, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &    fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,      !srv 13.10.06 25.12.06 23.10.09
     &    fhe, fhi, fnaPSch, fhePSch, fhiPSch, fch,
     &    fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert, !srv 21.10.02
     &    fchanml, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &    cdna, cdpa, cvsa_cl, cvla, chce, chve, chci, chvi, calf,       !xpb
     &    rlsa, rlra, rlqa, rlcx, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &    rza, rz2, rpt, rpi, sna, smo, smq, she, shi, sch, sne,
     &    wrong_flow, .false.)
        if(wrong_flow) write(*,*)
     &         'b2news [2]: wrong_flow returned from b2sral'
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

        if (ue_style.eq.2 .or. iout.ne.0)
*   ..compute ue, electron parallel velocity
     >   call b2xpve (nx, ny, ns, qe, rza, rz2, alfx_c, sigx_c, po,      !srv 05.07.17
     &              ne, te, na, ua, ue_ext, ne_ext, ue)
        if (ue_style.ne.2 .or. iout.ne.0) then
          do iy=-1,ny
            do ix=-1,nx
              if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
                wrk0(ix,iy)=-(fch_p(ix,iy,0)/pbs(ix,iy,0)+
     &           fch_p(rightix(ix,iy),rightiy(ix,iy),0)/
     &           pbs(rightix(ix,iy),rightiy(ix,iy),0))/qe/2.0_R8
              elseif(leftix(ix,iy).eq.-2) then
                wrk0(ix,iy)=-fch_p(rightix(ix,iy),rightiy(ix,iy),0)/
     &           pbs(rightix(ix,iy),rightiy(ix,iy),0)/qe
              elseif(rightix(ix,iy).eq.nx+1) then
                wrk0(ix,iy)=-fch_p(ix,iy,0)/pbs(ix,iy,0)/qe
              endif
!             wrk0(ix,iy)=wrk0(ix,iy)*prl_cur                            !srv 03.04.15
            enddo
          enddo
          do is=0,ns-1
            wrk0(:,:)=wrk0(:,:)+rza(:,:,is)*ua(:,:,is)*na(:,:,is)
          enddo
          wrk0(:,:)=wrk0(:,:)+ue_ext(:,:)*ne_ext(:,:)
          wrk0(:,:)=wrk0(:,:)/ne(:,:)                                    !srv 05.03.11
        endif
        if(iout.ne.0) then                                               !srv 05.07.17
          call b2xzef (nx, ny, ns, rz2, na, ne, zeff)
          call b2xpic (nx, ny, ns, qe, rza, bb, na, ua, zeff, fch_pi_c)
          call b2xpif (nx, ny, hx, hy, vol, fch_pi_c, fch_pi_f)
          call my_out(70,nx,ny,fch_pi_c,'b2news_fch_pi_c')
          call my_out(70,nx,ny,fch_pi_f,'b2news_fch_pi_f')
          call my_out(70,nx,ny,ue,'b2news_ve')
          call my_out(70,nx,ny,wrk0,'b2news_ue')
        endif                                                            !srv 05.07.17
        if (ue_style.ne.2) ue=wrk0

cxpb
*   ..iterate on parallel momentum balance
        call b2npmo (nx, ny, ns, ismain, ncall, solving(1), rxf1,        !srv 13.01.17
     &    qz, vol, hx, hy, hz, bb, pbshz, lnlam,                         !srv 01.07.09
     &    na, te, ti, po, ue, fna_fcor,                                  !srv 05.07.17
     &    cvsahz, cvsa, cvsa_cl, cvsahz_cl, kt_neo,                      !srv 02.07.08
     &    smo, smq, cthe, cthi, chci, ua, smof, wadia, vaecrb, resmo,
     &    resmt, corua, corut, pcca, pccm, ier0)
        call xertst (ier0.eq.0, 'error return from b2npmo')
*   ..set nep
*     (= ne before continuity update)
        call b2scopy (n2, ne, 1, nep, 1)
*   ..recompute currents
        call b2tanml (nx, ny, ns, vol, hx, hy, csig_an, po, fchanml)     !srv 20.09.06 }
        call b2tvspa (nx, ny, ns, vol, hx, hy, qz, bb, ua, cvsa_cl,      !srv 20.09.06 {
     &                       fac_vis, fchvispar)
        call b2tiner (nx, ny, ns, vol,
     &                       gs, hx, hy, bb, na, ua, facdrift, fchinert)
        call b2tvspr (nx, ny, ns, fac_vis, vol, hx, hy,  qz, bb,        !srv 10.10.17
     &       po, ti, na, avisper, bvisper, cvisper, fchvisper)
        call b2tvsq (nx, ny, ns, fac_vis, vol, hx, hy, qz, bb,          !srv 10.10.17
     &               ne2, lnlam, te, ti, ni, chci, fchvisq)

*     ..iterate on particle balance
#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
!!FR Possible loop (2) for use with OpenMP
!$OMP PARALLEL DO DEFAULT(none)
!$OMP& PRIVATE(ier0, solvedum, solvedm2, mach10, ix, iy, new_matrix)
!$OMP& PRIVATE(ireg)
!$OMP& FIRSTPRIVATE(pb, flob, conb)
#ifndef GFORTRAN
#ifndef PGF90
!$OMP& COPYIN(last_solve_5, last_solve_9) 
!$OMP& COPYIN(last_solve_7, last_solve_11)
#else
!$OMP& PRIVATE(last_solve_5, last_solve_7)
!$OMP& PRIVATE(last_solve_9, last_solve_11)
#endif
#else
!$OMP& PRIVATE(last_solve_5, last_solve_7)
!$OMP& PRIVATE(last_solve_9, last_solve_11)
#endif
!$OMP& SHARED(nsmin, nsmax, na_min, ismain, mdf_fnb)
!$OMP& SHARED(ns, nx, ny, region, leftix, rightix, bottomiy)
!$OMP& SHARED(topiy, na, rza, zn, zamax, te, ti, am, ua, ne)
!$OMP& SHARED(is_neutral, uadia, ncall, ncallout, resmo, corua)
!$OMP& SHARED(guard_flows, facdrift, fac_ExB, fac_vis, vol, pbs, pbshz)
!$OMP& SHARED(cssb, vaecrb, vadia)
!$OMP& SHARED(wadia, bb, hx, hy, hz, qz, gs, qc, po, cvla, cdna, cdpa)
!$OMP& SHARED(cddi, cvlahz, cvmahz, cdnahz, cdpahz, fchvispar)
!$OMP& SHARED(fchvisper, fchinert, fchanml, fchvisq)
!$OMP& SHARED(fna, fna_nodrift, fna_fcor, fna_he, fnaPSch, fna_mdf)
!$OMP& SHARED(fna_32, fna_52, fna_32nd, fna_52nd)
!$OMP& SHARED(kinrgy, fllim0fna, solving, nnreg, solveco, solvemo)
!$OMP& SHARED(rxf0, sna, pcca, resco, corpa, b2stbr_sna, b2stbm_sna) 
!$OMP& SHARED(ext_sna, use_b2npc9)
!        do is = 0, ns-1
        do is = nsmin, nsmax-1
*    ..DPC
          do iy=-1,ny
            do ix=-1,nx
              na(ix,iy,is)=max(na(ix,iy,is),na_min)
              if(region(ix,iy,0).ne.0 .and.
     &         leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
                mach10=10.0_R8*
     &           sqrt((rza(ix,iy,is)*te(ix,iy)+ti(ix,iy))/(am(is)*mp))
                if(ua(ix,iy,is).gt.mach10.and.(.not.is_neutral(is)))
     >           then
                  if(na(ix,iy,is)/ne(ix,iy).gt.1.0e-6_R8 .or.
     &               abs(b2stbr_sna(ix,iy,is)).gt.1.0e8_R8 .or.
     &               abs(b2stbm_sna(ix,iy,is)).gt.1.0e8_R8 .or.
     &               abs(ext_sna(ix,iy,is)).gt.1.0e8_R8) then
                    write(*,*)
     1               'Velocity limit active for ix,iy,is = ',ix,iy,is
                    write(*,'(a,1p,7g14.7)')
     .               'ua,uadia(x&y),te,ti,mach10,na/ne ',
     1                ua(ix,iy,is),(uadia(ix,iy,i,is),i=0,1),
     2                te(ix,iy)/ev,ti(ix,iy)/ev,mach10,
     3                na(ix,iy,is)/ne(ix,iy)
                    write(*,'(a,1p,3g14.7)') 'sources(r,m,ext): ',
     1               b2stbr_sna(ix,iy,is),b2stbm_sna(ix,iy,is),
     2               ext_sna(ix,iy,is)
                  endif
                  ua(ix,iy,is)=mach10
                elseif(ua(ix,iy,is).lt.-mach10.and.
     &                 (.not.is_neutral(is))) then
                  if(na(ix,iy,is)/ne(ix,iy).gt.1.0e-6_R8 .or.
     &               abs(b2stbr_sna(ix,iy,is)).gt.1.0e8_R8 .or.
     &               abs(b2stbm_sna(ix,iy,is)).gt.1.0e8_R8 .or.
     &               abs(ext_sna(ix,iy,is)).gt.1.0e8_R8) then
                    write(*,*)
     1               'Velocity limit active for ix,iy,is = ',ix,iy,is
                    write(*,'(a,1p,7g14.7)')
     .               'ua,uadia(x&y),te,ti,mach10,na/ne ',
     1                ua(ix,iy,is),(uadia(ix,iy,i,is),i=0,1),
     2                te(ix,iy)/ev,ti(ix,iy)/ev,mach10,
     3                na(ix,iy,is)/ne(ix,iy)
                    write(*,'(a,1p,3g14.7)') 'sources(r,m,ext): ',
     1                b2stbr_sna(ix,iy,is),b2stbm_sna(ix,iy,is),
     2                ext_sna(ix,iy,is)
                  endif
                  ua(ix,iy,is)=-mach10
                endif
              endif
            enddo
          enddo
          if (ncall.eq.ncallout) then
            write(*,*) 'ncall,is',ncall,is
            call hbout('b2ne_xppb   ',nx,ny,ns,te,'te',ti,'ti',
     &       na,'na',ua,'ua',resmo,'resmo',corua,'corua')
          endif
*    ..compute partial pressure
          call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)
*    ..compute particle flux
          call b2tfnb (nx, ny, is, ismain,                               !srv 08.09.01
     1     guard_flows, facdrift, fac_ExB, am(is)*mp, vol, cssb,         !xpb !iyv 06.02.14
     1     rza(-1,-1,is), pbs, pbshz, na(-1,-1,is), pb, ua(-1,-1,is),
     1     uadia(-1,-1,0,is), vaecrb(-1,-1,0,is), vadia(-1,-1,0,is),     !xpb
     &     wadia(-1,-1,0,is), ne, te, ti, bb, hx, hy, hz, qz, gs, qc,    !srv 09.01.01
     &     po, cvla(-1,-1,0,is), cdna(-1,-1,0,is), cdpa(-1,-1,0,is),
     &     cddi(-1,-1,0,is),                                             !srv 13.10.06
     &     cvlahz(-1,-1,0,is), cvmahz(-1,-1,0,is), cdnahz(-1,-1,0,is),   !srv 09.01.01 15.12.05
     &     cdpahz(-1,-1,0,is),
     &     fchvispar, fchvisper, fchinert, fchanml, fchvisq,             !srv 20.09.06
     1     fna(-1,-1,0,is), fna_nodrift(-1,-1,0,is),                     !xpb
     &     fna_fcor(-1,-1,0,is), fna_he(-1,-1,0,is),                     !srv 01.10.99 21.01.00
     &     fnaPSch(-1,-1,0,is), fna_mdf(-1,-1,0,is),                     !srv 13.10.06
     1     flob, conb, fna_32(-1,-1,0,is), fna_52(-1,-1,0,is),
     1     fna_32nd(-1,-1,0,is), fna_52nd(-1,-1,0,is),
     1     kinrgy(-1,-1,is), fllim0fna(-1,-1,0,is), .true.)              !djm Jan2017
          if (ncall.eq.ncallout) then
            write(*,*) 'ncall,is',ncall,is
            call hbout('b2ne_npco   ',nx,ny,ns,te,'te',ti,'ti',
     &       na,'na',ua,'ua',na,'na',ua,'ua')
          endif
!        if (solving(2)) then                           !srv 22.05.18
!        if (all(solveco(is,0:nnreg(0)))) then
*    ..perform iteration
            solvedum(0:nnreg(0))=solveco(is,0:nnreg(0))
            solvedm2(0:nnreg(0))=solvemo(is,0:nnreg(0))
            new_matrix = .false.
            if(mdf_fnb.eq.0) then                                        !srv 13.10.06
              if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then        !srv 02.04.18
                do ireg = 0, nnreg(0)
                  new_matrix = new_matrix.or.
     &            solveco(is,ireg).neqv.last_solve_7(ireg)
                end do
                if(new_matrix) call restart_ma28_for_7
                call b2npc7 (nx, ny, fac_vis, guard_flows,
     &           rza(-1,-1,is), nnreg(0), solving(2), solvedum,          !srv 22.05.18 
     &           solving(1), solvedm2, ncall, is, rxf0, vol, pbs, flob,  !srv 22.01.02 01.07.05
     &           conb, sna(-1,-1,0,is), pcca(-1,-1,0,is), ti, pb,
     &           na(-1,-1,is), ua(-1,-1,is), fna(-1,-1,0,is),
     &           resco(-1,-1,is), corpa(-1,-1,is), ier0)
                 call xertst (ier0.eq.0, 'error return from b2npc7')
                last_solve_7(0:nnreg(0)) = solveco(is,0:nnreg(0))
              else                                                       !srv 02.04.18
                do ireg = 0, nnreg(0)
                  new_matrix = new_matrix.or.
     &            solveco(is,ireg).neqv.last_solve_5(ireg)
                end do
                if(new_matrix) call restart_ma28
                call b2npco (nx, ny, guard_flows,                        !srv 05.07.08
     &           nnreg(0), solving(2), solvedum, solving(1), solvedm2,   !srv 22.05.18
     &           ncall, is, rxf0, vol, pbs, flob, conb,
     &           sna(-1,-1,0,is), pcca(-1,-1,0,is), pb, na(-1,-1,is),
     &           ua(-1,-1,is), fna(-1,-1,0,is), resco(-1,-1,is),
     &           corpa(-1,-1,is), ier0)
                call xertst (ier0.eq.0, 'error return from b2npco')
                last_solve_5(0:nnreg(0)) = solveco(is,0:nnreg(0))
              endif                                                      !srv 02.04.18
            else                                                         !srv 13.10.06 {
              if(.not.is_neutral(is) .and. use_b2npc9.eq.1) then         !srv 06.11.08
                if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then      !srv 02.04.18
                  do ireg = 0, nnreg(0)
                    new_matrix = new_matrix.or.
     &               solveco(is,ireg).neqv.last_solve_11(ireg)
                  end do
                  if(new_matrix) call restart_ma28_for_11
                  call b2npc11 (nx, ny, fac_vis, guard_flows,
     &             rza(-1,-1,is), nnreg(0), solving(2), solvedum,        !srv 22.05.18
     &             solving(1), solvedm2, ncall, is, rxf0, vol,
     &             pbs, flob, conb, sna(-1,-1,0,is), pcca(-1,-1,0,is),
     &             po, ti, pb, na(-1,-1,is), ua(-1,-1,is),
     &             fna_mdf(-1,-1,0,is), resco(-1,-1,is),
     &             corpa(-1,-1,is), ier0)
                  call xertst (ier0.eq.0, 'error return from b2npc11')
                  last_solve_11(0:nnreg(0)) = solveco(is,0:nnreg(0))
                else                                                     !srv 02.04.18
                  do ireg = 0, nnreg(0)
                    new_matrix = new_matrix.or.
     &               solveco(is,ireg).neqv.last_solve_9(ireg)
                  end do
                  if(new_matrix) call restart_ma28_for_9
                  call b2npc9 (nx, ny, guard_flows, nnreg(0), 
     &             solving(2), solvedum, solving(1), solvedm2, ncall,    !srv 22.05.18
     &             is, rxf0, vol, pbs,
     &             flob, conb, sna(-1,-1,0,is), pcca(-1,-1,0,is),
     &             pb, na(-1,-1,is), ua(-1,-1,is), fna_mdf(-1,-1,0,is),
     &             resco(-1,-1,is), corpa(-1,-1,is), ier0)
                  call xertst (ier0.eq.0, 'error return from b2npc9')
                  last_solve_9(0:nnreg(0)) = solveco(is,0:nnreg(0))
                endif
              else
                do ireg = 0, nnreg(0)
                  new_matrix = new_matrix.or.
     &             solveco(is,ireg).neqv.last_solve_5(ireg)
                end do
                if(new_matrix) call restart_ma28
                call b2npco (nx, ny, guard_flows, nnreg(0), solving(2),  !srv 22.05.18
     &           solvedum, solving(1), solvedm2, ncall, is, rxf0, vol, 
     &           pbs, flob, conb, sna(-1,-1,0,is), pcca(-1,-1,0,is),
     &           pb, na(-1,-1,is), ua(-1,-1,is), fna_mdf(-1,-1,0,is),
     &           resco(-1,-1,is), corpa(-1,-1,is), ier0)
                call xertst (ier0.eq.0, 'error return from b2npco')
                last_solve_5(0:nnreg(0)) = solveco(is,0:nnreg(0))
              endif
            endif                                                        !srv 13.10.06 }
!*    ..re-compute partial pressure                                      !srv 01.10.99
!            call b2xppb (nx, ny,
!     &                   rza(-1,-1,is), na(-1,-1,is), te, ti, pb) !srv 01.10.99  srv 12.11.08
*    ..re-compute particle flux
            call b2tfnb (nx, ny, is, ismain,                             !srv 08.09.01
     1       guard_flows, facdrift, fac_ExB, am(is)*mp, vol, cssb,       !xpb !iyv 06.02.14
     1       rza(-1,-1,is), pbs, pbshz, na(-1,-1,is), pb, ua(-1,-1,is),
     1       uadia(-1,-1,0,is), vaecrb(-1,-1,0,is), vadia(-1,-1,0,is),   !xpb
     &       wadia(-1,-1,0,is), ne, te, ti, bb, hx, hy, hz, qz, gs, qc,  !srv 09.01.01
     &       po, cvla(-1,-1,0,is), cdna(-1,-1,0,is), cdpa(-1,-1,0,is),   !xpb
     &       cddi(-1,-1,0,is),                                           !srv 13.10.06
     &       cvlahz(-1,-1,0,is), cvmahz(-1,-1,0,is), cdnahz(-1,-1,0,is), !srv 09.01.01 15.12.05
     &       cdpahz(-1,-1,0,is),
     &       fchvispar, fchvisper, fchinert, fchanml, fchvisq,           !srv 20.09.06
     1       fna(-1,-1,0,is), fna_nodrift(-1,-1,0,is),                   !xpb
     &       fna_fcor(-1,-1,0,is), fna_he(-1,-1,0,is),                   !srv 01.10.99 21.01.00
     &       fnaPSch(-1,-1,0,is),fna_mdf(-1,-1,0,is),                    !srv 13.10.06
     1       flob, conb, fna_32(-1,-1,0,is), fna_52(-1,-1,0,is),
     1       fna_32nd(-1,-1,0,is), fna_52nd(-1,-1,0,is),
     1       kinrgy(-1,-1,is), fllim0fna(-1,-1,0,is), .false.)           !djm Jan2017
!          endif                                                         !srv 22.05.18
        enddo
!$OMP END PARALLEL DO
!! FR End of loop (2)
#ifdef WITH_FTIMINGS
      call timer%stop("parallel")
#endif
        if(wrong_flow .or. ipass.ne.0) then                              !xpb
           call print_wrong_flows(
     &          'Before potential loop: fna(ions) at plates')
        endif
*      ..compute electron drift velocities
        call b2tfed (nx, ny, facdrift, fac_ExB, vol, hx, hy, qz, bb,
     &    te, po, vedia, veecrb)
*   ..re-compute ne, ni
        call b2xpne (nx, ny, ns, rza, na, ne_ext, ne)
        call nancheck(nx, ny, ne, 'b2news(ne): after call to b2xpne')
        call b2xpni (nx, ny, ns, na, ni)
*   ..check if we are solving for the potential
        if (pot_eq.eq.0) then
          po=0.0_R8
        else if (pot_eq.eq.2) then
          po(:,:)=3.1_R8*te(:,:)/qe
        else
*   ..correct po
          if (solving(3)) then
            call b2nxdp (nx, ny, nnreg(0), solvepo,
     .       qe, te, nep, ne, wrk0)
            call b2saxpy (n2, 1.0_R8, wrk0, 1, po, 1)
          endif
          if(.not.wrong_flow) then                                       !xpb !srv 22.05.18
*   ..iterate several times on the potential equation                    !xpb
            i=0                                                          !xpb
            t0=1.0e40_R8                                                 !xpb
            if(ncall.eq.0) then
              if(vis_per.eq.0.0_R8) then
                call restart_ma28
              else
                if (style.eq.2) call restart_ma28_for_7                  !srv 11.01.00 29.05.02
                if (style.eq.3) call restart_sdrv_for_7                  !srv 29.05.02
              endif
            endif
            do while(i.lt.potit.and.(t0.gt.potok.or.i.le.potitmin).and.
     1       .not.wrong_flow)                                            !xpb
*    ..compute current
              call b2xpfe (nx, ny, ns, ns_ext, qe,
     &         rza, za_ext, fna, fa_ext, fch, fne)
              call b2xpfe (nx, ny, ns, ns_ext, qe,
     1         rza, za_ext, fna_32nd, fa_ext, fch_32, fne_32)
              call b2xpfe (nx, ny, ns, 0, qe,
     1         rza, zero, fna_52nd, zero, fch_52, fne_52)
              call b2xpfi (nx, ny, ns, fna, fni)                         !xpb
              call b2xpfi (nx, ny, ns, fna_32nd, fni_32)                 !xpb
              call b2xpfi (nx, ny, ns, fna_52nd, fni_52)                 !xpb
              call b2xppz (nx, ny, ns, ne, na, te, ti, pz)               !srv 18.01.02
              if (ncall.eq.ncallout)
     >         write(*,*) 'b2news [a][32|52] ',
     1         maxval(abs(fna-fna_32-fna_52)), maxval(abs(fna))
              call b2tfed (nx, ny, facdrift, fac_ExB,
     &         vol, hx, hy, qz, bb,
     &         te, po, vedia, veecrb)
              call b2tfch (nx, ny, ns, ismain, ismain0,
     1         guard_flows, facdrift, fac_ExB, fac_vis, BoRiS,           !xpb
     1         vol, gs, hx, hy, qz, qc, bb, te, po, ne, ne2, lnlam, pz,  !xpb
     1         ti, ni(-1,-1,1), na, fna, uadia, vedia, wadia, vaecrb,    !xpb
     &         veecrb, ua, fne_32, fne_52,
     &         csig, csigin, chce, chve, chci,
     &         calf, chvemx, cvsa_cl, fch, conc, fchstoch, fchdia,       !srv 17.12.13 16.10.17
     &         fch_p, fchvispar, fchvisq, fchin, fchinert,
     &         fch_32, fch_52, avisper, bvisper, cvisper, fchvisper,     !srv 25.11.08
     1         csig_cl, calf_cl)
              call b2tfhe (nx, ny, ns, guard_flows, facdrift,
     1         BoRiS, vol, hx, hy, qz, te, po, ne, na, vedia, veecrb,    !xpb  !srv 13.01.17
     &         fne_32, fne_52, fch_p, csig, chce, chve, calf, chvemx,    !xpb
     &         cddi, fhe, fhe_mdf, fhePSch, fhet, floe, cone, fllime,    !xpb  !srv 13.10.06
     &         fchstoch, .false.)                                        !djm Jan2017 !srv 17.12.13
              call b2tfhi (nx, ny, ns, guard_flows, facdrift, BoRiS,     !xpb
     1         vol, hx, hy, gs, qz, qc, bb, ti, ni, na,
     2         fna, fni_32, fni_52, kinrgy,
     &         chci, cddi, chvi, chvimx, cvsa, vadia, vaecrb,
     &         fchvispar, fchvisper, fchvisq, fchinert, fchanml,         !lk 30.05.11
     3         fhi, fhi_mdf, fhiPSch, fhit, floi, coni, fllimi,          !srv 13.10.06
     &         .false.)                                                  !djm Jan2017
              call b2xpen (nx, ny, ns,
     1         fna, fne, fch, fhe, fhi, kinrgy, rpt, te, ti, po, BoRiS,
     2         fhm, fhp, fhj, fht)
              if(recalculate_contributions.ne.0) then
*     ..update sources                                                   !xpb
                call b2sral (nx, ny, ns, nxtl, nxtr,
     &           nscx, nscxmax, nsmin, nsmax, iscx, ismain, ismain0,
     &           dtim, BoRiS, facdrift, fac_ExB, fac_vis,                !xpb
     &           vol, crx, hx, hy, hz, qz, qc, gs, pbs, pbshz, bb,
     &           lnlam, na, ua,                                          !srv 09.01.01
     &           uadia, vedia, vadia, wadia, veecrb, vaecrb, ve, wedia,
     &           te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,         !srv 28.10.02
     &           fna, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &           fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift,       !srv 13.10.06 25.12.06 23.10.09
     &           fna_he, fhe, fhi, fnaPSch, fhePSch, fhiPSch, fch,
     &           fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq,    !srv 21.10.02
     &           fchinert, fchanml, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &           cdna, cdpa, cvsa_cl, cvla, chce,
     &           chve, chci, chvi, calf,
     &           rlsa, rlra, rlqa, rlcx, rlrd,
     &           rlbr, rlza, rlz2, rlpt, rlpi,
     &           rza, rz2, rpt, rpi,
     &           sna0, smo0, smq0, she0, shi0, sch0, sne0,
     &           wrong_flow, .false.)
                if(wrong_flow) write(*,*)
     &           'b2news [3]: wrong_flow returned from b2sral'
*     ..smooth out sources
                call b2srsm (nx, ny, ns, dtim,
     &           vol, na, ua, te, ti, ne, ni,
     &           sna0, smo0, she0, shi0, sne0, .false.)
*     ..include contributions due to timestep                            !xpb
                call b2srdt (nx, ny, ns, guard_flows, dtim, BoRiS, hz,   !xpb !srv 09.01.01 01.07.08
     &           vol, na0, ua0, te0, ti0, ne0, ni0, kinrgy0,
     &           na, ua, te, ti, ne, ni, kinrgy,                         !srv 16.07.10
     &           sna0, smo0, she0, shi0, sch0, sne0, .false.)            !xpb
*     ..stabilise the source coefficients                                !xpb
                call b2srst (nx, ny, ns, guard_flows, na, ua, te, ti,    !xpb !srv 01.07.08
     &           po, ne, ni, sna0, smo0, she0, shi0, sch0, sne0)         !xpb
                if(wrong_flow) then                                      !xpb
                  sna=sna0                                               !xpb
                  sne=sne0                                               !xpb
                  smo=smo0                                               !xpb
                  smq=smq0                                               !xpb
                  she=she0                                               !xpb
                  shi=shi0                                               !xpb
                  sch=sch0                                               !xpb
                endif                                                    !xpb
              endif  ! recalculate_contributions                         !dpc
*    ..perform iteration
              if(.not.wrong_flow) then                                   !xpb
                if(vis_per.eq.0.0_R8) then
                  new_matrix = .false.
                  do ireg = 0, nnreg(0)
                    new_matrix = new_matrix.or.
     &               solvepo(ireg).neqv.last_solve_5(ireg)
                  end do
                  if (new_matrix.and.ncall.ne.0) call restart_ma28
                  call b2nppo (nx, ny, nnreg(0), solving(3), solvepo,    !srv 22.05.18
     &             ncall, rxf2, ne, te, conc, sch, po, fch, wrk0, wrk1,
     &             ier0)
                  last_solve_5(0:nnreg(0)) = solvepo(0:nnreg(0))
                else
                  call b2npp7 (nx, ny, nnreg(0), solving(3), solvepo,    !srv 22.05.18
     &             ncall, rxf2, ne, te, conc,
     1             avisper, bvisper, cvisper,
     &             sch, po, fch, wrk0, wrk1, ier0)
                endif
                t0 = 0.0_R8
                tmax = -1.0_R8                                           !xpb
                ixmax = -2                                               !xpb
                iymax = -2                                               !xpb
                do iy = -1, ny
cxpb  Do not count residuals from dead or unsolved cells
                  do ix = -1, nx
                    ireg = region(ix,iy,0)
                    if(ireg.ne.0 .and. solvepo(ireg)) then               !xpb
                      t0 = t0 + abs(wrk0(ix,iy))
                      if (abs(wrk0(ix,iy)).gt.tmax) then                 !xpb
                        tmax = abs(wrk0(ix,iy))                          !xpb
                        ixmax = ix                                       !xpb
                        iymax = iy                                       !xpb
                      endif                                              !xpb
                    endif                                                !xpb
                  enddo
                enddo
                write(*,'(a,1p,2g14.7,2i4)')
     &           't0, tmax, ixmax, iymax = ', t0, tmax, ixmax, iymax
                if (i.eq.0) then
                  respo = wrk0
                  corpo = wrk1
                endif
                write (hlp_frm,'(a,i1,a)') '(i',len_of_digits(i),')'
                write (charnc(1:len_of_digits(i)),hlp_frm) i
c       call my_out(70,nx,ny,wrk0,'respo'//charnc(1:len_of_digits(i)))
c       call my_out(70,nx,ny,wrk1,'corpo'//charnc(1:len_of_digits(i)))
                if (ncall.eq.ncallout) then
                  write(*,*) 'ncall,i',ncall,i
                  call hbout('b2ne+nppo   ',nx,ny,ns,te,'te',ti,'ti',
     &             na,'na',ua,'ua',resmo,'resmo',corua,'corua')
                endif
                if(ion_vlct_restrict .ne. 0 ) then                       !iyv 06.02.14
                  call b2xprz (nx, ny, ns, mp, am, na, rz)
                  call b2xppz (nx, ny, ns, ne, na, te, ti, pz)
                  wrk0(:,:)=
     &                 sqrt(pz(:,:)/rz(:,:))*abs(bb(:,:,0)/bb(:,:,3))
                  call intfaceh(nx,ny,leftix,leftiy,vol,wrk0,
     &                        cssb(-1,-1,0))
                  call intfacev(nx,ny,bottomix,bottomiy,vol,wrk0,
     &                        cssb(-1,-1,1))
                else
                  cssb = 0.0_R8
                endif                                                    !iyv 06.02.14
#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
!! FR possible loop (3) for use with OpenMP
!$OMP PARALLEL DO DEFAULT(none) LASTPRIVATE(flob, conb, pb)
!$OMP& SHARED(nx, ny, rza, na, te, ti, ns, nsmin, nsmax, ismain)
!$OMP& SHARED(guard_flows, facdrift, fac_ExB, am, cssb)
!$OMP& SHARED(vol, pbs, pbshz, ua, uadia, vaecrb, vadia, wadia)
!$OMP& SHARED(ne, bb, hx, hy, hz, qz, gs, qc, po, cvla, cdna, cdpa)
!$OMP& SHARED(cddi, cvlahz, cvmahz, cdnahz, cdpahz, fchvispar)
!$OMP& SHARED(fchvisper, fchinert, fchanml, fchvisq)
!$OMP& SHARED(fna, fna_nodrift, fna_fcor, fna_he, fnaPSch, fna_mdf)
!$OMP& SHARED(fna_32, fna_52, fna_32nd, fna_52nd)
!$OMP& SHARED(kinrgy, fllim0fna) 
!                do is = 0, ns-1
                do is = nsmin, nsmax-1
*    ..compute partial pressure
                  call b2xppb (nx, ny,
     1             rza(-1,-1,is), na(-1,-1,is), te, ti, pb)
*    ..re-compute particle flux
                  call b2tfnb (nx, ny, is, ismain,                       !srv 08.09.01
     1             guard_flows, facdrift, fac_ExB, am(is)*mp, vol, cssb, !iyv 06.02.14
     1             rza(-1,-1,is), pbs, pbshz, na(-1,-1,is), pb,
     1             ua(-1,-1,is), uadia(-1,-1,0,is), vaecrb(-1,-1,0,is),  !xpb
     &             vadia(-1,-1,0,is), wadia(-1,-1,0,is),
     &             ne, te, ti, bb, hx, hy, hz, qz, gs, qc, po,           !srv 09.01.01
     &             cvla(-1,-1,0,is), cdna(-1,-1,0,is), cdpa(-1,-1,0,is),
     &             cddi(-1,-1,0,is), cvlahz(-1,-1,0,is),                 !srv 09.01.01 13.10.06
     &             cvmahz(-1,-1,0,is), cdnahz(-1,-1,0,is),               !srv 15.12.05
     &             cdpahz(-1,-1,0,is),
     &             fchvispar, fchvisper, fchinert, fchanml, fchvisq,     !srv 20.09.06
     &             fna(-1,-1,0,is), fna_nodrift(-1,-1,0,is),
     &             fna_fcor(-1,-1,0,is), fna_he(-1,-1,0,is),             !srv 01.10.99 21.01.00
     &             fnaPSch(-1,-1,0,is), fna_mdf(-1,-1,0,is),             !srv 13.10.06
     1             flob, conb, fna_32(-1,-1,0,is), fna_52(-1,-1,0,is),
     1             fna_32nd(-1,-1,0,is), fna_52nd(-1,-1,0,is),
     &             kinrgy(-1,-1,is), fllim0fna(-1,-1,0,is), .false.)     !djm Jan2017
                enddo
!$OMP END PARALLEL DO
!!FR End of loop (3)
#ifdef WITH_FTIMINGS
      call timer%stop("parallel")
#endif
                if(wrong_flow .or. ipass.ne.0) then                      !xpb
                   call print_wrong_flows(
     &                  'In potential loop: fna(ions) at plates')
                endif
                if (tmax.gt.0.0_R8) then
                  write(*,'(a,i6,1p,2g14.7,2i4)')
     &             'i,respo',i,t0,tmax,ixmax,iymax
                else
                  write(*,*) 'No potential residuals, respo = ', t0
                endif
                call xertst (ier0.eq.0, 'error return from b2nppo')
                if(.not.solving(3).or..not.any(solvepo(0:nnreg(0))))
     &           exit                                                    !srv 22.05.18
                if(t0.gt.potok.and.(i.eq.(potit/3).or.i.eq.(potit*2/3)
     &           .or.mod(i+1,50).eq.0)) then
                  if(vis_per.eq.0.0_R8) then
                    call restart_ma28
                  else
                    if (style.eq.2) call restart_ma28_for_7              !srv 11.01.00 29.05.02
                    if (style.eq.3) call restart_sdrv_for_7              !srv 29.05.02
                  endif
                endif
                i=i+1
              endif                                                      !xpb - wrong_flow
*   ..end of potential equation loop                                     !xpb
            enddo
            if(wrong_flow .or. ipass.ne.0) then                          !xpb
               call print_wrong_flows(
     &              'After potential loop: fna(ions) at plates')
            endif
            if(i.gt.potit) then
              ierr=10
cxpb
cxpb          ierr=0     !testing forcing potit iterations
cxpb
              write(*,*) 'Potential equation not converging'
              call b2news_facdrift_fac_ExB_dec(nx,ny,ramp_slow)
              goto 99999
            endif
          endif                         !xpb - wrong_flow .or. solving(3)
        endif                           ! if pot_eq
*   ..compute fne, fni
c     with diamagnetic drift part                                        !xpb
        call b2xpfe (nx, ny, ns, ns_ext, qe,
     .   rza, za_ext, fna, fa_ext, fch, fne)                             !xpb
        call b2xpfe (nx, ny, ns, ns_ext, qe,
     .   rza, za_ext, fna_32nd, fa_ext, fch_32, fne_32)
        call b2xpfe (nx, ny, ns, 0, qe,
     .   rza, zero, fna_52nd, zero, fch_52, fne_52)
        call b2xpfi (nx, ny, ns, fna, fni)                               !xpb
        call b2xpfi (nx, ny, ns, fna_32nd, fni_32)
        call b2xpfi (nx, ny, ns, fna_52nd, fni_52)
        call b2xppz (nx, ny, ns, ne, na, te, ti, pz)                     !srv 18.01.02
        if (ncall.eq.ncallout)
     >   write(*,*) 'b2news [a][32|52] ',
     1   maxval(abs(fna-fna_32-fna_52)), maxval(abs(fna))
*   ..compute heat fluxes
        call b2tfed (nx, ny, facdrift, fac_ExB, vol, hx, hy, qz, bb,
     &   te, po, vedia, veecrb)
        call b2tfch (nx, ny, ns, ismain, ismain0,
     1   guard_flows, facdrift, fac_ExB, fac_vis, BoRiS,                 !xpb
     1   vol, gs, hx, hy, qz, qc, bb, te, po, ne, ne2, lnlam, pz,        !xpb
     1   ti, ni(-1,-1,1), na, fna, uadia, vedia, wadia, vaecrb,          !xpb
     &   veecrb, ua, fne_32, fne_52, csig, csigin, chce, chve, chci,
     &   calf, chvemx, cvsa_cl, fch, conc, fchstoch, fchdia,             !srv 17.12.13 16.10.17
     &   fch_p, fchvispar, fchvisq, fchin, fchinert, fch_32, fch_52,     !srv 25.11.08
     1   avisper, bvisper, cvisper, fchvisper, csig_cl, calf_cl)
        call b2tfhe (nx, ny, ns, guard_flows, facdrift, BoRiS,           !xpb
     1   vol, hx, hy, qz, te, po, ne, na, vedia, veecrb,                 !srv 13.01.17
     1   fne_32, fne_52, fch_p,
     1   csig, chce, chve, calf, chvemx, cddi, fhe, fhe_mdf, fhePSch,    !srv 13.10.06
     &   fhet, floe, cone, fllime, fchstoch,                             !xpb !srv 17.12.13
     &   .true.)                                                         !djm Jan2017
        call b2tfhi (nx, ny, ns, guard_flows, facdrift, BoRiS,           !xpb
     1   vol, hx, hy, gs, qz, qc, bb, ti, ni, na,
     2   fna, fni_32, fni_52, kinrgy,
     &   chci, cddi, chvi, chvimx, cvsa, vadia, vaecrb,
     &   fchvispar, fchvisper, fchvisq, fchinert, fchanml,               !lk 30.05.11
     3   fhi, fhi_mdf, fhiPSch, fhit, floi, coni, fllimi,                !srv 13.10.06
     &   .true.)                                                         !djm Jan2017
        call b2xpen (nx, ny, ns,
     1   fna, fne, fch, fhe, fhi, kinrgy, rpt, te, ti, po, BoRiS,
     2   fhm, fhp, fhj, fht)
*   ..add contribution from cdpa to floe, cone, floi, coni
        call b2tfcc (nx, ny, ns, na, te, ti, cdpa, floe, cone,
     &   floi, coni)
*    ..add contributions from the divergent total diamagnetic flux
        call b2sdia (nx, ny, ns, na, ne, ti, te, bb, vol, hx, hy, qz,
     1   she, shi, facdrift)                                             !xpb
*   ..compute ue, electron parallel velocity
        ue=0.0_R8
        if(ue_style.ne.2) then                                           !srv 13.01.17
         do iy=-1,ny
           do ix=-1,nx
             if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
               ue(ix,iy)=-(fch_p(ix,iy,0)/pbs(ix,iy,0)+
     &          fch_p(rightix(ix,iy),rightiy(ix,iy),0)/
     &          pbs(rightix(ix,iy),rightiy(ix,iy),0))/qe/2.0_R8
             elseif(leftix(ix,iy).eq.-2) then
               ue(ix,iy)=-fch_p(rightix(ix,iy),rightiy(ix,iy),0)/
     &         pbs(rightix(ix,iy),rightiy(ix,iy),0)/qe
             elseif(rightix(ix,iy).eq.nx+1) then
               ue(ix,iy)=-fch_p(ix,iy,0)/pbs(ix,iy,0)/qe
             endif
           end do
         end do
#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
!!FR loop (4) over ns suitable for OpenMP parallelisation
!$OMP PARALLEL DO DEFAULT(none) SHARED(rza,ua,na,ns) REDUCTION(+:ue)
         do is=0,ns-1
           ue(:,:)=ue(:,:)+rza(:,:,is)*ua(:,:,is)*na(:,:,is)
         enddo
!$OMP END PARALLEL DO
!!FR end of loop (4)
#ifdef WITH_FTIMINGS
      call timer%stop("parallel")
#endif
         ue(:,:)=ue(:,:)+ue_ext(:,:)*ne_ext(:,:)
         ue(:,:)=ue(:,:)/ne(:,:)
        else
         call b2xpve (nx, ny, ns, qe, rza, rz2, alfx_c, sigx_c, po,      !srv 17.11.16 {
     &                ne, te, na, ua, ue_ext, ne_ext, ue)
        endif                                                            !srv 13.01.17
*   ..iterate on heat balance
        if(mdf_fhe.eq.0 .and. mdf_fhi.eq.0) then                         !srv 13.10.06
          call b2npht (nx, ny, ns, nnreg(0), ncall, ismain, solving,     !srv 05.07.17
     &     solveee, solveei, solveet, solvepo, solvemo,                  !srv 22.05.18
     &     rxf3, BoRiS, fac_ExB, hx, hy, qz, vol, bb,
     &     pbs, lnlam, floe, cone, floi, coni, ceqp, pccm, she, shi,     !srv 01.07.09
     &     fne, fch, cvsa, fchdia, uadia, wadia, vadia, vaecrb,          !srv 09.07.08 07.07.09
     &     na, ua, ue, te, ti, po, ne, ni, fhe, fhi,                     !srv 09.07.08
     &     reshe, reshi, wrk0, corte, corti, ier0)
        endif                                                            !srv 13.10.06 }
        if((mdf_fhe.ne.0 .and. mdf_fhi.ne.0) .and. use_b2nph9.eq.1) then !srv 13.10.06
          call b2nph9 (nx, ny, ns, nnreg(0), ncall, ismain, solving,     !srv 05.07.17
     &     solveee, solveei, solveet, solvepo, solvemo,                  !srv 22.05.18
     &     rxf3, BoRiS, fac_ExB, hx, hy, qz, vol, bb,                    !srv 07.06.00 02.06.99
     &     pbs, lnlam, floe, cone, floi, coni, ceqp, pccm, she, shi,     !srv 07.06.00 01.07.09
     &     fne, fch, cvsa, fchdia, uadia, wadia, vadia, vaecrb,          !srv 27.12.01 28.06.03 30.06.03 07.07.09
     &     na, ua, ue, te, ti, po, ne, ni, fhe_mdf, fhi_mdf, shef, shif, !srv 31.05.99
     &     reshe, reshi, wrk0, corte, corti, ier0)                       !srv 31.05.99 10.11.02
        endif                                                            !srv 13.10.06 }
        if((mdf_fhe.ne.0 .and. mdf_fhi.ne.0) .and. use_b2nph9.eq.0) then !srv 13.10.06
          call b2npht (nx, ny, ns, nnreg(0), ncall, ismain, solving,     !srv 05.07.17
     &     solveee, solveei, solveet, solvepo, solvemo,                  !srv 22.05.18
     &     rxf3, BoRiS, fac_ExB, hx, hy, qz, vol, bb,                    !srv 07.06.00 02.06.99
     &     pbs, lnlam, floe, cone, floi, coni, ceqp, pccm, she, shi,     !srv 07.06.00 01.07.09
     &     fne, fch, cvsa, fchdia, uadia, wadia, vadia, vaecrb,          !srv 27.12.01 28.06.03 30.06.03 07.07.09
     &     na, ua, ue, te, ti, po, ne, ni, fhe_mdf, fhi_mdf,             !srv 31.05.99
     &     reshe, reshi, wrk0, corte, corti, ier0)                       !srv 31.05.99 10.11.02
        endif                                                            !srv 13.10.06 }
        call xertst (ier0.eq.0, 'error return from b2npht')
        if(pot_eq.eq.0) po=0.0_R8     !xpb remove pot. correction if unneeded
        if(pot_eq.eq.2) po(:,:)=3.1_R8*te(:,:)/qe                        !xpb
        ipass = ipass + 1                                                !xpb
        if(ipass.gt.min(20,potit)) then                                  !xpb
          ierr=10                                                        !xpb
          write(*,*) 'Wrong flow correction inoperative'                 !xpb
          call b2news_facdrift_fac_ExB_dec(nx,ny,ramp_slow)
          goto 99999                                                     !xpb
        endif                                                            !xpb
*  ..end of wrong_flow loop                                              !xpb
c        write(*,*) 'Calling b2stbr_neutr_scl'
        call b2stbr_neutr_scl(nx, ny, ns)                                !srv 25.10.08
c        write(*,*) 'Called b2stbr_neutr_scl'
      enddo                                                              !xpb

* ..return
      ierr = 0
99999 continue
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2news

      contains

      subroutine print_wrong_flows(label)
      implicit none
      character (len=*) :: label
      integer iy
      logical inverted_ua

      inverted_ua = .true.
      write(*,*) label
      do iy = -1,ny
        if(topiy(-1,iy).eq.ny+1) cycle
        if(bottomiy(-1,iy).eq.-2) cycle
        if(region(-1,iy,0).eq.0) cycle
        if(.not.on_closed_surface(-1,iy))
     >        write(*,*) iy,'west =',
     >        fna(rightix(-1,iy),rightiy(-1,iy),0,ismain)
        inverted_ua = inverted_ua .and.
     .   ua(-1,iy,ismain)*pbs(rightix(-1,iy),rightiy(-1,iy),0).gt.0.0_R8
      enddo
      do iy = -1,ny
        if(topiy(nx,iy).eq.ny+1) cycle
        if(bottomiy(nx,iy).eq.-2) cycle
        if(region(nx,iy,0).eq.0) cycle
        if(.not.on_closed_surface(nx,iy))
     >        write(*,*) iy,'east =',fna(nx,iy,0,ismain)
        inverted_ua = inverted_ua .and.
     .   ua(nx,iy,ismain)*pbs(nx,iy,0).lt.0.0_R8
      enddo
      if (inverted_ua) then
        write(*,*) 'All target velocities have the wrong sign'
        write(*,*) 'Please consider using the b2mndr_inverse_ua switch'
        call xerrab ('Inverted velocity field!')
      end if
      end subroutine print_wrong_flows

      end subroutine b2news
c
      subroutine b2news_facdrift_fac_ExB_init(nx, ny, ramp_slow)
      use b2mod_types
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_geo
      use b2mod_work
      implicit none
      integer nx, ny, ramp_slow
      real (kind=R8) ::
     1 facdrift_target, facdrift_start, facdrift_inc, facdrift_dec,
     3 facExB_target, facExB_start, facExB_inc, facExB_dec,
     3 facvis_target, facvis_start, facvis_inc, facvis_dec,
     & facdrift_scalar, fac_ExB_scalar, fac_vis_scalar,
     & facdrift_tanh_a, facdrift_tanh_b, facExB_tanh_a, facExB_tanh_b,
     1 facvis_tanh_a, facvis_tanh_b
      integer ncall, iy_nocoreExB, fac_ref
      integer ix, iy, jxi, jxa, jsep
      real (kind=R8) cr, cz, dssep
      save
     1 facdrift_target, facdrift_start, facdrift_inc, facdrift_dec,
     1 facExB_target, facExB_start, facExB_inc, facExB_dec,
     1 facvis_target, facvis_start, facvis_inc, facvis_dec,
     3 facdrift_scalar, fac_ExB_scalar, ncall, iy_nocoreExB,
     & facdrift_tanh_a, facdrift_tanh_b, facExB_tanh_a, facExB_tanh_b,
     & fac_ref, fac_vis_scalar, facvis_tanh_a, facvis_tanh_b
      data
     2 facdrift_target/0.0_R8/, facdrift_start/0.0_R8/,
     3 facdrift_inc/1.0_R8/, facdrift_dec/0.0_R8/,
     2 facExB_target/0.0_R8/, facExB_start/0.0_R8/,
     3 facExB_inc/1.0_R8/, facExB_dec/0.0_R8/,
     2 facvis_target/0.0_R8/, facvis_start/0.0_R8/,
     3 facvis_inc/1.0_R8/, facvis_dec/0.0_R8/,
     & facdrift_scalar/0.0_R8/, fac_ExB_scalar/0.0_R8/,
     1 fac_vis_scalar/0.0_R8/,
     & ncall/0/, iy_nocoreExB/-2/,
     & facdrift_tanh_a/0.0_R8/, facdrift_tanh_b/0.0_R8/,
     & facExB_tanh_a/0.0_R8/, facExB_tanh_b/0.0_R8/,
     & facvis_tanh_a/0.0_R8/, facvis_tanh_b/0.0_R8/

      cr(ix,iy)=0.25_R8*
     1 (crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+crx(ix,iy,3))
      cz(ix,iy)=0.25_R8*
     1 (cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+cry(ix,iy,3))

      if(ncall.eq.0) then
       call get_jsep(nx,ny,jxi,jxa,jsep)
       fac_ref=jxa
       call ipgeti ('b2news_fac_ref',fac_ref)
       call xertst (-1.le.fac_ref.and.fac_ref.le.nx,
     2      'faulty internal parameter fac_ref')
       call ipgetr ('b2news_facdrift_start', facdrift_start)
       call ipgetr ('b2news_facdrift_target', facdrift_target)
       call ipgetr ('b2news_facdrift_inc', facdrift_inc)
       call ipgetr ('b2news_facdrift_dec', facdrift_dec)
       call xertst (0.0_R8.le.facdrift_start .and.                       !xpb
     &              1.0_R8.ge.facdrift_start,                            !xpb
     &   'faulty argument facdrift_start')                               !xpb
       call xertst (0.0_R8.le.facdrift_target .and.                      !xpb
     &              1.0_R8.ge.facdrift_target,                           !xpb
     &   'faulty argument facdrift_target')                              !xpb
       call xertst (1.0_R8.le.facdrift_inc,                              !xpb
     &   'faulty argument facdrift_inc')                                 !xpb
       call xertst (0.0_R8.le.facdrift_dec .and. 1.0_R8.ge.facdrift_dec,
     &   'faulty argument facdrift_dec')                                 !xpb
       facdrift_scalar=facdrift_start
       if(facdrift_start.eq.facdrift_target .and.
     &  facdrift_start.ne.0.0_R8) then
         write(*,*) 'facdrift set to a fixed value :', facdrift_start
       elseif(facdrift_start.ne.facdrift_target.and.ramp_slow.eq.0) then
         write(*,'(a,1p,1g9.2,a,1g9.2)')
     &    'facdrift ramped normally from ',facdrift_start,
     &    ' towards ',facdrift_target
       elseif(facdrift_start.ne.facdrift_target.and.ramp_slow.ne.0) then
         write(*,'(a,1p,1g9.2,a,1g9.2)')
     &    'facdrift ramped slowly from ',facdrift_start,
     &    ' towards ',facdrift_target
       else
         write(*,*) 'No diamagnetic drifts'
       endif
       call ipgetr ('b2news_facExB_start', facExB_start)
       call ipgetr ('b2news_facExB_target', facExB_target)
       call ipgetr ('b2news_facExB_inc', facExB_inc)
       call ipgetr ('b2news_facExB_dec', facExB_dec)
       call xertst (0.0_R8.le.facExB_start .and.
     &              1.0_R8.ge.facExB_start,
     &   'faulty argument facExB_start')
       call xertst (0.0_R8.le.facExB_target .and.
     &              1.0_R8.ge.facExB_target,
     &   'faulty argument facExB_target')
       call xertst (1.0_R8.le.facExB_inc,
     &   'faulty argument facExB_inc')
       call xertst (0.0_R8.le.facExB_dec .and. 1.0_R8.ge.facExB_dec,
     &   'faulty argument facExB_dec')
       call ipgetr ('b2news_ExB', fac_ExB_scalar)                        !xpb
       call xertst (0.0_R8.le.fac_ExB_scalar .and.
     &              1.0_R8.ge.fac_ExB_scalar,
     &  'faulty parameter b2news_ExB')                                   !xpb
       if(fac_ExB_scalar.ne.0.0_R8) then
c old style
         facExB_start=fac_ExB_scalar
         facExB_target=fac_ExB_scalar
         facExB_inc=1.0_R8
         facExB_dec=1.0_R8
         write(*,*) 'fac_ExB set to a fixed value: ', fac_ExB_scalar
       elseif(facExB_start.ne.0.0_R8.and.facExB_start.eq.facExB_target)
     >  then
         fac_ExB_scalar=facExB_start
         facExB_inc=1.0_R8
         facExB_dec=1.0_R8
         write(*,*) 'fac_ExB set to a fixed value: ', fac_ExB_scalar
       elseif(facExB_start.ne.0.0_R8.and.facExB_start.ne.facExB_target)
     >  then
         fac_ExB_scalar=facExB_start
         if (ramp_slow.eq.0) then
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_ExB ramped normally from ',facExB_start,
     &      ' towards ',facExB_target
         else
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_ExB ramped slowly from ',facExB_start,
     &      ' towards ',facExB_target
         endif
       else
         write(*,*) 'No ExB'
       endif
       call ipgeti ('b2news_iy_nocoreExB', iy_nocoreExB)                 !xpb
       call xertst (-2.le.iy_nocoreExB .and. ny.ge.iy_nocoreExB,         !xpb
     &   'faulty parameter iy_nocoreExB')                                !xpb
       do iy=-1, iy_nocoreExB
         do ix=-1,nx
           if(on_closed_surface(ix,iy)) then
             fac_ExB_profile(ix,iy)=0.0_R8
           endif
         enddo
       enddo
       call ipgetr ('b2news_facvis_start', facvis_start)
       call ipgetr ('b2news_facvis_target', facvis_target)
       call ipgetr ('b2news_facvis_inc', facvis_inc)
       call ipgetr ('b2news_facvis_dec', facvis_dec)
       call xertst (0.0_R8.le.facvis_start .and.
     &              1.0_R8.ge.facvis_start,
     &   'faulty argument facvis_start')
       call xertst (0.0_R8.le.facvis_target .and.
     &              1.0_R8.ge.facvis_target,
     &   'faulty argument facvis_target')
       call xertst (1.0_R8.le.facvis_inc,
     &   'faulty argument facvis_inc')
       call xertst (0.0_R8.le.facvis_dec .and. 1.0_R8.ge.facvis_dec,
     &   'faulty argument facvis_dec')
       call ipgetr ('b2news_vis', fac_vis_scalar)                        !xpb
       call xertst (0.0_R8.le.fac_vis_scalar .and.
     &              1.0_R8.ge.fac_vis_scalar,
     &  'faulty parameter b2news_vis')                                   !xpb
       if(fac_vis_scalar.ne.0.0_R8) then
c old style
         facvis_start=fac_vis_scalar
         facvis_target=fac_vis_scalar
         facvis_inc=1.0_R8
         facvis_dec=1.0_R8
         write(*,*) 'fac_vis set to a fixed value: ', fac_vis_scalar
       elseif(facvis_start.ne.0.0_R8.and.facvis_start.eq.facvis_target)
     >  then
         fac_vis_scalar=facvis_start
         facvis_inc=1.0_R8
         facvis_dec=1.0_R8
         write(*,*) 'fac_vis set to a fixed value: ', fac_vis_scalar
       elseif(facvis_start.ne.0.0_R8.and.facvis_start.ne.facvis_target)
     >  then
         fac_vis_scalar=facvis_start
         if (ramp_slow.eq.0) then
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_vis ramped normally from ',facvis_start,
     &      ' towards ',facvis_target
         else
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_vis ramped slowly from ',facvis_start,
     &      ' towards ',facvis_target
         endif
       else
         write(*,*) 'No vis'
       endif
       ds(-1)=
     1   sqrt((cr(fac_ref,-1)-
     1      0.5_R8*(crx(fac_ref,-1,0)+crx(fac_ref,-1,1)))**2+
     2      (cz(fac_ref,-1)-
     3      0.5_R8*(cry(fac_ref,-1,0)+cry(fac_ref,-1,1)))**2)
       do iy=0,ny
          ds(iy)=ds(iy-1)+
     1        sqrt((cr(fac_ref,iy)-cr(fac_ref,iy-1))**2+
     2        (cz(fac_ref,iy)-cz(fac_ref,iy-1))**2)
       enddo
       dssep=(ds(jsep)+ds(jsep+1))/2.0_R8
       do iy=-1,ny
         ds(iy)=ds(iy)-dssep
       enddo

       call ipgetr ('b2news_facdrift_tanh_a', facdrift_tanh_a)
       call ipgetr ('b2news_facdrift_tanh_b', facdrift_tanh_b)
       if(facdrift_tanh_b.ne.0.0_R8) then
          do iy=-1, ny
             do ix=-1, nx
                if(on_closed_surface(ix,iy)) then
                   facdrift_profile(ix,iy)=facdrift_profile(ix,iy)*
     &                  0.5_R8*(1.0_R8+
     &                    tanh((ds(iy)-facdrift_tanh_a)/
     &                   facdrift_tanh_b))
                endif
             enddo
          enddo
       endif
       call ipgetr ('b2news_fac_ExB_tanh_a', facExB_tanh_a)
       if (facExB_tanh_a.eq.0.0_R8)
     >  call ipgetr ('b2news_facExB_tanh_a', facExB_tanh_a)
       call ipgetr ('b2news_fac_ExB_tanh_b', facExB_tanh_b)
       if (facExB_tanh_b.eq.0.0_R8)
     >  call ipgetr ('b2news_facExB_tanh_b', facExB_tanh_b)
       if(facExB_tanh_b.ne.0.0_R8) then
          do iy=-1, ny
             do ix=-1, nx
                if(on_closed_surface(ix,iy)) then
                   fac_ExB_profile(ix,iy)=fac_ExB_profile(ix,iy)*
     &                  0.5_R8*(1.0_R8+
     &                    tanh((ds(iy)-facExB_tanh_a)/
     &                   facExB_tanh_b))
                endif
             enddo
          enddo
       endif
       call ipgetr ('b2news_fac_vis_tanh_a', facvis_tanh_a)
       if (facvis_tanh_a.eq.0.0_R8)
     >  call ipgetr ('b2news_facvis_tanh_a', facvis_tanh_a)
       call ipgetr ('b2news_fac_vis_tanh_b', facvis_tanh_b)
       if (facvis_tanh_b.eq.0.0_R8)
     >  call ipgetr ('b2news_facvis_tanh_b', facvis_tanh_b)
       if(facvis_tanh_b.ne.0.0_R8) then
          do iy=-1, ny
             do ix=-1, nx
                if(on_closed_surface(ix,iy)) then
                   fac_vis_profile(ix,iy)=fac_vis_profile(ix,iy)*
     &                  0.5_R8*(1.0_R8+
     &                    tanh((ds(iy)-facvis_tanh_a)/
     &                   facvis_tanh_b))
                endif
             enddo
          enddo
       endif
      endif
      goto 1000

      entry b2news_facdrift_fac_ExB_inc(nx,ny,ramp_slow)
      call xertst (0.lt.ncall,
     &     'b2news_facdrift_fac_ExB_inc called before initialisation')
      facdrift_scalar=min(facdrift_scalar*facdrift_inc,facdrift_target)
      fac_ExB_scalar=min(fac_ExB_scalar*facExB_inc,facExB_target)
      fac_vis_scalar=min(fac_vis_scalar*facvis_inc,facvis_target)
      goto 1000

      entry b2news_facdrift_fac_ExB_dec(nx,ny,ramp_slow)
      call xertst (0.lt.ncall,
     &     'b2news_facdrift_fac_ExB_dec called before initialisation')
      facdrift_scalar=max(facdrift_scalar*facdrift_dec,facdrift_start)
      fac_ExB_scalar=max(fac_ExB_scalar*facExB_dec,facExB_start)
      fac_vis_scalar=max(fac_vis_scalar*facvis_dec,facvis_start)
      goto 1000

 1000 continue
      facdrift=facdrift_scalar*facdrift_profile
      fac_ExB=fac_ExB_scalar*fac_ExB_profile
      fac_vis=fac_vis_scalar*fac_vis_profile
      write(*,'(1x,a,1p,6g12.5)') 'facdrift, fac_ExB, fac_vis',
     &   minval(facdrift),maxval(facdrift),
     &   minval(fac_ExB),maxval(fac_ExB),
     &   minval(fac_vis),maxval(fac_vis)
      ncall=ncall+1
      return

      end subroutine b2news_facdrift_fac_ExB_init
