*-----------------------------------------------------------------------
*.specification

      subroutine b2news_(nCv, nFc, nVx, ns, nscx, iscx,
     &  nscxmax, ismain, ismain0, dtim, switch, geo, mpg,
     &  st, st_ext, st_avg, ierr)
      use b2mod_types
      use b2mod_constants
      use b2mod_numerics_namelist
      use b2mod_b2cmpa
      use b2mod_time
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_ma28_for_us
     , , only : restart_ma28_for_us
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_tallies
     & , only : rrahireg, rranareg, rsahireg, rqradreg, rsanareg,
     &          b2bremreg, rqahereg, rramoreg, b2radreg, rqbrmreg,
     &          rsamoreg
      use b2mod_eirdiag
     & , only : tmb2, dib2, tab2, tib2, pefluxa, rfluxa, pfluxa,
     &          dmb2, lkindp, lkindm, lkindi, refluxa, dab2
      use b2mod_transport_namelist
     & , only : parm_hce, parm_hci, parm_dna, parm_dpa, parm_vla,
     &          parm_vsa, parm_alf, parm_sig
      use b2mod_recycle
     & , only : int0r, int0l, int3r, int3l, int6r, int6l, int2r, int2l,
     &          int5r, int5l, int1r, int1l, int4r, int4l, fna_mol
      use b2mod_boundary_namelist
     & , only : conpar, enepar, enipar, nbcd, enkpar, potpar, mompar,
     &          lfeedback
      use b2mod_neutrals_namelist
     & , only : maxw_eff, nstraid, userfluxparm
      use b2mod_b2cmpt
     & , only : cfvla, cflim, cfvsa, cfalf, cfdpa, cfsig, cfdna,
     &          cfhce, cfhci
      use b2mod_input_profile
     & , only : ndata, no_pflux, tr_ip_new_files, region_flags, tdata,
     &          tr_ip_elm_count, no_div, addspec, nss, nrr, nkind_coeff,
     &          nkind_data, poloidal_scaling, transport_ip_filename,
     &          nscale, transport_ip_time_switch, transport_ip_time_mod
      use b2us_feedback
     & , only : fb_target, fb_prev, fb_current, cum_volrec, fb_const,
     &          charge_frac, saved_fb_actuator, fb_rescale, dt_prev,
     &          feedback_namelist_used
      use b2mod_math
     & , only : cutlo, cutll, b2mod_math_initialised, small_r4_constant
#endif
      use b2mod_ad
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      integer nscx, nscxmax, iscx(0: nscxmax-1), ismain, ismain0
      real (kind=R8) :: dtim
      type (switches), intent (inout) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (inout) :: mpg
*   ..input/output arguments
      type (B2State), intent (inout) :: st
      type (B2StateExt), intent (inout) :: st_ext
      type (B2Average), intent (in) :: st_avg
*   ..output arguments                                                    !xpb
      integer ierr
      logical wrong_flow
*   ..common blocks

*-----------------------------------------------------------------------
*
*  1. purpose
*
*     B2NEWS_ performs one inner iteration for the B2 equations.
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
      integer is, ier0, i, iCv, ireg
      character charns*3                                                 !srv 30.06.03
      real (kind=R8) ::
     *  rxf0, rxf1, rxf2, rxf3, t0                                       !srv 09.07.99 06.07.00
      real (kind=R8) ::                                                  !xpb
     &  zeff(nCv),                                                       !srv 05.07.17 29.01.20
     &  rz(nCv),                                                         !iyv 06.02.14
     &  nep(nCv), wrk0(nCv),
     &  zero(nCv), zerof(nFc,0:1), wrkf(nFc),
     &  fni0(nFc,0:1)
      logical solvedum(0:CVREGMAX), solvedm2(0:CVREGMAX)
#ifdef WG_TODO
      logical inverted_ua
#endif
      logical new_matrix
*   ..procedures
      intrinsic mod
      real (kind=R8) ::
     *  smin, smax, b2sasum
      external xertst, ipgetr, smin, smax,
     &  b2scopy, b2saxpy, b2sasum, ipgeti
      external b2xvsg, b2xvps, b2xpne, b2xpni, b2xpnn, b2xpfe, b2xpfi,
     &  b2xppb, b2xpnm, b2tral, b2trql, b2treq, b2tfnb,
     &  b2npmo, b2npco, b2nppo, b2npht, b2nxdp, b2xppz  !srv 06.07.00 18.01.02
*   ..initialization
      integer, save :: no_neutr_scl = 0
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c     ingeger, save :: nstg2 = 0
!     2 facdrift_target, facdrift_start, facdrift_inc, facdrift_dec      !srv 09.07.99
!      real (kind=R8) ::                                                 !srv 04.06.99
!     *   fna_nodrift(-1:nx,-1:ny,0:1,0:1,0:ns-1)                        !srv 04.06.99
!     2 ,facdrift/0.0e0_R8/, facdrift_target/0.0e0_R8/, facdrift_start/0.0e0_R8/, !srv 09.07.99
!     3 facdrift_inc/0.0e0_R8/, facdrift_dec/0.0e0_R8/                 !srv 09.07.99
!     1  sna0(-1:nx,-1:ny,0:1,0:ns-1),                                    !xpb
!     1  sne0(-1:nx,-1:ny,0:1),                                           !xpb
!     2  shi0(-1:nx,-1:ny,0:3), she0(-1:nx,-1:ny,0:3),                    !xpb
!     3  smo0(-1:nx,-1:ny,0:3,0:ns-1),                                    !xpb
!     4  smq0(-1:nx,-1:ny,0:3,0:ns-1),                                    !xpb
!     5  sch0(-1:nx,-1:ny,0:3),                                           !xpb
!     2 facdrift_target, facdrift_start, facdrift_inc, facdrift_dec      !srv 09.07.99
!      integer, save :: re_compute_pb = 1                                !srv 15.02.10
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2news_')
*   ..set internal parameters on first call
      if (ncall_b2news_.eq.0) then
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

       b2news_solving(1:4)=switch%no_solve.eq.0
       do i=1,4
         b2news_solving(i)=b2news_solving(i).or.
     &    (switch%no_solve.lt.0.and.
     &     mod(abs(switch%no_solve),2**i).ge.2**(i-1))
       enddo
       if(.not.b2news_solving(1))
     >   write(*,*) 'Not solving parallel momentum equations!'
       if(.not.b2news_solving(2))
     >   write(*,*) 'Not solving continuity equations!'
       if(.not.b2news_solving(3))
     >   write(*,*) 'Not solving potential equation!'
       if(.not.b2news_solving(4))
     >   write(*,*) 'Not solving heat equations!'
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c       call ipgeti ('b2mndt_nstg2', nstg2)
c       if(nstg2 > 1 .and. no_neutr_scl .ne. 0) then
c         write(*,*) 'WARNING, b2mndt_nstg2>1 .and. no_neutr_scl.ne.0'
c       endif
c       call ipgeti ('b2news_no_neutr_scl', no_neutr_scl)
csw 26feb2014 some extra variables
c       call ipgeti ('eirene_neutr_avg', neutr_avg)
c       call ipgeti ('eirene_underrelax', nrelax)
c       call ipgeti ('b2stbr_eir_src_nhist', nhist)
c       call xertst (.not.(
c     .    (neutr_avg.gt.1 .and. (nhist.gt.1 .or. nrelax.gt.1)) .or.
c     .    (nhist.gt.1 .and. (neutr_avg.gt.1 .or. nrelax.gt.1)) .or.
c     .    (nrelax.gt.1 .and. (neutr_avg.gt.1 .or. nhist.gt.1))),
c     .   'Cannot have simultaneous multiple Eirene averaging schemes!')
!       call ipgetr ('b2news_facdrift_target', facdrift_target)          !srv 09.07.99
!       call ipgetr ('b2news_facdrift_start', facdrift_start)            !srv 09.07.99
!       call ipgetr ('b2news_facdrift_inc', facdrift_inc)                !srv 09.07.99
!       call ipgetr ('b2news_facdrift_dec', facdrift_dec)                !srv 09.07.99
!       facdrift=facdrift_start                                          !srv 09.07.99
!       call ipgeti ('b2news_re_compute_pb', re_compute_pb)              !srv 15.02.10
      endif
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFv')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
      call xertst (0.le.ismain0.and.ismain0.lt.ns,
     . 'invalid main neutral species index ismain0')
      call xertst (is_neutral(ismain0).or.ismain.eq.ismain0,
     . 'invalid main neutral species ismain0; must be neutral')
*   ..extensive tests on first few calls
      if (ncall_b2news_.lt.3) then
*    ..test sign of bb
       call b2xvsg (nCv, geo%cvBb(1,3), 1, 'bb3', '.gt.')
*    ..test cvVol, cvQgam
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
       call xertst (-1.0_R8.lt.smin(nCv,geo%cvQgam(1,1),1).and.
     &              smax(nCv,geo%cvQgam(1,1),1).lt.1.0_R8,
     .              'faulty argument range geo%cvQgam(,1)')
       call xertst (0.0_R8.lt.smin(nCv,geo%cvQgam(1,0),1).and.
     &              smax(nCv,geo%cvQgam(1,0),1).le.1.0_R8,
     .              'faulty argument range geo%cvQgam(,0)')
       call xertst (-1.0_R8.lt.smin(nFc,geo%fcQgam(1,1),1).and.
     &              smax(nFc,geo%fcQgam(1,1),1).lt.1.0_R8,
     .              'faulty argument range geo%fcQgam(,1)')
       call xertst (0.0_R8.lt.smin(nFc,geo%fcQgam(1,0),1).and.
     &              smax(nFc,geo%fcQgam(1,0),1).le.1.0_R8,
     .              'faulty argument range geo%fcQgam(,0)')
*    ..test present state
       call b2xvps (nCv, nFc, ns, st%pl, st%dv)
       call b2xvsg (nCv*ns, st%dv%kinrgy, 1, 'kinrgy', '.ge.')
      endif
*   ..set rxf0/1/2/3
*     Note: The potential equation needs an accurate solution on each
*     iteration; hence the different setting of rxf2.
      rxf0 = switch%b2mndt_rxf*switch%xfm0
      rxf1 = switch%b2mndt_rxf*switch%xfm1
      rxf2 = 1.0e0_R8 ! (rxf*xfm2)**0.2e0_R8                             !srv 15.01.03
      if (switch%b2mndt_rxf .eq. 0.0e0_R8) rxf2 = 0.0e0_R8
      rxf3 = switch%b2mndt_rxf*switch%xfm3

      zero  = 0.0_R8
      zerof = 0.0_R8
csw 01oct2012 initialize potential if necessary
      if (switch%pot_eq.eq.0) then                                       !srv 13.07.05 {
        st%pl%po = 0.0e0_R8
      else if (switch%pot_eq.eq.2) then
        st%pl%po = 3.1_R8*st%pl%te/qe
      endif

!      facdrift=min(facdrift*facdrift_inc,facdrift_target)               !srv 09.07.99
!      write(*,*) 'facdrift',facdrift                                    !srv 09.07.99
* ..main computation
*   ..compute ni, nn, ne, ne2
      call b2xpni (nCv, ns, st%pl%na, st%dv%ni)
      call b2xpnn (nCv, ns, st%pl%na, st%dv%nn)
      call b2xpne (nCv, ns, st%rt%rza, st%pl%na, st_ext%ne, st%dv%ne )
      call b2xpne (nCv, ns, st%rt%rz2, st%pl%na, st_ext%ne2, st%dv%ne2)
!*   ..compute ceqp                                                      !srv 03.07.03 {
!      call b2xpnm (nCv, ns, st%rt%rz2, am, st%pl%na, wrk0)
!      call b2treq (nCv, switch, geo%cvVol, st%pl%te, st%pl%ti,
!    &  st%pl%ne, st%dv%ni, wrk0, lnlam, st_ext, ceqp)                   !srv 03.07.03 }
*   ..compute flux limit coefficients
      call b2trql (nCv, nFc, ns, switch, geo, mpg, st%pl, st%dv, st_ext,
     &  st%co%chvemx, st%co%chvimx)
*   ..compute transport coefficients
      call b2tral (nCv, nFc, nVx, ns, nscx, nscxmax, iscx, ismain,
     1 switch, geo, mpg, st%pl, st%dv, st%rt, st_ext, st%co)
*   ..compute currents
      call b2tanml (nCv, nFc, ns, switch, geo, mpg, st%co%csig_an,
     &              st%pl%po, st%dv%fchanml)
      call b2tvspa (nCv, nFc, nVx, ns, switch, geo, mpg, st%pl%ua,
     &              st%co%vsaf_cl, st%dv%fac_vis, st%dv%fchvispar)
      call b2tiner (nCv, nFc, nVx, ns, switch, geo, mpg, st%pl%na,
     &              st%pl%ua, st%dv%facdrift, st%dv%fchinert)
*   ..iterate on momentum balance
!      if (ncall_b2news_.eq.0) then                                      !srv 01.10.99
      if(switch%ion_vlct_restrict .ne. 0 ) then                          !iyv 06.02.14
        call b2xprz (nCv, ns, mp, am, st%pl%na, rz, st_ext)
        call b2xppz (nCv, ns, st%dv%ne, st%pl%na, st%pl%te, st%pl%ti,
     &               st%dv%pz, st_ext)
        wrk0 = sqrt(st%dv%pz/rz)*geo%cvBb(:,0)/geo%cvBb(:,3)
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,wrk0,st%co%cssb)
      else
        st%co%cssb = 0.0_R8
      endif                                                              !iyv 06.02.14
!        do is=0,ns-1                                                    !srv 01.10.99
        do is = switch%nsmin, switch%nsmax-1
*    ..compute partial pressure
          call b2xppb (nCv, st%rt%rza(:,is), st%pl%na(:,is), st%pl%te,
     &      st%pl%ti, st%pl%tn, is, st%dv%pa(:,is))                      !srv 01.10.99
*    ..compute particle flux
          call b2tfnb (nCv, nFc, nVx, is, ismain,
     &      switch, geo, mpg, st%pl, st%dv, st%co, st%rt, .true.)
        enddo                                                            !srv 01.10.99
*      ..compute electron drift velocities
        call b2tfed (nCv, nFc, nVx, switch, geo, mpg, st%dv%facdrift,
     &    st%dv%fac_ExB, st%dv%ne, st%pl%te, st%pl%po, st%dv%vedia,
     &    st%dv%veecrb)
!      end if                                                            !srv 01.10.99

csw 10dec2012: added extra call if wanted/not-wanted (SOLPS4 has an
c   extra call at this place within internal iteration)
      if(switch%no_b2sral_call .eq. 0) then
        write(*,*) ' re-compute source '
*     ..re-compute source linearisation                                  !xpb
        call b2sral (nCv, nFc, nVx, ns,
     &    nscx, nscxmax, iscx, ismain, ismain0, dtim,
     &    switch, geo, mpg, st, st_ext, st_avg, wrong_flow, .false.)
#ifdef WG_TODO
        if(wrong_flow) then
          inverted_ua = .true.
          write(*,*) 'b2news_: wrong_flow returned from b2sral'
          do iy = -1,ny
            if(topiy(-1,iy).eq.ny+1) cycle
            if(bottomiy(-1,iy).eq.-2) cycle
            if(region(-1,iy,0).eq.0) cycle
            if(.not.on_closed_surface(-1,iy)) then
              write(*,*) iy,'west =',
     >           fna(rightix(-1,iy),rightiy(-1,iy),0,0,ismain)
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
              write(*,*) iy,'east =',fna(nx,iy,0,0,ismain)
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
#endif
*     ..smooth out sources
        call b2srsm (nCv, ns, dtim, switch, geo, mpg,
     &    st%pl%na, st%pl%ua, st%pl%te, st%pl%ti,
     %    st%dv%ne, st%dv%ni, st%dv%nn, st%sr, .false.)
*     ..include contributions due to timestep                            !xpb
        call b2srdt (nCv, ns, dtim, switch, geo, mpg,
     &    st%psnc%na, st%psnc%ua, st%psnc%te, st%psnc%ti, st%psnc%tn,
     &    st%psnc%ne, st%psnc%ni, st%psnc%nn, st%psnc%kinrgy,
     &    st%psnc%kt, st%psnc%zt,                                        !srv 23.07.21
     &    st%pl%na, st%pl%ua, st%pl%te, st%pl%ti, st%pl%tn, st%dv%ne,
     &    st%dv%ni, st%dv%nn, st%dv%kinrgy, st%pl%kt, st%pl%zt, st%sr,
     &    .false.)                                                       !srv 23.07.21
*     ..stabilise the source coefficients                                !xpb
        call b2srst (nCv, ns, switch, mpg,
     &    st%pl%na, st%pl%ua, st%pl%te, st%pl%ti, st%pl%tn, st%pl%po, !srv 01.07.08
     &    st%dv%ne, st%dv%ni, st%dv%nn, st%pl%kt, st%pl%zt, st%sr)    !srv 23.07.21
      endif
*   ..compute ue, electron parallel velocity
! lkw20.03.2023{
      wrkf = -st%dv%fch_p(:,0)/(qe*(geo%fcPbs + b2news_cutlo))
      call intcell(nFc,nCv,mpg,mpg%intcellP,wrkf,wrk0)
      do is=0,ns-1
        wrk0 = wrk0 + st%rt%rza(:,is)*st%pl%ua(:,is)*st%pl%na(:,is)
      enddo
      wrk0 = wrk0 + st_ext%ue*st_ext%ne
      wrk0 = wrk0/st%dv%ne                                                  !srv 05.03.11
      if(switch%b2news_iout.ne.0) then                                                 !srv 05.07.17
       call b2xzef(nCv, ns, st%rt%rz2, st%pl%na, st%dv%ne, zeff, st_ext)
       call b2xpic(nCv, ns, qe, st%rt%rza, geo%cvBb, st%pl%na,
     &   st%pl%ua, zeff, st_ext, st%dv%fch_pi_c)
       call b2xpif (nCv, nFc, switch, geo, mpg, st%dv%fch_pi_c,
     &   st%dv%fch_pi_f)
       call my_out_us(70,nCv,0,st%dv%fch_pi_c,'b2news__fch_pi_c')
       call my_out_us(70,nFc,1,st%dv%fch_pi_f,'b2news__fch_pi_f')
       call my_out_us(70,nCv,0,st%dv%ue,'b2news__ve')
       call my_out_us(70,nCv,0,wrk0,'b2news__ue')
      endif                                                              !srv 05.07.17
      st%dv%ue=wrk0
!lkw20.03.2023}
csw 11oct2012 added check on Coriolis force (taken from b2news)
      ier0=0
c Necessary for smooth restart with Rhie-Chow:
c The Rhie-Chow term in the particle fluxes uses pccb terms which are calculated in b2npmo.
c Hence to get a smooth restart after cp b2fstate bfstati we need to do a dummy call to
c b2npmo with rxf = 0 and recompute fluxes.
      if (ncall_b2news_.eq.0) then
        call b2npmo (nCv, nFc, nVx, ns, ismain,
     &    switch, geo, mpg,
     &    ncall_b2news_, b2news_solving(1), 0.0_R8,
     &    st%pl, st%dv, st%rt, st%co, st%sr, st_ext, ier0)
        do is = switch%nsmin, switch%nsmax-1
*    ..compute partial pressure
          call b2xppb (nCv, st%rt%rza(:,is), st%pl%na(:,is), st%pl%te,
     &      st%pl%ti, st%pl%tn, is, st%dv%pa(:,is))
*    ..compute particle flux
          call b2tfnb (nCv, nFc, nVx, is, ismain,
     &      switch, geo, mpg, st%pl, st%dv, st%co, st%rt, .true.)
        enddo
      endif
      call b2npmo (nCv, nFc, nVx, ns, ismain,
     &  switch, geo, mpg,
     &  ncall_b2news_, b2news_solving(1), rxf1,
     &  st%pl, st%dv, st%rt, st%co, st%sr, st_ext, ier0)
      call xertst (ier0.eq.0, 'error return from b2npmo')
*   ..set nep
*     (= ne before continuity update)
      call b2scopy (nCv, st%dv%ne, 1, nep, 1)
*   ..recompute currents
      call b2tanml (nCv, nFc, ns, switch, geo, mpg, st%co%csig_an,
     &              st%pl%po, st%dv%fchanml)
      call b2tvspa (nCv, nFc, nVx, ns, switch, geo, mpg, st%pl%ua,
     &              st%co%vsaf_cl, st%dv%fac_vis, st%dv%fchvispar)
      call b2tiner (nCv, nFc, nVx, ns, switch, geo, mpg, st%pl%na,
     &              st%pl%ua, st%dv%facdrift, st%dv%fchinert)
*     ..iterate on particle balance
!      do is = 0, ns-1
      do is = switch%nsmin, switch%nsmax-1
*    ..DPC
        do iCv = 1, nCv
          st%pl%na(iCv,is)=max(st%pl%na(iCv,is),switch%b2mndr_na_min) !srv 24.04.09
        enddo

!! option to impose a speed limit for fluid neutrals in case of spatially hybrid neutrals
        if (switch%spatial_hybrid.ne.0.and.is_neutral(is)) then
          do iCv = 1, nCv
            st%pl%ua(iCv,is)=sign(
     &                  min(abs(st%pl%ua(iCv,is)),switch%b2mndr_ua_max),
     &                          st%pl%ua(iCv,is))
          enddo
        endif

*    ..compute partial pressure
        call b2xppb (nCv, st%rt%rza(:,is), st%pl%na(:,is), st%pl%te,
     &    st%pl%ti, st%pl%tn, is, st%dv%pa(:,is))
*    ..compute particle flux
        call b2tfnb (nCv, nFc, nVx, is, ismain,
     &    switch, geo, mpg, st%pl, st%dv, st%co, st%rt, .true.)
*    ..perform iteration
csw 16oct2012 init ier0
        ier0=0
!       if (all(solveco(is,0:mpg%nnreg(0)))) then                        !srv 26.09.12 22.05.18
          solvedum(0:mpg%nnreg(0))=solveco(is,0:mpg%nnreg(0))
          solvedm2(0:mpg%nnreg(0))=solvemo(is,0:mpg%nnreg(0))
          new_matrix = .false.
          do ireg = 0, mpg%nnreg(0)
            new_matrix = new_matrix.or.
     &       solveco(is,ireg).neqv.last_solve_9(ireg)
          end do
          if(new_matrix.or..true.) call restart_ma28_for_us
          call b2npco (nCv, nFc, nVx, mpg%nnreg(0),                       !srv 22.05.18
     &      b2news_solving(2), solvedum, b2news_solving(1), solvedm2,
     &      ncall_b2news_, is, rxf0, switch, geo, mpg,
     &      st%pl, st%dv, st%sr, ier0)
          call xertst (ier0.eq.0, 'error return from b2npco')
          last_solve_9(0:mpg%nnreg(0))=solveco(is,0:mpg%nnreg(0))
!        endif                                                            !srv 13.10.06 } 22.05.18
        do iCv = 1, nCv
          st%pl%na(iCv,is)=max(st%pl%na(iCv,is),switch%b2mndr_na_min)    !srv 24.04.09
        enddo
!        if(re_compute_pb.eq.1) then                                     !srv 15.02.10
!*    ..re-compute partial pressure                                      !srv 01.10.99
!        call b2xppb (nCv, st%rt%rza(1,is), st%pl%na(1,is), st%pl%te,
!    &     st%pl%ti, st%pl%tn, is, pb)                                   !srv 01.10.99
!        endif                                                           !srv 15.02.10
*    ..re-compute particle flux
        call b2tfnb (nCv, nFc, nVx, is, ismain,
     &    switch, geo, mpg, st%pl, st%dv, st%co, st%rt, .false.)
      enddo
*      ..compute electron drift velocities
        call b2tfed (nCv, nFc, nVx, switch, geo, mpg, st%dv%facdrift,
     &    st%dv%fac_ExB, st%dv%ne, st%pl%te, st%pl%po,
     &    st%dv%vedia, st%dv%veecrb)

!WG_TODO      if(iout.ne.0 .or. iout_b2npco.ne.0 .or.                    !srv 16.07.10
!WG_TODO     &             .or. iout_b2npmo.ne.0 .or. iout_b2npht.ne.0
!WG_TODO     &             .or. iout_b2nph9.ne.0 .or. iout_b2nppo.ne.0
!WG_TODO     &             .or. iout_b2npp7.ne.0 .or. iout_b2wdat.eq.4)
!WG_TODO*    .. calculate contributions due to timestep for diagnostic
!WG_TODO        call b2srdt (nCv, ns, dtim, switch, geo, mpg,
!WG_TODO     &    st%psnc%na, st%psnc%ua, st%psnc%te, st%psnc%ti, st%psnc%tn,
!WG_TODO     &    st%psnc%ne, st%psnc%ni, st%psnc%nn, st%psnc%kinrgy,
!WG_TODO     &    st%psnc%kt, st%psnc%zt,                                        !srv 23.07.21
!WG_TODO     &    st%pl%na, st%pl%ua, st%pl%te, st%pl%ti, st%pl%tn, st%dv%ne,
!WG_TODO     &    st%dv%ni, st%dv%nn, st%dv%kinrgy, st%pl%kt, st%pl%zt, st%sr,
!WG_TODO     &    .true.)                                                       !srv 23.07.21

*   ..re-compute ne, ni, lnlam
      call b2xpne (nCv, ns, st%rt%rza, st%pl%na, st_ext%ne, st%dv%ne)
      call nancheck(nCv, st%dv%ne, 'b2news(ne): after call to b2xpne')
      call b2xpni (nCv, ns, st%pl%na, st%dv%ni)
      call b2xpnn (nCv, ns, st%pl%na, st%dv%nn)
      call b2xpne (nCv, ns, st%rt%rz2, st%pl%na, st_ext%ne2, st%dv%ne2)  !srv 11.07.99
      call b2xppz (nCv, ns, st%dv%ne, st%pl%na, st%pl%te, st%pl%ti,
     &             st%dv%pz, st_ext)                                     !srv 18.01.02

*   ..compute ceqp                                                       !srv 03.07.03 {
      call b2xpnm (nCv, ns, st%rt%rz2, am, st%pl%na, wrk0)
      call b2treq (nCv, switch, geo%cvVol, st%pl%te, st%pl%ti, st%dv%ne,
     &  st%dv%ni, wrk0, st%dv%lnlam, st_ext, st%co%ceqp)                 !srv 03.07.03 }
      if(switch%b2news_iout .ne. 0) then                                 !srv 12.07.05 {
        do is = 0,ns-1
          write(charns,'(i3.3)') is
          call my_out_us(70,nCv,0,st%pl%na(1,is),'b2news__na'//charns)
        enddo
        call my_out_us(70,nCv,0,st%dv%ni(1,0),'b2news__ni0')
        call my_out_us(70,nCv,0,st%dv%ni(1,1),'b2news__ni1')
        call my_out_us(70,nCv,0,st%dv%ne,'b2news__ne')
        call my_out_us(70,nCv,0,st%dv%pz,'b2news__pz')
        call my_out_us(70,nCv,0,st%co%ceqp,'b2news__ceqp')
      endif                                                              !srv 12.07.05 }

      if (switch%pot_eq.eq.0) then                                       !srv 13.07.05 {
        st%pl%po        = 0.0e0_R8
        !wdk call b2tfch_ to ensure all currents are properly set to 0
        call b2tfch_ (nCv, nFc, nVx, ns, ismain, ismain0, switch,
     &     geo, mpg, st%pl, st%dv, st%co, st%rt, st_ext)
      else if (switch%pot_eq.eq.2) then                                  !srv 11.01.12
        st%pl%po        = 3.1_R8*st%pl%te/qe
        !wdk call b2tfch_ to ensure all currents are properly set to 0
        call b2tfch_ (nCv, nFc, nVx, ns, ismain, ismain0, switch,
     &     geo, mpg, st%pl, st%dv, st%co, st%rt, st_ext)
      else                                                               !srv 13.07.05 }
*   ..correct po
        if (b2news_solving(3)) then
          call b2nxdp (nCv, mpg, solvepo, qe, st%pl%te, nep,
     &      st%dv%ne, wrk0)
          call b2saxpy (nCv, 1.0_R8, wrk0, 1, st%pl%po, 1)
          if(switch%b2news_iout .ne. 0) then                             !srv 12.07.05 {
            call my_out_us(70,nCv,0,st%pl%po,'b2nxdp_po')
          endif                                                          !srv 12.07.05 }
        endif
*   ..compute currents                                                   !srv 29.01.20 {
        call b2tfch_ (nCv, nFc, nVx, ns, ismain, ismain0, switch,
     &     geo, mpg, st%pl, st%dv, st%co, st%rt, st_ext)
        if (b2news_solving(3)) then
*   ..solve the potential equation
          if (switch%fhe_vis_per.eq.0.0e0_R8.or..true.) then             !srv 06.07.00
            new_matrix = .false.
            do ireg = 0, mpg%nnreg(0)
             new_matrix = new_matrix.or.
     &        solvepo(ireg).neqv.last_solve_9(ireg)
            end do
            if (new_matrix.or..true.) call restart_ma28_for_us
              call b2nppo (nCv, nFc, nVx, mpg%nnreg(0),                  !srv 22.05.18
     &        b2news_solving(3), solvepo, ncall_b2news_, rxf2, switch,
     &        geo, mpg, st%pl, st%dv, st%sr, ier0)
            last_solve_9(0:mpg%nnreg(0))=solvepo(0:mpg%nnreg(0))
          else
             call xerrab ('visper not yet available for WG')
#ifdef WG_TODO
            call b2npp7 (nx, ny, mpg%nnreg(0), b2news_solving(3),        !srv 27.01.00 22.05.18
     &       solvepo, ncall_b2news_, rxf2, ne, te, conc, avisper,        !srv 27.01.00
     &       bvisper, cvisper, sch, po, fch, respo, corpo, ier0)         !srv 27.01.00
#endif
          endif                                                          !srv 06.07.00
          t0 = b2sasum(nCv,st%dv%respo,1)
          write(*,*) 'respo before solving',t0
          call xertst (ier0.eq.0, 'error return from b2nppo')
*   ..recompute currents
          call b2tfch_ (nCv, nFc, nVx, ns, ismain, ismain0, switch,
     &     geo, mpg, st%pl, st%dv, st%co, st%rt, st_ext)
*   ..compute new residual following the tradition
          call b2ursd (nCv, nFc, mpg, st%dv%ne, st%pl%po,
     &     st%sr%sch, st%dv%fch, wrk0)
          t0 = b2sasum(nCv,wrk0,1)
          write(*,*) 'respo after solving',t0
        else
         if(switch%get_residuals.eq.1)
     &    call b2ursd (nCv, nFc, mpg,
     &                 st%dv%ne, st%pl%po, st%sr%sch, st%dv%fch,
     &                 st%dv%respo)
        endif
      endif           ! if pot_eq                                        !srv 13.07.05 29.01.20 }

      if(switch%b2news_iout .ne. 0) then                                 !srv 12.07.05 {
        write (charns,'(i3.3)') i
        call my_out_us(70,nCv,0,st%dv%respo,'b2news__respo'//charns)
        call my_out_us(70,nCv,0,st%pl%po,'b2news__po'//charns)
      endif                                                              !srv 12.07.05 }

      if(switch%ion_vlct_restrict .ne. 0 ) then                          !iyv 06.02.14
        call b2xprz (nCv, ns, mp, am, st%pl%na, rz, st_ext)
        call b2xppz (nCv, ns, st%dv%ne, st%pl%na, st%pl%te, st%pl%ti,
     &               st%dv%pz, st_ext)
        wrk0 = sqrt(st%dv%pz/rz)*geo%cvBb(:,0)/geo%cvBb(:,3)
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,wrk0,st%co%cssb)
      else
        st%co%cssb = 0.0_R8
      endif
      do is = 0, ns-1
*   ..compute partial pressure
        call b2xppb (nCv, st%rt%rza(:,is), st%pl%na(:,is), st%pl%te,
     &    st%pl%ti, st%pl%tn, is, st%dv%pa(:,is))
*   ..re-compute particle flux
#ifdef DBG
        write(*,*) 'is, za(is) ',is, (zamin(is)+zamax(is))/2.0_R8
#endif
        call b2tfnb (nCv, nFc, nVx, is, ismain,
     &   switch, geo, mpg, st%pl, st%dv, st%co, st%rt, .false.)
      enddo

*   ..compute fne_he, fni_he
c      with diamagnetic drift part                                       !srv 13.05.99
      call b2xpfe (nCv, nFc, ns, st_ext%ns, geo, mpg, qe,
     &  st%rt%rza, st_ext%za, st%dv%fna_he, st_ext%fa,
     &  st%dv%fch_p, st%dv%fne_he)                                       !srv 13.05.99 21.01.00 11.10.02
      call b2xpfe (nCv, nFc, ns, st_ext%ns, geo, mpg, qe,
     &  st%rt%rza, st_ext%za, st%dv%fna_53, st_ext%fa,
     &  st%dv%fch_p, st%dv%fne_53)
      call b2xpfe (nCv, nFc, ns, st_ext%ns, geo, mpg, qe,
     &  st%rt%rza, st_ext%za, st%dv%fna_32nd, st_ext%fa,
     &  st%dv%fch_32, st%dv%fne_32)
      call b2xpfe (nCv, nFc, ns, 0, geo, mpg, qe,
     &  st%rt%rza, zero, st%dv%fna_52nd, zerof,
     &  st%dv%fch_52, st%dv%fne_52)
      if (switch%tn_style.eq.0) then
        ! ti-equation: all species
        call b2xpfi (nFc, ns, st%dv%fna_he, st%dv%fni_he)                !srv 13.05.99 21.01.00
        call b2xpfi (nFc, ns, st%dv%fna_32nd, st%dv%fni_32)              !xpb
        call b2xpfi (nFc, ns, st%dv%fna_52nd, st%dv%fni_52)              !xpb
      elseif (switch%tn_style.eq.1) then
        ! ti-equation: ions only
        call b2xpfz (nFc, ns, st%dv%fna_he, st%dv%fni_he)
        call b2xpfz (nFc, ns, st%dv%fna_32nd, st%dv%fni_32)
        call b2xpfz (nFc, ns, st%dv%fna_52nd, st%dv%fni_52)
      else
        ! ti-equation: all species except hydrogenic neutrals
        ! tn-equation: hydrogenic neutrals only
        call b2xpfi (nFc, ns, st%dv%fna_he, fni0)
        call b2xpfn (nFc, ns, st%dv%fna_he, st%dv%fnn)
        st%dv%fni_he = fni0 - st%dv%fnn
        call b2xpfi (nFc, ns, st%dv%fna_32nd, fni0)
        call b2xpfn (nFc, ns, st%dv%fna_32nd, st%dv%fnn_32)
        st%dv%fni_32 = fni0 - st%dv%fnn_32
        call b2xpfi (nFc, ns, st%dv%fna_52nd, fni0)
        call b2xpfn (nFc, ns, st%dv%fna_52nd, st%dv%fnn_52)
        st%dv%fni_52 = fni0 - st%dv%fnn_52
      endif
*      ..compute electron drift velocities
      call b2tfed (nCv, nFc, nVx, switch, geo, mpg, st%dv%facdrift,
     &  st%dv%fac_ExB, st%dv%ne, st%pl%te, st%pl%po,
     &  st%dv%vedia, st%dv%veecrb)
      call b2tfhe_ (nCv, nFc, nVx, ns, switch, geo, mpg,                 !srv 29.01.20
     &  st%pl, st%dv, st%co, st%rt, st_ext, .true.)                      !djm Jan2017
csw 30apr2014 added fna_he, 17sep2014 added ua,cvsa
      call b2tfhi_(nCv, nFc, nVx, ns, ismain, switch, geo, mpg,          !srv 13.10.06
     &  st%pl, st%dv, st%co, .true.)                                     !djm Jan2017
       st%dv%floi_noc = st%dv%floi                                       !srv 24.10.02
       st%dv%floe_noc = st%dv%floe                                       !srv 24.10.02

      call b2xpen (nCv, nFc, ns, mpg,
     &  st%dv%fna, st%dv%fna_53, st%dv%fne_53,
     &  st%dv%fch, st%dv%fhe, st%dv%fhi, st%dv%fhn, st%rt%rpt, st%pl%te,
     &  st%pl%ti, st%pl%tn, st%pl%po, switch%BoRiS,
     &  st%dv%fhm, st%dv%fhp, st%dv%fhj, st%dv%fht)

*   ..add contribution from cdpa to floe, cone, floi, coni
      call b2tfcc (nCv, nFc, nVx, ns, switch, geo, mpg,
     &  st%pl, st%dv, st%co, st%rt)
*   ..compute ue, electron parallel velocity                             !srv 05.03.11 {
!lkw20.03.2023{
      wrkf = -st%dv%fch_p(:,0)/(qe*(geo%fcPbs + b2news_cutlo))
      call intcell(nFc,nCv,mpg,mpg%intcellP,wrkf,st%dv%ue)
      do is=0,ns-1
        st%dv%ue = st%dv%ue + 
     &           st%rt%rza(:,is)*st%pl%ua(:,is)*st%pl%na(:,is)
      enddo
      st%dv%ue = st%dv%ue + st_ext%ue*st_ext%ne
      st%dv%ue = st%dv%ue/st%dv%ne                                       !srv 05.03.11
!lkw20.03.2023}

*   ..iterate on heat balance
      call b2npht (nCv, nFc, nVx, ns, switch, geo, mpg,
     &  ncall_b2news_, ismain, b2news_solving, solveee, solveei,         !srv 22.05.18
     &  solveen, solveet, solvepo, solvemo, solvekt, solvezt, rxf3,
     &  st%pl, st%dv, st%co, st%rt, st%sr, st%srw, st_ext, ier0)
csw 16aug2011 init potential
      if(switch%pot_eq.eq.0) st%pl%po=0.0_R8                             !srv 13.07.05
      if(switch%pot_eq.eq.2) st%pl%po=3.1_R8*st%pl%te/qe
      call xertst (ier0.eq.0, 'error return from b2npht')

csw 12oct2012, extra call to relax particle balance if wanted (as in SOLPS4)
      if(switch%do_2nd_b2npco_call.ne.0) then
c   recompute particle fluxes
       do is=0,ns-1                                                      !srv 01.10.99
*    ..DPC
        do iCv = 1, nCv
          st%pl%na(iCv,is)=max(st%pl%na(iCv,is),switch%b2mndr_na_min) !srv 24.04.09
        enddo
*    ..compute partial pressure
        call b2xppb (nCv, st%rt%rza(:,is), st%pl%na(:,is), st%pl%te,
     &    st%pl%ti, st%pl%tn, is, st%dv%pa(:,is))
*    ..compute particle flux
        call b2tfnb (nCv, nFc, nVx, is, ismain,
     &    switch, geo, mpg, st%pl, st%dv, st%co, st%rt, .false.)
*    ..perform iteration
        new_matrix = .false.
        solvedum(0:mpg%nnreg(0))=solveco(is,0:mpg%nnreg(0))
        solvedm2(0:mpg%nnreg(0))=solvemo(is,0:mpg%nnreg(0))
        do ireg = 0, mpg%nnreg(0)
          new_matrix = new_matrix.or.
     &     solveco(is,ireg).neqv.last_solve_9(ireg)
        end do
        if(new_matrix) call restart_ma28_for_us
        call b2npco (nCv, nFc, nVx, mpg%nnreg(0),                     !srv 22.05.18
     &    b2news_solving(2), solvedum, b2news_solving(1), solvedm2,
     &    ncall_b2news_, is, rxf0, switch, geo, mpg,
     &    st%pl, st%dv, st%sr, ier0)
        last_solve_9(0:mpg%nnreg(0))=solveco(is,0:mpg%nnreg(0))
        call xertst (ier0.eq.0, 'error return from b2npco')
       enddo
*      ..compute electron drift velocities
       call b2tfed (nCv, nFc, nVx, switch, geo, mpg,
     &   st%dv%facdrift, st%dv%fac_ExB, st%dv%ne, st%pl%te, st%pl%po,
     &   st%dv%vedia, st%dv%veecrb)
      endif ! do_2nd_b2npco

csw 06dec2012 finally re-evaluate particle fluxes
c   (necessary to work out new rescaling factors sclrtio/sclalph in neutr_scl)
      if(switch%re_eval_prtls_fluxes.eq.1.or.
     &                                      switch%use_eirene.ne.0) then !srv 28.06.13 {
       do is=0,ns-1
        do iCv = 1, nCv
          st%pl%na(iCv,is)=max(st%pl%na(iCv,is),switch%b2mndr_na_min)    !srv 24.04.09
        enddo
*    ..compute partial pressure
        call b2xppb (nCv, st%rt%rza(:,is), st%pl%na(:,is), st%pl%te,     !srv 01.10.99
     &    st%pl%ti, st%pl%tn, is, st%dv%pa(:,is))
*    ..compute particle flux
        call b2tfnb (nCv, nFc, nVx, is, ismain,
     &    switch, geo, mpg, st%pl, st%dv, st%co, st%rt, .false.)
       enddo                                                             !srv 01.10.99
*      ..compute electron drift velocities
        call b2tfed (nCv, nFc, nVx, switch, geo, mpg, st%dv%facdrift,
     &    st%dv%fac_ExB, st%dv%ne, st%pl%te, st%pl%po,
     &    st%dv%vedia, st%dv%veecrb)
      endif                                                              !srv 28.06.13 }

csw 26oct2011 recalculate rescaling factors sclrtio/sclalph (SOLPS4)
      if(no_neutr_scl .eq. 0) then
        write(*,*) 'Calling b2stbr_neutr_scl'
        call b2stbr_neutr_scl(nCv, nFc, ns, mpg, geo, st%pl, st%dv,      !srv 25.10.08
     &    switch%use_eirene, switch%ank_mods)
        write(*,*) 'Called b2stbr_neutr_scl'
      endif

* ..return
      ierr = 0
      ncall_b2news_ = ncall_b2news_+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2news_

      end subroutine b2news_
