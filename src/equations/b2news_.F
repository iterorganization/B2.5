*-----------------------------------------------------------------------
*.specification

      subroutine b2news_(nCv, nFc, nVx, ns, nxtl, nxtr, nscx, iscx,
     &  nscxmax, ismain, ismain0, dtim, rxf, geo, mpg, pl, dv, co, rt, 
     &  rtw, sr, srw, st_ext, psnc, ierr)
      use b2mod_types
!      use b2mod_geo
!      use b2mod_plasma
!      use b2mod_rates
!      use b2mod_residuals
!      use b2mod_sources
!      use b2mod_transport
!      use b2mod_anomalous_transport            
!      use b2mod_work
!      use b2mod_indirect
      use b2mod_constants
      use b2mod_numerics_namelist
      use b2mod_b2cmpa
!      use b2mod_external
      use b2mod_time
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_ma28_for_us
     , , only : restart_ma28_for_us
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, nxtl, nxtr
      integer nscx, nscxmax, iscx(0: nscxmax-1), ismain, ismain0
      real (kind=R8) ::
     *   dtim, rxf
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
*   ..input/output arguments
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Coeff), intent (inout) :: co
      type (B2Rates), intent (inout) :: rt
      type (B2RatesWork), intent (inout) :: rtw
      type (B2Source), intent (inout) :: sr
      type (B2SourceWork), intent (inout) :: srw
      type (B2StateExt), intent (inout) :: st_ext
      type (B2PlasmaSnapshot), intent (inout) :: psnc
*   ..output arguments                                                    !xpb
      integer ierr
      logical wrong_flow
*   ..common blocks

*-----------------------------------------------------------------------
*
*  1. purpose
*
*     B2NEWS_ performs one inner iteration for the b2 equations.
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
      integer nsmin, nsmax, no_current
      integer, save :: style = 3, pot_eq = 1                             !srv 29.05.02 13.07.05
      integer, save :: mdf_fnb = 0, mdf_fhe = 0, mdf_fhi = 0             !srv 13.10.06
      integer, save :: no_solve = 0                                      !srv 02.07.08
      integer, save :: get_residuals = 0                                 !srv 22.05.18
!      integer, save :: re_compute_pb = 1                                !srv 15.02.10
      integer, save :: iout = 0,
     &                 iout_b2npc9 = 0, iout_b2npco = 0,                 !srv 11.09.09 {
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npmo = 0, iout_b2npp7 = 0,
     &                 iout_b2nppo = 0, iout_b2wdat = 0
      integer, save :: ion_vlct_restrict = 0                             !iyv 06.02.14
      integer, save :: ue_style = 2                                      !srv 17.11.16
      real (kind=R8), save :: cutlo = 1.0e-60_R8
      integer ncall, is, ier0, i, iCv, ireg
      character charns*3                                                 !srv 30.06.03
      real (kind=R8), save :: BoRiS = 0.0_R8
      real (kind=R8), save :: na_min = 1.0e4_R8, na_new = 1.0e14_R8
      real (kind=R8) :: 
     *  xfm0, xfm1, xfm2, xfm3, rxf0, rxf1, rxf2, rxf3,
     1  t0, vis_per                                                      !srv 09.07.99 06.07.00
!     2 facdrift_target, facdrift_start, facdrift_inc, facdrift_dec      !srv 09.07.99
      real (kind=R8) ::                                                  !xpb
!     1  sna0(-1:nx,-1:ny,0:1,0:ns-1),                                    !xpb
!     1  sne0(-1:nx,-1:ny,0:1),                                           !xpb
!     2  shi0(-1:nx,-1:ny,0:3), she0(-1:nx,-1:ny,0:3),                    !xpb
!     3  smo0(-1:nx,-1:ny,0:3,0:ns-1),                                    !xpb
!     4  smq0(-1:nx,-1:ny,0:3,0:ns-1),                                    !xpb
!     5  sch0(-1:nx,-1:ny,0:3),                                           !xpb
     &  ue(nCv), lnlam(nCv), zeff(nCv), ehx(nFc),                        !srv 05.07.17 29.01.20
     &  rz(nCv), cssb(nFc,0:1),                                          !iyv 06.02.14
     &  nep(nCv), wrk0(nCv), wrk1(nCv),
     &  zero(nCv), zerof(nFc,0:1), wrkf(nFc),
     &  respo0(nCv), corpo0(nCv)
      logical solving(4), solvedum(0:CVREGMAX), solvedm2(0:CVREGMAX)
      logical inverted_ua, new_matrix
      save ncall, xfm0, xfm1, xfm2, xfm3                                 !srv 09.07.99
     &     , vis_per, solving, no_current                                !srv 06.07.00
     &     , nsmin, nsmax
!     2 facdrift_target, facdrift_start, facdrift_inc, facdrift_dec      !srv 09.07.99
!      real (kind=R8) ::                                                 !srv 04.06.99
!     *   fna_nodrift(-1:nx,-1:ny,0:1,0:1,0:ns-1)                        !srv 04.06.99
*   ..procedures
      intrinsic mod
      real (kind=R8) ::
     *  smin, smax, b2sasum
      integer len_of_digits
      external len_of_digits
      external subini, subend, xertst, ipgetr, smin, smax,
     &  b2scopy, b2saxpy, b2sasum, ipgeti
      external b2xvsg, b2xvps, b2xpne, b2xpni, b2xpfe, b2xpfi, b2xppb,
     &  b2xpnm, b2tral, b2trql, b2treq, b2tfnb, b2tfhi,
     &  b2npmo, b2npco, b2nppo, b2npht, b2nxdp, b2xppz  !srv 06.07.00 18.01.02
*   ..initialization
      data ncall/0/, xfm0/1.0e0_R8/, xfm1/1.0e0_R8/, xfm2/1.0e0_R8/,
     1   xfm3/1.0e0_R8/,
     2   vis_per/0.0e0_R8/, no_current/0/                              !srv 06.07.00
!     2 ,facdrift/0.0e0_R8/, facdrift_target/0.0e0_R8/, facdrift_start/0.0e0_R8/, !srv 09.07.99
!     3 facdrift_inc/0.0e0_R8/, facdrift_dec/0.0e0_R8/                 !srv 09.07.99
      integer, save :: no_neutr_scl = 0, no_b2sral_call = 1
      integer, save :: do_2nd_b2npco_call = 0
      integer, save :: re_eval_prtls_fluxes = 0                        !srv 28.06.13
      integer, save :: nhist = 0, neutr_avg = 0, nrelax = 0, nstg2 = 0
      integer, save :: use_eirene = 0, tn_style = 0
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2news_')
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
!       call ipgetr ('b2news_facdrift_target', facdrift_target)          !srv 09.07.99
!       call ipgetr ('b2news_facdrift_start', facdrift_start)            !srv 09.07.99
!       call ipgetr ('b2news_facdrift_inc', facdrift_inc)                !srv 09.07.99
!       call ipgetr ('b2news_facdrift_dec', facdrift_dec)                !srv 09.07.99
!       facdrift=facdrift_start                                          !srv 09.07.99
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
       call ipgeti ('b2npmo_iout', iout_b2npmo)
       call ipgeti ('b2npp7_iout', iout_b2npp7)
       call ipgeti ('b2nppo_iout', iout_b2nppo)
       call ipgeti ('b2nph9_iout', iout_b2nph9)
       call ipgeti ('b2npht_iout', iout_b2npht)                          !srv 16.07.10 }
       call ipgeti ('b2wdat_iout', iout_b2wdat)                          !iys 17.04.17
!       call ipgeti ('b2news_re_compute_pb', re_compute_pb)              !srv 15.02.10
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
* SOLVEPO, SOLVEMT, SOLVEET array from b2.numerics.parameters.
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
       call ipgeti ('b2news_no_neutr_scl', no_neutr_scl)   
       call ipgeti ('b2news_no_b2sral_call', no_b2sral_call)
       call ipgeti ('b2news_do_2nd_b2npco_call', do_2nd_b2npco_call)
       call ipgeti ('b2news_re_eval_prtls_fluxes', re_eval_prtls_fluxes) !srv 28.06.13

       call ipgeti ('b2mn_tn_style', tn_style)

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
      endif                                                              !srv 02.07.08
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
      if (ncall.lt.3) then
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
       call b2xvps (nCv, nFc, ns, pl, dv)
       call b2xvsg (nCv*ns, dv % kinrgy, 1, 'kinrgy', '.ge.')
      endif
*   ..set rxf0/1/2/3
*     Note: The potential equation needs an accurate solution on each
*     iteration; hence the different setting of rxf2.
      rxf0 = rxf*xfm0
      rxf1 = rxf*xfm1
      rxf2 = 1.0e0_R8 ! (rxf*xfm2)**0.2e0_R8                             !srv 15.01.03
      if (rxf .eq. 0.0e0_R8) rxf2 = 0.0e0_R8
      rxf3 = rxf*xfm3

      zero  = 0.0_R8
      zerof = 0.0_R8
csw 01oct2012 initialize potential if necessary
      if (pot_eq.eq.0) then                                              !srv 13.07.05 {
        pl%po = 0.0e0_R8
      else if (pot_eq.eq.2) then
        pl%po = 3.1_R8*pl%te/qe
      endif

!      facdrift=min(facdrift*facdrift_inc,facdrift_target)               !srv 09.07.99
!      write(*,*) 'facdrift',facdrift                                    !srv 09.07.99
* ..main computation
*   ..compute ni, ne, ne2
      call b2xpni (nCv, ns, pl%na,  dv%ni)
      call b2xpne (nCv, ns, rt%rza, pl%na, st_ext%ne,  dv%ne )
      call b2xpne (nCv, ns, rt%rza, pl%na, st_ext%ne2, dv%ne2)
!*   ..compute ceqp                                                      !srv 03.07.03 {
!      call b2xpnm (nx, ny, ns, rz2, am, na, wrk0)
!      call b2treq (nx, ny, vol, te, ti, ne, ni, wrk0, lnlam, ceqp)      !srv 03.07.03 }
*   ..compute flux limit coefficients
      call b2trql (nCv, nFc, ns, geo, mpg, pl, dv, st_ext, 
     &  co%chvemx, co%chvimx)
*   ..compute transport coefficients
      call b2tral (nCv, nFc, nVx, ns, nscx, nscxmax, iscx, ismain,
     1 geo, mpg, pl, dv, rt, st_ext, co)

      call b2tanml (nCv, nFc, ns, geo, mpg, co%csig_an,pl%po,dv%fchanml)
      call b2tvspa (nCv, nFc, nVx, ns, geo, mpg, pl%ua, co%vsaf_cl, 
     &  dv%fac_vis, dv%fchvispar)
      call b2tiner (nCv, nFc, nVx, ns, geo, mpg, pl%na, pl%ua,
     &  dv%facdrift, dv%fchinert)
#ifdef WG_TODO
*   ..iterate on momentum balance
!      if (ncall.eq.0) then                                              !srv 01.10.99
      if(ion_vlct_restrict .ne. 0 ) then                                 !iyv 06.02.14
        call b2xprz (nx, ny, ns, mp, am, na, rz)
        call b2xppz (nx, ny, ns, ne, na, te, ti, pz)
        wrk0(:,:)=sqrt(pz(:,:)/rz(:,:))*abs(bb(:,:,0)/bb(:,:,3))
        call intfaceh(nx,ny,leftix,leftiy,vol,wrk0,cssb(-1,-1,0))
        call intfacev(nx,ny,bottomix,bottomiy,vol,wrk0,cssb(-1,-1,1))
      else
        cssb = 0.0_R8
      endif   
#else
      cssb = 0.0_R8
#endif
                                                           !iyv 06.02.14
!        do is=0,ns-1                                                    !srv 01.10.99
        do is = nsmin, nsmax-1
*    ..compute partial pressure
          call b2xppb (nCv, rt%rza(1,is), pl%na(1,is), pl%te,
     &      pl%ti, pl%tn, is, dv%pa(1,is))                                   !srv 01.10.99
*    ..compute particle flux
          call b2tfnb (nCv, nFc, nVx, is, ismain, geo, mpg, pl,
     &      dv, co, rt, .true.)
        enddo                                                            !srv 01.10.99
*      ..compute electron drift velocities
        call b2tfed (nCv, nFc, nVx, geo, mpg, dv%facdrift, dv%fac_ExB,
     &      dv%ne, pl%te, pl%po, dv%vedia, dv%veecrb)
!      end if                                                            !srv 01.10.99

csw 10dec2012: added extra call if wanted/not-wanted (SOLPS4 has an
c   extra call at this place within internal iteration)
      if(no_b2sral_call .eq. 0) then
        write(*,*) ' re-compute source '
*     ..re-compute source linearisation                                  !xpb
        call b2sral (nCv, nFc, nVx, ns, nxtl, nxtr,
     &    nscx, nscxmax, nsmin, nsmax, iscx, ismain, ismain0, dtim,
     &    BoRiS, geo, mpg, pl, dv, co, rt, rtw, st_ext, 
     &    sr, srw, wrong_flow, .false.)
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
        call b2srsm (nCv, ns, dtim, geo, mpg, pl%na, pl%ua, 
     &    pl%te, pl%ti, dv%ne, dv%ni, sr, .false.)
*     ..include contributions due to timestep                            !xpb
        call b2srdt (nCv, ns, dtim, BoRiS, geo, mpg,
     &    psnc%na, psnc%ua, psnc%te, psnc%ti, psnc%tn,
     &    psnc%ne, psnc%ni, psnc%kinrgy,
     &    pl%na, pl%ua, pl%te, pl%ti, pl%tn, dv%ne, 
     &    dv%ni, dv%kinrgy,  sr, .false.)
*     ..stabilise the source coefficients                                !xpb
        call b2srst (nCv, ns, mpg, pl%na, pl%ua, pl%te,   !srv 01.07.08
     &    pl%ti, pl%po, dv%ne, dv%ni, sr)
      endif
      if (ue_style.eq.2 .or. iout.ne.0)
*   ..compute ue, electron parallel velocity
     > call b2xpve (nCv, nFc, ns, geo, mpg, qe, rt%rza, rt%rz2, 
     &   co%alfx_c, co%sigx_c, dv%fch_pi_c, pl%po, dv%ne, 
     &   pl%te, pl%na, pl%ua, st_ext, dv%ue)
      if (ue_style.ne.2 .or. iout.ne.0) then
        wrkf = -dv%fch_p(:,0)/(qe*(geo%fcPbs + cutlo))
        call intcell(nFc,nCv,mpg,mpg%intcellP,wrkf,wrk0)
        do is=0,ns-1
          wrk0 = wrk0 + rt%rza(:,is)*pl%ua(:,is)*pl%na(:,is)
        enddo
        wrk0 = wrk0 + st_ext%ue*st_ext%ne
        wrk0 = wrk0/dv%ne                                                !srv 05.03.11
      endif
#ifdef WG_TODO
      if(iout.ne.0) then                                                 !srv 05.07.17
       call b2xzef (nx, ny, ns, rz2, na, ne, zeff)
       call b2xpic (nx, ny, ns, qe, rza, bb, na, ua, zeff, fch_pi_c)
       call b2xpif (nx, ny, hx, hy, vol, fch_pi_c, fch_pi_f)
       call my_out(70,nx,ny,fch_pi_c,'b2news__fch_pi_c')
       call my_out(70,nx,ny,fch_pi_f,'b2news__fch_pi_f')
       call my_out(70,nx,ny,ue,'b2news__ve')
       call my_out(70,nx,ny,wrk0,'b2news__ue')
      endif                                                              !srv 05.07.17
#endif
      if (ue_style.ne.2) dv%ue=wrk0

csw 11oct2012 added check on coriolis force (taken from b2news())
      ier0=0
c Necessary for smooth restart with Rhie-Chow:
c The Rhie-Chow term in the particle fluxes uses pccb terms which are calculated in b2npmo.
c Hence to get a smooth restart after cp b2fstate bfstati we need to do a dummy call to 
c b2npmo with rxf = 0 and recompute fluxes.
      if (ncall.eq.0) then
        call b2npmo (nCv, nFc, nVx, ns, ismain, geo, mpg,
     &    ncall, solving(1), 0.0_R8, pl, dv, rt, co, sr, st_ext,
     &    ier0)
        do is = nsmin, nsmax-1
*    ..compute partial pressure
          call b2xppb (nCv, rt%rza(1,is), pl%na(1,is), pl%te,
     &      pl%ti, pl%tn, is, dv%pa(1,is))                                  
*    ..compute particle flux
          call b2tfnb (nCv, nFc, nVx, is, ismain, geo, mpg, pl,
     &      dv, co, rt, .true.)
        enddo 
      endif
      call b2npmo (nCv, nFc, nVx, ns, ismain, geo, mpg,
     &  ncall, solving(1), rxf1, pl, dv, rt, co, sr, st_ext,
     &  ier0)
      call xertst (ier0.eq.0, 'error return from b2npmo') 
*   ..set nep
*     (= ne before continuity update)
      call b2scopy (nCv, dv%ne, 1, nep, 1)
*   ..recompute currents
      call b2tanml (nCv, nFc, ns, geo, mpg, co%csig_an,pl%po,dv%fchanml)
      call b2tvspa (nCv, nFc, nVx, ns, geo, mpg, pl%ua, co%vsaf_cl, 
     &  dv%fac_vis, dv%fchvispar)
      call b2tiner (nCv, nFc, nVx, ns, geo, mpg, pl%na, pl%ua,
     &  dv%facdrift, dv%fchinert)
*     ..iterate on particle balance
!      do is = 0, ns-1
      do is = nsmin, nsmax-1
*    ..DPC
        do iCv = 1, nCv
          pl%na(iCv,is)=max(pl%na(iCv,is),na_min) ! abs(na(ix,iy,is))    !srv 24.04.09
        enddo
*    ..compute partial pressure
        call b2xppb (nCv, rt%rza(1,is), pl%na(1,is), pl%te, 
     &    pl%ti, pl%tn, is, dv%pa(1,is))
*    ..compute particle flux
        call b2tfnb (nCv, nFc, nVx, is, ismain, geo, mpg, pl, dv,
     &    co, rt, .true.)
!        call b2tfco (nx, ny, is, na(-1,-1,is), te, ti, pb, 
!     &      cdpa(-1,-1,0,0,is), flob, conb)
*    ..perform iteration
csw 16oct2012 init ier0
        ier0=0
!        if (all(solveco(is,0:mpg%nnreg(0)))) then                        !srv 26.09.12 22.05.18
          solvedum(0:mpg%nnreg(0))=solveco(is,0:mpg%nnreg(0))
          solvedm2(0:mpg%nnreg(0))=solvemo(is,0:mpg%nnreg(0))
          new_matrix = .false.
          do ireg = 0, mpg%nnreg(0)
            new_matrix = new_matrix.or.
     &       solveco(is,ireg).neqv.last_solve_9(ireg)
          end do
          if(new_matrix.or..true.) call restart_ma28_for_us
          call b2npco (nCv, nFc, nVx, mpg%nnreg(0), solving(2),solvedum,  !srv 22.05.18
     &      solving(1),  solvedm2, ncall, is, rxf0, geo, mpg, pl, dv, 
     &      sr, ier0)
          call xertst (ier0.eq.0, 'error return from b2npco')
          last_solve_9(0:mpg%nnreg(0))=solveco(is,0:mpg%nnreg(0))
!        endif                                                            !srv 13.10.06 } 22.05.18
        do iCv = 1, nCv
          pl%na(iCv,is)=max(pl%na(iCv,is),na_min) ! abs(na(ix,iy,is))    !srv 24.04.09
        enddo
!        if(re_compute_pb.eq.1) then                                     !srv 15.02.10
!*    ..re-compute partial pressure                                      !srv 01.10.99
!        call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)   !srv 01.10.99
!        endif                                                           !srv 15.02.10
*    ..re-compute particle flux
        call b2tfnb (nCv, nFc, nVx, is, ismain, geo, mpg, pl, dv,
     &    co, rt, .false.)
      enddo 
*      ..compute electron drift velocities
        call b2tfed (nCv, nFc, nVx, geo, mpg, dv%facdrift, dv%fac_ExB,
     &      dv%ne, pl%te, pl%po, dv%vedia, dv%veecrb)

!WG_TODO      if(iout.ne.0 .or. iout_b2npco.ne.0 .or. iout_b2npc9.ne.0           !srv 16.07.10
!WG_TODO     &             .or. iout_b2npmo.ne.0 .or. iout_b2npht.ne.0 
!WG_TODO     &             .or. iout_b2nph9.ne.0 .or. iout_b2nppo.ne.0
!WG_TODO     &             .or. iout_b2npp7.ne.0 .or. iout_b2wdat.eq.4)
!WG_TODO*    .. calculate contributions due to timestep for diagnostic 
!WG_TODO     &  call b2srdt (nx, ny, ns, dtim, BoRiS, hz,
!WG_TODO     &       vol, na0, ua0, te0, ti0, ne0, ni0, kinrgy0,
!WG_TODO     &       na, ua, te, ti, ne, ni, kinrgy,
!WG_TODO     &       snadt, smodt, shedt, shidt, schdt, snedt, .true.)

*   ..re-compute ne, ni, lnlam
      call b2xpne (nCv, ns, rt%rza, pl%na, st_ext%ne, dv%ne)
      call nancheck(nCv, dv%ne, 'b2news(ne): after call to b2xpne')
      call b2xpni (nCv, ns, pl%na, dv%ni)
      call b2xpne (nCv, ns, rt%rz2, pl%na, st_ext%ne2, dv%ne2)  !srv 11.07.99
      call b2xppz (nCv, ns, dv%ne, pl%na, pl%te, pl%ti, dv%pz,st_ext)    !srv 18.01.02

*   ..compute ceqp                                                       !srv 03.07.03 {
      call b2xpnm (nCv, ns, rt%rz2, am, pl%na, wrk0)
      call b2treq (nCv, geo%cvVol, pl%te, pl%ti, dv%ne, dv%ni, wrk0, 
     %  dv%lnlam, st_ext, co%ceqp)                                       !srv 03.07.03 }
      if(iout .ne. 0) then                                               !srv 12.07.05 {
        do is = 0,ns-1
          write(charns,'(i3.3)') is
          call my_out_us(70,nCv,0,pl%na(1,is),'b2news__na'//charns)
        enddo
        call my_out_us(70,nCv,0,dv%ni(1,0),'b2news__ni0')
        call my_out_us(70,nCv,0,dv%ni(1,1),'b2news__ni1')
        call my_out_us(70,nCv,0,dv%ne,'b2news__ne')
        call my_out_us(70,nCv,0,dv%pz,'b2news__pz')
        call my_out_us(70,nCv,0,co%ceqp,'b2news__ceqp')
      endif                                                              !srv 12.07.05 }

      if (pot_eq.eq.0) then                                              !srv 13.07.05 {
        pl%po = 0.0e0_R8
      else if (pot_eq.eq.2) then                                         !srv 11.01.12
        pl%po = 3.1_R8*pl%te/qe
      else                                                               !srv 13.07.05 }
*   ..correct po
        if (solving(3)) then
          call b2nxdp (nCv, mpg, mpg%nnreg(0), solvepo, qe, pl%te, nep, 
     &       dv%ne, wrk0)
          call b2saxpy (nCv, 1.0_R8, wrk0, 1, pl%po, 1)
          if(iout .ne. 0) then                                           !srv 12.07.05 {
            call my_out_us(70,nCv,0,pl%po,'b2nxdp_po')
          endif                                                          !srv 12.07.05 }
        endif
*   ..compute currents                                                   !srv 29.01.20 {
        call b2tfch_ (nCv, nFc, nVx, ns, ismain, ismain0, geo, mpg, 
     &     pl, dv, co, rt)
        if (solving(3)) then
*   ..solve the potential equation
          if (vis_per.eq.0.0e0_R8) then                                 !srv 06.07.00
            new_matrix = .false.
            do ireg = 0, mpg%nnreg(0)
             new_matrix = new_matrix.or.
     &        solvepo(ireg).neqv.last_solve_9(ireg)
            end do
            if (new_matrix.or..true.) call restart_ma28_for_us
              call b2nppo (nCv, nFc, nVx, mpg%nnreg(0), solving(3),      !srv 22.05.18
     &        solvepo, ncall, rxf2, geo, mpg, pl, dv, sr, ier0)
            last_solve_9(0:mpg%nnreg(0))=solvepo(0:mpg%nnreg(0))
          else
             call xerrab ('visper not yet available for WG')
#ifdef WG_TODO
            call b2npp7 (nx, ny, mpg%nnreg(0), solving(3), solvepo,      !srv 27.01.00 22.05.18
     &       ncall, rxf2, ne, te, conc, avisper, bvisper, cvisper,       !srv 27.01.00
     &       sch, po, fch, respo, corpo, ier0)                           !srv 27.01.00
#endif
          endif                                                          !srv 06.07.00
          t0 = b2sasum(nCv,dv%respo,1)
          write(*,*) 'respo before solving',t0
          call xertst (ier0.eq.0, 'error return from b2nppo')
*   ..recompute currents
          call b2tfch_ (nCv, nFc, nVx, ns, ismain, ismain0, geo, mpg, 
     &     pl, dv, co, rt)
*   ..compute new residual following the tradition 
          call b2ursd (nCv, nFc, geo, mpg, dv%ne, pl%po, sr%sch, dv%fch, 
     &    wrk0)
          t0 = b2sasum(nCv,wrk0,1)
          write(*,*) 'respo after solving',t0
        else
         if(get_residuals.eq.1)
     &    call b2ursd (nCv, nFc, geo, mpg, dv%ne, pl%po, sr%sch, dv%fch, 
     &                 dv%respo)
        endif
      endif           ! if pot_eq                                        !srv 13.07.05 29.01.20 }

      if(iout .ne. 0) then                                               !srv 12.07.05 {
        write (charns,'(i3.3)') i
        call my_out_us(70,nCv,0,dv%respo,'b2news__respo'//charns)
        call my_out_us(70,nCv,0,pl%po,'b2news__po'//charns)
      endif                                                              !srv 12.07.05 }

#ifdef WG_TODO
          if(ion_vlct_restrict .ne. 0 ) then                             !iyv 06.02.14
            call b2xprz (nx, ny, ns, mp, am, na, rz)
            call b2xppz (nx, ny, ns, ne, na, te, ti, pz)
            wrk0(:,:)=sqrt(pz(:,:)/rz(:,:))*abs(bb(:,:,0)/bb(:,:,3))
            call intfaceh(nx,ny,leftix,leftiy,vol,wrk0,cssb(-1,-1,0))
            call intfacev(nx,ny,bottomix,bottomiy,vol,wrk0,
     &                    cssb(-1,-1,1))
          else
            cssb = 0.0_R8
          endif  
#else
          cssb = 0.0_R8
#endif
      do is = 0, ns-1
*   ..compute partial pressure
        call b2xppb (nCv, rt%rza(1,is), pl%na(1,is), pl%te,
     &      pl%ti, pl%tn, is, dv%pa(1,is))
*   ..re-compute particle flux
#ifdef DBG
            write(*,*) 'is, za(is) ',is, (zamin(is)+zamax(is))/2.0_R8
#endif
        call b2tfnb (nCv, nFc, nVx, is, ismain, geo, mpg, pl, dv,
     &   co, rt, .false.)
      enddo

*   ..compute fne_he, fni_he
c      with diamagnetic drift part                                       !srv 13.05.99
      call b2xpfe (nCv, nFc, ns, st_ext%ns, geo, mpg, qe, 
     &  rt%rza, st_ext%za, dv%fna_he, st_ext%fa, dv%fch_p, dv%fne_he)    !srv 13.05.99 21.01.00 11.10.02
      call b2xpfe (nCv, nFc, ns, st_ext%ns, geo, mpg, qe, 
     &  rt%rza, st_ext%za, dv%fna_53, st_ext%fa, dv%fch_p, dv%fne_53)
      call b2xpfe (nCv, nFc, ns, st_ext%ns, geo, mpg, qe, 
     &  rt%rza, st_ext%za, dv%fna_32nd, st_ext%fa, dv%fch_32, dv%fne_32)
      call b2xpfe (nCv, nFc, ns, 0, geo, mpg, qe, 
     &      rt%rza, zero, dv%fna_52nd, zerof, dv%fch_52, dv%fne_52)
      if (tn_style.eq.0) then
        call b2xpfi (nFc, ns, dv%fna_he, dv%fni_he)                      !srv 13.05.99 21.01.00
        call b2xpfi (nFc, ns, dv%fna_32nd, dv%fni_32)                    !xpb
        call b2xpfi (nFc, ns, dv%fna_52nd, dv%fni_52)                    !xpb
      elseif (tn_style.eq.1) then
        call b2xpfz (nFc, ns, dv%fna_he, dv%fni_he)
        call b2xpfz (nFc, ns, dv%fna_32nd, dv%fni_32)
        call b2xpfz (nFc, ns, dv%fna_52nd, dv%fni_52)
      else
        call b2xpfz (nFc, ns, dv%fna_he, dv%fni_he)
        call b2xpfz (nFc, ns, dv%fna_32nd, dv%fni_32)
        call b2xpfz (nFc, ns, dv%fna_52nd, dv%fni_52)
        call b2xpfn (nFc, ns, dv%fna_he, dv%fnn)
        call b2xpfn (nFc, ns, dv%fna_32nd, dv%fnn_32)
        call b2xpfn (nFc, ns, dv%fna_52nd, dv%fnn_52)
      endif
*      ..compute electron drift velocities
      call b2tfed (nCv, nFc, nVx, geo, mpg, dv%facdrift, dv%fac_ExB,
     &      dv%ne, pl%te, pl%po, dv%vedia, dv%veecrb)
      call b2tfhe_ (nCv, nFc, nVx, ns, geo, mpg,                         !srv 29.01.20
     &      pl, dv, co, .true.)                                          !djm Jan2017
csw 30apr2014 added fna_he, 17sep2014 added ua,cvsa
      call b2tfhi_(nCv, nFc, nVx, ns, geo, mpg, pl, dv, co, rt,          !srv 13.10.06
     &  .true.)                                                          !djm Jan2017
       dv%floi_noc = dv%floi                                             !srv 24.10.02
       dv%floe_noc = dv%floe                                             !srv 24.10.02

      call b2xpen (nCv, nFc, ns, mpg, dv%fna, dv%fna_53, dv%fne_53, 
     &  dv%fch, dv%fhe, dv%fhi, dv%fhn, rt%rpt, pl%te,
     &  pl%ti, pl%tn, pl%po, BoRiS, dv%fhm, dv%fhp, dv%fhj, dv%fht)

*   ..add contribution from cdpa to floe, cone, floi, coni
      call b2tfcc (nCv, nFc, nVx, ns, geo, mpg, pl, dv, co, rt)
      if(ue_style.ne.2) then
*   ..compute ue, electron parallel velocity                             !srv 05.03.11 {
       wrkf = -dv%fch_p(:,0)/(qe*(geo%fcPbs + cutlo))
       call intcell(nFc,nCv,mpg,mpg%intcellP,wrkf,dv%ue)
       do is=0,ns-1
         dv%ue = dv%ue + rt%rza(:,is)*pl%ua(:,is)*pl%na(:,is)
       enddo
       dv%ue = dv%ue + st_ext%ue*st_ext%ne
       dv%ue = dv%ue/dv%ne                                                !srv 05.03.11
      else
       call b2xpve (nCv, nFc, ns, geo, mpg, qe, rt%rza, rt%rz2, 
     &   co%alfx_c, co%sigx_c, dv%fch_pi_c, pl%po, dv%ne, 
     &   pl%te, pl%na, pl%ua, st_ext, dv%ue)
      endif

*   ..iterate on heat balance
      call b2npht (nCv, nFc, nVx, ns, geo, mpg, mpg%nnreg(0), ncall,
     &  ismain, solving, solveee, solveei, solveen, solveet, solvepo,       !srv 22.05.18
     &  solvemo, rxf3, BoRiS, pl, dv, co, sr, ier0)

csw 16aug2011 init potential 
      if(pot_eq.eq.0) pl%po=0.0_R8                                          !srv 13.07.05
      if(pot_eq.eq.2) pl%po=3.1_R8*pl%te/qe 
      call xertst (ier0.eq.0, 'error return from b2npht')

csw 12oct2012, extra call to relax particle balance if wanted (as in SOLPS4)
      if(do_2nd_b2npco_call.ne.0) then
        call xerrab ('b2news_ 2nd_no_call not yet for WG')
#ifdef WG_TODO
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
        call b2tfnb (nCv, nFc, nVx, is, ismain, geo, mpg, pl, dv,
     &  co, rt, .false.)

*    ..perform iteration
        new_matrix = .false.
        if(mdf_fnb.eq.0) then                                            !srv 13.10.06
!          if (all(solveco(is,0:mpg%nnreg(0)))) then                      !srv 22.05.18
            solvedum(0:mpg%nnreg(0))=solveco(is,0:mpg%nnreg(0))
            solvedm2(0:mpg%nnreg(0))=solvemo(is,0:mpg%nnreg(0))
            do ireg = 0, mpg%nnreg(0)
              new_matrix = new_matrix.or.
     &         solveco(is,ireg).neqv.last_solve_5(ireg)
            end do
            if(new_matrix) call restart_ma28
            call b2npco (nx, ny, mpg%nnreg(0), solving(2), solvedum,     !srv 22.05.18 
     &       solving(1), solvedm2, ncall, is, rxf0, vol, pbs, flob,conb, !srv 22.01.02 01.07.05
     &       sna(-1,-1,0,is), pcca(-1,-1,0,0,is), pb, na(-1,-1,is),
     &       ua(-1,-1,is), fna(-1,-1,0,0,is),
     &       resco(-1,-1,is), corpa(-1,-1,is), ier0)
            last_solve_5(0:mpg%nnreg(0))=solveco(is,0:mpg%nnreg(0))
!          endif                                                         !srv 22.05.18
        else                                                             !srv 13.10.06 {
 !         if (all(solveco(is,0:mpg%nnreg(0)))) then                      !srv 22.05.18
            solvedum(0:mpg%nnreg(0))=solveco(is,0:mpg%nnreg(0))
            solvedm2(0:mpg%nnreg(0))=solvemo(is,0:mpg%nnreg(0))
cwdk tmp: always use 9-point stencil, also for neutrals
cwdk      need a master-switch to decide to use 9-point stencil...
            if(.not.is_neutral(is).or. .true.) then
              do ireg = 0, mpg%nnreg(0)
                new_matrix = new_matrix.or.
     &           solveco(is,ireg).neqv.last_solve_9(ireg)
              end do
              if(new_matrix.or..true.) call restart_ma28_for_9
              call b2npc9 (nx, ny, mpg%nnreg(0), solvedum,
     &         solving(1), solvedm2, ncall, is, rxf0,
     &         vol, pbs, flob, conb,
     &         sna(-1,-1,0,is), pcca(-1,-1,0,0,is), pb, na(-1,-1,is),
     &         ua(-1,-1,is), fna_mdf(-1,-1,0,0,is),
     &         resco(-1,-1,is), corpa(-1,-1,is), ier0)
              last_solve_9(0:mpg%nnreg(0))=solveco(is,0:mpg%nnreg(0))
            else
              do ireg = 0, mpg%nnreg(0)
                new_matrix = new_matrix.or.
     &           solveco(is,ireg).neqv.last_solve_5(ireg)
              end do
              if(new_matrix) call restart_ma28
              call b2npco (nx, ny, mpg%nnreg(0),  solving(2), solvedum,  !srv 22.05.18  
     &         solving(1), solvedm2, ncall, is, rxf0,
     &         vol, pbs, flob, conb,
     &         sna(-1,-1,0,is), pcca(-1,-1,0,0,is), pb, na(-1,-1,is),
     &         ua(-1,-1,is), fna_mdf(-1,-1,0,0,is),
     &         resco(-1,-1,is), corpa(-1,-1,is), ier0)
              last_solve_5(0:mpg%nnreg(0))=solveco(is,0:mpg%nnreg(0))
            endif
!          endif                                                         !srv 22.05.18
        endif                                                            !srv 13.10.06 }
        call xertst (ier0.eq.0, 'error return from b2npco')
       enddo
*      ..compute electron drift velocities
        call b2tfed (nCv, nFc, nVx, geo, mpg, dv%facdrift, dv%fac_ExB,
     &      dv%ne, pl%te, pl%po, dv%vedia, dv%veecrb)
#endif
      endif ! do_2nd_b2npco

csw 06dec2012 finally re-evaluate particle fluxes
c   (necessary to work out new rescaling factors sclrtio/sclalph in neutr_scl)
      if(re_eval_prtls_fluxes.eq.1 .or. use_eirene.ne.0) then            !srv 28.06.13 {
       do is=0,ns-1
        do iCv = 1, nCv
          pl%na(iCv,is)=max(pl%na(iCv,is),na_min) ! abs(na(ix,iy,is))    !srv 24.04.09
        enddo
*    ..compute partial pressure
!pb        call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)    !srv 01.10.99
        call b2xppb (nCv, rt%rza(1,is), pl%na(1,is), pl%te,              !srv 01.10.99
     &    pl%ti, pl%tn, is, dv%pa(1,is))
*    ..compute particle flux
        call b2tfnb (nCv, nFc, nVx, is, ismain, geo, mpg, pl, dv,
     &  co, rt, .false.)
       enddo                                                             !srv 01.10.99
*      ..compute electron drift velocities
        call b2tfed (nCv, nFc, nVx, geo, mpg, dv%facdrift, dv%fac_ExB,
     &      dv%ne, pl%te, pl%po, dv%vedia, dv%veecrb)
      endif                                                              !srv 28.06.13 }

csw 26oct2011 recalculate rescaling factors sclrtio/sclalph (SOLPS4)
#ifdef WG_TODO
      if(no_neutr_scl .eq. 0) then
        write(*,*) 'Calling b2stbr_neutr_scl'
        call b2stbr_neutr_scl(nx, ny, ns)                                !srv 25.10.08
        write(*,*) 'Called b2stbr_neutr_scl'
      endif
#endif

* ..return
      ierr = 0
99999 continue
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2news_

      end
