*-----------------------------------------------------------------------
*.specification

      subroutine b2trcl (nCv, nFc, nVx, ns, nscx, iscx, ismain,       !srv 23.11.10
     1  switch, geo, mpg, pl, dv, rt, co, st_ext, cvsahz,             !srv 09.01.01 03.06.03
     3  cvsa, chce, chve, chci, chvi, csig, calf, cdkt, csigin,       !xpb
     4  cthe, cthi, vsaf_cl, cvsa_cl, cvsahz_cl, fllimvisc,           !srv 09.01.01
     5  csig_cl, calf_cl)                                             !xpb
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
#ifdef ASTRA
      use b2mod_b2_to_astra
#endif
      use b2mod_ad
     & , only : ncall_b2trcl, b2trcl_cutlo
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_math
     & , only : b2mod_math_initialised, cutlo, cutll, small_R4_constant
      use b2mod_ad
     & , only : my_out_folder,
#ifdef WG_DEBUG
     &          my_out_sf,                                            !srv 27.06.21
#endif
     &          ncall_b2tlnl,
     &          ncall_b2tlmv, ncall_b2tqin, ncall_b2tqce, ncall_b2tqca,
     &          ncall_b2ttia
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent (in) :: nCv, nFc, nVx, ns, nscx, ismain
      integer, intent (in) :: iscx(0:nscx-1)
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2StateExt), intent (in) :: st_ext
*   ..output arguments (unspecified on entry)
      type (B2Coeff), intent (inout) :: co
cwdk to consider: create an output struct co0 of type (B2Coeff) for this?
      real (kind=R8) ::
     &  cvsa(nFc,0:1,0:ns-1),
     &  chce(nFc,0:1), chve(nFc,0:1),
     &  chci(nFc,0:1), chvi(nFc,0:1),
     &  csig(nFc,0:1), calf(nFc,0:1),
     &  cdkt(nFc,0:1),
     &  csigin(nFc,0:1,0:ns-1,0:ns-1),                           !xpb
     &  cthe(nCv,0:ns-1), cthi(nCv,0:ns-1),
     &  vsaf_cl(nFc,0:1,0:ns-1),
     &  cvsa_cl(nFc,0:1,0:ns-1),
     &  cvsahz_cl(nFc,0:1,0:ns-1),                               !srv 09.01.01
     &  fllimvisc(nFc,0:ns-1),
     &  cvsahz(nFc,0:1,0:ns-1),   !srv 09.01.01
     &  csig_cl(nFc,0:1), calf_cl(nFc,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRCL computes the classical transport coefficients, including
*     geometric factors such as vol/hx**2 and vol/hy**2.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the following discrete classical transport
*     coefficients.
*     - cvsa: viscosity.
*     - chce: electron thermal conductivity.
*     - chve: electron thermal strange velocity.
*     - chci: atom thermal conductivity.
*     - chvi: atom thermal strange velocity.
*     - csig: electrical conductivity.
*     - csigin: ion-neutral conductivity.
*     - calf: thermo-electric coefficient.
*
*     The local differential transport coefficients, defined on cells,
*     are obtained by calls to b2tqca and b2tqce. To compute from these
*     the discrete transport coefficients involves multiplication by
*     geometric factors (vol, vol/hx**2, vol/hy**2), and, usually,
*     interpolation to cell faces.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer k, is, is0, iCv, iFc, iFt, mode
      character(3) chns, chns0                                        !srv 17.06.08
#ifdef ASTRA
      integer mintopcut                                               !lk 19.09.07
#endif
      real (kind=R8) ::
     *  cflmv, jemax                                                  !srv 13.01.17
      real (kind=R8) ::                                               !srv 03.06.03
     *   t0, dte(nFc), dti(nFc), flomx, connf, epsf
      real (kind=R8) :: tauia(nCv,0:ns-1),                            !lk 19.09.07
     &    collisnumf(mpg%nFc), collisnumef(mpg%nFc)
      real (kind=R8) :: inctran
#ifdef ASTRA
     & , tt0(mpg%nFt), tt1(mpg%nFt)                                   !lk 19.09.07
#endif
      real (kind=R8) :: tau(nCv), nu1(nCv),                           !lk 12.05.11
     &    k1(nCv), t1, t2                                             !lk 12.05.11
      real (kind=R8) ::
     *  hce(nCv,0:1), sig(nCv,0:1),
     1  alf(nCv,0:1), vsa(nCv,0:1,0:ns-1),
     2  hci(nCv,0:1), sigin(nCv,0:1,0:ns-1,0:nscx-1),
     3  wrk(nFc), wrk1(nFc,2), wrkv(nVx), wrkc(nCv),                  !srv 01.07.99 09.01.01
     4  tef(nFc), tif(nFc), ne2f(nFc), k1f, k2f, gtec(nCv),
     5  nu1f, nu2f, tauf(nFc), tauiaf(nFc)
*   ..procedures
      real (kind=R8) ::
     *  smin, smax                                                    !srv 22.04.16
      external xertst, xerrab, sfill, smin, smax
      external b2tqca, b2tqce, b2tqin, b2txcx, b2txcy, b2tlmv,        !xpb
     &  ipgeti, ipgetr, b2ttia                                        !lk 19.09.07
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2trcl')
*   ..test nCv, nFc, ns
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..set internal parameters on first call
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c        call ipgetr ('b2tfhe_vis_par',vis_par)
c        call ipgeti ('b2trcl_integral_modeling', integral_modeling)  !srv 01.11.11
c        call ipgeti ('b2trcl_lthf21',lthf21)
c       call ipgeti ('b2trcl_lvis21',lvis21)
c#ifdef NO_21
c        call xertst (lthf21.eq.0.and.lvis21.eq.0,
c     &   'Recompile without NO_21 option for 21-moment approach!')
c#endif
c        call xertst (.not.(lvis21.ne.0.and.vis_par.ne.0.0_R8),
c     &   '21-moment viscosity already includes parallel viscosity!')
*   ..extensive tests on first few calls
      if (ncall_b2trcl.lt.3) then
*    ..test bb
        call b2xvsg (nCv, geo%cvBb(1,3), 1, 'bb3', '.gt.')
*    ..test sign of vol
        call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
*    ..test range of cvQgam, fcQgam
       call xertst (-1.0_R8.lt.smin(nCv,geo%cvQgam(1,1),1).and.
     &   smax(nCv,geo%cvQgam(1,1),1).lt.1.0_R8,
     &   'faulty argument range cvQgam(,1)')
       call xertst (0.0_R8.lt.smin(nCv,geo%cvQgam(1,0),1).and.
     &   smax(nCv,geo%cvQgam(1,0),1).le.1.0_R8,
     &   'faulty argument range cvQgam(,0)')
       call xertst (-1.0_R8.lt.smin(nFc,geo%fcQgam(1,1),1).and.
     &   smax(nFc,geo%fcQgam(1,1),1).lt.1.0_R8,
     &   'faulty argument range fcQgam(,1)')
       call xertst (0.0_R8.lt.smin(nFc,geo%fcQgam(1,0),1).and.
     &   smax(nFc,geo%fcQgam(1,0),1).le.1.0_R8,
     &   'faulty argument range fcQgam(,0)')
*    ..test sign of na, te, ti
        call b2xvsg (nCv*ns, pl%na, 1, 'na', '.gt.')
        call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
        call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
      endif

      mode = switch%b2tral_mode

* ..compute transport coefficients
*   ..compute differential classical transport coefficients
      call b2tqca (nCv, ns, switch, geo, pl, dv, rt, st_ext, vsa, hci) !srv 01.07.09 20.09.11
      call b2tqce (nCv, switch, geo, mpg, pl, dv, hce, sig, alf)       !srv 01.07.09 20.09.11
      call b2tqin (nCv, ns, nscx, iscx, switch, geo, pl, rt, sigin)    !xpb
      call sfill (nCv*ns, switch%cthev, cthe,1)
      call sfill (nCv*ns, switch%cthiv, cthi,1)
      call sfill (nFc, 1.0_R8, co%fllim_ki, 1)                         !srv 13.01.17 {
      call sfill (nFc, 1.0_R8, co%fllim_ke, 1)
      call sfill (nFc, 1.0_R8, co%fllim_al, 1)
      call sfill (nCv, 1.0_R8, co%fllim_al_c, 1)
      co%alfx_c(:) = alf(:,0)  !/abs((bb(:,:,0))/bb(:,:,3))            !iys 19.11.13 !srv 08.10.14 { !srv 13.01.17
      co%sigx_c(:) = sig(:,0)  !/abs((bb(:,:,0))/bb(:,:,3))            !srv 13.01.17 }
      co%hcix_c(:) = hci(:,0)
*   ..interpolate vsa to faces
      !wdk need to write general routine for volume interpolation like this...
      !wdk intface_v? intface with mode switch?
      do is = 0, ns-1
        do iFc = 1, nFc
          vsaf_cl(iFc,0,is) =
     &      (geo%fcVol(iFc,1)*vsa(mpg%fcCv(iFc,1),0,is)+
     &       geo%fcVol(iFc,2)*vsa(mpg%fcCv(iFc,2),0,is))/
     &      (geo%fcVol(iFc,1)+geo%fcVol(iFc,2))
          vsaf_cl(iFc,1,is) =
     &      (geo%fcVol(iFc,1)*vsa(mpg%fcCv(iFc,1),1,is)+
     &       geo%fcVol(iFc,2)*vsa(mpg%fcCv(iFc,2),1,is))/
     &      (geo%fcVol(iFc,1)+geo%fcVol(iFc,2))
        enddo
      enddo
*   ..incorporate geometry with hz                                     !srv 09.01.01
      inctran = switch%b2txcx_increase_transp_coefs
      wrk = geo%fcS*geo%fcHz
      do iFc = 1, nFc
        wrk1(iFc,1) = geo%fcVol(iFc,1)*geo%cvHz(mpg%fcCv(iFc,1))
        wrk1(iFc,2) = geo%fcVol(iFc,2)*geo%cvHz(mpg%fcCv(iFc,2))
      enddo
      do is = 0, ns-1
       call b2txcx (nCv, nFc, inctran,
     &  mode, geo, mpg, wrk1, wrk, vsa(:,0,is), cvsahz(:,0,is))
       call b2txcy (nCv, nFc,
     &  mode, geo, mpg, wrk1, wrk, vsa(:,1,is), cvsahz(:,1,is))
      enddo                                                           !srv 09.01.01
*   ..incorporate geometry
      call sfill (nFc*2*ns*ns, 0.0_R8, csigin, 1)
#if (defined(_OPENMP) && !defined(NO_OPENMP_B2TRCL))
!svdk parallelize over species
!$omp parallel do default(none)
!$omp& private(is, k, is0)
!$omp& shared(ns, nCv, nFc, mode, geo, mpg, vsa, cvsa, nscx, iscx)
!$omp& shared(is_neutral, zn, sigin, csigin)
#endif
      do is = 0, ns-1
        call b2txcx (nCv, nFc, inctran,
     &   mode, geo, mpg, geo%fcVol, geo%fcS, vsa(:,0,is), cvsa(:,0,is))
        call b2txcy (nCv, nFc,
     &   mode, geo, mpg, geo%fcVol, geo%fcS, vsa(:,1,is), cvsa(:,1,is))
        do k = 0, nscx-1
          is0 = iscx(k)
          if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then
            call b2txcx (nCv, nFc, 0.0_R8,
     &       mode, geo, mpg, geo%fcVol, geo%fcS,
     &       sigin(:,0,is,k), csigin(:,0,is,is0))
            call b2txcy (nCv, nFc,
     &       mode, geo, mpg, geo%fcVol, geo%fcS,
     &       sigin(:,1,is,k), csigin(:,1,is,is0))
          endif
        enddo
      enddo
#if (defined(_OPENMP) && !defined(NO_OPENMP_B2TRCL))
!$omp end parallel do
#endif

      call b2txcx (nCv, nFc, inctran,
     & mode, geo, mpg, geo%fcVol, geo%fcS, hce(:,0), chce(:,0))
      call b2txcy (nCv, nFc,
     & mode, geo, mpg, geo%fcVol, geo%fcS, hce(:,1), chce(:,1))
      call b2txcx (nCv, nFc, inctran,
     & mode, geo, mpg, geo%fcVol, geo%fcS, hci(:,0), chci(:,0))
      call b2txcy (nCv, nFc,
     & mode, geo, mpg, geo%fcVol, geo%fcS, hci(:,1), chci(:,1))
      call b2txcx (nCv, nFc, 0.0_R8,
     & mode, geo, mpg, geo%fcVol, geo%fcS, sig(:,0), csig(:,0))
      call b2txcy (nCv, nFc,
     & mode, geo, mpg, geo%fcVol, geo%fcS, sig(:,1), csig(:,1))
      call b2txcx (nCv, nFc, 0.0_R8,
     &   -1, geo, mpg, geo%fcVol, geo%fcS, alf(:,0), calf(:,0))
      call b2txcy (nCv, nFc,
     &   -1, geo, mpg, geo%fcVol, geo%fcS, alf(:,1), calf(:,1))

      if(switch%b2npmo_iout.ne.0 .or. switch%iout_b2wdat.eq.4 ) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        if (.not. is_neutral(is)) then
          call my_out_us(70,nFc,1,cvsahz(1,0,is),'b2trcl_cvsahzx'//
     &      chns)
        endif
       enddo
      endif
      if(switch%b2npht_iout.ne.0 .or. switch%iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,chce(1,0),'b2trcl_chcex')
!        call my_out_us(70,nFc,1,chce(1,1),'b2trcl_chcey')
        call my_out_us(70,nFc,1,calf(1,0),'b2trcl_calfx')
!        call my_out_us(70,nFc,1,calf(1,1),'b2trcl_calfy')
        call my_out_us(70,nCv,0,co%alfx_c,'b2trcl_alfx_c')            !srv 05.07.17
        call my_out_us(70,nFc,1,chci(1,0),'b2trcl_chcix')
!        call my_out_us(70,nFc,1,chci(1,1),'b2trcl_chciy')
      endif
      if(switch%b2npp7_iout.ne.0 .or. switch%b2nppo_iout.ne.0
     &                   .or. switch%iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,csig(1,0),'b2trcl_csigx')
!        call my_out_us(70,nFc,1,csig(1,1),'b2trcl_csigy')
        call my_out_us(70,nCv,0,co%sigx_c,'b2trcl_sigx_c')            !srv 05.07.17
      endif

*   ..call 21 moment routines
#ifdef WG_TODO
#ifndef NO_21
      if (lthf21.ne.0) call b2t21t (nx, ny, ns, lnlam, bb, vol, hx,   !srv 01.07.09
     1 pbs, ne, na, ne2, te, ti, ua, fch_p, chce, chci, cthe, cthi)
      if (lvis21.ne.0) call b2t21v (nx, ny, ns, lnlam,                !srv 01.07.09
     1 bb, vol, hx, pbs, ne, na, te, ti, ua, cvsa)
#endif
#endif
*   ..luciani
      co%f_luc_ke = 1.0_R8                                            !srv 13.01.17 {
      co%f_luc_ki = 1.0_R8
      co%f_luc_et = 1.0_R8
      co%f_luc_sg = 1.0_R8
      co%f_luc_al = 1.0_R8                                            !srv 13.01.17 }
      co%alpha1f = 1.0_R8
      co%f_luc_ke_c = 1.0_R8
      co%f_luc_ki_c = 1.0_R8
      if(switch%lluciani.eq.1) then   ! connection length formulation
       !wdk There is implicit assumption here that flux tubes exist in the grid.
       !wdk Need to generalize to truly unstructured grids; using fcConn at each
       !wdk face?
*    ..interpolate te, ti and ne2 to cell faces
       wrk1 = 1.0_R8
       call intface(nCv,nFc,mpg%fcCv,wrk1,pl%te,tef)
       call intface(nCv,nFc,mpg%fcCv,wrk1,pl%ti,tif)
       call intface(nCv,nFc,mpg%fcCv,wrk1,dv%ne2,ne2f)
       do iFc = 1, nFc                                                !lkw 12.10.22 {
        if (mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .or.
     &      mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
          if (mpg%fcFs(iFc).gt.0) then
            connf = geo%fsConn(mpg%fcFs(iFc))
          else
            iFt = mpg%cvFt(mpg%fcCv(iFc,1))
            if (iFt.gt.0.and.iFt.le.mpg%nFt) connf = geo%ftConn(iFt)
          end if
          ! compute flux limits
          if (connf.gt.0.0_R8) then
            co%f_luc_ke(iFc) = 1.0_R8/(1.0_R8 +
     &       7.5e16_R8*(tef(iFc)/ev)**2/(connf*ne2f(iFc)))
            co%f_luc_ki(iFc) = 1.0_R8/(1.0_R8 +
     &       7.5e16_R8*(tif(iFc)/ev)**2/(connf*ne2f(iFc)))
          end if
          co%f_luc_sg(iFc) = co%f_luc_ke(iFc)
          co%f_luc_al(iFc) = co%f_luc_ke(iFc)
          co%f_luc_et(iFc) = co%f_luc_ki(iFc)
          do is=0,ns-1
            if(.not.is_neutral(is)) then
              cvsa(iFc,0,is)    = cvsa(iFc,0,is)*co%f_luc_et(iFc)
              cvsahz(iFc,0,is)  = cvsahz(iFc,0,is)*co%f_luc_et(iFc)
              vsaf_cl(iFc,0,is) = vsaf_cl(iFc,0,is)*co%f_luc_et(iFc)
            endif
          enddo
          chce(iFc,0) = chce(iFc,0)*co%f_luc_ke(iFc)
          csig(iFc,0) = csig(iFc,0)*co%f_luc_sg(iFc)
          calf(iFc,0) = calf(iFc,0)*co%f_luc_al(iFc)
          chci(iFc,0) = chci(iFc,0)*co%f_luc_ki(iFc)
        endif
       enddo
       do iCv = 1, nCv
        if (mpg%cvOnClosedSurface(iCv)) then
          iFt = mpg%cvFt(iCv)
          if (iFt.gt.0.and.iFt.le.mpg%nFt) then
            co%f_luc_ke_c(iCv) = 1.0_R8/(1.0_R8 +
     1       7.5e16_R8*(pl%te(iCv)/ev)**2/(geo%ftConn(iFt)*dv%ne2(iCv)))
            co%sigx_c(iCv) = co%sigx_c(iCv)*co%f_luc_ke_c(iCv)
            co%alfx_c(iCv) = co%alfx_c(iCv)*co%f_luc_ke_c(iCv)
            co%f_luc_ki_c(iCv) = 1.0_R8/(1.0_R8 +
     1       7.5e16_R8*(pl%ti(iCv)/ev)**2/(geo%ftConn(iFt)*dv%ne2(iCv)))
            co%hcix_c(iCv) = co%hcix_c(iCv)*co%f_luc_ki_c(iCv)
          endif
        endif
       enddo
      else if (switch%lluciani.eq.2) then  ! SOLPS4.0 formulation
*    ..interpolate te, ti and ne2 to cell faces
       wrk1 = 1.0_R8
       call intface(nCv,nFc,mpg%fcCv,wrk1,pl%te,tef)
       call intface(nCv,nFc,mpg%fcCv,wrk1,pl%ti,tif)
       call intface(nCv,nFc,mpg%fcCv,wrk1,dv%ne2,ne2f)
       do iFc = 1, nFc                                                !lkw 12.10.22 {
        if (mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .or.
     &      mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
          co%f_luc_ke(iFc) = 1.0_R8/(1.0_R8 +
     1      0.3e16_R8*(tef(iFc)/ev)**2/ne2f(iFc))
          co%f_luc_ki(iFc) = 1.0_R8/(1.0_R8 +
     1      0.3e16_R8*(tif(iFc)/ev)**2/ne2f(iFc))
          co%f_luc_sg(iFc) = co%f_luc_ke(iFc)
          co%f_luc_al(iFc) = co%f_luc_ke(iFc)
          co%f_luc_et(iFc) = co%f_luc_ki(iFc)
          do is = 0, ns-1
            if(.not.is_neutral(is)) then
              cvsa(iFc,0,is)    = cvsa(iFc,0,is)*co%f_luc_et(iFc)
              cvsahz(iFc,0,is)  = cvsahz(iFc,0,is)*co%f_luc_et(iFc)
              vsaf_cl(iFc,0,is) = vsaf_cl(iFc,0,is)*co%f_luc_et(iFc)
            endif
          enddo
          chce(iFc,0) = chce(iFc,0)*co%f_luc_ke(iFc)
          csig(iFc,0) = csig(iFc,0)*co%f_luc_sg(iFc)
          calf(iFc,0) = calf(iFc,0)*co%f_luc_al(iFc)
          chci(iFc,0) = chci(iFc,0)*co%f_luc_ki(iFc)
         endif
        enddo
       do iCv = 1, nCv
        if (mpg%cvOnClosedSurface(iCv)) then
          iFt = mpg%cvFt(iCv)
          if (iFt.gt.0.and.iFt.le.mpg%nFt) then
            co%f_luc_ke_c(iCv) = 1.0_R8/(1.0_R8 +
     1       0.3e16_R8*(pl%te(iCv)/ev)**2/dv%ne2(iCv))
            co%sigx_c(iCv) = co%sigx_c(iCv)*co%f_luc_ke_c(iCv)
            co%alfx_c(iCv) = co%alfx_c(iCv)*co%f_luc_ke_c(iCv)
            co%f_luc_ki_c(iCv) = 1.0_R8/(1.0_R8 +
     1       0.3e16_R8*(pl%ti(iCv)/ev)**2/dv%ne2(iCv))
            co%hcix_c(iCv) = co%hcix_c(iCv)*co%f_luc_ki_c(iCv)
          endif
        endif
       enddo
      else if (switch%lluciani.eq.3) then  ! new formulation
*   .. plateau and banana correction for ion heat conduction          !lk 19.09.07 {
       call b2ttia (nCv, ns, pl%ti, rt%rz2, dv%ne2, dv%lnlam, tauia)
*   .. calculate collisionality parameter nu_star
       nu1 = 0.0_R8                                                   !lk 12.05.11{
       k1 = 0.58_R8
       co%nu2 = 0.0_R8
       co%k2 = 1.6_R8
       t1 = (4.0_R8*pi*eps0/qe)/qe
       t2 = 3.0_R8/(4.0_R8*sqrt(2.0_R8*pi))*t1*t1*ev*sqrt(ev)*sqrt(me)
       tau = t2/dv%lnlam*pl%te/ev*sqrt(pl%te/ev)/dv%ne2
       do iCv = 1, nCv
         if (mpg%cvOnClosedSurface(iCv)) then
           iFt = mpg%cvFt(iCv)
           if (iFt.gt.0.and.iFt.le.mpg%nFt) then
             nu1(iCv) = geo%ftConn(iFt)/(tau(iCv)*
     &        sqrt(pl%te(iCv)/me)*2.0_R8*pi*geo%ftEps(iFt)**1.5_R8)
             k1(iCv) = 1.0_R8/(1.0_R8+2.01_R8*sqrt(nu1(iCv))
     &        +1.53_R8*nu1(iCv))
     &        +0.52_R8*(geo%ftEps(iFt)**3)*nu1(iCv)
     &        /(1.0_R8+0.89_R8*sqrt(geo%ftEps(iFt)**3)*nu1(iCv))
             co%nu2(iCv) = geo%ftConn(iFt)/(tauia(iCv,ismain)
     &        *sqrt(pl%ti(iCv)/(mp*am(ismain)))
     &        *2.0_R8*pi*geo%ftEps(iFt)**1.5_R8)
             co%k2(iCv) = (0.66_R8+1.88_R8*sqrt(geo%ftEps(iFt))
     &        -1.54_R8*geo%ftEps(iFt))
     &        /(1.0_R8+1.03_R8*sqrt(co%nu2(iCv))+0.31_R8*co%nu2(iCv))
     &        +1.17_R8*(geo%ftEps(iFt)**3)*co%nu2(iCv)
     &        /(1.0_R8+0.74_R8*sqrt(geo%ftEps(iFt)**3)*co%nu2(iCv))
           endif
         endif
       enddo                                                          !lk 19.09.07 }
*    ..interpolate te, ti and ne2 to cell faces
       wrk1 = 1.0_R8
       call intface(nCv,nFc,mpg%fcCv,wrk1,pl%te,tef)
       call intface(nCv,nFc,mpg%fcCv,wrk1,pl%ti,tif)
       call intface(nCv,nFc,mpg%fcCv,wrk1,dv%ne2,ne2f)
       call intface(nCv,nFc,mpg%fcCv,wrk1,tau,tauf)
       call intface(nCv,nFc,mpg%fcCv,wrk1,tauia(:,ismain),tauiaf)
*    ..compute flux limits at faces
       do iFc = 1, nFc                                                !lkw 12.10.22 {
        if (mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .or.
     &      mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
          ! compute conn and eps on the face
          if (mpg%fcFs(iFc).gt.0) then
            connf = geo%fsConn(mpg%fcFs(iFc))
            epsf  = 0.0_R8
            if (mpg%cvOnClosedSurface(mpg%fcCv(iFc,1))) then
              iFt = mpg%cvFt(mpg%fcCv(iFc,1))
              if (iFt.gt.0.and.iFt.le.mpg%nFt) then
                epsf  = epsf + geo%ftEps(iFt)
              end if
            end if
            if (mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
              iFt = mpg%cvFt(mpg%fcCv(iFc,2))
              if (iFt.gt.0.and.iFt.le.mpg%nFt) then
                epsf  = epsf + geo%ftEps(iFt)
              end if
            end if
            if (mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &          mpg%cvOnClosedSurface(mpg%fcCv(iFc,2)))
     >       epsf = 0.5_R8
          else
            iFt = mpg%cvFt(mpg%fcCv(iFc,1))
            if (iFt.gt.0.and.iFt.le.mpg%nFt) then
              connf = geo%ftConn(iFt)
              epsf  = geo%ftEps(iFt)
            end if
          end if
          ! parameters at faces
          nu1f = connf/(tauf(iFc)*
     &     sqrt(tef(iFc)/me)*2.0_R8*pi*epsf**1.5_R8)
          k1f = 1.0_R8/(1.0_R8+2.01_R8*sqrt(nu1f)
     &     +1.53_R8*nu1f)
     &     +0.52_R8*(epsf**3)*nu1f
     &     /(1.0_R8+0.89_R8*sqrt(epsf**3)*nu1f)
          nu2f = connf/(tauiaf(iFc)
     &     *sqrt(tif(iFc)/(mp*am(ismain)))
     &     *2.0_R8*pi*epsf**1.5_R8)
          k2f = (0.66_R8+1.88_R8*sqrt(epsf)
     &     -1.54_R8*epsf)
     &     /(1.0_R8+1.03_R8*sqrt(nu2f)+0.31_R8*nu2f)
     &     +1.17_R8*(epsf**3)*nu2f
     &     /(1.0_R8+0.74_R8*sqrt(epsf**3)*nu2f)
          ! compute flux limits
          co%f_luc_ke(iFc) = 1.0_R8/k1f*0.58_R8*sqrt(epsf**3)
          co%f_luc_ki(iFc) = 1.0_R8/k2f*1.60_R8*sqrt(epsf**3)
          if (connf.gt.0.0_R8) then
            co%f_luc_et(iFc)=switch%cvsa_mltpl/
     1       (1.0_R8+15.12e16_R8/connf*
     2       (tif(iFc)/ev)**2/ne2f(iFc))/
     4       (1.0_R8+15.12e16_R8*epsf**1.5_R8/
     6        connf*(tif(iFc)/ev)**2/ne2f(iFc))
          else
            co%f_luc_et(iFc)=switch%cvsa_mltpl
          end if
          co%f_luc_sg(iFc) = co%f_luc_ke(iFc)
          co%f_luc_al(iFc) = co%f_luc_ke(iFc)
          do is=0,ns-1
            if(.not.is_neutral(is)) then
              cvsa(iFc,0,is)    = cvsa(iFc,0,is)*co%f_luc_et(iFc)
              cvsahz(iFc,0,is)  = cvsahz(iFc,0,is)*co%f_luc_et(iFc)
              vsaf_cl(iFc,0,is) = vsaf_cl(iFc,0,is)*co%f_luc_et(iFc)
            endif
          enddo
          chce(iFc,0) = chce(iFc,0)*co%f_luc_ke(iFc)
          csig(iFc,0) = csig(iFc,0)*co%f_luc_sg(iFc)
          calf(iFc,0) = calf(iFc,0)*co%f_luc_al(iFc)
          chci(iFc,0) = chci(iFc,0)*co%f_luc_ki(iFc)
        endif                                                         !lkw 12.10.22 {
       enddo
       do iCv = 1, nCv
        if (mpg%cvOnClosedSurface(iCv)) then
          iFt = mpg%cvFt(iCv)
          if (iFt.gt.0.and.iFt.le.mpg%nFt) then
            co%f_luc_ke_c(iCv) = 1.0_R8/k1(iCv)*0.58_R8*
     &        sqrt(geo%ftEps(iFt)**3)
            co%sigx_c(iCv) = co%sigx_c(iCv)*co%f_luc_ke_c(iCv)
            co%alfx_c(iCv) = co%alfx_c(iCv)*co%f_luc_ke_c(iCv)
            co%f_luc_ki_c(iCv) = 1.0_R8/co%k2(iCv)*1.60_R8*
     &        sqrt(geo%ftEps(iFt)**3)
            co%hcix_c(iCv) = co%hcix_c(iCv)*co%f_luc_ki_c(iCv)
          endif
        endif
       enddo
      endif
#ifdef WG_TODO
#ifdef ASTRA
      bbaverage2 = 0.0_R8                                             !srv 28.07.11 {
      chci_koeff = 0.0_R8
      tt0 = 0.0_R8
      tt1 = 0.0_R8
      if(integral_modeling.ne.0) then
       do iy = 0,mintopcut-1    ! for inner surfaces
        inv_hi_astra(iy) =  -(Ti_astra(iy)-Ti_astra(iy-1))*1.0e3_R8/
     &   (qi_neo_astra(iy)*1.0e6_R8)/(Rmajc_astra(iy)-Rmajc_astra(iy-1))
        do ix=0,nx
         if(on_closed_surface(ix,iy)) then
          tt0(iy) = tt0(iy) +bb(ix,iy,3)*bb(ix,iy,3)*vol(ix,iy)
          tt1(iy) = tt1(iy) +vol(ix,iy)
         endif
        enddo
        bbaverage2(iy) = tt0(iy)/tt1(iy)
        tt0(iy) = 0.0_R8
        do ix=0,nx
         if(on_closed_surface(ix,iy)) then
          tt0(iy) = tt0(iy) +vol(ix,iy)/bb(ix,iy,3)/bb(ix,iy,3)*
     &     (1.0_R8-bb(ix,iy,3)*bb(ix,iy,3)/bbaverage2(iy))**2
         endif
        enddo
        !wdk For WG, need to add flux limits for hcix_c !!
        do ix=0,nx
         if(on_closed_surface(ix,iy)) then
          chci_koeff(ix,iy) = 25.0_R8/4.0_R8*hy1(nmdpl,iy)/qe*tt0(iy)* !srv 22.04.16
     &     0.5_R8*((bb(ix,iy,2)/bb(ix,iy,3)/hy1(ix,iy)/hx(ix,iy))**2  !srv 22.04.16
     &     *vol(ix,iy)+(bb(leftix(ix,iy),iy,2)/bb(leftix(ix,iy),iy,3)
     &     /hy1(leftix(ix,iy),iy)/hx(leftix(ix,iy),iy))**2            !srv 22.04.16
     &     *vol(leftix(ix,iy),iy))*inv_hi_astra(iy) !     chi_N_astra(iy)
         endif
        enddo
       enddo
!       do iy = 0, iy_Ratb    ! for the region of overlap with ASTRA
        !wdk For WG, need to add flux limits for hcix_c !!
       do iy=0,mintopcut-1
        do ix=0,nx
         if(on_closed_surface(ix,iy)) then
          chci(ix,iy,0,0) = chci_koeff(ix, iy)*
     &     0.5_R8*((ti(ix,iy)*na(ix,iy,ismain))**2
     &     +(ti(leftix(ix,iy),iy)*na(leftix(ix,iy),iy,ismain))**2)
         endif
        enddo
       enddo                                                          !srv 28.07.11 }
      endif
#endif
#endif
*   ..set strange velocities chve, chvi
      call sfill (nFc*2, 0.0_R8, chve, 1)
      call sfill (nFc*2, 0.0_R8, chvi, 1)

      if(switch%b2npmo_iout.ne.0 .or. switch%iout_b2wdat.eq.4) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        if(.not. is_neutral(is))
     &  call my_out_us(70,nFc,1,cvsahz(1,0,is),
     &              'b2trcl_luciani_cvsahzx'//chns)
!        call my_out_us(70,nFc,1,cvsahz(1,1,is),
!     &              'b2trcl_luciani_cvsahzy'//chns)
       enddo
      endif
      if(switch%b2npht_iout.ne.0 .or. switch%iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,chce(1,0),'b2trcl_luciani_chce_th')
!        call my_out_us(70,nFc,1,chce(1,1),'b2trcl_luciani_chcey')
        call my_out_us(70,nFc,1,calf(1,0),'b2trcl_luciani_calf_th')
        call my_out_us(70,nCv,0,co%alfx_c,'b2trcl_luciani_alfx_c')    !srv 05.07.17
!        call my_out_us(70,nFc,1,calf(1,1),'b2trcl_luciani_calfy')
        call my_out_us(70,nFc,1,chci(1,0),'b2trcl_luciani_chci_th')
!        call my_out_us(70,nFc,1,chci(1,1),'b2trcl_luciani_chciy')
      endif
      if(switch%b2npp7_iout.ne.0 .or. switch%b2nppo_iout.ne.0
     &                  .or. switch%iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,csig(1,0),'b2trcl_luciani_csig_th')
!        call my_out_us(70,nFc,1,csig(1,1),'b2trcl_luciani_csigy')
        call my_out_us(70,nCv,0,co%sigx_c,'b2trcl_luciani_sig_th_c')  !srv 05.07.17
      endif

*   .. save chce(:,0) after applying Luciani
      wrk = chce(:,0)                                                 !srv 13.01.17

      if (switch%min_collisions.gt.0.0_R8) then
        call b2txnci (mpg, geo, switch, nCv, ns, pl%ti, pl%te, rt%rz2,
     &   dv%ne, dv%ne2, pl%na, dv%lnlam, ismain, st_ext, collisnumf,
     &   wrkc)
        call b2txnce (mpg, geo, switch, nCv, pl%te, dv%ne2,
     &   dv%lnlam, collisnumef, wrkc)
      else
        collisnumf = 1.0_R8
        collisnumef = 1.0_R8
      endif

      if(switch%conductive_limit.ne.0) then                           !srv 16.06.08 {
*   .. pre-compute poloidal differences of te, ti
       call diff_p(nCv, nFc, nVx, 0, geo, mpg, pl%te, wrkv, dte)
       call diff_p(nCv, nFc, nVx, 0, geo, mpg, pl%ti, wrkv, dti)
*   .. apply the conductive limit of parallel electron heat flux      !srv 03.06.03 {
       if(cflim(0).ne.0.0_R8) then
        do iFc = 1, nFc
csw 19oct2011, also in core for SOLPS4.3 comparison!!!
c            if(.not.on_closed_surface(ix,iy)) then                   !lk 06.05.07
         if(switch%core_cond_limit.ne.0.or.
     &      .not.(mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)).and.
     &            mpg%cvOnClosedSurface(mpg%fcCv(iFc,2)))) then
          if (collisnumef(iFc).gt.switch%min_collisions) then
           flomx = co%chvemx(iFc)*(dv%ne(mpg%fcCv(iFc,1))+
     &                             dv%ne(mpg%fcCv(iFc,2)))/2.0_R8
           t0 = abs(chce(iFc,0)*dte(iFc))/(flomx*max(
     &      pl%te(mpg%fcCv(iFc,1)),pl%te(mpg%fcCv(iFc,2))) +
     &      b2trcl_cutlo)
          else
           t0 = (cflim(0)*max(collisnumef(iFc), 1.0e-16_R8))**(-1)
          endif
          co%fllim_ke(iFc)=1.0_R8/(1.0_R8+t0)
          chce(iFc,0)=chce(iFc,0)*co%fllim_ke(iFc)
         endif
        enddo
       endif
*   .. apply the conductive limit of parallel ion heat flux
       if(cflim(1).ne.0.0_R8) then
        do iFc = 1,nFc
csw 19oct2011, also in core for SOLPS4.3 comparison!!!
c          if(.not.on_closed_surface(ix,iy)) then                     !lk 06.05.07
         if(switch%core_cond_limit.ne.0.or.
     &      .not.(mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)).and.
     &            mpg%cvOnClosedSurface(mpg%fcCv(iFc,2)))) then
          if (collisnumf(iFc) .gt. switch%min_collisions) then
           flomx = co%chvimx(iFc)*(dv%ni(mpg%fcCv(iFc,1),0)
     &                            +dv%ni(mpg%fcCv(iFc,2),0))/2.0_R8
           t0 = abs(chci(iFc,0)*dti(iFc))/(flomx*max(
     &      pl%ti(mpg%fcCv(iFc,1)),pl%ti(mpg%fcCv(iFc,2))) +
     &      b2trcl_cutlo)
          else
           t0 = (cflim(1)*max(collisnumf(iFc), 1.0e-16_R8))**(-1)
          endif
          co%fllim_ki(iFc) = 1.0_R8/(1.0_R8+t0)
          chci(iFc,0)=chci(iFc,0)*co%fllim_ki(iFc)
         endif
        enddo
        !wdk pragmatic approximation for cell center flux limits: based on interpolation
        !wdk of face flux limit
        call intcell(nFc,nCv,mpg,mpg%intcellP,co%fllim_ki,co%fllim_ki_c)
        do iCv = 1, nCv
         if(switch%core_cond_limit.ne.0.or.
     &      .not.mpg%cvOnClosedSurface(iCv)) then
          co%hcix_c(iCv)=co%hcix_c(iCv)*co%fllim_ki_c(iCv)
         endif
        enddo
       endif
      endif                                                           !srv 16.06.08 }
*   ..apply flux limit to calf                                        !srv 13.01.17 {
      if(switch%b2sigp_style.ne.2) then
       if(cflim(0).ne.0.0_R8) then
        do iFc = 1, nFc
         if(.not.(mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &            mpg%cvOnClosedSurface(mpg%fcCv(iFc,2)))) then
            flomx = co%chvemx(iFc)*(dv%ne(mpg%fcCv(iFc,1))+
     &                              dv%ne(mpg%fcCv(iFc,2)))/2.0_R8
            t0 = abs(wrk(iFc)*dte(iFc))/(flomx*max(
     &        pl%te(mpg%fcCv(iFc,1)),pl%te(mpg%fcCv(iFc,2))) +
     &        b2trcl_cutlo)
            co%fllim_al(iFc)=1.0_R8/(1.0_R8+t0)
            calf(iFc,0)=calf(iFc,0)*co%fllim_al(iFc)                  !srv 13.01.17
         endif
        enddo
        !wdk pragmatic approximation for cell center flux limits: based on interpolation
        !wdk of face flux limit
        call intcell(nFc,nCv,mpg,mpg%intcellP,co%fllim_al,co%fllim_al_c)
        do iCv = 1, nCv
          if(.not.mpg%cvOnClosedSurface(iCv)) then
            co%alfx_c(iCv)=co%alfx_c(iCv)*co%fllim_al_c(iCv)
          endif
        enddo
       endif
      else
       if(cflim(3).ne.0.0_R8) then
*   ..apply flux limit to calf at cell faces
        do iFc = 1, nFc
         if(.not.(mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &            mpg%cvOnClosedSurface(mpg%fcCv(iFc,2)))) then
            jemax = sqrt(
     &       (pl%te(mpg%fcCv(iFc,1))+pl%te(mpg%fcCv(iFc,2)))/2.0_R8/me)*
     &        qe*(dv%ne(mpg%fcCv(iFc,1))+dv%ne(mpg%fcCv(iFc,2)))/2.0_R8*
     &        abs(geo%fcPbs(iFc))
            t0 = abs(calf(iFc,0)*dte(iFc))/(cflim(3)*jemax +
     &        b2trcl_cutlo)
            co%fllim_al(iFc)=1.0_R8/(1.0_R8+t0)
            calf(iFc,0)=calf(iFc,0)*co%fllim_al(iFc)                  !srv 13.01.17
         endif
        enddo
*   ..apply flux limit to alfx_c at cell centers
        call gradc_p(nCv, nFc, nVx, 0, geo, mpg, pl%te, wrkv, gtec)
        do iCv = 1, nCv
          if(.not.mpg%cvOnClosedSurface(iCv)) then
            jemax = sqrt(pl%te(iCv)/me)*
     &         qe*dv%ne(iCv)*geo%cvBb(iCv,0)/geo%cvBb(iCv,3)
            t0 = abs(co%alfx_c(iCv)*gtec(iCv)/
     &         (cflim(3)*jemax + b2trcl_cutlo))
            co%fllim_al_c(iCv)=1.0_R8/(1.0_R8+t0)
            co%alfx_c(iCv)=co%alfx_c(iCv)*co%fllim_al_c(iCv)          !iyv 19.11.13 !srv 08.10.14
          endif
        enddo
       endif
      endif
*   ..apply flux limit to cvsa
      cflmv = cflim(2)
      do is = 0, ns-1
       call b2tlmv (nCv, nFc, nVx, is, cflmv, switch, geo, mpg,
     &   pl%na(:,is),
     &   pl%ti, pl%ua(:,is), collisnumf, vsaf_cl(:,0,is),
     &   cvsa(:,0,is), cvsahz(:,0,is), fllimvisc(:,is))
      enddo
*   ..computation conductivity of kt due to classical parallel current
!wdk  current model: conductivity ~ sigx*rhol**2/(am(ismain)*mp)
      wrkc = sqrt(2.0_R8*pl%ti*(am(ismain)*mp))/
     &       dabs(rt%rza(:,ismain)*qe*geo%cvBb(:,3))
      wrkc = co%sigx_c*wrkc**2*geo%cvBb(:,3)**2/(am(ismain)*mp)
      call b2txcx (nCv, nFc, 0.0_R8,
     & mode, geo, mpg, geo%fcVol, geo%fcS, wrkc, cdkt(:,0))
      cdkt(:,0) = cdkt(:,0)*switch%b2tfhi_fsigkt
      cdkt(:,1) = 0.0_R8
      co%sigx_kt = wrkc
*   ..save the classical transport coefficients
      csig_cl = csig
      calf_cl = calf
      cvsa_cl = cvsa
      cvsahz_cl = cvsahz                                              !srv 09.01.01

      if(switch%b2npmo_iout.ne.0 .or. switch%iout_b2wdat.eq.4) then
       do is = 0, ns-1
       write (chns,'(i3.3)') is
        call my_out_us(70,nFc,1,cvsahz(1,0,is),
     &              'b2trcl_luciani_fllim_cvsahz_th'//chns)
!        call my_out_us(70,nFc,1,cvsahz(1,1,is),
!     &              'b2trcl_luciani_fllim_cvsahzy'//chns)
       enddo
      endif
      if(switch%b2npht_iout.ne.0 .or. switch%iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,chce(1,0),
     &             'b2trcl_luciani_fllim_chce_th')
!        call my_out_us(70,nFc,1,chce(1,1),'b2trcl_luciani_fllim_chcey')
        call my_out_us(70,nFc,1,calf(1,0),
     &             'b2trcl_luciani_fllim_calf_th')
!        call my_out_us(70,nFc,1,calf(1,1),'b2trcl_luciani_fllim_calfy')
        call my_out_us(70,nCv,0,co%alfx_c,
     &             'b2trcl_luciani_fllim_alf_th_c')                   !srv 05.07.17
        call my_out_us(70,nFc,1,chci(1,0),
     &             'b2trcl_luciani_fllim_chci_th')
!        call my_out_us(70,nFc,1,chci(1,1),'b2trcl_luciani_fllim_chciy')
      endif

      if(switch%b2trcl_iout.ne.0) then                                !srv 17.06.08 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out_us(70,nFc,1,cvsa(1,0,is),'cl_cvsa_th'//chns)
        call my_out_us(70,nFc,1,cvsa(1,1,is),'cl_cvsa_r'//chns)
        call my_out_us(70,nFc,1,cvsahz(1,0,is),'cl_cvsahz_th'//chns)
        call my_out_us(70,nFc,1,cvsahz(1,1,is),'cl_cvsahz_r'//chns)
        call my_out_us(70,nFc,1,cvsa_cl(1,0,is),'cl_cvsa_cl_th'//chns)
        call my_out_us(70,nFc,1,cvsa_cl(1,1,is),'cl_cvsa_cl_r'//chns)
        call my_out_us(70,nFc,1,cvsahz_cl(1,0,is),
     &              'cl_cvsahz_clx'//chns)
        call my_out_us(70,nFc,1,cvsahz_cl(1,1,is),
     &              'cl_cvsahz_cly'//chns)
        call my_out_us(70,nCv,0,cthe(1,is),'cl_cthe'//chns)
        call my_out_us(70,nCv,0,cthi(1,is),'cl_cthi'//chns)
        do is0=0,ns-1
         write (chns0,'(i3.3)') is0
         call my_out_us(70,nFc,1,csigin(1,0,is,is0),'cl_csigin_th'//
     &        chns//'_'//chns0)
         call my_out_us(70,nFc,1,csigin(1,1,is,is0),'cl_csigin_r'//
     &        chns//'_'//chns0)
        enddo
       enddo
       call my_out_us(70,nFc,1,chce(1,0),'cl_chce_th')
       call my_out_us(70,nFc,1,chce(1,1),'cl_chce_r')
       call my_out_us(70,nFc,1,chci(1,0),'cl_chci_th')
       call my_out_us(70,nFc,1,chci(1,1),'cl_chci_r')
       call my_out_us(70,nFc,1,chve(1,0),'cl_chve_th')
       call my_out_us(70,nFc,1,chve(1,1),'cl_chve_r')
       call my_out_us(70,nFc,1,chvi(1,0),'cl_chvi_th')
       call my_out_us(70,nFc,1,chvi(1,1),'cl_chvi_r')
       call my_out_us(70,nFc,1,csig(1,0),'cl_csig_th')
       call my_out_us(70,nFc,1,csig(1,1),'cl_csig_r')
       call my_out_us(70,nFc,1,calf(1,0),'cl_calf_th')
       call my_out_us(70,nFc,1,calf(1,1),'cl_calf_r')
       call my_out_us(70,nCv,0,co%sigx_c,'cl_sig_th_c')               !srv 05.07.17
       call my_out_us(70,nCv,0,co%alfx_c,'cl_alf_th_c')               !srv 05.07.17
      endif                                                           !srv 17.06.08 }

* ..return
      ncall_b2trcl = ncall_b2trcl+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2trcl

      end subroutine b2trcl
