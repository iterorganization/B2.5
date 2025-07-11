*-----------------------------------------------------------------------
*.specification

      subroutine b2sihs_ (nCv, nFc, nVx, ns, switch, geo, mpg,
     &  ismain, fac_ExB, lnlam, na, ua, ue, te, ti, tn, po, ne, ne2,
     &  ni, nn, fch, cvsa, cvsa_cl, cvsa_drho, cvsahz_eff,
     &  f_luc_sg, alfx_c, sigx_c, rza, rz2, st_ext,
     &  she0, shi0, shn0, srw)
      use b2mod_types
      use b2mod_tallies
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_zhfrtf
     & , only: avm_u,                                                 !som 21.03.21 16.08.21
     &         b2mod_zhfrtf_avm_u                                     !som 18.08.21
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
#ifndef FIXED_POINT
      use b2mod_balance !djm Jan2017
     & , only : b2sihs_divue0to3, b2sihs_exbe0to3,
     &          b2sihs_joule0to3, b2sihs_divua0to3,
     &          b2sihs_exba0to3, b2sihs_visa0to3, b2sihs_fraa0to3,
     &          balance_netcdf
#endif
      use b2mod_transport_fun                                         !srv 05.07.17
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2sihs_
#ifndef NAGFOR
     &        , ncall_b2xehy, ncall_b2xehx, ncall_b2tlnl,
#ifdef WG_DEBUG
     &          my_out_sf,                                              !srv 27.06.21
#endif
     &          my_out_folder, ncall_b2sifrtf
      use b2mod_zhfrtf
     & , only : rho_a_rel
#endif
      use b2mod_subsys
      use b2mod_geometry
      use b2mod_b2plot_wall_loading
     & , only: itrg
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain                               !srv 05.07.17
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2StateExt), intent (in) :: st_ext
      type (B2SourceWork), intent (inout) :: srw
      real (kind=R8) ::
     *  fac_ExB(nFc), na(nCv,0:ns-1), ua(nCv,0:ns-1), ue(nCv),        !srv 05.03.11
     &  te(nCv), ti(nCv), tn(nCv), po(nCv), ne(nCv), ni(nCv,0:1),   !srv 07.05.08 01.07.09
     &  lnlam(nCv), f_luc_sg(nFc), alfx_c(nCv), sigx_c(nCv),
     &  nn(nCv), ne2(nCv)
      real (kind=R8) ::
     *  fch(nFc,0:1),
     *  cvsa(nFc,0:1,0:ns-1), cvsahz_eff(nFc,0:1,0:ns-1),
     *  cvsa_cl(nFc,0:1,0:ns-1), cvsa_drho(nFc,0:1,0:ns-1),
     *  rza(nCv,0:ns-1), rz2(nCv,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   she0(nCv,0:3), shi0(nCv,0:3), shn0(nCv,0:3)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIHS computes the source coefficients for a set of terms in the
*     heat equations. These terms specify the contributions of the
*     divergence of the flow field, of viscous heating, and of the
*     electron-atom and atom-atom friction.
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      integer iCv, iFc, is, isb, k
      integer :: gridGeometry
      real (kind=R8) :: mm(0:ns-1), am1(0:ns-1)                         !som 13.11.20
      real (kind=R8), parameter :: eps = 1.0e-60_R8
      real (kind=R8) ::
     *  vti32(nCv),
     &  wrkc(nCv), wrkf(nFc,0:1), wrkf1(nFc,0:1),
     &  wrkv(nVx),                                                      !srv 02.06.99
     &  wrkc1(nCv),                                                     !srv 28.08.02
     &  dpoc(nCv,0:1), dOnedBsqc(nCv,0:1), visc(nCv,0:1)
#ifdef WG_TODO
      real (kind=R8) :: csig_stoch(nFc), fchstoch(nCv)                  !lk 31.01.08 }
#endif
      real (kind=R8) ::
     *  shidu(nCv,0:3), shidd(nCv,0:3),                                 !srv 11.09.09 {
     *  shiva(nCv,0:3), shivc(nCv,0:3),
     *  shivan(nCv,0:3), shivcn(nCv,0:3),
     *  shifr(nCv,0:3), shedu(nCv,0:3),
     *  shedd(nCv,0:3), shefr(nCv,0:3),
     *  shest(nCv,0:3),                                                 !srv 17.12.13
     *  shidun(nCv,0:3),
     *  wrks(nCv), ehx(nFc), ehy(nFc)                                   !srv 13.01.17i
      real (kind=R8) ::
     3  zeff(nCv), zmff(nCv), ce1(nCv),                                 !srv 13.01.17 05.07.17 {
     4  ce2(nCv), cimp1(nCv,0:ns-1), cimp2(nCv,0:ns-1),
     5  zetap(nCv), zetae(nCv),
     6  gti(nCv), gte(nCv),
     7  zmc(nCv), z11(nCv),                                             !som 13.11.20
     8  z12(nCv), z22(nCv),
     9  zf1(nCv), zf2(nCv),
     1  smbch(nCv,0:3), smbtf(nCv,0:3),
     *  smbfria(nCv,0:3), smbfrea(nCv,0:3),
     *  smbtfia(nCv,0:3), smbtfea(nCv,0:3)                              !srv 13.01.17 05.07.17 }
      character chk*1                                                   !srv 11.09.09 } !srv 05.07.17
      real (kind=R8) ::
     *   t0, t1, ctaup, coef, czetap, czetae, ft                        !srv 05.07.17
      real (kind=R8) ::
     *   fac_ExB_c
*   ..procedures
      intrinsic min, max, sqrt
      external xertst, ipgeti, ipgetr, sfill
      external b2xvfv, b2xvsg
*   ..initialization
*     Note 0: style determines the form of the strange electron-atom
*     energy transfer term.
*     Note 1: rf0-5 are numerical parameters, see the code. All must
*     be non-negative and of order unity; larger values are meant to
*     provide more stabilization. Some or all may be set to zero and
*     removed eventually. I have had some experience to indicate that
*     rf0, rf2, rf4 may have a useful function; however, later changes
*     in the code may have made all these parameters superfluous.

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sihs_')
*   ..test nCv, nFc, ns
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall_b2sihs_.lt.3) then
*    ..test sign of vol
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
*    ..test sign of na, te, ti, ne, ni
       call b2xvsg (nCv*ns, na, 1, 'na', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, tn, 1, 'tn', '.gt.')
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, ni, 1, 'ni', '.gt.')
      endif
      gridGeometry = geometryId( mpg, geo, 1 )
*   ..compute ctaup
*     Following Balescu; Braginskii has ctaup different by a factor
*     sqrt(2). The mathematical expression is:
*       t0 = (3/(4*sqrt(2*pi)))*(4*pi*eps0/qe**2)**2*ev*sqrt(ev)
*       ctaup = t0*sqrt(mp)/lnlam
*     The numerical value in SI units is: ctaup = 1.4743d13/lnlam.
*     The following calculation of t0 is designed to avoid overflow
*     even for IEEE single precision arithmetic.
      t1 = (4*pi*eps0/qe)/qe
      t0 = (3/(4*sqrt(2*pi)))*t1*(t1*ev*sqrt(ev))
      ctaup = t0*sqrt(mp) ! /lnlam it is taken below     !srv 01.07.09
*   ..compute vti32
      vti32 = geo%cvVol/((ti/ev)*sqrt(ti/ev))

* ..compute heat source terms
*   ..initialize sources to zero
      call sfill (nCv*4, 0.0_R8, she0, 1)
      call sfill (nCv*4, 0.0_R8, shi0, 1)
      call sfill (nCv*4, 0.0_R8, shn0, 1)
      call sfill (nCv*4, 0.0_R8, shedu, 1)                            !srv 11.09.09 {
      call sfill (nCv*4, 0.0_R8, shedd, 1)
      call sfill (nCv*4, 0.0_R8, shefr, 1)
      call sfill (nCv*4, 0.0_R8, shidu, 1)
      call sfill (nCv*4, 0.0_R8, shidun, 1)
      call sfill (nCv*4, 0.0_R8, shidd, 1)
      call sfill (nCv*4, 0.0_R8, shifr, 1)
      call sfill (nCv*4, 0.0_R8, shiva, 1)
      call sfill (nCv*4, 0.0_R8, shivc, 1)                            !srv 11.09.09 }
      call sfill (nCv*4, 0.0_R8, shivan, 1)
      call sfill (nCv*4, 0.0_R8, shivcn, 1)                           !srv 11.09.09 }
      call sfill (nCv*4, 0.0_R8, shest, 1)                            !srv 17.12.13
*   ..compute the Coulomb logarithm
      call b2tlnl (nCv, switch, switch%icase_ii, te, ti, ne, lnlam)   !srv 01.07.09  20.09.11
*   ..contribution from stochastic electron flow
      if(switch%fch_stochastic.ne.0.0_R8) then                        !lk 31.01.08 { !srv 17.12.13
       call xerrab ('b2sihs_: fch_stochastic not yet available for WG')
#ifdef WG_TODO
       call b2tstch (nx, ny, po, te, ne, csig_stoch, fchstoch)
       do iy = 0, ny-1
        do ix = 0, nx-1
         if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
          shest(ix,iy,0) = shest(ix,iy,0)-
     &    (ne(topix(ix,iy),topiy(ix,iy))*te(topix(ix,iy),topiy(ix,iy))
     &    -ne(bottomix(ix,iy),bottomiy(ix,iy))
     &    *te(bottomix(ix,iy),bottomiy(ix,iy)))
     &    /hy(ix,iy)/qe/ne(ix,iy)*vol(ix,iy)
     &    *(fchstoch(ix,iy)*(hy(ix,iy)
     &    +hy(bottomix(ix,iy),bottomiy(ix,iy)))
     &    /(vol(ix,iy) + vol(bottomix(ix,iy),bottomiy(ix,iy)))
     &    +fchstoch(topix(ix,iy),topiy(ix,iy))*(hy(ix,iy)
     &    +hy(topix(ix,iy),topiy(ix,iy)))
     &    /(vol(ix,iy) + vol(topix(ix,iy),topiy(ix,iy))))/4.0_R8
         endif
        enddo
       enddo
#endif
      endif                                                           !lk 31.01.08 }
      call b2saxpy (nCv*4, 1.0_R8, shest, 1, she0, 1)                 !srv 17.12.13
*   ..contribution from divergence ue
      if(switch%b2sihs_style_int_vel.eq.0) then
        call intface(nCv,nFc,mpg%fcCv,geo%fcHc,ue,wrkf(:,0))
      else
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,ue,wrkf(:,0))
      endif
      wrkf(:,0) = wrkf(:,0)*geo%fcPbs
      call div0(nCv, nFc, mpg, wrkf(:,0), wrkc)
      b2divue=0.0_R8
      srw%b2sihs_divue=0.0_R8
      do iCv = 1, mpg%nCi
         t0 = -switch%b2sihs_phm0*(1.0_R8-switch%BoRiS)*wrkc(iCv)
         shedu(iCv,0) = shedu(iCv,0)+
     &    (max(0.0_R8,t0)+switch%b2sihs__rf0*abs(t0))*ne(iCv)*te(iCv)
         shedu(iCv,3) = shedu(iCv,3)+min(0.0_R8,t0)-
     &                  switch%b2sihs__rf0*abs(t0)
         b2divue(mpg%cvReg(iCv))=b2divue(mpg%cvReg(iCv))+
     1     t0*ne(iCv)*te(iCv)
         srw%b2sihs_divue(iCv)=t0*ne(iCv)*te(iCv)
      enddo
      if (switch%b2sihs_shedu.ne.0) then
        do k=mpg%nCi+1,mpg%nCv
          iFc = mpg%cvFc(mpg%cvFcP(k,1))
          iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
          if (switch%b2sihs_shedu.eq.1) then
            shedu(iCv,:) = 0.0_R8
          elseif (switch%b2sihs_shedu.eq.2) then
            select case ( gridGeometry )
            case ( GEOMETRY_LINEAR )
              if (mpg%fcReg(iFc).eq.1 .and. itrg(1).ne.0 .or.
     &            mpg%fcReg(iFc).eq.2 .and. itrg(2).ne.0) cycle
            case ( GEOMETRY_LIMITER )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.2 .or.
     &            mpg%fcReg(iFc).eq.4) cycle
            case ( GEOMETRY_SN )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.8) cycle
            case ( GEOMETRY_CDN )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.5 .or. mpg%fcReg(iFc).eq.8 .or.
     &            mpg%fcReg(iFc).eq.14 .or. mpg%fcReg(iFc).eq.21) cycle
            case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.5 .or. mpg%fcReg(iFc).eq.8 .or.
     &            mpg%fcReg(iFc).eq.15 .or. mpg%fcReg(iFc).eq.22) cycle
            case ( GEOMETRY_LFS_SNOWFLAKE_MINUS,
     &             GEOMETRY_LFS_SNOWFLAKE_PLUS )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.8 .or. mpg%fcReg(iFc).eq.15) cycle
            end select
            shedu(iCv,:) = 0.0_R8
          elseif (switch%b2sihs_shedu.eq.3) then
            select case ( gridGeometry )
            case ( GEOMETRY_LINEAR )
              if (mpg%fcReg(iFc).ne.3) cycle
            case ( GEOMETRY_LIMITER )
              if (mpg%fcReg(iFc).ne.6) cycle
            case ( GEOMETRY_SN )
              if (mpg%fcReg(iFc).ne.12) cycle
            case ( GEOMETRY_CDN )
              if (mpg%fcReg(iFc).ne.18 .and. mpg%fcReg(iFc).ne.25) cycle
            case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
              if (mpg%fcReg(iFc).ne.19 .and. mpg%fcReg(iFc).ne.26) cycle
            case ( GEOMETRY_LFS_SNOWFLAKE_MINUS,
     &             GEOMETRY_LFS_SNOWFLAKE_PLUS )
              if (mpg%fcReg(iFc).ne.19) cycle
            end select
            shedu(iCv,:) = 0.0_R8
          endif
        enddo
      endif
#ifndef FIXED_POINT
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_divue0to3 = shedu
      endif
#endif
      call b2saxpy (nCv*4, 1.0_R8, shedu, 1, she0, 1)                 !srv 11.09.09
*   ..contribution from divergence ua
      b2divua=0.0_R8
      srw%b2sihs_divua=0.0_R8
      do is = 0, ns-1
        if(switch%b2sihs_style_int_vel.eq.0) then
          call intface(nCv,nFc,mpg%fcCv,geo%fcHc,ua(:,is),wrkf(:,0))
        else
          call intface(nCv,nFc,mpg%fcCv,geo%fcVol,ua(:,is),wrkf(:,0))
        endif
        wrkf(:,0) = wrkf(:,0)*geo%fcPbs
        call div0(nCv, nFc, mpg, wrkf(:,0), wrkc)
        do iCv = 1, mpg%nCi
          t0 = -switch%b2sihs_phm1*(na(iCv,is)/ni(iCv,0))
     &                            *(1.0_R8-switch%BoRiS)*wrkc(iCv)     !srv 21.07.05 {
          if (switch%tn_style.eq.0) then
            shidu(iCv,0) = shidu(iCv,0)+
     &        (max(0.0_R8,t0)+switch%b2sihs__rf1*abs(t0))*ni(iCv,0)*
     &         ti(iCv)
            shidu(iCv,3) = shidu(iCv,3)+min(0.0_R8,t0)
     &                     -switch%b2sihs__rf1*abs(t0)
            b2divua(mpg%cvReg(iCv))=b2divua(mpg%cvReg(iCv))
     1         +t0*ni(iCv,0)*ti(iCv)
          elseif (switch%tn_style.eq.1) then
           if (.not.is_neutral(is)) then
            shidu(iCv,0) = shidu(iCv,0)+
     &        (max(0.0_R8,t0)+switch%b2sihs__rf1*abs(t0))*ni(iCv,0)*
     &         ti(iCv)
            shidu(iCv,3) = shidu(iCv,3)+min(0.0_R8,t0)
     &                     -switch%b2sihs__rf1*abs(t0)
            b2divua(mpg%cvReg(iCv))=b2divua(mpg%cvReg(iCv))
     1         +t0*ni(iCv,0)*ti(iCv)
           endif
          else
           if (.not.is_neutral(is).or.nint(zn(is)).ne.1) then
            shidu(iCv,0) = shidu(iCv,0)+
     &        (max(0.0_R8,t0)+switch%b2sihs__rf1*abs(t0))*ni(iCv,0)*
     &         ti(iCv)
            shidu(iCv,3) = shidu(iCv,3)+min(0.0_R8,t0)
     &                     -switch%b2sihs__rf1*abs(t0)
            b2divua(mpg%cvReg(iCv))=b2divua(mpg%cvReg(iCv))
     1         +t0*ni(iCv,0)*ti(iCv)
           else
            t0 = t0*ni(iCv,0)/nn(iCv)
            shidun(iCv,0) = shidun(iCv,0)+
     &        (max(0.0_R8,t0)+switch%b2sihs__rf1*abs(t0))*nn(iCv)*
     &         tn(iCv)
            shidun(iCv,3) = shidun(iCv,3)+min(0.0_R8,t0)-
     &         switch%b2sihs__rf1*abs(t0)
            b2divua(mpg%cvReg(iCv))=b2divua(mpg%cvReg(iCv))
     1         +t0*nn(iCv)*tn(iCv)
           endif
          endif
          srw%b2sihs_divua(iCv)=srw%b2sihs_divua(iCv)+
     1         t0*ni(iCv,0)*ti(iCv)
        enddo
      enddo
      if (switch%b2sihs_shidu.ne.0) then
        do k=mpg%nCi+1,mpg%nCv
          iFc = mpg%cvFc(mpg%cvFcP(k,1))
          iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
          if (switch%b2sihs_shidu.eq.1) then
            shidu(iCv,:) = 0.0_R8
          elseif (switch%b2sihs_shidu.eq.2) then
            select case ( gridGeometry )
            case ( GEOMETRY_LINEAR )
              if (mpg%fcReg(iFc).eq.1 .and. itrg(1).ne.0 .or.
     &            mpg%fcReg(iFc).eq.2 .and. itrg(2).ne.0) cycle
            case ( GEOMETRY_LIMITER )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.2 .or.
     &            mpg%fcReg(iFc).eq.4) cycle
            case ( GEOMETRY_SN )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.8) cycle
            case ( GEOMETRY_CDN )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.5 .or. mpg%fcReg(iFc).eq.8 .or.
     &            mpg%fcReg(iFc).eq.14 .or. mpg%fcReg(iFc).eq.21) cycle
            case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.5 .or. mpg%fcReg(iFc).eq.8 .or.
     &            mpg%fcReg(iFc).eq.15 .or. mpg%fcReg(iFc).eq.22) cycle
            case ( GEOMETRY_LFS_SNOWFLAKE_MINUS,
     &             GEOMETRY_LFS_SNOWFLAKE_PLUS )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.8 .or. mpg%fcReg(iFc).eq.15) cycle
            end select
            shidu(iCv,:) = 0.0_R8
          elseif (switch%b2sihs_shidu.eq.3) then
            select case ( gridGeometry )
            case ( GEOMETRY_LINEAR )
              if (mpg%fcReg(iFc).ne.3) cycle
            case ( GEOMETRY_LIMITER )
              if (mpg%fcReg(iFc).ne.6) cycle
            case ( GEOMETRY_SN )
              if (mpg%fcReg(iFc).ne.12) cycle
            case ( GEOMETRY_CDN )
              if (mpg%fcReg(iFc).ne.18 .and. mpg%fcReg(iFc).ne.25) cycle
            case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
              if (mpg%fcReg(iFc).ne.19 .and. mpg%fcReg(iFc).ne.26) cycle
            case ( GEOMETRY_LFS_SNOWFLAKE_MINUS,
     &             GEOMETRY_LFS_SNOWFLAKE_PLUS )
              if (mpg%fcReg(iFc).ne.19) cycle
            end select
            shidu(iCv,:) = 0.0_R8
          endif
        enddo
      endif
#ifndef FIXED_POINT
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_divua0to3 = shidu
      endif
#endif
      if (switch%tn_style.lt.2) then
        call b2saxpy (nCv*4, 1.0_R8, shidu, 1, shi0, 1)
      else
        call b2saxpy (nCv*4, 1.0_R8, shidu , 1, shi0, 1)
        call b2saxpy (nCv*4, 1.0_R8, shidun, 1, shn0, 1)
      endif

*   ..contribution from divergence of the electrical drift            !srv 02.06.99
      call gradc(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, dpoc)
      call gradc(nCv, nFc, nVx, 1, geo, mpg, geo%cvOnedBsq,
     &  geo%vxOnedBsq, dOnedBsqc)
      wrkc = geo%cvBb(:,2)*geo%cvVol*                                 !srv 11.01.13
     &    (dpoc(:,1)*dOnedBsqc(:,0) - dpoc(:,0)*dOnedBsqc(:,1))
      b2exbe=0.0_R8
      srw%b2sihs_exbe=0.0_R8
      do iCv = 1, mpg%nCi
!WG_TODO need better solution here...
        fac_ExB_c = fac_ExB(mpg%cvFc(mpg%cvFcP(iCv,1)))
        t0 = switch%b2sihs_phm5*fac_ExB_c*wrkc(iCv)                   !srv 07.06.00 02.06.99
        shedd(iCv,0) = shedd(iCv,0)+
     &    (max(0.0_R8,t0)+switch%b2sihs__rf0*abs(t0))*ne(iCv)*te(iCv)
        shedd(iCv,3) = shedd(iCv,3)+min(0.0_R8,t0)-
     &                 switch%b2sihs__rf0*abs(t0)
        b2exbe(mpg%cvReg(iCv))=b2exbe(mpg%cvReg(iCv))+
     1      t0*ne(iCv)*te(iCv)
        srw%b2sihs_exbe(iCv)=t0*ne(iCv)*te(iCv)
      enddo
      if (switch%b2sihs_shedd.ne.0) then
        do k=mpg%nCi+1,mpg%nCv
          iFc = mpg%cvFc(mpg%cvFcP(k,1))
          iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
          if (switch%b2sihs_shedd.eq.1) then
            shedd(iCv,:) = 0.0_R8
          elseif (switch%b2sihs_shedd.eq.2) then
            select case ( gridGeometry )
            case ( GEOMETRY_LINEAR )
              if (mpg%fcReg(iFc).eq.1 .and. itrg(1).ne.0 .or.
     &            mpg%fcReg(iFc).eq.2 .and. itrg(2).ne.0) cycle
            case ( GEOMETRY_LIMITER )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.2 .or.
     &            mpg%fcReg(iFc).eq.4) cycle
            case ( GEOMETRY_SN )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.8) cycle
            case ( GEOMETRY_CDN )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.5 .or. mpg%fcReg(iFc).eq.8 .or.
     &            mpg%fcReg(iFc).eq.14 .or. mpg%fcReg(iFc).eq.21) cycle
            case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.5 .or. mpg%fcReg(iFc).eq.8 .or.
     &            mpg%fcReg(iFc).eq.15 .or. mpg%fcReg(iFc).eq.22) cycle
            case ( GEOMETRY_LFS_SNOWFLAKE_MINUS,
     &             GEOMETRY_LFS_SNOWFLAKE_PLUS )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.8 .or. mpg%fcReg(iFc).eq.15) cycle
            end select
            shedd(iCv,:) = 0.0_R8
          elseif (switch%b2sihs_shedd.eq.3) then
            select case ( gridGeometry )
            case ( GEOMETRY_LINEAR )
              if (mpg%fcReg(iFc).ne.3) cycle
            case ( GEOMETRY_LIMITER )
              if (mpg%fcReg(iFc).ne.6) cycle
            case ( GEOMETRY_SN )
              if (mpg%fcReg(iFc).ne.12) cycle
            case ( GEOMETRY_CDN )
              if (mpg%fcReg(iFc).ne.18 .and. mpg%fcReg(iFc).ne.25) cycle
            case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
              if (mpg%fcReg(iFc).ne.19 .and. mpg%fcReg(iFc).ne.26) cycle
            case ( GEOMETRY_LFS_SNOWFLAKE_MINUS,
     &             GEOMETRY_LFS_SNOWFLAKE_PLUS )
              if (mpg%fcReg(iFc).ne.19) cycle
            end select
            shedd(iCv,:) = 0.0_R8
          endif
        enddo
      endif
#ifndef FIXED_POINT
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_exbe0to3 = shedd
      endif
#endif
      call b2saxpy (nCv*4, 1.0_R8, shedd, 1, she0, 1)                 !srv 11.09.09
*   ..contribution from divergence of the electrical drift            !srv 02.06.99                                               !srv 02.06.99
      b2exba=0.0_R8
      srw%b2sihs_exba=0.0_R8
      do is = 0, ns-1
       if(.not.is_neutral(is)) then
        do iCv = 1, mpg%nCi
!WG_TODO need better solution here...
           fac_ExB_c = fac_ExB(mpg%cvFc(mpg%cvFcP(iCv,1)))
           t0=switch%b2sihs_phm6*fac_ExB_c*
     &      (na(iCv,is)/ni(iCv,0))*wrkc(iCv)                          !srv 07.06.00 02.06.99
           shidd(iCv,0) = shidd(iCv,0)+
     &      (max(0.0_R8,t0)+switch%b2sihs__rf1*abs(t0))*ni(iCv,0)*
     &       ti(iCv)
!          shi0(iCv,0) = shi0(iCv,0)+
!     &      (max(0.0_R8,t0)+switch%b2sihs__rf1*abs(t0))*ni(iCv,0)*ti(iCv)
           shidd(iCv,3) = shidd(iCv,3)+min(0.0_R8,t0)-
     &                    switch%b2sihs__rf1*abs(t0)
!          shi0(iCv,3) = shi0(iCv,3)+min(0.0_R8,t0)-switch%b2sihs__rf1*abs(t0)
           b2exba(mpg%cvReg(iCv))=b2exba(mpg%cvReg(iCv))+
     1          t0*ni(iCv,0)*ti(iCv)
           srw%b2sihs_exba(iCv)=srw%b2sihs_exba(iCv)+
     1          t0*ni(iCv,0)*ti(iCv)
        enddo
       endif
      enddo
      if (switch%b2sihs_shidd.ne.0) then
        do k=mpg%nCi+1,mpg%nCv
          iFc = mpg%cvFc(mpg%cvFcP(k,1))
          iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
          if (switch%b2sihs_shidd.eq.1) then
            shidd(iCv,:) = 0.0_R8
          elseif (switch%b2sihs_shidd.eq.2) then
            select case ( gridGeometry )
            case ( GEOMETRY_LINEAR )
              if (mpg%fcReg(iFc).eq.1 .and. itrg(1).ne.0 .or.
     &            mpg%fcReg(iFc).eq.2 .and. itrg(2).ne.0) cycle
            case ( GEOMETRY_LIMITER )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.2 .or.
     &            mpg%fcReg(iFc).eq.4) cycle
            case ( GEOMETRY_SN )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.8) cycle
            case ( GEOMETRY_CDN )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.5 .or. mpg%fcReg(iFc).eq.8 .or.
     &            mpg%fcReg(iFc).eq.14 .or. mpg%fcReg(iFc).eq.21) cycle
            case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.5 .or. mpg%fcReg(iFc).eq.8 .or.
     &            mpg%fcReg(iFc).eq.15 .or. mpg%fcReg(iFc).eq.22) cycle
            case ( GEOMETRY_LFS_SNOWFLAKE_MINUS,
     &             GEOMETRY_LFS_SNOWFLAKE_PLUS )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.8 .or. mpg%fcReg(iFc).eq.15) cycle
            end select
            shidd(iCv,:) = 0.0_R8
          elseif (switch%b2sihs_shidd.eq.3) then
            select case ( gridGeometry )
            case ( GEOMETRY_LINEAR )
              if (mpg%fcReg(iFc).ne.3) cycle
            case ( GEOMETRY_LIMITER )
              if (mpg%fcReg(iFc).ne.6) cycle
            case ( GEOMETRY_SN )
              if (mpg%fcReg(iFc).ne.12) cycle
            case ( GEOMETRY_CDN )
              if (mpg%fcReg(iFc).ne.18 .and. mpg%fcReg(iFc).ne.25) cycle
            case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
              if (mpg%fcReg(iFc).ne.19 .and. mpg%fcReg(iFc).ne.26) cycle
            case ( GEOMETRY_LFS_SNOWFLAKE_MINUS,
     &             GEOMETRY_LFS_SNOWFLAKE_PLUS )
              if (mpg%fcReg(iFc).ne.19) cycle
            end select
            shidd(iCv,:) = 0.0_R8
          endif
        enddo
      endif
#ifndef FIXED_POINT
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_exba0to3 = shidd
      endif
#endif
      call b2saxpy (nCv*4, 1.0_R8, shidd, 1, shi0, 1)                 !srv 11.09.09
*   ..contribution from viscous heating
      if ((switch%Zhdanov_closure.eq.1) .and. 
     &    (switch%Zhdanov_test.eq.0) .and.
     &    switch%Zhdanov_visu.eq.1) then                              !som 18.08.21
*   ..re-calculate mass-average velocity
        call b2mod_zhfrtf_avm_u(nCv, ns, na, ua)
      endif
      b2visa=0.0_R8
      srw%b2sihs_visa=0.0_R8
      if (switch%min_collisions .gt. 0.0_R8) then
        call b2txnci (mpg, geo, switch, nCv, ns, ti, te, rz2,
     &   ne, ne2, na, lnlam, ismain, st_ext, wrkf,
     &   wrkc)
      else
        wrkc = 1.0_R8
      endif
      do is = 0, ns-1
       if(.not.is_neutral(is).or.switch%use_eirene.eq.0.or.           !srv 29.08.02
     &                           switch%spatial_hybrid.eq.1) then
        if (switch%style_visc.eq.0) then
          wrkf1 = 1.0_R8
        else
          !wdk use effective viscous momentum terms from the hybrid scheme
          !wdk avoids excessive heating in convection-dominated regions
          !wdk (improved total energy conservation at discretized level)
          wrkf1(:,0) = cvsahz_eff(:,0,is)/((cvsa(:,0,is)+eps)*geo%fcHz)
          wrkf1(:,1) = cvsahz_eff(:,1,is)/((cvsa(:,1,is)+eps)*geo%fcHz)
        endif
        if ((switch%Zhdanov_closure.eq.0) .or.
     &      (switch%Zhdanov_test.eq.1) .or.
     &      (switch%Zhdanov_visu.eq.0) .or.                           !som 16.08.21
     &      (switch%zh_vis_type.eq.0) .or.
     &      is_neutral(is) ) then
          call diff(nCv, nFc, nVx, 0, geo, mpg, ua(:,is), wrkv, wrkf)
          wrkf = cvsa(:,:,is)*wrkf1*wrkf**2/(geo%fcQalf+eps)
          call intcell(nFc, nCv, mpg, mpg%intcellP,wrkf(:,0),visc(:,0))
          call intcell(nFc, nCv, mpg, mpg%intcellR,wrkf(:,1),visc(:,1))
        else
         !wdk Case zh_vis_type .eq. 0 is now treated above.
         !wdk Should only get here if zh_vis_type .eq. 1
         ft = 4.0e0_R8/3.0e0_R8
         call diff(nCv, nFc, nVx, 0, geo, mpg, ua(:,is), wrkv, wrkf)
         wrkf(:,0) = 
     &     (cvsa(:,0,is) - ft*cvsa_cl(:,0,is))*wrkf1(:,0)*wrkf(:,0)**2/
     &     (geo%fcQalf(:,0)+eps)
         wrkf(:,1) = cvsa(:,1,is)*wrkf1(:,1)*wrkf(:,1)**2/
     &     (geo%fcQalf(:,1)+eps)
         call intcell(nFc, nCv, mpg, mpg%intcellP,wrkf(:,0),visc(:,0))
         call intcell(nFc, nCv, mpg, mpg%intcellR,wrkf(:,1),visc(:,1))
         call diff(nCv, nFc, nVx, 0, geo, mpg, avm_u, wrkv, wrkf)
         wrkf(:,0) = ft*cvsa_drho(:,0,is)*wrkf1(:,0)*wrkf(:,0)**2/
     &     (geo%fcQalf(:,0)+eps)
         call intcell(nFc, nCv, mpg, mpg%intcellP,wrkf(:,0),wrkc)
         visc(:,0) = visc(:,0) + wrkc
        endif
        do iCv = 1, mpg%nCi
          if (wrkc(iCv) .le. switch%min_collisions .and.
     &     .not.mpg%cvOnClosedSurface(iCv)) cycle
          t0 = switch%b2sihs_phm2*(1.0_R8-switch%BoRiS)*visc(iCv,0)
          t1 = switch%b2sihs_phm8*(1.0_R8-switch%BoRiS)*visc(iCv,1)
          if (switch%tn_style.eq.0) then
            shivc(iCv,0) = shivc(iCv,0)+
     &                     (1.0_R8+switch%b2sihs__rf2)*t0
            shivc(iCv,3) = shivc(iCv,3)-switch%b2sihs__rf2*t0/
     &                     (ni(iCv,0)*ti(iCv))
            shiva(iCv,0) = shiva(iCv,0)+
     &                     (1.0_R8+switch%b2sihs__rf2)*t1
            shiva(iCv,3) = shiva(iCv,3)-switch%b2sihs__rf2*t1/
     &                     (ni(iCv,0)*ti(iCv))
          elseif (switch%tn_style.eq.1) then
           if (.not.is_neutral(is)) then
            shivc(iCv,0) = shivc(iCv,0)+
     &                     (1.0_R8+switch%b2sihs__rf2)*t0
            shivc(iCv,3) = shivc(iCv,3)-switch%b2sihs__rf2*t0/
     &                     (ni(iCv,0)*ti(iCv))
            shiva(iCv,0) = shiva(iCv,0)+
     &                     (1.0_R8+switch%b2sihs__rf2)*t1
            shiva(iCv,3) = shiva(iCv,3)-switch%b2sihs__rf2*t1/
     &                     (ni(iCv,0)*ti(iCv))
           endif
          else
           if (.not.is_neutral(is).or.nint(zn(is)).ne.1) then
            shivc(iCv,0) = shivc(iCv,0)+
     &                     (1.0_R8+switch%b2sihs__rf2)*t0
            shivc(iCv,3) = shivc(iCv,3)-switch%b2sihs__rf2*t0/
     &                     (ni(iCv,0)*ti(iCv))
            shiva(iCv,0) = shiva(iCv,0)+
     &                     (1.0_R8+switch%b2sihs__rf2)*t1
            shiva(iCv,3) = shiva(iCv,3)-switch%b2sihs__rf2*t1/
     &                     (ni(iCv,0)*ti(iCv))
           else
            shivcn(iCv,0) = shivcn(iCv,0) +
     &                      (1.0_R8+switch%b2sihs__rf2)*t0
            shivcn(iCv,3) = shivcn(iCv,3) - switch%b2sihs__rf2*t0/
     &                      (nn(iCv)*tn(iCv))
            shivan(iCv,0) = shivan(iCv,0) +
     &                      (1.0_R8+switch%b2sihs__rf2)*t1
            shivan(iCv,3) = shivan(iCv,3) - switch%b2sihs__rf2*t1/
     &                      (nn(iCv)*tn(iCv))
           endif
          endif
          b2visa(mpg%cvReg(iCv))=b2visa(mpg%cvReg(iCv))+t0+t1
          srw%b2sihs_visa(iCv)=srw%b2sihs_visa(iCv)+t0+t1
        enddo
       endif                                                           !srv 29.08.02
      enddo

      if (switch%b2sihs_shivis.eq.0) then
        do k=mpg%nCi+1,mpg%nCv
          iFc = mpg%cvFc(mpg%cvFcP(k,1))
          iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
          if (switch%b2sihs_shivis.eq.1) then
            shivc(iCv,:) = 0.0_R8
            shiva(iCv,:) = 0.0_R8
          elseif (switch%b2sihs_shivis.eq.2) then
            select case ( gridGeometry )
            case ( GEOMETRY_LINEAR )
              if (mpg%fcReg(iFc).eq.1 .and. itrg(1).ne.0 .or.
     &            mpg%fcReg(iFc).eq.2 .and. itrg(2).ne.0) cycle
            case ( GEOMETRY_LIMITER )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.2 .or.
     &            mpg%fcReg(iFc).eq.4) cycle
            case ( GEOMETRY_SN )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.8) cycle
            case ( GEOMETRY_CDN )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.5 .or. mpg%fcReg(iFc).eq.8 .or.
     &            mpg%fcReg(iFc).eq.14 .or. mpg%fcReg(iFc).eq.21) cycle
            case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.5 .or. mpg%fcReg(iFc).eq.8 .or.
     &            mpg%fcReg(iFc).eq.15 .or. mpg%fcReg(iFc).eq.22) cycle
            case ( GEOMETRY_LFS_SNOWFLAKE_MINUS,
     &             GEOMETRY_LFS_SNOWFLAKE_PLUS )
              if (mpg%fcReg(iFc).eq.1 .or. mpg%fcReg(iFc).eq.4 .or.
     &            mpg%fcReg(iFc).eq.8 .or. mpg%fcReg(iFc).eq.15) cycle
            end select
            shivc(iCv,:) = 0.0_R8
            shiva(iCv,:) = 0.0_R8
          elseif (switch%b2sihs_shivis.eq.3) then
            select case ( gridGeometry )
            case ( GEOMETRY_LINEAR )
              if (mpg%fcReg(iFc).ne.3) cycle
            case ( GEOMETRY_LIMITER )
              if (mpg%fcReg(iFc).ne.6) cycle
            case ( GEOMETRY_SN )
              if (mpg%fcReg(iFc).ne.12) cycle
            case ( GEOMETRY_CDN )
              if (mpg%fcReg(iFc).ne.18 .and. mpg%fcReg(iFc).ne.25) cycle
            case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
              if (mpg%fcReg(iFc).ne.19 .and. mpg%fcReg(iFc).ne.26) cycle
            case ( GEOMETRY_LFS_SNOWFLAKE_MINUS,
     &             GEOMETRY_LFS_SNOWFLAKE_PLUS )
              if (mpg%fcReg(iFc).ne.19) cycle
            end select
            shivc(iCv,:) = 0.0_R8
            shiva(iCv,:) = 0.0_R8
          endif
        enddo
      endif

#ifndef FIXED_POINT
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_visa0to3 = shivc+shiva
      endif
#endif
      if (switch%tn_style.lt.2) then
       call b2saxpy (nCv*4, 1.0_R8, shivc , 1, shi0, 1)
       call b2saxpy (nCv*4, 1.0_R8, shiva , 1, shi0, 1)
      elseif (switch%tn_style.eq.2) then
       call b2saxpy (nCv*4, 1.0_R8, shivc , 1, shi0, 1)
       call b2saxpy (nCv*4, 1.0_R8, shiva , 1, shi0, 1)
       call b2saxpy (nCv*4, 1.0_R8, shivcn, 1, shn0, 1)
       call b2saxpy (nCv*4, 1.0_R8, shivan, 1, shn0, 1)
      endif

*   ..contribution from electron-ion and ion-ion friction in case of arbitrary impurity amount
*     in case of small impurity amount see below for b2sigp_style.ne.2
      if(switch%b2sigp_style.eq.2) then                                !srv 05.07.17
       call b2xzef (nCv, ns, rz2, na, ne, zeff, st_ext)                !srv 13.01.17 {
       do iCv = 1, nCv
         ce1(iCv)=fce1(zeff(iCv))
         ce2(iCv)=fce2(zeff(iCv))
       enddo 
       if (switch%Zhdanov_closure.eq.0.or.switch%Zhdanov_test.eq.1 .or.
     &     switch%Zhdanov_tf_fr.eq.0) then
         mm(0:ns-1) = switch%mass_multiplicator
         mm(ismain) = 1.0_R8
         am1(0:ns-1)= mm(0:ns-1)*am(0:ns-1)
         call b2xzmf (nCv, ns, ismain, rz2, na, zmff)
         call b2xzmc (nCv, ns, ismain, rz2, na, am1, zmc)
         call b2xzf1 (nCv, ns, ismain, rz2, na, am1, zf1)
         call b2xzf2 (nCv, ns, ismain, rz2, na, am1, zf2)
         call b2xz11 (nCv, ns, ismain, rz2, na, am1, z11)
         call b2xz12 (nCv, ns, ismain, rz2, na, am1, z12)
         call b2xz22 (nCv, ns, ismain, rz2, na, am1, z22)
         if (switch%impr_form_fr.eq.0) then
           do iCv = 1, nCv
             do is = 0, ns-1
               cimp1(iCv,is)=fce1(zmff(iCv))
             end do
           enddo
         else
           do is = 0, ns-1
             do iCv = 1, nCv
               cimp1(iCv,is)=fcim1(ismain,is,ns,am1,z11(iCv),
     &           z12(iCv),z22(iCv),zf1(iCv),zf2(iCv))
             enddo
           enddo
         endif
         if (switch%impr_form_tf.eq.0) then
           do iCv = 1, nCv
             do is = 0, ns-1
               cimp2(iCv,is)=fce2(zmff(iCv))
             end do
           enddo
         else
           do is = 0, ns-1
             do iCv = 1, nCv
               cimp2(iCv,is)=fcim2(ismain,is,ns,am1,zmc(iCv),
     &           z11(iCv),z12(iCv),z22(iCv))
             enddo
           enddo
         endif
       endif
*   ..compute the Coulomb logarithm
       call b2tlnl (nCv, switch, switch%icase_ii, te, ti, ne, lnlam)
*   ..compute zetap and zetae
       t1 = (4.0_R8*pi*eps0/qe)/qe
       t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
       czetap = t0*sqrt(mp)
       czetae = t0*sqrt(me)
       zetap = czetap*(ti/ev)*sqrt(ti/ev)/lnlam
       zetae = czetae*(te/ev)*sqrt(te/ev)/lnlam
*   ..compute gti,gte at cell centers
       call gradc_p(nCv, nFc, nVx, 0, geo, mpg, ti, wrkv, gti)
       call gradc_p(nCv, nFc, nVx, 0, geo, mpg, te, wrkv, gte)
       shifr = 0.0_R8
       shefr = 0.0_R8
       do is = 0,ns-1
         if(.not. is_neutral(is)) then
           call b2sifrtf (nCv, nFc, nVx, ns, is, ismain, switch, geo,
     &      mpg, na, ua, ue, te, ti, ni, ne, ne2, rza, rz2, po,
     &      zetap, zetae, gti, gte, ce1, ce2, cimp1, cimp2, zeff,
     &      f_luc_sg, alfx_c, sigx_c, st_ext,
     &      smbfrea, smbfria, smbtfea, smbtfia, smbch, smbtf)
           wrks =
     &      (smbfrea(:,0)+
     &       smbfrea(:,1)*ua(:,is)+
     &       smbfrea(:,2)*na(:,is)*mp*am(is)+
     &       smbfrea(:,3)*ua(:,is)*na(:,is)*mp*am(is))+
     &      (smbtfea(:,0)+
     &       smbtfea(:,1)*ua(:,is)+
     &       smbtfea(:,2)*na(:,is)*mp*am(is)+
     &       smbtfea(:,3)*ua(:,is)*na(:,is)*mp*am(is))
           shefr(:,0) = shefr(:,0) + wrks
           shifr(:,0) = shifr(:,0) + ua(:,is)*(wrks+
     &      (smbfria(:,0)+
     &       smbfria(:,1)*ua(:,is)+
     &       smbfria(:,2)*na(:,is)*mp*am(is)+
     &       smbfria(:,3)*ua(:,is)*na(:,is)*mp*am(is))+
     &      (smbtfia(:,0)+
     &       smbtfia(:,1)*ua(:,is)+
     &       smbtfia(:,2)*na(:,is)*mp*am(is)+
     &       smbtfia(:,3)*ua(:,is)*na(:,is)*mp*am(is))
     &                                               )
         endif
       enddo
*   ..add contribution due to Zhdanov multi-species approach
       if(switch%Zhdanov_closure.eq.1.and.switch%Zhdanov_test.eq.0.and.
     &    switch%Zhdanov_tf_fr.eq.1) then                             !som 01.07.21
         shifr(:,0) = shifr(:,0)-switch%zhd_corr*avm_u*shefr(:,0)
         shefr(:,0) = shefr(:,0)-switch%zhd_corr*(avm_u/ue)*shefr(:,0)
       endif
       shefr(:,0) =  ue*shefr(:,0)

       b2joule=0.0_R8
       srw%b2sihs_joule=0.0_R8
       b2fraa=0.0_R8
       srw%b2sihs_fraa=0.0_R8
       do iCv = 1, mpg%nCi
         if(mpg%cvReg(iCv).ne.0) then
           t0 = switch%b2sihs_phm3*geo%cvVol(iCv)*shefr(iCv,0)
!          shefr(iCv,0) = t0
           shefr(iCv,0) = max(0.0_R8,t0)+switch%b2sihs__rf3*abs(t0)
           shefr(iCv,3) =
     &      (min(0.0_R8,t0)-switch%b2sihs__rf3*abs(t0))/
     &      (ne(iCv)*te(iCv))
           b2joule(mpg%cvReg(iCv))=b2joule(mpg%cvReg(iCv))+t0
           srw%b2sihs_joule(iCv)=t0
           t0 = -switch%b2sihs_phm4*(1.0_R8-switch%BoRiS)*geo%cvVol(iCv)
     &          *shifr(iCv,0)
!          shifr(iCv,0) = t0
!           shifr(iCv,0) = (2.5_R8+switch%b2sihs__rf4)*t0
!           shifr(iCv,3) = -(1.5_R8+switch%b2sihs__rf4)*t0/(ni(iCv,0)*ti(iCv))
           shifr(iCv,0) = max(0.0_R8,t0)+switch%b2sihs__rf4*abs(t0)
           shifr(iCv,3) =
     &      (min(0.0_R8,t0)-switch%b2sihs__rf4*abs(t0))/
     &      (ni(iCv,0)*ti(iCv))
           b2fraa(mpg%cvReg(iCv))=b2fraa(mpg%cvReg(iCv))+t0
           srw%b2sihs_fraa(iCv)=t0
         endif
       enddo
      endif                                                           !srv 05.07.17

*   ..contribution from electron-ion friction
*     (a.k.a. Joule heating) in case of small impurity amount;
*     in case of arbitrary impurity amount see above for b2sigp_style.eq.2
      if(switch%b2sigp_style.ne.2) then                               !srv 05.07.17
       call b2xehx (nCv, nFc, nVx, switch, geo, mpg, te, po, ne, ehx) !srv 13.01.17
       if (switch%istyle_Joule_heating.eq.0)
     &   call b2xehy (nCv, nFc, nVx, switch, geo, mpg, te, po, ne, ehy)       !srv 13.01.17
       b2joule=0.0_R8
       srw%b2sihs_joule=0.0_R8
       do iFc = 1, nFc
         wrkf(iFc,0) = fch(iFc,0)*ehx(iFc)/
     &    ((geo%fcS(iFc)*geo%fcQalf(iFc,0) + eps)*
     &     (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2))))
       enddo
       call intcell(nFc, nCv, mpg, mpg%intcellP, wrkf(:,0), wrkc)
       if (switch%istyle_Joule_heating.eq.0) then
        do iFc = 1, nFc
          wrkf(iFc,1) = fch(iFc,1)*ehy(iFc)/
     &     ((geo%fcS(iFc)*geo%fcQalf(iFc,1) + eps)*
     &      (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2))))
        enddo
        call intcell(nFc, nCv, mpg, mpg%intcellR, wrkf(:,1), wrkc1)
        wrkc = wrkc + wrkc1
       endif
       do iCv = 1, mpg%nCi
         if(mpg%cvReg(iCv).ne.0) then
          t0 = switch%b2sihs_phm3*geo%cvVol(iCv)*wrkc(iCv)
          shefr(iCv,0) = shefr(iCv,0)+max(0.0_R8,t0)+
     &                   switch%b2sihs__rf3*abs(t0)
          shefr(iCv,3) = shefr(iCv,3)+
     &    (min(0.0_R8,t0)-switch%b2sihs__rf3*abs(t0))/
     &    (ne(iCv)*te(iCv))
          b2joule(mpg%cvReg(iCv))=b2joule(mpg%cvReg(iCv))+t0
          srw%b2sihs_joule(iCv)=t0
        endif                                                         !srv 26.11.02
       enddo
      endif                                                           !srv 05.07.17
#ifndef FIXED_POINT
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_joule0to3 = shefr
      endif
#endif
      call b2saxpy (nCv*4, 1.0_R8, shefr, 1, she0, 1)                 !srv 11.09.09
*   ..contribution from atom-atom friction
      if(switch%b2sigp_style.ne.2) then                               !srv 05.07.17
       b2fraa=0.0_R8
       srw%b2sihs_fraa=0.0_R8
       do isb = 0, ns-1
        do is = isb+1, ns-1
*     ..consider distinct cases
         if (.not.is_neutral(is).and..not.is_neutral(isb)) then
*      ..contribution from friction between charged species in case of small impurity amount;
*       in case of arbitrary impurity amount see above for b2sigp_style.eq.2
          coef = switch%b2sihs_phm4/ctaup*mp*
     &      sqrt(am(isb)*am(is)/(am(isb)+am(is)))                     !srv 01.07.09 11.01.12
          do iCv = 1,mpg%nCi
            t0 = rz2(iCv,isb)*rz2(iCv,is)*
     &          (coef*lnlam(iCv)/4.0_R8)*(1.0_R8-switch%BoRiS)*       !srv 11.01.12
     &           na(iCv,is)*na(iCv,isb)*vti32(iCv)*
     &           ((ua(iCv,is)-ua(iCv,isb))**2)                        !srv 22.07.05
            shifr(iCv,0) = shifr(iCv,0)+
     &                        (2.5_R8+switch%b2sihs__rf4)*t0
            shifr(iCv,3) = shifr(iCv,3)-
     &        (1.5_R8+switch%b2sihs__rf4)*t0/(ni(iCv,0)*ti(iCv))
            b2fraa(mpg%cvReg(iCv))=b2fraa(mpg%cvReg(iCv))+t0
            srw%b2sihs_fraa(iCv)=srw%b2sihs_fraa(iCv)+t0
          enddo
         else
*      ..contribution from friction involving neutral species
*        (We ignore such a term, except for charge exchange friction
*        between H0 and H1.d0 That term is calculated via the atomic
*        physics routines.)
         endif
        enddo
       enddo
      endif                                                           !srv 05.07.17
#ifndef FIXED_POINT
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_fraa0to3 = shifr
      endif
#endif
      call b2saxpy (nCv*4, 1.0_R8, shifr, 1, shi0, 1)                 !srv 11.09.09
      b2str=0.0_R8
      srw%b2sihs_str=0.0_R8
      b2wrong1=0.0_R8
* ..return
      if(switch%b2sihs_iout.eq.1) then
        do k=0,3
          write (chk,'(i1)') k
          call my_out_us(70,nCv,0,shedu(1,k),'b2sihs__shedu'//chk)
          call my_out_us(70,nCv,0,shedd(1,k),'b2sihs__shedd'//chk)
          call my_out_us(70,nCv,0,shefr(1,k),'b2sihs__shefr'//chk)
!         call my_out_us(70,nCv,0,she0(1,k),'b2sihs__she'//chk)

          call my_out_us(70,nCv,0,shidu(1,k),'b2sihs__shidu'//chk)
          call my_out_us(70,nCv,0,shidd(1,k),'b2sihs__shidd'//chk)
          call my_out_us(70,nCv,0,shifr(1,k),'b2sihs__shifr'//chk)
          call my_out_us(70,nCv,0,shiva(1,k),'b2sihs__shiva'//chk)
          call my_out_us(70,nCv,0,shivc(1,k),'b2sihs__shivc'//chk)
          call my_out_us(70,nCv,0,shivan(1,k),'b2sihs__shivan'//chk)
          call my_out_us(70,nCv,0,shivcn(1,k),'b2sihs__shivcn'//chk)
!         call my_out_us(70,nCv,0,shi0(1,k),'b2sihs__shi'//chk)
        enddo
      endif
      if (switch%b2npht_iout.eq.1.or.switch%iout_b2wdat.eq.4) then
        wrks=she0(:,0)+she0(:,1)*te+
     &       she0(:,2)*ne + she0(:,3)*te*ne
        call my_out_us(70,nCv,0,wrks,'b2sihs__she')
        wrks=shest(:,0)+shest(:,1)*te+
     &       shest(:,2)*ne + shest(:,3)*te*ne
        call my_out_us(70,nCv,0,shest,'b2sihs__shest')                !srv 17.12.13
        wrks=shedu(:,0)+shedu(:,1)*te+
     &       shedu(:,2)*ne + shedu(:,3)*te*ne
        call my_out_us(70,nCv,0,wrks,'b2sihs__shedu')
        wrks=shedd(:,0)+shedd(:,1)*te+
     &       shedd(:,2)*ne + shedd(:,3)*te*ne
        call my_out_us(70,nCv,0,wrks,'b2sihs__shedd')
        wrks=shefr(:,0)+shefr(:,1)*te+
     &  shefr(:,2)*ne + shefr(:,3)*te*ne
        call my_out_us(70,nCv,0,wrks,'b2sihs__shefr')
        wrks=shi0(:,0)+shi0(:,1)*ti+
     &       shi0(:,2)*ni(:,0)+shi0(:,3)*ti*ni(:,0)
        call my_out_us(70,nCv,0,wrks,'b2sihs__shi')
        wrks=shidu(:,0)+shidu(:,1)*ti+
     &       shidu(:,2)*ni(:,0)+shidu(:,3)*ti*ni(:,0)
        call my_out_us(70,nCv,0,wrks,'b2sihs__shidu')
        wrks=shidd(:,0)+shidd(:,1)*ti+
     &       shidd(:,2)*ni(:,0)+shidd(:,3)*ti*ni(:,0)
        call my_out_us(70,nCv,0,wrks,'b2sihs__shidd')
        wrks=shifr(:,0)+shifr(:,1)*ti+
     &       shifr(:,2)*ni(:,0)+shifr(:,3)*ti*ni(:,0)
        call my_out_us(70,nCv,0,wrks,'b2sihs__shifr')
        if (switch%tn_style.le.1) then
          wrks=shiva(:,0)+shiva(:,1)*ti+
     &         shiva(:,2)*ni(:,0)+shiva(:,3)*ti*ni(:,0)
          call my_out_us(70,nCv,0,wrks,'b2sihs__shiva')
          wrks=shivc(:,0)+shivc(:,1)*ti+
     &         shivc(:,2)*ni(:,0)+shivc(:,3)*ti*ni(:,0)
          call my_out_us(70,nCv,0,wrks,'b2sihs__shivc')
          wrks=shivan(:,0)+shivan(:,1)*ti+
     &         shivan(:,2)*ni(:,0)+shivan(:,3)*ti*ni(:,0)
          call my_out_us(70,nCv,0,wrks,'b2sihs__shivan')
          wrks=shivcn(:,0)+shivcn(:,1)*ti+
     &         shivcn(:,2)*ni(:,0)+shivcn(:,3)*ti*ni(:,0)
          call my_out_us(70,nCv,0,wrks,'b2sihs__shivcn')
        else
          wrks=shiva(:,0)+shiva(:,1)*ti+
     &         shiva(:,2)*ni(:,0)+shiva(:,3)*ti*ni(:,0)
          call my_out_us(70,nCv,0,wrks,'b2sihs__shiva')
          wrks=shivc(:,0)+shivc(:,1)*ti+
     &         shivc(:,2)*ni(:,0)+shivc(:,3)*ti*ni(:,0)
          call my_out_us(70,nCv,0,wrks,'b2sihs__shivc')
          wrks=shivan(:,0)+shivan(:,1)*tn+
     &         shivan(:,2)*(ni(:,0)-ni(:,1))+
     &         shivan(:,3)*tn*(ni(:,0)-ni(:,1))
          call my_out_us(70,nCv,0,wrks,'b2sihs__shivan')
          wrks=shivcn(:,0)+shivcn(:,1)*tn+
     &         shivcn(:,2)*(ni(:,0)-ni(:,1))+
     &         shivcn(:,3)*tn*(ni(:,0)-ni(:,1))
          call my_out_us(70,nCv,0,wrks,'b2sihs__shivcn')
        endif
      endif
      ncall_b2sihs_ = ncall_b2sihs_+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sihs_

      end subroutine b2sihs_
