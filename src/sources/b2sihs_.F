*-----------------------------------------------------------------------
*.specification

      subroutine b2sihs_ (nCv, nFc, nVx, ns, geo, mpg, BoRiS, ismain, 
     &  fac_ExB, lnlam, na, ua, ue, te, ti, po, ne, ni,
     &  fne, fch, cvsa, she0, shi0)
      use b2mod_types
!      use b2mod_rates
!      use b2mod_sources
      use b2mod_tallies
!      use b2mod_geo_corner
      use b2mod_anomalous_transport
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
!WG_TODO      use b2mod_balance !djm Jan2017
!WG_TODO     & , only : b2sihs_diae0to3, b2sihs_divue0to3, b2sihs_exbe0to3,
!WG_TODO     &          b2sihs_joule0to3, b2sihs_diaa0to3, b2sihs_divua0to3,
!WG_TODO     &          b2sihs_exba0to3, b2sihs_visa0to3, b2sihs_fraa0to3,
!WG_TODO     &          balance_netcdf
      use b2mod_transport_fun                                         !srv 05.07.17
!WG_TODO      use b2mod_external                                              !srv 05.07.17
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain                               !srv 05.07.17
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) :: BoRiS,
     *  fac_ExB(nFc), na(nCv,0:ns-1), ua(nCv,0:ns-1), ue(nCv),        !srv 05.03.11
     &  te(nCv), ti(nCv), po(nCv), ne(nCv), ni(nCv,0:1), lnlam(nCv)   !srv 07.05.08 01.07.09
      real (kind=R8) ::
     *  fne(nFc,0:1), fch(nFc,0:1), cvsa(nFc,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   she0(nCv,0:3), shi0(nCv,0:3)
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
      integer ncall, style, iCv, iFc, is, isb, k
      integer, save :: iout = 0, iout_b2npht = 0,                        !srv 11.09.09
     &                 iout_b2wdat = 0
      integer, save :: icase = 0
      integer, save :: istyle_Joule_heating = 1                         !srv 17.01.14
      integer, save :: b2sigp_style = 2                                 !srv 05.07.17
      integer, save :: use_eirene = 0
      real (kind=R8), save :: eps = 1.0e-60_R8
      real (kind=R8) ::
     *  vti32(nCv)
      real (kind=R8) ::
     *  rf0, rf1, rf2, rf3, rf4, phm0, phm1, phm2, phm3,
     &  phm4, phm5, phm6, wrkc(nCv), wrkf(nFc,0:1), wrkv(nVx),      !srv 02.06.99
     &  phm7, phm8, wrkc1(nCv),                                  !srv 28.08.02
     &  dpoc(nCv,0:1), dOnedBsqc(nCv,0:1), visc(nCv,0:1)
      real (kind=R8), save :: fch_stochastic = 0.0e0_R8                 !lk 31.01.08 {
      real (kind=R8) :: csig_stoch(nFc), fchstoch(nCv)  !lk 31.01.08 }
      real (kind=R8) :: 
     *  shidu(nCv,0:3), shidd(nCv,0:3),               !srv 11.09.09 {
     *  shiva(nCv,0:3), shivc(nCv,0:3),
     *  shifr(nCv,0:3), shedu(nCv,0:3),
     *  shedd(nCv,0:3), shefr(nCv,0:3),
     *  shest(nCv,0:3),                          !srv 17.12.13
     *  wrks(nCv), ehx(nFc), ehy(nFc)   !srv 13.01.17
      real (kind=R8) ::
     3  zeff(nCv), zmff(nCv), ce1(nCv),                                 !srv 13.01.17 05.07.17 {
     4  ce2(nCv), cimp1(nCv), cimp2(nCv),
     5  zetap(nCv), zetae(nCv),
     6  gti(nFc), gte(nFc), ne2(nCv),
     7  smbch(nCv,0:3), smbtf(nCv,0:3),
     *  smbfria(nCv,0:3), smbfrea(nCv,0:3),
     *  smbtfia(nCv,0:3), smbtfea(nCv,0:3)                              !srv 13.01.17 05.07.17 }
      character chk*1                                                   !srv 11.09.09 } !srv 05.07.17
      save ncall, style, rf0, rf1, rf2, rf3, rf4, phm0, phm1, phm2,
     &  phm3, phm4, phm5, phm6, phm7, phm8                              !srv 02.06.99 27.12.01 28.08.02
      real (kind=R8) ::
     *   t0, t1, ctaup, coef, czetap, czetae                            !srv 05.07.17
*   ..procedures
      intrinsic min, max, sqrt
      external subini, subend, xertst, ipgeti, ipgetr, sfill
      external b2xvff, b2xvfx, b2xvfv, b2xvsg
*   ..initialization
      data ncall/0/, style/0/
      data rf0/1.0e0_R8/, rf1/1.0e0_R8/, rf2/1.0e0_R8/, rf3/1.0e0_R8/,
     &   rf4/1.0e0_R8/, phm0/1.0e0_R8/, phm1/1.0e0_R8/, phm2/1.0e0_R8/,
     &   phm3/1.0e0_R8/, phm4/1.0e0_R8/, phm5/1.0e0_R8/, phm6/1.0e0_R8/,!srv 02.06.99
     &   phm7/0.0e0_R8/, phm8/1.0e0_R8/                                 !srv 27.12.01 28.08.02
*     Note 0: style determines the form of the strange electron-atom
*     energy transfer term.
*     Note 1: rf0-5 are numerical parameters, see the code. All must
*     be non-negative and of order unity; larger values are meant to
*     provide more stabilization. Some or all may be set to zero and
*     removed eventually. I've had some experience to indicate that
*     rf0, rf2, rf4 may have a useful function; however, later changes
*     in the code may have made all these parameters superfluous.

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sihs_')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2sihs_style', style)
       call ipgetr ('b2sihs_rf0', rf0)
       call ipgetr ('b2sihs_rf1', rf1)
       call ipgetr ('b2sihs_rf2', rf2)
       call ipgetr ('b2sihs_rf3', rf3)
       call ipgetr ('b2sihs_rf4', rf4)
       call ipgetr ('b2sihs_phm0', phm0)
       call ipgetr ('b2sihs_phm1', phm1)
       call ipgetr ('b2sihs_phm2', phm2)
       call ipgetr ('b2sihs_phm3', phm3)
       call ipgetr ('b2sihs_phm4', phm4)
       call ipgetr ('b2sihs_phm5', phm5)                           !srv 02.06.99
       call ipgetr ('b2sihs_phm6', phm6)                           !srv 02.06.99
       call ipgetr ('b2sihs_phm7', phm7)                           !srv 27.12.01
       call ipgetr ('b2sihs_phm8', phm8)                           !srv 28.08.02
       call xertst (0.le.style, 'faulty internal parameter style')
       call xertst (0.le.min(rf0,rf1,rf2,rf3,rf4),
     &   'faulty internal parameter rf0/1/2/3/4')
       call xertst (0.le.min(phm0,phm1,phm2,phm3,phm4,phm5,phm6),  !srv 02.06.99
     &   'faulty internal parameter phm0/1/2/3/4/5/6')             !srv 02.06.99
       call ipgeti ('b2sihs_iout', iout)                           !srv 10.09.03
       call ipgeti ('b2npht_iout', iout_b2npht)                    !srv 11.09.09
       call ipgeti ('b2wdat_iout', iout_b2wdat)                    !iys 17.04.17
       call ipgeti ('b2tlnl_ii', icase)
       call xertst (0.eq.icase .or. 1.eq.icase, 
     &               'faulty input b2tlnl_ii')
       call ipgetr ('b2tfhe_stochastic', fch_stochastic)           !lk 31.01.08
       call ipgeti ('b2sihs_istyle_Joule_heating', 
     &                      istyle_Joule_heating)                  !srv 17.01.14
       call ipgeti ('b2sigp_style', b2sigp_style)                  !srv 05.07.17
       call ipgeti ('b2mndr_eirene', use_eirene)
      endif
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
*    ..test sign of na, te, ti, ne, ni
       call b2xvsg (nCv*ns, na, 1, 'na', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, ni, 1, 'ni', '.gt.')
      endif
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
      do iCv = 1, nCv
        vti32(iCv) = geo%cvVol(iCv)/((ti(iCv)/ev)*sqrt(ti(iCv)/ev))
      enddo

* ..compute heat source terms
*   ..initialize sources to zero
      call sfill (nCv*4, 0.0e0_R8, she0, 1)
      call sfill (nCv*4, 0.0e0_R8, shi0, 1)
      call sfill (nCv*4, 0.0e0_R8, shedu, 1)                         !srv 11.09.09 {
      call sfill (nCv*4, 0.0e0_R8, shedd, 1)
      call sfill (nCv*4, 0.0e0_R8, shefr, 1)
      call sfill (nCv*4, 0.0e0_R8, shidu, 1)
      call sfill (nCv*4, 0.0e0_R8, shidd, 1)
      call sfill (nCv*4, 0.0e0_R8, shifr, 1)
      call sfill (nCv*4, 0.0e0_R8, shiva, 1)
      call sfill (nCv*4, 0.0e0_R8, shivc, 1)                         !srv 11.09.09 }
      call sfill (nCv*4, 0.0e0_R8, shest, 1)                         !srv 17.12.13
*   ..compute the Coulomb logarithm
      call b2tlnl (nCv, icase, te, ti, ne, lnlam)                !srv 01.07.09  20.09.11
*   ..contribution from stochastic electron flow
      if(fch_stochastic.ne.0.0e0_R8) then                           !lk 31.01.08 { !srv 17.12.13
       call xerrab ('b2sihs_: fch_stochastic not yet available for WG')
#ifdef WG_TODO
       call b2tstch (nx, ny, vol, hy, qz, bb, po, te, ne, csig_stoch,    !srv 07.07.14
     &      fchstoch)
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
     &    /(vol(ix,iy) + vol(topix(ix,iy),topiy(ix,iy))))/4.0e0_R8
         endif
        enddo
       enddo
#endif
      endif                                                         !lk 31.01.08 }
      call b2saxpy (nCv*4, 1.0e0_R8, shest, 1, she0, 1)                   !srv 17.12.13
*   ..contribution from divergence ue
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,ue,wrkf(1,0))
      wrkf(:,0) = wrkf(:,0)*geo%fcPbs
      call div0(nCv, nFc, mpg, wrkf(1,0), wrkc)
!WG_TODO      b2divue=0.0_R8
!WG_TODO      b2sihs_divue=0.0_R8
      do iCv = 1, mpg%nCi
         t0 = -phm0*(1.0_R8-BoRiS)*wrkc(iCv)
         shedu(iCv,0) = shedu(iCv,0)+
     &    (max(0.0e0_R8,t0)+rf0*abs(t0))*ne(iCv)*te(iCv)
         shedu(iCv,3) = shedu(iCv,3)+min(0.0e0_R8,t0)-rf0*abs(t0)
!WG_TODO         b2divue(mpg%cvReg(iCv))=b2divue(mpg%cvReg(iCv))+
!WG_TODO     1     t0*ne(iCv)*te(iCv)
!WG_TODO         b2sihs_divue(iCv)=t0*ne(iCv)*te(iCv)
      enddo
!WG_TODO      if (balance_netcdf.ne.0) then !djm Jan2017
!WG_TODO        b2sihs_divue0to3 = shedu
!WG_TODO      endif
      call b2saxpy (nCv*4, 1.0e0_R8, shedu, 1, she0, 1)                !srv 11.09.09
*   ..contribution from divergence ua
!WG_TODO      b2divua=0.0_R8
!WG_TODO      b2sihs_divua=0.0_R8
      do is = 0, ns-1
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,ua(:,is),wrkf(1,0))
        wrkf(:,0) = wrkf(:,0)*geo%fcPbs
        call div0(nCv, nFc, mpg, wrkf(1,0), wrkc)
        do iCv = 1, mpg%nCi
          t0 = -phm1*(na(iCv,is)/ni(iCv,0))*(1.0_R8-BoRiS)*wrkc(iCv)  !srv 21.07.05 {
          shidu(iCv,0) = shidu(iCv,0)+
     &      (max(0.0e0_R8,t0)+rf1*abs(t0))*ni(iCv,0)*ti(iCv)
          shidu(iCv,3) = shidu(iCv,3)+min(0.0e0_R8,t0)-rf1*abs(t0)
!WG_TODO          b2divua(mpg%cvReg(iCv))=b2divua(mpg%cvReg(iCv))
!WG_TODO     1         +t0*ni(iCv,0)*ti(iCv)
!WG_TODO          b2sihs_divua(iCv)=b2sihs_divua(iCv)+
!WG_TODO     1         t0*ni(iCv,0)*ti(iCv)
        enddo
      enddo
!WG_TODO      if (balance_netcdf.ne.0) then !djm Jan2017
!WG_TODO        b2sihs_divua0to3 = shidu
!WG_TODO      endif
      call b2saxpy (nCv*4, 1.0e0_R8, shidu, 1, shi0, 1)                !srv 11.09.09
*   ..contribution from divergence of the electrical drift            !srv 02.06.99
      call grad(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, wrkf)
      call intcell(nFc, nCv, mpg, mpg%intcellP, wrkf(1,0),dpoc(1,0))
      call intcell(nFc, nCv, mpg, mpg%intcellR, wrkf(1,1),dpoc(1,1))
      call grad(nCv, nFc, nVx, 0, geo, mpg, geo%cvOnedBsq, wrkv, wrkf)
      call intcell(nFc, nCv, mpg, mpg%intcellP,wrkf(1,0),dOnedBsqc(1,0))
      call intcell(nFc, nCv, mpg, mpg%intcellR,wrkf(1,1),dOnedBsqc(1,1))
      wrkc = geo%cvBb(:,2)*geo%cvVol*                                  !srv 11.01.13
     &    (dpoc(:,1)*dOnedBsqc(:,0) - dpoc(:,0)*dOnedBsqc(:,1))
!WG_TODO      b2exbe=0.0_R8
!WG_TODO      b2sihs_exbe=0.0_R8
      do iCv = 1, mpg%nCi
        t0 = phm5*fac_ExB(iCv)*wrkc(iCv)                         !srv 07.06.00 02.06.99
        shedd(iCv,0) = shedd(iCv,0)+
     &    (max(0.0e0_R8,t0)+rf0*abs(t0))*ne(iCv)*te(iCv)
        shedd(iCv,3) = shedd(iCv,3)+min(0.0e0_R8,t0)-rf0*abs(t0)
!WG_TODO        b2exbe(mpg%cvReg(iCv))=b2exbe(mpg%cvReg(iCv))+
!WG_TODO     1      t0*ne(iCv)*te(iCv)
!WG_TODO        b2sihs_exbe(iCv)=t0*ne(iCv)*te(iCv)
      enddo
!WG_TODO      if (balance_netcdf.ne.0) then !djm Jan2017
!WG_TODO        b2sihs_exbe0to3 = shedd
!WG_TODO      endif
      call b2saxpy (nCv*4, 1.0e0_R8, shedd, 1, she0, 1)                !srv 11.09.09
*   ..contribution from divergence of the electrical drift            !srv 02.06.99                                               !srv 02.06.99
!WG_TODO      b2exba=0.0_R8
!WG_TODO      b2sihs_exba=0.0_R8
      do is = 0, ns-1
       if(.not.is_neutral(is)) then
        do iCv = 1, mpg%nCi
           t0=phm6*fac_ExB(iCv)*
     &      (na(iCv,is)/ni(iCv,0))*wrkc(iCv)                    !srv 07.06.00 02.06.99
           shidd(iCv,0) = shidd(iCv,0)+
     &      (max(0.0e0_R8,t0)+rf1*abs(t0))*ni(iCv,0)*ti(iCv)
!          shi0(iCv,0) = shi0(iCv,0)+
!     &      (max(0.0e0_R8,t0)+rf1*abs(t0))*ni(iCv,0)*ti(iCv)
           shidd(iCv,3) = shidd(iCv,3)+min(0.0e0_R8,t0)-rf1*abs(t0)
!          shi0(iCv,3) = shi0(iCv,3)+min(0.0e0_R8,t0)-rf1*abs(t0)
!WG_TODO           b2exba(mpg%cvReg(iCv))=b2exba(mpg%cvReg(iCv))+
!WG_TODO     1          t0*ni(iCv,0)*ti(iCv)
!WG_TODO           b2sihs_exba(iCv)=b2sihs_exba(iCv)+
!WG_TODO     1          t0*ni(iCv,0)*ti(iCv)
        enddo
       endif
      enddo
!WG_TODO      if (balance_netcdf.ne.0) then !djm Jan2017
!WG_TODO        b2sihs_exba0to3 = shidd
!WG_TODO      endif
      call b2saxpy (nCv*4, 1.0e0_R8, shidd, 1, shi0, 1)                !srv 11.09.09
*   ..contribution from viscous heating
!WG_TODO      b2visa=0.0_R8
!WG_TODO      b2sihs_visa=0.0_R8
      do is = 0, ns-1
       if(.not.is_neutral(is).or.use_eirene.eq.0) then                !srv 29.08.02
        call diff(nCv, nFc, nVx, 0, geo, mpg, ua(1,is), wrkv, wrkf)
        wrkf(:,0) = cvsa(:,0,is)/(geo%fcQalf(iFc,0)+eps)*wrkf(:,0)**2
        wrkf(:,1) = cvsa(:,1,is)/(geo%fcQalf(iFc,1)+eps)*wrkf(:,1)**2
        call intcell(nFc, nCv, mpg, mpg%intcellP,wrkf(1,0),visc(1,0))
        call intcell(nFc, nCv, mpg, mpg%intcellR,wrkf(1,1),visc(1,1))
        do iCv = 1, mpg%nCi
           t0 = phm2*(1.0_R8-BoRiS)*visc(iCv,0)
           t1 = phm8*(1.0_R8-BoRiS)*visc(iCv,1)
           shivc(iCv,0) = shivc(iCv,0)+(1.0e0_R8+rf2)*t0
           shivc(iCv,3) = shivc(iCv,3)-
     &                        rf2*t0/(ni(iCv,0)*ti(iCv))
           shiva(iCv,0) = shiva(iCv,0)+(1.0e0_R8+rf2)*t1
           shiva(iCv,3) = shiva(iCv,3)-
     &                        rf2*t1/(ni(iCv,0)*ti(iCv))
!WG_TODO           b2visa(mpg%cvReg(iCv))=b2visa(mpg%cvReg(iCv))+t0+t1
!WG_TODO           b2sihs_visa(iCv)=b2sihs_visa(iCv)+t0+t1
        enddo
       endif                                                          !srv 29.08.02
      enddo
!WG_TODO      if (balance_netcdf.ne.0) then !djm Jan2017
!WG_TODO        b2sihs_visa0to3 = shivc+shiva
!WG_TODO      endif
      call b2saxpy (nCv*4, 1.0e0_R8, shivc, 1, shi0, 1)                !srv 11.09.09
      call b2saxpy (nCv*4, 1.0e0_R8, shiva, 1, shi0, 1)                !srv 11.09.09

*   ..contribution from electron-ion and ion-ion friction in case of arbitrary impurity amount
*     in case of small impurity amount see below for b2sigp_style.ne.2
      if(b2sigp_style.eq.2) then                                      !srv 05.07.17
       call xerrab('b2sihs_ -- b2sigp_style=2 not yet supported')
#ifdef WG_TODO
       call b2xpne (nCv, ns, rz2, na, st_ext%ne2, ne2)
       call b2xzef (nCv, ns, rz2, na, ne, zeff, st_ext)                    !srv 13.01.17 {
       call b2xzmf (nCv, ns, ismain, rz2, na, zmff)
       do iCv = 1, nCv
         ce1(iCv)=fce1(zeff(iCv))
         ce2(iCv)=fce2(zeff(iCv))
         cimp1(iCv)=fce1(zmff(iCv))
         cimp2(iCv)=fce2(zmff(iCv))
       enddo                                                          !srv 13.01.17 }
*   ..compute the Coulomb logarithm
       call b2tlnl (nCv, icase, te, ti, ne, lnlam)
*   ..compute zetap and zetae
       t1 = (4.0_R8*pi*eps0/qe)/qe
       t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
       czetap = t0*sqrt(mp)
       czetae = t0*sqrt(me)
       zetap = czetap*(ti/ev)*sqrt(ti/ev)/lnlam
       zetae = czetae*(te/ev)*sqrt(te/ev)/lnlam
*   ..compute gti,gte at cell faces
       call grad_p(nCv, nFc, nVx, 0, geo, mpg, ti, wrkv, gti)
       call grad_p(nCv, nFc, nVx, 0, geo, mpg, te, wrkv, gte)
       shifr = 0.0_R8
       shefr = 0.0_R8
       do is = 0,ns-1
         if(.not. is_neutral(is)) then
           call b2sifrtf (nx, ny, ns, is, ismain, hx, bb, na, ua,
     &      ue, te, ti, ni, ne, ne2,
     &      zetap, zetae, gti, gte, ce1, ce2, cimp1, cimp2,           !srv 05.07.17
     &      zeff, smbfrea, smbfria, smbtfea, smbtfia,                 !srv 05.07.17
     &      smbch, smbtf)
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
       shefr(:,0) =  ue*shefr(:,0)

       b2joule=0.0_R8
       b2sihs_joule=0.0_R8
       b2fraa=0.0_R8
       b2sihs_fraa=0.0_R8
       do iCv = 1, mpg%nCi
         if(mpg%cvReg(iCv).ne.0) then
           t0 = phm3*geo%cvVol(iCv)*shefr(iCv,0)
!          shefr(iCv,0) = t0
           shefr(iCv,0) = max(0.0e0_R8,t0)+rf3*abs(t0)
           shefr(iCv,3) =
     &      (min(0.0e0_R8,t0)-rf3*abs(t0))/(ne(iCv)*te(iCv))
           b2joule(mpg%cvReg(iCv))=b2joule(mpg%cvReg(iCv))+t0
           b2sihs_joule(iCv)=t0
           t0 = - phm4*(1.0_R8-BoRiS)*geo%cvVol(iCv)*shifr(iCv,0)
!          shifr(iCv,0) = t0
!           shifr(iCv,0) = (2.5e0_R8+rf4)*t0
!           shifr(iCv,3) = -(1.5e0_R8+rf4)*t0/(ni(iCv,0)*ti(iCv))
           shifr(iCv,0) = max(0.0e0_R8,t0)+rf4*abs(t0)
           shifr(iCv,3) =
     &      (min(0.0e0_R8,t0)-rf4*abs(t0))/(ni(iCv,0)*ti(iCv))
           b2fraa(mpg%cvReg(iCv))=b2fraa(mpg%cvReg(iCv))+t0
           b2sihs_fraa(iCv)=t0
         endif
       enddo
#endif
      endif                                                           !srv 05.07.17

*   ..contribution from electron-ion friction
*     (a.k.a. Joule heating) in case of small impurity amount;
*     in case of arbitrary impurity amount see above for b2sigp_style.eq.2
      if(b2sigp_style.ne.2) then                                      !srv 05.07.17
       call b2xehx (nCv, nFc, nVx, geo, mpg, te, po, ne, ehx)         !srv 13.01.17
       if (istyle_Joule_heating.eq.0)
     &   call b2xehy (nCv, nFc, nVx, geo, mpg, te, po, ne, ehy)       !srv 13.01.17
!WG_TODO       b2joule=0.0_R8
!WG_TODO       b2sihs_joule=0.0_R8
       do iFc = 1, nFc
         wrkf(:,0) = fch(iFc,0)*ehx(iFc)/
     &    ((geo%fcS(iFc)*geo%fcQalf(iFc,0) + eps)*
     &     (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2))))
       enddo
       call intcell(nFc,nCv,mpg, mpg%intcellP, wrkf(1,0),wrkc)
       if (istyle_Joule_heating.eq.0) then
        do iFc = 1, nFc
          wrkf(:,1) = fch(iFc,1)*ehy(iFc)/
     &     ((geo%fcS(iFc)*geo%fcQalf(iFc,1) + eps)*
     &      (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2))))
        enddo
        call intcell(nFc,nCv,mpg, mpg%intcellR, wrkf(1,1),wrkc1)
        wrkc = wrkc + wrkc1
       endif
       do iCv = 1, mpg%nCi
         if(mpg%cvReg(iCv).ne.0) then
          t0 = phm3*geo%cvVol(iCv)*wrkc(iCv)
          shefr(iCv,0) = shefr(iCv,0)+max(0.0e0_R8,t0)+rf3*abs(t0)
          shefr(iCv,3) = shefr(iCv,3)+
     &    (min(0.0e0_R8,t0)-rf3*abs(t0))/(ne(iCv)*te(iCv))
!WG_TODO          b2joule(mpg%cvReg(iCv))=b2joule(mpg%cvReg(iCv))+t0
!WG_TODO          b2sihs_joule(iCv)=t0
        endif                                                        !srv 26.11.02
       enddo
      endif                                                           !srv 05.07.17
!WG_TODO      if (balance_netcdf.ne.0) then !djm Jan2017
!WG_TODO        b2sihs_joule0to3 = shefr
!WG_TODO      endif
      call b2saxpy (nCv*4, 1.0e0_R8, shefr, 1, she0, 1)                !srv 11.09.09
*   ..contribution from atom-atom friction
#ifdef WG_TODO
      if(b2sigp_style.ne.2) then                                      !srv 05.07.17
       b2fraa=0.0_R8
       b2sihs_fraa=0.0_R8
       do isb = 0, ns-1
        do is = isb+1, ns-1
*     ..consider distinct cases
         if (.not.is_neutral(is).and..not.is_neutral(isb)) then
*      ..contribution from friction between charged species in case of small impurity amount;
*       in case of arbitrary impurity amount see above for b2sigp_style.eq.2
          coef = phm4/ctaup*mp*sqrt(am(isb)*am(is)/(am(isb)+am(is)))  !srv 01.07.09 11.01.12
          do iy = 0, ny-1
           do ix = 0, nx-1
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then !srv 26.11.02
             t0 = rz2(ix,iy,isb)*rz2(ix,iy,is)*
     &          (coef*lnlam(ix,iy)/4.0_R8)*(1.0_R8-BoRiS)*            !srv 11.01.12
     &           na(ix,iy,is)*na(ix,iy,isb)*vti32(ix,iy)*
     &           ((ua(ix,iy,is)-ua(ix,iy,isb))**2)                    !srv 22.07.05
             shifr(ix,iy,0) = shifr(ix,iy,0)+(2.5e0_R8+rf4)*t0
             shifr(ix,iy,3) = shifr(ix,iy,3)-
     &        (1.5e0_R8+rf4)*t0/(ni(ix,iy,0)*ti(ix,iy))
             b2fraa(region(ix,iy,0))=b2fraa(region(ix,iy,0))+t0
             b2sihs_fraa(ix,iy)=b2sihs_fraa(ix,iy)+t0
            endif                                                     !srv 26.11.02
           enddo
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
!WG_TODO      if (balance_netcdf.ne.0) then !djm Jan2017
!WG_TODO        b2sihs_fraa0to3 = shifr
!WG_TODO      endif
      call b2saxpy (nCv*4, 1.0e0_R8, shifr, 1, shi0, 1)                !srv 11.09.09
#endif
!WG_TODO      b2str=0.0_R8
!WG_TODO      b2sihs_str=0.0_R8
      b2wrong1=0.0_R8
* ..return
      if(iout.eq.1) then
        do k=0,3
          write (chk,'(i1)') k
          call my_out_us(70,nCv,shedu(1,k),'b2sihs__shedu'//chk)
          call my_out_us(70,nCv,shedd(1,k),'b2sihs__shedd'//chk)
          call my_out_us(70,nCv,shefr(1,k),'b2sihs__shefr'//chk)
!         call my_out_us(70,nCv,she0(1,k),'b2sihs__she'//chk)

          call my_out_us(70,nCv,shidu(1,k),'b2sihs__shidu'//chk)
          call my_out_us(70,nCv,shidd(1,k),'b2sihs__shidd'//chk)
          call my_out_us(70,nCv,shifr(1,k),'b2sihs__shifr'//chk)
          call my_out_us(70,nCv,shiva(1,k),'b2sihs__shiva'//chk)
          call my_out_us(70,nCv,shivc(1,k),'b2sihs__shivc'//chk)
!         call my_out_us(70,nCv,shi0(1,k),'b2sihs__shi'//chk)
        enddo
      endif
      if (iout_b2npht.eq.1 .or. iout_b2npht.eq.1 .or. iout_b2wdat.eq.4)
     &                                            then
        wrks=she0(:,0)+she0(:,1)*te+
     &       she0(:,2)*ne + she0(:,3)*te*ne
        call my_out_us(70,nCv,wrks,'b2sihs__she')
        wrks=shest(:,0)+shest(:,1)*te+
     &       shest(:,2)*ne + shest(:,3)*te*ne
        call my_out_us(70,nCv,shest,'b2sihs__shest')        !srv 17.12.13
        wrks=shedu(:,0)+shedu(:,1)*te+
     &       shedu(:,2)*ne + shedu(:,3)*te*ne
        call my_out_us(70,nCv,wrks,'b2sihs__shedu')
        wrks=shedd(:,0)+shedd(:,1)*te+
     &       shedd(:,2)*ne + shedd(:,3)*te*ne
        call my_out_us(70,nCv,wrks,'b2sihs__shedd')
        wrks=shefr(:,0)+shefr(:,1)*te+
     &  shefr(:,2)*ne + shefr(:,3)*te*ne
        call my_out_us(70,nCv,wrks,'b2sihs__shefr')

        wrks=shi0(:,0)+shi0(:,1)*ti+
     &       shi0(:,2)*ni(:,0)+shi0(:,3)*ti*ni(:,0)
        call my_out_us(70,nCv,wrks,'b2sihs__shi')
        wrks=shidu(:,0)+shidu(:,1)*ti+
     &       shidu(:,2)*ni(:,0)+shidu(:,3)*ti*ni(:,0)
        call my_out_us(70,nCv,wrks,'b2sihs__shidu')
        wrks=shidd(:,0)+shidd(:,1)*ti+
     &       shidd(:,2)*ni(:,0)+shidd(:,3)*ti*ni(:,0)
        call my_out_us(70,nCv,wrks,'b2sihs__shidd')
        wrks=shifr(:,0)+shifr(:,1)*ti+
     &       shifr(:,2)*ni(:,0)+shifr(:,3)*ti*ni(:,0)
        call my_out_us(70,nCv,wrks,'b2sihs__shifr')
        wrks=shiva(:,0)+shiva(:,1)*ti+
     &       shiva(:,2)*ni(:,0)+shiva(:,3)*ti*ni(:,0)
        call my_out_us(70,nCv,wrks,'b2sihs__shiva')
        wrks=shivc(:,0)+shivc(:,1)*ti+
     &       shivc(:,2)*ni(:,0)+shivc(:,3)*ti*ni(:,0)
        call my_out_us(70,nCv,wrks,'b2sihs__shivc')
      endif
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sihs_

      end
