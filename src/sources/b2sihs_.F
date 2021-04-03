*-----------------------------------------------------------------------
*.specification

      subroutine b2sihs_ (nCv, nFc, nVx, ns, switch, geo, mpg, BoRiS, 
     &  ismain, fac_ExB, lnlam, na, ua, ue, te, ti, tn, po, ne, ni,
     &  fna_fcor, fne, fch, cvsa, cvsahz_eff, she0, shi0, shn0)
      use b2mod_types
      use b2mod_tallies
      use b2mod_anomalous_transport
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
!WG_TODO      use b2mod_balance !djm Jan2017
!WG_TODO     & , only : b2sihs_diae0to3, b2sihs_divue0to3, b2sihs_exbe0to3,
!WG_TODO     &          b2sihs_joule0to3, b2sihs_diaa0to3, b2sihs_divua0to3,
!WG_TODO     &          b2sihs_exba0to3, b2sihs_visa0to3, b2sihs_fraa0to3,
!WG_TODO     &          balance_netcdf
      use b2mod_transport_fun                                         !srv 05.07.17
!WG_TODO      use b2mod_external                                              !srv 05.07.17
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2sihs_, ncall_b2xehy, ncall_b2xehx, ncall_b2tlnl
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain                               !srv 05.07.17
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) :: BoRiS,
     *  fac_ExB(nFc), na(nCv,0:ns-1), ua(nCv,0:ns-1), ue(nCv),        !srv 05.03.11
     &  te(nCv), ti(nCv), tn(nCv), po(nCv), ne(nCv), ni(nCv,0:1),   !srv 07.05.08 01.07.09
     &  lnlam(nCv) 
      real (kind=R8) ::
     *  fna_fcor(nFc,0:1,0:ns-1), fne(nFc,0:1), fch(nFc,0:1),
     *  cvsa(nFc,0:1,0:ns-1), cvsahz_eff(nFc,0:1,0:ns-1)
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
      real (kind=R8), parameter :: eps = 1.0e-60_R8
      real (kind=R8) ::
     *  vti32(nCv),
     &  wrkc(nCv), wrkf(nFc,0:1),
     &  wrkv(nVx),                      !srv 02.06.99
     &  wrkc1(nCv),                                         !srv 28.08.02
     &  dpoc(nCv,0:1), dOnedBsqc(nCv,0:1), visc(nCv,0:1),
     &  dumm0(nFc,0:1), dumm1(nFc,0:1)
      real (kind=R8) :: csig_stoch(nFc), fchstoch(nCv)  !lk 31.01.08 }
      real (kind=R8) :: 
     *  shidu(nCv,0:3), shidd(nCv,0:3),               !srv 11.09.09 {
     *  shiva(nCv,0:3), shivc(nCv,0:3),
     *  shivan(nCv,0:3), shivcn(nCv,0:3),
     *  shifr(nCv,0:3), shedu(nCv,0:3),
     *  shedd(nCv,0:3), shefr(nCv,0:3),
     *  shest(nCv,0:3),                          !srv 17.12.13
     *  shidun(nCv,0:3),
     *  wrks(nCv), ehx(nFc), ehy(nFc)   !srv 13.01.17i
      real (kind=R8) ::
     3  zeff(nCv), zmff(nCv), ce1(nCv),                                 !srv 13.01.17 05.07.17 {
     4  ce2(nCv), cimp1(nCv), cimp2(nCv),
     5  zetap(nCv), zetae(nCv),
     6  gti(nFc), gte(nFc), ne2(nCv),
     7  smbch(nCv,0:3), smbtf(nCv,0:3),
     *  smbfria(nCv,0:3), smbfrea(nCv,0:3),
     *  smbtfia(nCv,0:3), smbtfea(nCv,0:3)                              !srv 13.01.17 05.07.17 }
      character chk*1                                                   !srv 11.09.09 } !srv 05.07.17
      real (kind=R8) ::
     *   t0, t1, ctaup, coef, czetap, czetae                            !srv 05.07.17
      real (kind=R8) ::
     *   fac_ExB_c
*   ..procedures
      intrinsic min, max, sqrt
      external subini, subend, xertst, ipgeti, ipgetr, sfill
      external b2xvff, b2xvfx, b2xvfv, b2xvsg
*   ..initialization
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
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
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
      call sfill (nCv*4, 0.0e0_R8, shn0, 1)
      call sfill (nCv*4, 0.0e0_R8, shedu, 1)                         !srv 11.09.09 {
      call sfill (nCv*4, 0.0e0_R8, shedd, 1)
      call sfill (nCv*4, 0.0e0_R8, shefr, 1)
      call sfill (nCv*4, 0.0e0_R8, shidu, 1)
      call sfill (nCv*4, 0.0e0_R8, shidun, 1)
      call sfill (nCv*4, 0.0e0_R8, shidd, 1)
      call sfill (nCv*4, 0.0e0_R8, shifr, 1)
      call sfill (nCv*4, 0.0e0_R8, shiva, 1)
      call sfill (nCv*4, 0.0e0_R8, shivc, 1)                         !srv 11.09.09 }
      call sfill (nCv*4, 0.0e0_R8, shivan, 1)
      call sfill (nCv*4, 0.0e0_R8, shivcn, 1)                        !srv 11.09.09 }
      call sfill (nCv*4, 0.0e0_R8, shest, 1)                         !srv 17.12.13
*   ..compute the Coulomb logarithm
      call b2tlnl (nCv, switch, switch%icase_ii, te, ti, ne, lnlam)   !srv 01.07.09  20.09.11
*   ..contribution from stochastic electron flow
      if(switch%fch_stochastic.ne.0.0e0_R8) then                      !lk 31.01.08 { !srv 17.12.13
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
      if(switch%style_int_vel.eq.0) then
        call intface(nCv,nFc,mpg%fcCv,geo%fcHc,ue,wrkf(:,0))
      else
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,ue,wrkf(:,0))
      endif
      wrkf(:,0) = wrkf(:,0)*geo%fcPbs
      call div0(nCv, nFc, mpg, wrkf(:,0), wrkc)
!WG_TODO      b2divue=0.0_R8
!WG_TODO      b2sihs_divue=0.0_R8
      do iCv = 1, mpg%nCi
         t0 = -switch%b2sihs_phm0*(1.0_R8-BoRiS)*wrkc(iCv)
         shedu(iCv,0) = shedu(iCv,0)+
     &    (max(0.0e0_R8,t0)+switch%b2sihs__rf0*abs(t0))*ne(iCv)*te(iCv)
         shedu(iCv,3) = shedu(iCv,3)+min(0.0e0_R8,t0)-
     &                  switch%b2sihs__rf0*abs(t0)
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
        if(switch%style_int_vel.eq.0) then
          call intface(nCv,nFc,mpg%fcCv,geo%fcHc,ua(:,is),wrkf(:,0))
        else
          call intface(nCv,nFc,mpg%fcCv,geo%fcVol,ua(:,is),wrkf(:,0))
        endif
        wrkf(:,0) = wrkf(:,0)*geo%fcPbs
        call div0(nCv, nFc, mpg, wrkf(:,0), wrkc)
        do iCv = 1, mpg%nCi
          t0 = -switch%b2sihs_phm1*(na(iCv,is)/ni(iCv,0))
     &                                      *(1.0_R8-BoRiS)*wrkc(iCv)  !srv 21.07.05 {
          if (.not.is_neutral(is).or.switch%tn_style.eq.0) then
            shidu(iCv,0) = shidu(iCv,0)+
     &        (max(0.0e0_R8,t0)+switch%b2sihs__rf1*abs(t0))*ni(iCv,0)*
     &         ti(iCv)
            shidu(iCv,3) = shidu(iCv,3)+min(0.0e0_R8,t0)
     &                     -switch%b2sihs__rf1*abs(t0)
          else
            t0 = t0*tn(iCv)/ti(iCv) 
            shidun(iCv,0) = shidun(iCv,0)+
     &        (max(0.0e0_R8,t0)+switch%b2sihs__rf1*abs(t0))*ni(iCv,0)*
     &         ti(iCv)
            shidun(iCv,3) = shidun(iCv,3)+min(0.0e0_R8,t0)-
     &         switch%b2sihs__rf1*abs(t0)
          endif


!WG_TODO          b2divua(mpg%cvReg(iCv))=b2divua(mpg%cvReg(iCv))
!WG_TODO     1         +t0*ni(iCv,0)*ti(iCv)
!WG_TODO          b2sihs_divua(iCv)=b2sihs_divua(iCv)+
!WG_TODO     1         t0*ni(iCv,0)*ti(iCv)
        enddo
      enddo
!WG_TODO      if (balance_netcdf.ne.0) then !djm Jan2017
!WG_TODO        b2sihs_divua0to3 = shidu
!WG_TODO      endif

      if (.not.switch%use_eirene) then
        if (switch%tn_style.eq.0) then
          call b2saxpy (nCv*4, 1.0e0_R8, shidun, 1, shidu, 1)                !srv 11.09.09
          call b2saxpy (nCv*4, 1.0e0_R8, shidu, 1, shi0, 1)    
        else
          shidun(:,3) = shidun(:,3)*ni(:,0)*ti/(ni(:,0)-ni(:,1))/tn
      	  call b2saxpy (nCv*4, 1.0e0_R8, shidu, 1, shi0, 1)  
          call b2saxpy (nCv*4, 1.0e0_R8, shidun, 1, shn0, 1) 
        endif
      else
        call b2saxpy (nCv*4, 1.0e0_R8, shidu, 1, shi0, 1)
      endif  


*   ..contribution from divergence of the electrical drift            !srv 02.06.99
      call gradc(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, dpoc)
      call gradc(nCv, nFc, nVx, 1, geo, mpg, geo%cvOnedBsq,
     &  geo%vxOnedBsq, dOnedBsqc)
      wrkc = geo%cvBb(:,2)*geo%cvVol*                                  !srv 11.01.13
     &    (dpoc(:,1)*dOnedBsqc(:,0) - dpoc(:,0)*dOnedBsqc(:,1))
!WG_TODO      b2exbe=0.0_R8
!WG_TODO      b2sihs_exbe=0.0_R8
      do iCv = 1, mpg%nCi
!WG_TODO need better solution here...
        fac_ExB_c = fac_ExB(mpg%cvFc(mpg%cvFcP(iCv,1)))
        t0 = switch%b2sihs_phm5*fac_ExB_c*wrkc(iCv)                    !srv 07.06.00 02.06.99
        shedd(iCv,0) = shedd(iCv,0)+
     &    (max(0.0e0_R8,t0)+switch%b2sihs__rf0*abs(t0))*ne(iCv)*te(iCv)
        shedd(iCv,3) = shedd(iCv,3)+min(0.0e0_R8,t0)-
     &                 switch%b2sihs__rf0*abs(t0)
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
!WG_TODO need better solution here...
           fac_ExB_c = fac_ExB(mpg%cvFc(mpg%cvFcP(iCv,1)))
           t0=switch%b2sihs_phm6*fac_ExB_c*
     &      (na(iCv,is)/ni(iCv,0))*wrkc(iCv)                    !srv 07.06.00 02.06.99
           shidd(iCv,0) = shidd(iCv,0)+
     &      (max(0.0e0_R8,t0)+switch%b2sihs__rf1*abs(t0))*ni(iCv,0)*
     &       ti(iCv)
!          shi0(iCv,0) = shi0(iCv,0)+
!     &      (max(0.0e0_R8,t0)+switch%b2sihs__rf1*abs(t0))*ni(iCv,0)*ti(iCv)
           shidd(iCv,3) = shidd(iCv,3)+min(0.0e0_R8,t0)-
     &                    switch%b2sihs__rf1*abs(t0)
!          shi0(iCv,3) = shi0(iCv,3)+min(0.0e0_R8,t0)-switch%b2sihs__rf1*abs(t0)
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
       if(.not.is_neutral(is).or.switch%use_eirene.eq.0.or.            !srv 29.08.02
     &                           switch%spatial_hybrid.eq.1) then               
        if (switch%style_visc.eq.0) then
         call diff(nCv, nFc, nVx, 0, geo, mpg, ua(:,is), wrkv, wrkf)
         wrkf = cvsa(:,:,is)/(geo%fcQalf+eps)*wrkf**2
         call intcell(nFc, nCv, mpg, mpg%intcellP,wrkf(:,0),visc(:,0))
         call intcell(nFc, nCv, mpg, mpg%intcellR,wrkf(:,1),visc(:,1))
        else
         !wdk use effective viscous momentum terms from the hybrid scheme
         !wdk avoids excessive heating in convection-dominated regions
         !wdk (improved total energy conservation at discretized level)
         call diff(nCv, nFc, nVx, 0, geo, mpg, ua(:,is), wrkv, wrkf)
         wrkf = cvsahz_eff(:,:,is)/(geo%fcQalf+eps)*wrkf**2
         wrkf(:,0) = wrkf(:,0)/geo%fcHz
         wrkf(:,1) = wrkf(:,1)/geo%fcHz
         call intcell(nFc, nCv, mpg, mpg%intcellP,wrkf(:,0),visc(:,0))
         call intcell(nFc, nCv, mpg, mpg%intcellR,wrkf(:,1),visc(:,1))
        endif
        do iCv = 1, mpg%nCi
           t0 = switch%b2sihs_phm2*(1.0_R8-BoRiS)*visc(iCv,0)
           t1 = switch%b2sihs_phm8*(1.0_R8-BoRiS)*visc(iCv,1)
           if (.not.is_neutral(is)) then
             shivc(iCv,0) = shivc(iCv,0)+
     &                      (1.0e0_R8+switch%b2sihs__rf2)*t0
             shivc(iCv,3) = shivc(iCv,3)-switch%b2sihs__rf2*t0/
     &                      (ni(iCv,0)*ti(iCv))
             shiva(iCv,0) = shiva(iCv,0)+
     &                      (1.0e0_R8+switch%b2sihs__rf2)*t1
             shiva(iCv,3) = shiva(iCv,3)-switch%b2sihs__rf2*t1/
     &                      (ni(iCv,0)*ti(iCv))
           else
             shivcn(iCv,0) = shivcn(iCv,0) + 
     &                       (1.0e0_R8+switch%b2sihs__rf2)*t0
             shivcn(iCv,3) = shivcn(iCv,3) - switch%b2sihs__rf2*t0/
     &                       (ni(iCv,0)*ti(iCv))
             shivan(iCv,0) = shivan(iCv,0) + 
     &                       (1.0e0_R8+switch%b2sihs__rf2)*t1
             shivan(iCv,3) = shivan(iCv,3) - switch%b2sihs__rf2*t1/
     &                       (ni(iCv,0)*ti(iCv))
           endif
!WG_TODO           b2visa(mpg%cvReg(iCv))=b2visa(mpg%cvReg(iCv))+t0+t1
!WG_TODO           b2sihs_visa(iCv)=b2sihs_visa(iCv)+t0+t1
        enddo
       endif                                                           !srv 29.08.02
      enddo
!WG_TODO      if (balance_netcdf.ne.0) then !djm Jan2017
!WG_TODO        b2sihs_visa0to3 = shivc+shiva
!WG_TODO      endif
      if (switch%tn_style.eq.0) then
       call b2saxpy (nCv*4, 1.0e0_R8, shivcn, 1, shivc, 1)
       call b2saxpy (nCv*4, 1.0e0_R8, shivan, 1, shiva, 1)
       call b2saxpy (nCv*4, 1.0e0_R8, shivc , 1, shi0, 1)
       call b2saxpy (nCv*4, 1.0e0_R8, shiva , 1, shi0, 1)
      elseif (switch%tn_style.eq.1) then
       call b2saxpy (nCv*4, 1.0e0_R8, shivc, 1, shi0, 1)
       call b2saxpy (nCv*4, 1.0e0_R8, shiva, 1, shi0, 1)
      elseif (switch%tn_style.eq.2) then
       shivcn(:,3) = shivcn(:,3)*ni(:,0)*ti/(ni(:,0)-ni(:,1))/tn
       shivan(:,3) = shivan(:,3)*ni(:,0)*ti/(ni(:,0)-ni(:,1))/tn
       call b2saxpy (nCv*4, 1.0e0_R8, shivc, 1, shi0, 1)
       call b2saxpy (nCv*4, 1.0e0_R8, shiva, 1, shi0, 1)
       call b2saxpy (nCv*4, 1.0e0_R8, shivcn, 1, shn0, 1)
       call b2saxpy (nCv*4, 1.0e0_R8, shivan, 1, shn0, 1)
      endif

*   ..contribution from electron-ion and ion-ion friction in case of arbitrary impurity amount
*     in case of small impurity amount see below for b2sigp_style.ne.2
      if(switch%b2sigp_style.eq.2) then                                      !srv 05.07.17
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
       call b2tlnl (nCv, switch, swtich%icase_ii, te, ti, ne, lnlam)
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
           t0 = switch%b2sihs_phm3*geo%cvVol(iCv)*shefr(iCv,0)
!          shefr(iCv,0) = t0
           shefr(iCv,0) = max(0.0e0_R8,t0)+switch%b2sihs__rf3*abs(t0)
           shefr(iCv,3) =
     &      (min(0.0e0_R8,t0)-switch%b2sihs__rf3*abs(t0))/
     &      (ne(iCv)*te(iCv))
           b2joule(mpg%cvReg(iCv))=b2joule(mpg%cvReg(iCv))+t0
           b2sihs_joule(iCv)=t0
           t0 = -switch%b2sihs_phm4*(1.0_R8-BoRiS)*geo%cvVol(iCv)*
     &           shifr(iCv,0)
!          shifr(iCv,0) = t0
!           shifr(iCv,0) = (2.5e0_R8+switch%b2sihs__rf4)*t0
!           shifr(iCv,3) = -(1.5e0_R8+switch%b2sihs__rf4)*t0/(ni(iCv,0)*ti(iCv))
           shifr(iCv,0) = max(0.0e0_R8,t0)+switch%b2sihs__rf4*abs(t0)
           shifr(iCv,3) =
     &      (min(0.0e0_R8,t0)-switch%b2sihs__rf4*abs(t0))/
     &      (ni(iCv,0)*ti(iCv))
           b2fraa(mpg%cvReg(iCv))=b2fraa(mpg%cvReg(iCv))+t0
           b2sihs_fraa(iCv)=t0
         endif
       enddo
#endif
      endif                                                           !srv 05.07.17

*   ..contribution from electron-ion friction
*     (a.k.a. Joule heating) in case of small impurity amount;
*     in case of arbitrary impurity amount see above for b2sigp_style.eq.2
      if(switch%b2sigp_style.ne.2) then                               !srv 05.07.17
       call b2xehx (nCv, nFc, nVx, switch, geo, mpg, te, po, ne, ehx)         !srv 13.01.17
       if (switch%istyle_Joule_heating.eq.0)
     &   call b2xehy (nCv, nFc, nVx, switch, geo, mpg, te, po, ne, ehy)       !srv 13.01.17
!WG_TODO       b2joule=0.0_R8
!WG_TODO       b2sihs_joule=0.0_R8
       do iFc = 1, nFc
         wrkf(iFc,0) = fch(iFc,0)*ehx(iFc)/
     &    ((geo%fcS(iFc)*geo%fcQalf(iFc,0) + eps)*
     &     (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2))))
       enddo
       call intcell(nFc,nCv,mpg, mpg%intcellP, wrkf(:,0),wrkc)
       if (switch%istyle_Joule_heating.eq.0) then
        do iFc = 1, nFc
          wrkf(iFc,1) = fch(iFc,1)*ehy(iFc)/
     &     ((geo%fcS(iFc)*geo%fcQalf(iFc,1) + eps)*
     &      (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2))))
        enddo
        call intcell(nFc,nCv,mpg, mpg%intcellR, wrkf(:,1),wrkc1)
        wrkc = wrkc + wrkc1
       endif
       do iCv = 1, mpg%nCi
         if(mpg%cvReg(iCv).ne.0) then
          t0 = switch%b2sihs_phm3*geo%cvVol(iCv)*wrkc(iCv)
          shefr(iCv,0) = shefr(iCv,0)+max(0.0e0_R8,t0)+
     &                   switch%b2sihs__rf3*abs(t0)
          shefr(iCv,3) = shefr(iCv,3)+
     &    (min(0.0e0_R8,t0)-switch%b2sihs__rf3*abs(t0))/
     &    (ne(iCv)*te(iCv))
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
      if(switch%b2sigp_style.ne.2) then                                      !srv 05.07.17
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
             shifr(ix,iy,0) = shifr(ix,iy,0)+
     &                        (2.5e0_R8+switch%b2sihs__rf4)*t0
             shifr(ix,iy,3) = shifr(ix,iy,3)-
     &        (1.5e0_R8+switch%b2sihs__rf4)*t0/(ni(ix,iy,0)*ti(ix,iy))
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
        call my_out_us(70,nCv,0,shest,'b2sihs__shest')        !srv 17.12.13
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

      end
