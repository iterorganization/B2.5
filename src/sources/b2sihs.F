*-----------------------------------------------------------------------
*.specification

      subroutine b2sihs (nx, ny, ns, BoRiS, ismain, fac_ExB, hx, hy,     !xpb !srv 05.07.17
     1  vol, bb, pbs, lnlam, na, ua, ue, te, ti, po, ne, ni,             !xpb !srv 01.07.09
     2  fch, cvsa, fchdia, uadia, wadia, vaecrb, she0, shi0)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_tallies
      use b2mod_sources
      use b2mod_b2cmpa
      use b2mod_balance !djm Jan2017
     & , only : b2sihs_divue0to3, b2sihs_exbe0to3,
     &          b2sihs_joule0to3, b2sihs_divua0to3,
     &          b2sihs_exba0to3, b2sihs_visa0to3, b2sihs_fraa0to3,
     &          balance_netcdf
      use b2mod_transport_fun                                            !srv 05.07.17
      use b2mod_external                                                 !srv 05.07.17
      use b2mod_geo
     , , only : hy1
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain                                         !xpb !srv 05.07.15
      real (kind=R8) ::
     *  BoRiS, fac_ExB(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), 
     &  vol(-1:nx,-1:ny), pbs(-1:nx,-1:ny,0:1), na(-1:nx,-1:ny,0:ns-1),
     &  ua(-1:nx,-1:ny,0:ns-1), te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     &  po(-1:nx,-1:ny), ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     1  ue(-1:nx,-1:ny), lnlam(-1:nx,-1:ny)                              !srv 01.07.09
      real (kind=R8) ::
     *  fch(-1:nx,-1:ny,0:1), cvsa(-1:nx,-1:ny,0:1,0:ns-1),
     2  fchdia(-1:nx,-1:ny,0:1), uadia(-1:nx,-1:ny,0:1,0:ns-1),
     3  bb(-1:nx,-1:ny,0:3), wadia(-1:nx,-1:ny,0:1,0:ns-1), 
     4  vaecrb(-1:nx,-1:ny,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3)
*   ..workspace argument (unspecified on entry and on exit)
      real (kind=R8) :: shetmp(-1:nx,-1:ny,0:3), shitmp(-1:nx,-1:ny,0:3) !djm Jan2017
*   ..common block
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

*   ..local variables
      integer ncall, style, n2, ix, iy, is, isb, nsmin, nsmax
      integer, save :: icase = 0
      integer, save :: istyle_Joule_heating = 1                         !srv 17.01.14
      integer, save :: b2sigp_style = 2                                 !srv 05.07.17
      integer, save :: use_eirene = 0
      real (kind=R8), save :: fch_stochastic = 0.0e0_R8                 !lk 31.01.08 {
      real (kind=R8) :: csig_stoch(-1:nx,-1:ny), fchstoch(-1:nx,-1:ny)  !lk 31.01.08 }
      real (kind=R8) ::
     *  vti32(-1:nx,-1:ny)
      real (kind=R8) ::
     *  rf0, rf1, rf2, rf3, rf4, phm0, phm1, phm2, phm3,
     &  phm4, phm5, phm6, phm7, work(-1:nx,-1:ny), 
     &  wpov(-1:nx,-1:ny), wpoh(-1:nx,-1:ny),
     *  ehx(-1:nx,-1:ny), ehy(-1:nx,-1:ny)                              !srv 13.01.17
      real (kind=R8) ::
     3  zeff(-1:nx,-1:ny), zmff(-1:nx,-1:ny), ce1(-1:nx,-1:ny),         !srv 13.01.17 05.07.17 {
     4  ce2(-1:nx,-1:ny), cimp1(-1:nx,-1:ny), cimp2(-1:nx,-1:ny),
     5  zetap(-1:nx,-1:ny), zetae(-1:nx,-1:ny),
     6  gti(-1:nx,-1:ny), gte(-1:nx,-1:ny), ne2(-1:nx,-1:ny),
     &  shefr(-1:nx,-1:ny,0:3), shifr(-1:nx,-1:ny,0:3),                 !srv 05.07.17
     &  smbch(-1:nx,-1:ny,0:3), smbtf(-1:nx,-1:ny,0:3),
     &  smbfria(-1:nx,-1:ny,0:3), smbfrea(-1:nx,-1:ny,0:3),
     &  smbtfia(-1:nx,-1:ny,0:3), smbtfea(-1:nx,-1:ny,0:3),             !srv 13.01.17 05.07.17 }
     &  shest(-1:nx,-1:ny,0:3)                                          !srv 17.12.13
      save ncall, style, rf0, rf1, rf2, rf3, rf4, 
     &  phm0, phm1, phm2, phm3, phm4, phm5, phm6, phm7, nsmin, nsmax
      real (kind=R8) ::
     *  t0, t1, ctaup, coef, czetap, czetae                             !srv 05.07.17
*   ..procedures
      intrinsic min, max, sqrt
      real (kind=R8) ::
     *  uasx, uesx, uasy                                                !srv 13.01.17
      external subini, subend, xertst, ipgeti, ipgetr, sfill
      external b2xvff, b2xvsg, intfacev, intfaceh
      uasx(ix,iy,is)=(hx(ix,iy)*ua(leftix(ix,iy),leftiy(ix,iy),is)+
     1  hx(leftix(ix,iy),leftiy(ix,iy))*ua(ix,iy,is))/
     2 (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*pbs(ix,iy,0)
      uesx(ix,iy)=(hx(ix,iy)*ue(leftix(ix,iy),leftiy(ix,iy))+
     1  hx(leftix(ix,iy),leftiy(ix,iy))*ue(ix,iy))/
     2 (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*pbs(ix,iy,0)
      uasy(ix,iy,is)=(hy1(ix,iy)*ua(bottomix(ix,iy),bottomiy(ix,iy),is)+
     1  hy1(bottomix(ix,iy),bottomiy(ix,iy))*ua(ix,iy,is))/
     2 (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))*pbs(ix,iy,1)
*   ..initialisation
      data ncall/0/, style/0/
      data rf0/1.0_R8/, rf1/1.0_R8/, rf2/1.0_R8/,
     1  rf3/1.0_R8/, rf4/1.0_R8/,
     &  phm0/1.0_R8/, phm1/1.0_R8/, phm2/1.0_R8/,
     1  phm3/1.0_R8/, phm4/1.0_R8/, phm5/1.0_R8/,
     2  phm6/1.0_R8/, phm7/0.0_R8/
*     Note 0: style determines the form of the strange electron-atom
*     energy transfer term.
*     Note 1: rf0-5 are numerical parameters, see the code. All must
*     be non-negative and of order unity; larger values are meant to
*     provide more stabilisation. Some or all may be set to zero and
*     removed eventually. I have had some experience to indicate that
*     rf0, rf2, rf4 may have a useful function; however, later changes
*     in the code may have made all these parameters superfluous.

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sihs')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       nsmin=0
       nsmax=ns
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
       call ipgetr ('b2sihs_phm5', phm5)
       call ipgetr ('b2sihs_phm6', phm6)
       call ipgetr ('b2sihs_phm7', phm7)
       call xertst (1.ge.abs(style), 'faulty internal parameter style')
       call xertst (0.0_R8.le.min(rf0,rf1,rf2,rf3,rf4),
     &   'faulty internal parameter rf0/1/2/3/4')
       call xertst 
     .  (0.0_R8.le.min(phm0,phm1,phm2,phm3,phm4,phm5,phm6,phm7),
     &   'faulty internal parameter phm0/1/2/3/4/5/6/7')
       call ipgeti ('b2tlnl_ii', icase)
       call xertst (0.eq.icase .or. 1.eq.icase, 
     &               'faulty input b2tlnl_ii')
       call ipgetr ('b2tfhe_stochastic', fch_stochastic)           !lk 31.01.08
       call ipgeti ('b2sihs_istyle_Joule_heating', istyle_Joule_heating)  !srv 17.01.14
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
       call ipgeti ('b2sigp_style', b2sigp_style)                  !srv 05.07.17
       call ipgeti ('b2mndr_eirene', use_eirene)
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of hx, hy, vol
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
*    ..test sign of na, te, ti, ne, ni
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
*    ..test edge values of fch
       call b2xvff (nx, ny, fch, 'fch')
*    ..test edge values of cvsa
       do is = 0, ns-1
        call b2xvff (nx, ny, cvsa(-1,-1,0,is), 'cvsa')
       enddo
*    ..test fac_ExB
       call xertst (0.0_R8.le.minval(fac_ExB) .and.
     &       1.0_R8.ge.maxval(fac_ExB), 'faulty argument fac_ExB')
      endif
*   ..compute ctaup
*     Following Balescu; Braginskii has ctaup different by a factor
*     sqrt(2). The mathematical expression is:
*       t0 = (3/(4*sqrt(2*pi)))*(4*pi*eps0/qe**2)**2*ev*sqrt(ev)
*       ctaup = t0*sqrt(mp)/lnlam
*     The numerical value in SI units is: ctaup = 1.4743e13_R8/lnlam.
*     The following calculation of t0 is designed to avoid overflow
*     even for IEEE single precision arithmetic.
!xpb The lnlam term has been removed and is put back in later so it can
!xpb vary with space and species.
      t1 = (4.0_R8*pi*eps0/qe)/qe
      t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
      ctaup = t0*sqrt(mp)                                                !xpb
*   ..compute vti32
      do iy = -1, ny
       do ix = -1, nx
        vti32(ix,iy) = vol(ix,iy)/((ti(ix,iy)/ev)*sqrt(ti(ix,iy)/ev))
       enddo
      enddo

* ..compute heat source terms
*   ..initialise sources to zero
      call sfill (n2*4, 0.0_R8, she0, 1)
      call sfill (n2*4, 0.0_R8, shi0, 1)
      call sfill (n2*4, 0.0e0_R8, shest, 1)                        !srv 17.12.13
*   ..compute the Coulomb logarithm
      call b2tlnl (nx, ny, icase, te, ti, ne, lnlam)               !srv 01.07.09  20.09.11
*   ..contribution from stochastic electron flow
      if(fch_stochastic.ne.0.0e0_R8) then                          !lk 31.01.08 { !srv 17.12.13
       call b2tstch (nx, ny, vol, bb, po, te, ne, csig_stoch,      !srv 07.07.14
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
      endif                                                               !lk 31.01.08 }
      call b2saxpy (n2*4, 1.0e0_R8, shest, 1, she0, 1)                   !srv 17.12.13
*   ..contribution from divergence(ue,ve)
!xpb !!! The radial derivative terms uesy are non-conservative !!!
cxpb    t0 = -phm0*((uesx(rightix(ix,iy),rightiy(ix,iy))-uesx(ix,iy))    !xpb
cxpb     &         +(uesy(topix(ix,iy),topiy(ix,iy))-uesy(ix,iy)))       !xpb
      b2divue=0.0_R8
      b2sihs_divue=0.0_R8
      if (balance_netcdf.ne.0) then !djm Jan2017
        shetmp=she0
      endif
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &     bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &     region(ix,iy,0).ne.0) then
          t0 = -phm0*(1.0_R8-BoRiS)*(
     1     uesx(rightix(ix,iy),rightiy(ix,iy))-uesx(ix,iy))
          she0(ix,iy,0) = she0(ix,iy,0)+
     &     (max(0.0_R8,t0)+rf0*abs(t0))*ne(ix,iy)*te(ix,iy)
          she0(ix,iy,3) = she0(ix,iy,3)+min(0.0_R8,t0)-rf0*abs(t0)
          b2divue(region(ix,iy,0))=b2divue(region(ix,iy,0))+
     1     t0*ne(ix,iy)*te(ix,iy)
          b2sihs_divue(ix,iy)=t0*ne(ix,iy)*te(ix,iy)
        endif
       enddo
      enddo
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_divue0to3=she0-shetmp
      endif
*   ..contribution from divergence(ua,va)
      b2divua=0.0_R8
      b2sihs_divua=0.0_R8
      if (balance_netcdf.ne.0) then !djm Jan2017
        shitmp=shi0
      endif
      do is = 0, ns-1
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &       bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &       region(ix,iy,0).ne.0) then
              t0 = -phm1*(na(ix,iy,is)/ni(ix,iy,0))*(1.0_R8-BoRiS)*
     1         (uasx(rightix(ix,iy),rightiy(ix,iy),is)-uasx(ix,iy,is))
              shi0(ix,iy,0) = shi0(ix,iy,0)+
     &         (max(0.0_R8,t0)+rf1*abs(t0))*ni(ix,iy,0)*ti(ix,iy)
              shi0(ix,iy,3) = shi0(ix,iy,3)+min(0.0_R8,t0)-rf1*abs(t0)
              b2divua(region(ix,iy,0))=b2divua(region(ix,iy,0))
     1         +t0*ni(ix,iy,0)*ti(ix,iy)
              b2sihs_divua(ix,iy)=b2sihs_divua(ix,iy)+
     1         t0*ni(ix,iy,0)*ti(ix,iy)
            endif
          enddo
        enddo
      enddo
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_divua0to3=shi0-shitmp
      endif
*   ..contribution from divergence of the electrical drift       !sv 02.06.99
      call intfacev(nx,ny,bottomix,bottomiy,vol,po,wpov)
      call intfaceh(nx,ny,leftix,leftiy,vol,po,wpoh)
      b2exbe=0.0_R8
      b2sihs_exbe=0.0_R8
      if (balance_netcdf.ne.0) then !djm Jan2017
        shetmp=she0
      endif
      do iy = -1, ny
        do ix = -1, nx
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &     bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &     region(ix,iy,0).ne.0) then
            work(ix,iy) = bb(ix,iy,3)*vol(ix,iy)/hx(ix,iy)/hy1(ix,iy)
     &       *((wpov(topix(ix,iy),topiy(ix,iy))-wpov(ix,iy))*
     &       (1.0_R8/wbbr(ix,iy,3)**2-1.0_R8/wbbl(ix,iy,3)**2)
     &       -(wpoh(rightix(ix,iy),rightiy(ix,iy))-wpoh(ix,iy))*
     &       (1.0_R8/wbbv(topix(ix,iy),topiy(ix,iy),3)**2-
     &        1.0_R8/wbbv(ix,iy,3)**2))
            t0 = phm5*fac_ExB(ix,iy)*work(ix,iy)                  	 !xpb
            she0(ix,iy,0) = she0(ix,iy,0)+
     &       (max(0.0_R8,t0)+rf0*abs(t0))*ne(ix,iy)*te(ix,iy)
            she0(ix,iy,3) = she0(ix,iy,3)+min(0.0_R8,t0)-rf0*abs(t0)
            b2exbe(region(ix,iy,0))=b2exbe(region(ix,iy,0))+
     1       t0*ne(ix,iy)*te(ix,iy)
            b2sihs_exbe(ix,iy)=t0*ne(ix,iy)*te(ix,iy)
          endif
        enddo
      enddo
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_exbe0to3=she0-shetmp
      endif
*   ..contribution from divergence of the electrical drift    !sv 02.06.99
      b2exba=0.0_R8
      b2sihs_exba=0.0_R8
      if (balance_netcdf.ne.0) then !djm Jan2017
        shitmp=shi0
      endif
      do is = 0, ns-1
        if(.not.is_neutral(is)) then
          do iy = -1, ny
            do ix = -1, nx
              if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
                t0 = phm6*fac_ExB(ix,iy)*
     &           (na(ix,iy,is)/ni(ix,iy,0))*work(ix,iy)  !xpb
                shi0(ix,iy,0) = shi0(ix,iy,0)+
     &           (max(0.0_R8,t0)+rf1*abs(t0))*ni(ix,iy,0)*ti(ix,iy)
                shi0(ix,iy,3) = shi0(ix,iy,3)+min(0.0_R8,t0)-rf1*abs(t0)
                b2exba(region(ix,iy,0))=b2exba(region(ix,iy,0))+
     1           t0*ni(ix,iy,0)*ti(ix,iy)
                b2sihs_exba(ix,iy)=b2sihs_exba(ix,iy)+
     1           t0*ni(ix,iy,0)*ti(ix,iy)
              endif
            enddo
          enddo
        endif
      enddo
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_exba0to3=shi0-shitmp
      endif
*   ..contribution from viscous heating
      b2visa=0.0_R8
      b2sihs_visa=0.0_R8
      if (balance_netcdf.ne.0) then !djm Jan2017
        shitmp=shi0
      endif
      do is = 0, ns-1
       if(.not.is_neutral(is).or.use_eirene.eq.0) then
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &       bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &       region(ix,iy,0).ne.0) then
           t0 = phm2*vol(ix,iy)*(1.0_R8-BoRiS)*     ! 4.0e0_R8/3.0e0_R8* cvsa already has 4/3 !srv 16.10.17
     &      (cvsa(ix,iy,0,is)*(ua(ix,iy,is)-
     &          ua(leftix(ix,iy),leftiy(ix,iy),is))**2/
     &       (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+
     &       cvsa(rightix(ix,iy),rightiy(ix,iy),0,is)*
     &  (ua(rightix(ix,iy),rightiy(ix,iy),is)-ua(ix,iy,is))**2/
     &       (vol(ix,iy)+vol(rightix(ix,iy),rightiy(ix,iy))))
           shi0(ix,iy,0) = shi0(ix,iy,0)+(1.0_R8+rf2)*t0
           shi0(ix,iy,3) = shi0(ix,iy,3)-rf2*t0/(ni(ix,iy,0)*ti(ix,iy))
           b2visa(region(ix,iy,0))=b2visa(region(ix,iy,0))+t0
           b2sihs_visa(ix,iy)=b2sihs_visa(ix,iy)+t0
*!!      (viscous heating from d/dy is temporarily ignored)
          endif
         enddo
        enddo
       endif
      enddo
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_visa0to3=shi0-shitmp
      endif

*   ..contribution from electron-ion and ion-ion friction in case of arbitrary impurity amount
*     in case of small impurity amount see below for b2sigp_style.ne.2
      if(b2sigp_style.eq.2) then                                         !srv 05.07.17
       call b2xpne (nx, ny, ns, rz2, na, ne2_ext, ne2)
       call b2xzef (nx, ny, ns, rz2, na, ne, zeff)                       !srv 13.01.17 {
       call b2xzmf (nx, ny, ns, ismain, rz2, na, zmff)
       do iy = -1, ny
         do ix = -1, nx
           ce1(ix,iy)=fce1(zeff(ix,iy))
           ce2(ix,iy)=fce2(zeff(ix,iy))
           cimp1(ix,iy)=fce1(zmff(ix,iy))
           cimp2(ix,iy)=fce2(zmff(ix,iy))
         enddo
       enddo                                                             !srv 13.01.17 }
*   ..compute the Coulomb logarithm
       call b2tlnl (nx, ny, icase, te, ti, ne, lnlam)
*   ..compute zetap and zetae
       t1 = (4.0_R8*pi*eps0/qe)/qe
       t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
       czetap = t0*sqrt(mp)
       czetae = t0*sqrt(me)
       zetap = czetap*(ti/ev)*sqrt(ti/ev)/lnlam
       zetae = czetae*(te/ev)*sqrt(te/ev)/lnlam
*   ..compute gti,gte
       do iy = -1, ny
         do ix = -1, nx
           if(leftix(ix,iy).eq.-2) then
             gti(ix,iy) = 0.0_R8
             gte(ix,iy) = 0.0_R8
           else
             gti(ix,iy) = (ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))/
     &        (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
             gte(ix,iy) = (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))/
     &        (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
           endif
         enddo
       enddo
       shifr = 0.0_R8
       shefr = 0.0_R8
       do is = 0,ns-1
         if(.not. is_neutral(is)) then
           call b2sifrtf (nx, ny, ns, is, ismain, bb, na, ua,
     &      ue, te, ti, ni, ne, ne2, zetap, zetae, gti, gte, ce1, ce2,   !srv 05.07.17
     &      cimp1, cimp2, zeff, smbfrea, smbfria, smbtfea, smbtfia,      !srv 05.07.17
     &      smbch, smbtf)
           work =
     &      (smbfrea(:,:,0)+
     &       smbfrea(:,:,1)*ua(:,:,is)+
     &       smbfrea(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smbfrea(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is))+
     &      (smbtfea(:,:,0)+
     &       smbtfea(:,:,1)*ua(:,:,is)+
     &       smbtfea(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smbtfea(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is))
           shefr(:,:,0) = shefr(:,:,0) + work
           shifr(:,:,0) = shifr(:,:,0) + ua(:,:,is)*(work+
     &      (smbfria(:,:,0)+
     &       smbfria(:,:,1)*ua(:,:,is)+
     &       smbfria(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smbfria(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is))+
     &      (smbtfia(:,:,0)+
     &       smbtfia(:,:,1)*ua(:,:,is)+
     &       smbtfia(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smbtfia(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is))
     &                                               )
         endif
       enddo
       shefr(:,:,0) =  ue*shefr(:,:,0)

       b2joule=0.0_R8
       b2sihs_joule=0.0_R8
       b2fraa=0.0_R8
       b2sihs_fraa=0.0_R8
       do iy = 0,ny-1
         do ix =0,nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &     bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &     region(ix,iy,0).ne.0) then
             t0 = phm3*vol(ix,iy)*shefr(ix,iy,0)
!            shefr(ix,iy,0) = t0
             shefr(ix,iy,0) = max(0.0e0_R8,t0)+rf3*abs(t0)
             shefr(ix,iy,3) =
     &        (min(0.0e0_R8,t0)-rf3*abs(t0))/(ne(ix,iy)*te(ix,iy))
             b2joule(region(ix,iy,0))=b2joule(region(ix,iy,0))+t0
             b2sihs_joule(ix,iy)=t0
             t0 = - phm4*vol(ix,iy)*shifr(ix,iy,0)
!            shifr(ix,iy,0) = t0
!             shifr(ix,iy,0) = (2.5e0_R8+rf4)*t0
!             shifr(ix,iy,3) = -(1.5e0_R8+rf4)*t0/(ni(ix,iy,0)*ti(ix,iy))
             shifr(ix,iy,0) = max(0.0e0_R8,t0)+rf4*abs(t0)
             shifr(ix,iy,3) =
     &        (min(0.0e0_R8,t0)-rf4*abs(t0))/(ni(ix,iy,0)*ti(ix,iy))
             b2fraa(region(ix,iy,0))=b2fraa(region(ix,iy,0))+t0
             b2sihs_fraa(ix,iy)=t0
           endif
         enddo
       enddo
       call b2saxpy (n2*4, 1.0e0_R8, shefr, 1, she0, 1)
       call b2saxpy (n2*4, 1.0e0_R8, shifr, 1, shi0, 1)
      else                                                               !srv 05.07.17

*    ..contribution from electron-atom friction
*      (a.k.a. Joule heating)
       call b2xehx (nx, ny, hx, te, po, ne, ehx)                         !srv 13.01.17
       if (istyle_Joule_heating.eq.0)
     &   call b2xehy (nx, ny, te, po, ne, ehy)                           !srv 13.01.17
       b2joule=0.0_R8
       b2sihs_joule=0.0_R8
       if (balance_netcdf.ne.0) then !djm Jan2017
        shetmp=she0
       endif
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &     bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &     region(ix,iy,0).ne.0) then
          t0 = fch(ix,iy,0)*ehx(ix,iy)/
     &     (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+
     &      fch(rightix(ix,iy),rightiy(ix,iy),0)*
     &      ehx(rightix(ix,iy),rightiy(ix,iy))/
     &     (vol(ix,iy)+vol(rightix(ix,iy),rightiy(ix,iy)))
          if (istyle_Joule_heating.eq.0) t0 = t0 +
     &      fch(ix,iy,1)*ehy(ix,iy)/
     &     (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))+
     &      fch(topix(ix,iy),topiy(ix,iy),1)*
     &      ehy(topix(ix,iy),topiy(ix,iy))/
     &     (vol(ix,iy)+vol(topix(ix,iy),topiy(ix,iy)))
          t0 = phm3*vol(ix,iy)*t0
          she0(ix,iy,0) = she0(ix,iy,0)+max(0.0_R8,t0)+rf3*abs(t0)
          she0(ix,iy,3) = she0(ix,iy,3)+
     &     (min(0.0_R8,t0)-rf3*abs(t0))/(ne(ix,iy)*te(ix,iy))
          b2joule(region(ix,iy,0))=b2joule(region(ix,iy,0))+t0
          b2sihs_joule(ix,iy)=t0
         endif
        enddo
       enddo
       if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_joule0to3=she0-shetmp
       endif
*   ..contribution from atom-atom friction
       b2fraa=0.0_R8
       b2sihs_fraa=0.0_R8
       if (balance_netcdf.ne.0) then !djm Jan2017
        shitmp=shi0
       endif
!      do isb = 0, ns-1
       do isb = nsmin, nsmax-1
        do is = isb+1, ns-1
*      ..consider distinct cases
         if (.not.is_neutral(is).and..not.is_neutral(isb)) then
*       ..contribution from friction between charged species
          coef = phm4*(1.0_R8/ctaup)*
     &     mp*sqrt(am(isb)*am(is)/(am(isb)+am(is)))
          do iy = -1, ny
           do ix = -1, nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &        bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &        region(ix,iy,0).ne.0) then
             t0 = rz2(ix,iy,isb)*rz2(ix,iy,is)*
     &        (coef*lnlam(ix,iy)/4.0_R8)*(1.0_R8-BoRiS)*                 !srv 03.05.10
     &         na(ix,iy,is)*na(ix,iy,isb)*vti32(ix,iy)*
     &        ((ua(ix,iy,is)-ua(ix,iy,isb))**2+
     &         (wadia(ix,iy,0,is )+vaecrb(ix,iy,0,is )-
     &          wadia(ix,iy,0,isb)-vaecrb(ix,iy,0,isb))**2+
     &         (wadia(ix,iy,1,is )+vaecrb(ix,iy,1,is )-
     &          wadia(ix,iy,1,isb)-vaecrb(ix,iy,1,isb))**2)
             shi0(ix,iy,0) = shi0(ix,iy,0)+(2.5_R8+rf4)*t0
             shi0(ix,iy,3) = shi0(ix,iy,3)-
     &        (1.5_R8+rf4)*t0/(ni(ix,iy,0)*ti(ix,iy))
             b2fraa(region(ix,iy,0))=b2fraa(region(ix,iy,0))+t0
             b2sihs_fraa(ix,iy)=b2sihs_fraa(ix,iy)+t0
            endif
           enddo
          enddo
         else
*       ..contribution from friction involving neutral species
*         (We ignore such a term, except for charge exchange friction
*         between H0 and H1. That term is calculated via the atomic
*         physics routines.)
         endif
        enddo
       enddo
       if (balance_netcdf.ne.0) then !djm Jan2017
         b2sihs_fraa0to3=shi0-shitmp
       endif
      endif                                                              !srv 05.07.17

      b2str=0.0_R8
      b2sihs_str=0.0_R8
      if (style.eq.-1) then
*    ..strange term, version -1
       if(ncall.eq.0) then
        write(*,*) 'Strange term active with sign -1'
       endif
       do is = 0, ns-1
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &       bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &       region(ix,iy,0).ne.0) then
           t0 = uasy(ix,iy,is)*
     &      (rza(ix,iy,is)*ne(ix,iy)*te(ix,iy)-
     &       rza(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &       ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &       te(bottomix(ix,iy),bottomiy(ix,iy)))*
     &       vol(ix,iy)/
     &      (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))+
     &       uasy(topix(ix,iy),topiy(ix,iy),is)*
     &      (rza(topix(ix,iy),topiy(ix,iy),is)*
     &       ne(topix(ix,iy),topiy(ix,iy))*
     &       te(topix(ix,iy),topiy(ix,iy))-
     &       rza(ix,iy,is)*ne(ix,iy)*te(ix,iy))*
     &       vol(ix,iy)/(vol(ix,iy)+vol(topix(ix,iy),topiy(ix,iy)))
           she0(ix,iy,0) = she0(ix,iy,0)+max(0.0_R8,t0)
           she0(ix,iy,3) = she0(ix,iy,3)+
     &      min(0.0_R8,t0)/(ne(ix,iy)*te(ix,iy))
           shi0(ix,iy,0) = shi0(ix,iy,0)+max(0.0_R8,-t0)
           shi0(ix,iy,3) = shi0(ix,iy,3)+
     &      min(0.0_R8,-t0)/(ni(ix,iy,0)*ti(ix,iy))
           b2str(region(ix,iy,0))=b2str(region(ix,iy,0))+t0
           b2sihs_str(ix,iy)=b2sihs_str(ix,iy)+t0
          endif
         enddo
        enddo
       enddo
       if(ncall.eq.0) then
        write(*,*) 'b2str: ',b2str
       endif
      else if (style.eq.0) then
*    ..strange term, version 0
*      (in this version the term is absent)
       if(ncall.eq.0) then
        write(*,*) 'Strange term not active'
       endif
      else if (style.eq.1) then
*    ..strange term, version 1
       if(ncall.eq.0) then
        write(*,*) 'Strange term active with sign +1'
       endif
       do is = 0, ns-1
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &       bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &       region(ix,iy,0).ne.0) then
           t0 = uasy(ix,iy,is)*
     &      (rza(ix,iy,is)*na(ix,iy,is)*ti(ix,iy)-
     &       rza(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &       na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &       ti(bottomix(ix,iy),bottomiy(ix,iy)))*
     &       vol(ix,iy)/
     &      (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))+
     &       uasy(topix(ix,iy),topiy(ix,iy),is)*
     &      (rza(topix(ix,iy),topiy(ix,iy),is)*
     &       na(topix(ix,iy),topiy(ix,iy),is)*
     &       ti(topix(ix,iy),topiy(ix,iy))-
     &       rza(ix,iy,is)*na(ix,iy,is)*ti(ix,iy))*
     &       vol(ix,iy)/(vol(ix,iy)+vol(topix(ix,iy),topiy(ix,iy)))
           if(ncall.eq.0.and.t0.eq.0.0_R8.and..not.is_neutral(is)) then
            write(*,*) 'str t0=0 ix,iy,is,uasy(ix,iy,is) ',
     1       ix,iy,is,uasy(ix,iy,is)
            write(*,*)
     &       hy1(ix,iy),ua(bottomix(ix,iy),bottomiy(ix,iy),is),
     1       hy1(bottomix(ix,iy),bottomiy(ix,iy)),ua(ix,iy,is),
     2       pbs(ix,iy,1)
            write(*,*) 'T1, T2 ',
     &       uasy(ix,iy,is)*
     &      (rza(ix,iy,is)*na(ix,iy,is)*ti(ix,iy)-
     &       rza(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &       na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &       ti(bottomix(ix,iy),bottomiy(ix,iy)))*
     &       vol(ix,iy)/
     &      (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy)),
     &       uasy(topix(ix,iy),topiy(ix,iy),is)*
     &      (rza(topix(ix,iy),topiy(ix,iy),is)*
     &       na(topix(ix,iy),topiy(ix,iy),is)*
     &       ti(topix(ix,iy),topiy(ix,iy))-
     &       rza(ix,iy,is)*na(ix,iy,is)*ti(ix,iy))*
     &       vol(ix,iy)/(vol(ix,iy)+vol(topix(ix,iy),topiy(ix,iy)))
           endif
           she0(ix,iy,0) = she0(ix,iy,0)+max(0.0_R8,-t0)
           she0(ix,iy,3) = she0(ix,iy,3)+
     &      min(0.0_R8,-t0)/(ne(ix,iy)*te(ix,iy))
           shi0(ix,iy,0) = shi0(ix,iy,0)+max(0.0_R8,t0)
           shi0(ix,iy,3) = shi0(ix,iy,3)+
     &      min(0.0_R8,t0)/(ni(ix,iy,0)*ti(ix,iy))
           b2str(region(ix,iy,0))=b2str(region(ix,iy,0))+t0
           b2sihs_str(ix,iy)=b2sihs_str(ix,iy)+t0
          endif
         enddo
        enddo
       enddo
       if(ncall.eq.0) then
        write(*,*) 'b2str: ',b2str
       endif
      else
        write(*,*) 'b2sihs: style should be -1, 0 or +1 --- not ',style
        call xerrab ('b2sihs--faulty style specified')
      endif
* ..drift
* ..term should not be present: in B2.5 form of final energy equations !xpb
* ..this is already included within the electron energy exchange terms !xpb
* ..i.e. it is counted within the 'phm3' electron-atom friction (Joule !xpb
* ..heating term)                                                      !xpb
* ..!!!!!!!!!!!!!!!!!!! CHECK WITH RUSSIANS !!!!!!!!!!!
      b2wrong1=0.0_R8
      do is = 0, ns-1
        do iy = -1, ny          ! was ny DPC 970415
          do ix = -1, nx        ! was nx DPC 970415
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then

*   ..standard case
              if(periodic_bc.ne.1 .or. rightix(ix,iy).eq.ix+1 .or.
     &           region(ix,iy,0).ne.1 .or. nncut.eq.0) then
                t0=phm7*0.5_R8*bb(ix,iy,3)*(
     1           (fchdia(rightix(ix,iy),rightiy(ix,iy),0)+
     1            fchdia(ix,iy,0))*uadia(ix,iy,1,is)-
     1           (fchdia(topix(ix,iy),topiy(ix,iy),1)+
     2            fchdia(ix,iy,1))*uadia(ix,iy,0,is))
*   ..case where the mesh is discontinuous
              else
                t0=phm7*bb(ix,iy,3)*(fchdia(ix,iy,0)*uadia(ix,iy,1,is)-
     1           0.5_R8*(fchdia(topix(ix,iy),topiy(ix,iy),1)+
     2                   fchdia(ix,iy,1))*uadia(ix,iy,0,is))
              endif
              she0(ix,iy,0) = she0(ix,iy,0)+max(0.0_R8,t0)
              she0(ix,iy,3) = she0(ix,iy,3)+
     &         min(0.0_R8,t0)/(ne(ix,iy)*te(ix,iy))
              shi0(ix,iy,0) = shi0(ix,iy,0)+max(0.0_R8,-t0)
              shi0(ix,iy,3) = shi0(ix,iy,3)+
     &         min(0.0_R8,-t0)/(ni(ix,iy,0)*ti(ix,iy))
              b2wrong1(region(ix,iy,0))=b2wrong1(region(ix,iy,0))+t0
            endif
          enddo
        enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sihs

      end subroutine b2sihs
c
      subroutine b2sdia (nx, ny, ns, na, ne, ti, te, bb, vol, hx,
     1 she, shi, facdrift)                                               !xpb
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_tallies
      use b2mod_b2cmpa
      use b2mod_balance                                                  !djm Jan2017
     & , only : b2sihs_diaa0to3, b2sihs_diae0to3, balance_netcdf
      use b2mod_geo
     , , only : hy1
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny)
      real (kind=R8) ::
     *  bb(-1:nx,-1:ny,0:3), hx(-1:nx,-1:ny)
      real (kind=R8) ::                                                  !xpb
     *  facdrift(-1:nx,-1:ny)                                            !xpb
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   she(-1:nx,-1:ny,0:3), shi(-1:nx,-1:ny,0:3)
*   ..common blocks
*   ..local variables
      integer ncall, n2, ix, iy, is
      integer ixl, iyl, ixr, iyr, ixt, iyt, ixb, iyb
      real (kind=R8) ::
     *  t0, facgt
      save ncall, facgt
      real (kind=R8) :: shetmp(-1:nx,-1:ny,0:3), shitmp(-1:nx,-1:ny,0:3) !djm Jan2017
*   ..procedures
      external subini, subend, xertst, b2xvsg, ipgetr
      intrinsic min, max
*   ..initialisation
      data ncall/0/, facgt/0.0_R8/
 
* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sdia')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test facdrift
      call xertst (0.0_R8.le.minval(facdrift) .and.
     &  1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
*   ..get facgt on first call
      if (ncall.eq.0) then
        call ipgetr ('b2sdia_facgt', facgt)
        call xertst (0.0_R8.le.facgt .and. 1.0_R8.ge.facgt,
     &   'faulty argument facgt')
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol
        call b2xvsg (n2, vol, 1, 'vol', '.gt.')
*    ..test sign of na, te, ti, ne
        call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
        call b2xvsg (n2, te, 1, 'te', '.gt.')
        call b2xvsg (n2, ti, 1, 'ti', '.gt.')
        call b2xvsg (n2, ne, 1, 'ne', '.gt.')
      endif

      if(maxval(facdrift).ne.0.0_R8 .and. facgt.ne.0.0_R8) then
        b2wrong2=0.0_R8
      if (balance_netcdf.ne.0) then !djm Jan2017
        shitmp=shi
      endif
        do is=0,ns-1
          if(.not.is_neutral(is)) then
            do iy=-1,ny
              do ix=-1,nx
                ixl=leftix(ix,iy)
                iyl=leftiy(ix,iy)
                ixr=rightix(ix,iy)
                iyr=rightiy(ix,iy)
                ixt=topix(ix,iy)
                iyt=topiy(ix,iy)
                ixb=bottomix(ix,iy)
                iyb=bottomiy(ix,iy)
                if(ixl.ne.-2.and.ixr.ne.nx+1.and.
     &             iyb.ne.-2.and.iyt.ne.ny+1.and.
     &             region(ix,iy,0).ne.0) then
                  t0=facdrift(ix,iy)*facgt*2.5_R8/qe*
     2             ((vol(ixr,iyr)*ti(ix,iy)**2*na(ix,iy,is)/
     .               rza(ix,iy,is)+
     3               vol(ix,iy)*ti(ixr,iyr)**2*na(ixr,iyr,is)/
     .               rza(ixr,iyr,is))/
     4             (vol(ix,iy)+vol(ixr,iyr))-
     5              (vol(ix,iy)*ti(ixl,iyl)**2*na(ixl,iyl,is)/
     .               rza(ixl,iyl,is)+
     6               vol(ixl,iyl)*ti(ix,iy)**2*na(ix,iy,is)/
     .               rza(ix,iy,is))/
     7             (vol(ixl,iyl)+vol(ix,iy)))/hx(ix,iy)*vol(ix,iy)*
     8             (1.0_R8/bb(ix,iy,3)**2*
     9             ((vol(ixt,iyt)*bb(ix,iy,2)+vol(ix,iy)*bb(ixt,iyt,2))/ !xpb
     1              (vol(ix,iy)+vol(ixt,iyt))-(vol(ix,iy)*bb(ixb,iyb,2)+ !xpb
     2               vol(ixb,iyb)*bb(ix,iy,2))/
     3              (vol(ixb,iyb)+vol(ix,iy)))-                          !xpb
     4             2.0_R8*bb(ix,iy,2)/bb(ix,iy,3)**3*
     5             ((vol(ixt,iyt)*bb(ix,iy,3)+vol(ix,iy)*bb(ixt,iyt,3))/ !xpb
     6              (vol(ix,iy)+vol(ixt,iyt))-(vol(ix,iy)*bb(ixb,iyb,3)+ !xpb
     7               vol(ixb,iyb)*bb(ix,iy,3))/                          !xpb
     8              (vol(ixb,iyb)+vol(ix,iy))))/hy1(ix,iy)
                  t0=t0+facdrift(ix,iy)*facgt*2.5_R8/qe*
     2             ((vol(ixt,iyt)*ti(ix,iy)**2*na(ix,iy,is)/
     .               rza(ix,iy,is)+
     3               vol(ix,iy)*ti(ixt,iyt)**2*na(ixt,iyt,is)/
     .               rza(ixt,iyt,is))/
     4             (vol(ix,iy)+vol(ixt,iyt))-
     5              (vol(ix,iy)*ti(ixb,iyb)**2*na(ixb,iyb,is)/
     .               rza(ixb,iyb,is)+
     6               vol(ixb,iyb)*ti(ix,iy)**2*na(ix,iy,is)/
     .               rza(ix,iy,is))/
     7             (vol(ixb,iyb)+vol(ix,iy)))/hy1(ix,iy)*vol(ix,iy)*
     8             (-1.0_R8/bb(ix,iy,3)**2*(wbbr(ix,iy,2)-wbbl(ix,iy,2)) !xpb
     4              +2.0_R8*bb(ix,iy,2)/bb(ix,iy,3)**3*
     5              (wbbr(ix,iy,3)-wbbl(ix,iy,3)))/hx(ix,iy)             !xpb
                  shi(ix,iy,0) = shi(ix,iy,0)+max(0.0_R8,-t0)            !xpb
                  shi(ix,iy,3) = shi(ix,iy,3)+                           !xpb
     &             min(0.0_R8,-t0)/(na(ix,iy,is)*ti(ix,iy))              !xpb
                  b2wrong2(region(ix,iy,0))=b2wrong2(region(ix,iy,0))+t0
                endif
              enddo
            enddo
          endif
        enddo
        if (balance_netcdf.ne.0) then !djm Jan2017
          b2sihs_diaa0to3=shi-shitmp
        endif
        b2wrong3=0.0_R8
        if (balance_netcdf.ne.0) then !djm Jan2017
          shetmp=she
        endif
        do iy=-1,ny
          do ix=-1,nx
            ixl=leftix(ix,iy)
            iyl=leftiy(ix,iy)
            ixr=rightix(ix,iy)
            iyr=rightiy(ix,iy)
            ixt=topix(ix,iy)
            iyt=topiy(ix,iy)
            ixb=bottomix(ix,iy)
            iyb=bottomiy(ix,iy)
            if(ixl.ne.-2 .and. ixr.ne.nx+1.and.
     &         iyb.ne.-2 .and. iyt.ne.ny+1.and.
     &         region(ix,iy,0).ne.0) then
              t0=facdrift(ix,iy)*facgt*2.5_R8/qe*
     2         ((vol(ixr,iyr)*te(ix,iy)**2*ne(ix,iy)+
     3           vol(ix,iy)*te(ixr,iyr)**2*ne(ixr,iyr))/                 !xpb
     4          (vol(ix,iy)+vol(ixr,iyr))-                               !xpb
     5          (vol(ix,iy)*te(ixl,iyl)**2*ne(ixl,iyl)+                  !xpb
     6           vol(ixl,iyl)*te(ix,iy)**2*ne(ix,iy))/
     7          (vol(ixl,iyl)+vol(ix,iy)))/hx(ix,iy)*vol(ix,iy)*
     8         (1.0_R8/bb(ix,iy,3)**2*
     9         ((vol(ixt,iyt)*bb(ix,iy,2)+vol(ix,iy)*bb(ixt,iyt,2))/     !xpb
     1          (vol(ix,iy)+vol(ixt,iyt))-                               !xpb
     2          (vol(ix,iy)*bb(ixb,iyb,2)+vol(ixb,iyb)*bb(ix,iy,2))/     !xpb
     3          (vol(ixb,iyb)+vol(ix,iy)))-
     4          2.0_R8*bb(ix,iy,2)/bb(ix,iy,3)**3*
     5         ((vol(ixt,iyt)*bb(ix,iy,3)+vol(ix,iy)*bb(ixt,iyt,3))/     !xpb
     6          (vol(ix,iy)+vol(ixt,iyt))-                               !xpb
     7          (vol(ix,iy)*bb(ixb,iyb,3)+vol(ixb,iyb)*bb(ix,iy,3))/     !xpb
     8          (vol(ixb,iyb)+vol(ix,iy))))/hy1(ix,iy)
              t0=t0+facdrift(ix,iy)*facgt*2.5_R8/qe*
     2         ((vol(ixt,iyt)*te(ix,iy)**2*ne(ix,iy)+                    !xpb
     3           vol(ix,iy)*te(ixt,iyt)**2*ne(ixt,iyt))/                 !xpb
     4          (vol(ix,iy)+vol(ixt,iyt))-
     5          (vol(ix,iy)*te(ixb,iyb)**2*ne(ixb,iyb)+                  !xpb
     6           vol(ixb,iyb)*te(ix,iy)**2*ne(ix,iy))/
     7          (vol(ixb,iyb)+vol(ix,iy)))/hy1(ix,iy)*vol(ix,iy)*
     8         (-1.0_R8/bb(ix,iy,3)**2*
     9         (wbbr(ix,iy,2)-wbbl(ix,iy,2))+                            !xpb
     4          2.0_R8*bb(ix,iy,2)/bb(ix,iy,3)**3*
     5         (wbbr(ix,iy,3)-wbbl(ix,iy,3)))/hx(ix,iy)                  !xpb
              she(ix,iy,0) = she(ix,iy,0)+max(0.0_R8,t0)                 !xpb
              she(ix,iy,3) = she(ix,iy,3)+                               !xpb
     &         min(0.0_R8,t0)/(ne(ix,iy)*te(ix,iy))                      !xpb
              b2wrong3(region(ix,iy,0))=b2wrong3(region(ix,iy,0))+t0
            endif
          enddo
        enddo
        if (balance_netcdf.ne.0) then !djm Jan2017
          b2sihs_diae0to3=she-shetmp
        endif
      endif

* ..return
      ncall = ncall+1
      call subend ()
      end subroutine b2sdia

