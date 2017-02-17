*-----------------------------------------------------------------------
*.specification

      subroutine b2sifr (nx, ny, ns, isb, vol, hx, hz, bb, na, ua,   !srv 09.01.01
     &   te, ti, ni, ne, ne2, ehxp, lnlam, cthe, cthi, smbch, smbtf, !srv 01.07.09 11.09.09
     &   smprb, smptb, smfrb, vti32, gti, gte)                       !srv 02.07.08
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_external
      use b2mod_b2cmpa
      use b2mod_balance                                              !djm Jan2017
     & , only : b2sifr_smoch0to3, b2sifr_smotf0to3, balance_netcdf
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, isb
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hz(-1:nx,-1:ny),          !srv 09.01.01
     &  bb(-1:nx,-1:ny,0:3),
     &  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  ti(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), ne(-1:nx,-1:ny),
     &  ne2(-1:nx,-1:ny), ehxp(-1:nx,-1:ny), lnlam(-1:nx,-1:ny),     !srv 01.07.09
     &  cthe(-1:nx,-1:ny,0:ns-1), cthi(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   smbch(-1:nx,-1:ny,0:3), smbtf(-1:nx,-1:ny,0:3),             !srv 11.09.09
     &   smprb(-1:nx,-1:ny), smptb(-1:nx,-1:ny), smfrb(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   vti32(-1:nx,-1:ny), gti(-1:nx,-1:ny), gte(-1:nx,-1:ny)      !srv 02.07.08
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIFR computes a linearised representation of the friction
*     force between one species and all others.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes a linear representation of the force
*     exerted on one species by all others. This includes the ordinary
*     friction force proportional to the velocity difference, as well
*     as electron-atom and atom-atom thermal forces. It is an important
*     property that these forces add up to 0 when summed over all
*     species. This property must be maintained numerically.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, styl0, styl1, n2, ix, iy, is, k
      integer, save :: iout = 0                                      !srv 17.06.02
      integer, save :: icase = 0
      character chns*3, chk*1                                        !srv 17.06.02
      real (kind=R8) ::
     *   phm0, phm1, phm2, phm3, ctaup, coef, t0, t1
      real (kind=R8) ::
     1 delti, delte, zlamii, zlamee, dtiapr, dteapr, limthii, limthee
      real (kind=R8) :: wrk0(-1:nx,-1:ny)
*   ..procedures
      intrinsic sqrt, log, sign, max, abs
      real (kind=R8) ::
     *   roxa
      external subini, subend, xertst, ipgeti, ipgetr, sfill
      external b2xvfx, b2xvsg, b2xxgs
      roxa(ix,iy,is) = am(is)*mp*na(ix,iy,is)
*   ..initialisation
      save ncall, styl0, styl1, 
     &  phm0, phm1, phm2, phm3, limthii, limthee
      data ncall/0/, styl0/0/, styl1/0/, 
     &  phm0/1.0_R8/, phm1/1.0_R8/, phm2/1.0_R8/, phm3/1.0_R8/, 
     &  limthii/0.3_R8/, limthee/0.3_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sifr')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2sifr_styl0', styl0)
       call ipgeti ('b2sifr_styl1', styl1)
       call ipgetr ('b2sifr_phm0', phm0)
       call ipgetr ('b2sifr_phm1', phm1)
       call ipgetr ('b2sifr_phm2', phm2)
       call ipgetr ('b2sifr_phm3', phm3)
       call ipgetr ('b2sifr_limthii', limthii)
       call ipgetr ('b2sifr_limthee', limthee)
       call xertst (0.le.styl0.and.styl0.lt.9.and.0.le.styl1.and.
     &  styl1.lt.2, 'faulty internal parameter styl0, styl1')
       call xertst (0.0_R8.le.phm0.and.0.0_R8.le.phm1.and.
     &  0.0_R8.le.phm2.and.0.0_R8.le.phm3,
     &  'faulty internal parameter phm0, phm1, phm2, phm3')
       call xertst (0.0_R8.le.limthii.and.0.0_R8.le.limthee,
     &  'faulty internal parameter limthii, limthee')
       call ipgeti ('b2sifr_iout', iout)                             !srv 17.06.02
       call ipgeti ('b2tlnl_ii', icase)
       call xertst (0.eq.icase .or. 1.eq.icase, 
     &               'faulty input b2tlnl_ii')
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
*    ..test sign of na, ni, ti, ne, ne2
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, ne2, 1, 'ne2', '.gt.')
      endif
*   ..compute ctaup
*     Following Balescu; Braginskii has ctaup different by a factor
*     sqrt(2). The mathematical expression is:
*       t0 = (3/(4*sqrt(2*pi)))*(4*pi*eps0/qe**2)**2*ev*sqrt(ev)
*       ctaup = t0*sqrt(mp)/lnlam
*     The numerical value in SI units is: ctaup = 1.4743e13/lnlam.
*     The following calculation of t0 is designed to avoid overflow
*     even for IEEE single precision arithmetic.
!xpb The lnlam term has been removed and is put back in later so it can
!xpb vary with space and species.
      t1 = (4.0_R8*pi*eps0/qe)/qe
      t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
      ctaup = t0*sqrt(mp)                                          !xpb
*   ..compute the Coulomb logarithm
      call b2tlnl (nx, ny, icase, te, ti, ne, lnlam)               !srv 01.07.09  20.09.11
*   ..compute vti32
      do iy = -1, ny
       do ix = -1, nx
        vti32(ix,iy) = vol(ix,iy)/((ti(ix,iy)/ev)*sqrt(ti(ix,iy)/ev))
     &                 *hz(ix,iy)                                     !srv 09.01.01
       enddo
      enddo
*   ..compute gti,gte
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         gti(ix,iy) = 0.0_R8
         gte(ix,iy) = 0.0_R8
        else
         gti(ix,iy) = ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy))
         gte(ix,iy) = te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))
        endif
       enddo
      enddo
      call b2xxgs (nx, ny, 0.0_R8, gti, 1)
      call b2xxgs (nx, ny, 0.0_R8, gte, 1)
* ..compute friction
*   ..initialise to 0
      call sfill (n2*4, 0.0_R8, smbch, 1)
      call sfill (n2*4, 0.0_R8, smbtf, 1)
      smprb = 0.0_R8
      smptb = 0.0_R8
      smfrb = 0.0_R8
*   ..loop over all other species
      do is = 0, ns-1
*    ..consider distinct cases
       if (is.ne.isb.and.
     &    (.not.is_neutral(is)).and.(.not.is_neutral(isb))) then
*     ..contribution from friction between charged species
        coef = phm0*(1.0_R8/ctaup)*
     &    mp*sqrt(am(isb)*am(is)/(am(isb)+am(is)))
        do iy = -1, ny
         if (styl0.eq.0) then
          do ix = -1, nx
           t0 = na(ix,iy,isb)*na(ix,iy,is)*vti32(ix,iy)*lnlam(ix,iy)*  !srv 01.07.09
     1          rz2(ix,iy,isb)*rz2(ix,iy,is)
           smbch(ix,iy,0) = smbch(ix,iy,0)+coef*t0*ua(ix,iy,is)
           smbch(ix,iy,1) = smbch(ix,iy,1)-coef*t0
           smfrb(ix,iy) = smfrb(ix,iy) + 
     &       coef*t0*(ua(ix,iy,is)-ua(ix,iy,isb))
          enddo
         else if (styl0.eq.1) then
          do ix = -1, nx
           t0 = na(ix,iy,isb)*na(ix,iy,is)*vti32(ix,iy)*lnlam(ix,iy)*  !srv 01.07.09
     1          rz2(ix,iy,isb)*rz2(ix,iy,is)
           smbch(ix,iy,0) = smbch(ix,iy,0)+coef*t0*ua(ix,iy,is)
           smbch(ix,iy,3) = smbch(ix,iy,3)-coef*t0/
     &       roxa(ix,iy,isb)
           smfrb(ix,iy) = smfrb(ix,iy) + 
     &       coef*t0*(ua(ix,iy,is)-ua(ix,iy,isb))
          enddo
         else if (styl0.eq.2) then
          do ix = -1, nx
           t0 = na(ix,iy,isb)*na(ix,iy,is)*vti32(ix,iy)*lnlam(ix,iy)*  !srv 01.07.09
     1          rz2(ix,iy,isb)*rz2(ix,iy,is)
           smbch(ix,iy,2) = smbch(ix,iy,2)+coef*t0*ua(ix,iy,is)/
     &       roxa(ix,iy,isb)
           smbch(ix,iy,3) = smbch(ix,iy,3)-coef*t0/
     &       roxa(ix,iy,isb)
           smfrb(ix,iy) = smfrb(ix,iy) + 
     &       coef*t0*(ua(ix,iy,is)-ua(ix,iy,isb))
          enddo
         endif
        enddo
       else
*     ..contribution from friction involving neutral species
*       (We ignore such a term, except for charge exchange friction
*       between H0 and H1. That term is calculated via the atomic
*       physics routines.)
       endif
      enddo
*   ..compute thermal force
      do iy = -1, ny
        if (styl1.eq.0) then
          do ix = -1, nx
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            t0 = na(ix,iy,isb)*hz(ix,iy)*                            !srv 09.01.01
     &       (rza(ix,iy,isb)/(ne(ix,iy)-ne_ext(ix,iy))-
     &        rz2(ix,iy,isb)/(ne2(ix,iy)-ne2_ext(ix,iy)))
            delte=vol(ix,iy)*
     &    (gte(ix,iy)/(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+
     &     gte(rightix(ix,iy),rightiy(ix,iy))/(vol(ix,iy)+
     &       vol(rightix(ix,iy),rightiy(ix,iy))))
            if(limthee.gt.0.0_R8) then
              zlamee=1.5e16_R8*(te(ix,iy)/ev)**2/ne(ix,iy)
              dteapr=limthee/zlamee*
     1         te(ix,iy)*hx(ix,iy)/(abs(bb(ix,iy,0))/bb(ix,iy,3))
              delte=sign(min(abs(delte),dteapr),delte)
            endif
            delti=vol(ix,iy)*
     &    (gti(ix,iy)/(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+
     &     gti(rightix(ix,iy),rightiy(ix,iy))/(vol(ix,iy)+
     &       vol(rightix(ix,iy),rightiy(ix,iy))))
            if(limthii.gt.0.0_R8) then
              zlamii=1.5e16_R8*(ti(ix,iy)/ev)**2/ne2(ix,iy)
              dtiapr=limthii/zlamii*
     1         ti(ix,iy)*hx(ix,iy)/(abs(bb(ix,iy,0))/bb(ix,iy,3))
              delti=sign(min(abs(delti),dtiapr),delti)
            endif
c           if(cthi(ix,iy,isb).ne.2.65_R8)
c     1      write(*,*) ix,iy,isb,cthi(ix,iy,isb)
            smbtf(ix,iy,0) = smbtf(ix,iy,0)+t0*
     &       (phm1*ehxp(ix,iy)-
     &       (phm2*cthe(ix,iy,isb)*delte*ne(ix,iy)
     &       +phm3*cthi(ix,iy,isb)*delti*(ne(ix,iy)-ne_ext(ix,iy)))*
     1       (vol(ix,iy)/hx(ix,iy))*bb(ix,iy,0)/bb(ix,iy,3))
            smprb(ix,iy) = smprb(ix,iy) + t0*phm1*ehxp(ix,iy)
            smptb(ix,iy) = smptb(ix,iy) - t0*
     &           (phm2*cthe(ix,iy,isb)*delte*ne(ix,iy)
     &           +phm3*cthi(ix,iy,isb)*delti*(ne(ix,iy)-ne_ext(ix,iy)))*
     1           (vol(ix,iy)/hx(ix,iy))*bb(ix,iy,0)/bb(ix,iy,3)
           endif
          enddo
        else if (styl1.eq.1) then
          do ix = -1, nx
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            t0 = vol(ix,iy)*
     &    (gti(ix,iy)/(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+
     &     gti(rightix(ix,iy),rightiy(ix,iy))/(vol(ix,iy)+
     &       vol(rightix(ix,iy),rightiy(ix,iy))))
            t1 = vol(ix,iy)*
     &    (gte(ix,iy)/(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+
     &     gte(rightix(ix,iy),rightiy(ix,iy))/(vol(ix,iy)+
     &       vol(rightix(ix,iy),rightiy(ix,iy))))
            smbtf(ix,iy,0) = smbtf(ix,iy,0)+
     &    na(ix,iy,isb)*hz(ix,iy)*                                  !srv 09.01.01
     &  (rza(ix,iy,isb)/(ne(ix,iy)-ne_ext(ix,iy))-
     &   rz2(ix,iy,isb)/(ne2(ix,iy)-ne2_ext(ix,iy)))*
     &  (phm1*ehxp(ix,iy)-
     &  (phm3*cthi(ix,iy,isb)*t0*(ne(ix,iy)-ne_ext(ix,iy))+
     &   phm2*cthe(ix,iy,isb)*t1*ne(ix,iy))*
     &  (vol(ix,iy)/hx(ix,iy))*bb(ix,iy,0)/bb(ix,iy,3))   !srv 12.04.10 11.01.13
            smprb(ix,iy) = smprb(ix,iy) + 
     &    na(ix,iy,isb)*
     &  (rza(ix,iy,isb)/(ne(ix,iy)-ne_ext(ix,iy))-
     &   rz2(ix,iy,isb)/(ne2(ix,iy)-ne2_ext(ix,iy)))*
     &  phm1*ehxp(ix,iy)
            smptb(ix,iy) = smptb(ix,iy) -
     &    na(ix,iy,isb)*
     &  (rza(ix,iy,isb)/(ne(ix,iy)-ne_ext(ix,iy))-
     &   rz2(ix,iy,isb)/(ne2(ix,iy)-ne2_ext(ix,iy)))*
     &  (phm3*cthi(ix,iy,isb)*t0*(ne(ix,iy)-ne_ext(ix,iy))+
     &   phm2*cthe(ix,iy,isb)*t1*ne(ix,iy))*
     &  (vol(ix,iy)/hx(ix,iy))*bb(ix,iy,0)/bb(ix,iy,3)
           endif
          enddo
        endif
      enddo
      if (iout.ne.0) then                        !srv 17.06.02 {
	write (chns,'(i3.3)') isb
	call my_out(70,nx,ny,ehxp,'b2sifr_ehxp'//chns)
	call my_out(70,nx,ny,cthe(-1,-1,isb),'b2sifr_cthe'//chns)
	call my_out(70,nx,ny,cthi(-1,-1,isb),'b2sifr_cthi'//chns)
	do k=0,3
	  write (chk,'(i1)') k
	  call my_out(70,nx,ny,smbtf(-1,-1,k),
     *                'b2sifr_smotf'//chk)
	  call my_out(70,nx,ny,smbch(-1,-1,k),
     *                'b2sifr_smoch'//chk)
	enddo
        do is=0,ns-1
	  write (chns,'(i3.3)') is
          wrk0=smbtf(:,:,0)+smbtf(:,:,1)*ua(:,:,is)+
     &         smbtf(:,:,2)*na(:,:,is)*mp*am(is)+
     &         smbtf(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is)
	  call my_out(70,nx,ny,wrk0,'b2sifr_smotf'//chns)
          wrk0=smbch(:,:,0)+smbch(:,:,1)*ua(:,:,is)+
     &         smbch(:,:,2)*na(:,:,is)*mp*am(is)+
     &         smbch(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is)
	  call my_out(70,nx,ny,wrk0,'b2sifr_smoch'//chns)
	enddo
      endif                                                            !srv 17.06.02 }

cdjm Jan2017 Store the linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2sifr_smoch0to3(-1:nx,-1:ny,0:3,isb)=smbch
        b2sifr_smotf0to3(-1:nx,-1:ny,0:3,isb)=smbtf
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sifr

      end subroutine b2sifr
