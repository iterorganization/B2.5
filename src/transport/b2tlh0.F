      subroutine b2tlh0 (nCv, nFc, nVx, ns, geo, mpg, na, tn,
     & chcb, fllim0fhi)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  na(nCv,0:ns-1), tn(nCv),
     1  chcb(nFc,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &	fllim0fhi(nFc,0:1,0:ns-1)
*   ..common blocks
*     ------------------------------------------------------------------
*     B2TLH0 computes the conductive neutral heat transport flux limit.
*     ------------------------------------------------------------------
*   ..local variables
      integer iFc, is, ncall, c3
      integer, save :: iout = 0                                   !srv 13.01.04 {
      real (kind=R8), allocatable :: w1(:)
      character charns*3                                          !srv 13.01.04 }
      real (kind=R8) ::
     *  alpha, gamma, t0, t1, t2, t3, dtn, vthi, maxfluxlimit(0:1),
     *  flux_limit_min_tn, tnv(nVx), dtnf(nFc,0:1)
      real (kind=R8), save :: cutlo = 1.0e-60_R8
      save ncall, alpha, gamma, flux_limit_min_tn
*   ..initialisation
      data ncall/0/, alpha/0.0_R8/, gamma/2.0_R8/,
     *     flux_limit_min_tn/0.0_R8/
*   ..procedures
      intrinsic max, abs, sqrt
      external subini, subend, xertst, ipgetr, ipgeti
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2tlh0')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv .and. 0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read input on first call
      if (ncall.eq.0) then
        call ipgetr ('b2tlh0_alpha', alpha)
        call xertst (0.0_R8.le.alpha, 'invalid input parameter alpha')
        call ipgetr ('b2tlh0_gamma', gamma)
        call xertst (0.0_R8.lt.gamma, 'invalid input parameter gamma')
        call ipgetr ('b2tlh0_flux_limit_min_ti', flux_limit_min_tn) !! do not change name of switch for compatibility with existing input files?
        flux_limit_min_tn=flux_limit_min_tn*ev
        call ipgeti ('b2tlh0_iout', iout)                         !srv 13.01.04
      endif

      call sfill(nFc*2*ns, 1.0_R8, fllim0fhi, 1)
      maxfluxlimit=0.0_R8

      if (alpha.ne.0.0_R8) then
*   .. compute differences of temperature
       call diff(nCv, nFc, nVx, 0, geo, mpg, tn, tnv, dtnf)
       do is = 0, ns-1
        if (is_neutral(is)) then
         do iFc = 1, nFc
          t0 = (tn(mpg%fcCv(iFc,1))    + tn(mpg%fcCv(iFc,2)))/2.0_R8
          t1 = (na(mpg%fcCv(iFc,1),is) + na(mpg%fcCv(iFc,2),is))/2.0_R8
*   .. compute neutral heat conductivity flux limit in poloidal direction
          dtn = abs(chcb(iFc,0,is)*dtnf(iFc,0))
          vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))*
     &      geo%fcS(iFc)*abs(geo%fcQalf(iFc,0))
          t2 = 1.0_R8/(1.0_R8 + 
     .      (dtn/(alpha*vthi*t1*t0 + cutlo))**gamma)**(1.0_R8/gamma)
          chcb(iFc,0,is) = chcb(iFc,0,is)*t2
	  fllim0fhi(iFc,0,is) = t2
          maxfluxlimit(0)=max(maxfluxlimit(0),t2)
*   .. compute neutral heat conductivity flux limit in radial diretion
          dtn = abs(chcb(iFc,1,is)*dtnf(iFc,1))
          vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))*
     &      geo%fcS(iFc)*abs(geo%fcQalf(iFc,1))
          t2 = 1.0_R8/(1.0_R8 +
     .      (dtn/(alpha*vthi*t1*t0 + cutlo))**gamma)**(1.0_R8/gamma)
          chcb(iFc,1,is) = chcb(iFc,1,is)*t2
	  fllim0fhi(iFc,1,is) = t2
          maxfluxlimit(1)=max(maxfluxlimit(1),t2)
         end do
        end if
       end do
      endif

 100  continue

       if (ncall.lt.100) write(*,'(1x,a,1p,2e10.2)')
     .  'b2tlh0: maxfluxlimit = ',maxfluxlimit

       if (iout.ne.0) then
        allocate(w1(nCv))
        do is=0,ns-1                                            !srv 13.01.04 {
         write (charns,'(i3.3)') is
         call my_out_us(70,nFc,1, fllim0fhi(1,0,is),
     &               'fllim0fhix'//charns)
         call my_out_us(70,nFc,1, fllim0fhi(1,1,is),
     &               'fllim0fhiy'//charns)
         call my_out_us(70,nFc,1, chcb(1,0,is), 'chcix'//charns)
         call my_out_us(70,nFc,1, chcb(1,1,is), 'chciy'//charns)
        enddo                                                   !srv 13.01.04 }
        deallocate(w1)
       endif

*   ..return
      ncall = ncall + 1
      call subend ()
      return
*     ------------------------------------------------------------------
      end
