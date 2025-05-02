      subroutine b2tlh0 (nCv, nFc, nVx, ns, switch, geo, mpg, na, ti,
     & tn, chcib, chcb, fllim0fhi)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2tlh0, b2tlh0_cutlo
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  na(nCv,0:ns-1), ti(nCv), tn(nCv), 
     1  chcib(nFc,0:ns-1), chcb(nFc,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fllim0fhi(nFc,0:1,0:ns-1)
*   ..common blocks
*     ------------------------------------------------------------------
*     B2TLH0 computes the conductive neutral heat transport flux limit.
*     ------------------------------------------------------------------
*   ..local variables
      integer iFc, is
      real (kind=R8), allocatable :: w1(:)                        !srv 13.01.04 {
      character charns*3                                          !srv 13.01.04 }
      real (kind=R8) ::
     *  alpha, gamma, t0, t1, t2, dtn, vthi, maxfluxlimit(0:1),
     *  tnv(nVx), dtnf(nFc,0:1)
*   ..initialisation
*   ..procedures
      intrinsic max, abs, sqrt
      external xertst, ipgetr, ipgeti
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2tlh0')
*   ..test nCv, nFc, ns
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read input on first call
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c        call ipgetr ('b2tlh0_flux_limit_min_ti', flux_limit_min_tn) !! do not change name of switch for compatibility with existing input files?
c        flux_limit_min_tn=flux_limit_min_tn*ev

      alpha = switch%b2tlh0_alpha
      gamma = switch%b2tlh0_gamma
      call sfill(nFc*2*ns, 1.0_R8, fllim0fhi, 1)
      maxfluxlimit=0.0_R8

      if (alpha.ne.0.0_R8) then
*   .. compute differences of temperature
       call diff(nCv, nFc, nVx, 0, geo, mpg, tn, tnv, dtnf)
       if (switch%b2tlh0_style.eq.0) then
#if (defined(_OPENMP) && !defined(NO_OPENMP_B2TLH0))
!$omp parallel do default(none)
!$omp& private(is, iFc, t0, t1, dtn, vthi, t2)
!$omp& shared(ns, nFc, geo, mpg, tn, dtnf, am, chcb, gamma, alpha)
!$omp& shared(fllim0fhi, is_neutral, b2tlh0_cutlo, na)
!$omp& reduction(max: maxfluxlimit)
#endif 
        do is = 0, ns-1
         if (is_neutral(is)) then
          do iFc = 1, nFc
           t0 = (tn(mpg%fcCv(iFc,1))    + tn(mpg%fcCv(iFc,2)))/2.0_R8
           t1 = (na(mpg%fcCv(iFc,1),is) + na(mpg%fcCv(iFc,2),is))/2.0_R8
*   .. compute neutral heat conductivity flux limit in poloidal direction
           dtn = abs(chcb(iFc,0,is)*dtnf(iFc,0))
           vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))*
     &       geo%fcS(iFc)*abs(geo%fcQalf(iFc,0))
           t2 = 1.0_R8/(1.0_R8 +
     .       (dtn/(alpha*vthi*t1*t0 + b2tlh0_cutlo))**gamma)**
     .        (1.0_R8/gamma)
           chcb(iFc,0,is) = chcb(iFc,0,is)*t2
           fllim0fhi(iFc,0,is) = t2
           maxfluxlimit(0)=max(maxfluxlimit(0),t2)
*   .. compute neutral heat conductivity flux limit in radial diretion
           dtn = abs(chcb(iFc,1,is)*dtnf(iFc,1))
           vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))*
     &      geo%fcS(iFc)*abs(geo%fcQalf(iFc,1))
           t2 = 1.0_R8/(1.0_R8 +
     .       (dtn/(alpha*vthi*t1*t0 + b2tlh0_cutlo))**gamma)**
     .        (1.0_R8/gamma)
           chcb(iFc,1,is) = chcb(iFc,1,is)*t2
           fllim0fhi(iFc,1,is) = t2
           maxfluxlimit(1)=max(maxfluxlimit(1),t2)
          end do
         end if
        end do
#if (defined(_OPENMP) && !defined(NO_OPENMP_B2TLH0))
!$omp end parallel do 
#endif 
       else
cwdk    flux limit total gradient-driven flow, not individual poloidal
cwdk    or radial contributions, since neutral transport is isotropic
#if (defined(_OPENMP) && !defined(NO_OPENMP_B2TLH0))
!$omp parallel do default(none)
!$omp& private(is, iFc, t0, t1, dtn, vthi, t2)
!$omp& shared(ns, nFc, geo, mpg, tn, dtnf, am, chcb, gamma, alpha)
!$omp& shared(fllim0fhi, is_neutral, b2tlh0_cutlo, chcib, na)
!$omp& reduction(max: maxfluxlimit)
#endif 
        do is = 0, ns-1
         if (is_neutral(is)) then
          do iFc = 1, nFc
           t0 = (tn(mpg%fcCv(iFc,1))    + tn(mpg%fcCv(iFc,2)))/2.0_R8
           t1 = (na(mpg%fcCv(iFc,1),is) + na(mpg%fcCv(iFc,2),is))/2.0_R8
*   .. compute neutral heat conductivity flux limit
           dtn = chcib(iFc,is)*sqrt(dtnf(iFc,0)**2+dtnf(iFc,1)**2)
           vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))*geo%fcS(iFc)
           t2 = 1.0_R8/(1.0_R8 + 
     .       (dtn/(alpha*vthi*t1*t0 + b2tlh0_cutlo))**gamma)**
     .        (1.0_R8/gamma)
*   .. limit in poloidal direction
           chcb(iFc,0,is) = chcb(iFc,0,is)*t2
           fllim0fhi(iFc,0,is) = t2
           maxfluxlimit(0)=max(maxfluxlimit(0),t2)
*   .. limit in radial diretion
           chcb(iFc,1,is) = chcb(iFc,1,is)*t2
           fllim0fhi(iFc,1,is) = t2
           maxfluxlimit(1)=max(maxfluxlimit(1),t2)
          end do
         end if
        end do
#if (defined(_OPENMP) && !defined(NO_OPENMP_B2TLH0))
!$omp end parallel do
#endif 
       endif
      endif

      if (ncall_b2tlh0.lt.100) write(*,'(1x,a,1p,2e10.2)')
     .  'b2tlh0: maxfluxlimit = ',maxfluxlimit

      if (switch%b2tlh0_iout.ne.0) then
        allocate(w1(nCv))
        do is=0,ns-1                                            !srv 13.01.04 {
         write (charns,'(i3.3)') is
         call my_out_us(70,nFc,1, fllim0fhi(1,0,is),
     &               'fllim0fhi_th'//charns)
         call my_out_us(70,nFc,1, fllim0fhi(1,1,is),
     &               'fllim0fhi_r'//charns)
         call my_out_us(70,nFc,1, chcb(1,0,is), 'chci_th'//charns)
         call my_out_us(70,nFc,1, chcb(1,1,is), 'chci_r'//charns)
        enddo                                                   !srv 13.01.04 }
        deallocate(w1)
      endif

*   ..return
      ncall_b2tlh0 = ncall_b2tlh0 + 1
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2tlh0
