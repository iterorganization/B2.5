      subroutine b2tlc0 (nCv, nFc, nVx, ns, switch, geo, mpg, na, te, 
     &  ti, tn, rza, cdpa0, cdpa, cdpahz, flc)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : b2tlc0_cutlo
#ifndef NAGFOR
      use b2mod_ad
     & , only : my_out_folder
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  na(nCv,0:ns-1), te(nCv), ti(nCv), tn(nCv), rza(nCv,0:ns-1),
     *  cdpa0(nFc,0:ns-1)
*   ..input and output arguments
      real (kind=R8) ::
     *  cdpa(nFc,0:1,0:ns-1), cdpahz(nFc,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  flc(nFc,0:1,0:ns-1)
*   ..common blocks
*     ------------------------------------------------------------------
*     B2TLC0 computes the neutral diffusion flux limit.              !srv 09.07.03
*     ------------------------------------------------------------------
*   ..local variables
      integer iFc, is
      real (kind=R8) ::
     *  alpha, gamma, t0, t1, t2, t3, dpb, vbar, maxfluxlimit(0:1),
     *  pb(nCv), dpbf(nFc,0:1), pbv(nVx)
      character charns*3
*   ..procedures
      intrinsic max, abs, sqrt
      external xertst
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2tlc0')
*   ..test nCv, nFc, ns
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')

      call sfill(nFc*2*ns, 1.0_R8, flc, 1)
      call sfill(2, 1.0_R8, maxfluxlimit, 1)

      alpha = switch%b2tlc0_alpha
      gamma = switch%b2tlc0_gamma

*   ..nothing done if alpha.eq.0
      if (alpha.ne.0.0_R8) then
       if (switch%b2tlc0_style.eq.0) then
        do is = 0, ns-1
         if (is_neutral(is)) then
*    ..compute partial pressure
         call b2xppb (nCv, rza(1,is), na(1,is), te, ti, tn, is, pb)
*    ..compute differences of partial pressure
         call diff (nCv, nFc, nVx, 0, geo, mpg, pb, pbv, dpbf)
         do iFc = 1, nFc
           t0 = (tn(mpg%fcCv(iFc,1)) + tn(mpg%fcCv(iFc,2)))/2.0_R8
           t1 = (na(mpg%fcCv(iFc,1),is) + na(mpg%fcCv(iFc,2),is))/2.0_R8
*    ..poloidal
            dpb = abs(cdpa(iFc,0,is)*dpbf(iFc,0))
            vbar = sqrt(8.0_R8 * t0 / (pi * (am(is)*mp)))
     &            *geo%fcS(iFc)*abs(geo%fcQalf(iFc,0))
            t2 = 1.0_R8/(1.0_R8 +
     &        (dpb/(alpha*vbar/4.0_R8*t1 + b2tlc0_cutlo))**gamma)
     &                                                 **(1.0_R8/gamma)
            cdpa(iFc,0,is) = cdpa(iFc,0,is)*t2
            cdpahz(iFc,0,is) = cdpahz(iFc,0,is)*t2
            flc(iFc,0,is) = t2
            maxfluxlimit(0)=max(maxfluxlimit(0),t2)
*    ..perpendicular
            dpb = abs(cdpa(iFc,1,is)*dpbf(iFc,1))
            vbar = sqrt(8.0_R8 * t0 / (pi * (am(is)*mp)))
     &        *geo%fcS(iFc)*abs(geo%fcQalf(iFc,1))
            t2 = 1.0_R8/(1.0_R8 +
     &        (dpb/(alpha*vbar/4.0_R8*t1 + b2tlc0_cutlo))**gamma)
     &                                                 **(1.0_R8/gamma)
            cdpa(iFc,1,is)   = cdpa(iFc,1,is)*t2
            cdpahz(iFc,1,is) = cdpahz(iFc,1,is)*t2
            flc(iFc,1,is) = t2
            maxfluxlimit(1)=max(maxfluxlimit(1),t2)
          enddo
         endif
        enddo
       else
cwdk    flux limit total gradient-driven flow, not individual poloidal
cwdk    or radial contributions, since neutral transport is isotropic
        do is = 0, ns-1
         if (is_neutral(is)) then
*    ..compute partial pressure
          call b2xppb (nCv, rza(1,is), na(1,is), te, ti, tn, is, pb)
*    ..compute differences of partial pressure
          call diff (nCv, nFc, nVx, 0, geo, mpg, pb, pbv, dpbf)
          do iFc = 1, nFc
            t0 = (tn(mpg%fcCv(iFc,1)) + tn(mpg%fcCv(iFc,2)))/2.0_R8
            t3 = dpbf(iFc,0)*geo%fcQalf(iFc,0)+                         ! lkw 30.09.2022{
     &           dpbf(iFc,1)*geo%fcQalf(iFc,1)
            t1 = na(mpg%fcCv(iFc,1),is)
            if (t3 .gt. 0) t1 = na(mpg%fcCv(iFc,2),is)                  ! lkw 30.09.2022}
*    ..limit total gradient
            dpb = cdpa0(iFc,is)*sqrt(dpbf(iFc,0)**2 + dpbf(iFc,1)**2)
            vbar = sqrt(8.0_R8 * t0 / (pi * (am(is)*mp)))
     &           *geo%fcS(iFc)
            t2 = 1.0_R8/(1.0_R8 + 
     &        (dpb/(alpha*vbar/4.0_R8*t1 + b2tlc0_cutlo))**gamma)
     &                                                  **(1.0_R8/gamma)
*    ..poloidal
            cdpa(iFc,0,is)   = cdpa(iFc,0,is)*t2
            cdpahz(iFc,0,is) = cdpahz(iFc,0,is)*t2
            flc(iFc,0,is) = t2
            maxfluxlimit(0)=max(maxfluxlimit(0),t2)
*    ..radial
            cdpa(iFc,1,is)   = cdpa(iFc,1,is)*t2
            cdpahz(iFc,1,is) = cdpahz(iFc,1,is)*t2
            flc(iFc,1,is) = t2
            maxfluxlimit(1)=max(maxfluxlimit(1),t2)
          enddo
         endif
        enddo
       endif
      endif

!      write(*,'(1x,a,1p,2e10.2)') 'b2tlc0: maxfluxlimit = ',maxfluxlimit

       if (switch%b2tlc0_iout.ne.0) then
         do is=0,ns-1          
           if (is_neutral(is)) then                            
             write (charns,'(i3.3)') is
             call my_out_us(70,nFc,1, flc(1,0,is),
     &               'fl_fn_th'//charns)
             call my_out_us(70,nFc,1, flc(1,1,is),
     &               'fl_fn_r'//charns)
           endif
         enddo                                            
       endif

*   ..return
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2tlc0
