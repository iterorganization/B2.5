      subroutine b2tvspr(nCv, nFc, nVx, ns, switch, geo, mpg,
     &   po, ti, na, rza, vsa0, fac_vis, fchvisper)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     1   fac_vis(nFc), po(nCv), ti(nCv),
     2   na(nCv,0:ns-1), rza(nCv,0:ns-1), vsa0(nCv,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     3   fchvisper(nFc,0:1)
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is
      integer, save :: ncall = 0, iout = 0                              !srv 06.12.17
      real (kind=R8) :: fac_vism
      real (kind=R8) :: va(nCv), dva(nCv), gpo(nFc), gpa(nFc),
     &  wrkc(nCv), wrkf(nFc), wrkv(nVx)
      real (kind=R8) :: weight(nFc,2)
*   ..procedures
      intrinsic max, abs
      external subini, subend, xertst
      external b2xvsg, b2xvff, b2xvfx
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tvspr')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..calculate the coefficients in difference approximation of viscosity perpendicular current
!      call b2txvspr(nx, ny, ns, fac_vis, vol, hx, hy,                   !srv 10.10.17
!     &              qz, bb, avisper, bvisper, cvisper)

*   ..calculate the contribution of viscosity perpendicular current       !srv 22.11.99
      fac_vism = maxval(fac_vis)
      if(fac_vism.ne.0.0_R8 .and. switch%fhe_vis_per.ne.0.0e0_R8.and.
     &   switch%no_current.eq.0) then
       do is = 0,ns-1                                                   !srv 16.10.17
        if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then             !srv 16.10.17

         weight = 1.0_R8

*      ..compute V_{perp,a} on cell faces
         call grad_r(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, gpo)
         call grad_r(nCv, nFc, nVx, 0, geo, mpg, na(:,is)*ti, wrkv, gpa)
         call intface(nCv,nFc,mpg%fcCv,weight,na(:,is)*rza(:,is),wrkf)
         wrkf = -(gpo + gpa/(qe*wrkf))/geo%fcBb(:,3)

*      ..interpolate V_{perp,a} to cell centers
         call intcell(nFc, nCv, mpg, mpg%intcellR, wrkf, va)

*      ..compute radial gradient of V_{perp,a}, in cell centers
         call gradc_div_r(nCv, nFc, nVx, 1, geo, mpg, va, wrkf, dva)
         dva = 0.25_R8*dva*vsa0(:,is)

*      ..set dva in guard cells (constant extrapolation from domain cell)
         call extend(dva,nCv,mpg,1)

*      ..compute current density in cell centers through divergence theorem
         call intface(nCv,nFc,mpg%fcCv,weight,dva,wrkf) 
         wrkf = wrkf*geo%fcS*geo%fcQalf(:,1)
         call div0(nCv, nFc, mpg, wrkf, wrkc)
         wrkc = wrkc/geo%cvVol

*      ..set zero in guard cells
         wrkc(mpg%nCi+1:mpg%nCv) = 0.0_R8

*      ..interpolate to face and multiply by radial area
         call intface(nCv, nFc, mpg%fcCv, weight, wrkc, fchvisper(:,1))
         fchvisper(:,1) = -fchvisper(:,1)*
     &     geo%fcS*geo%fcQalf(:,1)/geo%fcBb(:,3)

        endif                                                           !srv 16.10.17
       enddo                                                            !srv 16.10.17
      else
        fchvisper = 0.0e0_R8
      endif                                                             !srv 22.11.99

      if(switch%b2tvspr_iout.eq.1) then
       call my_out_us(70,nFc,1,fchvisper(1,0),'b2tvspr_fchvisperx')
       call my_out_us(70,nFc,1,fchvisper(1,1),'b2tvspr_fchvispery')
      endif
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvspr

      end subroutine b2tvspr                                            !srv 10.10.17

