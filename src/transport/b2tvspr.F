      subroutine b2tvspr(nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &   po, ti, na, rza, vsa0, fac_vis, fchvisper_a, fchvisper)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : my_out_folder
#ifdef DBG
     &        , my_out_sf                                               !srv 27.06.21
#endif
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     1   fac_vis(nFc), po(nCv), ti(nCv),
     2   na(nCv,0:ns-1), rza(nCv,0:ns-1), vsa0(nCv,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     3   fchvisper_a(nFc,0:1,0:ns-1), fchvisper(nFc,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TVSPR computes the poloidal and radial components of the
*     perpendicular viscosity current.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is
      integer, save :: ncall = 0                                         !srv 06.12.17
      real (kind=R8) :: fac_vism
      real (kind=R8) :: va(nCv), dva(nCv), gpo(nFc), gpa(nFc),
     &  wrkc(nCv), wrkf(nFc), wrkv(nVx)
      real (kind=R8) :: weight(nFc,2)
*   ..procedures
      external xertst
      external b2xvsg
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tvspr')
*   ..test nCv, nFc, ns
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..calculate the coefficients in difference approximation of viscosity perpendicular current
!      call b2txvspr(nx, ny, ns, fac_vis,                               !srv 10.10.17
!     &              avisper, bvisper, cvisper)

*   ..calculate the contribution of viscosity perpendicular current     !srv 22.11.99
      fac_vism = maxval(fac_vis)
      if(fac_vism.ne.0.0_R8 .and. switch%fhe_vis_per.ne.0.0_R8.and.
     &   switch%no_current.eq.0) then

       fchvisper = 0.0_R8
       do is = 0,ns-1                                                   !srv 16.10.17
        if((switch%visper_all_ions.ne.0 .and. .not.is_neutral(is)) .or.
     &     (switch%visper_all_ions.eq.0 .and. is.eq.ismain)) then       !srv 16.10.17

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
         call extend(dva,mpg,1)

*      ..compute current density in cell centers through divergence theorem
         call intface(nCv,nFc,mpg%fcCv,weight,dva,wrkf) 
         wrkf = wrkf*geo%fcS*geo%fcQalf(:,1)
         call div0(nCv, nFc, mpg, wrkf, wrkc)
         wrkc = wrkc/geo%cvVol

*      ..set zero in guard cells
         wrkc(mpg%nCi+1:mpg%nCv) = 0.0_R8

*      ..interpolate to face and multiply by radial area
         fchvisper_a(:,0,is) = 0.0_R8
         call intface(nCv, nFc, mpg%fcCv, weight, wrkc,
     &     fchvisper_a(:,1,is))
         fchvisper_a(:,1,is) = 
     &    -fchvisper_a(:,1,is)*geo%fcS*geo%fcQalf(:,1)/geo%fcBb(:,3)

        else

         fchvisper_a(:,:,is) = 0.0_R8

        endif                                                           !srv 16.10.17

*     ..accumulate individual contributions into fchvisper
        fchvisper = fchvisper + fchvisper_a(:,:,is)

       enddo                                                            !srv 16.10.17
      else
        fchvisper = 0.0_R8
        fchvisper_a = 0.e0_R8
      endif                                                             !srv 22.11.99

      if(switch%b2tvspr_iout.eq.1) then
       call my_out_us(70,nFc,1,fchvisper(1,0),'b2tvspr_fchvisper_th')
       call my_out_us(70,nFc,1,fchvisper(1,1),'b2tvspr_fchvisper_r')
      endif
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvspr

      end subroutine b2tvspr                                            !srv 10.10.17
