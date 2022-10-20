      subroutine b2tvskt(nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &   po, na, kt, vaecrb, dna_ExB, fchviskt)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : my_out_folder
#ifdef WG_DEBUG
     &        , my_out_sf                                               !srv 27.06.21
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     1   po(nCv), kt(nCv), na(nCv,0:ns-1), dna_ExB(nCv),
     2   vaecrb(nFc,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     3   fchviskt(nFc,0:1)
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer, save :: ncall = 0
      real (kind=R8) :: va(nCv), dva(nCv),
     &  wrk0(nCv), wrkf(nFc), wrkv(nVx), wrkg(nFc,0:1)
      real (kind=R8) :: weight(nFc,2)
*   ..procedures
      intrinsic max, abs
      external xertst
      external b2xvsg, b2xvff, b2xvfx
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tvskt')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')

*   ..calculate the contribution of perpendicular current due to RS in 
*     kt-model
      if(switch%b2tfhe_vis_kt.ne.0.0e0_R8.and.
     &   switch%no_current.eq.0) then

        weight = 1.0_R8

*     ..compute gradients of 'RS pressure' on cell faces
        wrk0 = 2.0_R8/3.0_R8*na(:,ismain)*kt
        call grad(nCv, nFc, nVx, 0, geo, mpg, wrk0, wrkv, wrkg)

*     ..contribution from trace RS tensor
        fchviskt(:,0) = -switch%b2tfhe_vis_kt*geo%fcS*geo%fcQalf(:,0)*
     &    wrkg(:,1)*geo%fcBb(:,2)/geo%fcBb(:,3)**2
        fchviskt(:,1) =  switch%b2tfhe_vis_kt*geo%fcS*geo%fcQalf(:,1)*
     &    wrkg(:,0)*geo%fcBb(:,2)/geo%fcBb(:,3)**2

*     ..compute V_{perp,a} on cell faces
        call grad(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, wrkg)
        wrkf = -wrkg(:,1)/geo%fcBb(:,3)

*     ..interpolate V_{perp,a} to cell centers
        call intcell(nFc, nCv, mpg, mpg%intcellR, wrkf, va)

*     ..compute radial gradient of V_{perp,a}, in cell centers
        call gradc_div_r(nCv, nFc, nVx, 1, geo, mpg, va, wrkf, dva)
        dva = switch%b2sikt_fac_vis_RS*dva*
     &    dna_ExB*am(ismain)*mp*na(:,ismain)

*     ..set dva in guard cells (constant extrapolation from domain cell)
        call extend(dva,nCv,mpg,1)

*     ..compute current density in cell centers through divergence theorem
        call intface(nCv,nFc,mpg%fcCv,weight,dva,wrkf) 
        wrkf = wrkf*geo%fcS*geo%fcQalf(:,1)
        call div0(nCv, nFc, mpg, wrkf, wrk0)
        wrk0 = wrk0/geo%cvVol

*     ..set zero in guard cells
        wrk0(mpg%nCi+1:mpg%nCv) = 0.0_R8

*     ..interpolate to face and multiply by radial area
        call intface(nCv, nFc, mpg%fcCv, weight, wrk0, wrkf)
        fchviskt(:,1) = fchviskt(:,1) + switch%b2tfhe_vis_kt*
     &    wrkf*geo%fcS*geo%fcQalf(:,1)/geo%fcBb(:,3)

      else

        fchviskt = 0.0_R8

      endif

      if(switch%b2tvskt_iout.eq.1) then
       call my_out_us(70,nFc,1,fchviskt(1,0),'b2tvskt_fchviskt_th')
       call my_out_us(70,nFc,1,fchviskt(1,1),'b2tvskt_fchviskt_r')
      endif
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvskt

      end subroutine b2tvskt

