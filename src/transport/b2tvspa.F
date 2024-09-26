*-----------------------------------------------------------------------
*.specification

      subroutine b2tvspa (nCv, nFc, nVx, ns, ismain, switch, geo, mpg,  !srv 20.09.06 { 16.10.17
     &   ua, vsaf_cl, fac_vis, fchvispar_a, fchvispar)                  !srv 22.06.08
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain                                 !lkw 20.06.2022
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &  ua(nCv,0:ns-1), vsaf_cl(nFc,0:1,0:ns-1), fac_vis(nFc)           !srv 22.06.08
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fchvispar_a(nFc,0:1,0:ns-1), fchvispar(nFc,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TVSPA computes the poloidal and radial components of the
*     parallel viscosity current.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      integer is, iFc
      real (kind=R8) ::
     *   dvpar(nFc),
     *   wrk(nCv), wrk0(nFc), wrk1(nFc), wrkVx(nVx), gOnedBsq(nFc,0:1)
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx
*   ..initialization
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tvspa')
*   ..set internal parameters on first call
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')

      dvpar     = 0.0_R8                                                 !srv 24.06.99
      if (switch%no_current.eq.0.and.switch%fhe_vis_par.ne.0.0_R8
     &                          .and.maxval(fac_vis).ne.0.0_R8) then     !srv

*     ..compute gradients of OnedBsq on cell faces
        call grad(nCv, nFc, nVx, 1, geo, mpg, geo%cvOnedBsq,
     &    geo%vxOnedBsq, gOnedBsq)

*     ..compute pre-factor
        ! Note: directly compute magnetic prefactor on faces
        !       (original code: averaging from cell centers)
        wrk0 = geo%fcBb(:,2)*sqrt(geo%fcBb(:,3))/(3.0_R8*geo%fcBb(:,0))

*     ..compute parallel viscosity current
        fchvispar = 0.0_R8
        do is = 0, ns-1                                                  !srv 24.06.99

*       ..compute dvpar
          if (is.eq.ismain .and. .not.is_neutral(is)) then
            wrk = ua(:,is)*sqrt(geo%cvBb(:,3))
            call grad_p(nCv, nFc, nVx, 0, geo, mpg, wrk, wrkVx, wrk1)
            dvpar = wrk0*vsaf_cl(:,0,is)*wrk1
          else
            dvpar = 0.0_R8
          endif                                                          !srv 24.06.99

*       ..components of the viscosity parallel current
          fchvispar_a(:,0,is) =
     &     -switch%fhe_vis_par*fac_vis*dvpar*gOnedBsq(:,1)
     &      *geo%fcS*geo%fcQalf(:,0)
          fchvispar_a(:,1,is) = 
     &      switch%fhe_vis_par*fac_vis*dvpar*gOnedBsq(:,0)
     &      *geo%fcS*geo%fcQalf(:,1)

*       ..accumulate individual contributions into fchvispar
          fchvispar = fchvispar + fchvispar_a(:,:,is)

        enddo                                                            !srv 24.06.99

*     ..set zero outside of core
        do iFc = 1, nFc                                                  !lkw 20.06.2022{
          if ( .not. mpg%cvOnClosedSurface(mpg%fcCv(iFc,1))  .and.
     &         .not. mpg%cvOnClosedSurface(mpg%fcCv(iFc,2)))  then
            fchvispar(iFc,0:1) = 0.0_R8
            fchvispar_a(iFc,0:1,0:ns-1) = 0.0_R8
          endif
        enddo                                                            !lkw 20.06.2022}

      else
        fchvispar = 0.0_R8
        fchvispar_a = 0.0_R8
      endif                                                              !srv

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvspa

      end subroutine b2tvspa                                            !srv 20.09.06 }
