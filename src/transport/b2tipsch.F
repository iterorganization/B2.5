      subroutine b2tipsch(nCv, nFc, nVx, ns, switch, geo, mpg,
     &    facdrift, cddi, ti, na, floi, coni, fhiPSch)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &  facdrift(nFc), cddi(nFc,0:1,0:ns-1), ti(nCv), na(nCv,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  floi(nFc,0:1), coni(nFc,0:1,0:2), fhiPSch(nFc,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TIPSCH computes the poloidal and radial components of the
*     Pfirsch-Schlueter all atom heat flux, and the coefficients that
*     specify its linearization.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc, is
      integer, save :: ncall = 0
      real (kind=R8) ::
     &  facdriftm, nati(nCv), nati2(nCv), nativ(nVx), nati2v(nVx),
     &  dnati2(nFc,0:1),
     &  floi0(nFc,0:1), coni0(nFc,0:1), cddi_eff(nFc,0:1),
     &  dumm0(nFc,0:1)
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tipsch')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')   

*   ..set dummy arrays
      dumm0 = 0.0_R8

      fhiPSch = 0.0e0_R8                                                 !srv 13.10.06
      coni    = 0.0e0_R8                                                 !srv 02.01.07
      floi    = 0.0e0_R8
      facdriftm = maxval(facdrift)
      if(facdriftm.ne.0.0e0_R8) then
       do is = 0, ns-1
        if(.not.is_neutral(is)) then

*      ..precompute face and vertex values of nati and nati2
         nati  = na(:,is)*ti
         nati2 = na(:,is)*ti*ti
         call intvertex(nCv,nVx,mpg,geo%vxVol,nati,nativ)
         call intvertex(nCv,nVx,mpg,geo%vxVol,nati2,nati2v)
         call diff(nCv, nFc, nVx, 1, geo, mpg, nati2, nati2v, dnati2)

*      ..compute P.Sch heat flux
         fhiPSch(:,0) =  2.5_R8*switch%fhiPSch*facdrift*cddi(:,0,is)*
     &     dnati2(:,1)
         fhiPSch(:,1) = -2.5_R8*switch%fhiPSch*facdrift*cddi(:,1,is)*
     &     dnati2(:,0)

*      ..linearization coefficients for P.Sch flows
         cddi_eff(:,0) = -5.0_R8*switch%fhiPSch*facdrift*cddi(:,1,is)
         cddi_eff(:,1) =  5.0_R8*switch%fhiPSch*facdrift*cddi(:,0,is)
         call calccoef(nCv,nFc,nVx,0,geo,dumm0,cddi_eff,floi0,coni0)
         do iFc = 1, nFc
          floi(iFc,0)   = floi(iFc,0) +
     &       0.5_R8*floi0(iFc,0)*
     &      (nati(mpg%fcCv(iFc,1)) + nati(mpg%fcCv(iFc,2))) -
     &       coni0(iFc,0)*
     &      (nati(mpg%fcCv(iFc,2)) - nati(mpg%fcCv(iFc,1)))
          floi(iFc,1)   = floi(iFc,1) +
     &       0.5_R8*floi0(iFc,1)*
     &      (nativ(mpg%fcVx(iFc,1))+nativ(mpg%fcVx(iFc,2))) -
     &       coni0(iFc,1)*
     &      (nativ(mpg%fcVx(iFc,2))-nativ(mpg%fcVx(iFc,1)))
!wdk      For now, directly add these contributions into coni(,,0),
!wdk      so they are directly added into the stencil.
!wdk      Might need to move back to coni(,,1) and dedicated treatment
!wdk      in b2ux9p if we decide to move away from inverse volume weighting
!wdk      to compute vertex values of drifts!
          coni(iFc,0,0) = coni(iFc,0,0) +
     &       0.5_R8*coni0(iFc,0)*
     &      (nati(mpg%fcCv(iFc,1))+nati(mpg%fcCv(iFc,2))) -
     &       0.25_R8*floi0(iFc,0)*
     &      (nati(mpg%fcCv(iFc,2))-nati(mpg%fcCv(iFc,1)))
          coni(iFc,1,0) = coni(iFc,1,0) +
     &       0.5_R8*coni0(iFc,1)*
     &      (nativ(mpg%fcVx(iFc,1))+nativ(mpg%fcVx(iFc,2))) -
     &       0.25_R8*floi0(iFc,1)*
     &      (nativ(mpg%fcVx(iFc,2))-nativ(mpg%fcVx(iFc,1)))
         enddo
        endif
       enddo
      endif
* ..return
      ncall = ncall+1
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2tipsch
      end subroutine b2tipsch
