      subroutine b2tepsch(nCv, nFc, nVx, switch, geo, mpg,
     &    facdrift, cdde, te, ne, floe, cone, fhePSch)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &  facdrift(nFc), cdde(nFc,0:1), te(nCv), ne(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  floe(nFc,0:1), cone(nFc,0:1,0:2), fhePSch(nFc,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2ETPSCH computes the poloidal and radial components of the
*     Pfirsch-Schlueter electron heat flux, and the coefficients that
*     specify its linearization.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc
      integer, save :: ncall = 0
      real (kind=R8) ::
     &  facdriftm, nete(nCv), nete2(nCv), netev(nVx), nete2v(nVx),
     &  dnete2(nFc,0:1),
     &  floe0(nFc,0:1), cone0(nFc,0:1), cdde_eff(nFc,0:1),
     &  dumm0(nFc,0:1)
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tepsch')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')

*   ..set dummy arrays
      dumm0 = 0.0_R8

      fhePSch = 0.0_R8                                                   !srv 13.10.06
      cone    = 0.0_R8                                                   !srv 02.01.07
      floe    = 0.0_R8
      facdriftm = maxval(facdrift)
      if(facdriftm.ne.0.0_R8) then

*   ..precompute face and vertex values of nete2
        nete  = ne*te
        nete2 = ne*te*te
        call intvertex(nCv,nVx,mpg,geo%vxVol,nete,netev)
        call intvertex(nCv,nVx,mpg,geo%vxVol,nete2,nete2v)
        call diff(nCv, nFc, nVx, 1, geo, mpg, nete2, nete2v, dnete2)

*   ..compute P.Sch heat flux
        fhePSch(:,0) = -2.5_R8*switch%fhePSch*facdrift*cdde(:,0)*
     &    geo%fcQalf(:,0)*dnete2(:,1)
        fhePSch(:,1) =  2.5_R8*switch%fhePSch*facdrift*cdde(:,1)*
     &    geo%fcQalf(:,1)*dnete2(:,0)
*     ..linearization coefficients for P.Sch flows
        cdde_eff(:,0) =  5.0_R8*switch%fhePSch*facdrift*cdde(:,1)*
     &    geo%fcQalf(:,1)
        cdde_eff(:,1) = -5.0_R8*switch%fhePSch*facdrift*cdde(:,0)*
     &    geo%fcQalf(:,0)
        call calccoef(nCv,nFc,nVx,0,geo,dumm0,cdde_eff,floe0,cone0)
        do iFc = 1, nFc
          floe(iFc,0)   = 0.5_R8*floe0(iFc,0)*
     &      (nete(mpg%fcCv(iFc,1))+nete(mpg%fcCv(iFc,2))) -
     &       cone0(iFc,0)*
     &      (nete(mpg%fcCv(iFc,2))-nete(mpg%fcCv(iFc,1)))
          floe(iFc,1)   = 0.5_R8*floe0(iFc,1)*
     &      (netev(mpg%fcVx(iFc,1))+netev(mpg%fcVx(iFc,2))) -
     &       cone0(iFc,1)*
     &      (netev(mpg%fcVx(iFc,2))-netev(mpg%fcVx(iFc,1)))
!wdk      For now, directly add these contributions into cone(,,0),
!wdk      so they are directly added into the stencil.
!wdk      Might need to move back to cone(,,1) and dedicated treatment
!wdk      in b2ux9p if we decide to move away from inverse volume weighting
!wdk      to compute vertex values of drifts!
          cone(iFc,0,0) = 0.5_R8*cone0(iFc,0)*
     &      (nete(mpg%fcCv(iFc,1))+nete(mpg%fcCv(iFc,2))) -
     &       0.25_R8*floe0(iFc,0)*
     &      (nete(mpg%fcCv(iFc,2))-nete(mpg%fcCv(iFc,1)))
          cone(iFc,1,0) = 0.5_R8*cone0(iFc,1)*
     &      (netev(mpg%fcVx(iFc,1))+netev(mpg%fcVx(iFc,2))) -
     &       0.25_R8*floe0(iFc,1)*
     &      (netev(mpg%fcVx(iFc,2))-netev(mpg%fcVx(iFc,1)))
        enddo
      endif
* ..return
      ncall = ncall+1
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2tepsch
      end subroutine b2tepsch
