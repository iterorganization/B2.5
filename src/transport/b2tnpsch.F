      subroutine b2tnpsch(nCv, nFc, nVx, isb, switch, geo, mpg,
     &    facdrift, cddb, ti, nb, kmprti, flob, conb, fnbPSch)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, isb
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &  facdrift(nFc), cddb(nFc,0:1), ti(nCv), nb(nCv), kmprti(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  flob(nFc,0:1), conb(nFc,0:1,0:2), fnbPSch(nFc,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TNPSCH computes the poloidal and radial components of the
*     Pfirsch-Schlueter particle flux of one atomic species, and the
*     coefficients that specify its linearization.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc
      integer, save :: ncall = 0
      real (kind=R8) ::
     &  facdriftm, nbti(nCv), nbtiv(nVx), kmprtiv(nVx),
     &  dnbti(nFc,0:1),
     &  flob0(nFc,0:1), conb0(nFc,0:1), cddb_eff(nFc,0:1),
     &  dumm0(nFc,0:1)
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tnpsch')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')   

*   ..set dummy arrays
      dumm0 = 0.0_R8

      fnbPSch = 0.0e0_R8                                                 !srv 13.10.06
      conb    = 0.0e0_R8                                                 !srv 02.01.07
      flob    = 0.0e0_R8
      facdriftm = maxval(facdrift)
      if(facdriftm.ne.0.0e0_R8) then

*     ..precompute face and vertex values of nbti and kmprti
        nbti  = nb*ti
        call intvertex(nCv,nVx,mpg,geo%vxVol,nbti,nbtiv)
        call intvertex(nCv,nVx,mpg,geo%vxVol,kmprti,kmprtiv)
        call diff(nCv, nFc, nVx, 1, geo, mpg, nbti, nbtiv, dnbti)

*     ..compute P.Sch heat flux
        fnbPSch(:,0) =  switch%fnbPSch*facdrift*cddb(:,0)*dnbti(:,1)
        fnbPSch(:,1) = -switch%fnbPSch*facdrift*cddb(:,1)*dnbti(:,0)

*     ..linearization coefficients for P.Sch flows
        cddb_eff(:,0) = -switch%fnbPSch*facdrift*cddb(:,1)
        cddb_eff(:,1) =  switch%fnbPSch*facdrift*cddb(:,0)
        call calccoef(nCv,nFc,nVx,0,geo,dumm0,cddb_eff,flob0,conb0)
        do iFc = 1, nFc
          flob(iFc,0)   = 0.5_R8*flob0(iFc,0)*
     &      (kmprti(mpg%fcCv(iFc,1)) + kmprti(mpg%fcCv(iFc,2))) -
     &       conb0(iFc,0)*
     &      (kmprti(mpg%fcCv(iFc,2)) - kmprti(mpg%fcCv(iFc,1)))
          flob(iFc,1)   = 0.5_R8*flob0(iFc,1)*
     &      (kmprtiv(mpg%fcVx(iFc,1))+kmprtiv(mpg%fcVx(iFc,2))) -
     &       conb0(iFc,1)*
     &      (kmprtiv(mpg%fcVx(iFc,2))-kmprtiv(mpg%fcVx(iFc,1)))
!wdk      For now, directly add these contributions into conb(,,0),
!wdk      so they are directly added into the stencil.
!wdk      Might need to move back to conb(,,3) and dedicated treatment
!wdk      in b2ux9p if we decide to move away from inverse volume weighting
!wdk      to compute vertex values of drifts!
          conb(iFc,0,0) = 0.5_R8*conb0(iFc,0)*
     &      (kmprti(mpg%fcCv(iFc,1))+kmprti(mpg%fcCv(iFc,2))) -
     &       0.25_R8*flob0(iFc,0)*
     &      (kmprti(mpg%fcCv(iFc,2))-kmprti(mpg%fcCv(iFc,1)))
          conb(iFc,1,0) = 0.5_R8*conb0(iFc,1)*
     &      (kmprtiv(mpg%fcVx(iFc,1))+kmprtiv(mpg%fcVx(iFc,2))) -
     &       0.25_R8*flob0(iFc,1)*
     &      (kmprtiv(mpg%fcVx(iFc,2))-kmprtiv(mpg%fcVx(iFc,1)))
        enddo
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2tnpsch
      end subroutine b2tnpsch
