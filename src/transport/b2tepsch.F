      subroutine b2tepsch(nCv, nFc, nVx, ns, geo, mpg, facdrift, cddi,
     &    te, ne, floe, cone, fhePSch)
      use b2mod_types
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &  facdrift(nFc), cddi(nFc,0:1,0:ns-1), te(nCv), ne(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  floe(nFc,0:1), cone(nFc,0:1,0:2), fhePSch(nFc,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2ETPSCH computes the X and Y components of the Pfirsch-Schlueter
*     heat flux
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc, is
      integer, save :: ncall = 0
      real (kind=R8), save :: PSch=1.0e0_R8                              !srv 17.01.07
      real (kind=R8) ::
     &  facdriftm, nete2(nCv), nete2v(nVx), dnete2(nFc,0:1), nev(nVx),
     &  floe0(nFc,0:1), cone0(nFc,0:1), cddi_eff(nFc,0:1),
     &  dumm0(nFc,0:1)
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tepsch')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tfhe_PSch', PSch)                                 !srv 17.01.07
      end if
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')   
      
*   ..set dummy arrays
      dumm0 = 0.0_R8

      fhePSch   = 0.0e0_R8                                               !srv 13.10.06
      cone      = 0.0e0_R8                                               !srv 02.01.07
      floe      = 0.0e0_R8
      facdriftm = maxval(facdrift)
      if(facdriftm.ne.0.0e0_R8) then
       do is = 0,ns-1                                                    !srv 12.05.11
        if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then              !srv 12.05.11

*   ..precompute face and vertex values of nete2
         nete2 = ne*te*te
         call intvertex(nCv,nVx,mpg,geo%vxVol,nete2,nete2v)
         call diff(nCv, nFc, nVx, 1, geo, mpg, nete2, nete2v, dnete2)
         call intvertex(nCv,nVx,mpg,geo%vxVol,ne,nev)

*   ..compute P.Sch heat flux
         fhePSch(:,0) =  2.5_R8*PSch*facdrift*cddi(:,0,is)*dnete2(:,1)
         fhePSch(:,1) = -2.5_R8*PSch*facdrift*cddi(:,1,is)*dnete2(:,0)
*     ..linearization coefficients for P.Sch flows
         cddi_eff(:,0) = -5.0_R8*PSch*facdrift*cddi(:,1,is)
         cddi_eff(:,1) =  5.0_R8*PSch*facdrift*cddi(:,0,is)
         call calccoef(nCv,nFc,nVx,0,geo,mpg,dumm0,cddi_eff,floe0,cone0)
         do iFc = 1, nFc
          floe(iFc,0)   = 0.5_R8*floe0(iFc,0)*
     &      (ne(mpg%fcCv(iFc,1))+ne(mpg%fcCv(iFc,2))) -
     &       cone0(iFc,0)*
     &      (ne(mpg%fcCv(iFc,2))-ne(mpg%fcCv(iFc,1)))
          floe(iFc,1)   = 0.5_R8*floe0(iFc,1)*
     &      (nev(mpg%fcVx(iFc,1))+nev(mpg%fcVx(iFc,2))) -
     &       cone0(iFc,1)*
     &      (nev(mpg%fcVx(iFc,2))-nev(mpg%fcVx(iFc,1)))
!wdk      For now, directly add these contributions into cone(,,0),
!wdk      so they are directly added into the stencil.
!wdk      Might need to move back to cone(,,1) and dedicated treatment
!wdk      in b2ux9p if we decide to move away from inverse volume weighting
!wdk      to compute vertex values of drifts!
          cone(iFc,0,0) = 0.5_R8*cone0(iFc,0)*
     &      (ne(mpg%fcCv(iFc,1))+ne(mpg%fcCv(iFc,2))) -
     &       0.25_R8*floe0(iFc,0)*
     &      (ne(mpg%fcCv(iFc,2))-ne(mpg%fcCv(iFc,1)))
          cone(iFc,1,0) = 0.5_R8*cone0(iFc,1)*
     &      (nev(mpg%fcVx(iFc,1))+nev(mpg%fcVx(iFc,2))) -
     &       0.25_R8*floe0(iFc,0)*
     &      (nev(mpg%fcVx(iFc,2))-nev(mpg%fcVx(iFc,1)))
         enddo
        endif
       enddo
      endif
* ..return
      ncall = ncall+1
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2tepsch
      end subroutine b2tepsch
