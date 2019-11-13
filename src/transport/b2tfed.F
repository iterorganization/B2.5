*-----------------------------------------------------------------------
*.specification

      subroutine b2tfed (nCv, nFc, nVx, geo, mpg, facdrift, fac_ExB,
     & ne, te, po, vedia, veecrb)
      use b2mod_types
      use b2mod_constants
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     * facdrift(nFc), fac_ExB(nFc), ne(nCv), te(nCv), po(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  vedia(nFc,0:1), veecrb(nFc,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFED computes the electron drift velocities.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall
      real (kind=R8) ::
     * poleldr, poldidr, radeldr, raddidr
      real (kind=R8) ::
     * pof(nFc), pov(nVx), dpo(nFc,0:1), dOnedBsq(nFc,0:1), nef(nFc),
     * wrk(nCv), wrk0(nFc), wght(nFc,2)
      save ncall, poleldr, poldidr, radeldr, raddidr
*   ..procedures
      intrinsic min
      external subini, subend, ipgeti, ipgetr, xertst
      external b2xvsg
*   ..initialisation
      data ncall/0/,
     * poleldr/1.0_R8/, poldidr/1.0_R8/,
     * radeldr/1.0_R8/, raddidr/1.0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfed')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgetr ('b2tfnb_poleldr', poleldr)                         !srv 25.04.00 {
        call ipgetr ('b2tfnb_radeldr', radeldr)
        call ipgetr ('b2tfnb_poldidr', poldidr)
        call ipgetr ('b2tfnb_raddidr', raddidr)                         !srv 25.04.00 }
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of te
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
      endif
*   ..test facdrift
      call xertst (0.0_R8.le.minval(facdrift) .and.
     &     1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')

c Drift terms
      vedia  = 0.0_R8                                                    !xpb
      veecrb = 0.0_R8                                                    !xpb
      if(maxval(facdrift).ne.0.0_R8 .or. maxval(fac_ExB).ne.0.0_R8) then !xpb
cwdk precomputation of gradients on faces
        call grad(nCv, nFc, nVx, 0, geo, mpg, po, pov, dpo)
        call grad(nCv, nFc, nVx, 1, geo, mpg, geo%cvOnedBsq, 
     &   geo%vxOnedBsq, dOnedBsq)
*     ..computation of ExB drift
cwdk    to be checked: now we use inverse volume weighting
cwdk                   (original code: simple averaging)
        wrk  = geo%cvBb(:,2)/geo%cvBb(:,3)**2
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,wrk,wrk0)
        veecrb(:,0) = -poleldr*fac_ExB*wrk0*dpo(:,1)
        veecrb(:,1) =  radeldr*fac_ExB*wrk0*dpo(:,0)
*     ..computation of grad-B drift
        wght = 1.0_R8
        wrk  = ne*te/qe
        call intface(nCv,nFc,mpg%fcCv,wght,wrk,wrk0)
        call intface(nCv,nFc,mpg%fcCv,wght,ne,nef)
cwdk    to be checked: interpolation of geometric factors (hy1)
cwdk                   now consistent between x and y faces
cwdk                   (and consistent with treatment of other drifts)
cwdk                   (implicitly through grad-routine)
        vedia(:,0) = -poldidr*facdrift*wrk0*
     &    geo%fcBb(:,2)*dOnedBsq(:,1)/nef
        vedia(:,1) =  raddidr*facdrift*wrk0*
     &    geo%fcBb(:,2)*dOnedBsq(:,0)/nef
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfed

      end subroutine b2tfed
