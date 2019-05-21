*-----------------------------------------------------------------------
*.specification

      subroutine b2tfcc (nCv, nFc, nVx, ns, geo, mpg, pl, dv, co, rt)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type(B2Plasma), intent (in) :: pl 
      type(B2Coeff), intent (in) :: co
      type(B2Rates), intent (in) :: rt
*   ..input/output arguments
      type(B2Derivatives), intent (inout) :: dv
*   ..common blocks
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, is
      real (kind=R8) ::
     *   xfac
      real (kind=R8) ::
     *   nwe(nCv), nwi(nCv), nwef(nFc), nwif(nFc),
     *   dnwe(nFc,0:1), dnwi(nFc,0:1), tef(nFc), tif(nFc),
     *   wght(nFc,2), wrkVx(nVx),
     *   floe(nFc,0:1), cone(nFc,0:1), floe0(nFc,0:1), cone0(nFc,0:1),
     *   floi(nFc,0:1), coni(nFc,0:1), floi0(nFc,0:1), coni0(nFc,0:1)
      save ncall, xfac
*   ..procedures
      external subini, subend, xertst, ipgetr
      external b2xvsg, b2xvff, diff, intface
*   ..initialisation
      data ncall/0/, xfac/1.0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfcc')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tfcc_xfac', xfac)
       call xertst (0.le.xfac, 'faulty internal parameter xfac')
      endif
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of cdpa
       do is = 0, ns-1
        call b2xvsg (nFc, co%cdpa(:,0,is)*geo%fcQalf(:,0),1,'cdpa0',
     .    '.ge.')
        call b2xvsg (nFc, co%cdpa(:,1,is)*geo%fcQalf(:,1),1,'cdpa1',
     .    '.ge.')
       enddo
      endif

*   ..interpolate te and ti to cell faces
      wght = 1.0_R8
      call intface(nCv,nFc,mpg%fcCv,wght,pl%te,tef)
      call intface(nCv,nFc,mpg%fcCv,wght,pl%ti,tif)

* ..correct floe, cone, floi, coni

      do is = 0, ns-1

cwdk   Note: could reduce memory by doing ion and electron contributions
cwdk   subsequently (use same work arrays nw, dnw, nwf, flo, con, flo0, con0)

*    ..pre-compute radial differences
       nwi = pl%na(:,is)*(rt%rza(:,is) + 1.0_R8)
       nwe = pl%na(:,is)*(rt%rza(:,is) + 1.0_R8)*rt%rza(:,is)
       call diff(nCv, nFc, nVx, 0, geo, mpg, nwi, wrkVx, dnwi)
       call diff(nCv, nFc, nVx, 0, geo, mpg, nwe, wrkVx, dnwe)

*    ..interpolate nwe and nwi to cell faces
       call intface(nCv,nFc,mpg%fcCv,wght,nwe,nwef)
       call intface(nCv,nFc,mpg%fcCv,wght,nwi,nwif)

*    ..compute poloidal/radial coefficients for pressure diffusion term
       floe(:,0) = -xfac*2.5_R8*co%cdpa(:,0,is)*dnwe(:,0)*tef
       floe(:,1) = -xfac*2.5_R8*co%cdpa(:,1,is)*dnwe(:,1)*tef
       cone(:,0) =  xfac*2.5_R8*co%cdpa(:,0,is)*nwef*tef
       cone(:,1) =  xfac*2.5_R8*co%cdpa(:,1,is)*nwef*tef

       floi(:,0) = -xfac*2.5_R8*co%cdpa(:,0,is)*dnwi(:,0)*tif
       floi(:,1) = -xfac*2.5_R8*co%cdpa(:,1,is)*dnwi(:,1)*tif
       coni(:,0) =  xfac*2.5_R8*co%cdpa(:,0,is)*nwif*tif
       coni(:,1) =  xfac*2.5_R8*co%cdpa(:,1,is)*nwif*tif

*    ..compute coefficients for stencil
       call calccoef(nCv,nFc,nVx,0,geo,mpg,floe,cone,floe0,cone0)
       call calccoef(nCv,nFc,nVx,0,geo,mpg,floi,coni,floi0,coni0)

*    ..add to previous contributions
       dv%floe        = dv%floe        + floe0
       dv%cone(:,:,0) = dv%cone(:,:,0) + cone0

       dv%floi        = dv%floi        + floi0
       dv%coni(:,:,0) = dv%coni(:,:,0) + coni0

      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfcc

      end
