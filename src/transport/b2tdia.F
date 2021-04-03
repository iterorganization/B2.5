      subroutine b2tdia(nCv, nFc, nVx, ns, switch, geo, mpg,
     &   pz, fna, uadia, facdrift, fchdia)       !srv 16.19.17 {
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_constants
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
     &  pz(nCv), facdrift(nFc), fna(nFc,0:1,0:ns-1), 
     &  uadia(nFc,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fchdia(nFc,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TDIA computes the X and Y components of the divergent part of
*     the diamagnetic current
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc
      integer, save :: ncall = 0
      real (kind=R8) ::
     &  facdriftm, wght(nFc,2), gOnedBsq(nFc,0:1), pzf(nFc)
*   ..procedures
      intrinsic max, abs
      real (kind=R8) ::
     *   smin, smax
      external subini, subend, xertst
      external b2xvsg, b2xvff, b2xvfx
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tdia')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*   ..test nx, ny, ns
       call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
       call xertst (1.le.ns, 'faulty argument ns')
      endif
      
      facdriftm = maxval(facdrift)
      if (facdriftm.ne.0.0_R8.and.switch%dia_cur.ne.0.0e0_R8) then            !srv 16.10.17

*     ..interpolate to faces
        wght = 1.0_R8
        call intface(nCv,nFc,mpg%fcCv,wght,pz,pzf)

*     ..compute gradients of OnedBsq on cell faces
        call grad(nCv, nFc, nVx, 1, geo, mpg, geo%cvOnedBsq,
     &    geo%vxOnedBsq, gOnedBsq)

*     ..DIVERGENT PART OF DIAMAGNETIC CURRENT
        fchdia(:,0) =  switch%dia_cur*facdrift*pzf*
     &    gOnedBsq(:,1)*geo%fcBb(:,2)*geo%fcS*geo%fcQalf(:,0)
        fchdia(:,1) = -switch%dia_cur*facdrift*pzf*
     &    gOnedBsq(:,0)*geo%fcBb(:,2)*geo%fcS*geo%fcQalf(:,1)

      else
        fchdia = 0.0_R8
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2tdia

      end subroutine b2tdia                                              !srv 16.10.17 }
