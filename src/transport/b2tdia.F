      subroutine b2tdia(nCv, nFc, nVx, ns, switch, geo, mpg, dv)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_constants
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
*   ..input/output arguments
      type (B2Derivatives), intent (inout) :: dv
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TDIA computes the poloidal and radial components of the
*     divergent part of the diamagnetic current.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer, save :: ncall = 0
      real (kind=R8) ::
     &  facdriftm, wght(nFc,2), gOnedBsq(nFc,0:1), pzf(nFc)
*   ..procedures
      external xertst
      external b2xvsg
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tdia')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*   ..test nCv, nFc, ns
       call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
       call xertst (1.le.ns, 'faulty argument ns')
      endif
      
      facdriftm = maxval(dv%facdrift)
      if (facdriftm.ne.0.0_R8 .and. switch%dia_cur.ne.0.0_R8             !srv 16.10.17
     &                        .and. switch%no_current.eq.0) then

*     ..interpolate to faces
        wght = 1.0_R8
        call intface(nCv,nFc,mpg%fcCv,wght,dv%pz,pzf)

*     ..compute gradients of OnedBsq on cell faces
        call grad(nCv, nFc, nVx, 1, geo, mpg, geo%cvOnedBsq,
     &    geo%vxOnedBsq, gOnedBsq)

*     ..DIVERGENT PART OF DIAMAGNETIC CURRENT
        dv%fchdia(:,0) =  switch%dia_cur*dv%facdrift*pzf*
     &    gOnedBsq(:,1)*geo%fcBb(:,2)*geo%fcS*geo%fcQalf(:,0)
        dv%fchdia(:,1) = -switch%dia_cur*dv%facdrift*pzf*
     &    gOnedBsq(:,0)*geo%fcBb(:,2)*geo%fcS*geo%fcQalf(:,1)

      else
        dv%fchdia = 0.0_R8
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2tdia

      end subroutine b2tdia                                              !srv 16.10.17 }
