*-----------------------------------------------------------------------
*.specification

      subroutine b2tiner (nCv, nFc, nVx, ns, switch, geo, mpg, na, ua,
     &  facdrift, fchinert)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmfs
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
     &  na(nCv,0:ns-1), ua(nCv,0:ns-1), facdrift(nFc)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fchinert(nFc,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TINER computes the X and Y components of the inertial current.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      integer is
      real (kind=R8) :: wght(nFc,2), gOnedBsq(nFc,0:1),
     &  wrk(nCv), wrk0(nFc), wrk1(nFc)
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx
*   ..initialization
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tiner')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
c Drift terms
      if(switch%fch_inert.ne.0.0_R8.and.maxval(facdrift).ne.0.0_R8
     .   .and.switch%no_current.eq.0)then                         !srv 07.07.99

*     ..compute gradients of OnedBsq on cell faces
        call grad(nCv, nFc, nVx, 1, geo, mpg, geo%cvOnedBsq,
     &    geo%vxOnedBsq, gOnedBsq)

*     ..sum of inertia terms
        wght = 1.0_R8
        wrk0 = 0.0_R8
        do is = 0, ns-1
          if (.not.is_neutral(is)) then
            wrk = na(:,is)*ua(:,is)**2
            call intface(nCv,nFc,mpg%fcCv,wght,wrk,wrk1)
            wrk0 = wrk0 + am(is)*wrk1
          endif
        enddo
        wrk0 = 0.5_R8*mp*wrk0

*     ..COMPONENTS OF INERTIA CURRENT
        fchinert(:,0) =  switch%fch_inert*facdrift*wrk0*
     &    gOnedBsq(:,1)*geo%fcBb(:,2)*geo%fcS*geo%fcQalf(:,0)
        fchinert(:,1) = -switch%fch_inert*facdrift*wrk0*
     &    gOnedBsq(:,0)*geo%fcBb(:,2)*geo%fcS*geo%fcQalf(:,1)

      else
        fchinert = 0.0_R8
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tiner

      end subroutine b2tiner                                            !srv 20.09.06 }
