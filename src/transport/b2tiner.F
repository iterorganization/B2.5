*-----------------------------------------------------------------------
*.specification

      subroutine b2tiner (nCv, nFc, nVx, ns, geo, mpg, na, ua,
     &  facdrift, fchinert)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
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
*     B2TINER computes the X,Y-COMPONENTS OF INERTIAL CURRENT
*
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants
      
*   ..local variables
      integer ncall, is
      integer, save :: no_current = 0
      integer, save :: guard_flows
      real (kind=R8) ::
     &  fch_inert, wght(nFc,2), gOnedBsq(nFc,0:1),
     &  wrk(nCv), wrk0(nFc), wrk1(nFc)
      save ncall, fch_inert
*   ..procedures
      intrinsic max, abs
      external subini, subend, xertst
      external b2xvsg, b2xvff, b2xvfx
*   ..initialization
      data ncall/0/, fch_inert/1.0e0_R8/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tiner')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tiner_inert', fch_inert)                           !srv 07.07.99
       call ipgeti ('b2tfhe_no_current', no_current)
      endif
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
      endif
c Drift terms
      if(fch_inert.ne.0.0e0_R8.and.maxval(facdrift).ne.0.0_R8.and.
     .   no_current.eq.0)then                                            !srv 07.07.99

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
        fchinert(:,0) =  fch_inert*facdrift*wrk0*
     &    gOnedBsq(:,1)*geo%fcBb(:,2)*geo%fcS*geo%fcQalf(:,0)
        fchinert(:,1) = -fch_inert*facdrift*wrk0*
     &    gOnedBsq(:,0)*geo%fcBb(:,2)*geo%fcS*geo%fcQalf(:,1)

      else
        fchinert =0.0e0_R8
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tiner

      end                                                              !srv 20.09.06 }
