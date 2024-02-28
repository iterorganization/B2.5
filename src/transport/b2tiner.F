*-----------------------------------------------------------------------
*.specification

      subroutine b2tiner (nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &  na, ua, facdrift, fchinert_a, fchinert)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2tiner
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &  na(nCv,0:ns-1), ua(nCv,0:ns-1), facdrift(nFc)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fchinert_a(nFc,0:1,0:ns-1), fchinert(nFc,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TINER computes the poloidal and radial components of the inertial
*     current.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      integer is
      real (kind=R8) :: wght(nFc,2), gOnedBsq(nFc,0:1),
     &  wrk(nCv), wrkf(nFc)
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
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall_b2tiner.lt.3) then
      endif
c Drift terms
      if(switch%fch_inert.ne.0.0e0_R8.and.maxval(facdrift).ne.0.0_R8
     .   .and.switch%no_current.eq.0)then                                       !srv 07.07.99

*     ..compute gradients of OnedBsq on cell faces
        call grad(nCv, nFc, nVx, 1, geo, mpg, geo%cvOnedBsq,
     &    geo%vxOnedBsq, gOnedBsq)

*     ..compute inertia current
        wght = 1.0_R8
        fchinert = 0.0_R8
        do is = 0, ns-1
          if (.not.is_neutral(is)) then

*         ..interpolate inertia factor to faces
            wrk = 0.5_R8*am(is)*mp*na(:,is)*ua(:,is)**2
            call intface(nCv,nFc,mpg%fcCv,wght,wrk,wrkf)

*         ..components of inertia current per species
            fchinert_a(:,0,is) =  switch%fch_inert*facdrift*wrkf*
     &        gOnedBsq(:,1)*geo%fcBb(:,2)*geo%fcS*geo%fcQalf(:,0)
            fchinert_a(:,1,is) = -switch%fch_inert*facdrift*wrkf*
     &        gOnedBsq(:,0)*geo%fcBb(:,2)*geo%fcS*geo%fcQalf(:,1)

*         ..accumulate individual contributions into fchinert
            fchinert = fchinert + fchinert_a(:,:,is)

          endif
        enddo

*     ..fully assign to main ion species if needed
        if (switch%inert_all_ions.eq.0) then
          do is = 0, ns-1
            if (is.eq.ismain) then
              fchinert_a(:,:,is) = fchinert
            else
              fchinert_a(:,:,is) = 0.0_R8
            endif
          enddo
        endif

      else
        fchinert = 0.0e0_R8
        fchinert_a = 0.0e0_R8
      endif

* ..return
      ncall_b2tiner = ncall_b2tiner+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tiner

      end subroutine b2tiner                                            !srv 20.09.06 }
