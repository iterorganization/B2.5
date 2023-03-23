      subroutine b2tbfl(nCv, nFc, lcore, mpg, flux, tbfl, div)     !srv 12.02.16 {
      use b2mod_types
      use b2mod_boundary_namelist
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc
      logical lcore
      type (mapping), intent (in) :: mpg
*   ..output arguments (unspecified on entry)
      real (kind=R8) :: flux(nFc,0:1), tbfl(1:nbc), div
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TBFL computes the total flux through the boundaries.
*            if lcore is .true. then all boundaries are included
*            if lcore is .false. then the core boundary is not included
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer :: ib, ibc, iCv, iFc
      integer, save :: ncall = 0
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tbfl')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      div = 0.0_R8
      do ib = 1, nbc
        tbfl(ib) = 0.0e0_R8
        ! compute flux through boundary; positive if out of the domain
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                       ! number of the guard cell
          iFc = mpg%cvFc(mpg%cvFcP(iCv,1))                              ! number of the guard cell face
          if (lcore.or.(.not.mpg%cvOnClosedSurface(iCv))) then
            tbfl(ib) = tbfl(ib) + mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &        (flux(iFc,0) + flux(iFc,1))
          endif
        enddo
        ! add boundary flux to total divergence
        div = div + tbfl(ib)
      enddo
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tbfl
      end subroutine b2tbfl                                !srv 12.02.16 }
