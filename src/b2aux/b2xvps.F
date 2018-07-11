*-----------------------------------------------------------------------
*.specification

      subroutine b2xvps (nCv, nFc, ns, pl, dv)
      use b2mod_types
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent(in) :: nCv, nFc, ns
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2XVPS tests the validity of the state of the plasma.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine tests the validity of the state of the plasma:
*     densities and temperatures must be positive.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ncv, nfc, ns - integer, input.
*     state - B2State, input
*     See the calling routine for description.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is
*   ..procedures
      external subini, subend, xertst
      external b2xvsg

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2xvps')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')

* ..test state
*   ..test sign of ne, na, te, ti
      call b2xvsg (nCv, dv % ne, 1, 'ne', '.gt.')
      call b2xvsg (nCv*ns, pl % na, 1, 'na', '.gt.')
      call b2xvsg (nCv, pl % te, 1, 'te', '.gt.')
      call b2xvsg (nCv, pl % ti, 1, 'ti', '.gt.')
      call b2xvsg (nCv*ns, dv%kinrgy, 1, 'kinrgy', '.ge.')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xvps

      end
