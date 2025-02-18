*-----------------------------------------------------------------------
*.specification

      subroutine b2mwit (nput, nCv, ns, itim, itcnt, dv)
      use b2mod_types
      use b2mod_subsys
      use b2us_plasma
      use b2mod_ad
     , , only : ncall_b2mwit
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent(in) :: nput, nCv, ns, itim, itcnt
      type(B2Derivatives), intent(in) :: dv
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MWIT writes to an un*formatted file detailed monitoring output
*     for the inner iterations of the B2 solution process.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit (nput) data describing the numerical
*     convergence. Unit nput must be connected to an un*formatted file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nput - integer, input.
*     nput specifies a unit number for un*formatted output.
*
*     nCv, ns - integer, input.
*     See the calling routine for description.
*
*     itim - integer, input.
*     itim specifies the timestep counter. It is initialised to 0 and
*     incremented by 1 following each successful call to b2mndt.
*
*     itcnt - integer, input.
*     itcnt counts the number of inner iterations. It is initialised to
*     0 and incremented by 1 following each call to b2news.
*
*     dv - B2Derivatives data structures, input.
*     Contains the residuals and corrections to be printed out.
*
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer idum(0:9)
*   ..procedures
      external xertst, cfwuin, cfwure

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mwit')
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test nCv, ns
      call xertst (1.le.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test itim, itcnt
      call xertst (0.le.itim.and.0.le.itcnt,
     &  'faulty argument itim, itcnt')

* ..save detailed convergence data
*   ..save itim, itcnt
      idum(0) = itim
      idum(1) = itcnt
      call cfwuin (nput, 2, idum, 'itim,itcnt')
*   ..save residuals and corrections
      call cfwure (nput, nCv*ns, dv%resco, 'resco')
      call cfwure (nput, nCv*ns, dv%corpa, 'corpa')
      call cfwure (nput, nCv*ns, dv%resmo, 'resmo')
      call cfwure (nput, nCv*ns, dv%corua, 'corua')
      call cfwure (nput, nCv, dv%resmt, 'resmt')
      call cfwure (nput, nCv, dv%corut, 'corut')
      call cfwure (nput, nCv, dv%reshe, 'reshe')
      call cfwure (nput, nCv, dv%corte, 'corte')
      call cfwure (nput, nCv, dv%reshi, 'reshi')
      call cfwure (nput, nCv, dv%corti, 'corti')
      call cfwure (nput, nCv, dv%respo, 'respo')
      call cfwure (nput, nCv, dv%corpo, 'corpo')

* ..return
      ncall_b2mwit = ncall_b2mwit+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mwit

      end subroutine b2mwit
