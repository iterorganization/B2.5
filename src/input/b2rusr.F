*-----------------------------------------------------------------------
*.specification

      subroutine b2rusr (nget)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nget
*   ..output arguments (unspecified on entry)
*     (in common)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2RUSR reads a table of source term data from an un*formatted
*     file.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine is intended to be part of the interface to a Monte
*     Carlo neutral particle code.
*
*     At present, the routine is empty.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nget - integer, input.
*     nget specifies a unit number for un*formatted input.
*
*     (parameters in common:)
*
*     /b2cmsr/ - output.
*     /b2cmsr/ specifies source terms from neutral gas model.
*     (Presently empty--to be used in conjunction with a Monte Carlo
*     neutral gas code.)
*
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
*     (none)
*   ..procedures
      external subini, subend, xertst

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2rusr')
*   ..test nget
      call xertst (1.le.nget, 'faulty argument nget')

* ..read source data
*     (to follow)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2rusr

      end subroutine b2rusr
