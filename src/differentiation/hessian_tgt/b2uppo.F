*-----------------------------------------------------------------------
*.specification

      subroutine b2uppo (nCv, rxf, corpo, po)
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv
      real (kind=R8) ::
     *   rxf, corpo(nCv)
*   ..input/output arguments
      real (kind=R8) ::
     *   po(nCv)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UPPO corrects the electric potential po.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
*   ..procedures
      external xertst
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2uppo')
*   ..test nCv
      call xertst (0.lt.nCv, 'faulty argument nCv')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     . 'faulty argument rxf')

* ..implement the correction
*   ..compute new po
      po = po + rxf*corpo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2uppo

      end subroutine b2uppo
