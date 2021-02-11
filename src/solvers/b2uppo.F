*-----------------------------------------------------------------------
*.specification

      subroutine b2uppo (nCv, rxf, corpo, po)
      use b2mod_types
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2uppo
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
      external subini, subend, xertst
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2uppo')
*   ..test nCv
      call xertst (0.le.nCv, 'faulty argument nCv')
*   ..test rxf
      call xertst (0.le.rxf.and.rxf.le.1, 'faulty argument rxf')
*   ..extensive tests on first few calls
      if (ncall_b2uppo.lt.3) then
      endif

* ..implement the correction
*   ..compute new po
      po = po + rxf*corpo

* ..return
      ncall_b2uppo = ncall_b2uppo+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2uppo

      end
