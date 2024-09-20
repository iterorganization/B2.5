*-----------------------------------------------------------------------
*.specification

      subroutine b2nxst (nCv, rob, smb, ro, smt)
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv
      real (kind=R8) ::
     *   rob(nCv), smb(nCv,0:3), ro(nCv)
*   ..input/output arguments
      real (kind=R8) ::
     *   smt(nCv,0:3)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXST accumulates the parallel momentum source linearisation over
*     the distinct atomic species.
*
*-----------------------------------------------------------------------
*.declarations

*   ..procedures
      external xertst
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2nxst')
*   ..test nCv
      call xertst (0.lt.nCv, 'faulty argument nCv')
*   ..update smt
      smt(:,1) = smt(:,1) + smb(:,1)
      smt(:,3) = smt(:,3) + smb(:,3)*rob/ro

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxst

      end subroutine b2nxst
