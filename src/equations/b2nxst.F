*-----------------------------------------------------------------------
*.specification

      subroutine b2nxst (nCv, rob, smb, ro, smt)
      use b2mod_types
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

*   ..local variables
      integer ncall, iCv
      save ncall
*   ..procedures
      external subini, subend, xertst
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2nxst')
*   ..test nCv
      call xertst (0.le.nCv, 'faulty argument nCv')
*   ..update smt
      smt(:,1) = smt(:,1) + smb(:,1)
      smt(:,3) = smt(:,3) + smb(:,3)*rob/ro

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxst

      end
