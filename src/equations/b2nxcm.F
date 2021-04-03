*-----------------------------------------------------------------------
*.specification

      subroutine b2nxcm (nCv, ns, pccm, pcca)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns
      real (kind=R8) ::
     *   pccm(nCv,0:1)
*   ..input/output arguments
      real (kind=R8) ::
     *   pcca(nCv,0:1,0:ns-1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXCM computes the adjustment to the coefficients pcca that
*     comes from the combined parallel momentum equation solver.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, is
*   ..procedures
      real (kind=R8) ::
     *   uxcm
      external subini, subend, xertst, uxcm
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2nxcm')
*   ..test nCv
      call xertst (0.le.nCv, 'faulty argument nCv')

*   ..update pcca
      do is = 0, ns-1
        do iCv = 1, nCv
          pcca(iCv,0,is) = uxcm(pccm(iCv,0),pcca(iCv,0,is))
          pcca(iCv,1,is) = uxcm(pccm(iCv,1),pcca(iCv,1,is))
        enddo
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxcm

      end
