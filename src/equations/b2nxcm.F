*-----------------------------------------------------------------------
*.specification

      subroutine b2nxcm (nx, ny, ns, pccm, pcca)
      use b2mod_types
      use b2mod_math
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *   pccm(-1:nx,-1:ny,0:1)
*   ..input/output arguments
      real (kind=R8) ::
     *   pcca(-1:nx,-1:ny,0:1,0:ns-1)

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
      integer ncall, ix, iy, is
      save ncall
*   ..procedures
      external subini, subend, xertst
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2nxcm')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

*   ..update pcca
      do is = 0, ns-1
        do iy = -1, ny
         do ix = -1, nx
          pcca(ix,iy,0,is) = uxcm(pccm(ix,iy,0),pcca(ix,iy,0,is))
          pcca(ix,iy,1,is) = uxcm(pccm(ix,iy,1),pcca(ix,iy,1,is))
         enddo
        enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxcm

      end subroutine b2nxcm
