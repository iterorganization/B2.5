*-----------------------------------------------------------------------
*.specification

      subroutine b2uppo (nx, ny, rxf, corpo, po)
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     *   rxf, corpo(-1:nx,-1:ny)
*   ..input/output arguments
      real (kind=R8) ::
     *   po(-1:nx,-1:ny)

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
      integer ncall, ix, iy
      save ncall
*   ..procedures
      external xertst
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2uppo')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test rxf
      call xertst (0.le.rxf.and.rxf.le.1, 'faulty argument rxf')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
      endif

* ..implement the correction
*   ..compute new po
      do iy = -1, ny
       do ix = -1, nx
        po(ix,iy) = po(ix,iy)+rxf*corpo(ix,iy)
       enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2uppo

      end
