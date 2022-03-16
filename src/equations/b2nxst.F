*-----------------------------------------------------------------------
*.specification

      subroutine b2nxst (nx, ny, rob, smb, ro, smt)
      use b2mod_types
      use b2mod_indirect
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     *   rob(-1:nx,-1:ny), smb(-1:nx,-1:ny,0:3), ro(-1:nx,-1:ny)
*   ..input/output arguments
      real (kind=R8) ::
     *   smt(-1:nx,-1:ny,0:3)

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
      integer ncall, ix, iy
      !$OMP THREADPRIVATE(ncall)
      save ncall
*   ..procedures
      external xertst
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2nxst')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..update smt
      do iy = -1, ny
       do ix = -1, nx
        smt(ix,iy,1) = smt(ix,iy,1)+smb(ix,iy,1)
        smt(ix,iy,3) = smt(ix,iy,3)+smb(ix,iy,3)*
     &     rob(ix,iy)/ro(ix,iy)
       enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxst

      end
