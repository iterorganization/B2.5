*-----------------------------------------------------------------------
*.specification

      subroutine b2mxzr (nx, ny, ns, resco, resmo, resmt, reshe, reshi,
     &  respo)
      use b2mod_types
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
*   ..input/output arguments
      real (kind=R8) ::
     *  resco(-1:nx,-1:ny,0:ns-1), resmo(-1:nx,-1:ny,0:ns-1),
     &  resmt(-1:nx,-1:ny), reshe(-1:nx,-1:ny), reshi(-1:nx,-1:ny),
     &  respo(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MXZR sets the residuals to zero along the boundaries of the
*     mesh and on the special regions.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine sets the residuals to zero along the boundaries
*     and on the special regions. The residuals on the boundary cells
*     and on the special regions often do not have any physical
*     meaning, and it is desirable to eliminate them from the
*     computation of the norms and from the output.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, is, ix, iy
      save ncall
*   ..procedures
      external subini, subend, xertst, sfill
      external b2xvfx, b2xxgs
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mxzr')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')

*   ..zero-out parts of (na) residual
      do is = 0, ns-1
       call b2xxgs (nx, ny, 0.0_R8, resco(-1,-1,is), 0)
       call sfill (nx+2, 0.0_R8, resco(-1,-1,is), 1)
       call sfill (nx+2, 0.0_R8, resco(-1,ny,is), 1)
       call sfill (ny, 0.0_R8, resco(-1,0,is), nx+2)
       call sfill (ny, 0.0_R8, resco(nx,0,is), nx+2)
      enddo
*   ..zero-out parts of (ua) residual
      do is = 0, ns-1
       call b2xxgs (nx, ny, 0.0_R8, resmo(-1,-1,is), 0)
       call sfill (nx+2, 0.0_R8, resmo(-1,-1,is), 1)
       call sfill (nx+2, 0.0_R8, resmo(-1,ny,is), 1)
       call sfill (ny, 0.0_R8, resmo(-1,0,is), nx+2)
       call sfill (ny, 0.0_R8, resmo(nx,0,is), nx+2)
      enddo
      call b2xxgs (nx, ny, 0.0_R8, resmt, 0)
      call sfill (nx+2, 0.0_R8, resmt(-1,-1), 1)
      call sfill (nx+2, 0.0_R8, resmt(-1,ny), 1)
      call sfill (ny, 0.0_R8, resmt(-1,0), nx+2)
      call sfill (ny, 0.0_R8, resmt(nx,0), nx+2)
*   ..zero-out parts of (te) residual
      call b2xxgs (nx, ny, 0.0_R8, reshe, 0)
      call sfill (nx+2, 0.0_R8, reshe(-1,-1), 1)
      call sfill (nx+2, 0.0_R8, reshe(-1,ny), 1)
      call sfill (ny, 0.0_R8, reshe(-1,0), nx+2)
      call sfill (ny, 0.0_R8, reshe(nx,0), nx+2)
*   ..zero-out parts of (ti) residual
      call b2xxgs (nx, ny, 0.0_R8, reshi, 0)
      call sfill (nx+2, 0.0_R8, reshi(-1,-1), 1)
      call sfill (nx+2, 0.0_R8, reshi(-1,ny), 1)
      call sfill (ny, 0.0_R8, reshi(-1,0), nx+2)
      call sfill (ny, 0.0_R8, reshi(nx,0), nx+2)
*   ..zero-out parts of (po) residual
      call b2xxgs (nx, ny, 0.0_R8, respo, 0)
      call sfill (nx+2, 0.0_R8, respo(-1,-1), 1)
      call sfill (nx+2, 0.0_R8, respo(-1,ny), 1)
      call sfill (ny, 0.0_R8, respo(-1,0), nx+2)
      call sfill (ny, 0.0_R8, respo(nx,0), nx+2)

      do ix = 0, nx-1
        do iy = 0, ny-1
*   ..zero out interior boundaries and resignore regions
          if(resignore(ix,iy,1).eq.0 .or. leftix(ix,iy).eq.-2 .or.
     &       rightix(ix,iy).eq.nx+1) then
            do is = 0, ns-1
              resco(ix,iy,is) = 0.0_R8
              resmo(ix,iy,is) = 0.0_R8
            enddo
            resmt(ix,iy) = 0.0_R8
            reshe(ix,iy) = 0.0_R8
            reshi(ix,iy) = 0.0_R8
            respo(ix,iy) = 0.0_R8
          endif
        enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mxzr

      end
