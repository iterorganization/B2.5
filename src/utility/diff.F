*-----------------------------------------------------------------------
*.specification

      subroutine diff(nx, ny, hx, hy, vol, qz, fun, dfun)
      use b2mod_types
      implicit none
*   ..input arguments
      integer nx,ny
      real (kind=R8) ::
     *   hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), vol(-1:nx,-1:ny),
     *   qz(-1:nx,-1:ny,0:1), fun(-1:nx,-1:ny)
*   ..output arguments
      real (kind=R8) ::
     *   dfun(-1:nx,-1:ny,0:1,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*     DIFF computes differences of a cell centered quantity fun on x and
*     y faces of a cell.
*     dfun(:,:,0,:) will hold differences at x-faces
*     dfun(:,:,1,:) will hold differences at y-faces
*     dfun(:,:,:,0) will hold poloidal differences
*     dfun(:,:,:,1) will hold radial differences
*
*     The routine is a simple interface to DIFF_P and DIFF_R, that will
*     compute the actual poloidal and radial differences
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      real (kind=R8) ::
     *  func(-1:nx,-1:ny)

*   ..procedures
      external diff_p, diff_r

*-----------------------------------------------------------------------
*.computation      

*   ..interpolate fun to cell corners
      call intcor(nx,ny,vol,fun,func)

*   ..compute poloidal and radial differences
      call diff_p(nx,ny, 1, hx, hy, vol, qz, fun, func, dfun(-1,-1,0,0))
      call diff_r(nx,ny, 1, hx, hy, vol, qz, fun, func, dfun(-1,-1,0,1))

*   ..return
      return

*-----------------------------------------------------------------------
*.end diff

      end subroutine diff
