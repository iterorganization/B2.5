*-----------------------------------------------------------------------
*.specification

      subroutine grad(nx, ny, hx, hy, hy1, vol, qs, fun, gfun)
      use b2mod_types
      implicit none
*   ..input arguments
      integer nx,ny
      real (kind=R8) ::
     *   hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), hy1(-1:nx,-1:ny),
     *   vol(-1:nx,-1:ny),
     *   qs(-1:nx,-1:ny,0:1), fun(-1:nx,-1:ny)
*   ..output arguments
      real (kind=R8) ::
     *   gfun(-1:nx,-1:ny,0:1,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*     GRAD computes gradients of a cell centered quantity fun on x and
*     y faces of a cell.
*     gfun(:,:,0,:) will hold gradients at x-faces
*     gfun(:,:,1,:) will hold gradients at y-faces
*     gfun(:,:,:,0) will hold poloidal gradients
*     gfun(:,:,:,1) will hold radial gradients
*
*     The routine is a simple interface to GRAD_P and GRAD_R, that will
*     compute the actual poloidal and radial gradients
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      real (kind=R8) ::
     *  func(-1:nx,-1:ny)

*   ..procedures
      external grad_p, grad_r

*-----------------------------------------------------------------------
*.computation

*   ..interpolate fun to cell corners
      call intcor(nx,ny,vol,fun,func)

*   ..compute poloidal and radial gradients
      call grad_p(nx, ny, 1, hx, vol, fun, func, gfun(-1,-1,0,0))
      call grad_r(nx, ny, 1, hx, hy, hy1, vol, qs, fun, func,
     &  gfun(-1,-1,0,1))

*   ..return
      return

*-----------------------------------------------------------------------
*.end grad

      end subroutine grad
