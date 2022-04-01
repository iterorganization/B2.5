      subroutine b2xtvx (nx, ny, nx1, ny1, fun, fun1, wrk)
      use b2mod_types
      use b2mod_math
      implicit none
      integer nx, ny, nx1, ny1
      real (kind=R8) ::
     *   fun(-1:nx,-1:ny), fun1(-1:nx1,-1:ny1), wrk(-1:nx,-1:ny1)
*     ------------------------------------------------------------------
*     B2XTVX transfers an x-face variable to a new mesh.
*     ------------------------------------------------------------------
      integer ix1, iy1, ix, iy
      real (kind=R8) ::
     *   t0, t1, t2
      external b2xvfx, rratio
      real (kind=R8) ::
     *   rratio
*     ------------------------------------------------------------------
*   ..test edge values of fun
      call b2xvfx (nx, ny, fun, 'b2xtvx--fun')
*   ..interpolate in the second index
      do ix = 0, nx
       wrk(ix,-1) = fun(ix,-1)
      enddo
      do iy1 = 0, ny1-1
       iy = (2*iy1*ny+ny)/(2*ny1)
*      (iy specifies the nearest point)
       if ((2*iy+1)*ny1.eq.(2*iy1+1)*ny) then
        do ix = 0, nx
         wrk(ix,iy1) = fun(ix,iy)
        enddo
       else if (2*iy*ny1.eq.(2*iy1+1)*ny) then
        do ix = 0, nx
           wrk(ix,iy1) = (fun(ix,iy-1)+fun(ix,iy))/2
        enddo
       else
        t0 = rratio((2*iy+1)*ny1-(2*iy1+1)*ny,2*ny1)
        do ix = 0, nx
         t1 = t0*fun(ix,iy-1)+(1-t0)*fun(ix,iy)
         t2 = (1+t0)*fun(ix,iy)-t0*fun(ix,iy+1)
         wrk(ix,iy1) = med(t1,t2,fun(ix,iy))
        enddo
       endif
      enddo
      do ix = 0, nx
       wrk(ix,ny1) = fun(ix,ny)
      enddo
*   ..interpolate in the first index
      do iy1 = -1, ny1
       fun1(-1,iy1) = 0
      enddo
      do ix1 = 0, nx1
       ix = (ix1*nx)/nx1
*      (ix specifies the lower neighbour)
       if (ix*nx1.eq.ix1*nx) then
        do iy1 = -1, ny1
         fun1(ix1,iy1) = wrk(ix,iy1)
        enddo
       else
        t0 = rratio(nx*ix1-ix*nx1,nx1)
        do iy1 = -1, ny1
         fun1(ix1,iy1) = (1-t0)*wrk(ix,iy1)+t0*wrk(ix+1,iy1)
        enddo
       endif
      enddo
      return
*     ------------------------------------------------------------------
      end subroutine b2xtvx
