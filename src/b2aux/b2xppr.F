      subroutine b2xppr (nx, ny, ns, ne, na, te, ti, pr)
      use b2mod_types
      implicit none
      integer nx, ny, ns
      real (kind=R8) ::
     *  ne(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1), 
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), pr(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XPPR computes the total pressure, pr:
*       pr(,) = ne(,)*te(,) + (sum is :: na(,,is)*ti(,)) .
*     ------------------------------------------------------------------
      integer n2, ix, iy, is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xppr')
      n2 = (nx+2)*(ny+2)
      call sfill (n2, 0.0_R8, pr, 1)
      pr = ne*te
      do is = 0, ns-1
       do iy = -1, ny
        do ix = -1, nx
         pr(ix,iy) = pr(ix,iy)+na(ix,iy,is)*ti(ix,iy)
        enddo
       enddo
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xppr
