      subroutine b2xpmc (nx, ny, ub, pr, ro, mcb)
      use b2mod_types
      use b2mod_indirect
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *   ub(-1:nx,-1:ny), pr(-1:nx,-1:ny), ro(-1:nx,-1:ny),
     &  mcb(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XPMC computes the mach number of the parallel flow.
*     ------------------------------------------------------------------
      integer ix, iy
      intrinsic sqrt
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpmc')
      do iy = -1, ny
       do ix = -1,nx
          mcb(ix,iy) = ub(ix,iy)/ sqrt(pr(ix,iy)/ro(ix,iy))
       enddo
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpmc
