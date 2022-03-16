      subroutine b2txsy (nx, ny, vol, hy, qz, fun, funsy)
      use b2mod_types
      use b2mod_indirect
      use b2mod_subsys
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *   vol(-1:nx,-1:ny), hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1),
     &  fun(-1:nx,-1:ny), funsy(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2TXSY computes the integrals over cell y-faces of a function
*     that is provided on the cells.
*!!   (Re-consider specification and use of this routine for the case
*     of a non-orthogonal grid.)
*     ------------------------------------------------------------------
      integer ix, iy
      real (kind=R8) ::
     *   hy1
      external b2xxgs
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
*     ------------------------------------------------------------------
      call subini ('b2txsy')
*   ..compute funsy
      do iy = -1, ny
       do ix = -1, nx
        if(bottomiy(ix,iy).eq.-2) then
         funsy(ix,iy) = 0
        else
         funsy(ix,iy) =
     &    (vol(ix,iy)*fun(bottomix(ix,iy),bottomiy(ix,iy))
     &    +vol(bottomix(ix,iy),bottomiy(ix,iy))*fun(ix,iy))/
     &    (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
        endif
       enddo
      enddo
*   ..reset null regions
      call b2xxgs (nx, ny, 0.0_R8, funsy, 2)
*   ..return
      call subend ()
      return
*     ------------------------------------------------------------------
      end
