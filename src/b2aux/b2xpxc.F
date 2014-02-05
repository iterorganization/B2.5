      subroutine b2xpxc (nx, ny, funx, fun)
      use b2mod_types
      use b2mod_indirect
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *   funx(-1:nx,-1:ny), fun(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*      B2XPXC is obsolete.
*      B2XPXC interpolates a function from x-faces to cell centers.
*     ------------------------------------------------------------------
      integer ix, iy
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpxc')
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         fun(ix,iy) = funx(rightix(ix,iy),rightiy(ix,iy))
        elseif(rightix(ix,iy).eq.nx+1) then
         fun(ix,iy) = funx(ix,iy)
        else
         fun(ix,iy) = 
     .    (funx(ix,iy)+funx(rightix(ix,iy),rightiy(ix,iy)))/2
        endif
       enddo
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end
