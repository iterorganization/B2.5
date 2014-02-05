      subroutine b2xppb (nx, ny, rzb, nb, te, ti, pb)
      use b2mod_types
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *  rzb(-1:nx,-1:ny), nb(-1:nx,-1:ny), te(-1:nx,-1:ny), 
     &  ti(-1:nx,-1:ny), pb(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XPPB computes the partial pressure, pb = nb*(rzb*te+ti).
*     ------------------------------------------------------------------
      integer ix, iy
*     ------------------------------------------------------------------
*$$$  call subini ('b2xppb')
      do iy = -1, ny
       do ix = -1, nx
        pb(ix,iy) = nb(ix,iy)*(rzb(ix,iy)*te(ix,iy)+ti(ix,iy))
       enddo
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end
