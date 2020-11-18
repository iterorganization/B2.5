      subroutine b2txsx (nx, ny, vol, hx, fun, funsx)
      use b2mod_types
      use b2mod_indirect
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny),
     &  fun(-1:nx,-1:ny), funsx(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2TXSX computes the integrals over cell x-faces of a function
*     that is provided on the cells.
*!!   (Reconsider specification and use of this routine for the case
*     of a non-orthogonal grid.)
*     ------------------------------------------------------------------
      integer n2, ix, iy
      external subini, subend, b2xxgs
*     ------------------------------------------------------------------
      call subini ('b2txsx')
*   ..test nx, ny
      call xertst (0.le.nx .and. 0.le.ny, 'faulty argument nx,ny')
*   ..set n2
      n2=(nx+2)*(ny+2)
*   ..test sign of vol, hx
      call b2xvsg (n2, vol, 1, 'vol', '.gt.')
      call b2xvsg (n2, hx, 1, 'hx', '.gt.')
*   ..compute funsx
      do iy = -1, ny
       do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
           funsx(ix,iy) = 0
         else
           funsx(ix,iy) =
     &      (vol(ix,iy)*fun(leftix(ix,iy),leftiy(ix,iy))
     &      +vol(leftix(ix,iy),leftiy(ix,iy))*fun(ix,iy))/
     &      (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
         endif
       enddo
      enddo
*   ..reset null regions
      call b2xxgs (nx, ny, 0.0_R8, funsx, 1)
*   ..return
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2txsx
