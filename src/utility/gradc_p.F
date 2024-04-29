*-----------------------------------------------------------------------
*.specification

      subroutine gradc_p(nx, ny, hx, fun, dfunp)
      use b2mod_types
      use b2mod_indirect
      implicit none
*   ..input arguments
      integer nx,ny
      real (kind=R8) ::
     *   hx(-1:nx,-1:ny), fun(-1:nx,-1:ny)
*   ..output arguments
      real (kind=R8) ::
     *   dfunp(-1:nx,-1:ny)
*-----------------------------------------------------------------------
*.documentation
*
*     GRADC_P computes gradients of a cell centered quantity fun in the 
*     poloidal direction, in cell centers. As long as 
*     the x-direction is aligned with the poloidal direction, these 
*     gradients are fairly trivial.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix,iy

*-----------------------------------------------------------------------
*.computation      

*   ..compute gradients
      do iy = -1,ny
        do ix = -1,nx
          if (leftix(ix,iy).eq.-2) then
cwdk for now: zero here; consider 1-sided differences if needed?
            dfunp(ix,iy) = 0.0_R8
          elseif (rightix(ix,iy).eq.nx+1) then
cwdk for now: zero here; consider 1-sided differences if needed?
            dfunp(ix,iy) = 0.0_R8
          else
cwdk treatment in N/S guard cells to be considered
            dfunp(ix,iy) = 
     &       (fun(rightix(ix,iy),rightiy(ix,iy)) - 
     &        fun(leftix(ix,iy),leftiy(ix,iy)))/
     &       (hx(ix,iy) + 0.5_R8*hx(rightix(ix,iy),rightiy(ix,iy)) + 
     &                    0.5_R8*hx(leftix(ix,iy),leftiy(ix,iy)))
          endif
        enddo
      enddo

*   ..return
      return

*-----------------------------------------------------------------------
*.end gradc_p

      end subroutine gradc_p
