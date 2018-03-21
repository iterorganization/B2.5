*-----------------------------------------------------------------------
*.specification

      subroutine gradc_r(nx, ny, hx, hy, vol, qz, fun, dfunr)
      use b2mod_types
      use b2mod_indirect
      implicit none
*   ..input arguments
      integer nx,ny,mode
      real (kind=R8) ::
     *   hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), vol(-1:nx,-1:ny),
     *   qz(-1:nx,-1:ny,0:1), fun(-1:nx,-1:ny)
*   ..output arguments
      real (kind=R8) ::
     *   dfunr(-1:nx,-1:ny)
*-----------------------------------------------------------------------
*.documentation
*
*     GRADC_R computes gradients of a cell centered quantity fun in the 
*     radial direction, in cell centers. It takes into 
*     account the corrections for the y-direction not perpendicular to
*     the poloidal field (9-point stencil).
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix,iy

*   ..procedures
      real (kind=R8) :: hy1
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)

*-----------------------------------------------------------------------
*.computation      

*   ..compute radial gradients
      do iy = -1,ny
        do ix = -1,nx
          if (bottomiy(ix,iy).eq.-2) then
cwdk for now:  zero here; consider 1-sided differences if needed?
            dfunr(ix,iy) = 0.0_R8
          elseif (topiy(ix,iy).eq.ny+1) then
cwdk for now:  zero here; consider 1-sided differences if needed?
            dfunr(ix,iy) = 0.0_R8
          else
            if (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
cwdk for now: no 9-point corrections in these guard cells;
cwdk          consider 1-sided differences if needed?
              dfunr(ix,iy) = 
     &         (fun(topix(ix,iy),topiy(ix,iy)) - 
     &          fun(bottomix(ix,iy),bottomiy(ix,iy)))/
     &         (hy1(ix,iy) + 0.5_R8*hy1(topix(ix,iy),topiy(ix,iy)) +
     &          0.5_R8*hy1(bottomix(ix,iy),bottomiy(ix,iy)))
            else
              dfunr(ix,iy) = 
     &         (fun(topix(ix,iy),topiy(ix,iy)) -
     &          fun(bottomix(ix,iy),bottomiy(ix,iy)))/
     &         (hy1(ix,iy) + 0.5_R8*hy1(topix(ix,iy),topiy(ix,iy)) +
     &          0.5_R8*hy1(bottomix(ix,iy),bottomiy(ix,iy))) - 
     &          qz(ix,iy,0)/qz(ix,iy,1)*
     &         (fun(rightix(ix,iy),rightiy(ix,iy)) - 
     &          fun(leftix(ix,iy),leftiy(ix,iy)))/
     &         (hx(ix,iy) + 0.5_R8*hx(rightix(ix,iy),rightiy(ix,iy)) + 
     &                      0.5_R8*hx(leftix(ix,iy),leftiy(ix,iy)))
            endif
          endif
        enddo
      enddo

*   ..return
      return

*-----------------------------------------------------------------------
*.end gradc_r

      end subroutine gradc_r
