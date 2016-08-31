*-----------------------------------------------------------------------
*.specification

      subroutine diff_r(nx, ny, hx, hy, vol, qz, fun, dfunr)
      use b2mod_types
      use b2mod_indirect
      implicit none
*   ..input arguments
      integer nx,ny
      real (kind=R8) ::
     *   hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), vol(-1:nx,-1:ny),
     *   qz(-1:nx,-1:ny,0:1), fun(-1:nx,-1:ny)
*   ..output arguments
      real (kind=R8) ::
     *   dfunr(-1:nx,-1:ny,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*     DIFF_R computes differences of a cell centered quantity fun in the 
*     radial direction, on south and west faces of a cell. It takes into 
*     account the corrections for x-faces not perpendicular to the 
*     poloidal field (9-point stencil).
*     dfunr(:,:,0) contains the radial differences at x-faces
*     dfunr(:,:,1) contains the radial differences at y-faces
*
*     Note: dfunr will not be computed between guard cells, because 
*           corner values in guard cells may not be (well) defined
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix,iy
      real (kind=R8) ::
     *   func(-1:nx,-1:ny), qzf(-1:nx,-1:ny)

*   ..procedures
      external intcor, intfaceh, intfacev

*-----------------------------------------------------------------------
*.computation      
      

*   ..initialize to zero
      dfunr = 0.0_R8

*   ..interpolate fun to cell corners, using volume weighted averaging
      call intcor(nx,ny,vol,fun,func)

*   ..interpolate qz(:,:,0) to x-faces
      call intfaceh(nx,ny,leftix,leftiy,vol,qz(:,:,0),qzf)

*   ..compute radial difference at x-faces
      do iy = -1,ny
        do ix = -1,nx
          if(leftix(ix,iy).ne.-2 .and.
     &       bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
            dfunr(ix,iy,0) = 
     &       func(topix(ix,iy),topiy(ix,iy)) - func(ix,iy) -            ! contribution from y-direction
     &       (fun(ix,iy) - fun(leftix(ix,iy),leftiy(ix,iy)))*           ! contribution from x-direction (non-orthogonal grids)
     &       (hy(ix,iy) + hy(leftix(ix,iy),leftiy(ix,iy)))/
     &       (hx(ix,iy) + hx(leftix(ix,iy),leftiy(ix,iy)))*
     &       qzf(ix,iy)
          endif
        enddo
      enddo

*   ..interpolate qz(:,:,0) to y-faces
      call intfacev(nx,ny,bottomix,bottomiy,vol,qz(:,:,0),qzf)

*   ..compute difference at y-faces
      do iy = -1,ny
        do ix = -1,nx
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and. 
     &       bottomiy(ix,iy).ne.-2) then
            dfunr(ix,iy,1) = 
     &       fun(ix,iy) - fun(bottomix(ix,iy),bottomiy(ix,iy)) -        ! contribution from y-direction
     &       (func(rightix(ix,iy),rightiy(ix,iy)) - func(ix,iy))*       ! contribution from x-direction (non-orthogonal grids)
     &       (hy(ix,iy) + hy(bottomix(ix,iy),bottomiy(ix,iy)))/
     &       (hx(ix,iy) + hx(bottomix(ix,iy),bottomiy(ix,iy)))*
     &       qzf(ix,iy)
          endif
        enddo
      enddo

*   ..reset null regions
      call b2xxgs (nx, ny, 0.0_R8, dfunr(-1,-1,0), 1)
      call b2xxgs (nx, ny, 0.0_R8, dfunr(-1,-1,1), 2)

*   ..return
      return

*-----------------------------------------------------------------------
*.end diff_r

      end subroutine diff_r
