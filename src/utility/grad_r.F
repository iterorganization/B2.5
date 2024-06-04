*-----------------------------------------------------------------------
*.specification

      subroutine grad_r(nx, ny, mode, hx, hy, hy1, vol, qs, fun, func,
     *   gfunr)
      use b2mod_types
      use b2mod_indirect
      implicit none
*   ..input arguments
      integer nx,ny,mode
      real (kind=R8) ::
     *   hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), hy1(-1:nx,-1:ny),
     *   vol(-1:nx,-1:ny), qs(-1:nx,-1:ny,0:1),
     *   fun(-1:nx,-1:ny), func(-1:nx,-1:ny)
*   ..output arguments
      real (kind=R8) ::
     *   gfunr(-1:nx,-1:ny,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*     GRAD_R computes gradients of a cell centered quantity fun in the
*     radial direction, on x and y-faces of a cell. It takes into
*     account the corrections for x-faces not perpendicular to the
*     poloidal field (9-point stencil).
*     gfunr(:,:,0) will hold the radial gradients at x-faces
*     gfunr(:,:,1) will hold the radial gradients at y-faces
*
*     To compute values of the field in cell corners, the routine relies
*     calls to intcor, using the input array 'vol' as weights. Different
*     interpolation schemes can be mimicked by using a different array
*     with weights in the calling routine.
*
*     Note: gfunr will not be computed between guard cells, because
*           corner values in guard cells may not be (well) defined
*
*     If mode.eq.0, the routine will interpolate fun to corners in func
*     If mode.eq.1, the routine assumess func has been pre-computed
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix,iy

*   ..procedures
      external intcor, intfaceh, intfacev

*-----------------------------------------------------------------------
*.computation

*   ..test mode
      call xertst (0.le.mode.and.mode.le.1,
     .  'grad_r--faulty argument mode')

*   ..interpolate fun to cell corners, using volume weighted averaging
      if (mode.eq.0) then
        call intcor(nx,ny,vol,fun,func)
      endif

*   ..compute radial gradients at x-faces
      do iy = -1,ny
        do ix = -1,nx
          if(leftix(ix,iy).ne.-2 .and.
     &       bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
            gfunr(ix,iy,0) = 2.0_R8*
     &       (func(topix(ix,iy),topiy(ix,iy)) - func(ix,iy) -           ! contribution from y-direction
     &       (fun(ix,iy) - fun(leftix(ix,iy),leftiy(ix,iy)))*           ! contribution from x-direction (non-orthogonal grids)
     &       (hy(ix,iy) + hy(leftix(ix,iy),leftiy(ix,iy)))/
     &       (hx(ix,iy) + hx(leftix(ix,iy),leftiy(ix,iy)))*
     &       qs(ix,iy,0))/
     &       (hy1(ix,iy) + hy1(leftix(ix,iy),leftiy(ix,iy)))
          else
            gfunr(ix,iy,0) = 0.0_R8
          endif
        enddo
      enddo

*   ..compute gradients at y-faces
      do iy = -1,ny
        do ix = -1,nx
          if (bottomiy(ix,iy).ne.-2) then
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
cwdk temporary: compute gradients between guard cells for backwards
cwdk            compatibility (avoid problems in corner cells)
cwdk later on : these faces will be removed altogether
              gfunr(ix,iy,1) = 2.0_R8*
     &          (fun(ix,iy) - fun(bottomix(ix,iy),bottomiy(ix,iy)))/    ! contribution from y-direction
     &          (hy1(ix,iy) + hy1(bottomix(ix,iy),bottomiy(ix,iy)))
            else
              gfunr(ix,iy,1) = 2.0_R8*
     &          (fun(ix,iy) - fun(bottomix(ix,iy),bottomiy(ix,iy)) -    ! contribution from y-direction
     &          (func(rightix(ix,iy),rightiy(ix,iy)) - func(ix,iy))*    ! contribution from x-direction (non-orthogonal grids)
     &          (hy(ix,iy) + hy(bottomix(ix,iy),bottomiy(ix,iy)))/
     &          (hx(ix,iy) + hx(bottomix(ix,iy),bottomiy(ix,iy)))*
     &          qs(ix,iy,1))/
     &          (hy1(ix,iy) + hy1(bottomix(ix,iy),bottomiy(ix,iy)))
            endif
          else
            gfunr(ix,iy,1) = 0.0_R8
          endif
        enddo
      enddo

*   ..reset null regions
      call b2xxgs (nx, ny, 0.0_R8, gfunr(-1,-1,0), 1)
      call b2xxgs (nx, ny, 0.0_R8, gfunr(-1,-1,1), 2)

*   ..return
      return

*-----------------------------------------------------------------------
*.end grad_r

      end subroutine grad_r
