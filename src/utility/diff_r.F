*-----------------------------------------------------------------------
*.specification

      subroutine diff_r(nx, ny, mode, hx, hy, vol, qs, fun, func, dfunr)
      use b2mod_types
      use b2mod_indirect
      implicit none
*   ..input arguments
      integer nx,ny,mode
      real (kind=R8) ::
     *   hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), vol(-1:nx,-1:ny),
     *   qs(-1:nx,-1:ny,0:1), fun(-1:nx,-1:ny), func(-1:nx,-1:ny)
*   ..output arguments
      real (kind=R8) ::
     *   dfunr(-1:nx,-1:ny,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*     DIFF_R computes differences of a cell centered quantity fun in the 
*     radial direction, on x and y-faces of a cell. It takes into 
*     account the corrections for x-faces not perpendicular to the 
*     poloidal field (9-point stencil).
*     dfunr(:,:,0) will hold the radial differences at x-faces
*     dfunr(:,:,1) will hold the radial differences at y-faces
*
*     To compute values of the field in cell corners, the routine relies 
*     calls to intcor, using the input array 'vol' as weights. Different
*     interpolation schemes can be mimicked by using a different array 
*     with weights in the calling routine.
*
*     Note: dfunr will not be computed between guard cells, because 
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
     .  'diff_p--faulty argument mode')

*   ..interpolate fun to cell corners, using volume weighted averaging
      if (mode.eq.0) then
        call intcor(nx,ny,vol,fun,func)
      endif

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
     &       qs(ix,iy,0)
          else
            dfunr(ix,iy,0) = 0.0_R8
          endif
        enddo
      enddo

*   ..compute difference at y-faces
      do iy = -1,ny
        do ix = -1,nx
          if (bottomiy(ix,iy).ne.-2) then
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
cwdk temporary: compute differences between guard cells for backwards
cwdk            compatibility (avoid problems in corner cells)
cwdk later on : these faces will be removed altogether
              dfunr(ix,iy,1) = 
     &          fun(ix,iy) - fun(bottomix(ix,iy),bottomiy(ix,iy))       ! contribution from y-direction
            else
              dfunr(ix,iy,1) = 
     &          fun(ix,iy) - fun(bottomix(ix,iy),bottomiy(ix,iy)) -     ! contribution from y-direction
     &          (func(rightix(ix,iy),rightiy(ix,iy)) - func(ix,iy))*    ! contribution from x-direction (non-orthogonal grids)
     &          (hy(ix,iy) + hy(bottomix(ix,iy),bottomiy(ix,iy)))/
     &          (hx(ix,iy) + hx(bottomix(ix,iy),bottomiy(ix,iy)))*
     &          qs(ix,iy,1)
            endif
          else
            dfunr(ix,iy,1) = 0.0_R8
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
