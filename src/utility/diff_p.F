*-----------------------------------------------------------------------
*.specification

      subroutine diff_p(nx, ny, mode, hx, hy, vol, qs, fun, func, dfunp)
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
     *   dfunp(-1:nx,-1:ny,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*     DIFF_P computes differences of a cell centered quantity fun in the 
*     poloidal direction, on x and y-faces of a cell. As long as 
*     the x-direction is aligned with the poloidal direction, these 
*     differences are fairly trivial.
*     dfunp(:,:,0) will hold the poloidal differences at x-faces
*     dfunp(:,:,1) will hold the poloidal differences at y-faces
*
*     Note: dfunp will not be computed between guard cells, because 
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

*   ..compute poloidal difference at x-faces
      do iy = -1,ny
        do ix = -1,nx
          if(leftix(ix,iy).ne.-2 .and.
     &       bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
            dfunp(ix,iy,0) = 
     &       fun(ix,iy) - fun(leftix(ix,iy),leftiy(ix,iy))
          else
            dfunp(ix,iy,0) = 0.0_R8
          endif
        enddo
      enddo

*   ..compute difference at y-faces
      do iy = -1,ny
        do ix = -1,nx
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and. 
     &       bottomiy(ix,iy).ne.-2) then
            dfunp(ix,iy,1) = 
     &       func(rightix(ix,iy),rightiy(ix,iy)) - func(ix,iy)
          else
            dfunp(ix,iy,1) = 0.0_R8
          endif
        enddo
      enddo

*   ..reset null regions
      call b2xxgs (nx, ny, 0.0_R8, dfunp(-1,-1,0), 1)
      call b2xxgs (nx, ny, 0.0_R8, dfunp(-1,-1,1), 2)

*   ..return
      return

*-----------------------------------------------------------------------
*.end diff_p

      end subroutine diff_p
