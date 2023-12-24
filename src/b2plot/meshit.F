c
      subroutine meshit(physical)
      use b2mod_types
      use b2mod_b2plot
      use b2us_data
      implicit none
      logical physical
      real (kind=R8) ::
     *   x1,x2,y1,y2
      integer iFc,iVx1,iVx2

      if(mcol .gt. 1) then
#ifndef PLPLOT
        call lastd
        call gsplci(mcol)
#endif
      endif
      if(.not.physical) then ! mesh computational domain
#ifdef WG_TODO
        do ix=-1,nx+1
          x1=ix-0.5_R8
          x2=ix-0.5_R8
          y1=-1-0.5_R8
          y2=ny+1-0.5_R8
          call lined_d(x1,y1,x2,y2)
        enddo
        do iy=-1,ny+1
          x1=-1-0.5_R8
          x2=nx+1-0.5_R8
          y1=iy-0.5_R8
          y2=iy-0.5_R8
          call lined_d(x1,y1,x2,y2)
        enddo
#endif
      endif
      if(physical) then ! mesh physical domain
        do iFc=1,mpg%nFc
          iVx1 = mpg%fcVx(iFc,1)
          iVx2 = mpg%fcVx(iFc,2)
          call lined_d(geo%vxX(iVx1),geo%vxY(iVx1),
     1                 geo%vxX(iVx2),geo%vxY(iVx2))
        enddo
      endif
      if(mcol .gt. 1) then
#ifndef PLPLOT
        call lastd
        call gsplci(1)
#endif
      endif
      return
      end subroutine meshit
