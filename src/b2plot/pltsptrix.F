c
      subroutine pltsptrix(physical)
      use b2mod_types
      use b2mod_b2plot
      use b2us_data
      implicit none
      logical physical
      integer i,iCv,iFc,iFs,iVx1,iVx2,iXpt

#ifndef PLPLOT
      call lastd
      if(scol .gt. 1) call gsplci(scol)
      call gslwsc_d(sepwidth)
#endif
      if(physical) then
        do iXpt = 1, mpg%nXpt
          iFs = mpg%vxFs(mpg%Xpt(iXpt))
          do i = 1, mpg%fsFcP(iFs,2)
            iFc = mpg%fsFc(mpg%fsFcP(iFs,1)+i-1)
            iVx1 = mpg%fcVx(iFc,1)
            iVx2 = mpg%fcVx(iFc,2)
            call lined_d(geo%vxX(iVx1), geo%vxY(iVx1),
     &                   geo%vxX(iVx2), geo%vxY(iVx2))
          end do
        end do
        if (mpg%nXpt.eq.0.and.mpg%iFssep.gt.0) then
          do i = 1, mpg%fsFcP(mpg%iFssep,2)
            iFc = mpg%fsFc(mpg%fsFcP(mpg%iFssep,1)+i-1)
            iVx1 = mpg%fcVx(iFc,1)
            iVx2 = mpg%fcVx(iFc,2)
            call lined_d(geo%vxX(iVx1), geo%vxY(iVx1),
     &                   geo%vxX(iVx2), geo%vxY(iVx2))
          end do
        end if
      else if(mpg%isClassicalGrid.eq.1) then
        do iXpt = 1, mpg%nXpt
          iFs = mpg%vxFs(mpg%Xpt(iXpt))
          do i = 1, mpg%fsFcP(iFs,2)
            iFc = mpg%fsFc(mpg%fsFcP(iFs,1)+i-1)
            iVx1 = mpg%fcVx(iFc,1)
            iVx2 = mpg%fcVx(iFc,2)
            call lined_d(real(vxIx(iVx1),R8)-0.5_R8,
     &                   real(vxIy(iVx1),R8)-0.5_R8,
     &                   real(vxIx(iVx2),R8)+0.5_R8,
     &                   real(vxIx(iVx2),R8)-0.5_R8)
          end do
        end do
        if (mpg%nXpt.eq.0.and.mpg%iFssep.gt.0) then
          do i = 1, mpg%fsFcP(mpg%iFssep,2)
            iFc = mpg%fsFc(mpg%fsFcP(mpg%iFssep,1)+i-1)
            iVx1 = mpg%fcVx(iFc,1)
            iVx2 = mpg%fcVx(iFc,2)
            call lined_d(real(vxIx(iVx1),R8)-0.5_R8,
     &                   real(vxIy(iVx1),R8)-0.5_R8,
     &                   real(vxIx(iVx2),R8)+0.5_R8,
     &                   real(vxIx(iVx2),R8)-0.5_R8)
          end do
        end if
#ifdef WG_TODO
        do ix = -1, nx
          do iy = -1, ny
            if(rightix(ix,iy).ne.ix+1) then
              x1=ix+0.5_R8
              y1=iy-0.5_R8
              x2=ix+0.5_R8
              y2=iy+0.5_R8
              call lined_d(x1,y1,x2,y2)
            end if
            if(leftix(ix,iy).ne.ix-1) then
              x1=ix-0.5_R8
              y1=iy-0.5_R8
              x2=ix-0.5_R8
              y2=iy+0.5_R8
              call lined_d(x1,y1,x2,y2)
            end if
            if(topiy(ix,iy).ne.iy+1) then
              x1=ix-0.5_R8
              y1=iy+0.5_R8
              x2=ix+0.5_R8
              y2=iy+0.5_R8
              call lined_d(x1,y1,x2,y2)
            end if
            if(bottomiy(ix,iy).ne.iy-1) then
              x1=ix-0.5_R8
              y1=iy-0.5_R8
              x2=ix+0.5_R8
              y2=iy-0.5_R8
              call lined_d(x1,y1,x2,y2)
            end if
          end do
        end do
#endif
      end if
#ifndef PLPLOT
      call gslwsc_d(lwidth)
      call lastd
      if(scol .gt. 1) call gsplci(1)
#endif
      return
      end subroutine pltsptrix
