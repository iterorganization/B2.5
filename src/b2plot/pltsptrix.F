c
      subroutine pltsptrix(physical)
      use b2mod_types
      use b2mod_b2plot
      use b2us_data
      implicit none
      logical physical
      integer i,iCv1,iCv2,iFc,iFs,iVx1,iVx2,iXpt

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
     &                   real(vxIx(iVx2),R8)-0.5_R8,
     &                   real(vxIy(iVx2),R8)-0.5_R8)
          end do
        end do
        if (mpg%nXpt.eq.0.and.mpg%iFssep.gt.0) then
          do i = 1, mpg%fsFcP(mpg%iFssep,2)
            iFc = mpg%fsFc(mpg%fsFcP(mpg%iFssep,1)+i-1)
            iVx1 = mpg%fcVx(iFc,1)
            iVx2 = mpg%fcVx(iFc,2)
            call lined_d(real(vxIx(iVx1),R8)-0.5_R8,
     &                   real(vxIy(iVx1),R8)-0.5_R8,
     &                   real(vxIx(iVx2),R8)-0.5_R8,
     &                   real(vxIy(iVx2),R8)-0.5_R8)
          end do
        end if
        do iFc = 1, mpg%nFc
          iCv1 = mpg%fcCv(iFc,1)
          iCv2 = mpg%fcCv(iFc,2)
          if (fcxIx(iFc).ne.-2.and.fcxIy(iFc).ne.-2) then
            if (cvIx(iCv1).eq.-2.or.cvIy(iCv1).eq.-2.or.
     &          cvIx(iCv2).eq.-2.or.cvIy(iCv2).eq.-2.or.
     &          cvIy(iCv1).ne.cvIy(iCv2).or.(.not.
     &         (cvIx(iCv1).eq.cvIx(iCv2)+1.or.
     &          cvIx(iCv1).eq.cvIx(iCv2)-1))) then
              call lined_d(real(fcxIx(iFc),R8)-0.5_R8,
     &                     real(fcxIy(iFc),R8)-0.5_R8,
     &                     real(fcxIx(iFc),R8)-0.5_R8,
     &                     real(fcxIy(iFc),R8)+0.5_R8)

            end if
          end if
          if (fcyIx(iFc).ne.-2.and.fcyIy(iFc).ne.-2) then
            if (cvIx(iCv1).eq.-2.or.cvIy(iCv1).eq.-2.or.
     &          cvIx(iCv2).eq.-2.or.cvIy(iCv2).eq.-2.or.
     &          cvIx(iCv1).ne.cvIx(iCv2).or.(.not.
     &         (cvIy(iCv1).eq.cvIy(iCv2)+1.or.
     &          cvIy(iCv1).eq.cvIy(iCv2)-1))) then
              call lined_d(real(fcyIx(iFc),R8)-0.5_R8,
     &                     real(fcyIy(iFc),R8)-0.5_R8,
     &                     real(fcyIx(iFc),R8)+0.5_R8,
     &                     real(fcyIy(iFc),R8)-0.5_R8)

            end if
          end if
        end do
      end if
#ifndef PLPLOT
      call gslwsc_d(lwidth)
      call lastd
      if(scol .gt. 1) call gsplci(1)
#endif
      return
      end subroutine pltsptrix
