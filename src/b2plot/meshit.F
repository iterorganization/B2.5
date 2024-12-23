c
      subroutine meshit(physical)
      use b2mod_types
      use b2mod_b2plot
      use b2us_data
      implicit none
      logical physical
      integer iFc,iVx1,iVx2

      if(mcol .gt. 1) then
#ifndef PLPLOT
        call lastd
        call gsplci(mcol)
#endif
      endif
      if(.not.physical.and.mpg%isClassicalGrid.eq.1) then ! mesh computational domain
        do iFc=1,mpg%nFc
          if(fcxIx(iFc).ne.-2.and.fcxIy(iFc).ne.-2) then
            call lined_d(real(fcxIx(iFc),R8)-0.5_R8,
     &                   real(fcxIy(iFc),R8)-0.5_R8,
     &                   real(fcxIx(iFc),R8)-0.5_R8,
     &                   real(fcxIy(iFc),R8)+0.5_R8)
          endif
          if(fcyIx(iFc).ne.-2.and.fcyIy(iFc).ne.-2) then
            call lined_d(real(fcyIx(iFc),R8)-0.5_R8,
     &                   real(fcyIy(iFc),R8)-0.5_R8,
     &                   real(fcyIx(iFc),R8)+0.5_R8,
     &                   real(fcyIy(iFc),R8)-0.5_R8)
          endif
        enddo
      else if(physical) then ! mesh physical domain
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
