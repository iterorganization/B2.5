c
      subroutine mesh
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2mod_boundary_namelist
      implicit none
      integer lnblnk,iFc
      character*500 plotlabel
      external lnblnk

      if(b2plot_debug.ge.7) write(*,*) 'typeofplot = ',typeofplot

      if(typeofplot .ge. 2 .and. mpg%isClassicalGrid .eq. 1) then
        if(b2plot_debug.ge.7) write(*,*) 'computational domain'
        plotlabel='Computational mesh'
        call init_page(plotlabel(1:lnblnk(plotlabel)),
     1   values(1),values(2),values(3),values(4),1,.false.,.false.,
     2   .false.)
        if(mcol .gt. 1) then
#ifndef PLPLOT
          call lastd
          call gsplci(mcol)
#endif
        endif
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
        if(mcol .gt. 1) then
#ifndef PLPLOT
          call lastd
          call gsplci(1)
#endif
        endif
        call data_copy(0)
        call additional_elements(mpg%nCv,.false.)
        call end_page
      endif
      if(mod(typeofplot,2) .eq. 1) then
        if(b2plot_debug.ge.7) write(*,*) 'physical domain'
        plotlabel='Physical mesh'
        call init_page(plotlabel(1:lnblnk(plotlabel)),
     1   values(5),values(6),values(7),values(8),1,.true.,.false.,
     2   .true.)
        if(mcol .gt. 1) then
#ifndef PLPLOT
          call lastd
          call gsplci(mcol)
#endif
        endif
        do iFc=1,mpg%nFc
          call lined_d(geo%vxX(mpg%fcVx(iFc,1)),
     .                 geo%vxY(mpg%fcVx(iFc,1)),
     .                 geo%vxX(mpg%fcVx(iFc,2)),
     .                 geo%vxY(mpg%fcVx(iFc,2)))
        enddo
        if(mcol .gt. 1) then
#ifndef PLPLOT
          call lastd
          call gsplci(1)
#endif
        endif
        call data_copy(0)
        call pltvessel
        call additional_elements(mpg%nCv,.true.)
        call end_page
      endif
      return
      end subroutine mesh
