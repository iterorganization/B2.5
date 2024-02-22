c
      subroutine mesh
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2mod_boundary_namelist
      implicit none
      real (kind=R8) ::
     *   x1,x2,y1,y2
      integer ix,iy
      character*500 plotlabel

      if(b2plot_debug.ge.7) write(*,*) 'typeofplot = ',typeofplot

      if(typeofplot .ge. 2) then
        if(b2plot_debug.ge.7) write(*,*) 'computational domain'
        plotlabel='Computational mesh'
        call init_page(trim(plotlabel),
     1   values(1),values(2),values(3),values(4),1,.false.,.false.,
     2   .false.)
        if(mcol .gt. 1) then
#ifndef PLPLOT
          call lastd
          call gsplci(mcol)
#endif
        endif
        do ix=-1,nx+1
          x1=ix-0.5_R8
          x2=ix-0.5_R8
          y1=-1-0.5_R8
          y2=ny+0.5_R8
          call lined_d(x1,y1,x2,y2)
        enddo
        do iy=-1,ny+1
          x1=-1-0.5_R8
          x2=nx+0.5_R8
          y1=iy-0.5_R8
          y2=iy-0.5_R8
          call lined_d(x1,y1,x2,y2)
        enddo
        if(mcol .gt. 1) then
#ifndef PLPLOT
          call lastd
          call gsplci(1)
#endif
        endif
        call data_copy(0)
        call additional_elements(nx,ny,.false.)
        call end_page
      endif
      if(mod(typeofplot,2) .eq. 1) then
        if(b2plot_debug.ge.7) write(*,*) 'physical domain'
        plotlabel='Physical mesh'
        call init_page(trim(plotlabel),
     1   values(5),values(6),values(7),values(8),1,.true.,.false.,
     2   .true.)
        if(mcol .gt. 1) then
#ifndef PLPLOT
          call lastd
          call gsplci(mcol)
#endif
        endif
        do iy=-1,ny
          do ix=-1,nx
            if(region(ix,iy,0).ne.0) then
              call lined_d(crx(ix,iy,0),cry(ix,iy,0),
     1         crx(ix,iy,1),cry(ix,iy,1))
              call lined_d(crx(ix,iy,0),cry(ix,iy,0),
     1         crx(ix,iy,2),cry(ix,iy,2))
              call lined_d(crx(ix,iy,1),cry(ix,iy,1),
     2         crx(ix,iy,3),cry(ix,iy,3))
              call lined_d(crx(ix,iy,2),cry(ix,iy,2),
     2         crx(ix,iy,3),cry(ix,iy,3))
            endif
          enddo
        enddo
        if(mcol .gt. 1) then
#ifndef PLPLOT
          call lastd
          call gsplci(1)
#endif
        endif
        call data_copy(0)
        call pltvessel
        call additional_elements(nx,ny,.true.)
        call end_page
      endif
      return
      end subroutine mesh
