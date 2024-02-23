c
      subroutine trigmesh
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2mod_boundary_namelist
#ifdef B25_EIRENE
      use eirmod_ctrig
      implicit none
      integer it
      character*500 plotlabel

      if(b2plot_debug.ge.7) write(*,*) 'physical domain'
      plotlabel='Triangular mesh'
      call init_page(trim(plotlabel),
     1   values(5),values(6),values(7),values(8),1,.true.,.false.,
     2   .true.)
      if(mcol .gt. 1) then
#ifndef PLPLOT
        call lastd
        call gsplci(mcol)
#endif
      endif
      do it=1,ntrii
        call lined_d(xtrian(necke(1,it)),ytrian(necke(1,it)),
     1               xtrian(necke(2,it)),ytrian(necke(2,it)))
        call lined_d(xtrian(necke(2,it)),ytrian(necke(2,it)),
     1               xtrian(necke(3,it)),ytrian(necke(3,it)))
        call lined_d(xtrian(necke(3,it)),ytrian(necke(3,it)),
     1               xtrian(necke(1,it)),ytrian(necke(1,it)))
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
#endif
      return
      end subroutine trigmesh
