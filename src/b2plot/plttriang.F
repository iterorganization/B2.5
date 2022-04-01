c
      subroutine plttriang()
      use b2mod_types
      use b2mod_b2plot
#ifdef B25_EIRENE
      use eirmod_ctrig
      implicit none
      integer it

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
#endif
      return
      end subroutine plttriang
