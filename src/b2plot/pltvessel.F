c
      subroutine pltvessel
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_debug
      implicit none
      logical exist
      real (kind=R8) ::
     * x1,y1,x2,y2,r1o,z1o,r2o,z2o,r1,z1,r2,z2,dist
      integer icnt

      dist(x1,y1,x2,y2)=sqrt((x1-x2)**2+(y1-y2)**2)

#ifndef PLPLOT
      call lastd
      if(vcol .gt. 1) call gsplci(vcol)
      call gslwsc_d(vwidth)
#endif
      icnt=0
      if(mod(typeofplot,2) .eq. 1) then ! mesh physical domain
        inquire(file='../mesh.extra',exist=exist)
        if(exist) then
          open(99,file='../mesh.extra')
          read(99,*,end=20) r1o,z1o,r2o,z2o
   10     continue
          read(99,*,end=20) r1,z1,r2,z2
          if(r2o .eq. r1 .and. z2o .eq. z1 .and. 
     1	   dist(r1o,z1o,r2,z2) .lt. 1.0e-2_R8) then
            r2o=r2
            z2o=z2
            icnt=icnt+1
          else
            call lined_d(r1o,z1o,r2o,z2o)
            r1o=r1
            z1o=z1
            r2o=r2
            z2o=z2
          endif
          goto 10
   20     continue
          call lined_d(r1o,z1o,r2,z2)
          close(99)
          if((icnt.gt.0.and.b2plot_debug.ge.9).or.b2plot_debug.ge.10)
     1     write(*,*) 'Saved plotting ',icnt,' line segments'
        endif
      endif
#ifndef PLPLOT
      call lastd
      if(vcol .gt. 1) call gsplci(1)
      call gslwsc_d(lwidth)
#endif
      return
      end subroutine pltvessel
