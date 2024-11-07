c
      subroutine nofill(filflg)
      use b2mod_types
      use b2mod_graphics
      implicit none
      logical filflg

      colprt=.not.filflg
      if(colprt) then
        write (*,*) 'Filled Contour Levels'
      else
        write (*,*) 'Empty Contour Levels'
      endif
      return
      end subroutine nofill
