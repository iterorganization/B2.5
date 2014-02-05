c
      subroutine nofill(filflg)
      use b2mod_types
#include "graphics_vars.h"
      logical filflg
      colprt=.not.filflg
      if(colprt) then
        write (*,*) 'Filled Contour Levels'
      else
        write (*,*) 'Empty Contour Levels'
      endif
      return
      end subroutine nofill
