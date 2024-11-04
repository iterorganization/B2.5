c
      subroutine clvl(onflag)
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      logical onflag
      if(onflag) then
        contours=.true.
        write (*,*) 'contour level lines on'
      else
        contours=.false.
        write (*,*) 'contour level lines off'
      endif
      return
      end subroutine clvl
