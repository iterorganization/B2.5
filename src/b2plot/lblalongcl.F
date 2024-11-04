c
      subroutine lblalongcl(along)
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      logical along

      if(along) then
        lalong=.true.
        write (*,*) 'Labels: along Contour Lines'
      else
        lalong=.false.
        write (*,*) 'Labels: normal'
      endif
      return
      end subroutine lblalongcl
