c
      subroutine adjust_axes(onflag)
      use b2mod_types
#include "graphics_vars.h"
      logical onflag
      if(onflag) then
        nice_axes=.not.nice_axes
        write (*,*) 'nice axes state changed'
      else
        nice_axes=.true.
        write (*,*) 'graphs plotted with nice axes'
      endif
      return
      end subroutine adjust_axes
