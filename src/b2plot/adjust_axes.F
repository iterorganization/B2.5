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
c
      subroutine adjust_bold(onflag)
      use b2mod_types
#include "graphics_vars.h"
      logical onflag
      if(onflag) then
        bold=.not.bold
        write (*,*) 'bold state changed'
      else
        bold=.false.
        write (*,*) 'graphs plotted with normal labels'
      endif
      return
      end subroutine adjust_bold
