c
      subroutine gridon(onflag)
      use b2mod_types
#include "graphics_vars.h"
      logical onflag
      if(onflag) then
        grid=.not.grid
        write (*,*) 'grid state changed'
      else
        grid=.false.
        write (*,*) 'grid off'
      endif
      return
      end subroutine gridon
