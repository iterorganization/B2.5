c
      subroutine trigon(onflag)
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      logical onflag

      if(onflag) then
        trig=.not.trig
        write (*,*) 'triangular grid state changed'
      else
        trig=.false.
        write (*,*) 'triangular grid off'
      endif
      return
      end subroutine trigon
