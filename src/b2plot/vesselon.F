c
      subroutine vesselon(onflag)
      use b2mod_types
#include "graphics_vars.h"
      logical onflag
      if(onflag) then
        vessel=.not.vessel
        write (*,*) 'vessel state changed'
      else
        vessel=.false.
        write (*,*) 'vessel off'
      endif
      return
      end subroutine vesselon
