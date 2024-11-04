c
      subroutine ipolate(onflag)
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      logical onflag

      if(onflag) then
        intpol=.true.
      else
        intpol=.false.
      endif
      return
      end subroutine ipolate
