c
      subroutine outlineon(onflag)
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      logical onflag

      if(onflag) then
        outline=.not.outline
        write (*,*) 'outline state changed'
      else
        outline=.false.
        write (*,*) 'outline off'
      endif
      return
      end subroutine outlineon
