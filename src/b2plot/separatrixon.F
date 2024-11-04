c
      subroutine separatrixon(onflag)
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      logical onflag

      if(onflag) then
        sptrix=.not.sptrix
        write (*,*) 'separatrix state changed'
      else
        sptrix=.false.
        write (*,*) 'separatrix off'
      endif
      return
      end subroutine separatrixon
