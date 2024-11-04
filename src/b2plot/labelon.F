c
      subroutine labelon(onflag)
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      logical onflag

      labels=onflag
      if(labels) then
        lalong=.true.
        write (*,*) 'Labeling on'
      else
        lalong=.false.
        write (*,*) 'Labeling off'
      endif
      return
      end subroutine labelon
