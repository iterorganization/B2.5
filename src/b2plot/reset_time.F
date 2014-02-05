c
      subroutine reset_time(time)
      use b2mod_types
      implicit none
      real (kind=R8) :: time
#include "graphics_vars.h"
      runtime=time
      return
      end subroutine reset_time
