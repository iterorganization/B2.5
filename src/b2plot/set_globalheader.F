
      subroutine set_globalheader
      use b2mod_types
      implicit none
#include <graphics_vars.h>

      globalheader=extralabel
      extralabel=' '
      return
      end subroutine set_globalheader
