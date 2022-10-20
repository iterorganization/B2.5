
      subroutine clear_extralabel(label)
      use b2mod_types
#include <graphics_vars.h>
      character*(*) label
      extralabel=label
      return
      end subroutine clear_extralabel
