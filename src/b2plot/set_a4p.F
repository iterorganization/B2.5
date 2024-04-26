
      subroutine set_a4p
      use b2mod_types
#include <graphics_vars.h>
      xwidth=1.0/sqrt(2.0)
      yheight=1.0
#ifdef PLPLOT
      call plsori(1)
#endif
      return
      end subroutine set_a4p
