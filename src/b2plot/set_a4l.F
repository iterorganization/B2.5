
      subroutine set_a4l
      use b2mod_types
#include <graphics_vars.h>
      xwidth=1.0
      yheight=1.0/sqrt(2.0)
#ifdef PLPLOT
      call plsori(0)
#endif
      return
      end
