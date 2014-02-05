c
      subroutine linewidth(width)
      use b2mod_types
#include "graphics_vars.h"
      real (kind=R8) ::
     * width
      lwidth=width		
      write (*,*) 'Line Width set to',width
      return
      end subroutine linewidth
c
      subroutine walllinewidth(width)
      use b2mod_types
#include "graphics_vars.h"
      real (kind=R8) ::
     * width
      wwidth=width		
      write (*,*) 'Wall Line Width set to',width
      return
      end subroutine walllinewidth
