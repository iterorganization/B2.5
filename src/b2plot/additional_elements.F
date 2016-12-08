c
      subroutine additional_elements(nx,ny,physical)
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      integer nx, ny
      logical physical

      if(grid) then
        call meshit(physical)
      endif
      if(sptrix) then
        call pltsptrix(physical)
      endif
      if(vessel.and.physical) then
        call pltvessel
      endif
      if(outline) then
        call pltoutline(physical)
      endif
#ifdef WOS
      call flux_add(nx,ny)
#endif
      if(wall) then
        call pltwall(nx,ny,physical,.false.)
      endif
      if(target.and.physical) then
        call plttarget
      endif
      if(trig.and.physical) then
        call plttriang
      endif
      return
      end subroutine additional_elements
