c
      subroutine additional_elements(nCv,physical)
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      integer nCv
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
      call flux_add(nCv)
#endif
      if(wall) then
        call pltwall(nCv,physical,.false.)
      endif
      if(target.and.physical) then
        call plttarget
      endif
      if(trig.and.physical) then
        call plttriang
      endif
      return
      end subroutine additional_elements
