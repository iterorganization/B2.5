
      subroutine terminate_plotting
      use b2mod_types
#include <graphics_vars.h>
#ifdef NCAR4      
      integer iw
#endif
      if(.not.initialised_plotting) return
      if(same_page.ne.0) then
        same_page=0
        call end_page
      endif
#ifdef PLPLOT
      call plend()
#else
      if(dirty_page) call frame
#ifdef NCAR4      
      do iw=1,nwork
        call gdawk(iw)
        call gclwk(iw)
      enddo
      call gclks
#else
      call clsgks
#endif
#endif
      return
      end subroutine terminate_plotting
