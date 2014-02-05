      subroutine wallon(onflag,inum,wmind,wmaxd,ilog)
c
c  version : 23.08.97 11:56
c
c======================================================================
      use b2mod_types
#include "graphics_vars.h"
      logical onflag
      integer inum, ilog
      real (kind=R8) ::
     * wmind,wmaxd
c======================================================================
c	 print *,'%% wallon: onflag,inum,wmind,wmaxd,wlinlog=',
c     , 	             onflag,inum,wmind,wmaxd,wlinlog
      wmin=wmind
      wmax=wmaxd
      wlinlog=ilog
      if(onflag) then
        wall=inum.ne.0
      else
        wall=.false.
      endif
      if(wall) call ank_interface(inum)
c	 if(onflag) call ank_interface(inum)
      return
c======================================================================
      entry set_wall_on(wmind,wmaxd,ilog)
      wmin=wmind
      wmax=wmaxd
      wlinlog=ilog
      wall=.true.
      return
c======================================================================
      entry set_wlinlog(ilog)
      wlinlog=ilog
      return
c======================================================================
      entry walloff
      wall=.false.
      return
c======================================================================
      end subroutine wallon
