      subroutine targeton(onflag,trgmind,trgmaxd,trglog)
c
c  version : 23.08.97 11:56
c
c======================================================================
      use b2mod_types
#include "graphics_vars.h"
      logical onflag
      integer trglog
      real (kind=R8) ::
     * trgmind,trgmaxd
c======================================================================
c	 print *,'%% targeton: onflag,trgmind,trgmaxd,trglog=',
c     , 	               onflag,trgmind,trgmaxd,trglog
      trgmin=trgmind
      trgmax=trgmaxd
      trglinlog=trglog
      target=onflag
      return
c======================================================================
      entry set_trglinlog(trglog)
      trglinlog=trglog
      return
c======================================================================
      entry targetoff
      target=.false.
      return
c======================================================================
      end subroutine targeton
