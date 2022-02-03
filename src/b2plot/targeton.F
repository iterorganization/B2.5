      subroutine targeton(onflag,trgmind,trgmaxd,trglog)
c
c  version : 23.08.97 11:56
c
c======================================================================
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      logical, intent(in) :: onflag
      integer, intent(in) :: trglog
      real (kind=R8), intent(in) :: trgmind,trgmaxd
c======================================================================
c     print *,'%% targeton: onflag,trgmind,trgmaxd,trglog=',
c    ,                      onflag,trgmind,trgmaxd,trglog
      trgmin=real(trgmind,R4)
      trgmax=real(trgmaxd,R4)
      trglinlog=trglog
      target=onflag
      return
      end subroutine targeton
c======================================================================
      subroutine set_trglinlog(trglog)
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      integer, intent(in) :: trglog
      trglinlog=trglog
      return
      end subroutine set_trglinlog
c======================================================================
      subroutine targetoff
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      target=.false.
      return
      end subroutine targetoff
c======================================================================
