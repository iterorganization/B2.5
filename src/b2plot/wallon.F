      subroutine wallon(onflag,inum,wmind,wmaxd,ilog)
c
c  version : 23.08.97 11:56
c
c======================================================================
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      logical, intent(in) :: onflag
      integer, intent(in) :: inum, ilog
      real (kind=R8), intent(in) :: wmind, wmaxd
c======================================================================
c	 print *,'%% wallon: onflag,inum,wmind,wmaxd,wlinlog=',
c     , 	             onflag,inum,wmind,wmaxd,wlinlog
      wmin=real(wmind,R4)
      wmax=real(wmaxd,R4)
      wlinlog=ilog
      if(onflag) then
        wall=inum.ne.0
      else
        wall=.false.
      endif
      if(wall) call ank_interface(inum)
c	 if(onflag) call ank_interface(inum)
      return
      end subroutine wallon
c======================================================================
      subroutine set_wall_on(wmind,wmaxd,ilog)
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      integer, intent(in) :: ilog
      real (kind=R8), intent(in) :: wmind, wmaxd
      wmin=real(wmind,R4)
      wmax=real(wmaxd,R4)
      wlinlog=ilog
      wall=.true.
      return
      end subroutine set_wall_on
c======================================================================
      subroutine set_wlinlog(ilog)
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      integer, intent(in) :: ilog
      wlinlog=ilog
      return
      end subroutine set_wlinlog
c======================================================================
      subroutine walloff
      use b2mod_types
      implicit none
#include "graphics_vars.h"
      wall=.false.
      return
      end subroutine walloff
c======================================================================
