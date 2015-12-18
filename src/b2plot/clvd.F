c
      subroutine clvd(clvmin,clvmax,delta)
      use b2mod_types
#include "graphics_vars.h"
      real (kind=R8) :: clvmin,clvmax,delta
      external subini, subend, xertst, xerrab

      call subini ('clvd')
      call xertst(delta.ne.0.0_R8, 'division by 0 in clvd')
      nofcls=(clvmax-clvmin)/delta
      if(nofcls.eq.0) then
        call xerrab ('Contour Level Max == Contour Level Min')
      else
        nofcls=nofcls+1
        write (*,*) 'User specified',nofcls,'contour levels for range'
        write (*,*) 'Using customised contour level range:'
        do iclv=1,nofcls
          clval(iclv)=clvmin+(iclv-1)*delta
          write (*,*) 'contour level',iclv,'at',clval(iclv)
        enddo
      endif
      call subend ()
      return
      end subroutine clvd
c
      subroutine wclvd(clvmin,clvmax,delta)
      use b2mod_types
#include "graphics_vars.h"
      real (kind=R8) :: clvmin,clvmax,delta
      external subini, subend, xertst, xerrab

      call subini ('wclvd')
      call xertst (delta.ne.0.0_R8, 'division by 0 in wclvd')
      nofwcls=(clvmax-clvmin)/delta
      if(nofwcls.eq.0) then
        call xerrab ('Wall Contour Level Max == Wall Contour Level Min')
      else
        nofwcls=nofwcls+1
        write (*,*) 'User specified',nofwcls,
     .   'wall contour levels for range'
        write (*,*) 'Using customised wall contour level range:'
        do iclv=1,nofwcls
          wclval(iclv)=clvmin+(iclv-1)*delta
          write (*,*) 'wall contour level',iclv,'at',wclval(iclv)
        enddo
      endif
      call subend ()
      return
      end subroutine wclvd
c
      subroutine tclvd(clvmin,clvmax,delta)
      use b2mod_types
#include "graphics_vars.h"
      real (kind=R8) :: clvmin,clvmax,delta
      external subini, subend, xertst, xerrab

      call subini ('tcvld')
      call xertst(delta.ne.0.0_R8, 'division by 0 in tclvd')
      noftcls=(clvmax-clvmin)/delta
      if(noftcls.eq.0) then
        call xerrab 
     .   ('Target Contour Level Max == Target Contour Level Min')
      else
        noftcls=noftcls+1
        write (*,*) 'User specified',noftcls,
     .   'target contour levels for range'
        write (*,*) 'Using customised target contour level range:'
        do iclv=1,noftcls
          tclval(iclv)=clvmin+(iclv-1)*delta
          write (*,*) 'target contour level',iclv,'at',tclval(iclv)
        enddo
      endif
      call subend ()
      return
      end subroutine tclvd
 
