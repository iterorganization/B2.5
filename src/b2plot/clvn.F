c
      subroutine clvn(clvmin,clvmax,ncls)
      use b2mod_types
      use b2mod_subsys
      implicit none
#include "graphics_vars.h"
      real (kind=R8) ::
     * clvmin,clvmax,delta
      integer :: iclv, ncls
      external xerrab

      call subini ('clvn')
      nofcls=ncls
      if(nofcls.eq.0) then
        call xerrab ('You set number of contour levels to 0')
      elseif(nofcls.eq.1) then
        clval(1)=real((clvmax+clvmin)*0.5_R8,R4)
        write (*,*) 'You specified 1 contour level at',clval(1)
      else
        delta=(clvmax-clvmin)/(nofcls-1)
        write (*,*) 'User specified',nofcls,'contour levels'
        write (*,*) 'Using customised contour levels:'
        do iclv=1,nofcls
          clval(iclv)=real(clvmin+(iclv-1)*delta,R4)
          write (*,*) 'contour level',iclv,'at',clval(iclv)
        enddo
      endif
      call subend ()
      return
      end subroutine clvn
c
      subroutine eclvn(clvmin,clvmax,ncls)
      use b2mod_types
      use b2mod_subsys
      implicit none
#include "graphics_vars.h"
      real (kind=R8) ::
     * clvmin,clvmax,delta
      integer :: iclv, ncls
      external xerrab

      call subini ('eclvn')
      nofecls=ncls
      if(nofecls.eq.0) then
        call xerrab ('You set number of Eirene contour levels to 0')
      elseif(nofecls.eq.1) then
        eclval(1)=real((clvmax+clvmin)*0.5_R8,R4)
        write (*,*) 'You specified 1 Eirene contour level at',eclval(1)
      else
        delta=(clvmax-clvmin)/(nofecls-1)
        write (*,*) 'User specified',nofecls,'Eirene contour levels'
        write (*,*) 'Using customised Eirene contour levels:'
        do iclv=1,nofecls
          eclval(iclv)=real(clvmin+(iclv-1)*delta,R4)
          write (*,*) 'Eirene contour level',iclv,'at',eclval(iclv)
        enddo
      endif
      call subend ()
      return
      end subroutine eclvn
c
      subroutine wclvn(clvmin,clvmax,ncls)
      use b2mod_types
      use b2mod_subsys
      implicit none
#include "graphics_vars.h"
      real (kind=R8) ::
     * clvmin,clvmax,delta
      integer :: iclv, ncls
      external xerrab

      call subini ('wclvn')
      nofwcls=ncls
      if(nofwcls.eq.0) then
        call xerrab ('You set number of wall contour levels to 0')
      elseif(nofwcls.eq.1) then
        wclval(1)=real((clvmax+clvmin)*0.5_R8,R4)
        write (*,*) 'You specified 1 wall contour level at',wclval(1)
      else
        delta=(clvmax-clvmin)/(nofwcls-1)
        write (*,*) 'User specified',nofwcls,'wall contour levels'
        write (*,*) 'Using customised wall contour levels:'
        do iclv=1,nofwcls
          wclval(iclv)=real(clvmin+(iclv-1)*delta,R4)
          write (*,*) 'wall contour level',iclv,'at',wclval(iclv)
        enddo
      endif
      call subend ()
      return
      end subroutine wclvn
c
      subroutine tclvn(clvmin,clvmax,ncls)
      use b2mod_types
      use b2mod_subsys
      implicit none
#include "graphics_vars.h"
      real (kind=R8) ::
     * clvmin,clvmax,delta
      integer :: iclv, ncls
      external xerrab

      call subini ('wclvn')
      noftcls=ncls
      if(noftcls.eq.0) then
        call xerrab ('You set number of target contour levels to 0')
      elseif(noftcls.eq.1) then
        tclval(1)=real((clvmax+clvmin)*0.5_R8,R4)
        write (*,*) 'You specified 1 target contour level at',tclval(1)
      else
        delta=(clvmax-clvmin)/(noftcls-1)
        write (*,*) 'User specified',noftcls,'target contour levels'
        write (*,*) 'Using customised target contour levels:'
        do iclv=1,noftcls
          tclval(iclv)=real(clvmin+(iclv-1)*delta,R4)
          write (*,*) 'target contour level',iclv,'at',tclval(iclv)
        enddo
      endif
      call subend ()
      return
      end subroutine tclvn
