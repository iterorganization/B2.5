c
      subroutine clvn(clvmin,clvmax,ncls)
      use b2mod_types
#include "graphics_vars.h"
      real (kind=R8) ::
     * clvmin,clvmax,delta
      nofcls=ncls
      if(nofcls.eq.0) then
        stop 'You set number of contour levels to 0'
      elseif(nofcls.eq.1) then
        clval(1)=(clvmax+clvmin)*0.5_R8
        write (*,*) 'You specified 1 contour level at',clval(1)
      else
        delta=(clvmax-clvmin)/(nofcls-1)
        write (*,*) 'User specified',nofcls,'contour levels'
        write (*,*) 'Using customised contour levels:'
        do iclv=1,nofcls
          clval(iclv)=clvmin+(iclv-1)*delta
          write (*,*) 'contour level',iclv,'at',clval(iclv)
        enddo
      endif
      return
      end subroutine clvn
c
      subroutine wclvn(clvmin,clvmax,ncls)
      use b2mod_types
#include "graphics_vars.h"
      real (kind=R8) ::
     * clvmin,clvmax,delta
      nofwcls=ncls
      if(nofwcls.eq.0) then
        stop 'You set number of wall contour levels to 0'
      elseif(nofwcls.eq.1) then
        wclval(1)=(clvmax+clvmin)*0.5_R8
        write (*,*) 'You specified 1 wall contour level at',wclval(1)
      else
        delta=(clvmax-clvmin)/(nofwcls-1)
        write (*,*) 'User specified',nofwcls,'wall contour levels'
        write (*,*) 'Using customised wall contour levels:'
        do iclv=1,nofwcls
          wclval(iclv)=clvmin+(iclv-1)*delta
          write (*,*) 'wall contour level',iclv,'at',wclval(iclv)
        enddo
      endif
      return
      end subroutine wclvn
c
      subroutine tclvn(clvmin,clvmax,ncls)
      use b2mod_types
#include "graphics_vars.h"
      real (kind=R8) ::
     * clvmin,clvmax,delta
      noftcls=ncls
      if(noftcls.eq.0) then
        stop 'You set number of target contour levels to 0'
      elseif(noftcls.eq.1) then
        tclval(1)=(clvmax+clvmin)*0.5_R8
        write (*,*) 'You specified 1 target contour level at',tclval(1)
      else
        delta=(clvmax-clvmin)/(noftcls-1)
        write (*,*) 'User specified',noftcls,'target contour levels'
        write (*,*) 'Using customised target contour levels:'
        do iclv=1,noftcls
          tclval(iclv)=clvmin+(iclv-1)*delta
          write (*,*) 'target contour level',iclv,'at',tclval(iclv)
        enddo
      endif
      return
      end subroutine tclvn
