c
      subroutine clvs(ncls,rstack,irstack)
      use b2mod_types
      use b2mod_graphics
      implicit none
      integer :: ncls, irstack, iclv
      real (kind=R8) :: rstack(irstack)
      nofcls=ncls
      if(nofcls.eq.0) then
        write (*,*) 'Using standard contour levels'
        return
      else
        write (*,*) 'User specified',nofcls,'contour levels'
        write (*,*) 'Using customised contour levels:'
        do iclv=1,nofcls
          clval(iclv)=real(rstack(irstack-nofcls+iclv),R4)
          write (*,*) 'contour level',iclv,'at',clval(iclv)
        enddo
      endif
      return
      end subroutine clvs
c
      subroutine eclvs(ncls,rstack,irstack)
      use b2mod_types
      use b2mod_graphics
      implicit none
      integer :: ncls, irstack, iclv
      real (kind=R8) :: rstack(irstack)
      nofecls=ncls
      if(nofecls.eq.0) then
        write (*,*) 'Using standard Eirene contour levels'
        return
      else
        write (*,*) 'User specified',nofecls,'Eirene contour levels'
        write (*,*) 'Using customised Eirene contour levels:'
        do iclv=1,nofecls
          eclval(iclv)=real(rstack(irstack-nofecls+iclv),R4)
          write (*,*) 'Eirene contour level',iclv,'at',eclval(iclv)
        enddo
      endif
      return
      end subroutine eclvs 
c
      subroutine wclvs(ncls,rstack,irstack)
      use b2mod_types
      use b2mod_graphics
      implicit none
      integer :: ncls, irstack, iclv
      real (kind=R8) :: rstack(irstack)
      nofwcls=ncls
      if(nofwcls.eq.0) then
        write (*,*) 'Using standard wall contour levels'
        return
      else
        write (*,*) 'User specified',nofwcls,'wall contour levels'
        write (*,*) 'Using customised wall contour levels:'
        do iclv=1,nofwcls
          wclval(iclv)=real(rstack(irstack-nofwcls+iclv),R4)
          write (*,*) 'wall contour level',iclv,'at',wclval(iclv)
        enddo
      endif
      return
      end subroutine wclvs 
c
      subroutine tclvs(ncls,rstack,irstack)
      use b2mod_types
      use b2mod_graphics
      implicit none
      integer :: ncls, irstack, iclv
      real (kind=R8) :: rstack(irstack)
      noftcls=ncls
      if(noftcls.eq.0) then
        write (*,*) 'Using standard target contour levels'
        return
      else
        write (*,*) 'User specified',noftcls,'target contour levels'
        write (*,*) 'Using customised target contour levels:'
        do iclv=1,noftcls
          tclval(iclv)=real(rstack(irstack-noftcls+iclv),R4)
          write (*,*) 'target contour level',iclv,'at',tclval(iclv)
        enddo
      endif
      return
      end subroutine tclvs
