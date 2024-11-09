c
      subroutine clvd(clvmin,clvmax,delta)
      use b2mod_types
      use b2mod_subsys
      use b2mod_graphics
      implicit none
      real (kind=R8) :: clvmin,clvmax,delta
      integer :: iclv
      external xertst, xerrab

      call subini ('clvd')
      call xertst(delta.ne.0.0_R8, 'division by 0 in clvd')
      nofcls=int((clvmax-clvmin)/delta)
      if(nofcls.eq.0) then
        call xerrab ('Contour Level Max == Contour Level Min')
      else
        nofcls=nofcls+1
        write (*,*) 'User specified',nofcls,'contour levels for range'
        write (*,*) 'Using customised contour level range:'
        do iclv=1,nofcls
          clval(iclv)=real(clvmin+(iclv-1)*delta,R4)
          write (*,*) 'contour level',iclv,'at',clval(iclv)
        enddo
      endif
      call subend ()
      return
      end subroutine clvd
c
      subroutine eclvd(clvmin,clvmax,delta)
      use b2mod_types
      use b2mod_subsys
      use b2mod_graphics
      implicit none
      real (kind=R8) :: clvmin,clvmax,delta
      integer :: iclv
      external xertst, xerrab

      call subini ('eclvd')
      call xertst (delta.ne.0.0_R8, 'division by 0 in eclvd')
      nofecls=int((clvmax-clvmin)/delta)
      if(nofecls.eq.0) then
        call xerrab
     .   ('Eirene Contour Level Max == Eirene Contour Level Min')
      else
        nofecls=nofecls+1
        write (*,*) 'User specified',nofecls,
     .   'Eirene contour levels for range'
        write (*,*) 'Using customised Eirene contour level range:'
        do iclv=1,nofecls
          eclval(iclv)=real(clvmin+(iclv-1)*delta,R4)
          write (*,*) 'Eirene contour level',iclv,'at',eclval(iclv)
        enddo
      endif
      call subend ()
      return
      end subroutine eclvd
c
      subroutine wclvd(clvmin,clvmax,delta)
      use b2mod_types
      use b2mod_subsys
      use b2mod_graphics
      implicit none
      real (kind=R8) :: clvmin,clvmax,delta
      integer :: iclv
      external xertst, xerrab

      call subini ('wclvd')
      call xertst (delta.ne.0.0_R8, 'division by 0 in wclvd')
      nofwcls=int((clvmax-clvmin)/delta)
      if(nofwcls.eq.0) then
        call xerrab ('Wall Contour Level Max == Wall Contour Level Min')
      else
        nofwcls=nofwcls+1
        write (*,*) 'User specified',nofwcls,
     .   'wall contour levels for range'
        write (*,*) 'Using customised wall contour level range:'
        do iclv=1,nofwcls
          wclval(iclv)=real(clvmin+(iclv-1)*delta,R4)
          write (*,*) 'wall contour level',iclv,'at',wclval(iclv)
        enddo
      endif
      call subend ()
      return
      end subroutine wclvd
c
      subroutine tclvd(clvmin,clvmax,delta)
      use b2mod_types
      use b2mod_subsys
      use b2mod_graphics
      implicit none
      real (kind=R8) :: clvmin,clvmax,delta
      integer :: iclv
      external xertst, xerrab

      call subini ('tcvld')
      call xertst(delta.ne.0.0_R8, 'division by 0 in tclvd')
      noftcls=int((clvmax-clvmin)/delta)
      if(noftcls.eq.0) then
        call xerrab 
     .   ('Target Contour Level Max == Target Contour Level Min')
      else
        noftcls=noftcls+1
        write (*,*) 'User specified',noftcls,
     .   'target contour levels for range'
        write (*,*) 'Using customised target contour level range:'
        do iclv=1,noftcls
          tclval(iclv)=real(clvmin+(iclv-1)*delta,R4)
          write (*,*) 'target contour level',iclv,'at',tclval(iclv)
        enddo
      endif
      call subend ()
      return
      end subroutine tclvd
 
