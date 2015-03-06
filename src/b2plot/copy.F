c

      subroutine copy(src,scale,nz,texttmp)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer nz
      real (kind=R8) ::
     *   src(-1:nx,-1:ny,0:nz-1),scale
      character*(*) texttmp(0:nz-1)
c      write(*,*) 'Copying array: nx, ny, nz = ',nx,ny,nz
c      write(*,*) 'nxd, nyd = ',nxd,nyd
c      write(*,*) nncut
      if(iwrk .eq. wklng) call dropstackbottom()
      iwrk=iwrk+1
      work(:,:,0:nz-1,iwrk)=src(:,:,0:nz-1)*scale
      values(9)=minval(work(:,:,0:nz-1,iwrk))
      values(10)=maxval(work(:,:,0:nz-1,iwrk))
      textwk(0:nz-1,iwrk)=texttmp(0:nz-1)
      ndim(1,iwrk)=nx
      ndim(2,iwrk)=ny
      ndim(3,iwrk)=nz
c      write(*,*) src(0,0,1)*scale,work(0,0,1,iwrk),
c     1 src(nx+1,ny+1,1)*scale,work(nx+1,ny+1,1,iwrk)
      return
      end subroutine copy

      subroutine wall_copy(src,scale,nz,texttmp)
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer nz
      real (kind=R8) ::
     *   src(1:nwall,0:nz-1),scale
      character*(*) texttmp(0:nz-1)
c      write(*,*) 'Copying wall array: nwall, nz = ',nwall,nz
      if(iwwrk .eq. wklng) call dropwallstackbottom()
      iwwrk=iwwrk+1
      wall_stack(1:nbnd,0:nz-1,iwwrk)=src(1:nbnd,0:nz-1)*scale
      values(14)=minval(wall_stack(1:nbnd,0:nz-1,iwwrk))
      values(15)=maxval(wall_stack(1:nbnd,0:nz-1,iwwrk))
      textwwk(0:nz-1,iwwrk)=texttmp(0:nz-1)
      nwdim(1,iwwrk)=nbnd
      nwdim(2,iwwrk)=nz
c      write(*,*) src(1,0)*scale,wall_stack(1,0,iwwrk)
      return
      end subroutine wall_copy

      subroutine target_copy(src,scale)
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      implicit none
      real (kind=R8) ::
     *   src(1:nwall,1:ndepth),scale
      if(b2plot_debug.gt.13) 
     > write(*,*) 'Copying target array: nbnd, ndepth_nml = ',
     . nbnd,ndepth_nml
      if(itwrk .eq. wklng) call droptargetstackbottom()
      itwrk=itwrk+1
      target_stack(1:nbnd,1:ndepth_nml,itwrk)=
     . src(1:nbnd,1:ndepth_nml)*scale
      values(32)=minval(target_stack(1:nbnd,1:ndepth_nml,itwrk))
      values(33)=maxval(target_stack(1:nbnd,1:ndepth_nml,itwrk))
      if(b2plot_debug.gt.13) 
     > write(*,*) src(1,1)*scale,target_stack(1,1,itwrk)
      return
      end subroutine target_copy 

      subroutine data_copy(is)
c  Copy the top of the wall stack to the data array
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_ank
      implicit none
      integer is
      external subini, subend, xertst

      call subini('data_copy')
      if (iwwrk.eq.0) then
        call subend ()
        return  ! wall stack empty : do nothing
      end if
      ndata = nwdim(1,iwwrk)
      call xertst (is.ge.0 .and. is.lt.nwdim(2,iwwrk),
     . 'Invalid IS input in data_copy')
      data(1:ndata) = wall_stack(1:ndata,is,iwwrk)
      call subend()

      return
      end subroutine data_copy
      
