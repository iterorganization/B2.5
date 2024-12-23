c

      subroutine copy(src,scale,nz,texttmp)
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer nz
      real (kind=R8) ::
     *   src(mpg%nCv,0:nz-1),scale
      character*(*) texttmp(0:nz-1)
c      write(*,*) 'Copying array: nCv, nz = ',mpg%nCv,nz
c      write(*,*) nncut
      if(iwrk .eq. wklng) call dropstackbottom()
      iwrk=iwrk+1
      work(:,0:nz-1,iwrk)=src(:,0:nz-1)*scale
      values(9)=minval(work(:,0:nz-1,iwrk))
      values(10)=maxval(work(:,0:nz-1,iwrk))
      textwk(0:nz-1,iwrk)=texttmp(0:nz-1)
      ndim(1,iwrk)=mpg%nCv
      ndim(2,iwrk)=nz
c      write(*,*) src(1,0)*scale,work(1,0,iwrk)
      return
      end subroutine copy

#ifdef B25_EIRENE
      subroutine eirene_copy(src,scale,nz,texttmp)
! careful here: Eirene arrays are sized src(nz,ntrii) !
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      use eirmod_ctrig
      integer nz
      real (kind=R8) ::
     *   src(nz,1:ntrii),scale
      character*(*) texttmp(0:nz-1)
      integer i,iz
c      write(*,*) 'Copying Eirene array: ntrii, nz = ',ntrii,nz
      if(iewrk .eq. wklng) call dropeirenestackbottom()
      iewrk=iewrk+1
      do iz = 1, nz
        do i = 1, ntrii
          eirene_stack(i,iz-1,iewrk)=src(iz,i)*scale
        end do
      end do
      values(38)=minval(eirene_stack(1:ntrii,0:nz-1,iewrk))
      values(39)=maxval(eirene_stack(1:ntrii,0:nz-1,iewrk))
      textewk(0:nz-1,iewrk)=texttmp(0:nz-1)
      nedim(1,iewrk)=ntrii
      nedim(2,iewrk)=nz
c      write(*,*) src(1,1)*scale,eirene_stack(1,0,iewrk)
      return
      end subroutine eirene_copy

      subroutine eirene_coarse_copy(src,scale,nCv,nz,texttmp)
! special mapping of an Eirene array only defined over the B2.5 control volumes
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      use eirmod_ctrig
      use eirmod_wneutrals
      integer nCv, nz
      real (kind=R8) ::
     *   src(1:nCv,0:nz-1), scale
      character*(*) texttmp(0:nz-1)
      integer iCv
c      write(*,*) 'Copying Eirene array: ntrii, nCv, nz = ',ntrii,nCv,nz
      if(iewrk .eq. wklng) call dropeirenestackbottom()
      iewrk=iewrk+1
      do i = 1, ntrii
        iCv = ixtri(i)
        if (iCv.gt.0) then
          eirene_stack(i,0:nz-1,iewrk)=src(iCv,0:nz-1)*scale
        else
          eirene_stack(i,0:nz-1,iewrk)=0.0_R8
        endif
      end do
      values(38)=minval(eirene_stack(1:ntrii,0:nz-1,iewrk))
      values(39)=maxval(eirene_stack(1:ntrii,0:nz-1,iewrk))
      textewk(0:nz-1,iewrk)=texttmp(0:nz-1)
      nedim(1,iewrk)=ntrii
      nedim(2,iewrk)=nz
      return
      end subroutine eirene_coarse_copy
#endif

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
      use b2mod_subsys
      implicit none
      integer is
      external xertst

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
      
