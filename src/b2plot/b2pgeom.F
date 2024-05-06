      subroutine b2pgeom(lun)
c
c  version : 07.09.2000 16:07
c=======================================================================
c***  Produces tables with coordinates of grid cell corners along the
c***  separatrix and along the magnetic surfaces
c***  for subsequent plotting with standalone software
c***
c*** Output files:
c***
c***  gm_sep.dat first separatrix
c***  gm_SOL.dat magnetic surfaces in SOL
c***  gm_PFR.dat magnetic surfaces in PFR
c***  gm_core.dat magnetic surfaces in the core part of the grid
c***
c=======================================================================
      use b2mod_types
      use b2us_data
      use b2mod_diag
      implicit none
      integer lun

c*** Local variables

      integer i,j,k,l,iCv,iVx
      integer, allocatable :: fcList(:), vxList(:,:)
      character hlp_frmr*100, tab, ofile*12
c=======================================================================

      tab=char(9)

      if(mpg%iFssep.ne.0) then !{

        ofile='gm_sep.dat'
        open(lun,file=ofile,err=999)
        rewind(lun)
        write(lun,' (a)',err=998) '# Separatrix shape'
        write(lun,'(2a)',err=998) '# ',trim(label)
        write(lun,'(3a)',err=998) '# R',tab,'Z'
        allocate(fcList(mpg%fsFcP(mpg%iFssep,2)))
        call sort_faces(
     ,   mpg%fsFcP(mpg%iFssep,2),fcList,mpg%nFc,mpg%fcFs,mpg%iFssep,mpg)
        do i=1,mpg%fsFcP(mpg%iFssep,2) !{
          iVx = mpg%fcVx(fcList(i),1)
          write(lun,'(1p,e15.7,a,e15.7)',err=998) geo%cvX(iVx),tab,
     ,                                                    geo%cvY(iVx)
        end do !}
        iVx = mpg%fcVx(fcList(mpg%fsFcP(mpg%iFssep,2)),2)
        write(lun,'(1p,e15.7,a,e15.7)',err=998) geo%cvX(iVx),tab,
     ,                                                    geo%cvY(iVx)
        deallocate(fcList)
        if (mpg%iFssep2.gt.0) then
          call sort_faces(
     ,     mpg%fsFcP(mpg%iFssep2,2),fcList,mpg%nFc,mpg%fcFs,
     ,     mpg%iFssep2,mpg)
          do i=1,mpg%fsFcP(mpg%iFssep2,2) !{
            iVx = mpg%fcVx(fcList(i),1)
            write(lun,'(1p,e15.7,a,e15.7)',err=998) geo%cvX(iVx),tab,
     ,                                                    geo%cvY(iVx)
          end do !}
          iVx = mpg%fcVx(fcList(mpg%fsFcP(mpg%iFssep2,2)),2)
          write(lun,'(1p,e15.7,a,e15.7)',err=998) geo%cvX(iVx),tab,
     ,                                                    geo%cvY(iVx)
          deallocate(fcList)
        end if
        close(lun)

        ofile='gm_SOL.dat'
        open(lun,file=ofile,err=999)
        rewind(lun)
        write(lun,' (a)',err=998) '# SOL magnetic surfaces'
        write(lun,'(2a)',err=998) '# ',trim(label)
        k=0
        l=0
        do i=1,mpg%nFs
          if (i.eq.mpg%iFssep.or.i.eq.mpg%iFssep2) cycle
          if (geo%fsPsi(i).lt.geo%fsPsi(mpg%iFssep).and.
     &        geo%psi_increasing) cycle
          if (geo%fsPsi(i).gt.geo%fsPsi(mpg%iFssep).and.
     &        .not.geo%psi_increasing) cycle
          k=k+1
          l=max(l,mpg%fsFcP(i,2))
        end do
        if (k.gt.0) then !{
          write(hlp_frmr,'(a,i3,a)')
     .     '(a,4x,',k,'(a1,i2.2,13x,a1,i2.2,13x))'
          write(lun,hlp_frmr,err=998) '#',('R',i,'Z',i,i=1,k)
          write(hlp_frmr,'(a,i3,a)') '(1p,',2*k,'(e14.7,2x))'
          allocate(vxList(l+1,k))
          vxList=-999.0_R8
          k=0
          do i=1,mpg%nFs !{
            if (i.eq.mpg%iFssep.or.i.eq.mpg%iFssep2) cycle
            if (geo%fsPsi(i).lt.geo%fsPsi(mpg%iFssep).and.
     &          geo%psi_increasing) cycle
            if (geo%fsPsi(i).gt.geo%fsPsi(mpg%iFssep).and.
     &          .not.geo%psi_increasing) cycle
            k=k+1
            allocate(fcList(mpg%fsFcP(i,2)))
            call sort_faces(
     .       mpg%fsFcP(i,2),fcList,mpg%nFc,mpg%fcFs,i,mpg)
            vxList(1:mpg%fsFcP(i,2),k)=
     .       mpg%fcVx(fcList(1:mpg%fsFcP(i,2)),1)
            vxList(mpg%fsFcP(i,2)+1,k)=
     .       mpg%fcVx(fcList(mpg%fsFcP(i,2)),2)
            deallocate(fcList)
          end do !}
          do i=1,l+1 !{
            write(lun,hlp_frmr,err=998) (geo%cvX(vxList(i,j)),
     .                                   geo%cvY(vxList(i,j)),j=1,k)
          end do !}
          deallocate(vxList)
        end if !}
        close(lun)

        ofile='gm_PFR.dat'
        open(lun,file=ofile,err=999)
        rewind(lun)
        write(lun,' (a)',err=998) '# PFR magnetic surfaces'
        write(lun,'(2a)',err=998) '# ',trim(label)
        k=0
        l=0
        do i=1,mpg%nFs
          if (i.eq.mpg%iFssep.or.i.eq.mpg%iFssep2) cycle
          if (geo%fsPsi(i).gt.geo%fsPsi(mpg%iFssep).and.
     &        geo%psi_increasing) cycle
          if (geo%fsPsi(i).lt.geo%fsPsi(mpg%iFssep).and.
     &        .not.geo%psi_increasing) cycle
          iCv = min(mpg%fcCv(mpg%fsFc(mpg%fsFcP(i,1)),1),
     &              mpg%fcCv(mpg%fsFc(mpg%fsFcP(i,1)),2))
          if (mpg%cvOnClosedSurface(iCv)) cycle
          k=k+1
          l=max(l,mpg%fsFcP(i,2))
        end do
        if (k.gt.0) then !{
          write(hlp_frmr,'(a,i3,a)')
     ,     '(a,4x,',k,'(a1,i2.2,13x,a1,i2.2,13x))'
          write(lun,hlp_frmr,err=998) '#',('R',i,'Z',i,i=1,k)
          write(hlp_frmr,'(a,i3,a)') '(1p,',2*k,'(e14.7,2x))'
          allocate(vxList(l+1,k))
          vxList=-999.0_R8
          k=0
          do i=1,mpg%nFs !{
            if (i.eq.mpg%iFssep.or.i.eq.mpg%iFssep2) cycle
            if (geo%fsPsi(i).gt.geo%fsPsi(mpg%iFssep).and.
     &          geo%psi_increasing) cycle
            if (geo%fsPsi(i).lt.geo%fsPsi(mpg%iFssep).and.
     &          .not.geo%psi_increasing) cycle
            iCv = min(mpg%fcCv(mpg%fsFc(mpg%fsFcP(i,1)),1),
     &                mpg%fcCv(mpg%fsFc(mpg%fsFcP(i,1)),2))
            if (mpg%cvOnClosedSurface(iCv)) cycle
            k=k+1
            allocate(fcList(mpg%fsFcP(i,2)))
            call sort_faces(
     .       mpg%fsFcP(i,2),fcList,mpg%nFc,mpg%fcFs,i,mpg)
            vxList(1:mpg%fsFcP(i,2),k)=
     .       mpg%fcVx(fcList(1:mpg%fsFcP(i,2)),1)
            vxList(mpg%fsFcP(i,2)+1,k)=
     .       mpg%fcVx(fcList(mpg%fsFcP(i,2)),2)
            deallocate(fcList)
          end do !}
          do i=1,l+1 !{
            write(lun,hlp_frmr,err=998) (geo%cvX(vxList(i,j)),
     .                                   geo%cvY(vxList(i,j)),j=1,k)
          end do !}
          deallocate(vxList)
        end if !}
        close(lun)

        ofile='gm_core.dat'
        open(lun,file=ofile,err=999)
        rewind(lun)
        write(lun,' (a)',err=998) '# Core magnetic surfaces'
        write(lun,'(2a)',err=998) '# ',trim(label)
        k=0
        l=0
        do i=1,mpg%nFs
          if (i.eq.mpg%iFssep.or.i.eq.mpg%iFssep2) cycle
          if (geo%fsPsi(i).gt.geo%fsPsi(mpg%iFssep).and.
     &        geo%psi_increasing) cycle
          if (geo%fsPsi(i).lt.geo%fsPsi(mpg%iFssep).and.
     &        .not.geo%psi_increasing) cycle
          iCv = min(mpg%fcCv(mpg%fsFc(mpg%fsFcP(i,1)),1),
     &              mpg%fcCv(mpg%fsFc(mpg%fsFcP(i,1)),2))
          if (.not.mpg%cvOnClosedSurface(iCv)) cycle
          k=k+1
          l=max(l,mpg%fsFcP(i,2))
        end do
        if (k.gt.0) then !{
          write(hlp_frmr,'(a,i3,a)')
     ,     '(a,4x,',k,'(a1,i2.2,13x,a1,i2.2,13x))'
          write(lun,hlp_frmr,err=998) '#',('R',i,'Z',i,i=1,k)
          write(hlp_frmr,'(a,i3,a)') '(1p,',2*k,'(e14.7,2x))'
          allocate(vxList(l,k))
          vxList=-999.0_R8
          k=0
          do i=1,mpg%nFs !{
            if (i.eq.mpg%iFssep.or.i.eq.mpg%iFssep2) cycle
            if (geo%fsPsi(i).gt.geo%fsPsi(mpg%iFssep).and.
     &          geo%psi_increasing) cycle
            if (geo%fsPsi(i).lt.geo%fsPsi(mpg%iFssep).and.
     &          .not.geo%psi_increasing) cycle
            iCv = min(mpg%fcCv(mpg%fsFc(mpg%fsFcP(i,1)),1),
     &                mpg%fcCv(mpg%fsFc(mpg%fsFcP(i,1)),2))
            if (.not.mpg%cvOnClosedSurface(iCv)) cycle
            k=k+1
            allocate(fcList(mpg%fsFcP(i,2)))
            call sort_faces(
     .       mpg%fsFcP(i,2),fcList,mpg%nFc,mpg%fcFs,i,mpg)
            vxList(1:mpg%fsFcP(i,2),k)=
     .       mpg%fcVx(fcList(1:mpg%fsFcP(i,2)),1)
            deallocate(fcList)
          end do !}
          do i=1,l !{
            write(lun,hlp_frmr,err=998) (geo%cvX(vxList(i,j)),
     .                                   geo%cvY(vxList(i,j)),j=1,k)
          end do !}
          deallocate(vxList)
        end if !}
        close (lun)

      else !}{
        write(0,*) 'b2pgeom: no separatrix flux surface declared! '
        return
      end if !}

      return
c=======================================================================
  998 write(0,*) 'b2pgeom: error writing to the file ',ofile
      return
c=======================================================================
  999 write(0,*) 'b2pgeom: unable to open file ',ofile
      return
c=======================================================================
      end subroutine b2pgeom
