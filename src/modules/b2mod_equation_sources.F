      module b2mod_equation_sources                                     !srv 26.02.18 {
      use b2mod_types
      implicit none
       real (kind=R8), allocatable, save ::
     *  art_sch(:,:,:),
     *  art_she(:,:,:),
     *  art_shi(:,:,:),
     *  art_smo(:,:,:,:),
     *  art_sna(:,:,:,:)

      contains

      subroutine alloc_b2mod_equation_sources(nxd,nyd,nsd)
       implicit none
       integer nxd, nyd, nsd

       allocate(art_sch(-1:nxd,-1:nyd,0:3))
       allocate(art_she(-1:nxd,-1:nyd,0:3))
       allocate(art_shi(-1:nxd,-1:nyd,0:3))
       allocate(art_smo(-1:nxd,-1:nyd,0:3,0:nsd-1))
       allocate(art_sna(-1:nxd,-1:nyd,0:1,0:nsd-1))

       art_sch =0.0_R8
       art_she =0.0_R8
       art_shi =0.0_R8
       art_smo =0.0_R8
       art_sna =0.0_R8

       write(*,*)'Allocated space for b2mod_equation_sources'
       return

      end subroutine alloc_b2mod_equation_sources

      subroutine dealloc_b2mod_equation_sources
       implicit none

       deallocate(art_sch)
       deallocate(art_she)
       deallocate(art_shi)
       deallocate(art_smo)
       deallocate(art_sna)

       return
      end subroutine dealloc_b2mod_equation_sources

      subroutine read_equation_sources(nx, ny, ns)
       use b2mod_types
       implicit none
       integer nx, ny, ns
       integer is, ix, iy
       character*24 filename, hlp_frm
       logical art_exists

       write(hlp_frm,'(a,i4,a)') '(',nx+2,'(2x,e23.16))'

       filename = 'art_sna.dat'
       call find_file(filename,art_exists)
       if(art_exists) then
         open(70,FILE=trim(filename))
         read (70,hlp_frm)
     &    (((art_sna(ix,iy,0,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
         close(70)
       else
         call xerrab ('There is no art_sna.dat file')
       endif

       filename='art_smo.dat'
       call find_file(filename,art_exists)
       if(art_exists) then
         open(70,FILE=trim(filename))
         read (70,hlp_frm)
     &    (((art_smo(ix,iy,0,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
         close(70)
       else
         call xerrab('There is no art_smo.dat file')
       endif

       filename='art_she.dat'
       call find_file(filename,art_exists)
       if(art_exists) then
         open(70,FILE=trim(filename))
         read (70,hlp_frm)
     &    ((art_she(ix,iy,0),ix=-1,nx),iy=-1,ny)
         close(70)
       else
         call xerrab('There is no art_she.dat file')
       endif

       filename='art_shi.dat'
       call find_file(filename,art_exists)
       if(art_exists) then
         open(70,FILE=trim(filename))
         read (70,hlp_frm)
     &    ((art_shi(ix,iy,0),ix=-1,nx),iy=-1,ny)
         close(70)
       else
         call xerrab('There is no art_shi.dat file')
       endif

       filename='art_sch.dat'
       call find_file(filename,art_exists)
       if(art_exists) then
         open(70,FILE=trim(filename))
         read (70,hlp_frm)
     &    ((art_sch(ix,iy,0),ix=-1,nx),iy=-1,ny)
         close(70)
       else
         call xerrab('There is no art_sch.dat file')
       endif

       return
      end subroutine read_equation_sources


      end module b2mod_equation_sources                                 !srv 26.02.18 }
