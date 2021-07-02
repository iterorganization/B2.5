      subroutine mds_put_string(node_name,string)
      implicit none
      character*(*) node_name,string
#ifdef MDSPLUS
#ifndef G95
      interface
        integer function descr2(idtype,iflag,ilen)
        integer, intent(in) :: idtype
        integer, intent(in) :: iflag
        integer, intent(in), optional :: ilen
        end function descr2

        integer function mdsput2(name,units,description,string,iflg)
        character*(*), intent(in) :: name
        character*(*), intent(in) :: units
        integer, intent(in) :: description
        character*(*), intent(in) :: string
        integer, intent(in) :: iflg
        end function mdsput2
      end interface
#else
      integer descr2, mdsput2
#endif

#     include "mdslib_b2.inc"
      integer :: status

      if(len_trim(string).gt.0) then
        write(*,'(a,a,a,a,a)') 'mds_put_string to ',trim(node_name),
     1   ' of "',trim(string),'"'
#ifndef G95
        status=mdsput2(trim(node_name)//char(0),'$'//char(0),
     1   descr2(idtype_cstring,iflag=0,ilen=len_trim(string)),
     1   trim(string),0)
#else
        status=mdsput2(trim(node_name)//char(0),'$'//char(0),
     1   descr2(idtype_cstring,0,len_trim(string)),
     1   trim(string),0)
#endif
        call mds_err_msg(status)
      else
        write(*,*) 'Skipped writing zero length string to ',
     1   trim(node_name)
      endif

#else
      call xerrab('mds_put_string: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_string

      subroutine mds_put_string_ns(node_name,ns,string)
      implicit none
      integer ns
      character*(*) node_name,string(ns)
#ifdef MDSPLUS
#ifndef G95
      interface
        integer function descr2(idtype,il1,iflag,ilen)
        integer, intent(in) :: idtype
        integer, intent(in) :: il1
        integer, intent(in) :: iflag
        integer, intent(in), optional :: ilen
        end function descr2
      end interface
#     include "mdslib_b2.inc"
      integer mdsput2
#else
#     include "mdslib.inc"
#endif
      integer is
      integer :: status
      character*18 hlp_frm

      write(hlp_frm,'(a,i3,a)') '(a,a,a,',ns,'(a,a),a)'
      write(*,hlp_frm) 'mds_put_string to ',trim(node_name),
     1 ' of "',(' ',trim(string(is)),is=1,ns),'"'
      do is=1,ns
        string(is)=trim(string(is)) // char(0)
      enddo
#ifndef G95
      status=mdsput2(trim(node_name)//char(0),'$'//char(0),
     1 descr2(idtype_cstring,ns,iflag=0,ilen=len(string(1))),string,0)
#else
      status=mdsput2(trim(node_name)//char(0),'$'//char(0),
     1 descr2(idtype_cstring,ns,0,len(string(1))),string,0)
#endif
      call mds_err_msg(status)

#else
      call xerrab('mds_put_string_ns: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_string_ns

      subroutine mds_put_int(node_name,int)
      implicit none
      character*(*) node_name
      integer int
#ifdef MDSPLUS
#ifndef G95
      interface
        integer function mdsput2(name,units,description,int,iflg)
        character*(*), intent(in) :: name
        character*(*), intent(in) :: units
        integer, intent(in) :: description
        integer, intent(in) :: int
        integer, intent(in) :: iflg
        end function mdsput2
      end interface
      integer descr2
#     include "mdslib_b2.inc"
#else
#     include "mdslib.inc"
#endif
      integer :: status

      write(*,'(a,a,a,i8)') 'mds_put_int to ',trim(node_name),
     1 ' of ',int
      status=mdsput2(trim(node_name)//char(0),'$'//char(0),
     1 descr2(idtype_long,0),int,0)
      call mds_err_msg(status)

#else
      call xerrab('mds_put_int: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_int

      subroutine mds_put_int_rz(node_name,nx,ny,data)
      use b2mod_types
      implicit none
      character*(*) :: node_name
      integer :: nx, ny, data(nx,ny)

#ifdef MDSPLUS
#ifndef G95
      interface
        integer function descr2(idtype,i1,i2,iflag)
        integer, intent(in) :: idtype
        integer, intent(in) :: i1
        integer, intent(in) :: i2
        integer, intent(in) :: iflag
        end function descr2

        integer function mdsput2(name,units,description,i2d,iflg)
        integer nx, ny
        character*(*), intent(in) :: name
        character*(*), intent(in) :: units
        integer, intent(in) :: description
        integer, intent(in), allocatable :: i2d(:,:)
        integer, intent(in) :: iflg
        end function mdsput2
      end interface
#     include "mdslib_b2.inc"
#else
#     include "mdslib.inc"
#endif
      integer :: status
      character :: tdi*(172)
      integer, allocatable :: i2d(:,:)

      allocate(i2d(nx,ny))
      i2d=data
      tdi='$'

      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_long,nx,ny,0),i2d,0)
      write(*,'(a,a,a,a)')
     * 'mdsput2 of ',trim(node_name),' with ',trim(tdi)
      call mds_err_msg(status)
      deallocate(i2d)

#else
      call xerrab('mds_put_int_rz: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_int_rz

      subroutine mds_put_int_ns(node_name,ns,data)
      use b2mod_types
      implicit none
      character*(*) :: node_name
      integer :: ns, data(ns)

#ifdef MDSPLUS
#ifndef G95
      interface
        integer function descr2(idtype,i1,iflag)
        integer, intent(in) :: idtype
        integer, intent(in) :: i1
        integer, intent(in) :: iflag
        end function descr2

        integer function mdsput2(name,units,description,i1d,iflg)
        character*(*), intent(in) :: name
        character*(*), intent(in) :: units
        integer, intent(in) :: description
        integer, intent(in), allocatable :: i1d(:)
        integer, intent(in) :: iflg
        end function mdsput2
      end interface
#     include "mdslib_b2.inc"
#else
#     include "mdslib.inc"
#endif
      integer :: status
      character :: tdi*(172)
      integer, allocatable :: i1d(:)

      allocate(i1d(ns))
      i1d=data
      tdi='$'

      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_long,ns,0),i1d,0)
      write(*,'(a,a,a,a)')
     * 'mdsput2 of ',trim(node_name),' with ',trim(tdi)
      call mds_err_msg(status)
      deallocate(i1d)

#else
      call xerrab('mds_put_int_ns: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_int_ns

      subroutine mds_put_double(node_name,units,value)
      use b2mod_types
      implicit none
      character*(*) node_name, units
      real (kind=R8) :: value
#ifdef MDSPLUS
#ifndef G95
      interface
        integer function descr2(idtype,iflag)
        integer, intent(in) :: idtype
        integer, intent(in) :: iflag
        end function descr2

        integer function mdsput2(name,units,description,value,iflg)
        use b2mod_types
        character*(*), intent(in) :: name
        character*(*), intent(in) :: units
        integer, intent(in) :: description
        real(kind=R8), intent(in) :: value
        integer, intent(in) :: iflg
        end function mdsput2
      end interface
#     include "mdslib_b2.inc"
#else
#     include "mdslib.inc"
#endif
      integer :: status
      character :: tdi*(172)

      tdi='build_with_units($,"'//trim(units)//'")'
      write(*,'(a,a,a,1p,1e22.15,a,a)')
     * 'mds_put_double to ',trim(node_name),' of ',value,
     * ' with ',trim(tdi)
      status=mdsput2(trim(node_name)//char(0),trim(tdi)//char(0),
     1 descr2(idtype_double,0),value,0)
      call mds_err_msg(status)

#else
      call xerrab('mds_put_double: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_double

      subroutine mds_put_double_ns(node_name,units,ns,data)
      use b2mod_types
      implicit none
      character*(*) :: node_name, units
      integer :: ns
      real (kind=R8) :: data(ns)

#ifdef MDSPLUS
#ifndef G95
      interface
        integer function descr2(idtype,i1,iflag)
        integer, intent(in) :: idtype
        integer, intent(in) :: i1
        integer, intent(in) :: iflag
        end function descr2

        integer function mdsput2(name,units,description,d1d,iflg)
        use b2mod_types
        character*(*), intent(in) :: name
        character*(*), intent(in) :: units
        integer, intent(in) :: description
        real(kind=R8), intent(in), allocatable :: d1d(:)
        integer, intent(in) :: iflg
        end function mdsput2
      end interface
#     include "mdslib_b2.inc"
#else
#     include "mdslib.inc"
#endif
      integer :: status
      character :: tdi*(172)
      real(kind=R8), allocatable :: d1d(:)

      allocate(d1d(ns))
      d1d=data
      tdi='build_with_units($,"'//trim(units)//'")'
      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,ns,0),d1d,0)
      write(*,'(a,a,a,a)')
     * 'mdsput2 of ',trim(node_name),' with ',trim(tdi)
      call mds_err_msg(status)
      deallocate(d1d)

#else
      call xerrab('mds_put_double_ns: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_double_ns

      subroutine mds_put_double_rz(node_name,units,nx,ny,data,dirsuffix)
      use b2mod_types
      implicit none
      character*(*) :: node_name, units, dirsuffix
      integer :: nx, ny
      real (kind=R8) :: data(nx,ny)

#ifdef MDSPLUS
#ifndef G95
      interface
        integer function descr2(idtype,i1,i2,iflag)
        integer, intent(in) :: idtype
        integer, intent(in) :: i1
        integer, intent(in) :: i2
        integer, intent(in) :: iflag
        end function descr2

        integer function mdsput2(name,units,description,d2d,iflg)
        use b2mod_types
        character*(*), intent(in) :: name
        character*(*), intent(in) :: units
        integer, intent(in) :: description
        real(kind=R8), intent(in), allocatable :: d2d(:,:)
        integer, intent(in) :: iflg
        end function mdsput2
      end interface
#     include "mdslib_b2.inc"
#else
#     include "mdslib.inc"
#endif
      integer :: status
      character :: tdi*(172)
      real(kind=R8), allocatable :: d2d(:,:)

      allocate(d2d(nx,ny))
      d2d=data
      tdi='build_signal(build_with_units($,"' //
     1 trim(units) //
     2 '"),,\TOP.SNAPSHOT.GRID:CR' // trim(dirsuffix) //
     3 ',\TOP.SNAPSHOT.GRID:CZ' // trim(dirsuffix) //
     4 ')'

      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nx,ny,0),d2d,0)
      write(*,'(a,a,a,a)')
     * 'mdsput2 of ',trim(node_name),' with ',trim(tdi)
      call mds_err_msg(status)
      deallocate(d2d)

#else
      call xerrab('mds_put_double_rz: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_double_rz

      subroutine mds_put_double_rzns(node_name,units,nx,ny,ns,
     1 data,nsindex,dirsuffix)
      use b2mod_types
      implicit none
      character*(*) :: node_name, units, nsindex, dirsuffix
      integer :: nx, ny, ns
      real (kind=R8) :: data(nx,ny,ns)

#ifdef MDSPLUS
#ifndef G95
      interface
        integer function descr2(idtype,i1,i2,i3,iflag)
        integer, intent(in) :: idtype
        integer, intent(in) :: i1
        integer, intent(in) :: i2
        integer, intent(in) :: i3
        integer, intent(in) :: iflag
        end function descr2

        integer function mdsput2(name,units,description,d3d,iflg)
        use b2mod_types
        character*(*), intent(in) :: name
        character*(*), intent(in) :: units
        integer, intent(in) :: description
        real(kind=R8), intent(in), allocatable :: d3d(:,:,:)
        integer, intent(in) :: iflg
        end function mdsput2
      end interface
#     include "mdslib_b2.inc"
#else
#     include "mdslib.inc"
#endif
      integer :: status
      character :: tdi*(172)
      real (kind=R8), allocatable :: d3d(:,:,:)

      allocate(d3d(nx,ny,ns))
      d3d=data
      tdi='build_signal(build_with_units($,"' //
     1 trim(units) //
     2 '"),,\TOP.SNAPSHOT.GRID:CR' // trim(dirsuffix) //
     3 ',\TOP.SNAPSHOT.GRID:CZ' // trim(dirsuffix) //
     4 ',' //
     5 '\TOP.SNAPSHOT:' // trim(nsindex) // ')'

      write(*,'(a,a,a,a)')
     * 'mdsput2 of ',trim(node_name),' with ',trim(tdi)
      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nx,ny,ns,0),d3d,0)
      call mds_err_msg(status)
      deallocate(d3d)

#else
      call xerrab('mds_put_double_rzns: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_double_rzns

      subroutine mds_put_double_rzdns(node_name,units,nx,ny,nd,ns,
     1 data,dirindex,nsindex,dirsuffix)
      use b2mod_types
      implicit none
      character*(*) :: node_name, units, dirindex, nsindex, dirsuffix
      integer :: nx, ny, nd, ns
      real (kind=R8) :: data(nx,ny,nd,ns)

#ifdef MDSPLUS
#ifndef G95
      interface
        integer function descr2(idtype,i1,i2,i3,i4,iflag)
        integer, intent(in) :: idtype
        integer, intent(in) :: i1
        integer, intent(in) :: i2
        integer, intent(in) :: i3
        integer, intent(in) :: i4
        integer, intent(in) :: iflag
        end function descr2

        integer function mdsput2(name,units,description,d4d,iflg)
        use b2mod_types
        character*(*), intent(in) :: name
        character*(*), intent(in) :: units
        integer, intent(in) :: description
        real(kind=R8), intent(in), allocatable :: d4d(:,:,:,:)
        integer, intent(in) :: iflg
        end function mdsput2
      end interface
#     include "mdslib_b2.inc"
#else
#     include "mdslib.inc"
#endif
      integer :: status
      character :: tdi*(172)
      real (kind=R8), allocatable :: d4d(:,:,:,:)

      allocate(d4d(nx,ny,nd,ns))
      d4d=data
      tdi='build_signal(build_with_units($,"' //
     1 trim(units) //
     2 '"),,\TOP.SNAPSHOT.GRID:CR' // trim(dirsuffix) //
     3 ',\TOP.SNAPSHOT.GRID:CZ' // trim(dirsuffix) //
     4 ',' //
     5 '\TOP.SNAPSHOT:' // trim(dirindex) //
     6 ',' //
     7 '\TOP.SNAPSHOT:' // trim(nsindex) // ')'

      write(*,'(a,a,a,a)')
     * 'mdsput2 of ',trim(node_name),' with ',trim(tdi)
      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nx,ny,nd,ns,0),d4d,0)
      call mds_err_msg(status)
      deallocate(d4d)

#else
      call xerrab('mds_put_double_rzdns: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_double_rzdns

      subroutine mds_err_msg(status)
      implicit none
      integer status
#ifdef MDSPLUS
#ifndef G95
      interface
        integer function descr2(idtype,il1,iflag,ilen)
        integer, intent(in) :: idtype
        integer, intent(in), optional :: il1
        integer, intent(in) :: iflag
        integer, intent(in), optional :: ilen
        end function descr2
      end interface
#     include "mdslib_b2.inc"
      integer mdsvalue2
#else
#     include "mdslib.inc"
#endif
      integer :: size
      character :: message*80

      if(mod(status,2).eq.0) then
#ifndef G95
        status=mdsvalue2('getmsg($)'//char(0),
     1   descr2(idtype_long,iflag=0),status,
     2   descr2(idtype_cstring,iflag=0,ilen=80),message,0,size)
#else
        status=mdsvalue2('getmsg($)'//char(0),
     1   descr2(idtype_long,0),status,
     2   descr2(idtype_cstring,0,80),message,0,size)
#endif
        write(*,*) trim(message)
      endif
#else
      call xerrab('mds_err_msg: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_err_msg
