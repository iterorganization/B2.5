      subroutine mds_put_string(node_name,string)
      use b2mod_mdsplus
      implicit none
      character*(*) node_name,string
#ifdef MDSPLUS
#     include "mdslib.inc"
      integer :: status

      if(len_trim(string).gt.0) then
        write(*,'(a,a,a,a,a)') 'mds_put_string to ',trim(node_name),
     1   ' of "',trim(string),'"'
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsput2_str(trim(node_name)//char(0),'$'//char(0),
     1   descr2_str0(idtype_cstring,0,len_trim(string)),
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
      use b2mod_mdsplus
      implicit none
      integer ns
      character*(*) node_name,string(ns)
#ifdef MDSPLUS
#     include "mdslib.inc"
      integer is
      integer :: status
      integer :: strlen
      character*18 hlp_frm
#ifdef F2003
      character(len=:), allocatable :: str1(:)
#else
      integer, parameter :: strmax = 1024
      character(len=strmax+1), allocatable :: str1(:)
#endif
#ifndef GFORTRAN
      integer lnblnk
      external lnblnk
#endif

      write(hlp_frm,'(a,i3,a)') '(a,a,a,',ns,'(a,a),a)'
      write(*,hlp_frm) 'mds_put_string to ',trim(node_name),
     1 ' of "',(' ',trim(string(is)),is=1,ns),'"'
      strlen=0
      do is=1,ns
        strlen = max(strlen,lnblnk(string(is)))
      enddo
#ifdef F2003
      allocate(character(len=strlen+1) :: str1(ns))
#else
      call xertst(strlen.le.strmax,
     &    'increase size of strmax in mds_routines')
      allocate(str1(ns))
#endif
      do is=1,ns
        str1(is)=trim(string(is)) // char(0)
      end do
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2_str1(trim(node_name)//char(0),'$'//char(0),
     1 descr2_str1(idtype_cstring,ns,0,strlen),str1,0)
#else
      status=mdsput2(trim(node_name)//char(0),'$'//char(0),
     1 descr2(idtype_cstring,ns,0,strlen),string,0)
#endif
      call mds_err_msg(status)
      deallocate(str1)
#else
      call xerrab('mds_put_string_ns: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_string_ns

      subroutine mds_put_int(node_name,int)
      use b2mod_mdsplus
      implicit none
      character*(*) node_name
      integer int
#ifdef MDSPLUS
#     include "mdslib.inc"
      integer :: status

      write(*,'(a,a,a,i8)') 'mds_put_int to ',trim(node_name),
     1 ' of ',int
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2_i0(trim(node_name)//char(0),'$'//char(0),
     1 descr2_i0(idtype_long,0),int,0)
#else
      status=mdsput2(trim(node_name)//char(0),'$'//char(0),
     1 descr2(idtype_long,0),int,0)
#endif
      call mds_err_msg(status)
#else
      call xerrab('mds_put_int: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_int

      subroutine mds_put_int_us(node_name,n,data)
      use b2mod_types
      use b2mod_mdsplus
      implicit none
      character*(*) :: node_name
      integer :: n, data(n)

#ifdef MDSPLUS
#     include "mdslib.inc"
      integer :: status
      character :: tdi*(172)
      integer, allocatable :: i1d(:)

      allocate(i1d(n))
      i1d=data
      tdi='$'

#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2_i1d(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2_i1(idtype_long,n,0),i1d,0)
#else
      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_long,n,0),i1d,0)
#endif
      write(*,'(a,a,a,a)')
     * 'mdsput2 of ',trim(node_name),' with ',trim(tdi)
      call mds_err_msg(status)
      deallocate(i1d)
#else
      call xerrab('mds_put_int_us: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_int_us

      subroutine mds_put_int_ns(node_name,ns,data)
      use b2mod_types
      use b2mod_mdsplus
      implicit none
      character*(*) :: node_name
      integer :: ns, data(ns)

#ifdef MDSPLUS
#     include "mdslib.inc"
      integer :: status
      character :: tdi*(172)
      integer, allocatable :: i1d(:)

      allocate(i1d(ns))
      i1d=data
      tdi='$'

#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2_i1d(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2_i1(idtype_long,ns,0),i1d,0)
#else
      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_long,ns,0),i1d,0)
#endif
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
      use b2mod_mdsplus
      implicit none
      character*(*) node_name, units
      real (kind=R8) :: value
#ifdef MDSPLUS
#     include "mdslib.inc"
      integer :: status
      character :: tdi*(172)

      tdi='build_with_units($,"'//trim(units)//'")'
      write(*,'(a,a,a,1p,1e22.15,a,a)')
     * 'mds_put_double to ',trim(node_name),' of ',value,
     * ' with ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2_d0(trim(node_name)//char(0),trim(tdi)//char(0),
     1 descr2_i0(idtype_double,0),value,0)
#else
      status=mdsput2(trim(node_name)//char(0),trim(tdi)//char(0),
     1 descr2(idtype_double,0),value,0)
#endif
      call mds_err_msg(status)
#else
      call xerrab('mds_put_double: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_double

      subroutine mds_put_double_ns(node_name,units,ns,data)
      use b2mod_types
      use b2mod_mdsplus
      implicit none
      character*(*) :: node_name, units
      integer :: ns
      real (kind=R8) :: data(ns)

#ifdef MDSPLUS
#     include "mdslib.inc"
      integer :: status
      character :: tdi*(172)
      real(kind=R8), allocatable :: d1d(:)

      allocate(d1d(ns))
      d1d=data
      tdi='build_with_units($,"'//trim(units)//'")'
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2_d1d(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2_i1(idtype_double,ns,0),d1d,0)
#else
      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,ns,0),d1d,0)
#endif
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

      subroutine mds_put_double_us(node_name,units,n,data,dirsuffix)
      use b2mod_types
      use b2mod_mdsplus
      implicit none
      character*(*) :: node_name, units, dirsuffix
      integer :: n
      real (kind=R8) :: data(n)

#ifdef MDSPLUS
#     include "mdslib.inc"
      integer :: status
      character :: tdi*(172)
      real(kind=R8), allocatable :: d1d(:)

      allocate(d1d(n))
      d1d=data
      tdi='build_signal(build_with_units($,"' //
     1 trim(units) //
     2 '"),,\TOP.SNAPSHOT.GRID:CR' // trim(dirsuffix) //
     3 ',\TOP.SNAPSHOT.GRID:CZ' // trim(dirsuffix) //
     4 ')'

#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2_d1d(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2_i1(idtype_double,n,0),d1d,0)
#else
      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,n,0),d1d,0)
#endif
      write(*,'(a,a,a,a)')
     * 'mdsput2 of ',trim(node_name),' with ',trim(tdi)
      call mds_err_msg(status)
      deallocate(d1d)
#else
      call xerrab('mds_put_double_us: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_double_us

      subroutine mds_put_double_usns(node_name,units,n,ns,
     1 data,nsindex,dirsuffix)
      use b2mod_types
      use b2mod_mdsplus
      implicit none
      character*(*) :: node_name, units, nsindex, dirsuffix
      integer :: n, ns
      real (kind=R8) :: data(n,ns)

#ifdef MDSPLUS
#     include "mdslib.inc"
      integer :: status
      character :: tdi*(172)
      real (kind=R8), allocatable :: d2d(:,:)

      allocate(d2d(n,ns))
      d2d=data
      tdi='build_signal(build_with_units($,"' //
     1 trim(units) //
     2 '"),,\TOP.SNAPSHOT.GRID:CR' // trim(dirsuffix) //
     3 ',\TOP.SNAPSHOT.GRID:CZ' // trim(dirsuffix) //
     4 ',' //
     5 '\TOP.SNAPSHOT:' // trim(nsindex) // ')'

      write(*,'(a,a,a,a)')
     * 'mdsput2 of ',trim(node_name),' with ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2_d2d(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2_i2(idtype_double,n,ns,0),d2d,0)
#else
      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,n,ns,0),d2d,0)
#endif
      call mds_err_msg(status)
      deallocate(d2d)
#else
      call xerrab('mds_put_double_usns: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_double_usns

      subroutine mds_put_double_usdns(node_name,units,n,nd,ns,
     1 data,dirindex,nsindex,dirsuffix)
      use b2mod_types
      use b2mod_mdsplus
      implicit none
      character*(*) :: node_name, units, dirindex, nsindex, dirsuffix
      integer :: n, nd, ns
      real (kind=R8) :: data(n,nd,ns)

#ifdef MDSPLUS
#     include "mdslib.inc"
      integer :: status
      character :: tdi*(172)
      real (kind=R8), allocatable :: d3d(:,:,:)

      allocate(d3d(n,nd,ns))
      d3d=data
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
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2_d3d(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2_i3(idtype_double,n,nd,ns,0),d3d,0)
#else
      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,n,nd,ns,0),d3d,0)
#endif
      call mds_err_msg(status)
      deallocate(d3d)
#else
      call xerrab('mds_put_double_usdns: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_put_double_usdns

      subroutine mds_err_msg(status)
      use b2mod_mdsplus
      implicit none
      integer status
#ifdef MDSPLUS
#     include "mdslib.inc"
      integer :: size
      character :: message*80

      if(mod(status,2).eq.0) then
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2('getmsg($)'//char(0),
     1   descr2_i0(idtype_long,0),status,
     2   descr2_str0(idtype_cstring,0,80),message,0,size)
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
