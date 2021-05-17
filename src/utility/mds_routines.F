      subroutine mds_put_string(node_name,string)
      implicit none
      character*(*) node_name,string
#ifdef MDSPLUS
#     include "mdslib.inc"
      integer :: status

      if(len_trim(string).gt.0) then
        write(*,'(a,a,a,a,a)') 'mds_put_string to ',trim(node_name),
     1   ' of "',trim(string),'"'
        status=mdsput2(trim(node_name)//char(0),'$'//char(0),
     1   descr2(idtype_cstring,0,len_trim(string)),
     1   trim(string),0)
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
#     include "mdslib.inc"
      integer is
      integer :: status
      character*18 hlp_frm

      write(hlp_frm,'(a,i3,a)') '(a,a,a,',ns,'(a,a),a)'
      write(*,hlp_frm) 'mds_put_string to ',trim(node_name),
     1 ' of "',(' ',trim(string(is)),is=1,ns),'"'
      do is=1,ns
        string(is)=trim(string(is)) // char(0)
      enddo
      status=mdsput2(trim(node_name)//char(0),'$'//char(0),
     1 descr2(idtype_cstring,ns,0,len(string(1))),string,0)
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
#     include "mdslib.inc"
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
#     include "mdslib.inc"
      integer :: status
      character :: tdi*(172)

      tdi='$'

      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_long,nx,ny,0),data,0)
      write(*,'(a,a,a,a)')
     * 'mdsput2 of ',trim(node_name),' with ',trim(tdi)
      call mds_err_msg(status)

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
#     include "mdslib.inc"
      integer :: status
      character :: tdi*(172)

      tdi='$'

      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_long,ns,0),data,0)
      write(*,'(a,a,a,a)')
     * 'mdsput2 of ',trim(node_name),' with ',trim(tdi)
      call mds_err_msg(status)

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
#     include "mdslib.inc"
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
#     include "mdslib.inc"
      integer :: status
      character :: tdi*(172)

      tdi='build_with_units($,"'//trim(units)//'")'
      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,ns,0),data,0)
      write(*,'(a,a,a,a)')
     * 'mdsput2 of ',trim(node_name),' with ',trim(tdi)
      call mds_err_msg(status)

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
#     include "mdslib.inc"
      integer :: status
      character :: tdi*(172)

      tdi='build_signal(build_with_units($,"' //
     1 trim(units) //
     2 '"),,\TOP.SNAPSHOT.GRID:CR' // trim(dirsuffix) //
     3 ',\TOP.SNAPSHOT.GRID:CZ' // trim(dirsuffix) //
     4 ')'

      status=mdsput2(
     * trim(node_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nx,ny,0),data,0)
      write(*,'(a,a,a,a)')
     * 'mdsput2 of ',trim(node_name),' with ',trim(tdi)
      call mds_err_msg(status)

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
#     include "mdslib.inc"
      integer :: status
      character :: tdi*(172)

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
     * descr2(idtype_double,nx,ny,ns,0),data,0)
      call mds_err_msg(status)

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
#     include "mdslib.inc"
      integer :: status
      character :: tdi*(172)

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
     * descr2(idtype_double,nx,ny,nd,ns,0),data,0)
      call mds_err_msg(status)

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
#     include "mdslib.inc"
      integer :: size
      character :: message*80

      if(mod(status,2).eq.0) then
        status=mdsvalue2('getmsg($)'//char(0),
     1   descr2(idtype_long,0),status,
     2   descr2(idtype_cstring,0,80),message,0,size)
        write(*,*) trim(message)
      endif
#else
      call xerrab('mds_err_msg: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      return
      end subroutine mds_err_msg
