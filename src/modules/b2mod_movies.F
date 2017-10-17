      module b2mod_movies
      use b2mod_types
      use b2mod_constants
#ifndef NO_CDF
#     include <netcdf.inc>
      integer, save :: ncid, nxdim, nydim, nsdim, natmdim, nmoldim, 
     & timedim, ntstep

      contains

      subroutine create_cdf_start(filename, nx, ny, ns, natm, nmol)
      implicit none
      character*(*) filename
      integer nx, ny, ns, natm, nmol
* variable shapes
      integer dims(3)
* dimensions and variables
      integer iret
      integer ntstepid, timesid
* enter define mode
      iret = nf_create (filename, NF_NOCLOBBER, ncid)
      call check_cdf_status(iret)
* define dimensions
      iret = nf_def_dim(ncid, 'nx', nx, nxdim)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'ny', ny, nydim)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'ns', ns, nsdim)
      call check_cdf_status(iret)
      if(natm.gt.0) then
        iret = nf_def_dim(ncid, 'natm', natm, natmdim)
      endif
      if(nmol.gt.0) then
        iret = nf_def_dim(ncid, 'nmol', nmol, nmoldim)
      endif
      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, timedim)
      call check_cdf_status(iret)
* define variables
      dims(1) = 1
      iret = nf_def_var (ncid, 'ntstep', NF_INT, 0, dims, ntstepid)
      call check_cdf_status(iret)
      dims(1) = timedim
      iret = nf_def_var (ncid, 'times', NF_DOUBLE, 1, dims, timesid)
      call check_cdf_status(iret)

      return
      end subroutine create_cdf_start

      subroutine create_cdf_field2(name, long_name, units)
      implicit none
      character*(*) name, long_name, units
* variable shapes
      integer dims(3)
* dimensions and variables
      integer iret
      integer z_xyt_id
* define variable
      dims(1) = nxdim
      dims(2) = nydim
      dims(3) = timedim
      iret = nf_def_var (ncid, name, NF_DOUBLE, 3, dims, z_xyt_id)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_xyt_id, 'long_name', 
     1 len_trim(long_name), trim(long_name))
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_xyt_id, 'units', 
     1 len_trim(units), trim(units))
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_field2

      subroutine create_cdf_field3(name, itype, long_name, units)
      implicit none
      integer itype
      character*(*) name, long_name, units
* variable shapes
      integer dims(4)
* dimensions and variables
      integer iret
      integer z_xyst_id
      external xerrab
* define variable
      dims(1) = nxdim
      dims(2) = nydim
      type: select case (itype)
      case(1)
        dims(3) = nsdim
      case(2)
        dims(3) = natmdim
      case(3)
        dims(3) = nmoldim
      case default type
        write(*,*) 'Undefined "itype" in create_cdf_field3 ',itype
        call xerrab ('Undefined itype value')
      end select type
      dims(4) = timedim
      iret = nf_def_var (ncid, name, NF_DOUBLE, 4, dims, z_xyst_id)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_xyst_id, 'long_name', 
     1 len_trim(long_name), trim(long_name))
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_xyst_id, 'units', 
     1 len_trim(units), trim(units))
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_field3

      subroutine create_cdf_end()
      implicit none
* dimensions and variables
      integer iret, varid, index(1)
* leave define mode
      iret = nf_enddef (ncid)
      call check_cdf_status(iret)
* set the initial timestep index (0)
      iret=nf_inq_varid(ncid,'ntstep',varid)
      call check_cdf_status(iret)
      index(1)=1
      ntstep=0
      write(*,*) 'put_cdf: ',ntstep
      iret=nf_put_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_end

      subroutine put_cdf_start(filename, time)
      implicit none
      character*(*) filename
      real (kind=R8) :: time(1)
      integer iret
      integer varid
      integer start(3), count(3), index(3)
c
      iret=nf_open(filename,NF_WRITE,ncid)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,'ntstep',varid)
      call check_cdf_status(iret)
      index(1)=1
      iret=nf_get_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
      ntstep=ntstep+1
      write(*,*) 'put_cdf: ',filename,ntstep
      iret=nf_put_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,'times',varid)
      call check_cdf_status(iret)
      start(1)=ntstep
      count(1)=1
      iret=nf_put_vara_double(ncid,varid,start,count,time)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_start

      subroutine put_cdf_field2(name, z_xyt, nx, ny)
      implicit none
      character*(*) name
      integer nx, ny
      real (kind=R8) :: z_xyt(nx, ny)
      integer iret
      integer varid
      integer start(3), count(3)
c
      iret=nf_inq_varid(ncid,name,varid)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=1
      start(3)=ntstep
      count(1)=nx
      count(2)=ny
      count(3)=1
      iret=nf_put_vara_double(ncid,varid,start,count,z_xyt)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_field2

      subroutine put_cdf_field3(name, z_xyst, nx, ny, ns)
      implicit none
      character*(*) name
      integer nx, ny, ns
      real (kind=R8) :: z_xyst(nx, ny, ns)
      integer iret
      integer varid
      integer start(4), count(4)
c
      iret=nf_inq_varid(ncid,name,varid)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=ntstep
      count(1)=nx
      count(2)=ny
      count(3)=ns
      count(4)=1
      iret=nf_put_vara_double(ncid,varid,start,count,z_xyst)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_field3

      subroutine put_cdf_end()
      implicit none
      integer iret
c
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_end

      subroutine get_cdf_ntstep(filename, ntstep)
      implicit none
      character*(*) filename
      integer ntstep
      integer iret
      integer ncid, varid
      integer index(1)
c
      iret=nf_open(filename,NF_NOWRITE,ncid)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,'ntstep',varid)
      call check_cdf_status(iret)
      index(1)=1
      iret=nf_get_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      end subroutine get_cdf_ntstep

      subroutine get_cdf(filename, ntstep, time, name, z_xyt, nx, ny)
      implicit none
      character*(*) filename, name
      integer ntstep, nx, ny
      real (kind=R8) :: time(1), z_xyt(-1:nx, -1:ny)
      integer iret
      integer ncid, varid
      integer start(3), count(3)
c
      iret=nf_open(filename,NF_NOWRITE,ncid)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,'times',varid)
      call check_cdf_status(iret)
      start(1)=ntstep
      count(1)=1
      iret=nf_get_vara_double(ncid,varid,start,count,time)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,name,varid)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=1
      start(3)=ntstep
      count(1)=nx+2
      count(2)=ny+2
      count(3)=1
      iret=nf_get_vara_double(ncid,varid,start,count,z_xyt)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      end subroutine get_cdf

c OLD interface

      subroutine create_cdf(filename, nx, ny, name, long_name, units)
      implicit none
      character*(*) filename, name, long_name, units
      integer nx, ny
      integer index(1), ntstep
      integer lnblnk
* variable shapes
      integer dims(3)
* dimensions and variables
      integer iret
      integer ncid, nxdim, nydim, timedim
      integer ntstepid, timesid, z_xyt_id
* enter define mode
      iret = nf_create (filename, NF_NOCLOBBER, ncid)
      call check_cdf_status(iret)
* define dimensions
      iret = nf_def_dim(ncid, 'nx', nx+2, nxdim)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'ny', ny+2, nydim)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, timedim)
      call check_cdf_status(iret)
* define variables
      dims(1) = 1
      iret = nf_def_var (ncid, 'ntstep', NF_INT, 0, dims, ntstepid)
      call check_cdf_status(iret)
      dims(1) = timedim
      iret = nf_def_var (ncid, 'times', NF_DOUBLE, 1, dims, timesid)
      call check_cdf_status(iret)
      dims(1) = nxdim
      dims(2) = nydim
      dims(3) = timedim
      iret = nf_def_var (ncid, name, NF_DOUBLE, 3, dims, z_xyt_id)
      call check_cdf_status(iret)
c      call ncaptc (ncid, z_xyt_id, 'long_name', NCCHAR,
c     1 lnblnk(long_name), long_name, iret)
c      call ncaptc(ncid, z_xyt_id, 'units', NCCHAR,
c     1 lnblnk(units), units, iret)
      iret = nf_put_att_text (ncid, z_xyt_id, 'long_name', 
     1 lnblnk(long_name), long_name)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_xyt_id, 'units', 
     1 lnblnk(units), units)
      call check_cdf_status(iret)

* leave define mode
      iret = nf_enddef (ncid)
      call check_cdf_status(iret)
* set the initial timestep index (0)
      ntstep=0
      write(*,*) 'put_cdf: ',filename,ntstep
      iret=nf_put_var1_int(ncid,ntstepid,index,ntstep)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine create_cdf

      subroutine put_cdf(filename, time, name, z_xyt, nx, ny)
      implicit none
      character*(*) filename, name
      integer nx, ny
      real (kind=R8) :: time(1), z_xyt(nx, ny)
      integer ntstep
      integer iret
      integer ncid, varid
      integer start(3), count(3), index(3)
c
      iret=nf_open(filename,NF_WRITE,ncid)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,'ntstep',varid)
      call check_cdf_status(iret)
      index(1)=1
      iret=nf_get_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
      ntstep=ntstep+1
      write(*,*) 'put_cdf: ',filename,ntstep
      iret=nf_put_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,'times',varid)
      call check_cdf_status(iret)
      start(1)=ntstep
      count(1)=1
      iret=nf_put_vara_double(ncid,varid,start,count,time)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,name,varid)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=1
      start(3)=ntstep
      count(1)=nx+2
      count(2)=ny+2
      count(3)=1
      iret=nf_put_vara_double(ncid,varid,start,count,z_xyt)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf
#endif

      end module b2mod_movies
