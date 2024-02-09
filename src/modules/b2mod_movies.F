      module b2mod_movies
      use b2mod_types
      use b2mod_constants
#ifndef NO_CDF
      integer, save :: ncvdim, nsdim, natmdim, nmoldim, timedim, ntstep

      contains

      subroutine create_cdf_start(ncid, filename, nCv, ns, natm, nmol)
      implicit none
#     include <netcdf.inc>
      character*(*) filename
      integer ncid, nCv, ns, natm, nmol
* variable shapes
      integer dims(3)
* dimensions and variables
      integer iret
      integer ntstepid, timesid
* enter define mode
      iret = nf_create (filename, NF_WRITE, ncid)
      call check_cdf_status(iret)
* define dimensions
      iret = nf_def_dim(ncid, 'ncv', nCv, ncvdim)
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

      subroutine create_cdf_field2(ncid, name, long_name, units)
      implicit none
#     include <netcdf.inc>
      character*(*) name, long_name, units
* file id
      integer ncid
* variable shapes
      integer dims(2)
* dimensions and variables
      integer iret
      integer z_cvt_id
* define variable
      dims(1) = ncvdim
      dims(2) = timedim
      iret = nf_def_var (ncid, name, NF_DOUBLE, 2, dims, z_cvt_id)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_cvt_id, 'long_name',
     1 len_trim(long_name), trim(long_name))
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_cvt_id, 'units',
     1 len_trim(units), trim(units))
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_field2

      subroutine create_cdf_field3(ncid, name, itype, long_name, units)
      implicit none
#     include <netcdf.inc>
      integer ncid, itype
      character*(*) name, long_name, units
* variable shapes
      integer dims(3)
* dimensions and variables
      integer iret
      integer z_cvst_id
      external xerrab
* define variable
      dims(1) = ncvdim
      type: select case (itype)
      case(1)
        dims(2) = nsdim
      case(2)
        dims(2) = natmdim
      case(3)
        dims(2) = nmoldim
      case default type
        write(*,*) 'Undefined "itype" in create_cdf_field3 ',itype
        call xerrab ('Undefined itype value')
      end select type
      dims(3) = timedim
      iret = nf_def_var (ncid, name, NF_DOUBLE, 3, dims, z_cvst_id)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_cvst_id, 'long_name',
     1 len_trim(long_name), trim(long_name))
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_cvst_id, 'units',
     1 len_trim(units), trim(units))
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_field3

      subroutine create_cdf_end(ncid)
      implicit none
#     include <netcdf.inc>
* file id
      integer ncid
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
      write(*,*) 'create_cdf_end: ',ntstep
      iret=nf_put_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_end

      subroutine put_cdf_start(ncid, filename, time)
      implicit none
#     include <netcdf.inc>
      integer ncid
      character*(*) filename
      real (kind=R8) :: time(1)
      integer iret
      integer varid
      integer start(1), count(1), index(1)
c
      iret=nf_open(filename,or(NF_WRITE,NF_SHARE),ncid)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,'ntstep',varid)
      call check_cdf_status(iret)
      index(1)=1
      iret=nf_get_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
      ntstep=ntstep+1
      write(*,*) 'put_cdf_start: ',trim(filename),ntstep
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

      subroutine put_cdf_field2(ncid, name, z_cvt, nCv)
      implicit none
#     include <netcdf.inc>
      integer ncid
      character*(*) name
      integer nCv
      real (kind=R8) :: z_cvt(nCv)
      integer iret
      integer varid
      integer start(2), count(2)
c
      iret=nf_inq_varid(ncid,name,varid)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=ntstep
      count(1)=nCv
      count(2)=1
      iret=nf_put_vara_double(ncid,varid,start,count,z_cvt)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_field2

      subroutine put_cdf_field3(ncid, name, z_cvst, nCv, ns)
      implicit none
#     include <netcdf.inc>
      integer ncid
      character*(*) name
      integer nCv, ns
      real (kind=R8) :: z_cvst(nCv, ns)
      integer iret
      integer varid
      integer start(3), count(3)
c
      iret=nf_inq_varid(ncid,name,varid)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=1
      start(3)=ntstep
      count(1)=nCv
      count(2)=ns
      count(3)=1
      iret=nf_put_vara_double(ncid,varid,start,count,z_cvst)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_field3

      subroutine put_cdf_end(ncid)
      implicit none
#     include <netcdf.inc>
      integer ncid
      integer iret
c
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_end

      subroutine get_cdf_ntstep(filename, ntstep, ncid)
      implicit none
#     include <netcdf.inc>
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

      subroutine get_cdf(filename, ncid, ntstep, time, name, z_cvt, nCv)
      implicit none
#     include <netcdf.inc>
      character*(*) filename, name
      integer ntstep, nCv
      real (kind=R8) :: time(1), z_cvt(1:nCv)
      integer iret
      integer ncid, varid
      integer start(2), count(2)
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
      start(2)=ntstep
      count(1)=nCv
      count(2)=1
      iret=nf_get_vara_double(ncid,varid,start,count,z_cvt)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      end subroutine get_cdf

c OLD interface

      subroutine create_cdf(filename, ncid, nCv, name, long_name, units)
      implicit none
#     include <netcdf.inc>
      character*(*) filename, name, long_name, units
      integer nCv
      integer index(1), ntstep
      integer lnblnk
* variable shapes
      integer dims(2)
* dimensions and variables
      integer iret
      integer ncid, ncvdim, timedim
      integer ntstepid, timesid, z_cvt_id
* enter define mode
      iret = nf_create (filename, NF_NOCLOBBER, ncid)
      call check_cdf_status(iret)
* define dimensions
      iret = nf_def_dim(ncid, 'ncv', nCv, ncvdim)
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
      dims(1) = ncvdim
      dims(2) = timedim
      iret = nf_def_var (ncid, name, NF_DOUBLE, 2, dims, z_cvt_id)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_cvt_id, 'long_name',
     1 lnblnk(long_name), long_name)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_cvt_id, 'units',
     1 lnblnk(units), units)
      call check_cdf_status(iret)

* leave define mode
      iret = nf_enddef (ncid)
      call check_cdf_status(iret)
* set the initial timestep index (0)
      ntstep=0
      write(*,*) 'create_cdf: ',filename,ntstep
      iret=nf_put_var1_int(ncid,ntstepid,index,ntstep)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine create_cdf

      subroutine put_cdf(ncid, filename, time, name, z_cvt, nCv)
      implicit none
#     include <netcdf.inc>
      character*(*) filename, name
      integer nCv
      real (kind=R8) :: time(1), z_cvt(nCv)
      integer ntstep
      integer iret
      integer ncid, varid
      integer start(2), count(2), index(2)
c
      iret=nf_open(filename,or(NF_WRITE,NF_SHARE),ncid)
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
      start(2)=ntstep
      count(1)=nCv
      count(2)=1
      iret=nf_put_vara_double(ncid,varid,start,count,z_cvt)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf
#endif

      end module b2mod_movies
