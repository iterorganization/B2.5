      module b2mod_chords
      use b2mod_types
      use b2mod_constants
#ifndef NO_CDF
c   ..netCDF-4 layout: shared dimensions/variables live in the root group
c   ..(ns/natm/nmol/time, ntstep, times, nchord, resolution); each chord
c   ..is a separate group 'chord_<k>' carrying its own npt dimension, its
c   ..coordinate variables (crx/cry/ds) and all the physics fields.
      integer, save :: nsdim, natmdim, nmoldim, timedim
      integer, save :: ntstepid, timesid, ntstep
      integer, save :: nptdim

      contains

c-----------------------------------------------------------------------
      subroutine chord_group_name(cnum, name)
c     Group name for chord number cnum (its original index in rzchords),
c     zero-padded so the groups sort numerically (chord_0001, ...).
      implicit none
      integer cnum
      character*(*) name
      write(name,'(a,i4.4)') 'chord_', cnum
      return
      end subroutine chord_group_name

c-----------------------------------------------------------------------
      subroutine create_cdf_chord_start(ncid, filename, ns, natm, nmol)
c     Create a netCDF-4 file and define the shared root dimensions and
c     variables. Each chord then gets its own group via
c     create_cdf_chord_group. Leaves the file in define mode.
      implicit none
#     include <netcdf.inc>
      character*(*) filename
      integer ncid, ns, natm, nmol
      integer dims(1)
      integer iret
      iret = nf_create(filename, or(NF_NETCDF4,NF_CLOBBER), ncid)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'ns', ns, nsdim)
      call check_cdf_status(iret)
      if (natm.gt.0) then
        iret = nf_def_dim(ncid, 'natm', natm, natmdim)
        call check_cdf_status(iret)
      endif
      if (nmol.gt.0) then
        iret = nf_def_dim(ncid, 'nmol', nmol, nmoldim)
        call check_cdf_status(iret)
      endif
      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, timedim)
      call check_cdf_status(iret)
      dims(1) = 1
      iret = nf_def_var(ncid, 'ntstep', NF_INT, 0, dims, ntstepid)
      call check_cdf_status(iret)
      dims(1) = timedim
      iret = nf_def_var(ncid, 'times', NF_DOUBLE, 1, dims, timesid)
      call check_cdf_status(iret)
      call def_coord_var(ncid, 'nchord', -1, NF_INT,
     1  'Number of chords (groups) in this file', '1')
      call def_coord_var(ncid, 'resolution', -1, NF_DOUBLE,
     1  'Chord sampling resolution', 'm')
      return
      end subroutine create_cdf_chord_start

c-----------------------------------------------------------------------
      subroutine create_cdf_chord_group(ncid, cnum, npt, grpid)
c     Define the group 'chord_<cnum>' with its own npt dimension and the
c     per-chord coordinate/metadata variables. Sets the module nptdim to
c     this group's npt dimension so the field helpers pick it up.
      implicit none
#     include <netcdf.inc>
      integer ncid, cnum, npt, grpid
      character*32 name
      integer iret
      call chord_group_name(cnum, name)
      iret = nf_def_grp(ncid, trim(name), grpid)
      call check_cdf_status(iret)
      iret = nf_def_dim(grpid, 'npt', npt, nptdim)
      call check_cdf_status(iret)
      call def_coord_var(grpid, 'crx', nptdim, NF_DOUBLE,
     1  'R coordinate of chord sample point', 'm')
      call def_coord_var(grpid, 'cry', nptdim, NF_DOUBLE,
     1  'Z coordinate of chord sample point', 'm')
      call def_coord_var(grpid, 'ds', nptdim, NF_DOUBLE,
     1  'Arc length along the chord from its start', 'm')
      call def_coord_var(grpid, 'chord_id', -1, NF_INT,
     1  'Original chord index (1..nchords in b2.user.parameters)', '1')
      call def_coord_var(grpid, 'rstart', -1, NF_DOUBLE,
     1  'R of chord start point', 'm')
      call def_coord_var(grpid, 'zstart', -1, NF_DOUBLE,
     1  'Z of chord start point', 'm')
      call def_coord_var(grpid, 'rend', -1, NF_DOUBLE,
     1  'R of chord end point', 'm')
      call def_coord_var(grpid, 'zend', -1, NF_DOUBLE,
     1  'Z of chord end point', 'm')
      return
      end subroutine create_cdf_chord_group

c-----------------------------------------------------------------------
      subroutine create_cdf_chord_field2(ncid, name, long_name, units)
c     Define a scalar (per-species-independent) field in group ncid.
      implicit none
#     include <netcdf.inc>
      character*(*) name, long_name, units
      integer ncid
      integer dims(2)
      integer iret
      integer z_id
      dims(1) = nptdim
      dims(2) = timedim
      iret = nf_def_var (ncid, name, NF_DOUBLE, 2, dims, z_id)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_id, 'long_name',
     1 len_trim(long_name), trim(long_name))
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_id, 'units',
     1 len_trim(units), trim(units))
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_chord_field2

c-----------------------------------------------------------------------
      subroutine create_cdf_chord_field3(ncid, name, itype, long_name,
     &      units)
c     Define a per-species field in group ncid. itype selects the species
c     dimension: 1=ns, 2=natm, 3=nmol (all inherited from the root group).
      implicit none
#     include <netcdf.inc>
      integer ncid, itype
      character*(*) name, long_name, units
      integer dims(3)
      integer iret
      integer z_id
      external xerrab
      dims(1) = nptdim
      select case (itype)
      case(1)
        dims(2) = nsdim
      case(2)
        dims(2) = natmdim
      case(3)
        dims(2) = nmoldim
      case default
        write(*,*) 'Undefined "itype" in create_cdf_chord_field3 ',itype
        call xerrab ('Undefined itype value')
      end select
      dims(3) = timedim
      iret = nf_def_var (ncid, name, NF_DOUBLE, 3, dims, z_id)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_id, 'long_name',
     1 len_trim(long_name), trim(long_name))
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_id, 'units',
     1 len_trim(units), trim(units))
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_chord_field3

c-----------------------------------------------------------------------
      subroutine create_cdf_chord_end(ncid, nchord, resolution)
c     Leave define mode, initialise ntstep=0 and write the root scalars
c     nchord and resolution, then close the file.
      implicit none
#     include <netcdf.inc>
      integer ncid, nchord
      real (kind=R8) :: resolution
      integer iret, varid, index(1)
      real (kind=R8) :: rtmp(1)
      iret = nf_enddef (ncid)
      call check_cdf_status(iret)
      index(1)=1
      ntstep=0
      iret=nf_inq_varid(ncid,'ntstep',varid)
      call check_cdf_status(iret)
      iret=nf_put_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
      iret=nf_inq_varid(ncid,'nchord',varid)
      call check_cdf_status(iret)
      iret=nf_put_var1_int(ncid,varid,index,nchord)
      call check_cdf_status(iret)
      rtmp(1)=resolution
      iret=nf_inq_varid(ncid,'resolution',varid)
      call check_cdf_status(iret)
      iret=nf_put_var_double(ncid,varid,rtmp)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_chord_end

c-----------------------------------------------------------------------
      subroutine put_cdf_chord_coords(filename, nchord, chordid, off,
     & cnt, nptt, crx, cry, ds, rstart, zstart, rend, zend)
c     Reopen the file and write, into each chord's group, the static
c     coordinate variables (crx/cry/ds) and the metadata scalars
c     (chord_id, rstart/zstart/rend/zend). crx/cry/ds are the
c     concatenated arrays; chord m occupies points off(m)+1..off(m)+cnt(m)
      implicit none
#     include <netcdf.inc>
      character*(*) filename
      integer nchord, nptt
      integer chordid(nchord), off(nchord), cnt(nchord)
      real (kind=R8) :: crx(nptt), cry(nptt), ds(nptt)
      real (kind=R8) :: rstart(nchord), zstart(nchord), rend(nchord),
     & zend(nchord)
      integer iret, ncid, grpid, m, varid, index(1)
      character*32 name
      real (kind=R8) :: rtmp(1)
      iret = nf_open(filename, NF_WRITE, ncid)
      call check_cdf_status(iret)
      index(1)=1
      do m = 1, nchord
        call chord_group_name(chordid(m), name)
        iret = nf_inq_grp_ncid(ncid, trim(name), grpid)
        call check_cdf_status(iret)
        call put_static_slice(grpid, 'crx', crx, off(m), cnt(m))
        call put_static_slice(grpid, 'cry', cry, off(m), cnt(m))
        call put_static_slice(grpid, 'ds',  ds,  off(m), cnt(m))
        iret = nf_inq_varid(grpid,'chord_id',varid)
        call check_cdf_status(iret)
        iret = nf_put_var1_int(grpid,varid,index,chordid(m))
        call check_cdf_status(iret)
        rtmp(1)=rstart(m)
        call put_scalar_double(grpid,'rstart',rtmp)
        rtmp(1)=zstart(m)
        call put_scalar_double(grpid,'zstart',rtmp)
        rtmp(1)=rend(m)
        call put_scalar_double(grpid,'rend',rtmp)
        rtmp(1)=zend(m)
        call put_scalar_double(grpid,'zend',rtmp)
      enddo
      iret = nf_close(ncid)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_chord_coords

c-----------------------------------------------------------------------
      subroutine put_cdf_chord_start(ncid, filename, time)
c     Reopen the file, advance and store ntstep and append the time value
c     (both live in the root group). Leaves the file open; the caller
c     then fetches the group ids and writes the frame fields.
      implicit none
#     include <netcdf.inc>
      integer ncid
      character*(*) filename
      real (kind=R8) :: time(1)
      integer iret
      integer varid
      integer start(1), count(1), index(1)
      iret=nf_open(filename,or(NF_WRITE,NF_SHARE),ncid)
      call check_cdf_status(iret)
      iret=nf_inq_varid(ncid,'ntstep',varid)
      call check_cdf_status(iret)
      index(1)=1
      iret=nf_get_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
      ntstep=ntstep+1
      iret=nf_put_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
      iret=nf_inq_varid(ncid,'times',varid)
      call check_cdf_status(iret)
      start(1)=ntstep
      count(1)=1
      iret=nf_put_vara_double(ncid,varid,start,count,time)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_chord_start

c-----------------------------------------------------------------------
      subroutine get_cdf_chord_group(ncid, cnum, grpid)
c     Return the group id of chord number cnum in the open file ncid.
      implicit none
#     include <netcdf.inc>
      integer ncid, cnum, grpid, iret
      character*32 name
      call chord_group_name(cnum, name)
      iret = nf_inq_grp_ncid(ncid, trim(name), grpid)
      call check_cdf_status(iret)
      return
      end subroutine get_cdf_chord_group

c-----------------------------------------------------------------------
      subroutine put_cdf_chord_field2(ncid, name, z, npt)
c     Write a scalar field frame into group ncid at the current ntstep.
      implicit none
#     include <netcdf.inc>
      integer ncid
      character*(*) name
      integer npt
      real (kind=R8) :: z(npt)
      integer iret
      integer varid
      integer start(2), count(2)
      iret=nf_inq_varid(ncid,name,varid)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=ntstep
      count(1)=npt
      count(2)=1
      iret=nf_put_vara_double(ncid,varid,start,count,z)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_chord_field2

c-----------------------------------------------------------------------
      subroutine put_cdf_chord_field3(ncid, name, z, npt, ns)
c     Write a per-species field frame into group ncid at the current
c     ntstep.
      implicit none
#     include <netcdf.inc>
      integer ncid
      character*(*) name
      integer npt, ns
      real (kind=R8) :: z(npt, ns)
      integer iret
      integer varid
      integer start(3), count(3)
      iret=nf_inq_varid(ncid,name,varid)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=1
      start(3)=ntstep
      count(1)=npt
      count(2)=ns
      count(3)=1
      iret=nf_put_vara_double(ncid,varid,start,count,z)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_chord_field3

c-----------------------------------------------------------------------
      subroutine put_cdf_chord_end(ncid)
      implicit none
#     include <netcdf.inc>
      integer ncid
      integer iret
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_chord_end

c-----------------------------------------------------------------------
      subroutine def_coord_var(ncid, name, dimid, xtype, long_name,
     & units)
c     Define a static (time-independent) coordinate/metadata variable.
c     dimid <= 0 defines a scalar; otherwise a rank-1 variable over dimid.
      implicit none
#     include <netcdf.inc>
      integer ncid, dimid, xtype
      character*(*) name, long_name, units
      integer iret, vid, dims(1)
      if (dimid.le.0) then
        iret = nf_def_var(ncid, name, xtype, 0, dims, vid)
      else
        dims(1) = dimid
        iret = nf_def_var(ncid, name, xtype, 1, dims, vid)
      endif
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vid, 'long_name',
     1  len_trim(long_name), trim(long_name))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vid, 'units',
     1  len_trim(units), trim(units))
      call check_cdf_status(iret)
      return
      end subroutine def_coord_var

c-----------------------------------------------------------------------
      subroutine put_static_slice(ncid, name, arr, off, n)
c     Write n values arr(off+1..off+n) into the rank-1 variable name.
      implicit none
#     include <netcdf.inc>
      integer ncid, off, n
      character*(*) name
      real (kind=R8) :: arr(*)
      integer iret, varid, start(1), count(1)
      iret = nf_inq_varid(ncid, name, varid)
      call check_cdf_status(iret)
      start(1)=1
      count(1)=n
      iret = nf_put_vara_double(ncid, varid, start, count, arr(off+1))
      call check_cdf_status(iret)
      return
      end subroutine put_static_slice

c-----------------------------------------------------------------------
      subroutine put_scalar_double(ncid, name, val)
      implicit none
#     include <netcdf.inc>
      integer ncid
      character*(*) name
      real (kind=R8) :: val(1)
      integer iret, varid
      iret = nf_inq_varid(ncid, name, varid)
      call check_cdf_status(iret)
      iret = nf_put_var_double(ncid, varid, val)
      call check_cdf_status(iret)
      return
      end subroutine put_scalar_double

c-----------------------------------------------------------------------
      subroutine get_cdf_chord_meta(filename, nmax, ok, nchord, chordid,
     & npt, rstart, zstart, rend, zend, resolution)
c     Read back the per-chord metadata (one group per chord) of an
c     existing file, to decide whether an append is compatible with the
c     current chords. Groups are returned in their stored order. ok is
c     .false. if the file is missing or lacks anything needed (then it is
c     rewritten instead of appended).
      implicit none
#     include <netcdf.inc>
      character*(*) filename
      integer nmax, nchord
      integer chordid(nmax), npt(nmax)
      real (kind=R8) :: rstart(nmax), zstart(nmax), rend(nmax),
     & zend(nmax), resolution
      logical ok
      integer status, ncid, varid, dimid, ng, i
      integer grpids(nmax)
      real (kind=R8) :: rtmp(1)
      ok = .false.
      nchord = 0
      resolution = 0.0_R8
      status = nf_open(filename, NF_NOWRITE, ncid)
      if (status.ne.NF_NOERR) return
      status = nf_inq_grps(ncid, ng, grpids)
      if (status.ne.NF_NOERR) goto 900
      if (ng.lt.1 .or. ng.gt.nmax) goto 900
      nchord = ng
      status = nf_inq_varid(ncid,'resolution',varid)
      if (status.ne.NF_NOERR) goto 900
      status = nf_get_var_double(ncid,varid,rtmp)
      if (status.ne.NF_NOERR) goto 900
      resolution = rtmp(1)
      do i = 1, ng
        status = nf_inq_dimid(grpids(i),'npt',dimid)
        if (status.ne.NF_NOERR) goto 900
        status = nf_inq_dimlen(grpids(i),dimid,npt(i))
        if (status.ne.NF_NOERR) goto 900
        status = nf_inq_varid(grpids(i),'chord_id',varid)
        if (status.ne.NF_NOERR) goto 900
        status = nf_get_var_int(grpids(i),varid,chordid(i))
        if (status.ne.NF_NOERR) goto 900
        status = nf_inq_varid(grpids(i),'rstart',varid)
        if (status.ne.NF_NOERR) goto 900
        status = nf_get_var_double(grpids(i),varid,rtmp)
        if (status.ne.NF_NOERR) goto 900
        rstart(i)=rtmp(1)
        status = nf_inq_varid(grpids(i),'zstart',varid)
        if (status.ne.NF_NOERR) goto 900
        status = nf_get_var_double(grpids(i),varid,rtmp)
        if (status.ne.NF_NOERR) goto 900
        zstart(i)=rtmp(1)
        status = nf_inq_varid(grpids(i),'rend',varid)
        if (status.ne.NF_NOERR) goto 900
        status = nf_get_var_double(grpids(i),varid,rtmp)
        if (status.ne.NF_NOERR) goto 900
        rend(i)=rtmp(1)
        status = nf_inq_varid(grpids(i),'zend',varid)
        if (status.ne.NF_NOERR) goto 900
        status = nf_get_var_double(grpids(i),varid,rtmp)
        if (status.ne.NF_NOERR) goto 900
        zend(i)=rtmp(1)
      enddo
      ok = .true.
 900  status = nf_close(ncid)
      return
      end subroutine get_cdf_chord_meta
#endif

      end module b2mod_chords
