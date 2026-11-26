      module b2mod_chords
      use b2mod_types
      use b2mod_constants
#ifndef NO_CDF
c   ..netCDF-4 layout: shared dimensions/variables live in the root group
c   ..(ns/natm/nmol/time, ntstep, times, nchord); each chord set (one
c   ..chord file) is a group named after its label, and each chord of
c   ..the set is a group <set>_<k> inside it, carrying its own npt
c   ..dimension, its coordinate variables (R/z/phi/ds) and all the
c   ..physics fields.
      integer, save :: nsdim, natmdim, nmoldim, timedim
      integer, save :: ntstepid, timesid, ntstep
      integer, save :: nptdim
c   ..largest number of chord files and of chords per file
      integer, parameter :: maxcgrp = 10000

      contains

c-----------------------------------------------------------------------
      subroutine chord_group_name(setname, cnum, name)
c     Group name of chord number cnum of the set setname, zero-padded
c     so that the groups sort numerically (<set>_0001, ...).
      implicit none
      character*(*) setname, name
      integer cnum
      if (cnum.lt.10000) then
        write(name,'(a,a,i4.4)') trim(setname), '_', cnum
      else
        write(name,'(a,a,i5)') trim(setname), '_', cnum
      endif
      return
      end subroutine chord_group_name

c-----------------------------------------------------------------------
      subroutine chord_set_name(label, name)
c     Group name of a chord set: its label with every character other
c     than a letter, a digit, _ - + or . (or a leading - + .) replaced
c     by an underscore.
      implicit none
      character*(*) label, name
      character c
      integer i
      name = adjustl(label)
      do i = 1, len_trim(name)
        c = name(i:i)
        if (.not.((lge(c,'a').and.lle(c,'z')) .or.
     &      (lge(c,'A').and.lle(c,'Z')) .or.
     &      (lge(c,'0').and.lle(c,'9')) .or.
     &      index('_-+.',c).gt.0)) name(i:i) = '_'
      enddo
      if (index('-+.',name(1:1)).gt.0) name(1:1) = '_'
      return
      end subroutine chord_set_name

c-----------------------------------------------------------------------
      subroutine read_chord_files(nset, setname, setlabel, setnp,
     & nch, cset, cnum, cxyz, cval)
c     Read the chords of every file *.chords in the working directory
c     and in the run directory (../), in the format of the b2plot chord
c     files: a first line with the label (in quotes) and npoints, then
c     one line per chord with x1 y1 z1 x2 y2 z2 value (Cartesian end
c     points in m, z vertical). Each file is a chord set; chord k is
c     line cnum(k) of set cset(k). A file with a format error is
c     skipped with a warning.
      implicit none
      integer nset, nch
      character*80, allocatable :: setname(:), setlabel(:)
      integer, allocatable :: setnp(:), cset(:), cnum(:)
      real (kind=R8), allocatable :: cxyz(:,:,:), cval(:)
      character*256 lstname, fname
      character*1024 line
      character*256, allocatable :: files(:), setfile(:)
      character*80 lbl, sname
      integer nfile, lu, ios, ios2, i, j, k, np, n, iline
      logical ex, ok, dup
      real (kind=R8) :: v(7)
      real (kind=R8), allocatable :: buf(:,:)
      external xerrab
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif

c   ..list the chord files, those of the working directory first
      lstname = 'cdfchord.lst'
#ifndef LEGACYCOMP
      call execute_command_line('(ls -1 *.chords; ls -1 ../*.chords) > '
     &  //trim(lstname)//' 2> /dev/null', exitstat=ios, cmdstat=ios2)
      if (ios2.ne.0) write(*,*) 'cdfchord: cannot list the chord files'
#else
      call system('(ls -1 *.chords; ls -1 ../*.chords) > '
     &  //trim(lstname)//' 2> /dev/null')
#endif
      allocate(files(maxcgrp))
      nfile = 0
      inquire(file=lstname, exist=ex)
      ios = 1
#ifndef LEGACYCOMP
      if (ex) open(newunit=lu, file=lstname, status='old', iostat=ios)
#else
      if (ex) then
        lu = newunit()
        open(unit=lu, file=lstname, status='old', iostat=ios)
      endif
#endif
      if (ex .and. ios.eq.0) then
        read(lu,'(a)',iostat=ios) fname
        do while (ios.eq.0)
c   ..a file of the run directory hidden by one of the same name here
          dup = .false.
          if (fname(1:3).eq.'../') then
            do i = 1, nfile
              if (files(i).eq.fname(4:)) dup = .true.
            enddo
          endif
          if (len_trim(fname).gt.0 .and. .not.dup) then
            if (nfile.lt.maxcgrp) then
              nfile = nfile + 1
              files(nfile) = fname
            else
              write(*,*) 'cdfchord: too many chord files, ignored ',
     &          trim(fname)
            endif
          endif
          read(lu,'(a)',iostat=ios) fname
        enddo
        close(lu, status='delete')
      endif

      allocate(setname(max(nfile,1)), setlabel(max(nfile,1)))
      allocate(setnp(max(nfile,1)), setfile(max(nfile,1)))
      allocate(cset(0), cnum(0), cxyz(3,2,0), cval(0))
      nset = 0
      nch = 0
      do i = 1, nfile
        call read_one_file(files(i), ok)
        if (.not.ok) cycle
        call chord_set_name(lbl, sname)
        if (len_trim(sname).eq.0) then
          write(*,*) 'cdfchord: empty label, chord file skipped: ',
     &      trim(files(i))
          cycle
        endif
        do j = 1, nset
          if (setname(j).eq.sname) then
            write(*,*) 'cdfchord: chord files ', trim(setfile(j)),
     &        ' and ', trim(files(i)), ' have the same label ',
     &        trim(sname)
            call xerrab ('two chord files with the same label')
          endif
        enddo
        nset = nset + 1
        setname(nset) = sname
        setlabel(nset) = lbl
        setnp(nset) = np
        setfile(nset) = files(i)
        call grow_chords(nch + n)
        do k = 1, n
          cset(nch+k) = nset
          cnum(nch+k) = k
          cxyz(1:3,1,nch+k) = buf(1:3,k)
          cxyz(1:3,2,nch+k) = buf(4:6,k)
          cval(nch+k) = buf(7,k)
        enddo
        nch = nch + n
        write(*,'(a,i6,a,a,a,a)') ' cdfchord: ', n,
     &    ' chords of set ', trim(sname), ' read from ', trim(files(i))
      enddo
      if (allocated(buf)) deallocate(buf)
      deallocate(files, setfile)
      return

      contains

      subroutine read_one_file(fn, good)
c     Read one chord file into lbl, np (npoints), n and buf(1:7,1:n).
      implicit none
      character*(*) fn
      logical good
      integer lf, ierr
      logical bad
      good = .false.
#ifndef LEGACYCOMP
      open(newunit=lf, file=fn, status='old', iostat=ierr)
#else
      lf = newunit()
      open(unit=lf, file=fn, status='old', iostat=ierr)
#endif
      if (ierr.ne.0) then
        write(*,*) 'cdfchord: cannot open chord file ', trim(fn)
        return
      endif
      if (.not.allocated(buf)) allocate(buf(7,1000))
      n = 0
      np = 0
      iline = 0
      bad = .false.
      read(lf,'(a)',iostat=ierr) line
      do while (ierr.eq.0 .and. .not.bad .and. n.le.maxcgrp)
        iline = iline + 1
        if (len_trim(line).gt.0) then
          if (np.eq.0) then
            read(line,*,iostat=ierr) lbl, np
            bad = ierr.ne.0 .or. np.eq.0
            np = abs(np)
          else
            read(line,*,iostat=ierr) v
            bad = ierr.ne.0
            if (.not.bad) then
              if (n.ge.size(buf,2)) call grow_buf(2*size(buf,2))
              n = n + 1
              buf(1:7,n) = v
            endif
          endif
        endif
        if (.not.bad) read(lf,'(a)',iostat=ierr) line
      enddo
      close(lf)
      if (bad) then
        write(*,'(a,a,a,i6,a)') ' cdfchord: format error in ',
     &    trim(fn), ' (line ', iline, '), chord file skipped'
      else if (n.gt.maxcgrp) then
        write(*,*) 'cdfchord: more than ', maxcgrp,
     &    ' chords, chord file skipped: ', trim(fn)
      else if (n.eq.0) then
        write(*,*) 'cdfchord: no chords, chord file skipped: ',trim(fn)
      else
        good = .true.
      endif
      return
      end subroutine read_one_file

      subroutine grow_buf(m)
c     Enlarge buf to m chords, keeping its contents.
      implicit none
      integer m
      real (kind=R8), allocatable :: tmp(:,:)
      allocate(tmp(7,m))
      tmp(:,1:size(buf,2)) = buf
      call move_alloc(tmp, buf)
      return
      end subroutine grow_buf

      subroutine grow_chords(m)
c     Enlarge the chord arrays to m chords, keeping their contents.
      implicit none
      integer m
      integer, allocatable :: itmp(:)
      real (kind=R8), allocatable :: rtmp(:), xtmp(:,:,:)
      allocate(itmp(m))
      itmp(1:nch) = cset(1:nch)
      call move_alloc(itmp, cset)
      allocate(itmp(m))
      itmp(1:nch) = cnum(1:nch)
      call move_alloc(itmp, cnum)
      allocate(rtmp(m))
      rtmp(1:nch) = cval(1:nch)
      call move_alloc(rtmp, cval)
      allocate(xtmp(3,2,m))
      xtmp(:,:,1:nch) = cxyz(:,:,1:nch)
      call move_alloc(xtmp, cxyz)
      return
      end subroutine grow_chords

      end subroutine read_chord_files

c-----------------------------------------------------------------------
      subroutine create_cdf_chord_start(ncid, filename, ns, natm, nmol)
c     Create a netCDF-4 file and define the shared root dimensions and
c     variables. Each chord set then gets its own group via
c     create_cdf_chord_set and each chord a group inside it via
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
      return
      end subroutine create_cdf_chord_start

c-----------------------------------------------------------------------
      subroutine create_cdf_chord_set(ncid, setname, label, npoints,
     & sgrpid)
c     Define the group of a chord set, with the attributes label (the
c     label of its chord file) and npoints (samples per chord).
      implicit none
#     include <netcdf.inc>
      integer ncid, npoints, sgrpid
      character*(*) setname, label
      integer iret, ival(1)
      iret = nf_def_grp(ncid, trim(setname), sgrpid)
      call check_cdf_status(iret)
      iret = nf_put_att_text(sgrpid, NF_GLOBAL, 'label',
     1  len_trim(label), trim(label))
      call check_cdf_status(iret)
      ival(1) = npoints
      iret = nf_put_att_int(sgrpid, NF_GLOBAL, 'npoints', NF_INT, 1,
     1  ival)
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_chord_set

c-----------------------------------------------------------------------
      subroutine create_cdf_chord_group(sgrpid, setname, cnum, npt,
     & grpid)
c     Define the group <set>_<cnum> inside the set group sgrpid with
c     its own npt dimension and the per-chord coordinate/metadata
c     variables. Sets the module nptdim to this group npt dimension so
c     the field helpers pick it up.
      implicit none
#     include <netcdf.inc>
      integer sgrpid, cnum, npt, grpid
      character*(*) setname
      character*(NF_MAX_NAME) name
      integer iret
      call chord_group_name(setname, cnum, name)
      iret = nf_def_grp(sgrpid, trim(name), grpid)
      call check_cdf_status(iret)
      iret = nf_def_dim(grpid, 'npt', npt, nptdim)
      call check_cdf_status(iret)
      call def_coord_var(grpid, 'R', nptdim, NF_DOUBLE,
     1  'R coordinate of chord sample point', 'm')
      call def_coord_var(grpid, 'z', nptdim, NF_DOUBLE,
     1  'Z coordinate of chord sample point', 'm')
      call def_coord_var(grpid, 'phi', nptdim, NF_DOUBLE,
     1  'Toroidal angle of chord sample point', 'deg')
      call def_coord_var(grpid, 'ds', nptdim, NF_DOUBLE,
     1  '3D arc length along the chord from its start', 'm')
      call def_coord_var(grpid, 'chord_id', -1, NF_INT,
     1  'Chord number (line of the chord in its chord file)', '1')
      call def_coord_var(grpid, 'chord_value', -1, NF_DOUBLE,
     1  'Chord value (last column of the chord file)', '1')
      call def_coord_var(grpid, 'Rstart', -1, NF_DOUBLE,
     1  'R of chord start point', 'm')
      call def_coord_var(grpid, 'zstart', -1, NF_DOUBLE,
     1  'Z of chord start point', 'm')
      call def_coord_var(grpid, 'phistart', -1, NF_DOUBLE,
     1  'Toroidal angle of chord start point', 'deg')
      call def_coord_var(grpid, 'Rend', -1, NF_DOUBLE,
     1  'R of chord end point', 'm')
      call def_coord_var(grpid, 'zend', -1, NF_DOUBLE,
     1  'Z of chord end point', 'm')
      call def_coord_var(grpid, 'phiend', -1, NF_DOUBLE,
     1  'Toroidal angle of chord end point', 'deg')
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
      subroutine create_cdf_chord_end(ncid, nchord)
c     Leave define mode, initialise ntstep=0 and write the root scalar
c     nchord, then close the file.
      implicit none
#     include <netcdf.inc>
      integer ncid, nchord
      integer iret, varid, index(1)
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
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_chord_end

c-----------------------------------------------------------------------
      subroutine put_cdf_chord_coords(filename, nchord, setname, cnum,
     & cval, off, cnt, nptt, crx, cry, cphi, ds, rstart, zstart,
     & phistart, rend, zend, phiend)
c     Reopen the file and write, into each chord group, the static
c     coordinate variables (R/z/phi/ds) and the metadata scalars
c     (chord_id, chord_value, Rstart/zstart/phistart/Rend/zend/phiend).
c     Chord m is number cnum(m) of set setname(m); crx/cry/cphi/ds are
c     the concatenated arrays, chord m occupying points
c     off(m)+1..off(m)+cnt(m).
      implicit none
#     include <netcdf.inc>
      character*(*) filename
      integer nchord, nptt
      character*(*) setname(nchord)
      integer cnum(nchord), off(nchord), cnt(nchord)
      real (kind=R8) :: cval(nchord)
      real (kind=R8) :: crx(nptt), cry(nptt), cphi(nptt), ds(nptt)
      real (kind=R8) :: rstart(nchord), zstart(nchord),
     & phistart(nchord), rend(nchord), zend(nchord), phiend(nchord)
      integer iret, ncid, grpid, m, varid, index(1)
      real (kind=R8) :: rtmp(1)
      iret = nf_open(filename, NF_WRITE, ncid)
      call check_cdf_status(iret)
      index(1)=1
      do m = 1, nchord
        call get_cdf_chord_group(ncid, setname(m), cnum(m), grpid)
        call put_static_slice(grpid, 'R', crx, off(m), cnt(m))
        call put_static_slice(grpid, 'z', cry, off(m), cnt(m))
        call put_static_slice(grpid, 'phi', cphi, off(m), cnt(m))
        call put_static_slice(grpid, 'ds',  ds,  off(m), cnt(m))
        iret = nf_inq_varid(grpid,'chord_id',varid)
        call check_cdf_status(iret)
        iret = nf_put_var1_int(grpid,varid,index,cnum(m))
        call check_cdf_status(iret)
        rtmp(1)=cval(m)
        call put_scalar_double(grpid,'chord_value',rtmp)
        rtmp(1)=rstart(m)
        call put_scalar_double(grpid,'Rstart',rtmp)
        rtmp(1)=zstart(m)
        call put_scalar_double(grpid,'zstart',rtmp)
        rtmp(1)=phistart(m)
        call put_scalar_double(grpid,'phistart',rtmp)
        rtmp(1)=rend(m)
        call put_scalar_double(grpid,'Rend',rtmp)
        rtmp(1)=zend(m)
        call put_scalar_double(grpid,'zend',rtmp)
        rtmp(1)=phiend(m)
        call put_scalar_double(grpid,'phiend',rtmp)
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
      subroutine get_cdf_chord_group(ncid, setname, cnum, grpid)
c     Return the group id of chord number cnum of set setname in the
c     open file ncid.
      implicit none
#     include <netcdf.inc>
      integer ncid, cnum, grpid, iret, sgrpid
      character*(*) setname
      character*(NF_MAX_NAME) name
      iret = nf_inq_grp_ncid(ncid, trim(setname), sgrpid)
      call check_cdf_status(iret)
      call chord_group_name(setname, cnum, name)
      iret = nf_inq_grp_ncid(sgrpid, trim(name), grpid)
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
      subroutine get_cdf_chord_meta(filename, ok, nchord, setname,
     & cnum, npt, npoints, cval, rstart, zstart, phistart, rend, zend,
     & phiend)
c     Read back the chord groups of an existing file (chord sets in the
c     root group, chords in each set), to decide whether an append is
c     compatible with the current chords. Chord m is number cnum(m) of
c     set setname(m). ok is .false. if the file is missing or lacks
c     anything needed, e.g. a file of an older layout (then it is
c     rewritten instead of appended).
      implicit none
#     include <netcdf.inc>
      character*(*) filename
      logical ok
      integer nchord
      character*80, allocatable :: setname(:)
      integer, allocatable :: cnum(:), npt(:), npoints(:)
      real (kind=R8), allocatable :: cval(:), rstart(:), zstart(:),
     & phistart(:), rend(:), zend(:), phiend(:)
      integer status, ncid, nset, nc, is, ic, m, np, dimid, varid
      integer, allocatable :: sgrp(:), cgrp(:)
      character*(NF_MAX_NAME) gname
      real (kind=R8) :: rtmp(1)
      ok = .false.
      nchord = 0
      allocate(setname(0), cnum(0), npt(0), npoints(0), cval(0))
      allocate(rstart(0), zstart(0), phistart(0))
      allocate(rend(0), zend(0), phiend(0))
      status = nf_open(filename, NF_NOWRITE, ncid)
      if (status.ne.NF_NOERR) return
      allocate(sgrp(maxcgrp), cgrp(maxcgrp))
      nset = 0
      status = nf_inq_grps(ncid, nset, sgrp)
      is = 0
      do while (status.eq.NF_NOERR .and. is.lt.nset)
        is = is + 1
        status = nf_inq_grps(sgrp(is), nc, cgrp)
        if (status.eq.NF_NOERR) nchord = nchord + nc
      enddo
      if (status.eq.NF_NOERR .and. nchord.ge.1) then
        deallocate(setname, cnum, npt, npoints, cval)
        deallocate(rstart, zstart, phistart, rend, zend, phiend)
        allocate(setname(nchord), cnum(nchord), npt(nchord))
        allocate(npoints(nchord), cval(nchord))
        allocate(rstart(nchord), zstart(nchord), phistart(nchord))
        allocate(rend(nchord), zend(nchord), phiend(nchord))
        m = 0
        is = 0
        do while (status.eq.NF_NOERR .and. is.lt.nset)
          is = is + 1
          status = nf_inq_grpname(sgrp(is), gname)
          if (status.eq.NF_NOERR)
     &      status = nf_get_att_int(sgrp(is), NF_GLOBAL, 'npoints', np)
          if (status.eq.NF_NOERR)
     &      status = nf_inq_grps(sgrp(is), nc, cgrp)
          ic = 0
          do while (status.eq.NF_NOERR .and. ic.lt.nc)
            ic = ic + 1
            m = m + 1
            setname(m) = gname(1:80)
            npoints(m) = np
            status = nf_inq_dimid(cgrp(ic), 'npt', dimid)
            if (status.eq.NF_NOERR)
     &        status = nf_inq_dimlen(cgrp(ic), dimid, npt(m))
            if (status.eq.NF_NOERR)
     &        status = nf_inq_varid(cgrp(ic), 'chord_id', varid)
            if (status.eq.NF_NOERR)
     &        status = nf_get_var_int(cgrp(ic), varid, cnum(m))
            if (status.eq.NF_NOERR)
     &        status = nf_inq_varid(cgrp(ic), 'R', varid)
            call get_dbl('chord_value', cval(m))
            call get_dbl('Rstart', rstart(m))
            call get_dbl('zstart', zstart(m))
            call get_dbl('phistart', phistart(m))
            call get_dbl('Rend', rend(m))
            call get_dbl('zend', zend(m))
            call get_dbl('phiend', phiend(m))
          enddo
        enddo
        ok = status.eq.NF_NOERR
      endif
      deallocate(sgrp, cgrp)
      status = nf_close(ncid)
      return

      contains

      subroutine get_dbl(vname, val)
c     Read the scalar vname of chord group cgrp(ic), unless an earlier
c     read failed.
      implicit none
      character*(*) vname
      real (kind=R8) :: val
      if (status.ne.NF_NOERR) return
      status = nf_inq_varid(cgrp(ic), vname, varid)
      if (status.eq.NF_NOERR)
     &  status = nf_get_var_double(cgrp(ic), varid, rtmp)
      if (status.eq.NF_NOERR) val = rtmp(1)
      return
      end subroutine get_dbl

      end subroutine get_cdf_chord_meta
#endif

      end module b2mod_chords
