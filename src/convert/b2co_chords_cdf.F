      subroutine b2co_chords_cdf
#ifndef NO_CDF
     . (ns, tim, switch)
#endif
      use b2mod_types
      use b2mod_switches
      use b2mod_subsys
      implicit none
#ifndef NO_CDF
      integer ns
      real (kind=R8) :: tim
      type (switches), intent(in) :: switch
#     include <netcdf.inc>
#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
      integer, parameter :: maxv = 4096, maxg = 1024
      character*4 :: extension
      character*256 :: plasmastate
#endif
      external xerrab, xertst

      call subini ('b2co_chords_cdf')
      write(*,*) 'b2co_chords_cdf ...'

#ifndef NO_CDF
      call xertst (0.lt.ns, 'faulty argument ns')

c   ..plasma-state extension for the truncated output file names
      extension='xxxx'
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
      CALL PXFGETENV ('PLASMASTATE', 0, plasmastate, lenval, ierror)
      if (ierror.eq.0) extension=plasmastate(lenval-3:lenval)
#else
      call get_environment_variable('PLASMASTATE',
     . status=ierror, length=lenval)
      if (ierror.eq.0) then
        call get_environment_variable('PLASMASTATE', value=plasmastate)
        extension=plasmastate(lenval-3:lenval)
      endif
#endif
#endif

      call truncate_chord_file('b2chords.nc')
#ifdef B25_EIRENE
      if (switch%write_chord_triangles.ge.1)
     &  call truncate_chord_file('eirenechords.nc')
#endif

#else
      write(*,*)
     1 'b2co_chords_cdf failed because compiled with NO_CDF option'
#endif
      call subend ()
      return

#ifndef NO_CDF
      contains

c-----------------------------------------------------------------------
      subroutine truncate_chord_file(fname)
c     Copy fname to fname.<extension>, truncating the time dimension to
c     the frames up to the restart time tim.
      implicit none
      character*(*) fname
      integer ncid, ncid2, ng, i, itim, status
      integer ib, it, ih, tlen, tvid, did
      integer gsrc(maxg), gdst(maxg)
      character*(NF_MAX_NAME) gname
      character*256 dstname
      real (kind=R8), allocatable :: tvals(:)
      logical exists

      call find_file(fname, exists)
      status = nf_open(fname, NF_NOWRITE, ncid)
      if (status.ne.NF_NOERR) then
        write(*,*) 'Skipped '//trim(fname)//' truncation: not found'
        return
      endif

c   ..truncation index itim from the root time vector
      itim = 0
      status = nf_inq_dimid(ncid, 'time', did)
      if (status.eq.NF_NOERR) then
        status = nf_inq_dimlen(ncid, did, tlen)
        if (status.eq.NF_NOERR .and. tlen.gt.0) then
          status = nf_inq_varid(ncid, 'times', tvid)
          call check_cdf_status(status)
          allocate(tvals(tlen))
          status = nf_get_var_double(ncid, tvid, tvals)
          call check_cdf_status(status)
          ib = 1
          it = tlen
          ih = (ib+it)/2
          do while (it-ib.gt.1)
            ih = (ib+it)/2
            if (tim.lt.tvals(ih)) then
              it = ih
            else if (tim.gt.tvals(ih)) then
              ib = ih
            else
              exit
            endif
          end do
          itim = ih
          deallocate(tvals)
        endif
      endif

c   ..create the truncated copy (netCDF-4)
      dstname = trim(fname)//'.'//extension
      status = nf_create(dstname, or(NF_NETCDF4,NF_CLOBBER), ncid2)
      call check_cdf_status(status)
c   ..definitions: root group, then one group per chord
      call copy_defs(ncid, ncid2)
      status = nf_inq_grps(ncid, ng, gsrc)
      call check_cdf_status(status)
      do i = 1, ng
        status = nf_inq_grpname(gsrc(i), gname)
        call check_cdf_status(status)
        status = nf_def_grp(ncid2, trim(gname), gdst(i))
        call check_cdf_status(status)
        call copy_defs(gsrc(i), gdst(i))
      enddo
      status = nf_enddef(ncid2)
      call check_cdf_status(status)
c   ..data: root group, then each chord group (truncating time to itim)
      call copy_data(ncid, ncid2, itim)
      do i = 1, ng
        call copy_data(gsrc(i), gdst(i), itim)
      enddo

      status = nf_close(ncid)
      call check_cdf_status(status)
      status = nf_close(ncid2)
      call check_cdf_status(status)
      write(*,*) trim(dstname)//' written'
      return
      end subroutine truncate_chord_file

c-----------------------------------------------------------------------
      subroutine copy_defs(src, dst)
      implicit none
      integer src, dst
      integer ndim, dimids(NF_MAX_VAR_DIMS), i, j, nvar
      integer varids(maxv), xtype, nvd, vdims(NF_MAX_VAR_DIMS), nva, k
      integer ddims(NF_MAX_VAR_DIMS), dvid, ddum, dlen, status
      character*(NF_MAX_NAME) dname, vname, aname
      status = nf_inq_dimids(src, ndim, dimids, 0)
      call check_cdf_status(status)
      do i = 1, ndim
        status = nf_inq_dim(src, dimids(i), dname, dlen)
        call check_cdf_status(status)
        if (dname.eq.'time') then
          status = nf_def_dim(dst, dname, NF_UNLIMITED, ddum)
        else
          status = nf_def_dim(dst, dname, dlen, ddum)
        endif
        call check_cdf_status(status)
      enddo
      status = nf_inq_varids(src, nvar, varids)
      call check_cdf_status(status)
      do i = 1, nvar
        status = nf_inq_var(src, varids(i), vname, xtype, nvd, vdims,
     &    nva)
        call check_cdf_status(status)
        do j = 1, nvd
          status = nf_inq_dim(src, vdims(j), dname, dlen)
          call check_cdf_status(status)
          status = nf_inq_dimid(dst, dname, ddims(j))
          call check_cdf_status(status)
        enddo
        status = nf_def_var(dst, vname, xtype, nvd, ddims, dvid)
        call check_cdf_status(status)
        do k = 1, nva
          status = nf_inq_attname(src, varids(i), k, aname)
          call check_cdf_status(status)
          status = nf_copy_att(src, varids(i), aname, dst, dvid)
          call check_cdf_status(status)
        enddo
      enddo
      return
      end subroutine copy_defs

c-----------------------------------------------------------------------
      subroutine copy_data(src, dst, itim)
      implicit none
      integer src, dst, itim
      integer nvar, varids(maxv), i, j, svid, dvid, xtype, nvd, nva
      integer vdims(NF_MAX_VAR_DIMS), dimsize(NF_MAX_VAR_DIMS)
      integer start(NF_MAX_VAR_DIMS), count(NF_MAX_VAR_DIMS), index(1)
      integer status, ival
      logical is_time
      character*(NF_MAX_NAME) vname, dname
      real (kind=R8) :: rtmp(1)
      real (kind=R8), allocatable :: data1(:), data2(:,:), data3(:,:,:)
      status = nf_inq_varids(src, nvar, varids)
      call check_cdf_status(status)
      do i = 1, nvar
        svid = varids(i)
        status = nf_inq_var(src, svid, vname, xtype, nvd, vdims, nva)
        call check_cdf_status(status)
        status = nf_inq_varid(dst, vname, dvid)
        call check_cdf_status(status)
        is_time = .false.
        do j = 1, nvd
          status = nf_inq_dim(src, vdims(j), dname, dimsize(j))
          call check_cdf_status(status)
          start(j) = 1
          count(j) = dimsize(j)
          if (dname.eq.'time') is_time = .true.
        enddo
        if (xtype.eq.NF_INT) then
c   ..scalar integer counters (ntstep, nchord, chord_id)
          index(1) = 1
          if (vname.eq.'ntstep') then
            status = nf_put_var1_int(dst, dvid, index, itim)
            call check_cdf_status(status)
          else
            status = nf_get_var1_int(src, svid, index, ival)
            call check_cdf_status(status)
            status = nf_put_var1_int(dst, dvid, index, ival)
            call check_cdf_status(status)
          endif
        else if (nvd.eq.0) then
          status = nf_get_var_double(src, svid, rtmp)
          call check_cdf_status(status)
          status = nf_put_var_double(dst, dvid, rtmp)
          call check_cdf_status(status)
        else if (nvd.eq.1) then
          allocate(data1(dimsize(1)))
          status = nf_get_vara_double(src, svid, start, count, data1)
          call check_cdf_status(status)
          if (is_time) count(1) = itim
          if (.not.is_time .or. itim.gt.0)
     &      status = nf_put_vara_double(dst, dvid, start, count, data1)
          call check_cdf_status(status)
          deallocate(data1)
        else if (nvd.eq.2) then
          allocate(data2(dimsize(1),dimsize(2)))
          status = nf_get_vara_double(src, svid, start, count, data2)
          call check_cdf_status(status)
          if (is_time) count(2) = itim
          if (.not.is_time .or. itim.gt.0)
     &      status = nf_put_vara_double(dst, dvid, start, count, data2)
          call check_cdf_status(status)
          deallocate(data2)
        else if (nvd.eq.3) then
          allocate(data3(dimsize(1),dimsize(2),dimsize(3)))
          status = nf_get_vara_double(src, svid, start, count, data3)
          call check_cdf_status(status)
          if (is_time) count(3) = itim
          if (.not.is_time .or. itim.gt.0)
     &      status = nf_put_vara_double(dst, dvid, start, count, data3)
          call check_cdf_status(status)
          deallocate(data3)
        else
          write(*,'(a,i2,a)')
     *     'Case for ',nvd,' dimensions not yet coded'
        endif
      enddo
      return
      end subroutine copy_data
#endif

      end subroutine b2co_chords_cdf
