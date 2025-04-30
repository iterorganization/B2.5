      subroutine b2co_movies_cdf
#ifndef NO_CDF
     . (nCv, ns, natm, nmol, tim)
#endif
      use b2mod_types
      use b2mod_movies
      use b2mod_subsys
      implicit none
#ifndef NO_CDF
      integer nCv, ns, natm, nmol
      real (kind=R8) :: tim
#     include <netcdf.inc>
#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
      integer :: status
      integer :: idim, ivar, iatt, ib, it, ih, idum,
     1 itim, ncid, ncid2, ndims, nvars, natts, unlimid, xtype,
     2 dimids(NF_MAX_VAR_DIMS), len, dimsize(NF_MAX_VAR_DIMS),
     3 start(NF_MAX_VAR_DIMS), count(NF_MAX_VAR_DIMS),
     4 index(NF_MAX_VAR_DIMS)
      real (kind=R8) :: scale
      real (kind=R8), allocatable :: data1(:), data2(:,:), data3(:,:,:),
     &                               data4(:,:,:,:)
      character*(NF_MAX_NAME) :: varname, dimname, attname
      character :: units*(24)
      character*4 :: extension
      character*256 filename, plasmastate
      logical :: exists
#endif
      integer, save :: fields_selection = 1
      external ipgeti, ipgetr, xerrab, xertst

      call ipgeti ('cdfmovie_fields',fields_selection)
      call subini ('b2co_movies_cdf')
      write(*,*) 'b2co_movies_cdf ...'

#ifndef NO_CDF
*   ..test input
      call xertst (1.le.nCv, 'faulty argument nCv')
      call xertst (0.lt.ns, 'faulty argument ns')
      call xertst (0.le.natm, 'faulty argument natm')
      call xertst (0.le.nmol, 'faulty argument nmol')
c access the b2movies.nc file
      filename='b2movies.nc'
      call find_file(filename,exists)
      status=nf_open(filename,NF_NOWRITE,ncid)
      call check_cdf_status(status)

      status=nf_inq(ncid,ndims,nvars,natts,unlimid)
      call check_cdf_status(status)

      itim=0
      timedim=-1

c loop over all variables in the file
      do ivar=1,nvars
        status=nf_inq_var(ncid,ivar,varname,xtype,
     1   ndims,dimids,natts)
        call check_cdf_status(status)
        do idim=1,ndims
          status=nf_inq_dim(ncid,dimids(idim),dimname,len)
          call check_cdf_status(status)
          dimsize(idim)=len
        enddo
        units=''
        scale=1
        do iatt=1,natts
          status=nf_inq_attname(ncid,ivar,iatt,attname)
          call check_cdf_status(status)
          status=nf_inq_att(ncid,ivar,attname,xtype,len)
          call check_cdf_status(status)
          if(attname.eq.'units') then
            if(len.le.24) then
              status=nf_get_att_text(ncid,ivar,attname,units)
              call check_cdf_status(status)
            else
              call xerrab ('increase size of units')
            endif
          endif
          if(attname.eq.'scale') then
            if(len.le.1) then
              status=nf_get_att_double(ncid,ivar,attname,scale)
              call check_cdf_status(status)
            else
              call xerrab ('increase size of scale')
            endif
          endif
        enddo
      enddo

c find the time vector size
      do ivar=1,nvars
        status=nf_inq_var(ncid,ivar,varname,xtype,
     1   ndims,dimids,natts)
        call check_cdf_status(status)
        do idim=1,ndims
          status=nf_inq_dim(ncid,dimids(idim),dimname,len)
          call check_cdf_status(status)
          dimsize(idim)=len
        enddo
        if(varname.eq.'times') then
          allocate(data1(dimsize(1)))
          status=nf_get_var_double(ncid,ivar,data1)
          timedim=dimsize(1)
          call check_cdf_status(status)
          ib=1
          it=timedim
          ih=(ib+it)/2
          do while (it-ib.gt.1)
            ih=(ib+it)/2
            if(tim.lt.data1(ih)) then
              it=ih
            else if(tim.gt.data1(ih)) then
              ib=ih
            else
              exit
            endif
          end do
          itim=ih
          deallocate(data1)
        end if
      end do

c prepare the b2movies.nc.xxxx file
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
      filename='b2movies.nc.'//extension
      idum=timedim
      timedim=itim
      call create_cdf_start(ncid2, filename, nCv, ns, natm, nmol)
      if (fields_selection.ge.1) then
        call create_cdf_field3(ncid2, 'na', 1,
     1   'Ion density', 'm^-3')
        call create_cdf_field2(ncid2, 'ne',
     1   'Electron density', 'm^-3')
        call create_cdf_field2(ncid2, 'te',
     1   'Electron temperature', 'eV')
        call create_cdf_field2(ncid2, 'ti',
     1   'Ion temperature', 'eV')
        call create_cdf_field2(ncid2, 'po',
     1   'Electric potential', 'V')
        if(natm.gt.0) then
          call create_cdf_field3(ncid2, 'dab2', 2,
     1     'Atom density', 'm^-3')
          call create_cdf_field3(ncid2, 'tab2', 2,
     1     'Atom temperature', 'eV')
        endif
        if(nmol.gt.0) then
          call create_cdf_field3(ncid2, 'dmb2', 3,
     1     'Molecular density', 'm^-3')
          call create_cdf_field3(ncid2, 'tmb2', 3,
     1     'Molecule temperature', 'eV')
        endif
        if (fields_selection.ge.2) then
          call create_cdf_field3(ncid2, 'rsana', 1,
     1      'Volumetric ionisation rate', 'm^-3 s^-1')
          call create_cdf_field3(ncid2, 'rsahi', 1,
     1      'Volumetric ionisation energy source/sink', 'W m^-3')
          call create_cdf_field2(ncid2, 'rsahisum',
     1      'Total volumetric ionisation energy source/sink', 'W m^-3')
          call create_cdf_field3(ncid2, 'rrana', 1,
     1      'Volumetric recombination rate', 'm^-3 s^-1')
          call create_cdf_field3(ncid2, 'rrahi', 1,
     1      'Volumetric recombination energy source/sink', 'W m^-3')
          call create_cdf_field2(ncid2, 'rrahisum',
     1      'Total volumetric recombination energy source/sink',
     2      'W m^-3')
          call create_cdf_field3(ncid2, 'rcxna', 1,
     1      'Volumetric charge exchange rate', 'm^-3 s^-1')
          call create_cdf_field3(ncid2, 'rcxhi', 1,
     1      'Volumetric charge exchange energy source/sink', 'W m^-3')
          call create_cdf_field2(ncid2, 'rcxhisum',
     1      'Total volumetric charge exchange energy source/sink',
     2      'W m^-3')
          call create_cdf_field3(ncid2, 'rqahe', 1,
     1      'Volumetric electron cooling rate', 'W m^-3')
          call create_cdf_field2(ncid2, 'rqahesum',
     1      'Total volumetric electron cooling rate', 'W m^-3')
          call create_cdf_field3(ncid2, 'rqrad', 1,
     1      'Volumetric line radiation rate', 'W m^-3')
          call create_cdf_field2(ncid2, 'rqradsum',
     1      'Total volumetric line radiation rate', 'W m^-3')
          if(natm.gt.0) then
            call create_cdf_field3(ncid2, 'eneutrad', 2,
     1       'Volumetric radiation rate due to atoms', 'W m^-3')
            call create_cdf_field2(ncid2, 'eneutradsum',
     1       'Total volumetric radiation rate due to atoms', 'W m^-3')
          endif
          if(nmol.gt.0) then
            call create_cdf_field3(ncid2, 'emolrad', 3,
     1       'Volumetric radiation rate due to molecules', 'W m^-3')
            call create_cdf_field2(ncid2, 'emolradsum',
     1       'Total volumetric radiation rate due to molecules',
     2       'W m^-3')
          endif
        endif
      endif
      call create_cdf_end(ncid2)
      timedim=idum

c read and rewrite data
      status=nf_open(filename,NF_WRITE,ncid2)
      call check_cdf_status(status)
      do ivar=1,nvars
        status=nf_inq_var(ncid,ivar,varname,xtype,
     1   ndims,dimids,natts)
        call check_cdf_status(status)
        count=0
        do idim=1,ndims
          status=nf_inq_dim(ncid,dimids(idim),dimname,len)
          call check_cdf_status(status)
          dimsize(idim)=len
          count(idim)=len
          start(idim)=1
        enddo
        if(varname.eq.'times') then
          allocate(data1(dimsize(1)))
          status=nf_get_vara_double(ncid,ivar,start,count,data1)
          call check_cdf_status(status)
          count(1)=itim
          status=nf_put_vara_double(ncid2,ivar,start,count,data1)
          deallocate(data1)
        elseif(varname.eq.'ntstep') then
          index(1)=1
          status=nf_get_var1_int(ncid,ivar,index,idum)
          call check_cdf_status(status)
          status=nf_put_var1_int(ncid2,ivar,index,itim)
          call check_cdf_status(status)
        else
          if(ndims.eq.4) then
            allocate(data4(dimsize(1),dimsize(2),dimsize(3),dimsize(4)))
            status=nf_get_vara_double(ncid,ivar,start,count,data4)
            call check_cdf_status(status)
            if (dimsize(4).eq.timedim) then
              count(4)=itim
              if (itim.gt.0) status=nf_put_vara_double(ncid2,ivar,
     .         start,count,data4)
            else
              status=nf_put_vara_double(ncid2,ivar,start,count,data4)
            end if
            call check_cdf_status(status)
            deallocate(data4)
          elseif(ndims.eq.3) then
            allocate(data3(dimsize(1),dimsize(2),dimsize(3)))
            status=nf_get_vara_double(ncid,ivar,start,count,data3)
            call check_cdf_status(status)
            if (dimsize(3).eq.timedim) then
              count(3)=itim
              if (itim.gt.0) status=nf_put_vara_double(ncid2,ivar,
     .                                            start,count,data3)
            else
              status=nf_put_vara_double(ncid2,ivar,start,count,data3)
            end if
            call check_cdf_status(status)
            deallocate(data3)
          elseif(ndims.eq.2) then
            allocate(data2(dimsize(1),dimsize(2)))
            status=nf_get_vara_double(ncid,ivar,start,count,data2)
            call check_cdf_status(status)
            if (dimsize(2).eq.timedim) then
              count(2)=itim
              if (itim.gt.0) status=nf_put_var_double(ncid2,ivar,
     .                                           start,count,data2)
            else
              status=nf_put_vara_double(ncid2,ivar,start,count,data2)
            end if
            call check_cdf_status(status)
            deallocate(data2)
          elseif(ndims.eq.1) then
            allocate(data1(dimsize(1)))
            status=nf_get_vara_double(ncid,ivar,start,count,data1)
            call check_cdf_status(status)
            if (dimsize(1).eq.timedim) then
              count(1)=itim
              status=nf_put_vara_double(ncid2,ivar,start,count,data1)
            else
              status=nf_put_vara_double(ncid2,ivar,start,count,data1)
            end if
            call check_cdf_status(status)
            deallocate(data1)
          else
            write(*,'(a,i2,a)')
     *       'Case for ',ndims,' dimensions not yet coded'
          endif
        endif
      enddo

      status=nf_close(ncid)
      call check_cdf_status(status)
      status=nf_close(ncid2)
      call check_cdf_status(status)

#else
      write(*,*)
     1 'b2co_movies_cdf failed because compiled with NO_CDF option'
#endif
      call subend ()
      return
      end subroutine b2co_movies_cdf
