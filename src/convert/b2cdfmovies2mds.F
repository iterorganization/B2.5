      subroutine b2cdfmovies2mds
      use b2mod_subsys
#ifdef MDSPLUS
      use b2mod_types
      implicit none
#     include "mdslib.inc"
#ifndef NO_CDF
#     include <netcdf.inc>
#endif
      logical :: want_transpose
#ifndef NO_CDF
      real (kind=R8), allocatable  :: data1(:),
     1 data2(:,:), data2T(:,:), data3(:,:,:), data3T(:,:,:),
     2 data4(:,:,:,:), data4T(:,:,:,:)
      integer :: i,ivar,idim,iatt,i1,i2,i3,i4
      integer :: status,
     1 ncid, ndims, nvars, natts, unlimid, xtype,
     2 dimids(NF_MAX_VAR_DIMS), length,
     3 dimsize(NF_MAX_VAR_DIMS), underscore_pos
      character*(NF_MAX_NAME) ::
     1 varname, dimname(NF_MAX_VAR_DIMS), attname
      character :: mds_signal_name*(64), units*(24),
     1 tdi*(172)
      character*256 filename
      logical exist
#endif
      external xerrab

      call subini ('b2cdfmovies2mds')
      write(*,*) 'b2cdfmovies2mds ...'

      want_transpose=.false.
      want_transpose=.true.

#ifndef NO_CDF
c access the b2movies.nc file
      filename='b2movies.nc'
      call find_file(filename,exist)
      status=nf_open(filename,NF_NOWRITE,ncid)
      call check_cdf_status(status)

      status=nf_inq(ncid,ndims,nvars,natts,unlimid)
      call check_cdf_status(status)

c loop over all variables in the file
      do ivar=1,nvars
        status=nf_inq_var(ncid,ivar,varname,xtype,
     1   ndims,dimids,natts)
        call check_cdf_status(status)
        write(*,'(a,a,a,i2,a,i2,a)')
     *   'variable: ',trim(varname),' has type ',xtype,
     1   ' and ',ndims,' dimensions'
        do idim=1,ndims
          status=nf_inq_dim(ncid,dimids(idim),dimname(idim),
     1     dimsize(idim))
          call check_cdf_status(status)
          write(*,'(a,a,a,i6)')
     *     '    dimension: ',trim(dimname(idim)),
     1     ' of length ',dimsize(idim)
        enddo
        units=''
        do iatt=1,natts
          status=nf_inq_attname(ncid,ivar,iatt,attname)
          call check_cdf_status(status)
          status=nf_inq_att(ncid,ivar,attname,xtype,length)
          call check_cdf_status(status)
          write(*,'(a,a,a,i2,a,i6)')
     *     '    attribute: ',trim(attname),
     1     ' of type ',xtype,' and length ',length
          if(attname.eq.'units') then
            if(length.le.24) then
              status=nf_get_att_text(ncid,ivar,attname,units)
              call check_cdf_status(status)
              if(length.gt.0) then
                if(units(length:length).eq.char(0))
     1           units(length:length)=' '
              endif
            else
              call xerrab ('increase size of units')
            endif
          endif
        enddo
        write(*,*) '  units = ',trim(units)

c first deal with the time vector
        if(varname.eq.'times') then
          allocate(data1(dimsize(1)))
          status=nf_get_var_double(ncid,ivar,data1)
          call check_cdf_status(status)
          do i=1,dimsize(1)
            if(data1(i).eq.fildoub) then
              if(i.gt.1) then
                if(i.lt.dimsize(1)) then
                  if(data1(i+1).ne.fildoub) then
                    data1(i)=(data1(i-1)+data1(i+1))/2
                    write(*,*) '[1] Fixed index ',i,' to ',data1(i)
                  else if(i.gt.2) then
                    data1(i)=data1(i-1)+(data1(i-1)-data1(i-2))
                    write(*,*) '[2] Fixed index ',i,' to ',data1(i)
                  else
                    write(*,*) '[3] Cannot fix index ',i
                  endif
                else if(i.gt.2) then
                  data1(i)=data1(i-1)+(data1(i-1)-data1(i-2))
                  write(*,*) '[4] Fixed index ',i,' to ',data1(i)
                else
                  write(*,*) '[5] Cannot fix index ',i
                endif
              else if(i+2.le.dimsize(1)) then
                if(data1(i+1).ne.fildoub.and.
     1           data1(i+2).ne.fildoub) then
                  data1(i)=data1(i+1)-(data1(i+2)-data1(i+1))
                  write(*,*) '[6] Fixed index ',i,' to ',data1(i)
                else
                  write(*,*) '[7] Cannot fix index ',i
                endif
              else
                write(*,*) '[8] Cannot fix index ',i
              endif
            endif
          enddo
          mds_signal_name='\TOP.MOVIES:TIME'
          write(*,'(a,a,a,a)')
     *     'We need to save ',trim(varname),
     1     ' to ',trim(mds_signal_name)
          status=mdsput2(
     *     trim(mds_signal_name)//char(0),
     *     'build_with_units($,"s")'//char(0),
     *     descr2(idtype_double,dimsize(1),0),data1,
     *     0)
          write(*,*) 'mdsput2'
          call mds_err_msg(status)
          deallocate(data1)
        elseif(ndims.eq.2) then     ! NOT USED AT THE MOMENT
          mds_signal_name='\TOP.MOVIES:' // trim(varname)
          underscore_pos=index(mds_signal_name,'_')
          do while(underscore_pos.gt.0)
            mds_signal_name(underscore_pos:len(mds_signal_name))=
     1       mds_signal_name(underscore_pos+1:len(mds_signal_name))
            underscore_pos=index(mds_signal_name,'_')
          enddo
          write(*,'(a,a,a,a)')
     *     'We need to save ',trim(varname),
     1     ' to ',trim(mds_signal_name)
          allocate(data2(dimsize(1),dimsize(2)))
          status=nf_get_var_double(ncid,ivar,data2)
          call check_cdf_status(status)
          if(.not.want_transpose) then
            tdi=
     *       'build_signal(build_with_units($,"' //
     *       trim(units) //
     *       '"),,' //
     *       '\TOP.MOVIES:' // trim(dimname(1)) //',' //
     *       '\TOP.MOVIES:TIME)'
            write(*,'(a,a)') 'Using TDI ',trim(tdi)
            status=mdsput2(
     *       trim(mds_signal_name)//char(0),
     *       tdi//char(0),
     *       descr2(idtype_double,dimsize(1),dimsize(2),0),data2,0)
            write(*,*) 'mdsput2'
            call mds_err_msg(status)
          else
            allocate(data2T(dimsize(2),dimsize(1)))
            data2T=transpose(data2)
            tdi=
     *       'build_signal(build_with_units($,"' //
     *       trim(units) //
     *       '"),,\TOP.MOVIES:TIME,' //
     *       '\TOP.MOVIES:' // trim(dimname(1)) // ')'
            write(*,'(a,a)') 'Using TDI ',trim(tdi)
            status=mdsput2(
     *       trim(mds_signal_name)//char(0),
     *       tdi//char(0),
     *       descr2(idtype_double,dimsize(2),dimsize(1),0),data2T,0)
            write(*,*) 'mdsput2'
            call mds_err_msg(status)
            deallocate(data2T)
          endif
          deallocate(data2)
        elseif(ndims.eq.3) then
          mds_signal_name='\TOP.MOVIES:' // trim(varname)
          underscore_pos=index(mds_signal_name,'_')
          do while(underscore_pos.gt.0)
            mds_signal_name(underscore_pos:len(mds_signal_name))=
     1       mds_signal_name(underscore_pos+1:len(mds_signal_name))
            underscore_pos=index(mds_signal_name,'_')
          enddo
          write(*,'(a,a,a,a)') 'We need to save ',trim(varname),
     1     ' to ',trim(mds_signal_name)
          allocate(data3(dimsize(1),dimsize(2),dimsize(3)))
          status=nf_get_var_double(ncid,ivar,data3)
          call check_cdf_status(status)
          if(.not.want_transpose) then
            tdi=
     *       'build_signal(build_with_units($,"' //
     *       trim(units) //
     *       '"),,' //
     *       '\TOP.MOVIES:' // trim(dimname(1)) // ',' //
     *       '\TOP.MOVIES:' // trim(dimname(2)) // ',' //
     *       '\TOP.MOVIES:TIME)'
            write(*,'(a,a)') 'Using TDI ',trim(tdi)
            status=mdsput2(
     *       trim(mds_signal_name)//char(0),
     *       tdi//char(0),
     *       descr2(idtype_double,dimsize(1),dimsize(2),
     1       dimsize(3),0),data3,0)
            write(*,*) 'mdsput2'
            call mds_err_msg(status)
          else
            allocate(data3T(dimsize(3),dimsize(1),dimsize(2)))
            forall(i1=1:dimsize(1),i2=1:dimsize(2),i3=1:dimsize(3))
     *       data3T(i3,i1,i2)=data3(i1,i2,i3)
            tdi=
     *       'build_signal(build_with_units($,"' //
     *       trim(units) //
     *       '"),,\TOP.MOVIES:TIME,' //
     *       '\TOP.MOVIES:' // trim(dimname(1)) // ',' //
     *       '\TOP.MOVIES:' // trim(dimname(2)) // ')'
            write(*,'(a,a)') 'Using TDI ',trim(tdi)
            status=mdsput2(
     *       trim(mds_signal_name)//char(0),
     *       tdi//char(0),
     *       descr2(idtype_double,dimsize(3),dimsize(1),
     1       dimsize(2),0),data3T,0)
            write(*,*) 'mdsput2'
            call mds_err_msg(status)
            deallocate(data3T)
          endif
          deallocate(data3)
        elseif(ndims.eq.4) then
          mds_signal_name='\TOP.MOVIES:' // trim(varname)
          underscore_pos=index(mds_signal_name,'_')
          do while(underscore_pos.gt.0)
            mds_signal_name(underscore_pos:len(mds_signal_name))=
     1       mds_signal_name(underscore_pos+1:len(mds_signal_name))
            underscore_pos=index(mds_signal_name,'_')
          enddo
          write(*,'(a,a,a,a)') 'We need to save ',trim(varname),
     1     ' to ',trim(mds_signal_name)
          allocate(data4(dimsize(1),dimsize(2),dimsize(3),dimsize(4)))
          status=nf_get_var_double(ncid,ivar,data4)
          call check_cdf_status(status)
          if(.not.want_transpose) then
            tdi=
     *       'build_signal(build_with_units($,"' //
     *       trim(units) //
     *       '"),,' //
     *       '\TOP.MOVIES:' // trim(dimname(1)) // ',' //
     *       '\TOP.MOVIES:' // trim(dimname(2)) // ',' //
     *       '\TOP.MOVIES:' // trim(dimname(3)) // ',' //
     *       '\TOP.MOVIES:TIME)'
            write(*,'(a,a)') 'Using TDI ',trim(tdi)
            status=mdsput2(
     *       trim(mds_signal_name)//char(0),
     *       tdi//char(0),
     *       descr2(idtype_double,dimsize(1),dimsize(2),
     1       dimsize(3),dimsize(4),0),data4,0)
            write(*,*) 'mdsput2'
            call mds_err_msg(status)
          else
            allocate(data4T(dimsize(4),dimsize(1),dimsize(2),
     1       dimsize(3)))
            forall(i1=1:dimsize(1),i2=1:dimsize(2),i3=1:dimsize(3),
     1       i4=1:dimsize(4))
     *       data4T(i4,i1,i2,i3)=data4(i1,i2,i3,i4)
            tdi=
     *       'build_signal(build_with_units($,"' //
     *       trim(units) //
     *       '"),,\TOP.MOVIES:TIME,' //
     *       '\TOP.MOVIES:' // trim(dimname(1)) // ',' //
     *       '\TOP.MOVIES:' // trim(dimname(2)) // ',' //
     *       '\TOP.MOVIES:' // trim(dimname(3)) // ')'
            write(*,'(a,a)') 'Using TDI ',trim(tdi)
            status=mdsput2(
     *       trim(mds_signal_name)//char(0),
     *       tdi//char(0),
     *       descr2(idtype_double,dimsize(4),dimsize(1),
     1       dimsize(2),dimsize(3),0),data4T,0)
            write(*,*) 'mdsput2'
            call mds_err_msg(status)
            deallocate(data4T)
          endif
          deallocate(data4)
        else
          write(*,'(a,i2,a)')
     *     'Case for ',ndims,' dimensions not yet coded'
        endif
      enddo

      status=nf_close(ncid)
      call check_cdf_status(status)
#else
      write(*,*)
     . 'b2cdfmovies2mds failed because compiled with NO_CDF option'
#endif
#else
      call xerrab('b2cdfmovies2mds: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      call subend ()
      return
      end subroutine b2cdfmovies2mds
