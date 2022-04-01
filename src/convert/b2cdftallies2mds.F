      subroutine b2cdftallies2mds
      use b2mod_subsys
#ifdef MDSPLUS
      use b2mod_types
      use b2mod_mdsplus
      implicit none
#     include "mdslib.inc"
#ifndef NO_CDF
#     include <netcdf.inc>
#endif
      logical :: want_transpose
#ifndef NO_CDF
      real (kind=R8), allocatable  :: data1(:),
     1 data2(:,:), data2T(:,:), data3(:,:,:), data3T(:,:,:)
      character, allocatable :: data2c(:,:)
#ifdef F2003
      character(len=:), allocatable :: str1(:)
#else
      integer, parameter :: strmax = 1024
      character(len=strmax+1), allocatable :: str1(:)
#endif
      integer :: i,ivar,idim,iatt,i1,i2,i3
      integer :: status,
     1 ncid, ndims, nvars, natts, unlimid, xtypev, xtypeu,
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

      call subini ('b2cdftallies2mds')
      write(*,*) 'b2cdftallies2mds ...'

      want_transpose=.false.
      want_transpose=.true.

#ifndef NO_CDF
c access the b2tallies.nc file
      filename='b2tallies.nc'
      call find_file(filename,exist)
      status=nf_open(filename,NF_NOWRITE,ncid)
      call check_cdf_status(status)

      status=nf_inq(ncid,ndims,nvars,natts,unlimid)
      call check_cdf_status(status)

c loop over all variables in the file
      do ivar=1,nvars
        status=nf_inq_var(ncid,ivar,varname,xtypev,
     1   ndims,dimids,natts)
        call check_cdf_status(status)
        if(ndims.ne.1) then
          write(*,'(a,a,a,i2,a,i6,a)')
     *     'variable: ',trim(varname),' has type ',xtypev,
     1     ' and ',ndims,' dimensions'
        else
          write(*,'(a,a,a,i2,a,i6,a)')
     *     'variable: ',trim(varname),' has type ',xtypev,
     1     ' and ',ndims,' dimension'
        endif
        do idim=1,ndims
          status=nf_inq_dim(ncid,dimids(idim),dimname(idim),
     1     dimsize(idim))
          call check_cdf_status(status)
          write(*,'(a,a,a,i6)') '    dimension: ',trim(dimname(idim)),
     1     ' of length ',dimsize(idim)
        enddo
        units=''
        do iatt=1,natts
          status=nf_inq_attname(ncid,ivar,iatt,attname)
          call check_cdf_status(status)
          status=nf_inq_att(ncid,ivar,attname,xtypeu,length)
          call check_cdf_status(status)
          write(*,'(a,a,a,i2,a,i6)') '    attribute: ',trim(attname),
     1     ' of type ',xtypeu,' and length ',length
          if(attname.eq.'units') then
            if(length.le.24) then
              status=nf_get_att_text(ncid,ivar,attname,units)
              call check_cdf_status(status)
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
          mds_signal_name='\TOP.TALLIES:TIME'
          write(*,'(a,a,a,a)') 'We need to save ',trim(varname),
     1     ' to ',trim(mds_signal_name)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 )
          status=mdsput2_d1d(
     *     trim(mds_signal_name)//char(0),
     *     'build_with_units($,"s")'//char(0),
     *     descr2_i1(idtype_double,dimsize(1),0),
     *     data1,0)
#else
          status=mdsput2(
     *     trim(mds_signal_name)//char(0),
     *     'build_with_units($,"s")'//char(0),
     *     descr2(idtype_double,dimsize(1),0),data1,0)
#endif
          write(*,*) 'mdsput2'
          call mds_err_msg(status)
          deallocate(data1)
        elseif(ndims.eq.3) then
          mds_signal_name='\TOP.TALLIES:' // trim(varname)
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
     *       '\TOP.TALLIES:' // trim(dimname(1)) // ',' //
     *       '\TOP.TALLIES:' // trim(dimname(2)) // ',' //
     *       '\TOP.TALLIES:TIME)'
            write(*,'(a,a)') 'Using TDI ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 )
            status=mdsput2_d3d(
     *       trim(mds_signal_name)//char(0),
     *       tdi//char(0),
     *       descr2_i3(idtype_double,dimsize(1),dimsize(2),
     1       dimsize(3),0),data3,0)
#else
            status=mdsput2(
     *       trim(mds_signal_name)//char(0),
     *       tdi//char(0),
     *       descr2(idtype_double,dimsize(1),dimsize(2),
     1       dimsize(3),0),data3,0)
#endif
            write(*,*) 'mdsput2'
            call mds_err_msg(status)
          else
            allocate(data3T(dimsize(3),dimsize(1),dimsize(2)))
            forall(i1=1:dimsize(1),i2=1:dimsize(2),i3=1:dimsize(3))
     *       data3T(i3,i1,i2)=data3(i1,i2,i3)
            tdi=
     *       'build_signal(build_with_units($,"' //
     *       trim(units) //
     *       '"),,\TOP.TALLIES:TIME,' //
     *       '\TOP.TALLIES:' // trim(dimname(1)) // ',' //
     *       '\TOP.TALLIES:' // trim(dimname(2)) // ')'
            write(*,'(a,a)') 'Using TDI ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 )
            status=mdsput2_d3d(
     *       trim(mds_signal_name)//char(0),
     *       tdi//char(0),
     *       descr2_i3(idtype_double,dimsize(3),dimsize(1),
     1       dimsize(2),0),data3T,0)
#else
            status=mdsput2(
     *       trim(mds_signal_name)//char(0),
     *       tdi//char(0),
     *       descr2(idtype_double,dimsize(3),dimsize(1),
     1       dimsize(2),0),data3T,0)
#endif
            write(*,*) 'mdsput2'
            call mds_err_msg(status)
            deallocate(data3T)
          endif
          deallocate(data3)
        elseif(ndims.eq.2) then
          mds_signal_name='\TOP.TALLIES:' // trim(varname)
          underscore_pos=index(mds_signal_name,'_')
          do while(underscore_pos.gt.0)
            mds_signal_name(underscore_pos:len(mds_signal_name))=
     1       mds_signal_name(underscore_pos+1:len(mds_signal_name))
            underscore_pos=index(mds_signal_name,'_')
          enddo
          write(*,'(a,a,a,a)') 'We need to save ',trim(varname),
     1     ' to ',trim(mds_signal_name)
          if(xtypev.eq.nf_char) then
! this is a 1d array of text strings and not time-dependent
            allocate(data2c(dimsize(1),dimsize(2)))
            status=nf_get_var_text(ncid,ivar,data2c)
            call check_cdf_status(status)
#ifdef F2003
            allocate(character(len=dimsize(1)+1) :: str1(dimsize(2)))
#else
            call xertst(dimsize(1).le.strmax,
     &           'increase size of strmax in b2cdftallies2mds')
            allocate(str1(dimsize(2)))
#endif
            do i2=1, dimsize(2)
              str1(i2)=''
              do i1=1,dimsize(1)
                str1(i2)(i1:i1) = data2c(i1,i2)
              enddo
              write(*,*) trim(str1(i2))
              str1(i2) = trim(str1(i2)) // char(0)
            enddo
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 )
            status=mdsput2_str1(
     *       trim(mds_signal_name)//char(0),
     *       '$'//char(0),
     *       descr2_str1(idtype_cstring,dimsize(2),0,len(str1(1))),
     *       str1,0)
#else
            status=mdsput2(
     *       trim(mds_signal_name)//char(0),
     *       '$'//char(0),
     *       descr2(idtype_cstring,dimsize(2),0,len(str1(1))),str1,0)
#endif
            write(*,*) 'mdsput2'
            call mds_err_msg(status)
            deallocate(data2c,str1)
          elseif(xtypev.eq.nf_double) then
            allocate(data2(dimsize(1),dimsize(2)))
            status=nf_get_var_double(ncid,ivar,data2)
            call check_cdf_status(status)
            if(.not.want_transpose) then
              tdi=
     *         'build_signal(build_with_units($,"' //
     *         trim(units) //
     *         '"),,' //
     *         '\TOP.TALLIES:' // trim(dimname(1)) //',' //
     *         '\TOP.TALLIES:TIME)'
              write(*,'(a,a)') 'Using TDI ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 )
              status=mdsput2_d2d(
     *         trim(mds_signal_name)//char(0),
     *         tdi//char(0),
     *         descr2_i2(idtype_double,dimsize(1),dimsize(2),0),
     *         data2,0)
#else
              status=mdsput2(
     *         trim(mds_signal_name)//char(0),
     *         tdi//char(0),
     *         descr2(idtype_double,dimsize(1),dimsize(2),0),data2,0)
#endif
              write(*,*) 'mdsput2'
              call mds_err_msg(status)
            else
              allocate(data2T(dimsize(2),dimsize(1)))
              data2T=transpose(data2)
              tdi=
     *         'build_signal(build_with_units($,"' //
     *         trim(units) //
     *         '"),,\TOP.TALLIES:TIME,' //
     *         '\TOP.TALLIES:' // trim(dimname(1)) // ')'
              write(*,'(a,a)') 'Using TDI ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 )
              status=mdsput2_d2d(
     *         trim(mds_signal_name)//char(0),
     *         tdi//char(0),
     *         descr2_i2(idtype_double,dimsize(2),dimsize(1),0),
     *         data2T,0)
#else
              status=mdsput2(
     *         trim(mds_signal_name)//char(0),
     *         tdi//char(0),
     *         descr2(idtype_double,dimsize(2),dimsize(1),0),data2T,0)
#endif
              write(*,*) 'mdsput2'
              call mds_err_msg(status)
              deallocate(data2T)
            endif
            deallocate(data2)
          else
            write(*,'(a,i2,a)')
     *       "Dim 2 case for type ",xtypev," not coded yet"
          endif
        else
          write(*,'(a,i2,a)')
     *     'Case for ',ndims,' dimensions not yet coded'
        endif
      enddo

      status=nf_close(ncid)
      call check_cdf_status(status)
#else
      write(*,*)
     . 'b2cdftallies2mds failed because compiled with NO_CDF option'
#endif
#else
      call xerrab('b2cdftallies2mds: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      call subend ()
      return
      end subroutine b2cdftallies2mds
