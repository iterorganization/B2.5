      subroutine b2cdf2mds
      use b2mod_subsys
#ifdef MDSPLUS
      use b2mod_types
      implicit none

      interface
        integer function descr2(idtype,il1,il2,iflag)
        integer, intent(in) :: idtype
        integer, intent(in) :: il1
        integer, intent(in), optional :: il2
        integer, intent(in) :: iflag
        end function descr2
      end interface

#     include "mdslib_b2.inc"
      integer mdsput2
#ifndef NO_CDF
#     include <netcdf.inc>
#endif
      logical :: want_transpose, found_time
      real (kind=R8) :: ds(2048), da(2048,7:10)
      integer :: status, len_da(7:10)
      integer :: i, il
      logical :: exists
#ifndef NO_CDF
      integer :: idim, ivar, i1, i2, iatt,
     1 ncid, ndims, nvars, natts, unlimid, xtype,
     2 dimids(NF_MAX_VAR_DIMS), len, index, dimsize(NF_MAX_VAR_DIMS)
      real (kind=R8) :: scale
      real (kind=R8), allocatable :: data1(:), data2(:,:), data2T(:,:)
      character*(NF_MAX_NAME) :: varname, dimname, attname
      character*8 :: mds_zone(6)
      character*2 :: cdf_suffix(6)
      character :: mds_signal_name*(64), units*(24),
     1 tdi*(172), mds_signal_prefix*(8), mds_signal_suffix*(8)
      character*256 filename
#endif

      character*4 :: ds_filename(11)
      character*24 :: ds_mds(11)
      data ds_filename /
     * 'dsl','dsi','dsa','dsr','dstl','dstr',
     1 'dsL','dsR','dsTL','dsTR','dsp'/
      data ds_mds/
     1 '\TOP.TIMEDEP.TARGET1:DS', '\TOP.TIMEDEP.IMP:DS',
     1 '\TOP.TIMEDEP.OMP:DS',     '\TOP.TIMEDEP.TARGET4:DS',
     1 '\TOP.TIMEDEP.TARGET2:DS', '\TOP.TIMEDEP.TARGET3:DS',
     1 '\TOP.TIMEDEP.TARGET1:DA', '\TOP.TIMEDEP.TARGET4:DA',
     1 '\TOP.TIMEDEP.TARGET2:DA', '\TOP.TIMEDEP.TARGET3:DA',
     1 '\TOP.TIMEDEP.INPUT:DS'/
#ifndef NO_CDF
      data cdf_suffix /'l','i','a','r','tl','tr'/
      data mds_zone /'TARGET1','IMP','OMP','TARGET4',
     1 'TARGET2','TARGET3'/
#endif
      external xerrab, xertst

      call subini ('b2cdf2mds')
      write(*,*) 'b2cdf2mds ...'

      want_transpose=.false.
      want_transpose=.true.

c first handle x-coordinate of profiles (and areas)
      len_da=0
      do i=1,11
        call find_file(ds_filename(i),exists)
        if(exists) then
          open(1,file=trim(ds_filename(i)),status='old')
          do il=1,2048
            read(1,*,end=100) ds(il)
          enddo
          call xerrab ('Increase length of ds in b2cdf2mds')
  100     continue
          close(1)
          il=il-1

          write(*,'(a,a,a,i6,a,a)') 'Putting ',trim(ds_filename(i)),
     1     ' with length ',il,' to ',trim(ds_mds(i))
          status=mdsput2(
     *     trim(ds_mds(i))//char(0),
     *     'build_with_units($,"m")'//char(0),
     *     descr2(idtype_double,il,iflag=0),ds,
     *     0)
          write(*,*) 'mdsput2'
          call mds_err_msg(status)
          if(i.gt.6.and.i.le.10) then
            da(1:il,i)=ds(1:il)
            len_da(i)=il
          endif
        endif
      enddo

#ifndef NO_CDF
c access the b2time.nc file
      filename='b2time.nc'
      call find_file(filename,exists)
      status=nf_open(filename,NF_NOWRITE,ncid)
      call check_cdf_status(status)

      status=nf_inq(ncid,ndims,nvars,natts,unlimid)
      call check_cdf_status(status)

c loop over all variables in the file
      do ivar=1,nvars
        found_time = .false.
        status=nf_inq_var(ncid,ivar,varname,xtype,
     1   ndims,dimids,natts)
        call check_cdf_status(status)
        write(*,'(a,a,a,i2,a,i2,a)')
     .   'variable: ',trim(varname),' has type ',xtype,
     1   ' and ',ndims,' dimensions'
        do idim=1,ndims
          status=nf_inq_dim(ncid,dimids(idim),dimname,len)
          call check_cdf_status(status)
          dimsize(idim)=len
          write(*,'(a,a,a,i6)')
     .     '    dimension: ',trim(dimname),' of length ',len
          found_time = found_time .or. (trim(dimname) .eq. 'time')
        enddo
        units=''
        scale=1
        do iatt=1,natts
          status=nf_inq_attname(ncid,ivar,iatt,attname)
          call check_cdf_status(status)
          status=nf_inq_att(ncid,ivar,attname,xtype,len)
          call check_cdf_status(status)
          write(*,'(a,a,a,i2,a,i6)')
     .     '    attribute: ',trim(attname),
     1     ' of type ',xtype,' and length ',len
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
        write(*,*) '  units = ',trim(units)
        write(*,*) '  scale = ',scale

        if(.not. found_time) then
          write(*,*) 'This variable did not have a ',
     &               'time dimension --- it has been skipped'
          cycle
        endif

c first deal with the time vector
        if(varname.eq.'timesa') then
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
          mds_signal_name='\TOP.TIMEDEP:TIME'
          write(*,'(a,a,a,a)') 'We need to save ',trim(varname),
     1     ' to ',trim(mds_signal_name)
          status=mdsput2(
     *     trim(mds_signal_name)//char(0),
     *     'build_with_units($,"s")'//char(0),
     *     descr2(idtype_double,dimsize(1),iflag=0),data1,
     *     0)
          write(*,*) 'mdsput2'
          call mds_err_msg(status)
          deallocate(data1)
        else
c conceptually two types of traces, f(t) and g(t,x); here handle g(t,x)
          if(varname(3:4).eq.'3d') then
            index=0
            do i=1,6
              if(varname(5:6).eq.cdf_suffix(i)) then
                index=i
                exit
              endif
            enddo
            call xertst (index.gt.0, '3d suffix not found')
            mds_signal_name='\TOP.TIMEDEP.' //
     1       trim(mds_zone(index)) // ':'
            write(*,'(a,a,a,a)') 'We need to save ',trim(varname),
     1       ' to ',trim(mds_signal_name) // varname(1:2)
            allocate(data2(dimsize(1),dimsize(2)))
            status=nf_get_var_double(ncid,ivar,data2)
            call check_cdf_status(status)
            data2=data2*scale
            if(varname(1:1).eq.'f') then
              call xertst (index.ne.2.and.index.ne.3,
     1         'Invalid fluxes at midplane')
              if(index.gt.3) index=index-2
              index=index+6
              call xertst(dimsize(1).eq.len_da(index),
     1         'dimsize/len_da mismatch')
              forall (i1=1:dimsize(1),i2=1:dimsize(2))
     1         data2(i1,i2)=data2(i1,i2)/
     1         da(i1,index)
              units=trim(units)  // '/m^2'
            endif
            if(.not.want_transpose) then
              tdi=
     *         'build_signal(build_with_units($,"' //
     *         trim(units) //
     *         '"),,' //
     *         trim(mds_signal_name) //
     *         'DS,' //
     *         '\TOP.TIMEDEP:TIME)'
              write(*,'(a,a)') 'Using TDI ',trim(tdi)
              status=mdsput2(
     *         trim(mds_signal_name)//varname(1:2)//char(0),
     *         tdi//char(0),
     *         descr2(idtype_double,dimsize(1),dimsize(2),0),data2,
     *         0)
              write(*,*) 'mdsput2'
              call mds_err_msg(status)
            else
              allocate(data2T(dimsize(2),dimsize(1)))
              data2T=transpose(data2)
              tdi=
     *         'build_signal(build_with_units($,"' //
     *         trim(units) //
     *         '"),,\TOP.TIMEDEP:TIME,' //
     *         trim(mds_signal_name) //
     *         'DS)'
              write(*,'(a,a)') 'Using TDI ',trim(tdi)
              status=mdsput2(
     *         trim(mds_signal_name)//varname(1:2)//char(0),
     *         tdi//char(0),
     *         descr2(idtype_double,dimsize(2),dimsize(1),0),data2T,
     *         0)
              write(*,*) 'mdsput2'
              call mds_err_msg(status)
              deallocate(data2T)
            endif
            deallocate(data2)
          else
c and now for f(t)
            mds_signal_name=''
            select case(varname(1:1))                       ! varname(1:1)
            case('f')                                       ! fluxes
              mds_signal_prefix=varname(1:3)
              select case(varname(4:7))                     ! varname(4:7)
              case('xip')
                mds_signal_suffix='TARG1INT'
              case('xap')
                mds_signal_suffix='TARG4INT'
              case('yip')
                mds_signal_suffix='DIVINT'
              case('yap')
                mds_signal_suffix='MCWINT'
              case('sip')
                mds_signal_suffix=''
              case('sap')
                mds_signal_suffix=''
              case('sipp')
                mds_signal_suffix=''
              case('sapp')
                mds_signal_suffix=''
              case default                                  ! varname(4:7)
                mds_signal_suffix=''
              end select                                    ! varname(4:7)
              if(mds_signal_suffix.ne.'') then
                mds_signal_name='\TOP.TIMEDEP.' //
     1           trim(mds_signal_prefix) // ':' //
     2           trim(mds_signal_suffix)
              else
                mds_signal_name=''
              endif
            case default                                    ! varname(1:1)
              select case(varname(1:2))                     ! varname(1:2)
              case('te','ti','tp','po','ne',
     1         'dn','dp','ke','ki','vx','vy','vs')
                select case(varname(3:6))                   ! varname(3:6)
                case('sepi')
                  mds_signal_suffix='TARG1SEP'
                case('sepm')
                  mds_signal_suffix='OMPSEP'
                case('sepa')
                  mds_signal_suffix='TARG4SEP'
                case('mxip')
                  mds_signal_suffix='TARG1MAX'
                case('mxap')
                  mds_signal_suffix='TARG4MAX'
                case default                                ! varname(3:6)
                  mds_signal_suffix=''
                end select                                  ! varname(3:6)
                if(mds_signal_suffix.ne.'') then
                  mds_signal_name='\TOP.TIMEDEP.' //
     1             varname(1:2) // ':' // mds_signal_suffix
                else
                  mds_signal_name=''
                endif
              case('tm')                                    ! varname(1:2)
                select case(varname(3:4))                   ! varname(3:4)
                case('ha')
                  mds_signal_name='\TOP.TIMEDEP:HAL' //
     1             varname(5:len_trim(varname))
                case default
                  mds_signal_name=''
                end select
              case default                                  ! varname(1:2)
                mds_signal_name=''
              end select                                    ! varname(1:2)
            end select                                      ! varname(1:1)
            if(mds_signal_name.ne.'') then
              if(ndims.eq.2) then
                write(*,*) '2d case'
                allocate(data2(dimsize(1),dimsize(2)))
                status=nf_get_var_double(ncid,ivar,data2)
                call check_cdf_status(status)
                allocate(data2T(dimsize(2),dimsize(1)))
                data2T=transpose(data2)
                write(*,'(a,a,a,a)')
     *           'We need to save ',trim(varname),
     1           ' to ',trim(mds_signal_name)
                tdi=
     *           'build_signal(build_with_units($,"' //
     *           trim(units) //
     *           '"),,\TOP.TIMEDEP:TIME)'
                write(*,'(a,a)') 'Using TDI ',trim(tdi)
                status=mdsput2(
     *           trim(mds_signal_name)//char(0),
     *           trim(tdi)//char(0),
     *           descr2(idtype_double,dimsize(2),iflag=0),data2T(:,1),
     *           0)
                write(*,*) 'mdsput2'
                call mds_err_msg(status)
                deallocate(data2T)
                deallocate(data2)
              elseif(ndims.eq.1) then
                write(*,*) '1d case'
                allocate(data1(dimsize(1)))
                status=nf_get_var_double(ncid,ivar,data1)
                call check_cdf_status(status)
                write(*,'(a,a,a,a)')
     *           'We need to save ',trim(varname),
     1           ' to ',trim(mds_signal_name)
                tdi=
     *           'build_signal(build_with_units($,"' //
     *           trim(units) //
     *           '"),,\TOP.TIMEDEP:TIME)'
                write(*,'(a,a)') 'Using TDI ',trim(tdi)
                status=mdsput2(
     *           trim(mds_signal_name)//char(0),
     *           trim(tdi)//char(0),
     *           descr2(idtype_double,dimsize(1),iflag=0),data1,
     *           0)
                write(*,*) 'mdsput2'
                call mds_err_msg(status)
                deallocate(data1)
              else
                write(*,'(a,i2,a)')
     *           'Case for ',ndims,' dimensions not yet coded'
              endif
            else
              write(*,'(a,a,a)')
     *         'Case for ',trim(varname),' not coded yet'
            endif
          endif
        endif
      enddo

      status=nf_close(ncid)
      call check_cdf_status(status)
#else
      write(*,*) 'b2cdf2mds failed because compiled with NO_CDF option'
#endif
#else
      call xerrab('b2cdf2mds: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      call subend ()
      return
      end subroutine b2cdf2mds
