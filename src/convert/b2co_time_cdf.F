      subroutine b2co_time_cdf(nx, ny, ns, tim)
      use b2mod_types
      use b2mod_geo
      use b2mod_mwti
      use b2mod_indirect
      use b2mod_cdf
      use b2mod_subsys
      implicit none
      integer nx, ny, ns
      real (kind=R8) :: tim
#ifdef WG_TODO
#ifndef NO_CDF
#     include <netcdf.inc>
#endif
      integer :: jxi, jxa, jsep, ixtl, ixtr
      integer :: iylstrt, iyistrt, iyastrt, iyrstrt, iytlstrt, iytrstrt,
     1           iylend,  iyiend,  iyaend,  iyrend,  iytlend,  iytrend,
     2           nybl, nybr, nytl, nytr, nya, nyi, nc, nbatch,
     3           timedim, batchdim, nastep, ntstep
      integer, save :: write_2d = 0, target_offset = 1, ntim_batch = 500
      real (kind=R8), save :: dtim = 1.0_R8
#ifndef NO_CDF
#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
      integer :: status
      integer :: idim, ivar, iatt, iret, ib, it, ih,
     1 itim, ibatch, ncid, ncid2, ndims, nvars, natts, unlimid, xtype,
     2 dimids(NF_MAX_VAR_DIMS), len, dimsize(NF_MAX_VAR_DIMS),
     3 start(NF_MAX_VAR_DIMS), count(NF_MAX_VAR_DIMS)
      real (kind=R8) :: scale
      real (kind=R8), allocatable :: data1(:), data2(:,:), data3(:,:,:)
      character*(NF_MAX_NAME) :: varname, dimname, attname
      character :: units*(24)
      character*4 :: extension
      character*32 hlp_frm
      character*256 filename, plasmastate, message
      logical :: exists
      integer :: lnblnk
      external lnblnk
#endif
      external ipgeti, ipgetr, xerrab, xertst
#endif

      call subini ('b2co_time_cdf')
      write(*,*) 'b2co_time_cdf ...'

#ifdef WG_TODO

      call get_jsep(nx,ny,jxi,jxa,jsep)
      call ipgetr ('b2mndr_dtim', dtim)
      call ipgeti ('b2mwti_2dwrite', write_2d)
      call ipgeti ('b2mndt_av_ntim_batch', ntim_batch)
      call ipgeti ('b2mwti_target_offset', target_offset)
      call output_ds(crx,cry,nx,ny, -1,+target_offset,jsep,
     .                                        iylstrt,iylend,'dsl')
      call output_ds(crx,cry,nx,ny,jxi,0,jsep,iyistrt,iyiend,'dsi')
      call output_ds(crx,cry,nx,ny,jxa,0,jsep,iyastrt,iyaend,'dsa')
      call output_ds(crx,cry,nx,ny, nx,-target_offset,jsep,
     .                                        iyrstrt,iyrend,'dsr')
#ifndef NO_CDF
      if (nnreg(0).eq.8) then
        ixtl = 0
        do while (rightix(ixtl,max(topcut(1),topcut(2))).ne.nx+1.and.
     .            ixtl.lt.nx)
          ixtl=ixtl+1
        enddo
        ixtr = ixtl
        do while (leftix(ixtr,max(topcut(1),topcut(2))).ne.-2 .and.
     .            ixtr.lt.nx)
          ixtr=ixtr+1
        enddo
        call output_ds(crx,cry,nx,ny,ixtl,-target_offset,jsep,
     .                                        iytlstrt,iytlend,'dstl')
        call output_ds(crx,cry,nx,ny,ixtr,+target_offset,jsep,
     .                                        iytrstrt,iytrend,'dstr')
        nytl = iytlend - iytlstrt + 1
        nytr = iytrend - iytrstrt + 1
      else
        nytl = 0
        nytr = 0
      endif
      nybl = iylend  - iylstrt  + 1
      nybr = iyrend  - iyrstrt  + 1
      nya  = iyaend  - iyastrt  + 1
      nyi  = iyiend  - iyistrt  + 1
      nc = max(nncut,1)
      nbatch = int((tim/dtim)/ntim_batch) + 1
c access the b2time.nc file
      filename='b2time.nc'
      call find_file(filename,exists)
      status=nf_open(filename,NF_NOWRITE,ncid)
      call check_cdf_status(status)
      status=nf_inq(ncid,ndims,nvars,natts,unlimid)
      call check_cdf_status(status)
      itim=0
      ibatch=0
      nastep=0
      ntstep=0
      timedim=-1
      batchdim=-1

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
        if(varname.eq.'ntstep') then
          allocate(data1(1))
          status=nf_get_var_double(ncid,ivar,data1)
          call check_cdf_status(status)
          ntstep = nint(data1(1))
          deallocate(data1)
        else if(varname.eq.'nastep') then
          allocate(data1(1))
          status=nf_get_var_double(ncid,ivar,data1)
          call check_cdf_status(status)
          nastep = nint(data1(1))
          deallocate(data1)
        endif
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

c find the time and batch vector sizes
      do ivar=1,nvars
        status=nf_inq_var(ncid,ivar,varname,xtype,
     1   ndims,dimids,natts)
        call check_cdf_status(status)
        do idim=1,ndims
          status=nf_inq_dim(ncid,dimids(idim),dimname,len)
          call check_cdf_status(status)
          dimsize(idim)=len
        enddo
        if(varname.eq.'timesa'.and.ntstep.gt.0) then
          allocate(data1(dimsize(1)))
          status=nf_get_var_double(ncid,ivar,data1)
          call check_cdf_status(status)
          timedim=dimsize(1)
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
        else if(varname.eq.'batchsa'.and.nastep.gt.0) then
          allocate(data1(dimsize(1)))
          status=nf_get_var_double(ncid,ivar,data1)
          call check_cdf_status(status)
          batchdim=dimsize(1)
          ib=1
          it=batchdim
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
          ibatch=ih
          deallocate(data1)
        end if
      end do

c prepare the b2time.nc.xxxx file
#ifdef NO_GETENV
      extension='xxxx'
#else
#ifdef USE_PXFGETENV
      CALL PXFGETENV ('PLASMASTATE', 0, plasmastate, lenval, ierror)
#else
      call get_environment_variable('PLASMASTATE',
     . status=ierror, length=lenval)
      if (ierror.eq.0) then
        call get_environment_variable('PLASMASTATE',value=plasmastate)
        extension=plasmastate(lenval-3:lenval)
      elseif (ierror.eq.1) then
        extension='xxxx'
      endif
#endif
      extension=plasmastate(lnblnk(plasmastate)-3:lnblnk(plasmastate))
#endif
      filename='b2time.nc.'//extension
      message = 'File '//trim(filename)//' already present! '//
     .          'Please remove first.'
      call find_file(filename,exists)
      call xertst(.not.exists,message)
      status=nf_create(trim(filename),NF_NOCLOBBER,ncid2)
      call check_cdf_status(status)
      call b2crtimecdf(filename, nx, ny, nybl, nytl, nytr, nybr,
     .                 nya, nyi, nc, ns, nbatch, write_2d, ncid2,
     .                 .false., iret)
      status=nf_open(trim(filename),NF_WRITE,ncid2)
      call check_cdf_status(status)

c read and rewrite data
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
        if (ndims.eq.0) then
          write(message,'(a)') 'Processing '//trim(varname)//' ...'
        elseif (ndims.eq.1) then
          write(message,'(a,i6,a)') 'Processing '//trim(varname)//
     .     ' with size ',len,' ...'
        else
          write(hlp_frm,'(a,i1,a)') '(3a,',ndims-1,'(i6,'',''),i6,a1)'
          write(message,hlp_frm) 'Processing ',trim(varname),
     .     ' with dimensions (',(count(idim),idim=1,ndims),')'
        end if
        write(*,'(a)') trim(message)
        if(varname.eq.'timesa') then
          if(itim.gt.0) then
            allocate(data1(dimsize(1)))
            status=nf_get_vara_double(ncid,ivar,start,count,data1)
            call check_cdf_status(status)
            count(1)=itim
            status=nf_put_vara_double(ncid2,ivar,start,count,data1)
            deallocate(data1)
          end if
        else if(varname.eq.'batchsa') then
          if(ibatch.gt.0) then
            allocate(data1(dimsize(1)))
            status=nf_get_vara_double(ncid,ivar,start,count,data1)
            call check_cdf_status(status)
            count(1)=ibatch
            status=nf_put_vara_double(ncid2,ivar,start,count,data1)
            deallocate(data1)
          endif
        else
          if(ndims.eq.3) then
            allocate(data3(dimsize(1),dimsize(2),dimsize(3)))
            status=nf_get_vara_double(ncid,ivar,start,count,data3)
            call check_cdf_status(status)
            if (dimsize(3).eq.timedim) then
              count(3)=itim
              if (itim.gt.0) status=nf_put_vara_double(ncid2,ivar,
     .                                                start,count,data3)
            else if (dimsize(3).eq.batchdim) then
              count(3)=ibatch
              if (ibatch.gt.0) status=nf_put_vara_double(ncid2,ivar,
     .                                                start,count,data3)
            else
              count(3)=dimsize(3)
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
              if (itim.gt.0) status=nf_put_vara_double(ncid2,ivar,
     .                                            start,count,data2)
            else if (dimsize(2).eq.batchdim) then
              count(2)=ibatch
              if (ibatch.gt.0) status=nf_put_vara_double(ncid2,ivar,
     .                                            start,count,data2)
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
              if (itim.gt.0) status=nf_put_vara_double(ncid2,ivar,
     .                                           start,count,data1)
            else if (dimsize(1).eq.batchdim) then
              count(1)=ibatch
              if (ibatch.gt.0) status=nf_put_vara_double(ncid2,ivar,
     .                                             start,count,data1)
            else
              status=nf_put_vara_double(ncid2,ivar,start,count,data1)
            end if
            call check_cdf_status(status)
            deallocate(data1)
          elseif(ndims.eq.0) then
            allocate(data1(1))
            status=nf_get_var_double(ncid,ivar,data1)
            call check_cdf_status(status)
            if (nint(data1(1)).eq.timedim) then
              data1(1) = real(itim)
            else if (nint(data1(1)).eq.batchdim) then
              data1(1) = real(ibatch)
            end if
            status=nf_put_var_double(ncid2,ivar,data1)
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
     1 'b2co_time_cdf failed because compiled with NO_CDF option'
#endif
#else
      write(*,*)
     1 'b2co_time_cdf failed, not yet converted to WG'
#endif
      call subend ()
      return
      end subroutine b2co_time_cdf
