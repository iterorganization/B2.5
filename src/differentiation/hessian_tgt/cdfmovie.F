      subroutine cdfmovie
#ifndef NO_CDF
     . (ncid,nCv,ns,geo,state,switch)
      use b2mod_types
      use b2mod_time
      use b2us_geo_diffv_diffv
      use b2us_plasma_diffv_diffv
      use b2mod_dimensions
      use b2mod_movies
      use b2mod_constants
#ifdef B25_EIRENE
      use eirmod_parmmod
     ,   , ncv_eir => ncv
      use eirmod_wneutrals
      use eirmod_extrab25
#else
      use b2mod_eirdiag
#endif
      use b2mod_switches_diffv_diffv
      use b2mod_subsys
      implicit none
      integer ncid, nCv, ns
      type (geometry), intent(in) :: geo
      type (B2State), intent(in) :: state
      type (switches), intent(in) :: switch
#     include <netcdf.inc>
#ifndef B25_EIRENE
#     include <COUPLE/KOPPLDIM.F>
#endif
      real (kind=R8) :: time(1), dummy(1:nCv), dummy3(1:nCv, 1:ns)
      integer iCv, is
      integer ncall
      integer status, ndims, nvars, natts, unlimid
      character*256 filename
      logical ex
      save ncall
      data ncall/0/

      call subini ('cdfmovie')
      if (ncall.eq.0) then
        filename='b2movies.nc'
        call find_file(filename,ex)
        if (ex.and.switch%b2mndr_stim.lt.0.0_R8) then
          write(6,'(a)') trim(filename)//' will be appended'
          status=nf_open(filename,or(NF_WRITE,NF_SHARE),ncid)
          call check_cdf_status(status)
          status=nf_inq(ncid,ndims,nvars,natts,unlimid)
          call check_cdf_status(status)
        else
          if (ex) then
            write(6,'(a)') trim(filename)//' will be replaced'
          else
            write(6,'(a)') trim(filename)//' will be created'
          end if
          call create_cdf_start(ncid, filename, nCv, ns, nnatmi, nnmoli)
          call create_cdf_field2(ncid, 'te',
     1     'Electron temperature', 'eV')
          call create_cdf_field2(ncid, 'ti',
     1     'Ion temperature', 'eV')
          call create_cdf_field2(ncid, 'ne',
     1     'Electron density', 'm^-3')
          call create_cdf_field3(ncid, 'na', 1,
     1     'Ion density', 'm^-3')
          call create_cdf_field2(ncid, 'rqahesum',
     1     'Electron cooling rate', 'W m^-3')
          call create_cdf_field2(ncid, 'rqradsum',
     1     'Radiation rate', 'W m^-3')
          if(nnatmi.gt.0) then
            call create_cdf_field3(ncid, 'dab2', 2,
     1       'Atom density', 'm^-3')
          endif
          if(nnmoli.gt.0) then
            call create_cdf_field3(ncid, 'dmb2', 3,
     1       'Molecular density', 'm^-3')
          endif
          call create_cdf_end(ncid)
        endif
      endif

      time(1) = tim

      write(*,*) 'put_cdf_start'
      call put_cdf_start(ncid,filename,time)
      write(*,*) 'put_cdf_field2: te'
      dummy(:) = state%pl%te(:)/qe
      call put_cdf_field2(ncid,'te',dummy,nCv)
      write(*,*) 'put_cdf_field2: ti'
      dummy(:) = state%pl%ti(:)/qe
      call put_cdf_field2(ncid,'ti',dummy,nCv)
      write(*,*) 'put_cdf_field2: ne'
      call put_cdf_field2(ncid,'ne',state%dv%ne,nCv)
      write(*,*) 'put_cdf_field3: na'
      call put_cdf_field3(ncid,'na',state%pl%na,nCv,ns)
      do iCv = 1, nCv
        dummy(iCv) = 0.0_R8
        do is = 0, ns-1
          dummy(iCv) = dummy(iCv) + state%srw%rqahe(iCv,is)
        enddo
        dummy(iCv) = dummy(iCv) / geo%cvVol(iCv)
      enddo
      write(*,*) 'put_cdf_field2: rqahesum'
      call put_cdf_field2(ncid,'rqahesum',dummy,nCv)
      do iCv = 1, nCv
        dummy(iCv) = 0.0_R8
        do is = 0, ns-1
          dummy(iCv) = dummy(iCv) + state%srw%rqrad(iCv,is)
        enddo
        dummy(iCv) = dummy(iCv) / geo%cvVol(iCv)
      enddo
      write(*,*) 'put_cdf_field2: rqradsum'
      call put_cdf_field2(ncid,'rqradsum',dummy,nCv)
      if(nnatmi.gt.0) then
        dummy3 = 0.0_R8
        dummy3(1:nCv,1:nnatmi) = dab2(1:nCv,1:nnatmi,1)
        write(*,*) 'put_cdf_field3: dab2'
        call put_cdf_field3(ncid,'dab2',dummy3,nCv,nnatmi)
      endif
      if(nnmoli.gt.0) then
        dummy3 = 0.0_R8
        dummy3(1:nCv,1:nnmoli) = dmb2(1:nCv,1:nnmoli,1)
        write(*,*) 'put_cdf_field3: dmb2'
        call put_cdf_field3(ncid,'dmb2',dummy3,nCv,nnmoli)
      endif
      write(*,*) 'put_cdf_end'
      call put_cdf_end(ncid)

* ..return
      ncall = ncall+1
      call subend ()
#endif
      return
      end subroutine cdfmovie
