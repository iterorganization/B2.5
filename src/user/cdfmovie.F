      subroutine cdfmovie(ncid,nx,ny,ns)
      use b2mod_types
      use b2mod_plasma
      use b2mod_sources
      use b2mod_geo
      use b2mod_time
#ifndef NO_CDF
      use b2mod_movies
#endif
      use b2mod_constants
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_extrab25
#else
      use b2mod_eirdiag
#endif
      implicit none
#ifndef B25_EIRENE
#     include <COUPLE/KOPPLDIM.F>
#endif
      integer ncid, nx, ny, ns
#ifndef NO_CDF
      real (kind=R8) :: time(1), dummy(-1:nx, -1:ny), 
     & dummy3(-1:nx, -1:ny, 1:ns)
      integer ix, iy, is
      integer ncall
      integer status, ndims, nvars, natts, unlimid
      real (kind=R8), save :: stim = 0.0_R8
      character*256 filename
      logical ex
      save ncall
      data ncall/0/

      call subini ('cdfmovie')
      if (ncall.eq.0) then
        filename='b2movies.nc'
        call find_file(filename,ex)
        call ipgetr ('b2mndr_stim', stim)
        if (ex.and.stim.lt.0.0_R8) then
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
          call create_cdf_start(ncid, filename, nx+2, ny+2, ns,
     &     nnatmi, nnmoli)
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
      dummy(:,:) = te(:,:)/qe
      call put_cdf_field2(ncid,'te',dummy,nx+2,ny+2)
      write(*,*) 'put_cdf_field2: ti'
      dummy(:,:) = ti(:,:)/qe
      call put_cdf_field2(ncid,'ti',dummy,nx+2,ny+2)
      write(*,*) 'put_cdf_field2: ne'
      call put_cdf_field2(ncid,'ne',ne,nx+2,ny+2)
      write(*,*) 'put_cdf_field3: na'
      call put_cdf_field3(ncid,'na',na,nx+2,ny+2,ns)
      do iy=-1,ny
        do ix=-1,nx
          dummy(ix,iy)=0.0_R8
          do is=0,ns-1
            dummy(ix,iy)=dummy(ix,iy)+rqahe(ix,iy,is)
          enddo
          dummy(ix,iy)=dummy(ix,iy)/vol(ix,iy)
        enddo
      enddo
      write(*,*) 'put_cdf_field2: rqahesum'
      call put_cdf_field2(ncid,'rqahesum',dummy,nx+2,ny+2)
      do iy=-1,ny
        do ix=-1,nx
          dummy(ix,iy)=0.0_R8
          do is=0,ns-1
            dummy(ix,iy)=dummy(ix,iy)+rqrad(ix,iy,is)
          enddo
          dummy(ix,iy)=dummy(ix,iy)/vol(ix,iy)
        enddo
      enddo
      write(*,*) 'put_cdf_field2: rqradsum'
      call put_cdf_field2(ncid,'rqradsum',dummy,nx+2,ny+2)
      if(nnatmi.gt.0) then
        dummy3=0.0_R8
        dummy3(-1:nx,-1:ny,1:nnatmi)=dab2(0:nx+1,0:ny+1,1:nnatmi,1)
        write(*,*) 'put_cdf_field3: dab2'
        call put_cdf_field3(ncid,'dab2',dummy3,nx+2,ny+2,nnatmi)
      endif
      if(nnmoli.gt.0) then
        dummy3=0.0_R8
        dummy3(-1:nx,-1:ny,1:nnmoli)=dmb2(0:nx+1,0:ny+1,1:nnmoli,1)
        write(*,*) 'put_cdf_field3: dmb2'
        call put_cdf_field3(ncid,'dmb2',dummy3,nx+2,ny+2,nnmoli)
      endif
      write(*,*) 'put_cdf_end'
      call put_cdf_end(ncid)

* ..return
      ncall = ncall+1
      call subend ()
#endif
      return
      end subroutine cdfmovie
