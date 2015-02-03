      subroutine cdfmovie(nx,ny,ns)
      use b2mod_types
      use b2mod_plasma
      use b2mod_sources
      use b2mod_geo
      use b2mod_time
#ifndef NO_CDF
      use b2mod_movies
#endif
      use b2mod_constants
#ifdef USE_EIRENE
      use eirmod_extrab25
#endif
      implicit none
#     include <COUPLE/KOPPLDIM.F>
      integer nx, ny, ns
#ifndef NO_CDF
      real (kind=R8) :: time(1), dummy(-1:nx, -1:ny), 
     & dummy3(-1:nx, -1:ny, 1:ns)
      integer ix, iy, is
      integer ncall
      save ncall
      data ncall/0/

      call subini ('cdfmovie')
      if (ncall.eq.0) then
	call create_cdf_start('b2movies.nc', nx+2, ny+2, ns, 
     &	 nnatmi, nnmoli)
	call create_cdf_field2('te',
     1   'Electron temperature', 'eV')
	call create_cdf_field2('ti',
     1   'Ion temperature', 'eV')
        call create_cdf_field2('ne',
     1   'Electron density', 'm^-3')
        call create_cdf_field3('na',1,
     1   'Ion density', 'm^-3')
        call create_cdf_field2('rqahesum',
     1   'Electron cooling rate', 'W m^-3')
        call create_cdf_field2('rqradsum',
     1   'Radiation rate', 'W m^-3')
	if(nnatmi.gt.0) then
	  call create_cdf_field3('dab2',2,
     1	   'Atom density', 'm^-3')
	endif
	if(nnmoli.gt.0) then
	  call create_cdf_field3('dmb2',3,
     1	   'Molecular density', 'm^-3')
	endif
	call create_cdf_end()
      endif

      time(1) = tim

      write(*,*) 'put_cdf_start'
      call put_cdf_start('b2movies.nc',time)
      write(*,*) 'put_cdf_field2: te'
      call put_cdf_field2('te',te/qe,nx+2,ny+2)
      write(*,*) 'put_cdf_field2: ti'
      call put_cdf_field2('ti',ti/qe,nx+2,ny+2)
      write(*,*) 'put_cdf_field2: ne'
      call put_cdf_field2('ne',ne,nx+2,ny+2)
      write(*,*) 'put_cdf_field3: na'
      call put_cdf_field3('na',na,nx+2,ny+2,ns)
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
      call put_cdf_field2('rqahesum',dummy,nx+2,ny+2)
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
      call put_cdf_field2('rqradsum',dummy,nx+2,ny+2)
      write(*,*) 'put_cdf_end'
      if(nnatmi.gt.0) then
	dummy3=0.0_R8
	dummy3(-1:nx,-1:ny,1:nnatmi)=dab2(0:nx+1,0:ny+1,1:nnatmi,1)
	write(*,*) 'put_cdf_field3: dab2'
	call put_cdf_field3('dab2',dummy3,nx+2,ny+2,nnatmi)
      endif
      if(nnmoli.gt.0) then
	dummy3=0.0_R8
	dummy3(-1:nx,-1:ny,1:nnmoli)=dmb2(0:nx+1,0:ny+1,1:nnmoli,1)
	write(*,*) 'put_cdf_field3: dmb2'
	call put_cdf_field3('dmb2',dummy3,nx+2,ny+2,nnmoli)
      endif
      call put_cdf_end()

* ..return
      ncall = ncall+1
      call subend ()
#endif
      return
      end subroutine cdfmovie
