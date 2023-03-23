      subroutine get_from_movie_cdf(ncidcdf,istep)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_movies
      use b2mod_plasma
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_extrab25
#else
      use b2mod_eirdiag
#endif
      use b2mod_subsys
      implicit none
#ifndef B25_EIRENE
#     include <COUPLE/KOPPLDIM.F>
#endif
      integer ncidcdf,istep
#ifndef NO_CDF
#     include <netcdf.inc>
      integer :: iret, start(4), count(4), varid
      real (kind=R8) :: time(1)
      real (kind=R8), allocatable :: dummy(:,:,:)
#endif
      external xertst

      call subini ('get_from_movie_cdf')
      call xertst (0.lt.ncidcdf, 'faulty input ncidcdf')
      call xertst (0.le.istep, 'faulty input istep')
#ifndef NO_CDF
      start(1)=istep
      count(1)=1
      write(*,*) 'Getting time for slice ',istep
      iret=nf_inq_varid(ncidcdf,'times',varid)
      call check_cdf_status(iret)
      iret=nf_get_vara_double(ncidcdf,varid,start,count,time)
      call check_cdf_status(iret)
      write(*,*) 'time = ',time(1)
      call reset_time(time(1))
      start(1)=1
      start(2)=1
      start(3)=istep
      count(1)=nx+2
      count(2)=ny+2
      count(3)=1
      write(*,*) 'Getting te for slice ',istep
      iret=nf_inq_varid(ncidcdf,'te',varid)
      call check_cdf_status(iret)
      iret=nf_get_vara_double(ncidcdf,varid,start,count,te)
      call check_cdf_status(iret)
      te=te*qe
      write(*,*) 'Getting ti for slice ',istep
      iret=nf_inq_varid(ncidcdf,'ti',varid)
      call check_cdf_status(iret)
      iret=nf_get_vara_double(ncidcdf,varid,start,count,ti)
      call check_cdf_status(iret)
      ti=ti*qe
      write(*,*) 'Getting ne for slice ',istep
      iret=nf_inq_varid(ncidcdf,'ne',varid)
      call check_cdf_status(iret)
      iret=nf_get_vara_double(ncidcdf,varid,start,count,ne)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=istep
      count(1)=nx+2
      count(2)=ny+2
      count(3)=ns
      count(4)=1
      write(*,*) 'Getting na for slice ',istep
      iret=nf_inq_varid(ncidcdf,'na',varid)
      call check_cdf_status(iret)
      iret=nf_get_vara_double(ncidcdf,varid,start,count,na)
      call check_cdf_status(iret)
      if(natmi.gt.0) then
        allocate(dummy(0:nx+1,0:ny+1,1:natmi))
        count(3)=natmi
        write(*,*) 'Getting dab2 for slice ',istep
        iret=nf_inq_varid(ncidcdf,'dab2',varid)
        call check_cdf_status(iret)
        iret=nf_get_vara_double(ncidcdf,varid,start,count,dab2)
        call check_cdf_status(iret)
!WG_TODOdab2(0:nx+1,0:ny+1,1:natmi,1)=dummy(0:nx+1,0:ny+1,1:natmi)
        deallocate(dummy)
      endif
      if(nmoli.gt.0) then
        allocate(dummy(0:nx+1,0:ny+1,1:nmoli))
        count(3)=nmoli
        write(*,*) 'Getting dmb2 for slice ',istep
        iret=nf_inq_varid(ncidcdf,'dmb2',varid)
        call check_cdf_status(iret)
        iret=nf_get_vara_double(ncidcdf,varid,start,count,dummy)
        call check_cdf_status(iret)
!WG_TODOdmb2(0:nx+1,0:ny+1,1:nmoli,1)=dummy(0:nx+1,0:ny+1,1:nmoli)
        deallocate(dummy)
      endif
#endif
      call subend ()
      return
      end subroutine get_from_movie_cdf
