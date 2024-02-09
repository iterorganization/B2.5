      subroutine get_from_movie_cdf(ncidcdf,istep,nCv,ns,state)
      use b2mod_types
      use b2mod_movies
      use b2us_plasma
#ifndef NO_CDF
      use b2mod_constants
#ifdef B25_EIRENE
      use eirmod_parmmod
     ,   , ncv_eir => ncv
      use eirmod_wneutrals
      use eirmod_extrab25
#else
      use b2mod_eirdiag
#endif
#endif
      use b2mod_subsys
      use b2mod_dimensions
      implicit none
      integer ncidcdf,istep,nCv,ns
      type (B2State) :: state
#ifndef NO_CDF
#ifndef B25_EIRENE
#     include <COUPLE/KOPPLDIM.F>
#endif
#     include <netcdf.inc>
      integer :: iret, start(3), count(3), varid
      real (kind=R8) :: time(1)
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
      start(2)=istep
      count(1)=nCv
      count(2)=1
      write(*,*) 'Getting te for slice ',istep
      iret=nf_inq_varid(ncidcdf,'te',varid)
      call check_cdf_status(iret)
      iret=nf_get_vara_double(ncidcdf,varid,start,count,state%pl%te)
      call check_cdf_status(iret)
      state%pl%te=state%pl%te*qe
      write(*,*) 'Getting ti for slice ',istep
      iret=nf_inq_varid(ncidcdf,'ti',varid)
      call check_cdf_status(iret)
      iret=nf_get_vara_double(ncidcdf,varid,start,count,state%pl%ti)
      call check_cdf_status(iret)
      state%pl%ti=state%pl%ti*qe
      write(*,*) 'Getting ne for slice ',istep
      iret=nf_inq_varid(ncidcdf,'ne',varid)
      call check_cdf_status(iret)
      iret=nf_get_vara_double(ncidcdf,varid,start,count,state%dv%ne)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=1
      start(3)=istep
      count(1)=nCv
      count(2)=ns
      count(3)=1
      write(*,*) 'Getting na for slice ',istep
      iret=nf_inq_varid(ncidcdf,'na',varid)
      call check_cdf_status(iret)
      iret=nf_get_vara_double(ncidcdf,varid,start,count,state%pl%na)
      call check_cdf_status(iret)
      if(nnatmi.gt.0) then
        count(2)=nnatmi
        write(*,*) 'Getting dab2 for slice ',istep
        iret=nf_inq_varid(ncidcdf,'dab2',varid)
        call check_cdf_status(iret)
        iret=nf_get_vara_double(ncidcdf,varid,start,count,dab2)
        call check_cdf_status(iret)
      endif
      if(nnmoli.gt.0) then
        count(2)=nnmoli
        write(*,*) 'Getting dmb2 for slice ',istep
        iret=nf_inq_varid(ncidcdf,'dmb2',varid)
        call check_cdf_status(iret)
        iret=nf_get_vara_double(ncidcdf,varid,start,count,dmb2)
        call check_cdf_status(iret)
      endif
#endif
      call subend ()
      return
      end subroutine get_from_movie_cdf
