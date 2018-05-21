
*     B2YE is a post-processor for b2. It calculates statistical errors
*     based on averaged quantities collected during the averaging phase
*     of a b2 calculation.

      program b2ye

      use b2mod_types
      use b2mod_version
      use b2mod_constants
      use b2mod_b2cmpa
#ifndef NO_CDF
      use b2mod_mwti
#endif
#ifdef B25_EIRENE
      use eirmod_infcop
#endif

      implicit none

#include "DIMENSIONS.F"
#ifndef NO_CDF
#     include <netcdf.inc>

      integer :: nastep, ntim_batch, ntotdt, naver
      integer, parameter :: nncutmax = DEF_NCUT
      integer :: imap(maxvdims), iret, ncid, iout
      real (kind=R8), allocatable, save ::
     *      batchsa(:)
      real (kind=R8), allocatable, save ::
     *      nemxip_av(:,:),  temxip_av(:,:),  timxip_av(:,:),
     *      nemxap_av(:,:),  temxap_av(:,:),  timxap_av(:,:),
     *      pomxip_av(:,:),  pomxap_av(:,:),
     *      nemxip_std(:,:), temxip_std(:,:), timxip_std(:,:),
     *      nemxap_std(:,:), temxap_std(:,:), timxap_std(:,:),
     *      pomxip_std(:,:), pomxap_std(:,:)
      real (kind=R8), allocatable, save ::
     *      nesepi_av(:,:),  tesepi_av(:,:),  tisepi_av(:,:),
     *      nesepm_av(:,:),  tesepm_av(:,:),  tisepm_av(:,:),
     *      nesepa_av(:,:),  tesepa_av(:,:),  tisepa_av(:,:),
     *      posepi_av(:,:),  posepm_av(:,:),  posepa_av(:,:),
     *      nesepi_std(:,:), tesepi_std(:,:), tisepi_std(:,:),
     *      nesepm_std(:,:), tesepm_std(:,:), tisepm_std(:,:),
     *      nesepa_std(:,:), tesepa_std(:,:), tisepa_std(:,:),
     *      posepi_std(:,:), posepm_std(:,:), posepa_std(:,:)
      real (kind=R8) :: tstepn(1), eps_s, tmax, meansample, t1, t2
      logical b2ye_tprf
      character(5) :: rw
      integer :: io, iex, iun, idum(3), ns
      character lble*120, lbli*120, check*10
      character*10, save ::
     *   b2faveri_version = '00.000.000'
      real (kind=R8) ::
     *   zpdum(0:nsdecl-1,4)
       
      naver = 0
      ntotdt = 0

#ifdef B25_EIRENE
      iun = newunit()
#else
      iun= 99
#endif
      open (unit=iun,file='b2favere',access='SEQUENTIAL',
     .        form='FORMATTED',status='OLD',iostat=iex)

      io = iex
      read (iun,'(A10)',iostat=io) check
*  read number of total averaged iterations during running average from b2favere
      if (io == 0) then
        backspace iun
        call cfverr (iun, b2faveri_version)
        call cfruin (iun, 3, idum, 'nx,ny,ns')
        ns = idum(3)
   	call cfruch (iun, 120, lbli, 'label')
        call b2ruzd (iun, b2faveri_version, ns, 
     &          zpdum(0,1), zpdum(0,2), zpdum(0,3), zpdum(0,4), .false.)
        call cfruin (iun, 1, naver, 'naver')
        call cfruin (iun, 1, ntotdt, 'ntotdt')  
	write (*,*) ' naver ', naver, ' ntotdt ', ntotdt
      else
        write (*,*) ' file b2favere not found '
	write (*,*) ' results are questionable! '
      end if
      close (iun)
 
*  prepare reading from NETCDF file b2time.nc
      rw='read'
      iret = nf_open('b2time.nc',ncwrite,ncid)

      imap(1)=1
      call rwcdf (rw, ncid, 'nastep', imap, tstepn, iret)
      nastep = nint(tstepn(1))

      ntim_batch = naver / nastep

      write (*,*) ' nastep ',nastep, 'nncutmax', nncutmax,
     .            ' ntim_batch ', ntim_batch

      allocate (batchsa(nastep))
      allocate (nemxip_av(nncutmax,nastep))
      allocate (temxip_av(nncutmax,nastep))
      allocate (timxip_av(nncutmax,nastep))
      allocate (pomxip_av(nncutmax,nastep))
      allocate (nemxap_av(nncutmax,nastep))
      allocate (temxap_av(nncutmax,nastep))
      allocate (timxap_av(nncutmax,nastep))
      allocate (pomxap_av(nncutmax,nastep))

      allocate (nesepi_av(nncutmax,nastep))
      allocate (tesepi_av(nncutmax,nastep))
      allocate (tisepi_av(nncutmax,nastep))
      allocate (posepi_av(nncutmax,nastep))
      allocate (nesepm_av(nncutmax,nastep))
      allocate (tesepm_av(nncutmax,nastep))
      allocate (tisepm_av(nncutmax,nastep))
      allocate (posepm_av(nncutmax,nastep))
      allocate (nesepa_av(nncutmax,nastep))
      allocate (tesepa_av(nncutmax,nastep))
      allocate (tisepa_av(nncutmax,nastep))
      allocate (posepa_av(nncutmax,nastep))

*  read averaged quantities from b2time.nc

      imap(1)=1
      imap(2)=nncutmax

      call rwcdf(rw,ncid,'batchsa',imap,batchsa,iret)

      call rwcdf(rw,ncid,'nemxip_av',imap,nemxip_av,iret)
      call rwcdf(rw,ncid,'temxip_av',imap,temxip_av,iret)
      call rwcdf(rw,ncid,'timxip_av',imap,timxip_av,iret)
      call rwcdf(rw,ncid,'pomxip_av',imap,pomxip_av,iret)

      call rwcdf(rw,ncid,'nemxap_av',imap,nemxap_av,iret)
      call rwcdf(rw,ncid,'temxap_av',imap,temxap_av,iret)
      call rwcdf(rw,ncid,'timxap_av',imap,timxap_av,iret)
      call rwcdf(rw,ncid,'pomxap_av',imap,pomxap_av,iret)

      call rwcdf(rw,ncid,'nesepi_av',imap,nesepi_av,iret)
      call rwcdf(rw,ncid,'tesepi_av',imap,tesepi_av,iret)
      call rwcdf(rw,ncid,'tisepi_av',imap,tisepi_av,iret)
      call rwcdf(rw,ncid,'posepi_av',imap,posepi_av,iret)

      call rwcdf(rw,ncid,'nesepm_av',imap,nesepm_av,iret)
      call rwcdf(rw,ncid,'tesepm_av',imap,tesepm_av,iret)
      call rwcdf(rw,ncid,'tisepm_av',imap,tisepm_av,iret)
      call rwcdf(rw,ncid,'posepm_av',imap,posepm_av,iret)

      call rwcdf(rw,ncid,'nesepa_av',imap,nesepa_av,iret)
      call rwcdf(rw,ncid,'tesepa_av',imap,tesepa_av,iret)
      call rwcdf(rw,ncid,'tisepa_av',imap,tisepa_av,iret)
      call rwcdf(rw,ncid,'posepa_av',imap,posepa_av,iret)

*  write table of statistical errors

#ifdef B25_EIRENE
      iout = newunit()
#else
      iout = 99
#endif
      open (unit=iout,file='b2ye.prt')
      write (iout,'(a16,5a12)') ' variable ','mean value','t',
     .             't1','t2','eps_s'

      if (b2ye_tprf(nemxip_av(1,:),nastep)) then
        call b2ye_sterr(nemxip_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' nemxip ', meansample, tmax, 
     .               t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant nemxip ', 
     .               nemxip_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(temxip_av(1,:),nastep)) then
        call b2ye_sterr(temxip_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' temxip ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant temxip ', 
     .               temxip_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(timxip_av(1,:),nastep)) then
        call b2ye_sterr(timxip_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' timxip ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant timxip ', 
     .               timxip_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(pomxip_av(1,:),nastep)) then
        call b2ye_sterr(pomxip_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' pomxip ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant pomxip ', 
     .               pomxip_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(nemxap_av(1,:),nastep)) then
        call b2ye_sterr(nemxap_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' nemxap ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant nemxap ', 
     .               nemxap_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(temxap_av(1,:),nastep)) then
        call b2ye_sterr(temxap_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' temxap ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant temxap ', 
     .               temxap_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(timxap_av(1,:),nastep)) then
        call b2ye_sterr(timxap_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' timxap ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant timxap ', 
     .               timxap_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(pomxap_av(1,:),nastep)) then
        call b2ye_sterr(pomxap_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' pomxap ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant pomxap ', 
     .               pomxap_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(nesepi_av(1,:),nastep)) then
        call b2ye_sterr(nesepi_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' nesepi ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant nesepi ', 
     .               nesepi_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(tesepi_av(1,:),nastep)) then
        call b2ye_sterr(tesepi_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' tesepi ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant tesepi ', 
     .               tesepi_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(tisepi_av(1,:),nastep)) then
        call b2ye_sterr(tisepi_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' tisepi ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant tisepi ', 
     .               tisepi_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(posepi_av(1,:),nastep)) then
        call b2ye_sterr(posepi_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' posepi ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant posepi ', 
     .               posepi_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(nesepm_av(1,:),nastep)) then
        call b2ye_sterr(nesepm_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' nesepm ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant nesepm ', 
     .               nesepm_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(tesepm_av(1,:),nastep)) then
        call b2ye_sterr(tesepm_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' tesepm ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant tesepm ', 
     .               tesepm_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(tisepm_av(1,:),nastep)) then
        call b2ye_sterr(tisepm_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' tisepm ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant tisepm ', 
     .               tisepm_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(posepm_av(1,:),nastep)) then
        call b2ye_sterr(posepm_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' posepm ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant posepm ', 
     .               posepm_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(nesepa_av(1,:),nastep)) then
        call b2ye_sterr(nesepa_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' nesepa ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant nesepa ', 
     .               nesepa_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(tesepa_av(1,:),nastep)) then
        call b2ye_sterr(tesepa_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' tesepa ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant tesepa ', 
     .               tesepa_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(tisepa_av(1,:),nastep)) then
        call b2ye_sterr(tisepa_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' tisepa ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant tisepa ', 
     .               tisepa_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(posepa_av(1,:),nastep)) then
        call b2ye_sterr(posepa_av(1,:), nastep, ntotdt, ntim_batch, 
     .                 eps_s, tmax, t1, t2, meansample)
        write (iout,'(a16,5es12.4)') ' posepa ', meansample, tmax, 
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') 'constant posepa ', 
     .               posepa_av(1,1),0.,0.,0.,0.
      end if

      close (iout)

      iret = nf_close(ncid)

      deallocate (batchsa)
      deallocate (nemxip_av)
      deallocate (temxip_av)
      deallocate (timxip_av)
      deallocate (pomxip_av)
      deallocate (nemxap_av)
      deallocate (temxap_av)
      deallocate (timxap_av)
      deallocate (pomxap_av)

      deallocate (nesepi_av)
      deallocate (tesepi_av)
      deallocate (tisepi_av)
      deallocate (posepi_av)
      deallocate (nesepm_av)
      deallocate (tesepm_av)
      deallocate (tisepm_av)
      deallocate (posepm_av)
      deallocate (nesepa_av)
      deallocate (tesepa_av)
      deallocate (tisepa_av)
      deallocate (posepa_av)
#else
      write(*,*) 'b2ye failed because compiled with NO_CDF option'
#endif

      stop
      end program b2ye
