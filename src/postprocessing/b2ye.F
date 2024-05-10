
*     B2YE is a post-processor for B2. It calculates statistical errors
*     based on averaged quantities collected during the averaging phase
*     of a B2 calculation.

      program b2ye

      use b2mod_types
      use b2mod_version
      use b2mod_constants
      use b2mod_b2cmpa
#ifndef NO_CDF
      use b2mod_mwti
     , , only : rwcdf
#endif
#ifdef B25_EIRENE
      use eirmod_infcop
#endif
      use b2mod_dimensions
      implicit none

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
     *      pomxip_av(:,:),  pomxap_av(:,:)
      real (kind=R8), allocatable, save ::
     *      nesepi_av(:,:),  tesepi_av(:,:),  tisepi_av(:,:),
     *      nesepm_av(:,:),  tesepm_av(:,:),  tisepm_av(:,:),
     *      nesepa_av(:,:),  tesepa_av(:,:),  tisepa_av(:,:),
     *      posepi_av(:,:),  posepm_av(:,:),  posepa_av(:,:)
      real (kind=R8), allocatable :: dum(:)
      real (kind=R8) :: tstepn(1), eps_s, tmax, meansample, t1, t2
      logical b2ye_tprf, ex
      character(5) :: rw
      integer :: io, iex, iun, idum(3), ns, varid, nbatch
      character lbli*120, check*10
      real (kind=R8) ::
     *   zpdum(0:nsdecl-1,4)
      intrinsic floor
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif

      naver = 0
      ntotdt = 0

#ifndef LEGACYCOMP
      open (newunit=iun,file='b2favere',access='SEQUENTIAL',
     .        form='FORMATTED',status='OLD',iostat=iex)
#else
      iun = newunit()
      open (unit=iun,file='b2favere',access='SEQUENTIAL',
     .        form='FORMATTED',status='OLD',iostat=iex)
#endif

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
     &    zpdum(0,1), zpdum(0,2), zpdum(0,3), zpdum(0,4), .false.)
        call cfruin (iun, 1, idum, 'naver')
        naver = idum(1)
        call cfruin (iun, 1, idum, 'ntotdt')
        ntotdt = idum(1)
      else
        write (*,*) ' file b2favere not found '
        write (*,*) ' results are questionable! '
      end if
      close (iun)

*  prepare reading from NETCDF file b2batch.nc or b2time.nc
      rw='read'
      call find_file('b2batch.nc',ex)
      if (.not. ex) then
        call find_file('b2time.nc',ex)
        if (.not. ex) then
          write(*,*) ' b2ye failed because no b2batch.nc or b2time.nc',
     .               ' files were available'
          goto 999
        endif
        iret = nf_open('b2time.nc',NF_NOWRITE,ncid)
        iret = nf_inq_varid(ncid,'batchsa',varid)
        if (varid.le.0) then
          write(*,*) ' b2ye failed because averaged quantities are',
     .               ' not present in b2batch.nc'
          goto 999
        endif
        write(*,*) ' reading batch averages from b2time.nc'
      else
        write(*,*) ' reading batch averages from b2batch.nc'
        iret = nf_open('b2batch.nc',NF_NOWRITE,ncid)
      endif

      imap(1)=1
      imap(2)=1
      call rwcdf (rw, ncid, 'nastep', imap, tstepn, iret)
      nastep = nint(tstepn(1))

      call rwcdf (rw, ncid, 'ntim_batch', imap, tstepn, iret)
      ntim_batch = nint(tstepn(1))
      nbatch = min(floor(real(naver/ntim_batch)), nastep)

      write (*,*) ' naver ', naver, ' ntotdt ', ntotdt
      write (*,*) ' nastep ',nastep, 'nncutmax', nncutmax,
     .            ' ntim_batch ', ntim_batch
      write (*,*) ' Calculating error on ',nbatch ,'batch averages'

      allocate (batchsa(nastep))
      allocate (nemxip_av(nastep,nncutmax))
      allocate (temxip_av(nastep,nncutmax))
      allocate (timxip_av(nastep,nncutmax))
      allocate (pomxip_av(nastep,nncutmax))
      allocate (nemxap_av(nastep,nncutmax))
      allocate (temxap_av(nastep,nncutmax))
      allocate (timxap_av(nastep,nncutmax))
      allocate (pomxap_av(nastep,nncutmax))

      allocate (nesepi_av(nastep,nncutmax))
      allocate (tesepi_av(nastep,nncutmax))
      allocate (tisepi_av(nastep,nncutmax))
      allocate (posepi_av(nastep,nncutmax))
      allocate (nesepm_av(nastep,nncutmax))
      allocate (tesepm_av(nastep,nncutmax))
      allocate (tisepm_av(nastep,nncutmax))
      allocate (posepm_av(nastep,nncutmax))
      allocate (nesepa_av(nastep,nncutmax))
      allocate (tesepa_av(nastep,nncutmax))
      allocate (tisepa_av(nastep,nncutmax))
      allocate (posepa_av(nastep,nncutmax))
      allocate (dum(nbatch))

*  read averaged quantities from b2time.nc

      call rwcdf(rw,ncid,'batchsa',imap,batchsa,iret)
      imap(1)=1
      imap(2)=nncutmax

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

      iret = nf_close(ncid)

*  write table of statistical errors

#ifndef LEGACYCOMP
      open (newunit=iout,file='b2ye.prt')
#else
      iout = newunit()
      open (unit=iout,file='b2ye.prt')
#endif
      write (iout,'(a16,5a12)') ' variable ','mean value','t',
     .             't1','t2','eps_s'
      if (b2ye_tprf(nemxip_av(:,1),nastep)) then
        dum = nemxip_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' nemxip ', meansample, tmax,
     .               t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant nemxip ',
     .               nemxip_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(temxip_av(:,1),nastep)) then
        dum = temxip_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' temxip ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant temxip ',
     .               temxip_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(timxip_av(:,1),nastep)) then
        dum = timxip_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' timxip ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant timxip ',
     .               timxip_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(pomxip_av(:,1),nastep)) then
        dum = pomxip_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' pomxip ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant pomxip ',
     .               pomxip_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(nemxap_av(:,1),nastep)) then
        dum = nemxap_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' nemxap ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant nemxap ',
     .               nemxap_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(temxap_av(:,1),nastep)) then
        dum = temxap_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' temxap ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant temxap ',
     .               temxap_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(timxap_av(:,1),nastep)) then
        dum = timxap_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' timxap ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant timxap ',
     .               timxap_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(pomxap_av(:,1),nastep)) then
        dum = pomxap_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' pomxap ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant pomxap ',
     .               pomxap_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(nesepi_av(:,1),nastep)) then
        dum = nesepi_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' nesepi ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant nesepi ',
     .               nesepi_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(tesepi_av(:,1),nastep)) then
        dum = tesepi_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' tesepi ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant tesepi ',
     .               tesepi_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(tisepi_av(:,1),nastep)) then
        dum = tisepi_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' tisepi ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant tisepi ',
     .               tisepi_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(posepi_av(:,1),nastep)) then
        dum = posepi_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' posepi ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant posepi ',
     .               posepi_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(nesepm_av(:,1),nastep)) then
        dum = nesepm_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' nesepm ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant nesepm ',
     .               nesepm_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(tesepm_av(:,1),nastep)) then
        dum = tesepm_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' tesepm ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant tesepm ',
     .               tesepm_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(tisepm_av(:,1),nastep)) then
        dum = tisepm_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' tisepm ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant tisepm ',
     .               tisepm_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(posepm_av(:,1),nastep)) then
        dum = posepm_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' posepm ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant posepm ',
     .               posepm_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(nesepa_av(:,1),nastep)) then
        dum = nesepa_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' nesepa ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant nesepa ',
     .               nesepa_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(tesepa_av(:,1),nastep)) then
        dum = tesepa_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' tesepa ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant tesepa ',
     .               tesepa_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(tisepa_av(:,1),nastep)) then
        dum = tisepa_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' tisepa ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant tisepa ',
     .               tisepa_av(1,1),0.,0.,0.,0.
      end if

      if (b2ye_tprf(posepa_av(:,1),nastep)) then
        dum = posepa_av(nastep-nbatch+1:nastep,1)
        call b2ye_sterr(dum, nbatch, eps_s, tmax, t1, t2,
     .                  meansample)
        write (iout,'(a16,5es12.4)') ' posepa ', meansample, tmax,
     .                 t1, t2, eps_s
      else
        write (iout,'(a16,5es12.4)') ' constant posepa ',
     .               posepa_av(1,1),0.,0.,0.,0.
      end if

      close (iout)

      deallocate (dum)
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

  999 stop 'b2ye'
      end program b2ye
