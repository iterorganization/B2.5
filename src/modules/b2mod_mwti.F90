module b2mod_mwti
  use b2mod_types , only : R8
  implicit none
  private
  public :: b2mwti, output_ds, dealloc_b2mod_mwti
  real (kind=R8), allocatable, save :: &
         nesepi_av(:), tesepi_av(:), tisepi_av(:), &
         nesepm_av(:), tesepm_av(:), tisepm_av(:), &
         nesepa_av(:), tesepa_av(:), tisepa_av(:), &
         posepi_av(:), posepm_av(:), posepa_av(:)
  real (kind=R8), allocatable, save :: &
         nemxip_av(:), temxip_av(:), timxip_av(:), &
         nemxap_av(:), temxap_av(:), timxap_av(:), &
         pomxip_av(:), pomxap_av(:)
  real (kind=R8), allocatable, save :: &
         nesepi_std(:), tesepi_std(:), tisepi_std(:), &
         nesepm_std(:), tesepm_std(:), tisepm_std(:), &
         nesepa_std(:), tesepa_std(:), tisepa_std(:), &
         posepi_std(:), posepm_std(:), posepa_std(:)
  real (kind=R8), allocatable, save :: &
         nemxip_std(:), temxip_std(:), timxip_std(:), &
         nemxap_std(:), temxap_std(:), timxap_std(:), &
         pomxip_std(:), pomxap_std(:)
#ifndef NO_CDF
  public :: rwcdf, rwcdf_settime, rwcdf_setbatch, b2crtimecdf
#endif
contains

  subroutine b2mwti (itim, tim, ntim, b2time, ntim_batch, &
                     nx, ny, ns, ismain, ismain0, BoRiS, &
                     lwti, lwav, luav)
    use b2mod_geo
    use b2mod_plasma
    use b2mod_rates
    use b2mod_residuals
    use b2mod_sources
    use b2mod_transport
    use b2mod_anomalous_transport
    use b2mod_neutrals_namelist
    use b2mod_work
    use b2mod_indirect
    use b2mod_constants
    use b2mod_tallies
    use b2mod_wall
    use b2mod_b2cmpa
    use b2mod_external
#ifndef SOLPS4_3
#ifdef B25_EIRENE
    use eirmod_extrab25
#endif
#endif
    implicit none
    !   ..input arguments (unchanged on exit)
    integer, Intent(In) :: itim, ntim, b2time, ntim_batch, &
                           nx, ny, ns, ismain, ismain0
    real (kind=R8), Intent(In) :: tim, BoRiS
    logical, Intent(In) :: lwti, lwav, luav
    !   ..output arguments (unspecified on entry)
    !     (none)
    !   ..common blocks
#ifndef NO_CDF
#     include <netcdf.inc>
#endif
    !-----------------------------------------------------------------------
    !.documentation
    !
    !  1. purpose
    !
    !     B2MWTI writes to an un*formatted file a selection of data
    !     describing the state of the plasma.
    !
    !
    !  3. description (see also routine b2cdca)
    !
    !     This routine writes to unit (nput) data that are suitable for
    !     producing movie output. Unit nput must be connected to an
    !     un*formatted file. Presently, only the changes in the plasma
    !     state are written out.
    !
    !
    !-----------------------------------------------------------------------
    !.declarations

    !   ..local variables
    integer ncall, ntstep, nastep
    real (kind=R8) :: &
         fnixip(nncutmax), feexip(nncutmax), feixip(nncutmax), &
         fnixap(nncutmax), feexap(nncutmax), feixap(nncutmax), &
         pwmxip(nncutmax), fchxip(nncutmax), fchxap(nncutmax), &
         fetxip(nncutmax), fetxap(nncutmax)
    real (kind=R8) :: &
         fniyip(nncutmax), feeyip(nncutmax), feiyip(nncutmax), &
         fniyap(nncutmax), feeyap(nncutmax), feiyap(nncutmax), &
         pwmxap(nncutmax), fchyip(nncutmax), fchyap(nncutmax), &
         fetyip(nncutmax), fetyap(nncutmax)
    real (kind=R8) :: &
         nemxip(nncutmax), temxip(nncutmax), timxip(nncutmax), &
         nemxap(nncutmax), temxap(nncutmax), timxap(nncutmax), &
         pomxip(nncutmax), pomxap(nncutmax), tpmxip(nncutmax), &
         tpmxap(nncutmax)
    real (kind=R8) :: &
         fnisip(nncutmax), feesip(nncutmax), feisip(nncutmax), &
         fnisap(nncutmax), feesap(nncutmax), feisap(nncutmax), &
         fchsip(nncutmax), fchsap(nncutmax), &
         fetsip(nncutmax), fetsap(nncutmax)
    real (kind=R8) :: &
         fnisipp(nncutmax), feesipp(nncutmax), feisipp(nncutmax), &
         fnisapp(nncutmax), feesapp(nncutmax), feisapp(nncutmax), &
         fchsipp(nncutmax), fchsapp(nncutmax), &
         fetsipp(nncutmax), fetsapp(nncutmax)
    real (kind=R8) :: &
         tmne(1),tmte(1),tmti(1),tmvol

    integer iy, ix, ic, ixtl, ixtr, jsep
    integer jxi, jxa, target_offset, ix_off
    integer iyastrt, iyistrt, iylstrt, iyrstrt, iytlstrt, iytrstrt, &
         iyaend, iyiend, iylend, iyrend, iytlend, iytrend, &
         nya, nyi, nybl, nybr, nytl, nytr, nc

    !   ..procedures
    external subini, subend, xertst, ipgeti, batch_average
    real(kind=R8) :: fnitmp, feetmp, feitmp, fchtmp, fettmp, pwrtmp
    integer, save :: write_2d = 0
#ifndef NO_CDF
    integer, save :: ncid, nbatch
    integer imap(maxvdims), dims(1), iret
    integer nvars, natts, ndims, unlimid, nastepid
    real (kind=R8) :: fac
    real (kind=R8) :: &
         nesepi(nncutmax), tesepi(nncutmax), tisepi(nncutmax), &
         nesepm(nncutmax), tesepm(nncutmax), tisepm(nncutmax), &
         nesepa(nncutmax), tesepa(nncutmax), tisepa(nncutmax), &
         posepi(nncutmax), posepm(nncutmax), posepa(nncutmax), &
         dnsepm(nncutmax), dpsepm(nncutmax), kesepm(nncutmax), &
         kisepm(nncutmax), vxsepm(nncutmax), vysepm(nncutmax), &
         vssepm(nncutmax), tpsepi(nncutmax), tpsepa(nncutmax)
    real (kind=R8) :: &
         tmhacore(1), tmhasol(1), tmhadiv(1), slice(-1:ny), tstepn(1)
    real (kind=R8) :: &
         timesa(1), batchsa(1)
    real (kind=R8), save :: stim = 0.0_R8
    logical ex
    character*5 rw
    character*256, save :: filename, filename_av
    real(kind=R8) :: rratio
    external rratio
#endif
    !   ..initialisation
    save ncall, ntstep, jxi, jxa, jsep, ixtl, ixtr, target_offset, &
         iyastrt, iyistrt, iylstrt, iyrstrt, iytlstrt, iytrstrt, &
         iyaend,  iyiend,  iylend,  iyrend,  iytlend,  iytrend, &
         nc, nya, nyi, nybl, nybr, nytl, nytr, nastep
    data ncall/0/, target_offset/1/

    !-----------------------------------------------------------------------
    !.computation

    ! ..preliminaries
    !   ..subprogram start-up calls
    call subini ('b2mwti')
    !     ..test nx, ny
    call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
    call xertst (1.le.ns, 'faulty argument ns')
    call xertst (0.le.ismain.and.ismain.lt.ns, &
         'invalid main plasma species index ismain')
    !   ..extensive tests on first few calls
    if (ncall.eq.0) then
      !   ..test state
      call get_jsep(nx,ny,jxi,jxa,jsep)
      call ipgeti ('b2mwti_2dwrite',write_2d)
      call xertst (0.le.write_2d.and.write_2d.le.2,'faulty internal parameter write_2d')
      call ipgeti ('b2mwti_target_offset',target_offset)
      call xertst (0.le.target_offset.and.target_offset.le.1,'faulty internal parameter target_offset')
      write(*,*) 'target_offset ', target_offset
      call output_ds(crx,cry,nx,ny, -1,+target_offset,jsep,iylstrt,iylend,'dsl')
      call output_ds(crx,cry,nx,ny,jxi,0,jsep,iyistrt,iyiend,'dsi')
      call output_ds(crx,cry,nx,ny,jxa,0,jsep,iyastrt,iyaend,'dsa')
      call output_ds(crx,cry,nx,ny, nx,-target_offset,jsep,iyrstrt,iyrend,'dsr')
      if (nnreg(0).eq.8) then
        ixtl = 0
        do while (rightix(ixtl,max(topcut(1),topcut(2))).ne.nx+1.and.ixtl.lt.nx)
          ixtl=ixtl+1
        enddo
        ixtr = ixtl
        do while (leftix(ixtr,max(topcut(1),topcut(2))).ne.-2 .and. ixtr.lt.nx)
          ixtr=ixtr+1
        enddo
        call output_ds(crx,cry,nx,ny,ixtl,-target_offset,jsep,iytlstrt,iytlend,'dstl')
        call output_ds(crx,cry,nx,ny,ixtr,+target_offset,jsep,iytrstrt,iytrend,'dstr')
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
      ! Target areas
      open(99,file='dsL')
      do iy=-1,ny
        if(region(-1,iy,0).ne.0) write(99,*) gs(rightix(-1,iy),rightiy(-1,iy),0)
      enddo
      close(99)
      open(99,file='dsR')
      do iy=-1,ny
        if(region(nx,iy,0).ne.0) write(99,*) gs(nx,iy,0)
      enddo
      close(99)
      if (nnreg(0).eq.8) then
        open(99,file='dsTL')
        do iy=iytlstrt,iytlend
          write(99,*) gs(ixtl,iy,0)
        enddo
        close(99)
        open(99,file='dsTR')
        do iy=iytrstrt,iytrend
          write(99,*) gs(rightix(ixtr,iy),rightiy(ixtr,iy),0)
        enddo
        close(99)
      endif
      ! Poloidal contact areas
      open(99,file='dsLP')
      do iy=-1,ny
        if(region(-1,iy,0).ne.0) write(99,*) gs(rightix(-1,iy),rightiy(-1,iy),0)*qc(rightix(-1,iy),rightiy(-1,iy))
      enddo
      close(99)
      open(99,file='dsRP')
      do iy=-1,ny
        if(region(nx,iy,0).ne.0) write(99,*) gs(nx,iy,0)*qc(nx,iy)
      enddo
      close(99)
      if (nnreg(0).eq.8) then
        open(99,file='dsTLP')
        do iy=iytlstrt,iytlend
          write(99,*) gs(ixtl,iy,0)*qc(ixtl,iy)
        enddo
        close(99)
        open(99,file='dsTRP')
        do iy=iytrstrt,iytrend
          write(99,*) gs(rightix(ixtr,iy),rightiy(ixtr,iy),0)*qc(rightix(ixtr,iy),rightiy(ixtr,iy))
        enddo
        close(99)
      endif
#ifndef NO_CDF
      if (b2time.gt.0) then
        filename='b2time.nc'
        call find_file(filename,ex)
        call ipgetr ('b2mndr_stim', stim)
        if (.not.ex.or.stim.ge.0.0_R8) then
          ntstep = 0
          write(6,'(a)') trim(filename)//' will be created'
          call b2crtimecdf(filename, &
            nx, ny, nybl, nytl, nytr, nybr, nya, nyi, nc, ns, write_2d, &
            ncid, .false., iret)
          call check_cdf_status(iret)
          iret = nf_open(trim(filename),or(NF_WRITE,NF_SHARE),ncid)
          call check_cdf_status(iret)
        else if (ex.and.stim.lt.0.0_R8) then
          rw='read'
          iret = nf_open(filename,NF_NOWRITE,ncid)
          call check_cdf_status(iret)
          iret = nf_inq(ncid,ndims,nvars,natts,unlimid)
          call check_cdf_status(iret)
          imap(1)=1
          call rwcdf (rw, ncid, 'ntstep', imap, tstepn, iret)
          call check_cdf_status(iret)
          ntstep = nint(tstepn(1))
          iret = nf_close(ncid)
          write(6,'(a)') trim(filename)//' will be appended'
          iret = nf_open(trim(filename),or(NF_WRITE,NF_SHARE),ncid)
          call check_cdf_status(iret)
        else
          ntstep = 0
          write(6,'(a)') trim(filename)//' will be replaced'
          call b2crtimecdf(filename, &
            nx, ny, nybl, nytl, nytr, nybr, nya, nyi, nc, ns, &
            write_2d, ncid, .false., iret)
          call check_cdf_status(iret)
          iret = nf_open(trim(filename),or(NF_WRITE,NF_SHARE),ncid)
          call check_cdf_status(iret)
        end if
        write(*,*) 'ntstep = ', ntstep
        rw='write'
        imap(1)=1
        tstepn(1) = ntstep
        call rwcdf (rw, ncid, 'ntstep', imap, tstepn, iret)
        call check_cdf_status(iret)
        iret = nf_close(ncid)
        call check_cdf_status(iret)
      end if

      if (ntim_batch.gt.0) then
        nbatch = ntim/ntim_batch
        if (nbatch.gt.0) then
          write(*,*) 'nbatch = ', nbatch
          filename_av='b2batch.nc'
          call find_file(filename_av,ex)
          if (.not.ex.or.stim.ge.0.0_R8) then
            nastep = 0
            write(6,'(a)') trim(filename_av)//' will be created'
            call b2crtimecdf(filename_av, &
              nx, ny, nybl, nytl, nytr, nybr, nya, nyi, nc, ns, write_2d, &
              ncid, .true., iret)
            call check_cdf_status(iret)
            iret = nf_open(trim(filename_av),or(NF_WRITE,NF_SHARE),ncid)
            call check_cdf_status(iret)
          else if (ex.and.stim.lt.0.0_R8) then
            rw='read'
            iret = nf_open(filename_av,NF_NOWRITE,ncid)
            call check_cdf_status(iret)
            iret = nf_inq(ncid,ndims,nvars,natts,unlimid)
            call check_cdf_status(iret)
            imap(1)=1
            call rwcdf (rw, ncid, 'nastep', imap, tstepn, iret)
            call check_cdf_status(iret)
            nastep = nint(tstepn(1))
            call rwcdf ('read', ncid, 'ntim_batch', imap, tstepn, iret)
            call check_cdf_status(iret)
            if (tstepn(1).eq.ntim_batch) then
              call rwcdf (rw, ncid, 'nastep', imap, tstepn, iret)
              call check_cdf_status(iret)
              nastep = nint(tstepn(1))
              iret = nf_close(ncid)
              call check_cdf_status(iret)
              write(6,'(a)') trim(filename_av)//' will be appended'
            else
              write(*,*)'WARNING: ntim_batch has been changed'
              write(*,*)'WARNING: statistical error assessment will not be reliable'
              write(*,*)'WARNING: restarting the batch averaging'
              write(6,'(a)') trim(filename_av)//' will be replaced'
              nastep = 0
              iret = nf_close(ncid)
              call check_cdf_status(iret)
              call b2crtimecdf(filename_av, &
               nx, ny, nybl, nytl, nytr, nybr, nya, nyi, nc, ns, &
               write_2d, ncid, .true., iret)
            endif
            iret = nf_open(trim(filename_av),or(NF_WRITE,NF_SHARE),ncid)
            call check_cdf_status(iret)
          else
            nastep = 0
            write(6,'(a)') trim(filename_av)//' will be replaced'
            call b2crtimecdf(filename_av, &
              nx, ny, nybl, nytl, nytr, nybr, nya, nyi, nc, ns, &
              write_2d, ncid, .true., iret)
            call check_cdf_status(iret)
            iret = nf_open(trim(filename_av),or(NF_WRITE,NF_SHARE),ncid)
            call check_cdf_status(iret)
          end if
          write(*,*) 'nastep = ', nastep
          rw='write'
          imap(1)=1
          tstepn(1) = nastep
          call rwcdf (rw, ncid, 'nastep', imap, tstepn, iret)
          call check_cdf_status(iret)
          tstepn(1) = ntim_batch
          call rwcdf (rw, ncid, 'ntim_batch', imap, tstepn, iret)
          call check_cdf_status(iret)
          iret = nf_close(ncid)
          call check_cdf_status(iret)
        end if
      end if
#endif
      allocate (nesepm_av(1:nncutmax))
      allocate (tesepm_av(1:nncutmax))
      allocate (tisepm_av(1:nncutmax))
      allocate (posepm_av(1:nncutmax))
      allocate (nesepi_av(1:nncutmax))
      allocate (tesepi_av(1:nncutmax))
      allocate (tisepi_av(1:nncutmax))
      allocate (posepi_av(1:nncutmax))
      allocate (nesepa_av(1:nncutmax))
      allocate (tesepa_av(1:nncutmax))
      allocate (tisepa_av(1:nncutmax))
      allocate (posepa_av(1:nncutmax))
      allocate (nemxip_av(1:nncutmax))
      allocate (temxip_av(1:nncutmax))
      allocate (timxip_av(1:nncutmax))
      allocate (pomxip_av(1:nncutmax))
      allocate (nemxap_av(1:nncutmax))
      allocate (temxap_av(1:nncutmax))
      allocate (timxap_av(1:nncutmax))
      allocate (pomxap_av(1:nncutmax))
      allocate (nesepm_std(1:nncutmax))
      allocate (tesepm_std(1:nncutmax))
      allocate (tisepm_std(1:nncutmax))
      allocate (posepm_std(1:nncutmax))
      allocate (nesepi_std(1:nncutmax))
      allocate (tesepi_std(1:nncutmax))
      allocate (tisepi_std(1:nncutmax))
      allocate (posepi_std(1:nncutmax))
      allocate (nesepa_std(1:nncutmax))
      allocate (tesepa_std(1:nncutmax))
      allocate (tisepa_std(1:nncutmax))
      allocate (posepa_std(1:nncutmax))
      allocate (nemxip_std(1:nncutmax))
      allocate (temxip_std(1:nncutmax))
      allocate (timxip_std(1:nncutmax))
      allocate (pomxip_std(1:nncutmax))
      allocate (nemxap_std(1:nncutmax))
      allocate (temxap_std(1:nncutmax))
      allocate (timxap_std(1:nncutmax))
      allocate (pomxap_std(1:nncutmax))
      nesepm_av = 0.0_R8
      tesepm_av = 0.0_R8
      tisepm_av = 0.0_R8
      posepm_av = 0.0_R8
      nesepi_av = 0.0_R8
      tesepi_av = 0.0_R8
      tisepi_av = 0.0_R8
      posepi_av = 0.0_R8
      nesepa_av = 0.0_R8
      tesepa_av = 0.0_R8
      tisepa_av = 0.0_R8
      posepa_av = 0.0_R8
      nemxip_av = 0.0_R8
      temxip_av = 0.0_R8
      timxip_av = 0.0_R8
      pomxip_av = 0.0_R8
      nemxap_av = 0.0_R8
      temxap_av = 0.0_R8
      timxap_av = 0.0_R8
      pomxap_av = 0.0_R8
      nesepm_std = 0.0_R8
      tesepm_std = 0.0_R8
      tisepm_std = 0.0_R8
      posepm_std = 0.0_R8
      nesepi_std = 0.0_R8
      tesepi_std = 0.0_R8
      tisepi_std = 0.0_R8
      posepi_std = 0.0_R8
      nesepa_std = 0.0_R8
      tesepa_std = 0.0_R8
      tisepa_std = 0.0_R8
      posepa_std = 0.0_R8
      nemxip_std = 0.0_R8
      temxip_std = 0.0_R8
      timxip_std = 0.0_R8
      pomxip_std = 0.0_R8
      nemxap_std = 0.0_R8
      temxap_std = 0.0_R8
      timxap_std = 0.0_R8
      pomxap_std = 0.0_R8
    endif! ncall == 0
    if(ncall.lt.3) then
      call xertst (0.le.itim, 'faulty parameter itim')
      call xertst (0.0_R8.le.tim, 'faulty parameter tim')
    endif
    !
    !   ..compute change in plasma state
    !
    if (lwti) then
      ntstep = ntstep + 1
    !     write(*,*) 'ntstep = ',ntstep
#ifndef NO_CDF
      call rwcdf_settime ('time', ntstep)
      timesa(1) = tim
#endif
    endif

    if (lwav) then
      nastep = nastep + 1
    !   write(*,*) 'nastep = ',nastep
#ifndef NO_CDF
      call rwcdf_setbatch ('batch', nastep)
      batchsa(1) = tim
#endif
    endif
    !
    !    total flows to the divertor plates
    !

    fnixip = 0.0_R8; feexip = 0.0_R8; feixip = 0.0_R8; fchxip = 0.0_R8; fetxip = 0.0_R8
    nemxip = 0.0_R8; temxip = 0.0_R8; timxip = 0.0_R8; pomxip = 0.0_R8; pwmxip = 0.0_R8; tpmxip = 0.0_R8
    ix = -1 ! 1
    ix_off  = ix + target_offset
    do iy = iylstrt,iylend
      call calc_fet(ix,iy,'L',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp,pwrtmp)
      fnixip(1) = fnixip(1) + fnitmp
      feexip(1) = feexip(1) + feetmp
      feixip(1) = feixip(1) + feitmp
      fchxip(1) = fchxip(1) + fchtmp
      fetxip(1) = fetxip(1) + fettmp
      nemxip(1) = max(nemxip(1), ne(ix_off,iy))
      temxip(1) = max(temxip(1), te(ix_off,iy))
      timxip(1) = max(timxip(1), ti(ix_off,iy))
      pomxip(1) = max(pomxip(1), po(ix_off,iy))
      pwmxip(1) = max(pwmxip(1), pwrtmp)
      if (bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and. xymap(ix,iy).ne.0) then
        tpmxip(1) = max(tpmxip(1), target_temp(xymap(ix,iy),1))
      endif
    enddo

    fnixap = 0.0_R8; feexap = 0.0_R8; feixap = 0.0_R8; fchxap = 0.0_R8; fetxap = 0.0_R8
    nemxap = 0.0_R8; temxap = 0.0_R8; timxap = 0.0_R8; pomxap = 0.0_R8; pwmxap = 0.0_R8; tpmxap = 0.0_R8
    ix = nx ! 2
    ix_off  = ix - target_offset
    do iy = iyrstrt,iyrend
      call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp,pwrtmp)
      fnixap(1) = fnixap(1) + fnitmp
      feexap(1) = feexap(1) + feetmp
      feixap(1) = feixap(1) + feitmp
      fchxap(1) = fchxap(1) + fchtmp
      fetxap(1) = fetxap(1) + fettmp
      nemxap(1) = max(nemxap(1), ne(ix_off,iy))
      temxap(1) = max(temxap(1), te(ix_off,iy))
      timxap(1) = max(timxap(1), ti(ix_off,iy))
      pomxap(1) = max(pomxap(1), po(ix_off,iy))
      pwmxap(1) = max(pwmxap(1), pwrtmp)
      if (bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and. xymap(ix,iy).ne.0) then
        tpmxap(1) = max(tpmxap(1), target_temp(xymap(ix,iy),1))
      endif
    enddo

    if(nncut.ge.2) then
      ix = ixtr ! 3
      ix_off  = ix + target_offset
      do iy = iytrstrt,iytrend
        call calc_fet(ix,iy,'L',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp,pwrtmp)
        fnixap(2) = fnixap(2) + fnitmp
        feexap(2) = feexap(2) + feetmp
        feixap(2) = feixap(2) + feitmp
        fchxap(2) = fchxap(2) + fchtmp
        fetxap(2) = fetxap(2) + fettmp
        nemxap(2) = max(nemxap(2), ne(ix_off,iy))
        temxap(2) = max(temxap(2), te(ix_off,iy))
        timxap(2) = max(timxap(2), ti(ix_off,iy))
        pomxap(2) = max(pomxap(2), po(ix_off,iy))
        pwmxap(2) = max(pwmxap(2), pwrtmp)
        if (bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and. xymap(ix,iy).ne.0) then
          tpmxap(2) = max(tpmxap(2), target_temp(xymap(ix,iy),1))
        endif
      enddo

      ix = ixtl ! 4
      ix_off  = ix - target_offset
      do iy = iytlstrt,iytlend
        call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp,pwrtmp)
        fnixip(2) = fnixip(2) + fnitmp
        feexip(2) = feexip(2) + feetmp
        feixip(2) = feixip(2) + feitmp
        fchxip(2) = fchxip(2) + fchtmp
        fetxip(2) = fetxip(2) + fettmp
        nemxip(2) = max(nemxip(2), ne(ix_off,iy))
        temxip(2) = max(temxip(2), te(ix_off,iy))
        timxip(2) = max(timxip(2), ti(ix_off,iy))
        pomxip(2) = max(pomxip(2), po(ix_off,iy))
        pwmxip(2) = max(pwmxip(2), pwrtmp)
        if (bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and. xymap(ix,iy).ne.0) then
          tpmxip(2) = max(tpmxip(2), target_temp(xymap(ix,iy),1))
        endif
      enddo
    endif

    fnisip = 0.0_R8; feesip = 0.0_R8; feisip = 0.0_R8; fchsip = 0.0_R8; fetsip = 0.0_R8
    fnisap = 0.0_R8; feesap = 0.0_R8; feisap = 0.0_R8; fchsap = 0.0_R8; fetsap = 0.0_R8
    fnisipp = 0.0_R8; feesipp = 0.0_R8; feisipp = 0.0_R8; fetsipp = 0.0_R8; fchsipp = 0.0_R8
    fnisapp = 0.0_R8; feesapp = 0.0_R8; feisapp = 0.0_R8; fetsapp = 0.0_R8; fchsapp = 0.0_R8
    if(nnreg(0).ge.3) then
      do ic = 1, nncut
        do iy = -1,jsep
          ix = leftcut(ic) ! 5
          call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fnisipp(ic) = fnisipp(ic) + fnitmp
          feesipp(ic) = feesipp(ic) + feetmp
          feisipp(ic) = feisipp(ic) + feitmp
          fchsipp(ic) = fchsipp(ic) + fchtmp
          fetsipp(ic) = fetsipp(ic) + fettmp
          ix = rightcut(ic) ! 6
          call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fnisapp(ic) = fnisapp(ic) + fnitmp
          feesapp(ic) = feesapp(ic) + feetmp
          feisapp(ic) = feisapp(ic) + feitmp
          fchsapp(ic) = fchsapp(ic) + fchtmp
          fetsapp(ic) = fetsapp(ic) + fettmp
        enddo
        do iy = jsep+1,ny
          ix = leftcut(ic) ! 7
          call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fnisip(ic) = fnisip(ic) + fnitmp
          feesip(ic) = feesip(ic) + feetmp
          feisip(ic) = feisip(ic) + feitmp
          fchsip(ic) = fchsip(ic) + fchtmp
          fetsip(ic) = fetsip(ic) + fettmp
          ix = rightcut(ic) ! 8
          call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fnisap(ic) = fnisap(ic) + fnitmp
          feesap(ic) = feesap(ic) + feetmp
          feisap(ic) = feisap(ic) + feitmp
          fchsap(ic) = fchsap(ic) + fchtmp
          fetsap(ic) = fetsap(ic) + fettmp
        enddo
      enddo
    endif

    fniyip = 0.0_R8; feeyip = 0.0_R8; feiyip = 0.0_R8; fetyip = 0.0_R8; fchyip = 0.0_R8
    fniyap = 0.0_R8; feeyap = 0.0_R8; feiyap = 0.0_R8; fetyap = 0.0_R8; fchyap = 0.0_R8

    if(nnreg(0).eq.4) then
      do ix = -1,nx
        if(region(ix,ny,0).eq.2) then
          iy = ny ! 9
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyip(1) = fniyip(1) + fnitmp
          feeyip(1) = feeyip(1) + feetmp
          feiyip(1) = feiyip(1) + feitmp
          fchyip(1) = fchyip(1) + fchtmp
          fetyip(1) = fetyip(1) + fettmp
        endif
        if(region(ix,ny,0).eq.3.or.region(ix,ny,0).eq.4) then
          iy = ny ! 10
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(1) = fniyap(1) + fnitmp
          feeyap(1) = feeyap(1) + feetmp
          feiyap(1) = feiyap(1) + feitmp
          fchyap(1) = fchyap(1) + fchtmp
          fetyap(1) = fetyap(1) + fettmp
        endif
        if(region(ix,-1,0).eq.3.or.region(ix,-1,0).eq.4) then
          iy = -1 ! 11
          call calc_fet(ix,iy,'B',-1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(1) = fniyap(1) + fnitmp
          feeyap(1) = feeyap(1) + feetmp
          feiyap(1) = feiyap(1) + feitmp
          fchyap(1) = fchyap(1) + fchtmp
          fetyap(1) = fetyap(1) + fettmp
        endif
      enddo
    elseif(nnreg(0).eq.5) then
      do ix = -1,nx
        if(region(ix,ny,0).eq.5) then
          iy = ny ! 12
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyip(1) = fniyip(1) + fnitmp
          feeyip(1) = feeyip(1) + feetmp
          feiyip(1) = feiyip(1) + feitmp
          fchyip(1) = fchyip(1) + fchtmp
          fetyip(1) = fetyip(1) + fettmp
        endif
        if(region(ix,-1,0).eq.3.or.region(ix,-1,0).eq.4) then
          iy = -1 ! 13
          call calc_fet(ix,iy,'B',-1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(1) = fniyap(1) + fnitmp
          feeyap(1) = feeyap(1) + feetmp
          feiyap(1) = feiyap(1) + feitmp
          fchyap(1) = fchyap(1) + fchtmp
          fetyap(1) = fetyap(1) + fettmp
        endif
      enddo
    elseif (nnreg(0).eq.8) then
      do ix = -1,nx
        if(mod(region(ix,ny,0),4).eq.2) then
          iy = ny ! 14
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyip(region(ix,ny,0)/4+1) = fniyip(region(ix,ny,0)/4+1) + fnitmp
          feeyip(region(ix,ny,0)/4+1) = feeyip(region(ix,ny,0)/4+1) + feetmp
          feiyip(region(ix,ny,0)/4+1) = feiyip(region(ix,ny,0)/4+1) + feitmp
          fchyip(region(ix,ny,0)/4+1) = fchyip(region(ix,ny,0)/4+1) + fchtmp
          fetyip(region(ix,ny,0)/4+1) = fetyip(region(ix,ny,0)/4+1) + fettmp
        endif
        if(region(ix,ny,0).eq.3 .or. region(ix,ny,0).eq.8) then
          iy = ny ! 15
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(1) = fniyap(1) + fnitmp
          feeyap(1) = feeyap(1) + feetmp
          feiyap(1) = feiyap(1) + feitmp
          fchyap(1) = fchyap(1) + fchtmp
          fetyap(1) = fetyap(1) + fettmp
        endif
        if(region(ix,-1,0).eq.3 .or. region(ix,-1,0).eq.8) then
          iy = -1 ! 16
          call calc_fet(ix,iy,'B',-1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(1) = fniyap(1) + fnitmp
          feeyap(1) = feeyap(1) + feetmp
          feiyap(1) = feiyap(1) + feitmp
          fchyap(1) = fchyap(1) + fchtmp
          fetyap(1) = fetyap(1) + fettmp
        endif
        if(region(ix,ny,0).eq.4 .or. region(ix,ny,0).eq.7) then
          iy = ny ! 17
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(2) = fniyap(2) + fnitmp
          feeyap(2) = feeyap(2) + feetmp
          feiyap(2) = feiyap(2) + feitmp
          fchyap(2) = fchyap(2) + fchtmp
          fetyap(2) = fetyap(2) + fettmp
        endif
        if(region(ix,-1,0).eq.4 .or. region(ix,-1,0).eq.7) then
          iy = -1 ! 18
          call calc_fet(ix,iy,'B',-1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(2) = fniyap(2) + fnitmp
          feeyap(2) = feeyap(2) + feetmp
          feiyap(2) = feiyap(2) + feitmp
          fchyap(2) = fchyap(2) + fchtmp
          fetyap(2) = fetyap(2) + fettmp
        endif
      enddo
    else if (nnreg(0).eq.2) then
      do ix = -1,nx
        if(region(ix,ny,0).eq.2) then
          iy = ny ! 19
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyip(1) = fniyip(1) + fnitmp
          feeyip(1) = feeyip(1) + feetmp
          feiyip(1) = feiyip(1) + feitmp
          fchyip(1) = fchyip(1) + fchtmp
          fetyip(1) = fetyip(1) + fettmp
        endif
      enddo
    else
      do ix = -1,nx
        iy = ny ! 20
        call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
        fniyip(1) = fniyip(1) + fnitmp
        feeyip(1) = feeyip(1) + feetmp
        feiyip(1) = feiyip(1) + feitmp
        fchyip(1) = fchyip(1) + fchtmp
        fetyip(1) = fetyip(1) + fettmp
      enddo
    endif ! nnreg check

    !
    !    other quantities related to the target plates
    !
#ifndef NO_CDF
    if(nnreg(0).ne.2) then
      nesepi(1) = 0.5_R8 * (ne(-1+target_offset,jsep)+ne(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep)))
      tesepi(1) = 0.5_R8/ev * (te(-1+target_offset,jsep) + te(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep)))
      tisepi(1) = 0.5_R8/ev * (ti(-1+target_offset,jsep) + ti(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep)))
      if(xymap(-1,jsep).gt.0 .and. xymap(topix(-1,jsep),topiy(-1,jsep)).gt.0) then
        tpsepi(1) = 0.5_R8 * (target_temp(xymap(-1,jsep),1) + target_temp(xymap(topix(-1,jsep),topiy(-1,jsep)),1))
      else
        tpsepi(1) = 0.0_R8
      endif
      posepi(1) = 0.5_R8 * (po(-1+target_offset,jsep) + po(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep)))
    else
      nesepi(1) = ne(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep))
      tesepi(1) = 1.0_R8/ev * te(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep))
      tisepi(1) = 1.0_R8/ev * ti(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep))
      if(xymap(topix(-1,jsep),topiy(-1,jsep)).gt.0) then
        tpsepi(1) = target_temp(xymap(topix(-1,jsep),topiy(-1,jsep)),1)
      else
        tpsepi(1) = 0.0_R8
      endif
      posepi(1) = po(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep))
    endif
    nesepm(1) = 0.5_R8 * (ne(jxa,jsep)+ ne(topix(jxa,jsep),topiy(jxa,jsep)))
    tesepm(1) = 0.5_R8 * (te(jxa,jsep)+ te(topix(jxa,jsep),topiy(jxa,jsep)))/ev
    tisepm(1) = 0.5_R8 * (ti(jxa,jsep)+ ti(topix(jxa,jsep),topiy(jxa,jsep)))/ev
    posepm(1) = 0.5_R8 * (po(jxa,jsep)+ po(topix(jxa,jsep),topiy(jxa,jsep)))
    dnsepm(1) = 0.5_R8 * (dna0(jxa,jsep,ismain)+ dna0(topix(jxa,jsep),topiy(jxa,jsep),ismain))
    dpsepm(1) = 0.5_R8 * ( &
         dpa0(jxa,jsep,ismain0)*( &
         rza(jxa,jsep,ismain0)*te(jxa,jsep)+ti(jxa,jsep))+ &
         dpa0(topix(jxa,jsep),topiy(jxa,jsep),ismain0)*( &
         rza(topix(jxa,jsep),topiy(jxa,jsep),ismain0)* &
         te(topix(jxa,jsep),topiy(jxa,jsep))+ &
         ti(topix(jxa,jsep),topiy(jxa,jsep))))
    kesepm(1) = 0.5_R8 * (hce0(jxa,jsep)/ne(jxa,jsep)+ hce0(topix(jxa,jsep),topiy(jxa,jsep))/ &
         ne(topix(jxa,jsep),topiy(jxa,jsep)))
    kisepm(1) = 0.5_R8 * (hci0(jxa,jsep)/ni(jxa,jsep,0) + hci0(topix(jxa,jsep),topiy(jxa,jsep))/ &
         ni(topix(jxa,jsep),topiy(jxa,jsep),0))
    vxsepm(1) = 0.5_R8 * (vla0(jxa,jsep,0,ismain)+ vla0(topix(jxa,jsep),topiy(jxa,jsep),0,ismain))
    vysepm(1) = 0.5_R8 * (vla0(jxa,jsep,1,ismain)+ vla0(topix(jxa,jsep),topiy(jxa,jsep),1,ismain))
    vssepm(1) = 0.5_R8 * ( &
         vsa0(jxa,jsep,ismain)/(mp*am(ismain)*na(jxa,jsep,ismain))+ &
         vsa0(topix(jxa,jsep),topiy(jxa,jsep),ismain)/ &
         (mp*am(ismain)*na(topix(jxa,jsep),topiy(jxa,jsep),ismain)))
    if(nnreg(0).ne.2) then
      nesepa(1) = 0.5_R8 * (ne(nx-target_offset,jsep)+ ne(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep)))
      tesepa(1) = 0.5_R8/ev * (te(nx-target_offset,jsep)+ te(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep)))
      tisepa(1) = 0.5_R8/ev * (ti(nx-target_offset,jsep)+ ti(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep)))
      if(xymap(nx,jsep).gt.0 .and. xymap(topix(nx,jsep),topiy(nx,jsep)).ge.0) then
        tpsepa(1) = 0.5_R8 * (target_temp(xymap(nx,jsep),1)+ target_temp(xymap(topix(nx,jsep),topiy(nx,jsep)),1))
      else
        tpsepa(1) = 0.0_R8
      endif
      posepa(1) = 0.5_R8 *(po(nx-target_offset,jsep)+po(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep)))
    else
      nesepa(1) = ne(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep))
      tesepa(1) = 1.0_R8/ev*te(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep))
      tisepa(1) = 1.0_R8/ev*ti(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep))
      if(xymap(topix(nx,jsep),topiy(nx,jsep)).gt.0) then
        tpsepa(1) = target_temp(xymap(topix(nx,jsep),topiy(nx,jsep)),1)
      else
        tpsepa(1) = 0.0
      endif
      posepa(1) = po(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep))
    endif
    if(nncut.eq.2) then
      nesepa(2) = 0.5_R8 * (ne(ixtr+target_offset,jsep)+ne(topix(ixtr+target_offset,jsep),topiy(ixtr+target_offset,jsep)))
      tesepa(2) = 0.5_R8 * (te(ixtr+target_offset,jsep)+te(topix(ixtr+target_offset,jsep),topiy(ixtr+target_offset,jsep)))/ev
      tisepa(2) = 0.5_R8 * (ti(ixtr+target_offset,jsep)+ti(topix(ixtr+target_offset,jsep),topiy(ixtr+target_offset,jsep)))/ev
      if(xymap(ixtr,jsep).gt.0 .and. xymap(topix(ixtr,jsep),topiy(ixtr,jsep)).gt.0) then
        tpsepa(2) = 0.5_R8 *(target_temp(xymap(ixtr,jsep),1)+target_temp(xymap(topix(ixtr,jsep),topiy(ixtr,jsep)),1))
      else
        tpsepa(2) = 0.0_R8
      endif
      posepa(2) = 0.5_R8 * (po(ixtr+target_offset,jsep)+ po(topix(ixtr+target_offset,jsep),topiy(ixtr+target_offset,jsep)))
      nesepm(2) = 0.5_R8 * (ne(jxi,jsep)+ ne(topix(jxi,jsep),topiy(jxi,jsep)))
      tesepm(2) = 0.5_R8 * (te(jxi,jsep)+ te(topix(jxi,jsep),topiy(jxi,jsep)))/ev
      tisepm(2) = 0.5_R8 * (ti(jxi,jsep)+ ti(topix(jxi,jsep),topiy(jxi,jsep)))/ev
      posepm(2) = 0.5_R8 * (po(jxi,jsep)+ po(topix(jxi,jsep),topiy(jxi,jsep)))
      dnsepm(2) = 0.5_R8 * (dna0(jxi,jsep,ismain)+ dna0(topix(jxi,jsep),topiy(jxi,jsep),ismain))
      dpsepm(2) = 0.5_R8 * ( &
           dpa0(jxi,jsep,ismain0)*( &
           rza(jxi,jsep,ismain0)*te(jxi,jsep)+ti(jxi,jsep))+ &
           dpa0(topix(jxi,jsep),topiy(jxi,jsep),ismain0)*( &
           rza(topix(jxi,jsep),topiy(jxi,jsep),ismain0)* &
           te(topix(jxi,jsep),topiy(jxi,jsep))+ &
           ti(topix(jxi,jsep),topiy(jxi,jsep))))
      kesepm(2) = 0.5_R8 * (hce0(jxi,jsep)/ne(jxi,jsep)+hce0(topix(jxi,jsep),topiy(jxi,jsep))/ne(topix(jxi,jsep),topiy(jxi,jsep)))
      kisepm(2) = 0.5_R8 * (hci0(jxi,jsep)/ni(jxi,jsep,0)+ hci0(topix(jxi,jsep),topiy(jxi,jsep))/ &
           ni(topix(jxi,jsep),topiy(jxi,jsep),0))
      vxsepm(2) = 0.5_R8 * (vla0(jxi,jsep,0,ismain) + vla0(topix(jxi,jsep),topiy(jxi,jsep),0,ismain))
      vysepm(2) = 0.5_R8 * (vla0(jxi,jsep,1,ismain) + vla0(topix(jxi,jsep),topiy(jxi,jsep),1,ismain))
      vssepm(2) = 0.5_R8 * ( &
           vsa0(jxi,jsep,ismain)/(mp*am(ismain)*na(jxi,jsep,ismain)) + vsa0(topix(jxi,jsep),topiy(jxi,jsep),ismain)/ &
           (mp*am(ismain)*na(topix(jxi,jsep),topiy(jxi,jsep),ismain)))
      nesepi(2) = 0.5_R8 * (ne(ixtl-target_offset,jsep)+ ne(topix(ixtl-target_offset,jsep), topiy(ixtl-target_offset,jsep)))
      tesepi(2) = 0.5_R8 * (te(ixtl-target_offset,jsep)+ te(topix(ixtl-target_offset,jsep), topiy(ixtl-target_offset,jsep)))/ev
      tisepi(2) = 0.5_R8 * (ti(ixtl-target_offset,jsep)+ ti(topix(ixtl-target_offset,jsep), topiy(ixtl-target_offset,jsep)))/ev
      if(xymap(ixtl,jsep).gt.0 .and. xymap(topix(ixtl,jsep),topiy(ixtl,jsep)).gt.0) then
        tpsepi(2) = 0.5_R8 * (target_temp(xymap(ixtl,jsep),1)+ target_temp(xymap(topix(ixtl,jsep),topiy(ixtl,jsep)),1))
      else
        tpsepi(2) = 0.0_R8
      endif
      posepi(2) = 0.5_R8 * (po(ixtl-target_offset,jsep) + po(topix(ixtl-target_offset,jsep),topiy(ixtl-target_offset,jsep)))
    endif
    nesepa(nc+1:nncutmax) = 0.0_R8
    tesepa(nc+1:nncutmax) = 0.0_R8
    tisepa(nc+1:nncutmax) = 0.0_R8
    tpsepa(nc+1:nncutmax) = 0.0_R8
    posepa(nc+1:nncutmax) = 0.0_R8
    nesepm(nc+1:nncutmax) = 0.0_R8
    tesepm(nc+1:nncutmax) = 0.0_R8
    tisepm(nc+1:nncutmax) = 0.0_R8
    posepm(nc+1:nncutmax) = 0.0_R8
    dnsepm(nc+1:nncutmax) = 0.0_R8
    dpsepm(nc+1:nncutmax) = 0.0_R8
    kesepm(nc+1:nncutmax) = 0.0_R8
    kisepm(nc+1:nncutmax) = 0.0_R8
    vxsepm(nc+1:nncutmax) = 0.0_R8
    vysepm(nc+1:nncutmax) = 0.0_R8
    vssepm(nc+1:nncutmax) = 0.0_R8
    nesepi(nc+1:nncutmax) = 0.0_R8
    tesepi(nc+1:nncutmax) = 0.0_R8
    tisepi(nc+1:nncutmax) = 0.0_R8
    tpsepi(nc+1:nncutmax) = 0.0_R8
    posepi(nc+1:nncutmax) = 0.0_R8
#endif
    !
    temxip(1:nc) = temxip(1:nc)/ev
    timxip(1:nc) = timxip(1:nc)/ev
    temxap(1:nc) = temxap(1:nc)/ev
    timxap(1:nc) = timxap(1:nc)/ev

    tmne(1)=0.0_R8
    tmte(1)=0.0_R8
    tmti(1)=0.0_R8
    tmvol=0.0_R8
    do iy=-1,ny
      do ix=-1,nx
        if(region(ix,iy,0).ne.0) then
          tmne(1)=tmne(1)+ne(ix,iy)*vol(ix,iy)
          tmte(1)=tmte(1)+te(ix,iy)*ne(ix,iy)*vol(ix,iy)
          tmti(1)=tmti(1)+ti(ix,iy)*ni(ix,iy,0)*vol(ix,iy)
          tmvol=tmvol+vol(ix,iy)
        endif
      enddo
    enddo
    tmne(1)=tmne(1)
    tmte(1)=tmte(1)/ev
    tmti(1)=tmti(1)/ev

#ifndef NO_CDF
    tmhacore(1)=0.0_R8
    tmhasol(1)=0.0_R8
    tmhadiv(1)=0.0_R8
#ifdef B25_EIRENE
    ! note no emission in guard cells and offset of array by 1
    do iy=-1,ny
      do ix=-1,nx
        if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and. bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 ) then
          if(on_closed_surface(ix,iy)) then
            tmhacore(1)=tmhacore(1)+(emiss(ix+1,iy+1,1,1)+emissmol(ix+1,iy+1,1,1))*vol(ix,iy)
          elseif(mod(region(ix,iy,0),4).eq.3 .or.(mod(region(ix,iy,0),4).eq.0 .and. region(ix,iy,0).ne.0)) then
            tmhadiv(1)=tmhadiv(1) + (emiss(ix+1,iy+1,1,1)+emissmol(ix+1,iy+1,1,1))*vol(ix,iy)
          elseif(mod(region(ix,iy,0),4).eq.2 .or. nnreg(0).eq.1) then
            tmhasol(1)=tmhasol(1)+ (emiss(ix+1,iy+1,1,1)+emissmol(ix+1,iy+1,1,1))*vol(ix,iy)
          elseif(region(ix,iy,0).ne.0) then
            write(*,*) 'b2mwti: unknown region @ ', ix,iy,region(ix,iy,0)
          endif
        endif
      enddo
    enddo
#endif

!wdk update batch averages
    if (luav) then
      if (ntim_batch .gt. 0 ) then
        call batch_average(nncutmax,nesepm,nesepm_av,itim,ntim_batch)
        call batch_average(nncutmax,tesepm,tesepm_av,itim,ntim_batch)
        call batch_average(nncutmax,tisepm,tisepm_av,itim,ntim_batch)
        call batch_average(nncutmax,posepm,posepm_av,itim,ntim_batch)
        call batch_average(nncutmax,nesepi,nesepi_av,itim,ntim_batch)
        call batch_average(nncutmax,tesepi,tesepi_av,itim,ntim_batch)
        call batch_average(nncutmax,tisepi,tisepi_av,itim,ntim_batch)
        call batch_average(nncutmax,posepi,posepi_av,itim,ntim_batch)
        call batch_average(nncutmax,nesepa,nesepa_av,itim,ntim_batch)
        call batch_average(nncutmax,tesepa,tesepa_av,itim,ntim_batch)
        call batch_average(nncutmax,tisepa,tisepa_av,itim,ntim_batch)
        call batch_average(nncutmax,posepa,posepa_av,itim,ntim_batch)
        call batch_average(nncutmax,nemxip,nemxip_av,itim,ntim_batch)
        call batch_average(nncutmax,temxip,temxip_av,itim,ntim_batch)
        call batch_average(nncutmax,timxip,timxip_av,itim,ntim_batch)
        call batch_average(nncutmax,pomxip,pomxip_av,itim,ntim_batch)
        call batch_average(nncutmax,nemxap,nemxap_av,itim,ntim_batch)
        call batch_average(nncutmax,temxap,temxap_av,itim,ntim_batch)
        call batch_average(nncutmax,timxap,timxap_av,itim,ntim_batch)
        call batch_average(nncutmax,pomxap,pomxap_av,itim,ntim_batch)
        call batch_average_sq(nncutmax,nesepm,nesepm_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,tesepm,tesepm_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,tisepm,tisepm_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,posepm,posepm_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,nesepi,nesepi_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,tesepi,tesepi_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,tisepi,tisepi_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,posepi,posepi_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,nesepa,nesepa_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,tesepa,tesepa_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,tisepa,tisepa_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,posepa,posepa_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,nemxip,nemxip_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,temxip,temxip_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,timxip,timxip_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,pomxip,pomxip_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,nemxap,nemxap_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,temxap,temxap_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,timxap,timxap_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,pomxap,pomxap_std,itim,ntim_batch)
      endif
    endif
!wdk end of batch averaging

!wdk only write time data if lwti is true
    if (lwti) then
      rw = 'write'
      iret = nf_open(filename, or(NF_WRITE,NF_SHARE), ncid)
      call check_cdf_status(iret)
      imap(1)=1
      tstepn(1) = ntstep
      call rwcdf(rw,ncid,'ntstep',imap,tstepn,iret)
      call rwcdf(rw,ncid,'timesa',imap,timesa,iret)
      if (write_2d .ge. 1) then
        call rwcdf(rw,ncid,'ne2d',(/1,1,1/),ne,iret)
        call rwcdf(rw,ncid,'te2d',(/1,1,1/),te,iret)
        call rwcdf(rw,ncid,'ti2d',(/1,1,1/),ti,iret)
        if (write_2d .ge. 2) then
          call rwcdf(rw,ncid,'po2d',(/1,1,1/),po,iret)
          call rwcdf(rw,ncid,'kin2d',(/1,1,1,1/),kinrgy,iret)
          call rwcdf(rw,ncid,'rsahi2d',(/1,1,1,1/),rsahi,iret)
          call rwcdf(rw,ncid,'rsana2d',(/1,1,1,1/),rsana,iret)
          call rwcdf(rw,ncid,'rrahi2d',(/1,1,1,1/),rrahi,iret)
          call rwcdf(rw,ncid,'rrana2d',(/1,1,1,1/),rrana,iret)
          call rwcdf(rw,ncid,'rcxhi2d',(/1,1,1,1/),rcxhi,iret)
          call rwcdf(rw,ncid,'rcxna2d',(/1,1,1,1/),rcxna,iret)
          call rwcdf(rw,ncid,'rqrad2d',(/1,1,1,1/),rqrad,iret)
          call rwcdf(rw,ncid,'fhe2d',(/1,1,1,1/),fhe,iret)
          call rwcdf(rw,ncid,'fhi2d',(/1,1,1,1/),fhi,iret)
          call rwcdf(rw,ncid,'fch2d',(/1,1,1,1/),fch,iret)
          call rwcdf(rw,ncid,'fna2d',(/1,1,1,1,1/),fna,iret)
        endif
      endif
      imap(1)=1
      imap(2)=1
      call rwcdf(rw,ncid,'fnixip',imap,fnixip,iret)
      call rwcdf(rw,ncid,'feexip',imap,feexip,iret)
      call rwcdf(rw,ncid,'feixip',imap,feixip,iret)
      call rwcdf(rw,ncid,'fetxip',imap,fetxip,iret)
      call rwcdf(rw,ncid,'fchxip',imap,fchxip,iret)
      call rwcdf(rw,ncid,'fnixap',imap,fnixap,iret)
      call rwcdf(rw,ncid,'feexap',imap,feexap,iret)
      call rwcdf(rw,ncid,'feixap',imap,feixap,iret)
      call rwcdf(rw,ncid,'fetxap',imap,fetxap,iret)
      call rwcdf(rw,ncid,'fchxap',imap,fchxap,iret)

      call rwcdf(rw,ncid,'nesepi',imap,nesepi,iret)
      call rwcdf(rw,ncid,'tesepi',imap,tesepi,iret)
      call rwcdf(rw,ncid,'tisepi',imap,tisepi,iret)
      call rwcdf(rw,ncid,'posepi',imap,posepi,iret)
      call rwcdf(rw,ncid,'nesepm',imap,nesepm,iret)
      call rwcdf(rw,ncid,'tesepm',imap,tesepm,iret)
      call rwcdf(rw,ncid,'tisepm',imap,tisepm,iret)
      call rwcdf(rw,ncid,'posepm',imap,posepm,iret)
      call rwcdf(rw,ncid,'dnsepm',imap,dnsepm,iret)
      call rwcdf(rw,ncid,'dpsepm',imap,dpsepm,iret)
      call rwcdf(rw,ncid,'kesepm',imap,kesepm,iret)
      call rwcdf(rw,ncid,'kisepm',imap,kisepm,iret)
      call rwcdf(rw,ncid,'vxsepm',imap,vxsepm,iret)
      call rwcdf(rw,ncid,'vysepm',imap,vysepm,iret)
      call rwcdf(rw,ncid,'vssepm',imap,vssepm,iret)
      call rwcdf(rw,ncid,'nesepa',imap,nesepa,iret)
      call rwcdf(rw,ncid,'tesepa',imap,tesepa,iret)
      call rwcdf(rw,ncid,'tisepa',imap,tisepa,iret)
      call rwcdf(rw,ncid,'posepa',imap,posepa,iret)
      call rwcdf(rw,ncid,'nemxip',imap,nemxip,iret)
      call rwcdf(rw,ncid,'temxip',imap,temxip,iret)
      call rwcdf(rw,ncid,'timxip',imap,timxip,iret)
      call rwcdf(rw,ncid,'pomxip',imap,pomxip,iret)
      call rwcdf(rw,ncid,'nemxap',imap,nemxap,iret)
      call rwcdf(rw,ncid,'temxap',imap,temxap,iret)
      call rwcdf(rw,ncid,'timxap',imap,timxap,iret)
      call rwcdf(rw,ncid,'pomxap',imap,pomxap,iret)
      call rwcdf(rw,ncid,'fniyip',imap,fniyip,iret)
      call rwcdf(rw,ncid,'feeyip',imap,feeyip,iret)
      call rwcdf(rw,ncid,'feiyip',imap,feiyip,iret)
      call rwcdf(rw,ncid,'fetyip',imap,fetyip,iret)
      call rwcdf(rw,ncid,'fchyip',imap,fchyip,iret)
      call rwcdf(rw,ncid,'fniyap',imap,fniyap,iret)
      call rwcdf(rw,ncid,'feeyap',imap,feeyap,iret)
      call rwcdf(rw,ncid,'feiyap',imap,feiyap,iret)
      call rwcdf(rw,ncid,'fetyap',imap,fetyap,iret)
      call rwcdf(rw,ncid,'fchyap',imap,fchyap,iret)
      call rwcdf(rw,ncid,'pwmxip',imap,pwmxip,iret)
      call rwcdf(rw,ncid,'pwmxap',imap,pwmxap,iret)

      imap(1)=1
      call rwcdf(rw,ncid,'tmne',imap,tmne,iret)
      call rwcdf(rw,ncid,'tmte',imap,tmte,iret)
      call rwcdf(rw,ncid,'tmti',imap,tmti,iret)
      call rwcdf(rw,ncid,'tmhacore',imap,tmhacore,iret)
      call rwcdf(rw,ncid,'tmhasol',imap,tmhasol,iret)
      call rwcdf(rw,ncid,'tmhadiv',imap,tmhadiv,iret)

      imap(1)=1
      imap(2)=1
      call rwcdf(rw,ncid,'fnisip',imap,fnisip,iret)
      call rwcdf(rw,ncid,'feesip',imap,feesip,iret)
      call rwcdf(rw,ncid,'feisip',imap,feisip,iret)
      call rwcdf(rw,ncid,'fetsip',imap,fetsip,iret)
      call rwcdf(rw,ncid,'fchsip',imap,fchsip,iret)
      call rwcdf(rw,ncid,'fnisap',imap,fnisap,iret)
      call rwcdf(rw,ncid,'feesap',imap,feesap,iret)
      call rwcdf(rw,ncid,'feisap',imap,feisap,iret)
      call rwcdf(rw,ncid,'fetsap',imap,fetsap,iret)
      call rwcdf(rw,ncid,'fchsap',imap,fchsap,iret)
      call rwcdf(rw,ncid,'fnisipp',imap,fnisipp,iret)
      call rwcdf(rw,ncid,'feesipp',imap,feesipp,iret)
      call rwcdf(rw,ncid,'feisipp',imap,feisipp,iret)
      call rwcdf(rw,ncid,'fetsipp',imap,fetsipp,iret)
      call rwcdf(rw,ncid,'fchsipp',imap,fchsipp,iret)
      call rwcdf(rw,ncid,'fnisapp',imap,fnisapp,iret)
      call rwcdf(rw,ncid,'feesapp',imap,feesapp,iret)
      call rwcdf(rw,ncid,'feisapp',imap,feisapp,iret)
      call rwcdf(rw,ncid,'fetsapp',imap,fetsapp,iret)
      call rwcdf(rw,ncid,'fchsapp',imap,fchsapp,iret)
    !
      imap(1)=nx+2     ! bl
      imap(2)=1
      call rwcdf(rw,ncid,'ne3dl',imap,ne(-1+target_offset,iylstrt),iret)
      call rwcdf(rw,ncid,'te3dl',imap,te(-1+target_offset,iylstrt),iret)
      call rwcdf(rw,ncid,'ti3dl',imap,ti(-1+target_offset,iylstrt),iret)
      call rwcdf(rw,ncid,'po3dl',imap,po(-1+target_offset,iylstrt),iret)
      call rwcdf(rw,ncid,'fn3dl',imap,fna(0,iylstrt,0,ismain),iret)
      call rwcdf(rw,ncid,'fe3dl',imap,fhe(0,iylstrt,0),iret)
      call rwcdf(rw,ncid,'fi3dl',imap,fhi(0,iylstrt,0),iret)
      call rwcdf(rw,ncid,'fc3dl',imap,fch(0,iylstrt,0),iret)
      call rwcdf(rw,ncid,'fl3dl',imap,fne(0,iylstrt,0),iret)
      call rwcdf(rw,ncid,'fo3dl',imap,fni(0,iylstrt,0),iret)
      imap(1)=nx+2     ! i
      imap(2)=1
      call rwcdf(rw,ncid,'ne3di',imap,ne(jxi,iyistrt),iret)
      call rwcdf(rw,ncid,'te3di',imap,te(jxi,iyistrt),iret)
      call rwcdf(rw,ncid,'ti3di',imap,ti(jxi,iyistrt),iret)
      call rwcdf(rw,ncid,'po3di',imap,po(jxi,iyistrt),iret)
      imap(1)=nx+2     ! a
      imap(2)=1
      call rwcdf(rw,ncid,'ne3da',imap,ne(jxa,iyastrt),iret)
      call rwcdf(rw,ncid,'te3da',imap,te(jxa,iyastrt),iret)
      call rwcdf(rw,ncid,'ti3da',imap,ti(jxa,iyastrt),iret)
      call rwcdf(rw,ncid,'po3da',imap,po(jxa,iyastrt),iret)
      imap(1)=nx+2     ! br
      imap(2)=1
      call rwcdf(rw,ncid,'ne3dr',imap,ne(nx-target_offset,iyrstrt),iret)
      call rwcdf(rw,ncid,'te3dr',imap,te(nx-target_offset,iyrstrt),iret)
      call rwcdf(rw,ncid,'ti3dr',imap,ti(nx-target_offset,iyrstrt),iret)
      call rwcdf(rw,ncid,'po3dr',imap,po(nx-target_offset,iyrstrt),iret)
      call rwcdf(rw,ncid,'fn3dr',imap,fna(nx,iyrstrt,0,ismain),iret)
      call rwcdf(rw,ncid,'fe3dr',imap,fhe(nx,iyrstrt,0),iret)
      call rwcdf(rw,ncid,'fi3dr',imap,fhi(nx,iyrstrt,0),iret)
      call rwcdf(rw,ncid,'fc3dr',imap,fch(nx,iyrstrt,0),iret)
      call rwcdf(rw,ncid,'fl3dr',imap,fne(nx,iyrstrt,0),iret)
      call rwcdf(rw,ncid,'fo3dr',imap,fni(nx,iyrstrt,0),iret)
      if (nnreg(0).ge.8) then
        imap(1)=nx+2     ! tr
        imap(2)=1
        call rwcdf(rw,ncid,'ne3dtr',imap,ne(ixtr+target_offset,iytrstrt),iret)
        call rwcdf(rw,ncid,'te3dtr',imap,te(ixtr+target_offset,iytrstrt),iret)
        call rwcdf(rw,ncid,'ti3dtr',imap,ti(ixtr+target_offset,iytrstrt),iret)
        call rwcdf(rw,ncid,'po3dtr',imap,po(ixtr+target_offset,iytrstrt),iret)
        call rwcdf(rw,ncid,'fn3dtr',imap,fna(ixtr+1,iytrstrt,0,ismain),iret)
        call rwcdf(rw,ncid,'fe3dtr',imap,fhe(ixtr+1,iytrstrt,0),iret)
        call rwcdf(rw,ncid,'fi3dtr',imap,fhi(ixtr+1,iytrstrt,0),iret)
        call rwcdf(rw,ncid,'fc3dtr',imap,fch(ixtr+1,iytrstrt,0),iret)
        call rwcdf(rw,ncid,'fl3dtr',imap,fne(ixtr+1,iytrstrt,0),iret)
        call rwcdf(rw,ncid,'fo3dtr',imap,fni(ixtr+1,iytrstrt,0),iret)
        imap(1)=nx+2     ! tl
        imap(2)=1
        call rwcdf(rw,ncid,'ne3dtl',imap,ne(ixtl-target_offset,iytlstrt),iret)
        call rwcdf(rw,ncid,'te3dtl',imap,te(ixtl-target_offset,iytlstrt),iret)
        call rwcdf(rw,ncid,'ti3dtl',imap,ti(ixtl-target_offset,iytlstrt),iret)
        call rwcdf(rw,ncid,'po3dtl',imap,po(ixtl-target_offset,iytlstrt),iret)
        call rwcdf(rw,ncid,'fn3dtl',imap,fna(ixtl,iytlstrt,0,ismain),iret)
        call rwcdf(rw,ncid,'fe3dtl',imap,fhe(ixtl,iytlstrt,0),iret)
        call rwcdf(rw,ncid,'fi3dtl',imap,fhi(ixtl,iytlstrt,0),iret)
        call rwcdf(rw,ncid,'fc3dtl',imap,fch(ixtl,iytlstrt,0),iret)
        call rwcdf(rw,ncid,'fl3dtl',imap,fne(ixtl,iytlstrt,0),iret)
        call rwcdf(rw,ncid,'fo3dtl',imap,fni(ixtl,iytlstrt,0),iret)
      endif
    !
      imap(1)=1
      imap(2)=1
      slice=0.0_R8
      if (ismain0.ne.ismain) then
        slice(iylstrt:iylend)=na(-1+target_offset,iylstrt:iylend,ismain0)
        slice(0:ny-1)=slice(0:ny-1)+dab2(1,1:ny,b2eatcr(ismain0),1)
      endif
      call rwcdf(rw,ncid,'an3dl',imap,slice(iylstrt),iret)
      slice=0.0_R8
      if (ismain0.ne.ismain) then
        slice(iyistrt:iyiend)=na(jxi,iyistrt:iyiend,ismain0)
        slice(0:ny-1)=slice(0:ny-1)+dab2(jxi+1,1:ny,b2eatcr(ismain0),1)
      endif
      call rwcdf(rw,ncid,'an3di',imap,slice(iyistrt),iret)
      slice=0.0_R8
      if (ismain0.ne.ismain) then
        slice(iyastrt:iyaend)=na(jxa,iyastrt:iyaend,ismain0)
        slice(0:ny-1)=slice(0:ny-1)+dab2(jxa+1,1:ny,b2eatcr(ismain0),1)
      endif
      call rwcdf(rw,ncid,'an3da',imap,slice(iyastrt),iret)
      slice=0.0_R8
      if (ismain0.ne.ismain) then
        slice(iyrstrt:iyrend)=na(nx-target_offset,iyrstrt:iyrend,ismain0)
        slice(0:ny-1)=slice(0:ny-1)+dab2(nx,1:ny,b2eatcr(ismain0),1)
      endif
      call rwcdf(rw,ncid,'an3dr',imap,slice(iyrstrt),iret)
      if (nnreg(0).ge.8) then
        slice=0.0_R8
        if (ismain0.ne.ismain) then
          slice(iytlstrt:iytlend)= na(ixtl-target_offset,iytlstrt:iytlend,ismain0)
          slice(0:ny-1)=slice(0:ny-1)+dab2(ixtl,1:ny,b2eatcr(ismain0),1)
        endif
        call rwcdf(rw,ncid,'an3dtl',imap,slice(iytlstrt),iret)
        slice=0.0_R8
        if (ismain0.ne.ismain) then
          slice(iytrstrt:iytrend)= na(ixtr+target_offset,iytrstrt:iytrend,ismain0)
          slice(0:ny-1)=slice(0:ny-1)+dab2(ixtr+1,1:ny,b2eatcr(ismain0),1)
        endif
        call rwcdf(rw,ncid,'an3dtr',imap,slice(iytrstrt),iret)
      endif
      if (nnmoli.gt.0) then
        slice=0.0_R8
        slice(0:ny-1)=dmb2(1,1:ny,1,1)
        call rwcdf(rw,ncid,'mn3dl',imap,slice,iret)
        slice(0:ny-1)=dmb2(jxi+1,1:ny,1,1)
        call rwcdf(rw,ncid,'mn3di',imap,slice,iret)
        slice(0:ny-1)=dmb2(jxa+1,1:ny,1,1)
        call rwcdf(rw,ncid,'mn3da',imap,slice,iret)
        slice(0:ny-1)=dmb2(nx,1:ny,1,1)
        call rwcdf(rw,ncid,'mn3dr',imap,slice,iret)
        if (nnreg(0).ge.8) then
          slice(0:ny-1)=dmb2(ixtl,1:ny,1,1)
          call rwcdf(rw,ncid,'mn3dtl',imap,slice,iret)
          slice(0:ny-1)=dmb2(ixtr+1,1:ny,1,1)
          call rwcdf(rw,ncid,'mn3dtr',imap,slice,iret)
        endif
      endif
    !
      slice=0.0_R8
      do iy = iylstrt, iylend
        call calc_fet(-1,iy,'L',1._R8,nx,ny,ns,ismain,BoRiS,slice(iy))
      enddo
      call rwcdf(rw,ncid,'ft3dl',imap,slice(iylstrt),iret)

      slice=0.0_R8
      do iy = iyrstrt, iyrend
        call calc_fet(nx,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,slice(iy))
      enddo
      call rwcdf(rw,ncid,'ft3dr',imap,slice(iyrstrt),iret)
      if (nnreg(0).ge.8) then
        slice=0.0_R8
        do iy = iytlstrt, iytlend
          call calc_fet(ixtl,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,slice(iy))
        enddo
        call rwcdf(rw,ncid,'ft3dtl',imap,slice(iytlstrt),iret)

        slice=0.0_R8
        do iy = iytrstrt, iytrend
          call calc_fet(ixtr,iy,'L',1._R8,nx,ny,ns,ismain,BoRiS,slice(iy))
        enddo
        call rwcdf(rw,ncid,'ft3dtr',imap,slice(iytrstrt),iret)
      endif
      slice(-1:ny)=dna0(jxi,-1:ny,ismain)
      call rwcdf(rw,ncid,'dn3di',imap,slice,iret)
      slice(-1:ny)=dna0(jxa,-1:ny,ismain)
      call rwcdf(rw,ncid,'dn3da',imap,slice,iret)
      slice(-1:ny)=dpa0(jxi,-1:ny,ismain0)* (rza(jxi,-1:ny,ismain0)*te(jxi,-1:ny)+ti(jxi,-1:ny))
      call rwcdf(rw,ncid,'dp3di',imap,slice,iret)
      slice(-1:ny)=dpa0(jxa,-1:ny,ismain0)* (rza(jxa,-1:ny,ismain0)*te(jxa,-1:ny)+ti(jxa,-1:ny))
      call rwcdf(rw,ncid,'dp3da',imap,slice,iret)
      slice(-1:ny)=fllim0fhi(jxi,-1:ny,1,ismain0)
      call rwcdf(rw,ncid,'lh3di',imap,slice,iret)
      slice(-1:ny)=fllim0fhi(jxa,-1:ny,1,ismain0)
      call rwcdf(rw,ncid,'lh3da',imap,slice,iret)
      slice(-1:ny)=fllim0fna(jxi,-1:ny,1,ismain0)
      call rwcdf(rw,ncid,'ln3di',imap,slice,iret)
      slice(-1:ny)=fllim0fna(jxa,-1:ny,1,ismain0)
      call rwcdf(rw,ncid,'ln3da',imap,slice,iret)
      slice(-1:ny)=hce0(jxi,-1:ny)/ne(jxi,-1:ny)
      call rwcdf(rw,ncid,'ke3di',imap,slice,iret)
      slice(-1:ny)=hce0(jxa,-1:ny)/ne(jxa,-1:ny)
      call rwcdf(rw,ncid,'ke3da',imap,slice,iret)
      slice(-1:ny)=hci0(jxi,-1:ny)/ni(jxi,-1:ny,0)
      call rwcdf(rw,ncid,'ki3di',imap,slice,iret)
      slice(-1:ny)=hci0(jxa,-1:ny)/ni(jxa,-1:ny,0)
      call rwcdf(rw,ncid,'ki3da',imap,slice,iret)
      slice(-1:ny)=vla0(jxi,-1:ny,0,ismain)
      call rwcdf(rw,ncid,'vx3di',imap,slice,iret)
      slice(-1:ny)=vla0(jxa,-1:ny,0,ismain)
      call rwcdf(rw,ncid,'vx3da',imap,slice,iret)
      slice(-1:ny)=vla0(jxi,-1:ny,1,ismain)
      call rwcdf(rw,ncid,'vy3di',imap,slice,iret)
      slice(-1:ny)=vla0(jxa,-1:ny,1,ismain)
      call rwcdf(rw,ncid,'vy3da',imap,slice,iret)
      slice(-1:ny)=vsa0(jxi,-1:ny,ismain)/(mp*am(ismain)*na(jxi,-1:ny,ismain))
      call rwcdf(rw,ncid,'vs3di',imap,slice,iret)
      slice(-1:ny)=vsa0(jxa,-1:ny,ismain)/(mp*am(ismain)*na(jxa,-1:ny,ismain))
      call rwcdf(rw,ncid,'vs3da',imap,slice,iret)
    !
      imap(1)=1
      imap(2)=1
      call rwcdf(rw,ncid,'tpsepi',imap,tpsepi,iret)
      call rwcdf(rw,ncid,'tpsepa',imap,tpsepa,iret)
      call rwcdf(rw,ncid,'tpmxip',imap,tpmxip,iret)
      call rwcdf(rw,ncid,'tpmxap',imap,tpmxap,iret)
      imap(1)=1
      imap(2)=1
      slice=0.0_R8
      if(minval(xymap(-1,0:ny-1)).gt.0) then
        slice(-1)=target_temp(xymap(-1,0),1)
        slice(0:ny-1)=target_temp(xymap(-1,0:ny-1),1)
        slice(ny)=target_temp(xymap(-1,ny-1),1)
      endif
      call rwcdf(rw,ncid,'tp3dl',imap,slice,iret)
      slice=0.0_R8
      if(minval(xymap(nx,0:ny-1)).gt.0) then
        slice(-1)=target_temp(xymap(nx,0),1)
        slice(0:ny-1)=target_temp(xymap(nx,0:ny-1),1)
        slice(ny)=target_temp(xymap(nx,ny-1),1)
      endif
      call rwcdf(rw,ncid,'tp3dr',imap,slice,iret)
      if (nnreg(0).ge.8) then
        slice=0.0_R8
        if(minval(xymap(ixtl,0:ny-1)).gt.0) then
          slice(-1)=target_temp(xymap(ixtl,0),1)
          slice(0:ny-1)=target_temp(xymap(ixtl,0:ny-1),1)
          slice(ny)=target_temp(xymap(ixtl,ny-1),1)
        endif
        call rwcdf(rw,ncid,'tp3dtl',imap,slice,iret)
        slice=0.0_R8
        if(minval(xymap(ixtr,0:ny-1)).gt.0) then
          slice(-1)=target_temp(xymap(ixtr,0),1)
          slice(0:ny-1)=target_temp(xymap(ixtr,0:ny-1),1)
          slice(ny)=target_temp(xymap(ixtr,ny-1),1)
        endif
        call rwcdf(rw,ncid,'tp3dtr',imap,slice,iret)
      endif

      iret = nf_close(ncid)
      call check_cdf_status(iret)
    endif

!wdk only write batch data if lwav is true
    if (lwav) then
      rw = 'write'
      iret = nf_open(filename_av, or(NF_WRITE,NF_SHARE), ncid)
      call check_cdf_status(iret)
!wdk compute the standard deviation from average and average of squares
      fac = rratio(ntim_batch,ntim_batch - 1)
      nesepm_std = ((nesepm_std - nesepm_av**2)*fac)**0.5
      tesepm_std = ((tesepm_std - tesepm_av**2)*fac)**0.5
      tisepm_std = ((tisepm_std - tisepm_av**2)*fac)**0.5
      posepm_std = ((posepm_std - posepm_av**2)*fac)**0.5
      nesepi_std = ((nesepi_std - nesepi_av**2)*fac)**0.5
      tesepi_std = ((tesepi_std - tesepi_av**2)*fac)**0.5
      tisepi_std = ((tisepi_std - tisepi_av**2)*fac)**0.5
      posepi_std = ((posepi_std - posepi_av**2)*fac)**0.5
      nesepa_std = ((nesepa_std - nesepa_av**2)*fac)**0.5
      tesepa_std = ((tesepa_std - tesepa_av**2)*fac)**0.5
      tisepa_std = ((tisepa_std - tisepa_av**2)*fac)**0.5
      posepa_std = ((posepa_std - posepa_av**2)*fac)**0.5
      nemxip_std = ((nemxip_std - nemxip_av**2)*fac)**0.5
      temxip_std = ((temxip_std - temxip_av**2)*fac)**0.5
      timxip_std = ((timxip_std - timxip_av**2)*fac)**0.5
      pomxip_std = ((pomxip_std - pomxip_av**2)*fac)**0.5
      nemxap_std = ((nemxap_std - nemxap_av**2)*fac)**0.5
      temxap_std = ((temxap_std - temxap_av**2)*fac)**0.5
      timxap_std = ((timxap_std - timxap_av**2)*fac)**0.5
      pomxap_std = ((pomxap_std - pomxap_av**2)*fac)**0.5

!wdk write into b2time.nc
      imap(1)=1
      tstepn(1) = nastep
      call rwcdf(rw,ncid,'nastep',imap,tstepn,iret)
      call rwcdf(rw,ncid,'batchsa',imap,batchsa,iret)

!wdk averages
      imap(1)=1
      imap(2)=nncutmax
      call rwcdf(rw,ncid,'nesepm_av',imap,nesepm_av,iret)
      call rwcdf(rw,ncid,'tesepm_av',imap,tesepm_av,iret)
      call rwcdf(rw,ncid,'tisepm_av',imap,tisepm_av,iret)
      call rwcdf(rw,ncid,'posepm_av',imap,posepm_av,iret)
      call rwcdf(rw,ncid,'nesepi_av',imap,nesepi_av,iret)
      call rwcdf(rw,ncid,'tesepi_av',imap,tesepi_av,iret)
      call rwcdf(rw,ncid,'tisepi_av',imap,tisepi_av,iret)
      call rwcdf(rw,ncid,'posepi_av',imap,posepi_av,iret)
      call rwcdf(rw,ncid,'nesepa_av',imap,nesepa_av,iret)
      call rwcdf(rw,ncid,'tesepa_av',imap,tesepa_av,iret)
      call rwcdf(rw,ncid,'tisepa_av',imap,tisepa_av,iret)
      call rwcdf(rw,ncid,'posepa_av',imap,posepa_av,iret)
      call rwcdf(rw,ncid,'nemxip_av',imap,nemxip_av,iret)
      call rwcdf(rw,ncid,'temxip_av',imap,temxip_av,iret)
      call rwcdf(rw,ncid,'timxip_av',imap,timxip_av,iret)
      call rwcdf(rw,ncid,'pomxip_av',imap,pomxip_av,iret)
      call rwcdf(rw,ncid,'nemxap_av',imap,nemxap_av,iret)
      call rwcdf(rw,ncid,'temxap_av',imap,temxap_av,iret)
      call rwcdf(rw,ncid,'timxap_av',imap,timxap_av,iret)
      call rwcdf(rw,ncid,'pomxap_av',imap,pomxap_av,iret)

!wdk standard deviations
      call rwcdf(rw,ncid,'nesepm_std',imap,nesepm_std,iret)
      call rwcdf(rw,ncid,'tesepm_std',imap,tesepm_std,iret)
      call rwcdf(rw,ncid,'tisepm_std',imap,tisepm_std,iret)
      call rwcdf(rw,ncid,'posepm_std',imap,posepm_std,iret)
      call rwcdf(rw,ncid,'nesepi_std',imap,nesepi_std,iret)
      call rwcdf(rw,ncid,'tesepi_std',imap,tesepi_std,iret)
      call rwcdf(rw,ncid,'tisepi_std',imap,tisepi_std,iret)
      call rwcdf(rw,ncid,'posepi_std',imap,posepi_std,iret)
      call rwcdf(rw,ncid,'nesepa_std',imap,nesepa_std,iret)
      call rwcdf(rw,ncid,'tesepa_std',imap,tesepa_std,iret)
      call rwcdf(rw,ncid,'tisepa_std',imap,tisepa_std,iret)
      call rwcdf(rw,ncid,'posepa_std',imap,posepa_std,iret)
      call rwcdf(rw,ncid,'nemxip_std',imap,nemxip_std,iret)
      call rwcdf(rw,ncid,'temxip_std',imap,temxip_std,iret)
      call rwcdf(rw,ncid,'timxip_std',imap,timxip_std,iret)
      call rwcdf(rw,ncid,'pomxip_std',imap,pomxip_std,iret)
      call rwcdf(rw,ncid,'nemxap_std',imap,nemxap_std,iret)
      call rwcdf(rw,ncid,'temxap_std',imap,temxap_std,iret)
      call rwcdf(rw,ncid,'timxap_std',imap,timxap_std,iret)
      call rwcdf(rw,ncid,'pomxap_std',imap,pomxap_std,iret)
      iret = nf_close(ncid)
      call check_cdf_status(iret)
    endif
#endif

    ! ..return
    ncall = ncall+1
    call subend ()
    return

    !-----------------------------------------------------------------------
    !.end b2mwti

  end subroutine b2mwti

  subroutine dealloc_b2mod_mwti

  if (.not.allocated(nesepi_av)) return

  deallocate(nesepi_av)
  deallocate(tesepi_av)
  deallocate(tisepi_av)
  deallocate(posepi_av)
  deallocate(nesepm_av)
  deallocate(tesepm_av)
  deallocate(tisepm_av)
  deallocate(posepm_av)
  deallocate(nesepa_av)
  deallocate(tesepa_av)
  deallocate(tisepa_av)
  deallocate(posepa_av)
  deallocate(nemxip_av)
  deallocate(temxip_av)
  deallocate(timxip_av)
  deallocate(pomxip_av)
  deallocate(nemxap_av)
  deallocate(temxap_av)
  deallocate(timxap_av)
  deallocate(pomxap_av)
  deallocate(nesepi_std)
  deallocate(tesepi_std)
  deallocate(tisepi_std)
  deallocate(posepi_std)
  deallocate(nesepm_std)
  deallocate(tesepm_std)
  deallocate(tisepm_std)
  deallocate(posepm_std)
  deallocate(nesepa_std)
  deallocate(tesepa_std)
  deallocate(tisepa_std)
  deallocate(posepa_std)
  deallocate(nemxip_std)
  deallocate(temxip_std)
  deallocate(timxip_std)
  deallocate(pomxip_std)
  deallocate(nemxap_std)
  deallocate(temxap_std)
  deallocate(timxap_std)
  deallocate(pomxap_std)

  return
  end subroutine dealloc_b2mod_mwti

#ifndef NO_CDF
  subroutine b2crtimecdf(filename, &
   nx, ny, nybl, nytl, nytr, nybr, nya, nyi, nc, ns, write_2d, &
   ncid, batch_only, iret)
    use b2mod_constants
#     include <netcdf.inc>
    integer nx, ny, nybl, nytl, nytr, nybr, nya, nyi, nc, ns, iret
    integer, intent(in) :: write_2d
    logical, intent(in) :: batch_only
    character*256 :: filename
    ! NetCDF id
    integer  ncid
    ! dimension ids
    integer :: nxdim, nydim, nsdim, timedim, batchdim, &
         nybldim, nytldim, nytrdim, nybrdim, nyadim, nyidim, ncdim, idirdim
    ! variable ids
    integer :: ntstepid, timesaid, fnixipid, feexipid, feixipid, &
         fnixapid, feexapid, feixapid, nesepiid, tesepiid, tisepiid, &
         nesepmid, tesepmid, tisepmid, nesepaid, tesepaid, tisepaid, &
         nemxipid, temxipid, timxipid, nemxapid, temxapid, timxapid, &
         fniyipid, feeyipid, feiyipid, fniyapid, feeyapid, feiyapid, &
         pwmxipid, pwmxapid, tmneid, tmteid, tmtiid, tmhacoreid, &
         tmhasolid, tmhadivid, fnisipid, feesipid, feisipid, fnisapid, &
         feesapid, feisapid, fnisippid, feesippid, feisippid, fnisappid, &
         feesappid, feisappid, ne3dlid, te3dlid, ti3dlid, ne3diid, &
         te3diid, ti3diid, ne3daid, te3daid, ti3daid, ne3drid, te3drid, &
         ti3drid, an3dlid, mn3dlid, an3diid, mn3diid, an3daid, mn3daid, &
         an3drid, mn3drid, fn3dlid, fe3dlid, fi3dlid, fn3drid, fe3drid, &
         fi3drid, ne3dtlid, te3dtlid, ti3dtlid, ne3dtrid, te3dtrid, &
         ti3dtrid, an3dtlid, mn3dtlid, &
         an3dtrid, mn3dtrid, fn3dtlid, fe3dtlid, fi3dtlid, fn3dtrid, &
         fe3dtrid, fi3dtrid, fetxipid, fetxapid, fetyipid, fetyapid, &
         fetsipid, fetsapid, fetsippid, fetsappid
    integer :: ne2did, te2did, ti2did, po2did, kin2did, rsahi2did, &
         rsana2did, rrahi2did, rrana2did, rcxhi2did, rcxna2did, rqrad2did, &
         rqahe2did, fch2did, fhe2did, fhi2did, fna2did
    integer :: fchxipid, fchxapid, posepiid, posepmid, posepaid, &
         pomxipid, pomxapid, fchyipid, fchyapid, &
         fchsipid, fchsapid, fchsippid, fchsappid, po3dlid, &
         po3diid, po3daid, po3drid, fc3dlid, fc3drid, &
         fl3dlid, fl3drid, fo3dlid, fo3drid, &
         ft3dlid, ft3drid, po3dtlid, po3dtrid, fc3dtlid, fc3dtrid, &
         fl3dtlid, fl3dtrid, fo3dtlid, fo3dtrid, &
         ft3dtlid, ft3dtrid, &
         dn3diid, dn3daid, dp3diid, dp3daid, ke3diid, ke3daid, &
         ki3diid, ki3daid, vx3diid, vx3daid, vy3diid, vy3daid, &
         vs3diid, vs3daid, lh3diid, lh3daid, ln3diid, ln3daid, &
         dnsepmid, dpsepmid, kesepmid, kisepmid, &
         vxsepmid, vysepmid, vssepmid, &
         tpmxipid, tpmxapid, tp3drid, tp3dlid, tp3dtlid, tp3dtrid, &
         tpsepiid, tpsepaid, &
         nastepid, ntimbatchid, batchsaid, &
         nesepm_avid, tesepm_avid, tisepm_avid, posepm_avid, &
         nesepi_avid, tesepi_avid, tisepi_avid, posepi_avid, &
         nesepa_avid, tesepa_avid, tisepa_avid, posepa_avid, &
         nemxip_avid, temxip_avid, timxip_avid, pomxip_avid, &
         nemxap_avid, temxap_avid, timxap_avid, pomxap_avid, &
         nesepm_stdid, tesepm_stdid, tisepm_stdid, posepm_stdid, &
         nesepi_stdid, tesepi_stdid, tisepi_stdid, posepi_stdid, &
         nesepa_stdid, tesepa_stdid, tisepa_stdid, posepa_stdid, &
         nemxip_stdid, temxip_stdid, timxip_stdid, pomxip_stdid, &
         nemxap_stdid, temxap_stdid, timxap_stdid, pomxap_stdid
    ! variable shapes
    integer :: dims(2)
    real (kind=R8) :: dvals(1)
    ! CDF format variable
    integer, save :: cdf_default = 0    ! used for setting a default NetCDF format
    ! Create and enter define mode
    call ipgeti ('b2mndr_cdf_default', cdf_default)
    if (cdf_default.eq.3 .or. cdf_default.eq.4) then
      iret = nf_create(trim(filename), or(ncclob,nf_netcdf4), ncid)
    else
      iret = nf_create(trim(filename), ncclob, ncid)
    end if
    call check_cdf_status(iret)
    ! define dimensions
    if (.not.batch_only) then
      iret = nf_def_dim(ncid, 'nx', nx+2, nxdim)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'ny', ny+2, nydim)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'nybl', nybl, nybldim)
      call check_cdf_status(iret)
      if(nytl.gt.0) then
        iret = nf_def_dim(ncid, 'nytl', nytl, nytldim)
        call check_cdf_status(iret)
      endif
      if(nytr.gt.0) then
        iret = nf_def_dim(ncid, 'nytr', nytr, nytrdim)
        call check_cdf_status(iret)
      endif
      iret = nf_def_dim(ncid, 'nybr', nybr, nybrdim)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'nyi', nyi, nyidim)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'nya', nya, nyadim)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'nc', nc, ncdim)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'ns', ns, nsdim)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'time', ncunlim, timedim)
      call check_cdf_status(iret)
    else
      iret = nf_def_dim(ncid, 'nc', nc, ncdim)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'batch', ncunlim, batchdim)
      call check_cdf_status(iret)
    end if
    ! define variables
    if (.not.batch_only) then
      dims(1) = 0
      iret = nf_def_var(ncid, 'ntstep', NCDOUBLE, 0, dims, ntstepid)
      call check_cdf_status(iret)
      dims(1) = timedim
      iret = nf_def_var(ncid, 'timesa', NCDOUBLE, 1, dims, timesaid)
      call check_cdf_status(iret)
      dvals(1) = 1.0_R8/ev
      if (write_2d .ge. 1) then
        iret = nf_def_var(ncid, 'ne2d', NCDOUBLE, 3, (/nxdim,nydim,timedim/), ne2did)
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ne2did, 'long_name', 2, 'ne')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ne2did, 'units', 4, 'm^-3')
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'te2d', NCDOUBLE, 3, (/nxdim,nydim,timedim/), te2did)
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, te2did, 'long_name', 2, 'Te')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, te2did, 'units', 2, 'eV')
        call check_cdf_status(iret)
        iret = nf_put_att_double(ncid, te2did, 'scale', NCDOUBLE, 1, dvals(1))
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'ti2d', NCDOUBLE, 3, (/nxdim,nydim,timedim/), ti2did)
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ti2did, 'long_name', 2, 'Ti')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ti2did, 'units', 2, 'eV')
        call check_cdf_status(iret)
        iret = nf_put_att_double(ncid, ti2did, 'scale', NCDOUBLE, 1, dvals(1))
        call check_cdf_status(iret)
        if (write_2d .ge. 2) then
          iret = nf_def_var(ncid, 'po2d', NCDOUBLE, 3, (/nxdim,nydim,timedim/), po2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, po2did, 'long_name', 9, 'potential')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, po2did, 'units', 1, 'V')
          call check_cdf_status(iret)
          iret = nf_def_var(ncid, 'kin2d', NCDOUBLE, 4, (/nxdim,nydim,nsdim,timedim/), kin2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, kin2did, 'long_name', 23, 'parallel kinetic energy')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, kin2did, 'units', 1, 'J')
          call check_cdf_status(iret)
          iret = nf_def_var(ncid, 'rsahi2d', NCDOUBLE, 4, (/nxdim,nydim,nsdim,timedim/), rsahi2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rsahi2did, 'long_name', 21, 'iz energy source/sink')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rsahi2did, 'units', 1, 'W')
          call check_cdf_status(iret)
          iret = nf_def_var(ncid, 'rsana2d', NCDOUBLE, 4, (/nxdim,nydim,nsdim,timedim/), rsana2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rsana2did, 'long_name', 7, 'iz rate')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rsana2did, 'units', 3, '1/s')
          call check_cdf_status(iret)
          iret = nf_def_var(ncid, 'rrahi2d', NCDOUBLE, 4, (/nxdim,nydim,nsdim,timedim/), rrahi2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rrahi2did, 'long_name', 21, 'rc energy source/sink')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rrahi2did, 'units', 1, 'W')
          call check_cdf_status(iret)
          iret = nf_def_var(ncid, 'rrana2d', NCDOUBLE, 4, (/nxdim,nydim,nsdim,timedim/), rrana2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rrana2did, 'long_name', 7, 'rc rate')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rrana2did, 'units', 3, '1/s')
          call check_cdf_status(iret)
          iret = nf_def_var(ncid, 'rcxhi2d', NCDOUBLE, 4, (/nxdim,nydim,nsdim,timedim/), rcxhi2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rcxhi2did, 'long_name', 21, 'cx energy source/sink')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rcxhi2did, 'units', 1, 'W')
          call check_cdf_status(iret)
          iret = nf_def_var(ncid, 'rcxna2d', NCDOUBLE, 4, (/nxdim,nydim,nsdim,timedim/), rcxna2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rcxna2did, 'long_name', 7, 'cx rate')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rcxna2did, 'units', 3, '1/s')
          call check_cdf_status(iret)
          iret = nf_def_var(ncid, 'rqrad2d', NCDOUBLE, 4, (/nxdim,nydim,nsdim,timedim/), rqrad2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rqrad2did, 'long_name', 19, 'Line radiation rate')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rqrad2did, 'units', 1, 'W')
          call check_cdf_status(iret)
          iret = nf_def_var(ncid, 'rqahe2d', NCDOUBLE, 4, (/nxdim,nydim,nsdim,timedim/), rqahe2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rqahe2did, 'long_name', 21, 'Electron cooling rate')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, rqahe2did, 'units', 1, 'W')
          call check_cdf_status(iret)

          iret = nf_def_dim(ncid, 'idir', 2, idirdim)  ! Needed for fluxes
          call check_cdf_status(iret)
          iret = nf_def_var(ncid, 'fhe2d', NCDOUBLE, 4, (/nxdim,nydim,idirdim,timedim/), fhe2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, fhe2did, 'long_name', 18, 'Electron heat flux')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, fhe2did, 'units', 1, 'W')
          call check_cdf_status(iret)
          iret = nf_def_var(ncid, 'fhi2d', NCDOUBLE, 4, (/nxdim,nydim,idirdim,timedim/), fhi2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, fhi2did, 'long_name', 13, 'Ion heat flux')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, fhi2did, 'units', 1, 'W')
          call check_cdf_status(iret)
          iret = nf_def_var(ncid, 'fch2d', NCDOUBLE, 4, (/nxdim,nydim,idirdim,timedim/), fch2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, fch2did, 'long_name', 7, 'Current')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, fch2did, 'units', 1, 'A')
          call check_cdf_status(iret)
          iret = nf_def_var(ncid, 'fna2d', NCDOUBLE, 5, (/nxdim,nydim,idirdim,nsdim,timedim/), fna2did)
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, fna2did, 'long_name', 13, 'Particle flux')
          call check_cdf_status(iret)
          iret = nf_put_att_text(ncid, fna2did, 'units', 3, '1/s')
          call check_cdf_status(iret)
        endif
      endif

      dims(1) = ncdim
      dims(2) = timedim
      iret = nf_def_var(ncid, 'fnixip', NCDOUBLE, 2, dims, fnixipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feexip', NCDOUBLE, 2, dims, feexipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feixip', NCDOUBLE, 2, dims, feixipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fetxip', NCDOUBLE, 2, dims, fetxipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fchxip', NCDOUBLE, 2, dims, fchxipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fnixap', NCDOUBLE, 2, dims, fnixapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feexap', NCDOUBLE, 2, dims, feexapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feixap', NCDOUBLE, 2, dims, feixapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fetxap', NCDOUBLE, 2, dims, fetxapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fchxap', NCDOUBLE, 2, dims, fchxapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'nesepi', NCDOUBLE, 2, dims, nesepiid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tesepi', NCDOUBLE, 2, dims, tesepiid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tisepi', NCDOUBLE, 2, dims, tisepiid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tpsepi', NCDOUBLE, 2, dims, tpsepiid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'posepi', NCDOUBLE, 2, dims, posepiid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'nesepm', NCDOUBLE, 2, dims, nesepmid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tesepm', NCDOUBLE, 2, dims, tesepmid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tisepm', NCDOUBLE, 2, dims, tisepmid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'posepm', NCDOUBLE, 2, dims, posepmid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'dnsepm', NCDOUBLE, 2, dims, dnsepmid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'dpsepm', NCDOUBLE, 2, dims, dpsepmid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'kesepm', NCDOUBLE, 2, dims, kesepmid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'kisepm', NCDOUBLE, 2, dims, kisepmid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'vxsepm', NCDOUBLE, 2, dims, vxsepmid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'vysepm', NCDOUBLE, 2, dims, vysepmid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'vssepm', NCDOUBLE, 2, dims, vssepmid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'nesepa', NCDOUBLE, 2, dims, nesepaid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tesepa', NCDOUBLE, 2, dims, tesepaid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tisepa', NCDOUBLE, 2, dims, tisepaid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tpsepa', NCDOUBLE, 2, dims, tpsepaid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'posepa', NCDOUBLE, 2, dims, posepaid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'nemxip', NCDOUBLE, 2, dims, nemxipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'temxip', NCDOUBLE, 2, dims, temxipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'timxip', NCDOUBLE, 2, dims, timxipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tpmxip', NCDOUBLE, 2, dims, tpmxipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'pomxip', NCDOUBLE, 2, dims, pomxipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'nemxap', NCDOUBLE, 2, dims, nemxapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'temxap', NCDOUBLE, 2, dims, temxapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'timxap', NCDOUBLE, 2, dims, timxapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tpmxap', NCDOUBLE, 2, dims, tpmxapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'pomxap', NCDOUBLE, 2, dims, pomxapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fniyip', NCDOUBLE, 2, dims, fniyipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feeyip', NCDOUBLE, 2, dims, feeyipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feiyip', NCDOUBLE, 2, dims, feiyipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fetyip', NCDOUBLE, 2, dims, fetyipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fchyip', NCDOUBLE, 2, dims, fchyipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fniyap', NCDOUBLE, 2, dims, fniyapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feeyap', NCDOUBLE, 2, dims, feeyapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feiyap', NCDOUBLE, 2, dims, feiyapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fetyap', NCDOUBLE, 2, dims, fetyapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fchyap', NCDOUBLE, 2, dims, fchyapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'pwmxip', NCDOUBLE, 2, dims, pwmxipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'pwmxap', NCDOUBLE, 2, dims, pwmxapid)
      call check_cdf_status(iret)
      dims(1) = timedim
      iret = nf_def_var(ncid, 'tmne', NCDOUBLE, 1, dims, tmneid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tmte', NCDOUBLE, 1, dims, tmteid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tmti', NCDOUBLE, 1, dims, tmtiid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tmhacore', NCDOUBLE, 1, dims, tmhacoreid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tmhasol', NCDOUBLE, 1, dims, tmhasolid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tmhadiv', NCDOUBLE, 1, dims, tmhadivid)
      call check_cdf_status(iret)
      dims(1) = ncdim
      dims(2) = timedim
      iret = nf_def_var(ncid, 'fnisip', NCDOUBLE, 2, dims, fnisipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feesip', NCDOUBLE, 2, dims, feesipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feisip', NCDOUBLE, 2, dims, feisipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fetsip', NCDOUBLE, 2, dims, fetsipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fchsip', NCDOUBLE, 2, dims, fchsipid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fnisap', NCDOUBLE, 2, dims, fnisapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feesap', NCDOUBLE, 2, dims, feesapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feisap', NCDOUBLE, 2, dims, feisapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fetsap', NCDOUBLE, 2, dims, fetsapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fchsap', NCDOUBLE, 2, dims, fchsapid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fnisipp', NCDOUBLE, 2, dims, fnisippid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feesipp', NCDOUBLE, 2, dims, feesippid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feisipp', NCDOUBLE, 2, dims, feisippid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fetsipp', NCDOUBLE, 2, dims, fetsippid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fchsipp', NCDOUBLE, 2, dims, fchsippid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fnisapp', NCDOUBLE, 2, dims, fnisappid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feesapp', NCDOUBLE, 2, dims, feesappid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'feisapp', NCDOUBLE, 2, dims, feisappid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fetsapp', NCDOUBLE, 2, dims, fetsappid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fchsapp', NCDOUBLE, 2, dims, fchsappid)
      call check_cdf_status(iret)
      dims(2) = timedim
      dims(1) = nybldim
      iret = nf_def_var(ncid, 'fn3dl', NCDOUBLE, 2, dims, fn3dlid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fe3dl', NCDOUBLE, 2, dims, fe3dlid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fi3dl', NCDOUBLE, 2, dims, fi3dlid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'ft3dl', NCDOUBLE, 2, dims, ft3dlid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fc3dl', NCDOUBLE, 2, dims, fc3dlid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fl3dl', NCDOUBLE, 2, dims, fl3dlid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fo3dl', NCDOUBLE, 2, dims, fo3dlid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'ne3dl', NCDOUBLE, 2, dims, ne3dlid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'te3dl', NCDOUBLE, 2, dims, te3dlid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'ti3dl', NCDOUBLE, 2, dims, ti3dlid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'po3dl', NCDOUBLE, 2, dims, po3dlid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'an3dl', NCDOUBLE, 2, dims, an3dlid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'mn3dl', NCDOUBLE, 2, dims, mn3dlid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tp3dl', NCDOUBLE, 2, dims, tp3dlid)
      call check_cdf_status(iret)
      if(nytl.gt.0) then
        dims(2) = timedim
        dims(1) = nytldim
        iret = nf_def_var(ncid, 'fn3dtl', NCDOUBLE, 2, dims, fn3dtlid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'fe3dtl', NCDOUBLE, 2, dims, fe3dtlid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'fi3dtl', NCDOUBLE, 2, dims, fi3dtlid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'ft3dtl', NCDOUBLE, 2, dims, ft3dtlid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'fc3dtl', NCDOUBLE, 2, dims, fc3dtlid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'fl3dtl', NCDOUBLE, 2, dims, fl3dtlid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'fo3dtl', NCDOUBLE, 2, dims, fo3dtlid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'ne3dtl', NCDOUBLE, 2, dims, ne3dtlid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'te3dtl', NCDOUBLE, 2, dims, te3dtlid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'ti3dtl', NCDOUBLE, 2, dims, ti3dtlid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'po3dtl', NCDOUBLE, 2, dims, po3dtlid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'an3dtl', NCDOUBLE, 2, dims, an3dtlid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'mn3dtl', NCDOUBLE, 2, dims, mn3dtlid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'tp3dtl', NCDOUBLE, 2, dims, tp3dtlid)
        call check_cdf_status(iret)
      endif
      dims(2) = timedim
      dims(1) = nyidim
      iret = nf_def_var(ncid, 'ne3di', NCDOUBLE, 2, dims, ne3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'te3di', NCDOUBLE, 2, dims, te3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'ti3di', NCDOUBLE, 2, dims, ti3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'po3di', NCDOUBLE, 2, dims, po3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'an3di', NCDOUBLE, 2, dims, an3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'mn3di', NCDOUBLE, 2, dims, mn3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'dn3di', NCDOUBLE, 2, dims, dn3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'dp3di', NCDOUBLE, 2, dims, dp3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'lh3di', NCDOUBLE, 2, dims, lh3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'ln3di', NCDOUBLE, 2, dims, ln3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'ke3di', NCDOUBLE, 2, dims, ke3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'ki3di', NCDOUBLE, 2, dims, ki3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'vx3di', NCDOUBLE, 2, dims, vx3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'vy3di', NCDOUBLE, 2, dims, vy3diid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'vs3di', NCDOUBLE, 2, dims, vs3diid)
      call check_cdf_status(iret)
      dims(2) = timedim
      dims(1) = nyadim
      iret = nf_def_var(ncid, 'ne3da', NCDOUBLE, 2, dims, ne3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'te3da', NCDOUBLE, 2, dims, te3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'ti3da', NCDOUBLE, 2, dims, ti3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'po3da', NCDOUBLE, 2, dims, po3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'an3da', NCDOUBLE, 2, dims, an3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'mn3da', NCDOUBLE, 2, dims, mn3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'dn3da', NCDOUBLE, 2, dims, dn3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'dp3da', NCDOUBLE, 2, dims, dp3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'lh3da', NCDOUBLE, 2, dims, lh3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'ln3da', NCDOUBLE, 2, dims, ln3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'ke3da', NCDOUBLE, 2, dims, ke3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'ki3da', NCDOUBLE, 2, dims, ki3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'vx3da', NCDOUBLE, 2, dims, vx3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'vy3da', NCDOUBLE, 2, dims, vy3daid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'vs3da', NCDOUBLE, 2, dims, vs3daid)
      call check_cdf_status(iret)
      dims(2) = timedim
      dims(1) = nybrdim
      iret = nf_def_var(ncid, 'ne3dr', NCDOUBLE, 2, dims, ne3drid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'te3dr', NCDOUBLE, 2, dims, te3drid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'ti3dr', NCDOUBLE, 2, dims, ti3drid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'po3dr', NCDOUBLE, 2, dims, po3drid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'an3dr', NCDOUBLE, 2, dims, an3drid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'mn3dr', NCDOUBLE, 2, dims, mn3drid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fn3dr', NCDOUBLE, 2, dims, fn3drid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fe3dr', NCDOUBLE, 2, dims, fe3drid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fi3dr', NCDOUBLE, 2, dims, fi3drid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'ft3dr', NCDOUBLE, 2, dims, ft3drid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fc3dr', NCDOUBLE, 2, dims, fc3drid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fl3dr', NCDOUBLE, 2, dims, fl3drid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'fo3dr', NCDOUBLE, 2, dims, fo3drid)
      call check_cdf_status(iret)
      iret = nf_def_var(ncid, 'tp3dr', NCDOUBLE, 2, dims, tp3drid)
      call check_cdf_status(iret)
      if(nytr.gt.0) then
        dims(2) = timedim
        dims(1) = nytrdim
        iret = nf_def_var(ncid, 'ne3dtr', NCDOUBLE, 2, dims, ne3dtrid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'te3dtr', NCDOUBLE, 2, dims, te3dtrid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'ti3dtr', NCDOUBLE, 2, dims, ti3dtrid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'po3dtr', NCDOUBLE, 2, dims, po3dtrid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'an3dtr', NCDOUBLE, 2, dims, an3dtrid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'mn3dtr', NCDOUBLE, 2, dims, mn3dtrid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'fn3dtr', NCDOUBLE, 2, dims, fn3dtrid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'fe3dtr', NCDOUBLE, 2, dims, fe3dtrid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'fi3dtr', NCDOUBLE, 2, dims, fi3dtrid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'ft3dtr', NCDOUBLE, 2, dims, ft3dtrid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'fc3dtr', NCDOUBLE, 2, dims, fc3dtrid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'fl3dtr', NCDOUBLE, 2, dims, fl3dtrid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'fo3dtr', NCDOUBLE, 2, dims, fo3dtrid)
        call check_cdf_status(iret)
        iret = nf_def_var(ncid, 'tp3dtr', NCDOUBLE, 2, dims, tp3dtrid)
        call check_cdf_status(iret)
      endif
    endif

    !wdk averages
    if (batch_only) then
      dims(1) = 0
      iret  = nf_def_var(ncid, 'nastep', NCDOUBLE, 0, dims, nastepid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'ntim_batch', NCDOUBLE, 0, dims, ntimbatchid)
      call check_cdf_status(iret)
      dims(1) = batchdim
      iret  = nf_def_var(ncid, 'batchsa', NCDOUBLE, 1, dims, batchsaid)
      call check_cdf_status(iret)
      dims(1) = ncdim
      dims(2) = batchdim
      iret  = nf_def_var(ncid, 'nesepm_av', NCDOUBLE, 2, dims, nesepm_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'tesepm_av', NCDOUBLE, 2, dims, tesepm_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'tisepm_av', NCDOUBLE, 2, dims, tisepm_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'posepm_av', NCDOUBLE, 2, dims, posepm_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'nesepi_av', NCDOUBLE, 2, dims, nesepi_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'tesepi_av', NCDOUBLE, 2, dims, tesepi_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'tisepi_av', NCDOUBLE, 2, dims, tisepi_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'posepi_av', NCDOUBLE, 2, dims, posepi_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'nesepa_av', NCDOUBLE, 2, dims, nesepa_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'tesepa_av', NCDOUBLE, 2, dims, tesepa_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'tisepa_av', NCDOUBLE, 2, dims, tisepa_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'posepa_av', NCDOUBLE, 2, dims, posepa_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'nemxip_av', NCDOUBLE, 2, dims, nemxip_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'temxip_av', NCDOUBLE, 2, dims, temxip_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'timxip_av', NCDOUBLE, 2, dims, timxip_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'pomxip_av', NCDOUBLE, 2, dims, pomxip_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'nemxap_av', NCDOUBLE, 2, dims, nemxap_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'temxap_av', NCDOUBLE, 2, dims, temxap_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'timxap_av', NCDOUBLE, 2, dims, timxap_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'pomxap_av', NCDOUBLE, 2, dims, pomxap_avid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'nesepm_std', NCDOUBLE, 2, dims, nesepm_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'tesepm_std', NCDOUBLE, 2, dims, tesepm_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'tisepm_std', NCDOUBLE, 2, dims, tisepm_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'posepm_std', NCDOUBLE, 2, dims, posepm_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'nesepi_std', NCDOUBLE, 2, dims, nesepi_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'tesepi_std', NCDOUBLE, 2, dims, tesepi_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'tisepi_std', NCDOUBLE, 2, dims, tisepi_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'posepi_std', NCDOUBLE, 2, dims, posepi_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'nesepa_std', NCDOUBLE, 2, dims, nesepa_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'tesepa_std', NCDOUBLE, 2, dims, tesepa_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'tisepa_std', NCDOUBLE, 2, dims, tisepa_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'posepa_std', NCDOUBLE, 2, dims, posepa_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'nemxip_std', NCDOUBLE, 2, dims, nemxip_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'temxip_std', NCDOUBLE, 2, dims, temxip_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'timxip_std', NCDOUBLE, 2, dims, timxip_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'pomxip_std', NCDOUBLE, 2, dims, pomxip_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'nemxap_std', NCDOUBLE, 2, dims, nemxap_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'temxap_std', NCDOUBLE, 2, dims, temxap_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'timxap_std', NCDOUBLE, 2, dims, timxap_stdid)
      call check_cdf_status(iret)
      iret  = nf_def_var(ncid, 'pomxap_std', NCDOUBLE, 2, dims, pomxap_stdid)
      call check_cdf_status(iret)
    endif
    ! assign attributes
    if (.not.batch_only) then
      iret = nf_put_att_text(ncid, timesaid, 'long_name', 4, 'time')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, timesaid, 'units', 2, 's ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fnixipid, 'long_name', 47, 'integrated poloidal particle flux, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fnixipid, 'units', 4, 's^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feexipid, 'long_name', 54, 'integrated poloidal electron energy flux, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feexipid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feixipid, 'long_name', 49, 'integrated poloidal ion energy flux, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feixipid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetxipid, 'long_name', 60, 'integrated poloidal total internal energy flux, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetxipid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchxipid, 'long_name', 41, 'integrated poloidal current, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchxipid, 'units', 2, 'A ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fnixapid, 'long_name', 47, 'integrated poloidal particle flux, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fnixapid, 'units', 4, 's^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feexapid, 'long_name', 54, 'integrated poloidal electron energy flux, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feexapid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feixapid, 'long_name', 49, 'integrated poloidal ion energy flux, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feixapid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetxapid, 'long_name', 60, 'integrated poloidal total internal energy flux, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetxapid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchxapid, 'long_name', 41, 'integrated poloidal current, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchxapid, 'units', 2, 'A ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepiid, 'long_name', 41, 'separatrix electron density, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepiid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepiid, 'long_name', 45, 'separatrix electron temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepiid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepiid, 'long_name', 40, 'separatrix ion temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepiid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tpsepiid, 'long_name', 42, 'separatrix plate temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tpsepiid, 'units', 2, 'K ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepiid, 'long_name', 34, 'separatrix potential, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepiid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepmid, 'long_name', 43, 'separatrix electron density, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepmid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepmid, 'long_name', 47, 'separatrix electron temperature, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepmid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepmid, 'long_name', 42, 'separatrix ion temperature, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepmid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepmid, 'long_name', 36, 'separatrix potential, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepmid, 'units', 2, 'V ')
      call check_cdf_status(iret)

      iret = nf_put_att_text(ncid, dnsepmid, 'long_name', 37, 'diffusion coefficient, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, dnsepmid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, dpsepmid, 'long_name', 46, 'pressure diffusion coefficient, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, dpsepmid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, kesepmid, 'long_name', 44, 'electron thermal diffusivity, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, kesepmid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, kisepmid, 'long_name', 39, 'ion thermal diffusivity, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, kisepmid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vxsepmid, 'long_name', 39, 'poloidal pinch velocity, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vxsepmid, 'units', 6, 'm.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vysepmid, 'long_name', 37, 'radial pinch velocity, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vysepmid, 'units', 6, 'm.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vssepmid, 'long_name', 37, 'viscosity coefficient, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vssepmid, 'units', 12, 'm.kg^-1.s^-1')
      call check_cdf_status(iret)

      iret = nf_put_att_text(ncid, nesepaid, 'long_name', 41, 'separatrix electron density, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepaid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepaid, 'long_name', 45, 'separatrix electron temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepaid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepaid, 'long_name', 40, 'separatrix ion temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepaid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tpsepaid, 'long_name', 42, 'separatrix plate temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tpsepaid, 'units', 2, 'K ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepaid, 'long_name', 34, 'separatrix potential, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepaid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nemxipid, 'long_name', 38, 'maximum electron density, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nemxipid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, temxipid, 'long_name', 42, 'maximum electron temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, temxipid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, timxipid, 'long_name', 37, 'maximum ion temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, timxipid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tpmxipid, 'long_name', 39, 'maximum plate temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tpmxipid, 'units', 2, 'K ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pomxipid, 'long_name', 31, 'maximum potential, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pomxipid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nemxapid, 'long_name', 38, 'maximum electron density, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nemxapid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, temxapid, 'long_name', 42, 'maximum electron temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, temxapid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, timxapid, 'long_name', 38, 'maximum ion temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, timxapid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tpmxapid, 'long_name', 40, 'maximum plate temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tpmxapid, 'units', 2, 'K ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pomxapid, 'long_name', 32, 'maximum potential, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pomxapid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fniyipid, 'long_name', 45, 'integrated radial particle flux, main chamber')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fniyipid, 'units', 4, 's^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feeyipid, 'long_name', 52, 'integrated radial electron energy flux, main chamber')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feeyipid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feiyipid, 'long_name', 47, 'integrated radial ion energy flux, main chamber')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feiyipid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetyipid, 'long_name', 58, 'integrated radial total internal energy flux, main chamber')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetyipid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchyipid, 'long_name', 39, 'integrated radial current, main chamber')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchyipid, 'units', 2, 'A ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fniyapid, 'long_name', 48, 'integrated radial particle flux, divertor region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fniyapid, 'units', 4, 's^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feeyapid, 'long_name', 55, 'integrated radial electron energy flux, divertor region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feeyapid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feiyapid, 'long_name', 50, 'integrated radial ion energy flux, divertor region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feiyapid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetyapid, 'long_name', 61, 'integrated radial total internal energy flux, divertor region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetyapid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchyapid, 'long_name', 42, 'integrated radial current, divertor region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchyapid, 'units', 2, 'A ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pwmxipid, 'long_name', 38, 'maximum total power flux, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pwmxipid, 'units', 6, 'W.m^-2')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pwmxapid, 'long_name', 38, 'maximum total power flux, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pwmxapid, 'units', 6, 'W.m^-2')
      call check_cdf_status(iret)

      ! global quantities
      iret = nf_put_att_text(ncid, tmneid, 'long_name', 25, 'total number of particles')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tmneid, 'units', 2, '  ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tmteid, 'long_name', 21, 'total electron energy')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tmteid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tmtiid, 'long_name', 16, 'total ion energy')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tmtiid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tmhacoreid, 'long_name', 24, 'H-alpha emissivity, core')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tmhacoreid, 'units', 19, 'photons.m^-2.sr^-1?')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tmhasolid, 'long_name', 41, 'H-alpha emissivity, SOL above the x-point')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tmhasolid, 'units', 19, 'photons.m^-2.sr^-1?')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tmhadivid, 'long_name', 35, 'H-alpha emissivity, divertor region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tmhadivid, 'units', 19, 'photons.m^-2.sr^-1?')
      call check_cdf_status(iret)

      ! internal flux quantities
      iret = nf_put_att_text(ncid, fnisipid, 'long_name', 54, 'poloidal particle flux, into Western separatrix throat')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fnisipid, 'units', 4, 's^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feesipid, 'long_name', 61, 'poloidal electron energy flux, into Western separatrix throat')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feesipid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feisipid, 'long_name', 56, 'poloidal ion energy flux, into Western separatrix throat')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feisipid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetsipid, 'long_name', 67, 'poloidal total internal energy flux, into Western separatrix throat')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetsipid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchsipid, 'long_name', 48, 'poloidal current, into Western separatrix throat')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchsipid, 'units', 2, 'A ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fnisapid, 'long_name', 54, 'poloidal particle flux, into Eastern separatrix throat')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fnisapid, 'units', 4, 's^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feesapid, 'long_name', 61, 'poloidal electron energy flux, into Eastern separatrix throat')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feesapid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feisapid, 'long_name', 56, 'poloidal ion energy flux, into Eastern separatrix throat')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feisapid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetsapid, 'long_name', 67, 'poloidal total internal energy flux, into Eastern separatrix throat')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetsapid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchsapid, 'long_name', 48, 'poloidal current, into Eastern separatrix throat')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchsapid, 'units', 2, 'A ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fnisippid, 'long_name', 40, 'poloidal particle flux, core x-pt region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fnisippid, 'units', 4, 's^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feesippid, 'long_name', 52, 'poloidal electron energy flux, core x-pt flux region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feesippid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feisippid, 'long_name', 47, 'poloidal ion energy flux, core x-pt flux region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feisippid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetsippid, 'long_name', 58, 'poloidal total internal energy flux, core x-pt flux region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetsippid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchsippid, 'long_name', 39, 'poloidal current, core x-pt flux region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchsippid, 'units', 2, 'A ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fnisappid, 'long_name', 48, 'poloidal particle flux, x-pt private flux region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fnisappid, 'units', 4, 's^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feesappid, 'long_name', 55, 'poloidal electron energy flux, x-pt private flux region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feesappid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feisappid, 'long_name', 50, 'poloidal ion energy flux, x-pt private flux region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, feisappid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetsappid, 'long_name', 61, 'poloidal total internal energy flux, x-pt private flux region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fetsappid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchsappid, 'long_name', 42, 'poloidal current, x-pt private flux region')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fchsappid, 'units', 2, 'A ')
      call check_cdf_status(iret)

      ! Western edge (inboard divertor for LSN, outboard divertor for USN) quantities
      iret = nf_put_att_text(ncid, ne3dlid, 'long_name', 30, 'electron density, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ne3dlid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, te3dlid, 'long_name', 34, 'electron temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, te3dlid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, te3dlid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ti3dlid, 'long_name', 29, 'ion temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ti3dlid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, ti3dlid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tp3dlid, 'long_name', 31, 'plate temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tp3dlid, 'units', 2, 'K ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, po3dlid, 'long_name', 23, 'potential, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, po3dlid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, an3dlid, 'long_name', 26, 'atom density, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, an3dlid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, mn3dlid, 'long_name', 30, 'molecule density, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, mn3dlid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      dvals(1) = -1.0_R8
      iret = nf_put_att_text(ncid, fn3dlid, 'long_name', 40, 'poloidal main species flux, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fn3dlid, 'units', 4, 's^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, fn3dlid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fl3dlid, 'long_name', 36, 'poloidal electron flux, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fl3dlid, 'units', 4, 's^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, fl3dlid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fo3dlid, 'long_name', 31, 'poloidal ion flux, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fo3dlid, 'units', 4, 's^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, fo3dlid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fe3dlid, 'long_name', 43, 'poloidal electron energy flux, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fe3dlid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, fe3dlid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fi3dlid, 'long_name', 38, 'poloidal ion energy flux, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fi3dlid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, fi3dlid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ft3dlid, 'long_name', 40, 'poloidal total energy flux, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ft3dlid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, ft3dlid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fc3dlid, 'long_name', 31, 'poloidal current, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fc3dlid, 'units', 2, 'A ')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, fc3dlid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)

      ! inboard midplane quantities
      iret = nf_put_att_text(ncid, ne3diid, 'long_name', 34, 'electron density, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ne3diid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, te3diid, 'long_name', 38, 'electron temperature, inboard midplane')
      call check_cdf_status(iret)
      dvals(1) = 1.0_R8/ev
      iret = nf_put_att_text(ncid, te3diid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, te3diid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ti3diid, 'long_name', 33, 'ion temperature, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ti3diid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, ti3diid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, po3diid, 'long_name', 27, 'potential, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, po3diid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, an3diid, 'long_name', 30, 'atom density, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, an3diid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, mn3diid, 'long_name', 34, 'molecule density, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, mn3diid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, dn3diid, 'long_name', 39, 'diffusion coefficient, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, dn3diid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, dp3diid, 'long_name', 48, 'pressure diffusion coefficient, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, dp3diid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, lh3diid, 'long_name', 50, 'radial neutral heat flux limiter, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, lh3diid, 'units', 1, ' ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ln3diid, 'long_name', 50, 'radial neutral part flux limiter, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ln3diid, 'units', 1, ' ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ke3diid, 'long_name', 46, 'electron thermal diffusivity, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ke3diid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ki3diid, 'long_name', 41, 'ion thermal diffusivity, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ki3diid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vx3diid, 'long_name', 41, 'poloidal pinch velocity, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vx3diid, 'units', 6, 'm.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vy3diid, 'long_name', 39, 'radial pinch velocity, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vy3diid, 'units', 6, 'm.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vs3diid, 'long_name', 39, 'viscosity coefficient, inboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vs3diid, 'units', 12, 'm.kg^-1.s^-1')
      call check_cdf_status(iret)

      ! upper inboard divertor quantities
      if(nytl.gt.0) then
        iret = nf_put_att_text(ncid, ne3dtlid, 'long_name', 40, 'electron density, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ne3dtlid, 'units', 4, 'm^-3')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, te3dtlid, 'long_name', 44, 'electron temperature, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, te3dtlid, 'units', 2, 'eV')
        call check_cdf_status(iret)
        iret = nf_put_att_double(ncid, te3dtlid, 'scale', NCDOUBLE, 1, dvals(1))
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ti3dtlid, 'long_name', 39, 'ion temperature, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ti3dtlid, 'units', 2, 'eV')
        call check_cdf_status(iret)
        iret = nf_put_att_double(ncid, ti3dtlid, 'scale', NCDOUBLE, 1, dvals(1))
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, tp3dtlid, 'long_name', 41, 'plate temperature, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, tp3dtlid, 'units', 2, 'K ')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, po3dtlid, 'long_name', 32, 'potential, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, po3dtlid, 'units', 2, 'V ')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, an3dtlid, 'long_name', 36, 'atom density, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, an3dtlid, 'units', 4, 'm^-3')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, mn3dtlid, 'long_name', 40, 'molecule density, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, mn3dtlid, 'units', 4, 'm^-3')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fn3dtlid, 'long_name', 50, 'poloidal main species flux, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fn3dtlid, 'units', 4, 's^-1')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fl3dtlid, 'long_name', 46, 'poloidal electron flux, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fl3dtlid, 'units', 4, 's^-1')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fo3dtlid, 'long_name', 41, 'poloidal ion flux, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fo3dtlid, 'units', 4, 's^-1')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fe3dtlid, 'long_name', 53, 'poloidal electron energy flux, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fe3dtlid, 'units', 2, 'W ')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fi3dtlid, 'long_name', 48, 'poloidal ion energy flux, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fi3dtlid, 'units', 2, 'W ')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ft3dtlid, 'long_name', 50, 'poloidal total energy flux, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ft3dtlid, 'units', 2, 'W ')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fc3dtlid, 'long_name', 41, 'poloidal current, upper inboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fc3dtlid, 'units', 2, 'A ')
        call check_cdf_status(iret)
      endif

      ! outboard midplane quantities
      iret = nf_put_att_text(ncid, ne3daid, 'long_name', 35, 'electron density, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ne3daid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, te3daid, 'long_name', 39, 'electron temperature, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, te3daid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, te3daid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ti3daid, 'long_name', 34, 'ion temperature, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ti3daid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, ti3daid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, po3daid, 'long_name', 28, 'potential, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, po3daid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, an3daid, 'long_name', 31, 'atom density, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, an3daid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, mn3daid, 'long_name', 35, 'molecule density, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, mn3daid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, dn3daid, 'long_name', 40, 'diffusion coefficient, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, dn3daid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, dp3daid, 'long_name', 49, 'pressure diffusion coefficient, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, dp3daid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, lh3daid, 'long_name', 51, 'radial neutral heat flux limiter, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, lh3daid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ln3daid, 'long_name', 51, 'radial neutral part flux limiter, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ln3daid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ke3daid, 'long_name', 47, 'electron thermal diffusivity, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ke3daid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ki3daid, 'long_name', 42, 'ion thermal diffusivity, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ki3daid, 'units', 8, 'm^2.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vx3daid, 'long_name', 42, 'poloidal pinch velocity, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vx3daid, 'units', 6, 'm.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vy3daid, 'long_name', 40, 'radial pinch velocity, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vy3daid, 'units', 6, 'm.s^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vs3daid, 'long_name', 40, 'viscosity coefficient, outboard midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, vs3daid, 'units', 12, 'm.kg^-1.s^-1')
      call check_cdf_status(iret)

      ! Eastern edge (outboard divertor for LSN, inboard divertor for USN) quantities
      iret = nf_put_att_text(ncid, ne3drid, 'long_name', 30, 'electron density, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ne3drid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, te3drid, 'long_name', 34, 'electron temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, te3drid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, te3drid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ti3drid, 'long_name', 29, 'ion temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ti3drid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_double(ncid, ti3drid, 'scale', NCDOUBLE, 1, dvals(1))
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tp3drid, 'long_name', 31, 'plate temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tp3drid, 'units', 2, 'K ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, po3drid, 'long_name', 23, 'potential, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, po3drid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, an3drid, 'long_name', 26, 'atom density, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, an3drid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, mn3drid, 'long_name', 30, 'molecule density, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, mn3drid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fn3drid, 'long_name', 40, 'poloidal main species flux, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fn3drid, 'units', 4, 's^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fl3drid, 'long_name', 36, 'poloidal electron flux, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fl3drid, 'units', 4, 's^-1')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fo3drid, 'long_name', 31, 'poloidal ion flux, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fo3drid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fe3drid, 'long_name', 43, 'poloidal electron energy flux, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fe3drid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fi3drid, 'long_name', 38, 'poloidal ion energy flux, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fi3drid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ft3drid, 'long_name', 40, 'poloidal total energy flux, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, ft3drid, 'units', 2, 'W ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fc3drid, 'long_name', 30, 'poloidal current, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, fc3drid, 'units', 2, 'A ')
      call check_cdf_status(iret)

      ! upper outboard divertor quantities
      if(nytr.gt.0) then
        iret = nf_put_att_text(ncid, ne3dtrid, 'long_name', 41, 'electron density, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ne3dtrid, 'units', 4, 'm^-3')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, te3dtrid, 'long_name', 45, 'electron temperature, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, te3dtrid, 'units', 2, 'eV')
        call check_cdf_status(iret)
        iret = nf_put_att_double(ncid, te3dtrid, 'scale', NCDOUBLE, 1, dvals(1))
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ti3dtrid, 'long_name', 40, 'ion temperature, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ti3dtrid, 'units', 2, 'eV')
        call check_cdf_status(iret)
        iret = nf_put_att_double(ncid, ti3dtrid, 'scale', NCDOUBLE, 1, dvals(1))
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, tp3dtrid, 'long_name', 42, 'plate temperature, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, tp3dtrid, 'units', 2, 'K ')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, po3dtrid, 'long_name', 34, 'potential, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, po3dtrid, 'units', 2, 'V ')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, an3dtrid, 'long_name', 37, 'atom density, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, an3dtrid, 'units', 4, 'm^-3')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, mn3dtrid, 'long_name', 41, 'molecule density, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, mn3dtrid, 'units', 4, 'm^-3')
        call check_cdf_status(iret)
        dvals(1) = -1.0_R8
        iret = nf_put_att_text(ncid, fn3dtrid, 'long_name', 51, 'poloidal main species flux, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fn3dtrid, 'units', 4, 's^-1')
        call check_cdf_status(iret)
        iret = nf_put_att_double(ncid, fn3dtrid, 'scale', NCDOUBLE, 1, dvals(1))
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fl3dtrid, 'long_name', 47, 'poloidal electron flux, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fl3dtrid, 'units', 4, 's^-1')
        call check_cdf_status(iret)
        iret = nf_put_att_double(ncid, fl3dtrid, 'scale', NCDOUBLE, 1, dvals(1))
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fo3dtrid, 'long_name', 42, 'poloidal ion flux, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fo3dtrid, 'units', 2, 'W ')
        call check_cdf_status(iret)
        iret = nf_put_att_double(ncid, fo3dtrid, 'scale', NCDOUBLE, 1, dvals(1))
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fe3dtrid, 'long_name', 54, 'poloidal electron energy flux, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fe3dtrid, 'units', 2, 'W ')
        call check_cdf_status(iret)
        iret = nf_put_att_double(ncid, fe3dtrid, 'scale', NCDOUBLE, 1, dvals(1))
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fi3dtrid, 'long_name', 49, 'poloidal ion energy flux, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fi3dtrid, 'units', 2, 'W ')
        call check_cdf_status(iret)
        iret = nf_put_att_double(ncid, fi3dtrid, 'scale', NCDOUBLE, 1, dvals(1))
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ft3dtrid, 'long_name', 51, 'poloidal total energy flux, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, ft3dtrid, 'units', 2, 'W ')
        call check_cdf_status(iret)
        iret = nf_put_att_double(ncid, ft3dtrid, 'scale', NCDOUBLE, 1, dvals(1))
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fc3dtrid, 'long_name', 42, 'poloidal current, upper outboard divertor')
        call check_cdf_status(iret)
        iret = nf_put_att_text(ncid, fc3dtrid, 'units', 2, 'A ')
        call check_cdf_status(iret)
        iret = nf_put_att_double(ncid, fc3dtrid, 'scale', NCDOUBLE, 1, dvals(1))
        call check_cdf_status(iret)
      endif
    else
    !wdk averaged quantities
      iret = nf_put_att_text(ncid, nesepm_avid, 'long_name', 52, 'averaged separatrix electron density, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepm_avid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepm_avid, 'long_name', 56, 'averaged separatrix electron temperature, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepm_avid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepm_avid, 'long_name', 51, 'averaged separatrix ion temperature, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepm_avid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepm_avid, 'long_name', 45, 'averaged separatrix potential, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepm_avid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepi_avid, 'long_name', 50, 'averaged separatrix electron density, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepi_avid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepi_avid, 'long_name', 54, 'averaged separatrix electron temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepi_avid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepi_avid, 'long_name', 49, 'averaged separatrix ion temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepi_avid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepi_avid, 'long_name', 43, 'averaged separatrix potential, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepi_avid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepa_avid, 'long_name', 50, 'averaged separatrix electron density, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepa_avid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepa_avid, 'long_name', 54, 'averaged separatrix electron temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepa_avid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepa_avid, 'long_name', 49, 'averaged separatrix ion temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepa_avid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepa_avid, 'long_name', 43, 'averaged separatrix potential, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepa_avid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nemxip_avid, 'long_name', 47, 'averaged maximum electron density, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nemxip_avid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, temxip_avid, 'long_name', 51, 'averaged maximum electron temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, temxip_avid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, timxip_avid, 'long_name', 46, 'averaged maximum ion temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, timxip_avid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pomxip_avid, 'long_name', 40, 'averaged maximum potential, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pomxip_avid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nemxap_avid, 'long_name', 47, 'averaged maximum electron density, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nemxap_avid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, temxap_avid, 'long_name', 51, 'averaged maximum electron temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, temxap_avid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, timxap_avid, 'long_name', 46, 'averaged maximum ion temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, timxap_avid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pomxap_avid, 'long_name', 40, 'averaged maximum potential, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pomxap_avid, 'units', 2, 'V ')
      call check_cdf_status(iret)

    !wdk standard deviation of averaged quantities
      iret = nf_put_att_text(ncid, nesepm_stdid, 'long_name', 55, 'variance of separatrix electron density, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepm_stdid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepm_stdid, 'long_name', 59, 'variance of separatrix electron temperature, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepm_stdid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepm_stdid, 'long_name', 54, 'variance of separatrix ion temperature, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepm_stdid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepm_stdid, 'long_name', 48, 'variance of separatrix potential, outer midplane')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepm_stdid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepi_stdid, 'long_name', 53, 'variance of separatrix electron density, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepi_stdid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepi_stdid, 'long_name', 57, 'variance of separatrix electron temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepi_stdid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepi_stdid, 'long_name', 52, 'variance of separatrix ion temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepi_stdid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepi_stdid, 'long_name', 46, 'variance of separatrix potential, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepi_stdid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepa_stdid, 'long_name', 53, 'variance of separatrix electron density, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nesepa_stdid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepa_stdid, 'long_name', 57, 'variance of separatrix electron temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tesepa_stdid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepa_stdid, 'long_name', 52, 'variance of separatrix ion temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, tisepa_stdid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepa_stdid, 'long_name', 46, 'variance of separatrix potential, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, posepa_stdid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nemxip_stdid, 'long_name', 50, 'variance of maximum electron density, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nemxip_stdid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, temxip_stdid, 'long_name', 54, 'variance of maximum electron temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, temxip_stdid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, timxip_stdid, 'long_name', 49, 'variance of maximum ion temperature, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, timxip_stdid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pomxip_stdid, 'long_name', 43, 'variance of maximum potential, Western edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pomxip_stdid, 'units', 2, 'V ')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nemxap_stdid, 'long_name', 50, 'variance of maximum electron density, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, nemxap_stdid, 'units', 4, 'm^-3')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, temxap_stdid, 'long_name', 54, 'variance of maximum electron temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, temxap_stdid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, timxap_stdid, 'long_name', 50, 'variance of maximum ion temperature, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, timxap_stdid, 'units', 2, 'eV')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pomxap_stdid, 'long_name', 43, 'variance of maximum potential, Eastern edge')
      call check_cdf_status(iret)
      iret = nf_put_att_text(ncid, pomxap_stdid, 'units', 2, 'V ')
      call check_cdf_status(iret)
    endif

    ! leave define mode
    iret = nf_enddef(ncid)
    call check_cdf_status(iret)
    iret = nf_close(ncid)
    call check_cdf_status(iret)
    return
  end subroutine b2crtimecdf

  subroutine rwcdf(rw,ncid,data_name,imap,data_set,iret)
#     include <netcdf.inc>

    character*(*) rw,data_name
    integer ncid,imap(*),iret,i,varid,dimlen
    real(kind=R8), Intent(InOut) :: data_set(*)
    character*(maxncnam) dimnam
    integer vartyp,nvdims,start(maxvdims),mycount(maxvdims),dimids(maxvdims)
    character*(*) timnam,batchnam
    character*(maxncnam) timsav,batchsav
    integer ntsav,ntstep,nasav,nastep
    integer :: istride, imax
    logical, parameter :: debug = .false.
    save timsav,ntsav,batchsav,nasav
    data timsav /'!!!! INVALID NAME !!!!'/
    external subini, subend, xerrab
    !
    call subini ('rwcdf')
    iret = nf_inq_varid(ncid,data_name,varid)
    call check_cdf_status(iret)
    if(iret.ne.0) then
      write(*,*) 'Error: Could not inquire data_name: ',trim(data_name)
      call xerrab('Data name '//trim(data_name)//' not declared')
    endif
    if (debug) write(*,*) "Subroutine rwcdf in mode: ", rw
    if (debug) write(*,*) "Working on variable: ", data_name
    iret = nf_inq_varndims(ncid,varid,nvdims)
    if (debug) write(*,*) "Variable has nvdims=", nvdims
    if (debug) write(*,*) "Input imap(:)=", imap(1:nvdims)
    iret = nf_inq_vardimid(ncid,varid,dimids)
    call check_cdf_status(iret)
    mycount(1) = 1 ! for scalars
    start(1) = 1
    do i=1,nvdims
      iret = nf_inq_dim(ncid,dimids(i),dimnam,dimlen)
      call check_cdf_status(iret)
      mycount(i) = dimlen
      if(dimnam.eq.timsav) then
        start(i)=ntsav
        mycount(i)=1
      elseif(dimnam.eq.batchsav) then
        start(i)=nasav
        mycount(i)=1
      else
        start(i)=1
      endif
    enddo
    if (debug) write(*,*) "start(:)=", start(1:nvdims)
    if (debug) write(*,*) "mycount(:)=", mycount(1:nvdims)
    istride = 1
    imax = 1
    do i=1,nvdims-1
      istride = istride*imap(i)
      imax = imax*imap(i)*mycount(i)
    enddo
    if (debug) write(*,*) "istride=", istride
    if (debug) write(*,*) "imax=", imax
    iret = nf_inq_vartype(ncid,varid,vartyp)
    call check_cdf_status(iret)
    if(rw.eq.'read') then
      select case (vartyp)
      case (NCDOUBLE)
        iret = nf_get_vara_double(ncid,varid,start,mycount,data_set(1:imax:istride))
        call check_cdf_status(iret)
      case default
        call xerrab ('Unknown data type in rwcdf read')
      end select
    elseif(rw.eq.'write') then
      select case (vartyp)
      case (NCDOUBLE)
        iret = nf_put_vara_double(ncid,varid,start,mycount,data_set(1:imax:istride))
        call check_cdf_status(iret)
      case default
        call xerrab ('Unknown data type in rwcdf write')
      end select
    else
      write(*,*) 'Either "read" or "write" must be chosen, not ', rw
    endif

    call subend ()
    return
    !
    entry rwcdf_settime(timnam,ntstep)
    write(*,*) 'Saving ',trim(timnam),' as the time dimension'
    write(*,*) 'ntstep = ',ntstep
    timsav=timnam
    ntsav=ntstep
    return
    !
    entry rwcdf_setbatch(batchnam,nastep)
    write(*,*) 'Saving ',trim(batchnam),' as the batch dimension'
    write(*,*) 'nastep = ',nastep
    batchsav=batchnam
    nasav=nastep
    return
  end subroutine rwcdf
#endif
!
  subroutine output_ds(crx,cry,nx,ny,iref,target_offset, &
       jsep,iystart,iyend,filename)
    use b2mod_indirect
    implicit none
    integer nx,ny,iref,jsep,iystart,iyend,target_offset
    real (kind=R8) :: &
         crx(-1:nx,-1:ny,0:3),cry(-1:nx,-1:ny,0:3)
    real (kind=R8) :: &
         ds(-1:ny), ds_offset
    character*(*) filename
    integer ix,iy
    external subini, subend, xertst
    intrinsic sqrt
    real (kind=R8) :: &
         cr,cz

    cr(ix,iy)=0.25_R8*(crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+crx(ix,iy,3))
    cz(ix,iy)=0.25_R8*(cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+cry(ix,iy,3))

    call subini ('output_ds')
    iystart=-1
    do while (region(iref,iystart,0).eq.0 .and. iystart.lt.ny)
      iystart=iystart+1
    enddo
    call xertst(iystart.le.ny, 'faulty parameter iystart')
    iyend=ny
    do while (region(iref,iyend,0).eq.0 .and. iyend.gt.-1)
      iyend=iyend-1
    enddo
    call xertst(iyend.ge.-1, 'faulty parameter iyend')
    if(iystart.eq.ny.and.iyend.eq.-1) then
      ! special case [DPC]
      iystart=-1
      iyend=ny
      write(*,*) 'special treatment for iystart, iyend for ',trim(filename)
    endif
    call xertst(iyend.ge.iystart, 'faulty parameter iystart & iyend')
    ds(iystart)= &
         sqrt((cr(iref+target_offset,iystart)-0.5_R8*(crx(iref+target_offset,iystart,0)+crx(iref+target_offset,iystart,1)))**2+ &
              (cz(iref+target_offset,iystart)-0.5_R8*(cry(iref+target_offset,iystart,0)+cry(iref+target_offset,iystart,1)))**2)
    do iy=iystart+1,iyend
      ds(iy)=ds(iy-1)+ &
           sqrt((cr(iref+target_offset,iy)-cr(iref+target_offset,iy-1))**2+ &
                (cz(iref+target_offset,iy)-cz(iref+target_offset,iy-1))**2)
    enddo
    if(iystart.le.jsep.and.iyend.gt.jsep) then
      ds_offset=(ds(jsep)+ds(jsep+1))/2.0_R8
      do iy=iystart,iyend
        ds(iy)=ds(iy)-ds_offset
      enddo
    endif
    open(99,file=filename)
    do iy=iystart,iyend
      write(99,*) ds(iy)
    enddo
    close(99)
    call subend ()
    return
  end subroutine output_ds


  subroutine calc_fet(ix,iy,side,fac_flux,nx,ny,ns,ismain,BoRiS,fet,fni0,fee0,fei0,fch0,pwr)
    use b2mod_plasma   , only : ti, te, fna, fhe, fhi, fch, fht, fhj
    use b2mod_indirect , only : rightix, rightiy, bottomix, bottomiy, topix, topiy, leftix, leftiy
    use b2mod_external , only : fhi_ext, pt_ext, ta_ext, ua_ext, am_ext, ns_ext, fa_ext
    use b2mod_constants , only : ev, mp
    use b2mod_geo , only : hx, hy, qz, gs
    implicit none
    integer, intent(in) :: ix, iy
    integer, intent(in) :: ismain, nx, ny, ns
    real(kind=R8), intent(in) :: BoRiS, fac_flux
    real(kind=R8), intent(out) :: fet
    real(kind=R8), intent(out), Optional :: fni0, fee0, fei0, fch0, pwr
    character(len=1) :: side
    ! Local vars
    integer :: ix_adj, iy_adj, is, ix_flux, iy_flux, idir
    real(kind=R8) :: kintmp, rpttmp, taf
    real(kind=R8) :: h(-1:nx,-1:ny)
    ! computation

    select case (side)
    case ('l','L')
      ix_flux = rightix(ix,iy) ! Index to cell with flux entering cell
      iy_flux = rightiy(ix,iy)
      ix_adj  = rightix(ix,iy) ! Index to cell adjacent
      iy_adj  = rightiy(ix,iy)
      idir = 0                 ! Index in flux variables (x vs y direction)
      h(-1:nx,-1:ny) = hx(-1:nx,-1:ny)
    case ('r','R')
      ix_flux = ix
      iy_flux = iy
      ix_adj = leftix(ix,iy)
      iy_adj = leftiy(ix,iy)
      idir = 0
      h(-1:nx,-1:ny) = hx(-1:nx,-1:ny)
    case ('t','T')
      ix_flux = ix
      iy_flux = iy
      ix_adj = bottomix(ix,iy)
      iy_adj = bottomiy(ix,iy)
      idir = 1
      h(-1:nx,-1:ny) = hy(-1:nx,-1:ny)*qz(-1:nx,-1:ny,1)
    case ('b','B')
      ix_flux = topix(ix,iy)
      iy_flux = topiy(ix,iy)
      ix_adj = topix(ix,iy)
      iy_adj = topiy(ix,iy)
      idir = 1
      h(-1:nx,-1:ny) = hy(-1:nx,-1:ny)*qz(-1:nx,-1:ny,1)
    case default
      call xerrab('Unknown side in calc_fet')
    end select
    if (present(fni0)) fni0 = fac_flux*fna(ix_flux,iy_flux,idir,ismain)
    if (present(fee0)) fee0 = fac_flux*fhe(ix_flux,iy_flux,idir)
    if (present(fei0)) fei0 = fac_flux*fhi(ix_flux,iy_flux,idir)
    if (present(fch0)) fch0 = fac_flux*fch(ix_flux,iy_flux,idir)
    fet = fac_flux*(fht(ix_flux,iy_flux,idir) - fhj(ix_flux,iy_flux,idir) + fhi_ext(ix_flux,iy_flux,idir))
    do is=0,ns_ext-1
      kintmp = 0.5_R8*am_ext(is)*mp*(ua_ext(ix,iy,is)**2 * h(ix_adj,iy_adj)+ &
           ua_ext(ix_adj,iy_adj,is)**2*h(ix,iy))/(h(ix_adj,iy_adj)+h(ix,iy))
      rpttmp = (pt_ext(ix,iy,is)*h(ix_adj,iy_adj)+pt_ext(ix_adj,iy_adj,is)*h(ix,iy))/(h(ix_adj,iy_adj)+h(ix,iy))
      taf = (ta_ext(ix_adj,iy_adj,is)*h(ix,iy)+ta_ext(ix,iy,is)*h(ix_adj,iy_adj))/(h(ix,iy)+h(ix_adj,iy_adj))
      fet = fet + fac_flux*(rpttmp*ev + (kintmp+taf)*(1.0_R8-BoRiS))*fa_ext(ix_flux,iy_flux,idir,is)
    enddo
    if (present(pwr)) pwr = Abs(fet)/gs(ix_flux,iy_flux,idir)
  end subroutine calc_fet

end module b2mod_mwti

!!!Local Variables:
!!! mode: f90
!!! End:
