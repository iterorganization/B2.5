module b2mod_mwti
  use b2mod_types , only : R8
  use b2mod_cdf
  use b2mod_subsys
  implicit none
  private
  public :: output_ds_cv, output_ds_fc
#ifndef SOLPS4_3
  public :: b2mwti, dealloc_b2mod_mwti
#endif
  real (kind=R8), allocatable, save, public :: &
         nesepi_av(:), tesepi_av(:), tisepi_av(:), &
         nesepm_av(:), tesepm_av(:), tisepm_av(:), &
         nesepa_av(:), tesepa_av(:), tisepa_av(:), &
         posepi_av(:), posepm_av(:), posepa_av(:), &
         ktsepi_av(:), ktsepm_av(:), ktsepa_av(:)
  real (kind=R8), allocatable, save, public :: &
         nemxip_av(:), temxip_av(:), timxip_av(:), &
         nemxap_av(:), temxap_av(:), timxap_av(:), &
         pomxip_av(:), pomxap_av(:)
  real (kind=R8), allocatable, save, public :: &
         nesepi_std(:), tesepi_std(:), tisepi_std(:), &
         nesepm_std(:), tesepm_std(:), tisepm_std(:), &
         nesepa_std(:), tesepa_std(:), tisepa_std(:), &
         posepi_std(:), posepm_std(:), posepa_std(:), &
         ktsepi_std(:), ktsepm_std(:), ktsepa_std(:)
  real (kind=R8), allocatable, save, public :: &
         nemxip_std(:), temxip_std(:), timxip_std(:), &
         nemxap_std(:), temxap_std(:), timxap_std(:), &
         pomxip_std(:), pomxap_std(:)
contains

#ifndef SOLPS4_3
  subroutine b2mwti (itim, tim, ntim, b2time, ntim_batch, &
                     nCv, ns, nncutmax, geo, mpg, switch, &
                     pl, dv, co, rt, srw, &
                     ismain, ismain0, lwti, lwav, luav)
!    use b2mod_geo
!    use b2mod_plasma
!    use b2mod_rates
!    use b2mod_residuals
!    use b2mod_sources
!    use b2mod_transport
!    use b2mod_anomalous_transport
    use b2mod_neutrals_namelist
!    use b2mod_work
!    use b2mod_indirect
    use b2mod_constants
!    use b2mod_tallies
!    use b2mod_wall
    use b2mod_b2cmpa
!    use b2mod_external
    use b2us_geo
    use b2us_map
    use b2us_plasma
    use b2mod_geometry &
    , only : geometryID
    use b2mod_user_namelist &
    , only : omp, imp, nimp, nomp, icsepimp
#ifndef NO_CDF
    use b2mod_geometry &
    , only : GEOMETRY_CDN, GEOMETRY_DDN_TOP, GEOMETRY_DDN_BOTTOM, &
             GEOMETRY_LFS_SNOWFLAKE_PLUS, GEOMETRY_LFS_SNOWFLAKE_MINUS
#endif
    use b2mod_user_namelist &
    , only : icsepomp
    use b2mod_switches
#ifndef SOLPS4_3
#ifdef B25_EIRENE
    use eirmod_extrab25
#endif
#endif
    implicit none
    !   ..input arguments (unchanged on exit)
    type (geometry), intent (in) :: geo
    type (mapping), intent (in) :: mpg
    type (switches), intent (in) :: switch
    type (B2Plasma), intent (inout) :: pl
    type (B2Derivatives), intent (inout) :: dv
    type (B2SourceWork), intent (inout) :: srw
    type (B2Coeff), intent (in) :: co
    type (B2Rates), intent (in) :: rt
    integer, intent(in) :: itim, ntim, b2time, ntim_batch, &
                           nCv, ns, nncutmax, ismain, ismain0
    real (kind=R8), intent(in) :: tim
    logical, intent(in) :: lwti, lwav, luav
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
    integer ncall
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
         pomxip(nncutmax), pomxap(nncutmax), &
         tpmxip(nncutmax), tpmxap(nncutmax)
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

    integer icv
    integer target_offset
    integer nc
#ifdef WG_TODO
    integer jxi, jxa, jsep
    integer ix, iy, ic, ixtl, ixtr, ix_off
    integer iyastrt, iyistrt, iylstrt, iyrstrt, iytlstrt, iytrstrt, &
         iyaend, iyiend, iylend, iyrend, iytlend, iytrend, &
         nya, nyi, nybl, nybr, nytl, nytr
#endif

    !   ..procedures
    external xertst, ipgeti, batch_average
#ifdef WG_TODO
    real(kind=R8) :: fnitmp, feetmp, feitmp, fchtmp, fettmp, pwrtmp
#endif
    integer, save :: write_2d = 0
    integer, save :: ntstep, nastep, geometryType
#ifndef NO_CDF
    integer, save :: ncid, nbatch
    integer imap(maxvdims), iret, ift, iatm, cvtrg
    integer nvars, natts, ndims, unlimid
    real (kind=R8) :: fac
    real (kind=R8) :: &
         nesepi(nncutmax), tesepi(nncutmax), tisepi(nncutmax), &
         nesepm(nncutmax), tesepm(nncutmax), tisepm(nncutmax), &
         nesepa(nncutmax), tesepa(nncutmax), tisepa(nncutmax), &
         posepi(nncutmax), posepm(nncutmax), posepa(nncutmax), &
         dnsepm(nncutmax), dpsepm(nncutmax), kesepm(nncutmax), &
         kisepm(nncutmax), vxsepm(nncutmax), vysepm(nncutmax), &
         vssepm(nncutmax), tpsepi(nncutmax), tpsepa(nncutmax), &
         ktsepm(nncutmax), ktsepi(nncutmax), ktsepa(nncutmax)
    real (kind=R8) :: &
         tmhacore(1), tmhasol(1), tmhadiv(1)
    real (kind=R8) :: &
         timesa(1), batchsa(1), tstepn(1)
    real (kind=R8), allocatable :: slice(:)
    logical ex
    character*5 rw
    character*256, save :: filename, filename_av
    real(kind=R8) :: rratio
    external rratio
#endif
    !   ..initialisation
#ifdef WG_TODO
    save jxi, jxa, jsep, ixtl, ixtr, &
         iylstrt, iyrstrt, iytlstrt, iytrstrt, &
         iylend,  iyrend,  iytlend,  iytrend, &
         nc, nya, nyi, nybl, nybr, nytl, nytr
#endif
    save ncall, target_offset
    data ncall/0/, target_offset/1/

    !-----------------------------------------------------------------------
    !.computation

    ! ..preliminaries
    !   ..subprogram start-up calls
    call subini ('b2mwti')
    !     ..test nx, ny
    call xertst (0.le.nCv,'faulty argument nCv')
    call xertst (1.le.ns, 'faulty argument ns')
    call xertst (0.le.ismain.and.ismain.lt.ns, &
         'invalid main plasma species index ismain')
    !   ..extensive tests on first few calls
    if (ncall.eq.0) then
      geometryType = geometryId ( mpg, geo )
      !   ..test state
      call ipgeti ('b2mwti_2dwrite',write_2d)
      call xertst (0.le.write_2d.and.write_2d.le.2,'faulty internal parameter write_2d')
      call ipgeti ('b2mwti_target_offset',target_offset)
      call xertst (0.le.target_offset.and.target_offset.le.1,'faulty internal parameter target_offset')
      write(*,*) 'target_offset ', target_offset
      call xertst(icsepomp.gt.0,'Invalid icsepomp value, check rzomp in b2.user.parameters')
      nc = max(mpg%nXpt,1)
      if (nimp.gt.0) call output_ds_cv(geo,nimp,imp,icsepimp-1,'dsi')
      if (nomp.gt.0) call output_ds_cv(geo,nomp,omp,icsepomp-1,'dsa')
#ifdef WG_TODO
      call get_jsep(nx,ny,jxi,jxa,jsep)
      call output_ds(ny, -1,+target_offset,jsep,iylstrt,iylend,'dsl')
      call output_ds(ny, nx,-target_offset,jsep,iyrstrt,iyrend,'dsr')
      if (nnreg(0).ge.7) then
        ixtl = 0
        do while (rightix(ixtl,max(topcut(1),topcut(2))).ne.nx+1.and.ixtl.lt.nx)
          ixtl=ixtl+1
        enddo
        ixtr = ixtl
        do while (leftix(ixtr,max(topcut(1),topcut(2))).ne.-2 .and. ixtr.lt.nx)
          ixtr=ixtr+1
        enddo
        call output_ds(ny,ixtl,-target_offset,jsep,iytlstrt,iytlend,'dstl')
        call output_ds(ny,ixtr,+target_offset,jsep,iytrstrt,iytrend,'dstr')
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
      nc = max(mpg%nXpt,1)
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
      if (nnreg(0).ge.7) then
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
        if(region(-1,iy,0).ne.0) write(99,*) gs(rightix(-1,iy),rightiy(-1,iy),0)*qc(rightix(-1,iy),rightiy(-1,iy),0)
      enddo
      close(99)
      open(99,file='dsRP')
      do iy=-1,ny
        if(region(nx,iy,0).ne.0) write(99,*) gs(nx,iy,0)*qc(nx,iy,0)
      enddo
      close(99)
      if (nnreg(0).ge.7) then
        open(99,file='dsTLP')
        do iy=iytlstrt,iytlend
          write(99,*) gs(ixtl,iy,0)*qc(ixtl,iy,0)
        enddo
        close(99)
        open(99,file='dsTRP')
        do iy=iytrstrt,iytrend
          write(99,*) gs(rightix(ixtr,iy),rightiy(ixtr,iy),0)*qc(rightix(ixtr,iy),rightiy(ixtr,iy),0)
        enddo
        close(99)
      endif
! WG_TODO
#endif
#ifndef NO_CDF
      if (b2time.gt.0) then
        filename='b2time.nc'
        call find_file(filename,ex)
        if (.not.ex.or.switch%b2mndr_stim.ge.0.0_R8) then
          ntstep = 0
          write(6,'(a)') trim(filename)//' will be created'
          call b2crtimecdf(filename, nCv, nc, ns, write_2d, &
            ncid, .false., iret)
          call check_cdf_status(iret)
          iret = nf_open(trim(filename),or(NF_WRITE,NF_SHARE),ncid)
          call check_cdf_status(iret)
        else if (ex.and.switch%b2mndr_stim.lt.0.0_R8) then
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
          call b2crtimecdf(filename, nCv, nc, ns, &
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
          if (.not.ex.or.switch%b2mndr_stim.ge.0.0_R8) then
            nastep = 0
            write(6,'(a)') trim(filename_av)//' will be created'
            call b2crtimecdf(filename_av, nCv, nc, ns, write_2d, &
              ncid, .true., iret)
            call check_cdf_status(iret)
            iret = nf_open(trim(filename_av),or(NF_WRITE,NF_SHARE),ncid)
            call check_cdf_status(iret)
          else if (ex.and.switch%b2mndr_stim.lt.0.0_R8) then
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
              call b2crtimecdf(filename_av, nCv, nc, ns, &
               write_2d, ncid, .true., iret)
            endif
            iret = nf_open(trim(filename_av),or(NF_WRITE,NF_SHARE),ncid)
            call check_cdf_status(iret)
          else
            nastep = 0
            write(6,'(a)') trim(filename_av)//' will be replaced'
            call b2crtimecdf(filename_av, nCv, nc, ns, &
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
      allocate (ktsepm_av(1:nncutmax))
      allocate (ktsepi_av(1:nncutmax))
      allocate (ktsepa_av(1:nncutmax))
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
      allocate (ktsepm_std(1:nncutmax))
      allocate (ktsepi_std(1:nncutmax))
      allocate (ktsepa_std(1:nncutmax))
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
      ktsepm_av = 0.0_R8
      ktsepi_av = 0.0_R8
      ktsepa_av = 0.0_R8
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
      ktsepm_std = 0.0_R8
      ktsepi_std = 0.0_R8
      ktsepa_std = 0.0_R8
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
#ifdef WG_TODO            
    ix = -1 ! 1
    ix_off  = ix + target_offset
    do iy = iylstrt,iylend
      call calc_fet(ix,iy,'L',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp,pwrtmp)
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
#endif

    fnixap = 0.0_R8; feexap = 0.0_R8; feixap = 0.0_R8; fchxap = 0.0_R8; fetxap = 0.0_R8
    nemxap = 0.0_R8; temxap = 0.0_R8; timxap = 0.0_R8; pomxap = 0.0_R8; pwmxap = 0.0_R8; tpmxap = 0.0_R8
#ifdef WG_TODO
    ix = nx ! 2
    ix_off  = ix - target_offset
    do iy = iyrstrt,iyrend
      call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp,pwrtmp)
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

    if(mpg%nXpt.ge.2) then
      ix = ixtr ! 3
      ix_off  = ix + target_offset
      do iy = iytrstrt,iytrend
        call calc_fet(ix,iy,'L',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp,pwrtmp)
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
        call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp,pwrtmp)
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
#endif

    fnisip = 0.0_R8; feesip = 0.0_R8; feisip = 0.0_R8; fchsip = 0.0_R8; fetsip = 0.0_R8
    fnisap = 0.0_R8; feesap = 0.0_R8; feisap = 0.0_R8; fchsap = 0.0_R8; fetsap = 0.0_R8
    fnisipp = 0.0_R8; feesipp = 0.0_R8; feisipp = 0.0_R8; fetsipp = 0.0_R8; fchsipp = 0.0_R8
    fnisapp = 0.0_R8; feesapp = 0.0_R8; feisapp = 0.0_R8; fetsapp = 0.0_R8; fchsapp = 0.0_R8
#ifdef WG_TODO
    if(nnreg(0).ge.3) then
      do ic = 1, mpg%nXpt
        do iy = -1,jsep
          ix = leftcut(ic) ! 5
          call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fnisipp(ic) = fnisipp(ic) + fnitmp
          feesipp(ic) = feesipp(ic) + feetmp
          feisipp(ic) = feisipp(ic) + feitmp
          fchsipp(ic) = fchsipp(ic) + fchtmp
          fetsipp(ic) = fetsipp(ic) + fettmp
          ix = rightcut(ic) ! 6
          call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fnisapp(ic) = fnisapp(ic) + fnitmp
          feesapp(ic) = feesapp(ic) + feetmp
          feisapp(ic) = feisapp(ic) + feitmp
          fchsapp(ic) = fchsapp(ic) + fchtmp
          fetsapp(ic) = fetsapp(ic) + fettmp
        enddo
        do iy = jsep+1,ny
          ix = leftcut(ic) ! 7
          call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fnisip(ic) = fnisip(ic) + fnitmp
          feesip(ic) = feesip(ic) + feetmp
          feisip(ic) = feisip(ic) + feitmp
          fchsip(ic) = fchsip(ic) + fchtmp
          fetsip(ic) = fetsip(ic) + fettmp
          ix = rightcut(ic) ! 8
          call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fnisap(ic) = fnisap(ic) + fnitmp
          feesap(ic) = feesap(ic) + feetmp
          feisap(ic) = feisap(ic) + feitmp
          fchsap(ic) = fchsap(ic) + fchtmp
          fetsap(ic) = fetsap(ic) + fettmp
        enddo
      enddo
    endif
#endif
    fniyip = 0.0_R8; feeyip = 0.0_R8; feiyip = 0.0_R8; fetyip = 0.0_R8; fchyip = 0.0_R8
    fniyap = 0.0_R8; feeyap = 0.0_R8; feiyap = 0.0_R8; fetyap = 0.0_R8; fchyap = 0.0_R8
#ifdef WG_TODO
    if(nnreg(0).eq.4 .or. nnreg(0).eq.7) then
      do ix = -1,nx
        if(region(ix,ny,0).eq.2) then
          iy = ny ! 9
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyip(1) = fniyip(1) + fnitmp
          feeyip(1) = feeyip(1) + feetmp
          feiyip(1) = feiyip(1) + feitmp
          fchyip(1) = fchyip(1) + fchtmp
          fetyip(1) = fetyip(1) + fettmp
        endif
        if(region(ix,ny,0).ge.3) then
          iy = ny ! 10
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(1) = fniyap(1) + fnitmp
          feeyap(1) = feeyap(1) + feetmp
          feiyap(1) = feiyap(1) + feitmp
          fchyap(1) = fchyap(1) + fchtmp
          fetyap(1) = fetyap(1) + fettmp
        endif
        if(region(ix,-1,0).ge.3) then
          iy = -1 ! 11
          call calc_fet(ix,iy,'B',-1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
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
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyip(1) = fniyip(1) + fnitmp
          feeyip(1) = feeyip(1) + feetmp
          feiyip(1) = feiyip(1) + feitmp
          fchyip(1) = fchyip(1) + fchtmp
          fetyip(1) = fetyip(1) + fettmp
        endif
        if(region(ix,-1,0).eq.3.or.region(ix,-1,0).eq.4) then
          iy = -1 ! 13
          call calc_fet(ix,iy,'B',-1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
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
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyip(region(ix,ny,0)/4+1) = fniyip(region(ix,ny,0)/4+1) + fnitmp
          feeyip(region(ix,ny,0)/4+1) = feeyip(region(ix,ny,0)/4+1) + feetmp
          feiyip(region(ix,ny,0)/4+1) = feiyip(region(ix,ny,0)/4+1) + feitmp
          fchyip(region(ix,ny,0)/4+1) = fchyip(region(ix,ny,0)/4+1) + fchtmp
          fetyip(region(ix,ny,0)/4+1) = fetyip(region(ix,ny,0)/4+1) + fettmp
        endif
        if(region(ix,ny,0).eq.3 .or. region(ix,ny,0).eq.8) then
          iy = ny ! 15
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(1) = fniyap(1) + fnitmp
          feeyap(1) = feeyap(1) + feetmp
          feiyap(1) = feiyap(1) + feitmp
          fchyap(1) = fchyap(1) + fchtmp
          fetyap(1) = fetyap(1) + fettmp
        endif
        if(region(ix,-1,0).eq.3 .or. region(ix,-1,0).eq.8) then
          iy = -1 ! 16
          call calc_fet(ix,iy,'B',-1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(1) = fniyap(1) + fnitmp
          feeyap(1) = feeyap(1) + feetmp
          feiyap(1) = feiyap(1) + feitmp
          fchyap(1) = fchyap(1) + fchtmp
          fetyap(1) = fetyap(1) + fettmp
        endif
        if(region(ix,ny,0).eq.4 .or. region(ix,ny,0).eq.7) then
          iy = ny ! 17
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(2) = fniyap(2) + fnitmp
          feeyap(2) = feeyap(2) + feetmp
          feiyap(2) = feiyap(2) + feitmp
          fchyap(2) = fchyap(2) + fchtmp
          fetyap(2) = fetyap(2) + fettmp
        endif
        if(region(ix,-1,0).eq.4 .or. region(ix,-1,0).eq.7) then
          iy = -1 ! 18
          call calc_fet(ix,iy,'B',-1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
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
          call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
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
        call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,switch%BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
        fniyip(1) = fniyip(1) + fnitmp
        feeyip(1) = feeyip(1) + feetmp
        feiyip(1) = feiyip(1) + feitmp
        fchyip(1) = fchyip(1) + fchtmp
        fetyip(1) = fetyip(1) + fettmp
      enddo
    endif ! nnreg check
!WG_TODO
#endif
    !
    !    other quantities related to the target plates
    !
#ifndef NO_CDF
    ! added #ifndef NO_CDF by srv 27.06.21
    nesepi = 0.0_R8; tesepi = 0.0_R8; tisepi = 0.0_R8; tpsepi = 0.0_R8; posepi = 0.0_R8
    nesepm = 0.0_R8; tesepm = 0.0_R8; tisepm = 0.0_R8; posepm = 0.0_R8; dnsepm = 0.0_R8
    dpsepm = 0.0_R8; kesepm = 0.0_R8; kisepm = 0.0_R8; vxsepm = 0.0_R8; vysepm = 0.0_R8
    vssepm = 0.0_R8; nesepa = 0.0_R8; tesepa = 0.0_R8; tisepa = 0.0_R8; tpsepa = 0.0_R8
    posepa = 0.0_R8; ktsepm = 0.0_R8; ktsepi = 0.0_R8; ktsepa = 0.0_R8;
! csc For now only identify outer midplane and targets separatrix values using flux tube concept
    if (nomp.gt.0) then
      nesepm(1) = 0.5_R8 * (dv%ne(omp(icsepomp-1))+dv%ne(omp(icsepomp)))
      tesepm(1) = 0.5_R8 * (pl%te(omp(icsepomp-1))+pl%te(omp(icsepomp)))/ev
      tisepm(1) = 0.5_R8 * (pl%ti(omp(icsepomp-1))+pl%ti(omp(icsepomp)))/ev
      posepm(1) = 0.5_R8 * (pl%po(omp(icsepomp-1))+pl%po(omp(icsepomp)))
      dnsepm(1) = 0.5_R8 * (co%dna0(omp(icsepomp-1),ismain0)+co%dna0(omp(icsepomp),ismain0))
      kesepm(1) = 0.5_R8 * (co%hce0(omp(icsepomp-1))/dv%ne(omp(icsepomp-1))+co%hce0(omp(icsepomp))/dv%ne(omp(icsepomp)))
      if (switch%tn_style.eq.0) then
        kisepm(1) = 0.5_R8 * (co%hci0(omp(icsepomp-1))/dv%ni(omp(icsepomp-1),0)+co%hci0(omp(icsepomp))/dv%ni(omp(icsepomp),0))
      else if (switch%tn_style.eq.1) then
        kisepm(1) = 0.5_R8 * (co%hci0(omp(icsepomp-1))/dv%ni(omp(icsepomp-1),1)+co%hci0(omp(icsepomp))/dv%ni(omp(icsepomp),1))
      else if (switch%tn_style.eq.2) then
        kisepm(1) = 0.5_R8 * (co%hci0(omp(icsepomp-1))/ &
                             (dv%ni(omp(icsepomp-1),0)-dv%nn(omp(icsepomp-1))) &
                             +co%hci0(omp(icsepomp))/ &
                             (dv%ni(omp(icsepomp),0)-dv%nn(omp(icsepomp-1))))
      end if
      if (switch%tn_style.eq.0.or..not.is_neutral(ismain0)) then
        dpsepm(1) = 0.5_R8 * ( &
         co%dpa0(omp(icsepomp-1),ismain0)*( &
         rt%rza(omp(icsepomp-1),ismain0)*pl%te(omp(icsepomp-1))+pl%ti(omp(icsepomp-1)))+ &
         co%dpa0(omp(icsepomp),ismain0)*( &
         rt%rza(omp(icsepomp),ismain0)*pl%te(omp(icsepomp))+pl%ti(omp(icsepomp))))
      else
        dpsepm(1) = 0.5_R8 * ( &
         co%dpa0(omp(icsepomp-1),ismain0)*pl%tn(omp(icsepomp-1))+ &
         co%dpa0(omp(icsepomp),ismain0)*pl%tn(omp(icsepomp)))
      endif
      ktsepm(1) = 0.5_R8 * (pl%kt(omp(icsepomp-1)) + pl%kt(omp(icsepomp)))/ev
      vxsepm(1) = 0.5_R8 * (co%vla0(omp(icsepomp-1),0,ismain)+ co%vla0(omp(icsepomp),0,ismain))
      vysepm(1) = 0.5_R8 * (co%vla0(omp(icsepomp-1),1,ismain)+ co%vla0(omp(icsepomp),1,ismain))
      vssepm(1) = 0.5_R8 * ( &
         co%vsa0(omp(icsepomp-1),ismain)/(mp*am(ismain)*pl%na(omp(icsepomp-1),ismain))+ &
         co%vsa0(omp(icsepomp),ismain)/(mp*am(ismain)*pl%na(omp(icsepomp),ismain)))
      iFt = mpg%cvFt(omp(icsepomp)) !separatrix flux tube
      cvtrg = mpg%ftCv(mpg%ftCvP(iFt,1)+target_offset) !inner
      if (mpg%nXpt.lt.2 .or. geometryType.eq.GEOMETRY_DDN_BOTTOM .or. &
        & geometryType.eq.GEOMETRY_LFS_SNOWFLAKE_MINUS .or. &
        & geometryType.eq.GEOMETRY_LFS_SNOWFLAKE_PLUS) then
        nesepi(1) = dv%ne(cvtrg)
        tesepi(1) = pl%te(cvtrg)/ev
        tisepi(1) = pl%ti(cvtrg)/ev
        posepi(1) = pl%po(cvtrg)
        ktsepi(1) = pl%kt(cvtrg)
      else
        nesepa(2) = dv%ne(cvtrg)
        tesepa(2) = pl%te(cvtrg)/ev
        tisepa(2) = pl%ti(cvtrg)/ev
        posepa(2) = pl%po(cvtrg)
        ktsepa(2) = pl%kt(cvtrg)
      end if
      cvtrg = mpg%ftCv(mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1-target_offset) !outer
      if (mpg%nXpt.lt.2 .or. geometryType.ne.GEOMETRY_DDN_TOP) then
        nesepa(1) = dv%ne(cvtrg)
        tesepa(1) = pl%te(cvtrg)/ev
        tisepa(1) = pl%ti(cvtrg)/ev
        posepa(1) = pl%po(cvtrg)
        ktsepa(1) = pl%kt(cvtrg)
      else
        nesepi(2) = dv%ne(cvtrg)
        tesepi(2) = pl%te(cvtrg)/ev
        tisepi(2) = pl%ti(cvtrg)/ev
        posepi(2) = pl%po(cvtrg)
        ktsepi(2) = pl%kt(cvtrg)
      end if
    end if
    if (nimp.gt.0.and.mpg%nXpt.eq.2) then
      nesepm(2) = 0.5_R8 * (dv%ne(imp(icsepimp-1))+dv%ne(imp(icsepimp)))
      tesepm(2) = 0.5_R8 * (pl%te(imp(icsepimp-1))+pl%te(imp(icsepimp)))/ev
      tisepm(2) = 0.5_R8 * (pl%ti(imp(icsepimp-1))+pl%ti(imp(icsepimp)))/ev
      posepm(2) = 0.5_R8 * (pl%po(imp(icsepimp-1))+pl%po(imp(icsepimp)))
      dnsepm(2) = 0.5_R8 * (co%dna0(imp(icsepimp-1),ismain0)+ &
        &                   co%dna0(imp(icsepimp),ismain0))
      kesepm(2) = 0.5_R8 * (co%hce0(imp(icsepimp-1))/dv%ne(imp(icsepimp-1))+ &
        &                   co%hce0(imp(icsepimp))/dv%ne(imp(icsepimp)))
      if (switch%tn_style.eq.0) then
        kisepm(2) = 0.5_R8 * &
          &  (co%hci0(imp(icsepimp-1))/dv%ni(imp(icsepimp-1),0)+ &
          &   co%hci0(imp(icsepimp))/dv%ni(imp(icsepimp),0))
      else if (switch%tn_style.eq.1) then
        kisepm(2) = 0.5_R8 * &
          &  (co%hci0(imp(icsepimp-1))/dv%ni(imp(icsepimp-1),1)+ &
          &   co%hci0(imp(icsepimp))/dv%ni(imp(icsepimp),1))
      else if (switch%tn_style.eq.2) then
        kisepm(2) = 0.5_R8 * (co%hci0(imp(icsepimp-1))/ &
          &                  (dv%ni(imp(icsepimp-1),0)-dv%nn(imp(icsepimp-1))) &
          &                  +co%hci0(imp(icsepimp))/ &
          &                  (dv%ni(imp(icsepimp),0)-dv%nn(imp(icsepimp-1))))
      end if
      if (switch%tn_style.eq.0.or..not.is_neutral(ismain0)) then
        dpsepm(2) = 0.5_R8 * ( &
         co%dpa0(imp(icsepimp-1),ismain0)*( &
         rt%rza(imp(icsepimp-1),ismain0)*pl%te(imp(icsepimp-1))+pl%ti(imp(icsepimp-1)))+ &
         co%dpa0(imp(icsepimp),ismain0)*( &
         rt%rza(imp(icsepimp),ismain0)*pl%te(imp(icsepimp))+pl%ti(imp(icsepimp))))
      else
        dpsepm(2) = 0.5_R8 * ( &
         co%dpa0(imp(icsepimp-1),ismain0)*pl%tn(imp(icsepimp-1))+ &
         co%dpa0(imp(icsepimp),ismain0)*pl%tn(imp(icsepimp)))
      endif
      ktsepm(2) = 0.5_R8 * (pl%kt(imp(icsepimp-1)) + pl%kt(imp(icsepimp)))/ev
      vxsepm(2) = 0.5_R8 * (co%vla0(imp(icsepimp-1),0,ismain)+ co%vla0(imp(icsepimp),0,ismain))
      vysepm(2) = 0.5_R8 * (co%vla0(imp(icsepimp-1),1,ismain)+ co%vla0(imp(icsepimp),1,ismain))
      vssepm(2) = 0.5_R8 * ( &
         co%vsa0(imp(icsepimp-1),ismain)/(mp*am(ismain)*pl%na(imp(icsepimp-1),ismain))+ &
         co%vsa0(imp(icsepimp),ismain)/(mp*am(ismain)*pl%na(imp(icsepimp),ismain)))
      if (geometryType.eq.GEOMETRY_CDN) then
        iFt = mpg%cvFt(imp(icsepimp)) ! HFS separatrix flux tube
        cvtrg = mpg%ftCv(mpg%ftCvP(ift,1)+target_offset) !inner
        nesepi(1) = dv%ne(cvtrg)
        tesepi(1) = pl%te(cvtrg)/ev
        tisepi(1) = pl%ti(cvtrg)/ev
        posepi(1) = pl%po(cvtrg)
        ktsepi(1) = pl%kt(cvtrg)
        cvtrg = mpg%ftCv(mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1-target_offset) !outer
        nesepi(2) = dv%ne(cvtrg)
        tesepi(2) = pl%te(cvtrg)/ev
        tisepi(2) = pl%ti(cvtrg)/ev
        posepi(2) = pl%po(cvtrg)
        ktsepi(2) = pl%kt(cvtrg)
      else
!! TODO: Must have a way to identify the flux tube outside the secondary separatrix
      end if
    end if
#endif

#ifdef WG_TODO            
#ifndef NO_CDF
    if(nnreg(0).ne.2) then
      if(xymap(-1,jsep).gt.0 .and. xymap(topix(-1,jsep),topiy(-1,jsep)).gt.0) then
        tpsepi(1) = 0.5_R8 * (target_temp(xymap(-1,jsep),1) + target_temp(xymap(topix(-1,jsep),topiy(-1,jsep)),1))
      else
        tpsepi(1) = 0.0_R8
      endif
    else
      if(xymap(topix(-1,jsep),topiy(-1,jsep)).gt.0) then
        tpsepi(1) = target_temp(xymap(topix(-1,jsep),topiy(-1,jsep)),1)
      else
        tpsepi(1) = 0.0_R8
      endif
    endif
    if(nnreg(0).ne.2) then
      if(xymap(nx,jsep).gt.0 .and. xymap(topix(nx,jsep),topiy(nx,jsep)).ge.0) then
        tpsepa(1) = 0.5_R8 * (target_temp(xymap(nx,jsep),1)+ target_temp(xymap(topix(nx,jsep),topiy(nx,jsep)),1))
      else
        tpsepa(1) = 0.0_R8
      endif
    else
      if(xymap(topix(nx,jsep),topiy(nx,jsep)).gt.0) then
        tpsepa(1) = target_temp(xymap(topix(nx,jsep),topiy(nx,jsep)),1)
      else
        tpsepa(1) = 0.0
      endif
    endif
    if(mpg%nXpt.eq.2) then
      if(xymap(ixtr,jsep).gt.0 .and. xymap(topix(ixtr,jsep),topiy(ixtr,jsep)).gt.0) then
        tpsepa(2) = 0.5_R8 *(target_temp(xymap(ixtr,jsep),1)+target_temp(xymap(topix(ixtr,jsep),topiy(ixtr,jsep)),1))
      else
        tpsepa(2) = 0.0_R8
      endif
      if(xymap(ixtl,jsep).gt.0 .and. xymap(topix(ixtl,jsep),topiy(ixtl,jsep)).gt.0) then
        tpsepi(2) = 0.5_R8 * (target_temp(xymap(ixtl,jsep),1)+ target_temp(xymap(topix(ixtl,jsep),topiy(ixtl,jsep)),1))
      else
        tpsepi(2) = 0.0_R8
      endif
    endif
!NO_CDF
#endif
    !
    temxip(1:nc) = temxip(1:nc)/ev
    timxip(1:nc) = timxip(1:nc)/ev
    temxap(1:nc) = temxap(1:nc)/ev
    timxap(1:nc) = timxap(1:nc)/ev
!WG_TODO
#endif

    tmne(1)=0.0_R8
    tmte(1)=0.0_R8
    tmti(1)=0.0_R8
    tmvol=0.0_R8
    do iCv = 1, nCv
      if(mpg%cvReg(iCv).ne.0) then
        tmne(1)=tmne(1)+dv%ne(iCv)*geo%cvVol(iCv)
        tmte(1)=tmte(1)+pl%te(iCv)*dv%ne(iCv)*geo%cvVol(iCv)
        tmti(1)=tmti(1)+pl%ti(iCv)*dv%ni(iCv,0)*geo%cvVol(iCv)
        tmvol=tmvol+geo%cvVol(iCv)
      endif
    enddo
    tmne(1)=tmne(1)
    tmte(1)=tmte(1)/ev
    tmti(1)=tmti(1)/ev

#ifndef NO_CDF
    tmhacore(1)=0.0_R8
    tmhasol(1)=0.0_R8
    tmhadiv(1)=0.0_R8
#ifdef WG_TODO
#ifdef B25_EIRENE
    ! note no emission in guard cells and offset of array by 1
    do iy=-1,ny
      do ix=-1,nx
        if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and. bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 ) then
          if(on_closed_surface(ix,iy)) then
            tmhacore(1)=tmhacore(1)+(emiss(ix+1,iy+1,1,1)+emissmol(ix+1,iy+1,1,1))*vol(ix,iy)
          elseif((region(ix,iy,0).eq.5 .or. region(ix,iy,0).eq.6) .and. nnreg(0).eq.7) then
            tmhadiv(1)=tmhadiv(1) + (emiss(ix+1,iy+1,1,1)+emissmol(ix+1,iy+1,1,1))*vol(ix,iy)
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
        call batch_average(nncutmax,ktsepm,ktsepm_av,itim,ntim_batch)
        call batch_average(nncutmax,ktsepi,ktsepi_av,itim,ntim_batch)
        call batch_average(nncutmax,ktsepa,ktsepa_av,itim,ntim_batch)
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
        call batch_average_sq(nncutmax,ktsepm,ktsepm_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,ktsepi,ktsepi_std,itim,ntim_batch)
        call batch_average_sq(nncutmax,ktsepa,ktsepa_std,itim,ntim_batch)
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
        call rwcdf(rw,ncid,'ne2d',(/1,1/),dv%ne,iret)
        call rwcdf(rw,ncid,'te2d',(/1,1/),pl%te,iret)
        call rwcdf(rw,ncid,'ti2d',(/1,1/),pl%ti,iret)
        if (write_2d .ge. 2) then
          call rwcdf(rw,ncid,'po2d',(/1,1/),pl%po,iret)
          call rwcdf(rw,ncid,'kin2d',(/1,1/),dv%kinrgy,iret)
          call rwcdf(rw,ncid,'rsahi2d',(/1,1/),srw%rsahi,iret)
          call rwcdf(rw,ncid,'rsana2d',(/1,1/),srw%rsana,iret)
          call rwcdf(rw,ncid,'rrahi2d',(/1,1/),srw%rrahi,iret)
          call rwcdf(rw,ncid,'rrana2d',(/1,1/),srw%rrana,iret)
          call rwcdf(rw,ncid,'rcxhi2d',(/1,1/),srw%rcxhi,iret)
          call rwcdf(rw,ncid,'rcxna2d',(/1,1/),srw%rcxna,iret)
          call rwcdf(rw,ncid,'rqrad2d',(/1,1/),srw%rqrad,iret)
          call rwcdf(rw,ncid,'fhe2d',(/1,1/),dv%fhe,iret)
          call rwcdf(rw,ncid,'fhi2d',(/1,1/),dv%fhi,iret)
          call rwcdf(rw,ncid,'fch2d',(/1,1/),dv%fch,iret)
          call rwcdf(rw,ncid,'fna2d',(/1,1,1/),dv%fna,iret)
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
      call rwcdf(rw,ncid,'ktsepi',imap,ktsepi,iret)
      call rwcdf(rw,ncid,'nesepm',imap,nesepm,iret)
      call rwcdf(rw,ncid,'tesepm',imap,tesepm,iret)
      call rwcdf(rw,ncid,'tisepm',imap,tisepm,iret)
      call rwcdf(rw,ncid,'posepm',imap,posepm,iret)
      call rwcdf(rw,ncid,'ktsepm',imap,ktsepm,iret)
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
      call rwcdf(rw,ncid,'ktsepa',imap,ktsepa,iret)
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
#ifdef WG_TODO
      imap(1)=nx+2     ! bl
      imap(2)=1
      call rwcdf(rw,ncid,'ne3dl',imap,ne(-1+target_offset,iylstrt),iret)
      call rwcdf(rw,ncid,'te3dl',imap,te(-1+target_offset,iylstrt),iret)
      call rwcdf(rw,ncid,'ti3dl',imap,ti(-1+target_offset,iylstrt),iret)
      call rwcdf(rw,ncid,'po3dl',imap,po(-1+target_offset,iylstrt),iret)
      call rwcdf(rw,ncid,'fn3dl',imap,fna(0,iylstrt,0,0,ismain),iret)
      call rwcdf(rw,ncid,'fe3dl',imap,fhe(0,iylstrt,0,0),iret)
      call rwcdf(rw,ncid,'fi3dl',imap,fhi(0,iylstrt,0,0),iret)
      call rwcdf(rw,ncid,'fc3dl',imap,fch(0,iylstrt,0,0),iret)
      call rwcdf(rw,ncid,'fl3dl',imap,fne(0,iylstrt,0,0),iret)
      call rwcdf(rw,ncid,'fo3dl',imap,fni(0,iylstrt,0,0),iret)
#endif
      if (nimp.gt.0) then
        imap(1)=1     ! i
        imap(2)=1
        allocate(slice(nimp))
        slice(1:nimp)=dv%ne(imp(1:nimp))
        call rwcdf(rw,ncid,'ne3di',imap,slice,iret)
        slice(1:nimp)=pl%te(imp(1:nimp))
        call rwcdf(rw,ncid,'te3di',imap,slice,iret)
        slice(1:nimp)=pl%ti(imp(1:nimp))
        call rwcdf(rw,ncid,'ti3di',imap,slice,iret)
        slice(1:nimp)=pl%po(imp(1:nimp))
        call rwcdf(rw,ncid,'po3di',imap,slice,iret)
        deallocate(slice)
      end if
      if (nomp.gt.0) then
        imap(1)=1     ! a
        imap(2)=1
        allocate(slice(nomp))
        slice(1:nomp)=dv%ne(omp(1:nomp))
        call rwcdf(rw,ncid,'ne3da',imap,slice,iret)
        slice(1:nomp)=pl%te(omp(1:nomp))
        call rwcdf(rw,ncid,'te3da',imap,slice,iret)
        slice(1:nomp)=pl%ti(omp(1:nomp))
        call rwcdf(rw,ncid,'ti3da',imap,slice,iret)
        slice(1:nomp)=pl%po(omp(1:nomp))
        call rwcdf(rw,ncid,'po3da',imap,slice,iret)
        deallocate(slice)
      end if
#ifdef WG_TODO
      imap(1)=nx+2     ! br
      imap(2)=1
      call rwcdf(rw,ncid,'ne3dr',imap,ne(nx-target_offset,iyrstrt),iret)
      call rwcdf(rw,ncid,'te3dr',imap,te(nx-target_offset,iyrstrt),iret)
      call rwcdf(rw,ncid,'ti3dr',imap,ti(nx-target_offset,iyrstrt),iret)
      call rwcdf(rw,ncid,'po3dr',imap,po(nx-target_offset,iyrstrt),iret)
      call rwcdf(rw,ncid,'fn3dr',imap,fna(nx,iyrstrt,0,0,ismain),iret)
      call rwcdf(rw,ncid,'fe3dr',imap,fhe(nx,iyrstrt,0,0),iret)
      call rwcdf(rw,ncid,'fi3dr',imap,fhi(nx,iyrstrt,0,0),iret)
      call rwcdf(rw,ncid,'fc3dr',imap,fch(nx,iyrstrt,0,0),iret)
      call rwcdf(rw,ncid,'fl3dr',imap,fne(nx,iyrstrt,0,0),iret)
      call rwcdf(rw,ncid,'fo3dr',imap,fni(nx,iyrstrt,0,0),iret)
      if (nnreg(0).ge.7) then
        imap(1)=nx+2     ! tr
        imap(2)=1
        call rwcdf(rw,ncid,'ne3dtr',imap,ne(ixtr+target_offset,iytrstrt),iret)
        call rwcdf(rw,ncid,'te3dtr',imap,te(ixtr+target_offset,iytrstrt),iret)
        call rwcdf(rw,ncid,'ti3dtr',imap,ti(ixtr+target_offset,iytrstrt),iret)
        call rwcdf(rw,ncid,'po3dtr',imap,po(ixtr+target_offset,iytrstrt),iret)
        call rwcdf(rw,ncid,'fn3dtr',imap,fna(ixtr+1,iytrstrt,0,0,ismain),iret)
        call rwcdf(rw,ncid,'fe3dtr',imap,fhe(ixtr+1,iytrstrt,0,0),iret)
        call rwcdf(rw,ncid,'fi3dtr',imap,fhi(ixtr+1,iytrstrt,0,0),iret)
        call rwcdf(rw,ncid,'fc3dtr',imap,fch(ixtr+1,iytrstrt,0,0),iret)
        call rwcdf(rw,ncid,'fl3dtr',imap,fne(ixtr+1,iytrstrt,0,0),iret)
        call rwcdf(rw,ncid,'fo3dtr',imap,fni(ixtr+1,iytrstrt,0,0),iret)
        imap(1)=nx+2     ! tl
        imap(2)=1
        call rwcdf(rw,ncid,'ne3dtl',imap,ne(ixtl-target_offset,iytlstrt),iret)
        call rwcdf(rw,ncid,'te3dtl',imap,te(ixtl-target_offset,iytlstrt),iret)
        call rwcdf(rw,ncid,'ti3dtl',imap,ti(ixtl-target_offset,iytlstrt),iret)
        call rwcdf(rw,ncid,'po3dtl',imap,po(ixtl-target_offset,iytlstrt),iret)
        call rwcdf(rw,ncid,'fn3dtl',imap,fna(ixtl,iytlstrt,0,0,ismain),iret)
        call rwcdf(rw,ncid,'fe3dtl',imap,fhe(ixtl,iytlstrt,0,0),iret)
        call rwcdf(rw,ncid,'fi3dtl',imap,fhi(ixtl,iytlstrt,0,0),iret)
        call rwcdf(rw,ncid,'fc3dtl',imap,fch(ixtl,iytlstrt,0,0),iret)
        call rwcdf(rw,ncid,'fl3dtl',imap,fne(ixtl,iytlstrt,0,0),iret)
        call rwcdf(rw,ncid,'fo3dtl',imap,fni(ixtl,iytlstrt,0,0),iret)
      endif
    !
      imap(1)=1
      imap(2)=1
      slice=0.0_R8
      if (ismain0.ne.ismain) then
        slice(iylstrt:iylend)=na(-1+target_offset,iylstrt:iylend,ismain0)
        do iatm=1,nnatmi
          if (b2espcr(ismain0).eq.latmscl(iatm)) &
           &  slice(0:ny-1)=slice(0:ny-1)+dab2(1,1:ny,iatm,1)
        enddo
      endif
      call rwcdf(rw,ncid,'an3dl',imap,slice(iylstrt),iret)
#endif
      imap(1)=1
      imap(2)=1
      if (nimp.gt.0) then
        allocate(slice(nimp))
        if (ismain0.ne.ismain) then
          slice(1:nimp)=pl%na(imp(1:nimp),ismain0)
          do iatm=1,nnatmi
            if (b2espcr(ismain0).eq.latmscl(iatm)) &
              &  slice(1:nimp)=slice(1:nimp)+dab2(imp(1:nimp),iatm,1)
          enddo
        endif
        call rwcdf(rw,ncid,'an3di',imap,slice,iret)
        deallocate(slice)
      end if
      if (nomp.gt.0) then
        allocate(slice(nomp))
        if (ismain0.ne.ismain) then
          slice(1:nomp)=pl%na(omp(1:nomp),ismain0)
          do iatm=1,nnatmi
            if (b2espcr(ismain0).eq.latmscl(iatm)) &
              &  slice(1:nomp)=slice(1:nomp)+dab2(omp(1:nomp),iatm,1)
          enddo
        endif
        call rwcdf(rw,ncid,'an3da',imap,slice,iret)
        deallocate(slice)
      end if
#ifdef WG_TODO
      slice=0.0_R8
      if (ismain0.ne.ismain) then
        slice(iyrstrt:iyrend)=na(nx-target_offset,iyrstrt:iyrend,ismain0)
        do iatm=1,nnatmi
          if (b2espcr(ismain0).eq.latmscl(iatm)) &
           &  slice(0:ny-1)=slice(0:ny-1)+dab2(nx,1:ny,iatm,1)
        enddo
      endif
      call rwcdf(rw,ncid,'an3dr',imap,slice(iyrstrt),iret)
      if (nnreg(0).ge.7) then
        slice=0.0_R8
        if (ismain0.ne.ismain) then
          slice(iytlstrt:iytlend)= na(ixtl-target_offset,iytlstrt:iytlend,ismain0)
          do iatm=1,nnatmi
            if (b2espcr(ismain0).eq.latmscl(iatm)) &
             &  slice(0:ny-1)=slice(0:ny-1)+dab2(ixtl,1:ny,iatm,1)
          enddo
        endif
        call rwcdf(rw,ncid,'an3dtl',imap,slice(iytlstrt),iret)
        slice=0.0_R8
        if (ismain0.ne.ismain) then
          slice(iytrstrt:iytrend)= na(ixtr+target_offset,iytrstrt:iytrend,ismain0)
          do iatm=1,nnatmi
            if (b2espcr(ismain0).eq.latmscl(iatm)) &
             &  slice(0:ny-1)=slice(0:ny-1)+dab2(ixtr+1,1:ny,iatm,1)
          enddo
        endif
        call rwcdf(rw,ncid,'an3dtr',imap,slice(iytrstrt),iret)
      endif
#endif
      if (nnmoli.gt.0) then
#ifdef WG_TODO
        slice(0:ny-1)=dmb2(1,1:ny,1,1)
        call rwcdf(rw,ncid,'mn3dl',imap,slice,iret)
#endif
        if (nimp.gt.0) then
          allocate(slice(nimp))
          slice(1:nimp)=dmb2(imp(1:nimp),1,1)
          call rwcdf(rw,ncid,'mn3di',imap,slice,iret)
          deallocate(slice)
        end if
        if (nomp.gt.0) then
          allocate(slice(nomp))
          slice(1:nomp)=dmb2(omp(1:nomp),1,1)
          call rwcdf(rw,ncid,'mn3da',imap,slice,iret)
          deallocate(slice)
        end if
#ifdef WG_TODO
        slice(0:ny-1)=dmb2(nx,1:ny,1,1)
        call rwcdf(rw,ncid,'mn3dr',imap,slice,iret)
        if (nnreg(0).ge.7) then
          slice(0:ny-1)=dmb2(ixtl,1:ny,1,1)
          call rwcdf(rw,ncid,'mn3dtl',imap,slice,iret)
          slice(0:ny-1)=dmb2(ixtr+1,1:ny,1,1)
          call rwcdf(rw,ncid,'mn3dtr',imap,slice,iret)
        endif
#endif
      endif
    !
#ifdef WG_TODO
      slice=0.0_R8
      do iy = iylstrt, iylend
        call calc_fet(-1,iy,'L',1._R8,nx,ny,ns,ismain,switch%BoRiS,slice(iy))
      enddo
      call rwcdf(rw,ncid,'ft3dl',imap,slice(iylstrt),iret)

      slice=0.0_R8
      do iy = iyrstrt, iyrend
        call calc_fet(nx,iy,'R',1._R8,nx,ny,ns,ismain,switch%BoRiS,slice(iy))
      enddo
      call rwcdf(rw,ncid,'ft3dr',imap,slice(iyrstrt),iret)
      if (nnreg(0).ge.7) then
        slice=0.0_R8
        do iy = iytlstrt, iytlend
          call calc_fet(ixtl,iy,'R',1._R8,nx,ny,ns,ismain,switch%BoRiS,slice(iy))
        enddo
        call rwcdf(rw,ncid,'ft3dtl',imap,slice(iytlstrt),iret)

        slice=0.0_R8
        do iy = iytrstrt, iytrend
          call calc_fet(ixtr,iy,'L',1._R8,nx,ny,ns,ismain,switch%BoRiS,slice(iy))
        enddo
        call rwcdf(rw,ncid,'ft3dtr',imap,slice(iytrstrt),iret)
      endif
#endif
      if (nimp.gt.0) then
        allocate(slice(nimp))
        slice(1:nimp)=co%dna0(imp(1:nimp),ismain)
        call rwcdf(rw,ncid,'dn3di',imap,slice,iret)
        if (switch%tn_style.eq.0.or..not.is_neutral(ismain0)) then
          slice(1:nimp)=co%dpa0(imp(1:nimp),ismain0)* &
            &  (rt%rza(imp(1:nimp),ismain0)*pl%te(imp(1:nimp))+pl%ti(imp(1:nimp)))
        else
          slice(1:nimp)=co%dpa0(imp(1:nimp),ismain0)*pl%tn(imp(1:nimp))
        end if
        call rwcdf(rw,ncid,'dp3di',imap,slice,iret)
        slice(1:nimp)=co%hce0(imp(1:nimp))/dv%ne(imp(1:nimp))
        call rwcdf(rw,ncid,'ke3di',imap,slice,iret)
        if (switch%tn_style.eq.0) then
          slice(1:nimp)=co%hci0(imp(1:nimp))/dv%ni(imp(1:nimp),0)
        else if (switch%tn_style.eq.1) then
          slice(1:nimp)=co%hci0(imp(1:nimp))/dv%ni(imp(1:nimp),1)
        else if (switch%tn_style.eq.2) then
          slice(1:nimp)=co%hci0(imp(1:nimp))/(dv%ni(imp(1:nimp),0)-dv%nn(imp(1:nimp)))
        end if
        call rwcdf(rw,ncid,'ki3di',imap,slice,iret)
        slice(1:nimp)=co%vla0(imp(1:nimp),0,ismain)
        call rwcdf(rw,ncid,'vx3di',imap,slice,iret)
        slice(1:nimp)=co%vla0(imp(1:nimp),1,ismain)
        call rwcdf(rw,ncid,'vy3di',imap,slice,iret)
        slice(1:nimp)=co%vsa0(imp(1:nimp),ismain)/(mp*am(ismain)*pl%na(imp(1:nimp),ismain))
        call rwcdf(rw,ncid,'vs3di',imap,slice,iret)
        deallocate(slice)
      end if
      if (nomp.gt.0) then
        allocate(slice(nomp))
        slice(1:nomp)=co%dna0(omp(1:nomp),ismain)
        call rwcdf(rw,ncid,'dn3da',imap,slice,iret)
        if (switch%tn_style.eq.0.or..not.is_neutral(ismain0)) then
          slice(1:nomp)=co%dpa0(omp(1:nomp),ismain0)* &
            &  (rt%rza(omp(1:nomp),ismain0)*pl%te(omp(1:nomp))+pl%ti(omp(1:nomp)))
        else
          slice(1:nomp)=co%dpa0(omp(1:nomp),ismain0)*pl%tn(omp(1:nomp))
        end if
        call rwcdf(rw,ncid,'dp3da',imap,slice,iret)
        slice(1:nomp)=co%hce0(omp(1:nomp))/dv%ne(omp(1:nomp))
        call rwcdf(rw,ncid,'ke3da',imap,slice,iret)
        if (switch%tn_style.eq.0) then
          slice(1:nomp)=co%hci0(omp(1:nomp))/dv%ni(omp(1:nomp),0)
        else if (switch%tn_style.eq.1) then
          slice(1:nomp)=co%hci0(omp(1:nomp))/dv%ni(omp(1:nomp),1)
        else if (switch%tn_style.eq.2) then
          slice(1:nomp)=co%hci0(omp(1:nomp))/(dv%ni(omp(1:nomp),0)-dv%nn(omp(1:nomp)))
        end if
        call rwcdf(rw,ncid,'ki3da',imap,slice,iret)
        slice(1:nomp)=co%vla0(omp(1:nomp),0,ismain)
        call rwcdf(rw,ncid,'vx3da',imap,slice,iret)
        slice(1:nomp)=co%vla0(omp(1:nomp),1,ismain)
        call rwcdf(rw,ncid,'vy3da',imap,slice,iret)
        slice(1:nomp)=co%vsa0(omp(1:nomp),ismain)/(mp*am(ismain)*pl%na(omp(1:nomp),ismain))
        call rwcdf(rw,ncid,'vs3da',imap,slice,iret)
        deallocate(slice)
      end if
#ifdef WG_TODO
      slice(-1:ny)=fllim0fhi(jxi,-1:ny,1,1,ismain0)
      call rwcdf(rw,ncid,'lh3di',imap,slice,iret)
      slice(-1:ny)=fllim0fhi(jxa,-1:ny,1,1,ismain0)
      call rwcdf(rw,ncid,'lh3da',imap,slice,iret)
      slice(-1:ny)=fllim0fna(jxi,-1:ny,1,1,ismain0)
      call rwcdf(rw,ncid,'ln3di',imap,slice,iret)
      slice(-1:ny)=fllim0fna(jxa,-1:ny,1,1,ismain0)
      call rwcdf(rw,ncid,'ln3da',imap,slice,iret)
#endif
    !
      imap(1)=1
      imap(2)=1
      call rwcdf(rw,ncid,'tpsepi',imap,tpsepi,iret)
      call rwcdf(rw,ncid,'tpsepa',imap,tpsepa,iret)
      call rwcdf(rw,ncid,'tpmxip',imap,tpmxip,iret)
      call rwcdf(rw,ncid,'tpmxap',imap,tpmxap,iret)
#ifdef WG_TODO
      imap(1)=1
      imap(2)=1
      slice=0.0_R8
      if(minval(xymap(-1,0:ny-1)).gt.0) then
        slice(-1)=target_temp(xymap(-1,0),1)
        do iy=0,ny-1
          slice(iy)=target_temp(xymap(-1,iy),1)
        end do
        slice(ny)=target_temp(xymap(-1,ny-1),1)
      endif
      call rwcdf(rw,ncid,'tp3dl',imap,slice,iret)
      slice=0.0_R8
      if(minval(xymap(nx,0:ny-1)).gt.0) then
        slice(-1)=target_temp(xymap(nx,0),1)
        do iy=0,ny-1
          slice(iy)=target_temp(xymap(nx,iy),1)
        end do
        slice(ny)=target_temp(xymap(nx,ny-1),1)
      endif
      call rwcdf(rw,ncid,'tp3dr',imap,slice,iret)
      if (nnreg(0).ge.7) then
        slice=0.0_R8
        if(minval(xymap(ixtl,0:ny-1)).gt.0) then
          slice(-1)=target_temp(xymap(ixtl,0),1)
          do iy=0,ny-1
            slice(iy)=target_temp(xymap(ixtl,iy),1)
          end do
          slice(ny)=target_temp(xymap(ixtl,ny-1),1)
        endif
        call rwcdf(rw,ncid,'tp3dtl',imap,slice,iret)
        slice=0.0_R8
        if(minval(xymap(ixtr,0:ny-1)).gt.0) then
          slice(-1)=target_temp(xymap(ixtr,0),1)
          do iy=0,ny-1
            slice(iy)=target_temp(xymap(ixtr,iy),1)
          end do
          slice(ny)=target_temp(xymap(ixtr,ny-1),1)
        endif
        call rwcdf(rw,ncid,'tp3dtr',imap,slice,iret)
      endif
#endif
      iret = nf_close(ncid)
      call check_cdf_status(iret)
    endif !lwti

!wdk only write batch data if lwav is true
    if (lwav) then
      rw = 'write'
      iret = nf_open(filename_av, or(NF_WRITE,NF_SHARE), ncid)
      call check_cdf_status(iret)
!wdk compute the standard deviation from average and average of squares
      fac = rratio(ntim_batch,ntim_batch - 1)
      nesepm_std(1:nc) = (abs(nesepm_std(1:nc) - nesepm_av(1:nc)**2)*fac)**0.5
      tesepm_std(1:nc) = (abs(tesepm_std(1:nc) - tesepm_av(1:nc)**2)*fac)**0.5
      tisepm_std(1:nc) = (abs(tisepm_std(1:nc) - tisepm_av(1:nc)**2)*fac)**0.5
      posepm_std(1:nc) = (abs(posepm_std(1:nc) - posepm_av(1:nc)**2)*fac)**0.5
      nesepi_std(1:nc) = (abs(nesepi_std(1:nc) - nesepi_av(1:nc)**2)*fac)**0.5
      tesepi_std(1:nc) = (abs(tesepi_std(1:nc) - tesepi_av(1:nc)**2)*fac)**0.5
      tisepi_std(1:nc) = (abs(tisepi_std(1:nc) - tisepi_av(1:nc)**2)*fac)**0.5
      posepi_std(1:nc) = (abs(posepi_std(1:nc) - posepi_av(1:nc)**2)*fac)**0.5
      nesepa_std(1:nc) = (abs(nesepa_std(1:nc) - nesepa_av(1:nc)**2)*fac)**0.5
      tesepa_std(1:nc) = (abs(tesepa_std(1:nc) - tesepa_av(1:nc)**2)*fac)**0.5
      tisepa_std(1:nc) = (abs(tisepa_std(1:nc) - tisepa_av(1:nc)**2)*fac)**0.5
      posepa_std(1:nc) = (abs(posepa_std(1:nc) - posepa_av(1:nc)**2)*fac)**0.5
      nemxip_std(1:nc) = (abs(nemxip_std(1:nc) - nemxip_av(1:nc)**2)*fac)**0.5
      temxip_std(1:nc) = (abs(temxip_std(1:nc) - temxip_av(1:nc)**2)*fac)**0.5
      timxip_std(1:nc) = (abs(timxip_std(1:nc) - timxip_av(1:nc)**2)*fac)**0.5
      pomxip_std(1:nc) = (abs(pomxip_std(1:nc) - pomxip_av(1:nc)**2)*fac)**0.5
      nemxap_std(1:nc) = (abs(nemxap_std(1:nc) - nemxap_av(1:nc)**2)*fac)**0.5
      temxap_std(1:nc) = (abs(temxap_std(1:nc) - temxap_av(1:nc)**2)*fac)**0.5
      timxap_std(1:nc) = (abs(timxap_std(1:nc) - timxap_av(1:nc)**2)*fac)**0.5
      pomxap_std(1:nc) = (abs(pomxap_std(1:nc) - pomxap_av(1:nc)**2)*fac)**0.5
      ktsepm_std(1:nc) = (abs(ktsepm_std(1:nc) - ktsepm_av(1:nc)**2)*fac)**0.5
      ktsepi_std(1:nc) = (abs(ktsepi_std(1:nc) - ktsepi_av(1:nc)**2)*fac)**0.5
      ktsepa_std(1:nc) = (abs(ktsepa_std(1:nc) - ktsepa_av(1:nc)**2)*fac)**0.5

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
      call rwcdf(rw,ncid,'ktsepm_av',imap,ktsepm_av,iret)
      call rwcdf(rw,ncid,'nesepi_av',imap,nesepi_av,iret)
      call rwcdf(rw,ncid,'tesepi_av',imap,tesepi_av,iret)
      call rwcdf(rw,ncid,'tisepi_av',imap,tisepi_av,iret)
      call rwcdf(rw,ncid,'posepi_av',imap,posepi_av,iret)
      call rwcdf(rw,ncid,'ktsepi_av',imap,ktsepi_av,iret)
      call rwcdf(rw,ncid,'nesepa_av',imap,nesepa_av,iret)
      call rwcdf(rw,ncid,'tesepa_av',imap,tesepa_av,iret)
      call rwcdf(rw,ncid,'tisepa_av',imap,tisepa_av,iret)
      call rwcdf(rw,ncid,'posepa_av',imap,posepa_av,iret)
      call rwcdf(rw,ncid,'ktsepa_av',imap,ktsepa_av,iret)
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
      call rwcdf(rw,ncid,'ktsepm_std',imap,ktsepm_std,iret)
      call rwcdf(rw,ncid,'nesepi_std',imap,nesepi_std,iret)
      call rwcdf(rw,ncid,'tesepi_std',imap,tesepi_std,iret)
      call rwcdf(rw,ncid,'tisepi_std',imap,tisepi_std,iret)
      call rwcdf(rw,ncid,'posepi_std',imap,posepi_std,iret)
      call rwcdf(rw,ncid,'ktsepi_std',imap,ktsepi_std,iret)
      call rwcdf(rw,ncid,'nesepa_std',imap,nesepa_std,iret)
      call rwcdf(rw,ncid,'tesepa_std',imap,tesepa_std,iret)
      call rwcdf(rw,ncid,'tisepa_std',imap,tisepa_std,iret)
      call rwcdf(rw,ncid,'posepa_std',imap,posepa_std,iret)
      call rwcdf(rw,ncid,'ktsepa_std',imap,ktsepa_std,iret)
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
    endif !lwav
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
  deallocate(ktsepm_av)
  deallocate(ktsepi_av)
  deallocate(ktsepa_av)
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
  deallocate(ktsepm_std)
  deallocate(ktsepi_std)
  deallocate(ktsepa_std)

  return
  end subroutine dealloc_b2mod_mwti
#endif
!
  subroutine output_ds_cv(geo,nlist,cvlist,isep,filename)
    use b2us_geo
    implicit none
    type (geometry), intent(in) :: geo
    integer nlist,isep
    integer cvlist(nlist)
    real (kind=R8) :: &
         ds(nlist), ds_offset
    character*(*) filename
    integer i
    intrinsic sqrt

    ds(1)= 0.0_R8
    do i=2,nlist
      ds(i)=ds(i-1)+ &
           sqrt((geo%cvX(cvlist(i))-geo%cvX(cvlist(i-1)))**2+ &
                (geo%cvY(cvlist(i))-geo%cvY(cvlist(i-1)))**2)
    enddo
    if(isep.ne.0) then
      ds_offset=(ds(isep)+ds(isep+1))/2.0_R8
      do i=1,nlist
        ds(i)=ds(i)-ds_offset
      enddo
    endif
    open(99,file=filename)
    do i=1,nlist
      write(99,*) ds(i)
    enddo
    close(99)
    return
  end subroutine output_ds_cv

  subroutine output_ds_fc(geo,nlist,fclist,isep,filename)
    use b2us_geo
    implicit none
    type (geometry), intent(in) :: geo
    integer nlist,isep
    integer fclist(nlist)
    real (kind=R8) :: &
         ds(nlist), ds_offset
    character*(*) filename
    integer i

    ds(1)= 0.5_R8 * geo%fcHt(fclist(1))
    do i=2,nlist
      ds(i)=ds(i-1)+ 0.5 * (geo%fcHt(fclist(i-1)) + geo%fcHt(fclist(i)))
    enddo
    if(isep.ne.0) then
      ds_offset=(ds(isep)+ds(isep+1))/2.0_R8
      do i=1,nlist
        ds(i)=ds(i)-ds_offset
      enddo
    endif
    open(99,file=filename)
    do i=1,nlist
      write(99,*) ds(i)
    enddo
    close(99)
    return
  end subroutine output_ds_fc

#ifdef WG_TODO
  subroutine calc_fet(ix,iy,side,fac_flux,nx,ny,ns,ismain,BoRiS,fet,fni0,fee0,fei0,fch0,pwr)
    use b2mod_plasma   , only : ti, te, fna, fne, fhe, fhi, fch, fhm, fhp
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
    real(kind=R8) :: kintmp, rpttmp, tif, tef, taf
    real(kind=R8) :: h(-1:nx,-1:ny)
    ! Procedures
    external xerrab

    ! Computation
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
    if (present(fni0)) fni0 = fac_flux*fna(ix_flux,iy_flux,idir,idir,ismain)
    if (present(fee0)) fee0 = fac_flux*fhe(ix_flux,iy_flux,idir,idir)
    if (present(fei0)) fei0 = fac_flux*fhi(ix_flux,iy_flux,idir,idir)
    if (present(fch0)) fch0 = fac_flux*fch(ix_flux,iy_flux,idir,idir)
    fet = fac_flux*(fhe(ix_flux,iy_flux,idir,idir) + fhi(ix_flux,iy_flux,idir,idir) + fhi_ext(ix_flux,iy_flux,idir,idir))
    tef = (te(ix_adj,iy_adj)*h(ix,iy)+te(ix,iy)*h(ix_adj,iy_adj))/(h(ix,iy)+h(ix_adj,iy_adj))
    tif = (ti(ix_adj,iy_adj)*h(ix,iy)+ti(ix,iy)*h(ix_adj,iy_adj))/(h(ix,iy)+h(ix_adj,iy_adj))
    fet = fet + fac_flux*fne(ix_flux,iy_flux,idir,idir)*tef*(1.0_R8-BoRiS)
    do is=0,ns-1
      fet = fet + fac_flux*(fhm(ix_flux,iy_flux,idir,idir,is)+tif)*(1.0_R8-BoRiS) + fac_flux*fhp(ix_flux,iy_flux,idir,idir,is)
    enddo
    do is=0,ns_ext-1
      kintmp = 0.5_R8*am_ext(is)*mp*(ua_ext(ix,iy,is)**2 * h(ix_adj,iy_adj)+  &
           ua_ext(ix_adj,iy_adj,is)**2*h(ix,iy))/(h(ix_adj,iy_adj)+h(ix,iy))
      rpttmp = (pt_ext(ix,iy,is)*h(ix_adj,iy_adj)+pt_ext(ix_adj,iy_adj,is)*h(ix,iy))/(h(ix_adj,iy_adj)+h(ix,iy))
      taf = (ta_ext(ix_adj,iy_adj,is)*h(ix,iy)+ta_ext(ix,iy,is)*h(ix_adj,iy_adj))/(h(ix,iy)+h(ix_adj,iy_adj))
      fet = fet + fac_flux*(rpttmp*ev + (kintmp+taf)*(1.0_R8-BoRiS))*fa_ext(ix_flux,iy_flux,idir,idir,is)
    enddo
    if (present(pwr)) pwr = Abs(fet)/gs(ix_flux,iy_flux,idir)
  end subroutine calc_fet
#endif

end module b2mod_mwti

!!!Local Variables:
!!! mode: f90
!!! End:
