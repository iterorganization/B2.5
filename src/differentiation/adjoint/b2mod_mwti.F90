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
  subroutine b2mwti (itim, tim, &
#ifndef NO_CDF
                     ntim, b2time, ntim_batch, &
#endif
                     nCv, nFc, ns, nncutmax, geo, mpg, switch, &
                     pl, dv, co, rt, srw, ext, &
                     ismain, ismain0, lwti, lwav, luav)
!    use b2mod_geo_diff
!    use b2mod_plasma_diff
!    use b2mod_rates
!    use b2mod_residuals
!    use b2mod_sources
!    use b2mod_transport
!    use b2mod_anomalous_transport
    use b2mod_neutrals_namelist_diff
!    use b2mod_work
!    use b2mod_indirect_diff
    use b2mod_constants
!    use b2mod_tallies
!    use b2mod_wall
    use b2mod_b2cmpa_diff
!    use b2mod_external_diff
    use b2us_geo_diff
    use b2us_map_diff
    use b2us_plasma_diff
    use b2mod_geometry &
    , only : geometryID, GEOMETRY_CDN
    use b2mod_user_namelist_diff &
    , only : omp, imp, nimp, nomp, icsepimp
#ifndef NO_CDF
    use b2mod_geometry &
    , only : GEOMETRY_DDN_TOP, GEOMETRY_DDN_BOTTOM, &
             GEOMETRY_LFS_SNOWFLAKE_PLUS, GEOMETRY_LFS_SNOWFLAKE_MINUS
#endif
    use b2mod_user_namelist_diff &
    , only : icsepomp
    use b2mod_switches_diff
#ifndef SOLPS4_3
#ifdef B25_EIRENE
    use eirmod_wneutrals
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
    type (B2StateExt), intent (in) :: ext
    integer, intent(in) :: itim, nCv, nFc, ns, nncutmax, ismain, ismain0
    real (kind=R8), intent(in) :: tim
    logical, intent(in) :: lwti, lwav, luav
    !   ..output arguments (unspecified on entry)
    !     (none)
    !   ..common blocks
#ifndef NO_CDF
    integer, intent(in) :: ntim, b2time, ntim_batch
#   include <netcdf.inc>
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
         pomxip(nncutmax), pomxap(nncutmax)
#ifdef WG_TODO
    real (kind=R8) :: &
         tpmxip(nncutmax), tpmxap(nncutmax)
#endif
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

    integer i, j, is, iCv, iFc, ireg
    integer target_offset

    real(kind=R8) :: fettmp, fdir
    real(kind=R8), allocatable :: ptf(:,:), taf(:,:), uaf(:,:)
    integer, save :: write_2d = 0
    integer, save :: nc, ntstep, nastep, geometryType
    integer, allocatable :: fclist(:)
#ifndef NO_CDF
    integer, save :: ncid, nbatch
    integer, allocatable :: cvlist(:), cnlist(:)
    integer imap(maxvdims), iret, iCv1, iCv2, iFt, iatm, cvtrg
    integer nvars, natts, ndims, unlimid
    real (kind=R8) :: fac
    real (kind=R8), allocatable :: fcOr(:)
    real (kind=R8) :: &
         nesepi(nncutmax), tesepi(nncutmax), tisepi(nncutmax), &
         nesepm(nncutmax), tesepm(nncutmax), tisepm(nncutmax), &
         nesepa(nncutmax), tesepa(nncutmax), tisepa(nncutmax), &
         posepi(nncutmax), posepm(nncutmax), posepa(nncutmax), &
         dnsepm(nncutmax), dpsepm(nncutmax), kesepm(nncutmax), &
         kisepm(nncutmax), vxsepm(nncutmax), vysepm(nncutmax), &
         vssepm(nncutmax), tpsepi(nncutmax), tpsepa(nncutmax), &
         ktsepm(nncutmax), ktsepi(nncutmax), ktsepa(nncutmax)
#ifdef WG_TODO
    real (kind=R8) :: &
         tmhacore(1), tmhasol(1), tmhadiv(1)
#endif
    real (kind=R8) :: &
         timesa(1), batchsa(1), tstepn(1)
    real (kind=R8), allocatable :: slice(:)
    logical ex
    character*5 rw
    character*256, save :: filename, filename_av
    real(kind=R8) :: rratio
    external rratio
#endif
    !   ..procedures
    external xertst, ipgeti, batch_average
    !   ..initialisation
    save ncall, target_offset
    data ncall/0/, target_offset/1/

    !-----------------------------------------------------------------------
    !.computation

    ! ..preliminaries
    !   ..subprogram start-up calls
    call subini ('b2mwti')
    !     ..test input
    call xertst (0.lt.nCv, 'faulty argument nCv')
    call xertst (0.lt.nFc, 'faulty argument nFc')
    call xertst (1.le.ns, 'faulty argument ns')
    call xertst (0.le.ismain.and.ismain.lt.ns, &
         'invalid main plasma species index ismain')
    call xertst (is_neutral(ismain0).or.ismain0.eq.ismain, &
         'invalid main neutral species index ismain0')
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
      if (nimp.gt.0) call output_ds_cv(mpg,geo,nimp,imp,icsepimp-1,'dsi')
      if (nomp.gt.0) call output_ds_cv(mpg,geo,nomp,omp,icsepomp-1,'dsa')
      do i = 1, maxval(mpg%strDiv)
        allocate(fclist(mpg%divFcP(i,2)))
        fclist(1:mpg%divFcp(i,2)) = &
     &   mpg%divFc(mpg%divFcP(i,1):mpg%divFcP(i,1)+mpg%divFcP(i,2)-1)
        select case (i)
        case (1)
          call output_ds_fc(geo,mpg%divFcP(i,2),fclist,mpg%ifdiv(i),'dsl')
          open(99,file='dsL')
          do j = 1, mpg%divFcP(i,2)
            write(99,*) geo%fcS(fclist(j))
          enddo
          close(99)
          open(99,file='dsLP')
          do j = 1, mpg%divFcP(i,2)
            write(99,*) geo%fcS(fclist(j))*abs(geo%fcQalf(fclist(j),0))
          enddo
          close(99)
        case (2)
          if (mpg%nnreg(0).ge.7) then
            call output_ds_fc(geo,mpg%divFcP(i,2),fclist,mpg%ifdiv(i),'dstl')
            open(99,file='dsTL')
            do j = 1, mpg%divFcP(i,2)
              write(99,*) geo%fcS(fclist(j))
            enddo
            close(99)
            open(99,file='dsTLP')
            do j = 1, mpg%divFcP(i,2)
              write(99,*) geo%fcS(fclist(j))*abs(geo%fcQalf(fclist(j),0))
            enddo
            close(99)
          else
            call output_ds_fc(geo,mpg%divFcP(i,2),fclist,mpg%ifdiv(i),'dsr')
            open(99,file='dsR')
            do j = 1, mpg%divFcP(i,2)
              write(99,*) geo%fcS(fclist(j))
            enddo
            close(99)
            open(99,file='dsRP')
            do j = 1, mpg%divFcP(i,2)
              write(99,*) geo%fcS(fclist(j))*abs(geo%fcQalf(fclist(j),0))
            enddo
            close(99)
          end if
        case (3)
          call output_ds_fc(geo,mpg%divFcP(i,2),fclist,mpg%ifdiv(i),'dstr')
          open(99,file='dsTR')
          do j = 1, mpg%divFcP(i,2)
            write(99,*) geo%fcS(fclist(j))
          enddo
          close(99)
          open(99,file='dsTRP')
          do j = 1, mpg%divFcP(i,2)
            write(99,*) geo%fcS(fclist(j))*abs(geo%fcQalf(fclist(j),0))
          enddo
          close(99)
        case (4)
          call output_ds_fc(geo,mpg%divFcP(i,2),fclist,mpg%ifdiv(i),'dsr')
          open(99,file='dsR')
          do j = 1, mpg%divFcP(i,2)
            write(99,*) geo%fcS(fclist(j))
          enddo
          close(99)
          open(99,file='dsRP')
          do j = 1, mpg%divFcP(i,2)
            write(99,*) geo%fcS(fclist(j))*abs(geo%fcQalf(fclist(j),0))
          enddo
          close(99)
        end select
        deallocate(fclist)
      end do
#ifndef NO_CDF
      if (b2time.gt.0) then
        filename='b2time.nc'
        call find_file(filename,ex)
        if (.not.ex.or.switch%b2mndr_stim.ge.0.0_R8) then
          ntstep = 0
          write(6,'(a)') trim(filename)//' will be created'
          call b2crtimecdf(filename, mpg, ns, ismain, ismain0, nnmoli, &
            write_2d, ncid, .false., iret)
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
          call b2crtimecdf(filename, mpg, ns, ismain, ismain0, nnmoli, &
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
            call b2crtimecdf(filename_av, mpg, ns, ismain, ismain0, nnmoli, &
              write_2d, ncid, .true., iret)
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
              call b2crtimecdf(filename_av, mpg, ns, ismain, ismain0, nnmoli, &
               write_2d, ncid, .true., iret)
            endif
            iret = nf_open(trim(filename_av),or(NF_WRITE,NF_SHARE),ncid)
            call check_cdf_status(iret)
          else
            nastep = 0
            write(6,'(a)') trim(filename_av)//' will be replaced'
            call b2crtimecdf(filename_av, mpg, ns, ismain, ismain0, nnmoli, &
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

    if (ext%ns.gt.0) then
      allocate(ptf(1:mpg%nFc,0:ext%ns-1))
      allocate(taf(1:mpg%nFc,0:ext%ns-1))
      allocate(uaf(1:mpg%nFc,0:ext%ns-1))
      call intface1(mpg%nCv,mpg%nFc,ns-1,mpg%fcCv,geo%fcVol,ext%pt,ptf)
      call intface1(mpg%nCv,mpg%nFc,ns-1,mpg%fcCv,geo%fcVol,ext%ta,taf)
      call intface1(mpg%nCv,mpg%nFc,ns-1,mpg%fcCv,geo%fcVol,ext%ua,uaf)
    end if
    fnixip = 0.0_R8; feexip = 0.0_R8; feixip = 0.0_R8; fchxip = 0.0_R8; fetxip = 0.0_R8
    nemxip = 0.0_R8; temxip = 0.0_R8; timxip = 0.0_R8; pomxip = 0.0_R8; pwmxip = 0.0_R8
    do i = mpg%divFcP(1,1), mpg%divFcP(1,1) + mpg%divFcP(1,2) - 1
      iFc = mpg%divFc(i)
      if (mpg%fcCv(iFc,1).le.mpg%nCi) then
        iCv = mpg%fcCv(iFc,1)
      else
        iCv = mpg%fcCv(iFc,2)
      end if
      fnixip(1) = fnixip(1) + &
        &  mpg%divFcOr(i)*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
      feexip(1) = feexip(1) + &
        &  mpg%divFcOr(i)*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
      feixip(1) = feixip(1) + &
        &  mpg%divFcOr(i)*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
      fchxip(1) = fchxip(1) + &
        &  mpg%divFcOr(i)*(dv%fch(iFc,0) + dv%fch(iFc,1))
      fettmp = mpg%divFcOr(i)*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
        &                      dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
        &                      ext%fhi(iFc,0) + ext%fhi(iFc,1) )
      do is = 0, ext%ns-1
        fettmp = fettmp + mpg%divFcOr(i)*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
      end do
      fetxip(1) = fetxip(1) + fettmp
      nemxip(1) = max(nemxip(1), dv%ne(iCv) )
      temxip(1) = max(temxip(1), pl%te(iCv) )
      timxip(1) = max(timxip(1), pl%ti(iCv) )
      pomxip(1) = max(pomxip(1), pl%po(iCv) )
      pwmxip(1) = max(pwmxip(1), abs(fettmp)/geo%fcS(iFc) )
    end do
#ifdef WG_TODO
    tpmxip = 0.0_R8
    ix = -1 ! 1
    ix_off  = ix + target_offset
    do iy = iylstrt,iylend
      if (bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and. xymap(ix,iy).ne.0) then
        tpmxip(1) = max(tpmxip(1), target_temp(xymap(ix,iy),1))
      endif
    enddo
#endif

    fnixap = 0.0_R8; feexap = 0.0_R8; feixap = 0.0_R8; fchxap = 0.0_R8; fetxap = 0.0_R8
    nemxap = 0.0_R8; temxap = 0.0_R8; timxap = 0.0_R8; pomxap = 0.0_R8; pwmxap = 0.0_R8
    do i = mpg%divFcP(maxval(mpg%strDiv),1), &
         & mpg%divFcP(maxval(mpg%strDiv),1) + mpg%divFcP(maxval(mpg%strDiv),2) - 1
      iFc = mpg%divFc(i)
      if (mpg%fcCv(iFc,1).le.mpg%nCi) then
        iCv = mpg%fcCv(iFc,1)
      else
        iCv = mpg%fcCv(iFc,2)
      end if
      fnixap(1) = fnixap(1) + &
        &  mpg%divFcOr(i)*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
      feexap(1) = feexap(1) + &
        &  mpg%divFcOr(i)*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
      feixap(1) = feixap(1) + &
        &  mpg%divFcOr(i)*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
      fchxap(1) = fchxap(1) + &
        &  mpg%divFcOr(i)*(dv%fch(iFc,0) + dv%fch(iFc,1))
      fettmp = mpg%divFcOr(i)*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
        &                      dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
        &                      ext%fhi(iFc,0) + ext%fhi(iFc,1) )
      do is = 0, ext%ns-1
        fettmp = fettmp + mpg%divFcOr(i)*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
      end do
      fetxap(1) = fetxap(1) + fettmp
      nemxap(1) = max(nemxap(1), dv%ne(iCv) )
      temxap(1) = max(temxap(1), pl%te(iCv) )
      timxap(1) = max(timxap(1), pl%ti(iCv) )
      pomxap(1) = max(pomxap(1), pl%po(iCv) )
      pwmxap(1) = max(pwmxap(1), abs(fettmp)/geo%fcS(iFc) )
    end do
#ifdef WG_TODO
    tpmxap = 0.0_R8
    ix = nx ! 2
    ix_off  = ix - target_offset
    do iy = iyrstrt,iyrend
      if (bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and. xymap(ix,iy).ne.0) then
        tpmxap(1) = max(tpmxap(1), target_temp(xymap(ix,iy),1))
      endif
    enddo
#endif

    if(mpg%nXpt.ge.2.and.maxval(mpg%strDiv).gt.2) then
      do i = mpg%divFcP(2,1), mpg%divFcP(2,1) + mpg%divFcP(2,2) - 1
        iFc = mpg%divFc(i)
        if (mpg%fcCv(iFc,1).le.mpg%nCi) then
          iCv = mpg%fcCv(iFc,1)
        else
          iCv = mpg%fcCv(iFc,2)
        end if
        fnixip(2) = fnixip(2) + &
          &  mpg%divFcOr(i)*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feexip(2) = feexip(2) + &
          &  mpg%divFcOr(i)*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feixip(2) = feixip(2) + &
          &  mpg%divFcOr(i)*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchxip(2) = fchxip(2) + &
          &  mpg%divFcOr(i)*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = mpg%divFcOr(i)*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
          &                      dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
          &                      ext%fhi(iFc,0) + ext%fhi(iFc,1) )
        do is = 0, ext%ns-1
          fettmp = fettmp + mpg%divFcOr(i)*(ptf(iFc,is)*ev + &
            & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
            & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetxip(2) = fetxip(2) + fettmp
        nemxip(2) = max(nemxip(2), dv%ne(iCv) )
        temxip(2) = max(temxip(2), pl%te(iCv) )
        timxip(2) = max(timxip(2), pl%ti(iCv) )
        pomxip(2) = max(pomxip(2), pl%po(iCv) )
        pwmxip(2) = max(pwmxip(2), abs(fettmp)/geo%fcS(iFc) )
      end do
      do i = mpg%divFcP(3,1), mpg%divFcP(3,1) + mpg%divFcP(3,2) - 1
        iFc = mpg%divFc(i)
        if (mpg%fcCv(iFc,1).le.mpg%nCi) then
          iCv = mpg%fcCv(iFc,1)
        else
          iCv = mpg%fcCv(iFc,2)
        end if
        fnixap(2) = fnixap(2) + &
          &  mpg%divFcOr(i)*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feexap(2) = feexap(2) + &
          &  mpg%divFcOr(i)*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feixap(2) = feixap(2) + &
          &  mpg%divFcOr(i)*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchxap(2) = fchxap(2) + &
          &  mpg%divFcOr(i)*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = mpg%divFcOr(i)*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
          &                      dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
          &                      ext%fhi(iFc,0) + ext%fhi(iFc,1) )
        do is = 0, ext%ns-1
          fettmp = fettmp + mpg%divFcOr(i)*(ptf(iFc,is)*ev + &
            & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
            & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetxap(2) = fetxap(2) + fettmp
        nemxap(2) = max(nemxap(2), dv%ne(iCv) )
        temxap(2) = max(temxap(2), pl%te(iCv) )
        timxap(2) = max(timxap(2), pl%ti(iCv) )
        pomxap(2) = max(pomxap(2), pl%po(iCv) )
        pwmxap(2) = max(pwmxap(2), abs(fettmp)/geo%fcS(iFc) )
      end do
#ifdef WG_TODO
      ix = ixtr ! 3
      ix_off  = ix + target_offset
      do iy = iytrstrt,iytrend
        if (bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and. xymap(ix,iy).ne.0) then
          tpmxap(2) = max(tpmxap(2), target_temp(xymap(ix,iy),1))
        endif
      enddo

      ix = ixtl ! 4
      ix_off  = ix - target_offset
      do iy = iytlstrt,iytlend
        if (bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and. xymap(ix,iy).ne.0) then
          tpmxip(2) = max(tpmxip(2), target_temp(xymap(ix,iy),1))
        endif
      enddo
#endif
    endif

    fnisip = 0.0_R8; feesip = 0.0_R8; feisip = 0.0_R8; fchsip = 0.0_R8; fetsip = 0.0_R8
    fnisap = 0.0_R8; feesap = 0.0_R8; feisap = 0.0_R8; fchsap = 0.0_R8; fetsap = 0.0_R8
    fnisipp = 0.0_R8; feesipp = 0.0_R8; feisipp = 0.0_R8; fetsipp = 0.0_R8; fchsipp = 0.0_R8
    fnisapp = 0.0_R8; feesapp = 0.0_R8; feisapp = 0.0_R8; fetsapp = 0.0_R8; fchsapp = 0.0_R8
    fniyip = 0.0_R8; feeyip = 0.0_R8; feiyip = 0.0_R8; fetyip = 0.0_R8; fchyip = 0.0_R8
    fniyap = 0.0_R8; feeyap = 0.0_R8; feiyap = 0.0_R8; fetyap = 0.0_R8; fchyap = 0.0_R8
    do iFc = 1, mpg%nFc
      ireg = 0
      if(mpg%nnreg(0).eq.2) then
        ireg = 3
      else if (mpg%nnreg(0).eq.4 .or. mpg%nnreg(0).eq.5) then
        ireg = 5
      else if (mpg%nnreg(0).ge.7) then
        ireg = 9
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fnisipp(1) = fnisipp(1) + dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain)
        feesipp(1) = feesipp(1) + dv%fhe(iFc,0) + dv%fhe(iFc,1)
        feisipp(1) = feisipp(1) + dv%fhi(iFc,0) + dv%fhi(iFc,1)
        fchsipp(1) = fchsipp(1) + dv%fch(iFc,0) + dv%fch(iFc,1)
        fettmp = dv%fht(iFc,0) + dv%fht(iFc,1) - &
               & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
               & ext%fhi(iFc,0) + ext%fhi(iFc,1)
        do is = 0, ext%ns-1
          fettmp = fettmp + (ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetsipp(1) = fetsipp(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.8.and.mpg%nXpt.ge.2) ireg = 11
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fnisipp(2) = fnisipp(2) + dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain)
        feesipp(2) = feesipp(2) + dv%fhe(iFc,0) + dv%fhe(iFc,1)
        feisipp(2) = feisipp(2) + dv%fhi(iFc,0) + dv%fhi(iFc,1)
        fchsipp(2) = fchsipp(2) + dv%fch(iFc,0) + dv%fch(iFc,1)
        fettmp = dv%fht(iFc,0) + dv%fht(iFc,1) - &
               & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
               & ext%fhi(iFc,0) + ext%fhi(iFc,1)
        do is = 0, ext%ns-1
          fettmp = fettmp + (ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetsipp(2) = fetsipp(2) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.4.or.mpg%nnreg(0).eq.5) then
        ireg = 6
      else if (mpg%nnreg(0).eq.7) then
        ireg = 10
      else if (mpg%nnreg(0).eq.8) then
        ireg = 12
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fnisapp(1) = fnisapp(1) + dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain)
        feesapp(1) = feesapp(1) + dv%fhe(iFc,0) + dv%fhe(iFc,1)
        feisapp(1) = feisapp(1) + dv%fhi(iFc,0) + dv%fhi(iFc,1)
        fchsapp(1) = fchsapp(1) + dv%fch(iFc,0) + dv%fch(iFc,1)
        fettmp = dv%fht(iFc,0) + dv%fht(iFc,1) - &
               & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
               & ext%fhi(iFc,0) + ext%fhi(iFc,1)
        do is = 0, ext%ns-1
          fettmp = fettmp + (ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetsapp(1) = fetsapp(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).ge.7.and.mpg%nXpt.ge.2) ireg = 12
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fnisapp(2) = fnisapp(2) + dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain)
        feesapp(2) = feesapp(2) + dv%fhe(iFc,0) + dv%fhe(iFc,1)
        feisapp(2) = feisapp(2) + dv%fhi(iFc,0) + dv%fhi(iFc,1)
        fchsapp(2) = fchsapp(2) + dv%fch(iFc,0) + dv%fch(iFc,1)
        fettmp = dv%fht(iFc,0) + dv%fht(iFc,1) - &
               & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
               & ext%fhi(iFc,0) + ext%fhi(iFc,1)
        do is = 0, ext%ns-1
          fettmp = fettmp + (ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetsapp(2) = fetsapp(2) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).ge.4) ireg = 2
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fnisip(1) = fnisip(1) + dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain)
        feesip(1) = feesip(1) + dv%fhe(iFc,0) + dv%fhe(iFc,1)
        feisip(1) = feisip(1) + dv%fhi(iFc,0) + dv%fhi(iFc,1)
        fchsip(1) = fchsip(1) + dv%fch(iFc,0) + dv%fch(iFc,1)
        fettmp = dv%fht(iFc,0) + dv%fht(iFc,1) - &
               & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
               & ext%fhi(iFc,0) + ext%fhi(iFc,1)
        do is = 0, ext%ns-1
          fettmp = fettmp + (ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetsip(1) = fetsip(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).ge.7.and.mpg%nXpt.ge.2) ireg = 7
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fnisip(2) = fnisip(2) + dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain)
        feesip(2) = feesip(2) + dv%fhe(iFc,0) + dv%fhe(iFc,1)
        feisip(2) = feisip(2) + dv%fhi(iFc,0) + dv%fhi(iFc,1)
        fchsip(2) = fchsip(2) + dv%fch(iFc,0) + dv%fch(iFc,1)
        fettmp = dv%fht(iFc,0) + dv%fht(iFc,1) - &
               & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
               & ext%fhi(iFc,0) + ext%fhi(iFc,1)
        do is = 0, ext%ns-1
          fettmp = fettmp + (ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetsip(2) = fetsip(2) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).ge.4) ireg = 3
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fnisap(1) = fnisap(1) + dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain)
        feesap(1) = feesap(1) + dv%fhe(iFc,0) + dv%fhe(iFc,1)
        feisap(1) = feisap(1) + dv%fhi(iFc,0) + dv%fhi(iFc,1)
        fchsap(1) = fchsap(1) + dv%fch(iFc,0) + dv%fch(iFc,1)
        fettmp = dv%fht(iFc,0) + dv%fht(iFc,1) - &
               & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
               & ext%fhi(iFc,0) + ext%fhi(iFc,1)
        do is = 0, ext%ns-1
          fettmp = fettmp + (ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetsap(1) = fetsap(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).ge.7.and.mpg%nXpt.ge.2) ireg = 6
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fnisap(2) = fnisap(2) + dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain)
        feesap(2) = feesap(2) + dv%fhe(iFc,0) + dv%fhe(iFc,1)
        feisap(2) = feisap(2) + dv%fhi(iFc,0) + dv%fhi(iFc,1)
        fchsap(2) = fchsap(2) + dv%fch(iFc,0) + dv%fch(iFc,1)
        fettmp = dv%fht(iFc,0) + dv%fht(iFc,1) - &
               & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
               & ext%fhi(iFc,0) + ext%fhi(iFc,1)
        do is = 0, ext%ns-1
          fettmp = fettmp + (ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetsap(2) = fetsap(2) + fettmp
      end if
      do i = 1, size(mpg%bcFc)
        if (mpg%bcFc(i).eq.iFc) fdir = mpg%bcFcOr(i)
      end do
      ireg = 0
      if (mpg%nnreg(0).eq.1) then
        ireg = 4
      else if (mpg%nnreg(0).eq.2) then
        ireg = 6
      else if (mpg%nnreg(0).eq.4) then
        ireg = 12
      else if (mpg%nnreg(0).eq.5) then
        ireg = 13
      else if (mpg%nnreg(0).eq.7) then
        ireg = 19
      else if (mpg%nnreg(0).eq.8) then
        if (geometryType.eq.GEOMETRY_CDN) then
          ireg = 18
        else
          ireg = 19
        end if
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyip(1) = fniyip(1) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyip(1) = feeyip(1) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyip(1) = feiyip(1) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyip(1) = fchyip(1) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyip(1) = fetyip(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.4) then
        ireg = 11
      else if (mpg%nnreg(0).eq.5) then
        ireg = 8
      else if (mpg%nnreg(0).eq.7) then
        ireg = 18
      else if (mpg%nnreg(0).eq.8) then
        if (geometryType.eq.GEOMETRY_CDN) then
          ireg = 17
        else
          ireg = 18
        end if
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(1) = fniyap(1) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(1) = feeyap(1) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(1) = feiyap(1) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(1) = fchyap(1) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(1) = fetyap(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.4) then
        ireg = 13
      else if (mpg%nnreg(0).eq.5) then
        ireg = 10
      else if (mpg%nnreg(0).eq.7) then
        ireg = 20
      else if (mpg%nnreg(0).eq.8) then
        if (geometryType.eq.GEOMETRY_CDN) then
          ireg = 26
        else
          ireg = 27
        end if
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(1) = fniyap(1) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(1) = feeyap(1) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(1) = feiyap(1) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(1) = fchyap(1) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(1) = fetyap(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.7) then
        ireg = 24
      else if (mpg%nnreg(0).eq.8) then
        if (geometryType.eq.GEOMETRY_CDN) then
          ireg = 13
        else
          ireg = 14
        end if
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(1) = fniyap(1) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(1) = feeyap(1) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(1) = feiyap(1) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(1) = fchyap(1) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(1) = fetyap(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.7) then
        ireg = 25
      else if (mpg%nnreg(0).eq.8) then
        if (geometryType.eq.GEOMETRY_CDN) then
          ireg = 22
        else
          ireg = 23
        end if
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(1) = fniyap(1) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(1) = feeyap(1) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(1) = feiyap(1) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(1) = fchyap(1) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(1) = fetyap(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.7) ireg = 26
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(1) = fniyap(1) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(1) = feeyap(1) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(1) = feiyap(1) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(1) = fchyap(1) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(1) = fetyap(1) + fettmp
      end if
      ireg = 0
      ireg = 0
      if (mpg%nnreg(0).eq.4) then
        ireg = 7
      else if (mpg%nnreg(0).eq.7) then
        ireg = 14
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(1) = fniyap(1) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(1) = feeyap(1) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(1) = feiyap(1) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(1) = fchyap(1) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(1) = fetyap(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.4) then
        ireg = 9
      else if (mpg%nnreg(0).eq.7) then
        ireg = 16
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(1) = fniyap(1) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(1) = feeyap(1) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(1) = feiyap(1) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(1) = fchyap(1) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(1) = fetyap(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.7) ireg = 21
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(1) = fniyap(1) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(1) = feeyap(1) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(1) = feiyap(1) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(1) = fchyap(1) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(1) = fetyap(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.7) ireg = 22
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(1) = fniyap(1) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(1) = feeyap(1) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(1) = feiyap(1) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(1) = fchyap(1) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(1) = fetyap(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.7) ireg = 23
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(1) = fniyap(1) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(1) = feeyap(1) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(1) = feiyap(1) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(1) = fchyap(1) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(1) = fetyap(1) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.8) then
        if (mpg%nXpt.eq.1) then
          ireg = 0
        elseif (geometryType.eq.GEOMETRY_CDN) then
          ireg = 25
        else
          ireg = 26
        end if
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyip(2) = fniyip(2) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyip(2) = feeyip(2) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyip(2) = feiyip(2) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyip(2) = fchyip(2) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyip(2) = fetyip(2) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.8) then
        if (mpg%nXpt.eq.1) then
          ireg = 0
        elseif (geometryType.eq.GEOMETRY_CDN) then
          ireg = 19
        else
          ireg = 20
        end if
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(2) = fniyap(2) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(2) = feeyap(2) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(2) = feiyap(2) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(2) = fchyap(2) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(2) = fetyap(2) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.8) then
        if (mpg%nXpt.eq.1) then
          ireg = 0
        elseif (geometryType.eq.GEOMETRY_CDN) then
          ireg = 24
        else
          ireg = 25
        end if
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(2) = fniyap(2) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(2) = feeyap(2) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(2) = feiyap(2) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(2) = fchyap(2) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(2) = fetyap(2) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.8) then
        if (mpg%nXpt.eq.1) then
          ireg = 0
        elseif (geometryType.eq.GEOMETRY_CDN) then
          ireg = 15
        else
          ireg = 16
        end if
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(2) = fniyap(2) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(2) = feeyap(2) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(2) = feiyap(2) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(2) = fchyap(2) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(2) = fetyap(2) + fettmp
      end if
      ireg = 0
      if (mpg%nnreg(0).eq.8) then
        if (mpg%nXpt.eq.1) then
          ireg = 0
        elseif (geometryType.eq.GEOMETRY_CDN) then
          ireg = 20
        else
          ireg = 21
        end if
      end if
      if (mpg%fcReg(iFc).eq.ireg.and.ireg.ne.0) then
        fniyap(2) = fniyap(2) + &
          & fdir*(dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))
        feeyap(2) = feeyap(2) + fdir*(dv%fhe(iFc,0) + dv%fhe(iFc,1))
        feiyap(2) = feiyap(2) + fdir*(dv%fhi(iFc,0) + dv%fhi(iFc,1))
        fchyap(2) = fchyap(2) + fdir*(dv%fch(iFc,0) + dv%fch(iFc,1))
        fettmp = fdir*(dv%fht(iFc,0) + dv%fht(iFc,1) - &
                     & dv%fhj(iFc,0) - dv%fhj(iFc,1) + &
                     & ext%fhi(iFc,0) + ext%fhi(iFc,1))
        do is = 0, ext%ns-1
          fettmp = fettmp + fdir*(ptf(iFc,is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(iFc,is)**2+taf(iFc,is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(iFc,0,is)+ext%fa(iFc,1,is))
        end do
        fetyap(2) = fetyap(2) + fettmp
      end if
      ireg = 0
    end do
    if (ext%ns.gt.0) deallocate(ptf,taf,uaf)
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
    if (nimp.gt.0.and.mpg%nXpt.ge.2) then
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
      else if (geometryType.eq.GEOMETRY_DDN_TOP) then
        iFc = mpg%divFc(mpg%divFcP(1,1)+mpg%ifdiv(1)-1)
        if ((mpg%fcCv(iFc,1).gt.mpg%nCi.and.target_offset.eq.0).or. &
          & (mpg%fcCv(iFc,1).le.mpg%nCi.and.target_offset.eq.1)) then
          cvtrg = mpg%fcCv(iFc,1)
        else
          cvtrg = mpg%fcCv(iFc,2)
        end if
        nesepi(1) = dv%ne(cvtrg)
        tesepi(1) = pl%te(cvtrg)/ev
        tisepi(1) = pl%ti(cvtrg)/ev
        posepi(1) = pl%po(cvtrg)
        ktsepi(1) = pl%kt(cvtrg)
        iFc = mpg%divFc(mpg%divFcP(4,1)+mpg%ifdiv(4)-1)
        if ((mpg%fcCv(iFc,1).gt.mpg%nCi.and.target_offset.eq.0).or. &
          & (mpg%fcCv(iFc,1).le.mpg%nCi.and.target_offset.eq.1)) then
          cvtrg = mpg%fcCv(iFc,1)
        else
          cvtrg = mpg%fcCv(iFc,2)
        end if
        nesepa(1) = dv%ne(cvtrg)
        tesepa(1) = pl%te(cvtrg)/ev
        tisepa(1) = pl%ti(cvtrg)/ev
        posepa(1) = pl%po(cvtrg)
        ktsepa(1) = pl%kt(cvtrg)
      else if (geometryType.eq.GEOMETRY_DDN_BOTTOM) then
        if (maxval(mpg%strDiv).gt.2) then
          iFc = mpg%divFc(mpg%divFcP(2,1)+mpg%ifdiv(2)-1)
          if ((mpg%fcCv(iFc,1).gt.mpg%nCi.and.target_offset.eq.0).or. &
            & (mpg%fcCv(iFc,1).le.mpg%nCi.and.target_offset.eq.1)) then
            cvtrg = mpg%fcCv(iFc,1)
          else
            cvtrg = mpg%fcCv(iFc,2)
          end if
          nesepi(2) = dv%ne(cvtrg)
          tesepi(2) = pl%te(cvtrg)/ev
          tisepi(2) = pl%ti(cvtrg)/ev
          posepi(2) = pl%po(cvtrg)
          ktsepi(2) = pl%kt(cvtrg)
          iFc = mpg%divFc(mpg%divFcP(3,1)+mpg%ifdiv(3)-1)
          if ((mpg%fcCv(iFc,1).gt.mpg%nCi.and.target_offset.eq.0).or. &
            & (mpg%fcCv(iFc,1).le.mpg%nCi.and.target_offset.eq.1)) then
            cvtrg = mpg%fcCv(iFc,1)
          else
            cvtrg = mpg%fcCv(iFc,2)
          end if
          nesepa(2) = dv%ne(cvtrg)
          tesepa(2) = pl%te(cvtrg)/ev
          tisepa(2) = pl%ti(cvtrg)/ev
          posepa(2) = pl%po(cvtrg)
          ktsepa(2) = pl%kt(cvtrg)
        end if
      end if
    end if
#endif

#ifndef NO_CDF
#ifdef WG_TODO
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
    if(mpg%nXpt.ge.2) then
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
!WG_TODO
#endif
!NO_CDF
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
#ifdef WG_TODO
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
#ifdef WG_TODO
      call rwcdf(rw,ncid,'tmhacore',imap,tmhacore,iret)
      call rwcdf(rw,ncid,'tmhasol',imap,tmhasol,iret)
      call rwcdf(rw,ncid,'tmhadiv',imap,tmhadiv,iret)
#endif

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
      if (maxval(mpg%strDiv).ge.1) then
        imap(1)=1     ! bl
        imap(2)=1
        allocate(slice(mpg%divFcP(1,2)))
        allocate(fclist(mpg%divFcP(1,2)))
        allocate(cvlist(mpg%divFcP(1,2)))
        allocate(cnlist(mpg%divFcP(1,2)))
        allocate(fcOr(mpg%divFcP(1,2)))
        fclist(1:mpg%divFcp(1,2)) = &
     &   mpg%divFc(mpg%divFcP(1,1):mpg%divFcP(1,1)+mpg%divFcP(1,2)-1)
        fcOr(1:mpg%divFcp(1,2)) = &
     &   mpg%divFcOr(mpg%divFcP(1,1):mpg%divFcP(1,1)+mpg%divFcP(1,2)-1)
        do i = mpg%divFcP(1,1), mpg%divFcP(1,1)+mpg%divFcP(1,2)-1
          j = i - mpg%divFcP(1,1) + 1
          iFc = mpg%divFc(i)
          iCv1 = mpg%fcCv(iFc,1)
          iCv2 = mpg%fcCv(iFc,2)
          if ((iCv1.le.mpg%nCi.and.target_offset.eq.1).or. &
            & (iCv1.gt.mpg%nCi.and.target_offset.eq.0)) then
            cvlist(j) = iCv1
          else
            cvlist(j) = iCv2
          end if
          if (iCv1.le.mpg%nCi) cnlist(j) = iCv1
          if (iCv2.le.mpg%nCi) cnlist(j) = iCv2
        end do
        slice = 0.0_R8
        slice(:) = dv%ne(cvlist(:))
        call rwcdf(rw,ncid,'ne3dl',imap,slice,iret)
        slice(:) = pl%te(cvlist(:))
        call rwcdf(rw,ncid,'te3dl',imap,slice,iret)
        slice(:) = pl%ti(cvlist(:))
        call rwcdf(rw,ncid,'ti3dl',imap,slice,iret)
        slice(:) = pl%po(cvlist(:))
        call rwcdf(rw,ncid,'po3dl',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fna(fclist(:),0,ismain) + &
                 &          dv%fna(fclist(:),1,ismain))
        call rwcdf(rw,ncid,'fn3dl',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fhe(fclist(:),0) + dv%fhe(fclist(:),1))
        call rwcdf(rw,ncid,'fe3dl',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fhi(fclist(:),0) + dv%fhi(fclist(:),1))
        call rwcdf(rw,ncid,'fi3dl',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fch(fclist(:),0) + dv%fch(fclist(:),1))
        call rwcdf(rw,ncid,'fc3dl',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fne(fclist(:),0) + dv%fne(fclist(:),1))
        call rwcdf(rw,ncid,'fl3dl',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fni(fclist(:),0) + dv%fni(fclist(:),1))
        call rwcdf(rw,ncid,'fo3dl',imap,slice,iret)
        if (ismain0.ne.ismain) then
          slice(:) = pl%na(cvlist(:),ismain0)
          do iatm = 1, nnatmi
            if (b2espcr(ismain0).eq.latmscl(iatm)) &
              & slice(:) = slice(:) + dab2(cnlist(:),iatm,1)
          end do
          call rwcdf(rw,ncid,'an3dl',imap,slice,iret)
        end if
        if (nnmoli.gt.0) then
          slice(:) = dmb2(cnlist(:),1,1)
          call rwcdf(rw,ncid,'mn3dl',imap,slice,iret)
        end if
        slice(:) = fcOr(:)*(dv%fht(fclist(:),0) + dv%fht(fclist(:),1) - &
                          & dv%fhj(fclist(:),0) - dv%fhj(fclist(:),1) + &
                          & ext%fhi(fclist(:),0) + ext%fhi(fclist(:),1))
        do is = 0, ext%ns-1
          slice(:) = slice(:) + fcOr(:)*(ptf(fclist(:),is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(fclist(:),is)**2+ &
          &  taf(fclist(:),is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(fclist(:),0,is)+ext%fa(fclist(:),1,is))
        end do
        call rwcdf(rw,ncid,'ft3dl',imap,slice,iret)
        deallocate(fclist,cvlist,cnlist,fcOr)
        deallocate(slice)
      end if
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
      if (maxval(mpg%strDiv).ge.2) then
        imap(1)=1     ! br
        imap(2)=1
        ireg=maxval(mpg%strDiv)
        allocate(slice(mpg%divFcP(ireg,2)))
        allocate(fclist(mpg%divFcP(ireg,2)))
        allocate(cvlist(mpg%divFcP(ireg,2)))
        allocate(cnlist(mpg%divFcP(ireg,2)))
        allocate(fcOr(mpg%divFcP(ireg,2)))
        fclist(1:mpg%divFcp(ireg,2)) = &
     &   mpg%divFc(mpg%divFcP(ireg,1):mpg%divFcP(ireg,1)+mpg%divFcP(ireg,2)-1)
        fcOr(1:mpg%divFcp(ireg,2)) = &
     &   mpg%divFcOr(mpg%divFcP(ireg,1):mpg%divFcP(ireg,1)+mpg%divFcP(ireg,2)-1)
        do i = mpg%divFcP(ireg,1), mpg%divFcP(ireg,1)+mpg%divFcP(ireg,2)-1
          j = i - mpg%divFcP(ireg,1) + 1
          iFc = mpg%divFc(i)
          iCv1 = mpg%fcCv(iFc,1)
          iCv2 = mpg%fcCv(iFc,2)
          if ((iCv1.le.mpg%nCi.and.target_offset.eq.1).or. &
            & (iCv1.gt.mpg%nCi.and.target_offset.eq.0)) then
            cvlist(j) = iCv1
          else
            cvlist(j) = iCv2
          end if
          if (iCv1.le.mpg%nCi) cnlist(j) = iCv1
          if (iCv2.le.mpg%nCi) cnlist(j) = iCv2
        end do
        slice = 0.0_R8
        slice(:) = dv%ne(cvlist(:))
        call rwcdf(rw,ncid,'ne3dr',imap,slice,iret)
        slice(:) = pl%te(cvlist(:))
        call rwcdf(rw,ncid,'te3dr',imap,slice,iret)
        slice(:) = pl%ti(cvlist(:))
        call rwcdf(rw,ncid,'ti3dr',imap,slice,iret)
        slice(:) = pl%po(cvlist(:))
        call rwcdf(rw,ncid,'po3dr',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fna(fclist(:),0,ismain) + &
                          & dv%fna(fclist(:),1,ismain))
        call rwcdf(rw,ncid,'fn3dr',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fhe(fclist(:),0) + dv%fhe(fclist(:),1))
        call rwcdf(rw,ncid,'fe3dr',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fhi(fclist(:),0) + dv%fhi(fclist(:),1))
        call rwcdf(rw,ncid,'fi3dr',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fch(fclist(:),0) + dv%fch(fclist(:),1))
        call rwcdf(rw,ncid,'fc3dr',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fne(fclist(:),0) + dv%fne(fclist(:),1))
        call rwcdf(rw,ncid,'fl3dr',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fni(fclist(:),0) + dv%fni(fclist(:),1))
        call rwcdf(rw,ncid,'fo3dr',imap,slice,iret)
        if (ismain0.ne.ismain) then
          slice(:) = pl%na(cvlist(:),ismain0)
          do iatm = 1, nnatmi
            if (b2espcr(ismain0).eq.latmscl(iatm)) &
              & slice(:) = slice(:) + dab2(cnlist(:),iatm,1)
          end do
          call rwcdf(rw,ncid,'an3dr',imap,slice,iret)
        end if
        if (nnmoli.gt.0) then
          slice(:) = dmb2(cnlist(:),1,1)
          call rwcdf(rw,ncid,'mn3dr',imap,slice,iret)
        end if
        slice(:) = fcOr(:)*(dv%fht(fclist(:),0) + dv%fht(fclist(:),1) - &
                          & dv%fhj(fclist(:),0) - dv%fhj(fclist(:),1) + &
                          & ext%fhi(fclist(:),0) + ext%fhi(fclist(:),1))
        do is = 0, ext%ns-1
          slice(:) = slice(:) + fcOr(:)*(ptf(fclist(:),is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(fclist(:),is)**2+ &
          &  taf(fclist(:),is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(fclist(:),0,is)+ext%fa(fclist(:),1,is))
        end do
        call rwcdf(rw,ncid,'ft3dr',imap,slice,iret)
        deallocate(fclist,cvlist,cnlist,fcOr)
        deallocate(slice)
      end if
      if (maxval(mpg%strDiv).ge.3) then
        imap(1)=1     ! tr
        imap(2)=1
        allocate(slice(mpg%divFcP(3,2)))
        allocate(fclist(mpg%divFcP(3,2)))
        allocate(cvlist(mpg%divFcP(3,2)))
        allocate(cnlist(mpg%divFcP(3,2)))
        allocate(fcOr(mpg%divFcP(3,2)))
        fclist(1:mpg%divFcp(3,2)) = &
     &   mpg%divFc(mpg%divFcP(3,1):mpg%divFcP(3,1)+mpg%divFcP(3,2)-1)
        fcOr(1:mpg%divFcp(3,2)) = &
     &   mpg%divFcOr(mpg%divFcP(3,1):mpg%divFcP(3,1)+mpg%divFcP(3,2)-1)
        do i = mpg%divFcP(3,1), mpg%divFcP(3,1)+mpg%divFcP(3,2)-1
          j = i - mpg%divFcP(3,1) + 1
          iFc = mpg%divFc(i)
          iCv1 = mpg%fcCv(iFc,1)
          iCv2 = mpg%fcCv(iFc,2)
          if ((iCv1.le.mpg%nCi.and.target_offset.eq.1).or. &
            & (iCv1.gt.mpg%nCi.and.target_offset.eq.0)) then
            cvlist(j) = iCv1
          else
            cvlist(j) = iCv2
          end if
          if (iCv1.le.mpg%nCi) cnlist(j) = iCv1
          if (iCv2.le.mpg%nCi) cnlist(j) = iCv2
        end do
        slice = 0.0_R8
        slice(:) = dv%ne(cvlist(:))
        call rwcdf(rw,ncid,'ne3dtr',imap,slice,iret)
        slice(:) = pl%te(cvlist(:))
        call rwcdf(rw,ncid,'te3dtr',imap,slice,iret)
        slice(:) = pl%ti(cvlist(:))
        call rwcdf(rw,ncid,'ti3dtr',imap,slice,iret)
        slice(:) = pl%po(cvlist(:))
        call rwcdf(rw,ncid,'po3dtr',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fna(fclist(:),0,ismain) + &
                          & dv%fna(fclist(:),1,ismain))
        call rwcdf(rw,ncid,'fn3dtr',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fhe(fclist(:),0) + dv%fhe(fclist(:),1))
        call rwcdf(rw,ncid,'fe3dtr',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fhi(fclist(:),0) + dv%fhi(fclist(:),1))
        call rwcdf(rw,ncid,'fi3dtr',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fch(fclist(:),0) + dv%fch(fclist(:),1))
        call rwcdf(rw,ncid,'fc3dtr',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fne(fclist(:),0) + dv%fne(fclist(:),1))
        call rwcdf(rw,ncid,'fl3dtr',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fni(fclist(:),0) + dv%fni(fclist(:),1))
        call rwcdf(rw,ncid,'fo3dtr',imap,slice,iret)
        if (ismain0.ne.ismain) then
          slice(:) = pl%na(cvlist(:),ismain0)
          do iatm = 1, nnatmi
            if (b2espcr(ismain0).eq.latmscl(iatm)) &
              & slice(:) = slice(:) + dab2(cnlist(:),iatm,1)
          end do
          call rwcdf(rw,ncid,'an3dtr',imap,slice,iret)
        end if
        if (nnmoli.gt.0) then
          slice(:) = dmb2(cnlist(:),1,1)
          call rwcdf(rw,ncid,'mn3dtr',imap,slice,iret)
        end if
        slice(:) = fcOr(:)*(dv%fht(fclist(:),0) + dv%fht(fclist(:),1) - &
                          & dv%fhj(fclist(:),0) - dv%fhj(fclist(:),1) + &
                          & ext%fhi(fclist(:),0) + ext%fhi(fclist(:),1))
        do is = 0, ext%ns-1
          slice(:) = slice(:) + fcOr(:)*(ptf(fclist(:),is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(fclist(:),is)**2+ &
          &  taf(fclist(:),is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(fclist(:),0,is)+ext%fa(fclist(:),1,is))
        end do
        call rwcdf(rw,ncid,'ft3dtr',imap,slice,iret)
        deallocate(fclist,cvlist,cnlist,fcOr)
        deallocate(slice)
      end if
      if (maxval(mpg%strDiv).ge.4) then
        imap(1)=1     ! tl
        imap(2)=1
        allocate(slice(mpg%divFcP(2,2)))
        allocate(fclist(mpg%divFcP(2,2)))
        allocate(cvlist(mpg%divFcP(2,2)))
        allocate(cnlist(mpg%divFcP(2,2)))
        allocate(fcOr(mpg%divFcP(2,2)))
        fclist(1:mpg%divFcp(2,2)) = &
     &   mpg%divFc(mpg%divFcP(2,1):mpg%divFcP(2,1)+mpg%divFcP(2,2)-1)
        fcOr(1:mpg%divFcp(2,2)) = &
     &   mpg%divFcOr(mpg%divFcP(2,1):mpg%divFcP(2,1)+mpg%divFcP(2,2)-1)
        do i = mpg%divFcP(2,1), mpg%divFcP(2,1)+mpg%divFcP(2,2)-1
          j = i - mpg%divFcP(2,1) + 1
          iFc = mpg%divFc(i)
          iCv1 = mpg%fcCv(iFc,1)
          iCv2 = mpg%fcCv(iFc,2)
          if ((iCv1.le.mpg%nCi.and.target_offset.eq.1).or. &
            & (iCv1.gt.mpg%nCi.and.target_offset.eq.0)) then
            cvlist(j) = iCv1
          else
            cvlist(j) = iCv2
          end if
          if (iCv1.le.mpg%nCi) cnlist(j) = iCv1
          if (iCv2.le.mpg%nCi) cnlist(j) = iCv2
        end do
        slice = 0.0_R8
        slice(:) = dv%ne(cvlist(:))
        call rwcdf(rw,ncid,'ne3dtl',imap,slice,iret)
        slice(:) = pl%te(cvlist(:))
        call rwcdf(rw,ncid,'te3dtl',imap,slice,iret)
        slice(:) = pl%ti(cvlist(:))
        call rwcdf(rw,ncid,'ti3dtl',imap,slice,iret)
        slice(:) = pl%po(cvlist(:))
        call rwcdf(rw,ncid,'po3dtl',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fna(fclist(:),0,ismain) + &
                          & dv%fna(fclist(:),1,ismain))
        call rwcdf(rw,ncid,'fn3dtl',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fhe(fclist(:),0) + dv%fhe(fclist(:),1))
        call rwcdf(rw,ncid,'fe3dtl',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fhi(fclist(:),0) + dv%fhi(fclist(:),1))
        call rwcdf(rw,ncid,'fi3dtl',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fch(fclist(:),0) + dv%fch(fclist(:),1))
        call rwcdf(rw,ncid,'fc3dtl',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fne(fclist(:),0) + dv%fne(fclist(:),1))
        call rwcdf(rw,ncid,'fl3dtl',imap,slice,iret)
        slice(:) = fcOr(:)*(dv%fni(fclist(:),0) + dv%fni(fclist(:),1))
        call rwcdf(rw,ncid,'fo3dtl',imap,slice,iret)
        if (ismain0.ne.ismain) then
          slice(:) = pl%na(cvlist(:),ismain0)
          do iatm = 1, nnatmi
            if (b2espcr(ismain0).eq.latmscl(iatm)) &
              & slice(:) = slice(:) + dab2(cnlist(:),iatm,1)
          end do
          call rwcdf(rw,ncid,'an3dtl',imap,slice,iret)
        end if
        if (nnmoli.gt.0) then
          slice(:) = dmb2(cnlist(:),1,1)
          call rwcdf(rw,ncid,'mn3dtl',imap,slice,iret)
        end if
        slice(:) = fcOr(:)*(dv%fht(fclist(:),0) + dv%fht(fclist(:),1) - &
                          & dv%fhj(fclist(:),0) - dv%fhj(fclist(:),1) + &
                          & ext%fhi(fclist(:),0) + ext%fhi(fclist(:),1))
        do is = 0, ext%ns-1
          slice(:) = slice(:) + fcOr(:)*(ptf(fclist(:),is)*ev + &
          & (0.5_R8*ext%am(is)*mp*uaf(fclist(:),is)**2+ &
          &  taf(fclist(:),is))*(1.0_R8-switch%BoRiS))* &
          & (ext%fa(fclist(:),0,is)+ext%fa(fclist(:),1,is))
        end do
        call rwcdf(rw,ncid,'ft3dtl',imap,slice,iret)
        deallocate(fclist,cvlist,cnlist,fcOr)
        deallocate(slice)
      endif
    !
      imap(1)=1
      imap(2)=1
      if (nimp.gt.0) then
        if (ismain0.ne.ismain) then
          allocate(slice(nimp))
          slice(1:nimp)=pl%na(imp(1:nimp),ismain0)
          do iatm=1,nnatmi
            if (b2espcr(ismain0).eq.latmscl(iatm)) &
              &  slice(1:nimp)=slice(1:nimp)+dab2(imp(1:nimp),iatm,1)
          enddo
          call rwcdf(rw,ncid,'an3di',imap,slice,iret)
          deallocate(slice)
        endif
      end if
      if (nomp.gt.0) then
        if (ismain0.ne.ismain) then
          allocate(slice(nomp))
          slice(1:nomp)=pl%na(omp(1:nomp),ismain0)
          do iatm=1,nnatmi
            if (b2espcr(ismain0).eq.latmscl(iatm)) &
              &  slice(1:nomp)=slice(1:nomp)+dab2(omp(1:nomp),iatm,1)
          enddo
          call rwcdf(rw,ncid,'an3da',imap,slice,iret)
          deallocate(slice)
        endif
      end if
      if (nnmoli.gt.0) then
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
      endif
    !
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
#ifdef WG_TODO
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
  subroutine output_ds_cv(mpg,geo,nlist,cvlist,isep,filename)
    use b2us_geo_diff
    use b2us_map_diff
    implicit none
    type (mapping), intent(in) :: mpg
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
    if(.not.mpg%cvOnClosedSurface(cvlist(1))) then ! flip sign
      ds(1:nlist)=-ds(1:nlist)
    endif
    open(99,file=filename)
    do i=1,nlist
      write(99,*) ds(i)
    enddo
    close(99)
    return
  end subroutine output_ds_cv

  subroutine output_ds_fc(geo,nlist,fclist,isep,filename)
    use b2us_geo_diff
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
    use b2mod_plasma_diff   , only : ti, te, fna, fne, fhe, fhi, fch, fhm, fhp
    use b2mod_indirect_diff , only : rightix, rightiy, bottomix, bottomiy, topix, topiy, leftix, leftiy
    use b2mod_external_diff , only : fhi_ext, pt_ext, ta_ext, ua_ext, am_ext, ns_ext, fa_ext
    use b2mod_constants , only : ev, mp
    use b2mod_geo_diff , only : hx, hy, qz, gs
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
