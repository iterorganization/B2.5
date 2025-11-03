      module b2mod_facdrift_ExB
      use b2mod_types
      use b2mod_switches
      use b2mod_user_namelist
      implicit none
      real (kind=R8), allocatable ::
     & fac_ExB_profile(:), facdrift_profile(:), fac_vis_profile(:)
      real (kind=R8) ::
     & facdrift_scalar, fac_ExB_scalar, fac_vis_scalar,
     & facdrift_tanh_a, facdrift_tanh_b, facExB_tanh_a, facExB_tanh_b,
     1 facvis_tanh_a, facvis_tanh_b
      integer ncall_drift, iy_nocoreExB
      save
     3 facdrift_scalar, fac_ExB_scalar, ncall_drift, iy_nocoreExB,
     & facdrift_tanh_a, facdrift_tanh_b, facExB_tanh_a, facExB_tanh_b,
     & fac_vis_scalar, facvis_tanh_a, facvis_tanh_b
      data
     & facdrift_scalar/0.0_R8/, fac_ExB_scalar/0.0_R8/,
     1 fac_vis_scalar/0.0_R8/,
     & ncall_drift/0/, iy_nocoreExB/0/,
     & facdrift_tanh_a/0.0_R8/, facdrift_tanh_b/0.0_R8/,
     & facExB_tanh_a/0.0_R8/, facExB_tanh_b/0.0_R8/,
     & facvis_tanh_a/0.0_R8/, facvis_tanh_b/0.0_R8/

      contains

      subroutine b2news_facdrift_fac_ExB_init(nFc, mpg, fac_ExB,
     &  facdrift, fac_vis, sw)
      use b2us_map
      implicit none
      integer, intent (in) :: nFc
      type (mapping), intent (in) :: mpg
      real (kind=R8), intent (inout) :: fac_ExB(nFc), facdrift(nFc),
     & fac_vis(nFc)
      type (switches), intent (inout) :: sw
      integer :: i, ix, iy, iCv, iFc, iFs, iFt

      if(ncall_drift.eq.0) then
       facdrift_scalar=sw%facdrift_start
       if(sw%facdrift_start.eq.sw%facdrift_target .and.
     &    sw%facdrift_start.ne.0.0_R8) then
         write(*,*) 'facdrift set to a fixed value :',sw%facdrift_start
       elseif(sw%facdrift_start.ne.sw%facdrift_target.and.
     &        sw%ramp_slow.eq.0) then
         write(*,'(a,1p,1g9.2,a,1g9.2)')
     &    'facdrift ramped normally from ',sw%facdrift_start,
     &    ' towards ',sw%facdrift_target
       elseif(sw%facdrift_start.ne.sw%facdrift_target.and.
     &        sw%ramp_slow.ne.0) then
         write(*,'(a,1p,1g9.2,a,1g9.2)')
     &    'facdrift ramped slowly from ',sw%facdrift_start,
     &    ' towards ',sw%facdrift_target
       else
         write(*,*) 'No diamagnetic drifts'
       endif
       fac_ExB_scalar = sw%b2news_ExB
       if(fac_ExB_scalar.ne.0.0_R8) then
c old style
         sw%facExB_start=fac_ExB_scalar
         sw%facExB_target=fac_ExB_scalar
         sw%facExB_inc=1.0_R8
         sw%facExB_dec=1.0_R8
         write(*,*) 'fac_ExB set to a fixed value: ', fac_ExB_scalar
       elseif(sw%facExB_start.ne.0.0_R8.and.
     &        sw%facExB_start.eq.sw%facExB_target) then
         fac_ExB_scalar=sw%facExB_start
         sw%facExB_inc=1.0_R8
         sw%facExB_dec=1.0_R8
         write(*,*) 'fac_ExB set to a fixed value: ', fac_ExB_scalar
       elseif(sw%facExB_start.ne.0.0_R8.and.
     &        sw%facExB_start.ne.sw%facExB_target) then
         fac_ExB_scalar=sw%facExB_start
         if (sw%ramp_slow.eq.0) then
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_ExB ramped normally from ',sw%facExB_start,
     &      ' towards ',sw%facExB_target
         else
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_ExB ramped slowly from ',sw%facExB_start,
     &      ' towards ',sw%facExB_target
         endif
       else
         write(*,*) 'No ExB'
       endif
       fac_vis_scalar = sw%b2news_vis
       if(fac_vis_scalar.ne.0.0_R8) then
c old style
         sw%facvis_start=fac_vis_scalar
         sw%facvis_target=fac_vis_scalar
         sw%facvis_inc=1.0_R8
         sw%facvis_dec=1.0_R8
         write(*,*) 'fac_vis set to a fixed value: ', fac_vis_scalar
       elseif(sw%facvis_start.ne.0.0_R8.and.
     &        sw%facvis_start.eq.sw%facvis_target) then
         fac_vis_scalar=sw%facvis_start
         sw%facvis_inc=1.0_R8
         sw%facvis_dec=1.0_R8
         write(*,*) 'fac_vis set to a fixed value: ', fac_vis_scalar
       elseif(sw%facvis_start.ne.0.0_R8.and.
     &        sw%facvis_start.ne.sw%facvis_target) then
         fac_vis_scalar=sw%facvis_start
         if (sw%ramp_slow.eq.0) then
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_vis ramped normally from ',sw%facvis_start,
     &      ' towards ',sw%facvis_target
         else
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_vis ramped slowly from ',sw%facvis_start,
     &      ' towards ',sw%facvis_target
         endif
       else
         write(*,*) 'No vis'
       endif

       allocate(fac_ExB_profile(1:nFc))
       allocate(facdrift_profile(1:nFc))
       allocate(fac_vis_profile(1:nFc))
       fac_ExB_profile = 1.0_R8
       facdrift_profile = 1.0_R8
       fac_vis_profile = 1.0_R8
       call ipgeti ('b2news_iy_nocoreExB', iy_nocoreExB)                 !xpb
       call xertst (0.le.iy_nocoreExB .and. nomp.ge.iy_nocoreExB,        !xpb
     &  'faulty parameter iy_nocoreExB')                                 !xpb
       do iy = 1, iy_nocoreExB
         iFt = mpg%cvFt(omp(iy))
         do ix = 1, mpg%ftCvP(iFt,2)
           iCv = mpg%ftCv(mpg%ftCvP(iFt,1)+ix-1)
           if(mpg%cvOnClosedSurface(iCv)) then
             do i = 1, mpg%cvFcP(iCv,2)
               iFc = mpg%cvFc(mpg%cvFcP(iCv,1)+i-1)
               fac_ExB_profile(iFc) = 0.0_R8
             end do
           end if
         enddo
       enddo
       call ipgetr ('b2news_facdrift_tanh_a', facdrift_tanh_a)
       call ipgetr ('b2news_facdrift_tanh_b', facdrift_tanh_b)
       if(facdrift_tanh_b.ne.0.0_R8) then
         do iy = 1, nomp
           iFt = mpg%cvFt(omp(iy))
           do ix = 1, mpg%ftCvP(iFt,2)
             iCv = mpg%ftCv(mpg%ftCvP(iFt,1)+ix-1)
             if(mpg%cvOnClosedSurface(iCv)) then
               do i = 1, mpg%cvFcP(iCv,2)
                 iFc = mpg%cvFc(mpg%cvFcP(iCv,1)+i-1)
                 iFs = mpg%fcFs(iFc)
                 if (iFs.eq.0) then
                   facdrift_profile(iFc)=facdrift_profile(iFc)*0.5_R8*
     &              (1.0_R8+tanh((ft_ds_omp(iFt)-facdrift_tanh_a)/
     &                      facdrift_tanh_b))
                 else
                   facdrift_profile(iFc)=facdrift_profile(iFc)*0.5_R8*
     &              (1.0_R8+tanh((fs_ds_omp(iFs)-facdrift_tanh_a)/
     &                      facdrift_tanh_b))
                 end if
               end do
             endif
           enddo
         enddo
       endif
       call ipgetr ('b2news_fac_ExB_tanh_a', facExB_tanh_a)
       if (facExB_tanh_a.eq.0.0_R8)
     >  call ipgetr ('b2news_facExB_tanh_a', facExB_tanh_a)
       call ipgetr ('b2news_fac_ExB_tanh_b', facExB_tanh_b)
       if (facExB_tanh_b.eq.0.0_R8)
     >  call ipgetr ('b2news_facExB_tanh_b', facExB_tanh_b)
       if(facExB_tanh_b.ne.0.0_R8) then
         do iy = 1, nomp
           iFt = mpg%cvFt(omp(iy))
           do ix = 1, mpg%ftCvP(iFt,2)
             iCv = mpg%ftCv(mpg%ftCvP(iFt,1)+ix-1)
             if(mpg%cvOnClosedSurface(iCv)) then
               do i = 1, mpg%cvFcP(iCv,2)
                 iFc = mpg%cvFc(mpg%cvFcP(iCv,1)+i-1)
                 iFs = mpg%fcFs(iFc)
                 if (iFs.eq.0) then
                   fac_ExB_profile(iFc)=fac_ExB_profile(iFc)*0.5_R8*
     &              (1.0_R8+tanh((ft_ds_omp(iFt)-facExB_tanh_a)/
     &                      facExB_tanh_b))
                 else
                   fac_ExB_profile(iFc)=fac_ExB_profile(iFc)*0.5_R8*
     &              (1.0_R8+tanh((fs_ds_omp(iFs)-facExB_tanh_a)/
     &                      facExB_tanh_b))
                 end if
               end do
             endif
           enddo
         enddo
       endif
       call ipgetr ('b2news_fac_vis_tanh_a', facvis_tanh_a)
       if (facvis_tanh_a.eq.0.0_R8)
     >  call ipgetr ('b2news_facvis_tanh_a', facvis_tanh_a)
       call ipgetr ('b2news_fac_vis_tanh_b', facvis_tanh_b)
       if (facvis_tanh_b.eq.0.0_R8)
     >  call ipgetr ('b2news_facvis_tanh_b', facvis_tanh_b)
       if(facvis_tanh_b.ne.0.0_R8) then
         do iy = 1, nomp
           iFt = mpg%cvFt(omp(iy))
           do ix = 1, mpg%ftCvP(iFt,2)
             iCv = mpg%ftCv(mpg%ftCvP(iFt,1)+ix-1)
             if(mpg%cvOnClosedSurface(iCv)) then
               do i = 1, mpg%cvFcP(iCv,2)
                 iFc = mpg%cvFc(mpg%cvFcP(iCv,1)+i-1)
                 iFs = mpg%fcFs(iFc)
                 if (iFs.eq.0) then
                   fac_vis_profile(iFc)=fac_vis_profile(iFc)*0.5_R8*
     &              (1.0_R8+tanh((ft_ds_omp(iFt)-facvis_tanh_a)/
     &                      facvis_tanh_b))
                 else
                   fac_vis_profile(iFc)=fac_vis_profile(iFc)*0.5_R8*
     &              (1.0_R8+tanh((fs_ds_omp(iFs)-facvis_tanh_a)/
     &                      facvis_tanh_b))
                 end if
               end do
             endif
           enddo
         enddo
       endif
      endif

      call b2news_facdrift_fac_ExB_update(nFc,fac_ExB,facdrift,fac_vis)
      return
      end subroutine b2news_facdrift_fac_ExB_init

      subroutine b2news_facdrift_fac_ExB_inc(nFc, fac_ExB, facdrift,
     &  fac_vis, sw)
      implicit none
      integer nFc
      real (kind=R8), intent (inout) :: fac_ExB(nFc), facdrift(nFc),
     & fac_vis(nFc)
      type (switches), intent (in) :: sw

      call xertst (0.lt.ncall_drift,
     &     'b2news_facdrift_fac_ExB_inc called before initialisation')
      facdrift_scalar=min(facdrift_scalar*sw%facdrift_inc,
     &                    sw%facdrift_target)
      fac_ExB_scalar=min(fac_ExB_scalar*sw%facExB_inc,sw%facExB_target)
      fac_vis_scalar=min(fac_vis_scalar*sw%facvis_inc,sw%facvis_target)
      call b2news_facdrift_fac_ExB_update(nFc,fac_ExB,facdrift,fac_vis)

      return
      end subroutine b2news_facdrift_fac_ExB_inc

      subroutine b2news_facdrift_fac_ExB_dec(nFc, fac_ExB, facdrift,
     &  fac_vis, sw)
      implicit none
      integer nFc
      real (kind=R8), intent (inout) :: fac_ExB(nFc), facdrift(nFc),
     & fac_vis(nFc)
      type (switches), intent (in) :: sw

      call xertst (0.lt.ncall_drift,
     &     'b2news_facdrift_fac_ExB_dec called before initialisation')
      facdrift_scalar=max(facdrift_scalar*sw%facdrift_dec,
     &                    sw%facdrift_start)
      fac_ExB_scalar=max(fac_ExB_scalar*sw%facExB_dec,sw%facExB_start)
      fac_vis_scalar=max(fac_vis_scalar*sw%facvis_dec,sw%facvis_start)
      call b2news_facdrift_fac_ExB_update(nFc,fac_ExB,facdrift,fac_vis)

      return
      end subroutine b2news_facdrift_fac_ExB_dec

      subroutine b2news_facdrift_fac_ExB_update(nFc, fac_ExB, facdrift,
     &  fac_vis)
      implicit none
      integer nFc
      real (kind=R8), intent (inout) :: fac_ExB(nFc), facdrift(nFc),
     & fac_vis(nFc)

      facdrift=facdrift_scalar*facdrift_profile
      fac_ExB=fac_ExB_scalar*fac_ExB_profile
      fac_vis=fac_vis_scalar*fac_vis_profile
      write(*,'(1x,a,1p,6g12.5)') 'facdrift, fac_ExB, fac_vis',
     &   minval(facdrift),maxval(facdrift),
     &   minval(fac_ExB),maxval(fac_ExB),
     &   minval(fac_vis),maxval(fac_vis)
      ncall_drift=ncall_drift+1

      return
      end subroutine b2news_facdrift_fac_ExB_update

      subroutine dealloc_b2mod_facdrift_ExB
      implicit none

      if (allocated(facdrift_profile)) then
        deallocate(facdrift_profile)
        deallocate(fac_ExB_profile)
        deallocate(fac_vis_profile)
      end if

      return
      end subroutine dealloc_b2mod_facdrift_ExB

      end module b2mod_facdrift_ExB
