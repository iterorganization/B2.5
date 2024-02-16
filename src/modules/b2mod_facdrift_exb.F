      module b2mod_facdrift_ExB
      use b2mod_types
      implicit none
      real (kind=R8) ::
     1 facdrift_target, facdrift_start, facdrift_inc, facdrift_dec,
     3 facExB_target, facExB_start, facExB_inc, facExB_dec,
     3 facvis_target, facvis_start, facvis_inc, facvis_dec,
     & facdrift_scalar, fac_ExB_scalar, fac_vis_scalar,
     & facdrift_tanh_a, facdrift_tanh_b, facExB_tanh_a, facExB_tanh_b,
     1 facvis_tanh_a, facvis_tanh_b
      integer ncall, iy_nocoreExB, fac_ref
      save
     1 facdrift_target, facdrift_start, facdrift_inc, facdrift_dec,
     1 facExB_target, facExB_start, facExB_inc, facExB_dec,
     1 facvis_target, facvis_start, facvis_inc, facvis_dec,
     3 facdrift_scalar, fac_ExB_scalar, ncall, iy_nocoreExB,
     & facdrift_tanh_a, facdrift_tanh_b, facExB_tanh_a, facExB_tanh_b,
     & fac_ref, fac_vis_scalar, facvis_tanh_a, facvis_tanh_b
      data
     2 facdrift_target/0.0_R8/, facdrift_start/0.0_R8/,
     3 facdrift_inc/1.0_R8/, facdrift_dec/0.0_R8/,
     2 facExB_target/0.0_R8/, facExB_start/0.0_R8/,
     3 facExB_inc/1.0_R8/, facExB_dec/0.0_R8/,
     2 facvis_target/0.0_R8/, facvis_start/0.0_R8/,
     3 facvis_inc/1.0_R8/, facvis_dec/0.0_R8/,
     & facdrift_scalar/0.0_R8/, fac_ExB_scalar/0.0_R8/,
     1 fac_vis_scalar/0.0_R8/,
     & ncall/0/, iy_nocoreExB/-2/,
     & facdrift_tanh_a/0.0_R8/, facdrift_tanh_b/0.0_R8/,
     & facExB_tanh_a/0.0_R8/, facExB_tanh_b/0.0_R8/,
     & facvis_tanh_a/0.0_R8/, facvis_tanh_b/0.0_R8/
      contains

      subroutine b2news_facdrift_fac_ExB_init(nx, ny, ramp_slow)
      use b2mod_indirect
      use b2mod_geo
     & , only: cr, cz, crx, cry 
      use b2mod_work
      implicit none
      integer nx, ny, ramp_slow
      integer ix, iy, jxi, jxa, jsep
      real (kind=R8) :: dssep

      if(ncall.eq.0) then
       call get_jsep(nx,ny,jxi,jxa,jsep)
       fac_ref=jxa
       call ipgeti ('b2news_fac_ref',fac_ref)
       call xertst (-1.le.fac_ref.and.fac_ref.le.nx,
     2      'faulty internal parameter fac_ref')
       call ipgetr ('b2news_facdrift_start', facdrift_start)
       call ipgetr ('b2news_facdrift_target', facdrift_target)
       call ipgetr ('b2news_facdrift_inc', facdrift_inc)
       call ipgetr ('b2news_facdrift_dec', facdrift_dec)
       call xertst (0.0_R8.le.facdrift_start .and.                       !xpb
     &              1.0_R8.ge.facdrift_start,                            !xpb
     &   'faulty argument facdrift_start')                               !xpb
       call xertst (0.0_R8.le.facdrift_target .and.                      !xpb
     &              1.0_R8.ge.facdrift_target,                           !xpb
     &   'faulty argument facdrift_target')                              !xpb
       call xertst (1.0_R8.le.facdrift_inc,                              !xpb
     &   'faulty argument facdrift_inc')                                 !xpb
       call xertst (0.0_R8.le.facdrift_dec .and. 1.0_R8.ge.facdrift_dec,
     &   'faulty argument facdrift_dec')                                 !xpb
       facdrift_scalar=facdrift_start
       if(facdrift_start.eq.facdrift_target .and.
     &  facdrift_start.ne.0.0_R8) then
         write(*,*) 'facdrift set to a fixed value :', facdrift_start
       elseif(facdrift_start.ne.facdrift_target.and.ramp_slow.eq.0) then
         write(*,'(a,1p,1g9.2,a,1g9.2)')
     &    'facdrift ramped normally from ',facdrift_start,
     &    ' towards ',facdrift_target
       elseif(facdrift_start.ne.facdrift_target.and.ramp_slow.ne.0) then
         write(*,'(a,1p,1g9.2,a,1g9.2)')
     &    'facdrift ramped slowly from ',facdrift_start,
     &    ' towards ',facdrift_target
       else
         write(*,*) 'No diamagnetic drifts'
       endif
       call ipgetr ('b2news_facExB_start', facExB_start)
       call ipgetr ('b2news_facExB_target', facExB_target)
       call ipgetr ('b2news_facExB_inc', facExB_inc)
       call ipgetr ('b2news_facExB_dec', facExB_dec)
       call xertst (0.0_R8.le.facExB_start .and.
     &              1.0_R8.ge.facExB_start,
     &   'faulty argument facExB_start')
       call xertst (0.0_R8.le.facExB_target .and.
     &              1.0_R8.ge.facExB_target,
     &   'faulty argument facExB_target')
       call xertst (1.0_R8.le.facExB_inc,
     &   'faulty argument facExB_inc')
       call xertst (0.0_R8.le.facExB_dec .and. 1.0_R8.ge.facExB_dec,
     &   'faulty argument facExB_dec')
       call ipgetr ('b2news_ExB', fac_ExB_scalar)                        !xpb
       call xertst (0.0_R8.le.fac_ExB_scalar .and.
     &              1.0_R8.ge.fac_ExB_scalar,
     &  'faulty parameter b2news_ExB')                                   !xpb
       if(fac_ExB_scalar.ne.0.0_R8) then
c old style
         facExB_start=fac_ExB_scalar
         facExB_target=fac_ExB_scalar
         facExB_inc=1.0_R8
         facExB_dec=1.0_R8
         write(*,*) 'fac_ExB set to a fixed value: ', fac_ExB_scalar
       elseif(facExB_start.ne.0.0_R8.and.facExB_start.eq.facExB_target)
     >  then
         fac_ExB_scalar=facExB_start
         facExB_inc=1.0_R8
         facExB_dec=1.0_R8
         write(*,*) 'fac_ExB set to a fixed value: ', fac_ExB_scalar
       elseif(facExB_start.ne.0.0_R8.and.facExB_start.ne.facExB_target)
     >  then
         fac_ExB_scalar=facExB_start
         if (ramp_slow.eq.0) then
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_ExB ramped normally from ',facExB_start,
     &      ' towards ',facExB_target
         else
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_ExB ramped slowly from ',facExB_start,
     &      ' towards ',facExB_target
         endif
       else
         write(*,*) 'No ExB'
       endif
       call ipgeti ('b2news_iy_nocoreExB', iy_nocoreExB)                 !xpb
       call xertst (-2.le.iy_nocoreExB .and. ny.ge.iy_nocoreExB,         !xpb
     &   'faulty parameter iy_nocoreExB')                                !xpb
       do iy=-1, iy_nocoreExB
         do ix=-1,nx
           if(on_closed_surface(ix,iy)) then
             fac_ExB_profile(ix,iy)=0.0_R8
           endif
         enddo
       enddo
       call ipgetr ('b2news_facvis_start', facvis_start)
       call ipgetr ('b2news_facvis_target', facvis_target)
       call ipgetr ('b2news_facvis_inc', facvis_inc)
       call ipgetr ('b2news_facvis_dec', facvis_dec)
       call xertst (0.0_R8.le.facvis_start .and.
     &              1.0_R8.ge.facvis_start,
     &   'faulty argument facvis_start')
       call xertst (0.0_R8.le.facvis_target .and.
     &              1.0_R8.ge.facvis_target,
     &   'faulty argument facvis_target')
       call xertst (1.0_R8.le.facvis_inc,
     &   'faulty argument facvis_inc')
       call xertst (0.0_R8.le.facvis_dec .and. 1.0_R8.ge.facvis_dec,
     &   'faulty argument facvis_dec')
       call ipgetr ('b2news_vis', fac_vis_scalar)                        !xpb
       call xertst (0.0_R8.le.fac_vis_scalar .and.
     &              1.0_R8.ge.fac_vis_scalar,
     &  'faulty parameter b2news_vis')                                   !xpb
       if(fac_vis_scalar.ne.0.0_R8) then
c old style
         facvis_start=fac_vis_scalar
         facvis_target=fac_vis_scalar
         facvis_inc=1.0_R8
         facvis_dec=1.0_R8
         write(*,*) 'fac_vis set to a fixed value: ', fac_vis_scalar
       elseif(facvis_start.ne.0.0_R8.and.facvis_start.eq.facvis_target)
     >  then
         fac_vis_scalar=facvis_start
         facvis_inc=1.0_R8
         facvis_dec=1.0_R8
         write(*,*) 'fac_vis set to a fixed value: ', fac_vis_scalar
       elseif(facvis_start.ne.0.0_R8.and.facvis_start.ne.facvis_target)
     >  then
         fac_vis_scalar=facvis_start
         if (ramp_slow.eq.0) then
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_vis ramped normally from ',facvis_start,
     &      ' towards ',facvis_target
         else
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_vis ramped slowly from ',facvis_start,
     &      ' towards ',facvis_target
         endif
       else
         write(*,*) 'No vis'
       endif
       ds(-1)=
     1   sqrt((cr(fac_ref,-1)-
     1      0.5_R8*(crx(fac_ref,-1,0)+crx(fac_ref,-1,1)))**2+
     2      (cz(fac_ref,-1)-
     3      0.5_R8*(cry(fac_ref,-1,0)+cry(fac_ref,-1,1)))**2)
       do iy=0,ny
          ds(iy)=ds(iy-1)+
     1        sqrt((cr(fac_ref,iy)-cr(fac_ref,iy-1))**2+
     2        (cz(fac_ref,iy)-cz(fac_ref,iy-1))**2)
       enddo
       dssep=(ds(jsep)+ds(jsep+1))/2.0_R8
       do iy=-1,ny
         ds(iy)=ds(iy)-dssep
       enddo

       call ipgetr ('b2news_facdrift_tanh_a', facdrift_tanh_a)
       call ipgetr ('b2news_facdrift_tanh_b', facdrift_tanh_b)
       if(facdrift_tanh_b.ne.0.0_R8) then
          do iy=-1, ny
             do ix=-1, nx
                if(on_closed_surface(ix,iy)) then
                   facdrift_profile(ix,iy)=facdrift_profile(ix,iy)*
     &                  0.5_R8*(1.0_R8+
     &                    tanh((ds(iy)-facdrift_tanh_a)/
     &                   facdrift_tanh_b))
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
          do iy=-1, ny
             do ix=-1, nx
                if(on_closed_surface(ix,iy)) then
                   fac_ExB_profile(ix,iy)=fac_ExB_profile(ix,iy)*
     &                  0.5_R8*(1.0_R8+
     &                    tanh((ds(iy)-facExB_tanh_a)/
     &                   facExB_tanh_b))
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
          do iy=-1, ny
             do ix=-1, nx
                if(on_closed_surface(ix,iy)) then
                   fac_vis_profile(ix,iy)=fac_vis_profile(ix,iy)*
     &                  0.5_R8*(1.0_R8+
     &                    tanh((ds(iy)-facvis_tanh_a)/
     &                   facvis_tanh_b))
                endif
             enddo
          enddo
       endif
      endif

      call b2news_facdrift_fac_ExB_update
      return

      end subroutine b2news_facdrift_fac_ExB_init

      subroutine b2news_facdrift_fac_ExB_inc
      implicit none
      call xertst (0.lt.ncall,
     &     'b2news_facdrift_fac_ExB_inc called before initialisation')
      facdrift_scalar=min(facdrift_scalar*facdrift_inc,facdrift_target)
      fac_ExB_scalar=min(fac_ExB_scalar*facExB_inc,facExB_target)
      fac_vis_scalar=min(fac_vis_scalar*facvis_inc,facvis_target)
      call b2news_facdrift_fac_ExB_update
      return
      end subroutine b2news_facdrift_fac_ExB_inc

      subroutine b2news_facdrift_fac_ExB_dec
      implicit none
      call xertst (0.lt.ncall,
     &     'b2news_facdrift_fac_ExB_dec called before initialisation')
      facdrift_scalar=max(facdrift_scalar*facdrift_dec,facdrift_start)
      fac_ExB_scalar=max(fac_ExB_scalar*facExB_dec,facExB_start)
      fac_vis_scalar=max(fac_vis_scalar*facvis_dec,facvis_start)
      call b2news_facdrift_fac_ExB_update
      return
      end subroutine b2news_facdrift_fac_ExB_dec

      subroutine b2news_facdrift_fac_ExB_update
      use b2mod_plasma
     , , only: facdrift, fac_ExB, fac_vis
      use b2mod_work
     , , only: facdrift_profile, fac_ExB_profile, fac_vis_profile
      implicit none
      facdrift=facdrift_scalar*facdrift_profile
      fac_ExB=fac_ExB_scalar*fac_ExB_profile
      fac_vis=fac_vis_scalar*fac_vis_profile
      write(*,'(1x,a,1p,6g12.5)') 'facdrift, fac_ExB, fac_vis',
     &   minval(facdrift),maxval(facdrift),
     &   minval(fac_ExB),maxval(fac_ExB),
     &   minval(fac_vis),maxval(fac_vis)
      ncall=ncall+1
      return
      end subroutine b2news_facdrift_fac_ExB_update

      end module b2mod_facdrift_ExB
