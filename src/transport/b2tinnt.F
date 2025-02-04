*-----------------------------------------------------------------------
*.specification

      subroutine b2tinnt (nx, ny, ns, nscx, iscx,
     &       facdrift, fac_ExB,
     &       vol, gs, hx, hy, hy1, hz, qz, qc, qs, bb,                  !srv 08.07.09 { !srv 20.10.17 !som 04.11.21
     &       te, po, ti, na, fna, csigin, fchin_a, fchin)               !srv 08.03.99 16.06.09
      use b2mod_types
      use b2mod_indirect
      use b2mod_rates
      use b2mod_geo_corner
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_neutrals_namelist
#ifdef B25_EIRENE
      use eirmod_wneutrals
      use b2mod_sources
     ,   , only: smo0_eir_tot, smr0_eir_tot, smd0_eir_tot
#endif
      use b2mod_math
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nscx, iscx(0:nscxmax-1)                       !som 04.11.21
      real (kind=R8) ::
     & facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny),                     !srv 20.10.17
     * vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2),
     1 hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     1 hy1(-1:nx,-1:ny), hz(-1:nx,-1:ny),
     1 qz(-1:nx,-1:ny,0:1), bb(-1:nx,-1:ny,0:3),
     1 te(-1:nx,-1:ny), po(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     & na(-1:nx,-1:ny,0:ns-1), fna(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     & csigin(-1:nx,-1:ny,0:1,0:1,0:ns-1,0:ns-1), qc(-1:nx,-1:ny,0:1),
     & qs(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     & fchin(-1:nx,-1:ny,0:1,0:1), fchin_a(-1:nx,-1:ny,0:1,0:1,0:ns-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TINNT computes the ion-neutral current.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer, save :: pot_eq = 1                                       !srv 26.07.06
      integer, save :: mode = 1
      integer, save :: iout = 0
      integer, save :: use_eirene, area_fix
      integer, save :: fchin_in_core = 0                                !srv 16.11.17
      integer, save :: no_current = 0
      integer, save :: use_kinetic = 0                                  !IYS 17.10.2023
      integer n2, ix, iy, is, is0, ia, k                                !som 04.11.21
#ifdef B25_EIRENE
      integer im
#endif
      integer, save :: ncall = 0, style = 1
      integer, save :: csigin_style = 1
      integer, save :: use_9pt_stencil = 1
      real (kind=R8), save :: phm0 = 1.0_R8
      real (kind=R8), save :: na_min
      real (kind=R8), save :: fch_ion_neutral = 0.0_R8
      real (kind=R8) ::
     *   fnac(-1:nx,-1:ny,0:1,0:ns-1), vntrc(-1:nx,-1:ny,0:1),          !srv 08.07.09
     *   vntrl(-1:nx,-1:ny), vntrb(-1:nx,-1:ny),
     *   til(-1:nx,-1:ny), tib(-1:nx,-1:ny),
     *   fnan(-1:nx,-1:ny,0:1,0:1,0:ns-1), work(-1:nx,-1:ny,0:1,0:ns-1),
     *   fnal(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     *   fnav(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     *   nal(-1:nx,-1:ny,0:ns-1), nab(-1:nx,-1:ny,0:ns-1),
     *   wght(-1:nx,-1:ny), wrk0(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     *   wrk1(-1:nx,-1:ny,0:1,0:1), dpo(-1:nx,-1:ny,0:1,0:1),
     *   wrk(-1:nx,-1:ny)
      real (kind=R8) ::
     &   sec(-1:nx,-1:ny), sigia(-1:nx,-1:ny),
     &   csigia(-1:nx,-1:ny,0:1,0:1),                                   !srv 15.02.11
     &   sigi(-1:nx,-1:ny), facdriftm, fac_ExBm                         !srv 20.10.17
      real (kind=R8) :: mu, t0                                          !som 05.11.21
      character chns*3, chns0*3                                         !som 04.11.21
      logical found
*   ..procedures
      external ipgeti, xertst
*   ..initialisation
      data use_eirene/0/, area_fix/3/, na_min/1.0e4_R8/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tinnt')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgeti ('b2tinnt_iout', iout)
        call ipgeti ('b2tinnt_style', style)
        call ipgeti ('b2tral_mode', mode)
        call ipgeti ('b2mndr_eirene', use_eirene)
        call ipgeti ('b2news_area_fix', area_fix)
        call ipgeti ('b2news_poteq', pot_eq)                            !srv 26.07.06
        call ipgetr ('b2tfhe_neutral', fch_ion_neutral)
        call ipgeti ('b2tinnt_fchin_in_core', fchin_in_core)            !srv 16.11.17
        call ipgeti ('b2tfhe_no_current', no_current)
        call ipgeti ('b2tqin_csigin_style', csigin_style)
        call ipgetr ('b2sqcx_phm0', phm0)
        call ipgetr ('b2mndr_na_min', na_min)
        call xertst (na_min.gt.0.0_R8,
     &   'faulty internal parameter na_min')
        call ipgeti ('b2mndr_use_9pt_stencil', use_9pt_stencil)
        call ipgeti ('b2tinnt_kinetic', use_kinetic)
        if (use_kinetic.eq.1) then
#ifdef B25_EIRENE
          if (use_eirene.eq.0) then
            write(*,*) ' USE_KINETIC == 1 while USE_EIRENE == 0!!!'
            call xerrab
     & ('Impossible to use kinetic momentum sources for fluid neutrals')
          endif
#else
          call xerrab
     &     ('b2tinnt not compiled with the B25_EIRENE option')
#endif
          if (style.eq.0)
     &     call xerrab ('b2tinnt_style = 0 is obsolete '//
     &      'and not compatible with b2tinnt_kinetic = 1')
        endif
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of hz
        call b2xvsg (n2, hz, 1, 'hz', '.gt.')
      end if

      wght      = 1.0_R8
      fchin     = 0.0_R8                                                !srv 09.07.99
      fchin_a   = 0.0_R8                                                !IYS 17.10.2023
      facdriftm = maxval(facdrift)                                      !srv 20.10.17
      fac_ExBm  = maxval(fac_ExB)
      if(pot_eq.ne.0 .and.(facdriftm.gt.0.0_R8 .or. fac_ExBm.gt.0.0_R8)
     &               .and. no_current.eq.0
     .               .and. fch_ion_neutral .ne. 0.0_R8) then            !srv 26.07.06 01.07.99
       if(style.eq.0) then
*   ..interpolate ti in the left and bottom cell face
         call intfaceh(nx, ny, leftix, leftiy, vol, ti, til)
         call intfacev(nx, ny, bottomix, bottomiy, vol, ti, tib)
*   ..pre-computation of differences
         do is = 0, ns-1
           wrk = log(na(:,:,is))
           call diff(nx, ny, hx, hy, wght, qs, wrk,
     &          wrk0(-1,-1,0,0,is))
         end do
         call diff(nx, ny, hx, hy, wght, qs, ti, wrk1)
       else
         do is = 0, ns-1
*   ..interpolate na in the left and bottom cell face
           call intfaceh(nx, ny, leftix, leftiy, vol,
     *                 na(-1,-1,is), nal(-1,-1,is))                     !som 04.11.21
           call intfacev(nx, ny, bottomix, bottomiy, vol,
     &                 na(-1,-1,is), nab(-1,-1,is))                     !som 04.11.21
*   ..pre-computation of differences
           wrk = na(:,:,is)*ti(:,:)
           call diff(nx, ny, hx, hy, wght, qs, wrk,
     &                 wrk0(-1,-1,0,0,is))
         enddo
       endif
       call diff(nx, ny, hx, hy, wght, qs, po, dpo)
       if (use_eirene.ne.0) then
!        fchin = 0.0_R8

        if (use_kinetic.eq.0) then                                      !srv 15.02.11 {
*   ..calculate section
         do is0 = 0, ns-1                                               !som 04.11.21 {
          if (is_neutral(is0).and.nint(zn(is0)).eq.1) then
*   ..identify matching CX rate
           k = 0
           found = .false.
           do while (.not.found .and. k.lt.nscx)
            found = iscx(k).eq.is0
            if (.not.found) k = k + 1
           enddo
*   ..get Eirene neutral fluxes divided by area
           fnan = 0.0_R8
           do ia = 1, nnatmi
            if (b2espcr(is0).ne.latmscl(ia)) cycle
            do iy=-1,ny
             do ix=-1,nx
              if(leftix(ix,iy).ne.-2) fnan(ix,iy,0,0,is0)=
     .         fnan(ix,iy,0,0,is0)+pfluxa(ix+1,iy+1,ia,1) !/gs(ix,iy,0) !som 04.11.21
              if(bottomiy(ix,iy).ne.-2) fnan(ix,iy,1,1,is0)=
     .         fnan(ix,iy,1,1,is0)+rfluxa(ix+1,iy+1,ia,1) !/gs(ix,iy,1) !som 04.11.21
             enddo
            enddo
           end do
*   ..get fluxes in the cell centre
           fnac(:,:,0,is0) = fnan(:,:,0,0,is0)                          !som 04.11.21
           fnac(:,:,1,is0) = fnan(:,:,1,1,is0)                          !som 04.11.21
*   ..interpolate fluxes in the cell face
           call intfaceh(nx,ny,leftix,leftiy,vol,fnac(-1,-1,1,is0),
     &                                           fnal(-1,-1,0,0,is0))   !som 04.11.21
           call intfaceh(nx,ny,leftix,leftiy,vol,fnac(-1,-1,0,is0),
     &                                           fnal(-1,-1,0,1,is0))
           call intfacev(nx, ny, bottomix, bottomiy, vol,
     &                       fnac(-1,-1,0,is0), fnav(-1,-1,1,1,is0))
           do is = 0, ns-1
            if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then
             sigia = 0.0_R8
             mu=am(is)*am(is0)/(am(is)+am(is0))                         !som 05.11.21
             if (csigin_style.eq.0.and.found) then
! From 5.0/5.1: Using the CX rate if available
               do iy = -1, ny
                do ix = -1, nx
                  t0 = log(ti(ix,iy)/(am(is0)*ev))
                  sec(ix,iy) = phm0*sqrt(1.0_R8/mu)*
     &             expu(rlcx(ix,iy,0,is,k)+rlcx(ix,iy,1,is,k)*t0)
                end do
               end do
             else
! From 5.2: Using the diffusivity of room temperature hydrogen gas
!         sec = 3.2e-15_R8*sqrt(ti/ev/0.026_R8)
               sec(:,:) = 3.2e-15_R8*sqrt(ti(:,:)/ev/0.026_R8)*
     &          sqrt(1.0_R8/mu)
             end if
             sigi(:,:) = mu*mp*na(:,:,is)*sec(:,:)/(bb(:,:,3)**2)       !som 05.11.21
             do ia = 1, nnatmi
              if (b2espcr(is0).ne.latmscl(ia)) cycle
              sigia(:,:) = sigia(:,:)+
     &         sigi(:,:)*dab2(0:nx+1,0:ny+1,ia,1)                       !som 04.11.21
             enddo
             call b2txcx (nx, ny, mode, vol, gs, hx, qz, sigi,
     &                                             csigia(-1,-1,0,0))
             call b2txcy (nx, ny, mode, vol, gs, hx, hy, hy1,
     &                               qc, qs, sigi, csigia(-1,-1,0,1))
             call b2txcx (nx, ny, mode, vol, gs, hx, qz, sigia,
     &                                 csigin(-1,-1,0,0,is,is0))        !som 04.11.21
             call b2txcy (nx, ny, mode, vol, gs, hx, hy, hy1,
     &                       qc, qs, sigia, csigin(-1,-1,0,1,is,is0))   !som 04.11.21
             if (iout.ne.0) then
               write (chns,'(i3.3)') is
               write (chns0,'(i3.3)') is0
               call my_out(70,nx,ny,sigia,'b2tinnt_sigia'
     &                                       //chns//'_'//chns0)
               call my_out(70,nx,ny,csigia(-1,-1,0,0),'b2tinnt_csigiax'
     &                                       //chns//'_'//chns0)
               call my_out(70,nx,ny,csigia(-1,-1,1,1),'b2tinnt_csigiay'
     &                                       //chns//'_'//chns0)
             end if

*       X-COMPONENT OF ION-NEUTRAL CURRENT

             if(style.eq.0) then
               do iy=0,ny-1
                do ix=-1,nx
                 if(leftix(ix,iy).ne.-2) then
                  fchin(ix,iy,0,0)=fchin(ix,iy,0,0)+
     &              csigin(ix,iy,0,0,is,is0)*(                          !som 04.11.21
     &            ( -dpo(ix,iy,0,0)
     &              -(wrk1(ix,iy,0,0)
     &               +til(ix,iy)*wrk0(ix,iy,0,0,is))/
     &                     (qe*rza(ix,iy,is)) )*
     &                        (wbbl(ix,iy,2)/wbbl(ix,iy,3))**2 )
     &         +csigia(ix,iy,0,0)*wbbl(ix,iy,2)*fnal(ix,iy,0,0,is0)*
     &              (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
                  if (use_9pt_stencil.ne.0)
     >             fchin(ix,iy,0,1)=fchin(ix,iy,0,1)+
     &              csigin(ix,iy,0,1,is,is0)*(
     &              -dpo(ix,iy,0,1)
     &              -(wrk1(ix,iy,0,1)
     &               +til(ix,iy)*wrk0(ix,iy,0,1,is))/
     &                     (qe*rza(ix,iy,is)) )
     &         -csigia(ix,iy,0,1)*wbbl(ix,iy,3)*fnal(ix,iy,0,1,is0)*
     &              (hy1(ix,iy)+hy1(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
                 endif
                enddo
               enddo
             else
               do iy=0,ny-1
                do ix=-1,nx
                 if(leftix(ix,iy).ne.-2) then
                  fchin(ix,iy,0,0)=fchin(ix,iy,0,0)+
     &              csigin(ix,iy,0,0,is,is0)*(                          !som 04.11.21
     &             ( -dpo(ix,iy,0,0)
     &               -wrk0(ix,iy,0,0,is)/
     &              (nal(ix,iy,is)*qe*rza(ix,iy,is)) )*
     &              (wbbl(ix,iy,2)/wbbl(ix,iy,3))**2
     &                                                      )
     &         +csigia(ix,iy,0,0)*wbbl(ix,iy,2)*fnal(ix,iy,0,0,is0)*
     &              (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
                  if (use_9pt_stencil.ne.0)
     >             fchin(ix,iy,0,1)=fchin(ix,iy,0,1)+
     &              csigin(ix,iy,0,1,is,is0)*(
     &              -dpo(ix,iy,0,1)
     &              -wrk0(ix,iy,0,1,is)/
     &              (nal(ix,iy,is)*qe*rza(ix,iy,is)) )
     &         -csigia(ix,iy,0,1)*wbbl(ix,iy,3)*fnal(ix,iy,0,1,is0)*
     &              (hy1(ix,iy)+hy1(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
                 endif
                enddo
               enddo
             endif

*      Y-COMPONENT OF ION-NEUTRAL CURRENT
             if(style.eq.0) then
               do iy=0,ny
                do ix=-1,nx
                 if(leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1) then
                  fchin(ix,iy,1,1)=fchin(ix,iy,1,1)+
     &              csigin(ix,iy,1,1,is,is0)*(                          !som 04.11.21
     &            - dpo(ix,iy,1,1)
     &            -(wrk1(ix,iy,1,1)+tib(ix,iy)*wrk0(ix,iy,1,1,is))/
     &              (qe*rza(ix,iy,is))
     &                                                      )
     &        - csigia(ix,iy,1,1)*wbbv(ix,iy,3)*fnav(ix,iy,1,1,is0)*    !srv 11.01.13
     &          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8
                 endif
                enddo
               enddo
             else
               do iy=0,ny
                do ix=-1,nx
                 if(leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1) then
                  fchin(ix,iy,1,1)=fchin(ix,iy,1,1)+
     &              csigin(ix,iy,1,1,is,is0)*(  ! csigia(ix,iy,1,1)*(   !som 04.11.21
     &            - dpo(ix,iy,1,1)
     &            - wrk0(ix,iy,1,1,is)/
     &             (nab(ix,iy,is)*qe*rza(ix,iy,is))
     &                                                      )
     &        - csigia(ix,iy,1,1)*wbbv(ix,iy,3)*fnav(ix,iy,1,1,is0)*    !srv 11.01.13
     &          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8
                 endif
                enddo
               enddo
             endif                                                      !srv 15.02.11 }

            endif
           enddo ! is
          endif
         enddo ! is0                                                    !som 04.11.21 }

#ifdef B25_EIRENE
        else !if (use_kinetic .eq. 1) then                              !IYS 17.10.2023

         do is = 0,ns-1
! interpolate kinetic momentum source densities to cell faces
          sec(-1:nx,-1:ny) = smr0_eir_tot(-1:nx,-1:ny,is)
     &                             / hz(-1:nx,-1:ny) / vol(-1:nx,-1:ny)
          call intfaceh(nx,ny,leftix,leftiy,vol,
     &                       sec(-1,-1),work(-1,-1,1,is))
          sec(-1:nx,-1:ny) = smd0_eir_tot(-1:nx,-1:ny,is)
     &                             / hz(-1:nx,-1:ny) / vol(-1:nx,-1:ny)
          call intfaceh(nx,ny,leftix,leftiy,vol,
     &                       sec(-1,-1),work(-1,-1,0,is))
! HERE:
!   work(:,:,1,is) is the y-component of momentum source density
!   work(:,:,0,is) is the DIAMAGNETIC component of momentum source density,
!                  and the sign convention should be EIRENE-compatible
          do ix=-1,nx
           do iy=-1,ny
            if (leftix(ix,iy).ne.-2 .and. bottomiy(ix,iy).ne.-2
     &         .and. topiy(ix,iy).ne.ny+1) then
              fchin_a(ix,iy,0,0,is) = gs(ix,iy,0) * qc(ix,iy,0) *
     &         wbbl(ix,iy,2)/(wbbl(ix,iy,3))**2 * work(ix,iy,1,is)
             if (use_9pt_stencil.ne.0) then
              fchin_a(ix,iy,0,1,is) = - gs(ix,iy,0) * qc(ix,iy,1) *
     &         work(ix,iy,0,is) / wbbl(ix,iy,3)
             else
              fchin_a(ix,iy,0,0,is) = fchin_a(ix,iy,0,0,is)
     &         - gs(ix,iy,0) * qc(ix,iy,1) *
     &           work(ix,iy,0,is) / wbbl(ix,iy,3)
             endif
            endif
           enddo
          enddo
          call intfacev(nx,ny,bottomix,bottomiy,vol,
     &                       sec(-1,-1),work(-1,-1,0,is))
          do ix=-1,nx
           do iy=-1,ny
            if (leftix(ix,iy).ne.-2 .and. bottomiy(ix,iy).ne.-2
     &         .and. rightix(ix,iy).ne.nx+1)
     &       fchin_a(ix,iy,1,1,is) = gs(ix,iy,1) *
     &                                work(ix,iy,0,is) / wbbv(ix,iy,3)  ! Mind the sign!!!!
            fchin(ix,iy,0,0) = fchin(ix,iy,0,0) + fchin_a(ix,iy,0,0,is)
            fchin(ix,iy,0,1) = fchin(ix,iy,0,1) + fchin_a(ix,iy,0,1,is)
            fchin(ix,iy,1,1) = fchin(ix,iy,1,1) + fchin_a(ix,iy,1,1,is)
           enddo
          enddo
         enddo
#endif
        endif ! use_kinetic

        if(fchin_in_core.eq.0) then                                     !srv 16.11.17 {
          do iy=-1,ny
            do ix=-1,nx
              if(on_closed_surface(ix,iy)) then
                fchin(ix,iy,:,:) = 0.0_R8
              endif
            enddo
          enddo
        endif                                                           !srv 16.11.17 }

       else ! use fluid neutrals
*   ..extensive tests on first few calls
        if (ncall.lt.3) then
*    ..test sign of ne, te, chce
         call b2xvsg (n2, te, 1, 'te', '.gt.')
*    ..test edge values and sign of csigin
         do is = 0, ns-1
          do is0 = 0, ns-1
           call b2xvsg (n2, csigin(-1,-1,0,0,is,is0), 1, 'csigin',
     .            '.ge.')
           call b2xvsg (n2, csigin(-1,-1,1,1,is,is0), 1, 'csigin',
     .            '.ge.')
           call b2xvff (nx, ny, csigin(-1,-1,0,0,is,is0), 'csigin')
          enddo
         enddo
        endif
*   ..calculate the contribution of ion-neutral current
!        fchin = 0.0_R8
        fnan = 0.0_R8
        fnac = 0.0_R8
        do is0 = 0, ns-1
         if (is_neutral(is0).and.nint(zn(is0)).eq.1) then               !srv 20.10.17
*       X-COMPONENT OF ION-NEUTRAL CURRENT
*      ..get fluxes divided by area
          do iy=-1,ny
           do ix=-1,nx
            if(leftix(ix,iy).ne.-2) then
             fnan(ix,iy,0,0,is0)=fna(ix,iy,0,0,is0)/gs(ix,iy,0)
             if (area_fix.ge.1)
     >        fnan(ix,iy,0,0,is0)=fnan(ix,iy,0,0,is0)/qc(ix,iy,0)
            else
             fnan(ix,iy,0,0:1,is0)=0.0_R8
            endif
            if(bottomiy(ix,iy).ne.-2) then
             fnan(ix,iy,1,1,is0)=fna(ix,iy,1,1,is0)/gs(ix,iy,1)
            else
             fnan(ix,iy,1,1,is0)=0.0_R8
            endif
           enddo
          enddo
*       ..get fluxes in the cell centre
          fnac(:,:,0,is0) = fnan(:,:,0,0,is0)
          fnac(:,:,1,is0) = fnan(:,:,1,1,is0)
          call mapx(nx,ny,1,fnac(-1,-1,0,is0),work(-1,-1,0,is0))
          call mapy(nx,ny,1,fnac(-1,-1,1,is0),work(-1,-1,1,is0))
*      ..get velocities in the cell centre
          do iy=-1,ny
           do ix=-1,nx
            if(na(ix,iy,is0).gt.na_min) then
             vntrc(ix,iy,0) = fnac(ix,iy,1,is0)/na(ix,iy,is0)
             vntrc(ix,iy,1) = fnac(ix,iy,0,is0)/na(ix,iy,is0)
            else
             vntrc(ix,iy,0) = 0.0_R8
             vntrc(ix,iy,1) = 0.0_R8
            endif
           enddo
          enddo
*      ..interpolate velocities in the cell face
          call intfaceh(nx,ny,leftix,leftiy,vol,vntrc(-1,-1,0),vntrl)
          call intfaceh(nx,ny,leftix,leftiy,vol,vntrc(-1,-1,1),vntrb)
          do is = 0, ns-1
           if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then          !srv 20.10.17
            if(style.eq.0) then
             do iy=0,ny-1
              do ix=-1,nx
               if(leftix(ix,iy).ne.-2) then
                fchin(ix,iy,0,0)=fchin(ix,iy,0,0)+
     &                        csigin(ix,iy,0,0,is,is0)*(                !srv 24.06.08
     &                    ( -dpo(ix,iy,0,0)
     &                      -(wrk1(ix,iy,0,0)
     &                      +til(ix,iy)*wrk0(ix,iy,0,0,is))/
     &                      (qe*rza(ix,iy,is)) )*
     &                     (wbbl(ix,iy,2)/wbbl(ix,iy,3))**2
     &                      +wbbl(ix,iy,2)*vntrl(ix,iy)*(hx(ix,iy)+
     &                       hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8 )
                if (use_9pt_stencil.ne.0)
     >           fchin(ix,iy,0,1)=fchin(ix,iy,0,1)+
     &                        csigin(ix,iy,0,1,is,is0)*(
     &                      -dpo(ix,iy,0,1)
     &                      -(wrk1(ix,iy,0,1) +
     &                       til(ix,iy)*wrk0(ix,iy,0,1,is))/
     &                      (qe*rza(ix,iy,is))
     &                      -wbbl(ix,iy,3)*vntrb(ix,iy)*(hy1(ix,iy)+
     &                       hy1(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8 )
               endif
              enddo
             enddo
            else
             do iy=0,ny-1
              do ix=-1,nx
               if(leftix(ix,iy).ne.-2) then
                fchin(ix,iy,0,0)=fchin(ix,iy,0,0)+
     &                        csigin(ix,iy,0,0,is,is0)*(                !srv 24.06.08
     &                    ( -dpo(ix,iy,0,0)
     &                      -wrk0(ix,iy,0,0,is)/
     &                       (nal(ix,iy,is)*qe*rza(ix,iy,is)) )*
     &                       (wbbl(ix,iy,2)/wbbl(ix,iy,3))**2
     &                       +wbbl(ix,iy,2)*vntrl(ix,iy)*(hx(ix,iy)+
     &                      hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8 )
                if (use_9pt_stencil.ne.0)
     >           fchin(ix,iy,0,1)=fchin(ix,iy,0,1)+
     &                        csigin(ix,iy,0,1,is,is0)*(
     &                      -dpo(ix,iy,0,1)
     &                      -wrk0(ix,iy,0,1,is)/
     &                      (nal(ix,iy,is)*qe*rza(ix,iy,is))
     &                      -wbbl(ix,iy,3)*vntrb(ix,iy)*(hy1(ix,iy)+
     &                       hy1(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8 )
               endif
              enddo
             enddo
            endif
           endif
          enddo
*      Y-COMPONENT OF ION-NEUTRAL CURRENT
*      ..interpolate velocities in the cell face
          call intfacev(nx,ny,bottomix,bottomiy,vol,
     .                    vntrc(-1,-1,1),vntrb)
          do is = 0, ns-1
           if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then          !srv 20.10.17
            if(style.eq.0) then
             do iy=0,ny
              do ix=-1,nx
               if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
                fchin(ix,iy,1,1)=fchin(ix,iy,1,1)+
     &                      csigin(ix,iy,1,1,is,is0)*(
     &                  -dpo(ix,iy,1,1)
     &                  -(wrk1(ix,iy,1,1)+
     &                    tib(ix,iy)*wrk0(ix,iy,1,1,is))/
     &                   (qe*rza(ix,iy,is))-wbbv(ix,iy,3)*vntrb(ix,iy)* !srv 11.01.13
     &                   (hy1(ix,iy)+
     &                    hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8)
               endif
              enddo
             enddo
            else
             do iy=0,ny
              do ix=-1,nx
               if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
                fchin(ix,iy,1,1)=fchin(ix,iy,1,1)+
     &                   csigin(ix,iy,1,1,is,is0)*(
     &                  -dpo(ix,iy,1,1)
     &                  -wrk0(ix,iy,1,1,is)/
     &                   (nab(ix,iy,is)*qe*rza(ix,iy,is))
     &                  - wbbv(ix,iy,3)*vntrb(ix,iy)*(hy1(ix,iy)        !srv 11.01.13
     &                   +hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8)
               endif
              enddo
             enddo
            endif
           endif
          enddo ! is
         endif
        enddo ! is0                                                     !srv 08.07.09 }
       endif ! use_eirene
      endif

#ifdef B25_EIRENE
      if(iout.ne.0.and.use_eirene.ne.0) then
       do ia = 1,natm
        write (chns,'(i3.3)') ia
        call my_out(70,nx,ny,rfluxa(0:nx+1,0:ny+1,ia,1),
     &                'b2tinnt_rfluxa_eir'//chns)
        call my_out(70,nx,ny,pfluxa(0:nx+1,0:ny+1,ia,1),
     &                'b2tinnt_pfluxa_eir'//chns)
        call my_out(70,nx,ny,dab2(0:nx+1,0:ny+1,ia,1),
     &                'b2tinnt_dab2_eir'//chns)
       enddo
       do im = 1,nmol
        write (chns,'(i3.3)') im
        call my_out(70,nx,ny,rfluxm(0:nx+1,0:ny+1,im,1),
     &              'b2tinnt_rfluxm_eir'//chns)
        call my_out(70,nx,ny,pfluxm(0:nx+1,0:ny+1,im,1),
     &              'b2tinnt_pfluxm_eir'//chns)
        call my_out(70,nx,ny,dmb2(0:nx+1,0:ny+1,im,1),
     &              'b2tinnt_dmb2_eir'//chns)
       enddo
       do is0 = 0, ns-1                                                 !som 04.11.21 {
        write (chns0,'(i3.3)') is0
        if (is_neutral(is0).and.nint(zn(is0)).eq.1) then
         do is = 0, ns-1
          if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then
           write (chns,'(i3.3)') is
           call my_out(70,nx,ny,csigin(-1,-1,0,0,is,is0),
     &                                      'b2tinnt_csiginx'
     &                                       //chns//'_'//chns0)
           call my_out(70,nx,ny,csigin(-1,-1,1,1,is,is0),
     &                                      'b2tinnt_csiginy'
     &                                       //chns//'_'//chns0)
          endif
         enddo
         call my_out(70,nx,ny,fnan(-1,-1,0,0,is0),'b2tinnt_fnanx'
     &                                       //chns0)
         call my_out(70,nx,ny,fnan(-1,-1,1,1,is0),'b2tinnt_fnany'
     &                                       //chns0)
        endif
        if (use_kinetic.eq.1) then
          call my_out(70,nx,ny,smo0_eir_tot(-1,-1,is0),
     &                                   'b2tinnt_smo'//chns0)
          call my_out(70,nx,ny,smr0_eir_tot(-1,-1,is0),
     &                                   'b2tinnt_smr'//chns0)
          call my_out(70,nx,ny,smd0_eir_tot(-1,-1,is0),
     &                                   'b2tinnt_smd'//chns0)
        endif
      call my_out(70,nx,ny,fchin_a(-1,-1,0,0,is0),
     &                                   'b2tinnt_fchinx'//chns0)
      call my_out(70,nx,ny,fchin_a(-1,-1,1,1,is0),
     &                                   'b2tinnt_fchiny'//chns0)
       enddo                                                            !som 04.11.21 }
       call my_out(70,nx,ny,fchin(-1,-1,0,0),'b2tinnt_fchinx')
       call my_out(70,nx,ny,fchin(-1,-1,1,1),'b2tinnt_fchiny')
      endif
#endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tinnt

      end subroutine b2tinnt
