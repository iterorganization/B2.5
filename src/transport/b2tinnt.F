*-----------------------------------------------------------------------
*.specification

      subroutine b2tinnt (nx, ny, ns, facdrift, fac_ExB,
     &       ismain, ismain0, vol, gs, hx, hy, qz, qc, bb,               !srv 08.07.09 { !srv 20.10.17
     &       te, po, ti, na, fna, csigin, fchin)                         !srv 08.03.99 16.06.09
      use b2mod_types
      use b2mod_indirect
      use b2mod_rates
      use b2mod_geo_corner
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_neutrals_namelist
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, ismain0
      real (kind=R8) ::
     & facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny),                       !srv 20.10.17
     * vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2), hx(-1:nx,-1:ny),
     1 hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), bb(-1:nx,-1:ny,0:3), 
     1 te(-1:nx,-1:ny), po(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     & na(-1:nx,-1:ny,0:ns-1), fna(-1:nx,-1:ny,0:1,0:ns-1), 
     & csigin(-1:nx,-1:ny,0:1,0:ns-1,0:ns-1), qc(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   fchin(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TINNT computes the ion-neutral current
*
*-----------------------------------------------------------------------
*.declarations
     
*   ..local variables
      integer, save :: pot_eq = 1                                         !srv 26.07.06
      integer, save :: mode = 1
      integer, save :: iout = 0
      integer, save :: use_eirene, area_fix
      integer, save :: fchin_in_core = 0                                  !srv 16.11.17
      integer n2, ix, iy, is, is0, iamain
#ifdef B25_EIRENE
      integer ia, im
      character chns*3
#endif
      integer, save :: ncall = 0, style = 1
      real (kind=R8), save :: fch_ion_neutral = 0.0_R8
      real (kind=R8) ::
     *   fnac(-1:nx,-1:ny,0:1,0:ns-1), vntrc(-1:nx,-1:ny,0:1),            !srv 08.07.09
     *   vntrl(-1:nx,-1:ny), vntrb(-1:nx,-1:ny),
     *   til(-1:nx,-1:ny), tib(-1:nx,-1:ny),
     *   fnan(-1:nx,-1:ny,0:1,0:ns-1), work(-1:nx,-1:ny,0:1,0:ns-1),
     *   nal(-1:nx,-1:ny,0:ns-1), nab(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     &   sec(-1:nx,-1:ny), sigia(-1:nx,-1:ny), csigia(-1:nx,-1:ny,0:1),   !srv 15.02.11
     &   sigi(-1:nx,-1:ny), facdriftm, fac_ExBm                           !srv 20.10.17
*   ..procedures
      real (kind=R8) ::
     *   hy1                                                              !srv 07.07.14
      external subini, subend, ipgeti, xertst
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)                                  !srv 07.07.14
*   ..initialisation
      data use_eirene/0/, area_fix/3/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tinnt')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgeti ('b2tinnt_iout', iout)
        call ipgeti ('b2tinnt_style', style)
        call ipgeti ('b2tral_mode', mode)
        call ipgeti ('b2mndr_eirene', use_eirene)
        call ipgeti ('b2news_area_fix', area_fix)
        call ipgeti ('b2news_poteq', pot_eq)                              !srv 26.07.06
        call ipgetr ('b2tfhe_neutral', fch_ion_neutral)
        call ipgeti ('b2tinnt_fchin_in_core', fchin_in_core)              !srv 16.11.17
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      fchin     = 0.0e0_R8                                                !srv 09.07.99
      facdriftm = maxval(facdrift)                                        !srv 20.10.17
      fac_ExBm  = maxval(fac_ExB)
      if(pot_eq.ne.0 .and.(facdriftm.gt.0.0_R8 .or. fac_ExBm.gt.0.0_R8)
     .               .and. fch_ion_neutral .ne. 0.0e0_R8) then            !srv 26.07.06 01.07.99
        if (use_eirene.ne.0) then
          fchin = 0.0e0_R8                                                !srv 15.02.11 {
*   ..calculate section
! From 5.2: Using the diffusivity of room temperature hydrogen gas
          sec = 3.2e-15_R8*sqrt(ti/ev/0.026_R8)
! From 5.0/5.1: Using the CX rate
!        k = 0
!        found = .false.
!        do while (.not.found .and. k.lt.nscx)
!          found = iscx(k).eq.ismain0
!        enddo
!        if (found) then
!          t0 = log(ti(ix,iy)/(am(ismain0)*ev))
!          sec = phm0*
!     &       expu(rlcx(ix,iy,0,ismain,k)+rlcx(ix,iy,1,ismain,k)*t0)
!        else go back to diffusivity relation
*   ..calculate sigin
          call xertst(is_neutral(ismain0),
     .   'ismain0 not a neutral species!')
          call xertst(nint(zn(ismain)).eq.1.and.nint(zn(ismain0)).eq.1,
     .   'ion-neutral current computation only supported for H/D/T !') 
          iamain = b2eatcr(ismain)     ! iamain should correspond to atom of deuterium
!*   ..get density in the guard cells (boundaries)
!      do ix = -1, nx
!       do iy = -1, ny
!        if(leftix(ix,iy).eq.-2) then
!         dab2_(ix,iy)=dab2(rightix(ix,iy)+1,rightiy(ix,iy)+1,iamain,1)
!        elseif(rightix(ix,iy).eq.nx+1) then
!         dab2_(ix,iy)=dab2(leftix(ix,iy)+1,leftiy(ix,iy)+1,iamain,1)
!        else
!         dab2_(ix,iy)=dab2(ix+1,iy+1,iamain,1)
!        endif
!        if(topiy(ix,iy).eq.ny+1) then
!         dab2_(ix,iy)=dab2(bottomix(ix,iy)+1,bottomiy(ix,iy)+1,iamain,1)
!        elseif(bottomiy(ix,iy).eq.-2) then
!         dab2_(ix,iy)=dab2(topix(ix,iy)+1,topiy(ix,iy)+1,iamain,1)
!        endif
!       enddo
!      enddo
!*   ..get density in the corner cells
!      do iy = -1,ny,ny+1
!       do ix =-1,nx
!        if(leftix(ix,iy).eq.-2) then
!         dab2_(ix,iy)=dab2_(rightix(ix,iy),rightiy(ix,iy))
!        elseif(rightix(ix,iy).eq.nx+1) then
!         dab2_(ix,iy)=dab2_(leftix(ix,iy),leftiy(ix,iy))
!        endif
!       enddo
!      enddo
          sigi(:,:) = am(ismain)*mp*na(:,:,ismain)*sec/
     .                              (2.0_R8*bb(:,:,3)**2)
!      sigia(:,:) = sigi(:,:)*dab2_(:,:)
          sigia(:,:) = sigi(:,:)*dab2(0:nx+1,0:ny+1,iamain,1)
          call b2txcx (nx, ny, mode, vol, gs, hx, qz, sigi,
     &                                             csigia(-1,-1,0))
          call b2txcy (nx, ny, mode, vol, gs, hy, qz, sigi,
     &                                             csigia(-1,-1,1))
          call b2txcx (nx, ny, mode, vol, gs, hx, qz, sigia,
     &                                 csigin(-1,-1,0,ismain,ismain0))
          call b2txcy (nx, ny, mode, vol, gs, hy, qz, sigia,
     &                                 csigin(-1,-1,1,ismain,ismain0))
          if(style.eq.0) then
*   ..interpolate ti in the left and bottom cell face
            call intfaceh(nx, ny, leftix, leftiy, vol, ti, til)
            call intfacev(nx, ny, bottomix, bottomiy, vol, ti, tib)
          else
*   ..interpolate na in the left and bottom cell face
           call intfaceh(nx, ny, leftix, leftiy, vol,
     *                 na(-1,-1,ismain), nal(-1,-1,ismain))
           call intfacev(nx, ny, bottomix, bottomiy, vol,
     &                 na(-1,-1,ismain), nab(-1,-1,ismain))
          endif
*       X-COMPONENT OF ION-NEUTRAL CURRENT
*      ..get fluxes divided by area
          do iy=-1,ny
            do ix=-1,nx
              if(leftix(ix,iy).ne.-2) then
                fnan(ix,iy,0,ismain0)=pfluxa(ix+1,iy+1,iamain,1) !/gs(ix,iy,0)
              else
                fnan(ix,iy,0,ismain0)=0.0_R8
              endif
              if(bottomiy(ix,iy).ne.-2) then
                fnan(ix,iy,1,ismain0)=rfluxa(ix+1,iy+1,iamain,1) !/gs(ix,iy,1)
              else
                fnan(ix,iy,1,ismain0)=0.0_R8
              endif
            enddo
          enddo
*       ..get fluxes in the cell centre
          fnac(:,:,:,ismain0) = fnan(:,:,:,ismain0)

!!! These mapping routines seems to be not necessary here, since [rp]flux[am]
!!! comes already in cell centers                         ! IYS 27.01.2017
!        call mapx(nx,ny,1,fnac(-1,-1,0,ismain0),work(-1,-1,0,ismain0))
!        call mapy(nx,ny,1,fnac(-1,-1,1,ismain0),work(-1,-1,1,ismain0))
!!! These mapping routines seems to be not necessary here, since [rp]flux[am]
!!! comes already in cell centers                         ! IYS 27.01.2017

!*      ..get velocities in the cell centre
!        vntrc(:,:,0) = fnac(:,:,1,ismain0)/dab2(0:nx+1,0:ny+1,iamain,1)
!*      ..interpolate velocities in the cell face
!        call intfaceh(nx,ny,leftix,leftiy,vol,vntrc(-1,-1,0),vntrl)

*      ..interpolate fluxes in the cell face
          call intfaceh(nx,ny,leftix,leftiy,vol,fnac(-1,-1,1,ismain0),
     &                                        fnan(-1,-1,0,ismain0))
          if(style.eq.0) then
            do iy=0,ny-1
              do ix=-1,nx
                if(leftix(ix,iy).ne.-2) then
                  fchin(ix,iy,0)=fchin(ix,iy,0)+
     &              csigin(ix,iy,0,ismain,ismain0)*(
     &            ( -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))
     &              -(ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy))
     &               +til(ix,iy)*(log(na(ix,iy,ismain))-
     &                    log(na(leftix(ix,iy),leftiy(ix,iy),ismain))))/
     &                     (qe*rza(ix,iy,ismain)) )*
     &                        (wbbl(ix,iy,2)/wbbl(ix,iy,3))**2 )
     &             +csigia(ix,iy,0)*wbbl(ix,iy,2)*fnan(ix,iy,0,ismain0)*
     &              (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
!     &           +csigia(ix,iy,0)*wbbl(ix,iy,2)*vntrl(ix,iy)* 
!     &            (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
                endif
             enddo
           enddo
          else
            do iy=0,ny-1
              do ix=-1,nx
                if(leftix(ix,iy).ne.-2) then
                  fchin(ix,iy,0)=fchin(ix,iy,0)+
     &              csigin(ix,iy,0,ismain,ismain0)*(
     &             ( -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))
     &               -(na(ix,iy,ismain)*ti(ix,iy)-
     &              na(leftix(ix,iy),leftiy(ix,iy),ismain)*
     &              ti(leftix(ix,iy),leftiy(ix,iy)))/
     &              (nal(ix,iy,ismain)*qe*rza(ix,iy,ismain)) )*
     &              (wbbl(ix,iy,2)/wbbl(ix,iy,3))**2
     &                                                      )
     &             +csigia(ix,iy,0)*wbbl(ix,iy,2)*fnan(ix,iy,0,ismain0)*
     &              (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
!     &           +csigia(ix,iy,0)*wbbl(ix,iy,2)*vntrl(ix,iy)* 
!     &            (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
                endif
              enddo
            enddo
          endif
*      Y-COMPONENT OF ION-NEUTRAL CURRENT
!*      ..get velocities in the cell centre
!        vntrc(:,:,1) = fnac(:,:,0,ismain0)/dab2(0:nx+1,0:ny+1,iamain,1)
!*      ..interpolate velocities in the cell face
!        call intfacev(nx,ny,bottomix,bottomiy,vol,vntrc(-1,-1,1),vntrb)

*      ..interpolate fluxes in the cell face
          call intfacev(nx, ny, bottomix, bottomiy, vol,
     &                   fnac(-1,-1,0,ismain0), fnan(-1,-1,1,ismain0))
          if(style.eq.0) then
            do iy=0,ny
              do ix=-1,nx
                if(leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1) then
                  fchin(ix,iy,1)=fchin(ix,iy,1)+
     &              csigin(ix,iy,1,ismain,ismain0)*(
     &            -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))
     &            -(ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &              tib(ix,iy)*(log(na(ix,iy,ismain))-
     &                  log(na(bottomix(ix,iy),bottomiy(ix,iy),ismain)))
     &              )/(qe*rza(ix,iy,ismain))
     &                                                      )
     &            - csigia(ix,iy,1)*wbbv(ix,iy,2)*fnan(ix,iy,1,ismain0)*  !srv 11.01.13
     &          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8
!     &       - csigia(ix,iy,1)*wbbv(ix,iy,3)*vntrb(ix,iy)*
!     &        (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8
                endif
              enddo
            enddo
          else
            do iy=0,ny
              do ix=-1,nx
                if(leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1) then
                  fchin(ix,iy,1)=fchin(ix,iy,1)+
     &              csigin(ix,iy,1,ismain,ismain0)*(  ! csigia(ix,iy,1)*(
     &            -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))
     &            -(na(ix,iy,ismain)*ti(ix,iy)-
     &                      na(bottomix(ix,iy),bottomiy(ix,iy),ismain)*
     &                      ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &             (nab(ix,iy,ismain)*qe*rza(ix,iy,ismain))
     &                                                      )
     &            - csigia(ix,iy,1)*wbbv(ix,iy,2)*fnan(ix,iy,1,ismain0)*  !srv 11.01.13
     &          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8
!     &        - csigia(ix,iy,1)*wbbv(ix,iy,3)*vntrb(ix,iy)*
!     &         (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8
                endif
              enddo
            enddo
          endif                                                           !srv 15.02.11 }

          if(fchin_in_core.eq.0) then                                     !srv 16.11.17 {
            do iy=-1,ny
              do ix=-1,nx
                if(on_closed_surface(ix,iy)) then
                  fchin(ix,iy,:) = 0.0_R8
                endif
              enddo
            enddo
          endif                                                           !srv 16.11.17 }

        else
*   ..extensive tests on first few calls
          if (ncall.lt.3) then
*    ..test sign of ne, te, csigin, chce
            call b2xvsg (n2, te, 1, 'te', '.gt.')
            call b2xvsg (n2*2*ns*ns, csigin, 1, 'csigin', '.ge.')
*    ..test edge values of csigin
            do is = 0, ns-1
              do is0 = 0, ns-1
                call b2xvff (nx, ny, csigin(-1,-1,0,is,is0), 'csigin')
              enddo
            enddo
          endif
          if(style.eq.0) then
*   ..interpolate ti in the left and bottom cell face
            call intfaceh(nx, ny, leftix, leftiy, vol, ti, til)
            call intfacev(nx, ny, bottomix, bottomiy, vol, ti, tib)
          else
*   ..interpolate na in the left and bottom cell face
            do is = 0, ns-1
              call intfaceh(nx, ny, leftix, leftiy, vol,
     *                    na(-1,-1,is), nal(-1,-1,is))
              call intfacev(nx, ny, bottomix, bottomiy, vol,
     &                    na(-1,-1,is), nab(-1,-1,is))
            enddo
          endif
*   ..calculate the contribution of ion-neutral current
          fchin = 0.0e0_R8
          do is0 = 0, ns-1
!            if (is_neutral(is0)) then
            if (is_neutral(is0).and.nint(zn(is0)).eq.1) then             !srv 20.10.17
*       X-COMPONENT OF ION-NEUTRAL CURRENT
*      ..get fluxes divided by area
              do iy=-1,ny
                do ix=-1,nx
                  if(leftix(ix,iy).ne.-2) then
                    fnan(ix,iy,0,is0)=fna(ix,iy,0,is0)/gs(ix,iy,0)
                    if (area_fix.ge.1)
     >               fnan(ix,iy,0,is0)=fnan(ix,iy,0,is0)/qc(ix,iy)
                  else
                    fnan(ix,iy,0,is0)=0.0_R8
                  endif
                  if(bottomiy(ix,iy).ne.-2) then
                    fnan(ix,iy,1,is0)=fna(ix,iy,1,is0)/gs(ix,iy,1)
                  else
                    fnan(ix,iy,1,is0)=0.0_R8
                  endif
                enddo
              enddo
*       ..get fluxes in the cell centre
              fnac(:,:,:,is0) = fnan(:,:,:,is0)
              call mapx(nx,ny,1, fnac(-1,-1,0,is0), work(-1,-1,0,is0))
              call mapy(nx,ny,1, fnac(-1,-1,1,is0), work(-1,-1,1,is0))
*      ..get velocities in the cell centre
              vntrc(:,:,0) = fnac(:,:,1,is0)/na(:,:,is0)
              vntrc(:,:,1) = fnac(:,:,0,is0)/na(:,:,is0)
*      ..interpolate velocities in the cell face
              call intfaceh(nx,ny,leftix,leftiy,vol,vntrc(-1,-1,0),
     &                                                         vntrl)
              do is = 0, ns-1
!                if (.not.is_neutral(is)) then
                if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then       !srv 20.10.17
                  if(style.eq.0) then
                    do iy=0,ny-1
                      do ix=-1,nx
                        if(leftix(ix,iy).ne.-2) then
                          fchin(ix,iy,0)=fchin(ix,iy,0)+
     &                        csigin(ix,iy,0,is,is0)*(                    !srv 24.06.08
     &                    ( -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))
     &                       -(ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy))
     &                      +til(ix,iy)*(log(na(ix,iy,is))-
     &                        log(na(leftix(ix,iy),leftiy(ix,iy),is))))/
     &                      (qe*rza(ix,iy,is)) )*
     &                     (wbbl(ix,iy,2)/wbbl(ix,iy,3))**2
     &                      +wbbl(ix,iy,2)*vntrl(ix,iy)*(hx(ix,iy)+
     &                       hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8 )
                        endif
                      enddo
                    enddo
                  else
                    do iy=0,ny-1
                      do ix=-1,nx
                        if(leftix(ix,iy).ne.-2) then
                          fchin(ix,iy,0)=fchin(ix,iy,0)+
     &                        csigin(ix,iy,0,is,is0)*(                    !srv 24.06.08
     &                    ( -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))
     &                      -(na(ix,iy,is)*ti(ix,iy)-
     &                        na(leftix(ix,iy),leftiy(ix,iy),is)*
     &                        ti(leftix(ix,iy),leftiy(ix,iy)))/
     &                       (nal(ix,iy,is)*qe*rza(ix,iy,is)) )*
     &                       (wbbl(ix,iy,2)/wbbl(ix,iy,3))**2
     &                       +wbbl(ix,iy,2)*vntrl(ix,iy)*(hx(ix,iy)+
     &                      hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8 )
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
!                if (.not.is_neutral(is)) then
                if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then      !srv 20.10.17
                  if(style.eq.0) then
                    do iy=0,ny
                      do ix=-1,nx
                        if(leftix(ix,iy).ne.-2 .and.
     &                     rightix(ix,iy).ne.nx+1) then
                          fchin(ix,iy,1)=fchin(ix,iy,1)+
     &                      csigin(ix,iy,1,is,is0)*(
     &                  -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))
     &                  -(ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                    tib(ix,iy)*(log(na(ix,iy,is))-
     &                    log(na(bottomix(ix,iy),bottomiy(ix,iy),is))))/
     &                   (qe*rza(ix,iy,is))-wbbv(ix,iy,2)*vntrb(ix,iy)*   !srv 11.01.13
     &                   (hy1(ix,iy)+
     &                    hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8)
                        endif
                      enddo
                    enddo
                  else
                    do iy=0,ny
                      do ix=-1,nx
                        if(leftix(ix,iy).ne.-2 .and.
     &                     rightix(ix,iy).ne.nx+1) then
                          fchin(ix,iy,1)=fchin(ix,iy,1)+
     &                   csigin(ix,iy,1,is,is0)*(
     &                  -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))
     &                  -(na(ix,iy,is)*ti(ix,iy)-
     &                      na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &                      ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &                   (nab(ix,iy,is)*qe*rza(ix,iy,is))
     &                  - wbbv(ix,iy,2)*vntrb(ix,iy)*(hy1(ix,iy)          !srv 11.01.13
     &                   +hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8)
                        endif
                      enddo
                    enddo
                  endif
                endif
              enddo
            endif
          enddo                                                           !srv 08.07.09 }
        endif
      endif  !
#ifdef B25_EIRENE
      if(iout.ne.0) then
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
        call my_out(70,nx,ny,pfluxa(0:nx+1,0:ny+1,im,1),
     &              'b2tinnt_pfluxm_eir'//chns)
        call my_out(70,nx,ny,dmb2(0:nx+1,0:ny+1,im,1),
     &              'b2tinnt_dmb2_eir'//chns)
       enddo
!       call my_out(70,nx,ny,dab2_,'b2tinnt_dab2_01')
       call my_out(70,nx,ny,sigia,'b2tinnt_sigia')
       call my_out(70,nx,ny,csigia(-1,-1,0),'b2tinnt_csigiax')
       call my_out(70,nx,ny,csigia(-1,-1,1),'b2tinnt_csigiay')
       call my_out(70,nx,ny,csigin(-1,-1,0,ismain,ismain0),
     &                                      'b2tinnt_csiginx')
       call my_out(70,nx,ny,csigin(-1,-1,1,ismain,ismain0),
     &                                      'b2tinnt_csiginy')
       call my_out(70,nx,ny,fnan(-1,-1,0,ismain0),'b2tinnt_fnanx')
       call my_out(70,nx,ny,fnan(-1,-1,1,ismain0),'b2tinnt_fnany')
       call my_out(70,nx,ny,fchin(-1,-1,0),'b2tinnt_fchinx')
       call my_out(70,nx,ny,fchin(-1,-1,1),'b2tinnt_fchiny')
!       call my_out(70,nx,ny,vntrc(-1,-1,0),'b2tinnt_vntrcx')
!       call my_out(70,nx,ny,vntrc(-1,-1,1),'b2tinnt_vntrcy')
      endif
#endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tinnt

      end
