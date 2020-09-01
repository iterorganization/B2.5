*-----------------------------------------------------------------------
*.specification

      subroutine b2tinnt (nCv, nFc, nVx, ns, switch, geo, mpg, rt, 
     &       facdrift, fac_ExB, ismain, ismain0, te, po, ti, na, fna, 
     &       csigin, fchin)                                              !srv 08.03.99 16.06.09
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_neutrals_namelist
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain, ismain0
      real (kind=R8) ::
     & facdrift(nFc), fac_ExB(nFc), te(nCv), po(nCv), ti(nCv),
     & na(nCv,0:ns-1), fna(nFc,0:1,0:ns-1), 
     & csigin(nFc,0:1,0:ns-1,0:ns-1)
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Rates), intent (in) :: rt
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   fchin(nFc,0:1)
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
      integer is, is0, iamain
#ifdef B25_EIRENE
      integer ia, im
      character chns*3
#endif
      integer, save :: ncall = 0, style = 1
      real (kind=R8), save :: eps = 1.0e-60_R8
      real (kind=R8) ::
     *   fnac(nCv,0:1,0:ns-1), vntrc(nCv,0:1),                            !srv 08.07.09
     *   vntrf(nFc,0:1), tif(nFc), naf(nFc,0:ns-1), dpo(nFc,0:1),
     *   fnan(nFc,0:1,0:ns-1),!WG_TODO work(-1:nx,-1:ny,0:1,0:ns-1),
     *   wrk(nCv), wrk0(nFc,0:1,0:ns-1), wrk1(nFc,0:1), wrk2(nFc),
     *   wrkVx(nVx)
      real (kind=R8) ::
     &   sec(nCv), sigia(nCv), 
     &   csigia(nFc,0:1),                                                 !srv 15.02.11
     &   sigi(nCv), facdriftm, fac_ExBm                                   !srv 20.10.17
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c      integer, save :: mode = 1
c      integer, save :: iout = 0
c      integer, save :: fchin_in_core = 0                                  !srv 16.11.17
*   ..procedures
*   ..initialisation
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tinnt')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c        call ipgeti ('b2tinnt_iout', iout)
c        call ipgeti ('b2tinnt_style', style)
c        call ipgeti ('b2tral_mode', mode)
c        call ipgeti ('b2tinnt_fchin_in_core', fchin_in_core)              !srv 16.11.17
      endif
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      fchin     = 0.0e0_R8                                                !srv 09.07.99
      facdriftm = maxval(facdrift)                                        !srv 20.10.17
      fac_ExBm  = maxval(fac_ExB)
      if(switch%pot_eq.ne.0 .and.(facdriftm.gt.0.0_R8 .or.
     . fac_ExBm.gt.0.0_R8)  
     .            .and. switch%fch_ion_neutral .ne. 0.0e0_R8) then         !srv 26.07.06 01.07.99
        if (switch%use_eirene.ne.0) then
          call xerrab ('b2tinnt not yet available for WG + EIRENE')
#ifdef WG_TODO
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
          call b2txcx (nx, ny, mode, vol, gs, hx, hy, qz, sigi,
     &                                             csigia(-1,-1,0,0))
          call b2txcy (nx, ny, mode, vol, gs, hx, hy, qz, qs, sigi,
     &                                             csigia(-1,-1,0,1))
          call b2txcx (nx, ny, mode, vol, gs, hx, hy, qz, sigia,
     &                                 csigin(-1,-1,0,0,ismain,ismain0))
        call b2txcy (nx, ny, mode, vol, gs, hx, hy, qz, qs,  sigia, 
     &                                 csigin(-1,-1,0,1,ismain,ismain0))
          if(style.eq.0) then
*   ..interpolate ti in the left and bottom cell face
            call intfaceh(nx, ny, leftix, leftiy, vol, ti, til)
            call intfacev(nx, ny, bottomix, bottomiy, vol, ti, tib)
*   ..pre-computation of differences
            wght = 1.0e0_R8
            wrk  = log(na(:,:,ismain))
            call diff(nx, ny, hx, hy, wght, qs, po,   dpo)
            call diff(nx, ny, hx, hy, wght, qs, wrk, 
     &        wrk0(-1,-1,0,0,ismain))
            call diff(nx, ny, hx, hy, wght, qs, ti, wrk1)
          else
*   ..interpolate na in the left and bottom cell face
           call intfaceh(nx, ny, leftix, leftiy, vol,
     *                 na(-1,-1,ismain), nal(-1,-1,ismain))
           call intfacev(nx, ny, bottomix, bottomiy, vol,
     &                 na(-1,-1,ismain), nab(-1,-1,ismain))
*   ..pre-computation of differences
            wght = 1.0e0_R8
            call diff(nx, ny, hx, hy, wght, qs, po,   dpo)
            call diff(nx, ny, hx, hy, wght, qs, na(-1,-1,ismain)*ti, 
     &        wrk0(-1,-1,0,0,ismain))
          endif
*       X-COMPONENT OF ION-NEUTRAL CURRENT
*      ..get fluxes divided by area, in cell centre
          do iy=-1,ny
            do ix=-1,nx
              if(leftix(ix,iy).ne.-2) then
                fnac(ix,iy,0,ismain0)=pfluxa(ix+1,iy+1,iamain,1) !/gs(ix,iy,0)
              else
                fnac(ix,iy,0,ismain0)=0.0_R8
              endif
              if(bottomiy(ix,iy).ne.-2) then
                fnac(ix,iy,1,ismain0)=rfluxa(ix+1,iy+1,iamain,1) !/gs(ix,iy,1)
              else
                fnac(ix,iy,1,ismain0)=0.0_R8
              endif
            enddo
          enddo

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
     &                                        fnan(-1,-1,0,0,ismain0))
          call intfaceh(nx,ny,leftix,leftiy,vol,fnac(-1,-1,0,ismain0),
     &                                        fnan(-1,-1,0,1,ismain0))
          if(style.eq.0) then
            do iy=0,ny-1
              do ix=-1,nx
                if(leftix(ix,iy).ne.-2) then
                  fchin(ix,iy,0,0)=fchin(ix,iy,0,0)+
     &              csigin(ix,iy,0,0,ismain,ismain0)*(
     &            ( -dpo(ix,iy,0,0)
     &              -(wrk1(ix,iy,0,0) + 
     &                til(ix,iy)*wrk0(ix,iy,0,0,ismain))/
     &                     (qe*rza(ix,iy,ismain)) )*
     &                        (wbbl(ix,iy,2)/wbbl(ix,iy,3))**2 )
     &         +csigia(ix,iy,0,0)*wbbl(ix,iy,2)*fnan(ix,iy,0,0,ismain0)*
     &              (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
!     &           +csigia(ix,iy,0,0)*wbbl(ix,iy,2)*vntrl(ix,iy)* 
!     &            (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
                  fchin(ix,iy,0,1)=fchin(ix,iy,0,1)+
     &              csigin(ix,iy,0,1,ismain,ismain0)*(
     &              -dpo(ix,iy,0,1)
     &              -(wrk1(ix,iy,0,1) + 
     &                til(ix,iy)*wrk0(ix,iy,0,1,ismain))/
     &                     (qe*rza(ix,iy,ismain)) )
     &         -csigia(ix,iy,0,1)*wbbl(ix,iy,2)*fnan(ix,iy,0,1,ismain0)*
     &              (hy1(ix,iy)+hy1(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
                endif
             enddo
           enddo
          else
            do iy=0,ny-1
              do ix=-1,nx
                if(leftix(ix,iy).ne.-2) then
                  fchin(ix,iy,0,0)=fchin(ix,iy,0,0)+
     &              csigin(ix,iy,0,0,ismain,ismain0)*(
     &             ( -dpo(ix,iy,0,0)
     &               -wrk0(ix,iy,0,0,ismain)/
     &              (nal(ix,iy,ismain)*qe*rza(ix,iy,ismain)) )*
     &              (wbbl(ix,iy,2)/wbbl(ix,iy,3))**2
     &                                                      )
     &         +csigia(ix,iy,0,0)*wbbl(ix,iy,2)*fnan(ix,iy,0,0,ismain0)*
     &              (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
!     &           +csigia(ix,iy,0,0)*wbbl(ix,iy,2)*vntrl(ix,iy)* 
!     &            (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
                  fchin(ix,iy,0,1)=fchin(ix,iy,0,1)+
     &              csigin(ix,iy,0,1,ismain,ismain0)*(
     &              -dpo(ix,iy,0,1)
     &              -wrk0(ix,iy,0,1,ismain)/
     &              (nal(ix,iy,ismain)*qe*rza(ix,iy,ismain)) )
     &         -csigia(ix,iy,0,1)*wbbl(ix,iy,2)*fnan(ix,iy,0,1,ismain0)*
     &              (hy1(ix,iy)+hy1(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
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
     &                   fnac(-1,-1,0,ismain0), fnan(-1,-1,1,1,ismain0))
          if(style.eq.0) then
            do iy=0,ny
              do ix=-1,nx
                if(leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1) then
                  fchin(ix,iy,1,1)=fchin(ix,iy,1,1)+
     &              csigin(ix,iy,1,1,ismain,ismain0)*(
     &            - dpo(ix,iy,1,1)
     &            -(wrk1(ix,iy,1,1)+tib(ix,iy)*wrk0(ix,iy,1,1,ismain))/
     &              (qe*rza(ix,iy,ismain))
     &                                                      )
     &        - csigia(ix,iy,1,1)*wbbv(ix,iy,2)*fnan(ix,iy,1,1,ismain0)*  !srv 11.01.13
     &          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8
!     &       - csigia(ix,iy,1,1)*wbbv(ix,iy,3)*vntrb(ix,iy)*
!     &        (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8
                endif
              enddo
            enddo
          else
            do iy=0,ny
              do ix=-1,nx
                if(leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1) then
                  fchin(ix,iy,1,1)=fchin(ix,iy,1,1)+
     &              csigin(ix,iy,1,1,ismain,ismain0)*(  ! csigia(ix,iy,1,1)*(
     &            - dpo(ix,iy,1,1)
     &            - wrk0(ix,iy,1,1,ismain)/
     &             (nab(ix,iy,ismain)*qe*rza(ix,iy,ismain))
     &                                                      )
     &        - csigia(ix,iy,1,1)*wbbv(ix,iy,2)*fnan(ix,iy,1,1,ismain0)*  !srv 11.01.13
     &          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8
!     &        - csigia(ix,iy,1,1)*wbbv(ix,iy,3)*vntrb(ix,iy)*
!     &         (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*0.5_R8
                endif
              enddo
            enddo
          endif                                                           !srv 15.02.11 }

          if(fchin_in_core.eq.0) then                                     !srv 16.11.17 {
            do iy=-1,ny
              do ix=-1,nx
                if(on_closed_surface(ix,iy)) then
                  fchin(ix,iy,:,:) = 0.0_R8
                endif
              enddo
            enddo
          endif                                                           !srv 16.11.17 }
#endif
        else
*   ..extensive tests on first few calls
          if (ncall.lt.3) then
*    ..test sign of te
            call b2xvsg (nCv, te, 1, 'te', '.gt.')
*    ..test sign of csigin
            do is = 0, ns-1
              do is0 = 0, ns-1
                call b2xvsg (nFc, csigin(:,0,is,is0)*geo%fcQalf(:,0), 1,
     .            'csigin0', '.ge.')
                call b2xvsg (nFc, csigin(:,1,is,is0)*geo%fcQalf(:,1), 1,
     .            'csigin1', '.ge.')
              enddo
            enddo
          endif
          if(style.eq.0) then
*   ..interpolate ti to cell faces
            call intface(nCv, nFc, mpg%fcCv, geo%fcVol, ti, tif)
*   ..pre-computation of differences
            call diff(nCv, nFc, nVx, 0, geo, mpg, po, wrkVx, dpo)
            do is = 0, ns-1
              wrk = log(na(:,is))                                      !srv 12.05.18
              call diff(nCv, nFc, nVx, 0, geo, mpg, wrk, wrkVx,
     &          wrk0(:,:,is))
            enddo
            call diff(nCv, nFc, nVx, 0, geo, mpg, ti, wrkVx, wrk1)
          else
            do is = 0, ns-1
*   ..interpolate na to cell faces
              call intface(nCv, nFc, mpg%fcCv, geo%fcVol, na(:,is),
     &          naf(:,is))
*   ..pre-computation of differences
              call diff(nCv, nFc, nVx, 0, geo, mpg, na(:,is)*ti, wrkVx,
     &          wrk0(:,:,is))
            enddo
            call diff(nCv, nFc, nVx, 0, geo, mpg, po, wrkVx, dpo)
          endif
*   ..calculate the contribution of ion-neutral current
          do is0 = 0, ns-1
            if (is_neutral(is0).and.nint(zn(is0)).eq.1) then             !srv 20.10.17
*      ..get fluxes divided by area
              fnan(:,0,is0) = fna(:,0,is0)/(geo%fcS*geo%fcQalf(:,0)+eps)
              fnan(:,1,is0) = fna(:,1,is0)/(geo%fcS*geo%fcQalf(:,1)+eps)
*       ..get fluxes in the cell centre
              call intcell(nFc,nCv,mpg,mpg%intcellP,fnac(:,0,is0),
     &          fnac(:,0,is0))
              call intcell(nFc,nCv,mpg,mpg%intcellR,fnac(:,1,is0),
     &          fnac(:,1,is0))
*      ..get velocities in the cell centre
              vntrc(:,0) = fnac(:,1,is0)/na(:,is0)
              vntrc(:,1) = fnac(:,0,is0)/na(:,is0)
*      ..interpolate velocities in the cell face
              call intface(nCv, nFc, mpg%fcCv, geo%fcVol, vntrc(:,0),
     &          vntrf(:,0))
              call intface(nCv, nFc, mpg%fcCv, geo%fcVol, vntrc(:,1),
     &          vntrf(:,1))
*      ..components of ion-neutral current
              do is = 0, ns-1
                if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then       !srv 20.10.17
                  wrk  = qe*rt%rza(:,is)
                  call intface(nCv,nFc,mpg%fcCv,geo%fcVol,wrk,wrk2)
                  if(style.eq.0) then
                    fchin(:,0) = fchin(:,0) + csigin(:,0,is,is0)*(
     &                (-dpo(:,0) - (wrk1(:,0) + tif*wrk0(:,0,is))/wrk2)*
     &                (geo%fcBb(:,2)/geo%fcBb(:,3))**2 + 
     &                geo%fcBb(:,2)*vntrf(:,1)*geo%fcQgam(:,0)*
     &                (geo%fcHc(:,1) + geo%fcHc(:,2)))
                    fchin(:,1) = fchin(:,1) + csigin(:,1,is,is0)*(
     &                (-dpo(:,1) - (wrk1(:,1) + tif*wrk0(:,1,is))/wrk2)- 
     &                geo%fcBb(:,2)*vntrf(:,0)*geo%fcQgam(:,0)*
     &                (geo%fcHc(:,1) + geo%fcHc(:,2)))
                  else
                    fchin(:,0) = fchin(:,0) + csigin(:,0,is,is0)*(
     &                (-dpo(:,0) - wrk0(:,0,is)/(naf(:,is)*wrk2))*
     &                (geo%fcBb(:,2)/geo%fcBb(:,3))**2 + 
     &                geo%fcBb(:,2)*vntrf(:,1)*geo%fcQgam(:,0)*
     &                (geo%fcHc(:,1) + geo%fcHc(:,2)))
                    fchin(:,1) = fchin(:,1) + csigin(:,1,is,is0)*(
     &                (-dpo(:,1) - wrk0(:,1,is)/(naf(:,is)*wrk2))- 
     &                geo%fcBb(:,2)*vntrf(:,0)*geo%fcQgam(:,0)*
     &                (geo%fcHc(:,1) + geo%fcHc(:,2)))
                  endif
                endif
              enddo
            endif
          enddo                                                           !srv 08.07.09 }
        endif
      endif  !
#ifdef B25_EIRENE
#ifdef WG_TODO
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
       call my_out(70,nx,ny,csigia(-1,-1,0,0),'b2tinnt_csigiax')
       call my_out(70,nx,ny,csigia(-1,-1,1,1),'b2tinnt_csigiay')
       call my_out(70,nx,ny,csigin(-1,-1,0,0,ismain,ismain0),
     &                                      'b2tinnt_csiginx')
       call my_out(70,nx,ny,csigin(-1,-1,1,1,ismain,ismain0),
     &                                      'b2tinnt_csiginy')
       call my_out(70,nx,ny,fnan(-1,-1,0,0,ismain0),'b2tinnt_fnanx')
       call my_out(70,nx,ny,fnan(-1,-1,1,1,ismain0),'b2tinnt_fnany')
       call my_out(70,nx,ny,fchin(-1,-1,0,0),'b2tinnt_fchinx')
       call my_out(70,nx,ny,fchin(-1,-1,1,1),'b2tinnt_fchiny')
!       call my_out(70,nx,ny,vntrc(-1,-1,0),'b2tinnt_vntrcx')
!       call my_out(70,nx,ny,vntrc(-1,-1,1),'b2tinnt_vntrcy')
      endif
#endif
#endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tinnt

      end
