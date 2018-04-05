*-----------------------------------------------------------------------
*.specification

      subroutine b2siav (nx, ny, ns, isb, mb, hx, hy, hz, vol,
     & qz, bb, ub, wbdia, vbecrb, te, ti, nb, po, ni, ne2, lnlam, chci,  !srv 01.07.09
     & cvsb, cvsb_cl, cvsbhz, cvsbhz_cl, smbvh, smbvv, kt_neo,
     & wrk0, wrk1, wrk2) !srv 23.07.99 04.10.99 17.01.01 09.01.01 31.07.02 01.07.05
      use b2mod_types
      use b2mod_constants
      use b2mod_indirect
      use b2mod_geo_corner
      use b2mod_b2cmpa
      use b2mod_balance                                               !djm Jan2017
     & , only : b2siav_smovh0to3, b2siav_smovv0to3, balance_netcdf
      use b2mod_openmp
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, isb
      real (kind=R8) ::
     & mb, hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), hz(-1:nx,-1:ny),
     & vol(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1),
     1 bb(-1:nx,-1:ny,0:3), lnlam(-1:nx,-1:ny),                       !srv 01.07.09
     2 ub(-1:nx,-1:ny), wbdia(-1:nx,-1:ny,0:1), vbecrb(-1:nx,-1:ny,0:1),
     & te(-1:nx,-1:ny), ti(-1:nx,-1:ny), nb(-1:nx,-1:ny),
     & po(-1:nx,-1:ny), ni(-1:nx,-1:ny), ne2(-1:nx,-1:ny),
     & chci(-1:nx,-1:ny,0:1), cvsb(-1:nx,-1:ny,0:1),                  !srv 23.07.99
     & cvsb_cl(-1:nx,-1:ny,0:1), cvsbhz(-1:nx,-1:ny,0:1),             !srv 04.10.99 09.01.01
     & cvsbhz_cl(-1:nx,-1:ny,0:1)                                     !srv 31.07.02
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   smbvh(-1:nx,-1:ny,0:3), smbvv(-1:nx,-1:ny,0:3),
     &   kt_neo(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny), wrk1(-1:nx,-1:ny), wrk2(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIAV computes a linearised representation of the additional
*     viscosity.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, k
      integer, save :: qstyle = 0, lluciani = 3                     !lk 03.12.09
      integer, save :: iout = 0                                     !srv 17.06.02
!$OMP THREADPRIVATE(ncall, qstyle, lluciani, iout)
      character chns*3, chk*1                                       !srv 17.06.02
      real (kind=R8) ::
     *     t0, ft                                                   !srv 16.10.17
      real (kind=R8), save ::
     &     addvis=0.0_R8, addvis1=1.0_R8, cqip1 = 0.5e0_R8,         !srv 28.03.06 !srv 19.10.08
     &     cvsa_mltpl = 1.0e0_R8                                    !iyv 06.12.13
!$OMP THREADPRIVATE(addvis, addvis1, cqip1, cvsa_mltpl)
      real (kind=R8) ::
     *   tauia(-1:nx,-1:ny,0:ns-1), qip0(-1:nx,-1:ny),              !srv 23.07.99
     *   qip1(-1:nx,-1:ny), qipgen(-1:nx,-1:ny),                    !srv 23.07.99 19.01.05
     &   tt0(-1:ny), tt1(-1:ny), bbaverage2(-1:ny),                 !lk 23.05.06 05.07.07 {
     &   conn(-1:ny), tt2(-1:ny), tt3(-1:ny), eps1(-1:ny)
      real (kind=R8) :: qip0a(-1:nx,-1:ny), qip0b(-1:nx,-1:ny)
     &                , tauia1(-1:nx,-1:ny), alpha1(-1:nx,-1:ny)
     &                , nu1(-1:nx,-1:ny), wrks(-1:nx,-1:ny), s2     !lk !iyv 02.12.13
      integer mintopcut                                             !lk 23.05.06 05.07.07 }
      save ncall
      data ncall/0/
*   ..procedures
      intrinsic abs, sqrt
      real (kind=R8) :: hy1
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2siav')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2siav_addvis', addvis)                        !srv 23.07.99
       call ipgetr ('b2siav_addvis1', addvis1)                      !srv 04.10.99
       call ipgetr ('b2siav_cqip1', cqip1)                          !srv 28.03.06
       call ipgeti ('b2siav_iout', iout)                            !lk 23.05.06 05.07.07
       if ((iout.ne.0) .and. in_parallel()) then
         !$OMP MASTER
         write(*,*) 'B2SIAV OpenMP warning: no file out in ',
     &              'parallel mode'
         !$OMP END MASTER
         iout = 0
       end if
       call ipgeti ('b2siav_style_qip', qstyle)                     !lk 23.05.06 05.07.07
       call ipgeti ('b2trcl_lluciani', lluciani)                    !lk 03.12.09
       call ipgetr ('b2trcl_cvsa_mltpl', cvsa_mltpl)	            !iyv 06.12.13
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)

      ft = 4.0e0_R8/3.0e0_R8                                        !srv 16.10.17
      s2 = sqrt(2.0_R8)                                             !iyv 02.12.13
      smbvh  = 0.0_R8
      smbvv = 0.0_R8

      if(nint(zn(isb)).eq.1 .and. nint(zamax(isb)).eq.1) then       !srv 26.10.08
        if (addvis.ne.0.0_R8) then                                  !srv 19.01.05 28.03.06 {
          call b2ttia (nx, ny, ns, ti, ne2, lnlam, tauia)           !srv 01.07.09
* ..calculate the parallel heat flux of ions in the center of the cell

          tt0 = 0.0e0_R8                                            !lk 23.05.06 05.07.07
          tt1 = 0.0e0_R8
          tt2 = 0.0e0_R8
          tt3 = 0.0e0_R8
          conn = 0.0e0_R8
          bbaverage2 = 0.0e0_R8
          eps1 = 0.0e0_R8
!pbwdk          mintopcut = ny
          mintopcut = -2
          if (nncut > 0) mintopcut = minval(topcut(1:nncut))
          do iy = 0,mintopcut-1    ! for inner surfaces
            do ix=0,nx
              if(on_closed_surface(ix,iy)) then
                tt0(iy) = tt0(iy) +bb(ix,iy,3)*bb(ix,iy,3)*vol(ix,iy)
                tt1(iy) = tt1(iy) +vol(ix,iy)
                conn(iy) = conn(iy) +
     .                     hx(ix,iy)*bb(ix,iy,3)/abs(bb(ix,iy,0))   !srv 11.01.13
                tt2(iy) = tt2(iy) +hx(ix,iy)
                tt3(iy) = tt3(iy) +vol(ix,iy)/hy1(ix,iy)
              endif
            enddo
          enddo

          do iy = 0,mintopcut-1
            bbaverage2(iy) = tt0(iy)/tt1(iy)
            eps1(iy) = tt2(iy)**2/tt3(iy)
          enddo

          qip0 = 0.0e0_R8

          if (qstyle .eq. 0) then
            qip0a = 0.0e0_R8
            qip0b = 0.0e0_R8

            do iy = 0,mintopcut-1    ! for inner surfaces
              do ix = 0,nx
                if(on_closed_surface(ix,iy)) then
                  qip0a(ix,iy) = (-2.5e0_R8)*(bb(ix,iy,3)*bb(ix,iy,2)
     &             /bb(ix,iy,0)/bbaverage2(iy))
     &             *(ni(ix,iy)*ti(ix,iy)/qe)*
     &              (ti(topix(ix,iy),topiy(ix,iy))
     &              -ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &              (hy1(ix,iy)+0.5e0_R8*hy1(topix(ix,iy),topiy(ix,iy))+
     &               0.5e0_R8*hy1(bottomix(ix,iy),bottomiy(ix,iy)))
                endif
              enddo
            enddo

            do iy = 0,mintopcut-1
              do ix = -1,nx
                if(on_closed_surface(ix,iy)) then
                  if (leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)
     >             qip0b(ix,iy) = 0.25e0_R8*(qip0a(ix,iy)*2.0e0_R8
     &              +qip0a(leftix(ix,iy),leftiy(ix,iy))
     &              +qip0a(rightix(ix,iy),rightiy(ix,iy)))
                endif
              enddo
            enddo

            do iy = 0,mintopcut-1
              do ix = -1,nx
                if(on_closed_surface(ix,iy)) then
                  if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)
     >             qip0(ix,iy) = 0.25e0_R8*(qip0b(ix,iy)*2.0e0_R8
     &              +qip0b(leftix(ix,iy),leftiy(ix,iy))
     &              +qip0b(rightix(ix,iy),rightiy(ix,iy)))
                endif
              enddo
            enddo

          endif

          if (qstyle.eq.1) then

            wrk0 = 0.0e0_R8
            do iy=0,ny-1
              do ix=-1,nx
                if(leftix(ix,iy).ne.-2) then                        !srv 26.11.02
                  wrk0(ix,iy)=-chci(ix,iy,0)
     &             *(ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))
                endif
              enddo
            enddo

            do iy=-1,ny
              do ix=-1,nx-1
                if(on_closed_surface(ix,iy)) then
                  if(rightix(ix,iy).ne.nx+1 .and.
     &               topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2)
     >             qip0(ix,iy)=0.5e0_R8*(wrk0(ix,iy)+
     &               wrk0(rightix(ix,iy),rightiy(ix,iy)))
     &               *hx(ix,iy)*bb(ix,iy,3)/bb(ix,iy,0)/vol(ix,iy)
     &               -2.5e0_R8*bb(ix,iy,2)
     &               /bb(ix,iy,0)/bb(ix,iy,3)
     &               *(ni(ix,iy)*ti(ix,iy)/qe)*
     &                (ti(topix(ix,iy),topiy(ix,iy))
     &                -ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &                (hy1(ix,iy)+
     &                 0.5e0_R8*hy1(topix(ix,iy),topiy(ix,iy))+
     &                 0.5e0_R8*hy1(bottomix(ix,iy),bottomiy(ix,iy)))
                endif
              enddo
            enddo

          endif

          tauia1 = 0.0e0_R8
          nu1 = 0.0e0_R8
          alpha1 = 0.0e0_R8
          kt_neo = 0.0e0_R8                                         !iyv 02.12.13
          do iy = 0,mintopcut-1
            do ix=0,nx
              if(on_closed_surface(ix,iy)) then
                if (lluciani.eq.3) then
                  tauia1(ix,iy) = cvsa_mltpl*tauia(ix,iy,isb)/      !lk 20.11.07 { iyv 06.12.13
     &             (1.0e0_R8+15.12e16_R8/conn(iy)*
     &             (ti(ix,iy)/ev)**2/ne2(ix,iy))
     &             /(1.0e0_R8+15.12e16_R8*eps1(iy)**1.5e0_R8
     &             /conn(iy)*(ti(ix,iy)/ev)**2/ne2(ix,iy))
                elseif (lluciani.eq.1) then
                  tauia1(ix,iy) = tauia(ix,iy,isb)/
     1             (1.0_R8+7.5e16_R8/conn(iy)*
     2             (ti(ix,iy)/ev)**2/ne2(ix,iy))
                else
                  tauia1(ix,iy) = tauia(ix,iy,isb)
                endif
                nu1(ix,iy) = conn(iy)/(tauia(ix,iy,isb)
     &           *sqrt(ti(ix,iy)/(mp*am(isb)))
     &           *2.0e0_R8*pi*sqrt(eps1(iy)**3))
                alpha1(ix,iy) = ((-0.17e0_R8+1.05e0_R8*sqrt(nu1(ix,iy))+
     &           2.7e0_R8*nu1(ix,iy)**2*eps1(iy)**3)/(1.0e0_R8+
     &           0.7e0_R8*sqrt(nu1(ix,iy))+
     &           nu1(ix,iy)**2*eps1(iy)**3)-1.0e0_R8)*0.533e0_R8    !lk 20.11.07 }
                kt_neo(ix,iy) =
     &           ((-0.17e0_R8+1.05e0_R8*sqrt(nu1(ix,iy)*s2)         !iyv 02.12.13
     &           +2.7e0_R8*(nu1(ix,iy)*s2)**2*eps1(iy)**3)/(1.0e0_R8+
     &           0.7e0_R8*sqrt(nu1(ix,iy)*s2)+
     &           (nu1(ix,iy)*s2)**2*eps1(iy)**3))
              endif
            enddo
          enddo

          qip1 = 0.0e0_R8
          do iy=0,ny-1
            do ix=0,nx-1
              qip1(ix,iy)=cqip1*
     &         ni(ix,iy)*ti(ix,iy)*sqrt(te(ix,iy)+ti(ix,iy))
     &                                /sqrt(am(isb)*mp)
            enddo
          enddo

          qipgen = 0.0e0_R8
          do iy=0,ny-1
            do ix=0,nx-1
              qipgen(ix,iy)=qip1(ix,iy)*atan(qip0(ix,iy)/qip1(ix,iy))
     &                   *sqrt(bb(ix,iy,3))
            enddo
          enddo

* ..calculate the difference from parallel heat flux of ion
*   at the face of cell

          wrk1=0.0e0_R8
          do iy=0,ny-1
            do ix=-1,nx
              if(on_closed_surface(ix,iy)) then
                if(leftix(ix,iy).ne.-2) then
                  wrk1(ix,iy) = 2.0e0_R8*(qipgen(ix,iy) -
     &             qipgen(leftix(ix,iy),leftiy(ix,iy)))*
     &             (bb(ix,iy,0)/bb(ix,iy,3)**3
     &             *tauia1(ix,iy)*alpha1(ix,iy)
     &             *hx(leftix(ix,iy),leftiy(ix,iy))+
     &              bb(leftix(ix,iy),leftiy(ix,iy),0)/
     &              bb(leftix(ix,iy),leftiy(ix,iy),3)**3
     &             *tauia1(leftix(ix,iy),leftiy(ix,iy))
     &             *alpha1(leftix(ix,iy),leftiy(ix,iy))
     &             *hx(ix,iy))/
     &             (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))**2
                endif
              endif
            enddo
          enddo

* finally calculate the r.h.s. contribution to the parallel balance
* corresponding to parallel viscosity connected with heat flow
          do iy=0,ny-1
            do ix=0,nx-1
              if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
                t0=addvis*vol(ix,iy)/hx(ix,iy)*bb(ix,iy,0)*
     &            sqrt(bb(ix,iy,3))*hz(ix,iy)
                smbvh(ix,iy,0) = smbvh(ix,iy,0)
     &            +t0*(wrk1(rightix(ix,iy),rightiy(ix,iy))-wrk1(ix,iy))
              endif
            enddo
          enddo
        endif                                                       !lk 23.05.06 05.07.07

        if (addvis1.ne.0.0_R8) then          !srv 19.01.05 30.01.10

          wrk0=0
          wrk2 = sqrt(bb(:,:,3))     ! wrk2 is used below twice
          wrk1 = log(hz*wrk2)
          do iy = 0,ny-1                                            !srv 01.07.05 {
            do ix = 0,nx
              if(leftix(ix,iy).eq.-2) then
                wrk0(ix,iy)=0.0e0_R8
              else
                wrk0(ix,iy)=2.0e0_R8*ft*cvsbhz_cl(ix,iy,0)/         !srv 16.10.17
     &            (wrk2(ix,iy)+wrk2(leftix(ix,iy),leftiy(ix,iy)))*
     &            (wrk1(ix,iy)-wrk1(leftix(ix,iy),leftiy(ix,iy)))
              endif
            enddo
          enddo

          do iy= 0,ny-1
            do ix= 0,nx-1
              if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)
     >         then
                t0=(wrk0(rightix(ix,iy),rightiy(ix,iy))-wrk0(ix,iy))
     *             *wrk2(ix,iy)                                     !srv 01.07.05 }
                if (t0.gt.0.0e0) then
                  smbvv(ix,iy,0)=t0*ub(ix,iy)
                else
                  smbvv(ix,iy,1)=t0
                endif
              else
                smbvv(ix,iy,0) = 0.0e0_R8
                smbvv(ix,iy,1) = 0.0e0_R8
              endif
            enddo
          enddo                                                     !srv 19.01.05 }
        endif
      endif                                                         !srv 26.10.08

      if (iout.ne.0) then                                           !srv 17.06.02 {
        call my_out(70,nx,ny,alpha1,'b2siav_alpha1')
        call my_out(70,nx,ny,kt_neo,'b2siav_kt_neo')
        call my_out(70,nx,ny,tauia1,'b2siav_tauia1')
        write (chns,'(i3.3)') isb
        do k=0,3
          write (chk,'(i1)') k
          call my_out(70,nx,ny,smbvv(-1,-1,k),'b2siav_smovv'//chk//chns)
          call my_out(70,nx,ny,smbvh(-1,-1,k),'b2siav_smovh'//chk//chns)
        enddo
        wrks=smbvv(:,:,0)+smbvv(:,:,1)*ub(:,:)+
     &       smbvv(:,:,2)*nb(:,:)*mb+smbvv(:,:,3)*ub(:,:)*nb(:,:)*mb
        call my_out(70,nx,ny,wrks,'b2siav_smovv'//chns)
        wrks=smbvh(:,:,0)+smbvh(:,:,1)*ub(:,:)+
     &       smbvh(:,:,2)*nb(:,:)*mb+smbvh(:,:,3)*ub(:,:)*nb(:,:)*mb
        call my_out(70,nx,ny,wrks,'b2siav_smovh'//chns)
        call my_out(70,nx,ny,cvsbhz(-1,-1,0),'b2siav_cvsbhzx'//chns)
        call my_out(70,nx,ny,cvsbhz(-1,-1,1),'b2siav_cvsbhzy'//chns)
      endif                                                         !srv 17.06.02 }

cdjm Jan2017 Store linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2siav_smovh0to3(-1:nx,-1:ny,0:3,isb)=smbvh
        b2siav_smovv0to3(-1:nx,-1:ny,0:3,isb)=smbvv
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2siav

      end subroutine b2siav
