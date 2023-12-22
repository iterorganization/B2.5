*-----------------------------------------------------------------------
*.specification

      subroutine b2siav_zh (nx, ny, ns, isb, ismain, mb, zb,          !srv 19.02.21
     & ub, wbdia, vbecrb, te, ti, nb, ne2, lnlam,                     !srv 01.07.09
     & cvsbhz, cvsbhz_drho, smbvh, smbvv, smbvi, kt_neo,
     & alpha1, chci_al, chci_a, cvsahz_hAdp_albe, cvsahz_hBdp_al,     !som 08.10.21
     * cvsahz_uAdp_albe, cvsahz_uBdp_al,                              !som 22.12.23
     & wrk0, wrk1, wrk2, rho_rel, avm_u,
     & na, ua)
      use b2mod_types
      use b2mod_geo
      use b2mod_constants
      use b2mod_indirect
      use b2mod_geo_corner
      use b2mod_b2cmpa
      use b2mod_balance                                               !djm Jan2017
     & , only : b2siav_smovh0to3, b2siav_smovv0to3, balance_netcdf,
     &          b2siav_smovi0to3
      use b2mod_openmp
      use b2mod_plasma
     & , only : fac_vis                                               !srv 22.05.18
      use b2mod_subsys
      use b2mod_zhfrtf
     & , only : nal, Z_to_m1, Ia
      use b2mod_b2mrel
     & , only : amnucl, s2nucl, nucl2s, znnucl_int
      use b2mod_b2zhco
     & , only : zhcscorr, zhcscorr_vq, zh_vis_type
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, isb, ismain
      real (kind=R8) ::
     &  mb, lnlam(-1:nx,-1:ny), zb(-1:nx,-1:ny), ub(-1:nx,-1:ny),     !srv 01.07.09
     &  wbdia(-1:nx,-1:ny,0:1), vbecrb(-1:nx,-1:ny,0:1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), nb(-1:nx,-1:ny),
     &  ne2(-1:nx,-1:ny),
     &  cvsbhz(-1:nx,-1:ny,0:1),                                      !srv 04.10.99 09.01.01
     &  cvsbhz_drho(-1:nx,-1:ny,0:1),                                 !som 18.08.21
     &  chci_al(-1:nx,-1:ny,0:1,0:nspecies-1),                        !som 08.10.21
     &  chci_a(-1:nx,-1:ny,0:1,0:ns-1),                               !som 08.10.21
     &  cvsahz_hAdp_albe(-1:nx,-1:ny,0:1,0:nspecies-1,0:nspecies-1),  !som 08.10.21
     &  cvsahz_hBdp_al(-1:nx,-1:ny,0:1,0:nspecies-1),                 !som 08.10.21
     *  cvsahz_uAdp_albe(-1:nx,-1:ny,0:1,0:nspecies-1,0:nspecies-1), 
     *  cvsahz_uBdp_al(-1:nx,-1:ny,0:1,0:nspecies-1),                 !som 22.12.23
     &  rho_rel(-1:nx,-1:ny), avm_u(-1:nx,-1:ny),                     !som 13.08.21
     &  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  smbvh(-1:nx,-1:ny,0:3), smbvv(-1:nx,-1:ny,0:3),
     *  smbvi(-1:nx,-1:ny,0:3),                                       !som 13.08.21
     &  kt_neo(-1:nx,-1:ny), alpha1(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk0(-1:nx,-1:ny), wrk1(-1:nx,-1:ny), wrk2(-1:nx,-1:ny)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIAV_ZH computes a linearised representation of the additional
*     viscosity if the Zhdanov closure is applied.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, ix, iy, k, inucl, iz
      integer, save :: hstyle = 0, lluciani = 3                     !lk 03.12.09
      integer, save :: iout = 0                                     !srv 17.06.02
      !$OMP THREADPRIVATE(ncall, hstyle, lluciani, iout)
      character chns*3, chns0*3, chk*1                              !srv 17.06.02
      real (kind=R8) ::
     *  t0, ft,                                                     !srv 16.10.17
     &  ft_h                                                        !som 11.10.21
      real (kind=R8), save ::
     &  addvis=1.0_R8, addvis1=1.0_R8, chip1 = 0.5e0_R8,            !srv 28.03.06 !srv 19.10.08
     &  cvsa_mltpl = 1.0e0_R8                                       !iyv 06.12.13
      !$OMP THREADPRIVATE(addvis, addvis1, chip1, cvsa_mltpl)
      real (kind=R8) ::
     *  tauia(-1:nx,-1:ny,0:ns-1), hip0(-1:nx,-1:ny),               !srv 23.07.99
     *  hip1(-1:nx,-1:ny), hipgen(-1:nx,-1:ny),                     !srv 23.07.99 19.01.05
     *  hipgen_al(-1:nx,-1:ny,0:nspecies-1),                        !som 11.10.21
     &  tt0(-1:ny), tt1(-1:ny), bbaverage2(-1:ny),                  !lk 23.05.06 05.07.07 {
     &  tt2(-1:ny), tt3(-1:ny)
      real (kind=R8) :: hip0a(-1:nx,-1:ny), hip0b(-1:nx,-1:ny),
     &                  tauia1(-1:nx,-1:ny),
     &                  wrks(-1:nx,-1:ny)                           !lk !iyv 02.12.13
      real (kind=R8) :: u_res(-1:nx,-1:ny)                          !som 13.08.21
      integer mintopcut                                             !lk 23.05.06 05.07.07 }
      save ncall
      data ncall/0/
*   ..procedures
      intrinsic sqrt
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2siav_zh')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2siav_addvis', addvis)                        !srv 23.07.99
       call ipgetr ('b2siav_addvis1', addvis1)                      !srv 04.10.99
       call ipgetr ('b2siav_cqip1', chip1)                          !srv 28.03.06
       call ipgeti ('b2siav_iout', iout)                            !lk 23.05.06 05.07.07
       if ((iout.ne.0) .and. in_parallel()) then
         !$OMP MASTER
         write(*,*) 'b2siav_zh OpenMP warning: no file output in ',
     &              'parallel mode'
         !$OMP END MASTER
         iout = 0
       end if
       call ipgeti ('b2siav_style_qip', hstyle)                     !lk 23.05.06 05.07.07
       call ipgeti ('b2trcl_lluciani', lluciani)                    !lk 03.12.09
       call ipgetr ('b2trcl_cvsa_mltpl', cvsa_mltpl)                !iyv 06.12.13
      endif

      ft = 4.0e0_R8/3.0e0_R8                                        !srv 16.10.17
      ft_h =  8.0_R8/15.0_R8                                        !som 11.10.21
      smbvh  = 0.0_R8
      smbvv = 0.0_R8
      smbvi = 0.0_R8                                                !som 13.08.21

!     if(nint(zn(isb)).eq.1 .and. nint(zamax(isb)).eq.1) then       !srv 26.10.08
!     if(isb.eq.ismain) then                                        !srv 19.02.21
      if (.not.is_neutral(isb)) then                                !som 11.10.21
        if (addvis.ne.0.0_R8.and.maxval(fac_vis).ne.0.0_R8) then    !srv 19.01.05 28.03.06 { 22.05.18
          hip0  = 0.0e0_R8
          hip1  = 0.0e0_R8
          hipgen= 0.0e0_R8
          hip0a = 0.0e0_R8
          hip0b = 0.0e0_R8
          wrk0 = 0.0e0_R8

          call b2ttia (nx, ny, ns, ti, ne2, lnlam, tauia)           !srv 01.07.09
* ..calculate the parallel heat flux of ions in the center of the cell

          tt0 = 0.0e0_R8                                            !lk 23.05.06 05.07.07
          tt1 = 0.0e0_R8
          tt2 = 0.0e0_R8
          tt3 = 0.0e0_R8
          bbaverage2 = 0.0e0_R8
!pbwdk          mintopcut = ny
          mintopcut = -2
          if (nncut > 0) then
            if (nnreg(0).eq.7) then
              mintopcut=topcut(1)
            else
              mintopcut=minval(topcut(1:nncut))
            endif
          endif
          do iy = 0, mintopcut-1    ! for inner surfaces
            do ix = 0, nx
              if(on_closed_surface(ix,iy)) then
                tt0(iy) = tt0(iy) + bb(ix,iy,3)*bb(ix,iy,3)*vol(ix,iy)
                tt1(iy) = tt1(iy) + vol(ix,iy)
                tt2(iy) = tt2(iy) + hx(ix,iy)
                tt3(iy) = tt3(iy) + vol(ix,iy)/hy1(ix,iy)
              endif
            enddo
          enddo

          do iy = 0, mintopcut-1
            bbaverage2(iy) = tt0(iy)/tt1(iy)
          enddo

* ..calculate the parallel heat flux of the isb species
          hip0 = 0.0e0_R8

          if (hstyle .eq. 0 .and. minval(zb).gt.0.0_R8) then
            hip0a = 0.0e0_R8
            hip0b = 0.0e0_R8

            do iy = 0, mintopcut-1    ! for inner surfaces
              do ix = 0, nx
                if(on_closed_surface(ix,iy)) then
                  hip0a(ix,iy) = (-2.5e0_R8)*(bb(ix,iy,3)*bb(ix,iy,2)
     &             /bb(ix,iy,0)/bbaverage2(iy))
     &             *(nb(ix,iy)*ti(ix,iy)/qe/zb(ix,iy))*
     &              (ti(topix(ix,iy),topiy(ix,iy))
     &              -ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &              (hy1(ix,iy)+0.5e0_R8*hy1(topix(ix,iy),topiy(ix,iy))+
     &               0.5e0_R8*hy1(bottomix(ix,iy),bottomiy(ix,iy)))
                endif
              enddo
            enddo

            do iy = 0, mintopcut-1
              do ix = -1, nx
                if(on_closed_surface(ix,iy)) then
                  if (leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)
     >             hip0b(ix,iy) = 0.25e0_R8*(hip0a(ix,iy)*2.0e0_R8
     &              +hip0a(leftix(ix,iy),leftiy(ix,iy))
     &              +hip0a(rightix(ix,iy),rightiy(ix,iy)))
                endif
              enddo
            enddo

            do iy = 0, mintopcut-1
              do ix = -1, nx
                if(on_closed_surface(ix,iy)) then
                  if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)
     >             hip0(ix,iy) = 0.25e0_R8*(hip0b(ix,iy)*2.0e0_R8
     &              +hip0b(leftix(ix,iy),leftiy(ix,iy))
     &              +hip0b(rightix(ix,iy),rightiy(ix,iy)))
                endif
              enddo
            enddo

          endif

          if (hstyle.eq.1 .and. minval(zb).gt.0.0_R8) then

            wrk0 = 0.0e0_R8
            do iy=0,ny-1
              do ix=-1,nx
                if(leftix(ix,iy).ne.-2) then                        !srv 26.11.02
                  wrk0(ix,iy)=-chci_a(ix,iy,0,isb)
     &             *(ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))
                endif
              enddo
            enddo

            do iy=-1,ny
              do ix=-1,nx-1
                if(on_closed_surface(ix,iy)) then
                  if(rightix(ix,iy).ne.nx+1 .and.
     &               topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2)
     >             hip0(ix,iy)=0.5e0_R8*(wrk0(ix,iy)+
     &               wrk0(rightix(ix,iy),rightiy(ix,iy)))
     &               *hx(ix,iy)*bb(ix,iy,3)/bb(ix,iy,0)/vol(ix,iy)
     &               -2.5e0_R8*bb(ix,iy,2)
     &               /bb(ix,iy,0)/bb(ix,iy,3)
     &               *(nb(ix,iy)*ti(ix,iy)/qe/zb(ix,iy))*
     &                (ti(topix(ix,iy),topiy(ix,iy))
     &                -ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &                (hy1(ix,iy)+
     &                 0.5e0_R8*hy1(topix(ix,iy),topiy(ix,iy))+
     &                 0.5e0_R8*hy1(bottomix(ix,iy),bottomiy(ix,iy)))
                endif
              enddo
            enddo

          endif

          hip1 = 0.0e0_R8
          do iy=0,ny-1
            do ix=0,nx-1
              hip1(ix,iy)=chip1*
     &         nb(ix,iy)*ti(ix,iy)*sqrt(te(ix,iy)+ti(ix,iy))
     &                                /sqrt(am(isb)*mp)
            enddo
          enddo

          hipgen = 0.0e0_R8
          do iy=0,ny-1
            do ix=0,nx-1
              hipgen(ix,iy)=hip1(ix,iy)*atan(hip0(ix,iy)/hip1(ix,iy))
     &                   *sqrt(bb(ix,iy,3))
            enddo
          enddo

* ..calculate the parallel charge state averaged heat flux of all species
          do inucl= 0, nspecies-1
            hip0 = 0.0e0_R8

            if (hstyle .eq. 0 .and. minval(zb).gt.0.0_R8) then
              hip0a = 0.0e0_R8
              hip0b = 0.0e0_R8

              do iy = 0, mintopcut-1    ! for inner surfaces
                do ix = 0, nx
                  if(on_closed_surface(ix,iy)) then
                    hip0a(ix,iy) = (-2.5e0_R8)*(bb(ix,iy,3)*bb(ix,iy,2)
     &               /bb(ix,iy,0)/bbaverage2(iy))
     &               *(nal(ix,iy,inucl)*ti(ix,iy)/
     &                 qe*Z_to_m1(ix,iy,inucl))
     &               *(ti(topix(ix,iy),topiy(ix,iy))
     &                -ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &                (hy1(ix,iy)+
     &                 0.5e0_R8*hy1(topix(ix,iy),topiy(ix,iy))+
     &                 0.5e0_R8*hy1(bottomix(ix,iy),bottomiy(ix,iy)))
                  endif
                enddo
              enddo

              do iy = 0, mintopcut-1
                do ix = -1, nx
                  if(on_closed_surface(ix,iy)) then
                    if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)
     >               hip0b(ix,iy) = 0.25e0_R8*(hip0a(ix,iy)*2.0e0_R8
     &                +hip0a(leftix(ix,iy),leftiy(ix,iy))
     &                +hip0a(rightix(ix,iy),rightiy(ix,iy)))
                  endif
                enddo
              enddo

              do iy = 0, mintopcut-1
                do ix = -1, nx
                  if(on_closed_surface(ix,iy)) then
                    if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)
     >               hip0(ix,iy) = 0.25e0_R8*(hip0b(ix,iy)*2.0e0_R8
     &                +hip0b(leftix(ix,iy),leftiy(ix,iy))
     &                +hip0b(rightix(ix,iy),rightiy(ix,iy)))
                  endif
                enddo
              enddo

            endif

            if (hstyle.eq.1 .and. minval(zb).gt.0.0_R8) then

              wrk0 = 0.0e0_R8
              do iy=0,ny-1
                do ix=-1,nx
                  if(leftix(ix,iy).ne.-2) then                        !srv 26.11.02
                    wrk0(ix,iy)=-chci_al(ix,iy,0,inucl)
     &               *(ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))
                  endif
                enddo
              enddo

              do iy=-1,ny
                do ix=-1,nx-1
                  if(on_closed_surface(ix,iy)) then
                    if(rightix(ix,iy).ne.nx+1 .and.
     &                 topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2)
     >               hip0(ix,iy)=0.5e0_R8*(wrk0(ix,iy)
     &                +wrk0(rightix(ix,iy),rightiy(ix,iy)))
     &                *hx(ix,iy)*bb(ix,iy,3)/bb(ix,iy,0)/vol(ix,iy)
     &                -2.5e0_R8*bb(ix,iy,2)
     &                /bb(ix,iy,0)/bb(ix,iy,3)
     &                *(nal(ix,iy,inucl)*ti(ix,iy)
     &                 /qe*Z_to_m1(ix,iy,inucl))
     &                *(ti(topix(ix,iy),topiy(ix,iy))
     &                 -ti(bottomix(ix,iy),bottomiy(ix,iy)))
     &                /(hy1(ix,iy)+
     &                  0.5e0_R8*hy1(topix(ix,iy),topiy(ix,iy))+
     &                  0.5e0_R8*hy1(bottomix(ix,iy),bottomiy(ix,iy)))
                  endif
                enddo
              enddo

            endif

            hip1 = 0.0e0_R8
            do iy=0,ny-1
              do ix=0,nx-1
                hip1(ix,iy)=chip1*
     &           nal(ix,iy,inucl)*ti(ix,iy)*sqrt(te(ix,iy)+ti(ix,iy))
     &                                  /sqrt(amnucl(inucl)*mp)
              enddo
            enddo

            hipgen_al(-1:nx,-1:ny,inucl) = 0.0e0_R8
            do iy=0,ny-1
              do ix=0,nx-1
                hipgen_al(ix,iy,inucl)=hip1(ix,iy)*
     &           atan(hip0(ix,iy)/hip1(ix,iy))*sqrt(bb(ix,iy,3))
              enddo
            enddo
          enddo

* ..main contribution
* ..calculate the difference from parallel heat flux of ion
*   at the face of cell

          wrk0=0.0_R8
          wrk2 = sqrt(bb(:,:,3))
          wrk1 = log(hz*wrk2)
          do iy=0,ny-1
            do ix=-1,nx
              if(on_closed_surface(ix,iy)) then
                if(leftix(ix,iy).ne.-2) then
                  do inucl=0,nspecies-1
                    wrk0(ix,iy) = wrk0(ix,iy) +
     &           (nb(ix,iy)/nal(ix,iy,s2nucl(isb)-1)+
     &           nb(leftix(ix,iy),leftiy(ix,iy))/
     &           nal(leftix(ix,iy),leftiy(ix,iy),s2nucl(isb)-1))*
     &           ft_h*cvsahz_hAdp_albe(ix,iy,0,s2nucl(isb)-1,inucl)/
     &           (bb(ix,iy,0)*bb(ix,iy,3)*
     &           vol(ix,iy)/hx(ix,iy)*hz(ix,iy)+
     &           bb(leftix(ix,iy),leftiy(ix,iy),0)*
     &           bb(leftix(ix,iy),leftiy(ix,iy),3)*
     &           vol(leftix(ix,iy),leftiy(ix,iy))/
     &           hx(leftix(ix,iy),leftiy(ix,iy))*
     &           hz(leftix(ix,iy),leftiy(ix,iy)))*
     &           (hipgen_al(ix,iy,inucl) -
     &           hipgen_al(leftix(ix,iy),leftiy(ix,iy),inucl))
                  enddo
                endif
              endif
            enddo
          enddo

* finally calculate the r.h.s. contribution to the parallel balance
* corresponding to parallel viscosity connected with heat flow
          do iy=0,ny-1
            do ix=0,nx-1
              if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
                t0=fac_vis(ix,iy)*addvis*                           !srv 22.05.18
     &            vol(ix,iy)/hx(ix,iy)*bb(ix,iy,0)*
     &            sqrt(bb(ix,iy,3))*hz(ix,iy)
                smbvh(ix,iy,0) = smbvh(ix,iy,0)
     &            +t0*(wrk0(rightix(ix,iy),rightiy(ix,iy))-wrk0(ix,iy))
              endif
            enddo
          enddo

* ..charge state correction
          if (zhcscorr.gt.0 .and. zhcscorr_vq.gt.0) then
* ..calculate the difference from parallel heat flux of ion
*   at the face of cell

            wrk0=0.0_R8
            wrk2 = sqrt(bb(:,:,3))
            wrk1 = log(hz*wrk2)
            do iy=0,ny-1
              do ix=-1,nx
                if(on_closed_surface(ix,iy)) then
                  if(leftix(ix,iy).ne.-2) then
                    wrk0(ix,iy) =
     &           (nb(ix,iy)/nal(ix,iy,s2nucl(isb)-1)+
     &           nb(leftix(ix,iy),leftiy(ix,iy))/
     &           nal(leftix(ix,iy),leftiy(ix,iy),s2nucl(isb)-1))*
     &           ft_h*cvsahz_hBdp_al(ix,iy,0,s2nucl(isb)-1)/
     &           (bb(ix,iy,0)*bb(ix,iy,3)*
     &           vol(ix,iy)/hx(ix,iy)*hz(ix,iy)+
     &           bb(leftix(ix,iy),leftiy(ix,iy),0)*
     &           bb(leftix(ix,iy),leftiy(ix,iy),3)*
     &           vol(leftix(ix,iy),leftiy(ix,iy))/
     &           hx(leftix(ix,iy),leftiy(ix,iy))*
     &           hz(leftix(ix,iy),leftiy(ix,iy)))*
     &           ((hipgen(ix,iy) -
     &           hipgen(leftix(ix,iy),leftiy(ix,iy)))/
     &           Ia(ix,iy,isb)-
     &           (hipgen_al(ix,iy,s2nucl(isb)-1) -
     &           hipgen_al(leftix(ix,iy),leftiy(ix,iy),s2nucl(isb)-1)))
                  endif
                endif
              enddo
            enddo

* finally calculate the r.h.s. contribution to the parallel balance
* corresponding to parallel viscosity connected with heat flow
            do iy=0,ny-1
              do ix=0,nx-1
                if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then
                  t0=fac_vis(ix,iy)*addvis*                         !srv 22.05.18
     &              vol(ix,iy)/hx(ix,iy)*bb(ix,iy,0)*
     &              sqrt(bb(ix,iy,3))*hz(ix,iy)
                  smbvh(ix,iy,0) = smbvh(ix,iy,0)
     &             +t0*(wrk0(rightix(ix,iy),rightiy(ix,iy))-wrk0(ix,iy))
                endif
              enddo
            enddo
          endif
        endif                                                       !lk 23.05.06 05.07.07
      endif                                                         !som 11.10.21

        if (addvis1.ne.0.0_R8) then                                 !srv 19.01.05 30.01.10

          wrk0=0
          wrk2 = sqrt(bb(:,:,3))     ! wrk2 is used below twice
          wrk1 = log(hz*wrk2)
          do iy = 0,ny-1                                            !srv 01.07.05 {
            do ix = 0,nx
              if(leftix(ix,iy).eq.-2) then
                wrk0(ix,iy)=0.0e0_R8
              else
                if(zh_vis_type.eq.2) then
                  wrk0(ix,iy)=2.0e0_R8*ft*cvsbhz(ix,iy,0)/          !som 21.12.23
     &            (wrk2(ix,iy)+wrk2(leftix(ix,iy),leftiy(ix,iy)))*
     &            (wrk1(ix,iy)-wrk1(leftix(ix,iy),leftiy(ix,iy)))
                else
                  wrk0(ix,iy)=2.0e0_R8*ft*cvsbhz_drho(ix,iy,0)/       !som 18.08.21
     &            (wrk2(ix,iy)+wrk2(leftix(ix,iy),leftiy(ix,iy)))*  !som 18.08.21
     &            (wrk1(ix,iy)-wrk1(leftix(ix,iy),leftiy(ix,iy)))
                endif
              endif
            enddo
          enddo

          do iy= 0,ny-1
            do ix= 0,nx-1
              if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)
     >         then
                if (zh_vis_type.eq.0 .or. zh_vis_type.eq.2) then    !som 13.01.22 21.12.23
                 t0=(wrk0(rightix(ix,iy),rightiy(ix,iy))-wrk0(ix,iy))
     *             *wrk2(ix,iy)
                else if (zh_vis_type.eq.1) then
                 t0=(wrk0(rightix(ix,iy),rightiy(ix,iy))-wrk0(ix,iy))
     *             *wrk2(ix,iy)*rho_rel(ix,iy)                      !srv 01.07.05 } !som 13.08.21
                endif
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

* viscosity source with contribution from other ion velocities      !som 13.08.21
          wrk0=0
* calculate part of the mass-averaged velocity which depends on the other
* ion velocities
          if (zh_vis_type.eq.1) then
           u_res(-1:nx,-1:ny)=
     &      avm_u(-1:nx,-1:ny)-rho_rel(-1:nx,-1:ny)*ub(-1:nx,-1:ny)
           do iy=0,ny-1
            do ix=0,nx
              if(leftix(ix,iy).eq.-2) then
                wrk0(ix,iy)=0.0e0_R8
              else
                wrk0(ix,iy) = 2.0e0_R8*ft*cvsbhz_drho(ix,iy,0)/     !som 21.10.21
     &           (bb(ix,iy,0)*bb(ix,iy,3)*
     &           vol(ix,iy)/hx(ix,iy)*hz(ix,iy)+
     &           bb(leftix(ix,iy),leftiy(ix,iy),0)*
     &           bb(leftix(ix,iy),leftiy(ix,iy),3)*
     &           vol(leftix(ix,iy),leftiy(ix,iy))/
     &           hx(leftix(ix,iy),leftiy(ix,iy))*
     &           hz(leftix(ix,iy),leftiy(ix,iy)))*
     &           (u_res(ix,iy)*wrk2(ix,iy) -
     &           u_res(leftix(ix,iy),leftiy(ix,iy))*
     &           wrk2(leftix(ix,iy),leftiy(ix,iy)))
              endif
            enddo
           enddo

           do iy=0,ny-1
            do ix=0,nx-1
              if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
                smbvi(ix,iy,0)=vol(ix,iy)/hx(ix,iy)*hz(ix,iy)*
     &            bb(ix,iy,0)*wrk2(ix,iy)*
     &            (wrk0(rightix(ix,iy),rightiy(ix,iy))-wrk0(ix,iy))
              else
                smbvi(ix,iy,0) = 0.0e0_R8
              endif
            enddo
           enddo
          else if (zh_vis_type.eq.2) then
           do iy=0,ny-1
            do ix=0,nx
              wrk0(ix,iy)=0.0e0_R8
              if(leftix(ix,iy).eq.-2) then
                wrk0(ix,iy)=0.0e0_R8
              else
               do inucl=0,nspecies-1
                if(inucl.ne.s2nucl(isb)-1) then
                 t0=0.0e0_R8
                 do iz= 1, znnucl_int(inucl)
                  t0 = t0 + 0.5e0_R8*(na(ix,iy,nucl2s(inucl,iz))*
     &            ti(ix,iy)+
     &            na(leftix(ix,iy),leftiy(ix,iy),nucl2s(inucl,iz))*
     &            ti(leftix(ix,iy),leftiy(ix,iy)))*
     &            (ua(ix,iy,nucl2s(inucl,iz))*wrk2(ix,iy) -
     &            ua(leftix(ix,iy),leftiy(ix,iy),
     &            nucl2s(inucl,iz))*
     &            wrk2(leftix(ix,iy),leftiy(ix,iy)))
                 enddo
                
                 wrk0(ix,iy) = wrk0(ix,iy) +
     &           (nb(ix,iy)/nal(ix,iy,s2nucl(isb)-1)+
     &           nb(leftix(ix,iy),leftiy(ix,iy))/
     &           nal(leftix(ix,iy),leftiy(ix,iy),s2nucl(isb)-1))*
     &           ft*cvsahz_uAdp_albe(ix,iy,0,s2nucl(isb)-1,inucl)/
     &           (bb(ix,iy,0)*bb(ix,iy,3)*
     &           vol(ix,iy)/hx(ix,iy)*hz(ix,iy)+
     &           bb(leftix(ix,iy),leftiy(ix,iy),0)*
     &           bb(leftix(ix,iy),leftiy(ix,iy),3)*
     &           vol(leftix(ix,iy),leftiy(ix,iy))/
     &           hx(leftix(ix,iy),leftiy(ix,iy))*
     &           hz(leftix(ix,iy),leftiy(ix,iy)))*
     &           t0
                endif
               enddo
              endif
            enddo
           enddo

           do iy=0,ny-1
            do ix=0,nx-1
              if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
                smbvi(ix,iy,0)=vol(ix,iy)/hx(ix,iy)*hz(ix,iy)*
     &            bb(ix,iy,0)*wrk2(ix,iy)*
     &            (wrk0(rightix(ix,iy),rightiy(ix,iy))-wrk0(ix,iy))
              else
                smbvi(ix,iy,0) = 0.0e0_R8
              endif
            enddo
           enddo
          endif
        endif
!     endif                                                         !srv 26.10.08

      if (iout.ne.0) then                                           !srv 17.06.02 {
        write (chns,'(i3.3)') isb
        call my_out(70,nx,ny,alpha1,'b2siav_zh_alpha1'//chns)
        call my_out(70,nx,ny,kt_neo,'b2siav_zh_kt_neo'//chns)
        call my_out(70,nx,ny,tauia1,'b2siav_zh_tauia1'//chns)
        do k=0,3
          write (chk,'(i1)') k
          call my_out(70,nx,ny,smbvv(-1,-1,k),'b2siav_zh_smovv'
     &      //chk//chns)
          call my_out(70,nx,ny,smbvv(-1,-1,k),'b2siav_zh_smovi'
     &      //chk//chns)
          call my_out(70,nx,ny,smbvh(-1,-1,k),'b2siav_zh_smovh'
     &      //chk//chns)
        enddo
        wrks=smbvv(:,:,0)+smbvv(:,:,1)*ub(:,:)+
     &       smbvv(:,:,2)*nb(:,:)*mb+smbvv(:,:,3)*ub(:,:)*nb(:,:)*mb
        call my_out(70,nx,ny,wrks,'b2siav_zh_smovv'//chns)
        wrks=smbvi(:,:,0)+smbvi(:,:,1)*ub(:,:)+
     &       smbvi(:,:,2)*nb(:,:)*mb+smbvi(:,:,3)*ub(:,:)*nb(:,:)*mb
        call my_out(70,nx,ny,wrks,'b2siav_zh_smovi'//chns)             !som 13.08.21
        wrks=smbvh(:,:,0)+smbvh(:,:,1)*ub(:,:)+
     &       smbvh(:,:,2)*nb(:,:)*mb+smbvh(:,:,3)*ub(:,:)*nb(:,:)*mb
        call my_out(70,nx,ny,wrks,'b2siav_zh_smovh'//chns)
        call my_out(70,nx,ny,cvsbhz(-1,-1,0),'b2siav_zh_cvsbhzx'//chns)
        call my_out(70,nx,ny,cvsbhz(-1,-1,1),'b2siav_zh_cvsbhzy'//chns)

        call my_out(70,nx,ny,hipgen(-1,-1),'b2siav_zh_hipgen'//chns)
        call my_out(70,nx,ny,nb(-1,-1),'b2siav_zh_nb'//chns)
        call my_out(70,nx,ny,1.0e0_R8/zb(-1:nx,-1:ny),
     &   'b2siav_zh_zb_to_m1'//chns)

        do k=0,ns-1
          write (chns0,'(i3.3)') k
          call my_out(70,nx,ny,Ia(-1,-1,k),
     &     'b2siav_zh_Ia'//chns//'_'//chns0)
        enddo
        do inucl=0,nspecies-1
          write (chns0,'(i3.3)') inucl
          call my_out(70,nx,ny,hipgen_al(-1,-1,inucl),
     &     'b2siav_zh_hipgen_al'//chns//'_'//chns0)
          call my_out(70,nx,ny,nal(-1,-1,inucl),
     &     'b2siav_zh_nal'//chns//'_'//chns0)
          call my_out(70,nx,ny,Z_to_m1(-1,-1,inucl),
     &     'b2siav_zh_Z_to_m1'//chns//'_'//chns0)
        enddo
      endif                                                         !srv 17.06.02 }

cdjm Jan2017 Store linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2siav_smovh0to3(-1:nx,-1:ny,0:3,isb)=smbvh
        b2siav_smovv0to3(-1:nx,-1:ny,0:3,isb)=smbvv
        b2siav_smovi0to3(-1:nx,-1:ny,0:3,isb)=smbvi                 !som 13.08.21
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2siav_zh

      end subroutine b2siav_zh
