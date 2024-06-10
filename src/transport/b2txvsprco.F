      subroutine b2txvsprco(nx, ny, isb, fac_vis, vol,                  !srv 02.04.18
     &    zb, ti, kmpr, nb, avisper_b, bvisper_b, cvisper_b, dvisper_b)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_anomalous_transport
      use b2mod_rates
      use b2mod_b2cmpa
      use b2mod_geo_corner
      use b2mod_openmp
      use b2mod_geo
     , , only : hy1
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb
      real (kind=R8) ::
     1   fac_vis(-1:nx,-1:ny),
     1   vol(-1:nx,-1:ny),
     &   ti(-1:nx,-1:ny), nb(-1:nx,-1:ny),
     &   kmpr(-1:nx,-1:ny), zb(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     3   avisper_b(-1:nx,-1:ny,0:1), bvisper_b(-1:nx,-1:ny,0:1),        !srv 27.01.00
     &   cvisper_b(-1:nx,-1:ny,0:1), dvisper_b(-1:nx,-1:ny,0:1)
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix, iy
      integer, save :: ncall = 0, iout = 0                              !srv 06.12.17
      real (kind=R8), save :: vis_per = 1.0_R8                          !srv 22.11.99
      real (kind=R8) ::
     &  fac_vism, t0, wrk0(-1:nx,-1:ny)
      character(3) :: chns
      !$OMP THREADPRIVATE(ncall, iout, vis_per)
*   ..procedures
      external b2xvsg, b2xvff, b2xvfx
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2txvsprco')
*   ..set internal parameters on first call
      if (ncall.eq.0) then                                              !srv
       call ipgetr ('b2tfhe_vis_per', vis_per)                          !srv 22.11.99
       call ipgeti ('b2txvsprco_iout', iout)                            !srv 06.12.17
       if (in_parallel() .and. iout.ne.0) then
         !$OMP MASTER
         write(*,*)
     .    'b2txvsprco OpenMP warning: no file IO in parallel mode'
         !$OMP END MASTER
         iout = 0
       endif
      end if
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..calculate the coefficients in difference approximation of viscosity perpendicular current       !srv 22.11.99
      avisper_b  = 0.0_R8; bvisper_b = 0.0_R8; cvisper_b = 0.0_R8;      !srv 27.01.00 22.03.02
      dvisper_b = 0.0_R8

      fac_vism = maxval(fac_vis)
      if(fac_vism.ne.0.0_R8 .and. vis_per.ne.0.0_R8) then               !srv 22.11.99

          do iy=-1,ny
            do ix=-1,nx
              if(bottomiy(ix,iy).eq.-2)then
                cvisper_b(ix,iy,1) = 0.0_R8
                dvisper_b(ix,iy,1) = 0.0_R8
              elseif(leftix(ix,iy).eq.-2.or.rightix(ix,iy).eq.nx+1) then
                cvisper_b(ix,iy,1) = 0.0_R8   ! we do not want to have fchvisper(:,:,1) in guard cells on West and East
                dvisper_b(ix,iy,1) = 0.0_R8   ! we do not want to have fchvisper(:,:,1) in guard cells on West and East
              else
                t0 = 2.0_R8/wbbv(ix,iy,3)/
     &            (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
                cvisper_b(ix,iy,1) = t0*vis_per*fac_vis(ix,iy)
                dvisper_b(ix,iy,1) = t0*vis_per*fac_vis(ix,iy)*
     &           (ti(ix,iy)*kmpr(ix,iy)+                              ! in order to avoid nonlinear entry of pb
     &            ti(bottomix(ix,iy),bottomiy(ix,iy))*                ! we change it by na/kmpr
     &            kmpr(bottomix(ix,iy),bottomiy(ix,iy)))/
     &            (qe*(zb(ix,iy)*nb(ix,iy)+
     &            zb(bottomix(ix,iy),bottomiy(ix,iy))*
     &            nb(bottomix(ix,iy),bottomiy(ix,iy))))
              endif
            enddo
          enddo

          do iy=-1,ny
            do ix=-1,nx
              if(leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1) then ! we do not want to have fchvisper(:,:,1) in guard cells on West and East
                t0 = 0.25_R8*vsa0(ix,iy,isb)*vol(ix,iy)/hy1(ix,iy)**2
                bvisper_b(ix,iy,1) = t0*vis_per*fac_vis(ix,iy)
              endif
            enddo
          enddo

          avisper_b(:,:,1) = cvisper_b(:,:,1)

      endif                                                             !srv 22.11.99

      if(iout.eq.1) then
       call my_out(70,nx,ny,avisper_b(-1,-1,1),'b2txvsprco_avisper_by')
       call my_out(70,nx,ny,bvisper_b(-1,-1,1),'b2txvsprco_bvisper_by')
       call my_out(70,nx,ny,cvisper_b(-1,-1,1),'b2txvsprco_cvisper_by')
       call my_out(70,nx,ny,dvisper_b(-1,-1,1),'b2txvsprco_dvisper_by')
       do iy=-1,ny
         do ix=-1,nx
           wrk0(ix,iy) = hy1(ix,iy)
         enddo
       enddo
       call my_out(70,nx,ny,wrk0,'b2tvsprco_hy1')
       call my_out(70,nx,ny,wbbv(-1,-1,3),'b2tvsprco_bbv')
       call my_out(70,nx,ny,vol,'b2tvspcor_vol')
       write (chns,'(i3.3)') isb
       call my_out(70,nx,ny,vsa0(-1,-1,isb),'b2txvsprco_vsa0_'//chns)
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2txvsprco

      end subroutine b2txvsprco                                         !srv 02.04.18

