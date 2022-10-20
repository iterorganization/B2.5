      subroutine b2txvspr(nx, ny, ns, fac_vis, vol, hx, hy,              !srv 10.10.17
     &   qz, bb, avisper, bvisper, cvisper)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_anomalous_transport
      use b2mod_rates
      use b2mod_b2cmpa
      use b2mod_geo_corner
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     1   fac_vis(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     1   vol(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), bb(-1:nx,-1:ny,0:3) 
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     3   avisper(-1:nx,-1:ny,0:1), bvisper(-1:nx,-1:ny,0:1),            !srv 27.01.00
     $   cvisper(-1:nx,-1:ny,0:1)
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, ix, iy, is
      integer, save :: ncall = 0, iout = 0                              !srv 06.12.17
      real (kind=R8), save :: vis_per = 0.0e0_R8                        !srv 22.11.99
      real (kind=R8) ::
     &  fac_vism, t0, wrk0(-1:nx,-1:ny)
      character(3) :: chns
*   ..procedures
      intrinsic max, abs
      real (kind=R8) ::
     *   hy1
      external xertst
      external b2xvsg, b2xvff, b2xvfx
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)                                !srv 07.05.08
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2txvspr')
*   ..set internal parameters on first call
      if (ncall.eq.0) then                                              !srv
       call ipgetr ('b2tfhe_vis_per', vis_per)                          !srv 22.11.99
       call ipgeti ('b2txvspr_iout',iout)                               !srv 06.12.17
      end if
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..calculate the coefficients in difference approximation of viscosity perpendicular current       !srv 22.11.99 
      avisper  = 0.0e0_R8; bvisper = 0.0e0_R8; cvisper = 0.0e0_R8       !srv 27.01.00 22.03.02
        
      fac_vism = maxval(fac_vis)
      if(fac_vism.ne.0.0_R8 .and. vis_per.ne.0.0e0_R8) then             !srv 22.11.99

       do is = 0,ns-1                                                   !srv 16.10.17
         if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then             !srv 16.10.17
          do iy=-1,ny
            do ix=-1,nx
              if(bottomiy(ix,iy).eq.-2)then
                cvisper(ix,iy,1) = 0.0_R8
              elseif(leftix(ix,iy).eq.-2.or.rightix(ix,iy).eq.nx+1) then
                cvisper(ix,iy,1) = 0.0_R8   ! we do not want to have fchvisper(:,:,1) in guard cells on West and East
              else
                t0 = 2.0_R8/wbbv(ix,iy,3)/
     &            (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
                cvisper(ix,iy,1) = cvisper(ix,iy,1) + t0
              endif
            enddo
          enddo
        
          do iy=-1,ny
            do ix=-1,nx
              if(leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1) then ! we do not want to have fchvisper(:,:,1) in guard cells on West and East
                t0 = 0.25_R8*vsa0(ix,iy,is)*vol(ix,iy)/hy1(ix,iy)**2
                bvisper(ix,iy,1) = bvisper(ix,iy,1) + t0
              endif
            enddo
          enddo
          
          avisper(:,:,1) = cvisper(:,:,1)

         endif                                                          !srv 16.10.17
       enddo                                                            !srv 16.10.17
      endif                                                             !srv 22.11.99
      
      if(iout.eq.1) then
       call my_out(70,nx,ny,avisper(-1,-1,1),'b2txvspr_avispery')
       call my_out(70,nx,ny,bvisper(-1,-1,1),'b2txvspr_bvispery')
       call my_out(70,nx,ny,cvisper(-1,-1,1),'b2txvspr_cvispery')
       do iy=-1,ny
        do ix=-1,nx
         wrk0(ix,iy) = hy1(ix,iy)
        enddo
       enddo
       call my_out(70,nx,ny,wrk0,'b2tvspr_hy1')
       call my_out(70,nx,ny,wbbv(-1,-1,3),'b2tvspr_bbv')
       call my_out(70,nx,ny,vol,'b2tvspr_vol')
       do is = 0,ns-1
         write (chns,'(i3.3)') is
         if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then
          call my_out(70,nx,ny,vsa0(-1,-1,is),'b2tvspr_vsa0_'//chns)
         endif
       enddo
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2txvspr

      end subroutine b2txvspr                                            !srv 10.10.17

