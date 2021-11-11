      subroutine b2txvspr(nx, ny, ns, ismain, fac_vis, vol,             !srv 10.10.17 19.02.21
     &   avisper_a, bvisper_a, cvisper_a)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_anomalous_transport
      use b2mod_rates
      use b2mod_b2cmpa
      use b2mod_geo_corner
      use b2mod_geo
     , , only : hy1
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain                                        !srv 19.02.21
      real (kind=R8) ::
     1   fac_vis(-1:nx,-1:ny),
     1   vol(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   avisper_a(-1:nx,-1:ny,0:1,0:ns-1),                               !srv 27.01.00
     &   bvisper_a(-1:nx,-1:ny,0:1,0:ns-1),
     &   cvisper_a(-1:nx,-1:ny,0:1,0:ns-1)
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix, iy, is
      integer, save :: ncall = 0, iout = 0                              !srv 06.12.17
      real (kind=R8), save :: vis_per = 1.0e0_R8                        !srv 22.11.99
      integer, save :: visper_all_ions = 0
      real (kind=R8) ::
     &  fac_vism, t0, wrk0(-1:nx,-1:ny)
      character(3) :: chns
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2txvspr')
*   ..set internal parameters on first call
      if (ncall.eq.0) then                                              !srv
       call ipgetr ('b2tfhe_vis_per', vis_per)                          !srv 22.11.99
       call ipgeti ('b2txvspr_iout',iout)                               !srv 06.12.17
       call ipgeti ('b2tvspr_visper_all_ions', visper_all_ions)
      end if
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..calculate the coefficients in difference approximation of viscosity perpendicular current       !srv 22.11.99 
      avisper_a  = 0.0e0_R8; bvisper_a = 0.0e0_R8; cvisper_a = 0.0e0_R8       !srv 27.01.00 22.03.02
        
      fac_vism = maxval(fac_vis)
      if(fac_vism.ne.0.0_R8 .and. vis_per.ne.0.0e0_R8) then             !srv 22.11.99

       do is = 0,ns-1                                                   !srv 16.10.17
         if((visper_all_ions.ne.0 .and. .not.is_neutral(is)) .or.
     &      (visper_all_ions.eq.0 .and. is.eq.ismain)) then
          do iy=-1,ny
            do ix=-1,nx
              if(bottomiy(ix,iy).eq.-2)then
                cvisper_a(ix,iy,1,is) = 0.0_R8
              elseif(leftix(ix,iy).eq.-2.or.rightix(ix,iy).eq.nx+1) then
                cvisper_a(ix,iy,1,is) = 0.0_R8   ! we do not want to have fchvisper_a(:,:,1) in guard cells on West and East
              else
                t0 = 2.0_R8/wbbv(ix,iy,3)/
     &            (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
                cvisper_a(ix,iy,1,is) = t0*vis_per*fac_vis(ix,iy)
              endif
            enddo
          enddo
        
          do iy=-1,ny
            do ix=-1,nx
              if(leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1) then ! we do not want to have fchvisper_a(:,:,1) in guard cells on West and East
                t0 = 0.25_R8*vsa0(ix,iy,is)*vol(ix,iy)/hy1(ix,iy)**2
                bvisper_a(ix,iy,1,is) = t0*vis_per*fac_vis(is,iy)
              endif
            enddo
          enddo
          
          avisper_a(:,:,1,is) = cvisper_a(:,:,1,is)

         endif                                                          !srv 16.10.17
       enddo                                                            !srv 16.10.17
      endif                                                             !srv 22.11.99
      
      if(iout.eq.1) then
       do iy=-1,ny
        do ix=-1,nx
          wrk0(ix,iy)=hy1(ix,iy)
        enddo
       enddo
       call my_out(70,nx,ny,wrk0,'hy1') 
       call my_out(70,nx,ny,wbbv(-1,-1,3),'b2tvspr_bbv')
       call my_out(70,nx,ny,vol,'b2tvspr_vol')
       do is = 0,ns-1
         write (chns,'(i3.3)') is
         if((visper_all_ions.ne.0 .and. .not.is_neutral(is)) .or.
     &      (visper_all_ions.eq.0 .and. is.eq.ismain)) then
          call my_out(70,nx,ny,avisper_a(-1,-1,1,is),
     &    'b2txvspr_avisper_ay'//chns)
          call my_out(70,nx,ny,bvisper_a(-1,-1,1,is),
     &    'b2txvspr_bvisper_ay'//chns)
          call my_out(70,nx,ny,cvisper_a(-1,-1,1,is),
     &    'b2txvspr_cvisper_ay'//chns)
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

