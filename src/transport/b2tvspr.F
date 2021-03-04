      subroutine b2tvspr(nx, ny, ns, fac_vis, vol,                      !srv 10.10.17
     &   po, ti, na, avisper, bvisper, cvisper, fchvisper)
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
     1   fac_vis(-1:nx,-1:ny), vol(-1:nx,-1:ny),
     1   po(-1:nx,-1:ny), ti(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     3   avisper(-1:nx,-1:ny,0:1), bvisper(-1:nx,-1:ny,0:1),            !srv 27.01.00
     3   cvisper(-1:nx,-1:ny,0:1), fchvisper(-1:nx,-1:ny,0:1)
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix, iy, is
      integer, save :: ncall = 0, iout = 0                              !srv 06.12.17
      real (kind=R8), save :: vis_per = 1.0e0_R8                        !srv 22.11.99
      real (kind=R8) :: va(-1:nx,-1:ny,0:1), dva(-1:nx,-1:ny)           !srv 22.11.99
      real (kind=R8) :: fac_vism
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tvspr')
*   ..set internal parameters on first call
      if (ncall.eq.0) then                                              !srv
       call ipgetr ('b2tfhe_vis_per', vis_per)                          !srv 22.11.99
       call ipgeti ('b2tvspr_iout', iout)                               !srv 06.12.17
      end if
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..calculate the coefficients in difference approximation of viscosity perpendicular current
      call b2txvspr(nx, ny, ns, fac_vis, vol,                           !srv 10.10.17
     &              avisper, bvisper, cvisper)
* ..calculate the contribution of viscosity perpendicular current       !srv 22.11.99 
      va=0.0e0_R8
      dva=0.0e0_R8
      fchvisper =0.0e0_R8                                               !srv 22.11.99
        
      fac_vism = maxval(fac_vis)
      if(fac_vism.ne.0.0_R8 .and. vis_per.ne.0.0e0_R8) then             !srv 22.11.99

      do is = 0,ns-1                                                    !srv 16.10.17
        if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then             !srv 16.10.17
          do iy=-1,ny
            do ix=-1,nx
              if(bottomiy(ix,iy).eq.-2)then
                va(ix,iy,1) = 0.0_R8
              elseif(leftix(ix,iy).eq.-2.or.rightix(ix,iy).eq.nx+1) then
                va(ix,iy,1) = 0.0_R8     ! we do not want to have va in guard cells on West and East
              else
                va(ix,iy,1) = -cvisper(ix,iy,1)*                      ! va is calculated on bottom cell face
     &            (po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy))+
     &            2.0_R8/qe*(na(ix,iy,is)*ti(ix,iy)-
     &            na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &            ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &            (rza(ix,iy,is)*na(ix,iy,is)+
     &            rza(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &            na(bottomix(ix,iy),bottomiy(ix,iy),is)))
              endif
            enddo
          enddo

          do iy=-1,ny
            do ix=-1,nx
              if(bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1) then
!               dva(ix,iy) will be calculated below in order to take into account boundary conditions
              elseif(leftix(ix,iy).eq.-2.or.rightix(ix,iy).eq.nx+1.or.
     &             bottomiy(ix,iy).eq.-2.or.topiy(ix,iy).eq.ny+1) then
                dva(ix,iy) = 0.0_R8                ! we do not want to have dva in guard cells on West and East
              else
                dva(ix,iy)=bvisper(ix,iy,1)*       ! dva is calculated at the cell centre
     &          (va(topix(ix,iy),topiy(ix,iy),1)-va(ix,iy,1))
              endif
            enddo
          enddo
*   .. take into account the consequence of the boundary conditions fchvisper = 0 at North and South
          do iy=-1,ny
            do ix=-1,nx
              if(bottomiy(ix,iy).eq.-2) then
                dva(ix,iy) = dva(topix(ix,iy),topiy(ix,iy))
              endif
              if(topiy(ix,iy).eq.ny+1) then
                dva(ix,iy) = dva(bottomix(ix,iy),bottomiy(ix,iy))
              endif
            enddo
          enddo

          do iy=-1,ny
            do ix=-1,nx
              if(bottomiy(ix,iy).eq.-2) then
                fchvisper(ix,iy,1) = 0.0_R8
              elseif(leftix(ix,iy).eq.-2.or.rightix(ix,iy).eq.nx+1) then
                fchvisper(ix,iy,1) = 0.0_R8   ! we do not want to have fchvisper(:,:,1) in guard cells on West and East
              else
                fchvisper(ix,iy,1) = fchvisper(ix,iy,1) -
     &            avisper(ix,iy,1)*
     &           (dva(ix,iy)-dva(bottomix(ix,iy),bottomiy(ix,iy)))
              endif
            enddo
          enddo

*          X-COMPONENT OF VISCOSITY PERPENDICULAR CURRENT               !srv 22.11.99
          fchvisper(:,:,0)=fchvisper(:,:,0)+0.0e0_R8
        endif                                                           !srv 16.10.17
      enddo                                                             !srv 16.10.17
      endif                                                             !srv 22.11.99

      if(iout.eq.1) then
       call my_out(70,nx,ny,fchvisper(-1,-1,0),'b2tvspr_fchvisperx')
       call my_out(70,nx,ny,fchvisper(-1,-1,1),'b2tvspr_fchvispery')
      endif
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvspr

      end subroutine b2tvspr                                            !srv 10.10.17

