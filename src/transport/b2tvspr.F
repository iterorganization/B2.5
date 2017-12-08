      subroutine b2tvspr(nx, ny, ns, fac_vis, vol, hx, hy,              !srv 10.10.17
     &   qz, bb, po, ti, na, avisper, bvisper, cvisper, fchvisper)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_anomalous_transport
      use b2mod_rates
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     1   fac_vis(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     1   vol(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), bb(-1:nx,-1:ny,0:3), 
     1   po(-1:nx,-1:ny), ti(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     3   avisper(-1:nx,-1:ny,0:1), bvisper(-1:nx,-1:ny,0:1),            !srv 27.01.00
     $   cvisper(-1:nx,-1:ny,0:1), fchvisper(-1:nx,-1:ny,0:1)
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, ix, iy, is
      integer, save :: ncall = 0
      real (kind=R8), save :: vis_per = 0.0e0_R8                        !srv 22.11.99
      real (kind=R8) :: 
     *    va(-1:nx,-1:ny,0:1), dva(-1:nx,-1:ny), wrk0(-1:nx,-1:ny)      !srv 22.11.99 16.10.17
      real (kind=R8) ::
     &  fac_vism
*   ..procedures
      intrinsic max, abs
      real (kind=R8) ::
     *   hy1
      external subini, subend, xertst
      external b2xvsg, b2xvff, b2xvfx
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)                                !srv 07.05.08
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tvspr')
*   ..set internal parameters on first call
      if (ncall.eq.0) then                                              !srv
       call ipgetr ('b2tfhe_vis_per', vis_per)                          !srv 22.11.99
      end if
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
* ..calculate the contribution of viscosity perpendicular current       !srv 22.11.99 
      va=0.0e0_R8
      dva=0.0e0_R8
      fchvisper =0.0e0_R8                                               !srv 22.11.99
      avisper  = 0.0e0_R8; bvisper = 0.0e0_R8; cvisper = 0.0e0_R8       !srv 27.01.00 22.03.02
        
      fac_vism = maxval(fac_vis)
      if(fac_vism.ne.0.0_R8 .and. vis_per.ne.0.0e0_R8) then             !srv 22.11.99

!      analogous part for avisper and va takes place in b2siav          !srv 17.01.01

        avisper(:,-1,1) = 0.0e0_R8                                      !srv 06.05.02
        cvisper(:,-1,1) = 0.0e0_R8
        do iy=0,ny
          do ix=-1,nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then    !srv 26.11.02
              cvisper(ix,iy,1)=
     &         4.0e0_R8/(bb(ix,iy,3)+
     1         bb(bottomix(ix,iy),bottomiy(ix,iy),3))/
     &         (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
             avisper(ix,iy,1) = cvisper(ix,iy,1)                        !srv 06.05.02
            endif                                                       !srv 26.11.02
          enddo
        enddo

      do is = 0,ns-1                                                    !srv 16.10.17
       if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then              !srv 16.10.17
         do iy=-1,ny
          do ix=-1,nx
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then     !srv 26.11.02
            wrk0(ix,iy) =
     &             0.25e0_R8*vsa0(ix,iy,is)*vol(ix,iy)/hy1(ix,iy)**2    !srv 16.10.17
            bvisper(ix,iy,1) = bvisper(ix,iy,1) + wrk0(ix,iy)           !srv 16.10.17
!            pra(ix,iy) = na(ix,iy,is)*ti(ix,iy) 
           endif                                                        !srv 26.11.02
          enddo
         enddo

!          call intcor(nx,ny,vol,pra,prac)
!            call intcor(nx,ny,vol,na(-1,-1,is),nac)

         do iy=0,ny
          do ix=0,nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then     !srv 26.11.02
              va(ix,iy,0)=-cvisper(ix,iy,1)
     &        *(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))
!     &        -cdna(ix,iy,1,is)*2.0e0_R8*
!     &        (nac(ix,iy)-nac(rightix(ix,iy),rightiy(ix,iy)))/
!     &        (hx(ix,iy)*na(ix,iy,is)+
!     &        hx(bottomix(ix,iy),bottomiy(ix,iy))*
!     &        na(bottomix(ix,iy),bottomiy(ix,iy),is))
!     &        -cdpa(ix,iy,1,is)*2.0e0_R8*
!     &        (prac(ix,iy)-prac(rightix(ix,iy),rightiy(ix,iy)))/
!     &        (hx(ix,iy)*na(ix,iy,is)+
!     &        hx(bottomix(ix,iy),bottomiy(ix,iy))*
!     &        na(bottomix(ix,iy),bottomiy(ix,iy),is))
     &         -cvisper(ix,iy,1)*(na(ix,iy,is)*ti(ix,iy)
     &         -na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &         ti(bottomix(ix,iy),bottomiy(ix,iy)))*
     &         (2.0e0_R8/(rza(ix,iy,is)*qe))
     &         /(na(bottomix(ix,iy),bottomiy(ix,iy),is)+na(ix,iy,is))
           endif                                                        !srv 26.11.02
          enddo
         enddo

         do iy=0,ny-1
          do ix=0,nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then     !srv 26.11.02
            dva(ix,iy)=wrk0(ix,iy)*  ! bvisper(ix,iy,1)*                !srv 16.10.17
     *        (va(topix(ix,iy),topiy(ix,iy),0)-va(ix,iy,0))
           endif                                                        !srv 26.11.02
          enddo
         enddo

*          Y-COMPONENT OF VISCOSITY PERPENDICULAR CURRENT               !srv 22.11.99
         fchvisper(:,0,1) = 0.0e0_R8    ! because of boundary condition !srv 28.01.00
         do iy=1,ny-1    
          do ix=0,nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then     !srv 26.11.02
            fchvisper(ix,iy,1)=fchvisper(ix,iy,1)
     *       -avisper(ix,iy,1)*vis_per*fac_vis(ix,iy)*
     *       (dva(ix,iy)-dva(bottomix(ix,iy),bottomiy(ix,iy)))
           endif                                                        !srv 26.11.02
          enddo
         enddo
         fchvisper(:,ny,1) = 0.0e0_R8   ! because of boundary condition !srv 28.01.00
*          X-COMPONENT OF VISCOSITY PERPENDICULAR CURRENT               !srv 22.11.99
         fchvisper(:,:,0)=fchvisper(:,:,0)+0.0e0_R8
          endif                                                         !srv 16.10.17
         enddo                                                          !srv 16.10.17
        endif                                                           !srv 22.11.99
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvspr

      end subroutine b2tvspr                                            !srv 10.10.17

