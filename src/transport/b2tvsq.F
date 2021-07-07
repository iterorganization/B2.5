      subroutine b2tvsq (nx, ny, ns, fac_vis, vol, hx, hy, qz, qs, bb,
     &       ne2, lnlam, te, ti, ni, chci, fchvisq)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_rates
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     &   fac_vis(-1:nx,-1:ny), vol(-1:nx,-1:ny), hx(-1:nx,-1:ny),
     &   hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), bb(-1:nx,-1:ny,0:3),
     &   ne2(-1:nx,-1:ny), lnlam(-1:nx,-1:ny), te(-1:nx,-1:ny),
     &   ti(-1:nx,-1:ny), ni(-1:nx,-1:ny), chci(-1:nx,-1:ny,0:1,0:1), 
     &   qs(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     $   fchvisq(-1:nx,-1:ny,0:1,0:1)                                   !srv 11.07.99
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, ix, iy, is
      integer mintopcut                                                 !lk
      integer, save :: ncall = 0
      integer, save :: qstyle = 0                                       !lk 23.05.06 05.07.07 {
      real (kind=R8), save :: vis_q = 1.0e0_R8                          !srv 11.07.99
      real (kind=R8), save :: cqip1 = 0.5e0_R8                          !lk
      real (kind=R8), save :: cvsa_mltpl = 1.0e0_R8                     !iyv 06.12.13
      real (kind=R8) :: tt0(-1:ny), tt1(-1:ny), bbaverage2(-1:ny)       !lk
     &                  , conn(-1:ny), tt2(-1:ny), tt3(-1:ny)
     &                  , eps1(-1:ny), fac_vism                         !srv 16.10.17
      real (kind=R8) :: qip0a(-1:nx,-1:ny), qip0b(-1:nx,-1:ny)
     &                  , tauia1(-1:nx,-1:ny), alpha1(-1:nx,-1:ny)
     &                  , nu1(-1:nx,-1:ny)                              !lk 23.05.06 05.07.07 }
      real (kind=R8) ::
     *    wrk0(-1:nx,-1:ny), qip0(-1:nx,-1:ny),                         !srv 11.07.99
     *    qip1(-1:nx,-1:ny), qipgen(-1:nx,-1:ny),                       !srv 20.12.99
     *    tauia(-1:nx,-1:ny,0:ns-1),                                    !srv 28.03.99 11.07.99
     *    dtir(-1:nx,-1:ny), wght(-1:nx,-1:ny),
     *    dOnedBsq(-1:nx,-1:ny,0:1,0:1)
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
      call subini ('b2tvsq')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tfhe_vis_q', vis_q)                              !srv 11.07.99
       call ipgeti ('b2siav_style_qip', qstyle)                         !lk 23.05.06 05.07.07
       call ipgetr ('b2trcl_cvsa_mltpl', cvsa_mltpl)                    !iyv 06.12.13
      end if
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      
      fchvisq   = 0.0e0_R8                                              !srv 11.07.99
      fac_vism = maxval(fac_vis)
      if(fac_vism.ne.0.0_R8 .and. vis_q.ne.0.0e0_R8) then               !srv 16.10.17

*     ..compute radial gradients of ti
#ifdef WG_TODO
         call gradc_r(nx, ny, hx, hy, vol, qz, ti, dtir)

*       compute differences of OnedBsq on cell faces
        call diff(nx, ny, hx, hy, vol, qs, OnedBsq, dOnedBsq)
#endif
!        wght = 1.0_R8
!        call diff_p(nx,ny, 1, hx, hy, wght, qs, OnedBsq, OnedBsqc, 
!     .    dOnedBsq(-1,-1,0,0))
!        call diff_r(nx,ny, 1, hx, hy, wght, qs, OnedBsq, OnedBsqc, 
!     .    dOnedBsq(-1,-1,0,1))

         do is = 0,ns-1                                                 !srv 12.05.11
          if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then           !srv 12.05.11
!           is = ismain    ! 1                                  
           qip0  = 0.0e0_R8
           qip1  = 0.0e0_R8
           qipgen= 0.0e0_R8
           qip0a = 0.0e0_R8
           qip0b = 0.0e0_R8
           wrk0 = 0.0e0_R8

!WG_TODO           call b2ttia (nCv, ns, ti, rz2, ne2, lnlam, tauia)
           
           tt0 = 0.0e0_R8                                                 
           tt1 = 0.0e0_R8
           tt2 = 0.0e0_R8                                                 
           tt3 = 0.0e0_R8
           conn = 0.0e0_R8
           bbaverage2 = 0.0e0_R8
           eps1 = 0.0e0_R8

!pbwdk           mintopcut = ny
           mintopcut = -2
           if (nncut > 0) mintopcut = minval(topcut(1:nncut))
        
           do iy = 0,mintopcut    ! for inner surfaces and first open surface
            do ix=0,nx                                                  
             if(on_closed_surface(bottomix(ix,iy),bottomiy(ix,iy))) then
              tt0(iy) = tt0(iy) +bb(ix,iy,3)*bb(ix,iy,3)*vol(ix,iy)
              tt1(iy) = tt1(iy) +vol(ix,iy)
              conn(iy) = conn(iy)+hx(ix,iy)*abs(bb(ix,iy,3)/bb(ix,iy,0)) !srv 11.01.13
              tt2(iy) = tt2(iy) +hx(ix,iy)
              tt3(iy) = tt3(iy) +vol(ix,iy)/hy1(ix,iy)
             endif
            enddo
           enddo

           do iy = 0,mintopcut
            bbaverage2(iy) = tt0(iy)/tt1(iy)
            eps1(iy) = tt2(iy)**2/tt3(iy)
           enddo                                                        
         
* ..calculate the parallel heat flux of ion at the center of cell       !lk 23.05.06 05.07.07 {
           if (qstyle .eq. 0) then

            do iy = 0,mintopcut-1    ! for inner surfaces
             do ix = 0,nx                                                  
              if(on_closed_surface(bottomix(ix,iy),
     .                             bottomiy(ix,iy))) then
               qip0a(ix,iy) = (-2.5e0_R8)*(bb(ix,iy,3)*bb(ix,iy,2)
     &          /bb(ix,iy,0)/bbaverage2(iy))
     &          *(ni(ix,iy)*ti(ix,iy)/qe)*dtir(ix,iy)    
              endif
             enddo
            enddo

            do iy = 0,mintopcut-1
             do ix = 0,nx                                                  
              if(on_closed_surface(bottomix(ix,iy),
     .                             bottomiy(ix,iy))) then
               if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
                 qip0b(ix,iy) = 0.25e0_R8*(qip0a(ix,iy)*2.0e0_R8
     &            +qip0a(leftix(ix,iy),leftiy(ix,iy))
     &            +qip0a(rightix(ix,iy),rightiy(ix,iy)))
               endif
              endif
             enddo
            enddo

            do iy = 0,mintopcut-1
             do ix = 0,nx
              if(on_closed_surface(bottomix(ix,iy),
     .                             bottomiy(ix,iy))) then
               if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
                 qip0(ix,iy) = 0.25e0_R8*(qip0b(ix,iy)*2.0e0_R8
     &            +qip0b(leftix(ix,iy),leftiy(ix,iy))
     &            +qip0b(rightix(ix,iy),rightiy(ix,iy)))
               endif
              endif
             enddo
            enddo
       
           endif

           if (qstyle.eq.1) then

            qip0a = 0.0e0_R8
            do iy=0,ny-1
             do ix=-1,nx
              if(leftix(ix,iy).ne.-2) then
               qip0a(ix,iy)=-chci(ix,iy,0,0)
     &          *(ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))
              endif                                                      
             enddo                                                       
            enddo 
         
            do iy=0,ny-1
             do ix=-1,nx-1
              if(on_closed_surface(bottomix(ix,iy),
     .                             bottomiy(ix,iy))) then
               qip0(ix,iy)=0.5e0_R8*(qip0a(ix,iy)+
     &          qip0a(rightix(ix,iy),rightiy(ix,iy)))
     &          *hx(ix,iy)*bb(ix,iy,3)/bb(ix,iy,0)/vol(ix,iy)
     &          -2.5e0_R8*bb(ix,iy,2)
     &          /bb(ix,iy,0)/bb(ix,iy,3)
     &          *(ni(ix,iy)*ti(ix,iy)/qe)*dtir(ix,iy) 
              endif
             enddo
            enddo
                                                                  
           endif
        
           do iy=0,ny-1
            do ix=0,nx-1
             qip1(ix,iy)=cqip1*ni(ix,iy)*ti(ix,iy)*
     &        sqrt(te(ix,iy)+ti(ix,iy))/sqrt(am(is)*mp)
            enddo
           enddo

           do iy=0,ny-1
            do ix=0,nx-1
             wrk0(ix,iy)=qip1(ix,iy)*datan(qip0(ix,iy)/qip1(ix,iy))
     &        *sqrt(bb(ix,iy,3))
            enddo
           enddo

           do iy = 0,mintopcut     
            do ix = 0,nx
             if(on_closed_surface(bottomix(ix,iy),
     .                            bottomiy(ix,iy))) then
              qipgen(ix,iy)=(wrk0(ix,iy)
     &         +wrk0(leftix(ix,iy),leftiy(ix,iy)))*0.5e0_R8
             endif                                                      
            enddo
           enddo                

           tauia1 = 0.0e0_R8
           nu1 = 0.0e0_R8
           alpha1 = 0.0e0_R8

           do iy = 0,mintopcut
            do ix=0,nx
             if(on_closed_surface(bottomix(ix,iy),
     .                            bottomiy(ix,iy))) then
	       tauia1(ix,iy) = cvsa_mltpl*tauia(ix,iy,is)/              !lk 20.11.07 { !iyv 06.12.13
     &         (1.0e0_R8+15.12e16_R8/conn(iy)*
     &         (ti(ix,iy)/ev)**2/ne2(ix,iy))
     &         /(1.0e0_R8+15.12e16_R8*eps1(iy)**1.5e0_R8
     &         /conn(iy)*(ti(ix,iy)/ev)**2/ne2(ix,iy))
              nu1(ix,iy) = conn(iy)/(tauia(ix,iy,is)
     &         *sqrt(ti(ix,iy)/(mp*am(is)))
     &         *2.0e0_R8*pi*eps1(iy)**1.5e0_R8) 
              alpha1(ix,iy) = ((-0.17e0_R8+1.05e0_R8*sqrt(nu1(ix,iy))+
     &         2.7e0_R8*nu1(ix,iy)**2*eps1(iy)**3)/(1.0e0_R8+
     &         0.7e0_R8*sqrt(nu1(ix,iy))+nu1(ix,iy)**2*eps1(iy)**3)
     &         -1.0e0_R8)
     &         *0.533e0_R8                                              !lk 20.11.07 }
             endif    
            enddo
           enddo                                                        !lk 23.05.06 05.07.07 }

* ..calculate the x-component of the parallel viscosity current
*                                  connected with the heat flux 

           do iy = 0,mintopcut-1                                        !lk 23.05.06 05.07.07 {
            do ix = 0,nx
             if(on_closed_surface(bottomix(ix,iy),
     .                            bottomiy(ix,iy))) then
              fchvisq(ix,iy,0,0) = fchvisq(ix,iy,0,0)-                  !srv 16.10.17
     &         fac_vis(ix,iy)*0.24e0_R8*vis_q*(             ! signmf*   !srv 11.01.13 16.10.17
     &         vol(ix,iy)/hx(ix,iy)**2/hy1(ix,iy)
     &         *tauia1(ix,iy)
     &         *alpha1(ix,iy)
     &         *bb(ix,iy,2)/bb(ix,iy,3)
     &         *bb(ix,iy,0)/sqrt(bb(ix,iy,3))
     &         +vol(leftix(ix,iy),leftiy(ix,iy))
     &         /hx(leftix(ix,iy),leftiy(ix,iy))**2
     &         /hy1(leftix(ix,iy),leftiy(ix,iy))
     &         *tauia1(leftix(ix,iy),leftiy(ix,iy))
     &         *alpha1(leftix(ix,iy),leftiy(ix,iy))
     &         *bb(leftix(ix,iy),leftiy(ix,iy),2)
     &         /bb(leftix(ix,iy),leftiy(ix,iy),3)
     &         *bb(leftix(ix,iy),leftiy(ix,iy),0)
     &         /sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3))
     &         )*0.5e0_R8
     &         *(wrk0(ix,iy)-wrk0(leftix(ix,iy),leftiy(ix,iy)))
     &         *dOnedBsq(ix,iy,0,1)
cwdk computation of radial current on 'poloidal' faces
cwdk somewhat messy interpolation of geometric factors, for exact
cwdk cancelation of some terms in poloidal and radial current components
              fchvisq(ix,iy,0,1) = fchvisq(ix,iy,0,1)-                  !srv 16.10.17
     &         fac_vis(ix,iy)*0.24e0_R8*vis_q*(             ! signmf*   !srv 11.01.13 16.10.17
     &         vol(ix,iy)/hx(ix,iy)**2/hy1(ix,iy)
     &         *tauia1(ix,iy)
     &         *alpha1(ix,iy)
     &         *bb(ix,iy,2)/bb(ix,iy,3)
     &         *bb(ix,iy,0)/sqrt(bb(ix,iy,3))
     &         +vol(leftix(ix,iy),leftiy(ix,iy))
     &         /hx(leftix(ix,iy),leftiy(ix,iy))**2
     &         /hy1(leftix(ix,iy),leftiy(ix,iy))
     &         *tauia1(leftix(ix,iy),leftiy(ix,iy))
     &         *alpha1(leftix(ix,iy),leftiy(ix,iy))
     &         *bb(leftix(ix,iy),leftiy(ix,iy),2)
     &         /bb(leftix(ix,iy),leftiy(ix,iy),3)
     &         *bb(leftix(ix,iy),leftiy(ix,iy),0)
     &         /sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3))
     &         )*0.5e0_R8
     &         *(wrk0(ix,iy)-wrk0(leftix(ix,iy),leftiy(ix,iy)))
     &         *dOnedBsq(ix,iy,0,0)*qs(ix,iy,0)*
     &         (hy(ix,iy) + hy(leftix(ix,iy),leftiy(ix,iy)))/
     &         (hy1(ix,iy) + hy1(leftix(ix,iy),leftiy(ix,iy)))
             endif                                                      !lk 23.05.06 05.07.07 }
            enddo
           enddo

* ..calculate the y-component of the parallel viscosity current
*                                  connected with the heat flux 
           do iy = 0,mintopcut                                          !lk 23.05.06 05.07.07 {
            do ix = 0,nx
             if(on_closed_surface(bottomix(ix,iy),
     .                            bottomiy(ix,iy))) then
              if(rightix(ix,iy).ne.nx+1) then
               fchvisq(ix,iy,1,1) = fchvisq(ix,iy,1,1)+                 !srv 16.10.17 
     &          fac_vis(ix,iy)*0.24e0_R8*vis_q*(            ! signmf*   !srv 11.01.13 16.10.17
     &          vol(ix,iy)/hx(ix,iy)**2/hy1(ix,iy)
     &          *tauia1(ix,iy)
     &          *alpha1(ix,iy)
     &          *bb(ix,iy,2)/bb(ix,iy,3)
     &          *bb(ix,iy,0)/sqrt(bb(ix,iy,3))
     &          +vol(bottomix(ix,iy),bottomiy(ix,iy))
     &          /hx(bottomix(ix,iy),bottomiy(ix,iy))**2
     &          /hy1(bottomix(ix,iy),bottomiy(ix,iy))
     &          *tauia1(bottomix(ix,iy),bottomiy(ix,iy))
     &          *alpha1(bottomix(ix,iy),bottomiy(ix,iy))
     &          *bb(bottomix(ix,iy),bottomiy(ix,iy),2)
     &          /bb(bottomix(ix,iy),bottomiy(ix,iy),3)
     &          *bb(bottomix(ix,iy),bottomiy(ix,iy),0)
     &          /sqrt(bb(bottomix(ix,iy),bottomiy(ix,iy),3))
     &                      )*0.5e0_R8
     &          *(qipgen(rightix(ix,iy),rightiy(ix,iy))
     &          +qipgen(rightix(bottomix(ix,iy),bottomiy(ix,iy)),
     &             rightiy(bottomix(ix,iy),bottomiy(ix,iy)))
     &          -qipgen(ix,iy)-qipgen(bottomix(ix,iy),bottomiy(ix,iy))
     &            )*0.5e0_R8
     &          *dOnedBsq(ix,iy,1,0)
              endif
             endif                                                      !lk 23.05.06 05.07.07 }
            enddo
           enddo
          endif
         enddo
      endif
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvsq

      end subroutine b2tvsq
