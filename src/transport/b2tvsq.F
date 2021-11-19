      subroutine b2tvsq (nx, ny, ns, ismain, fac_vis,                    !srv 19.02.21
     &           ne2, lnlam, te, ti, ni, chci, fchvisq_a, fchvisq)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_rates
      use b2mod_b2cmpa
      use b2mod_geo
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain                                         !srv 19.02.21
      real (kind=R8) ::
     &   fac_vis(-1:nx,-1:ny),
     &   ne2(-1:nx,-1:ny), lnlam(-1:nx,-1:ny), te(-1:nx,-1:ny),
     &   ti(-1:nx,-1:ny), ni(-1:nx,-1:ny), chci(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   fchvisq_a(-1:nx,-1:ny,0:1,0:ns-1), fchvisq(-1:nx,-1:ny,0:1)    !srv 11.07.99
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TVSQ computes the X and Y components of the heat viscosity
*     current.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix, iy, is
      integer mintopcut                                                 !lk
      integer, save :: ncall = 0
      integer, save :: qstyle = 0                                       !lk 23.05.06 05.07.07 {
      integer, save :: no_current = 0
      real (kind=R8), save :: vis_q = 1.0e0_R8                          !srv 11.07.99
      real (kind=R8), save :: cqip1 = 0.5e0_R8                          !lk
      real (kind=R8), save :: cvsa_mltpl = 1.0e0_R8                     !iyv 06.12.13
      integer, save :: drift_style                                      !IYS 26.05.2020
      real (kind=R8) :: tt0(-1:ny), tt1(-1:ny), bbaverage2(-1:ny)       !lk
     &                  , conn_(-1:ny), tt2(-1:ny), tt3(-1:ny)
     &                  , eps1(-1:ny), fac_vism                         !srv 16.10.17
      real (kind=R8) :: qip0a(-1:nx,-1:ny), qip0b(-1:nx,-1:ny)
     &                  , tauia1(-1:nx,-1:ny), alpha1(-1:nx,-1:ny)
     &                  , nu1(-1:nx,-1:ny)                              !lk 23.05.06 05.07.07 }
      real (kind=R8) ::
     *    wrk0(-1:nx,-1:ny), qip0(-1:nx,-1:ny),                         !srv 11.07.99
     *    qip1(-1:nx,-1:ny), qipgen(-1:nx,-1:ny),                       !srv 20.12.99
     *    tauia(-1:nx,-1:ny,0:ns-1)                                     !srv 28.03.99 11.07.99
*   ..procedures
      intrinsic abs
      external xertst
      external b2xvsg, b2xvff, b2xvfx
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
       if (redef_gmtry.eq.0) then
          drift_style = 1
       else
          drift_style = 2
       end if
       call ipgeti ('b2tfnb_drift_style', drift_style)                  !srv 21.06.08
!       write(*,*) 'redef_gmtry = ',redef_gmtry,
!     &               ' drfit_style = ',drift_style
       call ipgeti ('b2tfhe_no_current', no_current)
      end if
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

      fchvisq   = 0.0e0_R8                                              !srv 11.07.99
      fchvisq_a = 0.0e0_R8
      fac_vism = maxval(fac_vis)
      if(fac_vism.ne.0.0_R8 .and. vis_q.ne.0.0e0_R8 .and.               !srv 16.10.17
     &   no_current.eq.0) then
         do is = 0,ns-1                                                 !srv 12.05.11
!          if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then          !srv 12.05.11
          if(is.eq.ismain) then                                         !srv 19.02.21   !
           qip0  = 0.0e0_R8
           qip1  = 0.0e0_R8
           qipgen= 0.0e0_R8
           qip0a = 0.0e0_R8
           qip0b = 0.0e0_R8
           wrk0 = 0.0e0_R8

           call b2ttia (nx, ny, ns, ti, ne2, lnlam, tauia)              !srv 11.07.99

           tt0 = 0.0e0_R8
           tt1 = 0.0e0_R8
           tt2 = 0.0e0_R8
           tt3 = 0.0e0_R8
           conn_ = 0.0e0_R8
           bbaverage2 = 0.0e0_R8
           eps1 = 0.0e0_R8

!pbwdk           mintopcut = ny
           mintopcut = -2
           if (nncut > 0) mintopcut = minval(topcut(1:nncut))

           do iy = 0,mintopcut    ! for inner surfaces and first open surface
            do ix = 0,nx
             if(on_closed_surface(bottomix(ix,iy),bottomiy(ix,iy))) then
              tt0(iy) = tt0(iy) +bb(ix,iy,3)*bb(ix,iy,3)*vol(ix,iy)
              tt1(iy) = tt1(iy) +vol(ix,iy)
              conn_(iy) =
     &            conn_(iy)+hx(ix,iy)*abs(bb(ix,iy,3)/bb(ix,iy,0))      !srv 11.01.13
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
     &          *(ni(ix,iy)*ti(ix,iy)/qe)*
     &          (ti(topix(ix,iy),topiy(ix,iy))
     &          -ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &          (hy1(ix,iy)+0.5e0_R8*hy1(topix(ix,iy),topiy(ix,iy))+
     &          0.5e0_R8*hy1(bottomix(ix,iy),bottomiy(ix,iy)))
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
               qip0a(ix,iy)=-chci(ix,iy,0)
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
     &          *(ni(ix,iy)*ti(ix,iy)/qe)*
     &          (ti(topix(ix,iy),topiy(ix,iy))
     &          -ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &          (hy1(ix,iy)+0.5e0_R8*hy1(topix(ix,iy),topiy(ix,iy))+
     &          0.5e0_R8*hy1(bottomix(ix,iy),bottomiy(ix,iy)))
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
             wrk0(ix,iy)=qip1(ix,iy)*atan(qip0(ix,iy)/qip1(ix,iy))
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
            do ix = 0,nx
             if(on_closed_surface(bottomix(ix,iy),
     .                            bottomiy(ix,iy))) then
              tauia1(ix,iy) = cvsa_mltpl*tauia(ix,iy,is)/               !lk 20.11.07 { !iyv 06.12.13
     &         (1.0e0_R8+15.12e16_R8/conn_(iy)*
     &         (ti(ix,iy)/ev)**2/ne2(ix,iy))
     &         /(1.0e0_R8+15.12e16_R8*eps1(iy)**1.5e0_R8
     &         /conn_(iy)*(ti(ix,iy)/ev)**2/ne2(ix,iy))
              nu1(ix,iy) = conn_(iy)/(tauia(ix,iy,is)
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

           if (drift_style .le. 1) then                                 !IYS 26.05.2020
             do iy = 0,mintopcut-1                                      !lk 23.05.06 05.07.07 {
               do ix = 0,nx
                 if(on_closed_surface(bottomix(ix,iy),
     .                                bottomiy(ix,iy))) then
                   fchvisq_a(ix,iy,0,is) =                               !srv 16.10.17
     &         -fac_vis(ix,iy)*0.24e0_R8*vis_q*(             ! signmf*   !srv 11.01.13 16.10.17
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
     &         *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))
                 endif                                                   !lk 23.05.06 05.07.07 }
               enddo
             enddo

           elseif (drift_style .eq. 2) then                              !IYS 26.05.2020

             do iy = 0,mintopcut-1                                       !lk 23.05.06 05.07.07 {
               do ix = 0,nx
                 if(on_closed_surface(bottomix(ix,iy),
     .                                bottomiy(ix,iy))) then
                   fchvisq_a(ix,iy,0,is) =
     &         -fac_vis(ix,iy)*0.24e0_R8*vis_q*(            ! signmf*    !srv 11.01.13 16.10.17
     &         tauia1(ix,iy)*alpha1(ix,iy)/hx(ix,iy)/hy(ix,iy)
     &         +tauia1(leftix(ix,iy),leftiy(ix,iy))
     &         *alpha1(leftix(ix,iy),leftiy(ix,iy))
     &         /hx(leftix(ix,iy),leftiy(ix,iy))**2)*0.5e0_R8
     &         *(wrk0(ix,iy)-wrk0(leftix(ix,iy),leftiy(ix,iy)))
     &         *wbbl(ix,iy,2)/wbbl(ix,iy,3)*wbbl(ix,iy,0)
     &         /sqrt(wbbl(ix,iy,3))*gs(ix,iy,0)*qc(ix,iy)
     &         *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))
                 endif                                                   !lk 23.05.06 05.07.07 }
               enddo
             enddo
           endif


* ..calculate the y-component of the parallel viscosity current
*                                  connected with the heat flux
          if (drift_style .le. 1) then                                   !IYS 26.05.2020
            do iy = 0,mintopcut                                          !lk 23.05.06 05.07.07 {
              do ix = 0,nx
                if(on_closed_surface(bottomix(ix,iy),
     .                               bottomiy(ix,iy))) then
                  if(rightix(ix,iy).ne.nx+1) then
                    fchvisq_a(ix,iy,1,is) =
     &          fac_vis(ix,iy)*0.24e0_R8*vis_q*(             ! signmf*   !srv 11.01.13 16.10.17
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
     &          *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))
     &          -OnedBsqc(ix,iy))
                  endif
                endif                                                    !lk 23.05.06 05.07.07 }
              enddo
            enddo

           elseif (drift_style .eq. 2) then                              !IYS 26.05.2020

            do iy = 0,mintopcut                                          !lk 23.05.06 05.07.07 {
             do ix = 0,nx
              if(on_closed_surface(bottomix(ix,iy),
     .                               bottomiy(ix,iy))) then
               if(rightix(ix,iy).ne.nx+1) then
                fchvisq_a(ix,iy,1,is) =
     &          fac_vis(ix,iy)*0.24e0_R8*vis_q*(            ! signmf*    !srv 11.01.13 16.10.17
     &          tauia1(ix,iy)*alpha1(ix,iy)/hx(ix,iy)**2
     &          +tauia1(bottomix(ix,iy),bottomiy(ix,iy))
     &          *alpha1(bottomix(ix,iy),bottomiy(ix,iy))
     &          /hx(bottomix(ix,iy),bottomiy(ix,iy))**2
     &                      )*0.5e0_R8
     &          *(qipgen(rightix(ix,iy),rightiy(ix,iy))
     &          +qipgen(rightix(bottomix(ix,iy),bottomiy(ix,iy)),
     &             rightiy(bottomix(ix,iy),bottomiy(ix,iy)))
     &          -qipgen(ix,iy)-qipgen(bottomix(ix,iy),bottomiy(ix,iy))
     &            )*0.5e0_R8
     &          *wbbv(ix,iy,2)/wbbv(ix,iy,3)
     &          *wbbv(ix,iy,0)/sqrt(wbbv(ix,iy,3))*gs(ix,iy,1)
     &          *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))
     &          -OnedBsqc(ix,iy))
               endif
              endif                                                      !lk 23.05.06 05.07.07 }
             enddo
            enddo

           endif

          else
            fchvisq_a(:,:,:,is)=0.0_R8
          endif
          fchvisq(:,:,:) = fchvisq(:,:,:) + fchvisq_a(:,:,:,is)

         enddo
      endif
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvsq

      end subroutine b2tvsq
