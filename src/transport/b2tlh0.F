      subroutine b2tlh0 (nx, ny, ns, hx, hy, vol, qz, na, ti,
     & hcib, chcb, hcimx, fllim0fhi)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     *  vol(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1),
     *  na(-1:nx,-1:ny,0:ns-1), ti(-1:nx,-1:ny),
     *  hcib(-1:nx,-1:ny,0:ns-1), chcb(-1:nx,-1:ny,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  hcimx(-1:nx,-1:ny,0:1,0:ns-1),
     &  fllim0fhi(-1:nx,-1:ny,0:1,0:ns-1)
*   ..common blocks
*     ------------------------------------------------------------------
*     B2TLH0 computes the conductive neutral heat transport flux limit.
*     ------------------------------------------------------------------
*   ..local variables
      integer ix, iy, is, n2, ncall, c3, hcimx_flag
      integer, save :: iout = 0                                   !srv 13.01.04 {
      real (kind=R8), allocatable :: w1(:,:)
      character charns*3                                          !srv 13.01.04 }
      integer, save :: flux_limit_style = 2                       !srv 13.06.08
      real (kind=R8) ::
     *  alpha, gamma, t0, t1, t2, t3, dti, vthi, maxfluxlimit(0:1),
     *  flux_limit_min_ti
      save ncall, alpha, gamma, flux_limit_min_ti, hcimx_flag
*   ..initialisation
      data ncall/0/, alpha/1.0_R8/, gamma/1.0_R8/,
     *     flux_limit_min_ti/0.0_R8/, hcimx_flag/1/
*   ..procedures
      intrinsic max, abs, sqrt
      real (kind=R8) ::
     *   hy1
      external subini, subend, xertst, ipgetr, ipgeti
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2tlh0')
*   ..test nx, ny, ns
      call xertst (0.le.nx .and. 0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..read input on first call
      if (ncall.eq.0) then
        call ipgetr ('b2tlh0_alpha', alpha)
        call xertst (0.0_R8.le.alpha, 'invalid input parameter alpha')
        call ipgetr ('b2tlh0_gamma', gamma)
        call xertst (0.0_R8.lt.gamma, 'invalid input parameter gamma')
        call ipgetr ('b2tlh0_flux_limit_min_ti', flux_limit_min_ti)
        flux_limit_min_ti=flux_limit_min_ti*ev
        call ipgeti ('b2tlh0_hcimx_flag', hcimx_flag)
        call ipgeti ('b2tlh0_iout', iout)                         !srv 13.01.04
        call ipgeti ('b2tlh0_flux_limit_style', flux_limit_style) !srv 13.06.08
      endif

      call sfill(n2*2*ns, 0.0_R8, hcimx, 1)
      call sfill(n2*2*ns, 1.0_R8, fllim0fhi, 1)
      maxfluxlimit=0.0_R8

*   ..nothing done if alpha.eq.0
      if (alpha.eq.0.0_R8) then
        hcimx(:,:,0,:) = hcib(:,:,:)
        hcimx(:,:,1,:) = hcib(:,:,:)
        goto 100
      endif

      if(flux_limit_style.eq.0) then                              !srv 13.06.08

       do ix = -1, nx
         do iy = -1, ny
*   .. compute poloidal neutral heat conductivity flux limit
           do is = 0, ns-1
             if (is_neutral(is)) then
               t3=0.0_R8
               c3=0
               if (leftix(ix,iy).ne.-2) then
                 t0 = max(flux_limit_min_ti,
     *                 ti(ix,iy), ti(leftix(ix,iy),leftiy(ix,iy)))
                 t1 = (na(ix,iy,is) +
     .                 na(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8
                 dti = abs(ti(leftix(ix,iy),leftiy(ix,iy)) - ti(ix,iy))
     .            /(0.5_R8*(hx(ix,iy)+
     .                      hx(leftix(ix,iy),leftiy(ix,iy))))  ! dpc 2002-12-17
                 vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))
                 t2 = (1.0_R8 +
     .            ((hcib(ix,iy,is)*dti)/(alpha*vthi*t1*t0))**gamma)
     .            **(1.0_R8/gamma)
                 t3=t3+t2
                 c3=c3+1
               endif
               if (rightix(ix,iy).ne.nx+1.and.
     1            (hcimx_flag.gt.0.or.
     2            (leftix(ix,iy).eq.-2.and.hcimx_flag.ge.0))) then
                 t0 = max(flux_limit_min_ti,
     *                 ti(ix,iy), ti(rightix(ix,iy),rightiy(ix,iy)))
                 t1 = (na(ix,iy,is) +
     .                 na(rightix(ix,iy),rightiy(ix,iy),is)) /2.0_R8
                 dti = abs(ti(rightix(ix,iy),rightiy(ix,iy)) -
     .                     ti(ix,iy))
     .            /(0.5_R8*(hx(ix,iy)+
     .                      hx(rightix(ix,iy),rightiy(ix,iy))))  ! dpc 2002-12-17
                 vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))
                 t2 = (1.0_R8 +
     .            ((hcib(ix,iy,is)*dti)/(alpha*vthi*t1*t0))**gamma)
     .            **(1.0_R8/gamma)
                 t3=t3+t2
                 c3=c3+1
               endif
               if(c3.gt.0) then
                 t3=t3/c3
                 hcimx(ix,iy,0,is) = hcib(ix,iy,is)/t3
                 maxfluxlimit(0)=max(maxfluxlimit(0),t3)
                 fllim0fhi(ix,iy,0,is) = t3
               else
                 hcimx(ix,iy,0,is) = hcib(ix,iy,is)
               endif
             else
               hcimx(ix,iy,0,is) = hcib(ix,iy,is)
             endif
*   .. compute radial neutral heat conductivity flux limit
             if (is_neutral(is)) then
               t3=0.0_R8
               c3=0
               if (bottomiy(ix,iy).ne.-2) then
                 t0 = max(flux_limit_min_ti,
     *                 ti(ix,iy), ti(bottomix(ix,iy),bottomiy(ix,iy)))
                 t1 = (na(ix,iy,is) +
     .                 na(bottomix(ix,iy),bottomiy(ix,iy),is))/2.0_R8
                 dti = abs(ti(bottomix(ix,iy),bottomiy(ix,iy))-
     .                     ti(ix,iy))
     .            /(0.5_R8*(hy(ix,iy)+
     .                      hy(bottomix(ix,iy),bottomiy(ix,iy)))) ! dpc 2002-12-17
                 vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))
                 t2 = (1.0_R8 +
     .            ((hcib(ix,iy,is)*dti)/(alpha*vthi*t1*t0))**gamma)
     .            **(1.0_R8/gamma)
                 t3=t3+t2
                 c3=c3+1
               endif
               if (topiy(ix,iy).ne.ny+1.and.
     1            (hcimx_flag.gt.0.or.
     2            (bottomiy(ix,iy).eq.-2.and.hcimx_flag.ge.0))) then
                 t0 = max(flux_limit_min_ti,
     *                 ti(ix,iy), ti(topix(ix,iy),topiy(ix,iy)))
                 t1 = (na(ix,iy,is) +
     .                 na(topix(ix,iy),topiy(ix,iy),is))/2.0_R8
                 dti = abs(ti(topix(ix,iy),topiy(ix,iy))-ti(ix,iy))
     .            /(0.5_R8*(hy(ix,iy)+
     .                      hy(topix(ix,iy),topiy(ix,iy)))) ! dpc 2002-12-17
                 vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))
                 t2 = (1.0_R8 +
     .            ((hcib(ix,iy,is)*dti)/(alpha*vthi*t1*t0))**gamma)
     .            **(1.0_R8/gamma)
                 t3=t3+t2
                 c3=c3+1
               endif
               if(c3.gt.0) then
                 t3=t3/c3
                 hcimx(ix,iy,1,is) = hcib(ix,iy,is)/t3
                 maxfluxlimit(1)=max(maxfluxlimit(1),t3)
                 fllim0fhi(ix,iy,1,is) = t3
               else
                 hcimx(ix,iy,1,is) = hcib(ix,iy,is)
               endif
             else
               hcimx(ix,iy,1,is) = hcib(ix,iy,is)
             endif
           enddo
         enddo
       enddo

      else if (flux_limit_style.eq.1) then                    !srv 13.06.08 {

       do ix = -1, nx
         do iy = -1, ny
*   .. compute parallel neutral heat conductivity flux limit
           if (leftix(ix,iy).ne.-2) then
             do is = 0, ns-1
               if (is_neutral(is)) then
                 t0 = max(ti(ix,iy), ti(leftix(ix,iy),leftiy(ix,iy)))
                 t1 = (na(ix,iy,is) +
     .                 na(leftix(ix,iy),leftiy(ix,iy),is))
     .            /2.0_R8
                 dti = abs(ti(leftix(ix,iy),leftiy(ix,iy)) - ti(ix,iy))
     .            /(0.5_R8*(hx(ix,iy)+
     .                      hx(leftix(ix,iy),leftiy(ix,iy))))  ! dpc 2002-12-17
                 vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))
                 t2 = (1.0_R8 +
     .            ((hcib(ix,iy,is)*dti)/(alpha*vthi*t1*t0))**gamma)
     .            **(1.0_R8/gamma)
                 hcimx(ix,iy,0,is) = hcib(ix,iy,is)/t2
                 maxfluxlimit(0)=max(maxfluxlimit(0),t2)
               else
                 hcimx(ix,iy,0,is) = hcib(ix,iy,is)
               endif
             enddo
           else
             hcimx(ix,iy,0,:) = hcib(ix,iy,:)
           endif
*   .. compute perpendicular neutral heat conductivity flux limit
           if (bottomiy(ix,iy).ne.-2) then
             do is = 0, ns-1
               if (is_neutral(is)) then
                 t0 = max(ti(ix,iy),
     .                    ti(bottomix(ix,iy),bottomiy(ix,iy)))
                 t1 = (na(ix,iy,is) +
     .                 na(bottomix(ix,iy),bottomiy(ix,iy),is))/2.0_R8
                 dti = abs(ti(bottomix(ix,iy),bottomiy(ix,iy))-
     .                     ti(ix,iy))
     .            /(0.5_R8*(hy(ix,iy)+
     .                      hy(bottomix(ix,iy),bottomiy(ix,iy)))) ! dpc 2002-12-17
                 vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))
                 t2 = (1.0_R8 +
     .            ((hcib(ix,iy,is)*dti)/(alpha*vthi*t1*t0))**gamma)
     .            **(1.0_R8/gamma)
                 hcimx(ix,iy,1,is) = hcib(ix,iy,is)/t2
                 maxfluxlimit(1)=max(maxfluxlimit(1),t2)
                else
                 hcimx(ix,iy,1,is) = hcib(ix,iy,is)
                endif
              enddo
            else
              hcimx(ix,iy,1,:) = hcib(ix,iy,:)
           endif
         enddo
       enddo

      else
! new default treatment applying the flux limit to the transport coefficient directly
       do is = 0, ns-1
        if (is_neutral(is)) then
         do ix = -1, nx
          do iy = -1, ny
*   .. compute parallel neutral heat conductivity flux limit
           if (leftix(ix,iy).ne.-2) then
            t0 = (ti(ix,iy) + ti(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8
            t1 = (na(ix,iy,is) + na(leftix(ix,iy),leftiy(ix,iy),is))
     .            /2.0_R8
            dti = chcb(ix,iy,0,is)*
     &            abs(ti(leftix(ix,iy),leftiy(ix,iy)) - ti(ix,iy))
            vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))
     &            *(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))
     &            /(hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))
            t2 = 1.0_R8/(1.0_R8 + (dti/(alpha*vthi*t1*t0))**gamma)
     .                                           **(1.0_R8/gamma)
            chcb(ix,iy,0,is) = chcb(ix,iy,0,is)*t2
            fllim0fhi(ix,iy,0,is) = t2
            maxfluxlimit(0)=max(maxfluxlimit(0),t2)
           end if
*   .. compute perpendicular neutral heat conductivity flux limit
           if (bottomiy(ix,iy).ne.-2) then
            t0 = (ti(ix,iy)+ti(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8
            t1 = (na(ix,iy,is) + na(bottomix(ix,iy),bottomiy(ix,iy),is))
     &           /2.0_R8
            dti = chcb(ix,iy,1,is)*
     &            abs(ti(bottomix(ix,iy),bottomiy(ix,iy))-ti(ix,iy))
            vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))
     &           *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))
     &           /(hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
            t2 = 1.0_R8/(1.0_R8 +(dti/(alpha*vthi*t1*t0))**gamma)
     .                                          **(1.0_R8/gamma)
            chcb(ix,iy,1,is) = chcb(ix,iy,1,is)*t2
            fllim0fhi(ix,iy,1,is) = t2
            maxfluxlimit(1)=max(maxfluxlimit(1),t2)
           end if
          end do
         end do
        end if
       end do
      endif                                                       !srv 13.06.08 }

  100 continue

       if (ncall.lt.100) write(*,'(1x,a,1p,2e10.2)')
     .  'b2tlh0: maxfluxlimit = ',maxfluxlimit

      if (iout.ne.0) then
        allocate(w1(-1:nx,-1:ny))
        do is=0,ns-1                                            !srv 13.01.04 {
         write (charns,'(i3.3)') is
         if (flux_limit_style.eq.0 .or. flux_limit_style.eq.1) then
           call my_out(70, nx, ny, hcimx(-1,-1,0,is),'hcimx'//charns)
           call my_out(70, nx, ny, hcimx(-1,-1,1,is),'hcimy'//charns)
           w1 = 1.0e0_R8
           where (hcimx(:,:,0,is) /= 0.0e0_R8) w1 = hcib(:,:,is)
     &                                          /hcimx(:,:,0,is)
           call my_out(70, nx, ny, w1,'chcimx'//charns)
           w1 = 1.0e0_R8
           where (hcimx(:,:,1,is) /= 0.0e0_R8) w1 = hcib(:,:,is)
     &                                          /hcimx(:,:,1,is)
           call my_out(70, nx, ny, w1,'chcimy'//charns)
         else
           call my_out(70, nx, ny, fllim0fhi(-1,-1,0,is),
     &                 'fllim0fhix'//charns)
           call my_out(70, nx, ny, fllim0fhi(-1,-1,1,is),
     &                 'fllim0fhiy'//charns)
           call my_out(70, nx, ny, chcb(-1,-1,0,is), 'chcix'//charns)
           call my_out(70, nx, ny, chcb(-1,-1,1,is), 'chciy'//charns)
         end if
        enddo                                                   !srv 13.01.04 }
        deallocate(w1)
      endif

*   ..return
      ncall = ncall + 1
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2tlh0
