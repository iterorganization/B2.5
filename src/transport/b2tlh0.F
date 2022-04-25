      subroutine b2tlh0 (nx, ny, ns, hx, hy, gs, qc, qs, na, ti,
     & chcb, fllim0fhi)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2),
     *  qc(-1:nx,-1:ny,0:1), qs(-1:nx,-1:ny,0:1),
     *  na(-1:nx,-1:ny,0:ns-1), ti(-1:nx,-1:ny),
     *  chcb(-1:nx,-1:ny,0:1,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fllim0fhi(-1:nx,-1:ny,0:1,0:1,0:ns-1)
*   ..common blocks
*     ------------------------------------------------------------------
*     B2TLH0 computes the conductive neutral heat transport flux limit.
*     ------------------------------------------------------------------
*   ..local variables
      integer ix, iy, is, n2, ncall
      integer, save :: iout = 0                                   !srv 13.01.04 {
      real (kind=R8), allocatable :: w1(:,:)
      character charns*3                                          !srv 13.01.04 }
      real (kind=R8) ::
     *  alpha, gamma, t0, t1, t2, dti, vthi, maxfluxlimit(0:1),
     *  flux_limit_min_ti, dtir(-1:nx,-1:ny,0:1), wght(-1:nx,-1:ny),
     *  tic(-1:nx,-1:ny)
      save ncall, alpha, gamma, flux_limit_min_ti
*   ..initialisation
      data ncall/0/, alpha/0.0_R8/, gamma/2.0_R8/,
     *     flux_limit_min_ti/0.0_R8/
*   ..procedures
      intrinsic max, abs, sqrt
      external xertst, ipgetr, ipgeti
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
        call ipgeti ('b2tlh0_iout', iout)                         !srv 13.01.04
      endif

      call sfill(n2*4*ns, 1.0_R8, fllim0fhi, 1)
      maxfluxlimit=0.0_R8


      if (alpha.ne.0.0_R8) then
! new default treatment applying the flux limit to the transport coefficient directly
*    ..compute radial differences of temperature
       wght = 1.0_R8
       call diff_r (nx, ny, 0, hx, hy, wght, qs, ti, tic, dtir)
       do is = 0, ns-1
        if (is_neutral(is)) then
         do ix = -1, nx
          do iy = -1, ny
*   .. compute neutral heat conductivity flux limit on x-surfaces
           if (leftix(ix,iy).ne.-2) then
*   .. parallel
            t0 = (ti(ix,iy) + ti(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8
            t1 = (na(ix,iy,is) + na(leftix(ix,iy),leftiy(ix,iy),is))
     .            /2.0_R8
            dti = chcb(ix,iy,0,0,is)*
     &            abs(ti(leftix(ix,iy),leftiy(ix,iy)) - ti(ix,iy))
            vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))
     &            *gs(ix,iy,0)*qc(ix,iy,0)
            t2 = 1.0_R8/(1.0_R8 + (dti/(alpha*vthi*t1*t0))**gamma)
     .                                           **(1.0_R8/gamma)
            chcb(ix,iy,0,0,is) = chcb(ix,iy,0,0,is)*t2
            fllim0fhi(ix,iy,0,0,is) = t2
            maxfluxlimit(0)=max(maxfluxlimit(0),t2)
*   .. perpendicular
            if (chcb(ix,iy,0,1,is).ne.0.0_R8) then
            dti  = abs(chcb(ix,iy,0,1,is)*dtir(ix,iy,0))
            vthi = abs(sqrt(2.0_R8 * t0 / (am(is)*mp))*
     &        gs(ix,iy,0)*qs(ix,iy,0))
            t2 = 1.0_R8/(1.0_R8 + (dti/(alpha*vthi*t1*t0))**gamma)
     .                                           **(1.0_R8/gamma)
            chcb(ix,iy,0,1,is) = chcb(ix,iy,0,1,is)*t2
            fllim0fhi(ix,iy,0,1,is) = t2
            endif
           end if
*   .. compute perpendicular neutral heat conductivity flux limit
           if (bottomiy(ix,iy).ne.-2) then
            t0 = (ti(ix,iy)+ti(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8
            t1 = (na(ix,iy,is) + na(bottomix(ix,iy),bottomiy(ix,iy),is))
     &           /2.0_R8
            dti = chcb(ix,iy,1,1,is)*abs(dtir(ix,iy,1))
            vthi = sqrt(2.0_R8 * t0 / (am(is)*mp))
     &           *gs(ix,iy,1)
            t2 = 1.0_R8/(1.0_R8 +(dti/(alpha*vthi*t1*t0))**gamma)
     .                                          **(1.0_R8/gamma)
            chcb(ix,iy,1,1,is) = chcb(ix,iy,1,1,is)*t2
            fllim0fhi(ix,iy,1,1,is) = t2
            maxfluxlimit(1)=max(maxfluxlimit(1),t2)
           end if
          end do
         end do
        end if
       end do
      endif                                                       !srv 13.06.08 }

      if (ncall.lt.100) write(*,'(1x,a,1p,2e10.2)')
     .  'b2tlh0: maxfluxlimit = ',maxfluxlimit

      if (iout.ne.0) then
        allocate(w1(-1:nx,-1:ny))
        do is=0,ns-1                                            !srv 13.01.04 {
         write (charns,'(i3.3)') is
         call my_out(70, nx, ny, fllim0fhi(-1,-1,0,0,is),
     &                 'fllim0fhix'//charns)
         call my_out(70, nx, ny, fllim0fhi(-1,-1,1,1,is),
     &                 'fllim0fhiy'//charns)
         call my_out(70, nx, ny, chcb(-1,-1,0,0,is), 'chcix'//charns)
         call my_out(70, nx, ny, chcb(-1,-1,1,1,is), 'chciy'//charns)
        enddo                                                   !srv 13.01.04 }
        deallocate(w1)
      endif

*   ..return
      ncall = ncall + 1
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2tlh0
