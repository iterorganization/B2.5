      subroutine b2tlc0 (nx, ny, ns, hx, hy, vol, qz,
     &  na, te, ti, cdpa, flc)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_rates
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  hx(-1:nx,-1:ny), vol(-1:nx,-1:ny),
     *  hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1),
     *  na(-1:nx,-1:ny,0:ns-1), te(-1:nx,-1:ny), ti(-1:nx,-1:ny)
*   ..input and output arguments
      real (kind=R8) ::
     *  cdpa(-1:nx,-1:ny,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  flc(-1:nx,-1:ny,0:1,0:ns-1)
*   ..common blocks
*     ------------------------------------------------------------------
*     B2TLC0 computes the neutral diffusion flux limit.              !srv 09.07.03
*     ------------------------------------------------------------------
*   ..local variables
      integer ix, iy, is, n2, ncall
      real (kind=R8) ::
     *  alpha, gamma, t0, t1, t2, dpb, vbar, maxfluxlimit(0:1),
     *  pb(-1:nx,-1:ny)
      save ncall, alpha, gamma
*   ..initialization
      data ncall/0/, alpha/1.0_R8/, gamma/1.0_R8/
*   ..procedures
      intrinsic max, abs, sqrt
      real (kind=R8) ::
     *   hy1
      external subini, subend, xertst, ipgetr
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2tlc0')
*   ..test nx, ny, ns
      call xertst (0.le.nx .and. 0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..read input on first call
      if (ncall.eq.0) then
        call ipgetr ('b2tlc0_alpha', alpha)
        call xertst (0.le.alpha, 'invalid input parameter alpha')
        call ipgetr ('b2tlc0_gamma', gamma)
        call xertst (0.lt.gamma, 'invalid input parameter gamma')
      endif

      call sfill(n2*2*ns, 1.0_R8, flc, 1)
      call sfill(2, 1.0_R8, maxfluxlimit, 1)
*   ..nothing done if alpha.eq.0
      if (alpha.ne.0.0_R8) then
       do is = 0, ns-1
        if (is_neutral(is)) then
*    ..compute partial pressure
         call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)
         do ix = -1, nx
          do iy = -1, ny
*   ..compute parallel neutral diffusion flux limit
           if (leftix(ix,iy).ne.-2) then
            t0 = (ti(ix,iy) + ti(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8
            t1 = (na(ix,iy,is) + na(leftix(ix,iy),leftiy(ix,iy),is))
     &          /2.0_R8
            dpb = cdpa(ix,iy,0,is)*
     &            abs(pb(leftix(ix,iy),leftiy(ix,iy)) - pb(ix,iy))
            vbar = sqrt(8.0_R8 * t0 / (pi * (am(is)*mp)))
     &            *(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))
     &            /(hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))
            t2 = 1.0_R8/(1.0_R8 + (dpb/(alpha*vbar/4.0e0_R8*t1))**gamma)
     &                                                 **(1.0_R8/gamma)
            cdpa(ix,iy,0,is) = cdpa(ix,iy,0,is)*t2
            flc(ix,iy,0,is) = t2
            maxfluxlimit(0)=max(maxfluxlimit(0),t2)
           endif
          enddo
         enddo
*   ..compute perpendicular neutral diffusion flux limit
         do ix = -1, nx
          do iy = -1, ny
           if (bottomiy(ix,iy).ne.-2) then
            t0 = (ti(ix,iy)+ti(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8
            t1 = (na(ix,iy,is) + 
     &            na(bottomix(ix,iy),bottomiy(ix,iy),is))/2.0_R8
            dpb = cdpa(ix,iy,1,is)*
     &            abs(pb(bottomix(ix,iy),bottomiy(ix,iy))-pb(ix,iy))
            vbar = sqrt(8.0_R8 * t0 / (pi * (am(is)*mp)))
     &           *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))
     &           /(hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
            t2 = 1.0_R8/(1.0_R8 + (dpb/(alpha*vbar/4.0e0_R8*t1))**gamma)
     &                                                  **(1.0_R8/gamma)
            cdpa(ix,iy,1,is) = cdpa(ix,iy,1,is)*t2
            flc(ix,iy,1,is) = t2
            maxfluxlimit(1)=max(maxfluxlimit(1),t2)
           endif
          enddo
         enddo
        endif
       enddo
      endif

!      write(*,'(1x,a,1p,2e10.2)') 'b2tlc0: maxfluxlimit = ',maxfluxlimit

*   ..return
      ncall = ncall + 1
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2tlc0
