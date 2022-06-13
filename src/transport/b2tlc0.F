      subroutine b2tlc0 (nx, ny, ns, hx, hy, gs, qc, qs,
     &  na, te, ti, cdpa, flc)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_rates
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2),
     *  qc(-1:nx,-1:ny,0:1), qs(-1:nx,-1:ny,0:1),
     *  na(-1:nx,-1:ny,0:ns-1), te(-1:nx,-1:ny), ti(-1:nx,-1:ny)
*   ..input and output arguments
      real (kind=R8) ::
     *  cdpa(-1:nx,-1:ny,0:1,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  flc(-1:nx,-1:ny,0:1,0:1,0:ns-1)
*   ..common blocks
*     ------------------------------------------------------------------
*     B2TLC0 computes the neutral diffusion flux limit.              !srv 09.07.03
*     ------------------------------------------------------------------
*   ..local variables
      integer ix, iy, is, n2, ncall
      real (kind=R8) ::
     *  alpha, gamma, t0, t1, t2, dpb, vbar, maxfluxlimit(0:1),
     *  pb(-1:nx,-1:ny), dpbr(-1:nx,-1:ny,0:1), wght(-1:nx,-1:ny),
     *  pbc(-1:nx,-1:ny)
      save ncall, alpha, gamma
*   ..initialization
      data ncall/0/, alpha/1.0_R8/, gamma/1.0_R8/
*   ..procedures
      intrinsic max, abs, sqrt
      external xertst, ipgetr
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

      call sfill(n2*4*ns, 1.0_R8, flc, 1)
      call sfill(2, 1.0_R8, maxfluxlimit, 1)
*   ..nothing done if alpha.eq.0
      if (alpha.ne.0.0_R8) then
       wght=1.0_R8
       do is = 0, ns-1
        if (is_neutral(is)) then
*    ..compute partial pressure
         call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)
*    ..compute radial differences of partial pressure
         call diff_r (nx, ny, 0, hx, hy, wght, qs, pb, pbc, dpbr)
         do ix = -1, nx
          do iy = -1, ny
*   ..compute parallel neutral diffusion flux limit
           if (leftix(ix,iy).ne.-2) then
            t0 = (ti(ix,iy) + ti(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8
            t1 = (na(ix,iy,is) + na(leftix(ix,iy),leftiy(ix,iy),is))
     &          /2.0_R8
*    ..parallel
            dpb = cdpa(ix,iy,0,0,is)*
     &            abs(pb(leftix(ix,iy),leftiy(ix,iy)) - pb(ix,iy))
            vbar = sqrt(8.0_R8 * t0 / (pi * (am(is)*mp)))
     &            *gs(ix,iy,0)*qc(ix,iy,0)
            t2 = 1.0_R8/(1.0_R8 +
     &          (dpb/(alpha*vbar/4.0e0_R8*t1))**gamma)**(1.0_R8/gamma)
            cdpa(ix,iy,0,0,is) = cdpa(ix,iy,0,0,is)*t2
            flc(ix,iy,0,0,is) = t2
            maxfluxlimit(0)=max(maxfluxlimit(0),t2)
*    ..perpendicular
            if (cdpa(ix,iy,0,1,is).ne.0.0_R8) then
             dpb  = abs(cdpa(ix,iy,0,1,is)*dpbr(ix,iy,0))
             vbar = abs(sqrt(8.0_R8 * t0 / (pi * (am(is)*mp)))*
     &        *gs(ix,iy,0)*qs(ix,iy,0))
             t2 = 1.0_R8/(1.0_R8 +
     &           (dpb/(alpha*vbar/4.0e0_R8*t1))**gamma)**(1.0_R8/gamma)
             cdpa(ix,iy,0,1,is) = cdpa(ix,iy,0,1,is)*t2
             flc(ix,iy,0,1,is) = t2
            endif
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
            dpb = cdpa(ix,iy,1,1,is)*abs(dpbr(ix,iy,1))
            vbar = sqrt(8.0_R8 * t0 / (pi * (am(is)*mp)))
     &           *gs(ix,iy,1)
            t2 = 1.0_R8/(1.0_R8 + (dpb/(alpha*vbar/4.0e0_R8*t1))**gamma)
     &                                                  **(1.0_R8/gamma)
            cdpa(ix,iy,1,1,is) = cdpa(ix,iy,1,1,is)*t2
            flc(ix,iy,1,1,is) = t2
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
