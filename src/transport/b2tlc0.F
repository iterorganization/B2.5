      subroutine b2tlc0 (nx, ny, ns, hx, hy, na, te, ti, dpa0, dpa0mx)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
	use b2mod_b2cmpa
	use b2mod_rates
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     *  na(-1:nx,-1:ny,0:ns-1), te(-1:nx,-1:ny), ti(-1:nx,-1:ny), 
     1  dpa0(-1:nx,-1:ny,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  dpa0mx(-1:nx,-1:ny,0:1,0:ns-1)
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
      data ncall/0/, alpha/0.0_R8/, gamma/2.0_R8/
*   ..procedures
      intrinsic max, abs, sqrt
      external subini, subend, xertst, ipgetr
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

      call sfill(n2*2*ns, 0.0_R8, dpa0mx, 1)
*   ..nothing done if alpha.eq.0
      if (alpha.eq.0.0_R8) then
        dpa0mx(:,:,0,:) = dpa0(:,:,:)
        dpa0mx(:,:,1,:) = dpa0(:,:,:)
        goto 100
      endif 

      maxfluxlimit=0

      do is = 0, ns-1
       if (is_neutral(is)) then
*    ..compute partial pressure
        call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)
        do ix = -1, nx
         do iy = -1, ny
*   ..compute parallel neutral difusion flux limit
          if (leftix(ix,iy).ne.-2) then
           t0 = max(ti(ix,iy), ti(leftix(ix,iy),leftiy(ix,iy)))
           t1 = (na(ix,iy,is) + na(leftix(ix,iy),leftiy(ix,iy),is))
     .          /2.0_R8
           dpb = abs(pb(leftix(ix,iy),leftiy(ix,iy)) - pb(ix,iy))
     .           /(0.5_R8*(hx(ix,iy)+
     .                     hx(leftix(ix,iy),leftiy(ix,iy))))
           vbar = sqrt(8.0_R8 * t0 / (pi * (am(is)*mp)))
           t2 = (1.0_R8 + ((dpa0(ix,iy,is)*dpb)/
     .          (alpha*vbar/4.0e0_R8*t1))**gamma)**(1.0_R8/gamma)
           dpa0mx(ix,iy,0,is) = dpa0(ix,iy,is)/t2
           maxfluxlimit(0)=max(maxfluxlimit(0),t2)
          endif
	   enddo
	  enddo
*   ..compute perpendicular neutral difusion flux limit
        do ix = -1, nx
         do iy = -1, ny
          if (bottomiy(ix,iy).ne.-2) then
            t0 = max(ti(ix,iy), ti(bottomix(ix,iy),bottomiy(ix,iy)))
            t1 = (na(ix,iy,is) + 
     .            na(bottomix(ix,iy),bottomiy(ix,iy),is))/2.0_R8
            dpb = abs(pb(bottomix(ix,iy),bottomiy(ix,iy))-pb(ix,iy))
     .           /(0.5_R8*(hy(ix,iy)+
     .                     hy(bottomix(ix,iy),bottomiy(ix,iy))))
            vbar = sqrt(8.0_R8 * t0 / (pi * (am(is)*mp)))
            t2 = (1.0_R8 + ((dpa0(ix,iy,is)*dpb)/
     .            (alpha*vbar/4.0e0_R8*t1))**gamma)**(1.0_R8/gamma)
            dpa0mx(ix,iy,1,is) = dpa0(ix,iy,is)/t2
            maxfluxlimit(1)=max(maxfluxlimit(1),t2)
	    endif
	   enddo
	  enddo
       else
        dpa0mx(:,:,0,is) = dpa0(:,:,is)
        dpa0mx(:,:,1,is) = dpa0(:,:,is)
       endif
      enddo

!      write(*,"(1x,a,1p,2e10.2)") 'b2tlc0: maxfluxlimit = ',maxfluxlimit

 100  continue

*   ..return
      ncall = ncall + 1
      call subend ()
      return
*     ------------------------------------------------------------------
      end
