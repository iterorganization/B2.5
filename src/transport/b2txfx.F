      subroutine b2txfx (nx, ny, ns, mode, vol, gs, hx, qz, kapx, conx)
      use b2mod_types
      use b2mod_indirect
      implicit none
      integer nx, ny, ns, mode
      real (kind=R8) ::
     *   vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2), 
     &   qz(-1:nx,-1:ny,0:1),
     &   kapx(-1:nx,-1:ny,0:1,0:ns-1), conx(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2TXFX computes conductances between neighbouring cells in the
*     x-direction, given the flux-limited conductivities.
*     ------------------------------------------------------------------
      integer n2, ix, iy, is
      real (kind=R8) ::
     *   t0, t1, kaps, kapsn
      intrinsic max, sqrt
      real (kind=R8) ::
     *   harm, smin, smax
      external subini, subend, xertst, smin, smax
      external b2xvsg, b2xxgs
      harm(t0,t1) = 2.0_R8*t0*t1/(t0+t1)
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2txfx')
*   ..test nx, ny, ns
      call xertst (0.le.nx .and. 0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..set n2
      n2=(nx+2)*(ny+2)
*   ..test mode
      call xertst (-1.le.mode.and.mode.le.2,
     &  'b2txfx--faulty argument mode')
*   ..test sign of vol, hx
      call b2xvsg (n2, vol, 1, 'vol', '.gt.')
      call b2xvsg (n2, hx, 1, 'hx', '.gt.')
*   ..test range of qz
      call xertst (-1.0_R8.lt.smin(n2,qz(-1,-1,0),1).and.
     &  smax(n2,qz(-1,-1,0),1).lt.1.0_R8, 
     &  'faulty argument range qz(,,0)')
      call xertst (0.0_R8.lt.smin(n2,qz(-1,-1,1),1).and.
     &  smax(n2,qz(-1,-1,1),1).le.1.0_R8, 
     &  'faulty argument range qz(,,1)')

*   ..compute conx
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
          conx(ix,iy) = 0.0_R8
        else
         kaps = 0.0_R8
         kapsn= 0.0_R8
         do is = 0, ns-1
           kaps = kaps + kapx(ix,iy,0,is)
           kapsn= kapsn+ kapx(leftix(ix,iy),leftiy(ix,iy),0,is)
         enddo
         if (mode.eq.-1) then
*    ..use harmonic averaging
*      (needs somewhat messy code because the harmonic average is only
*      defined for two quantities of the same sign.)
           t0 = vol(leftix(ix,iy),leftiy(ix,iy))*kapsn/
     &        hx(leftix(ix,iy),leftiy(ix,iy))**2
           t1 = vol(ix,iy)*kaps/hx(ix,iy)**2
           if (t0*t1.le.0.0_R8) then
             conx(ix,iy) = 0.0_R8
           else
             conx(ix,iy) = harm(t0,t1)
           endif
         else if (mode.eq.0) then
*    ..use geometric averaging
           conx(ix,iy) = 2.0_R8*sqrt(max(0.0_R8,
     &      (vol(leftix(ix,iy),leftiy(ix,iy))*kapsn/
     &       hx(leftix(ix,iy),leftiy(ix,iy)))*
     &      (vol(ix,iy)*kaps/hx(ix,iy))))/
     &      (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
         else if (mode.eq.1) then
*    ..use arithmetic averaging (SOLPS5.0 formulation)
           conx(ix,iy) = 2.0_R8*
     &        (vol(leftix(ix,iy),leftiy(ix,iy))*kapsn
     &        +vol(ix,iy)*kaps)/
     &        (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))**2
         else if (mode.eq.2) then
*    ..use arithmetic averaging (SOLPS4.0 formulation)
           conx(ix,iy) = gs(ix,iy,0)*qz(ix,iy,1)*2.0_R8*
     &        (hx(leftix(ix,iy),leftiy(ix,iy))*kapsn
     &        +hx(ix,iy)*kaps)/
     &        (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))**2
         endif
        endif
       enddo
      enddo
*   ..reset null regions
      call b2xxgs (nx, ny, 0.0_R8, conx, 1)
*   ..return
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2txfx
