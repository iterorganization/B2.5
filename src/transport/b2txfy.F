      subroutine b2txfy (nx, ny, ns, mode, vol, gs, hy, qz, kapy, cony)
      use b2mod_types
      use b2mod_math
      use b2mod_indirect
      use b2mod_geo
     , , only : hy1
      use b2mod_subsys
      implicit none
      integer nx, ny, ns, mode
      real (kind=R8) ::
     *   vol(-1:nx,-1:ny), hy(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2), 
     &   qz(-1:nx,-1:ny,0:1),
     &   kapy(-1:nx,-1:ny,0:1,0:ns-1), cony(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2TXFY computes conductances between neighbouring cells in the
*     y-direction, given the flux-limited conductivities.
*     ------------------------------------------------------------------
      integer ix, iy, is
      real (kind=R8) ::
     *   t0, t1, kaps, kapsn
      intrinsic max, sqrt
      external xertst
      external b2xxgs
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2txfy')
*   ..test mode
      call xertst (-1.le.mode.and.mode.le.2,
     &  'b2txfy--faulty argument mode')

*   ..compute cony
      do iy = -1, ny
       do ix = -1, nx
        if(bottomiy(ix,iy).eq.-2) then
         cony(ix,iy) = 0.0_R8
        else
         kaps = 0.0_R8
         kapsn= 0.0_R8
         do is = 0, ns-1
           kaps = kaps + kapy(ix,iy,1,is)
           kapsn= kapsn+ kapy(bottomix(ix,iy),bottomiy(ix,iy),1,is)
         enddo
         if (mode.eq.-1) then
*    ..use harmonic averaging
           t0 = vol(bottomix(ix,iy),bottomiy(ix,iy))*kapsn/
     &        hy1(bottomix(ix,iy),bottomiy(ix,iy))**2
           t1 = vol(ix,iy)*kaps/hy1(ix,iy)**2
           if (t0*t1.le.0.0_R8) then
            cony(ix,iy) = 0.0_R8
           else
            cony(ix,iy) = harm(t0,t1)
           endif
         else if (mode.eq.0) then
*    ..use geometric averaging
           cony(ix,iy) = 2.0_R8*sqrt(max(0.0_R8,
     &      (vol(bottomix(ix,iy),bottomiy(ix,iy))*kapsn/
     &       hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &      (vol(ix,iy)*kaps/hy1(ix,iy))))/
     &      (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
         else if (mode.eq.1) then
*    ..use arithmetic averaging (SOLPS5.0 formulation)
           cony(ix,iy) = 2.0_R8*
     &       (vol(bottomix(ix,iy),bottomiy(ix,iy))*kapsn
     &       +vol(ix,iy)*kaps)/
     &     (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))**2
         else if (mode.eq.2) then
*    ..use arithmetic averaging (SOLPS4.0 formulation)
           cony(ix,iy) = gs(ix,iy,1)*2.0_R8*
     &       (hy1(bottomix(ix,iy),bottomiy(ix,iy))*kapsn
     &       +hy1(ix,iy)*kaps)/
     &     (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))**2
         endif
        endif
       enddo
      enddo
*     (Note: We used to have the loop over ix outside the one over iy
*     in the previous code, but then the AIX/XLF compiler version 2.2
*     produced wrong code when optimisation (-O) was on.)
*   ..reset null regions
      call b2xxgs (nx, ny, 0.0_R8, cony, 2)
*   ..return
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2txfy
