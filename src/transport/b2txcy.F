      subroutine b2txcy (nx, ny, mode, vol, gs, hy, qz, kapy, cony)
      use b2mod_types
      use b2mod_indirect
      implicit none
      integer nx, ny, mode
      real (kind=R8) ::
     *   vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2), hy(-1:nx,-1:ny), 
     &   qz(-1:nx,-1:ny,0:1), kapy(-1:nx,-1:ny), 
     &   cony(-1:nx,-1:ny,0:1)
*     ------------------------------------------------------------------
*     B2TXCY computes radial conductances between neighbouring, 
*     given the conductivities.
*     cony(:,:,0) contains radial conductivities across x-faces
*     cony(:,:,1) contains radial conductivities across y-faces
*     ------------------------------------------------------------------
      integer ix, iy
      real (kind=R8) ::
     *   t0, t1
      intrinsic max, sqrt
      real (kind=R8) ::
     *   hy1, harm
      external subini, subend, xertst
      external b2xxgs
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
      harm(t0,t1) = 2.0_R8*t0*t1/(t0+t1)
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2txcy')
*   ..test mode
      call xertst (-1.le.mode.and.mode.le.2,
     &  'b2txcy--faulty argument mode')

*   ..initialize cony
      cony = 0.0_R8
*   ..compute cony
      if (mode.eq.-1) then
*    ..use harmonic averaging
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).ne.-2) then
          t0 = vol(bottomix(ix,iy),bottomiy(ix,iy))*
     &         kapy(bottomix(ix,iy),bottomiy(ix,iy))/
     &         hy1(bottomix(ix,iy),bottomiy(ix,iy))**2
          t1 = vol(ix,iy)*kapy(ix,iy)/hy1(ix,iy)**2
          if (t0*t1.gt.0.0_R8) then
           cony(ix,iy,1) = harm(t0,t1)
          endif
         endif
        enddo
       enddo
      else if (mode.eq.0) then
*    ..use geometric averaging
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).ne.-2) then
          cony(ix,iy,1) = 2.0_R8*sqrt(max(0.0_R8,
     &     (vol(bottomix(ix,iy),bottomiy(ix,iy))*
     &      kapy(bottomix(ix,iy),bottomiy(ix,iy))/
     &      hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &     (vol(ix,iy)*kapy(ix,iy)/hy1(ix,iy))))/
     &     (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
         endif
        enddo
       enddo
      else if (mode.eq.1) then
*    ..use arithmetic averaging (SOLPS5.0 formulation)
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).ne.-2) then
          cony(ix,iy,1) = 2.0_R8*
     &     (vol(bottomix(ix,iy),bottomiy(ix,iy))*
     &        kapy(bottomix(ix,iy),bottomiy(ix,iy))
     &       +vol(ix,iy)*kapy(ix,iy))/
     &     (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))**2
         endif
        enddo
       enddo
      else if (mode.eq.2) then
*    ..use arithmetic averaging (SOLPS4.0 formulation)
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).ne.-2) then
          cony(ix,iy,1) = gs(ix,iy,1)*2.0_R8*
     &       (hy1(bottomix(ix,iy),bottomiy(ix,iy))*
     &        kapy(bottomix(ix,iy),bottomiy(ix,iy))
     &       +hy1(ix,iy)*kapy(ix,iy))/
     &     (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))**2
         endif
        enddo
       enddo
      endif
*     (Note: We used to have the loop over ix outside the one over iy
*     in the previous code, but then the AIX/XLF compiler version 2.2
*     produced wrong code when optimisation (-O) was on.)
*   ..reset null regions
      call b2xxgs (nx, ny, 0.0_R8, cony(-1,-1,0), 1)
      call b2xxgs (nx, ny, 0.0_R8, cony(-1,-1,1), 2)
*   ..return
      call subend ()
      return
*     ------------------------------------------------------------------
      end
