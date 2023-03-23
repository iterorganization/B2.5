      subroutine b2txcv (nx, ny, mode, vsb, cvsb)
      use b2mod_types
      use b2mod_geo
      use b2mod_indirect
      use b2mod_subsys
      implicit none
      integer nx, ny, mode
      real (kind=R8) ::
     &  vsb(-1:nx,-1:ny,0:1),
     &  cvsb(-1:nx,-1:ny,0:1)
*     ------------------------------------------------------------------
*     B2TXCV is obsolete.
*     ------------------------------------------------------------------
*     ------------------------------------------------------------------
*     B2TXCV computes viscous conductances between neighbouring
*     staggered cells, given the viscosities.
*     The viscosities vsb(,,0) and vsb(,,1) are defined on cells.
*     Because the parallel velocity is defined on cell faces staggered
*     in the x-direction, the viscous conductance cvsb(,,0) must be
*     defined on cells, and the viscous conductance cvsb(,,1) must be
*     defined on cell corners.
*     ------------------------------------------------------------------
      integer ix, iy
      real (kind=R8) ::
     *   t0, t1
      intrinsic max, sqrt
      real (kind=R8) ::
     *   harm
      external xertst
      external b2xxgs, harm
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2txcv')
*   ..test mode
      call xertst (-1.le.mode.and.mode.le.2,
     &  'b2txcv--faulty argument mode')
*   ..compute cvsb(,,0)
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         cvsb(ix,iy,0) = 0.0_R8
        elseif(rightix(ix,iy).eq.nx+1) then
         cvsb(ix,iy,0) = 0.0_R8
        else
         cvsb(ix,iy,0) = vol(ix,iy)*vsb(ix,iy,0)/hx(ix,iy)**2
        endif
       enddo
      enddo
*   ..reset null regions
      call b2xxgs (nx, ny, 0.0_R8, cvsb(-1,-1,0), 0)
*   ..compute cvsb(,,1)
*     (temporary assignment to cvsb(,,1))
      if (mode.eq.-1) then
*    ..use harmonic averaging
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).eq.-2) then
          cvsb(ix,iy,1) = 0.0_R8
         else 
          t0 = vol(bottomix(ix,iy),bottomiy(ix,iy))*
     &         vsb(bottomix(ix,iy),bottomiy(ix,iy),1)/
     &         hy1(bottomix(ix,iy),bottomiy(ix,iy))**2
          t1 = vol(ix,iy)*vsb(ix,iy,1)/hy1(ix,iy)**2
          if (t0*t1.le.0.0_R8) then
           cvsb(ix,iy,1) = 0.0_R8
          else
           cvsb(ix,iy,1) = harm(t0,t1)
          endif
         endif
        enddo
       enddo
      else if (mode.eq.0) then
*    ..use geometric averaging
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).eq.-2) then
          cvsb(ix,iy,1) = 0.0_R8
         else
          cvsb(ix,iy,1) = 2.0_R8*sqrt(max(0.0_R8,
     &     (vol(bottomix(ix,iy),bottomiy(ix,iy))*
     &      vsb(bottomix(ix,iy),bottomiy(ix,iy),1)/
     &      hy1(bottomix(ix,iy),bottomiy(ix,iy))*
     &     (vol(ix,iy)*vsb(ix,iy,1)/hy1(ix,iy)))))/
     &     (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
         endif
        enddo
       enddo
      else if (mode.eq.1) then
*    ..use arithmetic averaging (SOLPS5.0 formulation)
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).eq.-2) then
          cvsb(ix,iy,1) = 0.0_R8
         else
          cvsb(ix,iy,1) =
     &     2.0_R8*(vol(bottomix(ix,iy),bottomiy(ix,iy))*
     &        vsb(bottomix(ix,iy),bottomiy(ix,iy),1)
     &       +vol(ix,iy)*vsb(ix,iy,1))/
     &     (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))**2
         endif
        enddo
       enddo
      else if (mode.eq.2) then
*    ..use arithmetic averaging (SOLPS4.0 formulation)
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).eq.-2) then
          cvsb(ix,iy,1) = 0.0_R8
         else
          cvsb(ix,iy,1) =
     &     2.0_R8*(hy1(bottomix(ix,iy),bottomiy(ix,iy))*
     &             vsb(bottomix(ix,iy),bottomiy(ix,iy),1)
     &            +hy1(ix,iy)*vsb(ix,iy,1))/
     &     (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))**2
         endif
        enddo
       enddo
      endif
*     (transfer cvsb(,,1) to cell corners)
      do ix = nx, -1, -1
       do iy = -1, ny
        if(leftix(ix,iy).ne.-2) then
         cvsb(ix,iy,1) = (cvsb(ix,iy,1)
     &                   +cvsb(leftix(ix,iy),leftiy(ix,iy),1))/2
        else
         cvsb(ix,iy,1) = 0.0_R8
        endif
       enddo
      enddo
*   ..reset null regions
      call b2xxgs (nx, ny, 0.0_R8, cvsb(-1,-1,1), 3)
*   ..return
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2txcv
