      subroutine b2xgbs (nx, ny, wbbl, wbbv, gs, qc, crx, cry, pbs)
      use b2mod_types
      use b2mod_math
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *  wbbl(-1:nx,-1:ny,0:3), wbbv(-1:nx,-1:ny,0:3),
     &  gs(-1:nx,-1:ny,0:1), crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     &  pbs(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XGBS computes the product (magnetic field pitch)*(surface area)
*     at cell faces. pbs(,,0) will contain the left face value and
*     pbs(,,1) will contain the bottom face value.
*     pbs thus contains the parallel contact area.
*     ------------------------------------------------------------------
      integer n2, ix, iy, ncall, area_fix
      real (kind=R8) :: t1
      external subini, subend, xertst, b2xvsg, b2xvff, ipgeti
      save ncall, area_fix
      data ncall/0/, area_fix/3/
*-----------------------------------------------------------------------
*   ..preliminaries
      call subini ('b2xgbs')
*   ..read switches
      if (ncall.eq.0) then
        call ipgeti ('b2news_area_fix', area_fix)
      endif
*   ..test input
      call xertst (0.le.nx.and.0.le.ny, 'faulty input nx, ny')
*   ..test sign of input arrays
      n2 = (nx+2)*(ny+2)
      call b2xvsg (n2, wbbl(-1,-1,3), 1, 'wbbl', '.gt.')
      call b2xvsg (n2, wbbv(-1,-1,3), 1, 'wbbv', '.gt.')
*   ..test edge values of gs
      call b2xvff (nx, ny, gs, 'gs')
*   ..compute pbs(,,0:1)
      do iy = -1, ny
       do ix = -1, nx
        t1 = cosang (
     &     crx(ix,iy,1)-crx(ix,iy,0),  cry(ix,iy,1)-cry(ix,iy,0),
     &    (crx(ix,iy,2)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,1)),
     &    (cry(ix,iy,2)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,1)))
        pbs(ix,iy,0) = gs(ix,iy,0)*wbbl(ix,iy,0)/wbbl(ix,iy,3)        !srv 11.01.13
        if (area_fix.ge.3) pbs(ix,iy,0) = pbs(ix,iy,0)*qc(ix,iy)
        pbs(ix,iy,1) = (gs(ix,iy,1)/t1)*wbbv(ix,iy,1)/wbbv(ix,iy,3)   !srv 11.01.13
       enddo
      enddo
*   ..test sign of pbs
!      call b2xvsg (n2*2, pbs, 1, 'pbs', '.ge.')                      !srv 11.01.13
*   ..test edge values of pbs
      call b2xvff (nx, ny, pbs, 'pbs')
*   ..return
      ncall = ncall + 1
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xgbs

      subroutine b2xgbs1 (nx1, ny1, 
     . bottomiy1, wbbl1, wbbv1, gs1, qc1, crx1, cry1, pbs1)
      use b2mod_types
      use b2mod_math
      implicit none
      integer nx1, ny1
      integer bottomiy1(-1:nx1,-1:ny1)
      real (kind=R8) ::
     *  wbbl1(-1:nx1,-1:ny1,0:3), wbbv1(-1:nx1,-1:ny1,0:3),
     &  gs1(-1:nx1,-1:ny1,0:2), crx1(-1:nx1,-1:ny1,0:3),
     &  cry1(-1:nx1,-1:ny1,0:3), qc1(-1:nx1,-1:ny1),
     &  pbs1(-1:nx1,-1:ny1,0:1)
*     ------------------------------------------------------------------
*     B2XGBS1 computes the product (magnetic field pitch)*(surface area)
*     at cell faces on the basis mesh. 
*     pbs1(,,0) will contain the left face value and
*     pbs1(,,1) will contain the bottom face value.
*     pbs1 thus contains the parallel contact area.
*     ------------------------------------------------------------------
      integer n2, ix, iy, ncall, area_fix
      real (kind=R8) :: t1
      external subini, subend, xertst, b2xvsg, b2xvfx, b2xvfy1, ipgeti
      save ncall, area_fix
      data ncall/0/, area_fix/3/
*-----------------------------------------------------------------------
*   ..preliminaries
      call subini ('b2xgbs1')
*   ..read switches
      if (ncall.eq.0) then
        call ipgeti ('b2news_area_fix', area_fix)
      endif
*   ..test input
      call xertst (0.le.nx1.and.0.le.ny1, 'faulty input nx1, ny1')
*   ..test sign of input arrays
      n2 = (nx1+2)*(ny1+2)
      call b2xvsg (n2, wbbl1(-1,-1,3), 1, 'wbbl1', '.gt.')
      call b2xvsg (n2, wbbv1(-1,-1,3), 1, 'wbbv1', '.gt.')
*   ..test edge values of gs1
      call b2xvfx (nx1, ny1, gs1(-1,-1,0), 'gs1x')
      call b2xvfy1(nx1, ny1, gs1(-1,-1,1), 'gs1y', bottomiy1)
      call b2xvsg (n2, gs1(-1,-1,2), 1, 'gsz', '.gt.')
*   ..compute pbs1(,,0:1)
      do iy = -1, ny1
       do ix = -1, nx1
        t1 = cosang (
     &     crx1(ix,iy,1)-crx1(ix,iy,0),  cry1(ix,iy,1)-cry1(ix,iy,0),
     &    (crx1(ix,iy,2)+crx1(ix,iy,3))-(crx1(ix,iy,0)+crx1(ix,iy,1)),
     &    (cry1(ix,iy,2)+cry1(ix,iy,3))-(cry1(ix,iy,0)+cry1(ix,iy,1)))
        pbs1(ix,iy,0) = gs1(ix,iy,0)*wbbl1(ix,iy,0)/wbbl1(ix,iy,3)      !srv 11.01.13
        if (area_fix.ge.3) pbs1(ix,iy,0) = pbs1(ix,iy,0)*qc1(ix,iy)
        pbs1(ix,iy,1) = (gs1(ix,iy,1)/t1)*wbbv1(ix,iy,1)/wbbv1(ix,iy,3) !srv 11.01.13
       enddo
      enddo
*   ..test sign of pbs1
!      call b2xvsg (n2*2, pbs1, 1, 'pbs1', '.ge.')                      !srv 11.01.13
*   ..test edge values of pbs1
      call b2xvfx (nx1, ny1, pbs1(-1,-1,0), 'pbs1x')
      call b2xvfy1(nx1, ny1, pbs1(-1,-1,1), 'pbs1y', bottomiy1)
*   ..return
      ncall = ncall + 1
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xgbs1
