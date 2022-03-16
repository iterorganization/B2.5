      subroutine b2xgbs (nx, ny, wbbl, wbbv, gs, qc, crx, cry, pbs)
      use b2mod_types
      use b2mod_subsys
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *  wbbl(-1:nx,-1:ny,0:3), wbbv(-1:nx,-1:ny,0:3),
     &  gs(-1:nx,-1:ny,0:1), crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     &  pbs(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny,0:1)
*     ------------------------------------------------------------------
*     B2XGBS computes the product (magnetic field pitch)*(surface area)
*     at cell faces. pbs(,,0) will contain the left face value and
*     pbs(,,1) will contain the bottom face value.
*     pbs thus contains the parallel contact area.
*     ------------------------------------------------------------------
      integer n2, ix, iy, ncall, area_fix
      real (kind=R8) ::
     *   t1, p0x, p0y, p1x, p1y
      intrinsic abs, sqrt, min, max
      external xertst, b2xvsg, b2xvff, ipgeti
      real (kind=R8) ::
     *   trim1, norm, cosang
      norm(p0x,p0y) = sqrt(p0x**2+p0y**2)
      trim1(p0x) = max(-1.0_R8,min(1.0_R8,p0x))
      cosang(p0x,p0y,p1x,p1y) = trim1(
     &  (p0x*p1y-p0y*p1x)/(norm(p0x,p0y)*norm(p1x,p1y)))
      save ncall, area_fix
      data ncall/0/, area_fix/3/
*-----------------------------------------------------------------------
*.documentation-internal
*
*     Let (p0x,p0y) specify a vector p0, and (p1x,p1y) specify
*     a vector p1. Then norm(p0x,p0y) is the standard norm of p0,
*     and cosang(p0x,p0y,p1x,p1y) is cos(t), where t is the complementary 
*     angle between p0 and p1.
*     t=pi/2 when p0 and p1 are parallel, t=0 or t=pi when p0 and p1
*     are orthogonal, and t=-pi/2 when p0 and p1 are anti-parallel.
*
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
      call b2xvfx (nx, ny, gs(-1,-1,0), 'gs')
      call b2xvfy (nx, ny, gs(-1,-1,1), 'gs')
*   ..compute pbs(,,0:1)
      do iy = -1, ny
       do ix = -1, nx
        t1 = cosang (
     &     crx(ix,iy,1)-crx(ix,iy,0),  cry(ix,iy,1)-cry(ix,iy,0),
     &    (crx(ix,iy,2)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,1)),
     &    (cry(ix,iy,2)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,1)))
        pbs(ix,iy,0) = gs(ix,iy,0)*wbbl(ix,iy,0)/wbbl(ix,iy,3)        !srv 11.01.13
        if (area_fix.ge.3) pbs(ix,iy,0) = pbs(ix,iy,0)*qc(ix,iy,0)
        pbs(ix,iy,1) = (gs(ix,iy,1)/t1)*wbbv(ix,iy,1)/wbbv(ix,iy,3)   !srv 11.01.13
       enddo
      enddo
*   ..test sign of pbs
!      call b2xvsg (n2*2, pbs, 1, 'pbs', '.ge.')                      !srv 11.01.13
*   ..test edge values of pbs
      call b2xvfx (nx, ny, pbs(-1,-1,0), 'pbs')
      call b2xvfy (nx, ny, pbs(-1,-1,1), 'pbs')
*   ..return
      ncall = ncall + 1
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xgbs

      subroutine b2xgbs1 (nx1, ny1, 
     . bottomiy1, wbbl1, wbbv1, gs1, qc1, crx1, cry1, pbs1)
      use b2mod_types
      use b2mod_subsys
      implicit none
      integer nx1, ny1
      integer bottomiy1(-1:nx1,-1:ny1)
      real (kind=R8) ::
     *  wbbl1(-1:nx1,-1:ny1,0:3), wbbv1(-1:nx1,-1:ny1,0:3),
     &  gs1(-1:nx1,-1:ny1,0:1), crx1(-1:nx1,-1:ny1,0:3), 
     &  cry1(-1:nx1,-1:ny1,0:3), qc1(-1:nx1,-1:ny1,0:1),
     &  pbs1(-1:nx1,-1:ny1,0:1)
*     ------------------------------------------------------------------
*     B2XGBS1 computes the product (magnetic field pitch)*(surface area)
*     at cell faces on the basis mesh. 
*     pbs1(,,0) will contain the left face value and
*     pbs1(,,1) will contain the bottom face value.
*     pbs1 thus contains the parallel contact area.
*     ------------------------------------------------------------------
      integer n2, ix, iy, ncall, area_fix
      real (kind=R8) ::
     *   t1, p0x, p0y, p1x, p1y
      intrinsic abs, sqrt, min, max
      external xertst, b2xvsg, b2xvfx, b2xvfy1, ipgeti
      real (kind=R8) ::
     *   trim1, norm, cosang
      norm(p0x,p0y) = sqrt(p0x**2+p0y**2)
      trim1(p0x) = max(-1.0_R8,min(1.0_R8,p0x))
      cosang(p0x,p0y,p1x,p1y) = trim1(
     &  (p0x*p1y-p0y*p1x)/(norm(p0x,p0y)*norm(p1x,p1y)))
      save ncall, area_fix
      data ncall/0/, area_fix/3/
*-----------------------------------------------------------------------
*.documentation-internal
*
*     Let (p0x,p0y) specify a vector p0, and (p1x,p1y) specify
*     a vector p1. Then norm(p0x,p0y) is the standard norm of p0,
*     and cosang(p0x,p0y,p1x,p1y) is cos(t), where t is the complementary 
*     angle between p0 and p1.
*     t=pi/2 when p0 and p1 are parallel, t=0 or t=pi when p0 and p1
*     are orthogonal, and t=-pi/2 when p0 and p1 are anti-parallel.
*
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
*   ..compute pbs1(,,0:1)
      do iy = -1, ny1
       do ix = -1, nx1
        t1 = cosang (
     &     crx1(ix,iy,1)-crx1(ix,iy,0),  cry1(ix,iy,1)-cry1(ix,iy,0),
     &    (crx1(ix,iy,2)+crx1(ix,iy,3))-(crx1(ix,iy,0)+crx1(ix,iy,1)),
     &    (cry1(ix,iy,2)+cry1(ix,iy,3))-(cry1(ix,iy,0)+cry1(ix,iy,1)))
        pbs1(ix,iy,0) = gs1(ix,iy,0)*wbbl1(ix,iy,0)/wbbl1(ix,iy,3)      !srv 11.01.13
        if (area_fix.ge.3) pbs1(ix,iy,0) = pbs1(ix,iy,0)*qc1(ix,iy,0)
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
