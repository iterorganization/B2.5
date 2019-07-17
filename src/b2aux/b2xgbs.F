      subroutine b2xgbs (nx, ny, wbbl, wbbv, gs, qac, pbs,
     &  leftix, leftiy, bottomix, bottomiy, cflag)
      use b2mod_types
      use b2mod_cellhelper
      use b2mod_connectivity
      use carre_constants
      implicit none
      integer nx, ny
      integer cflag(-1:nx, -1:ny, CARREOUT_NCELLFLAGS),
     &  leftix(-1:nx,-1:ny), leftiy(-1:nx,-1:ny), 
     &  bottomix(-1:nx,-1:ny), bottomiy(-1:nx,-1:ny)
      real (kind=R8) ::
     *  wbbl(-1:nx,-1:ny,0:3), wbbv(-1:nx,-1:ny,0:3),
     &  gs(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1), qac(-1:nx,-1:ny,0:1)
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
      external subini, subend, xertst, b2xvsg, b2xvff, ipgeti
      save ncall
      data ncall/0/
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
*   ..test input
      call xertst (0.le.nx.and.0.le.ny, 'faulty input nx, ny')
*   ..test sign of input arrays
      n2 = (nx+2)*(ny+2)
cwdk  test for wbb.(:,:,3) switched to .ge.; could be zero in dead cells
      call b2xvsg (n2, wbbl(-1,-1,3), 1, 'wbbl', '.ge.')
      call b2xvsg (n2, wbbv(-1,-1,3), 1, 'wbbv', '.ge.')
*   ..test edge values of gs
      call b2xvfx (nx, ny, gs(-1,-1,0), 'gs')
      call b2xvfy (nx, ny, gs(-1,-1,1), 'gs')
*   ..compute pbs(,,0:1)
      do iy = -1, ny
       do ix = -1, nx
        if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
        if(isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) 
     &    pbs(ix,iy,0) = gs(ix,iy,0)*qac(ix,iy,0)*wbbl(ix,iy,0)/
     &      wbbl(ix,iy,3)
        if(isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) 
     &    pbs(ix,iy,1) = gs(ix,iy,1)*qac(ix,iy,1)*wbbv(ix,iy,0)/
     &      wbbv(ix,iy,3)
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


