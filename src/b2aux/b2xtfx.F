      subroutine b2xtfx (nx, ny, nx1, ny1, fx, fx1, 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hy, qc, hy1, qc1)
      use b2mod_types
      implicit none
      integer nx, ny, nx1, ny1, nncutmax
      real (kind=R8) ::
     *   fx(-1:nx,-1:ny), fx1(-1:nx1,-1:ny1),
     &   hy(-1:nx,-1:ny), hy1(-1:nx1,-1:ny1),
     &   qc(-1:nx,-1:ny), qc1(-1:nx1,-1:ny1)
      integer non_commensurate
      integer nnreg(0:2), nncut, periodic_bc, jsep, jsep1
      integer region(-1:nx,-1:ny,0:2), region1(-1:nx1,-1:ny1,0:2)
      integer leftcut(nncutmax), leftcut1(nncutmax)
      integer rightcut(nncutmax), rightcut1(nncutmax)
      integer bottomcut(nncutmax), bottomcut1(nncutmax)
      integer topcut(nncutmax), topcut1(nncutmax)
*     ------------------------------------------------------------------
*     B2XTFX transfers an x-face flux to a new mesh.
*     ------------------------------------------------------------------
      integer ix, iy, ix1, iy1, ixa, ixb, iya, iyb 
      integer ixc, ixd, iyc, iyd, iyp, iym, ixp, nd, ifail
      integer inseliy, inseliy1
      real (kind=R8), allocatable :: d(:)
      real (kind=R8) :: xx(-1:nx), xx1(-1:nx1), yy(-1:ny), yy1(-1:ny1) 
      real (kind=R8) :: fa(-1:nx,-1:ny), fa1(-1:nx1,-1:ny1)
      real (kind=R8) :: rratio
      external b2xvfx, xerrab, subini, subend, interp2d, rratio
      external E01BEF, E01BFF
      intrinsic min, max
*     ------------------------------------------------------------------
      call subini ('b2xtfx')
      if (non_commensurate.eq.0) then   ! commensurate mapping
*   ..transfer fx
       do iy1 = -1, ny1
        if (iy1.eq.-1) then
         iya = -1
        else
         iya = 2*ny-1-(4*ny1-1-2*iy1)*ny/(2*ny1)
        endif
        if (iy1.eq.ny1) then
         iyb = ny
        else
         iyb = (2*iy1+2*ny1+1)*ny/(2*ny1)-ny
        endif
        fx1(-1,iy1) = 0.0_R8
        do ix1 = 0, nx1
         ixa = nx-(2*(nx1-ix1)*nx+nx1)/(2*nx1)
         ixb = (2*ix1*nx+nx1)/(2*nx1)
         fx1(ix1,iy1) =
     &    ny*(fx(ixa,iya)+fx(ixb,iya)+fx(ixa,iyb)+fx(ixb,iyb))/
     &    (4.0_R8*ny1)
        enddo
       enddo
      else  ! non-commensurate conformal mapping
       nd=max(nx,nx1,ny,ny1)
       allocate(d(nd+2))
       ifail=0
*    ..all regions : all topologies : western edge
       fx1(-1,-1:ny1) = 0.0_R8
*   ..we proceed per region
*    ..region 1 : 
       if (nncut.eq.0 .and. nx.eq.1) then  ! 1-D radial case
        yy(-1) = 0.0_R8
        do iy = 0, ny-1
          yy(iy) = (iy+0.5_R8)/real(ny)
        enddo
        yy(ny) = 1.0_R8
        do iy1 = 0, ny1-1
          yy1(iy1) = (iy1+0.5_R8)/real(ny1)
        enddo
*     ..region 1 : 1-D radial case : central column (including top and bottom rows)
        fx1(0:nx1-1,-1) = fx(0,-1)
        fa(0,-1) = fx(0,0)/(hy(0,0)*qc(0,0))
        fa(0,0:ny-1) = fx(0,0:ny-1)/(hy(0,0:ny-1)*qc(0,0:ny-1))
        fa(0,ny) = fx(0,ny-1)/(hy(0,ny-1)*qc(0,ny-1))
        call E01BEF(ny+2,yy(-1:ny),fa(0,-1:ny),d,ifail)
        call E01BFF(ny+2,yy(-1:ny),fa(0,-1:ny),d,
     .              ny1,yy1(0:ny1-1),fa1(0,0:ny1-1),ifail)
        fx1(0,0:ny1-1) = fa1(0,0:ny1-1)*(hy1(0,0:ny1-1)*qc1(0,0:ny1-1))
        fx1(0,ny1) = fx(0,ny) 
        do ix1 = 1, nx1-1
         fx1(ix1,0:ny1) = fx1(0,0:ny1)
        enddo
*     ..region 1 : 1-D radial case : eastern edge
        fx1(nx1,-1) = fx(nx,-1)
        fa(nx,-1) = fx(nx,0)/(hy(nx,0)*qc(nx,0))
        fa(nx,0:ny-1) = fx(nx,0:ny-1)/(hy(nx,0:ny-1)*qc(nx,0:ny-1))
        fa(nx,ny) = fx(nx,ny)/(hy(nx,ny)*qc(nx,ny))
        call E01BEF(ny+2,yy(-1:ny),fa(nx,-1:ny),d,ifail)
        call E01BFF(ny+2,yy(-1:ny),fa(nx,-1:ny),d,
     .              ny1,yy1(0:ny1-1),fa1(nx1,0:ny1-1),ifail)
        fx1(nx1,0:ny1-1) = 
     .              fa1(nx1,0:ny1-1)*(hy1(nx1,0:ny1-1)*qc1(nx1,0:ny1-1))
        fx1(nx1,ny1) = fx(nx,ny) 
       else if (nncut.eq.0 .and. ny.eq.1) then  ! 1-D parallel case
        do ix = 0, nx
          xx(ix) = rratio(ix,nx)
        enddo
        do ix1 = 0, nx1
          xx1(ix1) = rratio(ix1,nx1)
        enddo
*     ..region 1 : 1-D parallel case : bottom row (including eastern edge)
        call E01BEF(nx+1,xx(0:nx),fx(0:nx,-1),d,ifail)
        call E01BFF(nx+1,xx(0:nx),fx(0:nx,-1),d,
     .              nx1+1,xx1(0:nx1),fx1(0:nx1,-1),ifail)
*     ..region 1 : 1-D parallel case : center row (including eastern edge)
        call E01BEF(nx+1,xx(0:nx),fx(0:nx,0),d,ifail)
        call E01BFF(nx+1,xx(0:nx),fx(0:nx,0),d,
     .              nx1+1,xx1(0:nx1),fa1(0:nx1,0),ifail)
        do iy1 = 0, ny1-1
          fx1(0:nx1,iy1) = fa1(0:nx1,0)/real(ny1)
        enddo
*     ..region 1 : 1-D parallel case : top row (including eastern edge)
        call E01BEF(nx+1,xx(0:nx),fx(0:nx,ny),d,ifail)
        call E01BFF(nx+1,xx(0:nx),fx(0:nx,ny),d,
     .              nx1+1,xx1(0:nx1),fx1(0:nx1,ny1),ifail)
       else if (nncut.eq.0 .and. periodic_bc.eq.0) then  ! 2-D slab case
        do ix = 0, nx
          xx(ix) = rratio(ix,nx)
        enddo
        do ix1 = 0, nx1
          xx1(ix1) = rratio(ix1,nx1)
        enddo
        yy(-1)=0.0_R8
        do iy = 0, ny-1
          yy(iy) = (iy+0.5_R8)/real(ny)
        enddo
        yy(ny)=1.0_R8
        do iy1 = 0, ny1-1
          yy1(iy1) = (iy1+0.5_R8)/real(ny1)
        enddo
*     ..region 1 : 2-D slab case : bottom row (including eastern edge)
        call E01BEF(nx+1,xx(0:nx),fx(0:nx,-1),d,ifail)
        call E01BFF(nx+1,xx(0:nx),fx(0:nx,-1),d,
     .              nx1+1,xx1(0:nx1),fx1(0:nx1,-1),ifail)
*     ..region 1 : 2-D slab case : central part (including eastern edge)
        fa(0:nx,-1) = fx(0:nx,0)/(hy(0:nx,0)*qc(0:nx,0))
        fa(0:nx,0:ny-1) = 
     .   fx(0:nx,0:ny-1)/(hy(0:nx,0:ny-1)*qc(0:nx,0:ny-1))
        fa(0:nx,ny) = fx(0:nx,ny-1)/(hy(0:nx,ny-1)*qc(0:nx,ny-1))
        call interp2d(nx+1,ny+2,nx1+1,ny1,
     .   xx(0:nx),yy(-1:ny),xx1(0:nx1),yy1(0:ny1-1),
     .   fa(0:nx,-1:ny),fa1(0:nx1,0:ny1-1))
        fx1(0:nx1,0:ny1-1) = 
     .   fa1(0:nx1,0:ny1-1)*(hy1(0:nx1,0:ny1-1)*qc1(0:nx1,0:ny1-1))
*     ..region 1 : 2-D slab case : top row (including eastern edge)
        call E01BEF(nx+1,xx(0:nx),fx(0:nx,ny),d,ifail)
        call E01BFF(nx+1,xx(0:nx),fx(0:nx,ny),d,
     .              nx1+1,xx1(0:nx1),fx1(0:nx1,ny1),ifail)
       else    ! X-point topologies
        if (nncut.eq.0 .and. periodic_bc.eq.1) then  ! limiter core
         ixa = 0
         iya = 0
         ixb = nx-1
         ixp = ixa
         iyb = topcut(1)-1
         ixc = 0
         iyc = 0
         ixd = nx1-1
         iyd = topcut1(1)-1
         fx1(nx1,-1:topcut1(1)) = 0.0_R8   ! eastern dead cells
        else if (nncut.eq.1) then  ! SN geometry (core)
         ixa = leftcut(1)
         iya = 0
         ixb = rightcut(1)-1
         ixp = ixa
         iyb = topcut(1)-1
         ixc = leftcut1(1)
         iyc = 0
         ixd = rightcut1(1)-1
         iyd = topcut1(1)-1
        else if (nncut.eq.2) then  ! DN geometry (HFS core)
         ixa = leftcut(1)
         iya = 0
         ixb = leftcut(2)-1
         ixp = rightcut(2)
         iyb = min(topcut(1),topcut(2))-1
         ixc = leftcut1(1)
         iyc = 0
         ixd = leftcut1(2)-1
         iyd = min(topcut1(1),topcut(2))-1
        else
         call xerrab ('Unexpected topology !')
        endif
        do ix = ixa, ixb
          xx(ix) = rratio(ix-ixa,ixb-ixa+1)
        enddo
        xx(ixb+1) = 1.0_R8   ! for periodicity
        do ix1 = ixc, ixd
          xx1(ix1) = rratio(ix1-ixc,ixd-ixc+1)
        enddo
        yy(iya-1)=0.0_R8
        do iy = iya, iyb
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1)=1.0_R8
        do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
*     ..region 1 : core : bottom row
        fa(ixa:ixb,-1) = fx(ixa:ixb,-1)
        fa(ixb+1,-1) = fx(ixp,-1)  ! periodicity condition
        call E01BEF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,-1),d,ifail)
        call E01BFF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,-1),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,-1),ifail)
*     ..region 1 : core : central part
        fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .   (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
        fa(ixb+1,iya-1) = fx(ixp,iya-1)/
     .   (hy(ixp,iya-1)*qc(ixp,iya-1))  ! periodicity condition
        fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .   (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
        fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .   (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! periodicity condition
        fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .   (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
        fa(ixb+1,iyb+1) = fx(ixp,iyb)/
     .   (hy(ixp,iyb)*qc(ixp,iyb))  ! periodicity condition
        call interp2d(ixb-ixa+2,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa:ixb+1,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
        fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .   (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
       endif
*    ..region 2 : 
       if (nncut.eq.0 .and. periodic_bc.eq.1) then  ! limiter SOL
        do ix = 0, nx
          xx(ix) = rratio(ix,nx)
        enddo
        do ix1 = 0, nx1
          xx1(ix1) = rratio(ix1,nx1)
        enddo
        yy(topcut(1)-1) = 0.0_R8
        do iy = topcut(1), ny-1
          yy(iy) = (iy+0.5_R8-topcut(1))/real(ny-topcut(1))
        enddo
        yy(ny) = 1.0_R8
        do iy1 = topcut1(1), ny1-1
          yy1(iy1) = (iy1+0.5_R8-topcut1(1))/real(ny1-topcut1(1))
        enddo
*     ..region 2 : limiter SOL : central part (including eastern edge)
        fa(0:nx,topcut(1)-1) = fx(0:nx,topcut(1))/
     .   (hy(0:nx,topcut(1))*qc(0:nx,topcut(1)))
        fa(0:nx,topcut(1):ny-1) = fx(0:nx,topcut(1):ny-1)/
     .   (hy(0:nx,topcut(1):ny-1)*qc(0:nx,topcut(1):ny-1))
        fa(0:nx,ny) = fx(0:nx,ny-1)/(hy(0:nx,ny-1)*qc(0:nx,ny-1))
        call interp2d(nx+1,ny-topcut(1)+2,nx1+1,ny1-topcut1(1),
     .   xx(0:nx),yy(topcut(1)-1:ny),xx1(0:nx1),yy1(topcut1(1):ny1-1),
     .   fa(0:nx,topcut(1)-1:ny),fa1(0:nx1,topcut1(1):ny1-1))
        fx1(0:nx1,topcut1(1):ny1-1) = fa1(0:nx1,topcut(1):ny1-1)*
     .   (hy1(0:nx1,topcut1(1):ny1-1)*qc1(0:nx1,topcut1(1):ny1-1))
*     ..region 2 : limiter SOL : top row (including eastern edge)
        call E01BEF(nx+1,xx(0:nx),fx(0:nx,ny),d,ifail)
        call E01BFF(nx+1,xx(0:nx),fx(0:nx,ny),d,
     .              nx1+1,xx1(0:nx1),fx1(0:nx1,ny1),ifail)
       else if (nncut.eq.1 .and. periodic_bc.eq.1) then  ! stellarator island SOL
        do ix = 0, nx
          xx(ix) = rratio(ix,nx)
        enddo
        do ix1 = 0, nx1
          xx1(ix1) = rratio(ix1,nx1)
        enddo
        iy = topcut(1)
        do while (region(nx/2,iy,0).eq.2)
          inseliy = iy+1
        enddo
        yy(topcut(1)-1) = 0.0_R8
        do iy = topcut(1), inseliy-1
          yy(iy) = (iy+0.5_R8-topcut(1))/real(inseliy-topcut(1))
        enddo
        yy(inseliy) = 1.0_R8
        iy1 = topcut1(1)
        do while (region1(nx1/2,iy1,0).eq.2)
          inseliy1 = iy1+1
        enddo
        do iy1 = topcut1(1), inseliy1-1
          yy1(iy1) = (iy1+0.5_R8-topcut1(1))/real(inseliy1-topcut1(1))
        enddo
*     ..region 2 : stellarator island SOL : central part (including eastern edge)
        fa(0:nx,topcut(1)-1) = fx(0:nx,topcut(1))/
     .   (hy(0:nx,topcut(1))*qc(0:nx,topcut(1)))
        fa(0:nx,topcut(1):inseliy-1) = fx(0:nx,topcut(1):inseliy-1)/
     .   (hy(0:nx,topcut(1):inseliy-1)*qc(0:nx,topcut(1):inseliy-1))
        fa(0:nx,inseliy) = fx(0:nx,inseliy-1)/
     .   (hy(0:nx,inseliy-1)*qc(0:nx,inseliy-1))
        call interp2d(nx+1,inseliy-topcut(1)+2,
     .   nx1+1,inseliy1-topcut1(1),
     .   xx(0:nx),yy(topcut(1)-1:inseliy),
     .   xx1(0:nx1),yy1(topcut1(1):inseliy1-1),
     .   fa(0:nx,topcut(1)-1:inseliy),fa1(0:nx1,topcut1(1):inseliy1-1))
        fx1(0:nx1,topcut1(1):inseliy1-1) = 
     .   fa1(0:nx1,topcut(1):inseliy1-1)*
     .   (hy1(0:nx1,topcut1(1):inseliy1-1)*
     .    qc1(0:nx1,topcut1(1):inseliy1-1))
       else if (nncut.eq.1 .or. 
     &  (nncut.eq.2 .and. topcut(1).eq.topcut(2))) then  ! SN or CDN HFS SOL
        ixa = leftcut(1)
        iya = topcut(1)
        iyb = ny-1
        ixc = leftcut1(1)
        iyc = topcut1(1)
        iyd = ny1-1
        if (nncut.eq.1) then
          ixb = rightcut(1)-1
          ixp = rightcut(1)   ! position past the X-pt
          ixd = rightcut1(1)-1
        else if (nncut.eq.2) then
          ixb = leftcut(2)-1
          ixp = leftcut(2)
          ixd = leftcut1(2)-1
        endif
        do ix = ixa, ixb
          xx(ix) = rratio(ix-ixa,ixb-ixa+1)
        enddo
        xx(ixb+1) = 1.0_R8   ! for continuity
        do ix1 = ixc, ixd
          xx1(ix1) = rratio(ix1-ixc,ixd-ixc+1)
        enddo
        yy(iya-1) = 0.0_R8
        do iy = iya, iyb
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
*     ..region 2 : SN or CDN HFS SOL : central part
        fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .   (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
        fa(ixb+1,iya-1) = fx(ixp,iya-1)/
     .   (hy(ixp,iya-1)*qc(ixp,iya-1))  ! periodicity condition
        fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .   (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
        fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .   (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! continuity condition
        fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .   (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
        fa(ixb+1,iyb+1) = fx(ixp,iyb)/
     .   (hy(ixp,iyb)*qc(ixp,iyb))  ! periodicity condition
        call interp2d(ixb-ixa+2,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa:ixb+1,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
        fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .   (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
*     ..region 2 : SN or CDN HFS SOL : top row
        fa(ixa:ixb,ny) = fx(ixa:ixb,ny)
        fa(ixb+1,ny) = fx(ixp,ny)  ! continuity condition
        call E01BEF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,ny),d,ifail)
        call E01BFF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,ny),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,ny1),ifail)
       else if (nncut.eq.2) then  ! DDN HFS SOL
*     ..region 2 : DDN HFS SOL : between separatrices
        ixa = leftcut(1)
        iya = min(topcut(1),topcut(2))
        ixb = leftcut(2)-1
        iyb = max(topcut(1),topcut(2))-1
        ixc = leftcut1(1)
        iyc = min(topcut1(1),topcut1(2))
        ixd = leftcut1(2)-1
        iyd = max(topcut1(1),topcut1(2))-1
        if (topcut(1).lt.topcut(2)) then    ! bottom X-point is active
          ixp = rightcut(2)                 ! position on LFS
        else if (topcut(1).gt.topcut(2)) then    ! top X-point is active
          ixp = leftcut(2)                  ! position past the X-point
        endif
        do ix = ixa, ixb
          xx(ix) = rratio(ix-ixa,ixb-ixa+1)
        enddo
        xx(ixb+1) = 1.0_R8   ! for continuity
        do ix1 = ixc, ixd
          xx1(ix1) = rratio(ix1-ixc,ixd-ixc+1)
        enddo
        yy(iya-1) = 0.0_R8
        do iy = iya, iyb
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
        fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .   (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
        fa(ixb+1,iya-1) = fx(ixp,iya-1)/
     .   (hy(ixp,iya-1)*qc(ixp,iya-1))  ! periodicity condition
        fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .   (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
        fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .   (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! continuity condition
        fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .   (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
        fa(ixb+1,iyb+1) = fx(ixp,iyb)/
     .   (hy(ixp,iyb)*qc(ixp,iyb))  ! periodicity condition
        call interp2d(ixb-ixa+2,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa:ixb+1,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
        fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .   (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
*     ..region 2 : DDN HFS SOL : above second separatrix
        iya = max(topcut(1),topcut(2))
        iyb = ny-1
        iyc = max(topcut1(1),topcut1(2))
        iyd = ny1-1
        ixp = leftcut(2)                  ! position past the X-point
        fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .   (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
        fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .   (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! continuity condition
        call interp2d(ixb-ixa+2,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
        fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .   (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
*     ..region 2 : DDN HFS SOL : top row
        fa(ixa:ixb,ny) = fx(ixa:ixb,ny)
        fa(ixb+1,ny) = fx(ixp,ny)  ! continuity condition
        call E01BEF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,ny),d,ifail)
        call E01BFF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,ny),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,ny1),ifail)
       endif
*    ..region 3 : 
       if (nncut.ge.1) then
        ixa = 0
        ixb = leftcut(1)-1
        ixp = rightcut(1)                 ! position on LFS
        ixc = 0
        ixd = leftcut1(1)-1
        do ix = ixa, ixb
         xx(ix) = rratio(ix-ixa,ixb-ixa+1)
        enddo
        xx(ixb+1) = 1.0_R8   ! for continuity
        do ix1 = ixc, ixd
         xx1(ix1) = rratio(ix1-ixc,ixd-ixc+1)
        enddo
*     ..region 3 : western PFR : bottom row
        fa(ixa:ixb,-1) = fx(ixa:ixb,-1)
        fa(ixb+1,-1) = fx(ixp,-1)  ! continuity condition
        call E01BEF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,-1),d,ifail)
        call E01BFF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,-1),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,-1),ifail)
*     ..region 3 : western PFR : below the first separatrix
        iya = 0
        iyc = 0
        if (nncut.eq.1) then
         iyb = topcut(1)-1
         iyd = topcut1(1)-1
        else if (nncut.eq.2) then
         iyb = min(topcut(1),topcut(2))
         iyd = min(topcut1(1),topcut1(2))
        endif
        yy(iya-1) = 0.0_R8
        do iy = iya, iyb
         yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
        fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .   (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
        fa(ixb+1,iya-1) = fx(ixp,iya-1)/
     .   (hy(ixp,iya-1)*qc(ixp,iya-1))  ! periodicity condition
        fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .    (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
        fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .    (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! continuity condition
        fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .   (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
        fa(ixb+1,iyb+1) = fx(ixp,iyb)/
     .   (hy(ixp,iyb)*qc(ixp,iyb))  ! periodicity condition
        call interp2d(ixb-ixa+2,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa:ixb+1,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
        fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .    (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
*     ..region 3 : DDN HFS lower divertor : between separatrices
        if (nncut.eq.2 .and. topcut(1).ne.topcut(2)) then
         if (topcut(1).lt.topcut(2)) ixp = leftcut(1)   ! position past the X-point
         iya = min(topcut(1),topcut(2))
         iyb = max(topcut(1),topcut(2))-1
         iyc = min(topcut1(1),topcut1(2))
         iyd = max(topcut1(1),topcut1(2))-1
         do iy = iya, iyb
           yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         do iy1 = iyc, iyd
           yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
         fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .    (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
         fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .    (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! continuity condition
         call interp2d(ixb-ixa+2,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
         fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .    (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
        endif
*     ..region 3 : western divertor SOL : above the separatrices
        ixp = leftcut(1)   ! position past the X-point
        if (nncut.eq.1 .and. periodic_bc.eq.0) then
         iya = topcut(1)
         iyb = ny-1
         iyc = topcut1(1)
         iyd = ny1-1
        else if (nncut.eq.2) then
         iya = max(topcut(1),topcut(2))
         iyb = ny-1
         iyc = max(topcut1(1),topcut1(2))
         iyd = ny1-1
        endif
        if (nncut.ge.1 .and. periodic_bc.eq.0) then
         yy(iya-1) = 0.0_R8
         do iy = iya, iyb
           yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         yy(iyb+1) = 1.0_R8
         do iy1 = iyc, iyd
           yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
         fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .    (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
         fa(ixb+1,iya-1) = fx(ixp,iya-1)/
     .    (hy(ixp,iya-1)*qc(ixp,iya-1))  ! periodicity condition
         fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .    (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
         fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .    (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! continuity condition
         fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .    (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
         fa(ixb+1,iyb+1) = fx(ixp,iyb)/
     .    (hy(ixp,iyb)*qc(ixp,iyb))  ! periodicity condition
         call interp2d(ixb-ixa+2,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa:ixb+1,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
         fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .    (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
*    ..region 3 : western divertor SOL : top row
         fa(ixa:ixb,ny) = fx(ixa:ixb,ny)
         fa(ixb+1,ny) = fx(ixp,ny)  ! continuity condition
         call E01BEF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,ny),d,ifail)
         call E01BFF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,ny),d,
     .               ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,ny1),ifail)
        endif       
       endif
*    ..region 4 : 
       if (nncut.ge.1) then
        if (nncut.eq.1) then
         ixa = rightcut(1)
         ixb = nx
         ixc = rightcut1(1)
         ixd = nx1
        else if (nncut.eq.2) then
         ixa = rightcut(2)
         ixb = ixa
         do while (region(ixb+1,jsep,0).eq.4)
          ixb = ixb+1
         enddo 
         ixc = rightcut1(2)
         ixd = ixc
         do while (region1(ixd+1,jsep1,0).eq.4)
          ixd = ixd + 1
         enddo
        endif
        do ix = ixa, ixb
          xx(ix) = rratio(ix-ixa,ixb-ixa)
        enddo
        do ix1 = ixc, ixd
          xx1(ix1) = rratio(ix1-ixc,ixd-ixc)
        enddo
*     ..region 4 : eastern SN PFR or HFS top DN PFR : bottom row (including eastern edge)
        call E01BEF(ixb-ixa+1,xx(ixa:ixb),fx(ixa:ixb,-1),d,ifail)
        call E01BFF(ixb-ixa+1,xx(ixa:ixb),fx(ixa:ixb,-1),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,-1),ifail)
*     ..region 4 : eastern SN PFR or HFS top DN PFR : below the separatrix (including eastern edge)
        iya = 0
        iyc = 0
        if (nncut.eq.1) then
         iyb = topcut(1)-1
         iyd = topcut1(1)-1
        else if (nncut.eq.2) then
         iyb = min(topcut(1),topcut(2))-1
         iyd = min(topcut1(1),topcut1(2))-1
        endif
        yy(iya-1) = 0.0_R8
        do iy = iya, iyb
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
        fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .   (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
        fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .   (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
        fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .   (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
        call interp2d(ixb-ixa+1,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa:ixb),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa:ixb,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
        fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .   (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
*     ..region 4 : HFS DDN top divertor : between the separatrices (including eastern edge)
        if (nncut.eq.2 .and. topcut(1).ne.topcut(2)) then
         iya = min(topcut(1),topcut(2))
         iyb = max(topcut(1),topcut(2))-1
         iyc = min(topcut1(1),topcut1(2))
         iyd = max(topcut1(1),topcut1(2))-1
         yy(iya-1) = 0.0_R8
         do iy = iya, iyb
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         yy(iyb+1) = 1.0_R8
         do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
         fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .    (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
         fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .    (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
         fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .    (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
         call interp2d(ixb-ixa+1,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa:ixb),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa:ixb,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
         fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .    (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
        endif
        if (nncut.ge.1 .and. periodic_bc.eq.0) then
*     ..region 4 : western divertor SN SOL or HFS DDN top divertor SOL : above the separatrices (including eastern edge)
         if (nncut.eq.1) then
          iya = topcut(1)
          iyc = topcut1(1)
         else if (nncut.eq.2) then
          iya = max(topcut(1),topcut(2))
          iyc = max(topcut1(1),topcut1(2))
         endif
         iyb = ny-1
         iyd = ny1-1
         yy(iya-1) = 0.0_R8
         do iy = iya, iyb
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         yy(iyb+1) = 1.0_R8
         do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
         fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .    (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
         fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .    (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
         fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .    (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
         call interp2d(ixb-ixa+1,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa:ixb),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa:ixb,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
         fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .    (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
*     ..region 4 : western divertor SN SOL or HFS DDN top divertor SOL : top row (including eastern edge)
         fa(ixa:ixb,ny) = fx(ixa:ixb,ny)
         call E01BEF(ixb-ixa+1,xx(ixa:ixb),fa(ixa:ixb,ny),d,ifail)
         call E01BFF(ixb-ixa+1,xx(ixa:ixb),fa(ixa:ixb,ny),d,
     .               ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,ny1),ifail)
        endif
       endif
*    ..region 5 : 
       if (nncut.eq.1 .and. periodic_bc.eq.1) then  ! stellarator island
        ixa = 0
        iya = inseliy
        ixb = nx-1
        iyb = ny-1
        ixp = 0   ! for periodicity
        ixc = 0
        iyc = inseliy1
        ixd = nx1-1
        iyd = ny1-1
        fx1(nx1,inseliy1:ny1) = 0.0_R8   ! eastern dead cells
        do ix = ixa, ixb
         xx(ix) = rratio(ix-ixa,ixb-ixa+1)
        enddo
        xx(ixb+1) = 1.0_R8   ! for periodicity
        do ix1 = ixc, ixd
         xx1(ix1) = rratio(ix1-ixc,ixd-ixc+1)
        enddo
        yy(iya-1) = 0.0_R8
        do iy = iya, iyb
         yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
*     ..region 5 : stellarator island : central part
        fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .   (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
        fa(ixb+1,iya-1) = fx(ixp,iya-1)/
     .   (hy(ixp,iya-1)*qc(ixp,iya-1))  ! periodicity condition
        fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .   (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
        fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .   (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! periodicity condition
        fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .   (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
        fa(ixb+1,iyb+1) = fx(ixp,iyb)/
     .   (hy(ixp,iyb)*qc(ixp,iyb))  ! periodicity condition
        call interp2d(ixb-ixa+2,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa:ixb+1,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
        fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .   (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
*     ..region 5 : stellarator island : top row
        fa(ixa:ixb,ny) = fx(ixa:ixb,ny)
        fa(ixb+1,ny) = fx(ixp,ny)  ! periodicity condition
        call E01BEF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,ny),d,ifail)
        call E01BFF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,ny),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,ny1),ifail)         
       else if (nncut.eq.2) then  ! DN LFS core
        ixa = rightcut(2)
        ixb = rightcut(1)-1
        ixp = leftcut(1)  ! for continuity
        ixc = rightcut1(2)
        ixd = rightcut1(1)-1
        do ix = ixa, ixb
         xx(ix) = rratio(ix-ixa,ixb-ixa+1)
        enddo
        xx(ixb+1) = 1.0_R8   ! for continuity
        do ix1 = ixc, ixd
         xx1(ix1) = rratio(ix1-ixc,ixd-ixc+1)
        enddo
*     ..region 5 : DN LFS core : bottom row
        fa(ixa:ixb,-1) = fx(ixa:ixb,-1)
        fa(ixb+1,-1) = fx(ixp,-1)  ! continuity condition
        call E01BEF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,-1),d,ifail)
        call E01BFF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,-1),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,-1),ifail)
*     ..region 5 : DN LFS core : central part
        iya = 0
        iyb = min(topcut(1),topcut(2))-1
        iyc = 0
        iyd = min(topcut1(1),topcut1(2))-1
        yy(iya-1) = 0.0_R8
        do iy = iya, iyb
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
        fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .   (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
        fa(ixb+1,iya-1) = fx(ixp,iya-1)/
     .   (hy(ixp,iya-1)*qc(ixp,iya-1))  ! periodicity condition
        fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .   (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
        fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .   (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! continuity condition
        fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .   (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
        fa(ixb+1,iyb+1) = fx(ixp,iyb)/
     .   (hy(ixp,iyb)*qc(ixp,iyb))  ! periodicity condition
        call interp2d(ixb-ixa+2,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa:ixb+1,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
        fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .   (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
       endif
*    ..region 6 :
       if (nncut.eq.2) then
         ixa = rightcut(2)
         iya = min(topcut(1),topcut(2))
         ixb = rightcut(1)-1
         ixc = rightcut1(2)
         iya = min(topcut1(1),topcut1(2))
         ixd = rightcut1(1)-1
         if (topcut(1).eq.topcut(2)) then   ! CDN Balanced X-points
           ixp = rightcut(1)  ! position past the X-point
           iyb = ny-1
           iyd = ny1-1
         else if (topcut(1).lt.topcut(2)) then   ! Bottom X-point is active or CDN
           ixp = rightcut(1)  ! position past the X-point
           iyb = max(topcut(1),topcut(2))-1
           iyd = max(topcut1(1),topcut1(2))-1
         else if (topcut(1).gt.topcut(2)) then  ! Top X-point is active
           ixp = leftcut(1)   ! continuity to HFS
           iyb = max(topcut(1),topcut(2))-1
           iyd = max(topcut1(1),topcut1(2))-1
         endif
         do ix = ixa, ixb
          xx(ix) = rratio(ix-ixa,ixb-ixa+1)
         enddo
         xx(ixb+1) = 1.0_R8   ! for periodicity
         do ix1 = ixc, ixd
          xx1(ix1) = rratio(ix1-ixc,ixd-ixc+1)
         enddo
         yy(iya-1) = 0.0_R8
         do iy = iya, iyb
           yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         yy(iyb+1) = 1.0_R8
         do iy1 = iyc, iyd
           yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
*     ..region 6 : CDN LFS SOL or DDN LFS SOL between the separatrices
         fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .    (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
         fa(ixb+1,iya-1) = fx(ixp,iya-1)/
     .    (hy(ixp,iya-1)*qc(ixp,iya-1))  ! periodicity condition
         fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .    (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
         fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .    (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! periodicity condition
         fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .    (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
         fa(ixb+1,iyb+1) = fx(ixp,iyb)/
     .    (hy(ixp,iyb)*qc(ixp,iyb))  ! periodicity condition
         call interp2d(ixb-ixa+2,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa:ixb+1,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
         fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .    (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
         if (topcut(1).ne.topcut(2)) then
*     ..region 6 : DDN LFS SOL above the second separatrix
           ixp = rightcut(1)  ! position past the X-point
           iya = max(topcut(1),topcut(2))
           iyb = ny-1
           iyc = max(topcut1(1),topcut1(2))
           iyd = ny1-1
           yy(iya-1) = 0.0_R8
           do iy = iya, iyb
            yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
           enddo
           yy(iyb+1) = 1.0_R8
           do iy1 = iyc, iyd
            yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
           enddo
           fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .      (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
           fa(ixb+1,iya-1) = fx(ixp,iya-1)/
     .      (hy(ixp,iya-1)*qc(ixp,iya-1))  ! periodicity condition
           fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .      (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
           fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .      (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! periodicity condition
           fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .      (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
           fa(ixb+1,iyb+1) = fx(ixp,iyb)/
     .      (hy(ixp,iyb)*qc(ixp,iyb))  ! periodicity condition
           call interp2d(ixb-ixa+2,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .      xx(ixa:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .      fa(ixa:ixb+1,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
           fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .      (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
         endif
*     ..region 6 : DN LFS SOL : top row
         fa(ixa:ixb,ny) = fx(ixa:ixb,ny)
         fa(ixb+1,ny) = fx(ixp,ny)  ! continuity condition
         call E01BEF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,ny),d,ifail)
         call E01BFF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,ny),d,
     .               ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,ny1),ifail)
       endif
*    ..region 7 : 
       if (nncut.eq.2) then
        ixb = rightcut(2)-1
        ixa = ixb
        do while (region(ixa-1,jsep,0).eq.7)
          ixa = ixa-1
        enddo
        ixp = leftcut(2)                 ! position on HFS
        ixd = rightcut1(2)-1
        ixc = ixd
        do while (region1(ixc-1,jsep1,0).eq.7)
          ixc = ixc-1
        enddo
        do ix = ixa, ixb
         xx(ix) = rratio(ix-ixa,ixb-ixa+1)
        enddo
        xx(ixb+1) = 1.0_R8   ! for continuity
        do ix1 = ixc, ixd
         xx1(ix1) = rratio(ix1-ixc,ixd-ixc+1)
        enddo
*     ..region 7 : top LFS PFR : bottom row
        fa(ixa:ixb,-1) = fx(ixa:ixb,-1)
        fa(ixb+1,-1) = fx(ixp,-1)  ! continuity condition
        call E01BEF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,-1),d,ifail)
        call E01BFF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,-1),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,-1),ifail)
*     ..region 7 : top LFS PFR : below the first separatrix
        iya = 0
        iyb = min(topcut(1),topcut(2))
        iyc = 0
        iyd = min(topcut1(1),topcut1(2))
        yy(iya-1) = 0.0_R8
        do iy = iya, iyb
         yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
        fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .   (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
        fa(ixb+1,iya-1) = fx(ixp,iya-1)/
     .   (hy(ixp,iya-1)*qc(ixp,iya-1))  ! periodicity condition
        fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .    (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
        fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .    (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! continuity condition
        fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .   (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
        fa(ixb+1,iyb+1) = fx(ixp,iyb)/
     .   (hy(ixp,iyb)*qc(ixp,iyb))  ! periodicity condition
        call interp2d(ixb-ixa+2,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa:ixb+1,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
        fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .    (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
*     ..region 7 : DDN LFS upper divertor : between separatrices
        if (nncut.eq.2 .and. topcut(1).ne.topcut(2)) then
         if (topcut(1).gt.topcut(2)) ixp = rightcut(2)   ! position past the X-point
         iya = min(topcut(1),topcut(2))
         iyb = max(topcut(1),topcut(2))-1
         iyc = min(topcut1(1),topcut1(2))
         iyd = max(topcut1(1),topcut1(2))-1
         yy(iya-1) = 0.0_R8
         do iy = iya, iyb
           yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         yy(iyb+1) = 1.0_R8
         do iy1 = iyc, iyd
           yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
         fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .    (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
         fa(ixb+1,iya-1) = fx(ixp,iya-1)/
     .    (hy(ixp,iya-1)*qc(ixp,iya-1))  ! periodicity condition
         fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .    (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
         fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .    (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! continuity condition
         fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .    (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
         fa(ixb+1,iyb+1) = fx(ixp,iyb)/
     .    (hy(ixp,iyb)*qc(ixp,iyb))  ! periodicity condition
         call interp2d(ixb-ixa+2,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa:ixb+1,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
         fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .    (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
        endif
*     ..region 7 : DDN LFS upper divertor : above the second separatrix
        ixp = rightcut(2)   ! position past the X-point
        iya = max(topcut(1),topcut(2))
        iyb = ny-1
        iyc = max(topcut1(1),topcut1(2))
        iyd = ny1-1
        yy(iya-1) = 0.0_R8
        do iy = iya, iyb
         yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
        fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .   (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
        fa(ixb+1,iya-1) = fx(ixp,iya-1)/
     .   (hy(ixp,iya-1)*qc(ixp,iya-1))  ! periodicity condition
        fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .   (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
        fa(ixb+1,iya:iyb) = fx(ixp,iya:iyb)/
     .   (hy(ixp,iya:iyb)*qc(ixp,iya:iyb))  ! continuity condition
        fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .   (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
        fa(ixb+1,iyb+1) = fx(ixp,iyb)/
     .   (hy(ixp,iyb)*qc(ixp,iyb))  ! periodicity condition
        call interp2d(ixb-ixa+2,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa:ixb+1,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
        fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .   (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
*    ..region 7 : LFS upper divertor SOL : top row
         fa(ixa:ixb,ny) = fx(ixa:ixb,ny)
         fa(ixb+1,ny) = fx(ixp,ny)  ! continuity condition
         call E01BEF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,ny),d,ifail)
         call E01BFF(ixb-ixa+2,xx(ixa:ixb+1),fa(ixa:ixb+1,ny),d,
     .               ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,ny1),ifail) 
       endif
*    ..region 8 : 
       if (nncut.eq.2) then
        ixa = rightcut(1)
        ixb = nx
        ixc = rightcut1(1)
        ixd = nx1
        do ix = ixa, ixb
          xx(ix) = rratio(ix-ixa,ixb-ixa)
        enddo
        do ix1 = ixc, ixd
          xx1(ix1) = rratio(ix1-ixc,ixd-ixc)
        enddo
*     ..region 8 : eastern DN LFS PFR : bottom row (including eastern edge)
        call E01BEF(ixb-ixa+1,xx(ixa:ixb),fx(ixa:ixb,-1),d,ifail)
        call E01BFF(ixb-ixa+1,xx(ixa:ixb),fx(ixa:ixb,-1),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,-1),ifail)
*     ..region 8 : eastern DN LFS PFR : below the separatrix (including eastern edge)
        iya = 0
        iyb = min(topcut(1),topcut(2))-1
        iyc = 0
        iyb = min(topcut1(1),topcut1(2))-1
        yy(iya-1) = 0.0_R8
        do iy = iya, iyb
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
        fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .   (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
        fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .   (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
        fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .   (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
        call interp2d(ixb-ixa+1,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa:ixb),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa:ixb,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
        fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .   (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
*     ..region 8 : eastern DN LFS divertor : between the separatrices (including eastern edge)
        if (topcut(1).ne.topcut(2)) then
         iya = min(topcut(1),topcut(2))
         iyb = max(topcut(1),topcut(2))-1
         iyc = min(topcut1(1),topcut1(2))
         iyd = max(topcut1(1),topcut1(2))-1
         yy(iya-1) = 0.0_R8
         do iy = iya, iyb
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         yy(iyb+1) = 1.0_R8
         do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
         fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .    (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
         fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .    (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
         fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .    (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
         call interp2d(ixb-ixa+1,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa:ixb),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa:ixb,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
         fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .    (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))  
        endif
*     ..region 8 : eastern DN LFS divertor : above the second separatrix (including eastern edge)
        iya = max(topcut(1),topcut(2))
        iyb = ny-1
        iyc = max(topcut1(1),topcut1(2))
        iyd = ny1-1
        yy(iya-1) = 0.0_R8
        do iy = iya, iyb
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
        fa(ixa:ixb,iya-1) = fx(ixa:ixb,iya)/
     .   (hy(ixa:ixb,iya)*qc(ixa:ixb,iya))
        fa(ixa:ixb,iya:iyb) = fx(ixa:ixb,iya:iyb)/
     .    (hy(ixa:ixb,iya:iyb)*qc(ixa:ixb,iya:iyb))
        fa(ixa:ixb,iyb+1) = fx(ixa:ixb,iyb)/
     .   (hy(ixa:ixb,iyb)*qc(ixa:ixb,iyb))
        call interp2d(ixb-ixa+1,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa:ixb),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa:ixb,iya-1:iyb+1),fa1(ixc:ixd,iyc:iyd))
        fx1(ixc:ixd,iyc:iyd) = fa1(ixc:ixd,iyc:iyd)*
     .    (hy1(ixc:ixd,iyc:iyd)*qc1(ixc:ixd,iyc:iyd))
*     ..region 8 : eastern DN LFS divertor : top row (including eastern edge)
        fa(ixa:ixb,ny) = fx(ixa:ixb,ny)
        call E01BEF(ixb-ixa+1,xx(ixa:ixb),fa(ixa:ixb,ny),d,ifail)
        call E01BFF(ixb-ixa+1,xx(ixa:ixb),fa(ixa:ixb,ny),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fx1(ixc:ixd,ny1),ifail)
       endif
*    ..clean-up
       deallocate(d)
      endif
*   ..test edge values of fx1
      call b2xvfx (nx1, ny1, fx1, 'b2xtfx--fx1')

      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xtfx
