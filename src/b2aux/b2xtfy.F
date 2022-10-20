      subroutine b2xtfy (nx, ny, nx1, ny1, fy, fy1,
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hx, hx1, bottomiy1)
      use b2mod_types
      use b2mod_subsys
      implicit none
      integer nx, ny, nx1, ny1, nncutmax
      real (kind=R8) ::
     *   fy(-1:nx,-1:ny), fy1(-1:nx1,-1:ny1),
     *   hx(-1:nx,-1:ny), hx1(-1:nx1,-1:ny1)
      integer non_commensurate
      integer nnreg(0:2), nncut, periodic_bc, jsep, jsep1
      integer region(-1:nx,-1:ny,0:2), region1(-1:nx1,-1:ny1,0:2)
      integer leftcut(nncutmax), leftcut1(nncutmax)
      integer rightcut(nncutmax), rightcut1(nncutmax)
      integer bottomcut(nncutmax), bottomcut1(nncutmax)
      integer topcut(nncutmax), topcut1(nncutmax)
      integer bottomiy1(-1:nx1,-1:ny1)
*     ------------------------------------------------------------------
*     B2XTFY transfers a y-face flux to a new mesh.
*     ------------------------------------------------------------------
      integer ix, iy, ix1, iy1, ixa, ixb, iya, iyb, nd, ifail
      integer ixm, ixp, ixc, ixd, iyc, iyd, inseliy, inseliy1
      real (kind=R8), allocatable :: d(:)
      real (kind=R8) :: xx(-1:nx), xx1(-1:nx1), yy(-1:ny), yy1(-1:ny1)
      real (kind=R8) :: fa(-1:nx,-1:ny), fa1(-1:nx1,-1:ny1)
      real (kind=R8) :: rratio
      external b2xvfy1, xerrab, interp2d, rratio
      external E01BEF, E01BFF
      intrinsic min, max
*     ------------------------------------------------------------------
      call subini ('b2xtfy')
      if (non_commensurate.eq.0) then   ! commensurate mapping
*   ..transfer fy
       do ix1 = -1, nx1
        fy1(ix1,-1) = 0.0_R8
       enddo
       do iy1 = 0, ny1
        iya = ny-(2*(ny1-iy1)*ny+ny1)/(2*ny1)
        iyb = (2*iy1*ny+ny1)/(2*ny1)
        do ix1 = -1, nx1
         if (ix1.eq.-1) then
          ixa = -1
         else if (ix1.eq.nx1) then
          ixa = nx
         else
          ixa = 2*nx-1-(4*nx1-1-2*ix1)*nx/(2*nx1)
         endif
         if (ix1.eq.nx1) then
          ixb = nx
         else if (ix1.eq.-1) then
          ixb = -1
         else
          ixb = (2*ix1+2*nx1+1)*nx/(2*nx1)-nx
         endif
         fy1(ix1,iy1) =
     &    nx*(fy(ixa,iya)+fy(ixb,iya)+fy(ixa,iyb)+fy(ixb,iyb))/
     &    (4.0_R8*nx1)
        enddo
       enddo
      else  ! non-commensurate conformal mapping
       nd=max(nx,nx1,ny,ny1)
       allocate(d(nd+2))
       ifail=0
*    ..all regions : all topologies : bottom row
       fy1(-1:nx1,-1) = 0.0_R8
*   ..we proceed per region
*    .. region 1 :
       if (nncut.eq.0 .and. nx.eq.1) then  ! 1-D radial case
        do iy = 0, ny
          yy(iy) = rratio(iy,ny)
        enddo
        do iy1 = 0, ny1
          yy1(iy1) = rratio(iy1,ny1)
        enddo
*     ..region 1 : 1-D radial case : western edge (including top row)
        fa(-1,0:ny) = fy(-1,0:ny)/hx(-1,0:ny)
        call E01BEF(ny+1,yy(0:ny),fa(-1,0:ny),d,ifail)
        call E01BFF(ny+1,yy(0:ny),fa(-1,0:ny),d,
     .              ny1+1,yy1(0:ny1),fa1(-1,0:ny1),ifail)
        fy1(-1,0:ny1) = fa1(-1,0:ny1)*hx(-1,0:ny1)
*     ..region 1 : 1-D radial case : central column (including top row)
        fa(0,0:ny) = fy(0,0:ny)/hx(0,0:ny)
        call E01BEF(ny+1,yy(0:ny),fa(0,0:ny),d,ifail)
        call E01BFF(ny+1,yy(0:ny),fa(0,0:ny),d,
     .              ny1+1,yy1(0:ny1),fa1(0,0:ny1),ifail)
        do ix1 = 0, nx1-1
          fy1(ix1,0:ny1) = fa1(0,0:ny1)/real(nx1)
        enddo
*     ..region 1 : 1-D radial case : eastern edge (including top row)
        fa(nx,0:ny) = fy(nx,0:ny)/hx(nx,0:ny)
        call E01BEF(ny+1,yy(0:ny),fa(nx,0:ny),d,ifail)
        call E01BFF(ny+1,yy(0:ny),fa(nx,0:ny),d,
     .              ny1+1,yy1(0:ny1),fa1(nx1,0:ny1),ifail)
        fy1(nx1,0:ny1) = fa1(nx1,0:ny1)*hx(nx1,0:ny1)
       else if (nncut.eq.0 .and. ny.eq.1) then  ! 1-D parallel case
        xx(-1) = 0.0_R8
        do ix = 0, nx-1
          xx(ix) = (ix+0.5_R8)/real(nx)
        enddo
        xx(nx) = 1.0_R8
        do ix1 = 0, nx1-1
          xx1(ix1) = (ix1+0.5_R8)/real(nx1)
        enddo
*     ..region 1 : 1-D parallel case : western edge (including top row)
        fy1(-1,0:ny1-1) = fy(-1,0)
        fy1(-1,ny1) = fy(-1,ny)
*     ..region 1 : 1-D parallel case : central part
        fa(-1,0) = fy(0,0)/hx(0,0)
        fa(0:nx-1,0) = fy(0:nx-1,0)/hx(0:nx-1,0)
        fy(nx,0) = fy(nx-1,0)/hx(nx-1,0)
        call E01BEF(nx+2,xx(-1:nx),fa(-1:nx,0),d,ifail)
        call E01BFF(nx+2,xx(-1:nx),fa(-1:nx,0),d,
     .              nx1,xx1(0:nx1-1),fa1(0:nx1-1,0),ifail)
        fy1(0:nx1-1,0) = fa1(0:nx1-1,0)*hx1(0:nx1-1,0)
        do iy1 = 1, ny1-1
          fy1(0:nx1-1,iy1) = fy1(0:nx1-1,0)
        enddo
*     ..region 1 : 1-D parallel case : top row
        fa(-1,ny) = fy(0,ny)/hx(0,ny)
        fa(0:nx-1,ny) = fy(0:nx-1,ny)/hx(0:nx-1,ny)
        fy(nx,ny) = fy(nx-1,ny)/hx(nx-1,ny)
        call E01BEF(nx+2,xx(-1:nx),fa(-1:nx,ny),d,ifail)
        call E01BFF(nx+2,xx(-1:nx),fa(-1:nx,ny),d,
     .              nx1,xx1(0:nx1-1),fa1(0:nx1-1,ny1),ifail)
        fy1(0:nx1-1,ny1) = fa1(0:nx1-1,ny1)*hx1(0:nx1-1,ny1)
*     ..region 1 : 1-D parallel case : eastern edge
        fy1(nx1,0:ny1-1) = fy(nx,0)
        fy1(nx1,ny1) = fy(nx,ny)
       else if (nncut.eq.0 .and. periodic_bc.eq.0) then  ! 2-D slab case
        do ix = 0, nx-1
          xx(ix) = (ix+0.5_R8)/real(nx)
        enddo
        xx(nx) = 1.0_R8
        do ix1 = 0, nx1-1
          xx1(ix1) = (ix1+0.5_R8)/real(nx1)
        enddo
        do iy = 0, ny
          yy(iy) = rratio(iy,ny)
        enddo
        do iy1 = 0, ny1
          yy1(iy1) = rratio(iy1,ny1)
        enddo
*     ..region 1 : 2-D slab case : western edge (includes top row)
        call E01BEF(ny+1,yy(0:ny),fy(-1,0:ny),d,ifail)
        call E01BFF(ny+1,yy(0:ny),fy(-1,0:ny),d,
     .              ny1+1,yy1(0:ny1),fy1(-1,0:ny1),ifail)
*     ..region 1 : 2-D slab case : central part (includes top row)
        fa(0:nx,0:ny) = fy(0:nx,0:ny)/hx(0:nx,0:ny)
        call interp2d(nx+1,ny+1,nx1+1,ny1+1,
     .   xx(0:nx),yy(0:ny),xx1(0:nx1),yy1(0:ny1),
     .   fa(0:nx,0:ny),fa1(0:nx1,0:ny1))
        fy1(0:nx1-1,0:ny1) = fa1(0:nx1-1,0:ny1)/hx1(0:nx1-1,0:ny1)
*     ..region 1 : 2-D slab case : eastern edge (includes top row)
        call E01BEF(ny+1,yy(0:ny),fy(nx,0:ny),d,ifail)
        call E01BFF(ny+1,yy(0:ny),fy(nx,0:ny),d,
     .              ny1+1,yy1(0:ny1),fy1(nx1,0:ny1),ifail)
       else  ! X-point topologies
        iya = 0
        iyb = jsep+1
        iyc = 0
        iyd = jsep1+1
        do iy = iya, iyb
          yy(iy) = rratio(iy-iya,iyb-iya)
        enddo
        do iy1 = iyc, iyd
          yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
        enddo
        if (nncut.eq.0 .and. periodic_bc.eq.1) then ! limiter core
          ixa = 0
          ixb = nx-1
          ixm = nx-1
          ixp = 0
          ixc = 0
          ixd = nx1-1
          fy1(-1,iyc:iyd) = 0.0_R8   ! western dead cells
          fy1(ny1,iyc:iyd) = 0.0_R8  ! eastern dead cells
        else if (nncut.eq.1 .or. nnreg(0).eq.7) then   ! SN core
          ixa = leftcut(1)
          ixb = rightcut(1)-1
          ixm = rightcut(1)-1
          ixp = leftcut(1)
          ixc = leftcut1(1)
          ixd = rightcut1(1)-1
        else if (nncut.eq.2) then   ! DN HFS core
          ixa = leftcut(1)
          ixb = leftcut(2)-1
          ixm = rightcut(1)-1
          ixp = rightcut(2)
          ixc = leftcut1(1)
          ixd = leftcut1(2)-1
        else
          call xerrab ('Unexpected topology !')
        endif
        do ix = ixa-1, ixb+1
          xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        do ix1 = ixc, ixd
          xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
*     ..region 1 : core : all
        fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
        fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
        fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
        fy1(ixc:ixd,iyc:iyd) =
     .   fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
       endif
*    ..region 2 :
       if (nncut.ge.1 .or. periodic_bc.eq.1) then
        iya = jsep+1
        iyc = jsep1+1
        if (nncut.eq.0 .and. periodic_bc.eq.1) then   ! limiter SOL
          ixa = 0
          ixb = nx-1
          ixc = 0
          ixd = nx1-1
          ixm = 0
          ixp = nx-1
          iyb = ny
          iyd = ny1
        else if (nncut.eq.2 .and. topcut(1).ne.topcut(2)) then   ! DDN HFS SOL
          if (nnreg(0).eq.7) then
            ixa = leftcut(1)
            ixm = leftcut(1)-1
            ixc = leftcut1(1)
            ixb = rightcut(1)-1
            ixp = rightcut(1)
            ixd = rightcut1(1)-1
          else
            ixa = leftcut(1)
            ixb = leftcut(2)-1
            if (topcut(1).lt.topcut(2)) then   ! Bottom X-point is active
              ixm = leftcut(1)-1
              ixp = rightcut(2)
            else if (topcut(1).gt.topcut(2)) then  ! Top X-point is active
              ixm = rightcut(1)-1
              ixp = leftcut(2)
            endif
            ixc = leftcut1(1)
            ixd = leftcut1(2)-1
          endif
          iyb = max(topcut(1),topcut(2))
          iyd = max(topcut1(1),topcut1(2))
        else if (nncut.eq.1 .and. periodic_bc.eq.1) then  ! stellarator island SOL
          ixa = 0
          ixb = nx-1
          ixm = 0
          ixp = nx-1
          ixc = 0
          ixd = nx1-1
          inseliy = iya
          do while (region(nx/2,inseliy+1,0).eq.2)
            inseliy = inseliy + 1
          enddo
          iyb = inseliy
          inseliy1 = iyc
          do while (region1(nx1/2,inseliy1+1,0).eq.2)
            inseliy1 = inseliy1 + 1
          enddo
          iyd = inseliy1
        else if (nncut.ge.1) then  ! CDN HFS SOL or SN SOL
          iyb = ny
          iyd = ny1
          ixa = leftcut(1)
          ixm = leftcut(1)-1
          ixc = leftcut1(1)
          if (nncut.eq.1) then   ! SN
            ixb = rightcut(1)-1
            ixp = rightcut(1)
            ixd = rightcut1(1)-1
          else if (nncut.eq.2) then  ! DN
            ixb = leftcut(2)-1
            ixp = leftcut(2)
            ixd = leftcut1(2)-1
          endif
        endif
        do iy = iya, iyb
          yy(iy) = rratio(iy-iya,iyb-iya)
        enddo
        do iy1 = iyc, iyd
          yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
        enddo
        do ix = ixa-1, ixb+1
          xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        do ix1 = ixc, ixd
          xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        if (periodic_bc.eq.1) then !
*     ..region 2 : SOL : western edge
         call E01BEF(iyb-iya+1,yy(iya:iyb),fy(-1,iya:iyb),d,ifail)
         call E01BFF(iyb-iya+1,yy(iya:iyb),fy(-1,iya:iyb),d,
     .               iyd-iyc+1,yy1(iyc:iyd),fy1(-1,iyc:iyd),ifail)
*     ..region 2 : SOL : eastern edge
         call E01BEF(iyb-iya+1,yy(iya:iyb),fy(nx,iya:iyb),d,ifail)
         call E01BFF(iyb-iya+1,yy(iya:iyb),fy(nx,iya:iyb),d,
     .               iyd-iyc+1,yy1(iyc:iyd),fy1(nx1,iyc:iyd),ifail)
        endif
*     ..region 2 : SOL : central part : until second separatrix
        fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
        fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
        fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
        fy1(ixc:ixd,iyc:iyd) =
     .   fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
        if (nncut.eq.2 .and. topcut(1).ne.topcut(2)) then   ! DDN HFS SOL
          iya = max(topcut(1),topcut(2))
          iyb = ny
          iyc = max(topcut1(1),topcut1(2))
          iyd = ny1
          do iy = iya, iyb
            yy(iy) = rratio(iy-iya,iyb-iya)
          enddo
          do iy1 = iyc, iyd
            yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
          enddo
*     ..region 2 : DDN HFS SOL : above second separatrix
          fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
          fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
          fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
          call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .     xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .     fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
          fy1(ixc:ixd,iyc:iyd) =
     .     fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
        endif
       endif
*    ..region 3 :
       if (nncut.ge.1) then
         ixa = 0
         ixb = leftcut(1)-1
         ixm = 0
         ixp = rightcut(1)
         ixc = 0
         ixd = leftcut1(1)-1
         iya = 0
         iyb = jsep+1
         iyc = 0
         iyd = jsep1+1
         xx(ixa-1) = 0.0_R8
         do ix = ixa, ixb+1
          xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
         enddo
         do ix1 = ixc, ixd
          xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
         enddo
         do iy = iya, iyb
          yy(iy) = rratio(iy-iya,iyb-iya)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
         enddo
*     ..region 3 : western PFR : below the first separatrix (western edge)
         call E01BEF(iyb-iya+1,yy(iya:iyb),fy(ixa-1,iya:iyb),d,ifail)
         call E01BFF(iyb-iya+1,yy(iya:iyb),fy(ixa-1,iya:iyb),d,
     .               iyd-iyc+1,yy1(iyc:iyd),fy1(ixc-1,iyc:iyd),ifail)         
*     ..region 3 : western PFR : below the first separatrix
         fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
         fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
         fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
         fy1(ixc:ixd,iyc:iyd) =
     .    fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
*     ..region 3 : DDN HFS lower divertor : between separatrices
        if (nncut.eq.2 .and. topcut(1).ne.topcut(2)) then   ! DDN HFS SOL
         if (topcut(1).lt.topcut(2) .and. nnreg(0).ne.7) then
          ixp = leftcut(1)   ! position past the X-point
         endif
         iya = min(topcut(1),topcut(2))
         iyb = max(topcut(1),topcut(2))
         iyc = min(topcut1(1),topcut1(2))
         iyd = max(topcut1(1),topcut1(2))
         do iy = iya, iyb
          yy(iy) = rratio(iy-iya,iyb-iya)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
         enddo
*     ..region 3 : western divertor : between separatrices (western edge)
         call E01BEF(iyb-iya+1,yy(iya:iyb),fy(ixa-1,iya:iyb),d,ifail)
         call E01BFF(iyb-iya+1,yy(iya:iyb),fy(ixa-1,iya:iyb),d,
     .               iyd-iyc+1,yy1(iyc:iyd),fy1(ixc-1,iyc:iyd),ifail)
*     ..region 3 : DDN HFS lower divertor : between separatrices
         fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
         fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
         fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
         fy1(ixc:ixd,iyc:iyd) =
     .    fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
        endif
*     ..region 3 : western divertor SOL : above the separatrices
        ixp = leftcut(1)   ! position past the X-point
        if (nncut.eq.1 .and. periodic_bc.eq.0) then
         iya = topcut(1)
         iyb = ny
         iyc = topcut1(1)
         iyd = ny1
        else if (nncut.eq.2) then
         iya = max(topcut(1),topcut(2))
         iyb = ny
         iyc = max(topcut1(1),topcut1(2))
         iyd = ny1
        endif
        if (nncut.ge.1 .and. periodic_bc.eq.0) then
         do iy = iya, iyb
          yy(iy) = rratio(iy-iya,iyb-iya)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
         enddo
*     ..region 3 : western divertor SOL : above the separatrices (western edge)
         call E01BEF(iyb-iya+1,yy(iya:iyb),fy(ixa-1,iya:iyb),d,ifail)
         call E01BFF(iyb-iya+1,yy(iya:iyb),fy(ixa-1,iya:iyb),d,
     .               iyd-iyc+1,yy1(iyc:iyd),fy1(ixc-1,iyc:iyd),ifail)         
*     ..region 3 : western divertor SOL : above the separatrices
         fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
         fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
         fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
         fy1(ixc:ixd,iyc:iyd) =
     .    fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
        endif
       endif
*    ..region 4 :
       if (nncut.ge.1) then
        if (nncut.eq.1) then
         ixa = rightcut(1)
         ixb = nx-1
         ixm = leftcut(1)-1
         ixp = nx-1
         ixc = rightcut1(1)
         ixd = nx1-1
        else if (nnreg(0).eq.7) then
         ixa = rightcut(1)
         ixm = leftcut(1)-1
         ixb = leftcut(2)-1
         ixp = leftcut(2)-1
         ixc = rightcut1(1)
         ixd = leftcut1(2)-1
        else if (nncut.eq.2) then
         ixa = rightcut(2)
         ixm = leftcut(2)-1
         ixb = ixa
         do while (region(ixb+1,jsep,0).eq.4)
          ixb = ixb+1
         enddo
         ixb = ixb-1  ! step back not to count the guard cell
         ixp = ixb
         ixc = rightcut1(2)
         ixd = ixc
         do while (region1(ixd+1,jsep1,0).eq.4)
          ixd = ixd + 1
         enddo
         ixd = ixd - 1 ! step back not to count the guard cell
        endif
        do ix = ixa-1, ixb
          xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        xx(ixb+1) = 1.0_R8
        do ix1 = ixc, ixd
          xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
*     ..region 4 : eastern SN PFR or HFS top DN PFR : below the separatrix
        iya = 0
        iyc = 0
        if (nncut.eq.1) then
         iyb = topcut(1)
         iyd = topcut1(1)
        else if (nncut.eq.2) then
         iyb = min(topcut(1),topcut(2))
         iyd = min(topcut1(1),topcut1(2))
        endif
        do iy = iya, iyb
         yy(iy) = rratio(iy-iya,iyb-iya)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
        enddo
        fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
        fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
        fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
        fy1(ixc:ixd,iyc:iyd) =
     .   fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
*     ..region 4 : eastern SN PFR or HFS top DN PFR : below the separatrix (eastern edge)
        call E01BEF(iyb-iya+1,yy(iya:iyb),fy(ixb+1,iya:iyb),d,ifail)
        call E01BFF(iyb-iya+1,yy(iya:iyb),fy(ixb+1,iya:iyb),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fy1(ixd+1,iyc:iyd),ifail)
*     ..region 4 : HFS DDN top divertor : between the separatrices
        if (nncut.eq.2 .and. topcut(1).ne.topcut(2)) then   ! DDN HFS SOL
         iya = min(topcut(1),topcut(2))
         iyb = max(topcut(1),topcut(2))
         iyc = min(topcut1(1),topcut1(2))
         iyd = max(topcut1(1),topcut1(2))
         do iy = iya, iyb
          yy(iy) = rratio(iy-iya,iyb-iya)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
         enddo
         fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
         fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
         fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
         fy1(ixc:ixd,iyc:iyd) =
     .    fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
*     ..region 4 : HFS DDN top divertor : between the separatrices (eastern edge)
         call E01BEF(iyb-iya+1,yy(iya:iyb),fy(ixb+1,iya:iyb),d,ifail)
         call E01BFF(iyb-iya+1,yy(iya:iyb),fy(ixb+1,iya:iyb),d,
     .               iyd-iyc+1,yy1(iyc:iyd),fy1(ixd+1,iyc:iyd),ifail)
        endif
        if (nncut.ge.1 .and. periodic_bc.eq.0) then
*     ..region 4 : western divertor SN SOL or HFS DDN top divertor SOL : above the separatrices
         if (nncut.eq.1) then
          iya = topcut(1)
          iyc = topcut1(1)
         else if (nncut.eq.2) then
          iya = max(topcut(1),topcut(2))
          iyc = max(topcut1(1),topcut1(2))
         endif
         iyb = ny
         iyd = ny1
         do iy = iya, iyb
          yy(iy) = rratio(iy-iya,iyb-iya)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
         enddo
         fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
         fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
         fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
         fy1(ixc:ixd,iyc:iyd) =
     .    fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
*     ..region 4 : western divertor SN SOL or HFS DDN top divertor SOL : above the separatrices (eastern edge)
         call E01BEF(iyb-iya+1,yy(iya:iyb),fy(ixb+1,iya:iyb),d,ifail)
         call E01BFF(iyb-iya+1,yy(iya:iyb),fy(ixb+1,iya:iyb),d,
     .               iyd-iyc+1,yy1(iyc:iyd),fy1(ixd+1,iyc:iyd),ifail)
        endif        
       endif
*    ..region 5 :
       if (nncut.eq.1 .and. periodic_bc.eq.1) then  ! stellarator island
        ixa = 0
        ixm = nx-1
        iya = inseliy
        ixb = nx-1
        iyb = ny
        ixp = 0   ! for periodicity
        ixc = 0
        iyc = inseliy1
        ixd = nx1-1
        iyd = ny1
        fy1(nx1,inseliy1:ny1) = 0.0_R8   ! eastern dead cells
        do ix = ixa-1, ixb+1
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        do iy = iya, iyb
         yy(iy) = rratio(iy-iya,iyb-iya)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
        enddo
*     ..region 5 : stellarator island : central part
        fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
        fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
        fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
        fy1(ixc:ixd,iyc:iyd) =
     .   fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
       else if (nncut.eq.2) then  ! DN LFS core
*     ..region 5 : DN LFS core : central part
        if (nnreg(0).eq.7) then
          iya = 0
          iyb = ny
          iyc = 0
          iyd = ny1
          ixa = leftcut(2)
          ixb = ixa
          ixc = leftcut1(2)
          ixd = ixc
          ixm = leftcut(2)-1
          ixp = leftcut(1)
          do while (region(ixb+1,jsep,0).eq.5)
           ixb = ixb + 1
          enddo
          ixb = ixb - 1   ! step back to avoid counting guard cells
          ixp = ixb
          ixd = ixc
          do while (region1(ixd+1,jsep1,0).eq.5)
           ixd = ixd + 1
          enddo
          ixd = ixd - 1   ! step back to avoid counting guard cells
        else
          ixa = leftcut(2)
          ixm = rightcut(2)-1  ! for continuity
          ixb = rightcut(1)-1
          ixp = leftcut(1)     ! for continuity
          ixc = leftcut1(2)
          ixd = rightcut1(1)-1
          iya = 0
          iyb = min(topcut(1),topcut(2))
          iyc = 0
          iyd = min(topcut1(1),topcut1(2))
        endif
        do ix = ixa-1, ixb+1
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        do iy = iya, iyb
         yy(iy) = rratio(iy-iya,iyb-iya)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
        enddo
        fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
        fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
        fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
        fy1(ixc:ixd,iyc:iyd) =
     .   fy1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
       endif
*    ..region 6 :
       if (nncut.eq.2) then
        if (nnreg(0).eq.7) then
          ixb = rightcut(2)-1
          ixd = rightcut1(2)-1
          do while (region(ixa-1,jsep,0).eq.6)
            ixa = ixa-1
          enddo
          ixa = ixa + 1  ! step back to not count guard cell
          ixm = ixa
          ixc = ixd
          do while (region1(ixc-1,jsep1,0).eq.6)
           ixc = ixc-1
          enddo
          ixc = ixc + 1  ! step back to not count guard cell
          iya = 0
          iyc = 0
        else
          ixa = rightcut(2)
          iya = min(topcut(1),topcut(2))
          ixb = rightcut(1)-1
          ixc = rightcut1(2)
          iya = min(topcut1(1),topcut1(2))
          ixd = rightcut1(1)-1
        endif
        if (topcut(1).eq.topcut(2)) then   ! CDN Balanced X-points
         ixp = rightcut(1)   ! position past the X-point
         ixm = rightcut(2)-1 ! position past the X-point
         iyb = ny
         iyd = ny1
        else if (nnreg(0).eq.7) then
         ixm = leftcut(2)-1
         ixp = rightcut(1)
         iyb = max(topcut(1),topcut(2))
         iyd = max(topcut1(1),topcut1(2))
        else if (topcut(1).lt.topcut(2)) then   ! Bottom X-point is active or CDN
         ixp = rightcut(1)   ! position past the X-point
         ixm = leftcut(2)-1  ! continuity to HFS
         iyb = max(topcut(1),topcut(2))
         iyd = max(topcut1(1),topcut1(2))
        else if (topcut(1).gt.topcut(2)) then  ! Top X-point is active
         ixp = leftcut(1)    ! continuity to HFS
         ixm = rightcut(2)-1 ! position past the X-point
         iyb = max(topcut(1),topcut(2))
         iyd = max(topcut1(1),topcut1(2))
        endif
        do ix = ixa-1, ixb+1
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        do iy = iya, iyb
         yy(iy) = rratio(iy-iya,iyb-iya)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
        enddo
*     ..region 6 : CDN LFS SOL or DDN LFS SOL between the separatrices
        fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
        fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
        fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
        fy1(ixc:ixd,iyc:iyd) =
     .   fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
        if (topcut(1).ne.topcut(2)) then   ! DDN HFS SOL
*     ..region 6 : DDN LFS SOL above the second separatrix
         if (nnreg(0).ne.7) then
          ixp = rightcut(1)   ! position past the X-point
          ixm = rightcut(2)-1 ! position past the X-point
         endif
         iya = max(topcut(1),topcut(2))
         iyb = ny
         iyc = max(topcut1(1),topcut1(2))
         iyd = ny1
         do iy = iya, iyb
          yy(iy) = rratio(iy-iya,iyb-iya)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
         enddo
         fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
         fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
         fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
         fy1(ixc:ixd,iyc:iyd) =
     .    fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
        endif
       endif
*    ..region 7 :
       if (nncut.eq.2) then
        if (nnreg(0).eq.7) then
          ixa = rightcut(2)
          ixm = leftcut(1)-1
          ixc = rightcut1(2)
          ixb = nx-1
          ixp = nx-1
          ixd = nx1-1
        else
          ixb = rightcut(2)-1
          ixa = ixb
          do while (region(ixa-1,jsep,0).eq.7)
            ixa = ixa-1
          enddo
          ixa = ixa + 1                    ! step back to avoid western edge
          ixp = leftcut(2)                 ! position on HFS
          ixd = rightcut1(2)-1
          ixc = ixd
          do while (region1(ixc-1,jsep1,0).eq.7)
            ixc = ixc-1
          enddo
          ixc = ixc + 1   ! step back to avoid western edge
        endif
        xx(ixa-1) = 0.0_R8
        do ix = ixa, ixb+1
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
*     ..region 7 : top LFS PFR : below the first separatrix
        iya = 0
        iyb = min(topcut(1),topcut(2))+1
        iyc = 0
        iyd = min(topcut1(1),topcut1(2))+1
        do iy = iya, iyb
         yy(iy) = rratio(iy-iya,iyb-iya)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
        enddo
*     ..region 7 : top LFS PFR : below the first separatrix (western edge)
        call E01BEF(iyb-iya+1,yy(iya:iyb),fy(ixa-1,iya:iyb),d,ifail)
        call E01BFF(iyb-iya+1,yy(iya:iyb),fy(ixa-1,iya:iyb),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fy1(ixc-1,iyc:iyd),ifail)
*     ..region 7 : top LFS PFR : below the first separatrix
        fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
        fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
        fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
        fy1(ixc:ixd,iyc:iyd) =
     .   fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
*     ..region 7 : DDN LFS upper divertor : between separatrices
        if (nncut.eq.2 .and. topcut(1).ne.topcut(2)) then   ! DDN HFS SOL
         if (topcut(1).gt.topcut(2) .and. nnreg(0).ne.7) then
          ixp = rightcut(2)   ! position past the X-point
         endif
         iya = min(topcut(1),topcut(2))
         iyb = max(topcut(1),topcut(2))
         iyc = min(topcut1(1),topcut1(2))
         iyd = max(topcut1(1),topcut1(2))
         do iy = iya, iyb
          yy(iy) = rratio(iy-iya,iyb-iya)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
         enddo
*     ..region 7 : DDN LFS upper divertor : between separatrices (western edge)
         call E01BEF(iyb-iya+1,yy(iya:iyb),fy(ixa-1,iya:iyb),d,ifail)
         call E01BFF(iyb-iya+1,yy(iya:iyb),fy(ixa-1,iya:iyb),d,
     .               iyd-iyc+1,yy1(iyc:iyd),fy1(ixc-1,iyc:iyd),ifail)
*     ..region 7 : DDN LFS upper divertor : between separatrices
         fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
         fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
         fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
         fy1(ixc:ixd,iyc:iyd) =
     .    fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
        endif
*     ..region 7 : DDN LFS upper divertor : above the second separatrix
        ixp = rightcut(2)   ! position past the X-point
        iya = max(topcut(1),topcut(2))
        iyb = ny
        iyc = max(topcut1(1),topcut1(2))
        iyd = ny1
        do iy = iya, iyb
         yy(iy) = rratio(iy-iya,iyb-iya)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
        enddo
*     ..region 7 : DDN LFS upper divertor : above the second separatrix (western edge)
        call E01BEF(iyb-iya+1,yy(iya:iyb),fy(ixa-1,iya:iyb),d,ifail)
        call E01BFF(iyb-iya+1,yy(iya:iyb),fy(ixa-1,iya:iyb),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fy1(ixc-1,iyc:iyd),ifail)
*     ..region 7 : DDN LFS upper divertor : above the second separatrix
        fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
        fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
        fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
        fy1(ixc:ixd,iyc:iyd) =
     .   fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
       endif
*    ..region 8 :
       if (nncut.eq.2 .and. nnreg(0).ne.7) then
        ixa = rightcut(1)
        ixm = leftcut(1)-1
        ixb = nx-1
        ixp = nx-1
        ixc = rightcut1(1)
        ixd = nx1
        do ix = ixa-1, ixb
          xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        xx(ixb+1) = 1.0_R8
        do ix1 = ixc, ixd
          xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
*     ..region 8 : eastern DN LFS PFR : below the separatrix (including eastern edge)
        iya = 0
        iyb = min(topcut(1),topcut(2))
        iyc = 0
        iyb = min(topcut1(1),topcut1(2))
        do iy = iya, iyb
         yy(iy) = rratio(iy-iya,iyb-iya)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
        enddo
        fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
        fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
        fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
        fy1(ixc:ixd,iyc:iyd) =
     .   fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
*     ..region 8 : eastern DN LFS PFR : below the separatrix (eastern edge)
        call E01BEF(iyb-iya+1,yy(iya:iyb),fy(ixb+1,iya:iyb),d,ifail)
        call E01BFF(iyb-iya+1,yy(iya:iyb),fy(ixb+1,iya:iyb),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fy1(ixd+1,iyc:iyd),ifail)
*     ..region 8 : eastern DN LFS divertor : between the separatrices
        if (topcut(1).ne.topcut(2)) then
         if (topcut(1).gt.topcut(1)) ixm = rightcut(1)-1
         iya = min(topcut(1),topcut(2))
         iyb = max(topcut(1),topcut(2))
         iyc = min(topcut1(1),topcut1(2))
         iyd = max(topcut1(1),topcut1(2))
         do iy = iya, iyb
          yy(iy) = rratio(iy-iya,iyb-iya)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
         enddo
         fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
         fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
         fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .    fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
         fy1(ixc:ixd,iyc:iyd) =
     .    fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
*     ..region 8 : eastern DN LFS divertor : between the separatrices (eastern edge)
         call E01BEF(iyb-iya+1,yy(iya:iyb),fy(ixb+1,iya:iyb),d,ifail)
         call E01BFF(iyb-iya+1,yy(iya:iyb),fy(ixb+1,iya:iyb),d,
     .               iyd-iyc+1,yy1(iyc:iyd),fy1(ixd+1,iyc:iyd),ifail)
        endif
*     ..region 8 : eastern DN LFS divertor : above the second separatrix
        ixm = rightcut(1)-1
        iya = max(topcut(1),topcut(2))
        iyb = ny
        iyc = max(topcut1(1),topcut1(2))
        iyd = ny1
        do iy = iya, iyb
         yy(iy) = rratio(iy-iya,iyb-iya)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = rratio(iy1-iyc,iyd-iyc)
        enddo
        fa(ixa-1,iya:iyb) = fy(ixm,iya:iyb)/hx(ixm,iya:iyb)
        fa(ixa:ixb,iya:iyb) = fy(ixa:ixb,iya:iyb)/hx(ixa:ixb,iya:iyb)
        fa(ixb+1,iya:iyb) = fy(ixp,iya:iyb)/hx(ixp,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+1,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya:iyb),xx1(ixc:ixd),yy1(iyc:iyd),
     .   fa(ixa-1:ixb+1,iya:iyb),fa1(ixc:ixd,iyc:iyd))
        fy1(ixc:ixd,iyc:iyd) =
     .   fa1(ixc:ixd,iyc:iyd)*hx1(ixc:ixd,iyc:iyd)
*     ..region 8 : eastern DN LFS divertor : above the second separatrix (eastern edge)
        call E01BEF(iyb-iya+1,yy(iya:iyb),fy(ixb+1,iya:iyb),d,ifail)
        call E01BFF(iyb-iya+1,yy(iya:iyb),fy(ixb+1,iya:iyb),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fy1(ixd+1,iyc:iyd),ifail)
       endif
*    ..clean-up
       deallocate(d)
      endif
*   ..test edge values of fy1
      call b2xvfy1 (nx1, ny1, fy1, 'b2xtfy--fy1', bottomiy1)

      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xtfy
