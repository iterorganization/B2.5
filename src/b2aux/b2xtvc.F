      subroutine b2xtvc (nx, ny, nx1, ny1, fun, fun1, 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk1)
      use b2mod_types
      use b2mod_math
      use b2mod_subsys
      implicit none
      integer nx, ny, nx1, ny1, nncutmax
      real (kind=R8) ::
     *   fun(-1:nx,-1:ny), fun1(-1:nx1,-1:ny1), wrk1(-1:nx,-1:ny1)
      integer non_commensurate
      integer nnreg(0:2), nncut, periodic_bc, jsep, jsep1
      integer region(-1:nx,-1:ny,0:2), region1(-1:nx1,-1:ny1,0:2)
      integer leftcut(nncutmax), leftcut1(nncutmax)
      integer rightcut(nncutmax), rightcut1(nncutmax)
      integer bottomcut(nncutmax), bottomcut1(nncutmax)
      integer topcut(nncutmax), topcut1(nncutmax)
      logical isSnowflake
*     ------------------------------------------------------------------
*     B2XTVC transfers a cell center variable to a new mesh.
*     ------------------------------------------------------------------
      integer ix1, iy1, ix, iy, nd, ifail
      integer ixa, ixb, iya, iyb, ixc, ixd, iyc, iyd
      integer ixm, ixp, iym, iyp
      integer inseliy, inseliy1
      real (kind=R8), allocatable :: d(:)
      real (kind=R8) :: xx(-1:nx), xx1(-1:nx1), yy(-1:ny), yy1(-1:ny1)
      external E01BEF, E01BFF
      real (kind=R8) ::
     *   t0, t1, t2
      external xerrab, interp2d
      intrinsic min, max, real
*     ------------------------------------------------------------------
      call subini ('b2xtvc')
      isSnowflake = .False.
      if (nncut.eq.2) then
         if (rightcut(1).le.rightcut(2)) then
            isSnowflake = .True.
         endif
      endif
      if (non_commensurate.eq.0) then   ! commensurate mapping
*   ..interpolate in the second index
       do ix = -1, nx
        wrk1(ix,-1) = fun(ix,-1)
       enddo
       do iy1 = 0, ny1-1
        iy = (2*iy1*ny+ny)/(2*ny1)
*      (iy specifies the nearest point)
        if ((2*iy+1)*ny1.eq.(2*iy1+1)*ny) then
         do ix = -1, nx
          wrk1(ix,iy1) = fun(ix,iy)
         enddo
        else if (2*iy*ny1.eq.(2*iy1+1)*ny) then
         do ix = -1, nx
          wrk1(ix,iy1) = (fun(ix,iy-1)+fun(ix,iy))/2.0_R8
         enddo
        else
         t0 = real((2*iy+1)*ny1-(2*iy1+1)*ny)/(2*ny1)
         do ix = -1, nx
          t1 = t0*fun(ix,iy-1)+(1.0_R8-t0)*fun(ix,iy)
          t2 = (1.0_R8+t0)*fun(ix,iy)-t0*fun(ix,iy+1)
          wrk1(ix,iy1) = med(t1,t2,fun(ix,iy))
         enddo
        endif
       enddo
       do ix = -1, nx
        wrk1(ix,ny1) = fun(ix,ny)
       enddo
*   ..interpolate in the first index
       do iy1 = -1, ny1
        fun1(-1,iy1) = wrk1(-1,iy1)
       enddo
       do ix1 = 0, nx1-1
        ix = (2*ix1*nx+nx)/(2*nx1)
*      (ix specifies the nearest point)
        if ((2*ix+1)*nx1.eq.(2*ix1+1)*nx) then
         do iy1 = -1, ny1
          fun1(ix1,iy1) = wrk1(ix,iy1)
         enddo
        else if (2*ix*nx1.eq.(2*ix1+1)*nx) then
         do iy1 = -1, ny1
          fun1(ix1,iy1) = (wrk1(ix-1,iy1)+wrk1(ix,iy1))/2.0_R8
         enddo
        else
         t0 = real((2*ix+1)*nx1-(2*ix1+1)*nx)/(2*nx1)
         do iy1 = -1, ny1
          t1 = t0*wrk1(ix-1,iy1)+(1.0_R8-t0)*wrk1(ix,iy1)
          t2 = (1.0_R8+t0)*wrk1(ix,iy1)-t0*wrk1(ix+1,iy1)
          fun1(ix1,iy1) = med(t1,t2,wrk1(ix,iy1))
         enddo
        endif
       enddo
       do iy1 = -1, ny1
        fun1(nx1,iy1) = wrk1(nx,iy1)
       enddo
      else  ! non-commensurate conformal mapping
       nd=max(nx,nx1,ny,ny1)
       allocate(d(nd+2))
       ifail=0
*   ..we proceed per region
*    .. all topologies : corners
       fun1( -1, -1) = fun(-1,-1)
       fun1( -1,ny1) = fun(-1,ny) 
       fun1(nx1, -1) = fun(nx,-1)
       fun1(nx1,ny1) = fun(nx,ny) 
*    .. region 1 : 
       if (nncut.eq.0 .and. periodic_bc.le.0) then ! slab cases
        ixa = 0
        ixm = 0
        ixb = nx-1
        ixp = nx-1
        ixc = 0
        ixd = nx1-1
        iya = 0
        iym = 0
        iyb = ny-1
        iyp = ny-1
        iyc = 0
        iyd = ny1-1
        xx(ixa-1) = 0.0_R8
        do ix = ixa, ixb
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        xx(ixb+1) = 1.0_R8
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        yy(iya-1) = 0.0_R8
        do iy = iya, iyb
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
*    .. region 1 : western edge
        wrk1(-1,iya-1) = fun(-1,iym)
        wrk1(-1,iya:iyb) = fun(-1,iya:iyb)
        wrk1(-1,iyb+1) = fun(-1,iyp)
        call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(-1,iya-1:iyb+1),d,
     .              ifail)
        call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(-1,iya-1:iyb+1),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fun1(-1,iyc:iyd),ifail)
*    .. region 1 : eastern edge
        wrk1(nx,iya-1) = fun(nx,iym)
        wrk1(nx,iya:iyb) = fun(nx,iya:iyb)
        wrk1(nx,iyb+1) = fun(nx,iyp)
        call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(nx,iya-1:iyb+1),d,
     .              ifail)
        call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(nx,iya-1:iyb+1),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fun1(nx1,iyc:iyd),ifail)
*    .. region 1 : top row
        wrk1(ixa-1,ny) = fun(ixm,ny)
        wrk1(ixa:ixb,ny) = fun(ixa:ixb,ny)
        wrk1(ixb+1,ny) = fun(ixp,ny)
        call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .              ifail)
        call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,ny1),ifail)
       else   ! core cases
        iya = 0
        iyb = jsep
        iyp = jsep+1
        iym = 0
        iyc = 0
        iyd = jsep1
        if (nncut.eq.0 .and. periodic_bc.eq.1) then    ! limiter case
          ixa = 0
          ixb = nx-1
          ixm = nx-1
          ixp = 0
          ixc = 0
          ixd = nx1-1
          fun1(-1,iyc:iyd) = fun(-1,(iya+iyb)/2)    ! western dead cells
          fun1(nx1,iyc:iyd) = fun(nx1,(iya+iyb)/2)  ! eastern dead cells
        else if (nncut.eq.1 .or. isSnowflake) then  ! SN and SF topologies
          ixa = leftcut(1)
          ixb = rightcut(1)-1
          ixp = leftcut(1)
          ixm = rightcut(1)-1
          ixc = leftcut1(1)
          ixd = rightcut1(1)-1
        else if (nncut.eq.2) then    ! DN topologies
          ixa = leftcut(1)
          ixb = leftcut(2)-1
          ixp = rightcut(2)
          ixm = rightcut(1)-1
          ixc = leftcut1(1)
          ixd = leftcut1(2)-1
        endif
        do ix = ixa-1, ixb+1
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        yy(iya-1) = 0.0_R8
        do iy = iya, iyb+1
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
       endif
*    .. region 1 : bottom row
       wrk1(ixa-1,-1) = fun(ixm,-1)
       wrk1(ixa:ixb,-1) = fun(ixa:ixb,-1)
       wrk1(ixb+1,-1) = fun(ixp,-1)
       call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,-1),d,
     .             ifail)
       call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,-1),d,
     .             ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,-1),ifail)
*    .. region 1 : central part
       wrk1(ixa-1,iya-1) = fun(ixm,iym)
       wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
       wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
       wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
       wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
       wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
       wrk1(ixb+1,iya-1) = fun(ixp,iym)
       wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
       wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
       call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .  xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .  wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
*    .. region 2 :
       if (periodic_bc.eq.1) then   ! limiter or stellarator island case
        ixa = 0
        ixm = 0
        ixb = nx-1
        ixp = nx-1
        ixc = 0
        ixd = nx1-1
        xx(ixa-1) = 0.0_R8
        do ix = ixa, ixb
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        xx(ixb+1) = 1.0_R8
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        iya = jsep+1
        iym = jsep
        iyc = jsep1+1
        if (nncut.eq.0) then  ! limiter case
         iyb = ny-1
         iyp = ny-1
         iyd = ny1-1
        else if (nncut.eq.1) then   ! stellarator island case
         inseliy = iya
         do while (region(inseliy+1,nx/2,0).eq.2)
          inseliy = inseliy+1
         enddo
         iyb = inseliy
         iyp = inseliy+1
         inseliy1 = iyc
         do while (region1(inseliy1+1,nx1/2,0).eq.2)
          inseliy1 = inseliy1+1
         enddo
         iyd = inseliy1
        endif
        do iy = iya-1, iyb+1
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
*    .. region 2 : periodic_bc.eq.1 : western edge
        wrk1(-1,iya-1) = fun(-1,iym)
        wrk1(-1,iya:iyb) = fun(-1,iya:iyb)
        wrk1(-1,iyb+1) = fun(-1,iyp)
        call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(-1,iya-1:iyb+1),d,
     .              ifail)
        call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(-1,iya-1:iyb+1),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fun1(-1,iyc:iyd),ifail)
*    .. region 2 : periodic_bc.eq.1 : eastern edge
        wrk1(nx,iya-1) = fun(nx,iym)
        wrk1(nx,iya:iyb) = fun(nx,iya:iyb)
        wrk1(nx,iyb+1) = fun(nx,iyp)
        call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(nx,iya-1:iyb+1),d,
     .              ifail)
        call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(nx,iya-1:iyb+1),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fun1(nx1,iyc:iyd),ifail)
       else
        iya = jsep+1
        iym = jsep
        iyc = jsep1+1
        ixa = leftcut(1)
        ixc = leftcut1(1)
        if (nncut.eq.1 .or. isSnowflake) then
         ixm = leftcut(1)-1
         ixb = rightcut(1)-1
         ixp = rightcut(1)
         ixd = rightcut1(1)-1
        else if (nncut.eq.2 .and. topcut(1).eq.topcut(2)) then
         ixm = leftcut(1)-1
         ixb = leftcut(2)-1
         ixp = leftcut(2)
         ixd = leftcut1(2)-1
        else if (nncut.eq.2 .and. topcut(1).ne.topcut(2)) then
         ixb = leftcut(2)-1
         ixd = leftcut1(2)-1
         if (topcut(1).lt.topcut(2)) then
          ixm = leftcut(1)-1
          ixp = rightcut(2)
         else if (topcut(1).gt.topcut(2)) then
          ixm = rightcut(1)-1
          ixp = leftcut(2)
         endif
         iyb = max(topcut(1),topcut(2))
         iyp = max(topcut(1),topcut(2))+1
         iyd = max(topcut1(1),topcut1(2))
         do ix = ixa-1, ixb+1
          xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
         enddo
         do ix1 = ixc, ixd
          xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
         enddo
         do iy = iya-1, iyb+1
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
*    .. region 2 : DN HFS SOL : between the separatrices
         wrk1(ixa-1,iya-1) = fun(ixm,iym)
         wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
         wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
         wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
         wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
         wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
         wrk1(ixb+1,iya-1) = fun(ixp,iym)
         wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
         wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
         ixm = leftcut(1)-1
         ixp = leftcut(2)
         iya = max(topcut(1),topcut(2))+1
         iym = max(topcut(1),topcut(2))
         iyc = max(topcut1(1),topcut1(2))+1
        endif
*    .. region 2 : SOL : above the separatrices
        iyb = ny-1
        iyp = ny-1
        iyd = ny1-1
        do ix = ixa-1, ixb+1
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        do iy = iya-1, iyb
         yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
        wrk1(ixa-1,iya-1) = fun(ixm,iym)
        wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
        wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
        wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
        wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
        wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
        wrk1(ixb+1,iya-1) = fun(ixp,iym)
        wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
        wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
        call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
        if (nnreg(0).ne.5) then
*    .. region 2 : SOL : top row
         wrk1(ixa-1,ny) = fun(ixm,ny)
         wrk1(ixa:ixb,ny) = fun(ixa:ixb,ny)
         wrk1(ixb+1,ny) = fun(ixp,ny)
         call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .               ifail)
         call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .               ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,ny1),ifail)
        endif
       endif
*    .. region 3 :
       if (nncut.ge.1) then
        ixa = 0
        ixm = 0
        ixb = leftcut(1)-1
        ixp = rightcut(1)
        ixc = 0
        ixd = leftcut1(1)-1
        iya = 0
        iyb = jsep
        iym = 0
        iyp = jsep+1
        iyc = 0
        iyd = jsep1
        xx(ixa-1) = 0.0_R8
        do ix = ixa, ixb+1
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        yy(iya) = 0.0_R8
        do iy = iya, iyb+1
         yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
*    .. region 3 : western PFR : bottom row
        wrk1(ixa-1,-1) = fun(ixm,-1)
        wrk1(ixa:ixb,-1) = fun(ixa:ixb,-1)
        wrk1(ixb+1,-1) = fun(ixp,-1)
        call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,-1),d,
     .              ifail)
        call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,-1),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,-1),ifail)
*    .. region 3 : western PFR : western edge
        wrk1(-1,iya-1) = fun(-1,iym)
        wrk1(-1,iya:iyb) = fun(-1,iya:iyb)
        wrk1(-1,iyb+1) = fun(-1,iyp)
        call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(-1,iya-1:iyb+1),d,
     .              ifail)
        call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(-1,iya-1:iyb+1),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fun1(-1,iyc:iyd),ifail)
*    .. region 3 : western PFR : central part
        wrk1(ixa-1,iya-1) = fun(ixm,iym)
        wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
        wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
        wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
        wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
        wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
        wrk1(ixb+1,iya-1) = fun(ixp,iym)
        wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
        wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
        if (nncut.eq.2 .and. topcut(1).ne.topcut(2)) then   ! DDN HFS SOL
         if (topcut(1).lt.topcut(2) .and. .not.isSnowflake) then
          ixp = leftcut(1)
         endif
         iya = min(topcut(1),topcut(2))
         iym = min(topcut(1),topcut(2))-1
         iyb = max(topcut(1),topcut(2))
         iyp = max(topcut(1),topcut(2))+1
         iyc = min(topcut1(1),topcut1(2))
         iyd = max(topcut1(1),topcut1(2))
         do iy = iya-1, iyb+1
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
*    .. region 3 : western divertor : between the separatrices : western edge
         wrk1(-1,iya-1) = fun(-1,iym)
         wrk1(-1,iya:iyb) = fun(-1,iya:iyb)
         wrk1(-1,iyb+1) = fun(-1,iyp)
         call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(-1,iya-1:iyb+1),d,
     .               ifail)
         call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(-1,iya-1:iyb+1),d,
     .               iyd-iyc+1,yy1(iyc:iyd),fun1(-1,iyc:iyd),ifail)
*    .. region 3 : western divertor : between the separatrices : central part
         wrk1(ixa-1,iya-1) = fun(ixm,iym)
         wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
         wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
         wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
         wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
         wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
         wrk1(ixb+1,iya-1) = fun(ixp,iym)
         wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
         wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
        endif
        if (periodic_bc.eq.0) then
         ixp = leftcut(1)
         iyb = ny-1
         iyp = ny-1
         iyd = ny1-1
         if (nncut.eq.1 .or. topcut(1).eq.topcut(2)) then
          iya = jsep+1
          iym = jsep
          iyc = jsep1
         else
          iya = max(topcut(1),topcut(2))+1
          iym = max(topcut(1),topcut(2))
          iyc = max(topcut1(1),topcut1(2))+1
         endif
         do iy = iya-1, iyb
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         yy(iyb+1) = 1.0_R8
         do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
*    .. region 3 : western divertor : above the separatrices : western edge
         wrk1(-1,iya-1) = fun(-1,iym)
         wrk1(-1,iya:iyb) = fun(-1,iya:iyb)
         wrk1(-1,iyb+1) = fun(-1,iyp)
         call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(-1,iya-1:iyb+1),d,
     .               ifail)
         call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(-1,iya-1:iyb+1),d,
     .               iyd-iyc+1,yy1(iyc:iyd),fun1(-1,iyc:iyd),ifail)
*    .. region 3 : western divertor : above the separatrices : central part
         wrk1(ixa-1,iya-1) = fun(ixm,iym)
         wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
         wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
         wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
         wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
         wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
         wrk1(ixb+1,iya-1) = fun(ixp,iym)
         wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
         wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
*    .. region 3 : western divertor : above the separatrices : top row
         wrk1(ixa-1,ny) = fun(ixm,ny)
         wrk1(ixa:ixb,ny) = fun(ixa:ixb,ny)
         wrk1(ixb+1,ny) = fun(ixp,ny)
         call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .               ifail)
         call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .               ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,ny1),ifail)
        endif
       endif
*    .. region 4 :
       if (nncut.ge.1) then
        iya = 0
        iyb = jsep
        iym = 0
        iyp = jsep+1
        iyc = 0
        iyd = jsep1
        if (nncut.eq.1) then
         ixa = rightcut(1)
         ixm = leftcut(1)-1
         ixb = nx-1
         ixp = nx-1
         ixc = rightcut1(1)
         ixd = nx1-1
        else if (isSnowflake) then
         ixa = rightcut(1)
         ixm = leftcut(1)-1
         ixb = leftcut(2)-1
         ixp = leftcut(2)-1
         ixc = rightcut1(1)
         ixd = leftcut1(2)-1
        else if (nncut.eq.2) then
         ixa = leftcut(2)
         ixm = rightcut(2)-1
         ixc = leftcut1(2)
         ixb = ixa
         do while (region(ixb+1,jsep,0).eq.4)
          ixb = ixb + 1
         enddo
         ixb = ixb - 1   ! step back to avoid counting guard cells
         ixp = ixb
         ixd = ixc
         do while (region1(ixd+1,jsep1,0).eq.4)
          ixd = ixd + 1
         enddo
         ixd = ixd - 1   ! step back to avoid counting guard cells
        endif
        do ix = ixa-1, ixb
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        xx(ixb+1) = 1.0_R8
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        yy(iya) = 0.0_R8
        do iy = iya, iyb+1
         yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
*    .. region 4 : eastern SN PFR or DN HFS top PFR : bottom row
        wrk1(ixa-1,-1) = fun(ixm,-1)
        wrk1(ixa:ixb,-1) = fun(ixa:ixb,-1)
        wrk1(ixb+1,-1) = fun(ixp,-1)
        call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,-1),d,
     .              ifail)
        call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,-1),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,-1),ifail)
*    .. region 4 : eastern SN PFR or DN HFS top PFR : eastern edge
        wrk1(ixb+1,iya-1) = fun(ixb+1,iym)
        wrk1(ixb+1,iya:iyb) = fun(ixb+1,iya:iyb)
        wrk1(ixb+1,iyb+1) = fun(ixb+1,iyp)
        call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),
     .              wrk1(ixb+1,iya-1:iyb+1),d,ifail)
        call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(ixb+1,iya-1:iyb+1),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fun1(ixd+1,iyc:iyd),ifail)
*    .. region 4 : eastern SN PFR or DN HFS top PFR : central part
        wrk1(ixa-1,iya-1) = fun(ixm,iym)
        wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
        wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
        wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
        wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
        wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
        wrk1(ixb+1,iya-1) = fun(ixp,iym)
        wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
        wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
        if (nncut.eq.2 .and. topcut(1).ne.topcut(2)) then   ! DDN HFS SOL
         if (topcut(1).gt.topcut(2) .and. .not.isSnowflake) then
          ixm = leftcut(2)-1
         endif
         iya = min(topcut(1),topcut(2))
         iym = min(topcut(1),topcut(2))-1
         iyb = max(topcut(1),topcut(2))
         iyp = max(topcut(1),topcut(2))+1
         iyc = min(topcut1(1),topcut1(2))
         iyd = max(topcut1(1),topcut1(2))
         do iy = iya-1, iyb+1
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
*    .. region 4 : eastern SN or DN HFS top divertor : between the separatrices : eastern edge
         wrk1(ixb+1,iya-1) = fun(ixb+1,iym)
         wrk1(ixb+1,iya:iyb) = fun(ixb+1,iya:iyb)
         wrk1(ixb+1,iyb+1) = fun(ixb+1,iyp)
         call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(ixb+1,iya-1:iyb+1),
     .               d,ifail)
         call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(ixb+1,iya-1:iyb+1),
     .               d,iyd-iyc+1,yy1(iyc:iyd),fun1(ixd+1,iyc:iyd),ifail)
*    .. region 4 : eastern SN PFR or DN HFS top PFR : between the separatrices : central part
         wrk1(ixa-1,iya-1) = fun(ixm,iym)
         wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
         wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
         wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
         wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
         wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
         wrk1(ixb+1,iya-1) = fun(ixp,iym)
         wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
         wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
        endif
        if (periodic_bc.eq.0) then
         ixm = ixa-1
         iyb = ny-1
         iyp = ny-1
         iyd = ny1-1
         if (nncut.eq.1 .or. topcut(1).eq.topcut(2)) then
          iya = jsep+1
          iym = jsep
          iyc = jsep1
         else
          iya = max(topcut(1),topcut(2))+1
          iym = max(topcut(1),topcut(2))
          iyc = max(topcut1(1),topcut1(2))+1
         endif
         do iy = iya-1, iyb
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         yy(iyb+1) = 1.0_R8
         do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
*    .. region 4 : eastern SN or DN HFS top divertor : above the separatrices : eastern edge
         wrk1(ixb+1,iya-1) = fun(ixb+1,iym)
         wrk1(ixb+1,iya:iyb) = fun(ixb+1,iya:iyb)
         wrk1(ixb+1,iyb+1) = fun(ixb+1,iyp)
         call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),
     .               wrk1(ixb+1,iya-1:iyb+1),d,ifail)
         call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),
     .               wrk1(ixb+1,iya-1:iyb+1),
     .               d,iyd-iyc+1,yy1(iyc:iyd),
     .               fun1(ixd+1,iyc:iyd),ifail)
*    .. region 4 : eastern SN or DN HFS top divertor : above the separatrices : central part
         wrk1(ixa-1,iya-1) = fun(ixm,iym)
         wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
         wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
         wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
         wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
         wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
         wrk1(ixb+1,iya-1) = fun(ixp,iym)
         wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
         wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
*    .. region 4 : eastern SN or DN HFS top divertor : above the separatrices : top row
         wrk1(ixa-1,ny) = fun(ixm,ny)
         wrk1(ixa:ixb,ny) = fun(ixa:ixb,ny)
         wrk1(ixb+1,ny) = fun(ixp,ny)
         call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .               ifail)
         call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .               ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,ny1),ifail)
        endif
       endif
*    .. region 5 :
       if (nncut.eq.1 .and. periodic_bc.eq.1) then
        ixa = 0
        ixb = nx-1
        ixm = nx-1
        ixp = 0
        ixc = 0
        ixd = nx1-1
        iya = inseliy+1
        iyb = ny-1
        iyc = inseliy1+1
        iyd = ny1-1
        iyp = ny-1
        iym = inseliy
        do ix = ixa-1, ixb+1
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        do iy = iya-1, iyb
         yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
*    .. region 5 : stellarator island : western and eastern edges : dead cells
        fun1(-1,iyc:iyd) = fun(-1,(iya+iyb)/2)
        fun1(nx1,iyc:iyd) = fun(nx,(iya+iyb)/2)
*    .. region 5 : stellarator island : top row
        wrk1(ixa-1,ny) = fun(ixm,ny)
        wrk1(ixa:ixb,ny) = fun(ixa:ixb,ny)
        wrk1(ixb+1,ny) = fun(ixp,ny)
        call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .              ifail)
        call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,ny1),ifail)
*    .. region 5 : stellarator island : central part
        wrk1(ixa-1,iya-1) = fun(ixm,iym)
        wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
        wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
        wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
        wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
        wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
        wrk1(ixb+1,iya-1) = fun(ixp,iym)
        wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
        wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
       else if (nncut.eq.2) then
        ixa = rightcut(2)
        ixb = rightcut(1)-1
        ixc = rightcut1(2)
        ixd = rightcut1(1)-1
        ixm = leftcut(2)-1
        ixp = leftcut(1)
        iya = 0
        iyb = jsep
        iyc = 0
        iyd = jsep1
        iym = 0
        iyp = jsep+1
        do ix = ixa-1, ixb+1
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        yy(iya-1) = 0.0_R8
        do iy = iya, iyb+1
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
*    .. region 5 : DN LFS core : bottom row
        wrk1(ixa-1,-1) = fun(ixm,-1)
        wrk1(ixa:ixb,-1) = fun(ixa:ixb,-1)
        wrk1(ixb+1,-1) = fun(ixp,-1)
        call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,-1),d,
     .              ifail)
        call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,-1),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,-1),ifail)
*    .. region 5 : DN LFS core : central part
        wrk1(ixa-1,iya-1) = fun(ixm,iym)
        wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
        wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
        wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
        wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
        wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
        wrk1(ixb+1,iya-1) = fun(ixp,iym)
        wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
        wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
       endif
*    .. region 6 :
       if (nncut.eq.2) then
        if (isSnowflake) then
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
          iym = 0
          iyc = 0
        else
          ixa = rightcut(2)
          ixb = rightcut(1)-1
          ixc = rightcut1(2)
          ixd = rightcut1(1)-1
          iya = jsep+1
          iym = jsep
          iyc = jsep1+1
          if (topcut(1).eq.topcut(2)) then
            ixm = rightcut(2)-1
            ixp = rightcut(1)
            iyb = ny-1
            iyp = ny-1
            iyd = ny1-1
          else if (isSnowflake) then
            ixm = leftcut(2)-1
            ixp = rightcut(1)
          else if (topcut(1).lt.topcut(2)) then
            ixm = leftcut(2)-1
            ixp = rightcut(1)
          else if (topcut(1).gt.topcut(2)) then
            ixm = rightcut(2)-1
            ixp = leftcut(1)
          endif
        endif
        if (topcut(1).ne.topcut(2)) then   ! DDN HFS SOL
         iyb = max(topcut(1),topcut(2))
         iyp = max(topcut(1),topcut(2))+1
         iyd = max(topcut1(1),topcut1(2))
         do ix = ixa-1, ixb+1
          xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
         enddo
         do ix1 = ixc, ixd
          xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
         enddo
         do iy = iya-1, iyb+1
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
*    .. region 6 : DN LFS SOL : between the separatrices
         wrk1(ixa-1,iya-1) = fun(ixm,iym)
         wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
         wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
         wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
         wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
         wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
         wrk1(ixb+1,iya-1) = fun(ixp,iym)
         wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
         wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
         ixm = rightcut(2)-1
         ixp = rightcut(1)
         iya = max(topcut(1),topcut(2))+1
         iym = max(topcut(1),topcut(2))
         iyc = max(topcut1(1),topcut1(2))+1
         iyb = ny-1
         iyp = ny-1
         iyd = ny1-1
        endif
*    .. region 6 : DN LFS SOL : above the separatrices
        do ix = ixa-1, ixb+1
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        do iy = iya-1, iyb
         yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
        wrk1(ixa-1,iya-1) = fun(ixm,iym)
        wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
        wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
        wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
        wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
        wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
        wrk1(ixb+1,iya-1) = fun(ixp,iym)
        wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
        wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
*    .. region 6 : DN LFS SOL : top row
        wrk1(ixa-1,ny) = fun(ixm,ny)
        wrk1(ixa:ixb,ny) = fun(ixa:ixb,ny)
        wrk1(ixb+1,ny) = fun(ixp,ny)
        call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .              ifail)
        call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,ny1),ifail)
       endif
*    .. region 7 :
       if (nncut.eq.2) then
        if (isSnowflake) then
          ixa = rightcut(2)
          ixm = leftcut(1)-1
          ixc = rightcut1(2)
          ixb = nx-1
          ixp = nx-1
          ixd = nx1-1
        else
          ixb = rightcut(2)-1
          ixp = leftcut(2)
          ixd = rightcut1(2)-1
          ixa = ixb
          do while (region(ixa-1,jsep,0).eq.7)
            ixa = ixa-1
          enddo
          ixa = ixa + 1  ! step back to not count guard cell
          ixm = ixa
          ixc = ixd
          do while (region1(ixc-1,jsep1,0).eq.7)
            ixc = ixc-1
          enddo
          ixc = ixc + 1  ! step back to not count guard cell
        endif
        xx(ixa-1) = 0.0_R8
        do ix = ixa, ixb+1
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        iya = 0
        iym = 0
        iyb = jsep
        iyp = jsep+1
        iyc = 0
        iyd = jsep1
        yy(iya) = 0.0_R8
        do iy = iya, iyb+1
         yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
*    .. region 7 : DN top LFS PFR : bottom row
        wrk1(ixa-1,-1) = fun(ixm,-1)
        wrk1(ixa:ixb,-1) = fun(ixa:ixb,-1)
        wrk1(ixb+1,-1) = fun(ixp,-1)
        call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,-1),d,
     .              ifail)
        call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,-1),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,-1),ifail)
*    .. region 7 : DN top LFS PFR : western edge
        wrk1(ixa-1,iya-1) = fun(ixa-1,iym)
        wrk1(ixa-1,iya:iyb) = fun(ixa-1,iya:iyb)
        wrk1(ixa-1,iyb+1) = fun(ixa-1,iyp)
        call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(ixa-1,iya-1:iyb+1),d,
     .              ifail)
        call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(ixa-1,iya-1:iyb+1),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fun1(ixc-1,iyc:iyd),ifail)
*    .. region 7 : DN top LFS PFR : central part
        wrk1(ixa-1,iya-1) = fun(ixm,iym)
        wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
        wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
        wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
        wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
        wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
        wrk1(ixb+1,iya-1) = fun(ixp,iym)
        wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
        wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
        if (topcut(1).ne.topcut(2)) then   ! DDN HFS SOL
         if (topcut(1).gt.topcut(2) .and. .not.isSnowflake) then
          ixp = rightcut(2)
         endif
         iya = min(topcut(1),topcut(2))
         iym = min(topcut(1),topcut(2))-1
         iyb = max(topcut(1),topcut(2))
         iyp = max(topcut(1),topcut(2))+1
         iyc = min(topcut1(1),topcut1(2))
         iyd = max(topcut1(1),topcut1(2))
         do iy = iya-1, iyb+1
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
*    .. region 7 : DN top LFS divertor : between the separatrices : western edge
         wrk1(ixa-1,iya-1) = fun(ixa-1,iym)
         wrk1(ixa-1,iya:iyb) = fun(ixa-1,iya:iyb)
         wrk1(ixa-1,iyb+1) = fun(ixa-1,iyp)
         call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(ixa-1,iya-1:iyb+1),
     .               d,ifail)
         call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(ixa-1,iya-1:iyb+1),
     .               d,iyd-iyc+1,yy1(iyc:iyd),fun1(ixc-1,iyc:iyd),ifail)
*    .. region 7 : DN top LFS divertor : between the separatrices : central part
         wrk1(ixa-1,iya-1) = fun(ixm,iym)
         wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
         wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
         wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
         wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
         wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
         wrk1(ixb+1,iya-1) = fun(ixp,iym)
         wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
         wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
        endif
        ixp = rightcut(2)
        iya = max(topcut(1),topcut(2))+1
        iym = max(topcut(1),topcut(2))
        iyc = max(topcut1(1),topcut1(2))+1
        iyb = ny-1
        iyp = ny-1
        iyd = ny1-1
        do iy = iya-1, iyb
         yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
*    .. region 7 : DN top LFS divertor : above the separatrices : western edge
        wrk1(-1,iya-1) = fun(-1,iym)
        wrk1(-1,iya:iyb) = fun(-1,iya:iyb)
        wrk1(-1,iyb+1) = fun(-1,iyp)
        call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(-1,iya-1:iyb+1),d,
     .              ifail)
        call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(-1,iya-1:iyb+1),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fun1(-1,iyc:iyd),ifail)
*    .. region 7 : DN top LFS divertor : above the separatrices : central part
        wrk1(ixa-1,iya-1) = fun(ixm,iym)
        wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
        wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
        wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
        wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
        wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
        wrk1(ixb+1,iya-1) = fun(ixp,iym)
        wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
        wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
*    .. region 7 : DN top LFS divertor : above the separatrices : top row
        wrk1(ixa-1,ny) = fun(ixm,ny)
        wrk1(ixa:ixb,ny) = fun(ixa:ixb,ny)
        wrk1(ixb+1,ny) = fun(ixp,ny)
        call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .              ifail)
        call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,ny1),ifail)
       endif
*    .. region 8 :
       if (nncut.eq.2 .and. .not.isSnowflake) then
        ixa = rightcut(1)
        ixm = leftcut(1)-1
        ixc = rightcut1(1)
        ixb = nx-1
        ixp = nx-1
        ixd = nx1-1
        do ix = ixa-1, ixb
         xx(ix) = (ix+0.5_R8-ixa)/real(ixb-ixa+1)
        enddo
        xx(ixb+1) = 1.0_R8
        do ix1 = ixc, ixd
         xx1(ix1) = (ix1+0.5_R8-ixc)/real(ixd-ixc+1)
        enddo
        iya = 0
        iym = 0
        iyb = jsep
        iyp = jsep+1
        iyc = 0
        iyd = jsep1
        yy(iya) = 0.0_R8
        do iy = iya, iyb+1
         yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
*    .. region 8 : eastern PFR : bottom row
        wrk1(ixa-1,-1) = fun(ixm,-1)
        wrk1(ixa:ixb,-1) = fun(ixa:ixb,-1)
        wrk1(ixb+1,-1) = fun(ixp,-1)
        call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,-1),d,
     .              ifail)
        call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,-1),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,-1),ifail)
*    .. region 8 : eastern PFR : eastern edge
        wrk1(ixb+1,iya-1) = fun(ixb+1,iym)
        wrk1(ixb+1,iya:iyb) = fun(ixb+1,iya:iyb)
        wrk1(ixb+1,iyb+1) = fun(ixb+1,iyp)
        call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(ixb+1,iya-1:iyb+1),d,
     .              ifail)
        call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(ixb+1,iya-1:iyb+1),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fun1(ixd+1,iyc:iyd),ifail)
*    .. region 8 : eastern PFR : central part
        wrk1(ixa-1,iya-1) = fun(ixm,iym)
        wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
        wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
        wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
        wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
        wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
        wrk1(ixb+1,iya-1) = fun(ixp,iym)
        wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
        wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
        if (topcut(1).ne.topcut(2)) then
         if (topcut(1).lt.topcut(2)) ixm = rightcut(1)-1
         iya = min(topcut(1),topcut(2))
         iym = min(topcut(1),topcut(2))-1
         iyb = max(topcut(1),topcut(2))
         iyp = max(topcut(1),topcut(2))+1
         iyc = min(topcut1(1),topcut1(2))
         iyd = max(topcut1(1),topcut1(2))
         do iy = iya-1, iyb+1
          yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
         enddo
         do iy1 = iyc, iyd
          yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
         enddo
*    .. region 8 : eastern divertor : between the separatrices : eastern edge
         wrk1(ixb+1,iya-1) = fun(ixb+1,iym)
         wrk1(ixb+1,iya:iyb) = fun(ixb+1,iya:iyb)
         wrk1(ixb+1,iyb+1) = fun(ixb+1,iyp)
         call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(ixb+1,iya-1:iyb+1),
     .               d,ifail)
         call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(ixb+1,iya-1:iyb+1),
     .               d,iyd-iyc+1,yy1(iyc:iyd),fun1(ixd+1,iyc:iyd),ifail)
*    .. region 8 : eastern divertor : between the separatrices : central part
         wrk1(ixa-1,iya-1) = fun(ixm,iym)
         wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
         wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
         wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
         wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
         wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
         wrk1(ixb+1,iya-1) = fun(ixp,iym)
         wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
         wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
         call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .    xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .    wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
        endif
        ixm = rightcut(1)-1
        iya = max(topcut(1),topcut(2))+1
        iym = max(topcut(1),topcut(2))
        iyc = max(topcut1(1),topcut1(2))+1
        iyb = ny-1
        iyp = ny-1
        iyd = ny1-1
        do iy = iya-1, iyb
         yy(iy) = (iy+0.5_R8-iya)/real(iyb-iya+1)
        enddo
        yy(iyb+1) = 1.0_R8
        do iy1 = iyc, iyd
         yy1(iy1) = (iy1+0.5_R8-iyc)/real(iyd-iyc+1)
        enddo
*    .. region 8 : eastern divertor : above the separatrices : eastern edge
        wrk1(ixb+1,iya-1) = fun(ixb+1,iym)
        wrk1(ixb+1,iya:iyb) = fun(ixb+1,iya:iyb)
        wrk1(ixb+1,iyb+1) = fun(ixb+1,iyp)
        call E01BEF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(ixb+1,iya-1:iyb+1),d,
     .              ifail)
        call E01BFF(iyb-iya+3,yy(iya-1:iyb+1),wrk1(ixb+1,iya-1:iyb+1),d,
     .              iyd-iyc+1,yy1(iyc:iyd),fun1(ixd+1,iyc:iyd),ifail)
*    .. region 8 : eastern divertor : above the separatrices : central part
        wrk1(ixa-1,iya-1) = fun(ixm,iym)
        wrk1(ixa-1,iya:iyb) = fun(ixm,iya:iyb)
        wrk1(ixa-1,iyb+1) = fun(ixm,iyp)
        wrk1(ixa:ixb,iyb+1) = fun(ixa:ixb,iyp)
        wrk1(ixb+1,iyb+1) = fun(ixp,iyp)
        wrk1(ixb+1,iya:iyb) = fun(ixp,iya:iyb)
        wrk1(ixb+1,iya-1) = fun(ixp,iym)
        wrk1(ixa:ixb,iya-1) = fun(ixa:ixb,iym)
        wrk1(ixa:ixb,iya:iyb) = fun(ixa:ixb,iya:iyb)
        call interp2d(ixb-ixa+3,iyb-iya+3,ixd-ixc+1,iyd-iyc+1,
     .   xx(ixa-1:ixb+1),yy(iya-1:iyb+1),xx1(ixc:ixd),yy1(iyc:iyd),
     .   wrk1(ixa-1:ixb+1,iya-1:iyb+1),fun1(ixc:ixd,iyc:iyd))
*    .. region 8 : eastern divertor : above the separatrices : top row
        wrk1(ixa-1,ny) = fun(ixm,ny)
        wrk1(ixa:ixb,ny) = fun(ixa:ixb,ny)
        wrk1(ixb+1,ny) = fun(ixp,ny)
        call E01BEF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .              ifail)
        call E01BFF(ixb-ixa+3,xx(ixa-1:ixb+1),wrk1(ixa-1:ixb+1,ny),d,
     .              ixd-ixc+1,xx1(ixc:ixd),fun1(ixc:ixd,ny1),ifail)
       endif
*    ..clean-up
       deallocate(d)
      endif

      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xtvc
