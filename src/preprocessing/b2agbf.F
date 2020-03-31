*-----------------------------------------------------------------------
*.specification

      subroutine b2agbf (nx, ny, isymm, periodic_bc, nncut,
     & crx, cry, fpsi, ffbz, vol, gs, qz, qc, qs, 
     & bb, wbbl, wbbr, wbbv, region, leftix, leftiy,
     & bottomix, bottomiy, rightix, rightiy, topix, topiy, cflag)
      use b2mod_types
      use b2mod_constants
      use tradui_constants
      use b2mod_connectivity
      use b2mod_cellhelper
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isymm, periodic_bc, nncut
      integer region(-1:nx,-1:ny,0:2), 
     &  leftix(-1:nx,-1:ny), leftiy(-1:nx,-1:ny), 
     &  bottomix(-1:nx,-1:ny), bottomiy(-1:nx,-1:ny),
     &  rightix(-1:nx,-1:ny), rightiy(-1:nx,-1:ny),
     &  topix(-1:nx,-1:ny), topiy(-1:nx,-1:ny),
     $  cflag(-1:nx,-1:ny,GRID_N_CELLFLAGS) 
      real (kind=R8) ::
     &  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     &  fpsi(-1:nx,-1:ny,0:3,0:2), ffbz(-1:nx,-1:ny,0:3),
     &  vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2),
     &  qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny,0:1), qs(-1:nx,-1:ny,0:1),
     &  bb(-1:nx,-1:ny,0:3)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  wbbl(-1:nx,-1:ny,0:3), wbbr(-1:nx,-1:ny,0:3), 
     &  wbbv(-1:nx,-1:ny,0:3)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AGBF computes the magnetic field at cell faces.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the three components of the magnetic field
*     and the total magnetic field strength, all at cell faces. The
*     components (bx,by) are obtained by differencing the potential fpsi
*     from its values at cell corners. The component bz is obtained
*     directly from the field function ffbz.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny - integer, input.
*     See the calling routine for description.
*
*     vol, crx, cry, qz - real*8 array, input.
*     See the calling routine for description.
*
*     fpsi - (-1:nx,-1:ny,0:3) real*8 array, input.
*     For (ix,iy,i) in (-1:nx,-1:ny,0:3), fpsi(ix,iy,i) specifies the
*     potential for the magnetic field components (bx,by) at the
*     (ix,iy,i) position.
*     (In the case of straight symmetry, bx=dpsi/dy and by=-dpsi/dx.
*     In the case of toroidal or cylindrical symmetry,
*       bx=(dpsi/dy)/(2*pi*R) and by=-(dpsi/dx)/(2*pi*R),
*     where R is the distance to the axis of symmetry.)
*
*     ffbz - (-1:nx,-1:ny,0:3) real*8 array, input.
*     For (ix,iy,i) in (-1:nx,-1:ny,0:3), ffbz(ix,iy,i) specifies the
*     field function for the magnetic field component bz at the
*     (ix,iy,i) position.
*     (In the case of straight symmetry, bz=ffbz. In the case of
*     toroidal or cylindrical symmetry, bz=ffbz/(2*pi*R), where R is
*     the distance to the axis of symmetry.)
*
*     bb - (-1:nx,-1:ny,0:3) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), bb(ix,iy,0:3) specifies the
*     magnetic field at the center of the (ix,iy) cell.
*     bb(ix,iy,0:2) are the (x,y,z)-components and bb(ix,iy,3) is
*     the absolute magnetic field strength.
*     (bb(,,3)=sqrt(bb(,,0)**2+bb(,,1)**2+bb(,,2)**2.)
*     It will hold that 0.lt.bb(,,3).
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix, iy, i, j
      real (kind=R8) ::
     *  t0, p0x, p0y, width, bb0, babs
      real (kind=R8) :: hzl(-1:nx,-1:ny),    ! hz at the left cell face
     *                  hzr(-1:nx,-1:ny),    ! hz at the right cell face
     *                  hzb(-1:nx,-1:ny)     ! hz at the bottom cell face
*   ..procedures
      intrinsic sqrt
      real (kind=R8) ::
     *  norm, cornerDist
      norm(p0x,p0y) = sqrt(p0x**2+p0y**2)
      cornerDist(ix,iy,i,j) = norm (
     &  crx(ix,iy,j)-crx(ix,iy,i), cry(ix,iy,j)-cry(ix,iy,i))
      external subini, subend, xertst
      data width /1.0_R8/
      save width

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agbf')
      call ipgetr ('b2agmt_1d_width', width)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

*   ..initialize
      wbbl = 0.0_R8
      wbbr = 0.0_R8
      wbbv = 0.0_R8

*    ..interpolate magnetic field to cell faces
      if (maxval(abs(fpsi)).eq.0.0_R8 .or.
     &    maxval(abs(ffbz)).eq.0.0_R8 ) then
        !wdk 'standard' treatment: inverse volume weighting
        call intfaceh1(nx,ny,2,leftix,leftiy,vol,bb,wbbl)
        call intfacer1(nx,ny,2,rightix,rightiy,vol,bb,wbbr)
        call intfacev1(nx,ny,2,bottomix,bottomiy,vol,bb,wbbv)
      else
        !wdk try to be more precise in case we know full fpsi/ffbz:
        !wdk direct computation by averaging from vertext values 
*     ..preparatory quantities
        hzl = INVALID_DOUBLE
        hzr = INVALID_DOUBLE
        hzb = INVALID_DOUBLE
        do iy = -1, ny
         do ix = -1, nx
          if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
          if(isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
           if (isymm.eq.0) then
             hzl(ix,iy) = width
           else if (isymm.eq.1 .or. isymm.eq.2) then
             hzl(ix,iy) = pi * ( crx(ix,iy,0) + crx(ix,iy,2) )
           else if (isymm.eq.3 .or. isymm.eq.4) then
             hzl(ix,iy) = pi * ( cry(ix,iy,0) + cry(ix,iy,2) )
           end if
          endif
          if(isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
           if (isymm.eq.0) then
             hzr(ix,iy) = width
           else if (isymm.eq.1 .or. isymm.eq.2) then
             hzr(ix,iy) = pi * ( crx(ix,iy,1) + crx(ix,iy,3) )
           else if (isymm.eq.3 .or. isymm.eq.4) then
             hzr(ix,iy) = pi * ( cry(ix,iy,1) + cry(ix,iy,3) )
           end if
          endif
          if(isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
           if (isymm.eq.0) then
             hzb(ix,iy) = width
           else if (isymm.eq.1 .or. isymm.eq.2) then
             hzb(ix,iy) = pi * ( crx(ix,iy,0) + crx(ix,iy,1) )
           else if (isymm.eq.3 .or. isymm.eq.4) then
             hzb(ix,iy) = pi * ( cry(ix,iy,0) + cry(ix,iy,1) )
           end if
          endif
         enddo
        enddo
*     ..computation of magnetic field components
        do iy = -1, ny
         do ix = -1, nx
          if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
          bb0 = (fpsi(ix,iy,2,0) - fpsi(ix,iy,0,0) + 
     .     fpsi(ix,iy,3,0)-fpsi(ix,iy,1,0))
          if(isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
           babs = 0.5_R8*(
     .       sqrt(fpsi(ix,iy,0,1)**2+fpsi(ix,iy,0,2)**2) + 
     .       sqrt(fpsi(ix,iy,2,1)**2+fpsi(ix,iy,2,2)**2) )/hzl(ix,iy)
           wbbl(ix,iy,0) = sign(babs,bb0)
           wbbl(ix,iy,1) = 0.0_R8
           wbbl(ix,iy,2) = (ffbz(ix,iy,0) + ffbz(ix,iy,2)) /
     .      (2.0_R8 * hzl(ix,iy))
          endif
          if(isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
           babs = 0.5_R8*(
     .       sqrt(fpsi(ix,iy,1,1)**2+fpsi(ix,iy,1,2)**2) + 
     .       sqrt(fpsi(ix,iy,3,1)**2+fpsi(ix,iy,3,2)**2) )/hzr(ix,iy)
           wbbr(ix,iy,0) = sign(babs,bb0)
           wbbr(ix,iy,1) = 0.0_R8
           wbbr(ix,iy,2) = (ffbz(ix,iy,1) + ffbz(ix,iy,3)) /
     .      (2.0_R8 * hzr(ix,iy))
          endif
          if(isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
           babs = 0.5_R8*(
     .       sqrt(fpsi(ix,iy,0,1)**2+fpsi(ix,iy,0,2)**2) + 
     .       sqrt(fpsi(ix,iy,1,1)**2+fpsi(ix,iy,1,2)**2) )/hzb(ix,iy)
           wbbv(ix,iy,0) = sign(babs,bb0)
           wbbv(ix,iy,1) = 0.0_R8
           wbbv(ix,iy,2) = (ffbz(ix,iy,0) + ffbz(ix,iy,1)) /
     .      (2.0_R8 * hzb(ix,iy))
          endif
         enddo
        enddo
      endif

*    ..extrapolate magnetic field to edges
      write (*,*) 'b2agfs -- not computing edges; not needed for US?'
#ifdef WG_TODO
      do iy = -1, ny
        do ix = -1, nx
          if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
          write(*,*) 'b2agbf ix, iy', ix, iy
          if (leftix(ix,iy).eq.-2) wbbl(ix,iy,0:2) = bb(ix,iy,0:2) -
     &      (wbbr(ix,iy,0:2)-bb(ix,iy,0:2))
          if (rightix(ix,iy).eq.nx+1) wbbr(ix,iy,0:2) = 
     &     bb(ix,iy,0:2) - (wbbl(ix,iy,0:2)-bb(ix,iy,0:2))
          if (bottomiy(ix,iy).eq.-2) wbbv(ix,iy,0:2) = bb(ix,iy,0:2)
     &     - (wbbv(topix(ix,iy),topiy(ix,iy),0:2)-bb(ix,iy,0:2))
        enddo
      enddo
*   .. correct interpolation of magnetic field for diverted island case
*      as B is not a periodic quantity
      if(periodic_bc.eq.1 .and. nncut.ge.1) then
         do iy = -1, ny
           do ix = -1, nx
             if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
             if(leftix(ix,iy).ne.ix-1) then
               if(region(ix,iy,0).eq.1 .or. region(ix,iy,0).eq.4) then
                 wbbl(ix,iy,0)=bb(ix,iy,0)-
     .            (wbbr(ix,iy,0)-bb(ix,iy,0))
                 wbbl(ix,iy,1)=bb(ix,iy,1)-
     .            (wbbr(ix,iy,1)-bb(ix,iy,1))
                 wbbl(ix,iy,2)=bb(ix,iy,2)-
     .            (wbbr(ix,iy,2)-bb(ix,iy,2))
               endif
             else if(rightix(ix,iy).ne.ix+1) then
               if(region(ix,iy,0).eq.1 .or. region(ix,iy,0).eq.3) then
                 wbbr(ix,iy,0)=bb(ix,iy,0)-
     .            (wbbl(ix,iy,0)-bb(ix,iy,0))
                 wbbr(ix,iy,1)=bb(ix,iy,1)-
     .            (wbbl(ix,iy,1)-bb(ix,iy,1))
                 wbbr(ix,iy,2)=bb(ix,iy,2)-
     .            (wbbl(ix,iy,2)-bb(ix,iy,2))
               endif
             endif
           enddo
         enddo
      endif
#endif
*    ..compute extrapolated wbb.(:,:,3) components
      wbbl(:,:,3)=sqrt(wbbl(:,:,0)**2+wbbl(:,:,1)**2+wbbl(:,:,2)**2)
      wbbr(:,:,3)=sqrt(wbbr(:,:,0)**2+wbbr(:,:,1)**2+wbbr(:,:,2)**2)
      wbbv(:,:,3)=sqrt(wbbv(:,:,0)**2+wbbv(:,:,1)**2+wbbv(:,:,2)**2)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agbf

      end subroutine b2agbf
