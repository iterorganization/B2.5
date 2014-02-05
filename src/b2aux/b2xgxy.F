*-----------------------------------------------------------------------
*.specification

      subroutine b2xgxy (nx, ny, crx, cry, ixr, iyr, xc, xf, yc, yf)
      use b2mod_types
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ixr, iyr
      real (kind=R8) ::
     *   crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   xc(-1:nx), xf(-1:nx), yc(-1:ny), yf(-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2XGXY computes cartesian coordinates suitable for graphical and
*     printed output.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes cartesian (x,y) coordinates on the b2
*     geometry, by integrating the line element along a pair of
*     reference lines.
*     These coordinates may be used as cartesian axes for graphical
*     or printed output. They have no role in the actual computation
*     performed in b2.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny - integer, input.
*     crx, cry - (-1:nx,-1:ny,0:3) real array, input.
*     See the calling routine for description.
*
*     ixr - integer, input.
*     ixr specifies the location of the reference line that is used to
*     define the y-coordinate (for output purposes only); this reference
*     line is the boundary of the cells between ix=ixr-1 and ix=ixr.
*     It will hold that 0.le.ixr.le.nx.
*
*     iyr - integer, input.
*     iyr specifies the location of the reference line that is used to
*     define the x-coordinate (for output purposes only); this reference
*     line is the boundary of the cells between iy=iyr-1 and iy=iyr.
*     It will hold that 0.le.iyr.le.ny.
*
*     xc - (-1:nx) real array, output.
*     For -1.le.ix.le.nx, xc(ix) specifies the cartesian x-coordinate
*     for the centers of cells (ix,iy), where iy is arbitrary.
*
*     xf - (-1:nx) real array, output.
*     For 0.le.ix.le.nx, xf(ix) specifies the cartesian x-coordinate
*     for the x-faces between cells (ix-1,iy) and cells (ix,iy), where
*     iy is arbitrary. xf(-1) is arbitrary.
*
*     yc - (-1:ny) real array, output.
*     For -1.le.iy.le.ny, yc(iy) specifies the cartesian y-coordinate
*     for the centers of cells (ix,iy), where ix is arbitrary.
*
*     yf - (-1:ny) real array, output.
*     For 0.le.iy.le.ny, yf(iy) specifies the cartesian y-coordinate
*     for the y-faces between cells (ix,iy-1) and cells (ix,iy), where
*     ix is arbitrary. yf(-1) is arbitrary.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix, iy, i, j
*   ..procedures
      intrinsic sqrt
      real (kind=R8) ::
     *   dist, disx, disy
      external subini, subend, xertst
      dist(ix,iy,i,j) = sqrt (
     &  (crx(ix,iy,j)-crx(ix,iy,i))**2+(cry(ix,iy,j)-cry(ix,iy,i))**2)
      disx(ix,iy) = (dist(bottomix(ix,iy),bottomiy(ix,iy),2,3)+
     1 dist(ix,iy,0,1))/2
      disy(ix,iy) = (dist(leftix(ix,iy),leftiy(ix,iy),1,3)+
     1 dist(ix,iy,0,2))/2

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2xgxy')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test ixr, iyr
      call xertst (0.le.ixr.and.ixr.le.nx, 'faulty argument ixr')
      call xertst (0.le.iyr.and.iyr.le.ny, 'faulty argument iyr')

* ..compute the coordinates
*   ..compute x-coordinates
*     (integrate distance across cell faces at iy=iyr)
      xf(ixr) = 0
      do ix = ixr-1, -1, -1
       xc(ix) = xf(ix+1)-disx(ix,iyr)/2
       xf(ix) = xf(ix+1)-disx(ix,iyr)
      enddo
      do ix = ixr+1, nx
       xc(ix-1) = xf(ix-1)+disx(ix-1,iyr)/2
       xf(ix) = xf(ix-1)+disx(ix-1,iyr)
      enddo
      xc(nx) = xf(nx)+disx(nx,iyr)/2
*   ..compute y-coordinates
*     (integrate distance across cell faces at ix=ixr)
      yf(iyr) = 0
      do iy = iyr-1, -1, -1
       yc(iy) = yf(iy+1)-disy(ixr,iy)/2
       yf(iy) = yf(iy+1)-disy(ixr,iy)
      enddo
      do iy = iyr+1, ny
       yc(iy-1) = yf(iy-1)+disy(ixr,iy-1)/2
       yf(iy) = yf(iy-1)+disy(ixr,iy-1)
      enddo
      yc(ny) = yf(ny)+disy(ixr,ny)/2

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xgxy

      end subroutine b2xgxy
