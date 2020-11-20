*-----------------------------------------------------------------------
*.specification

      subroutine b2xgxy (nx, ny, ixr, iyr, xc, xf, yc, yf)
      use b2mod_types
      use b2mod_geo
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ixr, iyr
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
*     This routine computes cartesian (x,y) coordinates on the B2
*     geometry, by integrating the line element along a pair of
*     reference lines.
*     These coordinates may be used as cartesian axes for graphical
*     or printed output. They have no role in the actual computation
*     performed in B2.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny - integer, input.
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
      integer ix, iy
*   ..procedures
      external subini, subend, xertst

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
      xf(ixr) = 0.0_R8
      do ix = ixr-1, -1, -1
       xc(ix) = xf(ix+1)-disx(ix,iyr)/2.0_R8
       xf(ix) = xf(ix+1)-disx(ix,iyr)
      enddo
      do ix = ixr+1, nx
       xc(ix-1) = xf(ix-1)+disx(ix-1,iyr)/2.0_R8
       xf(ix) = xf(ix-1)+disx(ix-1,iyr)
      enddo
      xc(nx) = xf(nx)+disx(nx,iyr)/2.0_R8
*   ..compute y-coordinates
*     (integrate distance across cell faces at ix=ixr)
      yf(iyr) = 0.0_R8
      do iy = iyr-1, -1, -1
       yc(iy) = yf(iy+1)-disy(ixr,iy)/2.0_R8
       yf(iy) = yf(iy+1)-disy(ixr,iy)
      enddo
      do iy = iyr+1, ny
       yc(iy-1) = yf(iy-1)+disy(ixr,iy-1)/2.0_R8
       yf(iy) = yf(iy-1)+disy(ixr,iy-1)
      enddo
      yc(ny) = yf(ny)+disy(ixr,ny)/2.0_R8

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xgxy

      end subroutine b2xgxy
