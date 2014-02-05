*-----------------------------------------------------------------------
*.specification

      subroutine b2wgxp (nx, ny, header, xc, xf, yc, yf, fun, iloc,
     &  scal, fa0, fa1, ncnta, fb0, fb1, ncntb)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, iloc, ncnta, ncntb
      character header*(*)
      real (kind=R8) ::
     *   xc(-1:nx), xf(-1:nx), yc(-1:ny), yf(-1:ny),
     &  fun(-1:nx,-1:ny), scal, fa0, fa1, fb0, fb1
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks


*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WGXP produces contour plots on a cartesian domain.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine produces a contour plot of some b2 quantity on a
*     rectangular mesh.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny - integer, input.
*     See the calling routine for description.
*
*     header - character*(*), input.
*     header specifies a label for the plots.
*
*     xc, xf - (-1:nx) real*8 array, input.
*     yc, yf - (-1:ny) real*8 array, input.
*     See the calling routine for description.
*
*     fun - (-1:nx,-1:ny) real*8 array, input.
*     fun specifies the function of which a contour plot is to be made.
*     See the description under iloc.
*
*     iloc - integer, input.
*     iloc specifies the class of locations on which fun is defined.
*     iloc.eq.0: fun(-1:nx,-1:ny) is defined at cell centers.
*     iloc.eq.1: fun( 0:nx,-1:ny) is defined at cell x-faces.
*     iloc.eq.2: fun(-1:nx, 0:ny) is defined at cell y-faces.
*     iloc.eq.3: fun( 0:nx, 0:ny) is defined at cell corners.
*     Other values of iloc are not allowed.
*
*     scal - real, input.
*     scal specifies a scaling factor for the function fun. The
*     contours will be labelled with the value of fun/scal.
*     It is required that 0.lt.scal.
*
*     fa0, fa1 - real, input.
*     fa0 and fa1 specify the extreme values of fun for the solid
*     contour region.
*
*     ncnta - integer, input.
*     ncnta specifies the number of solid contour intervals.
*     If 0.lt.ncnta a solid contour will be plotted at the values
*       fun = ((ncnta-i)*fa0+i*fa1)/ncnta, for 0.le.i.le.ncnta.
*     It is required that 0.le.ncnta.
*
*     fb0, fb1 - real, input.
*     fb0 and fb1 specify the extreme values of fun for the dashed
*     contour region.
*
*     ncntb - integer, input.
*     ncntb specifies the number of dashed contour intervals.
*     If 0.lt.ncntb a dashed contour will be plotted at the values
*       fun = ((ncntb-i)*fb0+i*fb1)/ncntb, for 0.le.i.le.ncntb.
*     It is required that 0.le.ncntb.
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix0, iy0, ix1, iy1, ixb, iyb
      character headr2*120
*   ..procedures
      external subini, subend, xertst, pagini, pagend, cpseti, gndfm2
      external b2xxmm, b2wgyp

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wgxp')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test iloc
      call xertst (0.le.iloc.and.iloc.lt.4, 'faulty argument iloc')
*   ..test scal
      call xertst (0.lt.scal, 'faulty argument scal')
*   ..test ncnta, ncntb
      call xertst (0.le.ncnta.and.0.le.ncntb,
     &  'faulty argument ncnta, ncntb')
*   ..find min and max
      call b2xxmm (nx, ny, fun, iloc, ix0, iy0, ix1, iy1)

* ..make plot
*   ..construct two-line headr2
      headr2(1:60) = header
      if (iloc.eq.0) then
        write(headr2(61:120),101)
     &    fun(ix0,iy0)/scal, xc(ix0), yc(iy0),
     &    fun(ix1,iy1)/scal, xc(ix1), yc(iy1)
      else if (iloc.eq.1) then
        write(headr2(61:120),101)
     &    fun(ix0,iy0)/scal, xf(ix0), yc(iy0),
     &    fun(ix1,iy1)/scal, xf(ix1), yc(iy1)
      else if (iloc.eq.2) then
        write(headr2(61:120),101)
     &    fun(ix0,iy0)/scal, xc(ix0), yf(iy0),
     &    fun(ix1,iy1)/scal, xc(ix1), yf(iy1)
      else if (iloc.eq.3) then
        write(headr2(61:120),101)
     &    fun(ix0,iy0)/scal, xf(ix0), yf(iy0),
     &    fun(ix1,iy1)/scal, xf(ix1), yf(iy1)
      endif
*   ..define range and mapping
      if (iloc.eq.0) then
       ixb = -1
       iyb = -1
       call gndfm2 (xc(-1), nx+1, yc(-1), ny+1)
      else if (iloc.eq.1) then
       ixb = 0
       iyb = -1
       call gndfm2 (xf(0), nx, yc(-1), ny+1)
      else if (iloc.eq.2) then
       ixb = -1
       iyb = 0
       call gndfm2 (xc(-1), nx+1, yf(0), ny)
      else if (iloc.eq.3) then
       ixb = 0
       iyb = 0
       call gndfm2 (xf(0), nx, yf(0), ny)
      endif
*   ..make plot
      call pagini (4, xc(-1), xc(nx), yc(-1), yc(ny), 'x', 'y',
     &  headr2)
      call cpseti ('MAP', 2)
      call b2wgyp (nx-ixb, ny-iyb, nx+1, fun(ixb,iyb), scal,
     &  fa0, fa1, ncnta, fb0, fb1, ncntb)
      call pagend ()

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.formats

  101 format (
     &  1pe9.2,' at (',0pf6.2,',',f6.2,'), ',
     &  1pe9.2,' at (',0pf6.2,',',f6.2,').')

*-----------------------------------------------------------------------
*.end b2wgxp

      end
