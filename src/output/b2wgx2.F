*-----------------------------------------------------------------------
*.specification

      subroutine b2wgx2 (nx, ny, header, fun, iloc, scal,
     &  fa0, fa1, ncnta, fb0, fb1, ncntb)
      use b2mod_types
      use b2mod_indirect
      use b2mod_b2cmwg
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, iloc, ncnta, ncntb
      character header*(*)
      real (kind=R8) ::
     *  fun(-1:nx,-1:ny)
      real (kind=R8) ::
     *  scal, fa0, fa1, fb0, fb1
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WGX2 produces contour plots on a sequence of cartesian domains.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine produces a contour plot of some quantity on a
*     sequence of cartesian representations of the B2 mesh.
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
*     /b2cmwg/ - input.
*     /b2cmwg/ contains some parameters to control graphical output.
*      nmbpmx - integer constant.
*     nmbpmx is used in dimension declarations instead of nmbpcc or
*     nmbpdm.
*      nmbpcc - integer.
*     nmbpcc specifies the number of cartesian domains on which any
*     quantity may be plotted.
*      nmbpdm - integer.
*     nmbpdm specifies the number of physical domains on which any
*     quantity may be plotted.
*      cwaxr, cwayr - (0:nmbpmx-1) real*8 array.
*     For each iplot in (0:nmbpcc-1), cwaxr(iplot) and cwayr(iplot)
*     specify coordinates for a cartesian domain plot. The coordinates
*     are obtained by integrating the grid metric along two reference
*     lines. The reference line for the x-coordinate is at position
*     nint(cwaxr*nx), and the reference line for the y-coordinate is
*     at nint(cwayr*ny). It will hold that 0.le.cwaxr.le.1 and
*     0.le.cwayr.le.1.
*      cwxmin, cwxmax, cwymin, cwymax - (0:nmbpmx-1) real*8 array.
*     For each iplot in (0:nmbpdm-1), the quadruple of cwxmin(iplot)
*     through cwymax(iplot) specifies a domain in the physical plane
*     for contour plots.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iplot, ixr, iyr
      real (kind=R8) ::
     *   xf(-1:nx),xc(-1:nx),yf(-1:ny),yc(-1:ny)
*   ..procedures
      intrinsic nint
      external xertst
      external b2wgxp, b2xgxy

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wgx2')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test iloc
      call xertst (0.le.iloc.and.iloc.lt.4, 'faulty argument iloc')
*   ..test scal
      call xertst (0.lt.scal, 'faulty argument scal')
*   ..test ncnta, ncntb
      call xertst (0.le.ncnta.and.0.le.ncntb,
     &  'faulty argument ncnta, ncntb')

* ..make plots
*   ..loop over plot index
      do iplot = 0, nmbpcc-1
*    ..obtain the coordinates
       ixr = nint(cwaxr(iplot)*nx)
       iyr = nint(cwayr(iplot)*ny)
       call b2xgxy (nx, ny, ixr, iyr, xc, xf, yc, yf)
*    ..make plot
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   fun, iloc, scal, fa0, fa1, ncnta, fb0, fb1, ncntb)
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wgx2

      end subroutine b2wgx2
