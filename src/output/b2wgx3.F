*-----------------------------------------------------------------------
*.specification

      subroutine b2wgx3 (nx, ny, header, crx, cry, fun, iloc, scal,
     &  fa0, fa1, ncnta, fb0, fb1, ncntb, wrkx, wrky)
      use b2mod_types
      use b2mod_b2cmwg
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, iloc, ncnta, ncntb
      character header*(*)
      real (kind=R8) ::
     *   crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3), 
     &  fun(-1:nx,-1:ny)
      real (kind=R8) ::
     *   scal, fa0, fa1, fb0, fb1
*   ..output arguments (unspecified on entry)
*     (none)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   wrkx(-1:nx,-1:ny), wrky(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WGX3 produces contour plots on the physical domain.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine produces a contour plot of some quantity on a
*     sequence of domains in the physical plane.
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
*     crx, cry - (-1:nx,-1:ny,0:3) real*8 array, input.
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
*     wrkx, wrky - (-1:nx,-1:ny) real*8 array, workspace.
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
      integer iplot, ix, iy, ix0, iy0, ix1, iy1, ixb, iyb
      character headr2*120
*   ..procedures
      external xertst, pagini, pagend, cpseti
      external b2xxmm, b2wgxb, b2wgyp

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wgx3')
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
*   ..compute coordinates
      if (iloc.eq.0) then
       do iy = -1, ny
        do ix = -1, nx
         wrkx(ix,iy) = (crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+
     &     crx(ix,iy,3))/4
         wrky(ix,iy) = (cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+
     &     cry(ix,iy,3))/4
        enddo
       enddo
      else if (iloc.eq.1) then
       do iy = -1, ny
        do ix = -1, nx
         wrkx(ix,iy) = (crx(ix,iy,0)+crx(ix,iy,2))/2
         wrky(ix,iy) = (cry(ix,iy,0)+cry(ix,iy,2))/2
        enddo
       enddo
      else if (iloc.eq.2) then
       do iy = -1, ny
        do ix = -1, nx
         wrkx(ix,iy) = (crx(ix,iy,0)+crx(ix,iy,1))/2
         wrky(ix,iy) = (cry(ix,iy,0)+cry(ix,iy,1))/2
        enddo
       enddo
      else if (iloc.eq.3) then
       do iy = -1, ny
        do ix = -1, nx
         wrkx(ix,iy) = crx(ix,iy,0)
         wrky(ix,iy) = cry(ix,iy,0)
        enddo
       enddo
      endif
*   ..construct two-line headr2
      headr2(1:60) = header
      write(headr2(61:120),101)
     &  fun(ix0,iy0)/scal, wrkx(ix0,iy0), wrky(ix0,iy0),
     &  fun(ix1,iy1)/scal, wrkx(ix1,iy1), wrky(ix1,iy1)
*   ..define range and mapping
      if (iloc.eq.0) then
       ixb = -1
       iyb = -1
      else if (iloc.eq.1) then
       ixb = 0
       iyb = -1
      else if (iloc.eq.2) then
       ixb = -1
       iyb = 0
      else if (iloc.eq.3) then
       ixb = 0
       iyb = 0
      endif
      call gndfm3 (nx+1, ny+1, nx-ixb, ny-iyb, wrkx, wrky)

* ..make plots
*   ..loop over plot index
      do iplot = 0, nmbpdm-1
*    ..make plot
       call pagini (4, cwxmin(iplot), cwxmax(iplot), cwymin(iplot),
     &   cwymax(iplot), 'crx', 'cry', headr2)
       call b2wgxb (nx, ny, crx, cry)
       call cpseti ('MAP', 3)
       call b2wgyp (nx-ixb, ny-iyb, nx+1, fun(ixb,iyb), scal,
     &   fa0, fa1, ncnta, fb0, fb1, ncntb)
       call pagend ()
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.formats

  101 format (
     &  1pe9.2,' at (',0pf6.2,',',f6.2,'), ',
     &  1pe9.2,' at (',0pf6.2,',',f6.2,').')

*-----------------------------------------------------------------------
*.end b2wgx3

      end
