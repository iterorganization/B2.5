*-----------------------------------------------------------------------
*.specification

      subroutine b2wgx0 (nx, ny, header)
      use b2mod_types
      use b2mod_b2cmwg
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      character header*(*)
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WGX0 displays the B2 mesh in a cartesian representation.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine displays the B2 mesh in a sequence of cartesian
*     representations.
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

*   ..dimensioning parameter
      integer nmax
      parameter (nmax=300)
*   ..local variables
      integer ix, iy, iplot, ixr, iyr
      real (kind=R8) ::
     *   xc(-1:nmax), xf(-1:nmax), yc(-1:nmax), yf(-1:nmax),
     &  xlohi(0:1), ylohi(0:1)
*   ..procedures
      intrinsic nint
      external subini, subend, xertst, pagini, pagend, gsln, gpl_d
      external b2xgxy

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wgx0')
*   ..test nx, ny, nmax
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (nx.le.nmax.and.ny.le.nmax,
     &  'faulty parameter nmax')

* ..display the cartesian geometry
*   ..loop over plot index
      do iplot = 0, nmbpcc-1
*    ..obtain the coordinates
       ixr = nint(cwaxr(iplot)*nx)
       iyr = nint(cwayr(iplot)*ny)
       call b2xgxy (nx, ny, ixr, iyr, xc, xf, yc, yf)
*    ..plot the grid
       call pagini (4, xc(-1), xc(nx), yc(-1), yc(ny), 'x', 'y',
     &   header)
       ylohi(0) = yc(-1)
       ylohi(1) = yc(ny)
       do ix = 0, nx
        xlohi(0) = xf(ix)
        xlohi(1) = xf(ix)
        if (mod(4*ix,nx).eq.0) then
         call gsln (1)
         call gpl_d (2, xlohi, ylohi)
        else
         call gsln (2)
         call gpl_d (2, xlohi, ylohi)
        endif
       enddo
       xlohi(0) = xc(-1)
       xlohi(1) = xc(nx)
       do iy = 0, ny
        ylohi(0) = yf(iy)
        ylohi(1) = yf(iy)
        if (mod(4*iy,ny).eq.0) then
         call gsln (1)
         call gpl_d (2, xlohi, ylohi)
        else
         call gsln (2)
         call gpl_d (2, xlohi, ylohi)
        endif
       enddo
       call pagend ()
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wgx0

      end subroutine b2wgx0
