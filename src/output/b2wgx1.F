*-----------------------------------------------------------------------
*.specification

      subroutine b2wgx1 (nx, ny, header, crx, cry)
      use b2mod_types
      use b2mod_b2cmwg
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      character header*(*)
      real (kind=R8) ::
     *   crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3)
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WGX1 displays the B2 mesh on the physical plane.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine displays the B2 mesh on a sequence of physical
*     domains.
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
      integer ix, iy, iplot
      real (kind=R8) ::
     *   xlohi(0:1), ylohi(0:1)
*   ..procedures
      external xertst, pagini, pagend, gpm_d, gpl_d, gsln
      external b2wgxb

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wgx1')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

* ..display the geometry in the physical plane
*   ..loop over plot index
      do iplot = 0, nmbpdm-1
*    ..start plot page
*      (first plot marks vertices)
       call pagini (4, cwxmin(iplot), cwxmax(iplot), cwymin(iplot),
     &   cwymax(iplot), 'crx', 'cry', header)
*    ..plot grid outline
       call b2wgxb (nx, ny, crx, cry)
*    ..plot complete grid
       do ix = -1, nx
        do iy = -1, ny
*      ..left face
         xlohi(0) = crx(ix,iy,0)
         ylohi(0) = cry(ix,iy,0)
         xlohi(1) = crx(ix,iy,2)
         ylohi(1) = cry(ix,iy,2)
         call gpm_d (2, xlohi, ylohi)
*      ..right face
         xlohi(0) = crx(ix,iy,1)
         ylohi(0) = cry(ix,iy,1)
         xlohi(1) = crx(ix,iy,3)
         ylohi(1) = cry(ix,iy,3)
         call gpm_d (2, xlohi, ylohi)
*      ..bottom face
         xlohi(0) = crx(ix,iy,0)
         ylohi(0) = cry(ix,iy,0)
         xlohi(1) = crx(ix,iy,1)
         ylohi(1) = cry(ix,iy,1)
         call gpm_d (2, xlohi, ylohi)
*      ..top face
         xlohi(0) = crx(ix,iy,2)
         ylohi(0) = cry(ix,iy,2)
         xlohi(1) = crx(ix,iy,3)
         ylohi(1) = cry(ix,iy,3)
         call gpm_d (2, xlohi, ylohi)
        enddo
       enddo
       call pagend ()
*    ..start plot page
*      (second plot draws edges)
       call pagini (4, cwxmin(iplot), cwxmax(iplot), cwymin(iplot),
     &   cwymax(iplot), 'crx', 'cry', header)
*    ..plot grid outline
       call b2wgxb (nx, ny, crx, cry)
*    ..change to dashed style
       call gsln (2)
*    ..plot complete grid
       do ix = -1, nx
        do iy = -1, ny
*      ..left face
         xlohi(0) = crx(ix,iy,0)
         ylohi(0) = cry(ix,iy,0)
         xlohi(1) = crx(ix,iy,2)
         ylohi(1) = cry(ix,iy,2)
         call gpl_d (2, xlohi, ylohi)
*      ..right face
         xlohi(0) = crx(ix,iy,1)
         ylohi(0) = cry(ix,iy,1)
         xlohi(1) = crx(ix,iy,3)
         ylohi(1) = cry(ix,iy,3)
         call gpl_d (2, xlohi, ylohi)
*      ..bottom face
         xlohi(0) = crx(ix,iy,0)
         ylohi(0) = cry(ix,iy,0)
         xlohi(1) = crx(ix,iy,1)
         ylohi(1) = cry(ix,iy,1)
         call gpl_d (2, xlohi, ylohi)
*      ..top face
         xlohi(0) = crx(ix,iy,2)
         ylohi(0) = cry(ix,iy,2)
         xlohi(1) = crx(ix,iy,3)
         ylohi(1) = cry(ix,iy,3)
         call gpl_d (2, xlohi, ylohi)
        enddo
       enddo
       call pagend ()
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wgx1

      end
