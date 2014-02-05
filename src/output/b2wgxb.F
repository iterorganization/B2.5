*-----------------------------------------------------------------------
*.specification

      subroutine b2wgxb (nx, ny, crx, cry)
      use b2mod_types
      use b2mod_indirect
      use b2mod_b2cmgs
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
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
*     B2WGXB plots an outline of the b2 mesh on the physical domain.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny - integer, input.
*     See the calling routine for description.
*
*     crx, cry - (-1:nx,-1:ny,0:3) real*8 array, input.
*     See the calling routine for description.
*
*     (parameters in common:)
*
*     /b2cmgs/ - input.
*     /b2cmgs/ specifies the singular regions in the geometry.
*      nlrmx - integer constant.
*     nlrmx is an upper bound for nlreg, and is used in dimension
*     specifications.
*      nlreg - integer.
*     nlreg specifies the number of singular regions. It will hold that
*     0.le.nlreg.le.nlrmx.
*      nlxlo, nlxhi, nlylo, nlyhi - (0:nlrmx-1) integer array.
*     For 0.le.i.lt.nlreg, nlxlo(i), nlxhi(i), nlylo(i) and nlyhi(i) are
*     the bounds of the i-th singular region. Mesh locations (ix,iy)
*     in (nlxlo(i):nlxhi(i),nlylo(i):nlyhi(i)) make up the region, which
*     will be included in (-1:nx,-1:ny).
*     For nlreg.le.i.lt.nlrmx, nlxlo(i), nlxhi(i), nlylo(i) and nlyhi(i)
*     are unspecified.
*      nlloc - (0:nlrmx-1) integer array.
*     For 0.le.i.lt.4, nlloc(i) specifies the class of quantities for
*     which the i-th region is special, as follows:
*       nlloc(i) = 0: cell centers, indices (-1:nx,-1:ny).
*       nlloc(i) = 1: cell x-faces, indices (0:nx,-1:ny).
*       nlloc(i) = 2: cell y-faces, indices (-1:nx,0:ny).
*       nlloc(i) = 3: cell corners, indices (0:nx,0:ny).
*     For nlreg.le.i.lt.nlrmx, nlloc(i) is unspecified.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ireg, ix, iy
      real (kind=R8) ::
     *   xlohi(0:1), ylohi(0:1)
*   ..procedures
      external subini, subend, xertst, gpl_d

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wgxb')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

* ..draw grid outline
*   ..south face
      do ix = -1, nx
       xlohi(0) = crx(ix,-1,0)
       ylohi(0) = cry(ix,-1,0)
       xlohi(1) = crx(ix,-1,1)
       ylohi(1) = cry(ix,-1,1)
       call gpl_d (2, xlohi, ylohi)
      enddo
*   ..north face
      do ix = -1, nx
       xlohi(0) = crx(ix,ny,2)
       ylohi(0) = cry(ix,ny,2)
       xlohi(1) = crx(ix,ny,3)
       ylohi(1) = cry(ix,ny,3)
       call gpl_d (2, xlohi, ylohi)
      enddo
*   ..west face
      do iy = -1, ny
       xlohi(0) = crx(-1,iy,0)
       ylohi(0) = cry(-1,iy,0)
       xlohi(1) = crx(-1,iy,2)
       ylohi(1) = cry(-1,iy,2)
       call gpl_d (2, xlohi, ylohi)
      enddo
*   ..east face
      do iy = -1, ny
       xlohi(0) = crx(nx,iy,1)
       ylohi(0) = cry(nx,iy,1)
       xlohi(1) = crx(nx,iy,3)
       ylohi(1) = cry(nx,iy,3)
       call gpl_d (2, xlohi, ylohi)
      enddo

* ..draw faces of singular regions
      do ireg = 0, nlreg-1
       if (nlloc(ireg).eq.1) then
*     ..test region boundaries
        call xertst (0.le.nlxlo(ireg).and.nlxhi(ireg).le.nx.and.
     &    -1.le.nlylo(ireg).and.nlyhi(ireg).le.ny,
     &    'faulty region bounds in /b2cmgs/ for nlloc.eq.1')
*     ..draw x-faces
        do iy = nlylo(ireg), nlyhi(ireg)
         do ix = nlxlo(ireg), nlxhi(ireg)
          xlohi(0) = crx(leftix(ix,iy),leftiy(ix,iy),1)
          ylohi(0) = cry(leftix(ix,iy),leftiy(ix,iy),1)
          xlohi(1) = crx(leftix(ix,iy),leftiy(ix,iy),3)
          ylohi(1) = cry(leftix(ix,iy),leftiy(ix,iy),3)
          call gpl_d (2, xlohi, ylohi)
          xlohi(0) = crx(ix,iy,0)
          ylohi(0) = cry(ix,iy,0)
          xlohi(1) = crx(ix,iy,2)
          ylohi(1) = cry(ix,iy,2)
          call gpl_d (2, xlohi, ylohi)
         enddo
        enddo
       else if (nlloc(ireg).eq.2) then
*     ..test region boundaries
        call xertst (-1.le.nlxlo(ireg).and.nlxhi(ireg).le.nx.and.
     &    0.le.nlylo(ireg).and.nlyhi(ireg).le.ny,
     &    'faulty region bounds in /b2cmgs/ for nlloc.eq.2')
*     ..draw y-faces
        do iy = nlylo(ireg), nlyhi(ireg)
         do ix = nlxlo(ireg), nlxhi(ireg)
          xlohi(0) = crx(bottomix(ix,iy),bottomiy(ix,iy),2)
          ylohi(0) = cry(bottomix(ix,iy),bottomiy(ix,iy),2)
          xlohi(1) = crx(bottomix(ix,iy),bottomiy(ix,iy),3)
          ylohi(1) = cry(bottomix(ix,iy),bottomiy(ix,iy),3)
          call gpl_d (2, xlohi, ylohi)
          xlohi(0) = crx(ix,iy,0)
          ylohi(0) = cry(ix,iy,0)
          xlohi(1) = crx(ix,iy,1)
          ylohi(1) = cry(ix,iy,1)
          call gpl_d (2, xlohi, ylohi)
         enddo
        enddo
       endif
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wgxb

      end
