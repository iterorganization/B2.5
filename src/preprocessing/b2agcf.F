*-----------------------------------------------------------------------
*.specification

      subroutine b2agcf (nx, ny, crx, cry, fpsi, ffbz)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     &  fpsi(-1:nx,-1:ny,0:3), ffbz(-1:nx,-1:ny,0:3)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AGCF computes the locations of the cell vertices, and computes
*     the flux functions at those locations.
*
*
*  3. description (see also routine b2cdca)
*
*     This subroutine obtains a description of the geometry and of the
*     magnetic field. The geometry is described by the arrays crx and
*     cry, which return horizontal (crx) and vertical (cry) coordinates
*     of the four vertices of all cells of the mesh. The magnetic field
*     is described by the arrays fpsi and ffbz, which return the
*     magnetic field potential and the field function at the same set of
*     locations.
*
*     The description of the geometry and the field is obtained by
*     a sequence of calls to the routine b2agsp. The content of b2agsp
*     depends upon the application at hand. The routine b2agcf is an
*     interface.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny - integer, input.
*     See the calling routine for description.
*
*     crx, cry - (-1:nx,-1:ny,0:3) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), crx(ix,iy,0:3) specify the
*     horizontal coordinates and cry(ix,iy,0:3) specify the vertical
*     coordinates (in the physical plane) of the four corners of the
*     (ix,iy) cell, mapped as follows:
*     (,,0): lower left corner,  (,,1): lower right corner,
*     (,,2): upper left corner,  (,,3): upper right corner.
*     One expects the following identities to hold away from the
*     boundaries:
*     crx(ix,iy,0)=crx(ix-1,iy,1)=crx(ix,iy-1,2)=crx(ix-1,iy-1,3),
*     cry(ix,iy,0)=cry(ix-1,iy,1)=cry(ix,iy-1,2)=cry(ix-1,iy-1,3),
*     However, there may be cuts in the geometry, in which case these
*     identities need not hold.
*
*     fpsi - (-1:nx,-1:ny,0:3) real*8 array, output.
*     For (ix,iy,i) in (-1:nx,-1:ny,0:3), fpsi(ix,iy,i) specifies the
*     potential for the magnetic field components (bx,by) at the
*     (ix,iy,i) position.
*     (In the case of straight symmetry, bx=dpsi/dy and by=-dpsi/dx.
*     In the case of toroidal or cylindrical symmetry,
*       bx=(dpsi/dy)/(2*pi*R) and by=-(dpsi/dx)/(2*pi*R),
*     where R is the distance to the axis of symmetry.)
*
*     ffbz - (-1:nx,-1:ny,0:3) real*8 array, output.
*     For (ix,iy,i) in (-1:nx,-1:ny,0:3), ffbz(ix,iy,i) specifies the
*     field function for the magnetic field component bz at the
*     (ix,iy,i) position.
*     (In the case of straight symmetry, bz=ffbz. In the case of
*     toroidal or cylindrical symmetry, bz=ffbz/(2*pi*R), where R is
*     the distance to the axis of symmetry.)
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix, iy
      real (kind=R8) ::
     *   prtst1, edg, ax, ay, t0
      parameter (prtst1=1.0_R8+1.0e-12_R8)
*     (prtst1.eq.1.0 will be used as test of the numerical precision)
*   ..procedures
      intrinsic min, max
      real (kind=R8) ::
     *   rratio, trunc
      external subini, subend, xertst, rratio
      external b2agsp
      trunc(t0,edg) = max(-edg,min(1+edg,t0))

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agcf')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..set edg
*     (An appropriate value for edg must be small compared to 1.0/nx
*     and 1.0/ny, and large compared to the machine precision. I prefer
*     not to invoke the precise machine parameter macheps here, because
*     that would make the results sensitive to the machine in a setting
*     where such a dependency is not really required. The test here
*     only distinguishes between single and double precision. I expect
*     that this code will normally be compiled for double precision.)
      if (prtst1.eq.1.0) then
*      (single precision)
       edg = 1.0e-4_R8
      else
*      (double precision)
       edg = 1.0e-7_R8
      endif

* ..compute physical coordinates
*   ..note
*     Each cell is indexed by a pair (ix,iy) in (-1:nx,-1:ny). For each
*     cell we compute the locations of the four cell corners; at all
*     these locations we also compute the magnetic potential, fpsi, and
*     the field function, ffbz. This is done through calls to the
*     application-dependent routine b2agsp.
*     The principal input argument to b2agsp is the specification
*     (ax,ay) of a point in the unit square, or just outside the unit
*     square in the case of the boundary cells. (ax,ay) are the
*     coordinates that make the mesh appear uniform on the unit square.
*     The routine b2agsp implements a transformation of this uniform
*     mesh to the actual mesh in the physical plane.
*     In many cases the transformation will have branch points,
*     corresponding to an x-point in the real*8 geometry. It is intended
*     that each branch point coincides with a cell vertex, and that
*     branch cuts follow cell boundaries. The arguments (kx,ky) passed
*     to b2agsp are intended to resolve the ambiguity in the case that
*     (ax,ay) lies on a cut. If (ax,ay) does not lie on a cut, then
*     (kx,ky) should be immaterial. If (ax,ay) does lie on a cut, then
*     (ax,ay) must be interpreted as (ax+kx*eps,ay+ky*eps), where eps
*     is a positive infinitesimal.
*     Legitimate values for kx and ky are -1, 0, and 1.
*     b2agsp returns the physical coordinates (crx,cry), and the values
*     of the magnetic field potential, fpsi, and of the magnetic field
*     function, ffbz, at that point.
*   ..loop over cells
      do iy = -1, ny
       do ix = -1, nx
*     ..compute quantities at cell corners
*       (corner 0)
        if (nx.eq.0) then
         ax = rratio(ix+1,2)
        else
         ax = trunc(rratio(ix,nx),edg)
        endif
        if (ny.eq.0) then
         ay = rratio(iy+1,2)
        else
         ay = trunc(rratio(iy,ny),edg)
        endif
        call b2agsp (ax, ay, 1, 1, crx(ix,iy,0), cry(ix,iy,0),
     &    fpsi(ix,iy,0), ffbz(ix,iy,0))
*       (corner 1)
        if (nx.eq.0) then
         ax = rratio(ix+2,2)
        else
         ax = trunc(rratio(ix+1,nx),edg)
        endif
        if (ny.eq.0) then
         ay = rratio(iy+1,2)
        else
         ay = trunc(rratio(iy,ny),edg)
        endif
        call b2agsp (ax, ay, -1, 1, crx(ix,iy,1), cry(ix,iy,1),
     &    fpsi(ix,iy,1), ffbz(ix,iy,1))
*       (corner 2)
        if (nx.eq.0) then
         ax = rratio(ix+1,2)
        else
         ax = trunc(rratio(ix,nx),edg)
        endif
        if (ny.eq.0) then
         ay = rratio(iy+2,2)
        else
         ay = trunc(rratio(iy+1,ny),edg)
        endif
        call b2agsp (ax, ay, 1, -1, crx(ix,iy,2), cry(ix,iy,2),
     &    fpsi(ix,iy,2), ffbz(ix,iy,2))
*       (corner 3)
        if (nx.eq.0) then
         ax = rratio(ix+2,2)
        else
         ax = trunc(rratio(ix+1,nx),edg)
        endif
        if (ny.eq.0) then
         ay = rratio(iy+2,2)
        else
         ay = trunc(rratio(iy+1,ny),edg)
        endif
        call b2agsp (ax, ay, -1, -1, crx(ix,iy,3), cry(ix,iy,3),
     &    fpsi(ix,iy,3), ffbz(ix,iy,3))
       enddo
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agcf

      end subroutine b2agcf
