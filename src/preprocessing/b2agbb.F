*-----------------------------------------------------------------------
*.specification

      subroutine b2agbb (nx, ny, vol, hx, hy, qz, fpsi, ffbz, bb)
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  qz(-1:nx,-1:ny,0:1), fpsi(-1:nx,-1:ny,0:3),
     &  ffbz(-1:nx,-1:ny,0:3)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   bb(-1:nx,-1:ny,0:3)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AGBB computes the magnetic field from the flux functions.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the three components of the magnetic field
*     and the total magnetic field strength, all at cell centers. The
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
*     vol, hx, hy, qz - real*8 array, input.
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
      integer ix, iy
      real (kind=R8) ::
     *   t0
*   ..procedures
      intrinsic sqrt
      external xertst

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agbb')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

* ..compute magnetic field
*   ..loop over cells
      do iy = -1, ny
       do ix = -1, nx
*     ..compute magnetic field at cell center
        t0 = vol(ix,iy)/(hx(ix,iy)*hy(ix,iy)*qz(ix,iy,1))
        bb(ix,iy,0) =
     &    (fpsi(ix,iy,2)-fpsi(ix,iy,0)+fpsi(ix,iy,3)-fpsi(ix,iy,1))*
     &    hx(ix,iy)/(2.0_R8*vol(ix,iy))
        bb(ix,iy,1) = 0.0_R8
*!!     (Note: Consider to introduce nonzero bb(,,1), as below)
*$$$    bb(ix,iy,1) =
*$$$ &    -(fpsi(ix,iy,1)-fpsi(ix,iy,0)+fpsi(ix,iy,3)-fpsi(ix,iy,2))*
*$$$ &    hy(ix,iy)*qz(ix,iy,1)/(2.0_R8*vol(ix,iy))
        bb(ix,iy,2) = (ffbz(ix,iy,0)+ffbz(ix,iy,1)+ffbz(ix,iy,2)+
     &    ffbz(ix,iy,3))/(4.0_R8*t0)
        bb(ix,iy,3) =
     &    sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2+bb(ix,iy,2)**2)
       enddo
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agbb

      end
