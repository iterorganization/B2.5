*-----------------------------------------------------------------------
*.specification

      subroutine b2cdcv ()
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2CDCV contains a description of some of the variables used in
*     the b2 code. No computation is performed in this routine.
*
*-----------------------------------------------------------------------
*.text
*
*  0. Introduction
*
*     ...
*
*
*  1. Principal code variables
*
*   ..variables describing the geometry
*
*     nx, ny - integer.
*     nx and ny specify the number of interior cells along the first
*     and the second coordinate, respectively. The total number of
*     cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
*     It will hold that 0.le.nx and 0.le.ny.
*
*     crx, cry - (-1:nx,-1:ny,0:3) real*8 array.
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
*     fpsi - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy,i) in (-1:nx,-1:ny,0:3), fpsi(ix,iy,i) specifies the
*     potential for the magnetic field components (bx,by) at the
*     (ix,iy,i) position.
*     (In the case of straight symmetry, bx=dpsi/dy and by=-dpsi/dx.
*     In the case of toroidal or cylindrical symmetry,
*       bx=(dpsi/dy)/(2*pi*R) and by=-(dpsi/dx)/(2*pi*R),
*     where R is the distance to the axis of symmetry.)
*
*     ffbz - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy,i) in (-1:nx,-1:ny,0:3), ffbz(ix,iy,i) specifies the
*     field function for the magnetic field component bz at the
*     (ix,iy,i) position.
*     (In the case of straight symmetry, bz=ffbz. In the case of
*     toroidal or cylindrical symmetry, bz=ffbz/(2*pi*R), where R is
*     the distance to the axis of symmetry.)
*
*     bb - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), bb(ix,iy,0:3) specifies the
*     magnetic field at the center of the (ix,iy) cell.
*     bb(ix,iy,0:2) are the (x,y,z)-components and bb(ix,iy,3) is
*     the absolute magnetic field strength.
*     (bb(,,3)=sqrt(bb(,,0)**2+bb(,,1)**2+bb(,,2)**2.)
*     It will hold that 0.lt.bb(,,3).
*
*     vol - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), vol(ix,iy) specifies the volume of
*     the (ix,iy) cell. (In the case of toroidal symmetry this volume
*     includes a factor 2*pi*crx, while in the case of cylindrical
*     symmetry it includes a factor 2*pi*cry.)
*     It will hold that 0.lt.vol(,).
*
*     hx, hy, hz - (-1:nx,-1:ny) real*8 array.                          !srv 09.06.08
*     For (ix,iy) in (-1:nx,-1:ny), hx(ix,iy) and hy(ix,iy) specify the
*     x-diameter, the y-diameter and z-diameter respectively, of the 
*     (ix,iy) cell.
*     It will hold that 0.lt.hx(,) and 0.lt.hy(,).
*
*     qz - (-1:nx,-1:ny,0:1) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), qz(ix,iy,0) specifies sin(t) and
*     qz(ix,iy,1) specifies cos(t), where t is the complementary angle
*     between the x-direction and the y-direction on the (ix,iy) cell.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. (qz(ix,iy,0) is the inner product of the x and y unit
*     vectors.)
*     It will hold that -1.lt.qz(,,0).lt.1 and 0.lt.qz(,,1).le.1.
*
*     qc - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qc(ix,iy) specifies cos(t), where 
*     t is the complementary angle between the x-direction of the 
*     (ix,iy) cell and its left face.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. 
*     It will hold that 0.lt.qc(,).le.1.
*
*     gs - (-1:nx,-1:ny,0:2) real*8 array.
*     (Let gsx(,)=gs(,,0), gsy(,)=gs(,,1) and gsz(,)=gs(,,2).)
*     For (ix,iy) in (0:nx,-1:ny), gsx(ix,iy) specifies the area of the
*     face between the (ix-1,iy) and the (ix,iy) cell. For ix.eq.-1,
*     gsx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), gsy(ix,iy) specifies the area of the
*     face between the (ix,iy-1) and the (ix,iy) cell. For iy.eq.-1,
*     gsy(-1:nx,iy) will hold 0.
*     For (ix,iy) in (-1:nx,-1:ny), gsz(ix,iy) specifies the area of the
*     (ix,iy) cell normal to the third (ignorable) coordinate.
*
*     isymm - integer.
*     isymm specifies the type of symmetry of the geometry: isymm.eq.0
*     implies a slab geometry, isymm.eq.1.or.isymm.eq.2 imply rotational 
*     symmetry about the crx=0 axis, and isymm.eq.3.or.isymm.eq.4 
*     indicate rotational symmetry about the cry=0 axis. 
*     Other values are not allowed. 
*     isymm.eq.1.or.isymm.eq.3 mean a toroidal geometry case, i.e. 
*     with a toroidal (out-of-plane) magnetic field component.
*     isymm.eq.2.or.isymm.eq.4 mean a cylindrical geometry case, i.e. 
*     no toroidal magnetic field component.
*
*     periodic_bc - integer.
*     periodic_bc indicates whether the geometry closes upon itself (i.e.
*     whether periodic boundary conditions are necessary). The default value
*     is 0. periodic_bc.eq.1 implies that periodic BC are being used. This
*     case is used for the simulation of diverted stellarator islands or
*     limiters.
*
*     /b2cmfs/ - common block.
*     /b2cmfs/ is reserved for numerical mhd equilibrium data.
*      isymm - integer.
*     See above for description. 
*
*     /b2cmgs/ - common block.
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
*   ..auxiliary variables describing the geometry
*
*     pbs - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let pbsx(,)=pbs(,,0) and pbsy(,)=pbs(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), pbsx(ix,iy) specifies the product
*     (bx/bb)*sx and pbsy(ix,iy) specifies (by/bb)*sy, or (magnetic 
*     field pitch)*(surface area), namely the parallel contact area, 
*     on the (ix,iy) cell left and bottom faces respectively.
*     In the case where the grid covers the full vacuum vessel, 
*     pbsy also contains a contribution from the narrowing of the
*     flux tube when impacting the wall structures.  
*     See routine b2xgbs for the form of interpolation employed.
*
*     pbshz - (-1:nx,-1:ny,0:1) real*8 array.                       !srv 09.06.08
*     (Let pbshzx(,)=pbshz(,,0) and pbshzy(,)=pbshz(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), pbshzx(ix,iy) specifies the product
*     hz*(bx/bb)*sx and pbshzy(ix,iy) specifies hz*(by/bb)*sy, or (magnetic 
*     field pitch)*(surface area)*hz
*
*     bzb - (-1:nx,-1:ny) real*8 array contains Bz/B**2 - Bz/<B**2> !srv 09.06.08
*     on the cell centre.
*
*     OnedBsq - (-1:nx,-1:ny) real*8 array.                       !sv
*     For (ix,iy) in (-1:nx,-1:ny), OnedBsq(ix,iy) specifies the
*     1/bb(,,3)**2 on the cell centre
*
*   ..basic physical parameters
*
*     /b2cmpa/ - common block
*     /b2cmpa/ specifies basic physics parameters.
*      me, mp, ev, qe, eps0, mu0 - real*8 constant.
*     These parameters are included with /b2cmpa/.
*       me = electron mass (kg)
*       mp = proton mass (kg)
*       ev = 1 electron-Volt (Joule)
*       qe = elementary charge (Coulomb)
*       eps0 = vacuum dielectric permittivity (SI)
*       mu0 = vacuum magnetic permeability (SI)
*      nsdecl - integer constant.
*     nsdecl is an upper bound for the number of species (ns), and is
*     used in dimension specifications.
*      zamin, zamax, zn, am - (0:nsdecl-1) real*8 array.
*     For (is) in (0:ns-1), zn(is) specifies the nuclear charge of
*     species (is), am(is) specifies its atomic mass, and zamin(is) and
*     zamax(is) specify its atomic charge range. All are dimensionless.
*     It will hold that 0.le.zamin(is).le.zamax(is), 0.lt.zn(is) and
*     0.lt.am(is). It is required that neutral species not be bundled
*     with charged species.
*     For ns.le.is.lt.nsdecl, zamin(is), zamax(is), zn(is), am(is) 
*     are undefined.
*
*   ..variables describing the plasma state
*
*     ns - integer.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*     na - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), na(ix,iy,is) specifies
*     the density of atomic species (is) on the (ix,iy) cell.
*     It will hold that 0.lt.na(,,).
*
*     ua - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), ua(ix,iy,is) specifies
*     the parallel velocity of atomic species (is) on the (ix,iy) cell.
*
*     te - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), te(ix,iy) specifies the electron
*     temperature on the (ix,iy) cell. It will hold that 0.lt.te(,).
*
*     ti - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), ti(ix,iy) specifies the all atom
*     temperature on the (ix,iy) cell. It will hold that 0.lt.ti(,).
*
*     po - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), po(ix,iy) specifies the electric
*     potential on the (ix,iy) cell.
*
*     fna - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let fnax(,,)=fna(,,0,) and fnay(,,)=fna(,,1,).)
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), fnax(ix,iy,is) specifies
*     the flux of atoms of species (is) through the face between the
*     (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fnax(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), fnay(ix,iy,is) specifies
*     the flux of atoms of species (is) through the face between the
*     (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fnay(-1:nx,iy,0:ns-1) will hold 0.
*
*     fna_mdf - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let fnax_mdf(,,)=fna_mdf(,,0,) and fnay_mdf(,,)=fna_mdf(,,1,).)
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), fnax_mdf(ix,iy,is) specifies
*     the modified flux of atoms of species (is) through the face between
*     the (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fnax_mdf(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), fnay_mdf(ix,iy,is) specifies
*     the modified flux of atoms of species (is) through the face between
*     the (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fnay_mdf(-1:nx,iy,0:ns-1) will hold 0.
*     fna_mdf is equal fna for runs without drifts. 
*
*     fna_fcor - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let fna_fcorx(,,)=fna_fcor(,,0,) and fna_fcory(,,)=fna_fcor(,,1,).)
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), fna_fcorx(ix,iy,is) specifies
*     the flux of atoms of species (is) for calculating momentum transport
*     through the face between the (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fna_fcorx(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), fna_fcory(ix,iy,is) specifies
*     the flux of atoms of species (is) for calculating momentum transport
*     through the face between the(ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fna_fcory(-1:nx,iy,0:ns-1) will hold 0.
*     fna_fcor takes into account Coriolis force.
*     fna_fcor is equal fna for runs without drifts. 
*
*     fhe - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhex(,)=fhe(,,0) and fhey(,)=fhe(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhex(ix,iy) specifies the electron
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For ix.eq.-1, fhex(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhey(ix,iy) specifies the electron
*     heat flux through the face between the (ix,iy-1) and (ix,iy)
*     cells. For iy.eq.-1, fhey(-1:nx,iy) will hold 0.
*
*     fhe_mdf - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhex_mdf(,)=fhe_mdf(,,0) and fhey_mdf(,)=fhe_mdf(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhex_mdf(ix,iy) specifies the modified
*     electron heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For ix.eq.-1, fhex_mdf(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhey_mdf(ix,iy) specifies the modified
*     electron heat flux through the face between the (ix,iy-1) and (ix,iy)
*     cells. For iy.eq.-1, fhey_mdf(-1:nx,iy) will hold 0.
*     fhe_mdf is equal fhe for runs without drifts. 
*
*     fhi - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhix(,)=fhi(,,0) and fhiy(,)=fhi(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhix(ix,iy) specifies the all atom
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For ix.eq.-1, fhix(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhiy(ix,iy) specifies the all atom
*     heat flux through the face between the (ix,iy-1) and (ix,iy)
*     cells. For iy.eq.-1, fhiy(-1:nx,iy) will hold 0.
*
*     fhi_mdf - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhix_mdf(,)=fhi_mdf(,,0) and fhiy_mdf(,)=fhi_mdf(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhix_mdf(ix,iy) specifies the 
*     modified all atom heat flux through the face between the (ix-1,iy)
*      and (ix,iy) cells. For ix.eq.-1, fhix_mdf(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhiy_mdf(ix,iy) specifies the modified
*     all atom heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For iy.eq.-1, fhiy_mdf(-1:nx,iy) will hold 0.
*
*     fch - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchx(,)=fch(,,0) and fchy(,)=fch(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchx(ix,iy) specifies the electric
*     current through the face between the (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fchx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchy(ix,iy) specifies the electric
*     current through the face between the (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fchy(-1:nx,iy) will hold 0.
*
*   ..auxiliary variables describing the plasma state
*
*     nb - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), nb(ix,iy) specifies the density of
*     one selected atomic species on the (ix,iy) cell. If that species
*     has index (is) then nb(,).eq.na(,,is).
*
*     ub - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), ub(ix,iy) specifies the parallel
*     velocity of one selected atomic species on the (ix,iy) cell.
*
*     fnb - (-1:nx,-1:ny,0:1) real*8 array.
*     fnb specifies the flux of one selected atomic species.
*     If that species has index (is) then fnb(,,)=fna(,,,is).
*     precisely: (Let fnbx(,)=fnb(,,0) and fnby(,)=fnb(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fnbx(ix,iy) specifies the flux of
*     the selected atomic species through the face between the
*     (ix-1,iy) and the (ix,iy) cell.
*     For ix.eq.-1, fnbx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fnby(ix,iy) specifies the flux of
*     the selected atomic species through the face between the
*     (ix,iy-1) and the (ix,iy) cell.
*     For iy.eq.-1, fnby(-1:nx,iy) will hold 0.
*     fnb_mdf is equal fnb for runs without drifts. 
*
*     fnb_mdf - (-1:nx,-1:ny,0:1) real*8 array.
*     fnb_mdf specifies the flux of one selected atomic species.
*     If that species has index (is) then fnb_mdf(,,)=fna_mdf(,,,is).
*     precisely: (Let fnbx_mdf(,)=fnb_mdf(,,0) and fnby_mdf(,)=fnb_mdf(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fnbx_mdf(ix,iy) specifies the modified
*     flux of the selected atomic species through the face between the
*     (ix-1,iy) and the (ix,iy) cell.
*     For ix.eq.-1, fnbx_mdf(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fnby_mdf(ix,iy) specifies the modified
*     flux of the selected atomic species through the face between the
*     (ix,iy-1) and the (ix,iy) cell.
*     For iy.eq.-1, fnby_mdf(-1:nx,iy) will hold 0.
*
*     fnb_fcor - (-1:nx,-1:ny,0:1) real*8 array.
*     fnb_fcor specifies the flux of one selected atomic species 
*     for calculating momentum transport.
*     If that species has index (is) then fnb_fcor(,,)=fna_fcor(,,,is).
*     precisely: (Let fnbx_fcor(,)=fnb_fcor(,,0) and fnby_fcor(,)=fnb_fcor(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fnb_fcorx(ix,iy) specifies the flux of
*     the selected atomic species through the face between the
*     (ix-1,iy) and the (ix,iy) cell.
*     For ix.eq.-1, fnb_fcorx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fnb_fcory(ix,iy) specifies the flux of
*     the selected atomic species through the face between the
*     (ix,iy-1) and the (ix,iy) cell.
*     For iy.eq.-1, fnb_fcory(-1:nx,iy) will hold 0.
*     fnb_fcor is equal fnb for runs without drifts. 
*
*     ni - (-1:nx,-1:ny,0:1) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), ni(ix,iy,0) specifies the total atom
*     density on the (ix,iy) cell. It will hold that 0.lt.ni(,0).
*     The defining expression for ni is
*       ni(,0) = (sum is : 0.le.is.lt.ns : na(,,is)).
*     For (ix,iy) in (-1:nx,-1:ny), ni(ix,iy,1) specifies the total ion
*     density on the (ix,iy) cell. It will hold that 0.lt.ni(,1).
*     The defining expression for ni is
*       ni(,1) = (sum is : 0.le.is.lt.ns : na(,,is)) when rza(,,is)>0.
*
*     ne - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), ne(ix,iy) specifies the electron
*     density on the (ix,iy) cell. It will hold that 0.lt.ne(,).
*     The defining expression for ne is
*       ne(,) = (sum is : 0.le.is.lt.ns : rza(,,is)*na(,,is)).
*
*     ne2 - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), ne2(ix,iy) specifies the atom
*     density weighted by rz2 on the (ix,iy) cell. It will hold that
*     0.lt.ne2(,).
*     The defining expression for ne2 is
*       ne2(,) = (sum is : 0.le.is.lt.ns : rz2(,,is)*na(,,is)).
*
*     pr - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), pr(ix,iy) specifies the total
*     pressure on the (ix,iy) cell. The defining expression for pr is
*       pr(,) = ne(,)*te(,) + (sum is :: na(,,is)*ti(,)).
*
*     pe - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), pe(ix,iy) specifies the electron
*     pressure on the (ix,iy) cell.
*     The defining expression for pe is  pe(,) = ne(,)*te(,).
*
*     ro - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), ro(ix,iy) specifies the mass
*     density on the (ix,iy) cell. The defining expression for ro is
*       ro(,) = (sum is : 0.le.is.lt.ns : am(is)*mp**na(,,is)).
*
*     ne2m - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), ne2m(ix,iy) specifies the atom
*     density weighted by rz2/am on the (ix,iy) cell. It will hold
*     that 0.lt.ne2m(,).
*     The defining expression for ne2m is
*       ne2m(,) = (sum is :: (rz2(,,is)/am(is))*na(,,is)).
*
*     fne - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fnex(,)=fne(,,0) and fney(,)=fne(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fnex(ix,iy) specifies the flux of
*     electrons through the face between the (ix-1,iy) and the (ix,iy)
*     cell. For ix.eq.-1, fnex(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fney(ix,iy) specifies the flux of
*     electrons through the face between the (ix,iy-1) and the (ix,iy)
*     cell. For iy.eq.-1, fney(-1:nx,iy) will hold 0.
*     The defining expression for fne is
*       fne(,,) = (sum is :: rza(,,is)*fna(,,,is))-(1/qe)*fch(,,)
*
*     fni - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fnix(,)=fni(,,0) and fniy(,)=fni(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fnix(ix,iy) specifies the flux of
*     atoms (all species) through the face between the (ix-1,iy) and
*     the (ix,iy) cell. For ix.eq.-1, fnix(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fniy(ix,iy) specifies the flux of
*     atoms (all species) through the face between the (ix,iy-1) and
*     the (ix,iy) cell. For iy.eq.-1, fniy(-1:nx,iy) will hold 0.
*     The defining expression for fni is
*       fni(,,) = (sum is :: fna(,,,is))
*
*     fchvispar - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchvisparx(,)=fch(,,0) and fchvispary(,)=fch(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchvisparx(ix,iy) specifies the
*     electric current driven by parallel viscosity through the face
*     between the (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fchx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchvispary(ix,iy) specifies the
*     electric current driven by parallel viscosity through the face
*     between the (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fchvispary(-1:nx,iy) will hold 0.
*
*     fchvisq - (-1:nx,-1:ny,0:1) real*8 array.
*     fchvisq specifies the electric current produced by the components
*     of the viscosity tensor which is connected with the heat fluxes.
*     (Let fchvisqx(,)=fch(,,0) and fchvisqy(,)=fch(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchvisqx(ix,iy) specifies the current
*     through the face between the (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fchx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchvisqy(ix,iy) specifies the current
*     through the face between the (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fchvisqy(-1:nx,iy) will hold 0.
*
*     fchinert - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchinertx(,)=fch(,,0) and fchinerty(,)=fch(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchinertx(ix,iy) specifies the
*     contributions from inertia and gyroviscosity current through
*     the face between the (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fchx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchinerty(ix,iy) specifies the
*     contributions from inertia and gyroviscosity current through
*     the face between the (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fchinerty(-1:nx,iy) will hold 0.
*
*     fchanml - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchanmlx(,)=fch(,,0) and fchanmly(,)=fch(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchanmlx(ix,iy) specifies the
*     anomalous current through the face between the (ix-1,iy)
*     and (ix,iy) cells.
*     For ix.eq.-1, fchx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchanmly(ix,iy) specifies the
*     amomalous current through the face between the (ix,iy-1)
*     and (ix,iy) cells.
*     For iy.eq.-1, fchanmly(-1:nx,iy) will hold 0.
*
*     fchdia - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchdiax(,)=fch(,,0) and fchdiay(,)=fch(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchdiax(ix,iy) specifies the
*     modified diamagnetic current through the face between 
*     the (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fchx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchdiay(ix,iy) specifies the
*     modified diamagnetic current through the face between
*     the (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fchdiay(-1:nx,iy) will hold 0.
*
*     fchin - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchinx(,)=fch(,,0) and fchiny(,)=fch(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchinx(ix,iy) specifies the
*     contribution from ion-neutral friction through the face between 
*     the (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fchx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchiny(ix,iy) specifies the
*     contribution from ion-neutral friction through the face between
*     the (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fchiny(-1:nx,iy) will hold 0.
*
*     fch_p - (-1:nx,-1:ny,0:1) real*8 array.
*     For (ix,iy) in (0:nx,-1:ny), fch_p(ix,iy,0) specifies the parallel 
*     electric current along the magnetic field through the face between 
*     the (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fch_p(ix,-1:ny,0) will hold 0. 
*     fch_p(-1:nx,-1:ny,1) = 0
*
*     fchvisper - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchvisperx(,)=fch(,,0) and fchvispery(,)=fch(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), fchvisperx(ix,iy) = 0.
*     For (ix,iy) in (-1:nx, 0:ny), fchvispery(ix,iy) specifies the
*     electric current connected with contribution from
*     the perpendicular viscositythrough the face between the (ix,iy-1) 
*     and (ix,iy) cells.
*     For iy.eq.-1, fchvispery(-1:nx,iy) will hold 0.
*
*     uadia - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let uadiax(,,)=uadia(,,0,) and uadiay(,)=uadia(,,1,).)
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), uadiax(ix,iy,is) specifies
*     the total effective drift velocity of species (is) in cell (ix,iy)
*     in the diamagnetic direction (i.e. within the flux surface but 
*     perpendicular to the magnetic field).
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), uadiay(ix,iy,is) specifies
*     the total effective drift velocity of species (is) in cell (ix,iy)
*     in the radial direction.
*     uadia refers to cell center if 'b2tfnb_drift_style' has '0'
*     uadia refers to cell faces if 'b2tfnb_drift_style' has '1'
*
*     vadia - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let vadiax(,,)=vadia(,,0,) and vadiay(,)=vadia(,,1,).)
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vadiax(ix,iy,is) specifies
*     the effective diamagnetic drift velocity of species (is) in cell (ix,iy)
*     in the diamagnetic direction (i.e. within the flux surface but 
*     perpendicular to the magnetic field).
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vadiay(ix,iy,is) specifies
*     the effective diamagnetic drift velocity of species (is) in cell (ix,iy)
*     in the radial direction.
*     vadia refers to cell center if 'b2tfnb_drift_style' has '0'
*     vadia refers to cell faces if 'b2tfnb_drift_style' has '1'
*
*     vedia - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let vediax(,)=vedia(,,0) and vediay(,)=vedia(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), vediax(ix,iy) specifies
*     the effective diamagnetic drift velocity of electrons in cell (ix,iy)
*     in the diamagnetic direction (i.e. within the flux surface but 
*     perpendicular to the magnetic field).
*     For (ix,iy) in (-1:nx,-1:ny), vediay(ix,iy) specifies
*     the effective diamagnetic drift velocity of electrons in cell (ix,iy)
*     in the radial direction.
*
*     wadia - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let wadiax(,,)=wadia(,,0,) and wadiay(,)=wadia(,,1,).)
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), wadiax(ix,iy,is) specifies
*     the total diamagnetic drift velocity of species (is) in cell (ix,iy)
*     in the diamagnetic direction (i.e. within the flux surface but 
*     perpendicular to the magnetic field).
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), wadiay(ix,iy,is) specifies
*     the total diamagnetic drift velocity of species (is) in cell (ix,iy)
*     in the radial direction.
*     wadia refers to cell center if 'b2tfnb_drift_style' has '0'
*     wadia refers to cell faces if 'b2tfnb_drift_style' has '1'
*
*     vaecrb - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let vaecrbx(,,)=vaecrb(,,0,) and vaecrby(,)=vaecrb(,,1,).)
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vaecrbx(ix,iy,is) specifies
*     the ExB drift velocity of species (is) in cell (ix,iy)
*     in the diamagnetic direction (i.e. within the flux surface but 
*     perpendicular to the magnetic field).
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vaecrby(ix,iy,is) specifies
*     the ExB drift velocity of species (is) in cell (ix,iy)
*     in the radial direction.
*     vaecrb refers to cell center if 'b2tfnb_drift_style' has '0'
*     vaecrb refers to cell faces if 'b2tfnb_drift_style' has '1'
*
*     kinrgy - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), kinrgy(ix,iy,is) specifies
*     the kinetic energy of a particle of species (is) in cell (ix,iy). 
*     The velocity used includes the parallel velocity as well
*     as the total diamagnetic and ExB drift velocities in both the radial
*     and diamagnetic directions.
*
*     veecrb - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let veecrbx(,)=veecrb(,,0) and veecrby(,)=veecrb(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), veecrbx(ix,iy) specifies
*     the ExB drift velocity of electrons in cell (ix,iy)
*     in the diamagnetic direction (i.e. within the flux surface but 
*     perpendicular to the magnetic field).
*     For (ix,iy) in (-1:nx,-1:ny), veecrby(ix,iy) specifies
*     the ExB drift velocity of electrons in cell (ix,iy)
*     in the radial direction.
*
*
*   ..residuals and corrections
*
*     resco - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), resco(ix,iy,is) specifies
*     the residual of the continuity equation on the (ix,iy) cell for
*     species (is).
*
*     resmo - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), resmo(ix,iy,is) specifies
*     the residual of the parallel momentum equation on the (ix,iy) 
*     cell for species (is).
*
*     resmb - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), resmb(ix,iy) specifies the residual
*     of the parallel momentum equation on the (ix,iy) cell for the 
*     selected species. If that species has index (is) then 
*     resmb(,).eq.resmo(,,is).
*
*     resmt - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), resmt(ix,iy) specifies the residual
*     of the total parallel momentum equation on the (ix,iy) cell.
*
*     reshe - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), reshe(ix,iy) specifies the residual
*     of the electron heat balance equation on the (ix,iy) cell.
*
*     reshi - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), reshe(ix,iy) specifies the residual
*     of the ion heat balance equation on the (ix,iy) cell.
*
*     resht - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), resht(ix,iy) specifies the residual
*     of the total heat balance equation on the (ix,iy) cell.
*
*     respo - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), respo(ix,iy) specifies the
*     residual of the electric potential equation on the (ix,iy) cell.
*
*     corpa - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), corpa(ix,iy,is) specifies
*     the last computed correction to the pressure of species (is) on
*     the (ix,iy) cell. corpa is computed by solving an equation that
*     has resco as its right hand side. It is not necessarily the same
*     as the pressure correction over one complete iteration.
*
*     corua - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), corua(ix,iy,is) specifies
*     the last computed correction to the parallel velocity ua of
*     species (is) on the (ix,iy) cell.
*     corua is computed by solving an equation that has resmo as its
*     right hand side. It is not necessarily the same as the change in
*     the parallel velocity over one complete iteration.
*
*     corut - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (0:nx,-1:ny), corut(ix,iy) specifies the last
*     computed correction to the parallel velocity ua of all species
*     simultaneously on the (ix,iy) cell.
*     corut is computed by solving an equation that has resmt as its
*     right hand side. It is not necessarily the same as the change in
*     the parallel velocity over one complete iteration.
*
*     corte - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), corte(ix,iy) specifies the last
*     computed correction to the electron temperature on the (ix,iy)
*     cell. corte is computed by solving an equation that has reshe as
*     its right hand side. It is not necessarily the same as the change
*     in the electron temperature over one complete iteration.
*
*     corti - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), corti(ix,iy) specifies the last
*     computed correction to the atom temperature on the (ix,iy) cell.
*     corti is computed by solving an equation that has reshi as its
*     right hand side. It is not necessarily the same as the change in
*     the atom temperature over one complete iteration.
*
*     cortt - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), cortt(ix,iy) specifies the last
*     computed correction to the combined temperature on the (ix,iy)
*     cell. cortt is computed by solving an equation that has resht as
*     its right hand side. It is not necessarily the same as the change
*     in the electron or atom temperature over one complete iteration.
*
*     corpo - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), corpo(ix,iy) specifies the last
*     computed correction to the electric potential on the (ix,iy)
*     cell. corpo is computed by solving an equation that has respo as
*     its right hand side. It is not necessarily the same as the change
*     in the electric potential over one complete iteration.
*
*     aresco - (0:ns-1) real*8 array.
*     For (is) in (0:ns-1), aresco(is) specifies the scaled norm of the
*     residual of the continuity equation for species (is).
*
*     aresmo - (0:ns-1) real*8 array.
*     For (is) in (0:ns-1), aresmo(is) specifies the scaled norm of the
*     residual of the parallel momentum equation for species (is).
*
*     aresmt - real.
*     aresmt specifies the scaled norm of the residual of the total
*     parallel momentum equation.
*
*     areshe - real.
*     areshe specifies the scaled norm of the residual of the electron
*     heat balance equation.
*
*     areshi - real.
*     areshi specifies the scaled norm of the residual of the all atom
*     heat balance equation.
*
*     arespo - real.
*     arespo specifies the scaled norm of the residual of the electric
*     potential equation.
*
*     acorpa - (0:ns-1) real*8 array.
*     For (is) in (0:ns-1), acorpa(is) specifies the norm of the last
*     computed correction for the continuity equation for species (is).
*
*     acorua - (0:ns-1) real*8 array.
*     For (is) in (0:ns-1), acorua(is) specifies the norm of the last
*     computed correction for the momentum balance equation for
*     species (is).
*
*     acorut - real.
*     acorut specifies the norm of the last computed correction for the
*     parallel momentum balance equation for all species simultaneously.
*
*     acorte - real.
*     acorte specifies the norm of the last computed correction for the
*     electron heat balance equation.
*
*     acorti - real.
*     acorti specifies the norm of the last computed correction for the
*     all atom heat balance equation.
*
*     acorpo - real.
*     acorpo specifies the norm of the last computed correction for the
*     electric potential equation.
*
*   ..rate coefficients
*
*     /b2cmrc/ - common block.
*     /b2cmrc/ contains a table of atomic rate coefficients.
*      rtnt - integer constant.
*     rtnt specifies the number of intervals in the discretisation
*     of the temperature range. The table will contain data for
*     (rtnt+1) values of the temperature variable.
*      rtnn - integer constant.
*     rtnn specifies the number of intervals in the discretisation
*     of the density range. The table will contain data for (rtnn+1)
*     values of the density variable.
*      rtns - integer.
*     rtns specifies the total number of atomic species; i.e.
*       (sum inuc : 0.le.inuc.lt.nnuc : izhi(inuc)-izlo(inuc)+1) ,
*     where inuc is a nuclear species index, izlo(inuc) is the lowest
*     and izhi(inuc) is the highest atomic charge considered in b2 for
*     that nuclear species.
*     rtns equals the variable ns used elsewhere in b2.
*      rtdumm - integer.
*     rtdumm is a dummy integer variable, introduced in order to obtain
*     an even-length sequence of integer variables in the common block.
*      rtzmin, rtzmax, rtzn - (0:rtns-1) real*8 array.
*     For 0.le.is.lt.rtns, rtzmin(is) and rtzmax(is) specify the atomic 
*     charge range, and rtzn(is) specifies the nuclear charge
*     of atomic species (is). 
*     The izhi(inuc)-izlo(inuc)+1 charge states for each nuclear
*     species are taken in natural order, starting at izlo(inuc); the
*     highest state of nuclear species (inuc) is followed by the lowest
*     state of (inuc+1). rtzmin, rtzmax and rtzn are equal to the arrays 
*     zamin, zamax and zn used elsewhere in b2.
*      rtt - (0:rtnt) real*8 array.
*     rtt(0:rtnt) specifies the sequence of values of the temperature
*     variable, expressed in ev, for which the rate coefficients are
*     evaluated. rtt(0:rtnt) will be positive and increasing.
*      rtn - (0:rtnn) real*8 array.
*     rtn(0:rtnn) specifies the sequence of values of the density
*     variable, expressed in m**-3, for which the rate coefficients
*     are evaluated. rtn(0:rtnn) will be positive and increasing.
*      rtlt - (0:rtnt) real*8 array.
*     rtlt(0:rtnt) holds log(rtt(0:rtnt)).
*      rtln - (0:rtnn) real*8 array.
*     rtln(0:rtnn) holds log(rtn(0:rtnn)).
*      rtlsa - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlsa(i,j,k) specifies
*     logu(rtrsa), where rtrsa is the ionisation rate coefficient at
*     electron temperature rtt(i), electron density rtn(j), for the
*     process k->k+1. rtrsa has dimension m**3/s. The function logu is
*     the same as (log) except that logu(0.0) is equal to some huge
*     negative number.
*      rtlra - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlra(i,j,k) specifies
*     logu(rtrra), where rtrra is the recombination rate coefficient at
*     electron temperature rtt(i) and electron density rtn(j), for the
*     process k->k-1. rtrra has dimension m**3/s.
*      rtlqa - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlqa(i,j,k) specifies
*     logu(rtrqa), where rtrqa is the electron heat loss rate
*     coefficient at electron temperature rtt(i) and electron density
*     rtn(j), for processes starting from atomic species k. rtrqa has
*     dimension ev*m**3/s.
*      rtlcx - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlcx(i,j,k) specifies
*     logu(rtrcx), where rtrcx is the charge exchange rate coefficient
*     at neutral hydrogen temperature/amu rtt(i) and electron density
*     rtn(j), for the process k->k-1. rtrcx has dimension m**3/s.
*      rtlza - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlza(i,j,k) specifies
*     the effective charge state at
*     electron temperature rtt(i), electron density rtn(j), for the
*     superstage k. rtlza has dimension e. 
*      rtlz2 - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlz2(i,j,k) specifies
*     the effective squared charge at
*     electron temperature rtt(i), electron density rtn(j), for the
*     superstage k. rtlz2 has dimension e**2.
*      rtlpt - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlpt(i,j,k) specifies
*     the cumulative ionisation energy at
*     electron temperature rtt(i), electron density rtn(j), 
*     to reach superstage k. rtlpt has dimension eV. 
*      rtlpi - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlpi(i,j,k) specifies
*     the effective ionisation energy at
*     electron temperature rtt(i), electron density rtn(j), 
*     to reach superstage k from superstage k-1. rtlpi has dimension eV. 
*
*     nscx - integer.
*     nscx is the number of charge-exchanging neutral hydrogen species
*     that are present in the calculation.  (So, usually nscx.eq.1).
*     It will hold that 0.le.nscx.lt.10.
*
*     iscx - (0:nscxmax-1) integer array.
*     iscx(0:nscx-1) will hold the species indices for the various kinds
*     of neutral hydrogen that are present in the calculation.
*     For 0.le.k.lt.nscx it will hold that 0.le.iscx(k).lt.ns.
*
*     rlsa - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), rlsa(ix,iy,0,is) and
*     rlsa(ix,iy,1,is) are the coefficients in the representation
*       rsa = exp(rlsa(,,0,)+rlsa(,,1,)*log(te(ix,iy)/ev))
*     that specifies the rate coefficient for ionisation from atomic
*     species (is) to species (is+1). (For is.eq.ns-1 or in case the
*     atomic species (is) and (is+1) belong to different nuclear
*     species this expression will underflow, and its value is then
*     understood as 0.)
*     rsa has dimension m**3/s.
*
*     rlra - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), rlra(ix,iy,0,is) and
*     rlra(ix,iy,1,is) are the coefficients in the representation
*       rra = exp(rlra(,,0,)+rlra(,,1,)*log(te(ix,iy)/ev))
*     that specifies the rate coefficient for recombination from atomic
*     species (is) to species (is-1). (For is.eq.0 or in case the
*     atomic species (is) and (is-1) belong to different nuclear
*     species this expression will underflow, and its value is then
*     understood as 0.)
*     rra has dimension m**3/s.
*
*     rlqa - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), rlqa(ix,iy,0,is) and
*     rlqa(ix,iy,1,is) are the coefficients in the representation
*       rqa = ev*exp(rlqa(,,0,)+rlqa(,,1,)*log(te(ix,iy)/ev))
*     that specifies the rate coefficient for electron heat loss terms
*     that are proportional to the density of the atomic species (is).
*     rqa has dimension J*m**3/s.
*
*     rlcx - (-1:nx,-1:ny,0:1,0:ns-1,0:nscxmax-1) real*8 array.
*     For (ix,iy,is,k) in (-1:nx,-1:ny,0:ns-1,0:nscx-1),
*     rlcx(ix,iy,0,is,k) and rlcx(ix,iy,1,is,k) are the coefficients in
*     the expression
*       rcx = exp(rlcx(,,0,,)+rlcx(,,1,,)*log(th0/(amh0*ev)))
*     that specifies the rate coefficient for charge exchange
*     recombination from state (is) to state (is-1) by collisions with
*     the neutral hydrogen species iscx(k); here, amh0 is the atomic
*     mass of that neutral species. (For is.eq.0 or in case the atomic
*     species (is) and (is-1) belong to different nuclear species this
*     expression will underflow, and its value is then understood as 0.)
*     rcx has dimension m**3/s.
*
*     rlza - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), rlza(ix,iy,0,is) and
*     rlza(ix,iy,1,is) are the coefficients in the representation
*       rza = rlza(,,0,)+rlza(,,1,)*log(te(ix,iy)/ev)
*     that specifies the effective charge of stage (is).
*     rza has dimension e.
*
*     rlz2 - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), rlz2(ix,iy,0,is) and
*     rlz2(ix,iy,1,is) are the coefficients in the representation
*       rz2 = rlz2(,,0,)+rlz2(,,1,)*log(te(ix,iy)/ev)
*     that specifies the effective square charge of stage (is).
*     rz2 has dimension e**2.
*
*     rlpt - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), rlpt(ix,iy,0,is) and
*     rlpt(ix,iy,1,is) are the coefficients in the representation
*       rpt = rlpt(,,0,)+rlpt(,,1,)*log(te(ix,iy)/ev)
*     that specifies the cumulative ionisation potential to stage (is).
*     rpt has dimension eV.
*
*     rlpi - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), rlpi(ix,iy,0,is) and
*     rlpi(ix,iy,1,is) are the coefficients in the representation
*       rpi = rlpi(,,0,)+rlpi(,,1,)*log(te(ix,iy)/ev)
*     that specifies the effective ionisation potential to stage (is)
*     from stage (is-1).
*     rpi has dimension eV.
*
*     rsa - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), rsa(ix,iy,is) specifies
*     the rate coefficient for ionisation from atomic species (is) to
*     species (is+1). For is.eq.ns-1 or in case the atomic species (is)
*     and (is+1) belong to different nuclear species, this coefficient
*     will be 0.
*     rsa has dimension m**3/s.
*
*     rra - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), rra(ix,iy,is) specifies
*     the rate coefficient for recombination from atomic species (is)
*     to species (is-1). For is.eq.0 or in case the atomic species (is)
*     and (is-1) belong to different nuclear species, this coefficient
*     will be 0.
*     rra has dimension m**3/s.
*
*     rqa - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), rqa(ix,iy,is) specifies
*     the rate coefficient for electron heat loss terms that are
*     proportional to the density of the atomic species (is).
*     rqa has dimension J*m**3/s.
*
*     rcx - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), rcx(ix,iy,is) specifies
*     the rate coefficient for charge exchange recombination from
*     atomic species (is) to species (is-1) by collisions with the
*     selected species of neutral hydrogen. For is.eq.0 or in case the
*     atomic species (is) and (is-1) belong to different nuclear
*     species, this coefficient will be 0.
*     rcx has dimension m**3/s.
*
*   ..boundary conditions, parameters
*
*     /b2cmpb/ - common block.
*     /b2cmpb/ specifies physics parameters for boundary conditions.
*     For this purpose the boundary is divided up into regions, and for
*     each region a set of coefficients defines the local conditions.
*     For detailed description of /b2cmpb/ see the code, especially the
*     routine b2stbc.
*      cbnrmx - integer constant.
*     cbnrmx is an upper bound for the number of boundary regions nreg,
*     (nreg.eq.cbnrso+cbnrno+cbnrwe+cbnrea) and is used in dimension
*     specifications.
*      cbirso, cbirno, cbirwe, cbirea - integer.
*     cbirso(-no,-we,-ea) is the index of the first region on the
*     south (north, west, east) boundary.
*      cbnrso, cbnrno, cbnrwe, cbnrea - integer.
*     cbnrso(-no,-we,-ea) is the number of regions on the south
*     (north, west, east) boundary.
*     It will hold that 1.le.cbnrso, 1.le.cbnrno, 1.le.cbnrwe, and
*     1.le.cbnrea; and 0.eq.cbirso, cbirso+cbnrso.eq.cbirno,
*     cbirno+cbnrno.eq.cbirwe, cbirwe+cbnrwe.eq.cbirea.
*      cbrbrk - (0:cbnrmx-1) real*8 array.
*     cbrbrk defines the boundaries between the boundary regions.
*      cbsna - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.is.lt.ns and 0.le.ireg.lt.nreg, the eight coefficients
*     cbsna(0:7,is,ireg) specify standard form boundary conditions for
*     the continuity equation for species (is) on the (ireg) boundary
*     region.
*      cbsmo - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.is.lt.ns and 0.le.ireg.lt.nreg, the eight coefficients
*     cbsmo(0:7,is,ireg) specify standard form boundary conditions for
*     the parallel momentum balance equation for species (is) on the
*     (ireg) boundary region.
*      cbshi - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.is.lt.ns and 0.le.ireg.lt.nreg, the eight coefficients
*     cbshi(0:7,is,ireg) specify standard form boundary conditions for
*     the contribution to that atom heat equation by species (is) on
*     the (ireg) boundary region.
*      cbshe - (0:7,0:cbnrmx-1) real*8 array.
*     For 0.le.ireg.lt.nreg, the eight coefficients cbshe(0:7,ireg)
*     specify standard form boundary conditions for the electron heat
*     equation on the (ireg) boundary region.
*      cbsch - (0:7,0:cbnrmx-1) real*8 array.
*     For 0.le.ireg.lt.nreg, the eight coefficients cbsch(0:7,ireg)
*     specify standard form boundary conditions for the electric
*     potential equation on the (ireg) boundary region.
*      cbrec - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.is.lt.ns-1 and 0.le.ireg.lt.nreg, the eight coefficients
*     cbrec(0:7,is,ireg) specify the recycling process for species (is)
*     on the boundary region ireg. The exact interpretation of these
*     coefficients is defined by the routine that computes the source
*     terms due to the recycling process (b2stbr). Routine b2stbr is
*     problem dependent.
*      cbmsa - (0:7,0:cbnrmx-1) real*8 array.
*      cbmsc - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.ireg.lt.nreg, the eight coefficients cbmsa(0:7,ireg)
*     specify miscellaneous processes on the boundary region ireg.
*     For 0.le.is.lt.ns-1 and 0.le.ireg.lt.nreg, the eight coefficients
*     cbmsc(0:7,is,ireg) specify miscellaneous processes for species
*     (is) on the boundary region ireg. The interpretation of these
*     coefficients is defined by the routine that computes the source
*     terms due to miscellaneous processes (b2stbm). Routine b2stbm is
*     application dependent.
*
*   ..source terms from neutral gas model
*
*     /b2cmsr/ - common block.
*     /b2cmsr/ specifies source terms from neutral gas model.
*     (Presently empty--to be used in conjunction with a Monte Carlo
*     neutral gas code.)
*
*   ..source coefficients
*
*     sna - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), sna(ix,iy,0,is) and
*     sna(ix,iy,1,is) are the coefficients in the expression
*       sna(ix,iy,0,is)+sna(ix,iy,1,is)*na(ix,iy,is)
*     that specifies a linear approximation to the assumed particle
*     source of species (is) integrated over the (ix,iy) cell.
*
*     sne - (-1:nx,-1:ny,0:1) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), sne(ix,iy,0) and
*     sne(ix,iy,1) are the coefficients in the expression
*       sne(ix,iy,0)+sne(ix,iy,1)*ne(ix,iy)
*     that specifies a linear approximation to the assumed non-ambipolar
*     electron source of integrated over the (ix,iy) cell.
*
*     smo - (-1:nx,-1:ny,0:3,0:ns-1) real*8 array.
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), smo(ix,iy,0:3,is) are the
*     coefficients in the expression
*       smo(ix,iy,0,is)+
*         smo(ix,iy,1,is)*ua(ix,iy,is)+
*         smo(ix,iy,2,is)*roxa(ix,iy,is)+
*         smo(ix,iy,3,is)*roxa(ix,iy,is)*ua(ix,iy,is)
*     that specifies a linear approximation to the parallel momentum
*     source for species (is) integrated over the (ix,iy) cell; here
*       roxa(ix,iy,is) = am(is)*mp*na(ix,iy,is) .
*
*     she - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), she(ix,iy,0:3) are the coefficients
*     in the expression
*       she(ix,iy,0)+
*         she(ix,iy,1)*te(ix,iy)+
*         she(ix,iy,2)*ne(ix,iy)+
*         she(ix,iy,3)*ne(ix,iy)*te(ix,iy)
*     that specifies a linear approximation to the assumed electron heat
*     source integrated over the (ix,iy) cell.
*
*     shi - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), shi(ix,iy,0:3) are the coefficients
*     in the expression
*       shi(ix,iy,0)+
*         shi(ix,iy,1)*ti(ix,iy)+
*         shi(ix,iy,2)*ni(ix,iy)+
*         shi(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
*     that specifies a linear approximation to the assumed all atom heat
*     source integrated over the (ix,iy) cell.
*
*     sch - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), sch(ix,iy,0:3) are the coefficients
*     in the expression
*       sch(ix,iy,0)+
*         sch(ix,iy,1)*po(ix,iy)+
*         sch(ix,iy,2)*ne(ix,iy)+
*         sch(ix,iy,3)*ne(ix,iy)*po(ix,iy)
*     that specifies a linear approximation to the assumed electric
*     charge source integrated over the (ix,iy) cell.
*
*   ..partial source coefficients
*
*     sna0 - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), sna0(ix,iy,0,is) and
*     sna0(ix,iy,1,is) are the coefficients in the expression
*       sna0(ix,iy,0,is)+sna0(ix,iy,1,is)*na(ix,iy,is)
*     that specifies a linear approximation to some specific
*     contribution to the particle source of species (is) integrated
*     over the (ix,iy) cell.
*
*     sne0 - (-1:nx,-1:ny,0:1) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), sne0(ix,iy,0) and
*     sne0(ix,iy,1) are the coefficients in the expression
*       sne0(ix,iy,0)+sne0(ix,iy,1)*ne(ix,iy)
*     that specifies a linear approximation to some specific
*     contribution to the particle non-ambipolar electron source 
*     integrated over the (ix,iy) cell.
*
*     smo0 - (-1:nx,-1:ny,0:3,0:ns-1) real*8 array.
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), smo0(ix,iy,0:3,is) are the
*     coefficients in the expression
*       smo0(ix,iy,0,is)+
*         smo0(ix,iy,1,is)*ua(ix,iy,is)+
*         smo0(ix,iy,2,is)*roxa(ix,iy,is)+
*         smo0(ix,iy,3,is)*roxa(ix,iy,is)*ua(ix,iy,is)
*     that specifies a linear approximation to some specific
*     contribution to the non-frictional parallel momentum source for
*     species (is) integrated over the (ix,iy) cell. Here,
*       roxa(ix,iy,is) = am(is)*mp*na(ix,iy,is) .
*
*     she0 - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), she0(ix,iy,0:3) are the coefficients
*     in the expression
*       she0(ix,iy,0)+
*         she0(ix,iy,1)*te(ix,iy)+
*         she0(ix,iy,2)*ne(ix,iy)+
*         she0(ix,iy,3)*ne(ix,iy)*te(ix,iy)
*     that specifies a linear approximation to some specific
*     contribution to the electron heat source integrated over the
*     (ix,iy) cell.
*
*     shi0 - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), shi0(ix,iy,0:3) are the coefficients
*     in the expression
*       shi0(ix,iy,0)+
*         shi0(ix,iy,1)*ti(ix,iy)+
*         shi0(ix,iy,2)*ni(ix,iy)+
*         shi0(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
*     that specifies a linear approximation to some specific
*     contribution to the all atom heat source integrated over the
*     (ix,iy) cell.
*
*     sch0 - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), sch0(ix,iy,0:3) are the coefficients
*     in the expression
*       sch0(ix,iy,0)+
*         sch0(ix,iy,1)*po(ix,iy)+
*         sch0(ix,iy,2)*ne(ix,iy)+
*         sch0(ix,iy,3)*ne(ix,iy)*po(ix,iy)
*     that specifies a linear approximation to some specific
*     contribution to the electric charge source integrated over the
*     (ix,iy) cell.
*
*   ..transport coefficients, parameters
*
*     /b2cmpt/ - common block.
*     (Note: the integer constant nsdecl, which comes with /b2cmpa/,
*     is used instead of (ns) in dimension specifications.)
*      cfdf0 - (0:7,0:ns-1) real*8 array.
*     cfdf0 contains some auxiliary coefficients.  For 0.le.is.lt.ns,
*     the local expression
*       df0(,,is) = cfdf0(0,is)+cfdf0(1,is)/(ni(,,0)/1.0e20)+
*         cfdf0(2,is)*sqrt(ti(,)/(am(is)*mp))/
*          (cfdf0(3,is)*ni(,,0)+cfdf0(4,is)*ne(,))
*     is used in the definition of the anomalous transport coefficients.
*     df0 is meant to represent a diffusivity.  If cfdf0(2,is).eq.0 then
*     the corresponding term in df0 is understood as 0 independent of
*     cfdf0(3:4,is).
*      cfdna - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfdna(,is) are the coefficients in the
*     local expression
*       difna(,,is) = cfdna(0,is)+cfdna(1,is)/(ne(,)/1.0e20)+
*         cfdna(2,is)*df0(,)+
*         cfdna(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous diffusivity with respect to
*     grad.na; the transport equation for species (is) contains a term
*     difna(,,is)*'grad.na'.
*     If cfdna(7,is) is not zero, the dna0 array is flux-scaled with
*     a scaling factor of cfdna(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfdna(4:6,) are not used in the default version of the code.
*      cfdpa - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfdpa(,is) are the coefficients in the
*     local expression
*       difpa(,,is) = cfdpa(0,is)+cfdpa(1,is)/(ne(,)/1.0e20)+
*         cfdpa(2,is)*df0(,,is)+
*         cfdpa(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous diffusivity with respect to
*     grad.pa; the transport equation for species (is) contains a term
*     (1/(rza(,,is)*te(,)+ti(,)))*difpa(,,is)*'grad.pa'.
*     If cfdpa(7,is) is not zero, the dpa0 array is flux-scaled with
*     a scaling factor of cfdpa(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfdpa(4:6,) are not used in the default version of the code.
*      cfvla - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfvla(,is) are the coefficients in the
*     local expression
*       vlay(,,is) = cfvla(0,is)+cfvla(1,is)/(ne(,)/1.0e20)+
*         cfvla(2,is)*df0(,,is)+
*         cfvla(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous y-velocity of species (is); the
*     transport equation for species (is) contains a term
*     vlay(,,is)*na(,,is).
*     If cfvla(7,is) is not zero, the vla0 array is flux-scaled with
*     a scaling factor of cfvla(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfvla(4:6,) are not used in the default version of the code.
*      cfvsa - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfvsa(,is) are the coefficients in the
*     local expression
*       dfvsa(,,is) = cfvsa(0,is)+cfvsa(1,is)/(ne(,)/1.0e20)+
*         cfvsa(2,is)*df0(,,is)+
*         cfvsa(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the diffusivity corresponding to the anomalous
*     transverse viscosity, so that this viscosity is
*       vsa(,,is) = mp*am(is)*na(,,is)*dfvsa
*     If cfvsa(7,is) is not zero, the vsa0 array is flux-scaled with
*     a scaling factor of cfvsa(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfvsa(4:6,) are not used in the default version of the code.
*      cfhci - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfhci(,is) are the coefficients in the
*     local expression
*       dfhca(,,is) = cfhci(0,is)+cfhci(1,is)/(ne(,)/1.0e20)+
*         cfhci(2,is)*df0(,,is)+
*         cfhci(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous atom heat conductivity associated
*     with species (is), so that the atom heat conduction coefficient
*     contains a term dfhca(,,is)*na(,,is).
*     If cfhci(7,is) is not zero, the hci arrays are flux-scaled with
*     a scaling factor of cfhci(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfhci(4:6,) are not used in the default version of the code.
*      cfhce - (0:7) real*8 array.
*     cfhce() are the coefficients in the local expression
*       dfhce(,) = cfhce(0)+cfhce(1)/(ne(,)/1.0e20)+
*         cfhce(3)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous electron heat conductivity, so that
*     the electron heat conduction coefficient is dfhce(,)*ne(,).
*     If cfhce(7) is not zero, the hce0 array is flux-scaled with
*     a scaling factor of cfhce(7)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfhce(2) and (4:6) are not used in the default version of the code.
*      cfsig - (0:7) real*8 array.
*     cfsig() are the coefficients in the local expression
*       dfsig(,) = cfsig(0)+cfsig(1)/(ne(,)/1.0e20)+
*         cfsig(3)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the diffusivity corresponding to the anomalous
*     electrical conductivity, which is
*       sig(,) = qe*ne(,)*dfsig(,).
*     If cfsig(7) is not zero, the sig0 array is flux-scaled with
*     a scaling factor of cfsig(7)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfsig(2) and (4:6) are not used in the default version of the code.
*      cfalf - (0:7) real*8 array.
*     cfalf() are the coefficients in the local expression
*       dfalf(,) = cfalf(0)+cfalf(1)/(ne(,)/1.0e20)+
*         cfalf(3)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the diffusivity corresponding to the anomalous
*     thermo-electric coefficient, which is
*       alf(,) = ne(,)*sqrt(qe/te(,))*dfalf(,).
*     If cfalf(7) is not zero, the alf0 array is flux-scaled with
*     a scaling factor of cfalf(7)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfalf(2) and (4:6) are not used in the default version of the code.
*      cflim - (0:7) real*8 array.
*     cflim contains coefficients for flux limiting. cflim(0)=cflme
*     specifies the flux limit factor for poloidal electron heat
*     transport. cflim(1)=cflmi specifies the flux limit factor for
*     poloidal all atom heat transport. cflim(2)=cflmv specifies the
*     flux limit factor for poloidal atom momentum transport.
*     Zero settings turn off these flux limits.
*     cflim(3:7) are not used in the default version of the code.
*
*   ..transport coefficients, differential
*
*     dna - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let dnax(,,)=dna(,,0,) and dnay(,,)=dna(,,1,).)
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), dnax(ix,iy,is) and
*     dnay(ix,iy,is) specify the x-component and the y-component of
*     the diffusivity of species (is) on the (ix,iy) cell. The
*     corresponding term in the particle flux density is
*     -dna*grad.na .
*
*     dpa - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let dpax(,,)=dpa(,,0,) and dpay(,,)=dpa(,,1,).)
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), dpax(ix,iy,is) and
*     dpay(ix,iy,is) specify the x-component and the y-component of
*     the diffusion coefficient (w.r.t. grad.pa) of species (is) on
*     the (ix,iy) cell. The corresponding term in the particle flux
*     density is  -dpa*grad.pa .
*
*     vla - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let vlax(,,)=vla(,,0,) and vlay(,,)=vla(,,1,).)
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vlax(ix,iy,is) specifies
*     the x-component of the anomalous velocity of species (is) on the
*     (ix,iy) cell.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vlay(ix,iy,is) specifies
*     the y-component of the anomalous velocity of species (is) on the
*     (ix,iy) cell.
*
*     vsa - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let vsax(,,)=vsa(,,0,) and vsay(,,)=vsa(,,1,).)
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vsax(ix,iy,is) specifies
*     the x-component and vsay(ix,iy,is) the y-component of the
*     viscosity of species (is) on the (ix,iy) cell.
*
*     hce - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let hcex(,)=hce(,,0) and hcey(,)=hce(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), hcex(ix,iy) specifies the
*     x-component and hcey(ix,iy) the y-component of the electron
*     thermal conductivity on the (ix,iy) cell.
*
*     hci - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let hcix(,)=hci(,,0) and hciy(,)=hci(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), hcix(ix,iy) specifies the
*     x-component and hciy(ix,iy) the y-component of the all atom
*     thermal conductivity on the (ix,iy) cell.
*
*     sig - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let sigx(,)=sig(,,0) and sigy(,)=sig(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), sigx(ix,iy) specifies the
*     x-component and sigy(ix,iy) the y-component of the electrical
*     conductivity on the (ix,iy) cell.
*
*     alf - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let alfx(,)=alf(,,0) and alfy(,)=alf(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), alfx(ix,iy) specifies the
*     x-component and alfy(ix,iy) the y-component of the thermo-
*     electric coefficient on the (ix,iy) cell.
*
*     dna0 - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), dna0(ix,iy,is) specifies
*     the anomalous diffusivity with respect to grad(na) of species
*     (is) on the (ix,iy) cell.
*
*     dpa0 - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), dpa0(ix,iy,is) specifies
*     the anomalous diffusivity with respect to grad(pa) of species
*     (is) on the (ix,iy) cell.
*
*     vla0 - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let vla0x(,,)=vla0(,,0,) and vla0y(,,)=vla0(,,1,).)
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vla0x(ix,iy,is) specifies
*     the x-component and vla0y(ix,iy,is) the y-component of the
*     anomalous velocity of species (is) on the (ix,iy) cell.
*
*     vsa0 - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vsa0(ix,iy,is) specifies
*     the anomalous viscosity of species (is) on the (ix,iy) cell.
*
*     hce0 - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), hce0(ix,iy) specifies the anomalous
*     electron thermal conductivity on the (ix,iy) cell.
*
*     hci0 - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), hci0(ix,iy) specifies the anomalous
*     atom thermal conductivity on the (ix,iy) cell.
*
*     sig0 - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), sig0(ix,iy) specifies the anomalous
*     electrical conductivity on the (ix,iy) cell.
*
*     alf0 - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), alf0(ix,iy) specifies the anomalous
*     thermo-electric coefficient on the (ix,iy) cell.
*
*   ..transport coefficients, discrete
*
*     ceqp - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), ceqp(ix,iy) specifies the energy
*     equipartition multiplier integrated over the (ix,iy) cell.
*     (The cell-integrated energy equipartition coefficient is
*     ceqp*ne*ni/te**1.5.)
*
*     cdna - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let cdnax=cdna(,,0,) and cdnay=cdna(,,1,).)
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), cdnax(ix,iy,is) specifies
*     the particle diffusivity w.r.t. grad(na) of species (is) between
*     the (ix-1,iy) and the (ix,iy) cell. This coefficient enters into
*     the expression for the particle flux in the form
*       -cdnax(ix,iy,is)*(na(ix,iy,is)-na(ix-1,iy,is))
*     For ix.eq.-1, cdnax(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), cdnay(ix,iy,is) specifies
*     the particle diffusivity w.r.t. grad(na) of species (is) between
*     the (ix,iy-1) and the (ix,iy) cell. This coefficient enters into
*     the expression for the particle flux in the form
*       -cdnay(ix,iy,is)*(na(ix,iy,is)-na(ix,iy-1,is))
*     For iy.eq.-1, cdnay(-1:nx,iy,0:ns-1) will hold 0.
*
*     cdpa - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let cdpax=cdpa(,,0,) and cdpay=cdpa(,,1,).)
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), cdpax(ix,iy,is) specifies
*     the particle diffusivity w.r.t. grad(pa) of species (is) between
*     the (ix-1,iy) and the (ix,iy) cell. This coefficient enters into
*     the expression for the particle flux in the form
*       -cdpax(ix,iy,is)*(pa(ix,iy,is)-pa(ix-1,iy,is))
*     For ix.eq.-1, cdpax(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), cdpay(ix,iy,is) specifies
*     the particle diffusivity w.r.t. grad(pa) of species (is) between
*     the (ix,iy-1) and the (ix,iy) cell. This coefficient enters into
*     the expression for the particle flux in the form
*       -cdpay(ix,iy,is)*(pa(ix,iy,is)-pa(ix,iy-1,is))
*     For iy.eq.-1, cdpay(-1:nx,iy,0:ns-1) will hold 0.
*
*     cvla - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let cvlax=cvla(,,0,) and cvlay=cvla(,,1,).)
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), cvlax(ix,iy,is) specifies
*     the surface-integrated anomalous velocity of species (is) between
*     the (ix-1,iy) and the (ix,iy) cell. This coefficient enters into
*     the expression for the particle flux in the form
*       cvlax(ix,iy,is)*(na)
*     in which (na) involves some interpolation to cell faces.
*     For ix.eq.-1, cvlax(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), cvlay(ix,iy,is) specifies
*     the surface-integrated anomalous velocity of species (is) between
*     the (ix,iy-1) and the (ix,iy) cell. This coefficient enters into
*     the expression for the particle flux in the form
*       cvlay(ix,iy,is)*(na)
*     in which (na) involves some interpolation to cell faces.
*     For iy.eq.-1, cvlay(-1:nx,iy,0:ns-1) will hold 0.
*
*     cvsa - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let cvsax=cvsa(,,0,) and cvsay=cvsa(,,1,).)
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), cvsax(ix,iy,is) specifies
*     the viscous conductance of species (is) between the (ix-1,iy) and
*     the (ix,iy) cell.
*     For ix.eq.-1, cvsax(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), cvsay(ix,iy,is) specifies
*     the viscous conductance of species (is) between the (ix,iy-1) and
*     the (ix,iy) cell.
*     For iy.eq.-1, cvsay(-1:nx,iy,0:ns-1) will hold 0.
*
*     chce - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let chcex=chce(,,0) and chcey=chce(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), chcex(ix,iy) specifies the electron
*     thermal conductance between the (ix-1,iy) and the (ix,iy) cell.
*     This coefficient enters into the expression for the electron
*     heat flux in the form
*       -chcex(ix,iy)*(te(ix,iy)-te(ix-1,iy)) .
*     For ix.eq.-1, chcex(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), chcey(ix,iy) specifies the electron
*     thermal conductance between the (ix,iy-1) and the (ix,iy) cell.
*     This coefficient enters into the expression for the electron
*     heat flux in the form
*       -chcey(ix,iy)*(te(ix,iy)-te(ix,iy-1)) .
*     For iy.eq.-1, chcey(-1:nx,iy) will hold 0.
*
*     chve - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let chvex=chve(,,0) and chvey=chve(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), chvex(ix,iy) specifies the electron
*     thermal strange velocity between the (ix-1,iy) and the (ix,iy)
*     cell. This coefficient enters into the expression for the
*     electron heat flux in the form  chvex*ne*te , interpolated to
*     cell faces in a suitable manner.
*     For ix.eq.-1, chvex(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), chvey(ix,iy) specifies the electron
*     thermal strange velocity between the (ix,iy-1) and the (ix,iy)
*     cell. This coefficient enters into the expression for the
*     electron heat flux in the form  chvey*ne*te , interpolated to
*     cell faces in a suitable manner.
*     For iy.eq.-1, chvey(-1:nx,iy) will hold 0.
*
*     chci - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let chcix=chci(,,0) and chciy=chci(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), chcix(ix,iy) specifies the all atom
*     thermal conductance between the (ix-1,iy) and the (ix,iy) cell.
*     This coefficient enters into the expression for the atom heat
*     flux in the form
*       -chcix(ix,iy)*(ti(ix,iy)-ti(ix-1,iy)) .
*     For ix.eq.-1, chcix(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), chciy(ix,iy) specifies the all atom
*     thermal conductance between the (ix,iy-1) and the (ix,iy) cell.
*     This coefficient enters into the expression for the atom heat
*     flux in the form
*       -chciy(ix,iy)*(ti(ix,iy)-ti(ix,iy-1)) .
*     For iy.eq.-1, chciy(-1:nx,iy) will hold 0.
*
*     chvi - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let chvix=chvi(,,0) and chviy=chvi(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), chvix(ix,iy) specifies the all atom
*     thermal strange velocity between the (ix-1,iy) and the (ix,iy)
*     cell. This coefficient enters into the expression for the
*     atom heat flux in the form  chvix*ni*ti , interpolated to cell
*     faces in a suitable manner.
*     For ix.eq.-1, chvix(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), chviy(ix,iy) specifies the all atom
*     thermal strange velocity between the (ix,iy-1) and the (ix,iy)
*     cell. This coefficient enters into the expression for the
*     atom heat flux in the form  chviy*ni*ti , interpolated to
*     cell faces in a suitable manner.
*     For iy.eq.-1, chviy(-1:nx,iy) will hold 0.
*
*     csig - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let csigx=csig(,,0) and csigy=csig(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), csigx(ix,iy) specifies the electrical
*     conductance between the (ix-1,iy) and the (ix,iy) cell.
*     For ix.eq.-1, csigx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), csigy(ix,iy) specifies the electrical
*     conductance between the (ix,iy-1) and the (ix,iy) cell.
*     For iy.eq.-1, csigx(-1:nx,iy) will hold 0.
*
*     calf - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let calfx=calf(,,0) and calfy=calf(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), calfx(ix,iy) specifies the thermo-
*     electric coefficient between the (ix-1,iy) and the (ix,iy) cell.
*     For ix.eq.-1, calfx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), calfy(ix,iy) specifies the thermo-
*     electric coefficient between the (ix,iy-1) and the (ix,iy) cell.
*     For iy.eq.-1, calfy(-1:nx,iy) will hold 0.
*
*     csigin - (-1:nx,-1:ny,0:1,0:ns-1,0:ns-1) real*8 array.		 !xpb
*     (Let csiginx=csigin(,,0,,) and csiginy=csigin(,,1,,).)
*     For (ix,iy,is,is0) in (0:nx,-1:ny,0:ns-1,0:ns-1), 
*     csiginx(ix,iy,is,is0) specifies the ion-neutral
*     conductivity between the (ix-1,iy) and the (ix,iy) cell, where 
*     the ion species is species is and the neutral species is is0.
*     For ix.eq.-1, csiginx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), csiginy(ix,iy) specifies the ion-neutral
*     conductivity between the (ix,iy-1) and the (ix,iy) cell, where 
*     the ion species is species is and the neutral species is is0.
*     For iy.eq.-1, csiginx(-1:nx,iy) will hold 0.
*     When is is a neutral species or is0 is an ionised species,
*     csigin(,,,is,is0) will be zero.
*
*     csig_an - (-1:nx,-1:ny,0:1) real*8 array.				 !xpb
*     (Let csig_anx=csig_an(,,0) and csig_any=csig_an(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), csig_anx(ix,iy) specifies the anomalous
*     conductivity between the (ix-1,iy) and the (ix,iy) cell.
*     For ix.eq.-1, csig_anx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), csig_any(ix,iy) specifies the anomalous
*     conductivity between the (ix,iy-1) and the (ix,iy) cell.
*     For iy.eq.-1, csig_anx(-1:nx,iy) will hold 0.
*
*     calf_an - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let calf_anx=calf_an(,,0) and calf_any=calf_an(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), calf_anx(ix,iy) specifies the 
*     anomalous part of the thermo-electric coefficient 
*     between the (ix-1,iy) and the (ix,iy) cell.
*     For ix.eq.-1, calfx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), calf_any(ix,iy) specifies the 
*     anomalous part of the thermo-electric coefficient 
*     between the (ix,iy-1) and the (ix,iy) cell.
*     For iy.eq.-1, calf_any(-1:nx,iy) will hold 0.
*
*   ..parameters for core numerical iteration
*
*     ntim - integer.
*     ntim specifies the number of timesteps to be performed.
*
*     itim - integer.
*     itim specifies the timestep counter. It is initialised to 0 and
*     incremented by 1 following each successful call to b2mndt.
*
*     dtim - real.
*     dtim specifies the timestep for the present call to b2mndt.
*     It will hold that 0.lt.dtim.
*
*     itcnt - integer.
*     itcnt counts the number of inner iterations. It is initialised to
*     0 and incremented by 1 following each call to b2news.
*
*     nstg - (0:2) integer array.
*     nstg(0:2) specifies the maximum number of iterations for each of 
*     the three main loops in the subroutine b2mndt.
*
*     rxf - real.
*     rxf specifies the main under-relaxation parameter. Computed
*     corrections to the unknowns are scaled by a factor rxf before
*     being applied. It will hold that 0.le.rxf.le.1.
*
*     rxg - real.
*     rxg specifies a special under-relaxation parameter. A term of the
*     form "abs(residual)/rxg" is added to the diagonal of the matrix
*     of the correction equation, with the effect of limiting the
*     computed correction. rxg is dimensionless and of order unity;
*     smaller values of rxg imply stronger damping.
*
*   ..variables for input-output
*
*     nrid - integer.
*     nrid specifies a unit number for formatted input.
*
*     nget - integer.
*     nget specifies a unit number for un*formatted input.
*
*     nprt - integer.
*     nprt specifies a unit number for formatted output.
*
*     nput - integer.
*     nput specifies a unit number for un*formatted output.
*
*     idout - character*(*).
*     idout specifies the desired selection of output segments for
*     printed and graphical output. Associated with each segment is a
*     four-character word, which is identified on the printed output.
*     (The same word may be associated with several segments.)
*     A particular segment will be produced only if its associated word
*     has a match, as defined by the routine strmas, in idout.
*
*     /b2cmwg/ - common block.
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
*   ..workspace
*
*     /stackw/ - common block.
*     /stackw/ contains two stacks, one for real*8 and one for integer
*     workspace.  They are used to provide local dynamic arrays within
*     the rules of Fortran 77.
*      nstckw - integer constant.
*     nstckw specifies the size of the real*8 workspace array w. If
*     nstckw is too small, the program will abort via a call to
*     xertst/xerrab.
*      nstcwi - integer constant.
*     nstcwi specifies the size of the integer workspace array wi. If
*     nstcwi is too small, the program will abort via a call to
*     xertst/xerrab.
*      w - (0:nstckw-1) real*8 array.
*     w is used as a stack for local dynamic real*8 arrays. Base
*     addresses for arrays in this stack are to be maintained by
*     appropriate calls to the subroutines subini, subend, and lwget.
*      wi - (0:nstcwi-1) integer array.
*     wi is used as a stack for local dynamic integer arrays. Base
*     addresses for arrays in this stack are to be maintained by
*     appropriate calls to the subroutines subini, subend, and lwiget.
*
*   ..some other variables
*
*     lblgm - character*120.
*     lblgm identifies the instance of the program employed to prepare
*     the geometry.
*
*     lblcp - character*120.
*     lblcp identifies the instance of the program employed to prepare
*     the table of physics parameters.
*
*     lblps - character*120.
*     lblps identifies the instance of the program employed to prepare
*     the initial state.
*
*     lblrc - character*120.
*     lblrc identifies the instance of the program employed to prepare
*     the table of atomic rate coefficients.
*
*     lblsr - character*120.
*     lblsr identifies the instance of the program employed to prepare
*     the neutral gas source data.
*
*     lblmn - character*120.
*     lblmn identifies the instance of the program employed for the main
*     calculation.
*
*     ncall - integer.
*     ncall counts the number of calls made to the present routine.
*
*
*  2. Notes
*
*     Variables local to a subroutine are not expected to be preserved
*     between calls unless they are declared in a 'save' statement.
*     Variables in common are not expected to be preserved between
*     calls to the outermost routine in which the common block appears.
*     (See ANSI X3.9-1978 Fortran 77, Section 17.3.6d.) The code may
*     be compiled and loaded under an option that initialises all
*     quantities not mentioned in a data statement to an invalid
*     numerical value.
*
*
*-----------------------------------------------------------------------
*.end b2cdcv

      end subroutine b2cdcv
