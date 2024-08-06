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
*     the B2.5 code. No computation is performed in this routine.
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
*     nCv, nCi, nCg - integer.
*     nCv = nCi + nCg specifies the number of finite volume cells in the grid.
*     nCi is the number of internal cells. nCg is the number of guard cells.
*
*   ....cell-centered variables
*
*     cvX, cvY - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), cvX(iCv) specifies the
*     Cartesian X-coordinate of the center of the cell and cvX(iCv) specifies
*     the Cartesian Y-coordinate of the center of the cell.
*
*     cvBb - (1:nCv,0:3) real*8 array.
*     For (iCv) in (1:nCv), cvBb(iCv,0:3) specifies the
*     magnetic field at the center of the (iCv) cell.
*     cvBb(iCv,0:2) are the (x,y,z)-components and cvBb(iCv,3) is
*     the absolute magnetic field strength.
*     (cvBb(,3)=sqrt(cvBb(,0)**2+cvBb(,1)**2+cvBb(,,2)**2.)
*     It will hold that 0.le.cvBb(,0:3), because the direction of the
*     magnetic field follows from cvEb (see below).
*
*     cvVol - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), cvVol(iCv) specifies the volume of
*     the (iCv) cell. (In the case of toroidal symmetry this volume
*     includes a factor 2*pi*R, while in the case of cylindrical
*     symmetry it includes a factor 2*pi*Z.)
*     It will hold that 0.lt.cvVol().
*
*     cvSz - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), cvSz(iCv) specifies the area of
*     the (iCv) cell perpendicular to the z-direction.
*
*     cvHz - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), cvHz(iCv) specifies the length of
*     the (iCv) cell in the z-direction, measured at the cell center.
*
*     cvOnedBsq - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), cvOnedBsq(iCv) specifies the quantity
*     1/B**2 in the iCv cell center.
*
*     cvBzb - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), cvBzb(iCv) specifies the quantity
*     Bz/B**2 - Bz/<B**2> in the iCv cell center, used for P-S/mdf flows.
*     Here, <B**2> is the magnetic field averaged over the core region.
*
*     cvEb - (1:nCv,0:2) real*8 array.
*     cvEb(iCv,0:2) specifies the direction of the magnetic field at
*     the center of the cell iCv in the (R,Z,toroidal) axis system.
*     (sqrt(cvEb(iCv,0)**2+cvEb(iCv,1)**2+cvEb(iCv,2)**2)=1)
*
*     signmf - real*8.
*     signfm specifies the direction of the poloidal field.
*     signmf = +1.0: the poloidal field goes clockwise around the plasma
*     signmf = -1.0: the poloidal field goes counterclockwise around the plasma
*
*   ....face-based variables
*
*     nFc - integer.
*     nFc is the number of cell faces in the grid. Each 2D cell is a polygon
*     bounded by N faces with N > 2. A face connects two internal
*     cells or an internal cell and a guard cell. There are no faces
*     between guard cells.
*
*     fcS - (1:nFc,1:2) real*8 array.
*     For iFc in (1:nFc), fcHc(iFc,:) specifies the total area of face iFc.
*
*     fcQgam - (1:nFc,0:1) real*8 array.
*     For iFc in (1:nFc), fcQgam(iFc,0) specifies cos(t) and
*     fcQgam(iFc,1) specifies sin(t), where t is the included angle
*     between the cell-connector vector and the normal vector
*     of face iFc. It will hold that 0.lt.fcQgam(,).le.1.
*
*     fcQalf - (1:nFc,0:1) real*8 array.
*     For iFc in (1:nFc), fcQalf(iFc,0) specifies cos(t) and
*     fcQalf(iFc,1) specifies sin(t), where t is the included angle
*     between the magnetic field vector and the face normal.
*
*     fcQbet - (1:nFc,0:1) real*8 array.
*     For iFc in (1:nFc), fcQbet(iFc,0) specifies cos(t) and
*     fcQbet(iFc,1) specifies sin(t), where t is the included angle
*     between the cell-connector and the magnetic field vector.
*
*     fcPbs - (1:nFc) real*8 array.
*     For iFc in (1:nFc), fcPbs(iFc) specifies the product
*     (fcBx/fcBb)*sx or (magnetic field pitch)*(surface area),
*     namely the parallel contact area.
*
*     fcPbshz - (1:nFc) real*8 array.
*     For iFc in (1:nFc), fcPbshz(iFc) specifies the product
*     (fcBx/fcBb)*sx*hz or (magnetic field pitch)*(surface area)*(hz).
*
*     fcHc - (1:nFc,1:2) real*8 array.
*     For iFc in (1:nFc), fcHc(iFc,:) are the distances from
*     the cell centers fcCv(iFc,:) to the face center of face iFc.
*
*     fcHt - (1:nFc) real*8 array.
*     For iFc in (1:nFc), fcHt is the face length of face iFc.
*
*     fcHz - (1:nFc) real*8 array.
*     For iFc in (1:nFc), fcHz is the length of the face in the z-direction.
*
*     fcBb - (1:nFc,0:3) real*8 array.
*     For iFc in (1:nFc), fcBb(iFc,0:3) specifies the magnetic field
*     at the cell faces, analogously to cvBb.
*
*     fcEb - (1:nFc,0:2) real*8 array.
*     For iFc in (1:nFc), fcEb(iFc,0:2) specifies the unit vector point along
*     the magnetic field at face center iFc. Analogous to cvEb.
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
*   ....vertex-based variables
*
*     vxX, vxY - (1:nVx) real*8 array.
*     vxX(iVx) specify the horizontal coordinate and vxY(iVx)
*     specify the vertical coordinates (in the physical plane) of
*     the vertex with number iVx.
*
*     vxFpsi - (1:nVx) real*8 array.
*     vxFpsi(iVx) specifies the potential for the magnetic field
*     components (bx,by) of the vertex with number iVx.
*     In the case of straight symmetry, bx=dpsi/dy and by=-dpsi/dx.
*     In the case of toroidal or cylindrical symmetry,
*     bx=(dpsi/dy)/(2*pi*R) and by=-(dpsi/dx)/(2*pi*R),
*     where R is the distance to the axis of symmetry.
*
*     vxFfbz - (1:nVx) real*8 array.
*     vxFfbz(iVx) specifies the field function for the magnetic
*     field component bz of the vertex with number iVx.
*     In the case of straight symmetry, bz=ffbz. In the case of
*     toroidal or cylindrical symmetry, bz=ffbz/(2*pi*R), where R is
*     the distance to the axis of symmetry.
*
*     vxEb - (1:nVx,0:2) real*8 array.
*     For iVx in (1:nVx), vxEb(iFc,0:2) specifies the unit vector point along
*     the magnetic field at the vertex iVx. Analogous to cvEb.
*
*     vxBb - (1:nVx,0:3) real*8 array.
*     For iVx in (1:nVx), vxBb(iFc,0:3) specifies the magnetic field at
*     the vertex iVx. Analogous to cvBb.
*
*     vxHz - (1:nVx) real*8 array.
*     For iVx in (1:nVx), vxHz(iFc) specifies the toroidal length of the vertex
*     iVx.
*
*     vxOnedBsq - (1:nVx) real*8 array.
*     For iVx in (1:nVx), vxOnedBsq(iFc) specifies the quantity 1/B**2 at the
*     vertex iVx.
*
*     vxBzb - (1:nVx) real*8 array.
*     For iVx in (1:nVx), vxBzb(iFc) specifies the quantity Bz/B**2 - Bz/<B**2>
*     at vertex iVx.
*
*   ....Flux-tube quantities
*
*     cvConn - (1:nCv) real*8 array.
*     For iCv in (1:nCv), cvConn(iCv) specifies the connection length for cell
*     iCv.
*
*     vxConn - (1:nVx) real*8 array.
*     For iVx in (1:nVx), vxConn(iVx) specifies the connection length for vertex
*     iVx.
*
*     ftConn - (1:nFt) real*8 array.
*     For iFt in (1:nFt), ftConn(iFt) specifies the connection length for
*     flux tube iFt.
*
*     fsConn - (1:nFs) real*8 array.
*     For iFs in (1:nFs), fsConn(iFs) specifies the connection length for
*     flux surface iFs.
*
*     ftEps - (1:nFt) real*8 array.
*     For iFt in (1:nFt), ftEps(iFt) specifies the inverse aspect ratio for
*     flux tube iFt.
*
*     ftBbav2 - (1:nFt) real*8 array.
*     For iFt in (1:nFt), ftBbav2(iFt) specifies the average |B|**2 on the flux
*     tube.
*
*     fsPsi - (1:nFs) real*8 array.
*     For iFs in (1:nFs), fsPsi(iFs) specifies Psi on each flux surface.
*
*     psi_increasing - logical
*     .true. if the Psi value increases from the O-point to the active X-point.
*
*     periodic_bc - integer
*     periodic_bc tells whether the cuts correspond to periodic boundary
*     conditions (for example in the case of limited plasmas or
*     a diverted stellarator island in a toroidally averaged geometry).
*     periodic_bc.eq.0 is the default.
*     periodic_bc.eq.1 indicates that periodic boundary conditions are
*     being used.
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
*     the bounds of the i-th singular region. Mesh locations (iCv)
*     in (nlxlo(i):nlxhi(i),nlylo(i):nlyhi(i)) make up the region, which
*     will be included in (1:nCv).
*     For nlreg.le.i.lt.nlrmx, nlxlo(i), nlxhi(i), nlylo(i) and nlyhi(i)
*     are unspecified.
*      nlloc - (0:nlrmx-1) integer array.
*     For 0.le.i.lt.4, nlloc(i) specifies the class of quantities for
*     which the i-th region is special, as follows:
*       nlloc(i) = 0: cell centers, indices (1:nCv).
*       nlloc(i) = 1: cell x-faces, indices (0:nx,-1:ny).
*       nlloc(i) = 2: cell y-faces, indices (-1:nx,0:ny).
*       nlloc(i) = 3: cell corners, indices (0:nx,0:ny).
*     For nlreg.le.i.lt.nlrmx, nlloc(i) is unspecified.
*
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
*     na - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), na(iCv,is) specifies
*     the density of atomic species (is) on the (iCv) cell.
*     It will hold that 0.lt.na(:).
*
*     ua - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), ua(iCv,is) specifies
*     the parallel velocity of atomic species (is) on the (iCv) cell.
*     The sign of the parallel velocity indicates whether it is flowing
*     in the same direction as the local poloidal magnetic field.
*
*     te - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), te(iCv) specifies the electron
*     temperature on the (iCv) cell. It will hold that 0.lt.te(:).
*
*     ti - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), ti(iCv) specifies the all atom
*     temperature on the (iCv) cell. It will hold that 0.lt.ti(:).
*
*     tn - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), tn(iCv) specifies temperature of the
*     hydrogenic neutrals on the (iCv) cell. It will hold that 0.lt.tn(:).
*     When no separate fluid neutral energy equation is solved, tn = ti.
*
*     po - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), po(iCv) specifies the electric
*     potential on the (iCv) cell.
*
*     kt - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), kt(iCv) specifies the turbulent kinetic
*     energy on the (iCv) cell.
*
*     zt - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), kt(iCv) specifies the enstrophy
*     on the (iCv) cell.
*
*     For all the fluxes below, the sign convention is as follows.
*     Fluxes are positive when they flow in the local x-direction,
*     where x is the connector line from C1 to C2.
*
*     fna - (1:nFc,0:1,0:ns-1) real*8 array.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), fna(iFc,0,is) specifies
*     the poloidal component and fna(iFc,1,is) specifies
*     the radial component of the flux of atoms of species (is) through the
*     face iFc.
*
*     fna_mdf - (1:nFc,0:1,0:ns-1) real*8 array.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), fna_mdf(iFc,0,is) specifies
*     the poloidal component and fna_mdf(iFc,1,is) specifies
*     the radial component of the modified flux of atoms of species (is)
*     through the face iFc.
*     fna_mdf is equal to fna in runs without drifts.
*
*     fna_fcor - (1:nFc,0:1,0:ns-1) real*8 array.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), fna_fcor(iFc,0,is) specifies
*     the poloidal component and fna_fcor(iFc,1,is) specifies
*     the radial component of the flux of atoms of species (is) for calculating
*     momentum transport through the face between the iFc.
*     fna_fcor takes into account the hz toroidal length derivative and
*     the Coriolis force.
*
*     fhe - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fhe(iFc,0) specifies the poloidal component and
*     fhe(iFc,1) specifies the radial component of the electron heat flux
*     through the face iFc.
*
*     fhe_mdf - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fhe_mdf(iFc,0) specifies the poloidal component and
*     fhe_mdf(iFc,1) specifies the radial component of the modified electron heat
*     flux through the face iFc.
*     fhe_mdf is equal to fhe in runs without drifts.
*
*     fhi - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fhi(iFc,0) specifies the poloidal component and
*     fhi(iFc,1) specifies the radial component of the all atom heat flux
*     through the face iFc.
*
*     fhi_mdf - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fhi_mdf(iFc,0) specifies the poloidal component and
*     fhi_mdf(iFc,1) specifies the radial component of the modified all atom
*     heat flux through the face iFc.
*     fhi_mdf is equal to fhi in runs without drifts.
*
*     fch - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fch(iFc,0) specifies the poloidal component and
*     fch(iFc,1) specifies the radial component of the electric
*     current through the face iFc.
*
*   ..auxiliary variables describing the plasma state
*
*     nb - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), nb(iCv) specifies the density of
*     one selected atomic species on the (iCv) cell. If that species
*     has index (is) then nb(:).eq.na(:,is).
*
*     ub - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), ub(iCv) specifies the parallel
*     velocity of one selected atomic species on the (iCv) cell.
*     The sign of the parallel velocity indicates whether it is flowing
*     in the same direction as the local poloidal magnetic field.
*
*     ni - (1:nCv,0:1) real*8 array.
*     For (iCv) in (1:nCv), ni(iCv,0) specifies the total atom
*     density on the (iCv) cell. It will hold that 0.lt.ni(:,0).
*     The defining expression for ni is
*       ni(:,0) = (sum is : 0.le.is.lt.ns : na(:,is)).
*     For (iCv) in (1:nCv), ni(iCv,1) specifies the total ion
*     density on the (iCv) cell. It will hold that 0.lt.ni(:,1).
*     The defining expression for ni is
*       ni(:,1) = (sum is : 0.le.is.lt.ns : na(:,is)) when rza(:,is)>0.
*
*     ne - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), ne(iCv) specifies the electron
*     density on the (iCv) cell. It will hold that 0.lt.ne(:).
*     The defining expression for ne is
*       ne(:) = (sum is : 0.le.is.lt.ns : rza(:,is)*na(:,is)).
*
*     ne2 - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), ne2(iCv) specifies the atom
*     density weighted by rz2 on the (iCv) cell. It will hold that
*     0.lt.ne2(:).
*     The defining expression for ne2 is
*       ne2(,) = (sum is : 0.le.is.lt.ns : rz2(:,is)*na(:,is)).
*
*     pr - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), pr(iCv) specifies the total
*     pressure on the (iCv) cell. The defining expression for pr is
*       pr(:) = ne(:)*te(:) + (sum is :: na(:,is)*ta(:)), where
*       ta = ti for ions and non-hydrogenic neutrals, and ta = tn for hydrogenic
*       neutrals.
*
*     pz - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), pz(iCv) specifies the plasma
*     pressure on the (iCv) cell. The defining expression for pr is
*       pr(:) = ne(:)*te(:) + (sum is :: na(:,is)*ti(:), not including neutrals).
*
*     pe - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), pe(iCv) specifies the electron
*     pressure on the (iCv) cell.
*     The defining expression for pe is pe(:) = ne(:)*te(:).
*
*     ro - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), ro(iCv) specifies the mass
*     density on the (iCv) cell. The defining expression for ro is
*       ro(:) = (sum is : 0.le.is.lt.ns : am(is)*mp**na(:,is)).
*
*     rz - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), rz(iCv) specifies the ion mass
*     density on the (iCv) cell. The defining expression for rz is
*       rz(:) = (sum is : 0.le.is.lt.ns .and. .not.is_neutral(is) :
*                         am(is)*mp**na(:,is)).
*
*     ne2m - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), ne2m(iCv) specifies the atom
*     density weighted by rz2/am on the (iCv) cell. It will hold
*     that 0.lt.ne2m(,).
*     The defining expression for ne2m is
*       ne2m(:) = (sum is :: (rz2(:,is)/am(is))*na(:,is)).
*
*     fne - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fne(iFc,0) specifies the poloidal component and
*     fne(iFc,1) specifies the radial component of the flux of
*     electrons through the face iFc.
*     The defining expression for fne is
*       fne(:,0) = (sum is :: rza(:,is)*fna(:,0,is))-(1/qe)*fch(:,0)
*       fne(:,1) = (sum is :: rza(:,is)*fna(:,1,is))-(1/qe)*fch(:,1)
*
*     fni - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fni(iFc,0) specifies the poloidal component
*     and fni(iFc,1) specifies the radial component of the flux of
*     atoms through the face iFc. If tn_style.eq.0 (default), fhi contains
*     all ions and fluid neutrals. If tn_style.eq.1, fhi contains only ions.
*     If tn_style.eq.2, fhi contains the ions and non-hydrogenic neutrals.
*
*     fnn - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fnn(iFc,0) specifies the poloidal component
*     and fnn(iFc,1) specifies the radial component of the flux of
*     hydrogenic neutral atoms through the face iFc if a separate neutral energy
*     equation is solved. If tn_style.eq.0 (default), or tn_style.eq.1, fhn = 0.
*
*     fchvispar - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fchvispar(iFc,0) specifies the poloidal component
*     and fchvispar(iFc,1) specifies the radial component of the
*     electric current driven by parallel viscosity through the face iFc.
*
*     fchvisq - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fchvisq(iFc,0) specifies the poloidal component
*     and fchvisq(iFc,1) specifies the radial component of the
*     electric current produced by the components
*     of the viscosity tensor which is connected with the heat fluxes through
*     face iFc.
*
*     fchinert - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fchinert(iFc,0) specifies the poloidal component
*     and fchinert(iFc,1) specifies the radial component of the
*     contributions from inertia and gyroviscosity current through face iFc.
*
*     fchanml - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fchanml(iFc,0) specifies the poloidal component
*     and fchanml(iFc,1) specifies the radial component of the
*     anomalous current through face iFc.
*
*     fchdia - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fchsia(iFc,0) specifies the poloidal component
*     and fchsia(iFc,1) specifies the radial component of the
*     modified diamagnetic current through face iFc.
*
*     fchin - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fchin(iFc,0) specifies the poloidal component
*     and fchin(iFc,1) specifies the radial component of the
*     current contribution from ion-neutral friction through face iFc.
*
*     fch_p - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fch_p(iFc,0) specifies the poloidal component
*     and fch_p(iFc,1) specifies the radial component of the product of
*     the parallel electric current and the poloidal magnetic field component
*     sign through the face iFc.
*
*     fchvisper - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fchvisper(iFc,0) specifies the poloidal component
*     and fchvisper(iFc,1) specifies the radial component of the electric
*     current connected with contribution from the perpendicular viscosity
*     through the face iFc.
*
*     fchvisper - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fchvisper(iFc,0) specifies the poloidal component
*     and fchvisper(iFc,1) specifies the radial component of the electric
*     current connected with contribution from the perpendicular viscosity
*     through the face iFc.
*
*     uadia - (1:nFc,0:1,0:ns-1) real*8 array.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), uadia(iFc,0,is) specifies
*     the total effective drift velocity of species (is) on the (iFc) face
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but projected onto the
*     poloidal direction.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), uadia(iFc,1,is) specifies
*     the total effective drift velocity of species (is) on the (iFc) face
*     in the radial direction.
*
*     vadia - (1:nFc,0:1,0:ns-1) real*8 array.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), vadia(iFc,0,is) specifies
*     the effective diamagnetic drift velocity of species (is) on the (iFc) face
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but projected onto the
*     poloidal direction.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), vadia(iFc,1,is) specifies
*     the effective diamagnetic drift velocity of species (is) on the (iFc) face
*     in the radial direction.
*
*     vedia - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), vedia(iFc,0) specifies
*     the effective diamagnetic drift velocity of electrons on the (iFc) face
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but projected onto the
*     poloidal direction.
*     For (iFc) in (1:nFc), vedia(iFc,1) specifies
*     the effective diamagnetic drift velocity of electrons on the (iFc) face
*     in the radial direction.
*
*     wadia - (1:nFc,0:1,0:ns-1) real*8 array.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), wadia(iFc,0,is) specifies
*     the total diamagnetic drift velocity of species (is) on the (iFc) face
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but projected onto the
*     poloidal direction.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), wadia(iFc,1,is) specifies
*     the total diamagnetic drift velocity of species (is) on the (iFc) face
*     in the radial direction.
*
*     vaecrb - (1:nFc,0:1,0:ns-1) real*8 array.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), vaecrb(iFc,0,is) specifies
*     the ExB drift velocity of species (is) on the (iFc) face
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but projected onto the
*     poloidal direction.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), vaecrb(iFc,1,is) specifies
*     the ExB drift velocity of species (is) on the (iFc) face
*     in the radial direction.
*
*     veecrb - (1:nFc,0:1) real*8 array.
*     For (iFc,is) in (1:nFc,0:1), veecrb(iFc,0) specifies
*     the ExB drift velocity of electrons in the (iFc) cell
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but projected onto the
*     poloidal direction.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), veecrb(iFc,1,is) specifies
*     the ExB drift velocity of electrons in the (iFc) cell
*     in the radial direction.
*
*     kinrgy - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), kinrgy(iCv,is) specifies
*     the kinetic energy of a particle of species (is) in cell (iCv).
*     The velocity used includes the parallel velocity as well
*     as the total diamagnetic and ExB drift velocities in both the radial
*     and diamagnetic directions.
*
*
*   ..residuals and corrections
*
*     resco - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), resco(iCv,is) specifies
*     the residual of the continuity equation on the (iCv) cell for
*     species (is).
*
*     resmo - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), resmo(iCv,is) specifies
*     the residual of the parallel momentum equation on the (iCv)
*     cell for species (is).
*
*     resmb - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), resmb(iCv) specifies the residual
*     of the parallel momentum equation on the (iCv) cell for the
*     selected species. If that species has index (is) then
*     resmb(:).eq.resmo(:,is).
*
*     resmt - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), resmt(iCv) specifies the residual
*     of the total parallel momentum equation on the (iCv) cell.
*
*     reshe - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), reshe(iCv) specifies the residual
*     of the electron heat balance equation on the (iCv) cell.
*
*     reshi - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), reshe(iCv) specifies the residual
*     of the ion heat balance equation on the (iCv) cell.
*
*     reshn - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), reshn(iCv) specifies the residual
*     of the fluid neutral heat balance equation on the (iCv) cell, if
*     a separate energy equation for the hydrogenic neutrals is used.
*     Otherwise, reshn is not used.
*
*     resht - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), resht(iCv) specifies the residual
*     of the total heat balance equation on the (iCv) cell.
*
*     respo - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), respo(iCv) specifies the
*     residual of the electric potential equation on the (iCv) cell.
*
*     corpa - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), corpa(iCv,is) specifies
*     the last computed correction to the pressure of species (is) on
*     the (iCv) cell. corpa is computed by solving an equation that
*     has resco as its right hand side. It is not necessarily the same
*     as the pressure correction over one complete iteration.
*
*     corua - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), corua(iCv,is) specifies
*     the last computed correction to the parallel velocity ua of
*     species (is) on the (iCv) cell.
*     corua is computed by solving an equation that has resmo as its
*     right hand side. It is not necessarily the same as the change in
*     the parallel velocity over one complete iteration.
*
*     corut - (1:nCv) real*8 array.
*     For (iCv) in (0:nx,-1:ny), corut(iCv) specifies the last
*     computed correction to the parallel velocity ua of all species
*     simultaneously on the (iCv) cell.
*     corut is computed by solving an equation that has resmt as its
*     right hand side. It is not necessarily the same as the change in
*     the parallel velocity over one complete iteration.
*
*     corte - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), corte(iCv) specifies the last
*     computed correction to the electron temperature on the (iCv)
*     cell. corte is computed by solving an equation that has reshe as
*     its right hand side. It is not necessarily the same as the change
*     in the electron temperature over one complete iteration.
*
*     corti - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), corti(iCv) specifies the last
*     computed correction to the atom temperature on the (iCv) cell.
*     corti is computed by solving an equation that has reshi as its
*     right hand side. It is not necessarily the same as the change in
*     the atom temperature over one complete iteration.
*
*     cortn - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), cortn(iCv) specifies the last
*     computed correction to the neutral temperature on the (iCv) cell.
*     cortn is computed by solving an equation that has reshn as its
*     right hand side. It is not necessarily the same as the change in
*     the neutral temperature over one complete iteration.
*
*     cortt - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), cortt(iCv) specifies the last
*     computed correction to the combined temperature on the (iCv)
*     cell. cortt is computed by solving an equation that has resht as
*     its right hand side. It is not necessarily the same as the change
*     in the electron or atom temperature over one complete iteration.
*
*     corpo - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), corpo(iCv) specifies the last
*     computed correction to the electric potential on the (iCv)
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
*     areshn - real.
*     areshn specifies the scaled norm of the residual of the fluid neutral
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
*     acortn - real.
*     acortn specifies the norm of the last computed correction for the
*     fluid neutral heat balance equation.
*
*     acorpo - real.
*     acorpo specifies the norm of the last computed correction for the
*     electric potential equation.
*
*   ..rate coefficients
*
*     /b2cmrc/ - common block.
*     /b2cmrc/ contains a table of atomic rate coefficients.
*
*      rtnt - integer constant.
*     rtnt specifies the number of intervals in the discretisation
*     of the temperature range. The table will contain data for
*     (rtnt+1) values of the temperature variable.
*
*      rtnn - integer constant.
*     rtnn specifies the number of intervals in the discretisation
*     of the density range. The table will contain data for (rtnn+1)
*     values of the density variable.
*
*      rtns - integer.
*     rtns specifies the total number of atomic species; i.e.
*       (sum inuc : 0.le.inuc.lt.nnuc : izhi(inuc)-izlo(inuc)+1) ,
*     where inuc is a nuclear species index, izlo(inuc) is the lowest
*     and izhi(inuc) is the highest atomic charge considered in B2.5 for
*     that nuclear species.
*     rtns equals the variable ns used elsewhere in B2.5.
*
*      rtdumm - integer.
*     rtdumm is a dummy integer variable, introduced in order to obtain
*     an even-length sequence of integer variables in the common block.
*
*      rtzmin, rtzmax, rtzn - (0:rtns-1) real*8 array.
*     For 0.le.is.lt.rtns, rtzmin(is) and rtzmax(is) specify the atomic
*     charge range, and rtzn(is) specifies the nuclear charge
*     of atomic species (is).
*     The izhi(inuc)-izlo(inuc)+1 charge states for each nuclear
*     species are taken in natural order, starting at izlo(inuc); the
*     highest state of nuclear species (inuc) is followed by the lowest
*     state of (inuc+1). rtzmin, rtzmax and rtzn are equal to the arrays
*     zamin, zamax and zn used elsewhere in B2.5.
*
*      rtt - (0:rtnt) real*8 array.
*     rtt(0:rtnt) specifies the sequence of values of the temperature
*     variable, expressed in ev, for which the rate coefficients are
*     evaluated. rtt(0:rtnt) will be positive and increasing.
*
*      rtn - (0:rtnn) real*8 array.
*     rtn(0:rtnn) specifies the sequence of values of the density
*     variable, expressed in m**-3, for which the rate coefficients
*     are evaluated. rtn(0:rtnn) will be positive and increasing.
*
*      rtlt - (0:rtnt) real*8 array.
*     rtlt(0:rtnt) holds log(rtt(0:rtnt)).
*
*      rtln - (0:rtnn) real*8 array.
*     rtln(0:rtnn) holds log(rtn(0:rtnn)).
*
*      rtlsa - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlsa(i,j,k) specifies
*     logu(rtrsa), where rtrsa is the ionisation rate coefficient at
*     electron temperature rtt(i), electron density rtn(j), for the
*     process k->k+1. rtrsa has dimension m**3/s. The function logu is
*     the same as (log) except that logu(0.0) is equal to some huge
*     negative number.
*
*      rtlra - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlra(i,j,k) specifies
*     logu(rtrra), where rtrra is the recombination rate coefficient at
*     electron temperature rtt(i) and electron density rtn(j), for the
*     process k->k-1. rtrra has dimension m**3/s.
*
*      rtlqa - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlqa(i,j,k) specifies
*     logu(rtrqa), where rtrqa is the electron heat loss rate
*     coefficient at electron temperature rtt(i) and electron density
*     rtn(j), for processes starting from atomic species k. rtrqa has
*     dimension ev*m**3/s.
*
*      rtlcx - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlcx(i,j,k) specifies
*     logu(rtrcx), where rtrcx is the charge exchange rate coefficient
*     at neutral hydrogen temperature/amu rtt(i) and electron density
*     rtn(j), for the process k->k-1. rtrcx has dimension m**3/s.
*
*      rtlza - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlza(i,j,k) specifies
*     the effective charge state at
*     electron temperature rtt(i), electron density rtn(j), for the
*     superstage k. rtlza has dimension e.
*
*      rtlz2 - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlz2(i,j,k) specifies
*     the effective squared charge at
*     electron temperature rtt(i), electron density rtn(j), for the
*     superstage k. rtlz2 has dimension e**2.
*
*      rtlpt - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlpt(i,j,k) specifies
*     the cumulative ionisation energy at
*     electron temperature rtt(i), electron density rtn(j),
*     to reach superstage k. rtlpt has dimension eV.
*
*      rtlpi - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlpi(i,j,k) specifies
*     the effective ionisation energy at
*     electron temperature rtt(i), electron density rtn(j),
*     to reach superstage k from superstage k-1. rtlpi has dimension eV.
*
*     nscx - integer.
*     nscx is the number of charge exchanging neutral hydrogen species
*     that are present in the calculation.  (So, usually nscx.eq.1).
*     It will hold that 0.le.nscx.lt.10.
*
*     iscx - (0:nscxmax-1) integer array.
*     iscx(0:nscx-1) will hold the species indices for the various kinds
*     of neutral hydrogen that are present in the calculation.
*     For 0.le.k.lt.nscx it will hold that 0.le.iscx(k).lt.ns.
*
*     rlsa - (1:nCv,0:1,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), rlsa(iCv,0,is) and
*     rlsa(iCv,1,is) are the coefficients in the representation
*       rsa = exp(rlsa(:,0,)+rlsa(:,1,)*log(te(iCv)/ev))
*     that specifies the rate coefficient for ionisation from atomic
*     species (is) to species (is+1). (For is.eq.ns-1 or in case the
*     atomic species (is) and (is+1) belong to different nuclear
*     species this expression will underflow, and its value is then
*     understood as 0.)
*     rsa has dimension m**3/s.
*
*     rlra - (1:nCv,0:1,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), rlra(iCv,0,is) and
*     rlra(iCv,1,is) are the coefficients in the representation
*       rra = exp(rlra(:,0,)+rlra(:,1,)*log(te(iCv)/ev))
*     that specifies the rate coefficient for recombination from atomic
*     species (is) to species (is-1). (For is.eq.0 or in case the
*     atomic species (is) and (is-1) belong to different nuclear
*     species this expression will underflow, and its value is then
*     understood as 0.)
*     rra has dimension m**3/s.
*
*     rlqa - (1:nCv,0:1,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), rlqa(iCv,0,is) and
*     rlqa(iCv,1,is) are the coefficients in the representation
*       rqa = ev*exp(rlqa(:,0,)+rlqa(:,1,)*log(te(iCv)/ev))
*     that specifies the rate coefficient for electron heat loss terms
*     that are proportional to the density of the atomic species (is).
*     rqa has dimension J*m**3/s.
*
*     rlcx - (1:nCv,0:1,0:ns-1,0:nscxmax-1) real*8 array.
*     For (iCv,is,k) in (1:nCv,0:ns-1,0:nscx-1),
*     rlcx(iCv,0,is,k) and rlcx(iCv,1,is,k) are the coefficients in
*     the expression
*       rcx = exp(rlcx(:,0,,)+rlcx(:,1,,)*log(th0/(amh0*ev)))
*     that specifies the rate coefficient for charge exchange
*     recombination from state (is) to state (is-1) by collisions with
*     the neutral hydrogen species iscx(k); here, amh0 is the atomic
*     mass of that neutral species. (For is.eq.0 or in case the atomic
*     species (is) and (is-1) belong to different nuclear species this
*     expression will underflow, and its value is then understood as 0.)
*     rcx has dimension m**3/s.
*
*     rlza - (1:nCv,0:1,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), rlza(iCv,0,is) and
*     rlza(iCv,1,is) are the coefficients in the representation
*       rza = rlza(,0,)+rlza(,1,)*log(te(iCv)/ev)
*     that specifies the effective charge of stage (is).
*     rza has dimension e.
*
*     rlz2 - (1:nCv,0:1,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), rlz2(iCv,0,is) and
*     rlz2(iCv,1,is) are the coefficients in the representation
*       rz2 = rlz2(,0,)+rlz2(,1,)*log(te(iCv)/ev)
*     that specifies the effective square charge of stage (is).
*     rz2 has dimension e**2.
*
*     rlpt - (1:nCv,0:1,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), rlpt(iCv,0,is) and
*     rlpt(iCv,1,is) are the coefficients in the representation
*       rpt = rlpt(,0,)+rlpt(,1,)*log(te(iCv)/ev)
*     that specifies the cumulative ionisation potential to stage (is).
*     rpt has dimension eV.
*
*     rlpi - (1:nCv,0:1,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), rlpi(iCv,0,is) and
*     rlpi(iCv,1,is) are the coefficients in the representation
*       rpi = rlpi(,0,)+rlpi(,1,)*log(te(iCv)/ev)
*     that specifies the effective ionisation potential to stage (is)
*     from stage (is-1).
*     rpi has dimension eV.
*
*     rsa - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), rsa(iCv,is) specifies
*     the rate coefficient for ionisation from atomic species (is) to
*     species (is+1). For is.eq.ns-1 or in case the atomic species (is)
*     and (is+1) belong to different nuclear species, this coefficient
*     will be 0.
*     rsa has dimension m**3/s.
*
*     rra - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), rra(iCv,is) specifies
*     the rate coefficient for recombination from atomic species (is)
*     to species (is-1). For is.eq.0 or in case the atomic species (is)
*     and (is-1) belong to different nuclear species, this coefficient
*     will be 0.
*     rra has dimension m**3/s.
*
*     rqa - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), rqa(iCv,is) specifies
*     the rate coefficient for electron heat loss terms that are
*     proportional to the density of the atomic species (is).
*     rqa has dimension J*m**3/s.
*
*     rcx - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), rcx(iCv,is) specifies
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
*     sna - (1:nCv,0:1,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), sna(iCv,0,is) and
*     sna(iCv,1,is) are the coefficients in the expression
*       sna(iCv,0,is)+sna(iCv,1,is)*na(iCv,is)
*     that specifies a linear approximation to the assumed particle
*     source of species (is) integrated over the (iCv) cell.
*
*     sne - (1:nCv,0:1) real*8 array.
*     For (iCv) in (1:nCv), sne(iCv,0) and
*     sne(iCv,1) are the coefficients in the expression
*       sne(iCv,0)+sne(iCv,1)*ne(iCv)
*     that specifies a linear approximation to the assumed non-ambipolar
*     electron source of integrated over the (iCv) cell.
*
*     smo - (1:nCv,0:3,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), smo(iCv,0:3,is) are the
*     coefficients in the expression
*       smo(iCv,0,is)+
*         smo(iCv,1,is)*ua(iCv,is)+
*         smo(iCv,2,is)*roxa(iCv,is)+
*         smo(iCv,3,is)*roxa(iCv,is)*ua(iCv,is)
*     that specifies a linear approximation to the parallel momentum
*     source for species (is) integrated over the (iCv) cell; here
*       roxa(iCv,is) = am(is)*mp*na(iCv,is) .
*
*     she - (1:nCv,0:3) real*8 array.
*     For (iCv) in (1:nCv), she(iCv,0:3) are the coefficients
*     in the expression
*       she(iCv,0)+
*         she(iCv,1)*te(iCv)+
*         she(iCv,2)*ne(iCv)+
*         she(iCv,3)*ne(iCv)*te(iCv)
*     that specifies a linear approximation to the assumed electron heat
*     source integrated over the (iCv) cell.
*
*     shi - (1:nCv,0:3) real*8 array.
*     For (iCv) in (1:nCv), shi(iCv,0:3) are the coefficients
*     in the expression
*       shi(iCv,0)+
*         shi(iCv,1)*ti(iCv)+
*         shi(iCv,2)*ni(iCv,0)+
*         shi(iCv,3)*ni(iCv,0)*ti(iCv)
*     that specifies a linear approximation to the assumed all atom heat
*     source integrated over the (iCv) cell.
*     Note: for now, shi(:,2) and shi(:,3) always scale with ni(:,0)
*     (all ions + fluid neutrals), even when tn_style.eq.2. It can be considered
*     for the future to let shi scale with ni(:,1) when tn_style.ne.0, but this
*     requires many code changes.
*
*     sch - (1:nCv,0:3) real*8 array.
*     For (iCv) in (1:nCv), sch(iCv,0:3) are the coefficients
*     in the expression
*       sch(iCv,0)+
*         sch(iCv,1)*po(iCv)+
*         sch(iCv,2)*ne(iCv)+
*         sch(iCv,3)*ne(iCv)*po(iCv)
*     that specifies a linear approximation to the assumed electric
*     charge source integrated over the (iCv) cell.
*
*   ..partial source coefficients
*
*     sna0 - (1:nCv,0:1,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), sna0(iCv,0,is) and
*     sna0(iCv,1,is) are the coefficients in the expression
*       sna0(iCv,0,is)+sna0(iCv,1,is)*na(iCv,is)
*     that specifies a linear approximation to some specific
*     contribution to the particle source of species (is) integrated
*     over the (iCv) cell.
*
*     sne0 - (1:nCv,0:1) real*8 array.
*     For (iCv) in (1:nCv), sne0(iCv,0) and
*     sne0(iCv,1) are the coefficients in the expression
*       sne0(iCv,0)+sne0(iCv,1)*ne(iCv)
*     that specifies a linear approximation to some specific
*     contribution to the particle non-ambipolar electron source
*     integrated over the (iCv) cell.
*
*     smo0 - (1:nCv,0:3,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), smo0(iCv,0:3,is) are the
*     coefficients in the expression
*       smo0(iCv,0,is)+
*         smo0(iCv,1,is)*ua(iCv,is)+
*         smo0(iCv,2,is)*roxa(iCv,is)+
*         smo0(iCv,3,is)*roxa(iCv,is)*ua(iCv,is)
*     that specifies a linear approximation to some specific
*     contribution to the non-frictional parallel momentum source for
*     species (is) integrated over the (iCv) cell. Here,
*       roxa(iCv,is) = am(is)*mp*na(iCv,is) .
*
*     she0 - (1:nCv,0:3) real*8 array.
*     For (iCv) in (1:nCv), she0(iCv,0:3) are the coefficients
*     in the expression
*       she0(iCv,0)+
*         she0(iCv,1)*te(iCv)+
*         she0(iCv,2)*ne(iCv)+
*         she0(iCv,3)*ne(iCv)*te(iCv)
*     that specifies a linear approximation to some specific
*     contribution to the electron heat source integrated over the
*     (iCv) cell.
*
*     shi0 - (1:nCv,0:3) real*8 array.
*     For (iCv) in (1:nCv), shi0(ix,iy,0:3) are the coefficients
*     in the expression
*       shi0(ix,iy,0)+
*         shi0(iCv,1)*ti(iCv)+
*         shi0(iCv,2)*ni(iCv)+
*         shi0(iCv,3)*ni(iCv)*ti(iCv)
*     that specifies a linear approximation to some specific
*     contribution to the all atom heat source integrated over the
*     (iCv) cell.
*
*     sch0 - (1:nCv,0:3) real*8 array.
*     For (iCv) in (1:nCv), sch0(ix,iy,0:3) are the coefficients
*     in the expression
*       sch0(iCv,0)+
*         sch0(iCv,1)*po(iCv)+
*         sch0(iCv,2)*ne(iCv)+
*         sch0(iCv,3)*ne(iCv)*po(iCv)
*     that specifies a linear approximation to some specific
*     contribution to the electric charge source integrated over the
*     (iCv) cell.
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
*         cfsig(3)*te(,)/ev/abs(bb(,,3))/16+
*         cfsig(4)*exp(-cfsig(5)*(y/L)**2)
*     that specifies the diffusivity corresponding to the anomalous
*     electrical conductivity, which is
*       sig(,) = qe*ne(,)*dfsig(,).
*     L is the local radial width of the plasma and
*     y is the radial distance from the inner (South) boundary.
*     If cfsig(7) is not zero, the sig0 array is flux-scaled with
*     a scaling factor of cfsig(7)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfsig(2) and (6) are not used in the default version of the code.
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
*     cflim(3)=cflal specifies the flux limit factor for the
*     thermo-electric coefficient.
*     cflim(4)=cflab specifies the flux limit factor for the
*     friction force.
*     Zero settings turn off these flux limits.
*     cflim(5:7) are not used in the default version of the code.
*
*   ..transport coefficients, differential
*
*     dna - (1:nCv,0:1,0:ns-1) real*8 array.
*     (Let dnax(,,)=dna(,,0,) and dnay(,,)=dna(,,1,).)
*     For (iCv,is) in (1:nCv,0:ns-1), dnax(iCv,is) and
*     dnay(iCv,is) specify the x-component and the y-component of
*     the diffusivity of species (is) on the (iCv) cell. The
*     corresponding term in the particle flux density is
*     -dna*grad.na .
*
*     dpa - (1:nCv,0:1,0:ns-1) real*8 array.
*     (Let dpax(,,)=dpa(,,0,) and dpay(,,)=dpa(,,1,).)
*     For (iCv,is) in (1:nCv,0:ns-1), dpax(iCv,is) and
*     dpay(iCv,is) specify the x-component and the y-component of
*     the diffusion coefficient (w.r.t. grad.pa) of species (is) on
*     the (iCv) cell. The corresponding term in the particle flux
*     density is  -dpa*grad.pa .
*
*     vla - (1:nCv,0:1,0:ns-1) real*8 array.
*     (Let vlax(,,)=vla(,,0,) and vlay(,,)=vla(,,1,).)
*     For (iCv,is) in (1:nCv,0:ns-1), vlax(iCv,is) specifies
*     the x-component of the anomalous velocity of species (is) on the
*     (iCv) cell.
*     For (iCv,is) in (1:nCv,0:ns-1), vlay(iCv,is) specifies
*     the y-component of the anomalous velocity of species (is) on the
*     (iCv) cell.
*
*     vsa - (1:nCv,0:1,0:ns-1) real*8 array.
*     (Let vsax(,,)=vsa(,,0,) and vsay(,,)=vsa(,,1,).)
*     For (iCv,is) in (1:nCv,0:ns-1), vsax(iCv,is) specifies
*     the x-component and vsay(iCv,is) the y-component of the
*     viscosity of species (is) on the (iCv) cell.
*
*     hce - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let hcex(,)=hce(,,0) and hcey(,)=hce(,,1).)
*     For (iCv) in (1:nCv), hcex(iCv) specifies the
*     x-component and hcey(iCv) the y-component of the electron
*     thermal conductivity on the (iCv) cell.
*
*     hci - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let hcix(,)=hci(,,0) and hciy(,)=hci(,,1).)
*     For (iCv) in (1:nCv), hcix(iCv) specifies the
*     x-component and hciy(iCv) the y-component of the all atom
*     thermal conductivity on the (iCv) cell.
*
*     sig - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let sigx(,)=sig(,,0) and sigy(,)=sig(,,1).)
*     For (iCv) in (1:nCv), sigx(iCv) specifies the
*     x-component and sigy(iCv) the y-component of the electrical
*     conductivity on the (iCv) cell.
*
*     alf - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let alfx(,)=alf(,,0) and alfy(,)=alf(,,1).)
*     For (iCv) in (1:nCv), alfx(iCv) specifies the
*     x-component and alfy(iCv) the y-component of the thermo-
*     electric coefficient on the (iCv) cell.
*
*     dna0 - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), dna0(iCv,is) specifies
*     the anomalous diffusivity with respect to grad(na) of species
*     (is) on the (iCv) cell.
*
*     dpa0 - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), dpa0(iCv,is) specifies
*     the anomalous diffusivity with respect to grad(pa) of species
*     (is) on the (iCv) cell.
*
*     vla0 - (1:nCv,0:1,0:ns-1) real*8 array.
*     (Let vla0x(,,)=vla0(,,0,) and vla0y(,,)=vla0(,,1,).)
*     For (iCv,is) in (1:nCv,0:ns-1), vla0x(iCv,is) specifies
*     the x-component and vla0y(iCv,is) the y-component of the
*     anomalous velocity of species (is) on the (iCv) cell.
*
*     vsa0 - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), vsa0(iCv,is) specifies
*     the anomalous viscosity of species (is) on the (iCv) cell.
*
*     hce0 - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), hce0(iCv) specifies the anomalous
*     electron thermal conductivity on the (iCv) cell.
*
*     hci0 - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), hci0(iCv) specifies the anomalous
*     atom thermal conductivity on the (iCv) cell.
*
*     sig0 - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), sig0(iCv) specifies the anomalous
*     electrical conductivity on the (iCv) cell.
*
*     alf0 - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), alf0(iCv) specifies the anomalous
*     thermo-electric coefficient on the (iCv) cell.
*
*   ..transport coefficients, discrete
*
*     ceqp - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), ceqp(iCv) specifies the energy
*     equipartition multiplier integrated over the (iCv) cell.
*     (The cell-integrated energy equipartition coefficient is
*     ceqp*ne*ni/te**1.5.)
*
*     cdna - (1:nCv,0:1,0:ns-1) real*8 array.
*     (Let cdnax=cdna(,,0,) and cdnay=cdna(,,1,).)
*     For (iCv,is) in (1:nCv,0:ns-1), cdnax(iCv,is) specifies
*     the particle diffusivity w.r.t. grad(na) of species (is) between
*     the (iCv) cell and its left neighbor. This coefficient enters into
*     the expression for the particle flux in the form
*       -cdnax(iCv,is)*(na(iCv,is)-na(ix-1,iy,is))
*     For ix.eq.-1, cdnax(ix,-1:ny,0:ns-1) will hold 0.
*     For (iCv,is) in (-1:nx,0:ny,0:ns-1), cdnay(iCv,is) specifies
*     the particle diffusivity w.r.t. grad(na) of species (is) between
*     the (iCv) cell and its bottom neighbor. This coefficient enters into
*     the expression for the particle flux in the form
*       -cdnay(iCv,is)*(na(iCv,is)-na(bottom(iCv),is))
*     For iy.eq.-1, cdnay(-1:nx,iy,0:ns-1) will hold 0.
*
*     cdpa - (1:nCv,0:1,0:ns-1) real*8 array.
*     (Let cdpax=cdpa(,,0,) and cdpay=cdpa(,,1,).)
*     For (iCv,is) in (1:nCv,0:ns-1), cdpax(iCv,is) specifies
*     the particle diffusivity w.r.t. grad(pa) of species (is) between
*     the (iCv) cell and its left neighbor. This coefficient enters into
*     the expression for the particle flux in the form
*       -cdpax(iCv,is)*(pa(iCv,is)-pa(ix-1,iy,is))
*     For ix.eq.-1, cdpax(ix,-1:ny,0:ns-1) will hold 0.
*     For (iCv,is) in (-1:nx,0:ny,0:ns-1), cdpay(iCv,is) specifies
*     the particle diffusivity w.r.t. grad(pa) of species (is) between
*     the (iCv) cell and its bottom neighbor. This coefficient enters into
*     the expression for the particle flux in the form
*       -cdpay(iCv,is)*(pa(iCv,is)-pa(bottom(iCv),is))
*     For iy.eq.-1, cdpay(-1:nx,iy,0:ns-1) will hold 0.
*
*     cvla - (1:nCv,0:1,0:ns-1) real*8 array.
*     (Let cvlax=cvla(,,0,) and cvlay=cvla(,,1,).)
*     For (iCv,is) in (1:nCv,0:ns-1), cvlax(iCv,is) specifies
*     the surface-integrated anomalous velocity of species (is) between
*     the (iCv) cell and its left neighbor. This coefficient enters into
*     the expression for the particle flux in the form
*       cvlax(iCv,is)*(na)
*     in which (na) involves some interpolation to cell faces.
*     For ix.eq.-1, cvlax(ix,-1:ny,0:ns-1) will hold 0.
*     For (iCv,is) in (-1:nx,0:ny,0:ns-1), cvlay(iCv,is) specifies
*     the surface-integrated anomalous velocity of species (is) between
*     the (iCv) cell and its bottom neighbor. This coefficient enters into
*     the expression for the particle flux in the form
*       cvlay(iCv,is)*(na)
*     in which (na) involves some interpolation to cell faces.
*     For iy.eq.-1, cvlay(-1:nx,iy,0:ns-1) will hold 0.
*
*     cvsa - (1:nCv,0:1,0:ns-1) real*8 array.
*     (Let cvsax=cvsa(,,0,) and cvsay=cvsa(,,1,).)
*     For (iCv,is) in (1:nCv,0:ns-1), cvsax(iCv,is) specifies
*     the viscous conductance of species (is) between the (iCv) cell
*     and its left neighbor.
*     For ix.eq.-1, cvsax(ix,-1:ny,0:ns-1) will hold 0.
*     For (iCv,is) in (-1:nx,0:ny,0:ns-1), cvsay(iCv,is) specifies
*     the viscous conductance of species (is) between the (iCv) cell
*     and its bottom neighbor.
*     For iy.eq.-1, cvsay(-1:nx,iy,0:ns-1) will hold 0.
*
*     chce - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let chcex=chce(,,0) and chcey=chce(,,1).)
*     For (iCv) in (0:nx,-1:ny), chcex(iCv) specifies the electron
*     thermal conductance between the (iCv) cell and its left neighbor.
*     This coefficient enters into the expression for the electron
*     heat flux in the form
*       -chcex(iCv)*(te(iCv)-te(left(iCv))) .
*     For ix.eq.-1, chcex(ix,-1:ny) will hold 0.
*     For (iCv) in (-1:nx,0:ny), chcey(iCv) specifies the electron
*     thermal conductance between the (iCv) cell and its bottom neighbor.
*     This coefficient enters into the expression for the electron
*     heat flux in the form
*       -chcey(iCv)*(te(iCv)-te(bottom(iCv))) .
*     For iy.eq.-1, chcey(-1:nx,iy) will hold 0.
*
*     chve - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let chvex=chve(,,0) and chvey=chve(,,1).)
*     For (iCv) in (0:nx,-1:ny), chvex(iCv) specifies the electron
*     thermal strange velocity between the (iCv) cell and its left
*     neighbor. This coefficient enters into the expression for the
*     electron heat flux in the form  chvex*ne*te , interpolated to
*     cell faces in a suitable manner.
*     For ix.eq.-1, chvex(ix,-1:ny) will hold 0.
*     For (iCv) in (-1:nx,0:ny), chvey(iCv) specifies the electron
*     thermal strange velocity between the (iCv) cell and its bottom
*     neighbor. This coefficient enters into the expression for the
*     electron heat flux in the form  chvey*ne*te , interpolated to
*     cell faces in a suitable manner.
*     For iy.eq.-1, chvey(-1:nx,iy) will hold 0.
*
*     chci - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let chcix=chci(,,0) and chciy=chci(,,1).)
*     For (iCv) in (0:nx,-1:ny), chcix(iCv) specifies the all atom
*     thermal conductance between the (iCv) cell and its left neighbor.
*     This coefficient enters into the expression for the atom heat
*     flux in the form
*       -chcix(iCv)*(ti(iCv)-ti(left(iCv))) .
*     For ix.eq.-1, chcix(ix,-1:ny) will hold 0.
*     For (iCv) in (-1:nx,0:ny), chciy(iCv) specifies the all atom
*     thermal conductance between the (iCv) cell and its bottom neighbor.
*     This coefficient enters into the expression for the atom heat
*     flux in the form
*       -chciy(iCv)*(ti(iCv)-ti(bottom(iCv))) .
*     For iy.eq.-1, chciy(-1:nx,iy) will hold 0.
*
*     chvi - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let chvix=chvi(,,0) and chviy=chvi(,,1).)
*     For (iCv) in (0:nx,-1:ny), chvix(iCv) specifies the all atom
*     thermal strange velocity between the (iCv) cell and its left
*     neighbor. This coefficient enters into the expression for the
*     atom heat flux in the form  chvix*ni*ti , interpolated to cell
*     faces in a suitable manner.
*     For ix.eq.-1, chvix(ix,-1:ny) will hold 0.
*     For (iCv) in (-1:nx,0:ny), chviy(iCv) specifies the all atom
*     thermal strange velocity between the (iCv) cell and its bottom
*     neighbor. This coefficient enters into the expression for the
*     atom heat flux in the form  chviy*ni*ti , interpolated to
*     cell faces in a suitable manner.
*     For iy.eq.-1, chviy(-1:nx,iy) will hold 0.
*
*     csig - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let csigx=csig(,,0) and csigy=csig(,,1).)
*     For (iCv) in (0:nx,-1:ny), csigx(iCv) specifies the electrical
*     conductance between the (iCv) cell and its left neighbor.
*     For ix.eq.-1, csigx(ix,-1:ny) will hold 0.
*     For (iCv) in (-1:nx,0:ny), csigy(iCv) specifies the electrical
*     conductance between the (iCv) cell and its bottom neighbor.
*     For iy.eq.-1, csigx(-1:nx,iy) will hold 0.
*
*     calf - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let calfx=calf(,,0) and calfy=calf(,,1).)
*     For (iCv) in (0:nx,-1:ny), calfx(iCv) specifies the thermo-
*     electric coefficient between the (iCv) cell and its left neighbor.
*     For ix.eq.-1, calfx(ix,-1:ny) will hold 0.
*     For (iCv) in (-1:nx,0:ny), calfy(iCv) specifies the thermo-
*     electric coefficient between the (iCv) cell and its bottom neighbor.
*     For iy.eq.-1, calfy(-1:nx,iy) will hold 0.
*
*     csigin - (-1:nx,-1:ny,0:1,0:ns-1,0:ns-1) real*8 array.		 !xpb
*     (Let csiginx=csigin(,,0,,) and csiginy=csigin(,,1,,).)
*     For (ix,iy,is,is0) in (0:nx,-1:ny,0:ns-1,0:ns-1),
*     csiginx(ix,iy,is,is0) specifies the ion-neutral
*     conductivity between the (iCv) cell and its left neighbor, where
*     the ion species is species is and the neutral species is is0.
*     For ix.eq.-1, csiginx(ix,-1:ny) will hold 0.
*     For (iCv) in (-1:nx,0:ny), csiginy(iCv) specifies the ion-neutral
*     conductivity between the (iCv) cell and its bottom neighbor, where
*     the ion species is species is and the neutral species is is0.
*     For iy.eq.-1, csiginx(-1:nx,iy) will hold 0.
*     When is is a neutral species or is0 is an ionised species,
*     csigin(,,,is,is0) will be zero.
*
*     csig_cl - (-1:nx,-1:ny,0:1) real*8 array.				 !xpb
*     (Let csig_clx=csig_cl(,,0) and csig_cly=csig_cl(,,1).)
*     For (iCv) in (0:nx,-1:ny), csig_clx(iCv) specifies the classical
*     conductivity between the (iCv) cell and its left neighbor.
*     For ix.eq.-1, csig_clx(ix,-1:ny) will hold 0.
*     For (iCv) in (-1:nx,0:ny), csig_cly(iCv) specifies the classical
*     conductivity between the (iCv) cell and its bottom neighbor.
*     For iy.eq.-1, csig_clx(-1:nx,iy) will hold 0.
*
*     calf_cl - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let calf_clx=calf_cl(,,0) and calf_cly=calf_cl(,,1).)
*     For (iCv) in (0:nx,-1:ny), calf_clx(iCv) specifies the
*     classical part of the thermo-electric coefficient
*     between the (iCv) cell and its left neighbor.
*     For ix.eq.-1, calf_clx(ix,-1:ny) will hold 0.
*     For (iCv) in (-1:nx,0:ny), calf_cly(iCv) specifies the
*     classical part of the thermo-electric coefficient
*     between the (iCv) cell and its bottom neighbor.
*     For iy.eq.-1, calf_cly(-1:nx,iy) will hold 0.
*
*     csig_an - (-1:nx,-1:ny,0:1) real*8 array.				 !xpb
*     (Let csig_anx=csig_an(,,0) and csig_any=csig_an(,,1).)
*     For (iCv) in (0:nx,-1:ny), csig_anx(iCv) specifies the anomalous
*     conductivity between the (iCv) cell and its left neighbor.
*     For ix.eq.-1, csig_anx(ix,-1:ny) will hold 0.
*     For (iCv) in (-1:nx,0:ny), csig_any(iCv) specifies the anomalous
*     conductivity between the (iCv) cell and its bottom neighbor.
*     For iy.eq.-1, csig_anx(-1:nx,iy) will hold 0.
*
*     calf_an - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let calf_anx=calf_an(,,0) and calf_any=calf_an(,,1).)
*     For (iCv) in (0:nx,-1:ny), calf_anx(iCv) specifies the
*     anomalous part of the thermo-electric coefficient
*     between the (iCv) cell and its left neighbor.
*     For ix.eq.-1, calf_anx(ix,-1:ny) will hold 0.
*     For (iCv) in (-1:nx,0:ny), calf_any(iCv) specifies the
*     anomalous part of the thermo-electric coefficient
*     between the (iCv) cell and its bottom neighbor.
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
