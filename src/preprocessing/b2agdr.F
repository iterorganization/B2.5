*-----------------------------------------------------------------------
*.specification

      subroutine b2agdr (ninp, nout, nx, ny, nx1, ny1)
      use b2mod_types
      use b2mod_geo
      use b2mod_geo2
      use b2mod_version
      use b2mod_indirect
      use b2mod_b2cmfs
      use b2mod_b2cmgs
      use b2mod_connectivity
      use b2ag_parameters
      use b2mod_interp
      use b2mod_grid_mapping
      use b2us_geo
      use b2us_map
      use b2us_prep
      use b2us_io
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:1), nout(0:1), nx, ny, nx1, ny1, nout_us, io
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AGDR computes and writes to an external file the specification
*     of a b2 geometry and magnetic field.
*
*
*  3. description (see also routine b2cdca)
*
*     This subroutine computes and writes out a specification of a b2
*     geometry and magnetic field. Summary information will be written
*     to unit nout(0), which must be connected to a formatted file.
*     A complete specification of the geometry and the magnetic field
*     will be written to unit nout(1), which must be connected to an
*     un*formatted file.
*
*     B2 employs a two-dimensional, curvilinear, topologically
*     rectangular mesh. The coordinates of the mesh are labelled 'x'
*     and 'y'. There are nx+2 cells along the x-coordinate and ny+2
*     cells along the y-coordinate. The cells of the mesh are indexed
*     by (ix,iy) in (-1:nx,-1:ny); the interior cells have index in
*     (0:nx-1,0:ny-1). The mesh may contain cuts, which typically
*     originate from an x-point in the real*8 geometry.
*
*     The primary data describing the mesh are the coordinates of the
*     four vertices of all the cells of the mesh. These data are stored
*     in the arrays crx and cry. The metric parameters (vol, hx, hy, qz,
*     qc and gs; described elsewhere) are derived from crx and cry.
*
*     The primary data describing the magnetic field are the flux
*     functions fpsi and ffbz at the cell vertices. The magnetic field
*     components at cell centers are derived from fpsi and ffbz.
*
*     Singular regions in the geometry are specified by the data in the
*     commmon block /b2cmgs/.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:1) integer array, input.
*     ninp specifies the input unit numbers, as follows.
*     ninp(0): formatted; provides geometric parameters.
*     ninp(1): formatted; provides mhd equilibrium data.
*
*     nout - (0:1) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides summary printed output.
*     nout(1): un*formatted; full specification of the geometry and
*      field.
*
*     nx, ny - integer, input.
*     nx and ny specify the number of interior cells along the first
*     and the second coordinate, respectively. The total number of
*     cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
*     It will hold that 0.le.nx and 0.le.ny.
*
*     nx1, ny1 - integer, input.
*     nx1 and ny1 specify the mesh dimensions for an auxiliary mesh
*     that will be used to evaluate accurately the metric coefficients.
*     This auxiliary mesh has no role outside the present program.
*     nx1 will be a multiple of nx and ny1 will be a multiple of ny.
*     The multiplicity test is put off till b2agmx, so we know if there
*     are any internal isolating cuts inside the domain.
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks
*   ..local variables
      integer idum(0:9), ix, iy, geometryType
      integer nxiso1, nxiso2, nyiso1, nyiso2
      character lblgm*120, idout*80, id*8, cnamip*80, cvalip*80
      character line*500
      integer istyle
      type(mapping) :: mpg
      type(geometry) :: gm
*   ..initialisation
      data nxiso1/-2/, nxiso2/-2/, nyiso1/-2/, nyiso2/-2/
*   ..procedures
      logical streql, is_comment
      external xertst, xerrab, streql, cfwuin, cfwuch, cfwure
      external ipgeti, intfaceh1, intfacer1, intfacev1, is_comment
      external b2xxid, b2agsi, b2agfs, b2agcf, b2agmt, b2agmx, b2agbb,
     &  b2xgbs1, b2wfgm, b2wugm

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The workspace array w is used as a stack for storage of dynamic
*     real*8 arrays. An array 'name' in w(0:*) will be stored starting at
*     location w(l'name').
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     /b2cmfs/ is reserved for numerical mhd equilibrium data.
*     (Presently empty--to be used in conjunction with grid generation
*     based on a numerical mhd equilibrium.)
*
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
*     Description of some local variables (including arrays in w(0:*):
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
*     hx, hy - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), hx(ix,iy) and hy(ix,iy) specify the
*     x-diameter and the y-diameter, respectively, of the (ix,iy) cell.
*     It will hold that 0.lt.hx(,) and 0.lt.hy(,).
*
*     cell_width, cell_height - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), cell_width(ix,iy) and cell_height(ix,iy)
*     specify the average width and height, respectively, of the (ix,iy) cell.
*     It will hold that 0.lt.cell_width(,) and 0.lt.cell_height(,).
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
*     face between the (ix,iy) cell and its left neighbor. For ix.eq.-1,
*     gsx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), gsy(ix,iy) specifies the area of the
*     face between the (ix,iy) cell and its bottom neighbor. For iy.eq.-1,
*     gsy(-1:nx,iy) will hold 0.
*     For (ix,iy) in (-1:nx,-1:ny), gsz(ix,iy) specifies the area of the
*     (ix,iy) cell normal to the third (ignorable) coordinate.
*
*     pbs - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let pbsx(,)=pbs(,,0) and pbsy(,)=pbs(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), pbsx(ix,iy) specifies the product
*     (bx/bb)*sx and pbsy(ix,iy) specifies (by/bb)*sy, or (magnetic 
*     field pitch)*(surface area), namely the parallel contact area, 
*     on the (ix,iy) cell left and bottom faces respectively.  
*     See routine b2xgbs for the form of interpolation employed.
*
*     crx1, cry1, fpsi1, ffbz1, vol1, 
*     hx1, hy1, qz1, qc1, gs1, pbs1 - real*8 arrays.
*     These arrays characterise the geometry on the auxiliary mesh.
*
*     lblgm - character*120.
*     lblgm identifies the instance of the program employed to prepare
*     the geometry.
*
*     idout - character*80.
*     idout specifies the desired selection of output segments for
*     printed and graphical output. Associated with each segment is a
*     four-character word, which is identified on the printed output.
*     (The same word may be associated with several segments.)
*     A particular segment will be produced only if its associated word
*     has a match, as defined by the routine strmas, in idout.
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
*     parg - (0:99) real*8 array.
*     parg(0:*) contains a set of numerical parameters that specify
*     details of the geometry. The interpretation of these parameters
*     is problem-dependent; see the code.
*
*     periodic_bc - integer
*     periodic_bc tells whether the cuts correspond to periodic boundary
*     conditions (for example in the case of limited plasmas or
*     a diverted stellarator island in a toroidally averaged geometry). 
*     periodic_bc.eq.0 is the default.
*     periodic_bc.eq.1 indicates that periodic boundary conditions are
*     being used.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agdr')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.
     &             1.le.nout(0).and.1.le.nout(1),
     &  'faulty argument ninp, nout')
*   ..test nx, ny
      call xertst (1.le.nx.and.1.le.ny, 'faulty argument nx, ny')
*   ..test nx1, ny1
      call xertst (1.le.nx1.and.1.le.ny1, 'faulty argument nx1, ny1')
*   ..set output control
      idout = '....'

* ..initialise the geometry
      redef_gmtry = 1 ! Since version 3.0.6, we apply the pbsx correction by default
      call ipgeti ('b2agdr_redef_gmtry', redef_gmtry)
*   ..use parg, obtain isymm and lblgm
      call b2xxid ('b2ag', lblgm(1:60))
      call b2agsi (nx, ny, parg, isymm, lblgm(61:120), istyle)
*!!   Be careful about the use of (nx,ny) in b2agsi: the mesh will also
*!!   be generated with dimensions (nx1,ny1).
*   ..obtain grid and mhd equilibrium data from file
      call b2agfs (ninp(1),nx,ny,nx1,ny1,istyle,
     1 crx1,cry1,fpsi1,ffbz1,bb1,wbbc1,cflags1,
     2 lblgm(61:120),leftix1,leftiy1,rightix1,rightiy1,
     3 topix1,topiy1,bottomix1,bottomiy1,
     4 leftcut1,rightcut1,bottomcut1,topcut1,
     5 periodic_bc,nncut,nncutmax)
      call init_region(nx1,ny1,nncut,nncutmax,
     1 leftcut1,rightcut1,topcut1,bottomcut1,
     2 leftix1,leftiy1,rightix1,rightiy1,
     3 topix1,topiy1,bottomix1,bottomiy1,
     4 region1,nnreg,resignore1,crx1,cry1,periodic_bc,cflags1)
      write(*,*) 'b2agdr0: Average bb1(:,:,0) ',
     . b2gavg1(nx1,ny1,bb1(-1:nx1,-1:ny1,0))
*   ..echo lblgm
      write (nout(0),'(/2x,a)') 'label for geometry'
      write (nout(0),'(/4x,a/4x,a)') lblgm(1:60), lblgm(61:120)
*   ..compute geometry and flux/field functions on auxiliary mesh
*   ..if not provided already by Carre grid
      if(istyle.ge.0) call b2agcf (nx1, ny1, crx1, cry1, fpsi1, ffbz1)
*   ..compute metric coefficients on auxiliary mesh
      call b2agmt (nx1, ny1, isymm, crx1, cry1, fpsi1, bb1, vol1,
     & hx1, hy1, hz1, hc1, ht1, qz1, qc1, qs1, qac1, qas1, ebc1, gs1,
     & cell_width1, cell_height1,
     & leftix1, leftiy1, bottomix1, bottomiy1, 
     & rightix1, rightiy1, topix1, topiy1, cflags1)
      call b2agbf (nx1, ny1, isymm, periodic_bc, nncut,
     & crx1, cry1, fpsi1, ffbz1,
     & vol1, gs1, qz1, qc1, qs1, bb1, wbbl1, wbbr1, wbbv1, region1,
     & leftix1, leftiy1, bottomix1, bottomiy1, rightix1, rightiy1,
     & topix1, topiy1, cflags1)
*    ..compute pbs1 (magnetic field pitch * surface area) at cell faces
      geometryType = geometryId(nnreg, periodic_bc, topcut1)
      call b2xgbs (nx1, ny1, wbbl1, wbbv1, gs1, qac1, pbs1,
     &  leftix1, leftiy1, bottomix1, bottomiy1, cflags1)
*    ..interpolate magnetic field to cell corners
      if (maxval(abs(wbbc1)).eq.0.0_R8) then
        !wdk if wbbc not yet provided by grid generator, compute here
        call my_intcor1(nx1,ny1,2,vol1,periodic_bc,nncut,
     &     leftix1,leftiy1,rightix1,rightiy1,
     &     bottomix1,bottomiy1,topix1,topiy1,region1,bb1,wbbc1)
        wbbc1(:,:,3)=
     &     sqrt(wbbc1(:,:,0)**2+wbbc1(:,:,1)**2+wbbc1(:,:,2)**2)
      endif
*   ..project metric coefficients to proper mesh
      call b2agmx (nx1, ny1, nx, ny, 
     & vol1, hx1, hy1, hz1, hc1, ht1, qz1, qc1, qs1, qac1, qas1, ebc1,
     & gs1, crx1, cry1, bb1, pbs1, fpsi1, ffbz1,
     & wbbl1, wbbr1, wbbv1, wbbc1,
     & cell_width1, cell_height1,
     & vol, hx, hy, hz, hc, ht, qz, qc, qs, qac, qas, ebc,
     & gs, crx, cry, bb, pbs, fpsi, ffbz,
     & wbbl, wbbr, wbbv, wbbc,
     & cell_width, cell_height,
     & leftix1, leftiy1, rightix1, rightiy1,
     & topix1, topiy1, bottomix1, bottomiy1,
     & leftcut1, rightcut1, bottomcut1, topcut1,
     & cflags1)
*   ..correct pbsx on the proper mesh
      if(redef_gmtry.ne.0) then
        call b2agbs (nx, ny, nncutmax,
     1   leftix, rightix, rightiy, region, crx, cry,
     1   ffbz, bb, wbbl, gs, qz, qc, hx, hy, vol,
     2   pbs, cell_width, cell_height)
      end if
*   ..establish an isolated region if needed
      call ipgeti('b2agdr_nxiso1', nxiso1)
      call ipgeti('b2agdr_nxiso2', nxiso2)
      call ipgeti('b2agdr_nyiso1', nyiso1)
      call ipgeti('b2agdr_nyiso2', nyiso2)
      if (nxiso1.ne.-2) then
        if (nxiso2.eq.-2) nxiso2 = nx
        if (nyiso1.eq.-2) nyiso1 = -1
        if (nyiso2.eq.-2) nyiso2 = ny
      endif
      if (nxiso2.ne.-2) then
        if (nxiso1.eq.-2) nxiso1 = -1
        if (nyiso1.eq.-2) nyiso1 = -1
        if (nyiso2.eq.-2) nyiso2 = ny
      endif
      if (nyiso1.ne.-2) then
        if (nxiso1.eq.-2) nxiso1 = -1
        if (nxiso2.eq.-2) nxiso2 = nx
        if (nyiso2.eq.-2) nyiso2 = ny
      endif
      if (nyiso2.ne.-2) then
        if (nxiso1.eq.-2) nxiso1 = -1
        if (nxiso2.eq.-2) nxiso2 = nx
        if (nyiso1.eq.-2) nyiso1 = -1
      endif
      if (nxiso1.ne.-2) then
        do iy = nyiso1, nyiso2
          if (isInDomain(nx,ny,leftix(nxiso1,iy),
     .                         leftiy(nxiso1,iy))) then
            rightix(leftix(nxiso1,iy),leftiy(nxiso1,iy)) = nx+1
            region(leftix(nxiso1,iy),leftiy(nxiso1,iy),1) =
     .       region(nxiso2,iy,1)
          end if
          if (isInDomain(nx,ny,rightix(nxiso2,iy),
     .                         rightiy(nxiso2,iy))) then
            leftix(rightix(nxiso2,iy),rightiy(nxiso2,iy)) = -2
            region(rightix(nxiso2,iy),rightiy(nxiso2,iy),1) =
     .       region(nxiso1,iy,1)
          end if
        end do
        do ix = nxiso1, nxiso2
          if (isInDomain(nx,ny,topix(ix,nyiso2),topiy(ix,nyiso2))) then
            bottomiy(topix(ix,nyiso2),topiy(ix,nyiso2)) = -2
            region(topix(ix,nyiso2),topiy(ix,nyiso2),2) =
     .       region(ix,nyiso1,2)
          endif
          if (isInDomain(nx,ny,bottomix(ix,nyiso1),
     .                         bottomiy(ix,nyiso1))) then
            topiy(bottomix(ix,nyiso1),bottomiy(ix,nyiso1)) = ny+1
            region(bottomix(ix,nyiso1),bottomiy(ix,nyiso1),2) =
     .       region(ix,nyiso2,2)
          end if
        end do
        write(*,'()')
        write(*,*) 'Established an isolated region in the range'
        write(*,*) 'NXISO1 = ', nxiso1,' , NXISO2 = ', nxiso2
        write(*,*) 'NYISO1 = ', nyiso1,' , NYISO2 = ', nyiso2
        write(*,'()')
        region(nxiso1:nxiso2,nyiso1:nyiso2,0:2) = 0
        resignore(nxiso1:nxiso2,nyiso1:nyiso2,1:2) = 0
        do iy = nyiso1,nyiso2
          if (isInDomain(nx,ny,leftix(nxiso1,iy),leftiy(nxiso1,iy)))
     >     leftix(nxiso1,iy) = nxiso2
          if (isInDomain(nx,ny,rightix(nxiso2,iy),rightiy(nxiso2,iy)))
     >     rightix(nxiso2,iy) = nxiso1
        end do
        do ix = nxiso1,nxiso2
          if (isInDomain(nx,ny,bottomix(ix,nyiso1),bottomiy(ix,nyiso1)))
     >     bottomiy(ix,nyiso1) = nyiso2
          if (isInDomain(nx,ny,topix(ix,nyiso2),topiy(ix,nyiso2)))
     >     topiy(ix,nyiso2) = nyiso1
        end do
      endif
cdpc-debug
      if(nx1.eq.nx.and.ny1.eq.ny) then
        write(*,*) 'b2agdr sanity check: vol ',
     1   maxval(abs(vol1-vol))
        write(*,*) 'b2agdr sanity check: hx  ',
     1   maxval(abs(hx1-hx))
        write(*,*) 'b2agdr sanity check: hy  ',
     1   maxval(abs(hy1-hy))
        write(*,*) 'b2agdr sanity check: qz  ',
     1   maxval(abs(qz1-qz))
        write(*,*) 'b2agdr sanity check: qc  ',
     1   maxval(abs(qc1-qc))
        write(*,*) 'b2agdr sanity check: qs ',
     1   maxval(abs(qs1-qs))
        write(*,*) 'b2agdr sanity check: gs  ',
     1   maxval(abs(gs1-gs))
        write(*,*) 'b2agdr sanity check: crx ',
     1   maxval(abs(crx1-crx))
        write(*,*) 'b2agdr sanity check: cry ',
     1   maxval(abs(cry1-cry))
        write(*,*) 'b2agdr sanity check: bb  ',
     1   maxval(abs(bb1-bb))
        write(*,*) 'b2agdr sanity check: pbs ',
     1   maxval(abs(pbs1-pbs))
      endif
      write(*,*) 'b2agdr1: bb(nx,ny,0) ',bb(nx,ny,0)
      if(istyle.ge.0) then
*   ..compute geometry and flux/field functions on proper mesh
        call b2agcf (nx, ny, crx, cry, fpsi, ffbz)
*   ..compute magnetic field on proper mesh
        call b2agbb (nx, ny, vol, gs, crx, cry, qz, 
     &   fpsi, ffbz, bb)
      endif
      write(*,*) 'b2agdr2: bb(nx,ny,0) ',bb(nx,ny,0)
      write(*,'()')
      if (bb(nx,ny,2).gt.0.0_R8) then
        write(*,'(a)')
     .   'This case has POSITIVE toroidal magnetic field.'
        if (isymm.ne.0) write(*,'(a)')
     .   'This means that the Bx(grad B) direction points DOWN.'
        write(*,'(a)')
     .   'Bz is out of the page towards you.'
        write(nout(0),'(/2x,a)')
     .   'This case has POSITIVE toroidal magnetic field.'
        if (isymm.ne.0) write(nout(0),'(2x,a)')
     .   'This means that the Bx(grad B) direction points DOWN.'
        write(nout(0),'(2x,a)')
     .   'Bz is out of the page towards you.'
      else if (bb(nx,ny,2).lt.0.0_R8) then
        write(*,'(a)')
     .   'This case has NEGATIVE toroidal magnetic field.'
        if (isymm.ne.0) write(*,'(a)')
     .   'This means that the Bx(grad B) direction points UP.'
        write(*,'(a)')
     .   'Bz is into the page away from you.'
        write(nout(0),'(/2x,a)')
     .   'This case has NEGATIVE toroidal magnetic field.'
        if (isymm.ne.0) write(nout(0),'(2x,a)')
     .   'This means that the Bx(grad B) direction points UP.'
        write(nout(0),'(2x,a)')
     .   'Bz is into the page away from you.'
      else
        write(*,'(a)')
     .   'This case has ZERO toroidal magnetic field.'
        write(nout(0),'(/2x,a)')
     .   'This case has ZERO toroidal magnetic field.'
      endif
      if (bb(nx,ny,2).ne.0.0_R8) then
        write(*,'(a)')
     .   'If this is not correct, please set :'
        write(*,'(a)')
     .   '"b2agfs_Bt_reversal" to "1" in "b2ag.dat" and rerun.'
        write(*,'(a)')
     .   'If already set to "1", set to "0" instead.'
      endif
      write(*,'()')
      if (bb(nx,ny,0).gt.0.0_R8) then
        write(*,'(a)')
     .   'This case has POSITIVE poloidal magnetic field.'
        write(*,'(a)')
     .   'This means the poloidal field points in the '//
     .   'direction of INCREASING poloidal index ix.'
        write(*,'(a)')
     .   'The parallel velocities ua and ue will be '//
     .   'positive when pointing in the direction of '//
     .   'INCREASING poloidal index ix.'
        write(*,'(a)')
     .   'This is the SAME convention as for poloidal '//
     .   'quantities.'
        write(nout(0),'(2x,a)')
     .   'This case has POSITIVE poloidal magnetic field.'
        write(nout(0),'(2x,a)')
     .   'This means the poloidal field points in the '//
     .   'direction of INCREASING poloidal index ix.'
        write(nout(0),'(2x,a)')
     .   'The parallel velocities ua and ue will be '//
     .   'positive when pointing in the direction of '//
     .   'INCREASING poloidal index ix.'
        write(nout(0),'(2x,a)')
     .   'This is the SAME convention as for poloidal '//
     .   'quantities.'
      else
        write(*,'(a)')
     .   'This case has NEGATIVE poloidal magnetic field.'
        write(*,'(a)')
     .   'This means the poloidal field points in the '//
     .   'direction of DECREASING poloidal index ix.'
        write(*,'(a)')
     .   'The parallel velocities ua and ue will be '//
     .   'positive when pointing in the direction of '//
     .   'DECREASING poloidal index ix.'
        write(*,'(a)')
     .   'This is the OPPOSITE convention as for '//
     .   'poloidal quantities.'
        write(nout(0),'(2x,a)')
     .   'This case has NEGATIVE poloidal magnetic field.'
        write(nout(0),'(2x,a)')
     .   'This means the poloidal field points in the '//
     .   'direction of DECREASING poloidal index ix.'
        write(nout(0),'(2x,a)')
     .   'The parallel velocities ua and ue will be '//
     .   'positive when pointing in the direction of '//
     .   'DECREASING poloidal index ix.'
        write(nout(0),'(2x,a)')
     .   'This is the OPPOSITE convention as for '//
     .   'poloidal quantities.'
      end if
      if (bb(nx,ny,0)*bb(nx,ny,2).lt.0.0_R8) then
        write(*,'(a)')
     .   'This case has CO-directed plasma current and '//
     .   'toroidal field.'
        write(nout(0),'(2x,a)')
     .   'This case has CO-directed plasma current and '//
     .   'toroidal field.'
      else if (bb(nx,ny,0)*bb(nx,ny,2).gt.0.0_R8) then
        write(*,'(a)')
     .   'This case has COUNTER-directed plasma current and '//
     .   'toroidal field.'
        write(nout(0),'(2x,a)')
     .   'This case has COUNTER-directed plasma current and '//
     .   'toroidal field.'
      endif
      write(*,'(a)')
     . 'If this is not correct, please set :'
      write(*,'(a)')
     . '"b2agfs_pit_rescale" to "-1.0" in "b2ag.dat" and rerun.'
      write(*,'(a)')
     . 'If already set to "-1.0", set to "1.0" instead.'
      write(*,'()')

* ..output
*   ..write summary of geometry
!      call b2wfgm (nout(0), idout, nx, ny)
*   ..write version number
      call cfverw (nout(1), newversion)
*   ..save dimensions
      idum(0) = nx
      idum(1) = ny
      call cfwuin (nout(1), 2, idum, 'nx,ny')
*   ..save label
      call cfwuch (nout(1), 120, lblgm, 'label')
*   ..save isymm
      idum(0) = isymm
      call cfwuin (nout(1), 1, idum, 'isymm')
*   ..save geometry
      call b2wugm (nout(1), nx, ny, crx, cry, fpsi, ffbz, bb,
     &  vol, hx, hy, hc, ht, qz, qc, qs, qac, qas, ebc, gs, pbs, wbbl,
     &  wbbc, cell_width, cell_height)
*   ..create B2GridMap data structure
      call b2CreateMap(nx, ny, crx, cry, cflags,
     &     leftix, leftiy, rightix, rightiy,
     &     topix, topiy, bottomix, bottomiy, .true., gmap)
      mapInitialized = .true.
*   ..save parameters
      call cfwure (nout(1), 100, parg, 'parg')
      write(*,*) 'b2agdr3: bb(nx,ny,0) ',bb(nx,ny,0)
*   ..in 5.1, and SOLPS-ITER we need a fort.30 file early to build triangular grids
      call write_f30

*   ..conversion to new, unstructured data-types
      call calc_mapping (nx, ny, mpg)
      write (*,*) 'Mapping created'
      call prep_geometry (nx, ny,
     .        mpg%nCv, mpg%nCi, mpg%nFc, mpg%nVx, mpg%nFs, mpg%nFt,
     .        mpg%nVmxCv, gm, mpg)
      write (*,*) 'Geometry created'
*     write unstructured b2fgmtry file
      open (newunit=nout_us,file='b2fgmtry_us',iostat=io)

*   ..write headers for output files
      call cfverw (nout_us, newversion)
      idum(0) = mpg%nCi
      idum(1) = mpg%nCg
      idum(2) = mpg%nCv
      idum(3) = mpg%nFc
      idum(4) = mpg%nVx
      idum(5) = mpg%nFs
      idum(6) = mpg%nFt
      call cfwuin (nout_us,7,idum,'nCi,nCg,nCv,nFc,nVx,nFs,nFt')
      idum(0) = mpg%nCmxVx
      idum(1) = mpg%nCmxFc
      idum(2) = mpg%nVmxCv
      idum(3) = mpg%nVmxFc
      idum(4) = mpg%nCmxNv
      call cfwuin (nout_us,5,idum,'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')
      idum(0) = mpg%isClassicalGrid
      call cfwuin (nout_us,1,idum,'isClassicalGrid')
      if (mpg%isClassicalGrid == 1) then
        idum(0) = mpg%nx
        idum(1) = mpg%ny
        idum(2) = mpg%nncut
      else
        ! dummy values for unstructured grid
        idum(0) = -1
        idum(1) = -1
        idum(2) = -1
      endif
      call cfwuin (nout_us, 3, idum, 'nx,ny,nncut')
      call cfwuch (nout_us, 120, lblgm, 'label')
      idum(0) = isymm
      call cfwuin (nout_us, 1, idum, 'isymm')
*   ..write geometry
      call write_b2fgmtry (nout_us, mpg, gm)
*   ..clean-up
      close (nout_us)
      call dealloc_geometry(gm)
      call dealloc_mapping(mpg)



* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read heading *param')
   92 call xerrab ('error trying to read parg')
   93 call xerrab
     1 ('error trying to read non-default internal parameters')

*-----------------------------------------------------------------------
*.end b2agdr

      end subroutine b2agdr
