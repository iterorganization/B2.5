*-----------------------------------------------------------------------
*.specification

      subroutine b2ygdr (ninp, nout, nx, ny)
      use b2mod_types
      use b2mod_geo
      use b2mod_b2cmgs
      use b2mod_b2cmwg
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:1), nout(0:1), nx, ny
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YGDR is a post-processor for b2. It reads and displays the
*     geometry and magnetic field.
*
*
*  3. description (see also routine b2cdca)
*
*     This subroutine displays a b2 geometry and magnetic flux
*     functions.
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
*     in the arrays crx and cry. The metric parameters (vol, hx, hy, qz, qc
*     and gs; described elsewhere) are derived from crx and cry.
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
*     ninp(0): formatted; provides output control parameters.
*     ninp(1): un*formatted; full specification of the geometry.
*
*     nout - (0:1) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides summary printed output.
*     nout(1): graphical; displays the geometry and the flux functions.
*
*     nx, ny - integer, input.
*     nx and ny specify the number of interior cells along the first
*     and the second coordinate, respectively. The total number of
*     cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
*     It will hold that 0.le.nx and 0.le.ny.
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks
*   ..local variables
      integer iftp, iplot, idum(0:9)
      character lblgm*120, idout*80, id*8
*   ..procedures
      logical streql
      external subini, subend, xertst, xerrab, streql, gacwk,
     &  gdawk, cfruin, cfruch, lbfilx
      external b2rugm, b2wfgm, b2wggm

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
*     lblgm - character*120
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
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ygdr')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.
     &             1.le.nout(0).and.1.le.nout(1),
     &  'faulty argument ninp, nout')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..read output control
      read (ninp(0),'(a8)',err=91) id
      call xertst (streql(id,'*idout'),
     &  'faulty input instead of *idout, id='//id)
      read (ninp(0),*,err=92) iftp, idout
      read (ninp(0),'(a8)',err=93) id
      call xertst (streql(id,'*nmbpcc'),
     &  'faulty input instead of *nmbpcc, id='//id)
      read (ninp(0),*,err=94) nmbpcc
      call xertst (0.le.nmbpcc.and.nmbpcc.le.nmbpmx,
     &  'faulty data for nmbpcc')
      do iplot = 0, nmbpcc-1
       read (ninp(0),*,err=94) cwaxr(iplot), cwayr(iplot)
       call xertst (0.0_R8.le.cwaxr(iplot).and.cwaxr(iplot).le.1.0_R8
     &         .and.0.0_R8.le.cwayr(iplot).and.cwayr(iplot).le.1.0_R8,
     &             'faulty data for cwaxr, cwayr')
      enddo
      read (ninp(0),'(a8)',err=95) id
      call xertst (streql(id,'*nmbpdm'),
     &  'faulty input instead of *nmbpdm, id='//id)
      read (ninp(0),*,err=96) nmbpdm
      call xertst (0.le.nmbpdm.and.nmbpdm.le.nmbpmx,
     &  'faulty data for nmbpdm')
      do iplot = 0, nmbpdm-1
       read (ninp(0),*,err=96) cwxmin(iplot), cwxmax(iplot),
     &                         cwymin(iplot), cwymax(iplot)
       call xertst (cwxmin(iplot).lt.cwxmax(iplot).and.
     &              cwymin(iplot).lt.cwymax(iplot),
     &   'faulty data for cwxmin, cwxmax, cwymin, cwymax')
      enddo
*   ..write output control
      write (nout(0),'(/2x,a,2x,i3/2x,a,2x,a)')
     &  'style of contour plots:', iftp,
     &  'idout:', idout
      write (nout(0),'(2x,a,2x,i3/(2x,0p,2f8.4))')
     &  'selection of cartesian plots', nmbpcc,
     &  (cwaxr(iplot), cwayr(iplot), iplot=0,nmbpcc-1)
      write (nout(0),'(2x,a,2x,i3/(2x,1p,4e12.3))')
     &  'selection of physical domain plots', nmbpdm,
     &  (cwxmin(iplot), cwxmax(iplot),cwymin(iplot), cwymax(iplot),
     &   iplot=0,nmbpdm-1)

* ..initialise the geometry
*   ..read label
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
*   ..echo lblgm
      write (nout(0),'(/2x,a)') 'label for geometry'
      write (nout(0),'(/4x,a/4x,a)') lblgm(1:60), lblgm(61:120)
*   ..read geometry
      call b2rugm (ninp(1), nx, ny, crx, cry, fpsi, ffbz,
     &  bb, vol, hx, hy, qz, qc, qcb, gs, pbs, wbbl, wbbr, wbbv, wbbc,
     &  cell_width, cell_height, gmap)

* ..output
*   ..write summary of geometry
      call b2wfgm (nout(0), idout, nx, ny)
*   ..activate graphics on nout(1)
      call gacwk (nout(1))
*   ..mark contour style
      call lbfilx (iftp)
*   ..plot geometry
      call b2wggm (nout(0), idout, nx, ny, lblgm(61:120), crx,
     &  cry, fpsi, ffbz, bb)
*   ..deactivate graphics on nout(1)
      call gdawk (nout(1))

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read heading *idout')
   92 call xerrab ('error trying to read output control')
   93 call xerrab ('error trying to read heading *nmbpcc')
   94 call xerrab ('error trying to read nmbpcc data')
   95 call xerrab ('error trying to read heading *nmbpdm')
   96 call xerrab ('error trying to read nmbpdm data')

*-----------------------------------------------------------------------
*.end b2ygdr

      end subroutine b2ygdr
