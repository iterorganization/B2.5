*-----------------------------------------------------------------------
*.specification

      subroutine b2rugm (nget, nx, ny, crx, cry, fpsi, ffbz,
     &  bb, vol, hx, hy, qz, qc, gs, pbs, wbbl, wbbr, wbbv)
      use b2mod_types
      use b2mod_indirect
      use b2mod_version
      use b2mod_b2cmgs
      implicit none
*   ..input arguments (unchanged on exit)
      integer nget, nx, ny
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     &  fpsi(-1:nx,-1:ny,0:3), ffbz(-1:nx,-1:ny,0:3),
     &  bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny), hx(-1:nx,-1:ny),
     &  hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny), 
     &  gs(-1:nx,-1:ny,0:2), pbs(-1:nx,-1:ny,0:1),
     &  wbbl(-1:nx,-1:ny,0:3), wbbr(-1:nx,-1:ny,0:3), 
     &  wbbv(-1:nx,-1:ny,0:3)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2RUGM reads from an un*formatted file the geometry and magnetic
*     field.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine reads from unit (nget) the specification of a b2
*     geometry and magnetic field. Unit nget must be connected to an
*     un*formatted file.
*     The data consists of the coordinates of the four corners of all
*     the cells in the mesh, the magnetic flux and field functions at
*     those positions, the components of the magnetic field at cell
*     centers, a collection of metric coefficients, and a specification
*     of singular regions in the geometry. See the code to discern the
*     required format of the input file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nget - integer, input.
*     nget specifies a unit number for un*formatted input.
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
*     bb - (-1:nx,-1:ny,0:3) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), bb(ix,iy,0:3) specifies the
*     magnetic field at the center of the (ix,iy) cell.
*     bb(ix,iy,0:2) are the (x,y,z)-components and bb(ix,iy,3) is
*     the absolute magnetic field strength.
*     (bb(,,3)=sqrt(bb(,,0)**2+bb(,,1)**2+bb(,,2)**2.)
*     It will hold that 0.lt.bb(,,3).
*
*     vol - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), vol(ix,iy) specifies the volume of
*     the (ix,iy) cell. (In the case of toroidal symmetry this volume
*     includes a factor 2*pi*crx, while in the case of cylindrical
*     symmetry it includes a factor 2*pi*cry.)
*     It will hold that 0.lt.vol(,).
*
*     hx, hy - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), hx(ix,iy) and hy(ix,iy) specify the
*     x-diameter and the y-diameter, respectively, of the (ix,iy) cell.
*     It will hold that 0.lt.hx(,) and 0.lt.hy(,).
*
*     qz - (-1:nx,-1:ny,0:1) real*8 array, output.
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
*     gs - (-1:nx,-1:ny,0:2) real*8 array, output.
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
*     pbs - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let pbsx(,)=pbs(,,0) and pbsy(,)=pbs(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), pbsx(ix,iy) specifies the product
*     (bx/bb)*sx and pbsy(ix,iy) specifies (by/bb)*sy, or (magnetic 
*     field pitch)*(surface area), namely the parallel contact area, 
*     on the (ix,iy) cell left and bottom faces respectively.  
*     See routine b2xgbs for the form of interpolation employed.
*
*     (parameters in common:)
*
*     /b2cmgs/ - output.
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
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, ix, iy, idum(0:9)
      real (kind=R8) ::
     *   p0x, p0y, p1x, p1y
*   ..procedures
      intrinsic min, max, sqrt
      real (kind=R8) ::
     *   trim1, norm, cosang
      norm(p0x,p0y) = sqrt(p0x**2+p0y**2)
      trim1(p0x) = max(-1.0_R8,min(1.0_R8,p0x))
      cosang(p0x,p0y,p1x,p1y) = trim1(
     &  (p0x*p1y-p0y*p1x)/(norm(p0x,p0y)*norm(p1x,p1y)))
      external subini, subend, xertst, cfrure, cfruin
      external intfaceh1, intfacer1, intfacev1
      external b2xgbs, b2xvsg, b2xvgm

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2rugm')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nget, nx, ny
      call xertst (1.le.nget, 'faulty argument nget')
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

* ..read geometry and field
*   ..read coordinates
      call cfrure (nget, n2*4, crx, 'crx')
      call cfrure (nget, n2*4, cry, 'cry')
*   ..read flux functions
      call cfrure (nget, n2*4, fpsi, 'fpsi')
      call cfrure (nget, n2*4, ffbz, 'ffbz')
*   ..read magnetic field
      call cfrure (nget, n2*4, bb, 'bb')
*   ..read vol, hx, hy
      call cfrure (nget, n2, vol, 'vol')
      call cfrure (nget, n2, hx, 'hx')
      call cfrure (nget, n2, hy, 'hy')
*   ..read qz, qc, gs
      call cfrure (nget, n2*2, qz, 'qz')
c-xpb
      if(b2fgmtry_version.ge.'01.001.022') then
	call cfrure (nget, n2, qc, 'qc')
      else
       do iy = -1, ny
        do ix = -1, nx
         qc(ix,iy) = cosang (
     &    (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &    (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &    crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))
         if(qc(ix,iy).le.0.0_R8) then
          write(*,*) 'Problem in b2rugm: qc(ix,iy) <= 0',
     1     qc(ix,iy),ix,iy
          write(*,'(a,4f12.6)') 'crx = ', crx(ix,iy,0:3)
          write(*,'(a,4f12.6)') 'cry = ', cry(ix,iy,0:3)
         endif
        enddo
       enddo
      endif
c-xpb-end
      call cfrure (nget, n2*3, gs, 'gs')
*   ..read /b2cmgs/
      call cfruin (nget, 1, idum, 'nlreg')
      nlreg = idum(0)
      call xertst (0.le.nlreg.and.nlreg.le.nlrmx, 'faulty input nlreg')
      call cfruin (nget, nlreg, nlxlo, 'nlxlo')
      call cfruin (nget, nlreg, nlxhi, 'nlxhi')
      call cfruin (nget, nlreg, nlylo, 'nlylo')
      call cfruin (nget, nlreg, nlyhi, 'nlyhi')
      call cfruin (nget, nlreg, nlloc, 'nlloc')
c-dpc
      call cfruin (nget, 1, idum, 'nncut')
      nncut=idum(0)
      call xertst (0.le.nncut.and.nncut.le.nncutmax,
     . 'faulty input nncut')
      call cfruin (nget, nncut, leftcut, 'leftcut')
      call cfruin (nget, nncut, rightcut, 'rightcut')
      call cfruin (nget, nncut, topcut, 'topcut')
      call cfruin (nget, nncut, bottomcut, 'bottomcut')
      call cfruin (nget, n2, leftix, 'leftix')
      call cfruin (nget, n2, rightix, 'rightix')
      call cfruin (nget, n2, topix, 'topix')
      call cfruin (nget, n2, bottomix, 'bottomix')
      call cfruin (nget, n2, leftiy, 'leftiy')
      call cfruin (nget, n2, rightiy, 'rightiy')
      call cfruin (nget, n2, topiy, 'topiy')
      call cfruin (nget, n2, bottomiy, 'bottomiy')
      call cfruin (nget, n2*3, region, 'region')
      call cfruin (nget, 3, nnreg, 'nnreg')
      call cfruin (nget, n2*2, resignore, 'resignore')
c-dpc-end 
c-xpb
      if(b2fgmtry_version.ge.'01.001.000') then
	call cfruin (nget, 1, idum, 'periodic_bc')
	periodic_bc=idum(0)
      else
        periodic_bc=0
      endif
c-xpb-end
*  ..interpolate magnetic field to cell faces
c!!!This interpolation scheme does not guarantee that bb(3)=sqrt(bb(0:2)**2)!!!
      call intfaceh1(nx,ny,2,leftix,leftiy,vol,bb,wbbl)
      call intfacer1(nx,ny,2,rightix,rightiy,vol,bb,wbbr)
      call intfacev1(nx,ny,2,bottomix,bottomiy,vol,bb,wbbv)
*   ..extrapolate magnetic field to edges
      do iy = -1, ny
        do ix = -1, nx
          if (leftix(ix,iy).eq.-2) wbbl(ix,iy,0:2) = 
     &     bb(ix,iy,0:2) - (wbbr(ix,iy,0:2)-bb(ix,iy,0:2))
          if (rightix(ix,iy).eq.nx+1) wbbr(ix,iy,0:2) = 
     &     bb(ix,iy,0:2) - (wbbl(ix,iy,0:2)-bb(ix,iy,0:2))
          if (bottomiy(ix,iy).eq.-2) wbbv(ix,iy,0:2) = bb(ix,iy,0:2) -
     &     (wbbv(topix(ix,iy),topiy(ix,iy),0:2)-bb(ix,iy,0:2))
        enddo
      enddo
*   .. correct interpolation of magnetic field for diverted island case
*      as B is not a periodic quantity
      if(periodic_bc.eq.1 .and. nncut.ge.1) then
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).ne.ix-1) then
              if(region(ix,iy,0).eq.1 .or. region(ix,iy,0).eq.4) then
                wbbl(ix,iy,0)=bb(ix,iy,0)-(wbbr(ix,iy,0)-bb(ix,iy,0))
                wbbl(ix,iy,1)=bb(ix,iy,1)-(wbbr(ix,iy,1)-bb(ix,iy,1))
                wbbl(ix,iy,2)=bb(ix,iy,2)-(wbbr(ix,iy,2)-bb(ix,iy,2))
              endif
            else if(rightix(ix,iy).ne.ix+1) then
              if(region(ix,iy,0).eq.1 .or. region(ix,iy,0).eq.3) then
                wbbr(ix,iy,0)=bb(ix,iy,0)-(wbbl(ix,iy,0)-bb(ix,iy,0))
                wbbr(ix,iy,1)=bb(ix,iy,1)-(wbbl(ix,iy,1)-bb(ix,iy,1))
                wbbr(ix,iy,2)=bb(ix,iy,2)-(wbbl(ix,iy,2)-bb(ix,iy,2))
              endif
            endif
          enddo
        enddo
      endif
*   ..compute extrapolated wbb.(:,:,3) components
      wbbl(:,:,3)=sqrt(wbbl(:,:,0)**2+wbbl(:,:,1)**2+wbbl(:,:,2)**2)
      wbbr(:,:,3)=sqrt(wbbr(:,:,0)**2+wbbr(:,:,1)**2+wbbr(:,:,2)**2)
      wbbv(:,:,3)=sqrt(wbbv(:,:,0)**2+wbbv(:,:,1)**2+wbbv(:,:,2)**2)
*   ..compute or read pbs (magnetic field pitch * surface area) at cell faces
      if(b2fgmtry_version.ge.'01.001.063') then                      !srv 15.02.11
        call cfrure (nget, n2*2, pbs, 'pbs')
      else
        call b2xgbs (nx, ny, wbbl, wbbv, gs, qc, crx, cry, pbs)      
      endif
*   ..test sign of bb
      call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
      call b2xvsg (n2, wbbl(-1,-1,3), 1, 'wbbl3', '.gt.')
      call b2xvsg (n2, wbbr(-1,-1,3), 1, 'wbbr3', '.gt.')
      call b2xvsg (n2, wbbv(-1,-1,3), 1, 'wbbv3', '.gt.')
*   ..test geometry
      call b2xvgm (nx, ny, vol, hx, hy, qz, qc, gs, pbs)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2rugm

      end subroutine b2rugm
