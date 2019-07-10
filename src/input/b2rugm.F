*-----------------------------------------------------------------------
*.specification

      subroutine b2rugm (nget, nx, ny, crx, cry, fpsi, ffbz,
     &  bb, vol, hx, hy, hz, qz, qc, qcb, gs, pbs,
     &  hyl, hxb, hzl, hzb, wbbl, wbbr, wbbv, wbbc,
     &  cell_width, cell_height, gmap)
      use b2mod_types
      use b2mod_indirect
      use b2mod_version
      use b2mod_b2cmgs
      use b2mod_b2cmfs
      use b2mod_constants
      use b2mod_grid_mapping
      implicit none
*   ..input arguments (unchanged on exit)
      integer nget, nx, ny
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     &  fpsi(-1:nx,-1:ny,0:3), ffbz(-1:nx,-1:ny,0:3),
     &  bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny), hx(-1:nx,-1:ny),
     &  hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny),
     &  hz(-1:nx,-1:ny), hyl(-1:nx,-1:ny), hxb(-1:nx,-1:ny),
     &  hzl(-1:nx,-1:ny), hzb(-1:nx,-1:ny),
     &  qcb(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2), pbs(-1:nx,-1:ny,0:1),
     &  wbbl(-1:nx,-1:ny,0:3), wbbr(-1:nx,-1:ny,0:3), 
     &  wbbv(-1:nx,-1:ny,0:3), wbbc(-1:nx+1,-1:ny+1,0:3),
     &  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny)
      type(B2GridMap) :: gmap
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
*     hxb, hyl - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), hxb(ix,iy) and hyl(ix,iy) specify the
*     x-diameter and the y-diameter, respectively, of the bottom and
*     left faces of the (ix,iy) cell.
*     It will hold that 0.lt.hxb(,) and 0.lt.hyl(,).
*
*     hz - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), hz(ix,iy) specifies the
*     length along the toroidal direction of the (ix,iy) cell.
*     It will hold that 0.lt.hz(,).
*
*     hzb, hzl - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), hzb(ix,iy) and hzl(ix,iy) specify the
*     length along the toroidal direction, respectively, of the bottom and
*     left faces of the (ix,iy) cell.
*     It will hold that 0.lt.hzb(,) and 0.lt.hzl(,).
*
*     cell_width, cell_height - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), cell_width(ix,iy) and cell_height(ix,iy)
*     specify the average width and height, respectively, of the (ix,iy) cell.
*     It will hold that 0.lt.cell_width(,) and 0.lt.cell_height(,).
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
*     qcb - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qcb(ix,iy) specifies cos(t), where
*     t is the complementary angle between the y-direction of the
*     (ix,iy) cell and its bottom face.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry.
*     It will hold that 0.lt.qcb(,).le.1.
*
*     gs - (-1:nx,-1:ny,0:2) real*8 array, output.
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
      integer iFace, offset
      real (kind=R8) ::
     *   ebc(-1:nx,-1:ny,0:1), ebf(-1:nx,-1:ny,0:1),
     *   hzc(-1:nx+1,-1:ny+1)
      real (kind=R8) ::
     *   p0x, p0y, p1x, p1y, t0, t1
      logical qc_to_do, pbs_to_do
      real (kind=R8), save :: fac_hz = 1.0e0_R8
      real (kind=R8) :: hzconst
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
*   ..read hz flags
      call ipgetr ('b2mndr_hz', fac_hz)
      call xertst (0.0_R8.le.fac_hz .and. 1.0_R8.ge.fac_hz,
     &  'faulty parameter b2mndr_hz')

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
cxpb
      qc_to_do = .false.
      if(b2fgmtry_version.ge.'01.001.022') then
       call cfrure (nget, n2, qc, 'qc')
      else
       qc_to_do = .true.
      endif
      qcb = 1.0_R8
      if(b2fgmtry_version.ge.'01.001.030') then
       call cfrure_opt (nget, n2, qcb, 'qcb')
      end if
      if (maxval(qcb).eq.1.0_R8 .and. minval(qcb).eq.1.0_R8) then
       do iy = -1, ny
        do ix = -1, nx
         qcb(ix,iy) = cosang (
     &     crx(ix,iy,1)-crx(ix,iy,0),  cry(ix,iy,1)-cry(ix,iy,0),
     &    (crx(ix,iy,2)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,1)),
     &    (cry(ix,iy,2)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,1)))
         if(qcb(ix,iy).le.0.0_R8) then
          write(*,'(a,1f12.6,2i5)')
     1     'Problem in b2rugm: qcb(ix,iy) <= 0',qcb(ix,iy),ix,iy
          write(*,'(a,4f12.6)') 'crx = ', crx(ix,iy,0:3)
          write(*,'(a,4f12.6)') 'cry = ', cry(ix,iy,0:3)
         endif
        enddo
       enddo
      endif
cxpb-end
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
cdpc
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
cdpc-end
cxpb
      if(b2fgmtry_version.ge.'01.001.000') then
        call cfruin (nget, 1, idum, 'periodic_bc')
        periodic_bc=idum(0)
      else
        periodic_bc=0
      endif
cxpb-end
*   ..read pbs (magnetic field pitch * surface area) at cell faces
      pbs_to_do = .false.
      if(b2fgmtry_version.ge.'01.001.063') then                      !srv 15.02.11
        call cfrure (nget, n2*2, pbs, 'pbs')
      else
        pbs_to_do = .true.
      endif
      if(b2fgmtry_version.eq.'03.000.006') then
        call cfruin (nget, 1, idum, 'redef_pbs')
        redef_gmtry=idum(0)
      elseif(b2fgmtry_version.ge.'03.000.007') then
        call cfruin (nget, 1, idum, 'redef_gmtry')
        redef_gmtry=idum(0)
        call cfrure (nget, n2*4, wbbl, 'wbbl')
        if (redef_gmtry.eq.0) write(*,*)
     .   'Warning: wbbl value read from b2fgmtry will be overwritten '//
     .   'because redef_gmtry.eq.0 !'
      else
        redef_gmtry=0
      endif

      do iy = -1, ny
        do ix = -1, nx
          if (nnreg(0).eq.1) then
            on_closed_surface(ix,iy) = periodic_bc.eq.1
          else if (nnreg(0).eq.2 .or. nnreg(0).eq.4) then
            on_closed_surface(ix,iy) = region(ix,iy,0).eq.1
          else if (nnreg(0).eq.5 .or. nnreg(0).eq.8) then
            on_closed_surface(ix,iy) = mod(region(ix,iy,0),4).eq.1
          else
            call xerrab ('Unknown number of regions!')
          end if
        end do
      end do

*  ..compute qc if not yet known
      if (qc_to_do) then
cwdk   poloidal unit vector in cell centers, ebc
       do iy = -1, ny
        do ix = -1, nx
         t0 = (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)) 
         t1 = (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2))
         ebc(ix,iy,0) = t0/norm(t0,t1)
         ebc(ix,iy,1) = t1/norm(t0,t1)
        enddo
       enddo
cwdk   poloidal unit vector at left faces, ebf
       call intfaceh(nx,ny,leftix,leftiy,vol,ebc(-1,-1,0),ebf(-1,-1,0))
       call intfaceh(nx,ny,leftix,leftiy,vol,ebc(-1,-1,1),ebf(-1,-1,1))
       do iy = -1, ny
        do ix = -1, nx
cwdk     note: qc in left ghosts unchanged
         if (leftix(ix,iy).ne.-2) then
          qc(ix,iy) = cosang ( ebf(ix,iy,0), ebf(ix,iy,1), 
     &      crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))
         else
          qc(ix,iy) = cosang (
     &     (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &     (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &      crx(ix,iy,2)-crx(ix,iy,0),  cry(ix,iy,2)-cry(ix,iy,0))
         end if
         if(qc(ix,iy).le.0.0_R8) then
          write(*,'(a,1f12.6,2i5)')
     1     'Problem in b2rugm: qc(ix,iy) <= 0',qc(ix,iy),ix,iy
          write(*,'(a,4f12.6)') 'crx = ', crx(ix,iy,0:3)
          write(*,'(a,4f12.6)') 'cry = ', cry(ix,iy,0:3)
         endif
        enddo
       enddo
      endif

*   ..compute hz                                                         !srv 09.01.01  09.06.08 {
!      hzconst = 2.0_R8*pi*6.0_R8   ! for ITER for testing
!      hzconst = 2.0_R8*pi*1.6_R8   ! for ASDEX-UPGRADE for testing
      hzconst = 1.0_R8                                                   !srv 30.01.10
      hz(:,:) = hzconst*(1.0e0_R8-fac_hz)+fac_hz*vol(:,:)/gs(:,:,2)      !srv 16.06.09
      hyl(:,:)= sqrt((crx(:,:,2)-crx(:,:,0))**2+
     &               (cry(:,:,2)-cry(:,:,0))**2)                         ! IYS 04.03.2019
      hxb(:,:)= sqrt((crx(:,:,1)-crx(:,:,0))**2+
     &               (cry(:,:,1)-cry(:,:,0))**2)                         ! IYS 04.03.2019
      if (isymm.eq.0) then
        hzl(:,:) = hz(:,:)
        hzb(:,:) = hz(:,:)
        hzc(-1:nx,-1:ny) = hz(-1:nx,-1:ny)
        hzc(nx+1,-1:ny) = hz(nx,-1:ny)
        hzc(-1:nx,ny+1) = hz(-1:nx,ny)
        hzc(nx+1,ny+1)  = hz(nx,ny)
      else if (isymm.eq.1 .or. isymm.eq.2) then
        hzl(:,:) = pi * (crx(:,:,0)+crx(:,:,2))
        hzb(:,:) = pi * (crx(:,:,0)+crx(:,:,1))
        hzc(-1:nx,-1:ny) = pi * crx(:,:,0)
        hzc(nx+1,-1:ny)  = pi * crx(nx,-1:ny,1)
        hzc(-1:nx,ny+1)  = pi * crx(-1:nx,ny,2)
        hzc(nx+1,ny+1)   = pi * crx(nx,ny,3)
      else if (isymm.eq.3 .or. isymm.eq.4) then
        hzl(:,:) = pi * (cry(:,:,0)+cry(:,:,2))
        hzb(:,:) = pi * (cry(:,:,0)+cry(:,:,1))
        hzc(-1:nx,-1:ny) = pi * cry(:,:,0)
        hzc(nx+1,-1:ny)  = pi * cry(nx,-1:ny,1)
        hzc(-1:nx,ny+1)  = pi * cry(-1:nx,ny,2)
        hzc(nx+1,ny+1)   = pi * cry(nx,ny,3)
      end if
      if (redef_gmtry.eq.0) then
*  ..old treatment: interpolate magnetic field to cell faces : overwrite wbbl
c!!!This interpolation scheme does not guarantee that bb(3)=sqrt(bb(0:2)**2)!!!
        call intfaceh1(nx,ny,2,leftix,leftiy,vol,bb,wbbl)
        call intfacer1(nx,ny,2,rightix,rightiy,vol,bb,wbbr)
        call intfacev1(nx,ny,2,bottomix,bottomiy,vol,bb,wbbv)
        call my_intcor1(nx,ny,2,vol,periodic_bc,nncut,
     1     leftix,leftiy,rightix,rightiy,
     1     bottomix,bottomiy,topix,topiy,region,bb,wbbc)
*   ..extrapolate magnetic field to edges
        do iy = -1, ny
          do ix = -1, nx
            if (leftix(ix,iy).eq.-2) wbbl(ix,iy,0:2) =
     &       bb(ix,iy,0:2) - (wbbr(ix,iy,0:2)-bb(ix,iy,0:2))
            if (rightix(ix,iy).eq.nx+1) wbbr(ix,iy,0:2) =
     &       bb(ix,iy,0:2) - (wbbl(ix,iy,0:2)-bb(ix,iy,0:2))
            if (bottomiy(ix,iy).eq.-2) wbbv(ix,iy,0:2) = bb(ix,iy,0:2) -
     &       (wbbv(topix(ix,iy),topiy(ix,iy),0:2)-bb(ix,iy,0:2))
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
                  wbbc(ix,iy,0)=
     &             wbbv(ix,iy,0)-(wbbc(ix+1,iy,0)-wbbv(ix,iy,0))
                  wbbc(ix,iy,1)=
     &             wbbv(ix,iy,1)-(wbbc(ix+1,iy,1)-wbbv(ix,iy,1))
                  wbbc(ix,iy,2)=
     &             wbbv(ix,iy,2)-(wbbc(ix+1,iy,2)-wbbv(ix,iy,2))
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
          do ix = -1, nx
            if(leftix(ix,ny).ne.ix-1) then
              if(region(ix,ny,0).eq.1 .or. region(ix,ny,0).eq.4) then
                wbbc(ix,ny,0)=
     &           wbbv(ix,ny,0)-(wbbc(ix+1,ny,0)-wbbv(ix,ny,0))
                wbbc(ix,ny,1)=
     &           wbbv(ix,ny,1)-(wbbc(ix+1,ny,1)-wbbv(ix,ny,1))
                wbbc(ix,ny,2)=
     &           wbbv(ix,ny,2)-(wbbc(ix+1,ny,2)-wbbv(ix,ny,2))
              endif
            endif
          enddo
        endif
*   ..compute extrapolated wbb.(:,:,3) components
        wbbl(:,:,3)=sqrt(wbbl(:,:,0)**2+wbbl(:,:,1)**2+wbbl(:,:,2)**2)
        wbbr(:,:,3)=sqrt(wbbr(:,:,0)**2+wbbr(:,:,1)**2+wbbr(:,:,2)**2)
        wbbv(:,:,3)=sqrt(wbbv(:,:,0)**2+wbbv(:,:,1)**2+wbbv(:,:,2)**2)
        wbbc(:,:,3)=sqrt(wbbc(:,:,0)**2+wbbc(:,:,1)**2+wbbc(:,:,2)**2)
      else ! redef_gmtry.eq.1 treatment
        wbbr = 0.0_R8
        wbbv = 0.0_R8
        wbbc = 0.0_R8
        do iy = -1, ny
          do ix = -1, nx
            if (rightix(ix,iy).lt.nx+1) then
              wbbr(ix,iy,0) = wbbl(rightix(ix,iy),rightiy(ix,iy),0)
              wbbr(ix,iy,2) = wbbl(rightix(ix,iy),rightiy(ix,iy),2)
            else
              wbbr(ix,iy,0) =
     &         bb(ix,iy,0) - (wbbl(ix,iy,0)-bb(ix,iy,0))
              wbbr(ix,iy,2) =
     &         bb(ix,iy,2) - (wbbl(ix,iy,2)-bb(ix,iy,2))
            end if
            if (bottomiy(ix,iy).gt.-2) wbbv(ix,iy,0) = (
     &          bb(ix,iy,0) * hy(ix,iy) * hz(ix,iy) * qz(ix,iy,1) +
     &          bb(bottomix(ix,iy),bottomiy(ix,iy),0) *
     &          hy(bottomix(ix,iy),bottomiy(ix,iy)) *
     &          hz(bottomix(ix,iy),bottomiy(ix,iy)) *
     &          qz(bottomix(ix,iy),bottomiy(ix,iy),1) ) /
     &        ( hy(ix,iy) * hz(ix,iy) * qz(ix,iy,1) +
     &          hy(bottomix(ix,iy),bottomiy(ix,iy)) *
     &          hz(bottomix(ix,iy),bottomiy(ix,iy)) *
     &          qz(bottomix(ix,iy),bottomiy(ix,iy),1) )
            wbbv(ix,iy,2) =(ffbz(ix,iy,0)+ffbz(ix,iy,1)) /
     &       (2.0_R8*hzb(ix,iy))
            if (bottomiy(ix,iy).gt.-2 .and. topiy(ix,iy).lt.ny+1) then
              if (leftix(ix,iy).gt.-2) then
               wbbc(ix,iy,0) =
     &         (wbbl(ix,iy,0) * gs(ix,iy,0) * qc(ix,iy) +
     &          wbbl(bottomix(ix,iy),bottomiy(ix,iy),0) *
     &            gs(bottomix(ix,iy),bottomiy(ix,iy),0) *
     &            qc(bottomix(ix,iy),bottomiy(ix,iy)) ) /
     &           (gs(ix,iy,0) * qc(ix,iy) +
     &            gs(bottomix(ix,iy),bottomiy(ix,iy),0) *
     &            qc(bottomix(ix,iy),bottomiy(ix,iy)) )
              else
               wbbc(ix,iy,0) =
     &         (wbbl(ix,iy,0) *
     &            gs(rightix(ix,iy),rightiy(ix,iy),0) *
     &            qc(rightix(ix,iy),rightiy(ix,iy)) +
     &          wbbl(rightix(bottomix(ix,iy),bottomiy(ix,iy)),
     &               rightiy(bottomix(ix,iy),bottomiy(ix,iy)),0) *
     &            gs(rightix(bottomix(ix,iy),bottomiy(ix,iy)),
     &               rightiy(bottomix(ix,iy),bottomiy(ix,iy)),0) *
     &            qc(rightix(bottomix(ix,iy),bottomiy(ix,iy)),
     &               rightiy(bottomix(ix,iy),bottomiy(ix,iy))) ) /
     &           (gs(rightix(ix,iy),rightiy(ix,iy),0) *
     &            qc(rightix(ix,iy),rightiy(ix,iy)) +
     &            gs(rightix(bottomix(ix,iy),bottomiy(ix,iy)),
     &               rightiy(bottomix(ix,iy),bottomiy(ix,iy)),0) *
     &            qc(rightix(bottomix(ix,iy),bottomiy(ix,iy)),
     &               rightiy(bottomix(ix,iy),bottomiy(ix,iy))) )
              endif
            endif
            wbbc(ix,iy,2) = ffbz(ix,iy,0) / hzc(ix,iy)
          end do
          if (bottomiy(nx,iy).gt.-2) wbbc(nx+1,iy,0) =
     &      (wbbr(nx,iy,0) * hy(nx,iy) * hz(nx,iy) * qz(nx,iy,1) +
     &       wbbr(bottomix(nx,iy),bottomiy(nx,iy),0) *
     &         hy(bottomix(nx,iy),bottomiy(nx,iy)) *
     &         hz(bottomix(nx,iy),bottomiy(nx,iy)) *
     &         qz(bottomix(nx,iy),bottomiy(nx,iy),1) ) /
     &      (hy(nx,iy) * hz(nx,iy) * qz(nx,iy,1) +
     &       hy(bottomix(nx,iy),bottomiy(nx,iy)) *
     &       hz(bottomix(nx,iy),bottomiy(nx,iy)) *
     &       qz(bottomix(nx,iy),bottomiy(nx,iy),1) )
          wbbc(nx+1,iy,2) = ffbz(nx,iy,1) / hzc(nx+1,iy)
        end do
        do iy = -1, ny
          do ix = -1, nx
            if (bottomiy(ix,iy).eq.-2) then
              wbbv(ix,iy,0) = bb(ix,iy,0) -
     &         (wbbv(topix(ix,iy),topiy(ix,iy),0)-bb(ix,iy,0))
              wbbc(ix,iy,0) = wbbl(ix,iy,0) -
     &         (wbbc(topix(ix,iy),topiy(ix,iy),0) - wbbl(ix,iy,0))
            end if
            if (topiy(ix,iy).eq.ny+1) then
              wbbc(ix,iy,0) = wbbl(bottomix(ix,iy),bottomiy(ix,iy),0) -
     &         (wbbc(bottomix(ix,iy),bottomiy(ix,iy),0) -
     &          wbbl(bottomix(bottomix(ix,iy),bottomiy(ix,iy)),
     &               bottomiy(bottomix(ix,iy),bottomiy(ix,iy)),0))
            end if
          end do
          if (bottomiy(nx,iy).eq.-2) then
            wbbc(nx+1,iy,0) = wbbl(nx,iy,0) -
     &       (wbbc(topix(nx,iy),topiy(nx,iy),0) - wbbl(nx,iy,0))
          end if
        end do
        do ix = -1, nx
          wbbc(ix,ny+1,0) =
     &     wbbl(ix,ny,0) - (wbbc(ix,ny,0) - wbbl(ix,ny,0))
          wbbc(ix,ny+1,2) = ffbz(ix,ny,2) / hzc(ix,ny+1)
        end do
        wbbc(nx+1,ny+1,0) =
     &   wbbr(nx,ny,0) - (wbbc(nx,ny,0) - wbbl(nx,ny,0))
        wbbc(nx+1,ny+1,2) = ffbz(nx,ny,3) / hzc(nx+1,ny+1)
*   .. correct interpolation of magnetic field for diverted island case
*      as B is not a periodic quantity
        if(periodic_bc.eq.1 .and. nncut.ge.1) then
          do iy = -1, ny
            do ix = -1, nx
              if(leftix(ix,iy).ne.ix-1) then
                if(region(ix,iy,0).eq.1 .or. region(ix,iy,0).eq.4) then
                  wbbc(ix,iy,0)=
     &             wbbv(ix,iy,0)-(wbbc(ix+1,iy,0)-wbbv(ix,iy,0))
                  wbbc(ix,iy,2)=
     &             wbbv(ix,iy,2)-(wbbc(ix+1,iy,2)-wbbv(ix,iy,2))
                endif
              elseif(rightix(ix,iy).ne.ix+1) then
                if(region(ix,iy,0).eq.1 .or. region(ix,iy,0).eq.3) then
                  wbbr(ix,iy,0)=bb(ix,iy,0)-(wbbl(ix,iy,0)-bb(ix,iy,0))
                  wbbr(ix,iy,2)=bb(ix,iy,2)-(wbbl(ix,iy,2)-bb(ix,iy,2))
                endif
              endif
            enddo
          enddo
          do ix = -1, nx
            if(leftix(ix,ny).ne.ix-1) then
              if(region(ix,ny,0).eq.1 .or. region(ix,ny,0).eq.4) then
                wbbc(ix,ny,0)=
     &           wbbv(ix,ny,0)-(wbbc(ix+1,ny,0)-wbbv(ix,ny,0))
                wbbc(ix,ny,2)=
     &           wbbv(ix,ny,2)-(wbbc(ix+1,ny,2)-wbbv(ix,ny,2))
              endif
            endif
          enddo
        endif
*   ..compute extrapolated wbb.(:,:,3) components
        wbbr(:,:,3)=sqrt(wbbr(:,:,0)**2+wbbr(:,:,1)**2+wbbr(:,:,2)**2)
        wbbv(:,:,3)=sqrt(wbbv(:,:,0)**2+wbbv(:,:,1)**2+wbbv(:,:,2)**2)
        wbbc(:,:,3)=sqrt(wbbc(:,:,0)**2+wbbc(:,:,1)**2+wbbc(:,:,2)**2)
      end if
*   ..compute pbs (magnetic field pitch * surface area) at cell faces if needed
      if(pbs_to_do)
     > call b2xgbs (nx, ny, wbbl, wbbv, gs, qc, crx, cry, pbs)

*   ..read extended grid cell flags
      cflags = 0
      if(b2fgmtry_version.ge.'01.001.028') then
        call cfruin_opt (nget, n2*GRID_N_CELLFLAGS, cflags, 'cflags')
      endif
      if(maxval(cflags).eq.0) then
        cflags(1:nx-2,1:ny-2,CELLFLAG_TYPE) = GRID_INTERNAL
        cflags(0,0:ny-1,CELLFLAG_TYPE) = GRID_BOUNDARY
        cflags(nx-1,0:ny-1,CELLFLAG_TYPE) = GRID_BOUNDARY
        cflags(0:nx-1,0,CELLFLAG_TYPE) = GRID_BOUNDARY
        cflags(0:nx-1,ny-1,CELLFLAG_TYPE) = GRID_BOUNDARY
        cflags(-1,-1:ny,CELLFLAG_TYPE) = GRID_GUARD
        cflags(nx,-1:ny,CELLFLAG_TYPE) = GRID_GUARD
        cflags(-1:nx,-1,CELLFLAG_TYPE) = GRID_GUARD
        cflags(-1:nx,ny,CELLFLAG_TYPE) = GRID_GUARD
        if (nncut.ge.2) then
          iy = ny/2
          do ix = 1, nx-2
            if (rightix(ix,iy).eq.nx+1) then
              cflags(ix-1,0:ny-1,CELLFLAG_TYPE) = GRID_BOUNDARY
              cflags(ix,-1:ny,CELLFLAG_TYPE) = GRID_GUARD
            end if
            if (leftix(ix,iy).eq.-2) then
              cflags(ix,-1:ny,CELLFLAG_TYPE) = GRID_GUARD
              cflags(ix+1,0:ny-1,CELLFLAG_TYPE) = GRID_BOUNDARY
            end if
          end do
        end if
        do iy=-1,ny
          do ix=-1,nx
            if (cflags(ix, iy, CELLFLAG_TYPE) /= GRID_BOUNDARY) cycle
            if (cflags(leftix(ix,iy),leftiy(ix,iy),CELLFLAG_TYPE) ==
     .          GRID_GUARD) then
              iFace = CELLFLAG_LEFTFACE
              offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
              cflags(ix, iy, iFace) = offset - region(ix, iy, 0)
              cflags(leftix(ix,iy),leftiy(ix,iy),CELLFLAG_RIGHTFACE) =
     .         cflags(ix, iy, iFace)
            end if
            if (cflags(rightix(ix,iy),rightiy(ix,iy),CELLFLAG_TYPE) ==
     .          GRID_GUARD) then
              iFace = CELLFLAG_RIGHTFACE
              offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
              cflags(ix, iy, iFace) = offset - region(ix, iy, 0)
              cflags(rightix(ix,iy),rightiy(ix,iy),CELLFLAG_LEFTFACE) =
     .         cflags(ix, iy, iFace)
            end if
            if (cflags(topix(ix,iy),topiy(ix,iy),CELLFLAG_TYPE) ==
     .          GRID_GUARD) then
              iFace = CELLFLAG_TOPFACE
              offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
              cflags(ix, iy, iFace) = offset - region(ix, iy, 0)
              cflags(topix(ix,iy),topiy(ix,iy),CELLFLAG_BOTTOMFACE) =
     .         cflags(ix, iy, iFace)
            end if
            if (cflags(bottomix(ix,iy),bottomiy(ix,iy),CELLFLAG_TYPE) ==
     .          GRID_GUARD) then
              iFace = CELLFLAG_BOTTOMFACE
              offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
              cflags(ix, iy, iFace) = offset - region(ix, iy, 0)
              cflags(bottomix(ix,iy),bottomiy(ix,iy),CELLFLAG_TOPFACE) =
     .         cflags(ix, iy, iFace)
            end if
          end do
        end do
      end if
      cell_width = hx
      cell_height(:,:) = hy(:,:)*qz(:,:,1)
      if(b2fgmtry_version.ge.'01.001.033') then
        call cfrure_opt (nget, n2, cell_width, 'cell_width')
        call cfrure_opt (nget, n2, cell_height, 'cell_height')
      endif
*   ..test sign of bb
      call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
      call b2xvsg (n2, wbbl(-1,-1,3), 1, 'wbbl3', '.gt.')
      call b2xvsg (n2, wbbr(-1,-1,3), 1, 'wbbr3', '.gt.')
      call b2xvsg (n2, wbbv(-1,-1,3), 1, 'wbbv3', '.gt.')
      call b2xvsg ((nx+3)*(ny+3), wbbc(-1,-1,3), 1, 'wbbc3', '.gt.')
*   ..test geometry
      call xertst(maxval(region(-1:nx,-1:ny,1)).gt.0 .or.
     &            minval(region(-1:nx,-1:ny,2)).lt.0,
     $            'X-regions not set, re-run b2ag!')
      call xertst(maxval(region(-1:nx,-1:ny,2)).gt.0 .or.
     &            minval(region(-1:nx,-1:ny,1)).lt.0,
     $            'Y-regions not set, re-run b2ag!')
      call b2xvgm (nx, ny, vol, hx, hy, qz, qc, qcb, gs, pbs,
     &     cell_width, cell_height)
*   ..create B2GridMap data structure
      if (nget.ne.70) then ! skip case when reading b2fgmtry.old from b2yt
        call b2CreateMap(nx, ny, crx, cry, cflags,
     &     leftix, leftiy, rightix, rightiy,
     &     topix, topiy, bottomix, bottomiy, .true., gmap)
        mapInitialized = .true.
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2rugm

      end subroutine b2rugm
