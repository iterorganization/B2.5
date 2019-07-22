*-----------------------------------------------------------------------
*.specification

      subroutine b2agmt (nx, ny, isymm, crx, cry,
     & fpsi, bb, vol, hx, hy, hc, ht, qz, qc, qs, qac, qas, ebc,
     & gs, cell_width, cell_height, leftix, leftiy, bottomix, bottomiy, 
     & rightix, rightiy, topix, topiy, cflag)
      use b2mod_types
      use b2mod_constants
      use tradui_constants
      use b2mod_connectivity
      use b2mod_cellhelper
      use b2mod_indirect
     , ,only : TriangleArea, TriangleCentroid, quadCentroid, quadArea
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isymm
      integer leftix(-1:nx,-1:ny), leftiy(-1:nx,-1:ny), 
     &  bottomix(-1:nx,-1:ny), bottomiy(-1:nx,-1:ny),
     &  rightix(-1:nx,-1:ny), rightiy(-1:nx,-1:ny),
     &  topix(-1:nx,-1:ny), topiy(-1:nx,-1:ny),
     $  cflag(-1:nx,-1:ny,GRID_N_CELLFLAGS) 
      real (kind=R8) ::
     *  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     &  fpsi(-1:nx,-1:ny,0:3,0:2), bb(-1:nx,-1:ny,0:2)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny,0:1), qs(-1:nx,-1:ny,0:1),
     &  gs(-1:nx,-1:ny,0:2), qac(-1:nx,-1:ny,0:1), qas(-1:nx,-1:ny,0:1),
     &  ebc(-1:nx,-1:ny,0:2), 
     &  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny),
     &  hc(-1:nx,-1:ny,0:3), ht(-1:nx,-1:ny,0:1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AGMT computes the metric coefficients given the locations of
*     cell vertices.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine takes as input the coordinates of the four corners
*     of all the cells in the mesh. It computes metric coefficients of
*     the mesh: vol, hx, hy, qz, and gs.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny - integer, input.
*     See the calling routine for description.
*
*     isymm - integer, input.
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
*     crx, cry - (-1:nx,-1:ny,0:3) real*8 array, input.
*     See the calling routine for description.
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
*     qc - (-1:nx,-1:ny,0:1) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qc(ix,iy) specifies cos(t), where 
*     t is the complementary angle between the x-direction of the 
*     (ix,iy) cell and its left face.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. 
*     It will hold that 0.lt.qc(,,).le.1.
*
*     qs - (-1:nx,-1:ny,0:1) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qs(ix,iy) specifies sin(t), where
*     t is the complementary angle between the x-direction of the
*     (ix,iy) cell and its left face.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry.
*     It will hold that -1.lt.qs(,,).le.1.
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
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix, iy, i, j
      integer i1, i2, ii
      integer volhmk, cgeo
      real (kind=R8) ::
     *   p0x, p0y, p1x, p1y, width, centroid(0:1), t0, t1, psic
      logical classical
      real (kind=R8) ::
     *   ecc(-1:nx,-1:ny,0:3,0:1), ecf(0:1), ebf(0:1),
     *   isx(0:1), isy(0:1)
      integer, save :: qc_style = 1
*   ..procedures
      intrinsic min, max, sqrt
      real (kind=R8) ::
     *   trim1, norm, dist, sinang, cosang
      real (kind=R8) ::
     &  sin_t, cos_t, sin_tl, cos_tl, tt
      external subini, subend, xertst, ipgetr, ipgeti
      external b2xvsg
      norm(p0x,p0y) = sqrt(p0x**2+p0y**2)
      trim1(p0x) = max(-1.0_R8,min(1.0_R8,p0x))
      dist(ix,iy,i,j) = norm (
     &  crx(ix,iy,j)-crx(ix,iy,i), cry(ix,iy,j)-cry(ix,iy,i))
      sinang(p0x,p0y,p1x,p1y) = trim1(
     &  (p0x*p1x+p0y*p1y)/(norm(p0x,p0y)*norm(p1x,p1y)))
      cosang(p0x,p0y,p1x,p1y) = trim1(
     &  (p0x*p1y-p0y*p1x)/(norm(p0x,p0y)*norm(p1x,p1y)))
*   ..initialisation
      data volhmk/1/, width /1.0_R8/
      save volhmk, width

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Let (p0x,p0y) specify a vector p0, and (p1x,p1y) specify
*     a vector p1. Then norm(p0x,p0y) is the standard norm of p0,
*     sinang(p0x,p0y,p1x,p1y) is sin(t), and cosang(p0x,p0y,p1x,p1y)
*     is cos(t), where t is the complementary angle between p0 and p1.
*     t=pi/2 when p0 and p1 are parallel, t=0 or t=pi when p0 and p1
*     are orthogonal, and t=-pi/2 when p0 and p1 are anti-parallel.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agmt')
*   ..read input data
      call ipgeti ('b2agmt_volhmk', volhmk)
      call ipgetr ('b2agmt_1d_width', width)
      call ipgeti ('b2agmt_qc_style', qc_style)
      call xertst (0.0_R8.lt.width, 'faulty input width')
*   ..test nx, ny, isymm
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (0.le.isymm.and.isymm.le.4,
     &  'faulty argument isymm')

* ..compute metric coefficients
      classical = isClassicalGrid(cflag)

*   ..compute hx, hy, hc, ht, ecc, ebc
      hx = INVALID_DOUBLE
      hy = INVALID_DOUBLE
      cell_width = INVALID_DOUBLE
      cell_height = INVALID_DOUBLE
      qz(:,:,0) = 0.0_R8
      qz(:,:,1) = 1.0_R8
      do iy = -1, ny
       do ix = -1, nx
        if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle

        ! Requirements for US solver:
        ! hx: not needed
        ! hy: not needed
        ! qz: not needed
        ! hc(:,:,0:3): distance between cell centroid and face centers
        ! ht(:,:,0:1): length of left and bottom faces
        ! ecc(:,:,0:3,0:1): unit vectors of the centroid-face connections
        ! ebc(:,:,0:2): unit vector magnetic field field in cell center

        cgeo = cellGeoType( crx(ix,iy,:), cry(ix,iy,:) )
        select case (cgeo)
        case(CGEO_BROKEN)
          write(*,*) 'Problem in b2agmt: broken cell', ix,iy
          write(*,'(a,4f12.6)') 'crx = ',crx(ix,iy,0:3)
          write(*,'(a,4f12.6)') 'cry = ',cry(ix,iy,0:3)
          call xerrab('Broken cell geometry')
        case(CGEO_QUAD)
          ! For quadrangular cells: assume cell center = geometric mean
          hx(ix,iy) = norm(0.5_R8*(crx(ix,iy,1)+crx(ix,iy,3))-
     &                     0.5_R8*(crx(ix,iy,0)+crx(ix,iy,2)),
     &                     0.5_R8*(cry(ix,iy,1)+cry(ix,iy,3))-
     &                     0.5_R8*(cry(ix,iy,0)+cry(ix,iy,2)))
          hy(ix,iy) = norm(0.5_R8*(crx(ix,iy,2)+crx(ix,iy,3))-
     &                     0.5_R8*(crx(ix,iy,0)+crx(ix,iy,1)),
     &                     0.5_R8*(cry(ix,iy,2)+cry(ix,iy,3))-
     &                     0.5_R8*(cry(ix,iy,0)+cry(ix,iy,1))) 
          cell_width(ix,iy)  = hx(ix,iy)
          cell_height(ix,iy) = hy(ix,iy)

          ! Distances centroid to faces
          hc(ix,iy,0) = 0.5*hx(ix,iy)
          hc(ix,iy,1) = 0.5*hy(ix,iy)
          hc(ix,iy,2) = 0.5*hx(ix,iy) 
          hc(ix,iy,3) = 0.5*hy(ix,iy) 

          ! Distances along faces (only W/S)
          ht(ix,iy,0) = dist(ix,iy,0,2)
          ht(ix,iy,1) = dist(ix,iy,0,1)

          ! Unit vectors of connector lines cell center to face centers
          t0 = (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)) 
          t1 = (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2))
          ecc(ix,iy,0,0) = t0/norm(t0,t1)
          ecc(ix,iy,0,1) = t1/norm(t0,t1)
          ecc(ix,iy,2,0) = t0/norm(t0,t1)
          ecc(ix,iy,2,1) = t1/norm(t0,t1)

          t0 = (crx(ix,iy,2)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,1)) 
          t1 = (cry(ix,iy,2)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,1))
          ecc(ix,iy,1,0) = t0/norm(t0,t1)
          ecc(ix,iy,1,1) = t1/norm(t0,t1)
          ecc(ix,iy,3,0) = t0/norm(t0,t1)
          ecc(ix,iy,3,1) = t1/norm(t0,t1)

          ! Unit vector of poloidal field in cell center
          if ((cflag(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL) .or.
     &         classical ) then 
            ! Assumptions: - internal cells have N and S faces aligned
            !              - for classical grids, the above is true for all cells
            ! Set ebc equal to ecc for numerical accuary
            ebc(ix,iy,0) = ecc(ix,iy,0,0)*bb(ix,iy,0)/abs(bb(ix,iy,0))
            ebc(ix,iy,1) = ecc(ix,iy,0,1)*bb(ix,iy,0)/abs(bb(ix,iy,0))
          else
            ! For arbitrary (quad) cells, compute as average between 
            ! poloidal directions at vertices
            t0 = 0.25_R8*(fpsi(ix,iy,0,1) + fpsi(ix,iy,1,1) + 
     &                    fpsi(ix,iy,2,1) + fpsi(ix,iy,3,1))
            t1 = 0.25_R8*(fpsi(ix,iy,0,2) + fpsi(ix,iy,1,2) + 
     &                    fpsi(ix,iy,2,2) + fpsi(ix,iy,3,2))
            ebc(ix,iy,0) = t0/norm(t0,t1)
            ebc(ix,iy,1) = t1/norm(t0,t1)
          endif

          ! complementary angle between x and y
          qz(ix,iy,0) = sinang (
     &      (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &      (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &      (crx(ix,iy,2)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,1)),
     &      (cry(ix,iy,2)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,1)))
          qz(ix,iy,1) = cosang (
     &      (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &      (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &      (crx(ix,iy,2)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,1)),
     &      (cry(ix,iy,2)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,1)))

        case(CGEO_TRIA_NOLEFT, CGEO_TRIA_NORIGHT,
     &       CGEO_TRIA_NOBOT,  CGEO_TRIA_NOTOP)
cwdk      General approach for arbitrary triangle cases:
cwdk      1) identify vertices with max and min psi values psi_min, psi_max
cwdk      2) assume the triangle lies on the flux surface with psic = 0.5*(psi_min + psi_max)
cwdk         => through middle of the corresponding face; defines 1st intersection
cwdk      3) identify 2nd intersection of this surface with one of the other two faces
cwdk      4) centroid: average between the intersections with faces

          ! Identify vertices with min and max psi values
          i1 = minloc(fpsi(ix,iy,0:3,0),1)-1
          i2 = maxloc(fpsi(ix,iy,0:3,0),1)-1
          psic = 0.5_R8*(fpsi(ix,iy,i1,0) + fpsi(ix,iy,i2,0))

          ! Loop over faces to find which ones are intersected by the psic surface
          ! (Loop includes degenerate face; should not cause problems)
          ii = 0
          if (psic.ge.min(fpsi(ix,iy,0,0),fpsi(ix,iy,2,0)) .and.
     &        psic.lt.max(fpsi(ix,iy,0,0),fpsi(ix,iy,2,0))) then ! intersection with W
            t0 = (psic - fpsi(ix,iy,0,0))/
     &           (fpsi(ix,iy,2,0) - fpsi(ix,iy,0,0))
            isx(ii) = crx(ix,iy,0) + t0*(crx(ix,iy,2) - crx(ix,iy,0)) ! x-co intersection
            isy(ii) = cry(ix,iy,0) + t0*(cry(ix,iy,2) - cry(ix,iy,0)) ! y-co intersection
            ii = ii + 1
          endif
          if (psic.ge.min(fpsi(ix,iy,0,0),fpsi(ix,iy,1,0)) .and.
     &        psic.lt.max(fpsi(ix,iy,0,0),fpsi(ix,iy,1,0))) then ! intersection with S
            t0 = (psic - fpsi(ix,iy,0,0))/
     &           (fpsi(ix,iy,1,0) - fpsi(ix,iy,0,0))
            isx(ii) = crx(ix,iy,0) + t0*(crx(ix,iy,1) - crx(ix,iy,0)) ! x-co intersection
            isy(ii) = cry(ix,iy,0) + t0*(cry(ix,iy,1) - cry(ix,iy,0)) ! y-co intersection
            ii = ii + 1
          endif
          if (psic.ge.min(fpsi(ix,iy,1,0),fpsi(ix,iy,3,0)) .and.
     &        psic.lt.max(fpsi(ix,iy,1,0),fpsi(ix,iy,3,0))) then ! intersection with E
            t0 = (psic - fpsi(ix,iy,1,0))/
     &           (fpsi(ix,iy,3,0) - fpsi(ix,iy,1,0))
            isx(ii) = crx(ix,iy,1) + t0*(crx(ix,iy,3) - crx(ix,iy,1)) ! x-co intersection
            isy(ii) = cry(ix,iy,1) + t0*(cry(ix,iy,3) - cry(ix,iy,1)) ! y-co intersection
            ii = ii + 1
          endif
          if (psic.ge.min(fpsi(ix,iy,2,0),fpsi(ix,iy,3,0)) .and.
     &        psic.lt.max(fpsi(ix,iy,2,0),fpsi(ix,iy,3,0))) then ! intersection with N
            t0 = (psic - fpsi(ix,iy,2,0))/
     &           (fpsi(ix,iy,3,0) - fpsi(ix,iy,2,0))
            isx(ii) = crx(ix,iy,2) + t0*(crx(ix,iy,3) - crx(ix,iy,2)) ! x-co intersection
            isy(ii) = cry(ix,iy,2) + t0*(cry(ix,iy,3) - cry(ix,iy,2)) ! y-co intersection
            ii = ii + 1
          endif
          if (ii.ne.2) then
            write(*,*) 'Problem with intersections in triangle  ',ix,iy
            write(*,*) 'Number of intersections found: ',ii
            call xerrab('Intersections in triangle not found!')
          endif

          ! Compute centroid
          centroid(0) = 0.5_R8*sum(isx)
          centroid(1) = 0.5_R8*sum(isy)

          ! hx:   distance between intersection points
          ! hy:   half distance between points of min/max psi
          ! Note: hx, hy not needed by unstructured solver!
          hx(ix,iy) = norm(isx(1)-isx(0),isy(1)-isy(0))
          hy(ix,iy) = 0.5_R8*dist(ix,iy,i1,i2)
          cell_width(ix,iy)  = hx(ix,iy)
          cell_height(ix,iy) = hy(ix,iy)

          ! Distances centroid to faces
          hc(ix,iy,0) = norm(
     &      0.5_R8*(crx(ix,iy,0)+crx(ix,iy,2)) - centroid(0),
     &      0.5_R8*(cry(ix,iy,0)+cry(ix,iy,2)) - centroid(1))          
          hc(ix,iy,1) = norm(
     &      0.5_R8*(crx(ix,iy,0)+crx(ix,iy,1)) - centroid(0),
     &      0.5_R8*(cry(ix,iy,0)+cry(ix,iy,1)) - centroid(1))
          hc(ix,iy,2) = norm(
     &      0.5_R8*(crx(ix,iy,1)+crx(ix,iy,3)) - centroid(0),
     &      0.5_R8*(cry(ix,iy,1)+cry(ix,iy,3)) - centroid(1))
          hc(ix,iy,3) = norm(
     &      0.5_R8*(crx(ix,iy,2)+crx(ix,iy,3)) - centroid(0),
     &      0.5_R8*(cry(ix,iy,2)+cry(ix,iy,3)) - centroid(1))

          ! Distances along faces (only W/S)
          ht(ix,iy,0) = dist(ix,iy,0,2)
          ht(ix,iy,1) = dist(ix,iy,0,1)

          ! Unit vectors of connector lines cell center to face centers
          t0 = centroid(0) - 0.5_R8*(crx(ix,iy,0)+crx(ix,iy,2)) 
          t1 = centroid(1) - 0.5_R8*(cry(ix,iy,0)+cry(ix,iy,2))
          ecc(ix,iy,0,0) = t0/norm(t0,t1)
          ecc(ix,iy,0,1) = t1/norm(t0,t1)

          t0 = centroid(0) - 0.5_R8*(crx(ix,iy,0)+crx(ix,iy,1)) 
          t1 = centroid(1) - 0.5_R8*(cry(ix,iy,0)+cry(ix,iy,1))
          ecc(ix,iy,1,0) = t0/norm(t0,t1)
          ecc(ix,iy,1,1) = t1/norm(t0,t1)

          t0 = 0.5_R8*(crx(ix,iy,1)+crx(ix,iy,3)) - centroid(0) 
          t1 = 0.5_R8*(cry(ix,iy,1)+cry(ix,iy,3)) - centroid(1)
          ecc(ix,iy,2,0) = t0/norm(t0,t1)
          ecc(ix,iy,2,1) = t1/norm(t0,t1)

          t0 = 0.5_R8*(crx(ix,iy,2)+crx(ix,iy,3)) - centroid(0) 
          t1 = 0.5_R8*(cry(ix,iy,2)+cry(ix,iy,3)) - centroid(1)
          ecc(ix,iy,3,0) = t0/norm(t0,t1)
          ecc(ix,iy,3,1) = t1/norm(t0,t1)

          ! Unit vector of poloidal field in cell center
          ! WG-TMP: based on simple geom. avg. (consistency with carre_b2agbb)
          t0 = 0.25_R8*(fpsi(ix,iy,0,1) + fpsi(ix,iy,1,1) + 
     &                  fpsi(ix,iy,2,1) + fpsi(ix,iy,3,1))
          t1 = 0.25_R8*(fpsi(ix,iy,0,2) + fpsi(ix,iy,1,2) + 
     &                  fpsi(ix,iy,2,2) + fpsi(ix,iy,3,2))
          ebc(ix,iy,0) = t0/norm(t0,t1)
          ebc(ix,iy,1) = t1/norm(t0,t1)

          ! Note: qz loses meaning for a triangle; not computed

        end select
        if(qz(ix,iy,1).le.0.0_R8) then
          write(*,'(a,1x,1f9.6,2i4)')
     1     'Problem in b2agmt: qz(ix,iy,1) <= 0',qz(ix,iy,1),ix,iy
          write(*,'(a,1p,4g13.6)') 'crx = ',(crx(ix,iy,i),i=0,3)
          write(*,'(a,1p,4g13.6)') 'cry = ',(cry(ix,iy,i),i=0,3)
        endif
       enddo
      enddo

*   ..compute qc, qs, qac, qas at cell faces
      ! Requirements for US solver:
      ! qc:  cosine of angle between connector line and face normal
      ! qs:    sine of angle between connector line and face normal
      !      => qc and qs will translate into fcQgam
      ! qac: cosine of angle between poloidal field and face normal
      ! qas:   sine of angle between poloidal field and face normal
      !     => qac and qas will translate into fcQalf
      qc  = 1.0_R8
      qs  = 0.0_R8
      qac = 0.0_R8
      qas = 1.0_R8
*   ..poloidal faces
cwdk  implementation based on b2agmt_qc_style.eq.3 (?): inverse distance
      do ix = -1, nx
       do iy = -1, ny
        if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
        if (leftix(ix,iy).eq.-2 .or. 
     &      leftix(ix,iy).eq.NO_CONNECTIVITY) cycle

        ! unit vectors at left faces (not re-normalized)
        if ((cflag(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL .and.
     &       cflag(leftix(ix,iy),leftiy(ix,iy),CELLFLAG_TYPE) == 
     &       GRID_INTERNAL) .or. classical ) then 
          
          ! Implicit assumption made here: internal cells are quads
          ! Standard treatment for magnetically rectangular cells
          ! Backwards compatible
          ebf(0) = 
     &       (hc(ix,iy,0)*ebc(leftix(ix,iy),leftiy(ix,iy),0)+
     &        hc(leftix(ix,iy),leftiy(ix,iy),2)*ebc(ix,iy,0))/
     &       (hc(leftix(ix,iy),leftiy(ix,iy),2)+hc(ix,iy,0))
          ebf(1) =
     &       (hc(ix,iy,0)*ebc(leftix(ix,iy),leftiy(ix,iy),1)+
     &        hc(leftix(ix,iy),leftiy(ix,iy),2)*ebc(ix,iy,1))/
     &       (hc(leftix(ix,iy),leftiy(ix,iy),2)+hc(ix,iy,0))

          ! Connector line is purely poloidal
          ecf(0) = 
     &       (hc(ix,iy,0)*ecc(leftix(ix,iy),leftiy(ix,iy),2,0)+
     &        hc(leftix(ix,iy),leftiy(ix,iy),2)*ecc(ix,iy,0,0))/
     &       (hc(leftix(ix,iy),leftiy(ix,iy),2)+hc(ix,iy,0))
          ecf(1) =
     &       (hc(ix,iy,0)*ecc(leftix(ix,iy),leftiy(ix,iy),2,1)+
     &        hc(leftix(ix,iy),leftiy(ix,iy),2)*ecc(ix,iy,0,1))/
     &       (hc(leftix(ix,iy),leftiy(ix,iy),2)+hc(ix,iy,0))

        else
          ! "Special" faces (not between aligned quads)
          ! Compute local field direction at face from fpsi
          t0 = 0.5_R8*(fpsi(ix,iy,0,1) + fpsi(ix,iy,2,1))
          t1 = 0.5_R8*(fpsi(ix,iy,0,2) + fpsi(ix,iy,2,2))
          ebf(0) = t0/norm(t0,t1)
          ebf(1) = t1/norm(t0,t1)

          ! Orientation of connector line at faces
          ecf(0) = 
     &       (hc(ix,iy,0)*ecc(ix,iy,0,0)+
     &        hc(leftix(ix,iy),leftiy(ix,iy),2)*
     &       ecc(leftix(ix,iy),leftiy(ix,iy),2,0))/
     &       (hc(leftix(ix,iy),leftiy(ix,iy),2)+hc(ix,iy,0))
          ecf(1) =
     &       (hc(ix,iy,0)*ecc(ix,iy,0,1)+
     &        hc(leftix(ix,iy),leftiy(ix,iy),2)*
     &       ecc(leftix(ix,iy),leftiy(ix,iy),2,1))/
     &       (hc(leftix(ix,iy),leftiy(ix,iy),2)+hc(ix,iy,0))

        endif

        ! angles at left face
        qc(ix,iy,0) = cosang ( ecf(0), ecf(1), 
     &    crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))
        qs(ix,iy,0) = sinang ( ecf(0), ecf(1), 
     &    crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))

        qac(ix,iy,0) = cosang ( ebf(0), ebf(1), 
     &    crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))
        qas(ix,iy,0) = sinang ( ebf(0), ebf(1), 
     &    crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))

        if(qc(ix,iy,0).le.0.0_R8) then
          write(*,'(a,1f12.6,2i5)') 
     1     'Problem in b2agmt: qc(ix,iy,0) <= 0',qc(ix,iy,0),ix,iy
          write(*,'(a,4f12.6)') 'crx = ',crx(ix,iy,0:3)
          write(*,'(a,4f12.6)') 'cry = ',cry(ix,iy,0:3)
          write(*,*) "Cell geometry is ", cgeo
        endif
       enddo
      enddo

*   ..radial faces
      do ix = -1, nx
       do iy = -1, ny
        if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
        if (bottomiy(ix,iy).eq.-2 .or. 
     &      bottomiy(ix,iy).eq.NO_CONNECTIVITY) cycle

        ! unit vectors at bottom faces (not re-normalized)
        if ((cflag(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL .and.
     &       cflag(bottomix(ix,iy),bottomiy(ix,iy),CELLFLAG_TYPE) == 
     &       GRID_INTERNAL) .or. classical ) then 

          ! Implicit assumption made here: internal cells are quads
          ! Standard treatment for magnetically rectangular cells
          ! First: estimate of ebf by interpolation
          ebf(0) = 
     &       (hc(ix,iy,1)*ebc(bottomix(ix,iy),bottomiy(ix,iy),0)+
     &        hc(bottomix(ix,iy),bottomiy(ix,iy),3)*ebc(ix,iy,0))/
     &       (hc(bottomix(ix,iy),bottomiy(ix,iy),2)+hc(ix,iy,0))
          ebf(1) = 
     &       (hc(ix,iy,1)*ebc(bottomix(ix,iy),bottomiy(ix,iy),1)+
     &        hc(bottomix(ix,iy),bottomiy(ix,iy),3)*ebc(ix,iy,1))/
     &       (hc(bottomix(ix,iy),bottomiy(ix,iy),2)+hc(ix,iy,0))

          ! Assume perfectly aligned face
          qac(ix,iy,1) = 0.0_R8
          qas(ix,iy,1) = sinang ( ebf(0), ebf(1), 
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0))
          qas(ix,iy,1) = sign(1.0_R8,qas(ix,iy,1))

          ! qc and qs for backwards compatibility
          qc(ix,iy,1) =
     .     (hy(ix,iy)*qz(ix,iy,1) + hy(bottomix(ix,iy),bottomiy(ix,iy))*
     .      qz(bottomix(ix,iy),bottomiy(ix,iy),1))/ 
     .     (hy(ix,iy)  + hy(bottomix(ix,iy),bottomiy(ix,iy)))
          qs(ix,iy,1) =
     .     (hy(ix,iy)*qz(ix,iy,0) + hy(bottomix(ix,iy),bottomiy(ix,iy))*
     .      qz(bottomix(ix,iy),bottomiy(ix,iy),0))/ 
     .     (hy(ix,iy)  + hy(bottomix(ix,iy),bottomiy(ix,iy)))

        else
          ! "Special" faces (not between aligned quads)
          ! Compute local field direction at face from fpsi
          t0 = 0.5_R8*(fpsi(ix,iy,0,1) + fpsi(ix,iy,1,1))
          t1 = 0.5_R8*(fpsi(ix,iy,0,2) + fpsi(ix,iy,1,2))
          ebf(0) = t0/norm(t0,t1)
          ebf(1) = t1/norm(t0,t1)

          ecf(0) = 
     &       (hc(ix,iy,1)*ecc(ix,iy,1,0)+
     &        hc(bottomix(ix,iy),bottomiy(ix,iy),3)*
     &       ecc(bottomix(ix,iy),bottomiy(ix,iy),3,0))/
     &       (hc(bottomix(ix,iy),bottomiy(ix,iy),3)+hc(ix,iy,1))
          ecf(1) = 
     &       (hc(ix,iy,1)*ecc(ix,iy,1,1)+
     &        hc(bottomix(ix,iy),bottomiy(ix,iy),3)*
     &       ecc(bottomix(ix,iy),bottomiy(ix,iy),3,1))/
     &       (hc(bottomix(ix,iy),bottomiy(ix,iy),3)+hc(ix,iy,1))

          ! angles at bottom face
          qc(ix,iy,1) = cosang ( ecf(0), ecf(1), 
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0))
          qs(ix,iy,1) = sinang ( ecf(0), ecf(1), 
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0))

          qac(ix,iy,1) = cosang ( ebf(0), ebf(1), 
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0))
          qas(ix,iy,1) = sinang ( ebf(0), ebf(1), 
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0))

        endif

        if(qc(ix,iy,1).le.0.0_R8) then
          write(*,'(a,1f12.6,2i4)')
     1     'Problem in b2agmt: qc(ix,iy,1) <= 0',qc(ix,iy,1),ix,iy
          write(*,'(a,4f13.6)') 'crx = ',(crx(ix,iy,i),i=0,3)
          write(*,'(a,4f13.6)') 'cry = ',(cry(ix,iy,i),i=0,3)
        endif
       enddo
      enddo
*   ..unit vector of total magnetic field
      ! until this point, ebc was unit vector of poloidal field
      ! now, correct for toroidal field
      do iy = -1, ny
        do ix = -1, nx
          if (bb(ix,iy,3) == 0.0_R8) then
            ebc(ix,iy,:) = 0.0_R8
          else
            ebc(ix,iy,0) = abs(bb(ix,iy,0))/bb(ix,iy,3)*ebc(ix,iy,0)
            ebc(ix,iy,1) = abs(bb(ix,iy,0))/bb(ix,iy,3)*ebc(ix,iy,1)
            ebc(ix,iy,2) = bb(ix,iy,2)/bb(ix,iy,3)
          endif
        enddo
      enddo
*   ..test sign of qz(,,1) and qc(,)
*     (A negative entry in qz(,,1) corresponds to a left-handed
*     orientation of the coordinate system, which is not allowed in
*     this code.)
      call b2xvsg1(nx, ny, qz(-1,-1,1), 1, 'qz1', '.gt.')
      call b2xvsg1(nx, ny, qc(-1,-1,0), 1, 'qc', '.gt.')
!      call b2xvsg1(nx, ny, qc(-1,-1,1), 1, 'qc', '.gt.')
*   ..compute vol, gs
      vol = INVALID_DOUBLE
      gs(:,:,0) = 0.0_R8
      gs(:,:,1) = 0.0_R8
      gs(:,:,2) = INVALID_DOUBLE
      if (isymm.eq.0) then
*    ..slab geometry
       do iy = -1, ny
        do ix = -1, nx
         if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
         gs(ix,iy,2) = quadArea(  ! the vertices must be in sequence
     &     crx(ix,iy,0), cry(ix,iy,0), crx(ix,iy,1), cry(ix,iy,1),
     &     crx(ix,iy,3), cry(ix,iy,3), crx(ix,iy,2), cry(ix,iy,2))
         vol(ix,iy) = gs(ix,iy,2)*width
        enddo
       enddo
       do iy = -1, ny
        do ix = -1, nx
         if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
         if(isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) 
     &    gs(ix,iy,0) = dist(ix,iy,0,2)*width
         if(isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) 
     &    gs(ix,iy,1) = dist(ix,iy,0,1)*width
        enddo
       enddo
      else if (isymm.eq.1 .or. isymm.eq.2) then
*    ..case of rotational symmetry about crx=0 axis
       do iy = -1, ny
        do ix = -1, nx
         if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
         gs(ix,iy,2) = quadArea(  ! the vertices must be in sequence
     &     crx(ix,iy,0), cry(ix,iy,0), crx(ix,iy,1), cry(ix,iy,1),
     &     crx(ix,iy,3), cry(ix,iy,3), crx(ix,iy,2), cry(ix,iy,2))
         if (volhmk.eq.0) then
          centroid(0) = 0.25_R8*
     &     (crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+crx(ix,iy,3))
          centroid(1) = 0.25_R8*
     &     (cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+cry(ix,iy,3))
         else
          centroid = quadCentroid(crx(ix,iy,0),cry(ix,iy,0),
     &                            crx(ix,iy,1),cry(ix,iy,1),
     &                            crx(ix,iy,2),cry(ix,iy,2),
     &                            crx(ix,iy,3),cry(ix,iy,3)) 
         endif
         vol(ix,iy) = gs(ix,iy,2)*2.0_R8*pi*centroid(0)
        enddo
       enddo
       do iy = -1, ny
        do ix = -1, nx
         if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
         if (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) 
     &    gs(ix,iy,0) = pi*dist(ix,iy,0,2)*(crx(ix,iy,0)+crx(ix,iy,2))
         if (isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) 
     &    gs(ix,iy,1) = pi*dist(ix,iy,0,1)*(crx(ix,iy,0)+crx(ix,iy,1))
        enddo
       enddo
      else if (isymm.eq.3 .or. isymm.eq.4) then
*    ..case of rotational symmetry about cry=0 axis
       do iy = -1, ny
        do ix = -1, nx
         if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
         gs(ix,iy,2) = quadArea(  ! the vertices must be in sequence
     &     crx(ix,iy,0), cry(ix,iy,0), crx(ix,iy,1), cry(ix,iy,1),
     &     crx(ix,iy,3), cry(ix,iy,3), crx(ix,iy,2), cry(ix,iy,2))
         if (volhmk.eq.0) then
          centroid(0) = 0.25_R8*
     &     (crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+crx(ix,iy,3))
          centroid(1) = 0.25_R8*
     &     (cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+cry(ix,iy,3))
         else
          centroid = quadCentroid(crx(ix,iy,0),cry(ix,iy,0),
     &                            crx(ix,iy,1),cry(ix,iy,1),
     &                            crx(ix,iy,2),cry(ix,iy,2),
     &                            crx(ix,iy,3),cry(ix,iy,3)) 
         endif
         vol(ix,iy) = gs(ix,iy,2)*2.0_R8*pi*centroid(1)
        enddo
       enddo
       do iy = -1, ny
        do ix = -1, nx
         if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
         if (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) 
     &    gs(ix,iy,0) = pi*dist(ix,iy,0,2)*(cry(ix,iy,0)+cry(ix,iy,2))
         if (isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) 
     &    gs(ix,iy,1) = pi*dist(ix,iy,0,1)*(cry(ix,iy,0)+cry(ix,iy,1))
        enddo
       enddo
      endif
*   ..test sign of vol, gs
      call b2xvsg1(nx, ny, vol(-1,-1),  1, 'vol', '.gt.')
!      call b2xvsg1(nx, ny, gs(-1,-1,0), 1, 'gsx', '.gt.')
!      call b2xvsg1(nx, ny, gs(-1,-1,1), 1, 'gsy', '.gt.')
      call b2xvsg1(nx, ny, gs(-1,-1,2), 1, 'gsz', '.gt.')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agmt

      end subroutine b2agmt
