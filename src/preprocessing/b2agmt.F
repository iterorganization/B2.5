*-----------------------------------------------------------------------
*.specification

      subroutine b2agmt (nx, ny, isymm, crx, cry,
     & fpsi, bb, vol, hx, hy, hz, hc, ht, qz, qc, qs, qac, qas, ebc,
     & gs, cell_width, cell_height,
     & leftix, leftiy, bottomix, bottomiy, cflag)
      use b2mod_types
      use b2mod_math
      use b2mod_constants
      use tradui_constants
      use b2mod_connectivity
      use b2mod_cellhelper
      use b2mod_indirect
     , ,only : quadCentroid, quadArea
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isymm
      integer leftix(-1:nx,-1:ny), leftiy(-1:nx,-1:ny),
     &  bottomix(-1:nx,-1:ny), bottomiy(-1:nx,-1:ny),
     &  cflag(-1:nx,-1:ny,GRID_N_CELLFLAGS)
      real (kind=R8) ::
     *  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     &  fpsi(-1:nx,-1:ny,0:3,0:2), bb(-1:nx,-1:ny,0:3)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  hz(-1:nx,-1:ny),
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
      integer ix, iy, i
      integer i1, i2, ii
      integer volhmk, cgeo, cgeon
      real (kind=R8) ::
     *   width, centroid(0:1), t0, t1, psic
      logical classical
      real (kind=R8) ::
     *   ecc(-1:nx,-1:ny,0:3,0:1), ecf(0:1), ebf(0:1),
     *   isx(0:1), isy(0:1),
     *   Xvertices(0:3), Yvertices(0:3)
      integer, save :: style = 0
*   ..procedures
      intrinsic min, max
      external xertst, ipgetr, ipgeti
      external b2xvsg1
*   ..initialisation
      data volhmk/1/, width /1.0_R8/
      save volhmk, width

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agmt')
*   ..read input data
      call ipgeti ('b2agmt_volhmk', volhmk)
      call ipgetr ('b2agmt_1d_width', width)
      call ipgeti ('b2agmt_style', style)
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

        Xvertices(:) = crx(ix,iy,:)
        Yvertices(:) = cry(ix,iy,:)
        cgeo = cellGeoType(Xvertices, Yvertices)
        select case (cgeo)
        case(CGEO_BROKEN)
          write(*,*) 'Problem in b2agmt: broken cell', ix,iy
          write(*,'(a,4f12.6)') 'crx = ',(crx(ix,iy,i),i=0,3)
          write(*,'(a,4f12.6)') 'cry = ',(cry(ix,iy,i),i=0,3)
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
            ! TODO: extend this treatment for boundary cells with both N/S faces aligned
            ebc(ix,iy,0) = ecc(ix,iy,0,0)*bb(ix,iy,0)/abs(bb(ix,iy,0))
            ebc(ix,iy,1) = ecc(ix,iy,0,1)*bb(ix,iy,0)/abs(bb(ix,iy,0))
          else
            ! For arbitrary (quad) cells, compute as average between
            ! poloidal directions at vertices
            t0 =  0.25_R8*(fpsi(ix,iy,0,2) + fpsi(ix,iy,1,2) +
     &                     fpsi(ix,iy,2,2) + fpsi(ix,iy,3,2))
            t1 = -0.25_R8*(fpsi(ix,iy,0,1) + fpsi(ix,iy,1,1) +
     &                     fpsi(ix,iy,2,1) + fpsi(ix,iy,3,1))
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
          t0 =  0.25_R8*(fpsi(ix,iy,0,2) + fpsi(ix,iy,1,2) +
     &                   fpsi(ix,iy,2,2) + fpsi(ix,iy,3,2))
          t1 = -0.25_R8*(fpsi(ix,iy,0,1) + fpsi(ix,iy,1,1) +
     &                   fpsi(ix,iy,2,1) + fpsi(ix,iy,3,1))
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
        Xvertices(:) = crx(ix,iy,:)
        Yvertices(:) = cry(iy,iy,:)
        cgeo = cellGeoType(Xvertices, Yvertices)
        Xvertices(:) = crx(leftix(ix,iy),leftiy(ix,iy),:)
        Yvertices(:) = cry(leftix(ix,iy),leftiy(ix,iy),:)
        cgeon = cellGeoType(Xvertices, Yvertices)

        ! unit vectors at left faces (not re-normalized)
        if (classical.and.(style.eq.0)) then 
          
          ! Inverse distance interpolation
          ! Backwards compatible
          ebf(0) = 
     &       (hc(ix,iy,0)*ebc(leftix(ix,iy),leftiy(ix,iy),0)+
     &        hc(leftix(ix,iy),leftiy(ix,iy),2)*ebc(ix,iy,0))/
     &       (hc(leftix(ix,iy),leftiy(ix,iy),2)+hc(ix,iy,0))
          ebf(1) =
     &       (hc(ix,iy,0)*ebc(leftix(ix,iy),leftiy(ix,iy),1)+
     &        hc(leftix(ix,iy),leftiy(ix,iy),2)*ebc(ix,iy,1))/
     &       (hc(leftix(ix,iy),leftiy(ix,iy),2)+hc(ix,iy,0))

          ! (Same) inverse distance interpolation;
          ! Connector line is purely poloidal: ecc = +- ebc
          ! for classical grids!
          ecf(0) = 
     &       (hc(ix,iy,0)*ecc(leftix(ix,iy),leftiy(ix,iy),2,0)+
     &        hc(leftix(ix,iy),leftiy(ix,iy),2)*ecc(ix,iy,0,0))/
     &       (hc(leftix(ix,iy),leftiy(ix,iy),2)+hc(ix,iy,0))
          ecf(1) =
     &       (hc(ix,iy,0)*ecc(leftix(ix,iy),leftiy(ix,iy),2,1)+
     &        hc(leftix(ix,iy),leftiy(ix,iy),2)*ecc(ix,iy,0,1))/
     &       (hc(leftix(ix,iy),leftiy(ix,iy),2)+hc(ix,iy,0))

        elseif (
     &       (cflag(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL .and.         ! 1) between 2 internal cells
     &        cflag(leftix(ix,iy),leftiy(ix,iy),CELLFLAG_TYPE) == 
     &        GRID_INTERNAL) .or.
     &       (cgeo == CGEO_QUAD .and. cgeon == CGEO_QUAD .and.          ! 2) between 2 quads with aligned N/S faces
     &        cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(ix,iy,CELLFLAG_TOPFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_TOPFACE) .le. 0) .or.
     &       (cgeo == CGEO_QUAD .and. cgeon == CGEO_TRIA_NOTOP .and.    ! 3) between quad and tria with aligned N/S faces
     &        cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(ix,iy,CELLFLAG_TOPFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_BOTTOMFACE) .le. 0) .or.
     &       (cgeo == CGEO_QUAD .and. cgeon == CGEO_TRIA_NOBOT .and.
     &        cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(ix,iy,CELLFLAG_TOPFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_TOPFACE) .le. 0) .or.
     &       (cgeo == CGEO_TRIA_NOTOP .and. cgeon == CGEO_QUAD .and.    ! 4) between tria and quad with aligned N/S faces
     &        cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_TOPFACE) .le. 0) .or.
     &       (cgeo == CGEO_TRIA_NOBOT .and. cgeon == CGEO_QUAD .and.
     &        cflag(ix,iy,CELLFLAG_TOPFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_TOPFACE) .le. 0) .or.
     &       (cgeo == CGEO_QUAD .and. cgeon == CGEO_TRIA_NOLEFT .and.   ! 5) between quad and tria with aligned N/S faces
     &        cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(ix,iy,CELLFLAG_TOPFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_BOTTOMFACE) .le. 0) .or.
     &       (cgeo == CGEO_QUAD .and. cgeon == CGEO_TRIA_NOLEFT .and.
     &        cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(ix,iy,CELLFLAG_TOPFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_TOPFACE) .le. 0) .or.
     &       (cgeo == CGEO_TRIA_NORIGHT .and. cgeon == CGEO_QUAD .and.  ! 6) between tria and quad with aligned N/S faces
     &        cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_TOPFACE) .le. 0) .or.
     &       (cgeo == CGEO_TRIA_NORIGHT .and. cgeon == CGEO_QUAD .and.
     &        cflag(ix,iy,CELLFLAG_TOPFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_TOPFACE) .le. 0) .or.
     &       (cgeo == CGEO_TRIA_NORIGHT .and.
     &        cgeon == CGEO_TRIA_NOLEFT .and.                           ! 7) between 2 trias with aligned S faces
     &        cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_BOTTOMFACE) .le. 0) .or.
     &       (cgeo == CGEO_TRIA_NORIGHT .and.
     &        cgeon == CGEO_TRIA_NOTOP .and.
     &        cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_BOTTOMFACE) .le. 0) .or.
     &       (cgeo == CGEO_TRIA_NOTOP .and.
     &        cgeon == CGEO_TRIA_NOLEFT .and.
     &        cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_BOTTOMFACE) .le. 0) .or.
     &       (cgeo == CGEO_TRIA_NOTOP .and.
     &        cgeon == CGEO_TRIA_NOTOP .and.
     &        cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_BOTTOMFACE) .le. 0) .or.
     &       (cgeo == CGEO_TRIA_NORIGHT .and.
     &        cgeon == CGEO_TRIA_NOLEFT .and.                           ! 8) between 2 trias with aligned N faces
     &        cflag(ix,iy,CELLFLAG_TOPFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_TOPFACE) .le. 0) .or.
     &       (cgeo == CGEO_TRIA_NORIGHT .and.
     &        cgeon == CGEO_TRIA_NOBOT .and.
     &        cflag(ix,iy,CELLFLAG_TOPFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_TOPFACE) .le. 0) .or.
     &       (cgeo == CGEO_TRIA_NOBOT .and.
     &        cgeon == CGEO_TRIA_NOLEFT .and.
     &        cflag(ix,iy,CELLFLAG_TOPFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_TOPFACE) .le. 0) .or.
     &       (cgeo == CGEO_TRIA_NOBOT .and.
     &        cgeon == CGEO_TRIA_NOBOT .and.
     &        cflag(ix,iy,CELLFLAG_TOPFACE) .le. 0 .and.
     &        cflag(leftix(ix,iy),leftiy(ix,iy),
     &        CELLFLAG_TOPFACE) .le. 0) ) then
          ! Faces between "regular" internal cells, assumed to 
          ! be quads with N/S faces aligned to magnetic field.
          ! Also: 
          ! * faces between any 2 quad cells where all N/S faces are aligned
          ! * faces between a guard cell and a triangle cell with aligned N or S

          ! Compute local field direction at face from fpsi
          t0 =  0.5_R8*(fpsi(ix,iy,0,2) + fpsi(ix,iy,2,2))
          t1 = -0.5_R8*(fpsi(ix,iy,0,1) + fpsi(ix,iy,2,1))
          ebf(0) = t0/norm(t0,t1)
          ebf(1) = t1/norm(t0,t1)

          ! Connector line is purely poloidal: equal to ebf,
          ! but possibly with sign difference
          t0 = ebf(0)*(cry(ix,iy,2) - cry(ix,iy,0)) -
     &         ebf(1)*(crx(ix,iy,2) - crx(ix,iy,0))
          if (t0.ge.0.0_R8) then
            ecf = ebf
          else
            ecf = -ebf
          endif
        else
          ! "Special" faces (not between aligned quads)
          ! Compute local field direction at face from fpsi
          t0 =  0.5_R8*(fpsi(ix,iy,0,2) + fpsi(ix,iy,2,2))
          t1 = -0.5_R8*(fpsi(ix,iy,0,1) + fpsi(ix,iy,2,1))
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
        qas(ix,iy,0) = -sinang ( ebf(0), ebf(1), 
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
        Xvertices(:) = crx(ix,iy,:)
        Yvertices(:) = cry(iy,iy,:)
        cgeo = cellGeoType(Xvertices, Yvertices)
        Xvertices(:) = crx(bottomix(ix,iy),bottomiy(ix,iy),:)
        Yvertices(:) = cry(bottomix(ix,iy),bottomiy(ix,iy),:)
        cgeon = cellGeoType(Xvertices, Yvertices)

        ! unit vectors at bottom faces (not re-normalized)
        if (classical.and.(style.eq.0)) then

          ! Backwards compatibility option.
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
          qas(ix,iy,1) = sinang (
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0),
     &      ebf(0), ebf(1))
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
          qs(ix,iy,1) = sign(sqrt(1.0_R8 - qc(ix,iy,1)**2),qs(ix,iy,1))

        elseif (cflag(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL .or.
     &          cflag(bottomix(ix,iy),bottomiy(ix,iy),CELLFLAG_TYPE) ==
     &          GRID_INTERNAL .or.
     &          (cgeo  == CGEO_QUAD .and. 
     &           cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0) .or.
     &          (cgeon == CGEO_QUAD .and. 
     &           cflag(bottomix(ix,iy),bottomiy(ix,iy),
     &           CELLFLAG_TOPFACE) .le. 0) .or.
     &          (cgeo  == CGEO_TRIA_NOLEFT  .and.
     &           cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0) .or.
     &          (cgeo  == CGEO_TRIA_NORIGHT .and.
     &           cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0) .or.
     &          (cgeo  == CGEO_TRIA_NOTOP   .and.
     &           cflag(ix,iy,CELLFLAG_BOTTOMFACE) .le. 0) .or.
     &          (cgeon == CGEO_TRIA_NOLEFT  .and.
     &           cflag(bottomix(ix,iy),bottomiy(ix,iy),
     &           CELLFLAG_TOPFACE) .le. 0) .or.
     &          (cgeon == CGEO_TRIA_NORIGHT .and.
     &           cflag(bottomix(ix,iy),bottomiy(ix,iy),
     &           CELLFLAG_TOPFACE) .le. 0) .or.
     &          (cgeon == CGEO_TRIA_NOBOT   .and.
     &           cflag(bottomix(ix,iy),bottomiy(ix,iy),
     &           CELLFLAG_TOPFACE) .le. 0)) then 

          ! Treatment for aligned faces

          ! Compute local field direction at face from fpsi
          t0 =  0.5_R8*(fpsi(ix,iy,0,2) + fpsi(ix,iy,1,2))
          t1 = -0.5_R8*(fpsi(ix,iy,0,1) + fpsi(ix,iy,1,1))
          ebf(0) = t0/norm(t0,t1)
          ebf(1) = t1/norm(t0,t1)

          ! Assume perfectly aligned face for qac and qas
          qac(ix,iy,1) = 0.0_R8
          qas(ix,iy,1) = sinang (
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0),
     &      ebf(0), ebf(1))
          qas(ix,iy,1) = sign(1.0_R8,qas(ix,iy,1))

          ! Orientation of connector line at faces
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
          qc(ix,iy,1) = cosang (
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0),
     &      ecf(0), ecf(1))
          qs(ix,iy,1) = sinang (
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0),
     &      ecf(0), ecf(1))

        elseif ((cgeo  == CGEO_TRIA_NOLEFT  .and.
     &           cflag(ix,iy,CELLFLAG_TOPFACE) == 0) .or.
     &          (cgeo  == CGEO_TRIA_NORIGHT .and.
     &           cflag(ix,iy,CELLFLAG_TOPFACE) == 0) .or.
     &          (cgeon == CGEO_TRIA_NOLEFT  .and.
     &           cflag(bottomix(ix,iy),bottomiy(ix,iy),
     &           CELLFLAG_BOTTOMFACE) == 0) .or.
     &          (cgeon == CGEO_TRIA_NORIGHT .and.
     &           cflag(bottomix(ix,iy),bottomiy(ix,iy),
     &           CELLFLAG_BOTTOMFACE) == 0) ) then

          ! Special treatment for certain triangle faces connecting to boundary

          ! Compute local field direction at face from fpsi
          t0 =  0.5_R8*(fpsi(ix,iy,0,2) + fpsi(ix,iy,1,2))
          t1 = -0.5_R8*(fpsi(ix,iy,0,1) + fpsi(ix,iy,1,1))
          ebf(0) = t0/norm(t0,t1)
          ebf(1) = t1/norm(t0,t1)

          ! Connector line is purely poloidal: equal to ebf,
          ! but possibly with sign difference
          t0 = ebf(0)*(cry(ix,iy,0) - cry(ix,iy,1)) -
     &         ebf(1)*(crx(ix,iy,0) - crx(ix,iy,1))
          if (t0.ge.0.0_R8) then
            ecf = ebf
          else
            ecf = -ebf
          endif

          ! angles at bottom face
          qc(ix,iy,1) = cosang (
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0),
     &      ecf(0), ecf(1))
          qs(ix,iy,1) = sinang (
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0),
     &      ecf(0), ecf(1))

          qac(ix,iy,1) = cosang (
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0),
     &      ebf(0), ebf(1))
          qas(ix,iy,1) = sinang (
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0),
     &      ebf(0), ebf(1))

        else
          ! "Special" faces (not between aligned quads)
          ! Compute local field direction at face from fpsi
          t0 =  0.5_R8*(fpsi(ix,iy,0,2) + fpsi(ix,iy,1,2))
          t1 = -0.5_R8*(fpsi(ix,iy,0,1) + fpsi(ix,iy,1,1))
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
          qc(ix,iy,1) = cosang (
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0),
     &      ecf(0), ecf(1))
          qs(ix,iy,1) = sinang (
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0),
     &      ecf(0), ecf(1))

          qac(ix,iy,1) = cosang (
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0),
     &      ebf(0), ebf(1))
          qas(ix,iy,1) = sinang (
     &      crx(ix,iy,1)-crx(ix,iy,0), cry(ix,iy,1)-cry(ix,iy,0),
     &      ebf(0), ebf(1))

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
      call b2xvsg1(nx, ny, qc(-1,-1,1), 1, 'qc', '.gt.')
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
*   ..compute hz
      hz = INVALID_DOUBLE
      do iy = -1, ny
        do ix = -1, nx
          if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
          hz(ix,iy) = vol(ix,iy)/gs(ix,iy,2)
        enddo
      enddo

*   ..test sign of vol, gs
      call b2xvsg1(nx, ny, vol(-1,-1),  1, 'vol', '.gt.')
!      call b2xvsg1(nx, ny, ny, gs(-1,-1,0), 1, 'gsx', '.gt.')
!      call b2xvsg1(nx, ny, gs(-1,-1,1), 1, 'gsy', '.gt.')
      call b2xvsg1(nx, ny, gs(-1,-1,2), 1, 'gsz', '.gt.')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agmt

      contains

      function dist(ix,iy,i,j)
      use b2mod_math
      implicit none
      integer, intent(in) :: ix, iy, i, j
      real(kind=R8) :: dist
      dist = norm(crx(ix,iy,j)-crx(ix,iy,i), cry(ix,iy,j)-cry(ix,iy,i))
      return
      end function dist

      end subroutine b2agmt
