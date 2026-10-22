      module b2mod_chord_interp
!=======================================================================
!     Interpolation of B2.5 fields onto arbitrary (R,Z) points and along
!     straight chords.
!
!     Method:
!
!     * Point location. A point is inside a control volume if it
!       lies in one of the triangles formed by the cell center and
!       each of the cell's faces (fan triangulation from the center,
!       which is interior for a convex cell). This is independent of the
!       vertex ordering. Walking a chord, the previous point's cell is
!       used as a hint and only its 9-point stencil is searched before
!       falling back to a global scan.
!
!     * Volume-centered interpolation. Inverse-distance weighting (IDW)
!       over the containing cell and its face-adjacent neighbours, using
!       the cell-center coordinates geo%cvX/cvY.
!       The weight is 1/distance**power (power=2 by default). Guard cells
!       are excluded as neighbours unless include_guard is set.
!
!     * Face-centered fluxes.
!       chord_profile_fc always returns flux densities: the fluxes are
!       first divided by the projected face area with
!       divide_by_area (geo%fcS*geo%fcQalf). Because of the 9-point
!       stencil the components are not simple scalars per face, so
!       the densities are collapsed to cell centers with
!       intcell(...,intcellP,...) for the poloidal component and
!       intcell(...,intcellR,...) for the radial component, then
!       interpolated like any cell-centered quantity.
!       The two profiles returned are the poloidal and
!       radial flux densities.
!=======================================================================
      use b2mod_types, only : R8
      use b2us_geo,    only : geometry, divide_by_area
      use b2us_map,    only : mapping
      implicit none
      private

      public :: locate_cv, interp_cv_at_rz
      public :: chord_profile_cv, chord_profile_fc
      public :: chords_total_points, chords_geometry
      public :: sample_cv_at_points, facefield_to_cellcenters
      public :: mask_inside_cv, compact_chords

!     Points closer than this (in metres) to a cell center are treated
!     as coincident with it (return the cell value exactly).
      real (kind=R8), parameter :: dist_eps = 1.0e-12_R8

      contains

!-----------------------------------------------------------------------
      subroutine locate_cv (geo, mpg, r, z, iCv, hint)
!     Return in iCv the internal control volume (1..nCi) containing
!     (r,z), or 0 if the point is outside the internal mesh. If present,
!     hint is a cell to try first (together with its 9-point stencil)
!     before a global scan.
      implicit none
      type (geometry), intent (in)  :: geo
      type (mapping),  intent (in)  :: mpg
      real (kind=R8),  intent (in)  :: r, z
      integer,         intent (out) :: iCv
      integer, intent (in), optional :: hint
      integer :: k, c, p0, np, ic
      logical :: inside

      iCv = 0

!     1) the hint cell and its 9-point stencil
      if (present(hint)) then
        if (hint.ge.1 .and. hint.le.mpg%nCi) then
          call point_in_cell(geo, mpg, hint, r, z, inside)
          if (inside) then
            iCv = hint
            return
          end if
          p0 = mpg%cvNvP(hint,1)
          np = mpg%cvNvP(hint,2)
          do k = 0, np-1
            c = mpg%cvNv(p0+k)
            if (c.ge.1 .and. c.le.mpg%nCi) then
              call point_in_cell(geo, mpg, c, r, z, inside)
              if (inside) then
                iCv = c
                return
              end if
            end if
          end do
        end if
      end if

!     2) global scan over the internal cells (cold start / fallback)
      do ic = 1, mpg%nCi
        call point_in_cell(geo, mpg, ic, r, z, inside)
        if (inside) then
          iCv = ic
          return
        end if
      end do

      return
      end subroutine locate_cv

!-----------------------------------------------------------------------
      subroutine interp_cv_at_rz (geo, mpg, fld, r, z, val, iCv, ok,     &
     &                            hint, power, include_guard)
!     Interpolate the cell-centered field fld at the point (r,z) by
!     inverse-distance weighting over the containing cell and its
!     face-adjacent neighbours. On exit ok=.true. and val holds the
!     interpolated value if the point is inside the mesh; otherwise
!     ok=.false., val=0 and iCv=0. iCv returns the containing cell (0 if
!     none), useful as the hint for the next call.
      implicit none
      type (geometry), intent (in)  :: geo
      type (mapping),  intent (in)  :: mpg
      real (kind=R8),  intent (in)  :: fld(mpg%nCv)
      real (kind=R8),  intent (in)  :: r, z
      real (kind=R8),  intent (out) :: val
      integer,         intent (out) :: iCv
      logical,         intent (out) :: ok
      integer,       intent (in), optional :: hint
      real (kind=R8),intent (in), optional :: power
      logical,       intent (in), optional :: include_guard

      real (kind=R8) :: p, wsum, vsum, dist, w
      logical        :: guard
      integer        :: i, iFc, p0, nf, c

      val = 0.0_R8
      ok  = .false.
      iCv = 0

      p = 2.0_R8
      if (present(power)) p = power
      guard = .false.
      if (present(include_guard)) guard = include_guard

      if (present(hint)) then
        call locate_cv(geo, mpg, r, z, iCv, hint)
      else
        call locate_cv(geo, mpg, r, z, iCv)
      end if
      if (iCv.eq.0) return

      wsum = 0.0_R8
      vsum = 0.0_R8

!     containing cell
      dist = sqrt((geo%cvX(iCv)-r)**2 + (geo%cvY(iCv)-z)**2)
      if (dist.lt.dist_eps) then
        val = fld(iCv)
        ok  = .true.
        return
      end if
      w    = 1.0_R8/dist**p
      wsum = wsum + w
      vsum = vsum + w*fld(iCv)

!     face-adjacent neighbours
      p0 = mpg%cvFcP(iCv,1)
      nf = mpg%cvFcP(iCv,2)
      do i = 0, nf-1
        iFc = mpg%cvFc(p0+i)
        c = mpg%fcCv(iFc,1)
        if (c.eq.iCv) c = mpg%fcCv(iFc,2)
        if (c.lt.1) cycle
        if (.not.guard .and. c.gt.mpg%nCi) cycle
        dist = sqrt((geo%cvX(c)-r)**2 + (geo%cvY(c)-z)**2)
        if (dist.lt.dist_eps) then
          val = fld(c)
          ok  = .true.
          return
        end if
        w    = 1.0_R8/dist**p
        wsum = wsum + w
        vsum = vsum + w*fld(c)
      end do

      if (wsum.gt.0.0_R8) then
        val = vsum/wsum
        ok  = .true.
      end if

      return
      end subroutine interp_cv_at_rz

!-----------------------------------------------------------------------
      subroutine chord_profile_cv (geo, mpg, fld, r1, z1, r2, z2, dl,    &
     &                             npts, s, rr, zz, val, ok,             &
     &                             power, include_guard)
!     Sample the cell-centered field fld along the chord from (r1,z1) to
!     (r2,z2) at a spacing of about dl (metres). The chord is divided
!     into n = ceiling(length/dl) equal segments, giving npts=n+1 sample
!     points (both endpoints included). On exit the allocatable arrays
!     hold, for each sample point: s = arc length from (r1,z1), (rr,zz) =
!     coordinates, val = interpolated value, ok = whether the point was
!     inside the mesh (val=0 where ok is .false.).
      implicit none
      type (geometry), intent (in) :: geo
      type (mapping),  intent (in) :: mpg
      real (kind=R8),  intent (in) :: fld(mpg%nCv)
      real (kind=R8),  intent (in) :: r1, z1, r2, z2, dl
      integer,         intent (out) :: npts
      real (kind=R8), allocatable, intent (out) :: s(:), rr(:), zz(:),   &
     &                                             val(:)
      logical,        allocatable, intent (out) :: ok(:)
      real (kind=R8), intent (in), optional :: power
      logical,        intent (in), optional :: include_guard

      real (kind=R8) :: length, t, pw
      logical        :: guard
      integer        :: n, k, hint, iCv

      pw = 2.0_R8
      if (present(power)) pw = power
      guard = .false.
      if (present(include_guard)) guard = include_guard

      length = sqrt((r2-r1)**2 + (z2-z1)**2)
      if (dl.le.0.0_R8 .or. length.le.0.0_R8) then
        npts = 0
        allocate(s(0), rr(0), zz(0), val(0), ok(0))
        return
      end if

      n    = max(1, ceiling(length/dl))
      npts = n + 1
      allocate(s(npts), rr(npts), zz(npts), val(npts), ok(npts))

      hint = 0
      do k = 1, npts
        t     = real(k-1, R8)/real(n, R8)
        s(k)  = t*length
        rr(k) = r1 + t*(r2-r1)
        zz(k) = z1 + t*(z2-z1)
        call interp_cv_at_rz(geo, mpg, fld, rr(k), zz(k), val(k), iCv,   &
     &                       ok(k), hint=hint, power=pw,                 &
     &                       include_guard=guard)
!       keep the last valid cell as the hint so the search re-homes
!       quickly if the chord leaves and re-enters the mesh
        if (iCv.gt.0) hint = iCv
      end do

      return
      end subroutine chord_profile_cv

!-----------------------------------------------------------------------
      subroutine chord_profile_fc (geo, mpg, face2, r1, z1, r2, z2, dl,  &
     &                             npts, s, rr, zz, val_pol, val_rad, ok,&
     &                             power, include_guard)
!     Sample a face-centered flux field along a chord. face2(nFc,0:1)
!     holds the poloidal component in (:,0) and the radial component in
!     (:,1), as for dv%fna(:,:,is), dv%fhe, dv%fhi, etc. The fluxes are
!     converted to flux densities with divide_by_area (projected face
!     area geo%fcS*geo%fcQalf, the b2mod_mwti target convention), then
!     reconstructed to cell centers with the code's stencil-consistent
!     operator (intcell with intcellP / intcellR) and interpolated. On
!     exit val_pol is the poloidal flux-density profile and
!     val_rad the radial flux-density profile; the sampling
!     (npts, s, rr, zz, ok) is common to both.
      implicit none
      type (geometry), intent (in) :: geo
      type (mapping),  intent (in) :: mpg
      real (kind=R8),  intent (in) :: face2(mpg%nFc,0:1)
      real (kind=R8),  intent (in) :: r1, z1, r2, z2, dl
      integer,         intent (out) :: npts
      real (kind=R8), allocatable, intent (out) :: s(:), rr(:), zz(:),   &
     &                                             val_pol(:), val_rad(:)
      logical,        allocatable, intent (out) :: ok(:)
      real (kind=R8), intent (in), optional :: power
      logical,        intent (in), optional :: include_guard

      real (kind=R8), allocatable :: cpol(:), crad(:)
      real (kind=R8), allocatable :: s2(:), rr2(:), zz2(:)
      logical,        allocatable :: ok2(:)
      real (kind=R8) :: pw
      logical        :: guard
      integer        :: npts2

      pw = 2.0_R8
      if (present(power)) pw = power
      guard = .false.
      if (present(include_guard)) guard = include_guard

      allocate(cpol(mpg%nCv), crad(mpg%nCv))

!     face flux -> cell-centered poloidal/radial flux densities
      call facefield_to_cellcenters(geo, mpg, face2, cpol, crad)

      call chord_profile_cv(geo, mpg, cpol, r1, z1, r2, z2, dl,          &
     &                      npts, s, rr, zz, val_pol, ok,                &
     &                      power=pw, include_guard=guard)
      call chord_profile_cv(geo, mpg, crad, r1, z1, r2, z2, dl,          &
     &                      npts2, s2, rr2, zz2, val_rad, ok2,           &
     &                      power=pw, include_guard=guard)

      deallocate(cpol, crad, s2, rr2, zz2, ok2)

      return
      end subroutine chord_profile_fc

!-----------------------------------------------------------------------
      subroutine facefield_to_cellcenters (geo, mpg, face2, cpol, crad)
!     Convert a face-centered flux field face2(nFc,0:1) (component 0 =
!     poloidal, 1 = radial) into cell-centered poloidal (cpol) and radial
!     (crad) flux densities: divide by the projected face area with
!     divide_by_area (geo%fcS*geo%fcQalf, the b2mod_mwti target
!     convention), then collapse to cell centers with intcell.
      implicit none
      type (geometry), intent (in)  :: geo
      type (mapping),  intent (in)  :: mpg
      real (kind=R8),  intent (in)  :: face2(mpg%nFc,0:1)
      real (kind=R8),  intent (out) :: cpol(mpg%nCv), crad(mpg%nCv)
      real (kind=R8), allocatable :: dens2(:,:)
      external :: intcell
      allocate(dens2(mpg%nFc,0:1))
      call divide_by_area(mpg%nFc, geo, face2, dens2)
      call intcell(mpg%nFc, mpg%nCv, mpg, mpg%intcellP, dens2(1,0), cpol)
      call intcell(mpg%nFc, mpg%nCv, mpg, mpg%intcellR, dens2(1,1), crad)
      deallocate(dens2)
      return
      end subroutine facefield_to_cellcenters

!-----------------------------------------------------------------------
      subroutine chord_npoints (r1, z1, r2, z2, dl, np)
!     Number of sample points on one chord (endpoints inclusive), for a
!     spacing of about dl. A degenerate/zero chord returns 1.
      implicit none
      real (kind=R8), intent (in)  :: r1, z1, r2, z2, dl
      integer,        intent (out) :: np
      real (kind=R8) :: length
      length = sqrt((r2-r1)**2 + (z2-z1)**2)
      if (dl.le.0.0_R8 .or. length.le.0.0_R8) then
        np = 1
      else
        np = max(1, ceiling(length/dl)) + 1
      end if
      return
      end subroutine chord_npoints

!-----------------------------------------------------------------------
      subroutine chords_total_points (rzchords, nchords, dl, npt)
!     Total number of sample points over all nchords chords.
      implicit none
      integer,        intent (in)  :: nchords
      real (kind=R8), intent (in)  :: rzchords(2,2,nchords)
      real (kind=R8), intent (in)  :: dl
      integer,        intent (out) :: npt
      integer :: k, np
      npt = 0
      do k = 1, nchords
        call chord_npoints(rzchords(1,1,k), rzchords(2,1,k),            &
     &                     rzchords(1,2,k), rzchords(2,2,k), dl, np)
        npt = npt + np
      end do
      return
      end subroutine chords_total_points

!-----------------------------------------------------------------------
      subroutine chords_geometry (rzchords, nchords, dl, npt,           &
     &                            crx, cry, ds, ichord)
!     Fill the concatenated sample-point coordinates for all chords:
!     crx/cry are the (R,Z) of each point, ds the arc length from the
!     start of its chord, ichord the chord index (1..nchords). The point
!     ordering matches chords_total_points / sample_cv_at_points.
      implicit none
      integer,        intent (in)  :: nchords, npt
      real (kind=R8), intent (in)  :: rzchords(2,2,nchords), dl
      real (kind=R8), intent (out) :: crx(npt), cry(npt), ds(npt)
      integer,        intent (out) :: ichord(npt)
      integer :: k, i, np, off
      real (kind=R8) :: r1, z1, r2, z2, length, t
      off = 0
      do k = 1, nchords
        r1 = rzchords(1,1,k)
        z1 = rzchords(2,1,k)
        r2 = rzchords(1,2,k)
        z2 = rzchords(2,2,k)
        call chord_npoints(r1, z1, r2, z2, dl, np)
        length = sqrt((r2-r1)**2 + (z2-z1)**2)
        do i = 1, np
          if (np.gt.1) then
            t = real(i-1, R8)/real(np-1, R8)
          else
            t = 0.0_R8
          end if
          crx(off+i)    = r1 + t*(r2-r1)
          cry(off+i)    = z1 + t*(z2-z1)
          ds(off+i)     = t*length
          ichord(off+i) = k
        end do
        off = off + np
      end do
      return
      end subroutine chords_geometry

!-----------------------------------------------------------------------
      subroutine sample_cv_at_points (geo, mpg, fld, crx, cry, ichord,  &
     &                                npt, val, include_guard)
!     Interpolate the cell-centered field fld at the npt sample points
!     (crx,cry). The cell hint is chained along a chord and reset at each
!     chord boundary (ichord change). Points outside the mesh get 0.
      implicit none
      type (geometry), intent (in)  :: geo
      type (mapping),  intent (in)  :: mpg
      real (kind=R8),  intent (in)  :: fld(mpg%nCv)
      integer,         intent (in)  :: npt
      real (kind=R8),  intent (in)  :: crx(npt), cry(npt)
      integer,         intent (in)  :: ichord(npt)
      real (kind=R8),  intent (out) :: val(npt)
      logical, intent (in), optional :: include_guard
      logical :: guard, ok
      integer :: i, iCv, hint
      guard = .false.
      if (present(include_guard)) guard = include_guard
      hint = 0
      do i = 1, npt
        if (i.gt.1) then
          if (ichord(i).ne.ichord(i-1)) hint = 0
        end if
        call interp_cv_at_rz(geo, mpg, fld, crx(i), cry(i), val(i),     &
     &                       iCv, ok, hint=hint, include_guard=guard)
        if (iCv.gt.0) hint = iCv
        if (.not.ok) val(i) = 0.0_R8
      end do
      return
      end subroutine sample_cv_at_points

!-----------------------------------------------------------------------
      subroutine mask_inside_cv (geo, mpg, crx, cry, ichord, npt, inside)
!     Flag each of the npt sample points that lies inside an internal
!     control volume (inside=.true.) or outside the B2.5 grid
!     (inside=.false.). The cell hint is chained along a chord and reset
!     at each chord boundary (ichord change), exactly as in
!     sample_cv_at_points, so the mask is consistent with what the
!     sampler will later find at those same points.
      implicit none
      type (geometry), intent (in)  :: geo
      type (mapping),  intent (in)  :: mpg
      integer,         intent (in)  :: npt
      real (kind=R8),  intent (in)  :: crx(npt), cry(npt)
      integer,         intent (in)  :: ichord(npt)
      logical,         intent (out) :: inside(npt)
      integer :: i, iCv, hint
      hint = 0
      do i = 1, npt
        if (i.gt.1) then
          if (ichord(i).ne.ichord(i-1)) hint = 0
        end if
        call locate_cv(geo, mpg, crx(i), cry(i), iCv, hint)
        inside(i) = (iCv.gt.0)
        if (iCv.gt.0) hint = iCv
      end do
      return
      end subroutine mask_inside_cv

!-----------------------------------------------------------------------
      subroutine compact_chords (crx0, cry0, ds0, ichord0, npt0,         &
     &                           nchords, inside, crx, cry, ds, ichord,   &
     &                           npt, chordid, nchord, rstart, zstart,    &
     &                           rend, zend)
!     Keep only the sample points flagged inside(:) and drop any chord
!     with no inside point. The kept points are packed into
!     crx/cry/ds/ichord (ichord renumbered 1..nchord over the surviving
!     chords, in their original order); npt is the kept-point count and
!     nchord the surviving-chord count. chordid(m) is the original chord
!     index (1..nchords) of surviving chord m. rstart/zstart and
!     rend/zend are the (R,Z) of the first and last kept point of each
!     surviving chord (its in-domain extent). ds keeps its original
!     arc-length-from-chord-start value, so a chord that leaves and
!     re-enters the domain shows a gap in ds. All output arrays are
!     allocated here.
      implicit none
      integer,        intent (in)  :: npt0, nchords
      real (kind=R8), intent (in)  :: crx0(npt0), cry0(npt0), ds0(npt0)
      integer,        intent (in)  :: ichord0(npt0)
      logical,        intent (in)  :: inside(npt0)
      real (kind=R8), allocatable, intent (out) :: crx(:), cry(:), ds(:)
      integer,        allocatable, intent (out) :: ichord(:)
      integer,        intent (out) :: npt
      integer,        allocatable, intent (out) :: chordid(:)
      integer,        intent (out) :: nchord
      real (kind=R8), allocatable, intent (out) :: rstart(:), zstart(:), &
     &                                             rend(:), zend(:)
      integer :: i, k, m, j
      integer, allocatable :: newid(:)
      logical, allocatable :: seen(:)

!     which chords survive, and their new contiguous numbering
      allocate(newid(nchords))
      newid = 0
      do i = 1, npt0
        if (inside(i)) newid(ichord0(i)) = 1
      end do
      nchord = 0
      do k = 1, nchords
        if (newid(k).gt.0) then
          nchord = nchord + 1
          newid(k) = nchord
        end if
      end do

!     number of kept points
      npt = 0
      do i = 1, npt0
        if (inside(i)) npt = npt + 1
      end do

      allocate(crx(npt), cry(npt), ds(npt), ichord(npt))
      allocate(chordid(nchord))
      allocate(rstart(nchord), zstart(nchord), rend(nchord), zend(nchord))
      allocate(seen(nchord))
      seen = .false.

      do k = 1, nchords
        if (newid(k).gt.0) chordid(newid(k)) = k
      end do

      j = 0
      do i = 1, npt0
        if (.not.inside(i)) cycle
        j         = j + 1
        crx(j)    = crx0(i)
        cry(j)    = cry0(i)
        ds(j)     = ds0(i)
        m         = newid(ichord0(i))
        ichord(j) = m
        if (.not.seen(m)) then
          rstart(m) = crx0(i)
          zstart(m) = cry0(i)
          seen(m)   = .true.
        end if
        rend(m) = crx0(i)
        zend(m) = cry0(i)
      end do

      deallocate(newid, seen)
      return
      end subroutine compact_chords

!-----------------------------------------------------------------------
      subroutine point_in_cell (geo, mpg, iCv, r, z, inside)
!     Return inside=.true. if (r,z) lies in control volume iCv, tested by
!     fanning the cell into triangles (cell center, face vertex 1, face
!     vertex 2) over the cell's faces. Independent of vertex ordering;
!     exact for convex cells.
      implicit none
      type (geometry), intent (in)  :: geo
      type (mapping),  intent (in)  :: mpg
      integer,         intent (in)  :: iCv
      real (kind=R8),  intent (in)  :: r, z
      logical,         intent (out) :: inside
      integer :: i, iFc, p0, nf, v1, v2
      real (kind=R8) :: cx, cy

      inside = .false.
      cx = geo%cvX(iCv)
      cy = geo%cvY(iCv)
      p0 = mpg%cvFcP(iCv,1)
      nf = mpg%cvFcP(iCv,2)
      do i = 0, nf-1
        iFc = mpg%cvFc(p0+i)
        v1  = mpg%fcVx(iFc,1)
        v2  = mpg%fcVx(iFc,2)
        call point_in_triangle(r, z, cx, cy,                            &
     &       geo%vxX(v1), geo%vxY(v1), geo%vxX(v2), geo%vxY(v2), inside)
        if (inside) return
      end do

      return
      end subroutine point_in_cell

!-----------------------------------------------------------------------
      subroutine point_in_triangle (px, py, ax, ay, bx, by, cx, cy,     &
     &                              inside)
!     Standard half-plane test: inside=.true. if the point is in triangle
!     (a,b,c), i.e. the three edge signs do not straddle zero (points on
!     an edge count as inside).
      implicit none
      real (kind=R8), intent (in)  :: px, py, ax, ay, bx, by, cx, cy
      logical,        intent (out) :: inside
      real (kind=R8) :: d1, d2, d3
      logical :: has_neg, has_pos

      call edge_sign(px, py, ax, ay, bx, by, d1)
      call edge_sign(px, py, bx, by, cx, cy, d2)
      call edge_sign(px, py, cx, cy, ax, ay, d3)
      has_neg = (d1.lt.0.0_R8) .or. (d2.lt.0.0_R8) .or. (d3.lt.0.0_R8)
      has_pos = (d1.gt.0.0_R8) .or. (d2.gt.0.0_R8) .or. (d3.gt.0.0_R8)
      inside = .not. (has_neg .and. has_pos)

      return
      end subroutine point_in_triangle

!-----------------------------------------------------------------------
      subroutine edge_sign (px, py, ax, ay, bx, by, sgn)
!     Twice the signed area of triangle (a,b,point); its sign tells which
!     side of the directed edge a->b the point lies on.
      implicit none
      real (kind=R8), intent (in)  :: px, py, ax, ay, bx, by
      real (kind=R8), intent (out) :: sgn
      sgn = (px - bx)*(ay - by) - (ax - bx)*(py - by)
      return
      end subroutine edge_sign

      end module b2mod_chord_interp
