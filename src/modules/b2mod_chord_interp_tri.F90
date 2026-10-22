      module b2mod_chord_interp_tri
!=======================================================================
!     Interpolation of EIRENE triangular-grid quantities onto arbitrary
!     (R,Z) points and along chords.
!
!     ncltal(itri) maps triangle itri (1..NTRII) to its cell, so a
!     per-triangle field is built with build_tri_field and then
!     interpolated. The triangle geometry is the LEVGEO=4 mesh from
!     eirmod_ctrig: vertices XTRIAN/YTRIAN, triangle->vertex map NECKE,
!     triangle->neighbour map NCHBAR (0 on a boundary side). Location is
!     a point-in-triangle test walked through NCHBAR from a hint, with a
!     global scan fallback; interpolation is inverse-distance weighting
!     over the containing triangle and its edge neighbours (triangle
!     centroids).
!=======================================================================
#ifdef B25_EIRENE
      use b2mod_types,  only : R8
      use eirmod_ctrig, only : XTRIAN, YTRIAN, NECKE, NCHBAR, NTRII
      use eirmod_cgeom, only : NCLTAL
      implicit none
      private

      public :: build_tri_field, sample_tri_at_points, mask_inside_tri

      real (kind=R8), parameter :: m_to_cm = 100.0_R8
      real (kind=R8), parameter :: dist_eps = 1.0e-10_R8   ! cm

      contains

!-----------------------------------------------------------------------
      subroutine build_tri_field (fcell, ncell, ftri)
!     Map a per-EIRENE-cell field fcell(ncell) to a per-triangle field
!     ftri(NTRII) using the triangle->cell map ncltal.
      implicit none
      integer,        intent (in)  :: ncell
      real (kind=R8), intent (in)  :: fcell(ncell)
      real (kind=R8), intent (out) :: ftri(NTRII)
      integer :: itri, ic
      do itri = 1, NTRII
        ic = NCLTAL(itri)
        if (ic.ge.1 .and. ic.le.ncell) then
          ftri(itri) = fcell(ic)
        else
          ftri(itri) = 0.0_R8
        end if
      end do
      return
      end subroutine build_tri_field

!-----------------------------------------------------------------------
      subroutine sample_tri_at_points (ftri, crx, cry, ichord, npt, val)
!     Interpolate the per-triangle field ftri at the npt chord points
!     (crx,cry) given in metres. The triangle hint is chained along a
!     chord and reset at each chord boundary. Points outside the mesh
!     get 0.
      implicit none
      real (kind=R8), intent (in)  :: ftri(NTRII)
      integer,        intent (in)  :: npt
      real (kind=R8), intent (in)  :: crx(npt), cry(npt)
      integer,        intent (in)  :: ichord(npt)
      real (kind=R8), intent (out) :: val(npt)
      integer :: i, itri, hint
      logical :: ok
      hint = 0
      do i = 1, npt
        if (i.gt.1) then
          if (ichord(i).ne.ichord(i-1)) hint = 0
        end if
        call interp_tri(ftri, crx(i)*m_to_cm, cry(i)*m_to_cm,           &
     &                  val(i), itri, ok, hint)
        if (itri.gt.0) hint = itri
        if (.not.ok) val(i) = 0.0_R8
      end do
      return
      end subroutine sample_tri_at_points

!-----------------------------------------------------------------------
      subroutine mask_inside_tri (crx, cry, ichord, npt, inside)
!     Flag each of the npt sample points (given in metres) that lies
!     inside a triangle (inside=.true.) or outside the triangular mesh
!     (inside=.false.). The triangle hint is chained along a chord and
!     reset at each chord boundary, consistent with sample_tri_at_points.
      implicit none
      integer,        intent (in)  :: npt
      real (kind=R8), intent (in)  :: crx(npt), cry(npt)
      integer,        intent (in)  :: ichord(npt)
      logical,        intent (out) :: inside(npt)
      integer :: i, itri, hint
      hint = 0
      do i = 1, npt
        if (i.gt.1) then
          if (ichord(i).ne.ichord(i-1)) hint = 0
        end if
        call locate_tri(crx(i)*m_to_cm, cry(i)*m_to_cm, hint, itri)
        inside(i) = (itri.gt.0)
        if (itri.gt.0) hint = itri
      end do
      return
      end subroutine mask_inside_tri

!-----------------------------------------------------------------------
      subroutine interp_tri (ftri, xc, yc, val, itri, ok, hint)
!     Interpolate ftri at (xc,yc) [cm] by inverse-distance weighting over
!     the containing triangle and its edge neighbours (triangle
!     centroids). itri returns the containing triangle (0 if outside).
      implicit none
      real (kind=R8), intent (in)  :: ftri(NTRII)
      real (kind=R8), intent (in)  :: xc, yc
      real (kind=R8), intent (out) :: val
      integer,        intent (out) :: itri
      logical,        intent (out) :: ok
      integer,        intent (in)  :: hint
      real (kind=R8) :: wsum, vsum, dist, w, gx, gy
      integer :: is, c

      val = 0.0_R8
      ok  = .false.
      itri = 0

      call locate_tri(xc, yc, hint, itri)
      if (itri.eq.0) return

      wsum = 0.0_R8
      vsum = 0.0_R8
!     containing triangle
      call tri_centroid(itri, gx, gy)
      dist = sqrt((gx-xc)**2 + (gy-yc)**2)
      if (dist.lt.dist_eps) then
        val = ftri(itri)
        ok  = .true.
        return
      end if
      w    = 1.0_R8/(dist*dist)
      wsum = wsum + w
      vsum = vsum + w*ftri(itri)
!     edge neighbours
      do is = 1, 3
        c = NCHBAR(is, itri)
        if (c.lt.1 .or. c.gt.NTRII) cycle
        call tri_centroid(c, gx, gy)
        dist = sqrt((gx-xc)**2 + (gy-yc)**2)
        if (dist.lt.dist_eps) then
          val = ftri(c)
          ok  = .true.
          return
        end if
        w    = 1.0_R8/(dist*dist)
        wsum = wsum + w
        vsum = vsum + w*ftri(c)
      end do

      if (wsum.gt.0.0_R8) then
        val = vsum/wsum
        ok  = .true.
      end if
      return
      end subroutine interp_tri

!-----------------------------------------------------------------------
      subroutine locate_tri (xc, yc, hint, itri)
!     Find the triangle containing (xc,yc) [cm]: try the hint and its
!     three edge neighbours, then scan all triangles.
      implicit none
      real (kind=R8), intent (in)  :: xc, yc
      integer,        intent (in)  :: hint
      integer,        intent (out) :: itri
      integer :: is, c, k
      logical :: inside

      itri = 0
      if (hint.ge.1 .and. hint.le.NTRII) then
        call point_in_tri(hint, xc, yc, inside)
        if (inside) then
          itri = hint
          return
        end if
        do is = 1, 3
          c = NCHBAR(is, hint)
          if (c.ge.1 .and. c.le.NTRII) then
            call point_in_tri(c, xc, yc, inside)
            if (inside) then
              itri = c
              return
            end if
          end if
        end do
      end if
      do k = 1, NTRII
        call point_in_tri(k, xc, yc, inside)
        if (inside) then
          itri = k
          return
        end if
      end do
      return
      end subroutine locate_tri

!-----------------------------------------------------------------------
      subroutine tri_centroid (itri, gx, gy)
      implicit none
      integer,        intent (in)  :: itri
      real (kind=R8), intent (out) :: gx, gy
      gx = (XTRIAN(NECKE(1,itri)) + XTRIAN(NECKE(2,itri)) +             &
     &      XTRIAN(NECKE(3,itri))) / 3.0_R8
      gy = (YTRIAN(NECKE(1,itri)) + YTRIAN(NECKE(2,itri)) +             &
     &      YTRIAN(NECKE(3,itri))) / 3.0_R8
      return
      end subroutine tri_centroid

!-----------------------------------------------------------------------
      subroutine point_in_tri (itri, px, py, inside)
!     Half-plane test for the triangle itri (vertices from NECKE).
      implicit none
      integer,        intent (in)  :: itri
      real (kind=R8), intent (in)  :: px, py
      logical,        intent (out) :: inside
      real (kind=R8) :: ax, ay, bx, by, cx, cy, d1, d2, d3
      logical :: has_neg, has_pos
      ax = XTRIAN(NECKE(1,itri)); ay = YTRIAN(NECKE(1,itri))
      bx = XTRIAN(NECKE(2,itri)); by = YTRIAN(NECKE(2,itri))
      cx = XTRIAN(NECKE(3,itri)); cy = YTRIAN(NECKE(3,itri))
      d1 = (px-bx)*(ay-by) - (ax-bx)*(py-by)
      d2 = (px-cx)*(by-cy) - (bx-cx)*(py-cy)
      d3 = (px-ax)*(cy-ay) - (cx-ax)*(py-ay)
      has_neg = (d1.lt.0.0_R8) .or. (d2.lt.0.0_R8) .or. (d3.lt.0.0_R8)
      has_pos = (d1.gt.0.0_R8) .or. (d2.gt.0.0_R8) .or. (d3.gt.0.0_R8)
      inside = .not. (has_neg .and. has_pos)
      return
      end subroutine point_in_tri
#endif
      end module b2mod_chord_interp_tri
