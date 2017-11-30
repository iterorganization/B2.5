!!-----------------------------------------------------------------------------
!! DOCUMENTATION (doxygen 1.8.8):
!>      @section b2mod_ch_desc   Description
!!      Module providing routines for point distance measure and cell geometry
!!      type recognition.
!!
!!      @subsection b2mod_ch_pv    Parameters/variables
!!      @note   see also routine b2cdcv
!!
!!      @param  crx, cry - For (ix,iy) in (-1:nx,-1:ny), crx(ix,iy,0:3) specify
!!              the horizontal coordinates and cry(ix,iy,0:3) specify the
!!              vertical coordinates (in the physical plane) of the four
!!              corners of the (ix,iy) cell, mapped as follows:
!!              (,,0): lower left corner, (,,1): lower right corner,
!!              (,,2): upper left corner, (,,3): upper right corner.
!!              One expects the following identities to hold away from the
!!              boundaries:
!!              crx(ix,iy,0)=crx(ix-1,iy,1)=crx(ix,iy-1,2)=crx(ix-1,iy-1,3),
!!              cry(ix,iy,0)=cry(ix-1,iy,1)=cry(ix,iy-1,2)=cry(ix-1,iy-1,3),
!!              However, there may be cuts in the geometry, in which case these
!!              identities need not hold.
!!      @param  x1 - First coordinate of the first point
!!      @param  y1 - Second coordinate of the first point
!!      @param  x2 - First coordinate of the second point
!!      @param  y2 - Second coordinate of the second point
!!
!!-----------------------------------------------------------------------------
module b2mod_cellhelper
    use b2mod_types

    implicit none

    !! Direction numbering
    integer, parameter :: LEFT = 0, BOTTOM = 1, RIGHT = 2, TOP = 3
    integer, parameter :: NODIRECTION = -1

    !! Directed volumes numbering
    integer, parameter :: TO_SELF = 0
    integer, parameter :: TO_LEFT = 1
    integer, parameter :: TO_BOTTOM = 2
    integer, parameter :: TO_RIGHT = 3
    integer, parameter :: TO_TOP = 4

    !! Cell offset when moving to neighbour cell in on of the above directions
    integer, parameter :: dxDir(0:3) = (/-1,  0, 1,  0/)
    integer, parameter :: dyDir(0:3) = (/ 0, -1, 0,  1/)

    !! Cell corner vertex indices
    integer, parameter :: VX_LOWER = 0, VX_UPPER = 2, VX_LEFT = 0, VX_RIGHT = 1
    integer, parameter :: VX_LOWERLEFT = VX_LOWER + VX_LEFT    !! 0
    integer, parameter :: VX_LOWERRIGHT = VX_LOWER + VX_RIGHT  !! 1
    integer, parameter :: VX_UPPERLEFT = VX_UPPER + VX_LEFT    !! 2
    integer, parameter :: VX_UPPERRIGHT = VX_UPPER + VX_RIGHT  !! 3
    integer, parameter :: VX_UNDEFINED = -1

    integer, parameter :: B2_GRID_UNDEFINED = -2

    !! Circling around vertices, following connectivity in the B2 grid
    integer, parameter :: CLOCKWISE = 1
    integer, parameter :: COUNTERCLOCKWISE = 2

    !! Directions for clockwise circle around corner vertex
    !! To circle around cell corner vertex iCorner in clockwise direction,
    !! do istep = 1, 4; VXCIRCLE_STEPDIR(istep, iCorner, CLOCKWISE)
    !! To circle around cell corner vertex iCorner in counterclockwise direction,
    !! do istep = 1, 4; VXCIRCLE_CLOCKWISE(istep, iCorner, COUNTERCLOCKWISE)
    integer, parameter :: VXCIRCLE_STEPDIR(1:4, 0:3, 2) = reshape(  &
      & (/ BOTTOM, LEFT, TOP, RIGHT, &  ! Vertex 0/Lower left, clockwise
      &    RIGHT, BOTTOM, LEFT, TOP, &  !
      &    LEFT, TOP, RIGHT, BOTTOM, &
      &    TOP, RIGHT, BOTTOM, LEFT, &
      &    LEFT, BOTTOM, RIGHT, TOP, &
      &    BOTTOM, RIGHT, TOP, LEFT, &
      &    TOP, LEFT, BOTTOM, RIGHT, &
      &    RIGHT, TOP, LEFT, BOTTOM /), &
      & (/4, 4, 2/) )

    !! Lookup table stating what corner vertex index a vertex has
    !! in the neighbour cell when stepping in a given direction.
    !! Usage: corner iCorner=VX_..., step direction iDir = LEFT, BOTTOM, ...
    !! newICorner = VXCORNER_NEXTINDEX(iDir, iCorner)
    integer, parameter :: VXCORNER_NEXTINDEX(0:3, 0:3) = reshape(   &
      & (/ VX_LOWERRIGHT, VX_UPPERLEFT, VX_UNDEFINED, VX_UNDEFINED, &
      &    VX_UNDEFINED, VX_UPPERRIGHT, VX_LOWERLEFT, VX_UNDEFINED, &
      &    VX_UPPERRIGHT, VX_UNDEFINED, VX_UNDEFINED, VX_LOWERLEFT, &
      &    VX_UNDEFINED, VX_UNDEFINED, VX_UPPERLEFT, VX_LOWERRIGHT /),  &
      & (/4, 4/) )

    !! Start and end vertices of faces
    integer, parameter :: VX_START = 0
    integer, parameter :: VX_END = 1

    integer, parameter :: FACE_VERTEX(0:1, 0:3) = reshape( &
      & (/ VX_LOWERLEFT, VX_UPPERLEFT, &
      &    VX_LOWERLEFT, VX_LOWERRIGHT, &
      &    VX_LOWERRIGHT, VX_UPPERRIGHT, &
      &    VX_UPPERLEFT, VX_UPPERRIGHT /), &
      & (/ 2, 4 /) )

    !! Value marking a coordinate position to be invalid
    real(R8), parameter :: INVALID_POSITION = 1.0e3_R8

    !! Value marking a coordinate position to be invalid
    real(R8), parameter :: INVALID_DOUBLE = 9.99999e99_R8


    !! Distance between two points at which the points are declared to be equal
    real(R8), parameter, private :: geom_match_dist = 1.0e-9_R8

    !! Cell geometry types for cellGeoType and isTriangle
    integer, parameter :: CGEO_BROKEN = 0
    integer, parameter :: CGEO_QUAD = 1
    integer, parameter :: CGEO_TRIA_NOLEFT = 2
    integer, parameter :: CGEO_TRIA_NOBOT = 3
    integer, parameter :: CGEO_TRIA_NORIGHT = 4
    integer, parameter :: CGEO_TRIA_NOTOP = 5

contains

    !> Compute distance between points
    !! @param[in]   x1 - First coordinate of the first point
    !! @param[in]   y1 - Second coordinate of the first point
    !! @param[in]   x2 - First coordinate of the second point
    !! @param[in]   y2 - Second coordinate of the second point
    real(R8) function points_dist( x1, y1, x2, y2 )
        !! Internal variables
        real(R8), intent(in) :: x1, y1, x2, y2

        points_dist = sqrt( (x1-x2)**2 + (y1-y2)**2 )
    end function points_dist

    !> Check whether to points coincide
    !! It tests whether their distance is closer than geom_match_dist.
    !! @param[in]   x1 - First coordinate of the first point
    !! @param[in]   y1 - Second coordinate of the first point
    !! @param[in]   x2 - First coordinate of the second point
    !! @param[in]   y2 - Second coordinate of the second point
    logical function points_match( x1, y1, x2, y2 )
        !! Internal variables
        real(R8), intent(in) :: x1, y1, x2, y2

        points_match = points_dist( x1, y1, x2, y2 ).lt.geom_match_dist
    end function points_match


    !> Determine the geometry type of a cell
    !! @param[in]   crx - horizontal coordinate of the four corners of
    !!              the (ix,iy) cell
    !! @param[in]   cry - vertical coordinate of the four corners of
    !!              the (ix,iy) cell
    integer function cellGeoType( crx, cry )
        real(R8), dimension(0:3), intent(in) :: crx, cry

        !! Internal variable
        logical :: leftFace, botFace, rightFace, topFace
        integer :: fcount

        !! check which faces are present
        leftFace =  .not. points_match(crx(0), cry(0), crx(2), cry(2))
        botFace =   .not. points_match(crx(0), cry(0), crx(1), cry(1))
        rightFace = .not. points_match(crx(1), cry(1), crx(3), cry(3))
        topFace =   .not. points_match(crx(2), cry(2), crx(3), cry(3))

        fcount = 0
        if( leftFace ) fcount = fcount + 1
        if( botFace ) fcount = fcount + 1
        if( rightFace ) fcount = fcount + 1
        if( topFace ) fcount = fcount + 1

        select case (fcount)
        case (0:2)
            cellGeoType = CGEO_BROKEN
        case (3)
            if( .not. leftFace ) cellGeoType = CGEO_TRIA_NOLEFT
            if( .not. botFace ) cellGeoType = CGEO_TRIA_NOBOT
            if( .not. rightFace ) cellGeoType = CGEO_TRIA_NORIGHT
            if( .not. topFace ) cellGeoType = CGEO_TRIA_NOTOP
        case (4)
            cellGeoType = CGEO_QUAD
        end select

    end function cellGeoType

    !> Check whether a cell geometry type returned by cellGeoType
    !! indicates a triangular cell
    !! @param[in]   cellGeoType - Cell type ID (e.g. CGEO_TRIA_NOLEFT)
    logical function isTriangleCell( cellGeoType )
        !! Internal variables
        integer, intent(in) :: cellGeoType

        isTriangleCell = ( cellGeoType == CGEO_TRIA_NOLEFT ) .or.   &
            &   ( cellGeoType == CGEO_TRIA_NOBOT ) .or.             &
            &   ( cellGeoType == CGEO_TRIA_NORIGHT ) .or.           &
            &   ( cellGeoType == CGEO_TRIA_NOTOP )

    end function isTriangleCell

end module b2mod_cellhelper

!!!Local Variables:
!!! mode: f90
!!! End:
