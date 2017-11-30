!!-----------------------------------------------------------------------------
!! DOCUMENTATION:
!>      @section b2mod_gmap_desc  Description
!!      Module providing a mechanism to map from the B2 data structure to the
!!      IDS and CPO data structures, consisting of a routine to set up the map
!!      (b2CreateMap), a data structure to hold the map information (B2GridMap)
!!      and some service routines to handle this data structure.
!!
!!      @subsection b2mod_gmap_pv  Parameters/variables
!!      @param  geom_match_dist - Distance between two points at which the
!!              points are declared to be equal
!!      @param  ALIGNX, ALIGNY - Alignment index (for example in B2 flux arrays)
!!      @param  MAX_SPECIAL_VERTICES - Maximum number of special vertices
!!              expected in the grid
!!
!!-----------------------------------------------------------------------------

module b2mod_grid_mapping

    use b2mod_types , B2_R8 => R8, B2_R4 => R4
    use helper
    use logging , only: logmsg, LOGDEBUG
    use b2mod_connectivity , REMOVED_B2_R8 => R8
    use carre_constants
    use b2mod_cellhelper

    implicit none

    !! Distance between two points at which the points are declared to be equal
    real(R8), private, parameter :: geom_match_dist = 1.0e-6_R8

    !! Alignment index (for example in B2 flux arrays)
    integer, parameter :: ALIGNX = 1
    integer, parameter :: ALIGNY = 0

    !! Maximum number of special vertices expected in the grid
    integer, parameter :: MAX_SPECIAL_VERTICES = 10

    !> Data structure holding an intermediate grid description to be
    !! transferred into a CPO or IDS
    type B2GridMap
        integer :: ncv  !< Number of all cells in the domain (2D objects)
        integer :: nfcx !< Number of x-aligned faces/edges in the domain
                        !< (1D objects)
        integer :: nfcy !< Number of y-aligned faces/edges in the domain
                        !< (1D objects)
        integer :: nvx  !< Number of all vertices/nodes in the domain
                        !<(0D objects)
        integer :: b2nx
        integer :: b2ny

        integer, dimension(:), allocatable :: mapCvix !< Array of horizontal
            !< (x-aligned) cell position in computational space
        integer, dimension(:), allocatable :: mapCviy !< Array of vertical
            !< (y-aligned) cell position in computational space
        integer, dimension(:), allocatable :: mapFcix !< Array of horizontal
            !< (x-aligned) face/edge position in computational space
        integer, dimension(:), allocatable :: mapFciy !< Array of vertical
            !< (y-aligned) face/edge position in computational space
        integer, dimension(:), allocatable :: mapFcIFace
        integer, dimension(:), allocatable :: mapVxix !< Array of horizontal
            !< (x-aligned) vertex/node position in computational space
        integer, dimension(:), allocatable :: mapVxiy !< Array of vertical
            !< (y-aligned) vertex/node position in computational space
        integer, dimension(:), allocatable :: mapVxIVx
        integer, dimension(:,:), allocatable :: mapCvI
        integer, dimension(:,:,:), allocatable :: mapFcI !< 3D array of faces
            !< composing the quadliateral in the list.
            !< gmap%mapFcI( ix, iy, POSITION ), where POSITION is LEFT, BOTTOM,
            !< RIGHT OR TOP
        integer, dimension(:,:,:), allocatable ::mapVxI

        !! Special vertices (x-points)
        integer :: nsv !< Number of special vertices
        !! svix, sviy : positions of special vertices in B2 data structure
        !! (lower left corner of svix, sviy cell)
        integer, dimension(:), allocatable :: svix !< Array of horizontal
            !< (x-aligned) positions of special vertices in B2 data structure
        integer, dimension(:), allocatable :: sviy !< Array of vertical
            !< (y-aligned) positions of special vertices in B2 data structure
        integer, dimension(:), allocatable :: svi !< Array of indices of
            !< special vertices in the CPO/IDS data structure

        integer, dimension(:,:), allocatable :: mapCvixVx !< Correspondence
            !< between vertices and cells
        integer, dimension(:,:), allocatable :: mapCviyVx !< Correspondence
            !< between vertices and cells
    end type B2GridMap

    logical, save :: mapInitialized = .false.
    logical, save :: map1Initialized = .false.

    private :: R8

!!$        !! Mapping arrays:
!!$        !! 1d lists -> 2d B2 data structure ( i -> (ix, iy) )
!!$        ! b2cv( mapCvix(i), mapCviy(i) ) = cpocv(i)
!!$        ! b2fc( mapFcix(i), mapFciy(i), mapFciFace(i) ) = cpofc(i)
!!$        !! for "normal" faces, mapFcIFace will be LEFT or BOTTOM

!!$        ! b2vx( mapVxix(i), mapVxiy(i), mapVxIVx(i) ) = cpovx(i)
!!$        !! for "normal" vertices, mapVxIVx(i) will be 1 (lower left vertex)

!!$        !! 2d B2 data structure -> 1d CPO lists ( (ix, iy) -> i )
!!$        ! cpocv( mapCvI(ix, iy) ) = b2cv(ix, iy)
!!$        ! cpofc( mapFcI(ix, iy, iFace) ) = b2fc(ix, iy, iFace)
!!$        ! cpovx( mapVxI(ix, iy, iVertex) ) = b2vx(ix, iy, iVertex)

contains

    !! service routines for B2GridData

    !> Set B2GridMap type, intended to be filled with grid geometry information
    !! using b2CreateMap subroutine
    subroutine allocateB2GridMap( gd, nx, ny, ncv, nfcx, nfcy, nvx )
        type(B2GridMap), intent(inout) :: gd    !< The grid mapping as computed
            !< by b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS
        integer, intent(in) ::  nx  !< Specifies the number of interior cells
                                    !< along the first coordinate
        integer, intent(in) ::  ny  !< Specifies the number of interior cells
                                    !< along the second coordinate
        integer, intent(in) ::  ncv !< Number of all cells (2D objects)
        integer, intent(in) ::  nfcx    !< Number of x-aligned faces/edges
                                        !< (1D objects)
        integer, intent(in) ::  nfcy    !< Number of y-aligned faces/edges
                                        !< (1D objects)
        integer, intent(in) ::  nvx     !< Number of all vertices/nodes
                                        !< (0D objects)

        gd%ncv = ncv
        gd%nfcx = nfcx
        gd%nfcy = nfcy
        gd%nvx = nvx

        gd%b2nx = nx
        gd%b2ny = ny

        allocate( gd%mapCvI(-1:nx+1, -1:ny+1) )
        gd%mapCvI = GRID_UNDEFINED
        allocate( gd%mapFcI(-1:nx+1, -1:ny+1, 0:3) )
        gd%mapFcI = GRID_UNDEFINED
        allocate( gd%mapVxI(-1:nx+1, -1:ny+1, 0:3) )
        gd%mapVxI = GRID_UNDEFINED

        allocate( gd%mapCvix(ncv), gd%mapCviy(ncv) )
        gd%mapCvix = GRID_UNDEFINED
        gd%mapCviy = GRID_UNDEFINED
        allocate( gd%mapFcix(nfcx+nfcy), gd%mapFciy(nfcx+nfcy),     &
            &   gd%mapFcIFace(nfcx+nfcy) )
        gd%mapFcix = GRID_UNDEFINED
        gd%mapFciy = GRID_UNDEFINED
        gd%mapFcIFace = GRID_UNDEFINED
        allocate( gd%mapVxix(nvx), gd%mapVxiy(nvx), gd%mapVxIVx(nvx) )
        gd%mapVxix = GRID_UNDEFINED
        gd%mapVxiy = GRID_UNDEFINED
        gd%mapVxIVx = GRID_UNDEFINED

        gd%nsv = 0
        allocate( gd%svix(MAX_SPECIAL_VERTICES),    &
            &   gd%sviy(MAX_SPECIAL_VERTICES), gd%svi(MAX_SPECIAL_VERTICES) )
        gd%svix = GRID_UNDEFINED
        gd%sviy = GRID_UNDEFINED
        gd%svi = GRID_UNDEFINED

        allocate( gd%mapCvixVx(nvx, 8), gd%mapCviyVx(nvx, 8))
        gd%mapCvixVx = B2_GRID_UNDEFINED
        gd%mapCviyVx = B2_GRID_UNDEFINED

    end subroutine allocateB2GridMap

    !> Deallocate B2GridMap
    subroutine deallocateB2GridMap( gd )
        type(B2GridMap), intent(inout) :: gd   !< The grid mapping as computed
            !< by b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS

        deallocate( gd%mapCvI )
        deallocate( gd%mapFcI )
        deallocate( gd%mapVxI )

        deallocate( gd%mapCvix, gd%mapCviy )
        deallocate( gd%mapFcix, gd%mapFciy, gd%mapFcIFace )
        deallocate( gd%mapVxix, gd%mapVxiy )

        deallocate( gd%svix, gd%sviy, gd%svi )

        deallocate( gd%mapCvixVx, gd%mapCviyVx )

    end subroutine deallocateB2GridMap

    !> Create B2GridMap, containing grid geometry information
    subroutine b2CreateMap( nx, ny, crx, cry, cflag, leftix, leftiy,    &
            &   rightix, rightiy, topix, topiy, bottomix, bottomiy,     &
            &   includeGhostCells, gd)

        use b2mod_cellhelper

        !! Input arguments (unchanged on exit)
        !! Size of grid arrays: (-1:nx, -1:ny)
        integer :: nx   !< Specifies the number of interior cells
                        !< along the first coordinate (poloidal)
        integer :: ny   !< Specifies the number of interior cells
                        !< along the second coordinate (radial)

        !! Output arguments
        !! vertex coordinates
        real(R8), intent(in) :: crx( -1:nx, -1:ny, 0:3 ) !< Horizontal vertex
            !< coordinates of the four corners of the (ix, iy) cell
        real(R8), intent(in) :: cry( -1:nx, -1:ny, 0:3 ) !< Vertical vertex
            !< coordinates of the four corners of the (ix, iy) cell
        integer cflag( -1:nx, -1:ny, CARREOUT_NCELLFLAGS )
        integer, intent(in) :: leftix( -1:nx, -1:ny)    !< Left neighbour
            !< poloidal (first coordinate) index array
        integer, intent(in) :: leftiy( -1:nx, -1:ny )   !< Left neighbour radial
            !< (second coordinate) index
        integer, intent(in) :: rightix( -1:nx, -1:ny )  !< Right neighbour
            !< poloidal (first coordinate) index array
        integer, intent(in) :: rightiy( -1:nx, -1:ny )  !< Right neighbour
            !< radial (second coordinate) index
        integer, intent(in) :: topix( -1:nx, -1:ny )    !< Top neighbour
            !< poloidal (first coordinate) index array
        integer, intent(in) :: topiy( -1:nx, -1:ny )    !< Top neighbour radial
            !< (second coordinate) index
        integer, intent(in) :: bottomix( -1:nx, -1:ny ) !< Bottom neighbour
            !< poloidal (first coordinate) index array
        integer, intent(in) :: bottomiy( -1:nx, -1:ny ) !< Bottom neighbour
            !< radial (second coordinate) index
        logical, intent(in) :: includeGhostCells    !< Include "fake" cells

        type(B2GridMap), intent(inout) :: gd    !< The grid mapping as computed
            !< by b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS

        !! internal variables
        integer :: ix   !< x-aligned (poloidal) cell index
        integer :: iy   !< y-aligned (radial) cell index
        integer :: ic   !< Cell index
        integer :: ifcx !< x-aligned face index
        integer :: ifcy !< y-aligned face index
        integer :: ivx  !< Vertex/node index
        integer :: i1   !< Iterator
        integer :: i2   !< Iterator
        integer :: nbix
        integer :: nbiy
        integer :: nbix2
        integer :: nbiy2
        integer :: iFace !< Face/edge index
        integer :: iCorner
        integer :: index
        integer :: iPass

        !! There is an additional strip of "fake" cells (even faker than
        !! ghost cells) at the top and right to be able to also write out the
        !! ghost cells
        integer, dimension(-1:nx, -1:ny) :: cvi !< Control volume index
        integer, dimension(-1:nx, -1:ny, 0:3) :: fcxi   !< x-aligned face index
            !< (only BOTTOM and TOP used)
        integer, dimension(-1:nx, -1:ny, 0:3) :: fcyi   !< Y-aligned face index
            !< (only LEFT and RIGHT used)
        integer, dimension(-1:nx+1, -1:ny+1, 0:3) :: vxi !< vertex index

        logical :: cvNeeded( (nx + 2) * (ny + 2) )
        logical :: vxNeeded( (nx + 2) * (ny + 2) * 4)
        logical :: fcXNeeded( (nx + 2) * (ny + 2) * 2)
        logical :: fcYNeeded( (nx + 2) * (ny + 2) * 2)

        integer :: fcxiReduce((nx + 2) * (ny + 2) * 2)
        integer :: fcyiReduce((nx + 2) * (ny + 2) * 2)
        integer :: vxiReduce((nx + 2) * (ny + 2) * 4)
        integer :: nsector((nx + 2) * (ny + 2) * 4)

        logical :: check, cell_done

        !! list of identified special vertices
        !! vertex indices (ix, iy)
        integer, dimension(MAX_SPECIAL_VERTICES) :: svix !< Array of horizontal
            !< (x-aligned) positions of special vertices in B2 data structure
        integer, dimension(MAX_SPECIAL_VERTICES) :: sviy !< Array of vertical
            !< (y-aligned) positions of special vertices in B2 data structure
        integer, dimension(MAX_SPECIAL_VERTICES) :: svixAlias
        integer, dimension(MAX_SPECIAL_VERTICES) :: sviyAlias
        integer :: svc  !< Number of special vertices
        integer :: svcDuplicates    !< number of duplicate special vertices

        !! numbers of unique objects
        integer :: ncv  !< Number of all cells (2D objects)
        integer :: nfcx !< Number of x-aligned faces/edges (1D objects)
        integer :: nfcy !< Number of y-aligned faces/edges (1D objects)
        integer :: nvx  !< Number of all vertices/nodes (0D objects)

        call logmsg( LOGDEBUG, "b2CreateMap: create map for a nx="  &
            &   //int2str(nx)//", ny="//int2str(ny)//" b2 grid" )

        cvi = GRID_UNDEFINED
        fcxi = GRID_UNDEFINED
        fcyi = GRID_UNDEFINED
        vxi = GRID_UNDEFINED

        !! set up initial lexicographic indices
        ic = 0 !! cvs
        ifcx = 0 !! x-aligned faces
        ifcy = 0 !! y-aligned faces
        ivx = 0 !! vertices
        do ix = -1, nx
            do iy = -1, ny
                !! Cell
                ic = ic + 1
                cvi( ix, iy ) = ic

                !! Faces: associate left and bottom face with every cell
                !! where possible
                ifcy = ifcy + 1
                fcyi( ix, iy, LEFT ) = ifcy !! left face
                ifcx = ifcx + 1
                fcxi( ix, iy, BOTTOM ) = ifcx !! bottom face

                !! Vertices: associate bottom left vertex with every cell
                !! where possible
                call findExistingVertexIndex( ix, iy, VX_LOWERLEFT, index )
                if(index == GRID_UNDEFINED) then
                    ivx = ivx + 1
                    vxi( ix, iy, VX_LOWERLEFT ) = ivx
                else
                    vxi( ix, iy, VX_LOWERLEFT ) = index
                end if
            end do
        end do

        !! fill in index numbers for remaining faces
        do ix = -1, nx
            do iy = -1, ny

                !! Right face: left face of left neighbour
                call getNeighbour(nx, ny, leftix, leftiy, rightix, rightiy, &
                    &   topix, topiy, bottomix, bottomiy,                   &
                    &   ix, iy, RIGHT, nbix, nbiy)

                if( isCellInDomain(nx, ny, nbix, nbiy) ) then
                    fcyi( ix, iy, RIGHT ) = fcyi( nbix, nbiy, LEFT )
                else
                    ifcy = ifcy + 1
                    fcyi( ix, iy, RIGHT ) = ifcy
                end if

                !! Top face: bottom face of top neighbour
                !! also top-left vertex
                call getNeighbour(nx, ny, leftix, leftiy, rightix, rightiy, &
                    & topix, topiy, bottomix, bottomiy,                        &
                    & ix, iy, TOP, nbix, nbiy)

                if( isCellInDomain(nx, ny, nbix, nbiy) ) then
                    fcxi( ix, iy, TOP ) = fcxi( nbix, nbiy, BOTTOM )
                else
                    ifcx = ifcx + 1
                    fcxi( ix, iy, TOP ) = ifcx
                end if
            end do
        end do

        call xertst (count(fcxi(:,:,BOTTOM) == GRID_UNDEFINED) == 0,    &
            &  "b2CreateMap: there are unnumbered x-aligned faces" )
        call xertst (count(fcxi(:,:,TOP) == GRID_UNDEFINED) == 0,       &
            &  "b2CreateMap: there are unnumbered x-aligned faces" )
        call xertst (count(fcyi(:,:,LEFT) == GRID_UNDEFINED) == 0,      &
            &   "b2CreateMap: there are unnumbered y-aligned faces" )
        call xertst (count(fcyi(:,:,RIGHT) == GRID_UNDEFINED) == 0,     &
            &   "b2CreateMap: there are unnumbered y-aligned faces" )

        !! Fill in vertex numbers for remaining vertices
        !! A vertex can be shared among 4 cells (possibly more for special
        !! vertices, but they are assumed to have a proper cell associated
        !! with them)
        do ix = -1, nx
            do iy = -1, ny

                do iCorner = VX_LOWERRIGHT, VX_UPPERRIGHT  ! 1, 3
                    call findExistingVertexIndex( ix, iy, iCorner, index )
                    if( index /= GRID_UNDEFINED ) then
                        vxi( ix, iy, iCorner ) = index
                    else
                        ivx = ivx + 1
                        vxi( ix, iy, iCorner ) = ivx
                    end if
                end do

            end do
        end do

        call xertst( count( vxi( -1:nx, -1:ny, 0:3 ) == GRID_UNDEFINED) == 0,   &
            & "b2CreateMap: there are unnumbered vertices" )

        !! Mark which cells, vertices and faces are needed in 1d lists
        cvNeeded = .false.
        fcXNeeded = .false.
        fcYNeeded = .false.
        vxNeeded = .false.
        do ix = -1, nx
            do iy = -1, ny
                if( .not. isUnneededCell( nx, ny, cflag,   &
                        &   includeGhostCells, ix, iy ) ) then
                    cvNeeded(cvi(ix, iy)) = .true.
                    do iCorner = 0, 3
                        vxNeeded(vxi( ix, iy, iCorner ) ) = .true.
                    end do
                    fcYNeeded( fcyi( ix, iy, LEFT ) ) = .true.
                    fcYNeeded( fcyi( ix, iy, RIGHT ) ) = .true.
                    fcXNeeded( fcxi( ix, iy, BOTTOM ) ) = .true.
                    fcXNeeded( fcxi( ix, iy, TOP ) ) = .true.
                end if
            end do
        end do

        !! search x-points.
        svc = 0
        do ix = -1, nx
            do iy = -1, ny
                !! do not do this for unneeded cells, might be not initialized
                !! This is not perfect, but special vertices are usually inside
                !! the domain, and there it should be ok.
                if( .not. cvNeeded( cvi( ix, iy ) ) ) cycle

                !! test whether vertex is special vertex
                if( isSpecialVertex( ix, iy ) ) then
                    svc = svc + 1
                    svix( svc ) = ix
                    sviy( svc ) = iy
                end if

            end do
        end do

        !! identify duplicate special vertices

        svixAlias = 0
        sviyAlias = 0
        svcDuplicates = 0

        do i1 = 1, svc - 1

            ix = svix( i1 )
            iy = sviy( i1 )

            !! has the point already been identified as unneeded?
            if( .not. vxNeeded( vxi( ix, iy, VX_LOWERLEFT ) ) ) cycle

                !! compare with remaining vertices
                do i2 = i1 + 1, svc

                !! has the point already been identified as unneeded?
                if( .not. vxNeeded( vxi( svix(i2), sviy(i2),   &
                    &   VX_LOWERLEFT ) ) ) cycle

                if( points_match( crx( ix, iy, 0 ), cry( ix, iy, 0 ),   &
                    &   crx( svix(i2), sviy(i2), 0 ), cry( svix( i2 ),  &
                    &   sviy( i2 ), 0 ) ) ) then

                    !! The special vertex with index i2 is a duplicate of the
                    !! one with index i1.
                    !! Mark as unneeded and set up all references to this point
                    !! as aliased to the first one.
                    !vxNeeded( vxi(svix(i2), sviy(i2), VX_LOWERLEFT) ) = .false.
                    ! where (vxi == vxi( svix(i2), sviy(i2), VX_LOWERLEFT ))    &
                    !    &   vxi = vxi( ix, iy, VX_LOWERLEFT )

                    svixAlias( i2 ) = ix
                    sviyAlias( i2 ) = iy
                    !! Bookkeeping for diagnostics
                    svcDuplicates = svcDuplicates + 1
                end if
            end do
        end do

        if( svc - svcDuplicates .eq. 1 ) then
            call logmsg( LOGDEBUG, "b2CreateMap: found "            &
                &   //int2str( svc - svcDuplicates )//              &
                &   " special vertex (x-point), "//int2str( svc )   &
                &   //" including duplicates." )
        else
            call logmsg( LOGDEBUG, "b2CreateMap: found "            &
                &   //int2str( svc - svcDuplicates )                &
                &   //" special vertices (x-points), "              &
                &   //int2str( svc )//" including duplicates." )
        end if

        !! number of unique cells
        !ncv = ( nx + 2 ) * ( ny + 2 ) - count( .not. isNeeded( cvi ) )
        ncv = count( cvNeeded )
        !! number of unique faces (x and y direction)
        nfcx = count( fcXNeeded )
        nfcy = count( fcYNeeded )
        !! number of unique vertices
        nvx = count( vxNeeded )

        !! allocate the mapping structure
        call allocateB2GridMap( gd, nx, ny, ncv, nfcx, nfcy, nvx )

        !! build the mappings

        !! ...for cells
        ic = 0
        gd%mapCvI = GRID_UNDEFINED
        do ix = -1, nx
            do iy = -1, ny
                if( cvNeeded( cvi( ix, iy ) ) ) then
                    ic = ic + 1
                    call xertst ( ic .le. ncv , &
                        & 'b2CreateMap: found more cells than expected' )
                    !! Map CPO -> B2
                    gd%mapCvix( ic ) = ix
                    gd%mapCviy( ic ) = iy
                    !! Map B2 -> CPO
                    gd%mapCvI( ix, iy ) = ic
                else
                    !! not needed
                    gd%mapCvI( ix, iy ) = GRID_UNDEFINED
                end if
            end do
        end do
        call xertst ( ic == ncv , 'b2CreateMap: found less cells than expected' )

        !! ...for x faces
        !! first set up numbering
        ic = 0
        fcxiReduce = GRID_UNDEFINED
        do i1 = 1, size( fcXNeeded )
            if( fcXNeeded( i1 ) ) then
                ic = ic + 1
                call xertst ( ic .le. nfcx , &
                    & 'b2CreateMap: found more x-aligned faces than expected' )
                fcxiReduce( i1 ) = ic
            end if
        end do
        !! sanity check for number of x faces
        call xertst ( ic == nfcx , &
            & 'b2CreateMap: found less x-aligned faces than expected' )

        !! We want a CPO face to point to a LEFT or BOTTOM cell face if possible.
        !! If no LEFT or BOTTOM face exists, alternatively point to a TOP or
        !! RIGHT face. For this we need two passes.

        gd%mapFcI = B2_GRID_UNDEFINED

        !! first x-aligned faces
        gd%mapFcix = B2_GRID_UNDEFINED !! cannot use GRID_UNDEFINED here...

        do iPass = 1, 2
            do ix = -1, nx
                do iy = -1, ny

                    if(iPass == 1) then
                        iFace = BOTTOM
                    else
                        iFace = TOP
                    end if

                    if( fcXNeeded( fcxi( ix, iy, iFace ) ) ) then
                        !if(isUnneededCell(nx, ny, cflag,  &
                        !   &   includeGhostCells, ix, iy)) cycle

                        !! Get the CPO index for this face
                        ic = fcxiReduce( fcxi( ix, iy, iFace ) )
                        !! Map CPO -> B2
                        if( gd%mapFcix( ic ) == B2_GRID_UNDEFINED) then
                            gd%mapFcix( ic ) = ix
                            gd%mapFciy( ic ) = iy
                            gd%mapFcIFace( ic ) = iFace
                        end if
                        !! Map B2 -> CPO
                        gd%mapFcI( ix, iy, iFace ) = ic
                    end if
                end do
            end do
        end do

        !! ...for y faces, continue counting the total faces in ic
        !! again, first set up numbering
        ic = nfcx
        fcyiReduce = GRID_UNDEFINED
        do i1 = 1, size( fcYNeeded )
            if( fcYNeeded( i1 ) ) then
                ic = ic + 1
                fcyiReduce( i1 ) = ic
            end if
        end do
        call xertst ( ic == nfcx + nfcy,    &
            & 'b2CreateMap: found less y-aligned faces than expected' )

        do iPass = 1, 2
            do ix = -1, nx
                do iy = -1, ny
                    if(iPass == 1) then
                        iFace = LEFT
                    else
                        iFace = RIGHT
                    end if
                    if( fcYNeeded( fcyi( ix, iy, iFace ) ) ) then
                        !! do not associate with unused cell
                        !if(isUnneededCell(nx, ny, cflag,  &
                        !   &   includeGhostCells, ix, iy)) cycle

                        !! Get the CPO index for this face
                        ic = fcyiReduce( fcyi( ix, iy, iFace ) )
                        !! Map CPO -> B2
                        if( gd%mapFcix( ic ) == B2_GRID_UNDEFINED) then
                            gd%mapFcix( ic ) = ix
                            gd%mapFciy( ic ) = iy
                            gd%mapFcIFace( ic ) = iFace
                        end if
                        !! Map B2 -> CPO
                        gd%mapFcI( ix, iy, iFace ) = ic
                    end if
                end do
            end do
        end do

        if( count( gd%mapFcI == B2_GRID_UNDEFINED) > 0) then
            call logmsg( LOGKNOWNWARNING, "b2CreateMap: have faces with "//  &
                &  "missing mapping data (ignored)" )
            ! call xerrab ( "b2CreateMap: have faces with missing mapping data" )
        endif

        !! vertices
        !! Like for faces, first set up unique numbering
        ic = 0
        vxiReduce = GRID_UNDEFINED
        do i1 = 1, size( vxNeeded )
            if( vxNeeded( i1 ) ) then
                ic = ic + 1
                vxiReduce( i1 ) = ic
            end if
        end do
        call xertst ( ic == nvx , &
            & 'b2CreateMap: did not find expected number of vertices' )

        gd%mapVxI = GRID_UNDEFINED
        gd%mapVxix = B2_GRID_UNDEFINED
        !! First looping over the corners makes sure that all vertices for which
        !! this is possible are associated with a lower-left vertex
        do iCorner = VX_LOWERLEFT, VX_UPPERRIGHT ! 0, 3
            do ix = -1, nx
                do iy = -1, ny
                    if( vxNeeded( vxi( ix, iy, iCorner ) ) ) then
                        ic = vxiReduce( vxi( ix, iy, iCorner ) )
                        !! Map CPO -> B2
                        if( gd%mapVxix( ic ) == B2_GRID_UNDEFINED ) then
                            gd%mapVxix( ic ) = ix
                            gd%mapVxiy( ic ) = iy
                            gd%mapVxIVx( ic ) = iCorner
                        end if
                        !! Map B2 -> CPO
                        call xertst ( ic /= GRID_UNDEFINED,     &
                            & "b2CreateMap: found vertex pointing to ic = 0" )
                        gd%mapVxI( ix, iy, iCorner ) = ic
                    end if
                end do
            end do
        end do

        !! sanity check: are all vertices initialized?
        do ic = 1, nvx
            if( gd%mapVxix( ic ) == B2_GRID_UNDEFINED ) then
                call logmsg( LOGWARNING, "b2CreateMap: vertex " &
                    &   //int2str( ic )//" not properly initialized (1)")
            end if
            if( gd%mapVxIVx( ic ) == B2_GRID_UNDEFINED ) then
                call logmsg( LOGWARNING, "b2CreateMap: vertex " &
                    &   //int2str( ic )//" not properly initialized (2)")
            end if
        end do

        call xertst ( count( gd%mapVxix == B2_GRID_UNDEFINED ) == 0,    &
            & "b2CreateMap: mapVxix broken" )

        !! After completing the vertex map, we can complete the special
        !! vertex (x-point) map

        !! The special vertex with index i1 is unique and needed. Copy it to
        !! the map structure identify duplicate special vertices
        gd%nsv = 0
        do ic = 1, svc
            ix = svix( ic )
            iy = sviy( ic )

            !! has the point already been identified as unneeded?
            if( .not. vxNeeded( vxi( ix, iy, VX_LOWERLEFT ) ) ) cycle

            !! if it is needed, copy it to the mapping structure
            gd%nsv = gd%nsv + 1
            gd%svix(gd%nsv) = ix
            gd%sviy(gd%nsv) = iy
            gd%svi(ic) = gd%mapVxI( gd%svix(ic), gd%sviy(ic), VX_LOWERLEFT )
        end do

        !! Now fill in missing vertex numbers where possible
        !! FIXME: do the same for the faces?
        vxi = gd%mapVxI
        do ix = -1, nx
            do iy = -1, ny
                do iCorner = VX_LOWERRIGHT, VX_UPPERRIGHT  ! 1, 3
                    if( gd%mapVxI( ix, iy, iCorner) == GRID_UNDEFINED ) then
                        call findExistingVertexIndex( ix, iy, iCorner,  &
                            &   index, stopOnUnneededCells = .true.)
                        if( index /= GRID_UNDEFINED) then
                            gd%mapVxI( ix, iy, iCorner ) = index
                        end if
                    end if
                end do
            end do
        end do

        call logmsg( LOGDEBUG, "b2CreateMap: finished setting up map" )
        call logmsg( LOGDEBUG, "b2CreateMap: map contains  "    &
            &   //int2str(gd%ncv)//" unique cells")
        call logmsg( LOGDEBUG, "b2CreateMap: map contains  "    &
            &   //int2str(gd%nfcx)//" unique x-aligned faces")
        call logmsg( LOGDEBUG, "b2CreateMap: map contains  "    &
            &   //int2str(gd%nfcy)//" unique y-aligned faces")
        call logmsg( LOGDEBUG, "b2CreateMap: map contains  "    &
            &   //int2str(gd%nvx)//" unique vertices")

    !! Find out how many and which control volumes touch any given vertex
        nsector = 0
        do ix = -1, nx
            do iy = -1, ny
                do iCorner = VX_LOWERLEFT, VX_UPPERRIGHT
                    if( vxi( ix, iy, iCorner ) == GRID_UNDEFINED ) cycle
                    nsector ( gd%mapVxI( ix, iy, iCorner) ) = 1
                    gd%mapCvixVx( gd%mapVxI( ix, iy, iCorner), 1 ) = ix
                    gd%mapCviyVx( gd%mapVxI( ix, iy, iCorner), 1 ) = iy
                    do i1 = -1, nx
                        do i2 = -1, ny
                            cell_done = i1 .eq. ix .and. i2 .eq. iy
                            do ic = VX_LOWERLEFT, VX_UPPERRIGHT
                                if( vxi( i1, i2, ic ) == GRID_UNDEFINED ) cycle
                                if( cell_done ) cycle
                                if( points_match( crx( ix, iy, iCorner ),  &
                                    &   cry( ix, iy, iCorner ),             &
                                    &   crx( i1, i2, ic), cry( i1, i2, ic ))) then
                                    nsector (                                   &
                                        &   gd%mapVxI( ix, iy, iCorner )) =     &
                                        &   nsector(                            &
                                        &   gd%mapVxI( ix, iy, iCorner )) + 1

                                    gd%mapCvixVx(                               &
                                        &   gd%mapVxI( ix, iy, iCorner ),       &
                                        &   nsector(                            &
                                        &   gd%mapVxI( ix, iy, iCorner ))) = ix

                                    gd%mapCviyVx(                               &
                                        &   gd%mapVxI( ix, iy, iCorner ),       &
                                        &   nsector(                            &
                                        &   gd%mapVxI( ix, iy, iCorner ) ) ) = iy
                                    cell_done = .true.
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do

        return

contains

    !> For a given corner vertex of a cell, check whether in any connected
    !! cell a vertex index was already assigned to this vertex
    subroutine findExistingVertexIndex( ix, iy, iCorner, index, &
            &   stopOnUnneededCells)
        integer, intent(in) :: ix   !< Specifies index of interior cell along
                                    !< the first coordinate
        integer, intent(in) :: iy   !< Specifies index of interior cell along
                                    !< the second coordinate
        integer, intent(in) :: iCorner
        integer, intent(out) :: index
        logical, intent(in), optional :: stopOnUnneededCells

        !! internal
        integer :: iRot
        integer :: iStep
        integer :: nix
        integer :: niy
        integer :: nix2
        integer :: niy2
        integer :: iDir
        integer :: nICorner
        logical :: lStopOnUnneededCells

        lStopOnUnneededCells = .false.
        if( present(stopOnUnneededCells) ) lStopOnUnneededCells =  &
            &   stopOnUnneededCells

        !! already have index?
        index = vxi( ix, iy, iCorner )
        if(index /= GRID_UNDEFINED) return

        !! Circle through the cells connected to the corner vertex with
        !! index iCorner in clockwise and counterclockwise direction
        do iRot = CLOCKWISE, COUNTERCLOCKWISE
            nix = ix
            niy = iy
            nICorner = iCorner

            do iStep = 1, 4
                !! take step
                iDir = VXCIRCLE_STEPDIR( iStep, iCorner, iRot )

                call getNeighbour( nx, ny, leftix, leftiy, rightix, &
                    &   rightiy, topix, topiy, bottomix, bottomiy,  &
                    &   nix, niy, iDir, nix2, niy2 )

                if(.not. isCellInDomain( nx, ny, nix2, niy2 )) exit
                if( lStopOnUnneededCells .and. isUnneededCell( nx, ny,     &
                    &   cflag, includeGhostCells, nix2, niy2) ) exit

                nix = nix2
                niy = niy2
                nICorner = VXCORNER_NEXTINDEX( iDir, nICorner )

                if( points_match( crx( ix, iy, iCorner ),                   &
                    &   cry( ix, iy, iCorner ), crx(nix, niy, nICorner ),   &
                    &   cry( nix, niy, nICorner ) ) ) then
                    index = vxi( nix, niy, nICorner)
                    if( index /= GRID_UNDEFINED ) return
                end if
            end do
        end do

        index = GRID_UNDEFINED
    end subroutine findExistingVertexIndex

    !> Test whether the vertex associated with the cell (ix, iy) is special,
    !! i.e. a x-point
    !! This steps around the vertex (left-bottom-right-top) and checks
    !! whether the resulting position is equal to the starting position
    logical function isSpecialVertex( ix, iy )
        integer, intent(in) :: ix   !< Specifies index of interior cell along
                                    !< the first coordinate
        integer, intent(in) :: iy   !< Specifies index of interior cell along
                                    !< the second coordinate

        !! internal
        integer :: x, y, xn, yn

        x = ix
        y = iy
        !! TODO: this can be simplified with the new structures in
        !! b2mod_cellhelper

        !! step left
        call getNeighbour( nx, ny, leftix, leftiy, rightix, rightiy,    &
            &   topix, topiy, bottomix, bottomiy, x, y, LEFT, xn, yn )
        if( .not. isCellInDomain( nx, ny, xn, yn, extended =   &
            &   includeGhostCells )) then
            isSpecialVertex = .false.
            return
        end if
        if( isGhostCell( cflag( xn, yn, CELLFLAG_TYPE ) ) ) then
            isSpecialVertex = .false.
            return
        end if
        x = xn
        y = yn

        !! step bottom
        call getNeighbour( nx, ny, leftix, leftiy, rightix, rightiy,    &
            &   topix, topiy, bottomix, bottomiy, x, y, BOTTOM, xn, yn )
        if( .not. isCellInDomain( nx, ny, xn, yn, extended =   &
            &   includeGhostCells ) ) then
            isSpecialVertex = .false.
            return
        end if
        if( isGhostCell(cflag(xn,yn,CELLFLAG_TYPE)) ) then
            isSpecialVertex = .false.
            return
        end if
        x = xn
        y = yn

        !! step right
        call getNeighbour(nx, ny, leftix, leftiy, rightix, rightiy, &
            &   topix,topiy,bottomix,bottomiy, x, y, RIGHT, xn, yn )
        if( .not. isCellInDomain( nx, ny, xn, yn, extended =   &
            &   includeGhostCells )) then
            isSpecialVertex = .false.
            return
        end if
        if( isGhostCell(cflag(xn,yn,CELLFLAG_TYPE)) ) then
            isSpecialVertex = .false.
            return
        end if
        x = xn
        y = yn

        !! step top
        call getNeighbour(nx, ny, leftix, leftiy, rightix, rightiy, &
            &   topix, topiy, bottomix, bottomiy, x, y, TOP, xn, yn )
        if( .not. isCellInDomain( nx, ny, xn, yn, extended =   &
            &   includeGhostCells)) then
            isSpecialVertex = .false.
            return
        end if
        if( isGhostCell(cflag(xn,yn,CELLFLAG_TYPE)) ) then
            isSpecialVertex = .false.
            return
        end if
        x = xn
        y = yn

        !! do we end up where we left?
        isSpecialVertex = .not. ( ( x == ix ) .and. ( y == iy ) )

    end function isSpecialVertex

  end subroutine b2CreateMap

    !> test whether cell (ix,iy) is actually used
    function isUnneededCell( nx, ny, cflag, includeGhostCells, ix, iy )
        integer, intent(in) :: nx   !< Specifies the number of interior cells
                                    !< along the first coordinate
        integer, intent(in) :: ny   !< Specifies the number of interior cells
                                    !< along the second coordinate
        integer, intent(in) :: ix   !< Specifies index of interior cell along
                                    !< the first coordinate
        integer, intent(in) :: iy   !< Specifies index of interior cell along
                                    !< the second coordinate
        integer, intent(in) :: cflag(-1:nx,-1:ny, CARREOUT_NCELLFLAGS)
        logical, intent(in) :: includeGhostCells    !< Include "fake" cells
        logical isUnneededCell

        !! Only cells inside the "normal" B2 domain can be needed
        !! (this catches fake cells and connectivity pointing outside the domain)
        isUnneededCell = .not. isCellInDomain(nx, ny, ix, iy, extended = .true.)
        if(isUnneededCell) return

        if(includeGhostCells) then
            isUnneededCell = isUnusedCell( cflag(ix,iy,CELLFLAG_TYPE) )
        else
            isUnneededCell = isUnusedCell( cflag(ix,iy,CELLFLAG_TYPE) ) &
                &   .or. isGhostCell( cflag(ix,iy,CELLFLAG_TYPE) )
        end if

        !! Classical treatment (without cflag, no extended grid) - for reference
        !if(.not. includeGhostCells) then
        !    isUnneededCell = &
        !        & ( leftix( ix, iy ) == -2 ) &
        !        & .or. ( rightix( ix, iy ) == ( nx + 1 ) ) &
        !        & .or. ( bottomiy( ix, iy ) == -2 ) &
        !        & .or. ( topiy( ix, iy ) == ( ny + 1 ) )
        !end if

    end function isUnneededCell


    !> Check whether the cell at position (ix,iy) is inside the 'classical'
    !! B2 grid.
    !! The default is to check whether it is inside  the extended grid
    !! (including the ghost cells).
    !! If the optional parameter extended is given, extended = .false. will
    !! check whether the position is inside the actual physical domain of the
    !! grid.
    function isCellInDomain( nx, ny, ix, iy, extended )
        integer, intent(in) :: nx   !< Specifies the number of interior cells
                                    !< along the first coordinate
        integer, intent(in) :: ny   !< Specifies the number of interior cells
                                    !< along the second coordinate
        integer, intent(in) :: ix   !< Specifies index of interior cell along
                                    !< the first coordinate
        integer, intent(in) :: iy   !< Specifies index of interior cell along
                                    !< the second coordinate
        logical, intent(in), optional :: extended
        logical :: isCellInDomain

        !! internal
        logical :: lExtended

        lExtended = .true.
        if( present( extended ) ) lExtended = extended

        if( lExtended ) then
            !! in extended domain (including ghost cells)?
            isCellInDomain = ( ix >= -1 ) .and. (ix <= nx) .and.    &
                &   ( iy >= -1 ) .and. ( iy <= ny )
        else
            isCellInDomain = ( ix > -1 ) .and. (ix < nx) .and.      &
                &   ( iy > -1 ) .and. ( iy < ny )
        end if
    end function isCellInDomain


    !> Check whether the node associate with the cell at position (ix,iy) is
    !! included in the grid.
    function isNodeInDomain( nx, ny, ix, iy, extended )
        integer, intent(in) :: nx   !< Specifies the number of interior cells
                                    !< along the first coordinate
        integer, intent(in) :: ny   !< Specifies the number of interior cells
                                    !< along the second coordinate
        integer, intent(in) :: ix   !< Specifies index of interior cell along
                                    !< the first coordinate
        integer, intent(in) :: iy   !< Specifies index of interior cell along
                                    !< the second coordinate
        logical, intent(in), optional :: extended
        logical :: isNodeInDomain

        !! internal
        logical :: lExtended

        lExtended = .true.
        if( present( extended ) ) lExtended = extended

        if( lExtended ) then
            !! in extended domain (including ghost cells)?
            isNodeInDomain = ( ix >= -1 ) .and. (ix <= nx + 1) .and.    &
                &   ( iy >= -1 ) .and. ( iy <= ny + 1 )
        else
            isNodeInDomain = ( ix > -1 ) .and. (ix < nx + 1) .and.      &
                &   ( iy > -1 ) .and. ( iy < ny + 1)
        end if
    end function isNodeInDomain

    !> extended neighbourhood mappings
    subroutine getNeighbour(nx, ny, leftix, leftiy, rightix, rightiy,   &
            &   topix, topiy, bottomix,bottomiy, ix, iy, dir, nbix, nbiy )
        integer, intent(in) :: nx   !< Specifies the number of interior cells
                                    !< along the first coordinate
        integer, intent(in) :: ny   !< Specifies the number of interior cells
                                    !< along the second coordinate
        integer, intent(in) :: ix   !< Specifies index of interior cell along
                                    !< the first coordinate
        integer, intent(in) :: iy   !< Specifies index of interior cell along
                                    !< the second coordinate
        integer, intent(in) :: dir
        integer, intent(in) :: leftix( -1:nx, -1:ny)    !< Left neighbour
            !< poloidal (first coordinate) index array
        integer, intent(in) :: leftiy( -1:nx, -1:ny )   !< Left neighbour radial
            !< (second coordinate) index
        integer, intent(in) :: rightix( -1:nx, -1:ny )  !< Right neighbour
            !< poloidal (first coordinate) index array
        integer, intent(in) :: rightiy( -1:nx, -1:ny )  !< Right neighbour
            !< radial (second coordinate) index
        integer, intent(in) :: topix( -1:nx, -1:ny )    !< Top neighbour
            !< poloidal (first coordinate) index array
        integer, intent(in) :: topiy( -1:nx, -1:ny )    !< Top neighbour radial
            !< (second coordinate) index
        integer, intent(in) :: bottomix( -1:nx, -1:ny ) !< Bottom neighbour
            !< poloidal (first coordinate) index array
        integer, intent(in) :: bottomiy( -1:nx, -1:ny ) !< Bottom neighbour
            !< radial (second coordinate) index
        integer, intent(out) :: nbix
        integer, intent(out) :: nbiy

        select case(dir)
        case(LEFT)
            nbix = leftix(ix, iy)
            nbiy = leftiy(ix, iy)
        case(BOTTOM)
            nbix = bottomix(ix, iy)
            nbiy = bottomiy(ix, iy)
        case(RIGHT)
            nbix = rightix(ix, iy)
            nbiy = rightiy(ix, iy)
        case(TOP)
            nbix = topix(ix, iy)
            nbiy = topiy(ix, iy)
        end select

    end subroutine getNeighbour

end module b2mod_grid_mapping

!!!Local Variables:
!!! mode: f90
!!! End:
