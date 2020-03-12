!!-----------------------------------------------------------------------------
!! DOCUMENTATION:
!>      @section b2mod_gmap_desc  Description
!!      Module providing a mechanism to map from the B2 data structure to the
!!      IDS and CPO data structures, consisting of a routine to set up the map
!!      (b2CreateMap), a data structure to hold the map information (B2GridMap)
!!      and some service routines to handle this data structure.
!!
!!      @subsection b2mod_gmap_pv  Parameters/variables
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

    !! Alignment index (for example in B2 flux arrays)
    integer, parameter :: ALIGNX = 1
    integer, parameter :: ALIGNY = 0

    !! Maximum number of special vertices expected in the grid
    integer, parameter :: MAX_SPECIAL_VERTICES = 10

    !> Data structure holding an intermediate grid description to be
    !! transferred into a CPO or IDS
    type B2GridMap
        integer :: nCv  !< Number of all cells in the domain (2D objects)
        integer :: nFcx !< Number of x-aligned faces/edges in the domain
                        !< (1D objects)
        integer :: nFcy !< Number of y-aligned faces/edges in the domain
                        !< (1D objects)
        integer :: nVx  !< Number of all vertices/nodes in the domain
                        !<(0D objects) (nVx = ( nx+1 )*( ny+1 ) - 1 )
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
            !< composing the quadrilateral in the list.
            !< gmap%mapFcI( ix, iy, POSITION ), where POSITION is LEFT, BOTTOM,
            !< RIGHT or TOP
        integer, dimension(:,:,:), allocatable :: mapVxI

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
!!$        ! b2Cv( mapCvix(i), mapCviy(i) ) = cpoCv(i)
!!$        ! b2Fc( mapFcix(i), mapFciy(i), mapFciFace(i) ) = cpoFc(i)
!!$        !! for "normal" faces, mapFcIFace will be LEFT or BOTTOM

!!$        ! b2Vx( mapVxix(i), mapVxiy(i), mapVxIVx(i) ) = cpoVx(i)
!!$        !! for "normal" vertices, mapVxIVx(i) will be 1 (lower left vertex)

!!$        !! 2d B2 data structure -> 1d CPO lists ( (ix, iy) -> i )
!!$        ! cpoCv( mapCvI(ix, iy) ) = b2Cv(ix, iy)
!!$        ! cpoFc( mapFcI(ix, iy, iFace) ) = b2Fc(ix, iy, iFace)
!!$        ! cpoVx( mapVxI(ix, iy, iVertex) ) = b2Vx(ix, iy, iVertex)

contains

    !! service routines for B2GridData

    !> Set B2GridMap type, intended to be filled with grid geometry information
    !! using b2CreateMap subroutine
    subroutine allocateB2GridMap( gd, nx, ny, nCv, nFcx, nFcy, nVx )
        type(B2GridMap), intent(inout) :: gd    !< The grid mapping as computed
            !< by b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS
        integer, intent(in) :: nx  !< Specifies the number of interior cells
                                   !< along the first coordinate
        integer, intent(in) :: ny  !< Specifies the number of interior cells
                                   !< along the second coordinate
        integer, intent(in) :: nCv     !< Number of all cells (2D objects)
        integer, intent(in) :: nFcx    !< Number of x-aligned faces/edges
                                       !< (1D objects)
        integer, intent(in) :: nFcy    !< Number of y-aligned faces/edges
                                       !< (1D objects)
        integer, intent(in) :: nVx     !< Number of all vertices/nodes
                                       !< (0D objects)

        gd%nCv = nCv
        gd%nFcx = nFcx
        gd%nFcy = nFcy
        gd%nVx = nVx

        gd%b2nx = nx
        gd%b2ny = ny

        allocate( gd%mapCvI(-1:nx+1, -1:ny+1) )
        gd%mapCvI = GRID_UNDEFINED
        allocate( gd%mapFcI(-1:nx+1, -1:ny+1, 0:3) )
        gd%mapFcI = GRID_UNDEFINED
        allocate( gd%mapVxI(-1:nx+1, -1:ny+1, 0:3) )
        gd%mapVxI = GRID_UNDEFINED

        allocate( gd%mapCvix(nCv), gd%mapCviy(nCv) )
        gd%mapCvix = GRID_UNDEFINED
        gd%mapCviy = GRID_UNDEFINED
        allocate( gd%mapFcix(nFcx+nFcy), gd%mapFciy(nFcx+nFcy),     &
            &     gd%mapFcIFace(nFcx+nFcy) )
        gd%mapFcix = GRID_UNDEFINED
        gd%mapFciy = GRID_UNDEFINED
        gd%mapFcIFace = GRID_UNDEFINED
        allocate( gd%mapVxix(nVx), gd%mapVxiy(nVx), gd%mapVxIVx(nVx) )
        gd%mapVxix = GRID_UNDEFINED
        gd%mapVxiy = GRID_UNDEFINED
        gd%mapVxIVx = GRID_UNDEFINED

        gd%nsv = 0
        allocate( gd%svix(MAX_SPECIAL_VERTICES),    &
            &   gd%sviy(MAX_SPECIAL_VERTICES), gd%svi(MAX_SPECIAL_VERTICES) )
        gd%svix = GRID_UNDEFINED
        gd%sviy = GRID_UNDEFINED
        gd%svi = GRID_UNDEFINED

        allocate( gd%mapCvixVx(nVx, 8), gd%mapCviyVx(nVx, 8))
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
        logical, intent(in) :: includeGhostCells        !< Include "fake" cells

        type(B2GridMap), intent(inout) :: gd    !< The grid mapping as computed
            !< by b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS

        !! internal variables
        integer :: ix    !< x-aligned (poloidal) cell index
        integer :: iy    !< y-aligned (radial) cell index
        integer :: ic    !< Cell index
        integer :: iFcx  !< x-aligned face index
        integer :: iFcy  !< y-aligned face index
        integer :: iVx   !< Vertex/node index
        integer :: i1    !< Iterator
        integer :: i2    !< Iterator
        integer :: nbix
        integer :: nbiy
        integer :: iFace !< Face/edge index
        integer :: iCorner
        integer :: index
        integer :: iPass

        !! There is an additional strip of "fake" cells (even faker than
        !! ghost cells) at the top and right to be able to also write out the
        !! ghost cells
        integer, dimension(-1:nx, -1:ny) :: Cvi          !< Control volume index
        integer, dimension(-1:nx, -1:ny, 0:3) :: Fcxi    !< x-aligned face index
            !< (only BOTTOM and TOP used)
        integer, dimension(-1:nx, -1:ny, 0:3) :: Fcyi    !< y-aligned face index
            !< (only LEFT and RIGHT used)
        integer, dimension(-1:nx+1, -1:ny+1, 0:3) :: Vxi !< vertex index

        logical :: CvNeeded( (nx + 2) * (ny + 2) )
        logical :: VxNeeded( (nx + 2) * (ny + 2) * 4)
        logical :: FcXNeeded( (nx + 2) * (ny + 2) * 2)
        logical :: FcYNeeded( (nx + 2) * (ny + 2) * 2)

        integer :: FcxiReduce( (nx + 2) * (ny + 2) * 2)
        integer :: FcyiReduce( (nx + 2) * (ny + 2) * 2)
        integer :: VxiReduce( (nx + 2) * (ny + 2) * 4)
        integer :: nsector( (nx + 2) * (ny + 2) * 4)

        logical :: cell_done

        !! list of identified special vertices
        !! vertex indices (ix, iy)
        integer, dimension(MAX_SPECIAL_VERTICES) :: svix !< Array of horizontal
            !< (x-aligned) positions of special vertices in B2 data structure
        integer, dimension(MAX_SPECIAL_VERTICES) :: sviy !< Array of vertical
            !< (y-aligned) positions of special vertices in B2 data structure
        integer, dimension(MAX_SPECIAL_VERTICES) :: svixAlias
        integer, dimension(MAX_SPECIAL_VERTICES) :: sviyAlias
        integer :: sVc  !< Number of special vertices
        integer :: sVcDuplicates    !< number of duplicate special vertices

        !! numbers of unique objects
        integer :: nCv  !< Number of all cells (2D objects)
        integer :: nFcx !< Number of x-aligned faces/edges (1D objects)
        integer :: nFcy !< Number of y-aligned faces/edges (1D objects)
        integer :: nVx  !< Number of all vertices/nodes (0D objects)

        character(*), parameter :: VERTEX_FILE = 'vertex_data.out'
        character(256) :: VERTEX_FILE_TEMP
        integer, parameter :: VERTEX_UNIT = 99
        logical :: vertexfileExists

        call logmsg( LOGDEBUG, "b2CreateMap: create map for a nx="  &
            &   //int2str(nx)//", ny="//int2str(ny)//" b2 grid" )

        Cvi = GRID_UNDEFINED
        Fcxi = GRID_UNDEFINED
        Fcyi = GRID_UNDEFINED
        Vxi = GRID_UNDEFINED

        !! set up initial lexicographic indices
        ic = 0   !! cells
        iFcx = 0 !! x-aligned faces
        iFcy = 0 !! y-aligned faces
        iVx = 0  !! vertices
        do ix = -1, nx
            do iy = -1, ny
                !! Cell
                ic = ic + 1
                Cvi( ix, iy ) = ic

                !! Faces: associate left and bottom face with every cell
                !! where possible
                iFcy = iFcy + 1
                Fcyi( ix, iy, LEFT ) = iFcy !! left face
                iFcx = iFcx + 1
                Fcxi( ix, iy, BOTTOM ) = iFcx !! bottom face

                !! Vertices: associate bottom left vertex with every cell
                !! where possible
                call find_Existing_Vertex_Index( ix, iy, VX_LOWERLEFT, index )
                if(index == GRID_UNDEFINED) then
                    iVx = iVx + 1
                    Vxi( ix, iy, VX_LOWERLEFT ) = iVx
                else
                    Vxi( ix, iy, VX_LOWERLEFT ) = index
                end if
            end do
        end do

        !! fill in index numbers for remaining faces
        do ix = -1, nx
            do iy = -1, ny

                !! Right face: left face of left neighbour
                call get_Neighbour( nx, ny, leftix, leftiy, rightix,    &
                    &   rightiy, topix, topiy, bottomix, bottomiy,      &
                    &   ix, iy, RIGHT, nbix, nbiy )

                if( is_Cell_In_Domain(nx, ny, nbix, nbiy) ) then
                    Fcyi( ix, iy, RIGHT ) = Fcyi( nbix, nbiy, LEFT )
                else
                    iFcy = iFcy + 1
                    Fcyi( ix, iy, RIGHT ) = iFcy
                end if

                !! Top face: bottom face of top neighbour
                !! also top-left vertex
                call get_Neighbour( nx, ny, leftix, leftiy, rightix,    &
                    &   rightiy, topix, topiy, bottomix, bottomiy,      &
                    & ix, iy, TOP, nbix, nbiy)

                if( is_Cell_In_Domain(nx, ny, nbix, nbiy) ) then
                    Fcxi( ix, iy, TOP ) = Fcxi( nbix, nbiy, BOTTOM )
                else
                    iFcx = iFcx + 1
                    Fcxi( ix, iy, TOP ) = iFcx
                end if
            end do
        end do

        call xertst (count(Fcxi(:,:,BOTTOM) == GRID_UNDEFINED) == 0,    &
            &  "b2CreateMap: there are unnumbered x-aligned faces" )
        call xertst (count(Fcxi(:,:,TOP) == GRID_UNDEFINED) == 0,       &
            &  "b2CreateMap: there are unnumbered x-aligned faces" )
        call xertst (count(Fcyi(:,:,LEFT) == GRID_UNDEFINED) == 0,      &
            &   "b2CreateMap: there are unnumbered y-aligned faces" )
        call xertst (count(Fcyi(:,:,RIGHT) == GRID_UNDEFINED) == 0,     &
            &   "b2CreateMap: there are unnumbered y-aligned faces" )

        !! Fill in vertex numbers for remaining vertices
        !! A vertex can be shared among 4 cells (possibly more for special
        !! vertices, but they are assumed to have a proper cell associated
        !! with them)
        do ix = -1, nx
            do iy = -1, ny

                do iCorner = VX_LOWERRIGHT, VX_UPPERRIGHT  ! 1, 3
                    call find_Existing_Vertex_Index( ix, iy, iCorner, index )
                    if( index /= GRID_UNDEFINED ) then
                        Vxi( ix, iy, iCorner ) = index
                    else
                        iVx = iVx + 1
                        Vxi( ix, iy, iCorner ) = iVx
                    end if
                end do

            end do
        end do

        call xertst( count( Vxi( -1:nx, -1:ny, 0:3 ) == GRID_UNDEFINED) == 0,   &
            & "b2CreateMap: there are unnumbered vertices" )

        !! Mark which cells, vertices and faces are needed in 1d lists
        CvNeeded = .false.
        FcXNeeded = .false.
        FcYNeeded = .false.
        VxNeeded = .false.
        do ix = -1, nx
            do iy = -1, ny
                if( .not. is_Unneeded_Cell( nx, ny, cflag,  &
                        &   includeGhostCells, ix, iy ) ) then
                    CvNeeded(Cvi(ix, iy)) = .true.
                    do iCorner = 0, 3
                        VxNeeded(Vxi( ix, iy, iCorner ) ) = .true.
                    end do
                    FcYNeeded( Fcyi( ix, iy, LEFT ) ) = .true.
                    FcYNeeded( Fcyi( ix, iy, RIGHT ) ) = .true.
                    FcXNeeded( Fcxi( ix, iy, BOTTOM ) ) = .true.
                    FcXNeeded( Fcxi( ix, iy, TOP ) ) = .true.
                end if
            end do
        end do

        !! search x-points.
        sVc = 0
        do ix = -1, nx
            do iy = -1, ny
                !! do not do this for unneeded cells, might be not initialized
                !! This is not perfect, but special vertices are usually inside
                !! the domain, and there it should be ok.
                if( .not. CvNeeded( Cvi( ix, iy ) ) ) cycle

                !! test whether vertex is special vertex
                if( is_Special_Vertex( ix, iy ) ) then
                    sVc = sVc + 1
                    sVix( sVc ) = ix
                    sViy( sVc ) = iy
                end if

            end do
        end do

        !! identify duplicate special vertices

        sVixAlias = 0
        sViyAlias = 0
        sVcDuplicates = 0

        do i1 = 1, sVc - 1

            ix = sVix( i1 )
            iy = sViy( i1 )

            !! has the point already been identified as unneeded?
            if( .not. VxNeeded( Vxi( ix, iy, VX_LOWERLEFT ) ) ) cycle

                !! compare with remaining vertices
                do i2 = i1 + 1, sVc

                !! has the point already been identified as unneeded?
                if( .not. VxNeeded( Vxi( sVix(i2), sViy(i2),   &
                    &   VX_LOWERLEFT ) ) ) cycle

                if( points_match( crx( ix, iy, 0 ), cry( ix, iy, 0 ),   &
                    &             crx( sVix(i2), sViy(i2), 0 ), &
                    &             cry( sVix(i2), sViy(i2), 0 ) ) ) then

                    !! The special vertex with index i2 is a duplicate of the
                    !! one with index i1.
                    !! Mark as unneeded and set up all references to this point
                    !! as aliased to the first one.
                    ! VxNeeded( Vxi(sVix(i2), sViy(i2), VX_LOWERLEFT) ) = .false.
                    ! where (Vxi == Vxi( sVix(i2), sViy(i2), VX_LOWERLEFT ))    &
                    !    &   Vxi = Vxi( ix, iy, VX_LOWERLEFT )

                    sVixAlias( i2 ) = ix
                    sViyAlias( i2 ) = iy
                    !! Bookkeeping for diagnostics
                    sVcDuplicates = sVcDuplicates + 1
                end if
            end do
        end do

        if( sVc - sVcDuplicates .eq. 1 ) then
            call logmsg( LOGDEBUG, "b2CreateMap: found "            &
                &   //int2str( sVc - sVcDuplicates )//              &
                &   " special vertex (x-point), "//int2str( sVc )   &
                &   //" including duplicates." )
        else
            call logmsg( LOGDEBUG, "b2CreateMap: found "            &
                &   //int2str( sVc - sVcDuplicates )                &
                &   //" special vertices (x-points), "              &
                &   //int2str( sVc )//" including duplicates." )
        end if

        !! number of unique cells
        !nCv = ( nx + 2 ) * ( ny + 2 ) - count( .not. isNeeded( Cvi ) )
        nCv = count( CvNeeded )
        !! number of unique faces (x and y direction)
        nFcx = count( FcXNeeded )
        nFcy = count( FcYNeeded )
        !! number of unique vertices
        nVx = count( VxNeeded )

        !! allocate the mapping structure
        call allocateB2GridMap( gd, nx, ny, nCv, nFcx, nFcy, nVx )

        !! build the mappings

        !! ...for cells
        ic = 0
        gd%mapCvI = GRID_UNDEFINED
        do ix = -1, nx
            do iy = -1, ny
                if( cvNeeded( Cvi( ix, iy ) ) ) then
                    ic = ic + 1
                    call xertst ( ic .le. nCv , &
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
        call xertst ( ic == nCv , 'b2CreateMap: found less cells than expected' )

        !! ...for x faces
        !! first set up numbering
        ic = 0
        FcxiReduce = GRID_UNDEFINED
        do i1 = 1, size( FcXNeeded )
            if( FcXNeeded( i1 ) ) then
                ic = ic + 1
                call xertst ( ic .le. nFcx , &
                    & 'b2CreateMap: found more x-aligned faces than expected' )
                FcxiReduce( i1 ) = ic
            end if
        end do
        !! sanity check for number of x faces
        call xertst ( ic == nFcx , &
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

                    if( FcXNeeded( Fcxi( ix, iy, iFace ) ) ) then
                        !if(is_Unneeded_Cell(nx, ny, cflag,  &
                        !   &   includeGhostCells, ix, iy)) cycle

                        !! Get the CPO index for this face
                        ic = FcxiReduce( Fcxi( ix, iy, iFace ) )
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
        ic = nFcx
        FcyiReduce = GRID_UNDEFINED
        do i1 = 1, size( FcYNeeded )
            if( FcYNeeded( i1 ) ) then
                ic = ic + 1
                FcyiReduce( i1 ) = ic
            end if
        end do
        call xertst( ic == nFcx + nFcy, &
            & 'b2CreateMap: found less y-aligned faces than expected' )

        do iPass = 1, 2
            do ix = -1, nx
                do iy = -1, ny
                    if(iPass == 1) then
                        iFace = LEFT
                    else
                        iFace = RIGHT
                    end if
                    if( FcYNeeded( Fcyi( ix, iy, iFace ) ) ) then
                        !! do not associate with unused cell
                        !if(is_Unneeded_Cell(nx, ny, cflag,  &
                        !   &   includeGhostCells, ix, iy)) cycle

                        !! Get the CPO index for this face
                        ic = FcyiReduce( Fcyi( ix, iy, iFace ) )
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
        VxiReduce = GRID_UNDEFINED
        do i1 = 1, size( VxNeeded )
            if( VxNeeded( i1 ) ) then
                ic = ic + 1
                VxiReduce( i1 ) = ic
            end if
        end do
        call xertst ( ic == nVx , &
            & 'b2CreateMap: did not find expected number of vertices' )

        gd%mapVxI = GRID_UNDEFINED
        gd%mapVxix = B2_GRID_UNDEFINED
        !! First looping over the corners makes sure that all vertices for which
        !! this is possible are associated with a lower-left vertex
        do iCorner = VX_LOWERLEFT, VX_UPPERRIGHT ! 0, 3
            do ix = -1, nx
                do iy = -1, ny
                    if( VxNeeded( vxi( ix, iy, iCorner ) ) ) then
                        ic = VxiReduce( Vxi( ix, iy, iCorner ) )
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
        do ic = 1, nVx
            if( gd%mapVxix( ic ) == B2_GRID_UNDEFINED ) then
                call logmsg( LOGWARNING, "b2CreateMap: vertex " &
                    &   //int2str( ic )//" not properly initialized (1)")
            end if
            if( gd%mapVxIVx( ic ) == B2_GRID_UNDEFINED ) then
                call logmsg( LOGWARNING, "b2CreateMap: vertex " &
                    &   //int2str( ic )//" not properly initialized (2)")
            end if
        end do

        call xertst( count( gd%mapVxix == B2_GRID_UNDEFINED ) == 0, &
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
            if( .not. VxNeeded( Vxi( ix, iy, VX_LOWERLEFT ) ) ) cycle

            !! if it is needed, copy it to the mapping structure
            gd%nsv = gd%nsv + 1
            gd%svix(gd%nsv) = ix
            gd%sviy(gd%nsv) = iy
            gd%svi(ic) = gd%mapVxI( gd%svix(ic), gd%sviy(ic), VX_LOWERLEFT )
        end do

        !! Now fill in missing vertex numbers where possible
        !! FIXME: do the same for the faces?
        Vxi = gd%mapVxI
        do ix = -1, nx
            do iy = -1, ny
                do iCorner = VX_LOWERRIGHT, VX_UPPERRIGHT  ! 1, 3
                    if( gd%mapVxI( ix, iy, iCorner) == GRID_UNDEFINED ) then
                        call find_Existing_Vertex_Index( ix, iy, iCorner,   &
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
            &   //int2str(gd%nCv)//" unique cells")
        call logmsg( LOGDEBUG, "b2CreateMap: map contains  "    &
            &   //int2str(gd%nFcx)//" unique x-aligned faces")
        call logmsg( LOGDEBUG, "b2CreateMap: map contains  "    &
            &   //int2str(gd%nFcy)//" unique y-aligned faces")
        call logmsg( LOGDEBUG, "b2CreateMap: map contains  "    &
            &   //int2str(gd%nVx)//" unique vertices")

        !! Find out how many and which control volumes touch any given vertex
        VERTEX_FILE_TEMP = trim(VERTEX_FILE)
        call find_file(VERTEX_FILE_TEMP,vertexfileExists)
        if (vertexfileExists) then
          ! read vertex file
          open(unit=VERTEX_UNIT, file=VERTEX_FILE_TEMP)
          read(VERTEX_UNIT,*) gd%mapCvixVx
          read(VERTEX_UNIT,*) gd%mapCviyVx
          close(VERTEX_UNIT)
        else
          nsector = 0
          do ix = -1, nx
            do iy = -1, ny
                do iCorner = VX_LOWERLEFT, VX_UPPERRIGHT
                    if( Vxi( ix, iy, iCorner ) == GRID_UNDEFINED ) cycle
                    nsector ( gd%mapVxI( ix, iy, iCorner) ) = 1
                    gd%mapCvixVx( gd%mapVxI( ix, iy, iCorner), 1 ) = ix
                    gd%mapCviyVx( gd%mapVxI( ix, iy, iCorner), 1 ) = iy
                    do i1 = -1, nx
                        do i2 = -1, ny
                            cell_done = i1 .eq. ix .and. i2 .eq. iy
                            do ic = VX_LOWERLEFT, VX_UPPERRIGHT
                                if( Vxi( i1, i2, ic ) == GRID_UNDEFINED ) cycle
                                if( cell_done ) cycle
                                if( points_match( crx( ix, iy, iCorner ),   &
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
#ifndef BUILDING_CARRE
          VERTEX_FILE_TEMP = trim("../"//VERTEX_FILE)
#else
          VERTEX_FILE_TEMP = trim(VERTEX_FILE)
#endif
          open(unit=VERTEX_UNIT, file=trim(VERTEX_FILE_TEMP))
          write(VERTEX_UNIT,*) gd%mapCvixVx
          write(VERTEX_UNIT,*) gd%mapCviyVx
          close(VERTEX_UNIT)
        endif
        return

contains

    !> For a given corner vertex of a cell, check whether in any connected
    !! cell a vertex index was already assigned to this vertex
    subroutine find_Existing_Vertex_Index( ix, iy, iCorner, index,  &
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
        index = Vxi( ix, iy, iCorner )
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

                call get_Neighbour( nx, ny, leftix, leftiy, rightix, &
                    &   rightiy, topix, topiy, bottomix, bottomiy,  &
                    &   nix, niy, iDir, nix2, niy2 )

                if(.not. is_Cell_In_Domain( nx, ny, nix2, niy2 )) exit
                if( lStopOnUnneededCells .and. is_Unneeded_Cell( nx, ny,    &
                    &   cflag, includeGhostCells, nix2, niy2) ) exit

                nix = nix2
                niy = niy2
                nICorner = VXCORNER_NEXTINDEX( iDir, nICorner )

                if( points_match( crx( ix, iy, iCorner ),                   &
                    &   cry( ix, iy, iCorner ), crx(nix, niy, nICorner ),   &
                    &   cry( nix, niy, nICorner ) ) ) then
                    index = Vxi( nix, niy, nICorner)
                    if( index /= GRID_UNDEFINED ) return
                end if
            end do
        end do

        index = GRID_UNDEFINED
    end subroutine find_Existing_Vertex_Index

    !> Test whether the vertex associated with the cell (ix, iy) is special,
    !! i.e. a x-point
    !! This steps around the vertex (left-bottom-right-top) and checks
    !! whether the resulting position is equal to the starting position
    logical function is_Special_Vertex( ix, iy )
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
        call get_Neighbour( nx, ny, leftix, leftiy, rightix, rightiy,   &
            &   topix, topiy, bottomix, bottomiy, x, y, LEFT, xn, yn )
        if( .not. is_Cell_In_Domain( nx, ny, xn, yn, extended =   &
            &   includeGhostCells )) then
            is_Special_Vertex = .false.
            return
        end if
        if( isGhostCell( cflag( xn, yn, CELLFLAG_TYPE ) ) ) then
            is_Special_Vertex = .false.
            return
        end if
        x = xn
        y = yn

        !! step bottom
        call get_Neighbour( nx, ny, leftix, leftiy, rightix, rightiy,   &
            &   topix, topiy, bottomix, bottomiy, x, y, BOTTOM, xn, yn )
        if( .not. is_Cell_In_Domain( nx, ny, xn, yn, extended =   &
            &   includeGhostCells ) ) then
            is_Special_Vertex = .false.
            return
        end if
        if( isGhostCell(cflag(xn,yn,CELLFLAG_TYPE)) ) then
            is_Special_Vertex = .false.
            return
        end if
        x = xn
        y = yn

        !! step right
        call get_Neighbour(nx, ny, leftix, leftiy, rightix, rightiy,    &
            &   topix,topiy,bottomix,bottomiy, x, y, RIGHT, xn, yn )
        if( .not. is_Cell_In_Domain( nx, ny, xn, yn, extended =   &
            &   includeGhostCells )) then
            is_Special_Vertex = .false.
            return
        end if
        if( isGhostCell(cflag(xn,yn,CELLFLAG_TYPE)) ) then
            is_Special_Vertex = .false.
            return
        end if
        x = xn
        y = yn

        !! step top
        call get_Neighbour(nx, ny, leftix, leftiy, rightix, rightiy,    &
            &   topix, topiy, bottomix, bottomiy, x, y, TOP, xn, yn )
        if( .not. is_Cell_In_Domain( nx, ny, xn, yn, extended =   &
            &   includeGhostCells)) then
            is_Special_Vertex = .false.
            return
        end if
        if( isGhostCell(cflag(xn,yn,CELLFLAG_TYPE)) ) then
            is_Special_Vertex = .false.
            return
        end if
        x = xn
        y = yn

        !! do we end up where we left?
        is_Special_Vertex = .not. ( ( x == ix ) .and. ( y == iy ) )

    end function is_Special_Vertex

  end subroutine b2CreateMap

    !> test whether cell (ix,iy) is actually used
    function is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, ix, iy )
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
        logical is_Unneeded_Cell

        !! Only cells inside the "normal" B2 domain can be needed
        !! (this catches fake cells and connectivity pointing outside the domain)
        is_Unneeded_Cell =  &
            &   .not. is_Cell_In_Domain(nx, ny, ix, iy, extended = .true.)
        if(is_Unneeded_Cell) return

        if(includeGhostCells) then
            is_Unneeded_Cell = isUnusedCell( cflag(ix,iy,CELLFLAG_TYPE) )
        else
            is_Unneeded_Cell = isUnusedCell( cflag(ix,iy,CELLFLAG_TYPE) )   &
                &   .or. isGhostCell( cflag(ix,iy,CELLFLAG_TYPE) )
        end if

        !! Classical treatment (without cflag, no extended grid) - for reference
        !if(.not. includeGhostCells) then
        !    is_Unneeded_Cell = &
        !        & ( leftix( ix, iy ) == -2 ) &
        !        & .or. ( rightix( ix, iy ) == ( nx + 1 ) ) &
        !        & .or. ( bottomiy( ix, iy ) == -2 ) &
        !        & .or. ( topiy( ix, iy ) == ( ny + 1 ) )
        !end if

    end function is_Unneeded_Cell


    !> Check whether the cell at position (ix,iy) is inside the 'classical'
    !! B2 grid.
    !! The default is to check whether it is inside  the extended grid
    !! (including the ghost cells).
    !! If the optional parameter extended is given, extended = .false. will
    !! check whether the position is inside the actual physical domain of the
    !! grid.
    function is_Cell_In_Domain( nx, ny, ix, iy, extended )
        integer, intent(in) :: nx   !< Specifies the number of interior cells
                                    !< along the first coordinate
        integer, intent(in) :: ny   !< Specifies the number of interior cells
                                    !< along the second coordinate
        integer, intent(in) :: ix   !< Specifies index of interior cell along
                                    !< the first coordinate
        integer, intent(in) :: iy   !< Specifies index of interior cell along
                                    !< the second coordinate
        logical, intent(in), optional :: extended
        logical :: is_Cell_In_Domain

        !! internal
        logical :: lExtended

        lExtended = .true.
        if( present( extended ) ) lExtended = extended

        if( lExtended ) then
            !! in extended domain (including ghost cells)?
            is_Cell_In_Domain = ( ix >= -1 ) .and. (ix <= nx) .and. &
                &   ( iy >= -1 ) .and. ( iy <= ny )
        else
            is_Cell_In_Domain = ( ix > -1 ) .and. (ix < nx) .and.   &
                &   ( iy > -1 ) .and. ( iy < ny )
        end if
    end function is_Cell_In_Domain


    !> Check whether the node associate with the cell at position (ix,iy) is
    !! included in the grid.
    function is_Node_In_Domain( nx, ny, ix, iy, extended )
        integer, intent(in) :: nx   !< Specifies the number of interior cells
                                    !< along the first coordinate
        integer, intent(in) :: ny   !< Specifies the number of interior cells
                                    !< along the second coordinate
        integer, intent(in) :: ix   !< Specifies index of interior cell along
                                    !< the first coordinate
        integer, intent(in) :: iy   !< Specifies index of interior cell along
                                    !< the second coordinate
        logical, intent(in), optional :: extended
        logical :: is_Node_In_Domain

        !! internal
        logical :: lExtended

        lExtended = .true.
        if( present( extended ) ) lExtended = extended

        if( lExtended ) then
            !! in extended domain (including ghost cells)?
            is_Node_In_Domain = ( ix >= -1 ) .and. (ix <= nx + 1) .and. &
                &   ( iy >= -1 ) .and. ( iy <= ny + 1 )
        else
            is_Node_In_Domain = ( ix > -1 ) .and. (ix < nx + 1) .and.   &
                &   ( iy > -1 ) .and. ( iy < ny + 1)
        end if
    end function is_Node_In_Domain

    !> extended neighbourhood mappings
    subroutine get_Neighbour(nx, ny, leftix, leftiy, rightix, rightiy,   &
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

    end subroutine get_Neighbour

end module b2mod_grid_mapping

!!!Local Variables:
!!! mode: f90
!!! End:
