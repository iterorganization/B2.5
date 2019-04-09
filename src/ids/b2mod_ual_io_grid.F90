!!-----------------------------------------------------------------------------
!! DOCUMENTATION:
!>      @section b2uw_ualio_grid_desc Description
!!      Module providing routines to set the B2 grid geometry, including grid
!!      subsets (subgrids), to ITM CPO or IMAS IDS grid description data
!!      structure.
!!
!!      The main two routines are:
!!          - b2_IMAS_Fill_Grid_Desc (for IMAS edge_profiles, edge_sources
!!            and edge_transport IDSs)
!!          - b2ITMFillGridDescription (for ITM edge CPO)
!!
!!      @subsection b2uw_ualio_grid_syx     Exceptional syntax explanation
!!      @code
!!          ! IGNORE    !! syntax used to ignore this module in list
!!                      !! dependency when compiling the code
!!      @endcode
!!
!!-----------------------------------------------------------------------------

module b2mod_ual_io_grid
    use b2mod_types , B2_R8 => R8, B2_R4 => R4
    use b2mod_constants , B2_PI => PI
#ifdef IMAS
    use ids_string        & ! IGNORE
     & , only : idsInt2str
    use ids_assert        & ! IGNORE
     & , only : assert
    use ids_grid_subgrid  & ! IGNORE
     & , only : getGridSubsetSize, getGridSubsetObject, findGridSubsetByName, &
     &          CreateGridSubsetForClass, CreateEmptyGridSubset, &
     &          CreateExplicitObjectListSingleSpace
    use ids_grid_object   & ! IGNORE
     & , only : ids_generic_grid_aos3_root, IDS_real, &
     &          ids_generic_grid_dynamic_grid_subset, &
     &          GRID_SUBSET_NODES, GRID_SUBSET_FACES, &
     &          GRID_SUBSET_X_ALIGNED_FACES, GRID_SUBSET_X_POINTS, &
     &          GRID_SUBSET_Y_ALIGNED_FACES, GRID_SUBSET_CELLS, &
     &          GridObject
    use ids_grid_structured & ! IGNORE
     & , only : GridWriteData, GridSetupStruct1dSpace
    use ids_grid_common   & ! IGNORE
     & , only : COORDTYPE_R, COORDTYPE_Z, COORDTYPE_PHI
#else
# ifdef ITM
    use itm_types , ITM_R8 => R8, ITM_R4 => R4 ! IGNORE
    use euITM_schemas ! IGNORE
    use itm_constants , pi => itm_pi ! IGNORE
    use itm_grid_object , only : GridObject ! IGNORE
    use itm_grid_common & ! IGNORE
     & , only : COORDTYPE_R, COORDTYPE_Z, COORDTYPE_PHI
    use itm_grid_objectlist & ! IGNORE
     & , only : createIndexListForRange
    use itm_grid_subgrid & ! IGNORE
     & , only : gridFindSubGridByName, gridSubGridSize, gridSubGridSize, &
     &          subGridGetObject
# endif
#endif
    use helper
    use logging , only: logmsg, LOGDEBUG
    use b2mod_connectivity , REMOVED_B2_R8 => R8
    use carre_constants
    use b2mod_cellhelper

    use b2mod_grid_mapping
    use b2mod_indirect

    implicit none

    !! Constants for use with the ITM grid description

    !! Space indices
    integer, parameter :: SPACE_POLOIDALPLANE = 1   !< Space indice
                                                    !< (polodal plane)
    integer, parameter :: SPACE_TOROIDALANGLE = 2   !< Space indice
                                                    !< (toroidal angle)
    integer, parameter :: SPACE_COUNT = SPACE_POLOIDALPLANE !< Space count
        !< set up:
        !< SPACE_COUNT = SPACE_POLOIDALPLANE: do only the poloidal plane space;
        !< SPACE_COUNT = SPACE_TOROIDALANGLE: will do the full 3d grid with two
        !< spaces
    !! Will have at maximum so many spaces
    integer, parameter :: SPACE_COUNT_MAX = 2

    !! Number of points in toroidal direction (only 1 makes sense here, this is
    !! for playing around)
    integer, parameter :: NNODES_TOROIDAL = 1   !< Number of points in toroidal
        !< direction (only 1 makes sense here, this is for playing around)

    logical, parameter :: TOROIDAL_PERIODIC = .false.   !< Flag controlling
        !< whether the toroidal space is set up periodic or not.
        !< If periodic, the last node is connect to the first by an edge.
        !< If not periodic, an additional node at 2 pi is added

    !! Object class tuples (ITM class definition)
    integer, dimension( SPACE_COUNT_MAX ), parameter :: CLASS_NODE =        &
        &   (/ 0, 0 /)  !< Object class tuple (ITM class definition): Node
    integer, dimension( SPACE_COUNT_MAX ), parameter :: CLASS_RZ_EDGE =     &
        &   (/ 1, 0 /)  !< Object class tuple (ITM class definition): Edge (R, Z)
    integer, dimension( SPACE_COUNT_MAX ), parameter :: CLASS_PHI_EDGE =    &
        &   (/ 0, 1 /)  !< Object class tuple (ITM class definition): Edge (phi)
    integer, dimension( SPACE_COUNT_MAX ), parameter ::   &
        &   CLASS_POLOIDALRADIAL_FACE = (/ 1, 1 /)  !< Object class tuple
        !< (ITM class definition): Poloidalradial face (in observed existing
        !< CPO cases this class refers to 2D cells)
    integer, dimension( SPACE_COUNT_MAX ), parameter ::   &
        &   CLASS_TOROIDAL_FACE = (/ 2, 0 /)    !< Object class tuple
        !< (ITM class definition): Toroidal face
    integer, dimension( SPACE_COUNT_MAX ), parameter :: CLASS_CELL = (/ 2, 1 /)
        !< Object class tuple (ITM class definition): Cell

    !! Object class tuples (one digit IDS classes, transformed from the above
    !! ITM classes )
    !! Primary IDS classes are:
    !!      Class 1 - nodes/vertices (0D objects)
    !!      Class 2 - edges/faces (1D objects)
    !!      Class 2 - 2D cells (2D objects)
    !! Fortran90 does not allow initialization of constants using SUM. This is
    !! permitted in newer Fortran 2003. Current workaround is to directly
    !! specify the primary IDS class constants

    ! integer, parameter :: IDS_CLASS_NODE = sum(CLASS_NODE) + 1
    integer, parameter :: IDS_CLASS_NODE = 1    !< Object class tuple (IMAS
        !< class definition): Node
    ! integer, parameter :: IDS_CLASS_RZ_EDGE = sum(CLASS_RZ_EDGE) + 1
    integer, parameter :: IDS_CLASS_RZ_EDGE = 2 !< Object class tuple (IMAS
        !< class definition): Edge
    ! integer, parameter :: IDS_CLASS_PHI_EDGE = sum(CLASS_PHI_EDGE) + 1
    integer, parameter :: IDS_CLASS_PHI_EDGE = 2    !< Object class tuple (IMAS
        !< class definition): Edge
    ! integer, parameter :: IDS_CLASS_POLOIDALRADIAL_FACE =   &
        ! &   sum(CLASS_POLOIDALRADIAL_FACE)
    integer, parameter :: IDS_CLASS_POLOIDALRADIAL_FACE = 2 !< Object class
        !< tuple (IMAS class definition): Edge
    ! integer, parameter :: IDS_CLASS_TOROIDAL_FACE =         &
        ! &   sum(CLASS_TOROIDAL_FACE)
    integer, parameter :: IDS_CLASS_TOROIDAL_FACE = 2 !< Object class tuple
        !< (IMAS !< class definition): Edge
    ! integer, parameter :: IDS_CLASS_CELL = sum(IDS_CLASS_CELL)
    integer, parameter :: IDS_CLASS_CELL = 3 !< Object class tuple (IMAS class
        !< definition): Cell (2D)

    !! Subgrid/Grid subset name constants

    integer, parameter :: B2_GENERIC_GSUBSET_COUNT = 6  !< Total number of
        !< generic grid subsets

    !! Generic grid subsets (all cells, all faces)
    !! Note: special grid subsets (given by region ids) do not have specific
    !! constants (see also b2mod_connectivity.f90)

    !! IMAS uses GGD grid subset identifier definitions defined in GSL
    !! (in ids_grid_common)
#ifdef ITM
    !! For ITM duplicates were made (old and new variable) in case of ITM code
    !! requiring old variables
    integer, parameter :: B2_GENERIC_SUBGRID_COUNT = 6  !< Total number of
        !< generic grid subsets

    integer, parameter :: B2_SUBGRID_UNSPECIFIED = 0    !< Unspecified grid
                                                        !< subset
    integer, parameter :: B2_GSUBSET_UNSPECIFIED = 0    !< Unspecified grid
                                                        !< subset
    integer, parameter :: B2_SUBGRID_NODES = 1  !< Grid subset containing all
        !< nodes (0D objects) belonging to associated space
    integer, parameter :: B2_GSUBSET_NODES = 1  !< Grid subset containing all
        !< nodes (0D objects) belonging to associated space
    integer, parameter :: B2_SUBGRID_FACES = 2  !< Grid subset containing all
        !< faces (first x-aligned, then y-aligned) belonging to the
        !< associated space (order given by grid map)
    integer, parameter :: B2_GSUBSET_FACES = 2  !< Grid subset containing all
        !< faces (first x-aligned, then y-aligned) belonging to the
        !< associated space (order given by grid map)
    integer, parameter :: B2_SUBGRID_FACES_X  = 3   !< Grid subset containing
        !< all x-aligned (poloidally) aligned faces belonging to the associated
        !< space (order given by grid map)
    integer, parameter :: B2_GSUBSET_X_ALIGNED_FACES = 3    !< Grid subset
        !< containing  all x-aligned (poloidally) aligned faces belonging to
        !< the associated  space (order given by grid map)
    integer, parameter :: B2_SUBGRID_FACES_Y = 4    !< Grid subset containing
        !< all y-aligned (radially) aligned faces belonging to the associated
        !< space (order given by grid map)
    integer, parameter :: B2_GSUBSET_Y_ALIGNED_FACES = 4 !< Grid subset
        !< containing  all y-aligned (radially) aligned faces belonging to
        !< the associated  space (order given by grid map)
    integer, parameter :: B2_SUBGRID_CELLS = 5 !< Grid subset containing all
        !< cells belonging to the associated space
    integer, parameter :: B2_GSUBSET_CELLS = 5 !< Grid subset containing all
        !< cells belonging to the associated space
    integer, parameter :: B2_SUBGRID_XPOINTS = 6    !< Grid subset containing
                                                    !< nodes defining x-points
    integer, parameter :: B2_GSUBSET_X_POINTS = 6   !< Grid subset containing
                                                    !< nodes defining x-points
#endif

contains

#ifdef IMAS

    !> Routine that fills in a grid description which is part of a IMAS IDS
    !! using the given grid data and prepared mappings
    subroutine b2_IMAS_Fill_Grid_Desc( gmap, ggd_grid, nx, ny, crx, cry,    &
        &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,       &
        &   bottomiy, nnreg, topcut, region, cflag, includeGhostCells, vol, &
        &   gs, qc )
        type(B2GridMap), intent(in) :: gmap !< The grid mapping as computed
            !< by b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS
#if IMAS_MINOR_VERSION < 15
        type(ids_generic_grid_dynamic), intent(out) :: ggd_grid !< Type of IDS
            !< data structure, designed for handling grid geometry data
#else
        type(ids_generic_grid_aos3_root), intent(out) :: ggd_grid !< Type of IDS
            !< data structure, designed for handling grid geometry data
#endif
        integer, intent(in) :: nx   !< Number of interior cells
            !< along the first coordinate (used to define size of grid arrays:
            !< (-1:nx, -1:ny)
        integer, intent(in) :: ny   !< Number of interior cells
            !< along the second coordinate (used to define size of grid arrays:
            !< (-1:nx, -1:ny)

        !! Output arguments
        real(R8), intent(in) :: crx( -1:nx, -1:ny, 0:3 )    !< Horizontal vertex
            !< coordinates of the four corners of the (ix, iy) cell
        real(R8), intent(in) :: cry( -1:nx, -1:ny, 0:3 )    !< Vertical vertex
            !< coordinates of the four corners of the (ix, iy) cell

        !! B2 connectivity array
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

        !! B2 region & cut information
        integer, intent(in) :: nnreg(0:2)
        integer, intent(in) :: topcut(:)
        integer, intent(in) :: region( -1:nx, -1:ny, 0:2 )
        !! Cell flags
        integer cflag( -1:nx, -1:ny,  CARREOUT_NCELLFLAGS ) !< Cell flag
        logical, intent(in) :: includeGhostCells    !< Include "fake" cells
        !! Optional B2 measure information
        real(R8), intent(in), optional :: vol( -1:nx, -1:ny, 0:4) !< Cell volume
        real(R8), intent(in), optional :: gs( -1:nx, -1:ny, 0:2)
        real(R8), intent(in), optional :: qc(-1:nx,-1:ny)   !< Cosine of the
            !< angle between flux line direction and left cell face

        !! Internal variables
        integer, parameter :: NDIM = 2  !< Dimension of the space

        call assert( present( gs ) .EQV. present( qc ) )

        !! Set GGD grid geometry
        call fill_In_Grid_Desc()
        !! Set grid subsets
        call fill_In_GridSubset_Desc()

contains
    !> Fill in the general grid description
    subroutine fill_In_Grid_Desc()
        !! Internal variables
        integer :: ivx  !< Vertex/node index
        integer :: ifc  !< Face/edge index
        integer :: icv  !< Cell index
        integer :: ix   !< x-aligned (poloidal) cell index
        integer :: iy   !< y-aligned (radial) cell index
        integer :: nix
        integer :: niy
        integer :: i    !< Iterator
        integer :: j    !< Iterator
        integer :: dir
        integer :: nfc  !< Number of all faces/edges (x + y aligned)
        integer :: geometryType  !< Geometry identifier index

        geometryType = geometryId(nnreg, periodic_bc, topcut)

        allocate( ggd_grid%identifier%name(1) )
        ggd_grid%identifier%name = geometryName(geometryType)
        ggd_grid%identifier%index = geometryType
        allocate( ggd_grid%identifier%description(1) )
        ggd_grid%identifier%description = geometryDescription(geometryType)

        allocate( ggd_grid%space( SPACE_COUNT ) )

        !! Coordinate types
        !! (dimension of space = NDIM = size( coordtype )

        allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%geometry_type%name(1) )
        ggd_grid%space(SPACE_POLOIDALPLANE)%geometry_type%name = 'Poloidal'

        !! Set the space coordinates, also defining the dimension of the space
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%coordinates_type(NDIM) )
        ggd_grid%space( SPACE_POLOIDALPLANE )%coordinates_type(1) = COORDTYPE_R
        ggd_grid%space( SPACE_POLOIDALPLANE )%coordinates_type(2) = COORDTYPE_Z

        !! Have two types of objects: 0d nodes, 1d edges, 2d cells
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension( NDIM + 1 ) )

        !! Allocate the number of objects of each type
        !! nodes
        !! gmap%nvx = ( nx+1 )*( ny+1 ) - 1)
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(1)%object( gmap%nvx ) )
        !! edges
        !! gmap%nfcx + gmap%nfcy = nx*( ny+1 ) + ( nx+1 ) * ny
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(2)%object( gmap%nfcx + gmap%nfcy ) )
        !! cells
        !! gmap%ncv = nx*ny
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(3)%object( gmap%ncv ) )

        !! Fill in node information
        do ivx = 1, gmap%nvx
            !! Allocate geometry leaf for each node
            allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(1)%object( ivx )%geometry(2) )
            !! Set geometry (R and Z coordinates) of each node
            !! The way of writing nodes data identically to the to IDS converted
            !! CPO 16151/1000 case, available on (written in 11 October 2017)
            !! ITER HPC custer in directory
            !! /home/ITER/penkod/public/imasdb/solps-iter/3/0
            !! or
            !! /home/ITER/kosl/public/imasdb/solps-iter/3/0
            !! while on GateWay Marconi
            !! /marconi_work/eufus_gw/work/g2penkod/imasdb/solps-iter/3/0
            !! or
            !! /marconi_work/eufus_gw/work/g2kosl/imasdb/solps-iter/3/0
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%   &
                &   object( ivx )%geometry(1) = crx(  gmap%mapVxix( ivx ),    &
                                                    & gmap%mapVxiy( ivx ),    &
                                                    & gmap%mapVxIVx( ivx ))
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%   &
                &   object( ivx )%geometry(2) = cry(  gmap%mapVxix( ivx ),    &
                                                    & gmap%mapVxiy( ivx ),    &
                                                    & gmap%mapVxIVx( ivx ))

            !! Set additional node index (REQUIRED!)
            allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(1)%object( ivx )%nodes(1))
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%   &
                &   object( ivx )%nodes(1) = ivx
        end do

#if 0
        !! Set geometry (R and Z coordinates) of each node
        !! (original b2mod. This node writing order does not work well with the
        !! code used to write 1D objects data (faces/edges). Thats why it is
        !! currently not used (#if 0) and the the method above is used instead)
        ivx = 0
        do iy = 0, ny-1
            do ix = 0, nx-1
                ivx = ivx + 1 !! Lower left corners
                ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                    &   objects_per_dimension(1)%object(ivx)%geometry(1) =  &
                    &   crx(ix,iy,0)
                ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                    &   objects_per_dimension(1)%object(ivx)%geometry(2) =  &
                    &   cry(ix,iy,0)
                if( ivx .eq. ( ( nx+1 )*( ny+1 ) - 1) ) exit
                if( ix.eq.nx-1) then
                    ivx = ivx + 1  !! Lower right corners
                    ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                        &   objects_per_dimension(1)%object(ivx)%geometry(1) =  &
                        &   crx( ix, iy, 1)
                    ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                        &   objects_per_dimension(1)%object(ivx)%geometry(2) =  &
                        &   cry( ix, iy, 1)
                    if( ivx .eq. ((nx+1)*(ny+1) - 1) ) exit
                end if
                if( iy.eq.ny-1) then
                    ivx = ivx + 1  !! Upper left corners
                    ggd_grid%space( SPACE_POLOIDALPLANE )%                      &
                        &   objects_per_dimension(1)%object( ivx )%geometry(1) =&
                        &   crx( ix, iy, 2 )
                    ggd_grid%space( SPACE_POLOIDALPLANE )%                      &
                        &   objects_per_dimension(1)%object( ivx )%geometry(2) =&
                        &   cry( ix, iy, 2 )
                    if( ivx .eq. ((nx+1)*(ny+1) - 1) ) exit
                    if( ix.eq.nx-1) then
                        ivx = ivx + 1  !! Upper right corners
                        ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                            &   objects_per_dimension(1)%object( ivx )% &
                            &   geometry(1) = crx( ix, iy, 3 )
                        ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                            &   objects_per_dimension(1)%object( ivx )% &
                            &   geometry(2) = cry( ix, iy, 3 )
                        if( ivx .eq. ( ( nx+1 )*( ny+1 ) - 1) ) exit
                    end if
                end if
            end do
        end do
#endif

        !! Fill in object definitions (i.e. what objects compose an object)
        !! 1D objects: faces/edges
        nfc = gmap%nfcx + gmap%nfcy !! Number of all faces/edges
        !! Allocate 1D objects
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(2)%object( gmap%nfcx + gmap%nfcy) )
        !! Each 1D object has two boundaries and two neighbours, in positive
        !! and negative coordinate direction, one on each side (for x-aligned
        !! faces: along flux surface, for y-aligned faces: orthogonal to flux
        !! surface)
        do ifc = 1, nfc
            !! Allocate and set all boundary & connectivity information to
            !! undefined
            allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(2)%object( ifc )%boundary(2) )
            !! Allocate list of 0D objects forming the 1D object
            !! Two 0D objects (vertices/nodes) form one 1D object (face/edge)
            allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(2)%object( ifc )%nodes(2) )
            do i = 1, 2
                !! Boundary to undefined
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(i)%index = B2_GRID_UNDEFINED
                !! Neighbours to undefined
                allocate(ggd_grid%space( SPACE_POLOIDALPLANE )% &
                    &   objects_per_dimension(2)%object( ifc )%boundary(i)% &
                    &   neighbours(2))
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(i)%neighbours(1) =   &
                    &   B2_GRID_UNDEFINED
                !! Nodes to undefined
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(i) = B2_GRID_UNDEFINED
            end do
        end do

        !! x-aligned faces
        do ifc = 1, gmap%nfcx
            !! get position of this face in the b2 grid
            ix = gmap%mapFcix( ifc )
            iy = gmap%mapFciy( ifc )
            !! get index of start vertex
            !! objdef dims: index of face, 1=start node, 1=one-dimensional
            !! object
            select case ( gmap%mapFcIFace( ifc ) )
            case( BOTTOM )
                !! start index: 1=start node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%index =     &
                    &   gmap%mapVxI( ix, iy, VX_LOWERLEFT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(1) = gmap%mapVxI( ix, iy, VX_LOWERLEFT )
                if (gmap%mapVxI( ix, iy, VX_LOWERLEFT ) ==  &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: BOTTOM face at " &
                        &   //idsInt2str(ix)//","//idsInt2str(iy)//         &
                        &   " has no start node")
                end if
                !! end vertex: 2=end node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(2)%index =     &
                    &   gmap%mapVxI( ix, iy, VX_LOWERRIGHT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(2) = gmap%mapVxI( ix, iy, VX_LOWERRIGHT )
                if (gmap%mapVxI( ix, iy, VX_LOWERRIGHT ) == &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: BOTTOM face at " &
                        &   //idsInt2str(ix)//","//idsInt2str(iy)//         &
                        &   " has no end node")
                end if
            case( TOP )
                !! start index: 1=start node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%index =     &
                    &   gmap%mapVxI( ix, iy, VX_UPPERLEFT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(1) = gmap%mapVxI( ix, iy, VX_UPPERLEFT )
                if (gmap%mapVxI( ix, iy, VX_UPPERLEFT ) ==  &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: TOP face at "    &
                        &   //idsInt2str(ix)//","//idsInt2str(iy)//         &
                        &   " has no start node")
                end if
                !! end vertex: 2=end node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(2)%index =     &
                    &   gmap%mapVxI( ix, iy, VX_UPPERRIGHT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(2) = gmap%mapVxI( ix, iy, VX_UPPERRIGHT )
                if (gmap%mapVxI( ix, iy, VX_UPPERRIGHT ) == &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: TOP face at "    &
                        &   //idsInt2str(ix)//","//idsInt2str(iy)//         &
                        &   " has no end node")
                end if
            end select

            !! Neighbour faces of this face
            !! Left neighbour: face continuing to the left of this face
            nix = leftix( ix, iy )
            niy = leftiy( ix, iy )
            if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%neighbours(1) =    &
                    &   gmap%mapFcI( nix, niy, gmap%mapFcIFace( ifc ) )
            end if
            !! Right neighbour: face continuing to the right of this face
            nix = rightix( ix, iy )
            niy = rightiy( ix, iy )
            if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(2)%neighbours(2) =     &
                    &   gmap%mapFcI( nix, niy, gmap%mapFcIFace( ifc ) )
            end if

            !! 1d object measure: face area
            if (present(gs)) then
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%measure = &
                    &   gs(ix, iy, ALIGNX)
            end if
        end do

        !! y-aligned faces
        do ifc = gmap%nfcx + 1, gmap % nfcx + gmap % nfcy
            !! get position of this face in the b2 grid
            ix = gmap % mapFcix( ifc )
            iy = gmap % mapFciy( ifc )
!!$            if (gmap%mapCvI(ix, iy) == B2_GRID_UNDEFINED) then
!!$                call logmsg(LOGWARNING,
!!$                "b2IMASFillGD: writing out faces for unused cell "//
!!$                idsint2str(ix)//","//idsint2str(iy))
!!$            end if

            select case ( gmap % mapFcIFace( ifc ) )
            case( LEFT )
                !! start index: 1=start node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%index =                       &
                    &   gmap%mapVxI( ix, iy, VX_LOWERLEFT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(1) = gmap%mapVxI( ix, iy, VX_LOWERLEFT )
                if (gmap%mapVxI( ix, iy, VX_LOWERLEFT ) ==  &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: LEFT face at "   &
                        &   //idsint2str(ix)//","//idsint2str(iy)//         &
                        &   " has no start node")
                end if
                !! end vertex: 2=end node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(2)%index =                       &
                    &   gmap%mapVxI( ix, iy, VX_UPPERLEFT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(2) = gmap%mapVxI( ix, iy, VX_UPPERLEFT )
                if (gmap%mapVxI( ix, iy, VX_UPPERLEFT ) ==  &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: LEFT face at "   &
                        &   //idsint2str(ix)//","//idsint2str(iy)//         &
                        &   " has no end node")
                end if

            case( RIGHT )
                !! start index: 1=start node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%index =                       &
                    &   gmap%mapVxI( ix, iy, VX_LOWERRIGHT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(1) = gmap%mapVxI( ix, iy, VX_LOWERRIGHT )
                if (gmap%mapVxI( ix, iy, VX_LOWERRIGHT ) == &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: RIGHT face at "  &
                        &   //idsint2str(ix)//","//idsint2str(iy)//         &
                        &   " has no start node")
                end if
                !! end vertex: 2=end node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(2)%index =                       &
                    &   gmap%mapVxI( ix, iy, VX_UPPERRIGHT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(2) = gmap%mapVxI( ix, iy, VX_UPPERRIGHT )
                if (gmap%mapVxI( ix, iy, VX_UPPERRIGHT ) == &
                    &    B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: RIGHT face at "  &
                        &   //idsint2str(ix)//","//idsint2str(iy)//         &
                        &   " has no end node")
                end if
            end select

            !! Neighbour faces of this face
            !! Bottom neighbour: face continuing to the bottom of this face
            nix = bottomix( ix, iy )
            niy = bottomiy( ix, iy )
            if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%neighbours(1) =    &
                    &   gmap%mapFcI( nix, niy, gmap%mapFcIFace( ifc ) )
            end if
            !! Top neighbour: face continuing to the top of this face
            nix = topix( ix, iy )
            niy = topiy( ix, iy )
            if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%neighbours(1) =    &
                    &   gmap%mapFcI( nix, niy, gmap%mapFcIFace( ifc ) )
            end if
            !! 1d object measure: edge length
            if (present(gs)) then
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%measure = &
                    &   gs(ix, iy, ALIGNY)*qc(ix, iy)
            end if
        end do

        !! Fill in object definitions (i.e. what objects compose an object)
        !! 2D objects: Cells
        ! write(0,*) "num_obj_2D: ", gmap%ncv
        !! Allocate 2D objects
        allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)% &
            &   object( gmap%ncv ) )

        !! Each 2D object has four boundaries
        !! Each boundary has one neighbour
        do icv = 1, gmap%ncv
            !! Allocate and set all boundary & connectivity information to
            !! undefined
            allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%  &
                &   objects_per_dimension(3)%object( icv )%boundary(4) )
            !! Allocate list of 0D objects forming the 2D object
            !! Four 0D objects (vertices/nodes) form one 2D object (2D cell)
            allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(3)%object( icv )%nodes(4) )
            !! Also store additional geometry information: position in
            !! computational space
            !! FIXME:   this should go into alternate geometry, which is not
            !!          available yet for grid objects
            allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%  &
                &   objects_per_dimension(3)%object( icv )%geometry(2) )
            do i = 1, 4
                !! Boundary to undefined
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                    &   object( icv )%boundary(i)%index = B2_GRID_UNDEFINED
                !! Neighbours to undefined
                allocate(ggd_grid%space( SPACE_POLOIDALPLANE )% &
                    &   objects_per_dimension(3)%object( icv )%boundary(i)% &
                    &   neighbours(2))
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                    &   object( icv )%boundary(i)%neighbours(1) =   &
                    &   B2_GRID_UNDEFINED
                !! Nodes to undefined
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                    &   object( icv )%nodes(i) = B2_GRID_UNDEFINED
                !! Geometry to undefined
                if (i < 3) then
                    ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                        &   object( icv )%geometry(i) = B2_GRID_UNDEFINED
                end if
            end do
        end do

        do icv = 1, gmap%ncv
            !! Set position in computational space
            ix = gmap%mapCvix( icv )
            iy = gmap%mapCviy( icv )
            ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                &   object( icv )%geometry(1) = ix
            ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                &   object( icv )%geometry(2) = iy

            !! Set faces composing the quadliateral in the list:
            !! left face (y-aligned)
            ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                &   object( icv )%boundary(1)%index = gmap%mapFcI( ix, iy, LEFT )
            !! bottom face (x-aligned
            ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                &   object( icv )%boundary(2)%index = gmap%mapFcI( ix, iy, BOTTOM )
            !! right face (y-aligned)
            ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                &   object( icv )%boundary(3)%index = gmap%mapFcI( ix, iy, RIGHT )
            !! top face (x-aligned)
            ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                &   object( icv )%boundary(4)%index = gmap%mapFcI( ix, iy, TOP )
            do dir = LEFT, TOP
                call get_Neighbour(nx, ny, leftix, leftiy, rightix, rightiy,     &
                    &   topix, topiy, bottomix, bottomiy, ix, iy, dir, nix, niy)
                if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells,    &
                    &   nix, niy ) ) then
                    ggd_grid%space(SPACE_POLOIDALPLANE)%            &
                        &   objects_per_dimension(3)%object (icv )% &
                        &   boundary( dir + 1 )%neighbours(1) =     &
                        &   gmap%mapCvI( nix, niy )
                end if
            end do
            !! 2d object measure: cell area
            if (present(vol)) then
                ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                    &   object( icv )%measure = vol(ix, iy, 1)
            end if
        end do

        !! Set nodes list, composing the 2D objects - Cells, using a subroutine
        call set_Cells_Conn_Array_Nodes(ggd_grid)

#if 0
        !! TODO
        !! Fill in x-point indices
        !! In edge_profiles no node for data on x-points was found. There is
        !! hovewer one in equilibrium%boundary%x_point
        allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%xpoints( gmap%nsv ) )
        ggd_grid%space(SPACE_POLOIDALPLANE)%xpoints = gmap%svi(1:gmap%nsv)
#endif

        !! If requested, add a second space for the toroidal angle
        if (SPACE_COUNT == SPACE_TOROIDALANGLE) then
            if ( TOROIDAL_PERIODIC ) then
                call gridSetupStruct1dSpace(                                    &
                    &   ggd_grid%space(SPACE_TOROIDALANGLE), COORDTYPE_PHI,     &
                    &   (/                                                      &
                    &   ( ( 2*B2_PI/NNODES_TOROIDAL )*i, i=0, NNODES_TOROIDAL-1 )   &
                    &   /),                                                     &
                    &   .true. ) !! periodic = .true.
            else
                call gridSetupStruct1dSpace(                                    &
                    &   ggd_grid%space( SPACE_TOROIDALANGLE ), COORDTYPE_PHI,   &
                    &   (/                                                      &
                    &   ( ( 2*B2_PI/NNODES_TOROIDAL )*i, i=0, NNODES_TOROIDAL ) &
                    &   /),                                                     &
                    &   .false. ) !! periodic = .false.
            end if
        end if

    end subroutine fill_In_Grid_Desc

    !> Set connectivity array for cells by defining nodes that form each cell
    subroutine set_Cells_Conn_Array_Nodes(ggd_grid)
#if IMAS_MINOR_VERSION < 15
        type(ids_generic_grid_dynamic), intent(inout) :: ggd_grid !< Type of IDS
            !< data structure, designed for handling grid geometry data
#else
        type(ids_generic_grid_aos3_root), intent(inout) :: ggd_grid !< Type of IDS
            !< data structure, designed for handling grid geometry data
#endif
        !! Internal variables
        integer, allocatable :: objects2Darray(:,:)
        integer :: num_nodes_2D     !< Total number of nodes forming one cell
        integer :: num_boundary_2D  !< Total number of boundary edges forming
            !< one cell
        integer :: node1
        integer :: node2
        integer, allocatable :: node_idx(:)
        integer, allocatable :: free_edge(:)
        integer :: edge_idx
        integer :: last_idx
        integer :: icv  !< Cell index
        integer :: ix   !< x-aligned (poloidal) cell index
        integer :: iy   !< y-aligned (radial) cell index
        integer :: m    !< Iterator
        integer :: loop_count   !< Loop counter


        !! Get the list of 0D objects forming the 2D objects and wirte it to IDS
        allocate( objects2Darray( gmap%ncv, 4 ) )
        ! objects2Darray = numpy.array([], dtype='int')

        ! ids_dim_2D.object.resize(ncv)

        ! Already done
        allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%  &
            &   objects_per_dimension(3)%object( gmap%ncv ) )

        !! Get 2D objects geometry (nodes) data from CPO, sort them into more
        !! orderly form  and put them into the IDS
        num_nodes_2D    = 4
        num_boundary_2D = 4
        allocate( node_idx(4) )
        allocate( free_edge(3) )

        do icv = 1, gmap%ncv
            ix = gmap%mapCvix( icv )
            iy = gmap%mapCviy( icv )
            allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%  &
                &   objects_per_dimension(3)%object(icv)%nodes(num_nodes_2D) )

            edge_idx = gmap%mapFcI( ix, iy, LEFT )
            free_edge(1) = gmap%mapFcI( ix, iy, BOTTOM )
            free_edge(2) = gmap%mapFcI( ix, iy, RIGHT )
            free_edge(3) = gmap%mapFcI( ix, iy, TOP )

            node_idx(1) = ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(2)%object( edge_idx )%    &
                &   boundary(1)%index
            last_idx = 2
            node_idx(last_idx) = ggd_grid%space( SPACE_POLOIDALPLANE )% &
                &   objects_per_dimension(2)%object( edge_idx )%        &
                &   boundary(2)%index

            do loop_count = 1, 4
                if (last_idx < 4) then
                    do m = 1, 3
                        edge_idx = free_edge(m)
                        if (edge_idx < 1) then
                            continue
                        end if
                        node1 = ggd_grid%space( SPACE_POLOIDALPLANE )%   &
                            &   objects_per_dimension(2)%object( edge_idx )% &
                            &   boundary(1)%index
                        node2 = ggd_grid%space( SPACE_POLOIDALPLANE )%   &
                            &   objects_per_dimension(2)%object( edge_idx )% &
                            &   boundary(2)%index

                        if (node_idx(last_idx) == node1) then
                            free_edge(m) = -1
                            last_idx = last_idx + 1
                            node_idx(last_idx)= node2
                            exit
                        end if
                        if (node_idx(last_idx) == node2) then
                            free_edge(m) = -1
                            last_idx = last_idx + 1
                            node_idx(last_idx) = node1
                            exit
                        end if
                    end do
                end if
            end do

            !! Set nodes list, composing the 2D objects - Cells
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                &   object( icv )%nodes(1) = node_idx(1)
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                &   object( icv )%nodes(2) = node_idx(4)
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                &   object( icv )%nodes(3) = node_idx(3)
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                &   object( icv )%nodes(4) = node_idx(2)
        end do

    end subroutine set_Cells_Conn_Array_Nodes

    !> Define grid subsets
    subroutine fill_In_GridSubset_Desc
        !! Internal variables
        integer :: geoId
        integer :: iRegion
        integer :: GSubsetCount
        integer :: iType
        integer :: nGSubset !< Total number of grid subsets
        integer :: xIn
        integer :: yIn
        integer :: xOut
        integer :: yOut
        integer :: iCoreGS !< Core grid subset ID
        integer :: cls(SPACE_COUNT_MAX)
        integer, allocatable :: xpoints(:,:)
        integer, allocatable :: indexList1d(:)
        integer, dimension(:,:), allocatable :: indexList2d
        integer :: i    !< Iterator

        geoId = geometryId(nnreg, periodic_bc, topcut)

        !! Figure out total number of grid subsets
        !! Do generic grid subsets + grid subsets
        nGSubset = B2_GENERIC_GSUBSET_COUNT + regionCountTotal(geoId)
        !! Inner/outer midplane grid subsets
        nGSubset = nGSubset + 2

        call logmsg( LOGDEBUG, "b2_IMAS_Fill_Grid_Desc: expecting total of " &
            &//idsInt2str(nGSubset)//" grid subsets" )
        allocate( ggd_grid%grid_subset( nGSubset ) )

        !! Set up generic grid subsets
        !! GRID_SUBSET_NODES: all nodes, one implicit object list
        !! The commented codes use CLASS_NODE, CLASS_POLOIDALRADIAL_FACE and
        !! CLASS_CELL parameters, using CPO class definitions, which might not
        !! work with the current IMAS GSL

        call createGridSubsetForClass(  ggd_grid,                   &
            &   ggd_grid%grid_subset( GRID_SUBSET_NODES ),          &
            &   IDS_CLASS_NODE, 1, GRID_SUBSET_NODES, "Nodes",      &
            &   "All nodes (0D objects) in the domain."  )

        !! GRID_SUBSET_FACES: all faces, one implicit object list
        call createGridSubsetForClass(  ggd_grid,                       &
            &   ggd_grid%grid_subset( GRID_SUBSET_FACES ),              &
            &   IDS_CLASS_POLOIDALRADIAL_FACE, 1, GRID_SUBSET_FACES,    &
            &   "Faces", "All faces (1D objects) in the domain."  )

        !! GRID_SUBSET_X_ALIGNED_FACES: x-aligned faces. One implicit object
        !! list, range over x faces
        !! Create grid subset with one object list
        call createEmptyGridSubset(                                     &
            &   ggd_grid%grid_subset( GRID_SUBSET_X_ALIGNED_FACES ),    &
            &   GRID_SUBSET_X_ALIGNED_FACES, 'x-aligned faces' )
        !! Initialize implicit object list for faces (class (/2/) )
        allocate(indexList1d(gmap%nfcx))
        indexList1d = (/ (i, i = 1, gmap%nfcx) /)
        call createExplicitObjectListSingleSpace( ggd_grid,            &
            &   ggd_grid%grid_subset( GRID_SUBSET_X_ALIGNED_FACES),    &
            &   IDS_CLASS_POLOIDALRADIAL_FACE, indexList1d,            &
            &   IDS_CLASS_POLOIDALRADIAL_FACE, 1)

        if ( SPACE_COUNT == SPACE_TOROIDALANGLE ) then
            deallocate(indexList1d)
            allocate(indexList1d(1))
            indexList1d = (/ 1 /)
            call createExplicitObjectListSingleSpace( ggd_grid,         &
                &   ggd_grid%grid_subset( GRID_SUBSET_X_ALIGNED_FACES), &
                &   IDS_CLASS_POLOIDALRADIAL_FACE, indexList1d,         &
                &   IDS_CLASS_POLOIDALRADIAL_FACE, 1)
        end if

        !! GRID_SUBSET_Y_ALIGNED_FACES: y-aligned faces. One implicit object
        !! list, range over y faces
        !! Create grid subset with one object list
        call createEmptyGridSubset(                                     &
            &   ggd_grid%grid_subset( GRID_SUBSET_Y_ALIGNED_FACES ),    &
            &   GRID_SUBSET_Y_ALIGNED_FACES, 'y-aligned faces' )
        !! Initialize implicit object list for faces (class (/2/) )
        deallocate(indexList1d)
        allocate(indexList1d(gmap%nfcy))
        indexList1d = (/ (i, i = gmap%nfcx + 1, gmap%nfcx + gmap%nfcy) /)
        call createExplicitObjectListSingleSpace( ggd_grid,             &
            &   ggd_grid%grid_subset( GRID_SUBSET_Y_ALIGNED_FACES),     &
            &   IDS_CLASS_POLOIDALRADIAL_FACE,                          &
            &   indexList1d, IDS_CLASS_POLOIDALRADIAL_FACE, 1)

        if ( SPACE_COUNT == SPACE_TOROIDALANGLE ) then
        deallocate(indexList1d)
        allocate(indexList1d(1))
        indexList1d = (/ 1 /)
        call createExplicitObjectListSingleSpace( ggd_grid,           &
            &   ggd_grid%grid_subset( GRID_SUBSET_Y_ALIGNED_FACES),   &
            &   IDS_CLASS_POLOIDALRADIAL_FACE, indexList1d,           &
            &   IDS_CLASS_POLOIDALRADIAL_FACE, 1)
        end if

        !! GRID_SUBSET_CELLS: all 2d cells, one implicit object list
        call createGridSubsetForClass(  ggd_grid,               &
            &   ggd_grid%grid_subset( GRID_SUBSET_CELLS ),      &
            &   IDS_CLASS_CELL, 1, GRID_SUBSET_CELLS, "Cells",  &
            &   "All faces (1D objects) in the domain."  )

        !! Grid subset of all x-points (in one poloidal plane at toroidal
        !! index 1)
        !! Assemble object descriptor for x-points
        allocate( xpoints(gmap%nsv, SPACE_COUNT) )
        xpoints = 1
        xpoints(:, SPACE_POLOIDALPLANE) = gmap%svi(1:gmap%nsv)
        !! Create grid subset with one object list
        call createEmptyGridSubset(                                     &
            &   ggd_grid%grid_subset( GRID_SUBSET_X_POINTS ),           &
            &   GRID_SUBSET_X_POINTS, 'x-points' )
        !! Initialize explicit object list for faces (class (/1/) )
        !! TODO: xpoints(:, 1 ) -> taking values for first space only. Set for
        !! all spaces.
        call createExplicitObjectListSingleSpace( ggd_grid,         &
                &   ggd_grid%grid_subset( GRID_SUBSET_X_POINTS),    &
                &   IDS_CLASS_NODE , xpoints(:, 1), IDS_CLASS_NODE, 1)

        !! Set up specific grid subset by collection faces for regions

        !! Start counting from end of generic grid subset
        GSubsetCount = B2_GENERIC_GSUBSET_COUNT

        !! Cell + face grid subset
        do iType = REGIONTYPE_CELL, REGIONTYPE_YFACE

            select case(iType)
            case( REGIONTYPE_CELL )
                cls = CLASS_CELL
            case( REGIONTYPE_YFACE, REGIONTYPE_XFACE )
                cls = CLASS_POLOIDALRADIAL_FACE
            end select

            do iRegion = 1, regionCount(geoId, iType)
                GSubsetCount = GSubsetCount + 1

                call logmsg( LOGDEBUG,                                      &
                    &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//       &
                    &   idsInt2str(GSubsetCount)//                          &
                    &   " for iType "//idsInt2str( iType )//                &
                    &   ", iRegion "//idsInt2str( iRegion )//": "//         &
                    &   regionName(geoId, iType, iRegion) )

                !! Create grid subset with one object list
                call createEmptyGridSubset(                     &
                    &   ggd_grid%grid_subset( GSubsetCount ),   &
                    &   GSubsetCount, regionName( geoId, iType, iRegion ) )

                !! Get explicit object list of the grid subset using
                !! subroutine collectIndexListForRegionSubroutine
                !! (function collectIndexListForRegion transferred to subroutine,
                !! as array of certain dimension is required as an output)
                call collectIndexListForRegionSubroutine( gmap, region,   &
                    &   iType, iRegion, indexList2d )

                !! Initialize explicit object list for grid subset
                !! TODO: Currently taking object indices only from space 1
                !!      ( %space(1) ). Set
                !!       to search all spaces
                call createExplicitObjectListSingleSpace( ggd_grid,           &
                    &   ggd_grid%grid_subset( GSubsetCount ), sum( cls ) - 1, &
                    &   indexList2d(:,1), sum(cls), 1)
            end do
        end do

        deallocate(indexList2d)
        !! Add midplane node grid subsets
        !! Find the core boundary grid subset by looking for its name as
        !! defined in b2mod_connectivity
        iCoreGS = findGridSubsetByName(ggd_grid, "Core boundary")
        !! For double null, we need the outer half of the core boundary
        if (iCoreGS == B2_GRID_UNDEFINED) then
            iCoreGS = findGridSubsetByName(ggd_grid, "Outer core boundary")
        end if
        if (iCoreGS == B2_GRID_UNDEFINED) stop "fill_In_GridSubset_Desc: "// &
            & "did not find core boundary grid subset for assembling " // &
            & " outer midplane grid subset"

        !! Figure out starting points for inner and outer midplane on core
        !! boundary
        call find_Midplane_Cells(ggd_grid%grid_subset( iCoreGS ), gmap, crx,  &
            &   xIn, yIn, xOut, yOut)

        GSubsetCount = GSubsetCount + 1
        !! Create grid subset with one object list
        call createEmptyGridSubset(                     &
            &   ggd_grid%grid_subset( GSubsetCount ),   &
            &       GSubsetCount, "Inner Midplane" )

        !! Get explicit object list of the grid subset using
        !! subroutine collectRadialVertexIndexListSubroutine
        !! (function collectRadialVertexIndexList transferred to subroutine,
        !! as array of certain dimension is required as an output)
        call collectRadialVertexIndexListSubroutine(gmap, cflag, xIn, yIn,  &
            &   topix, topiy, indexList2d)

        !! Initialize explicit object list for grid subset
        !! TODO: Currently taking object indices only from space 1
        !!      ( %space(1) ). Set
        !!       to search all spaces
        call createExplicitObjectListSingleSpace( ggd_grid,                 &
            &   ggd_grid%grid_subset( GSubsetCount ), IDS_CLASS_NODE - 1,   &
            &   indexList2d(:,1), IDS_CLASS_NODE, 1)

        GSubsetCount = GSubsetCount + 1

        !! Create grid subset with one object list
        call createEmptyGridSubset(                     &
            &   ggd_grid%grid_subset( GSubsetCount ),   &
            &       GSubsetCount, "Outer Midplane" )

        !! Get explicit object list of the grid subset using
        !! subroutine collectRadialVertexIndexListSubroutine
        !! (function collectRadialVertexIndexList transferred to subroutine,
        !! as array of certain dimension is required as an output)
        call collectRadialVertexIndexListSubroutine(gmap, cflag, xOut, yOut,  &
            &   topix, topiy, indexList2d)

        !! Initialize explicit object list for grid subset
        !! TODO: Currently taking object indices only from space 1
        !!      ( %space(1) ). Set
        !!       to search all spaces
        call createExplicitObjectListSingleSpace( ggd_grid,                 &
            &   ggd_grid%grid_subset( GSubsetCount ), IDS_CLASS_NODE - 1,   &
            &   indexList2d(:,1), IDS_CLASS_NODE, 1)

        call logmsg( LOGDEBUG, "b2_IMAS_Fill_Grid_Desc: wrote total of "    &
            &   //idsInt2str(GSubsetCount)//" grid subsets (expected was "  &
            &   //idsInt2str(size(ggd_grid%grid_subset))//")" )

        call assert( GSubsetCount == size(ggd_grid%grid_subset) )
    end subroutine fill_In_GridSubset_Desc

    end subroutine b2_IMAS_Fill_Grid_Desc

    !> Figure out starting cells for inner and outer midplane on core boundary
    !! by finding the points on the core boundary with minimum and maximum r
    !! positions
    subroutine find_Midplane_Cells( GridSubset, gmap, crx, xIn, yIn,  &
            &   xOut, yOut )
        type(ids_generic_grid_dynamic_grid_subset), intent(in) :: GridSubset
            !< Type of IDS data structure, designed for handling grid subset
            !< definitions
        type(B2GridMap), intent(in) :: gmap !< The grid mapping as computed by
            !< b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS
        real(IDS_real), intent(in) :: crx( -1:gmap%b2nx, -1:gmap%b2ny, 0:3 )
            !< Horizontal (x/radial) vertex coordinates of the four corners
            !< of the (ix, iy) cell
        integer, intent(out) :: xIn
        integer, intent(out) :: yIn
        integer, intent(out) :: xOut
        integer, intent(out) :: yOut

        !! Internal variables
        real(IDS_real) :: rMin
        real(IDS_real) :: rMax
        type(GridObject) :: obj !< GGD grid object
        integer :: ix   !< x-aligned cell index
        integer :: iy   !< y-aligned cell index
        integer :: iObj !< Object index

        rMin = huge(rMin)
        rMax = -huge(rMax)

        xIn = huge(xIn)
        xOut = huge(xOut)

        !! Loop over all faces in core boundary grid subset
        do iObj = 1, getGridSubsetSize(GridSubset)
            obj = getGridSubsetObject(GridSubset, iObj)
            !! Expect a face
            call xertst( all( obj%cls( 1:SPACE_COUNT ) ==               &
                &   IDS_CLASS_POLOIDALRADIAL_FACE ), &
                & "b2mod_ual_io_grid find_Midplane_Cells: assertion failure." )
            !! ...which is aligned along the x-direction
            call xertst( gmap%mapFcIFace( obj%ind( SPACE_POLOIDALPLANE ) ) ==   &
                &   BOTTOM, &
                & "b2mod_ual_io_grid find_Midplane_Cells: assertion failure." )
            ix = gmap % mapFcix( obj%ind( SPACE_POLOIDALPLANE ) )
            iy = gmap % mapFciy( obj%ind( SPACE_POLOIDALPLANE ) )

            !! We want the vertex associated with the cell at ix, iy, which is number 0
            if ( crx(ix, iy, 0) < rMin ) then
                rMin = crx(ix, iy, 0)
                xIn = ix
                yIn = iy
            end if
            if ( crx(ix, iy, 0) > rMax ) then
                rMax = crx(ix, iy, 0)
                xOut = ix
                yOut = iy
            end if
        end do

        call xertst( xIn /= huge( xIn ),    &
            &   "find_Midplane_Cells: did not find inner midplane position")
        call xertst( xOut /= huge( xOut ),  &
            &   "find_Midplane_Cells: did not find outer midplane position")
    end subroutine find_Midplane_Cells

#else
#ifdef ITM

  !> Routine that fills in a grid description which is part of a CPO
  !> using the given grid data and prepared mappings
  subroutine b2ITMFillGridDescription( gmap,itmgrid, &
      & nx,ny,crx,cry, &
      & leftix,leftiy,rightix,rightiy, &
      & topix,topiy,bottomix,bottomiy,&
      & nnreg,topcut,region,cflag,includeGhostCells,vol,gs,qc )

    type(B2GridMap), intent(in) :: gmap
    type(type_complexgrid), intent(out) :: itmgrid

    ! Size of grid arrays: (-1:nx, -1:ny)
    integer, intent(in) :: nx, ny
    !   .. output arguments
    ! vertex coordinates
    real (ITM_R8), intent(in) :: &
        & crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3)
    ! B2 connectivity array
    integer, intent(in) :: &
        & leftix(-1:nx,-1:ny),leftiy(-1:nx,-1:ny),&
        & rightix(-1:nx,-1:ny),rightiy(-1:nx,-1:ny),&
        & topix(-1:nx,-1:ny),topiy(-1:nx,-1:ny),&
        & bottomix(-1:nx,-1:ny),bottomiy(-1:nx,-1:ny)
    ! B2 region & cut information
    integer, intent(in) :: &
        & nnreg(0:2), topcut(:), &
        & region(-1:nx,-1:gmap%b2ny,0:2)
    ! Cell flags
    integer cflag(-1:nx,-1:ny, CARREOUT_NCELLFLAGS)
    logical, intent(in) :: includeGhostCells
    ! Optional B2 measure information
    real(ITM_R8), intent(in), optional :: vol(-1:nx,-1:ny,0:4), gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny)

    ! internal
    integer, parameter :: NDIM = 2

    call assert( present(gs) .EQV. present(qc) )

    call fill_In_Grid_Desc()
    call fillInSubGridDescription()

  contains

    !! Part 1: fill in grid description
    subroutine fill_In_Grid_Desc()

      !! internal
      integer :: ivx, ifc, icv, ix, iy, nix, niy, i, dir

      allocate( itmgrid % spaces(SPACE_COUNT) )

      !! Coordinate types
      ! (dimension of space = NDIM = size( coordtype )
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % coordtype(NDIM, 1) )
      itmgrid % spaces(SPACE_POLOIDALPLANE) % coordtype(:, 1) &
           & = (/ COORDTYPE_R, COORDTYPE_Z /)

      !! Have two types of objects: 1d edges, 2d cells
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(NDIM + 1) )

      !! Fill in node information
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(1) % geo(gmap%nvx, NDIM, 1, 1) )
      do ivx = 1, gmap % nvx
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(1) % geo(ivx, 1, 1, 1) = &
              & crx( gmap % mapVxix( ivx ), gmap % mapVxiy( ivx ), gmap % mapVxIVx( ivx ) )
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(1) % geo(ivx, 2, 1, 1) = &
              & cry( gmap % mapVxix( ivx ), gmap % mapVxiy( ivx ), gmap % mapVxIVx( ivx ) )
      end do

      !! Fill in object definitions (i.e. what objects compose an object)

      !! 1d objects: faces
      !! ...have two boundaries
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( gmap%nfcx + gmap%nfcy, 2) )
      !! ...have two neighbours, in positive and negative coordinate direction, one on each side
      !! (for x-aligned faces: along flux surface, for y-aligned faces: orthogonal to flux surface)
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour(gmap%nfcx + gmap%nfcy, 2, 1) )
      !! 1d object measure: face area
      if (present(gs)) allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % measure( gmap % nfcx + gmap % nfcy, 1 ) )
      !! first set all boundary & connectivity information to undefined
      itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary = GRID_UNDEFINED
      itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour = GRID_UNDEFINED
      !! x-aligned faces
      do ifc = 1, gmap % nfcx
          !! get position of this face in the b2 grid
          ix = gmap % mapFcix( ifc )
          iy = gmap % mapFciy( ifc )
          !! get index of start vertex
          !! objdef dims: index of face, 1=start node, 1=one-dimensional object
          select case ( gmap % mapFcIFace( ifc ) )
          case( BOTTOM )
             !! start index: 1=start node
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 ) = gmap % mapVxI( ix, iy, VX_LOWERLEFT )
             if (gmap % mapVxI( ix, iy, VX_LOWERLEFT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2ITMFillGD: BOTTOM face at "//int2str(ix)//","//int2str(iy)//" has no start node")
             end if
             !! end vertex: 2=end node
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 2 ) = gmap % mapVxI( ix, iy, VX_LOWERRIGHT )
             if (gmap % mapVxI( ix, iy, VX_LOWERRIGHT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2ITMFillGD: BOTTOM face at "//int2str(ix)//","//int2str(iy)//" has no end node")
             end if
          case( TOP )
             !! start index: 1=start node
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 ) = gmap % mapVxI( ix, iy, VX_UPPERLEFT )
             if (gmap % mapVxI( ix, iy, VX_UPPERLEFT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2ITMFillGD: TOP face at "//int2str(ix)//","//int2str(iy)//" has no start node")
             end if
             !! end vertex: 2=end node
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 2 ) = gmap % mapVxI( ix, iy, VX_UPPERRIGHT )
             if (gmap % mapVxI( ix, iy, VX_UPPERRIGHT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2ITMFillGD: TOP face at "//int2str(ix)//","//int2str(iy)//" has no end node")
             end if
          end select

          !! Neighbour faces of this face
          !! Left neighbour: face continuing to the left of this face
          nix = leftix( ix, iy )
          niy = leftiy( ix, iy )
          if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour( ifc, 1, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( ifc ) )
          end if
          !! Right neighbour: face continuing to the right of this face
          nix = rightix( ix, iy )
          niy = rightiy( ix, iy )
          if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour( ifc, 2, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( ifc ) )
          end if

          !! measure: area
          if (present(gs)) itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % measure( ifc, 1 ) = gs(ix, iy, ALIGNX)
      end do

      !! y-aligned faces
      do ifc = gmap % nfcx + 1, gmap % nfcx + gmap % nfcy
          !! get position of this face in the b2 grid
          ix = gmap % mapFcix( ifc )
          iy = gmap % mapFciy( ifc )
!!$          if (gmap%mapCvI(ix, iy) == GRID_UNDEFINED) then
!!$                  call logmsg(LOGWARNING, "b2ITMFillGD: writing out faces for unused cell "//int2str(ix)//","//int2str(iy))
!!$          end if

          select case ( gmap % mapFcIFace( ifc ) )
          case( LEFT )
              !! start index: 1=start node
              itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 ) = gmap % mapVxI( ix, iy, VX_LOWERLEFT )
              if (gmap % mapVxI( ix, iy, VX_LOWERLEFT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2ITMFillGD: LEFT face at "//int2str(ix)//","//int2str(iy)//" has no start node")
              end if
          !! end vertex: 2=end node
              itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 2 ) = gmap % mapVxI( ix, iy, VX_UPPERLEFT )
              if (gmap % mapVxI( ix, iy, VX_UPPERLEFT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2ITMFillGD: LEFT face at "//int2str(ix)//","//int2str(iy)//" has no end node")
          end if
              !if (itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 )


          case( RIGHT )
              !! start index: 1=start node
              itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 ) = gmap % mapVxI( ix, iy, VX_LOWERRIGHT )
              if (gmap % mapVxI( ix, iy, VX_LOWERRIGHT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2ITMFillGD: RIGHT face at "//int2str(ix)//","//int2str(iy)//" has no start node")
          end if
              !! end vertex: 2=end node
              itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 2 ) = gmap % mapVxI( ix, iy, VX_UPPERRIGHT )
              if (gmap % mapVxI( ix, iy, VX_UPPERRIGHT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2ITMFillGD: RIGHT face at "//int2str(ix)//","//int2str(iy)//" has no end node")
              end if
          end select


          !! Neighbour faces of this face
          !! Bottom neighbour: face continuing to the bottom of this face
          nix = bottomix( ix, iy )
          niy = bottomiy( ix, iy )
          if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour( ifc, 1, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( ifc ) )
          end if
          !! Top neighbour: face continuing to the top of this face
          nix = topix( ix, iy )
          niy = topiy( ix, iy )
          if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour( ifc, 2, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( ifc ) )
          end if
          !! measure: area
          if (present(gs)) itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % measure( ifc, 1 ) = gs(ix, iy, ALIGNY)*qc(ix, iy)
      end do

      !! 2d objects: cells
      !! ...have four boundaries
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( gmap%ncv, 4) )
      !! 2d object measure: cell volume
      if (present(vol)) allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % measure( gmap % ncv, 1 ) )
      !! Also store additional geometry information: position in computational space
      !! FIXME: this should go into alternate geometry, which is not available yet for grid objects
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % geo(gmap%ncv, 2, 1, 1) )

      !! first set all boundary information to undefined
      itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary = GRID_UNDEFINED

      do icv = 1, gmap % ncv
          ix = gmap % mapCvix( icv )
          iy = gmap % mapCviy( icv )

          !! Set position in computational space
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % geo(icv, 1, 1, 1) = ix
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % geo(icv, 2, 1, 1) = iy
          !! put faces composing the quadliateral in the list: left face (y-aligned)
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( icv, 1 ) = gmap % mapFcI( ix, iy, LEFT )
          !! bottom face (x-aligned
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( icv, 2 ) = gmap % mapFcI( ix, iy, BOTTOM )
          !! right face (y-aligned)
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( icv, 3 ) = gmap % mapFcI( ix, iy, RIGHT )
          !! top face (x-aligned)
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( icv, 4 ) = gmap % mapFcI( ix, iy, TOP )
      end do

      !! Fill in connectivity information
      !! ...have one neighbour per boundary
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % neighbour( gmap%ncv, 4, 1) )
      !! first set all to undefined
      itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % neighbour = GRID_UNDEFINED

      do icv = 1, gmap % ncv
          ix = gmap % mapCvix( icv )
          iy = gmap % mapCviy( icv )

          do dir = LEFT, TOP
             call get_Neighbour(nx, ny, leftix, leftiy, rightix, rightiy, topix, topiy, bottomix, bottomiy, &
                  & ix, iy, dir, nix, niy)
             if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % neighbour(icv, dir+1, 1) = gmap % mapCvI( nix, niy )
             end if
          end do

      end do

      !! Fill in x-point indices
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % xpoints( gmap % nsv ) )
      itmgrid % spaces(SPACE_POLOIDALPLANE) % xpoints = gmap % svi(1:gmap % nsv)

      !! If requested, add a second space for the toroidal angle
      if (SPACE_COUNT == SPACE_TOROIDALANGLE) then

          if ( TOROIDAL_PERIODIC ) then
              call gridSetupStruct1dSpace( itmgrid%spaces(SPACE_TOROIDALANGLE), &
                  & COORDTYPE_PHI, &
                  & (/  ( ( 2*pi / NNODES_TOROIDAL ) * i, i = 0, NNODES_TOROIDAL - 1 ) /), &
                  & .true. ) !! periodic = .true.
          else
              call gridSetupStruct1dSpace( itmgrid%spaces(SPACE_TOROIDALANGLE), &
                  & COORDTYPE_PHI, &
                  & (/  ( ( 2*pi / NNODES_TOROIDAL ) * i, i = 0, NNODES_TOROIDAL ) /) , &
                  & .false. ) !! periodic = .false.
          end if

      end if

    end subroutine fill_In_Grid_Desc


    !! Part 2: define subgrids
    subroutine fillInSubGridDescription

      !! internal
      integer :: geoId, iRegion, subgridCount, iType, nSubgrid
      integer :: xIn, yIn, xOut, yOut, iCoreSg
      integer :: cls(SPACE_COUNT_MAX)
      integer, allocatable :: xpoints(:,:)

      geoId = geometryId(nnreg, periodic_bc, topcut)

      !! Figure out total number of subgrids
      !! Do generic subgrids + subgrids
      nSubgrid = B2_GENERIC_SUBGRID_COUNT + regionCountTotal(geoId)
      !! Inner/outer midplane subgrids
      nSubgrid = nSubgrid + 2

      call logmsg( LOGDEBUG, "b2ITMFillGridDescription: expecting total of "&
          &//Int2str(nSubgrid)//" subgrids" )
      allocate( itmgrid % subgrids( nSubgrid ) )

      !! Set up generic subgrids

      !! B2_SUBGRID_CELLS: all 2d cells, one implicit object list
      call createSubGridForClass( itmgrid, itmgrid % subgrids( B2_SUBGRID_CELLS ), &
          & CLASS_CELL(1:SPACE_COUNT), 'Cells' )

      !! B2_SUBGRID_NODES: all nodes, one implicit object list
      call createSubGridForClass( itmgrid, itmgrid % subgrids( B2_SUBGRID_NODES ), &
          & CLASS_NODE(1:SPACE_COUNT), 'Nodes' )

      !! B2_SUBGRID_FACES: all faces, one implicit object list
      call createSubGridForClass( itmgrid, itmgrid % subgrids( B2_SUBGRID_FACES ), &
          & CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT), 'Faces' )

      !! B2_SUBGRID_FACES_X: x-aligned faces. One implicit object list, range over x faces
      !! Create subgrid with one object list
      call createSubGrid( itmgrid % subgrids( B2_SUBGRID_FACES_X ), 1, 'x-aligned faces' )
      !! Initialize implicit object list for faces (class (/1/) )
      call createImplicitObjectList( itmgrid, itmgrid % subgrids( B2_SUBGRID_FACES_X ) % list(1), &
          & CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT) )
      itmgrid % subgrids( B2_SUBGRID_FACES_X ) % list(1) % indset(1) &
          & = createIndexListForRange( 1, gmap%nfcx )
      if ( SPACE_COUNT == SPACE_TOROIDALANGLE ) then
          itmgrid % subgrids( B2_SUBGRID_FACES_X ) % list(1) % indset(2) &
              & = createIndexListForRange( 1, 1 )
      end if

      !! B2_SUBGRID_FACES_Y: y-aligned faces. One implicit object list, range over y faces. Same procedure.
      call createSubGrid( itmgrid % subgrids( B2_SUBGRID_FACES_Y ), 1, 'y-aligned faces' )
      call createImplicitObjectList( itmgrid, itmgrid % subgrids( B2_SUBGRID_FACES_Y ) % list(1)&
          & , CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT) )
      itmgrid % subgrids( B2_SUBGRID_FACES_Y ) % list(1) % indset(1) &
          & = createIndexListForRange( gmap%nfcx + 1, gmap%nfcx + gmap%nfcy )
      if ( SPACE_COUNT == SPACE_TOROIDALANGLE ) then
          itmgrid % subgrids( B2_SUBGRID_FACES_Y ) % list(1) % indset(2) &
              & = createIndexListForRange( 1, 1 )
      end if

      !! Subgrid of all x-points (in one poloidal plane at toroidal index 1)
      !! Assemble object descriptor for x-points
      allocate( xpoints(gmap%nsv, SPACE_COUNT) )
      xpoints = 1
      xpoints(:, SPACE_POLOIDALPLANE) = gmap%svi(1:gmap%nsv)
      call createSubGridForExplicitList( itmgrid, itmgrid % subgrids( B2_SUBGRID_XPOINTS ), &
          & CLASS_NODE(1:SPACE_COUNT), xpoints, 'x-points' )

      !! Set up specific subgrids by collection faces for regions

      !! Start counting from end of generic subgrids
      subgridCount = B2_GENERIC_SUBGRID_COUNT

      !! Cell + face subgrids
      do iType = REGIONTYPE_CELL, REGIONTYPE_YFACE

          select case(iType)
          case( REGIONTYPE_CELL )
              cls = CLASS_CELL
          case( REGIONTYPE_YFACE, REGIONTYPE_XFACE )
              cls = CLASS_POLOIDALRADIAL_FACE
          end select

          do iRegion = 1, regionCount(geoId, iType)
              subgridCount = subgridCount + 1

              call logmsg( LOGDEBUG, "b2ITMFillGridDescription: add subgrid #"//int2str(subgridCount)//&
                  & " for iType "//int2str(iType)//&
                  &", iRegion "//int2str(iRegion)//": "//regionName(geoId, iType, iRegion) )

              call createSubGridForExplicitList( itmgrid, itmgrid % subgrids( subgridCount ), &
                  & cls(1:SPACE_COUNT), &
                  & collectIndexListForRegion(gmap, region, iType, iRegion), &
                  & regionName(geoId, iType, iRegion) )

          end do
      end do

      !! Add midplane node subgrids

      !! Find the core boundary subgrid by looking for its name as defined in b2mod_connectivity
      iCoreSg = gridFindSubGridByName(itmgrid, "Core boundary")
      !! For double null, we need the outer half of the core boundary
      if (iCoreSg == GRID_UNDEFINED) then
          iCoreSg = gridFindSubGridByName(itmgrid, "Outer core boundary")
      end if
      if (iCoreSg == GRID_UNDEFINED) stop "fillInSubGridDescription: "// &
          & "did not find core boundary subgrid for assembling outer midplane subgrid"

      !! Figure out starting points for inner and outer midplane on core boundary
      call find_Midplane_Cells(itmgrid%subgrids(iCoreSg), gmap, crx, xIn, yIn, xOut, yOut)

      subgridCount = subgridCount + 1
      call createSubGridForExplicitList( itmgrid, itmgrid % subgrids( subgridCount ), &
          & CLASS_NODE(1:SPACE_COUNT), &
          & collectRadialVertexIndexList(gmap, cflag, xIn, yIn, topix, topiy), &
          & "Inner midplane" )

      subgridCount = subgridCount + 1
      call createSubGridForExplicitList( itmgrid, itmgrid % subgrids( subgridCount ), &
          & CLASS_NODE(1:SPACE_COUNT), &
          & collectRadialVertexIndexList(gmap, cflag, xOut, yOut, topix, topiy), &
          & "Outer midplane" )

      call logmsg( LOGDEBUG, "b2ITMFillGridDescription: wrote total of "&
          &//int2str(subgridCount)//" subgrids (expected was "//int2str(size(itmgrid%subgrids))//")" )

      call assert( subgridCount == size(itmgrid%subgrids) )
    end subroutine fillInSubGridDescription

  end subroutine b2ITMFillGridDescription

  !> Figure out starting cells for inner and outer midplane on core boundary
  !> by finding the points on the core boundary with minimum and maximum r positions
  subroutine find_Midplane_Cells(coreBndSubgrid, gmap, crx, xIn, yIn, xOut, yOut)
    type(type_complexgrid_subgrid), intent(in) :: coreBndSubgrid
    type(B2GridMap), intent(in) :: gmap
    !! x/radial vertex coordinates
    real (ITM_R8), intent(in) :: crx(-1:gmap%b2nx,-1:gmap%b2ny,0:3)
    integer, intent(out) :: xIn, yIn, xOut, yOut


    !! internal
    real(ITM_R8) :: rMin, rMax
    type(GridObject) :: obj
    integer :: ix, iy, iObj

    rMin = huge(rMin)
    rMax = -huge(rMax)

    xIn = huge(xIn)
    xOut = huge(xOut)

    !! Loop over all faces in core boundary subgrid
    do iObj = 1, gridSubGridSize(coreBndSubgrid)
        obj = subGridGetObject(coreBndSubgrid, iObj)
        !! Expect a face
        call assert( all(obj%cls(1:SPACE_COUNT) == CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT)) )
        !! ...which is aligned along the x-direction
        call assert( gmap % mapFcIFace(obj%ind(SPACE_POLOIDALPLANE)) == BOTTOM )
        ix = gmap % mapFcix( obj%ind(SPACE_POLOIDALPLANE) )
        iy = gmap % mapFciy( obj%ind(SPACE_POLOIDALPLANE) )

        !! We want the vertex associated with the cell at ix, iy, which is number 0
        if ( crx(ix, iy, 0) < rMin ) then
            rMin = crx(ix, iy, 0)
            xIn = ix
            yIn = iy
        end if
        if ( crx(ix, iy, 0) > rMax ) then
            rMax = crx(ix, iy, 0)
            xOut = ix
            yOut = iy
        end if
    end do

    call xertst(xIn /= huge(xIn),   &
        &   "find_Midplane_Cells: did not find inner midplane position")
    call xertst(xOut /= huge(xOut), &
        &   "find_Midplane_Cells: did not find outer midplane position")
  end subroutine

#endif
#endif

    !> Collect the grid indices of all vertices on the radial grid line outward
    !! starting at the vertex at position six,siy in computational space.
    function collectRadialVertexIndexList( gmap, cflag, six, siy, topix,    &
        &   topiy ) result( indexList )
        type(B2GridMap), intent(in) :: gmap !< The grid mapping as computed by
            !< b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS
        integer, allocatable, dimension(:,:) :: indexList   !< 2D index list

        integer, intent(in) :: six  !< First vertex position
        integer, intent(in) :: siy  !< Second vertex position
        integer, intent(in) ::  &
            &   cflag( -1:gmap%b2nx, -1:gmap%b2ny, CARREOUT_NCELLFLAGS )

        !! B2 connectivity array
        integer, intent(in) :: topix( -1:gmap%b2nx, -1:gmap%b2ny ) !< Top
            !< neighbour poloidal (first coordinate) index array
        integer, intent(in) :: topiy( -1:gmap%b2nx, -1:gmap%b2ny ) !< Top
            !< neighbour radial (second coordinate) index

        !! Internal variables
        integer :: ix   !< x-aligned cell index
        integer :: iy   !< y-aligned cell index
        integer :: nix
        integer :: niy
        integer :: nVx
        integer :: iVx

        !! First figure out how many points we have: start at six, siy,
        !! go towards top until running out of physical domain
        nVx = 1
        ix = six
        iy = siy
        do
            !! Take a step upwards
            nix = topix(ix, iy)
            niy = topiy(ix, iy)
            !! Stepped outside grid or into ghost cell?
            if (is_Unneeded_Cell( gmap%b2nx, gmap%b2ny, cflag, .true., nix, niy) ) then
                exit
            else
                nVx = nVx + 1
                ix = nix
                iy = niy
            end if
        end do

        allocate( indexList(nVx, SPACE_COUNT) )
        if (SPACE_COUNT == SPACE_TOROIDALANGLE)&
            & indexList(:,SPACE_TOROIDALANGLE) = 1

        !! collect indices: repeat above loop with storing indices
        ix = six
        iy = siy
        !! Store starting point index
        iVx = 1
        indexList(iVx, SPACE_POLOIDALPLANE) = gmap%mapVxI(ix, iy, VX_LOWERLEFT)
        do
            !! take a step
            nix = topix(ix, iy)
            niy = topiy(ix, iy)

            !! Stepped outside grid?
            if (is_Unneeded_Cell( gmap%b2nx, gmap%b2ny, cflag, .true., nix, niy) ) then
                exit
            end if

            !! Store index for new point
            iVx = iVx + 1
            if (gmap%mapVxI(nix, niy, VX_LOWERLEFT) /= B2_GRID_UNDEFINED) then
               indexList(iVx, SPACE_POLOIDALPLANE) =    &
                &   gmap%mapVxI(nix, niy, VX_LOWERLEFT)
            else if (gmap%mapVxI(ix, iy, VX_UPPERLEFT) /= B2_GRID_UNDEFINED &
                 & .and. iVx == nVx) then
               indexList(iVx, SPACE_POLOIDALPLANE) =    &
                &   gmap%mapVxI(ix, iy, VX_UPPERLEFT)
            else
               stop "collectRadialVertexIndexList: cannot find expected vertex index"
            end if

            ix = nix
            iy = niy
        end do

        call assert( iVx == nVx )

    end function collectRadialVertexIndexList

    !> Collect the grid indices of all vertices on the radial grid line outward
    !! starting at the vertex at position six,siy in computational space.
    subroutine collectRadialVertexIndexListSubroutine(gmap, cflag, six, siy, &
            &   topix, topiy, indexList)
        type(B2GridMap), intent(in) :: gmap !< The grid mapping as computed by
            !< b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS
        integer, allocatable, dimension(:,:), intent(out) :: indexList !< 2D
            !< index list

        integer, intent(in) :: six  !< First vertex position
        integer, intent(in) :: siy  !< Second vertex position
        integer, intent(in) ::  &
            &   cflag( -1:gmap%b2nx, -1:gmap%b2ny, CARREOUT_NCELLFLAGS )

        !! B2 connectivity array
        integer, intent(in) :: topix( -1:gmap%b2nx, -1:gmap%b2ny ) !< Top
            !< neighbour poloidal (first coordinate) index array
        integer, intent(in) :: topiy( -1:gmap%b2nx, -1:gmap%b2ny ) !< Top
            !< neighbour radial (second coordinate) index

        !! Internal variables
        integer :: ix   !< x-aligned cell index
        integer :: iy   !< y-aligned cell index
        integer :: nix
        integer :: niy
        integer :: nVx
        integer :: iVx

        !! First figure out how many points we have: start at six, siy,
        !! go towards top until running out of physical domain
        nVx = 1
        ix = six
        iy = siy
        do
            !! Take a step upwards
            nix = topix(ix, iy)
            niy = topiy(ix, iy)
            !! Stepped outside grid or into ghost cell?
            if (is_Unneeded_Cell( gmap%b2nx, gmap%b2ny, cflag, .true., nix, niy) ) then
                exit
            else
                nVx = nVx + 1
                ix = nix
                iy = niy
            end if
        end do

        allocate( indexList(nVx, SPACE_COUNT) )
        if (SPACE_COUNT == SPACE_TOROIDALANGLE)&
            & indexList(:,SPACE_TOROIDALANGLE) = 1

        !! collect indices: repeat above loop with storing indices
        ix = six
        iy = siy
        !! Store starting point index
        iVx = 1
        indexList(iVx, SPACE_POLOIDALPLANE) = gmap%mapVxI(ix, iy, VX_LOWERLEFT)
        do
            !! take a step
            nix = topix(ix, iy)
            niy = topiy(ix, iy)

            !! Stepped outside grid?
            if (is_Unneeded_Cell( gmap%b2nx, gmap%b2ny, cflag, .true., nix, niy) ) then
                exit
            end if

            !! Store index for new point
            iVx = iVx + 1
            if (gmap%mapVxI(nix, niy, VX_LOWERLEFT) /= B2_GRID_UNDEFINED) then
                indexList(iVx, SPACE_POLOIDALPLANE) =   &
                    &   gmap%mapVxI(nix, niy, VX_LOWERLEFT)
            else if (gmap%mapVxI(ix, iy, VX_UPPERLEFT) /= B2_GRID_UNDEFINED &
                    & .and. iVx == nVx) then
                indexList(iVx, SPACE_POLOIDALPLANE) =   &
                    &   gmap%mapVxI(ix, iy, VX_UPPERLEFT)
            else
                stop "collectRadialVertexIndexListSubroutine: cannot " // &
                   & "find expected vertex index"
            end if

            ix = nix
            iy = niy
        end do

        call assert( iVx == nVx )

    end subroutine collectRadialVertexIndexListSubroutine


    !> Build an index list of all objects of a given region type
    !! (b2mod_connectivity.REGIONTYPE_*) for a given region id.
    !! @result The list of indices for all objects that constitute this grid
    !! region. The array has two dimensions because it is given as a list of
    !! object descriptors.
    function collectIndexListForRegion(gmap, region, iRegionType,   &
        &   iRegion) result( indexList )
        integer, allocatable, dimension(:,:) :: indexList

        type(B2GridMap), intent(in) :: gmap
        integer, intent(in) :: region(-1:gmap%b2nx,-1:gmap%b2ny,0:2)
        integer, intent(in) :: iRegionType, iRegion

        !! internal
        integer :: ix, iy, nInd, iInd, ind

        !! Figure out how many indices to expect. A simple count of the form
        !! nInd = count( region(:,:,iRegionType) == iRegion )
        !! will not do, because we have to account for removed objects (ghost cells/faces).

        !! search the relevant objects and count them
        nInd = 0
        do ix = -1, gmap%b2nx
            do iy = -1, gmap%b2ny

                if ( region(ix, iy, iRegionType) == iRegion ) then
                  !! Get index depending on what object type we are looking at
                    select case (iRegionType)
                    case (REGIONTYPE_CELL)
                        ind = gmap%mapCvI(ix, iy)
                    case (REGIONTYPE_XFACE)
                        ind = gmap%mapFcI(ix, iy, LEFT)
                    case (REGIONTYPE_YFACE)
                        ind = gmap%mapFcI(ix, iy, BOTTOM)
                    end select

                    !! Only count this index if not undefined
                    if ( ind /= B2_GRID_UNDEFINED ) nInd = nInd + 1
                end if

            end do
        end do

        allocate( indexList(nInd, SPACE_COUNT) )
        indexList = 1

        !! search the relevant objects and store their index consecutively
        iInd = 0
        do ix = -1, gmap%b2nx
            do iy = -1, gmap%b2ny

                if ( region(ix, iy, iRegionType) == iRegion ) then
                  !! Get index depending on what object type we are looking at
                    select case (iRegionType)
                    case (REGIONTYPE_CELL)
                        ind = gmap%mapCvI(ix, iy)
                    case (REGIONTYPE_XFACE)
                        ind = gmap%mapFcI(ix, iy, LEFT)
                    case (REGIONTYPE_YFACE)
                        ind = gmap%mapFcI(ix, iy, BOTTOM)
                    end select

                    if ( ind /= B2_GRID_UNDEFINED ) then
                        iInd = iInd + 1
                        call assert(iInd <= nInd)
                        indexList( iInd, SPACE_POLOIDALPLANE ) = ind
                    end if
                end if

            end do
        end do

        call assert( iInd == nInd )

    end function collectIndexListForRegion

    !> Build an index list of all objects of a given region type
    !! (b2mod_connectivity.REGIONTYPE_*) for a given region id.
    !! @result The list of indices for all objects that constitute this grid
    !! region. The array has two dimensions because it is given as a list of
    !! object descriptors.
    subroutine collectIndexListForRegionSubroutine(gmap, region, iRegionType, &
        &   iRegion, indexlist)
        type(B2GridMap), intent(in) :: gmap !< The grid mapping as computed by
            !< b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS
        integer, allocatable, dimension(:,:), intent(out) :: indexList !< 2D
            !< index list
        integer, intent(in) :: region( -1:gmap%b2nx, -1:gmap%b2ny, 0:2 )
        integer, intent(in) :: iRegionType
        integer, intent(in) :: iRegion

        !! Internal variables
        integer :: ix   !< x-aligned cell index
        integer :: iy   !< y-aligned cell index
        integer :: nInd
        integer :: iInd
        integer :: ind

        !! Figure out how many indices to expect. A simple count of the form
        !! nInd = count( region(:,:,iRegionType) == iRegion )
        !! will not do, because we have to account for removed objects
        !! (ghost cells/faces).

        !! search the relevant objects and count them
        nInd = 0
        do ix = -1, gmap%b2nx
            do iy = -1, gmap%b2ny

                if ( region(ix, iy, iRegionType) == iRegion ) then
                  !! Get index depending on what object type we are looking at
                    select case (iRegionType)
                    case (REGIONTYPE_CELL)
                        ind = gmap%mapCvI(ix, iy)
                    case (REGIONTYPE_XFACE)
                        ind = gmap%mapFcI(ix, iy, LEFT)
                    case (REGIONTYPE_YFACE)
                        ind = gmap%mapFcI(ix, iy, BOTTOM)
                    end select

                    !! Only count this index if not undefined
                    if ( ind /= B2_GRID_UNDEFINED ) nInd = nInd + 1
                end if

            end do
        end do

        allocate( indexList(nInd, SPACE_COUNT) )
        indexList = 1

        !! search the relevant objects and store their index consecutively
        iInd = 0
        do ix = -1, gmap%b2nx
            do iy = -1, gmap%b2ny

                if ( region(ix, iy, iRegionType) == iRegion ) then
                  !! Get index depending on what object type we are looking at
                    select case (iRegionType)
                    case (REGIONTYPE_CELL)
                        ind = gmap%mapCvI(ix, iy)
                    case (REGIONTYPE_XFACE)
                        ind = gmap%mapFcI(ix, iy, LEFT)
                    case (REGIONTYPE_YFACE)
                        ind = gmap%mapFcI(ix, iy, BOTTOM)
                    end select

                    if ( ind /= B2_GRID_UNDEFINED ) then
                        iInd = iInd + 1
                        call assert(iInd <= nInd)
                        indexList( iInd, SPACE_POLOIDALPLANE ) = ind
                    end if
                end if

            end do
        end do

        call assert( iInd == nInd )

    end subroutine collectIndexListForRegionSubroutine

end module b2mod_ual_io_grid

!!!Local Variables:
!!! mode: f90
!!! End:
