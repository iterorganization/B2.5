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
#if IMAS_MINOR_VERSION > 8
    use ids_schemas  & ! IGNORE
     & , only : IDS_real
#endif
#if IMAS_MINOR_VERSION > 11
    use ids_grid_subgrid  & ! IGNORE
     & , only : getGridSubsetSize, getGridSubsetObject, findGridSubsetByName, &
     &          CreateGridSubsetForClass, CreateEmptyGridSubset, &
     &          CreateExplicitObjectListSingleSpace
#if IMAS_MINOR_VERSION < 15
   use ids_grid_object    & ! IGNORE
     & , only : ids_generic_grid_dynamic
#else
   use ids_grid_object    & ! IGNORE
     & , only : ids_generic_grid_aos3_root
#endif
    use ids_grid_object   & ! IGNORE
     & , only : ids_generic_grid_dynamic_grid_subset, &
     &          GRID_SUBSET_NODES, GRID_SUBSET_X_POINTS, GRID_SUBSET_CELLS, &
     &          GridObject
#if GGD_MINOR_VERSION > 9
    use ids_grid_object   & ! IGNORE
     & , only : GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
     &          GRID_SUBSET_EDGES
#else
    use ids_grid_object   & ! IGNORE
     & , only : GRID_SUBSET_X_ALIGNED_FACES, GRID_SUBSET_Y_ALIGNED_FACES, &
     &          GRID_SUBSET_FACES
#endif
    use ids_grid_structured & ! IGNORE
     & , only : GridWriteData, GridSetupStruct1dSpace
    use ids_grid_common     & ! IGNORE
     & , only : COORDTYPE_R, COORDTYPE_Y, COORDTYPE_Z, COORDTYPE_PHI,         &
     &          gridSubsetName, gridSubsetDescription,                        &
     &          GRID_SUBSET_CORE, GRID_SUBSET_SOL,                            &
     &          GRID_SUBSET_CORE_BOUNDARY,                                    &
     &          GRID_SUBSET_INNER_DIVERTOR, GRID_SUBSET_OUTER_DIVERTOR,       &
     &          GRID_SUBSET_INNER_DIVERTOR_INACTIVE,                          &
     &          GRID_SUBSET_OUTER_DIVERTOR_INACTIVE,                          &
     &          GRID_SUBSET_CORE_CUT, GRID_SUBSET_PFR_CUT,                    &
     &          GRID_SUBSET_CORE_CUT_INACTIVE, GRID_SUBSET_PFR_CUT_INACTIVE,  &
     &          GRID_SUBSET_OUTER_THROAT, GRID_SUBSET_INNER_THROAT,           &
     &          GRID_SUBSET_OUTER_THROAT_INACTIVE,                            &
     &          GRID_SUBSET_INNER_THROAT_INACTIVE,                            &
     &          GRID_SUBSET_INNER_TARGET, GRID_SUBSET_OUTER_TARGET,           &
     &          GRID_SUBSET_INNER_TARGET_INACTIVE,                            &
     &          GRID_SUBSET_OUTER_TARGET_INACTIVE,                            &
     &          GRID_SUBSET_MAIN_CHAMBER_WALL, GRID_SUBSET_MAIN_WALL,         &
     &          GRID_SUBSET_PFR_WALL,                                         &
     &          GRID_SUBSET_OUTER_PFR_WALL, GRID_SUBSET_INNER_PFR_WALL,       &
     &          GRID_SUBSET_OUTER_PFR_WALL_INACTIVE,                          &
     &          GRID_SUBSET_INNER_PFR_WALL_INACTIVE,                          &
     &          GRID_SUBSET_OUTER_BAFFLE, GRID_SUBSET_INNER_BAFFLE,           &
     &          GRID_SUBSET_OUTER_BAFFLE_INACTIVE,                            &
     &          GRID_SUBSET_INNER_BAFFLE_INACTIVE,                            &
     &          GRID_SUBSET_ACTIVE_SEPARATRIX, GRID_SUBSET_SEPARATRIX,        &
     &          GRID_SUBSET_SECOND_SEPARATRIX,                                &
     &          GRID_SUBSET_BETWEEN_SEPARATRICES,                             &
     &          GRID_SUBSET_OUTER_MIDPLANE, GRID_SUBSET_INNER_MIDPLANE,       &
     &          GRID_SUBSET_OUTER_MIDPLANE_SEPARATRIX,                        &
     &          GRID_SUBSET_INNER_MIDPLANE_SEPARATRIX,                        &
     &          GRID_SUBSET_INNER_STRIKEPOINT, GRID_SUBSET_OUTER_STRIKEPOINT, &
     &          GRID_SUBSET_INNER_STRIKEPOINT_INACTIVE,                       &
     &          GRID_SUBSET_OUTER_STRIKEPOINT_INACTIVE,                       &
     &          IDS_GRID_UNDEFINED => GRID_UNDEFINED
#endif
#else
# ifdef ITM_ENVIRONMENT_LOADED
    use itm_types , ITM_R8 => R8, ITM_R4 => R4 ! IGNORE
    use euITM_schemas ! IGNORE
    use itm_constants , pi => itm_pi ! IGNORE
    use itm_grid_object , only : GridObject ! IGNORE
    use itm_grid_common & ! IGNORE
     & , only : COORDTYPE_R, COORDTYPE_Y, COORDTYPE_Z, COORDTYPE_PHI
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
    use b2mod_b2cmfs
    use b2mod_ppout

    implicit none

    !! Constants for use with the ITM grid description

    logical, parameter :: INCLUDE_GHOST_CELLS = .false.

    !! Space indices
    integer, parameter :: SPACE_POLOIDALPLANE = 1   !< Space index
                                                    !< (polodal plane)
    integer, parameter :: SPACE_TOROIDALANGLE = 2   !< Space index
                                                    !< (toroidal angle)
    integer, parameter :: SPACE_COUNT = SPACE_TOROIDALANGLE !< Space count
        !< set up:
        !< SPACE_COUNT = SPACE_POLOIDALPLANE: do only the poloidal plane space;
        !< SPACE_COUNT = SPACE_TOROIDALANGLE: will do the full 3d grid with two
        !< spaces
    !! Will have at maximum so many spaces
    integer, parameter :: SPACE_COUNT_MAX = 2

    !! Number of points in toroidal direction
    !! (only 1 makes sense here, this is for playing around)
    integer, parameter :: NNODES_TOROIDAL = 1   !< Number of points in toroidal
        !< direction (only 1 makes sense here, this is for playing around)

    logical, parameter :: TOROIDAL_PERIODIC = .false.   !< Flag controlling
        !< whether the toroidal space is set up periodic or not.
        !< If periodic, the last node is connected to the first by an edge.
        !< If not periodic, an additional node at 2 pi is added

    integer :: jxi, jxa !< B2.5 Inner and Outer midplane default locations

    !! Object class tuples (ITM class definition)
    integer, dimension( SPACE_COUNT_MAX ), parameter :: CLASS_NODE =        &
        &   (/ 0, 0 /)  !< Object class tuple (ITM class definition): Node
    integer, dimension( SPACE_COUNT_MAX ), parameter :: CLASS_RZ_EDGE =     &
        &   (/ 1, 0 /)  !< Object class tuple (ITM class definition): Edge (R, Z)
    integer, dimension( SPACE_COUNT_MAX ), parameter :: CLASS_PHI_EDGE =    &
        &   (/ 0, 1 /)  !< Object class tuple (ITM class definition): Edge (phi)
    integer, dimension( SPACE_COUNT_MAX ), parameter :: &
        &   CLASS_POLOIDALRADIAL_EDGE = (/ 1, 1 /)  !< Object class tuple
        !< (ITM class definition): Poloidal/radial edge (in observed existing
        !< CPO cases this class refers to 2D cells)
    integer, dimension( SPACE_COUNT_MAX ), parameter :: &
        &   CLASS_TOROIDAL_EDGE = (/ 2, 0 /)    !< Object class tuple
        !< (ITM class definition): Toroidal edge
    integer, dimension( SPACE_COUNT_MAX ), parameter :: CLASS_CELL = (/ 2, 1 /)
        !< Object class tuple (ITM class definition): Cell

    !! Object class tuples
    !! (one digit IDS classes, transformed from the above ITM classes )
    !! Primary IDS classes are:
    !!      Class 1 - nodes/vertices (0D objects)
    !!      Class 2 - edges/faces (1D objects)
    !!      Class 3 - 2D cells (2D objects)
    !! Fortran90 does not allow initialization of constants using SUM. This is
    !! permitted in newer Fortran 2003. Current workaround is to directly
    !! specify the primary IDS class constants

    ! integer, parameter :: IDS_CLASS_NODE = sum(CLASS_NODE) + 1
    integer, parameter :: IDS_CLASS_NODE = 1    !< Object class tuple
        !< (IMAS class definition): Node
    ! integer, parameter :: IDS_CLASS_RZ_EDGE = sum(CLASS_RZ_EDGE) + 1
    integer, parameter :: IDS_CLASS_RZ_EDGE = 2 !< Object class tuple
        !< (IMAS class definition): Edge
    ! integer, parameter :: IDS_CLASS_PHI_EDGE = sum(CLASS_PHI_EDGE) + 1
    integer, parameter :: IDS_CLASS_PHI_EDGE = 2    !< Object class tuple
        !< (IMAS class definition): Edge
    ! integer, parameter :: IDS_CLASS_POLOIDALRADIAL_EDGE =   &
        ! &   sum(CLASS_POLOIDALRADIAL_EDGE)
    integer, parameter :: IDS_CLASS_POLOIDALRADIAL_EDGE = 2 !< Object class tuple
        !< (IMAS class definition): Edge
    ! integer, parameter :: IDS_CLASS_TOROIDAL_EDGE =         &
        ! &   sum(CLASS_TOROIDAL_EDGE)
    integer, parameter :: IDS_CLASS_TOROIDAL_EDGE = 2 !< Object class tuple
        !< (IMAS class definition): Edge
    ! integer, parameter :: IDS_CLASS_CELL = sum(IDS_CLASS_CELL)
    integer, parameter :: IDS_CLASS_CELL = 3 !< Object class tuple
        !< (IMAS class definition): Cell (2D)

    !! Subgrid/Grid subset name constants

    integer, parameter :: B2_GENERIC_GSUBSET_COUNT = 6  !< Total number of
        !< generic grid subsets

    !! Generic grid subsets (all cells, all edges)
    !! Note: special grid subsets (given by region ids) do not have specific
    !! constants (see also b2mod_connectivity.f90)

    !! IMAS uses GGD grid subset identifier definitions defined in GSL
    !! (in ids_grid_common)
#if GGD_MINOR_VERSION < 9
    !! IMAS GGD grid subset identifier definitions
    integer, parameter :: GRID_SUBSET_TYPES = 106

    !> x-aligned edges defining part of active separatrix separating core and SOL
    integer, parameter :: GRID_SUBSET_ACTIVE_SEPARATRIX = 26
    !> main_chamber_wall + outer_baffle(s) + inner_baffle(s)
    integer, parameter :: GRID_SUBSET_MAIN_WALL = 27
    !> outer_PFR_wall(s) + inner_PFR_wall(s)
    integer, parameter :: GRID_SUBSET_PFR_WALL = 28
    !> y-aligned edges inside the separatrix connecting to the non-active x-point
    integer, parameter :: GRID_SUBSET_CORE_CUT_INACTIVE = 29
    !> y-aligned edges in the private flux region connecting to the non-active x-point
    integer, parameter :: GRID_SUBSET_PFR_CUT_INACTIVE = 30
    !> y-aligned edges in the outer SOL connecting to the non-active x-point
    integer, parameter :: GRID_SUBSET_OUTER_THROAT_INACTIVE = 31
    !> y-aligned edges in the inner SOL connecting to the non-active x-point
    integer, parameter :: GRID_SUBSET_INNER_THROAT_INACTIVE = 32
    !> x-aligned edges defining the non-active separatrix
    integer, parameter :: GRID_SUBSET_SECOND_SEPARATRIX = 33
    !> x-aligned edges defining the chamber wall of the outer non-active divertor region
    integer, parameter :: GRID_SUBSET_OUTER_BAFFLE_INACTIVE = 34
    !> x-aligned edges defining the chamber wall of the inner non-active divertor region
    integer, parameter :: GRID_SUBSET_INNER_BAFFLE_INACTIVE = 35
    !> x-aligned edges defining the private flux region wall of the outer non-active
    !> divertor region
    integer, parameter :: GRID_SUBSET_OUTER_PFR_WALL_INACTIVE = 36
    !> x-aligned edges defining the private flux region wall of the inner non-active
    !> divertor region
    integer, parameter :: GRID_SUBSET_INNER_PFR_WALL_INACTIVE = 37
    !> Cells between the two separatrices
    integer, parameter :: GRID_SUBSET_BETWEEN_SEPARATRICES = 38
    !> Cells defining the outer inactive divertor region
    integer, parameter :: GRID_SUBSET_OUTER_DIVERTOR_INACTIVE = 39
    !> Cells defining the inner inactive divertor region
    integer, parameter :: GRID_SUBSET_INNER_DIVERTOR_INACTIVE = 40
    !> y-aligned edges defining the outer inactive target
    integer, parameter :: GRID_SUBSET_OUTER_TARGET_INACTIVE = 41
    !> y-aligned edges defining the inner inactive target
    integer, parameter :: GRID_SUBSET_INNER_TARGET_INACTIVE = 42
    !> Point on active separatrix at outer midplane
    integer, parameter :: GRID_SUBSET_OUTER_MIDPLANE_SEPARATRIX = 101
    !> Point on active separatrix at inner midplane
    integer, parameter :: GRID_SUBSET_INNER_MIDPLANE_SEPARATRIX = 102
    !> Point on active separatrix at outer active target
    integer, parameter :: GRID_SUBSET_OUTER_STRIKEPOINT = 103
    !> Point on active separatrix at inner active target
    integer, parameter :: GRID_SUBSET_INNER_STRIKEPOINT = 104
    !> Point on non-active separatrix at outer active target
    integer, parameter :: GRID_SUBSET_OUTER_STRIKEPOINT_INACTIVE = 105
    !> Point on non-active separatrix at inner active target
    integer, parameter :: GRID_SUBSET_INNER_STRIKEPOINT_INACTIVE = 106

    character*(26), parameter, private :: UU = &
       &    'UNSPECIFIED               ' !< Unspecified string

    character*(26), dimension (0:GRID_SUBSET_TYPES), parameter :: gridSubsetName = &
       &  (/   &
       &    'UNSPECIFIED               ' , &
       &    'NODES                     ' , &
       &    'EDGES                     ' , &
       &    'X_ALIGNED_EDGES           ' , &
       &    'Y_ALIGNED_EDGES           ' , &
       &    'CELLS                     ' , &
       &    'X_POINTS                  ' , &
       &    'CORE_CUT                  ' , &
       &    'PFR_CUT                   ' , &
       &    'OUTER_THROAT              ' , &
       &    'INNER_THROAT              ' , &
       &    'OUTER_MIDPLANE            ' , &
       &    'INNER_MIDPLANE            ' , &
       &    'OUTER_TARGET              ' , &
       &    'INNER_TARGET              ' , &
       &    'CORE_BOUNDARY             ' , &
       &    'SEPARATRIX                ' , &
       &    'MAIN_CHAMBER_WALL         ' , &
       &    'OUTER_BAFFLE              ' , &
       &    'INNER_BAFFLE              ' , &
       &    'OUTER_PFR_WALL            ' , &
       &    'INNER_PFR_WALL            ' , &
       &    'CORE                      ' , &
       &    'SOL                       ' , &
       &    'OUTER_DIVERTOR            ' , &
       &    'INNER_DIVERTOR            ' , &
       &    'ACTIVE_SEPARATRIX         ' , &
       &    'MAIN_WALL                 ' , &
       &    'PFR_WALL                  ' , &
       &    'CORE_CUT_INACTIVE         ' , &
       &    'PFR_CUT_INACTIVE          ' , &
       &    'OUTER_THROAT_INACTIVE     ' , &
       &    'INNER_THROAT_INACTIVE     ' , &
       &    'SECOND_SEPARATRIX         ' , &
       &    'OUTER_BAFFLE_INACTIVE     ' , &
       &    'INNER_BAFFLE_INACTIVE     ' , &
       &    'OUTER_PFR_WALL_INACTIVE   ' , &
       &    'INNER_PFR_WALL_INACTIVE   ' , &
       &    'BETWEEN_SEPARATRICES      ' , &
       &    'OUTER_DIVERTOR_INACTIVE   ' , &
       &    'INNER_DIVERTOR_INACTIVE   ' , &
       &    'OUTER_TARGET_INACTIVE     ' , &
       &    'INNER_TARGET_INACTIVE     ' , &
       &     UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, &
       &     UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, &
       &     UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, &
       &     UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, &
       &     UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, &
       &     UU, UU, UU, UU, UU, UU, UU, UU,         &
       &    'OUTER_MIDPLANE_SEPARATRIX ' , &
       &    'INNER_MIDPLANE_SEPARATRIX ' , &
       &    'OUTER_STRIKEPOINT         ' , &
       &    'INNER_STRIKEPOINT         ' , &
       &    'OUTER_STRIKEPOINT_INACTIVE' , &
       &    'INNER_STRIKEPOINT_INACTIVE'   &
       &   /)

    character*(93), parameter, private :: US = & !< Unspecified subset string
       &    'Unspecified grid subset                                                                      '

    character*(93), dimension (0:GRID_SUBSET_TYPES), parameter :: gridSubsetDescription = &
       &  (/   &
       &    'Unspecified grid subset                                                                      ' , &
       &    'All nodes (0D objects)                                                                       ' , &
       &    'All edges (1D objects)                                                                       ' , &
       &    'All x-aligned (poloidally) aligned edges                                                     ' , &
       &    'All y-aligned (radially) aligned edges                                                       ' , &
       &    'All cells (2D objects)                                                                       ' , &
       &    'X-points                                                                                     ' , &
       &    'y-aligned edges inside the separatrix connecting to the x-point                              ' , &
       &    'y-aligned edges in the private flux region connecting to the x-point                         ' , &
       &    'y-aligned edges in the outer SOL connecting to the x-point                                   ' , &
       &    'y-aligned edges in the inner SOL connecting to the x-point                                   ' , &
       &    'y-aligned edges connecting to the node closest to outer midplane on the separatrix           ' , &
       &    'y-aligned edges connecting to the node closest to inner midplane on the separatrix           ' , &
       &    'y-aligned edges defining the outer target                                                    ' , &
       &    'y-aligned edges defining the inner target                                                    ' , &
       &    'Innermost x-aligned edges                                                                    ' , &
       &    'x-aligned edges defining the separatrix                                                      ' , &
       &    'x-aligned edges defining main chamber wall outside of the divertor regions                   ' , &
       &    'x-aligned edges defining the chamber wall of the outer divertor region                       ' , &
       &    'x-aligned edges defining the chamber wall of the inner divertor region                       ' , &
       &    'x-aligned edges defining the private flux region wall of the outer divertor region           ' , &
       &    'x-aligned edges defining the private flux region wall of the inner divertor region           ' , &
       &    'Cells inside the separatrix                                                                  ' , &
       &    'Cells defining the main SOL outside of the divertor regions                                  ' , &
       &    'Cells defining the outer divertor region                                                     ' , &
       &    'Cells defining the inner divertor region                                                     ' , &
       &    'x-aligned edges defining part of active separatrix separating core and SOL                   ' , &
       &    'main_chamber_wall + outer_baffle(s) + inner_baffle(s)                                        ' , &
       &    'outer_PFR_wall(s) + inner_PFR_wall(s)                                                        ' , &
       &    'y-aligned edges inside the separatrix connecting to the non-active x-point                   ' , &
       &    'y-aligned edges in the private flux region connecting to the non-active x-point              ' , &
       &    'y-aligned edges in the outer SOL connecting to the non-active x-point                        ' , &
       &    'y-aligned edges in the inner SOL connecting to the non-active x-point                        ' , &
       &    'x-aligned edges defining the non-active separatrix                                           ' , &
       &    'x-aligned edges defining the chamber wall of the outer non-active divertor region            ' , &
       &    'x-aligned edges defining the chamber wall of the inner non-active divertor region            ' , &
       &    'x-aligned edges defining the private flux region wall of the outer non-active divertor region' , &
       &    'x-aligned edges defining the private flux region wall of the inner non-active divertor region' , &
       &    'Cells between the two separatrices                                                           ' , &
       &    'Cells defining the outer inactive divertor region                                            ' , &
       &    'Cells defining the inner inactive divertor region                                            ' , &
       &    'y-aligned edges defining the outer inactive target                                           ' , &
       &    'y-aligned edges defining the inner inactive target                                           ' , &
       &     US, US, US, US, US, US, US, US, US, US,                                                          &
       &     US, US, US, US, US, US, US, US, US, US,                                                          &
       &     US, US, US, US, US, US, US, US, US, US,                                                          &
       &     US, US, US, US, US, US, US, US, US, US,                                                          &
       &     US, US, US, US, US, US, US, US, US, US,                                                          &
       &     US, US, US, US, US, US, US, US,                                                                  &
       &    'Point on active separatrix at outer midplane                                                 ' , &
       &    'Point on active separatrix at inner midplane                                                 ' , &
       &    'Point on active separatrix at outer active target                                            ' , &
       &    'Point on active separatrix at inner active target                                            ' , &
       &    'Point on non-active separatrix at outer active target                                        ' , &
       &    'Point on non-active separatrix at inner active target                                        '   &
       &   /)
#endif
#ifdef IMAS
#if GGD_MINOR_VERSION < 10
    integer, parameter :: GRID_SUBSET_X_ALIGNED_EDGES = GRID_SUBSET_X_ALIGNED_FACES
    integer, parameter :: GRID_SUBSET_Y_ALIGNED_EDGES = GRID_SUBSET_Y_ALIGNED_FACES
    integer, parameter :: GRID_SUBSET_EDGES = GRID_SUBSET_FACES
#endif
#endif
#ifdef ITM_ENVIRONMENT_LOADED
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
    integer, parameter :: B2_SUBGRID_EDGES = 2  !< Grid subset containing all
        !< edges (first x-aligned, then y-aligned) belonging to the
        !< associated space (order given by grid map)
    integer, parameter :: B2_GSUBSET_EDGES = 2  !< Grid subset containing all
        !< edges (first x-aligned, then y-aligned) belonging to the
        !< associated space (order given by grid map)
    integer, parameter :: B2_SUBGRID_EDGES_X  = 3   !< Grid subset containing
        !< all x-aligned (poloidally) aligned edges belonging to the associated
        !< space (order given by grid map)
    integer, parameter :: B2_GSUBSET_X_ALIGNED_EDGES = 3    !< Grid subset
        !< containing  all x-aligned (poloidally) aligned edges belonging to
        !< the associated  space (order given by grid map)
    integer, parameter :: B2_SUBGRID_EDGES_Y = 4    !< Grid subset containing
        !< all y-aligned (radially) aligned edges belonging to the associated
        !< space (order given by grid map)
    integer, parameter :: B2_GSUBSET_Y_ALIGNED_EDGES = 4 !< Grid subset
        !< containing  all y-aligned (radially) aligned edges belonging to
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

#if IMAS_MINOR_VERSION > 11
    !> Routine that fills in a grid description which is part of a IMAS IDS
    !! using the given grid data and prepared mappings
    subroutine b2_IMAS_Fill_Grid_Desc( gmap, grid_ggd, nx, ny, crx, cry,    &
        &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,       &
        &   bottomiy, nnreg, topcut, region, cflag, includeGhostCells, vol, &
        &   gs, qc )
        type(B2GridMap), intent(in) :: gmap !< The grid mapping as computed
            !< by b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS
#if IMAS_MINOR_VERSION < 15
        type(ids_generic_grid_dynamic), intent(out) :: grid_ggd !< Type of IDS
            !< data structure, designed for handling grid geometry data
#else
        type(ids_generic_grid_aos3_root), intent(out) :: grid_ggd !< Type of IDS
            !< data structure, designed for handling grid geometry data
#endif
        integer, intent(in) :: nx   !< Number of interior cells
            !< along the first coordinate
            !< Used to define size of grid arrays: (-1:nx, -1:ny)
        integer, intent(in) :: ny   !< Number of interior cells
            !< along the second coordinate
            !< Used to define size of grid arrays: (-1:nx, -1:ny)

        !! Output arguments
        real(R8), intent(in) :: crx( -1:nx, -1:ny, 0:3 )    !< Horizontal vertex
            !< coordinates of the four corners of the (ix, iy) cell
        real(R8), intent(in) :: cry( -1:nx, -1:ny, 0:3 )    !< Vertical vertex
            !< coordinates of the four corners of the (ix, iy) cell

        !! B2 connectivity array
        integer, intent(in) :: leftix( -1:nx, -1:ny )   !< Left neighbour
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
        integer cflag( -1:nx, -1:ny, CARREOUT_NCELLFLAGS ) !< Cell flag
        logical, intent(in) :: includeGhostCells    !< Include "fake" cells
        !! Optional B2 measure information
        real(R8), intent(in), optional :: vol( -1:nx, -1:ny) !< Cell volume
        real(R8), intent(in), optional :: gs( -1:nx, -1:ny, 0:2)
        real(R8), intent(in), optional :: qc(-1:nx,-1:ny)   !< Cosine of the
            !< angle between flux line direction and left cell face
        real(R8), save :: width = 1.0_R8

        !! Internal variables
        integer, parameter :: NDIM = 2  !< Dimension of the space

        call ipgetr ('b2agmt_1d_width', width)
        call xertst (0.0_R8.lt.width, 'faulty input width')
        call xertst( present( gs ) .EQV. present( qc ) , &
            "Assert error ( gz or qc missing ) in b2_IMAS_Fill_Grid_Desc" )

        !! Set GGD grid geometry
        call fill_In_Grid_Desc()
        !! Set grid subsets
        call fill_In_GridSubset_Desc()

contains
    !> Fill in the general grid description
    subroutine fill_In_Grid_Desc()
        !! Internal variables
        integer :: iVx  !< Vertex/node index
        integer :: iFc  !< Face/edge index
        integer :: iCv  !< Cell index
        integer :: ix   !< x-aligned (poloidal) cell index
        integer :: iy   !< y-aligned (radial) cell index
        integer :: nix
        integer :: niy
        integer :: i    !< Iterator
        integer :: dir
        integer :: nFc  !< Number of all faces/edges (x + y aligned)
        integer :: geometryType  !< Geometry identifier index

        geometryType = geometryId(nnreg, isymm, periodic_bc, topcut)

        allocate( grid_ggd%identifier%name(1) )
        grid_ggd%identifier%name = geometryName(geometryType)
        grid_ggd%identifier%index = geometryType
        allocate( grid_ggd%identifier%description(1) )
        grid_ggd%identifier%description = geometryDescription(geometryType)

        allocate( grid_ggd%space( SPACE_COUNT ) )

#if IMAS_MINOR_VERSION > 19
        allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%identifier%name(1) )
        allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%identifier%description(1) )
        grid_ggd%space( SPACE_POLOIDALPLANE )%identifier%name = "Standard grid"
        grid_ggd%space( SPACE_POLOIDALPLANE )%identifier%index = 1
        grid_ggd%space( SPACE_POLOIDALPLANE )%identifier%description = labgeo
#endif

        !! Coordinate types
        !! dimension of space = NDIM = size( coordtype )

        grid_ggd%space( SPACE_POLOIDALPLANE )%geometry_type%index = 0
        allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%geometry_type%name(1) )
        grid_ggd%space( SPACE_POLOIDALPLANE )%geometry_type%name = 'Poloidal'
        allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%geometry_type%description(1) )
        grid_ggd%space( SPACE_POLOIDALPLANE )%geometry_type%description = &
            &   'Poloidal plane cross-section'

        !! Set the space coordinates, also defining the dimension of the space
        allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%coordinates_type(NDIM) )
        grid_ggd%space( SPACE_POLOIDALPLANE )%coordinates_type(1) = COORDTYPE_R
        grid_ggd%space( SPACE_POLOIDALPLANE )%coordinates_type(2) = COORDTYPE_Z

        !! Three types of objects: 0D vertices/nodes, 1D faces/edges, 2D cells
        allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension( NDIM + 1 ) )

        !! Allocate the number of objects of each type:
        !! 0D vertices/nodes
        !! For SN: gmap%nVx = ( nx+1 )*( ny+1 ) - 1
        allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(1)%object( gmap%nVx ) )
        !! 1D faces/edges
        !! For SN: gmap%nFcx + gmap%nFcy = nx*( ny+1 ) + ( nx+1 )*ny
        allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(2)%object( gmap%nFcx + gmap%nFcy ) )
        !! 2D cells
        !! For SN: gmap%nCv = nx*ny
        allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(3)%object( gmap%nCv ) )

        !! Fill in vertex/node information
        do iVx = 1, gmap%nVx
            !! Allocate geometry leaf for each vertex/node
            allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(1)%object( iVx )%geometry(2) )
            !! Set geometry (R and Z coordinates) of each vertex/node
            !! The way of writing vertices/nodes data identically to the
            !! conversion to IDS of
            !! CPO 16151/1000 case, available on (written in 11 October 2017)
            !! ITER HPC cluster in directory
            !! /home/ITER/penkod/public/imasdb/solps-iter/3/0
            !! or
            !! /home/ITER/kosl/public/imasdb/solps-iter/3/0
            !! while on Marconi Gateway
            !! /marconi_work/eufus_gw/work/g2penkod/imasdb/solps-iter/3/0
            !! or
            !! /marconi_work/eufus_gw/work/g2kosl/imasdb/solps-iter/3/0
            if (isymm.eq.0 .or. isymm.eq.1 .or. isymm.eq.2) then
              grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%  &
                  &   object( iVx )%geometry(1) = crx( gmap%mapVxix( iVx ),    &
                                                   & gmap%mapVxiy( iVx ),      &
                                                   & gmap%mapVxIVx( iVx ) )
              grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%  &
                  &   object( iVx )%geometry(2) = cry( gmap%mapVxix( iVx ),    &
                                                   & gmap%mapVxiy( iVx ),      &
                                                   & gmap%mapVxIVx( iVx ) )
            else if (isymm.eq.3 .or. isymm.eq.4) then
              grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%  &
                  &   object( iVx )%geometry(1) = cry( gmap%mapVxix( iVx ),    &
                                                   & gmap%mapVxiy( iVx ),      &
                                                   & gmap%mapVxIVx( iVx ) )
              grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%  &
                  &   object( iVx )%geometry(2) = crx( gmap%mapVxix( iVx ),    &
                                                   & gmap%mapVxiy( iVx ),      &
                                                   & gmap%mapVxIVx( iVx ) )
            end if

            !! Set additional vertex/node index (REQUIRED!)
            allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(1)%object( iVx )%nodes(1) )
            grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%   &
                &   object( iVx )%nodes(1) = iVx
        end do

#if 0
        !! Set geometry (R and Z coordinates) of each node
        !! (original b2mod. This node writing order does not work well with the
        !! code used to write 1D objects data (faces/edges). Thats why it is
        !! currently not used (#if 0) and the method above is used instead)
        iVx = 0
        do iy = 0, ny-1
            do ix = 0, nx-1
                iVx = iVx + 1 !! Lower left corners
                if (isymm.eq.0 .or. isymm.eq.1 .or. isymm.eq.2) then
                    grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                        &   objects_per_dimension(1)%object(iVx)%geometry(1) =  &
                        &   crx(ix,iy,0)
                    grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                        &   objects_per_dimension(1)%object(iVx)%geometry(2) =  &
                        &   cry(ix,iy,0)
                else if (isymm.eq.3 .or. isymm.eq.4) then
                    grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                        &   objects_per_dimension(1)%object(iVx)%geometry(1) =  &
                        &   cry(ix,iy,0)
                    grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                        &   objects_per_dimension(1)%object(iVx)%geometry(2) =  &
                        &   crx(ix,iy,0)
                end if
                if( iVx .eq. gmap%nVx ) exit
                if( ix.eq.nx-1) then
                    iVx = iVx + 1  !! Lower right corners
                    if (isymm.eq.0 .or. isymm.eq.1 .or. isymm.eq.2) then
                        grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                            &   objects_per_dimension(1)%object(iVx)%geometry(1) =  &
                            &   crx( ix, iy, 1)
                        grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                            &   objects_per_dimension(1)%object(iVx)%geometry(2) =  &
                            &   cry( ix, iy, 1)
                    else if (isymm.eq.3 .or. isymm.eq.4) then
                        grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                            &   objects_per_dimension(1)%object(iVx)%geometry(1) =  &
                            &   cry( ix, iy, 1)
                        grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                            &   objects_per_dimension(1)%object(iVx)%geometry(2) =  &
                            &   crx( ix, iy, 1)
                    end if
                    if( iVx .eq. gmap%nVx ) exit
                end if
                if( iy.eq.ny-1) then
                    iVx = iVx + 1  !! Upper left corners
                    if (isymm.eq.0 .or. isymm.eq.1 .or. isymm.eq.2) then
                        grid_ggd%space( SPACE_POLOIDALPLANE )%                       &
                            &   objects_per_dimension(1)%object( iVx )%geometry(1) = &
                            &   crx( ix, iy, 2 )
                        grid_ggd%space( SPACE_POLOIDALPLANE )%                       &
                            &   objects_per_dimension(1)%object( iVx )%geometry(2) = &
                            &   cry( ix, iy, 2 )
                    else if (isymm.eq.3 .or. isymm.eq.4) then
                        grid_ggd%space( SPACE_POLOIDALPLANE )%                       &
                            &   objects_per_dimension(1)%object( iVx )%geometry(1) = &
                            &   cry( ix, iy, 2 )
                        grid_ggd%space( SPACE_POLOIDALPLANE )%                       &
                            &   objects_per_dimension(1)%object( iVx )%geometry(2) = &
                            &   crx( ix, iy, 2 )
                    end if
                    if( iVx .eq. gmap%nVx ) exit
                    if( ix.eq.nx-1) then
                        iVx = iVx + 1  !! Upper right corners
                        if (isymm.eq.0 .or. isymm.eq.1 .or. isymm.eq.2) then
                            grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                                &   objects_per_dimension(1)%object( iVx )% &
                                &   geometry(1) = crx( ix, iy, 3 )
                            grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                                &   objects_per_dimension(1)%object( iVx )% &
                                &   geometry(2) = cry( ix, iy, 3 )
                        else if (isymm.eq.3 .or. isymm.eq.4) then
                            grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                                &   objects_per_dimension(1)%object( iVx )% &
                                &   geometry(1) = cry( ix, iy, 3 )
                            grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                                &   objects_per_dimension(1)%object( iVx )% &
                                &   geometry(2) = crx( ix, iy, 3 )
                        end if
                        if( iVx .eq. gmap%nVx ) exit
                    end if
                end if
            end do
        end do
#endif

        !! Fill in object definitions (i.e. what objects compose an object)
        !! 1D objects: faces/edges
        nFc = gmap%nFcx + gmap%nFcy !! Number of all faces/edges
        !! Allocate 1D objects
        allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(2)%object( gmap%nFcx + gmap%nFcy ) )
        !! Each 1D object has two boundaries and two neighbours, in positive
        !! and negative coordinate direction, one on each side (for x-aligned
        !! edges: along flux surface, for y-aligned edges: orthogonal to flux
        !! surface)
        do iFc = 1, nFc
            !! Allocate and set all boundary & connectivity information to
            !! undefined
            allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(2)%object( iFc )%boundary(2) )
            !! Allocate list of 0D objects forming the 1D object
            !! Two 0D objects (vertices/nodes) form one 1D object (face/edge)
            allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(2)%object( iFc )%nodes(2) )
            do i = 1, 2
                !! Boundary to undefined
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(i)%index = B2_GRID_UNDEFINED
                !! Neighbours to undefined
                allocate(grid_ggd%space( SPACE_POLOIDALPLANE )% &
                    &   objects_per_dimension(2)%object( iFc )%boundary(i)% &
                    &   neighbours(2))
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(i)%neighbours(1) =   &
                    &   B2_GRID_UNDEFINED
                !! Nodes to undefined
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%nodes(i) = B2_GRID_UNDEFINED
            end do
        end do

        !! x-aligned edges
        do iFc = 1, gmap%nFcx
            !! get position of this face in the B2 grid
            ix = gmap%mapFcix( iFc )
            iy = gmap%mapFciy( iFc )
            !! get index of start vertex
            !! objdef dims: index of edge, 1=start node, 1=one-dimensional
            !! object
            select case ( gmap%mapFcIFace( iFc ) )
            case( BOTTOM )
                !! start index: 1=start node
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(1)%index =     &
                    &   gmap%mapVxI( ix, iy, VX_LOWERLEFT )
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%nodes(1) = gmap%mapVxI( ix, iy, VX_LOWERLEFT )
                if (gmap%mapVxI( ix, iy, VX_LOWERLEFT ) ==  &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: BOTTOM edge at (" &
                        &   //int2str(ix)//","//int2str(iy)//                &
                        &   ") has no start node")
                end if
                !! end vertex: 2=end node
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(2)%index =     &
                    &   gmap%mapVxI( ix, iy, VX_LOWERRIGHT )
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%nodes(2) = gmap%mapVxI( ix, iy, VX_LOWERRIGHT )
                if (gmap%mapVxI( ix, iy, VX_LOWERRIGHT ) == &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: BOTTOM edge at (" &
                        &   //int2str(ix)//","//int2str(iy)//                &
                        &   ") has no end node")
                end if
            case( TOP )
                !! start index: 1=start node
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(1)%index =     &
                    &   gmap%mapVxI( ix, iy, VX_UPPERLEFT )
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%nodes(1) = gmap%mapVxI( ix, iy, VX_UPPERLEFT )
                if (gmap%mapVxI( ix, iy, VX_UPPERLEFT ) ==  &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: TOP edge at ("    &
                        &   //int2str(ix)//","//int2str(iy)//                &
                        &   ") has no start node")
                end if
                !! end vertex: 2=end node
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(2)%index =     &
                    &   gmap%mapVxI( ix, iy, VX_UPPERRIGHT )
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%nodes(2) = gmap%mapVxI( ix, iy, VX_UPPERRIGHT )
                if (gmap%mapVxI( ix, iy, VX_UPPERRIGHT ) == &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: TOP edge at ("    &
                        &   //int2str(ix)//","//int2str(iy)//                &
                        &   ") has no end node")
                end if
            end select

            !! Neighbour edges of this edge
            !! Left neighbour: edge continuing to the left of this edge
            nix = leftix( ix, iy )
            niy = leftiy( ix, iy )
            if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(1)%neighbours(1) =    &
                    &   gmap%mapFcI( nix, niy, gmap%mapFcIFace( iFc ) )
            end if
            !! Right neighbour: edge continuing to the right of this edge
            nix = rightix( ix, iy )
            niy = rightiy( ix, iy )
            if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(2)%neighbours(2) =     &
                    &   gmap%mapFcI( nix, niy, gmap%mapFcIFace( iFc ) )
            end if

            !! 1d object measure: face area
            if (present(gs)) then
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%measure = &
                    &   gs(ix, iy, ALIGNX)
            end if
        end do

        !! y-aligned edges
        do iFc = gmap%nFcx + 1, gmap % nFcx + gmap % nFcy
            !! get position of this face in the B2 grid
            ix = gmap % mapFcix( iFc )
            iy = gmap % mapFciy( iFc )
!!$            if (gmap%mapCvI(ix, iy) == B2_GRID_UNDEFINED) then
!!$                call logmsg(LOGWARNING,
!!$                "b2IMASFillGD: writing out edges for unused cell ("//
!!$                int2str(ix)//","//int2str(iy)//")")
!!$            end if

            select case ( gmap % mapFcIFace( iFc ) )
            case( LEFT )
                !! start index: 1=start node
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(1)%index =                       &
                    &   gmap%mapVxI( ix, iy, VX_LOWERLEFT )
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%nodes(1) = gmap%mapVxI( ix, iy, VX_LOWERLEFT )
                if (gmap%mapVxI( ix, iy, VX_LOWERLEFT ) ==  &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: LEFT edge at ("   &
                        &   //int2str(ix)//","//int2str(iy)//                &
                        &   ") has no start node")
                end if
                !! end vertex: 2=end node
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(2)%index =                       &
                    &   gmap%mapVxI( ix, iy, VX_UPPERLEFT )
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%nodes(2) = gmap%mapVxI( ix, iy, VX_UPPERLEFT )
                if (gmap%mapVxI( ix, iy, VX_UPPERLEFT ) ==  &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: LEFT edge at ("   &
                        &   //int2str(ix)//","//int2str(iy)//                &
                        &   ") has no end node")
                end if

            case( RIGHT )
                !! start index: 1=start node
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(1)%index =                       &
                    &   gmap%mapVxI( ix, iy, VX_LOWERRIGHT )
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%nodes(1) = gmap%mapVxI( ix, iy, VX_LOWERRIGHT )
                if (gmap%mapVxI( ix, iy, VX_LOWERRIGHT ) == &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: RIGHT edge at ("  &
                        &   //int2str(ix)//","//int2str(iy)//                &
                        &   ") has no start node")
                end if
                !! end vertex: 2=end node
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(2)%index =                       &
                    &   gmap%mapVxI( ix, iy, VX_UPPERRIGHT )
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%nodes(2) = gmap%mapVxI( ix, iy, VX_UPPERRIGHT )
                if (gmap%mapVxI( ix, iy, VX_UPPERRIGHT ) == &
                    &    B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: RIGHT edge at ("  &
                        &   //int2str(ix)//","//int2str(iy)//                &
                        &   ") has no end node")
                end if
            end select

            !! Neighbour edges of this edge
            !! Bottom neighbour: edge continuing to the bottom of this edge
            nix = bottomix( ix, iy )
            niy = bottomiy( ix, iy )
            if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(1)%neighbours(1) =    &
                    &   gmap%mapFcI( nix, niy, gmap%mapFcIFace( iFc ) )
            end if
            !! Top neighbour: edge continuing to the top of this edge
            nix = topix( ix, iy )
            niy = topiy( ix, iy )
            if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%boundary(1)%neighbours(1) =    &
                    &   gmap%mapFcI( nix, niy, gmap%mapFcIFace( iFc ) )
            end if
            !! 1d object measure: edge length
            if (present(gs)) then
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( iFc )%measure = &
                    &   gs(ix, iy, ALIGNY)*qc(ix, iy)
            end if
        end do

        !! Fill in object definitions (i.e. what objects compose an object)
        !! 2D objects: Cells
        ! write(0,*) "num_obj_2D: ", gmap%nCv
        !! Allocate 2D objects
        allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
            &   object( gmap%nCv ) )

        !! Each 2D object has four boundaries
        !! Each boundary has one neighbour
        do iCv = 1, gmap%nCv
            !! Allocate and set all boundary & connectivity information to
            !! undefined
            allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%  &
                &   objects_per_dimension(3)%object( iCv )%boundary(4) )
            !! Allocate list of 0D objects forming the 2D object
            !! Four 0D objects (vertices/nodes) form one 2D object (2D cell)
            allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(3)%object( iCv )%nodes(4) )
            !! Also store additional geometry information: position in
            !! computational space
            !! FIXME: this should go into alternate geometry, which is not
            !!        available yet for grid objects
            allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%  &
                &   objects_per_dimension(3)%object( iCv )%geometry(2) )
            do i = 1, 4
                !! Boundary to undefined
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                    &   object( iCv )%boundary(i)%index = B2_GRID_UNDEFINED
                !! Neighbours to undefined
                allocate(grid_ggd%space( SPACE_POLOIDALPLANE )% &
                    &   objects_per_dimension(3)%object( iCv )%boundary(i)% &
                    &   neighbours(2))
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                    &   object( iCv )%boundary(i)%neighbours(1) =   &
                    &   B2_GRID_UNDEFINED
                !! Nodes to undefined
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                    &   object( iCv )%nodes(i) = B2_GRID_UNDEFINED
                !! Geometry to undefined
                if (i < 3) then
                    grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                        &   object( iCv )%geometry(i) = B2_GRID_UNDEFINED
                end if
            end do
        end do

        do iCv = 1, gmap%nCv
            !! Set position in computational space
            ix = gmap%mapCvix( iCv )
            iy = gmap%mapCviy( iCv )
            grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)%   &
                &   object( iCv )%geometry(1) = ix
            grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)%   &
                &   object( iCv )%geometry(2) = iy

            !! Set edges composing the quadrilateral in the list:
            !! left edge (y-aligned)
            grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)%   &
                &   object( iCv )%boundary(1)%index = gmap%mapFcI( ix, iy, LEFT )
            !! bottom edge (x-aligned)
            grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)%   &
                &   object( iCv )%boundary(2)%index = gmap%mapFcI( ix, iy, BOTTOM )
            !! right edge (y-aligned)
            grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)%   &
                &   object( iCv )%boundary(3)%index = gmap%mapFcI( ix, iy, RIGHT )
            !! top edge (x-aligned)
            grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)%   &
                &   object( iCv )%boundary(4)%index = gmap%mapFcI( ix, iy, TOP )
            do dir = LEFT, TOP
                call get_Neighbour(nx, ny, leftix, leftiy, rightix, rightiy,     &
                    &   topix, topiy, bottomix, bottomiy, ix, iy, dir, nix, niy)
                if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells,   &
                    &   nix, niy ) ) then
                    grid_ggd%space( SPACE_POLOIDALPLANE )%          &
                        &   objects_per_dimension(3)%object( iCv )% &
                        &   boundary( dir + 1 )%neighbours(1) =     &
                        &   gmap%mapCvI( nix, niy )
                end if
            end do
            !! 2d object measure: cell area
            if (present(vol)) then
                grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)%   &
                    &   object( iCv )%measure = vol(ix, iy)
            end if
        end do

        !! Set nodes list, composing the 2D objects - Cells, using a subroutine
        call set_Cells_Conn_Array_Nodes(grid_ggd)

#if 0
        !! TODO
        !! Fill in x-point indices
        !! In edge_profiles no node for data on x-points was found. There is
        !! hovewer one in equilibrium%boundary%x_point
        allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%xpoints( gmap%nsv ) )
        grid_ggd%space( SPACE_POLOIDALPLANE )%xpoints = gmap%svi(1:gmap%nsv)
#endif

        !! If requested, add a second space for the toroidal angle
        if (SPACE_COUNT == SPACE_TOROIDALANGLE) then
#if IMAS_MINOR_VERSION > 19
          allocate( grid_ggd%space( SPACE_TOROIDALANGLE )%identifier%name(1) )
          allocate( grid_ggd%space( SPACE_TOROIDALANGLE )%identifier%description(1) )
          grid_ggd%space( SPACE_TOROIDALANGLE )%identifier%index = 1
#endif
          grid_ggd%space( SPACE_TOROIDALANGLE )%geometry_type%index = 0
          allocate( grid_ggd%space( SPACE_TOROIDALANGLE )%geometry_type%name(1) )
          allocate( grid_ggd%space( SPACE_TOROIDALANGLE )%geometry_type%description(1) )
          if (isymm.eq.0) then
            call ipgetr ('b2agmt_1d_width', width)
            call gridSetupStruct1dSpace(                                      &
                &   grid_ggd%space( SPACE_TOROIDALANGLE ), COORDTYPE_Y,       &
                &   (/                                                        &
                &   ( ( width/NNODES_TOROIDAL )*i, i=0, NNODES_TOROIDAL )     &
                &   /),                                                       &
                &   .false. ) !! periodic = .false.
#if IMAS_MINOR_VERSION > 19
            grid_ggd%space( SPACE_TOROIDALANGLE )%identifier%name = "Z-direction"
            grid_ggd%space( SPACE_TOROIDALANGLE )%identifier%description =    &
                &   "Cylindrical symmetry"
#endif
            grid_ggd%space( SPACE_TOROIDALANGLE )%geometry_type%name = "Cylindrical"
            grid_ggd%space( SPACE_TOROIDALANGLE )%geometry_type%description = &
                &   "Length along cylindrical direction"
          else
#if IMAS_MINOR_VERSION > 19
            grid_ggd%space( SPACE_TOROIDALANGLE )%identifier%name = "Toroidal direction"
#endif
            grid_ggd%space( SPACE_TOROIDALANGLE )%geometry_type%name = "Toroidal angle"
            if ( TOROIDAL_PERIODIC ) then
              call gridSetupStruct1dSpace(                                      &
                  &   grid_ggd%space( SPACE_TOROIDALANGLE ), COORDTYPE_PHI,     &
                  &   (/                                                        &
                  &   ( ( 2*B2_PI/NNODES_TOROIDAL )*i, i=0, NNODES_TOROIDAL-1 ) &
                  &   /),                                                       &
                  &   .true. ) !! periodic = .true.
#if IMAS_MINOR_VERSION > 19
              grid_ggd%space( SPACE_TOROIDALANGLE )%identifier%description =    &
                  &   "Toroidally symmetric and periodic"
#endif
              grid_ggd%space( SPACE_TOROIDALANGLE )%geometry_type%description = &
                  &   "Toroidal angle"
            else
              call gridSetupStruct1dSpace(                                      &
                  &   grid_ggd%space( SPACE_TOROIDALANGLE ), COORDTYPE_PHI,     &
                  &   (/                                                        &
                  &   ( ( 2*B2_PI/NNODES_TOROIDAL )*i, i=0, NNODES_TOROIDAL )   &
                  &   /),                                                       &
                  &   .false. ) !! periodic = .false.
#if IMAS_MINOR_VERSION > 19
              grid_ggd%space( SPACE_TOROIDALANGLE )%identifier%description =    &
                  &   "Toroidally symmetric"
#endif
              grid_ggd%space( SPACE_TOROIDALANGLE )%geometry_type%description = &
                  &   "Toroidal angle, full circle"
            end if
          end if
        end if

    end subroutine fill_In_Grid_Desc

    !> Set connectivity array for cells by defining nodes that form each cell
    subroutine set_Cells_Conn_Array_Nodes(grid_ggd)
#if IMAS_MINOR_VERSION < 15
        type(ids_generic_grid_dynamic), intent(inout) :: grid_ggd !< Type of IDS
            !< data structure, designed for handling grid geometry data
#else
        type(ids_generic_grid_aos3_root), intent(inout) :: grid_ggd !< Type of IDS
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
        integer :: iCv  !< Cell index
        integer :: ix   !< x-aligned (poloidal) cell index
        integer :: iy   !< y-aligned (radial) cell index
        integer :: m    !< Iterator
        integer :: loop_count   !< Loop counter


        !! Get the list of 0D objects forming the 2D objects and write it to IDS
        allocate( objects2Darray( gmap%nCv, 4 ) )
        ! objects2Darray = numpy.array([], dtype='int')

        ! ids_dim_2D.object.resize(nCv)

        ! Already done
        if (.not.associated( grid_ggd%space( SPACE_POLOIDALPLANE )%  &
            &   objects_per_dimension(3)%object ) )                  &
            & allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%       &
            &   objects_per_dimension(3)%object( gmap%nCv ) )

        !! Get 2D objects geometry (nodes) data from CPO, sort them into more
        !! orderly form  and put them into the IDS
        num_nodes_2D    = 4
        num_boundary_2D = 4
        allocate( node_idx(4) )
        allocate( free_edge(3) )

        do iCv = 1, gmap%nCv
            ix = gmap%mapCvix( iCv )
            iy = gmap%mapCviy( iCv )
            if (.not.associated( grid_ggd%space( SPACE_POLOIDALPLANE )%  &
                &   objects_per_dimension(3)%object(iCv)%nodes ) )       &
                & allocate( grid_ggd%space( SPACE_POLOIDALPLANE )%       &
                &   objects_per_dimension(3)%object(iCv)%nodes(num_nodes_2D) )

            edge_idx = gmap%mapFcI( ix, iy, LEFT )
            free_edge(1) = gmap%mapFcI( ix, iy, BOTTOM )
            free_edge(2) = gmap%mapFcI( ix, iy, RIGHT )
            free_edge(3) = gmap%mapFcI( ix, iy, TOP )

            node_idx(1) = grid_ggd%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(2)%object( edge_idx )%    &
                &   boundary(1)%index
            last_idx = 2
            node_idx(last_idx) = grid_ggd%space( SPACE_POLOIDALPLANE )% &
                &   objects_per_dimension(2)%object( edge_idx )%        &
                &   boundary(2)%index

            do loop_count = 1, 4
                if (last_idx < 4) then
                    do m = 1, 3
                        edge_idx = free_edge(m)
                        if (edge_idx < 1) then
                            continue
                        end if
                        node1 = grid_ggd%space( SPACE_POLOIDALPLANE )%   &
                            &   objects_per_dimension(2)%object( edge_idx )% &
                            &   boundary(1)%index
                        node2 = grid_ggd%space( SPACE_POLOIDALPLANE )%   &
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
            grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                &   object( iCv )%nodes(1) = node_idx(1)
            grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                &   object( iCv )%nodes(2) = node_idx(4)
            grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                &   object( iCv )%nodes(3) = node_idx(3)
            grid_ggd%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                &   object( iCv )%nodes(4) = node_idx(2)
        end do

    end subroutine set_Cells_Conn_Array_Nodes

    !> Define grid subsets
    subroutine fill_In_GridSubset_Desc
        !! Internal variables
        integer :: geoId
#if GGD_MINOR_VERSION > 8
        integer :: iRegion
        integer :: iPrivateB2
#endif
        integer :: GSubsetCount
        integer :: iType
        integer :: RegionsInSubset(6)
        integer :: nGSubset !< Total number of grid subsets
        integer :: nInd     !< Size of grid subset element list
        integer :: xIn
        integer :: yIn
        integer :: xOut
        integer :: yOut
        integer :: iCoreGS !< Core grid subset ID
        integer :: cls(SPACE_COUNT_MAX)
        integer, allocatable :: xpoints(:,:)
        integer, allocatable :: indexList1d(:)
        integer, dimension(:,:), allocatable :: indexList2d, indexPart2d, indextmp2d
        integer :: i, j    !< Iterators
        integer :: ix, iy  !< Grid coordinates
        integer :: iVx     !< Vertex/node index
        integer :: iSubset !< Iterator on standard subsets
        integer :: ireg    !< Iterator on regions included in subset
        integer :: ind     !< indexList2d start index
        integer :: iInd    !< indexList iterator
        integer :: isize
        integer :: jsep, nxtl, nxtr
        character*128 RegionDescription

        geoId = geometryId(nnreg, isymm, periodic_bc, topcut)
        call get_jsep( nx, ny, jxi, jxa, jsep )
        call get_nxt ( nx, nxtl, nxtr )

        !! Figure out total number of grid subsets
        !! Do generic + private grid subsets
#if GGD_MINOR_VERSION > 8
        nGSubset = B2_GENERIC_GSUBSET_COUNT + regionCountTotal(geoId)
#else
        nGSubset = B2_GENERIC_GSUBSET_COUNT
#endif
        !! Add pre-defined grid subsets (regions + points)
        select case ( geoId )
        case ( GEOMETRY_LINEAR )
            nGSubset = nGSubset + 5 + 2
            if ( jsep .ne. B2_GRID_UNDEFINED ) nGSubset = nGSubset + 2
        case ( GEOMETRY_CYLINDER, GEOMETRY_ANNULUS )
            nGSubset = nGSubset + 3
        case ( GEOMETRY_LIMITER )
            nGSubset = nGSubset + 10 + 2
        case ( GEOMETRY_SN )
            nGSubset = nGSubset + 20 + 2
        case ( GEOMETRY_STELLARATORISLAND )
            nGSubset = nGSubset + 14 + 2
        case ( GEOMETRY_CDN)
            nGSubset = nGSubset + 32 + 4
        case ( GEOMETRY_DDN_TOP, GEOMETRY_DDN_BOTTOM )
            nGSubset = nGSubset + 33 + 4
        end select
        !! Inner/outer midplane grid subsets
        nGSubset = nGSubset + 2
        !! Inner/outer midplane separatrix
        if (jsep /= B2_GRID_UNDEFINED) nGSubset = nGSubset + 2

        call logmsg( LOGDEBUG, "b2_IMAS_Fill_Grid_Desc: expecting total of " &
            &//int2str(nGSubset)//" grid subsets" )
        allocate( grid_ggd%grid_subset( nGSubset ) )

        !! Set up generic grid subsets
        !! The 6 generic grid subsets MUST be declared in that order
        !! matching the order in ids_grid_common
        !! GRID_SUBSET_NODES: all nodes, one implicit object list

        call createGridSubsetForClass( grid_ggd,                    &
            &   grid_ggd%grid_subset( GRID_SUBSET_NODES ),          &
            &   IDS_CLASS_NODE, SPACE_POLOIDALPLANE,                &
            &   GRID_SUBSET_NODES,                                  &
            &   "Nodes", "All nodes (0D objects) in the domain." )

        !! GRID_SUBSET_EDGES: all edges, one implicit object list
        call createGridSubsetForClass( grid_ggd,                    &
            &   grid_ggd%grid_subset( GRID_SUBSET_EDGES ),          &
            &   IDS_CLASS_POLOIDALRADIAL_EDGE, SPACE_POLOIDALPLANE, &
            &   GRID_SUBSET_EDGES,                                  &
            &   "Edges", "All edges (1D objects) in the domain." )

        !! GRID_SUBSET_X_ALIGNED_EDGES: x-aligned edges. One implicit object
        !! list, range over x edges
        !! Create grid subset with one object list
        call createEmptyGridSubset(                                    &
            &   grid_ggd%grid_subset( GRID_SUBSET_X_ALIGNED_EDGES ),   &
            &   GRID_SUBSET_X_ALIGNED_EDGES, 'x-aligned edges',        &
            &   "All X-aligned edges (1D objects) in the domain." )
        !! Initialize implicit object list for edges (class (/2/) )
        allocate(indexList1d(gmap%nFcx))
        indexList1d = (/ (i, i = 1, gmap%nFcx) /)
        call createExplicitObjectListSingleSpace( grid_ggd,            &
            &   grid_ggd%grid_subset( GRID_SUBSET_X_ALIGNED_EDGES ),   &
            &   IDS_CLASS_POLOIDALRADIAL_EDGE, indexList1d,            &
            &   IDS_CLASS_POLOIDALRADIAL_EDGE, SPACE_POLOIDALPLANE )

        !! GRID_SUBSET_Y_ALIGNED_EDGES: y-aligned edges. One implicit object
        !! list, range over y edges
        !! Create grid subset with one object list
        call createEmptyGridSubset(                                    &
            &   grid_ggd%grid_subset( GRID_SUBSET_Y_ALIGNED_EDGES ),   &
            &   GRID_SUBSET_Y_ALIGNED_EDGES, 'y-aligned edges',        &
            &   "All Y-aligned edges (1D objects) in the domain." )
        !! Initialize implicit object list for edges (class (/2/) )
        deallocate(indexList1d)
        allocate(indexList1d(gmap%nFcy))
        indexList1d = (/ (i, i = gmap%nFcx + 1, gmap%nFcx + gmap%nFcy) /)
        call createExplicitObjectListSingleSpace( grid_ggd,             &
            &   grid_ggd%grid_subset( GRID_SUBSET_Y_ALIGNED_EDGES ),    &
            &   IDS_CLASS_POLOIDALRADIAL_EDGE, indexList1d,             &
            &   IDS_CLASS_POLOIDALRADIAL_EDGE, SPACE_POLOIDALPLANE )

        !! GRID_SUBSET_CELLS: all 2d cells, one implicit object list
        call createGridSubsetForClass( grid_ggd,                &
            &   grid_ggd%grid_subset( GRID_SUBSET_CELLS ),      &
            &   IDS_CLASS_CELL, 1, GRID_SUBSET_CELLS, "Cells",  &
            &   "All cells (2D objects) in the domain." )

        !! Grid subset of all x-points
        !! (in one poloidal plane at toroidal index 1)
        !! Assemble object descriptor for x-points
        allocate( xpoints(gmap%nsv, SPACE_COUNT) )
        xpoints = 1
        xpoints(:, SPACE_POLOIDALPLANE) = gmap%svi(1:gmap%nsv)
        !! Create grid subset with one object list
        call createEmptyGridSubset(                                  &
            &   grid_ggd%grid_subset( GRID_SUBSET_X_POINTS ),        &
            &   GRID_SUBSET_X_POINTS, 'x-points',                    &
            &   "All X-points (0D objects) in the domain." )
        !! Initialize explicit object list for edges (class (/1/) )
        call createExplicitObjectListSingleSpace( grid_ggd,          &
                &   grid_ggd%grid_subset( GRID_SUBSET_X_POINTS ),    &
                &   IDS_CLASS_NODE, xpoints(:, SPACE_POLOIDALPLANE), &
                &   IDS_CLASS_NODE, SPACE_POLOIDALPLANE )

        !! Set up specific grid subset by collecting edges for regions

        !! Start counting from end of generic grid subset
        GSubsetCount = B2_GENERIC_GSUBSET_COUNT

#if GGD_MINOR_VERSION > 8
        iPrivateB2 = 0
        !! Cell + edge grid subset
        !! These are the "private" B2 regions, so will be given negative
        !! grid subset identifiers
        do iType = REGIONTYPE_CELL, REGIONTYPE_YEDGE

            select case(iType)
            case( REGIONTYPE_CELL )
                cls = CLASS_CELL
            case( REGIONTYPE_YEDGE, REGIONTYPE_XEDGE )
                cls = CLASS_POLOIDALRADIAL_EDGE
            end select

            do iRegion = 1, regionCount(geoId, iType)
                iPrivateB2 = iPrivateB2 - 1
                GSubsetCount = GSubsetCount + 1
                select case(iType)
                case( REGIONTYPE_CELL )
                  RegionDescription = "Volumetric B2.5 internal region #"// &
                    &   int2str(iRegion)
                case( REGIONTYPE_XEDGE )
                  RegionDescription = "Y-aligned B2.5 internal region #"//  &
                    &   int2str(iRegion)
                case( REGIONTYPE_YEDGE )
                  RegionDescription = "X-aligned B2.5 internal region #"//  &
                    &   int2str(iRegion)
                end select

                call logmsg( LOGDEBUG, "b2_IMAS_Fill_Grid_Desc:"// &
                    &   " add (private) grid subset #"//           &
                    &   int2str(GSubsetCount)//                    &
                    &   " for iType "//int2str( iType )//          &
                    &   ", iRegion "//int2str( iRegion )//": "//   &
                    &   regionName(geoId, iType, iRegion) )

                !! Create grid subset with one object list
                call createEmptyGridSubset(                              &
                    &   grid_ggd%grid_subset( GSubsetCount ),            &
                    &   iPrivateB2, regionName( geoId, iType, iRegion ), &
                    &   RegionDescription )

                !! Get explicit object list of the grid subset using
                !! subroutine collectIndexListForRegionSubroutine
                !! (function collectIndexListForRegion transferred to subroutine,
                !! as array of certain dimension is required as an output)
                call collectIndexListForRegionSubroutine( gmap, cflag, region, &
                    &   iType, iRegion, indexList2d )

                !! Initialize explicit object list for grid subset
                call createExplicitObjectListSingleSpace( grid_ggd,     &
                    &   grid_ggd%grid_subset( GSubsetCount ), sum(cls), &
                    &   indexList2d(:,SPACE_POLOIDALPLANE), sum(cls),   &
                    &   SPACE_POLOIDALPLANE )

            end do
        end do
        deallocate(indexList2d)
#endif

!! Do the grid subsets that map directly to B2 regions
        do iSubset = GRID_SUBSET_CORE_CUT, GRID_SUBSET_INNER_TARGET_INACTIVE
            if (iSubset == GRID_SUBSET_OUTER_MIDPLANE ) cycle ! Already handled below
            if (iSubset == GRID_SUBSET_INNER_MIDPLANE ) cycle ! Already handled below
            if (iSubset == GRID_SUBSET_BETWEEN_SEPARATRICES ) cycle ! Handled below

            select case( iSubset )
            case ( GRID_SUBSET_CORE, GRID_SUBSET_SOL, &
                & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
                & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
                iType = REGIONTYPE_CELL
            case ( GRID_SUBSET_CORE_BOUNDARY, GRID_SUBSET_SEPARATRIX, &
                & GRID_SUBSET_ACTIVE_SEPARATRIX, GRID_SUBSET_MAIN_CHAMBER_WALL, &
                & GRID_SUBSET_OUTER_BAFFLE, GRID_SUBSET_INNER_BAFFLE, &
                & GRID_SUBSET_OUTER_PFR_WALL, GRID_SUBSET_INNER_PFR_WALL, &
                & GRID_SUBSET_MAIN_WALL, GRID_SUBSET_PFR_WALL, &
                & GRID_SUBSET_SECOND_SEPARATRIX, &
                & GRID_SUBSET_OUTER_BAFFLE_INACTIVE, GRID_SUBSET_INNER_BAFFLE_INACTIVE, &
                & GRID_SUBSET_OUTER_PFR_WALL_INACTIVE, GRID_SUBSET_INNER_PFR_WALL_INACTIVE )
                iType = REGIONTYPE_YEDGE
            case ( GRID_SUBSET_CORE_CUT, GRID_SUBSET_PFR_CUT, &
                & GRID_SUBSET_OUTER_THROAT, GRID_SUBSET_INNER_THROAT, &
                & GRID_SUBSET_OUTER_TARGET, GRID_SUBSET_INNER_TARGET, &
                & GRID_SUBSET_CORE_CUT_INACTIVE, GRID_SUBSET_PFR_CUT_INACTIVE, &
                & GRID_SUBSET_OUTER_THROAT_INACTIVE, GRID_SUBSET_INNER_THROAT_INACTIVE, &
                & GRID_SUBSET_OUTER_TARGET_INACTIVE, GRID_SUBSET_INNER_TARGET_INACTIVE )
                iType = REGIONTYPE_XEDGE
            end select

            select case(iType)
            case( REGIONTYPE_CELL )
                cls = CLASS_CELL
            case( REGIONTYPE_YEDGE, REGIONTYPE_XEDGE )
                cls = CLASS_POLOIDALRADIAL_EDGE
            end select

            RegionsinSubset = 0
            select case( geoId )
            case ( GEOMETRY_LINEAR )
                select case( iSubset )
                case ( GRID_SUBSET_CORE_BOUNDARY )
                    RegionsinSubset(1) = 1
                case ( GRID_SUBSET_MAIN_CHAMBER_WALL, GRID_SUBSET_MAIN_WALL )
                    RegionsinSubset(1) = 2
                case ( GRID_SUBSET_OUTER_TARGET )
                    RegionsinSubset(1) = 2
                case ( GRID_SUBSET_INNER_TARGET )
                    RegionsinSubset(1) = 1
                end select
            case ( GEOMETRY_CYLINDER, GEOMETRY_ANNULUS )
                select case( iSubset )
                case ( GRID_SUBSET_CORE )
                    RegionsinSubset(1) = 1
                case ( GRID_SUBSET_CORE_BOUNDARY )
                    RegionsinSubset(1) = 1
                case ( GRID_SUBSET_CORE_CUT )
                    RegionsinSubset(1) = 1
                end select
            case ( GEOMETRY_LIMITER )
                select case( iSubset )
                case ( GRID_SUBSET_CORE )
                    RegionsinSubset(1) = 1
                case ( GRID_SUBSET_SOL )
                    RegionsinSubset(1) = 2
                case ( GRID_SUBSET_CORE_BOUNDARY )
                    RegionsinSubset(1) = 1
                case ( GRID_SUBSET_SEPARATRIX, GRID_SUBSET_ACTIVE_SEPARATRIX )
                    RegionsinSubset(1) = 2
                case ( GRID_SUBSET_MAIN_CHAMBER_WALL, GRID_SUBSET_MAIN_WALL )
                    RegionsinSubset(1) = 3
                case ( GRID_SUBSET_OUTER_TARGET )
                    if (LSN) then
                      RegionsinSubset(1) = 2
                    else
                      RegionsinSubset(1) = 1
                    end if
                case ( GRID_SUBSET_INNER_TARGET )
                    if (LSN) then
                      RegionsinSubset(1) = 1
                    else
                      RegionsinSubset(1) = 2
                    end if
                case ( GRID_SUBSET_CORE_CUT )
                    RegionsinSubset(1) = 3
                end select
            case ( GEOMETRY_SN )
                select case( iSubset )
                case ( GRID_SUBSET_CORE )
                    RegionsinSubset(1) = 1
                case ( GRID_SUBSET_SOL )
                    RegionsinSubset(1) = 2
                case ( GRID_SUBSET_OUTER_DIVERTOR )
                    if (LSN) then
                        RegionsinSubset(1) = 4
                    else
                        RegionsinSubset(1) = 3
                    end if
                case ( GRID_SUBSET_INNER_DIVERTOR )
                    if (LSN) then
                        RegionsinSubset(1) = 3
                    else
                        RegionsinSubset(1) = 4
                    end if
                case ( GRID_SUBSET_CORE_BOUNDARY )
                    RegionsinSubset(1) = 2
                case ( GRID_SUBSET_ACTIVE_SEPARATRIX )
                    RegionsinSubset(1) = 4
                case ( GRID_SUBSET_MAIN_CHAMBER_WALL )
                    RegionsinSubset(1) = 6
                case ( GRID_SUBSET_OUTER_BAFFLE )
                    if (LSN) then
                        RegionsinSubset(1) = 7
                    else
                        RegionsinSubset(1) = 5
                    end if
                case ( GRID_SUBSET_INNER_BAFFLE )
                    if (LSN) then
                        RegionsinSubset(1) = 5
                    else
                        RegionsinSubset(1) = 7
                    end if
                case ( GRID_SUBSET_OUTER_PFR_WALL )
                    if (LSN) then
                        RegionsinSubset(1) = 3
                    else
                        RegionsinSubset(1) = 1
                    end if
                case ( GRID_SUBSET_INNER_PFR_WALL )
                    if (LSN) then
                        RegionsinSubset(1) = 1
                    else
                        RegionsinSubset(1) = 3
                    end if
                case ( GRID_SUBSET_MAIN_WALL )
                    RegionsinSubset(1) = 5
                    RegionsinSubset(2) = 6
                    RegionsinSubset(3) = 7
                case ( GRID_SUBSET_PFR_WALL )
                    RegionsinSubset(1) = 1
                    RegionsinSubset(2) = 3
                case ( GRID_SUBSET_CORE_CUT )
                    RegionsinSubset(1) = 5
                case ( GRID_SUBSET_PFR_CUT )
                    RegionsinSubset(1) = 6
                case ( GRID_SUBSET_OUTER_THROAT )
                    if (LSN) then
                        RegionsinSubset(1) = 3
                    else
                        RegionsinSubset(1) = 2
                    end if
                case ( GRID_SUBSET_INNER_THROAT )
                    if (LSN) then
                        RegionsinSubset(1) = 2
                    else
                        RegionsinSubset(1) = 3
                    end if
                case ( GRID_SUBSET_OUTER_TARGET )
                    if (LSN) then
                        RegionsinSubset(1) = 4
                    else
                        RegionsinSubset(1) = 1
                    end if
                case ( GRID_SUBSET_INNER_TARGET )
                    if (LSN) then
                        RegionsinSubset(1) = 1
                    else
                        RegionsinSubset(1) = 4
                    end if
                end select
            case ( GEOMETRY_STELLARATORISLAND )
                select case( iSubset )
                case ( GRID_SUBSET_CORE )
                    RegionsinSubset(1) = 1
                case ( GRID_SUBSET_SOL )
                    RegionsinSubset(1) = 2
                case ( GRID_SUBSET_CORE_BOUNDARY )
                    RegionsinSubset(1) = 2
                case ( GRID_SUBSET_SEPARATRIX )
                    RegionsinSubset(1) = 5
                    RegionsinSubset(2) = 4
                    RegionsinSubset(3) = 7
                case ( GRID_SUBSET_ACTIVE_SEPARATRIX )
                    RegionsinSubset(1) = 4
                case ( GRID_SUBSET_OUTER_PFR_WALL )
                    RegionsinSubset(1) = 3
                case ( GRID_SUBSET_INNER_PFR_WALL )
                    RegionsinSubset(1) = 1
                case ( GRID_SUBSET_PFR_WALL )
                    RegionsinSubset(1) = 1
                    RegionsinSubset(2) = 3
                case ( GRID_SUBSET_CORE_CUT )
                    RegionsinSubset(1) = 5
                case ( GRID_SUBSET_PFR_CUT )
                    RegionsinSubset(1) = 6
                case ( GRID_SUBSET_OUTER_THROAT )
                    RegionsinSubset(1) = 3
                case ( GRID_SUBSET_INNER_THROAT )
                    RegionsinSubset(1) = 2
                case ( GRID_SUBSET_OUTER_TARGET )
                    RegionsinSubset(1) = 4
                case ( GRID_SUBSET_INNER_TARGET )
                    RegionsinSubset(1) = 1
                end select
            case ( GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM )
                select case( iSubset )
                case ( GRID_SUBSET_CORE )
                    RegionsinSubset(1) = 1
                    RegionsinSubset(2) = 5
                case ( GRID_SUBSET_SOL )
                    RegionsinSubset(1) = 2
                    RegionsinSubset(2) = 6
                case ( GRID_SUBSET_OUTER_DIVERTOR )
                    RegionsinSubset(1) = 8
                case ( GRID_SUBSET_INNER_DIVERTOR )
                    RegionsinSubset(1) = 3
                case ( GRID_SUBSET_CORE_BOUNDARY )
                    RegionsinSubset(1) = 2
                    RegionsinSubset(2) = 9
                case ( GRID_SUBSET_ACTIVE_SEPARATRIX )
                    RegionsinSubset(1) = 4
                    RegionsinSubset(2) = 11
                case ( GRID_SUBSET_MAIN_CHAMBER_WALL )
                    RegionsinSubset(1) = 6
                    RegionsinSubset(2) = 13
                case ( GRID_SUBSET_OUTER_BAFFLE )
                    RegionsinSubset(1) = 14
                case ( GRID_SUBSET_INNER_BAFFLE )
                    RegionsinSubset(1) = 5
                case ( GRID_SUBSET_OUTER_PFR_WALL )
                    RegionsinSubset(1) = 10
                case ( GRID_SUBSET_INNER_PFR_WALL )
                    RegionsinSubset(1) = 1
                case ( GRID_SUBSET_MAIN_WALL )
                    RegionsinSubset(1) = 5
                    RegionsinSubset(2) = 6
                    RegionsinSubset(3) = 7
                    RegionsinSubset(4) = 12
                    RegionsinSubset(5) = 13
                    RegionsinSubset(6) = 14
                case ( GRID_SUBSET_PFR_WALL )
                    RegionsinSubset(1) = 1
                    RegionsinSubset(2) = 3
                    RegionsinSubset(3) = 8
                    RegionsinSubset(4) = 10
                case ( GRID_SUBSET_OUTER_BAFFLE_INACTIVE )
                    RegionsinSubset(1) = 12
                case ( GRID_SUBSET_INNER_BAFFLE_INACTIVE )
                    RegionsinSubset(1) = 7
                case ( GRID_SUBSET_OUTER_PFR_WALL_INACTIVE )
                    RegionsinSubset(1) = 8
                case ( GRID_SUBSET_INNER_PFR_WALL_INACTIVE )
                    RegionsinSubset(1) = 3
                case ( GRID_SUBSET_CORE_CUT )
                    RegionsinSubset(1) = 9
                case ( GRID_SUBSET_PFR_CUT )
                    RegionsinSubset(1) = 12
                case ( GRID_SUBSET_OUTER_THROAT )
                    RegionsinSubset(1) = 7
                case ( GRID_SUBSET_INNER_THROAT )
                    RegionsinSubset(1) = 2
                case ( GRID_SUBSET_OUTER_TARGET )
                    RegionsinSubset(1) = 8
                case ( GRID_SUBSET_INNER_TARGET )
                    RegionsinSubset(1) = 1
                case ( GRID_SUBSET_CORE_CUT_INACTIVE )
                    RegionsinSubset(1) = 11
                case ( GRID_SUBSET_PFR_CUT_INACTIVE )
                    RegionsinSubset(1) = 10
                case ( GRID_SUBSET_OUTER_THROAT_INACTIVE )
                    RegionsinSubset(1) = 6
                case ( GRID_SUBSET_INNER_THROAT_INACTIVE )
                    RegionsinSubset(1) = 3
                case ( GRID_SUBSET_OUTER_DIVERTOR_INACTIVE )
                    RegionsinSubset(1) = 7
                case ( GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
                    RegionsinSubset(1) = 4
                case ( GRID_SUBSET_OUTER_TARGET_INACTIVE )
                    RegionsinSubset(1) = 5
                case ( GRID_SUBSET_INNER_TARGET_INACTIVE )
                    RegionsinSubset(1) = 4
                end select
            case ( GEOMETRY_DDN_TOP )
                select case( iSubset )
                case ( GRID_SUBSET_CORE )
                    RegionsinSubset(1) = 1
                    RegionsinSubset(2) = 5
                case ( GRID_SUBSET_SOL )
                    RegionsinSubset(1) = 2
                    RegionsinSubset(2) = 6
                case ( GRID_SUBSET_OUTER_DIVERTOR )
                    RegionsinSubset(1) = 7
                case ( GRID_SUBSET_INNER_DIVERTOR )
                    RegionsinSubset(1) = 4
                case ( GRID_SUBSET_CORE_BOUNDARY )
                    RegionsinSubset(1) = 2
                    RegionsinSubset(2) = 9
                case ( GRID_SUBSET_ACTIVE_SEPARATRIX )
                    RegionsinSubset(1) = 4
                    RegionsinSubset(2) = 11
                case ( GRID_SUBSET_MAIN_CHAMBER_WALL )
                    RegionsinSubset(1) = 6
                    RegionsinSubset(2) = 13
                case ( GRID_SUBSET_OUTER_BAFFLE )
                    RegionsinSubset(1) = 12
                case ( GRID_SUBSET_INNER_BAFFLE )
                    RegionsinSubset(1) = 7
                case ( GRID_SUBSET_OUTER_PFR_WALL )
                    RegionsinSubset(1) = 8
                case ( GRID_SUBSET_INNER_PFR_WALL )
                    RegionsinSubset(1) = 3
                case ( GRID_SUBSET_MAIN_WALL )
                    RegionsinSubset(1) = 5
                    RegionsinSubset(2) = 6
                    RegionsinSubset(3) = 7
                    RegionsinSubset(4) = 12
                    RegionsinSubset(5) = 13
                    RegionsinSubset(6) = 14
                case ( GRID_SUBSET_PFR_WALL )
                    RegionsinSubset(1) = 1
                    RegionsinSubset(2) = 3
                    RegionsinSubset(3) = 8
                    RegionsinSubset(4) = 10
                case ( GRID_SUBSET_OUTER_BAFFLE_INACTIVE )
                    RegionsinSubset(1) = 14
                case ( GRID_SUBSET_INNER_BAFFLE_INACTIVE )
                    RegionsinSubset(1) = 5
                case ( GRID_SUBSET_OUTER_PFR_WALL_INACTIVE )
                    RegionsinSubset(1) = 10
                case ( GRID_SUBSET_INNER_PFR_WALL_INACTIVE )
                    RegionsinSubset(1) = 1
                case ( GRID_SUBSET_CORE_CUT )
                    RegionsinSubset(1) = 11
                case ( GRID_SUBSET_PFR_CUT )
                    RegionsinSubset(1) = 10
                case ( GRID_SUBSET_OUTER_THROAT )
                    RegionsinSubset(1) = 6
                case ( GRID_SUBSET_INNER_THROAT )
                    RegionsinSubset(1) = 3
                case ( GRID_SUBSET_OUTER_TARGET )
                    RegionsinSubset(1) = 5
                case ( GRID_SUBSET_INNER_TARGET )
                    RegionsinSubset(1) = 4
                case ( GRID_SUBSET_CORE_CUT_INACTIVE )
                    RegionsinSubset(1) = 9
                case ( GRID_SUBSET_PFR_CUT_INACTIVE )
                    RegionsinSubset(1) = 12
                case ( GRID_SUBSET_OUTER_THROAT_INACTIVE )
                    RegionsinSubset(1) = 7
                case ( GRID_SUBSET_INNER_THROAT_INACTIVE )
                    RegionsinSubset(1) = 2
                case ( GRID_SUBSET_OUTER_DIVERTOR_INACTIVE )
                    RegionsinSubset(1) = 8
                case ( GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
                    RegionsinSubset(1) = 3
                case ( GRID_SUBSET_OUTER_TARGET_INACTIVE )
                    RegionsinSubset(1) = 1
                case ( GRID_SUBSET_INNER_TARGET_INACTIVE )
                    RegionsinSubset(1) = 8
                end select
            case ( GEOMETRY_UNSPECIFIED )
                continue
            end select
            if (RegionsinSubset(1) == 0) cycle
            GSubsetCount = GSubsetCount + 1

            call logmsg( LOGDEBUG,                                     &
               &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//       &
               &   int2str(GSubsetCount)//": "//                       &
               &   trim(gridSubsetName ( iSubset ))//                  &
               &   ", iType "//int2str(iType) )

            !! Create grid subset with one object list
            call createEmptyGridSubset(                     &
               &   grid_ggd%grid_subset( GSubsetCount ),    &
               &   iSubset, gridSubsetName ( iSubset ),     &
               &   gridSubsetDescription( iSubset )  )

            !! Get explicit object list of the grid subset using
            !! subroutine collectIndexListForRegionSubroutine
            !! (function collectIndexListForRegion transferred to subroutine,
            !! as array of certain dimension is required as an output)
            select case ( iType )
            case ( REGIONTYPE_CELL )
                allocate( indextmp2d ( gmap%nCv , SPACE_COUNT ) )
            case ( REGIONTYPE_XEDGE )
                allocate( indextmp2d ( gmap%nFcx , SPACE_COUNT ) )
            case ( REGIONTYPE_YEDGE )
                allocate( indextmp2d ( gmap%nFcy , SPACE_COUNT ) )
            end select
            indextmp2d = 1
            isize = 0
            do ireg = 1, size(RegionsinSubset)
                if (RegionsinSubset(ireg) == 0) cycle
                call collectIndexListForRegionSubroutine( gmap, cflag, region,   &
                   &   iType, RegionsinSubset( ireg ), indexPart2d )
                indextmp2d( isize+1 : isize+size(indexPart2d,1),:) = indexPart2d(:,:)
                isize = isize + size(indexPart2d,1)
            end do
            allocate( indexList2d ( isize, SPACE_COUNT ) )
            indexList2d(1:isize,:) = indextmp2d(1:isize,:)

            !! Initialize explicit object list for grid subset
            call createExplicitObjectListSingleSpace( grid_ggd,     &
               &   grid_ggd%grid_subset( GSubsetCount ), sum(cls),  &
               &   indexList2d(:,SPACE_POLOIDALPLANE), sum(cls),    &
               &   SPACE_POLOIDALPLANE )
            deallocate(indexList2d,indexPart2d,indextmp2d)

        end do

        !! Add midplane node grid subsets
        !! Find the core boundary grid subset by looking for its name as
        !! defined in b2mod_connectivity
        iCoreGS = findGridSubsetByName(grid_ggd, gridSubsetName( GRID_SUBSET_CORE_BOUNDARY ) )
        !! For double null, we need the outer half of the core boundary
        if (iCoreGS == B2_GRID_UNDEFINED) then
            iCoreGS = findGridSubsetByName(grid_ggd, "Outer core boundary")
        end if
        if (iCoreGS == B2_GRID_UNDEFINED) stop "fill_In_GridSubset_Desc: "// &
            & "did not find core boundary grid subset for assembling " // &
            & " outer midplane grid subset"

        !! Figure out starting points for inner and outer midplane on core
        !! boundary
        call find_Midplane_Cells(grid_ggd%grid_subset( iCoreGS ), gmap, crx,  &
            &   xIn, yIn, xOut, yOut)

        GSubsetCount = GSubsetCount + 1
        !! Create grid subset with one object list
        call createEmptyGridSubset(                           &
            &   grid_ggd%grid_subset( GSubsetCount ),         &
            &   GRID_SUBSET_INNER_MIDPLANE, "Inner Midplane", &
            &   "All nodes (0D objects) along the inner midplane." )

        !! Get explicit object list of the grid subset using
        !! subroutine collectRadialVertexIndexListSubroutine
        !! (function collectRadialVertexIndexList transferred to subroutine,
        !! as array of certain dimension is required as an output)
        call collectRadialVertexIndexListSubroutine(gmap, cflag, xIn, yIn,  &
            &   topix, topiy, indexList2d)

        !! Initialize explicit object list for grid subset
        call createExplicitObjectListSingleSpace( grid_ggd,             &
            &   grid_ggd%grid_subset( GSubsetCount ), IDS_CLASS_NODE,   &
            &   indexList2d(:,SPACE_POLOIDALPLANE), IDS_CLASS_NODE,     &
            &   SPACE_POLOIDALPLANE )

        GSubsetCount = GSubsetCount + 1

        !! Create grid subset with one object list
        call createEmptyGridSubset(                           &
            &   grid_ggd%grid_subset( GSubsetCount ),         &
            &   GRID_SUBSET_OUTER_MIDPLANE, "Outer Midplane", &
            &   "All nodes (0D objects) along the outer midplane." )

        !! Get explicit object list of the grid subset using
        !! subroutine collectRadialVertexIndexListSubroutine
        !! (function collectRadialVertexIndexList transferred to subroutine,
        !! as array of certain dimension is required as an output)
        call collectRadialVertexIndexListSubroutine(gmap, cflag, xOut, yOut,  &
            &   topix, topiy, indexList2d)

        !! Initialize explicit object list for grid subset
        call createExplicitObjectListSingleSpace( grid_ggd,             &
            &   grid_ggd%grid_subset( GSubsetCount ), IDS_CLASS_NODE,   &
            &   indexList2d(:,SPACE_POLOIDALPLANE), IDS_CLASS_NODE,     &
            &   SPACE_POLOIDALPLANE )

        !! Adding other special cases
        deallocate(indexList2d)
        select case ( geoId )
        case ( GEOMETRY_LINEAR )
            iType = REGIONTYPE_YEDGE
            cls = CLASS_POLOIDALRADIAL_EDGE
            do j = 1, 2
                if ( j.eq.1 ) iSubset = GRID_SUBSET_SEPARATRIX
                if ( j.eq.2 ) iSubset = GRID_SUBSET_ACTIVE_SEPARATRIX
                if ( jsep .ne. B2_GRID_UNDEFINED ) then
                    GSubsetCount = GSubsetCount + 1

                    call logmsg( LOGDEBUG,                                &
                        &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"// &
                        &   int2str(GSubsetCount)//": "//                 &
                        &   gridSubsetName ( iSubset ) )

                    !! Create grid subset with one object list
                    call createEmptyGridSubset(                     &
                        &   grid_ggd%grid_subset( GSubsetCount ),   &
                        &   iSubset, gridSubsetName ( iSubset ),    &
                        &   gridSubsetDescription( iSubset ) )

                    nInd = gmap%b2nx
                    allocate( indexList2d(nInd, SPACE_COUNT) )
                    iInd = 0
                    do ix = 0, gmap%b2nx-1
                        iInd = iInd + 1
                        ind = gmap%mapFcI(ix, jsep, TOP)
                        indexList2d( iInd, SPACE_POLOIDALPLANE ) = ind
                    end do

                    !! Initialize explicit object list for grid subset
                    call createExplicitObjectListSingleSpace( grid_ggd,     &
                        &   grid_ggd%grid_subset( GSubsetCount ), sum(cls), &
                        &   indexList2d(:,SPACE_POLOIDALPLANE), sum(cls),   &
                        &   SPACE_POLOIDALPLANE )
                    deallocate(IndexList2d)

                end if
            end do
        case ( GEOMETRY_SN )
            iType = REGIONTYPE_YEDGE
            cls = CLASS_POLOIDALRADIAL_EDGE
            iSubset = GRID_SUBSET_SEPARATRIX
            GSubsetCount = GSubsetCount + 1

            call logmsg( LOGDEBUG,                                      &
                &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//       &
                &   int2str(GSubsetCount)//": "//                       &
                &   gridSubsetName ( iSubset ) )

            !! Create grid subset with one object list
            call createEmptyGridSubset(                     &
                &   grid_ggd%grid_subset( GSubsetCount ),   &
                &   iSubset, gridSubsetName ( iSubset ),    &
                &   gridSubsetDescription ( iSubset ) )

            nInd = gmap%b2nx
            allocate( indexList2d(nInd, SPACE_COUNT) )
            iInd = 0
            do ix = 0, gmap%b2nx-1
                iInd = iInd + 1
                ind = gmap%mapFcI(ix, jsep, TOP)
                indexList2d( iInd, SPACE_POLOIDALPLANE ) = ind
            end do

            !! Initialize explicit object list for grid subset
            call createExplicitObjectListSingleSpace( grid_ggd,     &
                &   grid_ggd%grid_subset( GSubsetCount ), sum(cls), &
                &   indexList2d(:,SPACE_POLOIDALPLANE), sum(cls),   &
                &   SPACE_POLOIDALPLANE )
            deallocate(IndexList2d)

        case ( GEOMETRY_CDN )
            iType = REGIONTYPE_YEDGE
            cls = CLASS_POLOIDALRADIAL_EDGE
            iSubset = GRID_SUBSET_SEPARATRIX
            GSubsetCount = GSubsetCount + 1

            call logmsg( LOGDEBUG,                                      &
                &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//       &
                &   int2str(GSubsetCount)//": "//                       &
                &   gridSubsetName ( iSubset ) )

            !! Create grid subset with one object list
            call createEmptyGridSubset(                     &
                &   grid_ggd%grid_subset( GSubsetCount ),   &
                &   iSubset, gridSubsetName ( iSubset ),    &
                &   gridSubsetDescription ( iSubset ) )

            nInd = 0
            do ix = 0, gmap%b2nx-1
                if (isRealCell(cflags(ix,jsep,CELLFLAG_TYPE))) nInd = nInd + 1
            end do
            allocate( indexList2d(nInd, SPACE_COUNT) )
            iInd = 0
            do ix = 0, gmap%b2nx-1
                if (.not.isRealCell(cflags(ix,jsep,CELLFLAG_TYPE))) cycle
                iInd = iInd + 1
                ind = gmap%mapFcI(ix, jsep, TOP)
                indexList2d( iInd, SPACE_POLOIDALPLANE ) = ind
            end do

            !! Initialize explicit object list for grid subset
            call createExplicitObjectListSingleSpace( grid_ggd,     &
                &   grid_ggd%grid_subset( GSubsetCount ), sum(cls), &
                &   indexList2d(:,SPACE_POLOIDALPLANE), sum(cls),   &
                &   SPACE_POLOIDALPLANE )
            deallocate(IndexList2d)

        case ( GEOMETRY_DDN_BOTTOM )
            iType = REGIONTYPE_YEDGE
            cls = CLASS_POLOIDALRADIAL_EDGE
            iSubset = GRID_SUBSET_SEPARATRIX
            GSubsetCount = GSubsetCount + 1

            call logmsg( LOGDEBUG,                                      &
                &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//       &
                &   int2str(GSubsetCount)//": "//                       &
                &   gridSubsetName ( iSubset ) )

            !! Create grid subset with one object list
            call createEmptyGridSubset(                     &
                &   grid_ggd%grid_subset( GSubsetCount ),   &
                &   iSubset, gridSubsetName ( iSubset ),    &
                &   gridSubsetDescription ( iSubset ) )

            nInd = 0
            do ix = 0, gmap%b2nx-1
                if (isRealCell(cflags(ix,jsep,CELLFLAG_TYPE))) then
                    ireg = region(ix,jsep,0)
                    if (ireg.eq.1 .or. ireg.eq.3 .or. &
                      & ireg.eq.5 .or. ireg.eq.8) then
                        nInd = nInd + 1
                    end if
                end if
            end do
            allocate( indexList2d(nInd, SPACE_COUNT) )
            iInd = 0
            do ix = 0, gmap%b2nx-1
                if (.not.isRealCell(cflags(ix,jsep,CELLFLAG_TYPE))) cycle
                ireg = region(ix,jsep,0)
                if (ireg.eq.1 .or. ireg.eq.3 .or. &
                  & ireg.eq.5 .or. ireg.eq.8) then
                    iInd = iInd + 1
                    ind = gmap%mapFcI(ix, jsep, TOP)
                    indexList2d( iInd, SPACE_POLOIDALPLANE ) = ind
                end if
            end do

            !! Initialize explicit object list for grid subset
            call createExplicitObjectListSingleSpace( grid_ggd,     &
                &   grid_ggd%grid_subset( GSubsetCount ), sum(cls), &
                &   indexList2d(:,SPACE_POLOIDALPLANE), sum(cls),   &
                &   SPACE_POLOIDALPLANE )
            deallocate(IndexList2d)

            iType = REGIONTYPE_CELL
            cls = CLASS_CELL
            iSubset = GRID_SUBSET_BETWEEN_SEPARATRICES
            GSubsetCount = GSubsetCount + 1

            call logmsg( LOGDEBUG,                                      &
                &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//       &
                &   int2str(GSubsetCount)//": "//                       &
                &   gridSubsetName ( iSubset ) )

            nInd = (topcut(2)-topcut(1))*((leftcut(2)+1)+(gmap%b2nx-rightcut(2)))
            allocate( indexList2d(nInd, SPACE_COUNT) )
            iInd = 0
            do iy = topcut(1)+1, topcut(2)
                do ix = 0, leftcut(2)
                    iInd = iInd + 1
                    ind = gmap%mapCvI(ix, iy)
                    indexList2d( iInd, SPACE_POLOIDALPLANE ) = ind
                end do
                do ix = rightcut(2)+1, gmap%b2nx-1
                    iInd = iInd + 1
                    ind = gmap%mapCvI(ix, iy)
                    indexList2d( iInd, SPACE_POLOIDALPLANE ) = ind
                end do
            end do

            !! Create grid subset with one object list
            call createEmptyGridSubset(                     &
                &   grid_ggd%grid_subset( GSubsetCount ),   &
                &   iSubset, gridSubsetName ( iSubset ),    &
                &   gridSubsetDescription ( iSubset ) )

            !! Initialize explicit object list for grid subset
            call createExplicitObjectListSingleSpace( grid_ggd,     &
                &   grid_ggd%grid_subset( GSubsetCount ), sum(cls), &
                &   indexList2d(:,SPACE_POLOIDALPLANE), sum(cls),   &
                &   SPACE_POLOIDALPLANE )
            deallocate(IndexList2d)

        case ( GEOMETRY_DDN_TOP )
            iType = REGIONTYPE_YEDGE
            cls = CLASS_POLOIDALRADIAL_EDGE
            iSubset = GRID_SUBSET_SEPARATRIX
            GSubsetCount = GSubsetCount + 1

            call logmsg( LOGDEBUG,                                      &
                &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//       &
                &   int2str(GSubsetCount)//": "//                       &
                &   gridSubsetName ( iSubset ) )

            !! Create grid subset with one object list
            call createEmptyGridSubset(                     &
                &   grid_ggd%grid_subset( GSubsetCount ),   &
                &   iSubset, gridSubsetName ( iSubset ),    &
                &   gridSubsetDescription ( iSubset ) )

            nInd = 0
            do ix = 0, gmap%b2nx-1
                if (isRealCell(cflags(ix,jsep,CELLFLAG_TYPE))) then
                    ireg = region(ix,jsep,0)
                    if (ireg.eq.1 .or. ireg.eq.4 .or. &
                      & ireg.eq.5 .or. ireg.eq.7) then
                        nInd = nInd + 1
                    end if
                end if
            end do
            allocate( indexList2d(nInd, SPACE_COUNT) )
            iInd = 0
            do ix = 0, gmap%b2nx-1
                if (.not.isRealCell(cflags(ix,jsep,CELLFLAG_TYPE))) cycle
                ireg = region(ix,jsep,0)
                if (ireg.eq.1 .or. ireg.eq.4 .or. &
                  & ireg.eq.5 .or. ireg.eq.7) then
                    iInd = iInd + 1
                    ind = gmap%mapFcI(ix, jsep, TOP)
                    indexList2d( iInd, SPACE_POLOIDALPLANE ) = ind
                end if
            end do

            !! Initialize explicit object list for grid subset
            call createExplicitObjectListSingleSpace( grid_ggd,     &
                &   grid_ggd%grid_subset( GSubsetCount ), sum(cls), &
                &   indexList2d(:,SPACE_POLOIDALPLANE), sum(cls),   &
                &   SPACE_POLOIDALPLANE )
            deallocate(IndexList2d)

            iType = REGIONTYPE_CELL
            cls = CLASS_CELL
            iSubset = GRID_SUBSET_BETWEEN_SEPARATRICES
            GSubsetCount = GSubsetCount + 1

            call logmsg( LOGDEBUG,                                      &
                &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//       &
                &   int2str(GSubsetCount)//": "//                       &
                &   gridSubsetName ( iSubset ) )

            nInd = (topcut(1)-topcut(2))*((nxtl-leftcut(1)-1)+(rightcut(1)-nxtr))
            allocate( indexList2d(nInd, SPACE_COUNT) )
            iInd = 0
            do iy = topcut(2)+1, topcut(1)
                do ix = leftcut(1)+1, nxtl-1
                    iInd = iInd + 1
                    ind = gmap%mapCvI(ix, iy)
                    indexList2d( iInd, SPACE_POLOIDALPLANE ) = ind
                end do
                do ix = nxtr+1, rightcut(1)
                    iInd = iInd + 1
                    ind = gmap%mapCvI(ix, iy)
                    indexList2d( iInd, SPACE_POLOIDALPLANE ) = ind
                end do
            end do

            !! Create grid subset with one object list
            call createEmptyGridSubset(                     &
                &   grid_ggd%grid_subset( GSubsetCount ),   &
                &   iSubset, gridSubsetName ( iSubset ),    &
                &   gridSubsetDescription ( iSubset ) )

            !! Initialize explicit object list for grid subset
            call createExplicitObjectListSingleSpace( grid_ggd,     &
                &   grid_ggd%grid_subset( GSubsetCount ), sum(cls), &
                &   indexList2d(:,SPACE_POLOIDALPLANE), sum(cls),   &
                &   SPACE_POLOIDALPLANE )
            deallocate(IndexList2d)

        end select

        !! Do grid subsets for separatrix points of interest
        !! Outer midplane separatrix
        if (jsep .ne. B2_GRID_UNDEFINED) then
            GSubsetCount = GSubsetCount + 1

            call logmsg( LOGDEBUG,                                      &
                &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//       &
                &   int2str(GSubsetCount)//": "//                       &
                &   gridSubsetName ( GRID_SUBSET_OUTER_MIDPLANE_SEPARATRIX ) )

            !! Create grid subset with one object list
            call createEmptyGridSubset(                                 &
                &   grid_ggd%grid_subset( GSubsetCount ),               &
                &   GRID_SUBSET_OUTER_MIDPLANE_SEPARATRIX,              &
                &   gridSubsetName ( GRID_SUBSET_OUTER_MIDPLANE_SEPARATRIX ), &
                &   gridSubsetDescription ( GRID_SUBSET_OUTER_MIDPLANE_SEPARATRIX ) )

            nInd = 1
            allocate( indexList2d(nInd, SPACE_COUNT) )
            indexList2d( 1, SPACE_POLOIDALPLANE ) =  &
                &   gmap%mapVxI( xOut, jsep, VX_UPPERLEFT )

            !! Initialize explicit object list for grid subset
            call createExplicitObjectListSingleSpace( grid_ggd,           &
                &   grid_ggd%grid_subset( GSubsetCount ), IDS_CLASS_NODE, &
                &   indexList2d(:,SPACE_POLOIDALPLANE), IDS_CLASS_NODE,   &
                &   SPACE_POLOIDALPLANE )
            deallocate(IndexList2d)

        !! Inner midplane separatrix
            GSubsetCount = GSubsetCount + 1

            call logmsg( LOGDEBUG,                                        &
                &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//         &
                &   int2str(GSubsetCount)//": "//                         &
                &   gridSubsetName ( GRID_SUBSET_INNER_MIDPLANE_SEPARATRIX ) )

            !! Create grid subset with one object list
            call createEmptyGridSubset(                                   &
                &   grid_ggd%grid_subset( GSubsetCount ),                 &
                &   GRID_SUBSET_INNER_MIDPLANE_SEPARATRIX,                &
                &   gridSubsetName ( GRID_SUBSET_INNER_MIDPLANE_SEPARATRIX ), &
                &   gridSubsetDescription ( GRID_SUBSET_INNER_MIDPLANE_SEPARATRIX ) )

            nInd = 1
            allocate( indexList2d(nInd, SPACE_COUNT) )
            indexList2d( 1, SPACE_POLOIDALPLANE ) =    &
                &   gmap%mapVxI( xIn, jsep, VX_UPPERLEFT )

            !! Initialize explicit object list for grid subset
            call createExplicitObjectListSingleSpace( grid_ggd,           &
                &   grid_ggd%grid_subset( GSubsetCount ), IDS_CLASS_NODE, &
                &   indexList2d(:,SPACE_POLOIDALPLANE), IDS_CLASS_NODE,   &
                &   SPACE_POLOIDALPLANE )
            deallocate(IndexList2d)

        end if

        !! Outer strikepoint
        select case ( geoId )
           case ( GEOMETRY_SN, GEOMETRY_LIMITER )
              if (LSN) then
                  ix = gmap%b2nx - 1
                  iVx = VX_LOWERRIGHT
              else
                  ix = 0
                  iVx = VX_LOWERLEFT
              end if
           case ( GEOMETRY_LINEAR, GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM, &
               &  GEOMETRY_STELLARATORISLAND )
              ix = gmap%b2nx - 1
              iVx = VX_LOWERRIGHT
           case ( GEOMETRY_DDN_TOP )
              ix = nxtr + 1
              iVx = VX_LOWERLEFT
           case default
              ix = B2_GRID_UNDEFINED
              iVx = GRID_UNDEFINED
        end select
        if (ix .ne. B2_GRID_UNDEFINED) then
            GSubsetCount = GSubsetCount + 1

            call logmsg( LOGDEBUG,                                      &
                &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//       &
                &   int2str(GSubsetCount)//": "//                       &
                &   gridSubsetName ( GRID_SUBSET_OUTER_STRIKEPOINT ) )

            !! Create grid subset with one object list
            call createEmptyGridSubset(                                 &
                &   grid_ggd%grid_subset( GSubsetCount ),               &
                &   GRID_SUBSET_OUTER_STRIKEPOINT,                      &
                &   gridSubsetName ( GRID_SUBSET_OUTER_STRIKEPOINT ),   &
                &   gridSubsetDescription ( GRID_SUBSET_OUTER_STRIKEPOINT ) )

            nInd = 1
            allocate( indexList2d(nInd, SPACE_COUNT) )
            indexList2d( 1, SPACE_POLOIDALPLANE ) = gmap%mapVxI( ix, jsep, iVx )

            !! Initialize explicit object list for grid subset
            call createExplicitObjectListSingleSpace( grid_ggd,           &
                &   grid_ggd%grid_subset( GSubsetCount ), IDS_CLASS_NODE, &
                &   indexList2d(:,SPACE_POLOIDALPLANE), IDS_CLASS_NODE,   &
                &   SPACE_POLOIDALPLANE )
            deallocate(IndexList2d)

        end if

        !! Inner strikepoint
        select case ( geoId )
           case ( GEOMETRY_SN, GEOMETRY_LIMITER )
              if (LSN) then
                  ix = 0
                  iVx = VX_LOWERLEFT
              else
                  ix = gmap%b2nx - 1
                  iVx = VX_LOWERRIGHT
              end if
           case ( GEOMETRY_LINEAR, GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM, &
               &  GEOMETRY_STELLARATORISLAND )
              ix = 0
              iVx = VX_LOWERLEFT
           case ( GEOMETRY_DDN_TOP )
              ix = nxtl - 1
              iVx = VX_LOWERRIGHT
           case default
              ix = B2_GRID_UNDEFINED
              iVx = GRID_UNDEFINED
        end select
        if (ix .ne. B2_GRID_UNDEFINED) then
            GSubsetCount = GSubsetCount + 1

            call logmsg( LOGDEBUG,                                      &
                &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//       &
                &   int2str(GSubsetCount)//": "//                       &
                &   gridSubsetName ( GRID_SUBSET_INNER_STRIKEPOINT ) )

            !! Create grid subset with one object list
            call createEmptyGridSubset(                                 &
                &   grid_ggd%grid_subset( GSubsetCount ),               &
                &   GRID_SUBSET_INNER_STRIKEPOINT,                      &
                &   gridSubsetName ( GRID_SUBSET_INNER_STRIKEPOINT ),   &
                &   gridSubsetDescription ( GRID_SUBSET_INNER_STRIKEPOINT ) )

            nInd = 1
            allocate( indexList2d(nInd, SPACE_COUNT) )
            indexList2d( 1, SPACE_POLOIDALPLANE ) = gmap%mapVxI( ix, jsep, iVx )

            !! Initialize explicit object list for grid subset
            call createExplicitObjectListSingleSpace( grid_ggd,           &
                &   grid_ggd%grid_subset( GSubsetCount ), IDS_CLASS_NODE, &
                &   indexList2d(:,SPACE_POLOIDALPLANE), IDS_CLASS_NODE,   &
                &   SPACE_POLOIDALPLANE )
            deallocate(IndexList2d)

        end if

        !! Outer strikepoint inactive
        select case ( geoId )
           case ( GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM )
              ix = nxtr + 1
              iVx = VX_LOWERLEFT
           case ( GEOMETRY_DDN_TOP )
              ix = gmap%b2nx - 1
              iVx = VX_LOWERRIGHT
           case default
              ix = B2_GRID_UNDEFINED
              iVx = GRID_UNDEFINED
        end select
        if (ix .ne. B2_GRID_UNDEFINED) then
            GSubsetCount = GSubsetCount + 1

            call logmsg( LOGDEBUG,                                      &
                &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//       &
                &   int2str(GSubsetCount)//": "//                       &
                &   gridSubsetName ( GRID_SUBSET_OUTER_STRIKEPOINT_INACTIVE ) )

            !! Create grid subset with one object list
            call createEmptyGridSubset(                                 &
                &   grid_ggd%grid_subset( GSubsetCount ),               &
                &   GRID_SUBSET_OUTER_STRIKEPOINT_INACTIVE,             &
                &   gridSubsetName ( GRID_SUBSET_OUTER_STRIKEPOINT_INACTIVE ), &
                &   gridSubsetDescription ( GRID_SUBSET_OUTER_STRIKEPOINT_INACTIVE ) )

            nInd = 1
            allocate( indexList2d(nInd, SPACE_COUNT) )
            indexList2d( 1, SPACE_POLOIDALPLANE ) = gmap%mapVxI( ix, jsep, iVx )

            !! Initialize explicit object list for grid subset
            call createExplicitObjectListSingleSpace( grid_ggd,           &
                &   grid_ggd%grid_subset( GSubsetCount ), IDS_CLASS_NODE, &
                &   indexList2d(:,SPACE_POLOIDALPLANE), IDS_CLASS_NODE,   &
                &   SPACE_POLOIDALPLANE )
            deallocate(IndexList2d)

        end if

        !! Inner strikepoint inactive
        select case ( geoId )
           case ( GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM )
              ix = nxtl - 1
              iVx = VX_LOWERRIGHT
           case ( GEOMETRY_DDN_TOP )
              ix = 0
              iVx = VX_LOWERLEFT
           case default
              ix = B2_GRID_UNDEFINED
              iVx = GRID_UNDEFINED
        end select
        if (ix .ne. B2_GRID_UNDEFINED) then
            GSubsetCount = GSubsetCount + 1

            call logmsg( LOGDEBUG,                                      &
                &   "b2_IMAS_Fill_Grid_Desc: add grid subset #"//       &
                &   int2str(GSubsetCount)//": "//                       &
                &   gridSubsetName ( GRID_SUBSET_INNER_STRIKEPOINT_INACTIVE ) )

            !! Create grid subset with one object list
            call createEmptyGridSubset(                                 &
                &   grid_ggd%grid_subset( GSubsetCount ),               &
                &   GRID_SUBSET_INNER_STRIKEPOINT_INACTIVE,             &
                &   gridSubsetName ( GRID_SUBSET_INNER_STRIKEPOINT_INACTIVE ), &
                &   gridSubsetDescription ( GRID_SUBSET_INNER_STRIKEPOINT_INACTIVE ) )

            nInd = 1
            allocate( indexList2d(nInd, SPACE_COUNT) )
            indexList2d( 1, SPACE_POLOIDALPLANE ) = gmap%mapVxI( ix, jsep, iVx )

            !! Initialize explicit object list for grid subset
            call createExplicitObjectListSingleSpace( grid_ggd,           &
                &   grid_ggd%grid_subset( GSubsetCount ), IDS_CLASS_NODE, &
                &   indexList2d(:,SPACE_POLOIDALPLANE), IDS_CLASS_NODE,   &
                &   SPACE_POLOIDALPLANE )
            deallocate(IndexList2d)

        end if

        deallocate(indexList1d)

        call logmsg( LOGDEBUG, "b2_IMAS_Fill_Grid_Desc: wrote total of " &
            &   //int2str(GSubsetCount)//" grid subsets (expected was "  &
            &   //int2str(size(grid_ggd%grid_subset))//")" )

        call xertst( GSubsetCount == size(grid_ggd%grid_subset), &
            &  "Assert error (grid subset count) in fill_In_GridSubset_Desc" )

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

        !! Loop over all edges in core boundary grid subset
        do iObj = 1, getGridSubsetSize(GridSubset)
            obj = getGridSubsetObject(GridSubset, iObj)
            !! Expect an edge
            call xertst( obj%cls( SPACE_POLOIDALPLANE ) == &
                &   IDS_CLASS_POLOIDALRADIAL_EDGE, &
                & "b2mod_ual_io_grid find_Midplane_Cells: assertion failure." )
            !! ...which is aligned along the x-direction
            call xertst( gmap%mapFcIFace( obj%ind( SPACE_POLOIDALPLANE ) ) ==  &
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

        if ( xIn == huge( xIn ) ) then
          call logmsg( LOGWARNING, "find_Midplane_Cells: "   &
              &   //"did not find inner midplane position. " &
              &   //"Assigning it to jxi index." )
          xIn = jxi
        end if
        if ( xOut == huge( xOut ) ) then
          call logmsg( LOGWARNING, "find_Midplane_Cells: "   &
              &   //"did not find outer midplane position. " &
              &   //"Assigning it to jxa index." )
          xOut = jxa
        end if

    end subroutine find_Midplane_Cells

#endif
#else
#ifdef ITM_ENVIRONMENT_LOADED

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
    real(ITM_R8), intent(in), optional :: vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny)
    real(ITM_R8), save :: width = 1.0_ITM_R8

    ! internal
    integer, parameter :: NDIM = 2

    call ipgetr('b2agmt_1d_width', width)
    call xertst(0.0_ITM_R8.lt.width, 'faulty input width')
    call xertst( present( gs ) .EQV. present( qc ) , &
        "Assert error ( gz or qc missing ) in b2ITMFillGridDescription" )

    call fill_In_Grid_Desc()
    call fillInSubGridDescription()

  contains

    !! Part 1: fill in grid description
    subroutine fill_In_Grid_Desc()

      !! internal
      integer :: iVx, iFc, iCv, ix, iy, nix, niy, i, dir

      allocate( itmgrid % spaces(SPACE_COUNT) )

      !! Coordinate types
      ! (dimension of space = NDIM = size( coordtype )
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % coordtype(NDIM, 1) )
      itmgrid % spaces(SPACE_POLOIDALPLANE) % coordtype(:, 1) &
           & = (/ COORDTYPE_R, COORDTYPE_Z /)

      !! Have two types of objects: 1d edges, 2d cells
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(NDIM + 1) )

      !! Fill in node information
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(1) % geo(gmap%nVx, NDIM, 1, 1) )
      if (isymm.eq.0 .or. isymm.eq.1 .or. isymm.eq.2) then
        do iVx = 1, gmap % nVx
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(1) % geo(iVx, 1, 1, 1) = &
              & crx( gmap % mapVxix( iVx ), gmap % mapVxiy( iVx ), gmap % mapVxIVx( iVx ) )
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(1) % geo(iVx, 2, 1, 1) = &
              & cry( gmap % mapVxix( iVx ), gmap % mapVxiy( iVx ), gmap % mapVxIVx( iVx ) )
        end do
      else if (isymm.eq.3 .or. isymm.eq.4) then
        do iVx = 1, gmap % nVx
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(1) % geo(iVx, 1, 1, 1) = &
              & cry( gmap % mapVxix( iVx ), gmap % mapVxiy( iVx ), gmap % mapVxIVx( iVx ) )
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(1) % geo(iVx, 2, 1, 1) = &
              & crx( gmap % mapVxix( iVx ), gmap % mapVxiy( iVx ), gmap % mapVxIVx( iVx ) )
        end do
      end if

      !! Fill in object definitions (i.e. what objects compose an object)

      !! 1d objects: edges
      !! ...have two boundaries
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( gmap%nFcx + gmap%nFcy, 2) )
      !! ...have two neighbours, in positive and negative coordinate direction, one on each side
      !! (for x-aligned edges: along flux surface, for y-aligned edges: orthogonal to flux surface)
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour(gmap%nFcx + gmap%nFcy, 2, 1) )
      !! 1d object measure: face area
      if (present(gs)) allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % measure( gmap % nFcx + gmap % nFcy, 1 ) )
      !! first set all boundary & connectivity information to undefined
      itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary = GRID_UNDEFINED
      itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour = GRID_UNDEFINED
      !! x-aligned edges
      do iFc = 1, gmap % nFcx
          !! get position of this face in the B2 grid
          ix = gmap % mapFcix( iFc )
          iy = gmap % mapFciy( iFc )
          !! get index of start vertex
          !! objdef dims: index of edge, 1=start node, 1=one-dimensional object
          select case ( gmap % mapFcIFace( iFc ) )
          case( BOTTOM )
             !! start index: 1=start node
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( iFc, 1 ) = gmap % mapVxI( ix, iy, VX_LOWERLEFT )
             if (gmap % mapVxI( ix, iy, VX_LOWERLEFT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2ITMFillGD: BOTTOM edge at ("//int2str(ix)//","//int2str(iy)//") has no start node")
             end if
             !! end vertex: 2=end node
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( iFc, 2 ) = gmap % mapVxI( ix, iy, VX_LOWERRIGHT )
             if (gmap % mapVxI( ix, iy, VX_LOWERRIGHT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2ITMFillGD: BOTTOM edge at ("//int2str(ix)//","//int2str(iy)//") has no end node")
             end if
          case( TOP )
             !! start index: 1=start node
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( iFc, 1 ) = gmap % mapVxI( ix, iy, VX_UPPERLEFT )
             if (gmap % mapVxI( ix, iy, VX_UPPERLEFT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2ITMFillGD: TOP edge at ("//int2str(ix)//","//int2str(iy)//") has no start node")
             end if
             !! end vertex: 2=end node
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( iFc, 2 ) = gmap % mapVxI( ix, iy, VX_UPPERRIGHT )
             if (gmap % mapVxI( ix, iy, VX_UPPERRIGHT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2ITMFillGD: TOP edge at ("//int2str(ix)//","//int2str(iy)//") has no end node")
             end if
          end select

          !! Neighbour edges of this edge
          !! Left neighbour: edge continuing to the left of this edge
          nix = leftix( ix, iy )
          niy = leftiy( ix, iy )
          if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour( iFc, 1, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( iFc ) )
          end if
          !! Right neighbour: edge continuing to the right of this edge
          nix = rightix( ix, iy )
          niy = rightiy( ix, iy )
          if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour( iFc, 2, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( iFc ) )
          end if

          !! measure: area
          if (present(gs)) itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % measure( iFc, 1 ) = gs(ix, iy, ALIGNX)
      end do

      !! y-aligned edges
      do iFc = gmap % nFcx + 1, gmap % nFcx + gmap % nFcy
          !! get position of this face in the B2 grid
          ix = gmap % mapFcix( iFc )
          iy = gmap % mapFciy( iFc )
!!$          if (gmap%mapCvI(ix, iy) == GRID_UNDEFINED) then
!!$                  call logmsg(LOGWARNING, "b2ITMFillGD: writing out edges for unused cell ("//int2str(ix)//","//int2str(iy)//")")
!!$          end if

          select case ( gmap % mapFcIFace( iFc ) )
          case( LEFT )
              !! start index: 1=start node
              itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( iFc, 1 ) = gmap % mapVxI( ix, iy, VX_LOWERLEFT )
              if (gmap % mapVxI( ix, iy, VX_LOWERLEFT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2ITMFillGD: LEFT edge at ("//int2str(ix)//","//int2str(iy)//") has no start node")
              end if
          !! end vertex: 2=end node
              itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( iFc, 2 ) = gmap % mapVxI( ix, iy, VX_UPPERLEFT )
              if (gmap % mapVxI( ix, iy, VX_UPPERLEFT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2ITMFillGD: LEFT edge at ("//int2str(ix)//","//int2str(iy)//") has no end node")
          end if
              !if (itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( iFc, 1 )


          case( RIGHT )
              !! start index: 1=start node
              itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( iFc, 1 ) = gmap % mapVxI( ix, iy, VX_LOWERRIGHT )
              if (gmap % mapVxI( ix, iy, VX_LOWERRIGHT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2ITMFillGD: RIGHT edge at ("//int2str(ix)//","//int2str(iy)//") has no start node")
          end if
              !! end vertex: 2=end node
              itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( iFc, 2 ) = gmap % mapVxI( ix, iy, VX_UPPERRIGHT )
              if (gmap % mapVxI( ix, iy, VX_UPPERRIGHT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2ITMFillGD: RIGHT edge at ("//int2str(ix)//","//int2str(iy)//") has no end node")
              end if
          end select


          !! Neighbour edges of this edge
          !! Bottom neighbour: edge continuing to the bottom of this edge
          nix = bottomix( ix, iy )
          niy = bottomiy( ix, iy )
          if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour( iFc, 1, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( iFc ) )
          end if
          !! Top neighbour: edge continuing to the top of this edge
          nix = topix( ix, iy )
          niy = topiy( ix, iy )
          if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour( iFc, 2, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( iFc ) )
          end if
          !! measure: area
          if (present(gs)) itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % measure( iFc, 1 ) = gs(ix, iy, ALIGNY)*qc(ix, iy)
      end do

      !! 2d objects: cells
      !! ...have four boundaries
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( gmap%nCv, 4) )
      !! 2d object measure: cell volume
      if (present(vol)) allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % measure( gmap % nCv, 1 ) )
      !! Also store additional geometry information: position in computational space
      !! FIXME: this should go into alternate geometry, which is not available yet for grid objects
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % geo(gmap%nCv, 2, 1, 1) )

      !! first set all boundary information to undefined
      itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary = GRID_UNDEFINED

      do iCv = 1, gmap % nCv
          ix = gmap % mapCvix( iCv )
          iy = gmap % mapCviy( iCv )

          !! Set position in computational space
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % geo(iCv, 1, 1, 1) = ix
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % geo(iCv, 2, 1, 1) = iy
          !! put edges composing the quadrilateral in the list: left edge (y-aligned)
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( iCv, 1 ) = gmap % mapFcI( ix, iy, LEFT )
          !! bottom edge (x-aligned)
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( iCv, 2 ) = gmap % mapFcI( ix, iy, BOTTOM )
          !! right edge (y-aligned)
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( iCv, 3 ) = gmap % mapFcI( ix, iy, RIGHT )
          !! top edge (x-aligned)
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( iCv, 4 ) = gmap % mapFcI( ix, iy, TOP )
      end do

      !! Fill in connectivity information
      !! ...have one neighbour per boundary
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % neighbour( gmap%nCv, 4, 1) )
      !! first set all to undefined
      itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % neighbour = GRID_UNDEFINED

      do iCv = 1, gmap % nCv
          ix = gmap % mapCvix( iCv )
          iy = gmap % mapCviy( iCv )

          do dir = LEFT, TOP
             call get_Neighbour(nx, ny, leftix, leftiy, rightix, rightiy, topix, topiy, bottomix, bottomiy, &
                  & ix, iy, dir, nix, niy)
             if ( .not. is_Unneeded_Cell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % neighbour(iCv, dir+1, 1) = gmap % mapCvI( nix, niy )
             end if
          end do

      end do

      !! Fill in x-point indices
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % xpoints( gmap % nsv ) )
      itmgrid % spaces(SPACE_POLOIDALPLANE) % xpoints = gmap % svi(1:gmap % nsv)

      !! If requested, add a second space for the toroidal angle
      if (SPACE_COUNT == SPACE_TOROIDALANGLE) then
        if (isymm.eq.0) then
          call gridSetupStruct1dSpace( itmgrid%spaces(SPACE_TOROIDALANGLE), &
              & COORDTYPE_Y, &
              & (/  ( ( width / NNODES_TOROIDAL ) * i, i = 0, NNODES_TOROIDAL ) /) , &
              & .false. ) !! periodic = .false.
        else
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
      end if

    end subroutine fill_In_Grid_Desc


    !! Part 2: define subgrids
    subroutine fillInSubGridDescription

      !! internal
      integer :: geoId, iRegion, subgridCount, iType, nSubgrid
      integer :: xIn, yIn, xOut, yOut, iCoreSg
      integer :: cls(SPACE_COUNT_MAX)
      integer, allocatable :: xpoints(:,:)

      geoId = geometryId(nnreg, isymm, periodic_bc, topcut)

      !! Figure out total number of subgrids
      !! Do generic subgrids + subgrids
      nSubgrid = B2_GENERIC_SUBGRID_COUNT + regionCountTotal(geoId)
      !! Inner/outer midplane subgrids
      nSubgrid = nSubgrid + 2

      call logmsg( LOGDEBUG, "b2ITMFillGridDescription: expecting total of "&
          &//int2str(nSubgrid)//" subgrids" )
      allocate( itmgrid % subgrids( nSubgrid ) )

      !! Set up generic subgrids

      !! B2_SUBGRID_CELLS: all 2d cells, one implicit object list
      call createSubGridForClass( itmgrid, itmgrid % subgrids( B2_SUBGRID_CELLS ), &
          & CLASS_CELL(1:SPACE_COUNT), 'Cells' )

      !! B2_SUBGRID_NODES: all nodes, one implicit object list
      call createSubGridForClass( itmgrid, itmgrid % subgrids( B2_SUBGRID_NODES ), &
          & CLASS_NODE(1:SPACE_COUNT), 'Nodes' )

      !! B2_SUBGRID_EDGES: all edges, one implicit object list
      call createSubGridForClass( itmgrid, itmgrid % subgrids( B2_SUBGRID_EDGES ), &
          & CLASS_POLOIDALRADIAL_EDGE(1:SPACE_COUNT), 'Edges' )

      !! B2_SUBGRID_EDGES_X: x-aligned edges. One implicit object list, range over x edges
      !! Create subgrid with one object list
      call createSubGrid( itmgrid % subgrids( B2_SUBGRID_EDGES_X ), 1, 'x-aligned edges' )
      !! Initialize implicit object list for edges (class (/1/) )
      call createImplicitObjectList( itmgrid, itmgrid % subgrids( B2_SUBGRID_EDGES_X ) % list(1), &
          & CLASS_POLOIDALRADIAL_EDGE(1:SPACE_COUNT) )
      itmgrid % subgrids( B2_SUBGRID_EDGES_X ) % list(1) % indset(1) &
          & = createIndexListForRange( 1, gmap%nFcx )
      if ( SPACE_COUNT == SPACE_TOROIDALANGLE ) then
          itmgrid % subgrids( B2_SUBGRID_EDGES_X ) % list(1) % indset(2) &
              & = createIndexListForRange( 1, 1 )
      end if

      !! B2_SUBGRID_EDGES_Y: y-aligned edges.
      !! One implicit object list, range over y edges. Same procedure.
      call createSubGrid( itmgrid % subgrids( B2_SUBGRID_EDGES_Y ),  &
          &   1, 'y-aligned edges' )
      call createImplicitObjectList( itmgrid,                        &
          &   itmgrid % subgrids( B2_SUBGRID_EDGES_Y ) % list(1),    &
          &   CLASS_POLOIDALRADIAL_EDGE(1:SPACE_COUNT) )
      itmgrid % subgrids( B2_SUBGRID_EDGES_Y ) % list(1) % indset(1) &
          & = createIndexListForRange( gmap%nFcx + 1, gmap%nFcx + gmap%nFcy )
      if ( SPACE_COUNT == SPACE_TOROIDALANGLE ) then
          itmgrid % subgrids( B2_SUBGRID_EDGES_Y ) % list(1) % indset(2) &
              & = createIndexListForRange( 1, 1 )
      end if

      !! Subgrid of all x-points (in one poloidal plane at toroidal index 1)
      !! Assemble object descriptor for x-points
      allocate( xpoints(gmap%nsv, SPACE_COUNT) )
      xpoints = 1
      xpoints(:, SPACE_POLOIDALPLANE) = gmap%svi(1:gmap%nsv)
      call createSubGridForExplicitList( itmgrid,     &
          & itmgrid % subgrids( B2_SUBGRID_XPOINTS ), &
          & CLASS_NODE(1:SPACE_COUNT), xpoints, 'x-points' )

      !! Set up specific subgrids by collecting edges for regions

      !! Start counting from end of generic subgrids
      subgridCount = B2_GENERIC_SUBGRID_COUNT

      !! Cell + edge subgrids
      do iType = REGIONTYPE_CELL, REGIONTYPE_YEDGE

          select case(iType)
          case( REGIONTYPE_CELL )
              cls = CLASS_CELL
          case( REGIONTYPE_YEDGE, REGIONTYPE_XEDGE )
              cls = CLASS_POLOIDALRADIAL_EDGE
          end select

          do iRegion = 1, regionCount(geoId, iType)
              subgridCount = subgridCount + 1

              call logmsg( LOGDEBUG, "b2ITMFillGridDescription: add subgrid #"//&
                  & int2str(subgridCount)//" for iType "//int2str(iType)//&
                  &", iRegion "//int2str(iRegion)//": "//regionName(geoId, iType, iRegion) )

              call createSubGridForExplicitList( itmgrid, itmgrid % subgrids( subgridCount ), &
                  & cls(1:SPACE_COUNT), &
                  & collectIndexListForRegion(gmap, cflag, region, iType, iRegion), &
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

      call xertst( subgridCount == size(itmgrid%subgrids), &
       &  "Assert error (subgrid count) in fillInSubGridDescription" )

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

    !! Loop over all edges in core boundary subgrid
    do iObj = 1, gridSubGridSize(coreBndSubgrid)
        obj = subGridGetObject(coreBndSubgrid, iObj)
        !! Expect an edge
        call xertst( all(obj%cls(1:SPACE_COUNT) == CLASS_POLOIDALRADIAL_EDGE(1:SPACE_COUNT)), &
        &   "Assert error 1 (edge test) in find_Midplane_Cells" )

        !! ...which is aligned along the x-direction
        call xertst( gmap % mapFcIFace(obj%ind(SPACE_POLOIDALPLANE)) == BOTTOM, &
        &   "Assert error 2 (bottom edge) in find_Midplane_Cells" )
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
            !! Stepped outside grid?
            if (is_Unneeded_Cell( gmap%b2nx, gmap%b2ny, cflag, .true., nix, niy) ) then
                exit
            else
                nVx = nVx + 1
                ix = nix
                iy = niy
            end if
        end do

        allocate( indexList(nVx, SPACE_COUNT) )
        if (SPACE_COUNT == SPACE_TOROIDALANGLE) &
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

        call xertst( iVx == nVx , &
        &   "Assert error (vertex count) in collectRadialVertexIndexList" )

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
            !! Stepped outside grid ?
            if (is_Unneeded_Cell( gmap%b2nx, gmap%b2ny, cflag, .true., nix, niy) ) then
                exit
            else
                nVx = nVx + 1
                ix = nix
                iy = niy
            end if
        end do

        allocate( indexList(nVx, SPACE_COUNT) )
        if (SPACE_COUNT == SPACE_TOROIDALANGLE) &
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

        call xertst( iVx == nVx , &
        &   "Assert error (vertex count) in collectRadialVertexIndexListSubroutine" )

    end subroutine collectRadialVertexIndexListSubroutine


    !> Build an index list of all objects of a given region type
    !! (b2mod_connectivity.REGIONTYPE_*) for a given region id.
    !! @result The list of indices for all objects that constitute this grid
    !! region. The array has two dimensions because it is given as a list of
    !! object descriptors.
    function collectIndexListForRegion(gmap, cflag, region, iRegionType,   &
        &   iRegion) result( indexList )
        integer, allocatable, dimension(:,:) :: indexList

        type(B2GridMap), intent(in) :: gmap
        integer, intent(in) :: region(-1:gmap%b2nx,-1:gmap%b2ny,0:2)
        integer, intent(in) :: iRegionType, iRegion
        integer, intent(in) ::  &
            &   cflag( -1:gmap%b2nx, -1:gmap%b2ny, CARREOUT_NCELLFLAGS )

        !! internal
        integer :: ix, iy, nInd, iInd, ind

        !! Figure out how many indices to expect. A simple count of the form
        !! nInd = count( region(:,:,iRegionType) == iRegion )
        !! will not do, because we have to account for removed objects (ghost cells/edges).

        !! search the relevant objects and count them
        nInd = 0
        do ix = -1, gmap%b2nx
            do iy = -1, gmap%b2ny

                if ( region(ix, iy, iRegionType) == iRegion ) then
                  !! Get index depending on what object type we are looking at
                    select case (iRegionType)
                    case (REGIONTYPE_CELL)
                        if ( is_Unneeded_Cell( gmap%b2nx, gmap%b2ny, cflag, &
                           & INCLUDE_GHOST_CELLS, ix, iy ) ) cycle
                        ind = gmap%mapCvI(ix, iy)
                    case (REGIONTYPE_XEDGE)
                        ind = gmap%mapFcI(ix, iy, LEFT)
                    case (REGIONTYPE_YEDGE)
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
                        if ( is_Unneeded_Cell( gmap%b2nx, gmap%b2ny, cflag, &
                           & INCLUDE_GHOST_CELLS, ix, iy ) ) cycle
                        ind = gmap%mapCvI(ix, iy)
                    case (REGIONTYPE_XEDGE)
                        ind = gmap%mapFcI(ix, iy, LEFT)
                    case (REGIONTYPE_YEDGE)
                        ind = gmap%mapFcI(ix, iy, BOTTOM)
                    end select

                    if ( ind /= B2_GRID_UNDEFINED ) then
                        iInd = iInd + 1
                        call xertst(iInd <= nInd, &
                        &   "Assert error 1 (index) in collectIndexListForRegion" )
                        indexList( iInd, SPACE_POLOIDALPLANE ) = ind
                    end if
                end if

            end do
        end do

        call xertst( iInd == nInd, &
        &   "Assert error 2 (index) in collectIndexListForRegion" )

    end function collectIndexListForRegion

    !> Build an index list of all objects of a given region type
    !! (b2mod_connectivity.REGIONTYPE_*) for a given region id.
    !! @result The list of indices for all objects that constitute this grid
    !! region. The array has two dimensions because it is given as a list of
    !! object descriptors.
    subroutine collectIndexListForRegionSubroutine(gmap, cflag, region, iRegionType, &
        &   iRegion, indexlist)
        type(B2GridMap), intent(in) :: gmap !< The grid mapping as computed by
            !< b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS
        integer, allocatable, dimension(:,:), intent(out) :: indexList !< 2D
            !< index list
        integer, intent(in) :: region( -1:gmap%b2nx, -1:gmap%b2ny, 0:2 )
        integer, intent(in) :: iRegionType
        integer, intent(in) :: iRegion
        integer, intent(in) ::  &
            &   cflag( -1:gmap%b2nx, -1:gmap%b2ny, CARREOUT_NCELLFLAGS )

        !! Internal variables
        integer :: ix   !< x-aligned cell index
        integer :: iy   !< y-aligned cell index
        integer :: nInd
        integer :: iInd
        integer :: ind

        !! Figure out how many indices to expect. A simple count of the form
        !! nInd = count( region(:,:,iRegionType) == iRegion )
        !! will not do, because we have to account for removed objects
        !! (ghost cells/edges).

        !! search the relevant objects and count them
        nInd = 0
        do ix = -1, gmap%b2nx
            do iy = -1, gmap%b2ny

                if ( region(ix, iy, iRegionType) == iRegion ) then
                  !! Get index depending on what object type we are looking at
                    select case (iRegionType)
                    case (REGIONTYPE_CELL)
                        if ( is_Unneeded_Cell( gmap%b2nx, gmap%b2ny, cflag, &
                           & INCLUDE_GHOST_CELLS, ix, iy ) ) cycle
                        ind = gmap%mapCvI(ix, iy)
                    case (REGIONTYPE_XEDGE)
                        ind = gmap%mapFcI(ix, iy, LEFT)
                    case (REGIONTYPE_YEDGE)
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
                        if ( is_Unneeded_Cell( gmap%b2nx, gmap%b2ny, cflag, &
                           & INCLUDE_GHOST_CELLS, ix, iy ) ) cycle
                        ind = gmap%mapCvI(ix, iy)
                    case (REGIONTYPE_XEDGE)
                        ind = gmap%mapFcI(ix, iy, LEFT)
                    case (REGIONTYPE_YEDGE)
                        ind = gmap%mapFcI(ix, iy, BOTTOM)
                    end select

                    if ( ind /= B2_GRID_UNDEFINED ) then
                        iInd = iInd + 1
                        call xertst(iInd <= nInd, &
                        &   "Assert error 1 (index) in collectIndexListForRegionSubroutine" )
                        indexList( iInd, SPACE_POLOIDALPLANE ) = ind
                    end if
                end if

            end do
        end do

        call xertst( iInd == nInd, &
        &   "Assert error 2 (index) in collectIndexListForRegionSubroutine" )

    end subroutine collectIndexListForRegionSubroutine

end module b2mod_ual_io_grid

!!!Local Variables:
!!! mode: f90
!!! End:
