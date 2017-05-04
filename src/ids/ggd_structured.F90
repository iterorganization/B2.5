module ggd_structured

  !> Service routines for accessing structured grids and associated data

  use b2mod_types
  use ggd_assert
  use combinations

  use ggd_common
  use ggd_object
  use ggd_access
  use ggd_subgrid
!  use ggd_data
  use ggd_transform

  implicit none

  !> Definition of default subgrids
  integer, parameter :: GRID_STRUCT_SUBGRID_0D = 1 ! all 0d objects
  integer, parameter :: GRID_STRUCT_SUBGRID_1D = 2 ! all 1d objects
  integer, parameter :: GRID_STRUCT_SUBGRID_2D = 3 ! all 2d objects
  integer, parameter :: GRID_STRUCT_SUBGRID_3D = 4 ! all 3d objects
  integer, parameter :: GRID_STRUCT_SUBGRID_4D = 5 ! all 4d objects
  integer, parameter :: GRID_STRUCT_SUBGRID_5D = 6 ! all 5d objects
  integer, parameter :: GRID_STRUCT_SUBGRID_6D = 7 ! all 6d objects

  ! Same as above, with human-readable names
  integer, parameter :: GRID_STRUCT_NODES = GRID_STRUCT_SUBGRID_0D ! all 0d objects
  integer, parameter :: GRID_STRUCT_EDGES = GRID_STRUCT_SUBGRID_1D ! all 1d objects
  integer, parameter :: GRID_STRUCT_FACES = GRID_STRUCT_SUBGRID_2D ! all 2d objects  
  integer, parameter :: GRID_STRUCT_CELLS = GRID_STRUCT_SUBGRID_3D ! all 3d objects 

#ifdef ITM

  interface gridStructWriteData
     module procedure gridStructWriteData1d, gridStructWriteData2d , gridStructWriteData3d ,&
          & gridStructWriteData4d, gridStructWriteData5d, gridStructWriteData6d, &
          & gridStructWriteData1dComplex, gridStructWriteData2dComplex
  end interface

  interface gridStructReadData
     module procedure gridStructReadData1d, gridStructReadData2d, gridStructReadData3d , &
          & gridStructReadData4d, gridStructReadData5d, gridStructReadData6d, &
          & gridStructReadData1dComplex, gridStructReadData2dComplex
  end interface

  private computeMeasure

#endif

  !> Write a n-dimensional structured grid 
  !> into a grid descriptor, as well as the default subgrids for objects of all dimensions.
  !>
  !> The dimension n of the grid is taken as size(coordtype).
  !> 
  !> @param grid Grid descriptor to fill
  !> \param coordtype Dimension(n). Defines coordinate types / labels for the 
  !>                  individual axes. See the 
  !>                  constants defined in ggd_common.f90 (COORDTYPE_*)
  !> \param gshape    Dimension(n). Shape of the grid. In dimension i the grid 
  !>                  has shape(i) grid points.
  !> @param x         Dimension( maxval( gshape(n) ), n ). 
  !>                  Grid node coordinates in the individual dimensions. 
  !>                  The node positions in  space i are given by 
  !>                  x( 1 : gshape( i ), id ).
  !> @param createSubgrids Optional flag controlling whether default subgrids
  !>                  are created. Default is .true.
  !> @param periodicSpaces Optional integer array containing the indices
  !>                  of the coordinate directions that are periodic. This will
  !>                  result in the last node in these coordinate directions to
  !>                  be connected to the first node by an edge. Note that if periodic
  !>                  spaces are present, no metric information is computed.
  !> @param uid A unique identifier number for the          
  
  !> Set up a 1d structured space. 
  !>
  !> Helper routine used by gridSetupStructured. Sets up a space descriptor for the case
  !> of a simple 1d structured grid with standard connectivity

#ifdef IMAS

contains

  subroutine gridSetupStruct1dSpace( space, coordtype, nodes, periodic )
    type(ids_generic_grid_dynamic_space), intent(inout) :: space !> The space descriptor to fill
    integer, intent(in) :: coordtype !> The coordinate type of the space
    real(R8), intent(in), dimension(:) :: nodes !> The node positions in the space (assumed to be in increasing order)
    logical, intent(in), optional :: periodic

    ! internal
    integer, parameter :: NDIM = 1 ! this is a 1d grid

    logical :: lperiodic = .false.
    integer :: i, nobjects, iNb
    
    if (present(periodic)) lperiodic = periodic

#if 0
    ! Set coordinate types
    ! (dimension of space = NDIM = size( coordtype )
    allocate( space % coordtype(NDIM, 1) )    
    space % coordtype(:,1) = (/ coordtype /)

    ! Allocate object definition arrays     
    allocate( space % objects(NDIM + 1) )

    ! Set up node information
    ! Geometry: normal geometry representation, third dimension unused.
    ! TODO: support multiple geometries
    allocate( space % objects(1) % geo(size(nodes), NDIM, 1, 1) )
    space % objects(1) % geo(:, 1, 1, 1) = nodes
    ! space % nodes % xpoints, altgeo, alias: unused for this grid

    ! Allocate arrays for 1d-objects (edges)
    nobjects = size(nodes) - 1
    ! In the periodic case we have one additional edge
    if (lperiodic) nobjects = nobjects + 1

    ! Boundary array (2 boundaries per object)
    allocate( space % objects(2) % boundary( nobjects, 2 ) )
    ! connectivity array (2 boundaries per object, maximum of 1 neighbour per object)
    allocate( space % objects(2) % neighbour( nobjects, 2, 1 ) )
    
    ! Additional object geometry information (space % objects % geo) currently unused for this grid
    
    ! Fill in object definitions (i.e. what objects compose an object)
    if (nobjects > 0) then
       ! First set all to undefined   
       space % objects(2) % boundary = GRID_UNDEFINED
       ! Now fill in the nodes defining the edges
       do i = 1, nobjects
          ! left point of edge
          space % objects(2) % boundary( i, 1 ) = i
          ! right point of edge
          space % objects(2) % boundary( i, 2 ) = i + 1
          ! wrap around right edge if exceed node count in periodic case
          if (space % objects(2) % boundary( i, 2 ) > size(nodes)) then
             space % objects(2) % boundary( i, 2 ) = 1
          end if
       end do

       ! Fill in connectivity information
       ! first set all to undefined
       space%objects(2)%neighbour = GRID_UNDEFINED

       ! first edge: 

       if (.not. lperiodic) then
          ! standard case: no left neighbour
          space%objects(2)%neighbour( 1, 1, 1 ) = GRID_UNDEFINED
       else
          ! periodic case: left neighbour is last edge in space
          space%objects(2)%neighbour( 1, 1, 1 ) = nobjects
       end if
       space%objects(2)%neighbour( 1, 2, 1 ) = 2

       ! edges: left + right neighbour
       do i = 1, nobjects
          ! Left neighbour
          iNb = i - 1
          if ( iNb < 1 ) then
             iNb = GRID_UNDEFINED
             if (lperiodic) iNb = nobjects
          end if
          space%objects(2)%neighbour( i, 1, 1 ) =iNb

          ! Right neighbour
          iNb = i + 1
          if ( iNb > nobjects ) then
             iNb = GRID_UNDEFINED
             if (lperiodic) iNb = 1
          end if
          space%objects(2)%neighbour( i, 2, 1 ) = iNb
       end do
    end if


    ! Measures: store length of edges - only in non-periodic case
    if (.not. lperiodic) then
        allocate(space%objects(2)%measure( nobjects, 1 ))
        ! TODO: support multiple geometries
        if (nobjects > 0) then
           space%objects(2)%measure(:,1) = nodes(2:ubound(nodes,1)) - nodes(1:ubound(nodes,1)-1)
        end if
    end if

#endif

  end subroutine gridSetupStruct1dSpace

#else
#ifdef ITM

contains

  !> @see gridSetupStructuredSep
  subroutine gridSetupStructured( grid, coordtype, gshape, x, id, createSubgrids, periodicSpaces, uid, computeMeasures )
    type(type_complexgrid), intent(out) :: grid 
    integer, dimension(:), intent(in) :: coordtype
    integer, dimension(size(coordtype)), intent(in) :: gshape
    real(R8), dimension(:, :), intent(in) :: x
    character(*), intent(in), optional :: id
    logical, intent(in), optional :: createSubgrids
    integer, intent(in), optional :: periodicSpaces(:)
    integer, intent(in), optional :: uid
    logical, intent(in), optional :: computeMeasures 
    
    ! internal
    integer :: idim, ndim, iSg, iObj
    type(GridObject) :: obj
    integer, allocatable :: classes(:,:)
    logical :: lCreateSubgrids = .true., periodic
    logical :: lComputeMeasures = .true.

    if (present(createSubgrids)) lCreateSubgrids = createSubgrids
    if (present(computeMeasures)) lcomputeMeasures = computeMeasures

    call assert( size( coordtype ) == size( gshape ), &
         & "gridWriteStructured: size of coordtype and gshape do not match" )
    call assert( maxval( shape( x ) ) == maxval( gshape ), &
         & "gridWriteStructured: shape of x seems to be inconsistent with gshape" )

    if (present(id)) then
        allocate( grid%id(1) )
        grid%id(1) = id
    end if

    ndim = size(coordtype)
    allocate( grid % spaces(ndim) )

    ! ... fill in the grid data
    do idim = 1, ndim
        periodic = .false.
        if (present(periodicSpaces)) then
            periodic = any(periodicSpaces == idim)
        end if
        call gridSetupStruct1dSpace( grid % spaces(idim), &
            & coordtype( idim ), x( 1 : gshape(idim) ,idim ), periodic )   
    end do

    ! Create subgrids & store metric information for implicitly defined objects
    if (lCreateSubgrids) then

        call gridCreateDefaultSubGrids(grid, id)

        ! ...also set up measures for subgrids with objects of dimension > 1
        ! ...but only if there are no periodic spaces
        if (lComputeMeasures .and. .not. present(periodicSpaces) ) then

            allocate(grid%metric%measure(ndim - 1))
            
            do idim = 0, ndim
                iSg = idim + 1
                ! Store measures for implicitly defined objects (dim > 1) 
                if (idim > 1) then                
                    grid%metric%measure(idim-1)%subgrid = iSg
                    allocate( grid%metric%measure(idim-1)%scalar( gridSubGridSize(grid%subgrids(iSg)) ) )
                    ! loop over all objects in subgrid
                    do iObj = 1, gridSubGridSize( grid%subgrids(iSg) )
                        obj = subGridGetObject( grid%subgrids(iSg), iObj )
                        grid%metric%measure(idim-1)%scalar(iObj) = computeMeasure( grid, obj )
                    end do
                end if
            end do
        end if

    end if

  end subroutine gridSetupStructured

  ! Compute measure (length, area, volume, ...) for an object in a structured grid
  real(R8) function computeMeasure( grid, obj ) result( m )
    type(type_complexgrid), intent(in) :: grid
    type(GridObject), intent(in) :: obj

    ! internal
    integer :: iSp

    call assert( objectDim(obj) > 0 )

    m = 1.0_R8
    do iSp = 1,gridNSpace( grid ) 
        if (obj%cls(iSp) == 1) then
            ! TODO: add support for multiple geometries
            m = m * grid%spaces(iSp)%objects(2)%measure(obj%ind(iSp),1)
        end if
    end do
  end function computeMeasure

  !> Write a n-dimensional structured grid into a grid descriptor (alternate version 
  !> with separate dimension vectors)

  !> Alternate wrapper for gridSetupStructured, which makes it easier
  !> to give the node positions as individual arrays

  !> @see gridSetupStructured 
  subroutine gridSetupStructuredSep( grid, ndim, c1, x1, c2, x2, c3, x3, c4, x4, c5, x5, c6, x6, id, createSubgrids, periodicSpaces, uid, computeMeasures )
    type(type_complexgrid), intent(out) :: grid
    integer, intent(in) :: ndim
    real(R8), intent(in), dimension(:) :: x1 ! have to have at least one dimension
    integer, intent(in) :: c1 ! have to have at least one dimension
    real(R8), intent(in), dimension(:), optional :: x2, x3, x4, x5, x6
    integer, intent(in),optional :: c2, c3, c4, c5, c6
    character(*), intent(in), optional :: id
    logical, intent(in), optional :: createSubgrids
    integer, intent(in), optional :: periodicSpaces(:)
    integer, intent(in), optional :: uid
    logical, intent(in), optional :: computeMeasures


    ! internal
    integer :: lndim, nmax, i
    real(R8), dimension(:,:), allocatable :: x
    integer, dimension(:), allocatable :: gshape, coordtype
    

    lndim = 1
    nmax = size( x1 )

    if ( present( x2 ) ) then
            lndim = 2
            nmax = max( nmax, size( x2 ) )
            call assert( present( c2 ), "gridSetupStructuredSep: x2 given, but not c2" )
    end if
    if ( present( x3 ) ) then
            lndim = 3
            nmax = max( nmax, size( x3 ) )
            call assert( present( c3 ), "gridSetupStructuredSep: x3 given, but not c3" )
    end if
    if ( present( x4 ) ) then
            lndim = 4
            nmax = max( nmax, size( x4 ) )
            call assert( present( c4 ), "gridSetupStructuredSep: x4 given, but not c4" )
    end if
    if ( present( x5 ) ) then
            lndim = 5
            nmax = max( nmax, size( x5 ) )
            call assert( present( c5 ), "gridSetupStructuredSep: x5 given, but not c5" )
    end if
    if ( present( x6 ) ) then
            lndim = 6
            nmax = max( nmax, size( x6 ) )
            call assert( present( c6 ), "gridSetupStructuredSep: x6 given, but not c6" )
    end if

    call assert( lndim == ndim, "gridWriteStructured: error in call, ndim does not match number of arguments" )

    ! allocate and assemble temporary data structure
    allocate( x( nmax, ndim ) )
    allocate( gshape( ndim ) )
    allocate( coordtype( ndim ) )

    do i = 1, ndim

            select case( i )
            case( 1 )
                    x( 1:size( x1 ), 1 ) = x1
                    gshape(1) = size( x1 )
                    coordtype(1) = c1
            case( 2 )
                    x( 1:size( x2 ), 2 ) = x2
                    gshape(2) = size( x2 )
                    coordtype(2) = c2
            case( 3 )
                    x( 1:size( x3 ), 3 ) = x3
                    gshape(3) = size( x3 )
                    coordtype(3) = c3
            case( 4 )
                    x( 1:size( x4 ), 4 ) = x4
                    gshape(4) = size( x4 )
                    coordtype(4) = c4
            case( 5 )
                    x( 1:size( x5 ), 5 ) = x5
                    gshape(5) = size( x5 )
                    coordtype(5) = c5
            case( 6 )
                    x( 1:size( x6 ), 6 ) = x6
                    gshape(6) = size( x6 )
                    coordtype(6) = c6
            end select

    end do

    ! write
    call gridSetupStructured( grid, coordtype, gshape, x, id, createSubgrids, periodicSpaces, uid, computeMeasures )
        
  end subroutine gridSetupStructuredSep

    
  !> Set up a 1d structured space. 
  !>
  !> Helper routine used by gridSetupStructured. Sets up a space descriptor for the case
  !> of a simple 1d structured grid with standard connectivity

  subroutine gridSetupStruct1dSpace( space, coordtype, nodes, periodic )
    type(type_complexgrid_space), intent(inout) :: space !> The space descriptor to fill
    integer, intent(in) :: coordtype !> The coordinate type of the space
    real(R8), intent(in), dimension(:) :: nodes !> The node positions in the space (assumed to be in increasing order)
    logical, intent(in), optional :: periodic

    ! internal
    integer, parameter :: NDIM = 1 ! this is a 1d grid

    logical :: lperiodic = .false.
    integer :: i, nobjects, iNb
    
    if (present(periodic)) lperiodic = periodic

    ! Set coordinate types
    ! (dimension of space = NDIM = size( coordtype )
    allocate( space % coordtype(NDIM, 1) )    
    space % coordtype(:,1) = (/ coordtype /)

    ! Allocate object definition arrays     
    allocate( space % objects(NDIM + 1) )

    ! Set up node information
    ! Geometry: normal geometry representation, third dimension unused.
    ! TODO: support multiple geometries
    allocate( space % objects(1) % geo(size(nodes), NDIM, 1, 1) )
    space % objects(1) % geo(:, 1, 1, 1) = nodes
    ! space % nodes % xpoints, altgeo, alias: unused for this grid

    ! Allocate arrays for 1d-objects (edges)
    nobjects = size(nodes) - 1
    ! In the periodic case we have one additional edge
    if (lperiodic) nobjects = nobjects + 1

    ! Boundary array (2 boundaries per object)
    allocate( space % objects(2) % boundary( nobjects, 2 ) )
    ! connectivity array (2 boundaries per object, maximum of 1 neighbour per object)
    allocate( space % objects(2) % neighbour( nobjects, 2, 1 ) )
    
    ! Additional object geometry information (space % objects % geo) currently unused for this grid
    
    ! Fill in object definitions (i.e. what objects compose an object)
    if (nobjects > 0) then
       ! First set all to undefined   
       space % objects(2) % boundary = GRID_UNDEFINED
       ! Now fill in the nodes defining the edges
       do i = 1, nobjects
          ! left point of edge
          space % objects(2) % boundary( i, 1 ) = i
          ! right point of edge
          space % objects(2) % boundary( i, 2 ) = i + 1
          ! wrap around right edge if exceed node count in periodic case
          if (space % objects(2) % boundary( i, 2 ) > size(nodes)) then
             space % objects(2) % boundary( i, 2 ) = 1
          end if
       end do

       ! Fill in connectivity information
       ! first set all to undefined
       space%objects(2)%neighbour = GRID_UNDEFINED

       ! first edge: 

       if (.not. lperiodic) then
          ! standard case: no left neighbour
          space%objects(2)%neighbour( 1, 1, 1 ) = GRID_UNDEFINED
       else
          ! periodic case: left neighbour is last edge in space
          space%objects(2)%neighbour( 1, 1, 1 ) = nobjects
       end if
       space%objects(2)%neighbour( 1, 2, 1 ) = 2

       ! edges: left + right neighbour
       do i = 1, nobjects
          ! Left neighbour
          iNb = i - 1
          if ( iNb < 1 ) then
             iNb = GRID_UNDEFINED
             if (lperiodic) iNb = nobjects
          end if
          space%objects(2)%neighbour( i, 1, 1 ) =iNb

          ! Right neighbour
          iNb = i + 1
          if ( iNb > nobjects ) then
             iNb = GRID_UNDEFINED
             if (lperiodic) iNb = 1
          end if
          space%objects(2)%neighbour( i, 2, 1 ) = iNb
       end do
    end if


    ! Measures: store length of edges - only in non-periodic case
    if (.not. lperiodic) then
        allocate(space%objects(2)%measure( nobjects, 1 ))
        ! TODO: support multiple geometries
        if (nobjects > 0) then
           space%objects(2)%measure(:,1) = nodes(2:ubound(nodes,1)) - nodes(1:ubound(nodes,1)-1)
        end if
    end if

  end subroutine gridSetupStruct1dSpace

  
  !> Test whether the given grid descriptor contains a structured
  !> grid in the sense of this service module. 
  logical function gridIsStructured( grid ) result ( isStruct )
    type(type_complexgrid), intent(in) :: grid

    ! internal 
    integer :: is

    isStruct = .true.

    ! Consists of 1d spaces?
!!$    do is = 1, size( grid % spaces )
!!$       if ( gridSpaceNDim( grid%spaces(is) ) /= 1 ) isStruct = .false.
!!$    end do
    isStruct = all(gridSpaceNDims(grid) == 1)


    ! 1d spaces have default simple connectivity?
    ! TODO: implement connectivity test

  end function gridIsStructured


  !> Return the axes description of a structured grid. Essentially the equivalent read routine to gridSetupStructured.
  !> 
  !> @param grid The grid descriptor to read from
  !> @param coordtype The coordinate types of the individual axes/spaces
  !> @param gshape Number of grid nodes on the individual axes/spaces
  !> @param x The position of the grid nodes. x(i,s) is the position of node i in space s. 
  !>          All nodes in space s are given by x( 1:gshape(s), s )
  !> @see gridSetupStructured
  subroutine gridStructGetAxes( grid, coordtype, gshape, x )
    type(type_complexgrid),  intent(in) :: grid 
    integer, dimension(:), allocatable, intent(out) :: coordtype, gshape
    real(R8), dimension(:,:), allocatable, intent(out) :: x

    ! internal
    integer :: id, ndim

    call assert( gridIsStructured( grid ), "gridStructGetAxes: not a structured grid" )
    
    ndim = gridNDim( grid )

    allocate( coordtype( ndim ) )
    allocate( gshape( ndim ) )

    coordtype = gridCoordTypes( grid )
    do id = 1, ndim
       gshape( id ) = gridSpaceNNodes( grid%spaces(id) )
    end do

    allocate( x( 1 : maxval( gshape ), ndim ) )
    
    x = 0.0_R8
    do id = 1, ndim
       ! TODO: support multiple geometries
       x( 1 : gshape( id ), id ) = grid % spaces( id ) % objects(1) % geo(:, 1, 1, 1)
    end do

  end subroutine gridStructGetAxes


  !> Return the shape (number of points in every dimension) of a structured grid. 
  !> 
  !> @param grid The grid descriptor to read from
  !> @param gshape Number of grid nodes on the individual axes/spaces
  !> @see gridSetupStructured
  subroutine gridStructGetShape( grid, gshape )
    type(type_complexgrid),  intent(in) :: grid 
    integer, dimension(:), allocatable, intent(out) ::  gshape

    ! internal
    integer :: id, ndim

    call assert( gridIsStructured( grid ), "gridStructGetShape: not a structured grid" )
    
    ndim = gridNDim( grid )

    allocate( gshape( ndim ) )

    do id = 1, ndim
       gshape( id ) = gridSpaceNNodes( grid%spaces(id) )
    end do

  end subroutine gridStructGetShape


  ! Below here are some service routines to assemble object descriptors for structured grids

  !> Return an object descriptor for a cell in an n-dimensional
  !> structured grid for the canonical coordinates given in index.
  !> 
  !> @param Index of the cell (=indices of the composing faces)
  !> @return The object descriptor

  function gridStructGetCell( index ) result( object )
    type(GridObject) :: object
    integer, dimension(:), intent(in) :: index

    allocate( object % cls( size( index ) ) )
    allocate( object % ind( size( index ) ) )

    object % cls = 1
    object % ind = index

  end function gridStructGetCell

  !> Return a descriptor for a grid node with the given index. 
  !>
  !> @param index Index of the grid node
  !> @return The object descriptor

  function gridStructGetNode( index ) result( object )
    type(GridObject) :: object
    integer, dimension(:), intent(in) :: index

    allocate( object % cls( size( index ) ) )
    allocate( object % ind( size( index ) ) )

    object % cls = 0
    object % ind = index

  end function gridStructGetNode

  !> Return a descriptor for an (one-dimensional) edge, the start
  !> point of which is given by node with the given index and which
  !> extends to the next point in the increasing coordinate direction
  !> of dimension dim. 
  !> 
  !> @param index Index of starting node of the edge
  !> @param dim Index of the dimension along which the edge is aligned
  !> @return The object descriptor

  function gridStructGetEdge( index, dim ) result( object )
    type(GridObject) :: object
    integer, dimension(:), intent(in) :: index
    integer, intent(in) :: dim

    allocate( object % cls( size( index ) ) )
    allocate( object % ind( size( index ) ) )

    object % cls = 0
    object % cls( dim ) = 1
    object % ind = index

  end function gridStructGetEdge

  !> Return a descriptor for a (two-dimensional) face. The lower-right
  !> node of the face is given by index, and the face extends along the 
  !> two coordinate directions given in dims.
  !>
  !> @param index Index of lower-right node.
  !> @param dims The dimensions along which the face extends
  !> @return The object descriptor

  function gridStructGetFace( index, dims ) result( object )
    type(GridObject) :: object
    integer, dimension(:), intent(in) :: index
    integer, intent(in), dimension(2) :: dims

    allocate( object % cls( size( index ) ) )
    allocate( object % ind( size( index ) ) )

    object % cls = 0
    object % cls( dims(1) ) = 1
    object % cls( dims(2) ) = 1
    object % ind = index

  end function gridStructGetFace

  ! Note: in the gridStructWrite routines, grid is currently unused but might be used in the future

  subroutine gridStructWriteData1d( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar), intent(inout) :: cpofield
    real(R8), dimension(:), intent(in) :: data

    call gridWriteDataScalar( cpofield, subgrid, reshape(data, (/ size( data ) /)) )   
  end subroutine gridStructWriteData1d

  subroutine gridStructWriteData2d( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar), intent(inout) :: cpofield
    real(R8), dimension(:,:), intent(in) :: data

    call gridWriteDataScalar( cpofield, subgrid, reshape(data, (/ size( data ) /)) )  
  end subroutine gridStructWriteData2d

  subroutine gridStructWriteData3d( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar), intent(inout) :: cpofield
    real(R8), dimension(:,:,:), intent(in) :: data

    call gridWriteDataScalar( cpofield, subgrid, reshape(data, (/ size( data ) /)) )   
  end subroutine gridStructWriteData3d

  subroutine gridStructWriteData4d( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar), intent(inout) :: cpofield
    real(R8), dimension(:,:,:,:), intent(in) :: data

    call gridWriteDataScalar( cpofield, subgrid, reshape(data, (/ size( data ) /)) )   
  end subroutine gridStructWriteData4d

  subroutine gridStructWriteData5d( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar), intent(inout) :: cpofield
    real(R8), dimension(:,:,:,:,:), intent(in) :: data

    call gridWriteDataScalar( cpofield, subgrid, reshape(data, (/ size( data ) /)) )   
  end subroutine gridStructWriteData5d

  subroutine gridStructWriteData6d( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar), intent(inout) :: cpofield
    real(R8), dimension(:,:,:,:,:,:), intent(in) :: data

    call gridWriteDataScalar( cpofield, subgrid, reshape(data, (/ size( data ) /)) )   
  end subroutine gridStructWriteData6d

  subroutine gridStructWriteData1dComplex( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar_cplx), intent(inout) :: cpofield
    complex(R8), dimension(:), intent(in) :: data

    call gridWriteDataScalarComplex( cpofield, subgrid, reshape(data, (/ size( data ) /)) )   
  end subroutine gridStructWriteData1dComplex

  subroutine gridStructWriteData2dComplex( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar_cplx), intent(inout) :: cpofield
    complex(R8), dimension(:,:), intent(in) :: data

    call gridWriteDataScalarComplex( cpofield, subgrid, reshape(data, (/ size( data ) /)) )   
  end subroutine gridStructWriteData2dComplex

  !> Body of the data read routine for data arrays with arbitrary rank
  subroutine gridStructReadDataBody( grid, cpofield, subgrid, gshape )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar), intent(in) :: cpofield
    integer, dimension(:),  intent(in) :: gshape

    ! internal
    integer :: ndata, i

    ! check whether inputs make basic sense
    call assert( gridIsStructured( grid ), "gridStructReadDataBody: grid not structured" )
    call assert( cpofield%subgrid == subgrid, "gridStructReadDataBody: subgrid does not match" )
    call assert( associated( cpofield%scalar ), "gridStructReadDataBody: cpofield%scalar not associated" )
    
    ndata = product( gshape )
    call assert( ndata == gridSubGridSize(grid%subgrids(subgrid)), &
        & "gridStructReadDataBody: data size does not match subgrid size" )

    ! Actual reading the data is deferred to the rank-specific routines, to avoid unnecessary copying
    
  end subroutine gridStructReadDataBody

  ! FIXME: cpofield -> cpodata
  subroutine gridStructReadData1d( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar), intent(in) :: cpofield
    real(R8), dimension(:), intent(out) :: data

    call gridStructReadDataBody( grid, cpofield, subgrid, shape(data) )
    data = reshape( cpofield%scalar, shape(data) )   
  end subroutine gridStructReadData1d

  subroutine gridStructReadData2d( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar), intent(in) :: cpofield
    real(R8), dimension(:,:), intent(out) :: data

    call gridStructReadDataBody( grid, cpofield, subgrid, shape(data) )
    data = reshape( cpofield%scalar, shape(data) )   
  end subroutine gridStructReadData2d

  subroutine gridStructReadData3d( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar), intent(in) :: cpofield
    real(R8), dimension(:,:,:), intent(out) :: data

    call gridStructReadDataBody( grid, cpofield, subgrid, shape(data) )     
    data = reshape( cpofield%scalar, shape(data) )   
  end subroutine gridStructReadData3d

  subroutine gridStructReadData4d( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar), intent(in) :: cpofield
    real(R8), dimension(:,:,:,:), intent(out) :: data

    call gridStructReadDataBody( grid, cpofield, subgrid, shape(data) )     
    data = reshape( cpofield%scalar, shape(data) )   
  end subroutine gridStructReadData4d

  subroutine gridStructReadData5d( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar), intent(in) :: cpofield
    real(R8), dimension(:,:,:,:,:), intent(out) :: data

    call gridStructReadDataBody( grid, cpofield, subgrid, shape(data) )     
    data = reshape( cpofield%scalar, shape(data) )   
  end subroutine gridStructReadData5d

  subroutine gridStructReadData6d( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar), intent(in) :: cpofield
    real(R8), dimension(:,:,:,:,:,:), intent(out) :: data

    call gridStructReadDataBody( grid, cpofield, subgrid, shape(data) )     
    data = reshape( cpofield%scalar, shape(data) )   
  end subroutine gridStructReadData6d


  !> Body of the data read routine for data arrays with arbitrary rank
  subroutine gridStructReadDataBodyComplex( grid, cpofield, subgrid, gshape )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar_cplx), intent(in) :: cpofield
    integer, dimension(:),  intent(in) :: gshape

    ! internal
    integer :: ndata, i

    ! check whether inputs make basic sense
    call assert( gridIsStructured( grid ), "gridStructReadDataBodyComplex: grid not structured" )
    call assert( cpofield%subgrid == subgrid, "gridStructReadDataBodyComplex: subgrid does not match" )
    call assert( associated( cpofield%scalar ), "gridStructReadDataBodyComplex: cpofield%scalar not associated" )
    
    ndata = product( gshape )
    call assert( ndata == gridSubGridSize(grid%subgrids(subgrid)), &
        & "gridStructReadDataBodyComplex: data size does not match subgrid size" )

    ! Actual reading the data is deferred to the rank-specific routines, to avoid unnecessary copying
    
  end subroutine gridStructReadDataBodyComplex

  ! FIXME: cpofield -> cpodata
  subroutine gridStructReadData1dComplex( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar_cplx), intent(in) :: cpofield
    complex(R8), dimension(:), intent(out) :: data

    call gridStructReadDataBodyComplex( grid, cpofield, subgrid, shape(data) )
    data = reshape( cpofield%scalar, shape(data))
  end subroutine gridStructReadData1dComplex

  subroutine gridStructReadData2dComplex( grid, cpofield, subgrid, data )
    type(type_complexgrid),  intent(in) :: grid
    integer, intent(in) :: subgrid
    type(type_complexgrid_scalar_cplx), intent(in) :: cpofield
    complex(R8), dimension(:,:), intent(out) :: data

    call gridStructReadDataBodyComplex( grid, cpofield, subgrid, shape(data) )
    data = reshape( cpofield%scalar, shape(data) )
  end subroutine gridStructReadData2dComplex


  ! Routines for structured curvilinear grids

  subroutine gridSetupStructuredCurvilinear3d( grid, coordtype, x, y, z, id, createSubgrids)
    type(type_complexgrid), intent(out) :: grid 
    integer, dimension(3), intent(in) :: coordtype
    real(R8), dimension(:,:,:), intent(in) :: x, y ,z
    character(*), intent(in), optional :: id
    logical, intent(in), optional :: createSubgrids

    ! internal
    type(type_complexgrid) :: structgrid
    !real(R8), dimension(maxval(shape(x)), 3) :: structx
    real(R8), dimension(:, :), allocatable :: structx
    integer :: ix, iy, iz, iNode
    logical :: lCreateSubgrids = .true., periodic

    if (present(createSubgrids)) lCreateSubgrids = createSubgrids

    !write (*,*) shape(x)
    !write (*,*) maxval(shape(x))
    allocate( structx(maxval(shape(x)), 3) )

    ! First set up a 3d structured grid
    call gridSetupStructured( structgrid, coordtype, shape(x), structx, id, createSubgrids=.true. )

    ! transform it into an equivalent grid with one space
    call transformGridToSingleSpace( structgrid, grid )

    ! Now set the node positions properly. The nodes are in "fortran order", i.e. lowest dimension 
    ! iterating the fastest
    iNode = 0
    do iz = 1, size(x, 3)
        do iy = 1, size(x, 2)
            do ix = 1, size(x, 1)
                iNode = iNode + 1
                ! TODO: support multiple geometry types
                grid%spaces(1)%objects(1)%geo(iNode, 1, 1, 1) = x(ix, iy, iz)
                grid%spaces(1)%objects(1)%geo(iNode, 2, 1, 1) = y(ix, iy, iz)
                grid%spaces(1)%objects(1)%geo(iNode, 3, 1, 1) = z(ix, iy, iz)
            end do
        end do
    end do

    ! Create subgrids 
    if (lCreateSubgrids) then
        call gridCreateDefaultSubGrids(grid, id)
        ! TODO: what about the metrics?
    end if

  end subroutine gridSetupStructuredCurvilinear3d

#endif
#endif

end module ggd_structured

!!!Local Variables:
!!! mode: f90
!!! End:
