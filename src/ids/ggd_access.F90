module ggd_access

  use b2mod_types
  use ggd_assert
#ifdef IMAS
  use ids_schemas ! IGNORE
#else
#ifdef ITM
  use euITM_schemas ! IGNORE
#endif
#endif

  implicit none

  ! TODO: getSpaceNodeCount to replace size( node_value )

#ifdef ITM

contains

  !> Return the UID number of this grid
  pure integer function gridUid( grid )
    type(type_complexgrid), intent(in) :: grid

    gridUid = grid%uid
  end function gridUid

  !> Return the ID string of this grid
  character(132) function gridId(grid) 
    type(type_complexgrid), intent(in) :: grid

    gridId = repeat(' ', 132)
    if (associated( grid%id )) then
        if (len(grid%id) > 0) then
            gridId = grid%id(1)
        end if
    end if
  end function gridId

  !> Return total dimension of the grid described by a given grid descriptor.
  pure integer function gridNDim( grid )
    type(type_complexgrid), intent(in) :: grid

    ! internal
    integer :: i

    gridNDim = 0
    do i = 1, size( grid % spaces ) 
            gridNDim = gridNDim + size( grid % spaces( i ) % coordtype )
    end do

  end function gridNDim


  !> Return the number of spaces in a grid description.
  pure integer function gridNSpace( grid )
    type(type_complexgrid), intent(in) :: grid
    
    gridNSpace = size( grid % spaces ) 
  end function gridNSpace


  !> Return the dimension of an individual space.
  pure integer function gridSpaceNDim( space )
    type(type_complexgrid_space), intent(in) :: space

    gridSpaceNDim = size( space % coordtype )
  end function gridSpaceNDim


  !> Returns the dimension of all individual spaces
  pure function gridSpaceNDims( grid ) result( dims )
    type(type_complexgrid), intent(in) :: grid
    integer, dimension( size( grid % spaces ) ) :: dims

    ! internal
    integer :: i

    do i = 1, size( grid % spaces ) 
            dims( i ) = gridSpaceNDim( grid % spaces( i ) ) 
    end do

  end function gridSpaceNDims

  !> Returns the highest dimension for which objects are defined in the space
  pure function gridSpaceMaxObjDim( space ) result( dim )
    type(type_complexgrid_space), intent(in) :: space
    integer :: dim
   
    dim = 0
    if (associated(space%objects)) dim = size(space%objects) - 1
  end function gridSpaceMaxObjDim

  !> Return number of nodes (0d-objects) in the space.
  pure integer function gridSpaceNNodes( space )
    type(type_complexgrid_space), intent(in) :: space

    gridSpaceNNodes = 0
    if (associated( space % objects(1) % geo )) gridSpaceNNodes = size( space % objects(1) % geo, 1 )
  end function gridSpaceNNodes  


  !> Get the total number of objects 
  !> of the given dimension in the given space
  integer function gridSpaceNObject( space, dim ) result( objcount )
    type(type_complexgrid_space), intent(in) :: space
    integer, intent(in) :: dim

    call assert( ( dim >= 0 ) .and. ( dim <= gridSpaceNdim( space ) ), &
         & "gridSpaceNObject: dim out of bounds" )

    if ( dim == 0 ) then
        objcount = gridSpaceNNodes( space )
    else
        objcount = 0
        if (associated(space%objects)) then
            if ( dim <= size(space%objects) ) then
                objcount = size( space % objects(dim+1) % boundary, 1 )
            end if
        end if
    end if
  end function gridSpaceNObject


  !> Return maximum number of boundaries an object of dimension dim can have in the space.
  integer function gridSpaceMaxNBoundaries( space, dim ) 
    type(type_complexgrid_space), intent(in) :: space
    integer, intent(in) :: dim

    gridSpaceMaxNBoundaries = size( space % objects(dim+1) % boundary, 2 )
  end function gridSpaceMaxNBoundaries


  !> Returns the coordinate types for the individual dimensions of the grid
  !> @note Can be made pure by replacing gridSpaceNDim
  function gridCoordTypes( grid ) result( coordtype ) 
    type(type_complexgrid), intent(in) :: grid
    integer, dimension( gridNDim( grid ) ) :: coordtype

    ! internal
    integer :: is, ic, sdim

    ic = 0
    do is = 1, gridNSpace( grid )
       sdim = gridSpaceNDim( grid % spaces(is) )
       ! TODO: add treatement of multiple geometries here 
       coordtype( ic + 1 : ic + sdim ) &
            & = grid % spaces(is) % coordtype( 1 : sdim, 1 )
       ic = ic + sdim
    end do

  end function gridCoordTypes


  !> Returns index of a node according to the implicit ordering rules for the grid descriptor
  integer function gridNodeIndex( grid, nodeind ) result( index )
    type(type_complexgrid), intent(in) :: grid
    integer, dimension(:), intent(in) :: nodeind
    
    ! internal
    integer :: i, s

    call assert( size( nodeind ) == size( grid % spaces ), &
         & "gridNodeIndex: size of nodeind does not match the grid description" )
    
    index = nodeind(1)
    s = gridSpaceNNodes( grid % spaces(1) )

    do i = 2, size( grid % spaces )
            index = index + s * ( nodeind( i ) - 1 )
            s = s * gridSpaceNNodes( grid % spaces(i) ) 
    end do
    
  end function gridNodeIndex


  !> Get the coordinates of a node according to its index-tuple (assuming
  !> the default geometry representation.
  function gridNodeCoord( grid, nodeind ) result ( coord )
    type(type_complexgrid), intent(in) :: grid
    integer, dimension(gridNSpace(grid)), intent(in) :: nodeind
    real(R8), dimension( gridNDim( grid ) ) :: coord
    
    ! internal   
    integer :: is, id, nd

    call assert( size( nodeind ) == size( grid % spaces ), &
         & "gridNodeCoord: size of nodeind does not match the grid description" )

    ! FIXME: add test for default geometry representation.

    coord = 0.0_R8

    id = 0 ! coordinate counter
    do is = 1, gridNSpace( grid )
            ! get dimension of current space
            nd = gridSpaceNDim( grid % spaces(is) )
            ! copy coordinates
            ! TODO: add handling of multiple geometries here
            coord(id + 1 : id + nd ) = grid % spaces( is ) % objects(1) % geo( nodeind( is ), 1:nd, 1, 1 )
            ! increase coordinate counter
            id = id + nd
    end do

  end function gridNodeCoord
#endif

end module ggd_access

!!!Local Variables:
!!! mode: f90
!!! End:
