module ggd_object

  use b2mod_types
  use ggd_assert
#ifdef IMAS
  use ids_schemas ! IGNORE
#else
#ifdef ITM
  use euitm_schemas ! IGNORE
  use ggd_common
  use ggd_access
#endif
#endif

  implicit none

  ! TODO: rename GridObject -> GridObject
  type GridObject
     ! cls: class of object. Gives Dimension of subobjects / elemental objects in the individual spaces
     ! ind: indices of subobjects in the individual spaces
     ! Wanted to call the components class and index, but both are Fortran keywords/intrinsics. 
     ! Therefore the short names. Also shorter to type.
     integer, dimension(:), allocatable :: cls, ind
  end type GridObject

   INTERFACE OPERATOR (.eq.)
       module procedure objectsEqual
   END INTERFACE

   interface assignment (=)
       module procedure objectAssign
   end interface

contains

  !> Test whether two GridObject strucutres are equal, 
  !> i.e. they describe the same grid object
  logical function objectsEqual( objA, objB )
    type(GridObject), intent(in) :: objA, objB

    objectsEqual = all(objA%cls == objB%cls) &
        & .and. all(objA%ind == objB%ind)
  end function objectsEqual

  !> Assign one object variable to another object variable.
  !> Used to extend the assignment operator (=).
  !>
  !> Note: one should not have to define this explicitly,
  !> as this should be the default behaviour of the compiler for 
  !> allocatable components. However, the PGI compiler (v10) fails at this.
  subroutine objectAssign(out, in)
    type(GridObject), intent(out) :: out
    type(GridObject), intent(in) :: in
    
    if (allocated(in%cls)) then
        allocate(out%cls(size(in%cls)))
        allocate(out%ind(size(in%ind)))
        out%cls = in%cls
        out%ind = in%ind
    end if
  end subroutine objectAssign

#ifdef ITM

  !> Get class descriptor for the highest-dimensional objects defined in this grid.
  !> (Depending on how the grid is defined, this is not necessarily the grid dimension.)
  function gridGetMaxDimClass(grid) result (cls)
    type(type_complexgrid), intent(in) :: grid
    integer, dimension( gridNSpace( grid ) ) :: cls

    ! internal
    integer :: iSp

    do iSp = 1, gridNSpace(grid)
        cls(iSp) = gridSpaceMaxObjDim(grid%spaces(iSp))
    end do
  end function gridGetMaxDimClass

  !> Return the dimension of the objects of class cls
  integer function classDim(cls) result(dim)
    integer, dimension(:), intent(in) :: cls
    
    dim = sum(cls)
  end function classDim

  !> Return GridObject for given class and index tuple
  function getObject( cls, ind  ) result( object )
    type(GridObject) :: object
    integer, dimension(:), intent(in) :: cls, ind

    call assert( size( cls ) == size( ind ), &
         & "getObject: size of cls does not match size of ind" )

    allocate( object % cls( size( cls ) ) )
    allocate( object % ind( size( ind ) ) )

    object % cls = cls
    object % ind = ind
  end function getObject


  !> Test whether the given object is a grid node (i.e. a zero-dimensional object)
  !> @return True if this is the case, false if not.
  logical function objectIsNode( obj ) 
    type(GridObject), intent(in) :: obj

    objectIsNode = ( sum( abs( obj % cls ) ) == 0 )
  end function objectIsNode

  !> Return the dimension of an object.
  integer function objectDim( obj )
    type(GridObject), intent(in) :: obj

    objectDim = sum(obj%cls)
  end function objectDim


  !> Get the total number of objects of the given cls in the grid
  integer function gridNObject( grid, cls ) result( objcount )
    type(type_complexgrid), intent(in) :: grid
    integer, dimension( gridNSpace( grid ) ), intent(in) :: cls

    ! internal
    integer :: is

    objcount = 1

    do is = 1, gridNSpace( grid )
       objcount = objcount * gridSpaceNObject( grid % spaces( is ), cls( is ) )
    end do

  end function gridNObject


  !> Return the global index of an object
  integer function objectGlobalIndex( grid, object ) result( ind )
    type(type_complexgrid), intent(in) :: grid
    type(GridObject), intent(in) :: object   

    ! internal
    integer :: is, s

    ind = object % ind( 1 )
    s = gridSpaceNObject( grid % spaces( 1 ), object % cls( 1 ) )

    do is = 2, gridNSpace( grid )
       ind = ind + s * ( object % ind( is ) - 1 )
       s = s * gridSpaceNObject( grid % spaces( is ), object % cls( is ) )
    end do

  end function objectGlobalIndex


  !> TODO: rename to gridGetObjectByGlobalIndex

  type(GridObject) function getObjectByGlobalIndex( grid, cls, ind ) result( object )
    type(type_complexgrid), intent(in) :: grid
    integer, dimension(gridNSpace( grid )), intent(in) :: cls
    integer, intent(in) :: ind

    ! internal
    integer :: is, s, tind, irem, oc( gridNSpace( grid ) )

    call assert( size( cls ) == gridNSpace( grid ), &
         & "getObjectByGlobalIndex: size of cls does not match grid description" )

    allocate( object % cls( gridNSpace( grid ) ) )
    allocate( object % ind( gridNSpace( grid ) ) )

    object % cls = cls

    do is = 1, gridNSpace( grid )
       oc( is ) = gridSpaceNObject( grid % spaces( is ), cls( is ) )
    end do

    irem = ind
    do is = gridNSpace( grid ), 2, -1
       s = product( oc( 1 : is - 1 ) )
       tind = ( ( irem - 1 ) / s )
       irem = irem - tind * s
       object % ind( is ) = tind + 1
    end do
    object % ind( 1 ) = irem

  end function getObjectByGlobalIndex


  !> For a given object descriptor, retrieves a list that holds the object descriptors
  !> describing the lower-dimensional objects composing the given object
  subroutine getComposingObjects( grid, object, objlist )
    type(type_complexgrid), intent(in) :: grid
    type(GridObject), intent(in) :: object
    type(GridObject), dimension(:), allocatable, intent(out) :: objlist
    
    integer :: is, cobj, iobj, nobjs, nbounds( size( grid % spaces ) )

    ! Corner case: the composing object of a node (0d object) is the node itself
    if ( objectIsNode( object ) ) then
       allocate( objlist(1) )
       objlist(1) = object
       return
    end if   

    ! General case: figure out how many objects are returned
    nobjs = 0

    do is = 1, size( grid % spaces )
       ! every space contributes objects according to the boundaries of the subobject in the spaces

       ! if subobject in the space is a node, skip
       if ( object % cls( is ) == 0 ) cycle

       ! get number of boundaries of subobject
       do iobj = 1, gridSpaceMaxNBoundaries( grid % spaces( is ), object % cls( is ) )
          if ( grid % spaces( is ) % objects( object % cls( is ) + 1 ) % &
              & boundary( object % ind( is ), iobj ) == GRID_UNDEFINED ) then
             exit
          end if
          nbounds( is ) = iobj
       end do

       nobjs = nobjs + nbounds( is )
    end do

    allocate( objlist( nobjs ) )

    ! For every space: get composing objects (boundaries) of the subobject in this space
    ! and transform them into the requested composing objects.
    ! If the subobject in the space is a node, skip the space

    cobj = 0

    do is = 1, size( grid % spaces )

       ! if subobject in this space is a node, cycle
       if ( object % cls( is ) == 0 ) cycle

       ! for all boundary subobjects...
       do iobj = 1, nbounds( is )
             
          cobj = cobj + 1
          ! assert: cobj < size( objlist )
          
          ! copy current object (automatic allocation)...
          objlist( cobj ) = object 
          ! ...reduce dimension...
          objlist( cobj ) % cls( is ) = objlist( cobj ) % cls( is ) - 1
          ! ...and set index
          objlist( cobj ) % ind( is ) =  grid % spaces( is ) % &
              & objects( object % cls( is ) + 1 ) % boundary( object % ind( is ), iobj )

       end do

    end do

    call assert( cobj == nobjs, 'getComposingObjects: found less objects than anticipated' )

  end subroutine getComposingObjects


  ! Output routines

  !> Write a list of object descriptors to stdout. 
  subroutine writeObjectList( objs )
    type(GridObject), dimension(:), intent(in) :: objs

    ! internal
    integer :: i

    do i = 1, size( objs )
      call gridWriteGridObject( objs(i) )
    end do

  end subroutine writeObjectList

  !> Write an object descriptor to stdout.
  subroutine gridWriteGridObject( obj )
    type(GridObject), intent(in) :: obj

    write (*,*) '( (', obj % cls, ') (', obj % ind, ') )'

  end subroutine gridWriteGridObject

#endif

end module ggd_object

!!!Local Variables:
!!! mode: f90
!!! End:
