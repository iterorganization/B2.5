module ggd_objectlist

  use b2mod_types
  use ggd_assert
  use string
#ifdef IMAS
  use ids_schemas ! IGNORE
#else
#ifdef ITM
  use euitm_schemas ! IGNORE
  use ggd_common
  use ggd_access
  use ggd_object
#endif
#endif

  implicit none

#ifdef ITM

contains

  ! Routines for index sets

  !> Basic setup of structure (required because all arrays in 
  !> the ITM CPO data structures are pointers, we cannot have fixed size arrays.
  subroutine allocateIndexList( iList )
    type(type_complexgrid_indexlist), intent(inout) :: iList

    allocate(ilist%range(2))
  end subroutine allocateIndexList

  !> Create an index set: single object
  type(type_complexgrid_indexlist) function createIndexListForSingle( index ) result ( ilist )
    integer, intent(in) :: index

    call allocateIndexList(ilist)
    ilist%range(1) = index
    ilist%range(2) = index
  end function createIndexListForSingle


  !> Create an index set: all objects of given class
  type(type_complexgrid_indexlist) function createIndexListForAll( grid, cls, ispace ) result ( ilist )
    type(type_complexgrid), intent(in) :: grid
    integer, intent(in) :: cls(gridNSpace(grid)), ispace

    call allocateIndexList(ilist)
    ilist%range(1) = 1
    ilist%range(2) = gridSpaceNObject( grid%spaces(ispace), cls(ispace) )
  end function createIndexListForAll


  !> Create an index list: index range
  type(type_complexgrid_indexlist) function createIndexListForRange( ifrom, ito ) result ( ilist )
    integer, intent(in) ::  ifrom, ito

    call allocateIndexList(ilist)
    ilist%range(1) = ifrom
    ilist%range(2) = ito    
  end function createIndexListForRange

  !> Create an index list: explicit list of indices
  type(type_complexgrid_indexlist) function createIndexListForList( grid, nind, ind ) result ( ilist )
    type(type_complexgrid), intent(in) :: grid
    integer, intent(in) :: nind
    integer, intent(in), optional :: ind(nind)

    ilist%range = GRID_UNDEFINED

    allocate( ilist%ind(nind) ) 
    if ( present(ind) ) then
        ilist%ind = ind
    else
        ilist%ind = GRID_UNDEFINED
    end if
  end function createIndexListForList

  !> Get size of the index list
  integer function indexListSize( ilist ) result( icount )
    type(type_complexgrid_indexlist), intent(in) :: ilist

    if ( associated( ilist%ind ) ) then
        icount = size( ilist%ind )
    else
        icount = ilist%range(2) - ilist%range(1) + 1
    end if

  end function indexListSize

  !> Get the index specified in the index list for the local index localind
  integer function indexListGetIndex( ilist, localind ) result (ind)
    type(type_complexgrid_indexlist), intent(in) :: ilist
    integer, intent(in) :: localind

    call assert( (localind >= 1) .and. (localind <= indexListSize(ilist) ), &
        & "indexListGetIndex: index localind out of range" )

    if ( associated( ilist%ind ) ) then
        ind = ilist%ind(localind)
    else
        ind = ilist%range(1) + localind - 1
    end if

  end function indexListGetIndex

  !> Get the local index of the index ind in this index list. The the index
  !> is not contained in the index list, returns GRID_UNDEFINED
  integer function indexListGetLocalIndex( ilist, ind ) result (localind)
    type(type_complexgrid_indexlist), intent(in) :: ilist
    integer, intent(in) :: ind

    ! internal
    integer :: i

    if ( associated( ilist%ind ) ) then
        ! do a search
        do i = 1, size(ilist%ind)
            localind = i
            if (ilist%ind(i) == ind) return
        end do
    else
        ! index is in range
        if ( ( ilist%range(1) <= ind ) .and. ( ilist%range(2) >= ind )  ) then
            localind = ind - ilist%range(1) + 1
            return
        end if
    end if

    localind = GRID_UNDEFINED

  end function indexListGetLocalIndex


  ! ObjectList service routines

  !> Create an explicit object list for a given class and given objects (optional: directly set the objects)
  subroutine createExplicitObjectList( grid, olist, cls, ind )
    type(type_complexgrid), intent(in) :: grid
    type(type_complexgrid_objectlist), intent(out) :: olist
    integer, intent(in) :: cls(gridNSpace(grid))
    integer, intent(in), optional :: ind(:,:)

    call assert( size(ind, 2) == gridNSpace(grid) )

    allocate( olist%cls(size(cls)) )
    olist%cls = cls

    allocate( olist%ind(size(ind, 1), size(cls)) )    
    if ( present(ind) ) then
        olist%ind = ind
    else
        olist%ind = GRID_UNDEFINED
    end if
  end subroutine createExplicitObjectList


  !> Create an implicit object list for a given class. This routine 
  !> does not fill in any object indices, this has to be done in a separate step.
  subroutine createImplicitObjectList( grid, olist, cls )
    type(type_complexgrid), intent(in) :: grid
    type(type_complexgrid_objectlist), intent(out) :: olist
    integer, intent(in) :: cls(gridNSpace(grid))

    ! internal
    integer :: i

    allocate( olist%cls(size(cls)) )
    olist%cls = cls

    allocate( olist%indset(gridNSpace(grid)) )
    
  end subroutine createImplicitObjectList

  !> Create an implicit object list for a given class, all objects
  subroutine createImplicitObjectListForAll( grid, olist, cls )
    type(type_complexgrid), intent(in) :: grid
    type(type_complexgrid_objectlist), intent(out) :: olist
    integer, intent(in) :: cls(gridNSpace(grid))

    ! internal
    integer :: ispace

    call createImplicitObjectList( grid, olist, cls )
    do ispace = 1, size(cls)
        olist%indset(ispace) = createIndexListForAll( grid, cls, ispace )
    end do

  end subroutine createImplicitObjectListForAll

  !> Check whether an object list is an explicit list of objects.
  logical function objectListIsExplicit( olist )
    type(type_complexgrid_objectlist), intent(in) :: olist

    if ( .not. associated(olist%ind) ) then
        objectListIsExplicit = .false.
        return
    end if
    objectListIsExplicit = (size(olist%ind) > 0)

  end function objectListIsExplicit

  !> Get size of the index set
  integer function objectListSize( olist ) result( nobj )
    type(type_complexgrid_objectlist), intent(in) :: olist

    ! internal
    integer :: is

    if ( objectListIsExplicit(olist) ) then
        nobj = size(olist%ind, 1)
    else
        nobj = 1
        do is = 1, size(olist%cls) ! 1, nspace
            nobj = nobj * indexListSize( olist%indset(is) )
        end do
    end if
  end function objectListSize

  !> Get index tuple for the object with local index lind in this list
  function objectListGetIndex( olist, localind ) result( ind )
    type(type_complexgrid_objectlist), intent(in) :: olist
    integer, intent(in) :: localind
    integer, dimension(size(olist%cls)) :: ind

    ! internal
    integer :: is, irem, s, tind, nspace
    integer, dimension(size(olist%cls)) :: objcount

    nspace = size(olist%cls)

    if ( objectListIsExplicit(olist) ) then
        ind = olist%ind(localind,:)
    else
        ! set up array with list-local object counts for the individual spaces
        do is = 1, nspace
            objcount( is ) = indexListSize( olist%indset(is) )
        end do

        ! compute the local index tuple in this implicit object list
        irem = localind ! remaining index
        do is = nspace, 2, -1
            s = product( objcount( 1 : is - 1 ) )
            tind = ( ( irem - 1 ) / s )
            irem = irem - tind * s
            ind( is ) = tind + 1
        end do
        ind( 1 ) = irem

        ! Convert to grid-global index tuple
        do is = 1, nspace
            ind(is) = indexListGetIndex( olist%indset(is), ind(is) )
        end do
    end if

  end function objectListGetIndex

  !> Get ith object specified in the object list
  type(GridObject) function objectListGetObject( olist, localind ) result( obj )
    type(type_complexgrid_objectlist), intent(in) :: olist
    integer, intent(in) :: localind

    obj = getObject( olist%cls, objectListGetIndex( olist, localind ) )
  end function objectListGetObject


  !> Return the local index of the given object obj in the object list olist,
  !> according to the object ordering in the object list.
  !> If the object is not contained in the object list, GRID_UNDEFINED is returned.
  integer function objectListGetIndexForObject( olist, obj ) result( ind )
    type(type_complexgrid_objectlist), intent(in) :: olist
    type(GridObject), intent(in) :: obj

    ! internal
    integer :: is, s, iobj, nspace
    integer, dimension(size(olist%cls)) :: localInds, objcount

    ind = GRID_UNDEFINED

    ! Does the class match?   
    if (.not. all(olist%cls == obj%cls)) return

    nspace = size(olist%cls)

    if ( objectListIsExplicit(olist) ) then
        ! Do a linear search
        do iobj = 1, size(olist%ind, 1)
            if (all(olist%ind(iobj, :) == obj%ind)) then
                ind = iobj
                return
            end if
        end do
    else
        ! Look up list-local index tuple 
        do is = 1, nspace
            localInds(is) = indexListGetLocalIndex( olist%indset(is), obj%ind(is) )
        end do     
        ! If any of the local indices is GRID_UNDEFINED, the object is not in this object list
        if (any(localInds == GRID_UNDEFINED)) return

        ! Now compute the local index of the object in this object list

        ! set up array with list-local object counts for the individual spaces
        do is = 1, nspace
            objcount( is ) = indexListSize( olist%indset(is) )
        end do

        ! As always, the leftmost index in the list-local index tuple changes the fastest
        ind = localInds(1)
        s = objcount(1)
        do is = 2, nspace
            ind = ind + s * (localInds(is) - 1)
            s = s * objcount(is)
        end do
    end if

  end function objectListGetIndexForObject

#endif

end module ggd_objectlist

!!!Local Variables:
!!! mode: f90
!!! End:
