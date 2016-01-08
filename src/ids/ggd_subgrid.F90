module ggd_subgrid

  !> @author H.-J. Klingshirn

  use b2mod_types
  use ggd_assert
#ifdef IMAS
  use ids_schemas ! IGNORE
#else
#ifdef ITM
  use euitm_schemas ! IGNORE
#endif
#endif

  use combinations

  use ggd_common
  use ggd_access
  use ggd_object
  use ggd_objectlist

  implicit none

contains

#ifdef ITM

  !> Create a subgrid for a given number of object lists
  subroutine createSubGrid( sg, nobjlist, id )
    type(type_complexgrid_subgrid), intent(out) :: sg
    integer, intent(in) :: nobjlist
    character(*), intent(in), optional :: id


    allocate( sg%list(nobjlist) )

    allocate( sg%id(1) )
    if ( present( id ) ) then
        ! use given subgrid id
        sg % id(1) = id(1:min(len(sg%id(1)), len_trim(id)))
    else
        sg % id(1) = 'UNSPECIFIED'
    end if

  end subroutine createSubGrid


  !> Convenience routine: create a subgrid for one specific object class
  subroutine createSubGridForClass( grid, sg, cls, id )
    type(type_complexgrid), intent(in) :: grid
    type(type_complexgrid_subgrid), intent(out) :: sg
    integer, intent(in) :: cls(gridNSpace(grid))
    character(*), intent(in), optional :: id

    call createSubGrid( sg, 1, id )
    call createImplicitObjectListForAll( grid, sg%list(1), cls )
  end subroutine createSubGridForClass


  !> Convenience routine: create a subgrid for a list of object classes.
  !> @param classes The list of object class. First index: class index, 
  !>  second index: space index, i.e. classes(i,:) is the ith object class tuple.
  subroutine createSubGridForClasses( grid, sg, classes, id )
    type(type_complexgrid), intent(in) :: grid
    type(type_complexgrid_subgrid), intent(out) :: sg
    integer, intent(in) :: classes(:,:)
    character(*), intent(in), optional :: id

    ! internal
    integer :: iCls

    call assert(size(classes, 2) == gridNSpace(grid))

    ! Create one implicit object list for every object class
    call createSubGrid( sg, size(classes, 1), id )
    do iCls = 1, size(classes, 1)
        call createImplicitObjectListForAll( grid, sg%list(iCls), classes(iCls, :) )
    end do
  end subroutine createSubGridForClasses


  !> Convenience routine: create a subgrid for an explicit object list
  subroutine createSubGridForExplicitList( grid, sg, cls, indlist, id )
    type(type_complexgrid), intent(in) :: grid
    type(type_complexgrid_subgrid), intent(out) :: sg
    integer, intent(in) :: cls(gridNSpace(grid))
    integer, intent(in) :: indlist(:,:)
    character(*), intent(in), optional :: id

    call assert( size(indlist, 2) == gridNSpace(grid) )

    call createSubGrid( sg, 1, id )
    call createExplicitObjectList( grid, sg%list(1), cls, indlist )
  end subroutine createSubGridForExplicitList

  !> Return number of subgrids in the grid.
  integer function gridNSubGrid(grid)
    type(type_complexgrid), intent(in) :: grid

    gridNSubGrid = 0
    if (associated(grid%subgrids)) &
        & gridNSubGrid = size(grid%subgrids)        
  end function gridNSubGrid


  !> Returns the subgrid index for the subgrid with this name
  !> If none found, returns GRID_UNDEFINED
  integer function gridFindSubGridByName(grid, name) result (sgInd)
    type(type_complexgrid), intent(in) :: grid
    character(*), intent(in) :: name

    ! internal
    integer :: iSg
    
    do iSg = 1, gridNSubGrid(grid)       
        if (.not. associated(grid%subgrids(iSg)%id) ) cycle
        if (grid%subgrids(iSg)%id(1)(1:len(name)) == name) then
            sgInd = iSg
            return
        end if
    end do

    sgInd = GRID_UNDEFINED
  end function gridFindSubGridByName


  !> Return the number of objects in the subgrid
  ! TODO: rename to subGridSize
  integer function gridSubGridSize(sg) result( nobj )
    type(type_complexgrid_subgrid), intent(in) :: sg    

    ! internal
    integer :: ilist

    nobj = 0
    do ilist = 1, size(sg%list)
        nobj = nobj + objectListSize( sg%list(ilist) )
    end do
  end function gridSubGridSize


  !> Return the object with index iobj according to the implicit object ordering of the subgrid  
  type(GridObject) function subGridGetObject(sg, iobj) result( obj )
    type(type_complexgrid_subgrid), intent(in) :: sg  
    integer, intent(in) :: iobj

    ! internal
    integer :: ilist, listsize, offset    

    offset = 0
    do ilist = 1, size(sg%list)
        listsize = objectListSize( sg%list(ilist) )
        if ( (offset < iobj) .and. (iobj <= offset + listsize) ) then
            obj = objectListGetObject( sg%list(ilist), iobj - offset )
            return
        end if
        offset = offset + listsize
    end do

    stop "subGridGetObject: index iobj is out of range"
  end function subGridGetObject


  !> Return the local index of the given object in the subgrid, according to the
  !> implicit object ordering of the subgrid
  integer function subGridGetIndexForObject(sg, obj) result( index )
    type(type_complexgrid_subgrid), intent(in) :: sg  
    type(GridObject), intent(in) :: obj

    ! internal
    integer :: ilist, listsize, offset, localIndex

    offset = 0
    do ilist = 1, size(sg%list)
        listsize = objectListSize( sg%list(ilist) )

        localIndex = objectListGetIndexForObject(sg%list(ilist), obj)
        if (localIndex /= GRID_UNDEFINED) then
            index = offset + localIndex
            return
        end if

        offset = offset + listsize
    end do

    ! didn't find anything
    index = GRID_UNDEFINED
  end function subGridGetIndexForObject

  
  !> Add a default set of subgrids for a grid. One subgrid is added
  !> for every dimension for which objects exist in the grid. This subgrid
  !> will contain all objects of that dimension in the canonical implicit ordering.
  subroutine gridCreateDefaultSubGrids(grid, id)
    type(type_complexgrid), intent(inout) :: grid
    character(*), intent(in), optional :: id

    ! internal
    integer :: iDim, iSg
    integer, dimension(gridNSpace(grid)) :: maxCls
    integer, allocatable :: classes(:,:)

    character(7), parameter :: genId(0:6) = (/ 'nodes  ', 'edges  ', 'faces  ', 'volumes', 'obj4d  ', 'obj5d  ', 'obj6d  ' /)

    maxCls = gridGetMaxDimClass( grid )
    
    ! ...set up default subgrids: one subgrid for all objects of a given dimension
    allocate(grid%subgrids(classDim(maxCls)+1))

    do iDim = 0, classDim(maxCls)
        ! Create subgrid
        iSg = idim + 1
        call allocate_combinations(gridSpaceNDims(grid), iDim, classes)
        if (present(id)) then
            call createSubGridForClasses(grid, grid%subgrids(iSg), classes, trim(id//'_'//genId(idim)) )
        else
            call createSubGridForClasses(grid, grid%subgrids(iSg), classes, trim(genId(idim)) )
        end if
        deallocate(classes)
    end do

  end subroutine gridCreateDefaultSubGrids

#endif

end module ggd_subgrid

!!!Local Variables:
!!! mode: f90
!!! End:
