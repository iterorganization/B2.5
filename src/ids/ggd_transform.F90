module ggd_transform

  use b2mod_types
  use ggd_assert
#ifdef IMAS
  use ids_schemas ! IGNORE
#else
#ifdef ITM
  use euITM_schemas ! IGNORE
#endif
#endif
  use ggd_object
  use ggd_subgrid
  use combinations

  implicit none

#ifdef ITM

contains

  !> Transform a grid with possibly multiple spaces into an equivalent grid
  !> with only one space
  subroutine transformGridToSingleSpace( gIn, gOut  )
    type(type_complexgrid), intent(in) :: gin
    type(type_complexgrid), intent(out) :: gout

    ! internal
    integer :: i, iDim, iSp, nNodes, nObj, iObj, maxCompObjs, iCompObj
    type(type_complexgrid_space) :: s

    integer, dimension(gridNSpace(gIn)) :: cls, maxCls
    type(type_complexgrid_subgrid), dimension(:), allocatable :: subgrids
    integer, allocatable :: classes(:,:)
    type(GridObject) :: obj, node
    type(GridObject), allocatable :: compObjs(:)

    ! Some basic checks
    do iSp = 1, gridNSpace(gIn)
        ! Can only do simple geometry representations
        call assert(size(gIn%spaces(iSp)%objects(1)%geo,3) == 1)
    end do

    ! Set up new grid with one space
    allocate(gOut%spaces(1))
    
    ! Basic space information
    allocate(s%coordtype(gridNDim(gIn), 1))
    s%coordtype(:,1) = gridCoordTypes(gIn)

    ! Fill in object descriptions in new space s
    maxCls = gridGetMaxDimClass( gIn )
    allocate(s%objects(classDim(maxCls)+1))

    ! Nodes (object class (0,...,0))
    cls = (/ (0, i = 1, gridNSpace(gIn) ) /)
    nNodes = gridNObject( gin, cls )
    allocate(s%objects(1)%geo(nNodes, gridNDim(gIn), 1, 1))
    
    do i = 1, nNodes
        node = getObjectByGlobalIndex(gIn, cls, i)
        s%objects(1)%geo(i,:,1,1) = gridNodeCoord(gIn, node%ind)
    end do
    ! TODO for nodes: altgeo, xpoints, alias

    ! Higher-dimensional objects
    
    ! Build subgrids containing all objects of every dimension present grid
    allocate(subgrids(0:classDim(maxCls)))
    do iDim = 0, classDim(maxCls)
        call allocate_combinations(maxCls, iDim, classes)
        call createSubGridForClasses( gIn, subgrids(iDim), classes )
        deallocate(classes)
    end do

    ! Fill in the object descriptions

    do iDim = 1, classDim(maxCls)
        nObj = gridSubGridSize( subgrids(iDim) )

        ! Figure out maximum number of boundaries for objects
        maxCompObjs = 0
        do iObj = 1, nObj
            obj = subGridGetObject( subgrids(iDim), iObj )
            call getComposingObjects( gIn, obj, compObjs )
            maxCompObjs = max( maxCompObjs, size(compObjs) )
        end do

        ! Fill in boundary
        allocate(s%objects(iDim + 1)%boundary(nObj, maxCompObjs))
        s%objects(iDim + 1)%boundary = GRID_UNDEFINED

        do iObj = 1, nObj
            obj = subGridGetObject( subgrids(iDim), iObj )
            call getComposingObjects( gIn, obj, compObjs )
            do iCompObj = 1, size(compObjs)
                s%objects(iDim + 1)%boundary(iObj, iCompObj) = &
                    & subGridGetIndexForObject(subgrids(iDim-1), compObjs(iCompObj) )
            end do            
        end do

        ! TODO: fill in connectivity
    end do
   
    gOut%spaces(1) = s

  end subroutine transformGridToSingleSpace

#endif

end module ggd_transform

!!!Local Variables:
!!! mode: f90
!!! End:
