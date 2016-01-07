module b2mod_ual_io_data

  ! This module provides:
  !
  ! 2. A routine (b2ITMFillGridDescription) to write the B2 grid into an ITM 
  ! grid description data structure (which usually is part of a CPO). It also
  ! sets up the default subgrids for the B2 grid.
  ! 
  ! 3. Routines to transform variables stored in the B2 data structure into
  ! the form expected CPO data structure.

  use b2mod_types , B2_R8 => R8, B2_R4 => R4
#ifdef ITM
  use itm_types , ITM_R8 => R8, ITM_R4 => R4 ! IGNORE
  use euITM_schemas ! IGNORE
  use itm_grid , ITM_GRID_UNDEFINED => GRID_UNDEFINED ! IGNORE
  use itm_string ! IGNORE
  use itm_assert ! IGNORE
  use itm_grid_structured , only: gridSetupStruct1dSpace ! IGNORE
  use itm_constants , pi => itm_pi ! IGNORE
#endif
  use helper
  use logging , only: logmsg, LOGDEBUG
  use b2mod_connectivity , REMOVED_B2_R8 => R8
  use carre_constants
  use b2mod_cellhelper

!xpb  use b2mod_grid_mapping
  use b2mod_ual_io_grid

#ifdef ITM

  implicit none

  interface b2ITMTransformDataB2ToCpo
      module procedure b2ITMTransformDataB2ToCPOCell, b2ITMTransformDataB2ToCPOFace
  end interface

contains


  ! Below here: service routines to transform data from B2 to CPO

  function b2ITMTransformDataB2ToCPOCell( grid, subgridId, gmap, b2CellData ) result( cpodata )
    real(ITM_R8), dimension(:), pointer :: cpodata

    type(type_complexgrid), intent(in) :: grid
    integer, intent(in) :: subgridId
    type(B2ITMGridMap), intent(in) :: gmap
    real(ITM_R8), intent(in) :: b2CellData(-1:gmap%b2nx, -1:gmap%b2ny)

    cpodata => b2ITMTransformDataB2ToCPOGeneral( grid, subgridId, gmap, b2CellData = b2CellData )
  end function b2ITMTransformDataB2ToCPOCell

  function b2ITMTransformDataB2ToCPOFace( grid, subgridId, gmap, b2FaceData ) result( cpodata )
    real(ITM_R8), dimension(:), pointer :: cpodata

    type(type_complexgrid), intent(in) :: grid
    integer, intent(in) :: subgridId
    type(B2ITMGridMap), intent(in) :: gmap
    real(ITM_R8), intent(in) :: b2FaceData(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)

    cpodata => b2ITMTransformDataB2ToCPOGeneral( grid, subgridId, gmap, b2FaceData = b2FaceData )
  end function b2ITMTransformDataB2ToCPOFace

! TODO: find a way to include this subroutine in the b2ITMTransformDataB2ToCPO interface

  function b2ITMTransformDataB2ToCPOVertex( grid, subgridId, gmap, b2VertexData ) result( cpodata )
    real(ITM_R8), dimension(:), pointer :: cpodata

    type(type_complexgrid), intent(in) :: grid
    integer, intent(in) :: subgridId
    type(B2ITMGridMap), intent(in) :: gmap
    real(ITM_R8), intent(in) :: b2VertexData(-1:gmap%b2nx, -1:gmap%b2ny)

    cpodata => b2ITMTransformDataB2ToCPOGeneral( grid, subgridId, gmap, b2VertexData = b2VertexData )
  end function b2ITMTransformDataB2ToCPOVertex

!> Transform a quantity stored on faces from a 2d B2 array into a 1d CPO array for a given
!> subgrid id. Either b2CellData or b2FaceData must be given. Do not use this directly,
!> use the provided general interface b2ITMTransformDataB2ToCPO instead.
!>
!> @param grid The CPO grid description
!> @param subgridId Id of the subgrid the data is to be stored on.
!> @param gmap The grid mapping as computed by b2ITMCreateMap
!> @param b2CellData Cell data given on the 2d b2 data structure
!> @param b2FaceData Face data given on the 2d b2 data structure

  function b2ITMTransformDataB2ToCPOGeneral( grid, subgridId, gmap, b2CellData, b2FaceData, b2VertexData ) result( cpodata )
    real(ITM_R8), dimension(:), pointer :: cpodata

    type(type_complexgrid), intent(in) :: grid
    integer, intent(in) :: subgridId
    type(B2ITMGridMap), intent(in) :: gmap
    real(ITM_R8), intent(in), optional :: b2CellData(-1:gmap%b2nx, -1:gmap%b2ny)
    real(ITM_R8), intent(in), optional :: b2FaceData(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)
    real(ITM_R8), intent(in), optional :: b2VertexData(-1:gmap%b2nx, -1:gmap%b2ny)

    ! internal
    integer :: nobjs, iobj, ifc, icv, ivx
    type(GridObject) :: curObj

    ! .neqv. is xor (exclusive or) 
    ! (http://de.wikibooks.org/wiki/Fortran:_Fortran_95:_Logische_Ausdr%C3%BCcke)
    ! TODO: FIX
    !call assert( present(b2CellData) .neqv. present(b2FaceData) )

    ! Allocate result vector according to subgrid size
    nobjs = gridSubGridSize( grid%subgrids(subgridId) )
    allocate( cpodata(nobjs) )

    ! Collect all data items for the subgrid objects
    do iobj = 1, nobjs
        ! Get the object descriptor
        curObj = subGridGetObject( grid%subgrids(subgridId), iobj )

        if (present(b2CellData)) then
            ! Cell data case
            ! check that it's a cell
            call assert( all( curObj%cls == CLASS_CELL(1:SPACE_COUNT) ) )
            ! get the subobject index for the face in the 2d poloidal plane space
            icv = curObj%ind(SPACE_POLOIDALPLANE)
            ! copy data
            cpodata(iobj) = b2CellData( gmap%mapCvix(icv), gmap%mapCviy(icv) )
        else if (present(b2FaceData)) then
            ! Face data case
            ! check that it's a face
            call assert( all( curObj%cls == CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT) ) )
            ! get the subobject index for the face in the 2d poloidal plane space
            ifc = curObj%ind(SPACE_POLOIDALPLANE)
            ! copy data
            cpodata(iobj) = b2FaceData( gmap%mapFcix(ifc), gmap%mapFciy(ifc), gmap%mapFcIFace(ifc) )
        else if (present(b2VertexData)) then
            ! Vertex/Node data case
            ! check that it's a vertex
            call assert( all( curObj%cls == CLASS_NODE(1:SPACE_COUNT) ) )
            ! get the subobject index for the face in the 2d poloidal plane space
            ivx = curObj%ind(SPACE_POLOIDALPLANE)
            ! copy data
            cpodata(iobj) = b2VertexData( gmap%mapVxix(ivx), gmap%mapVxiy(ivx) )
        end if

    end do

  end function b2ITMTransformDataB2ToCPOGeneral

#endif

end module b2mod_ual_io_data

