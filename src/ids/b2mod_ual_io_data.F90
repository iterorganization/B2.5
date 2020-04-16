!!-----------------------------------------------------------------------------
!! DOCUMENTATION:
!>      @section b2uw_ualio_data_desc Description
!!      Module providing routines to transform variables stored in the B2
!!      data structure into the form expected by CPO/IDS data structure.
!!
!!      @subsection b2uw_ualio_data_syx     Exceptional syntax explanation
!!      @code
!!          ! IGNORE    !! syntax used to ignore this module in list
!!                      !! dependency when compiling the code
!!      @endcode
!!
!!-----------------------------------------------------------------------------
module b2mod_ual_io_data

    use b2mod_types , B2_R8 => R8, B2_R4 => R4
    use helper
    use logging , only: logmsg, LOGDEBUG
    use b2mod_connectivity , REMOVED_B2_R8 => R8
    use carre_constants
    use b2mod_cellhelper

    use b2mod_grid_mapping
#ifdef IMAS
#if IMAS_MINOR_VERSION > 11
#if IMAS_MINOR_VERSION > 14
    use b2mod_ual_io_grid &
     & , only : ids_generic_grid_aos3_root
#else
    use b2mod_ual_io_grid &
     & , only : ids_generic_grid_dynamic
#endif
    use b2mod_ual_io_grid &
     & , only : IDS_real, GridObject, &
     &          getGridSubsetObject, SPACE_POLOIDALPLANE, &
     &          GridWriteData, IDS_CLASS_CELL, IDS_CLASS_NODE, &
     &          IDS_CLASS_POLOIDALRADIAL_FACE
    use ids_grid_subgrid  & ! IGNORE
     & , only : getGridSubsetSize

    implicit none

    !> Provides service routines to transform data from B2 to IMAS IDS
    !! (data in form of vertex, face or cell)
    interface b2_IMAS_Transform_Data_B2_To_IDS
        module procedure b2_IMAS_Transform_Data_B2_To_IDS_Cell, &
            &   b2_IMAS_Transform_Data_B2_To_IDS_Face
    end interface

contains

    !! Below here: service routines to transform data from B2 to IMAS IDS

    !> Transform data from B2 to IDS cell
    function b2_IMAS_Transform_Data_B2_To_IDS_Cell( grid, gridSubsetInd, gmap,  &
            &   b2CellData ) result( idsdata )
#if IMAS_MINOR_VERSION > 14
        type(ids_generic_grid_aos3_root), intent(in) :: grid !< Type of IDS
            !< data structure, designed for handling grid geometry data
#else
        type(ids_generic_grid_dynamic), intent(in) :: grid !< Type of IDS
            !< data structure, designed for handling grid geometry data
#endif
        integer, intent(in) :: gridSubsetInd !< Index of the
            !< grid subset the data is to be stored on
        type(B2GridMap), intent(in) :: gmap !< The grid mapping as computed
            !< by b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS
        real(IDS_real), intent(in) :: b2CellData( -1:gmap%b2nx, -1:gmap%b2ny )
        real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
                !< handing data field values

        idsdata => b2_IMAS_Transform_Data_B2_To_IDS_General( grid,  &
            &   gridSubsetInd, gmap, b2CellData = b2CellData )
    end function b2_IMAS_Transform_Data_B2_To_IDS_Cell

    !> Transform data from B2 to IDS face
    function b2_IMAS_Transform_Data_B2_To_IDS_Face( grid, gridSubsetInd, gmap,  &
            &   b2FaceData ) result( idsdata )
#if IMAS_MINOR_VERSION > 14
        type(ids_generic_grid_aos3_root), intent(in)  :: grid !< Type of IDS
            !< data structure, designed for handling grid geometry data
#else
        type(ids_generic_grid_dynamic), intent(in) :: grid !< Type of IDS
            !< data structure, designed for handling grid geometry data
#endif
        integer, intent(in) :: gridSubsetInd !< Index of the
            !< grid subset the data is to be stored on
        type(B2GridMap), intent(in) :: gmap !< The grid mapping as computed
            !< by b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS
        real(IDS_real), intent(in) ::  &
            &   b2FaceData( -1:gmap%b2nx, -1:gmap%b2ny, 0:1 )   !< Face data
            !< given on the 2D B2 data structure
        real(IDS_real), dimension(:), pointer :: idsdata !< Array for
            !< handing data field values

        idsdata => b2_IMAS_Transform_Data_B2_To_IDS_General( grid,  &
            &   gridSubsetInd, gmap, b2FaceData = b2FaceData )
    end function b2_IMAS_Transform_Data_B2_To_IDS_Face

    !! TODO: find a way to include this subroutine in the
    !! b2_IMAS_Transform_Data_B2_To_IDS interface

    !> Transform data from B2 to IDS vertex
    function b2_IMAS_Transform_Data_B2_To_IDS_Vertex( grid, gridSubsetInd,   &
            &   gmap, b2VertexData ) result( idsdata )
#if IMAS_MINOR_VERSION > 14
        type(ids_generic_grid_aos3_root), intent(in)  :: grid !< Type of IDS
            !< data structure, designed for handling grid geometry data
#else
        type(ids_generic_grid_dynamic), intent(in) :: grid !< Type of IDS
            !< data structure, designed for handling grid geometry data
#endif
        integer, intent(in)         :: gridSubsetInd !< Index of the
            !< grid subset the data is to be stored on
        type(B2GridMap), intent(in) :: gmap !< The grid mapping as computed by
            !< b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS
        real(IDS_real), intent(in)  :: b2VertexData( -1:gmap%b2nx, -1:gmap%b2ny )
            !< Array holding vertex coordinates (2D space)
        real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
            !< handing data field values

        idsdata => b2_IMAS_Transform_Data_B2_To_IDS_General( grid, gridSubsetInd,   &
            &   gmap, b2VertexData = b2VertexData )
    end function b2_IMAS_Transform_Data_B2_To_IDS_Vertex

    !> Transform a quantity stored on faces from a 2d B2 array into a 1d IMAS
    !! IDS array for a given grid subset id. Either b2CellData or b2FaceData
    !! must be given. Do not use this directly, use the provided general
    !! interface b2_IMAS_Transform_Data_B2_To_IDS instead.
    function b2_IMAS_Transform_Data_B2_To_IDS_General( grid, gridSubsetInd,  &
            &   gmap, b2CellData, b2FaceData, b2VertexData ) result( idsdata )
#if IMAS_MINOR_VERSION > 14
        type(ids_generic_grid_aos3_root), intent(in)  :: grid !< Type of IDS
            !< data structure, designed for handling grid geometry data
#else
        type(ids_generic_grid_dynamic), intent(in) :: grid !< Type of IDS
            !< data structure, designed for handling grid geometry data
#endif
        integer, intent(in) :: gridSubsetInd !< Base grid subset index
        type(B2GridMap), intent(in) :: gmap !< The grid mapping as computed
            !< by b2CreateMap holding an intermediate grid description to be
            !< transferred into a CPO or IDS
        real(IDS_real), intent(in), optional :: &
            &   b2CellData( -1:gmap%b2nx, -1:gmap%b2ny )      !< Cell data given
            !< on the 2D B2 data structure
        real(IDS_real), intent(in), optional :: &
            &   b2FaceData( -1:gmap%b2nx, -1:gmap%b2ny, 0:1 ) !< Face data
            !< given on the 2D B2 data structure
        real(IDS_real), intent(in), optional :: &
            &   b2VertexData( -1:gmap%b2nx, -1:gmap%b2ny )    !< Vertex data
            !< given on the 2D B2 data structure
        real(IDS_real), dimension(:), pointer :: idsdata      !< Array for
            !< handing data field values

        !! Internal variables
        integer :: nobjs    !< Total number of objects
        integer :: iobj     !< Object index (iterator)
        integer :: ifc      !< Face index
        integer :: icv      !< Cell index
        integer :: ivx      !< vertex index
        type(GridObject) :: curObj

        !! .neqv. is xor (exclusive or)
        !! (http://de.wikibooks.org/wiki/Fortran:_Fortran_95:_Logische_Ausdr%C3%BCcke)
        !! TODO: FIX
        !call assert( present(b2CellData) .neqv. present(b2FaceData) )
        !! Allocate result vector according to grid subset size
        nobjs = getGridSubsetSize( grid%grid_subset( gridSubsetInd ) )
        allocate( idsdata( nobjs ) )
        !! Collect all data items for the grid subset objects
        do iobj = 1, nobjs
            !! Get the object descriptor
            curObj = getGridSubsetObject( grid%grid_subset( gridSubsetInd ), iobj )
            if( present( b2CellData ) ) then
                !! Cell data case
                !! check that it is a cell
                call xertst( all( curObj%cls == IDS_CLASS_CELL ),   &
                    &   "Assert error 1 in b2_IMAS_Transform_Data_B2_To_IDS_General!" )
                !! get the subobject index for the face in the 2d poloidal
                !! plane space
                icv = curObj%ind(SPACE_POLOIDALPLANE)
                if( icv .eq. 0 ) then
                   idsdata( iobj ) = 0.0 !! TODO DP skip when no index present
                else
                                !! copy data
                   idsdata( iobj ) =   &
                      & b2CellData( gmap%mapCvix( icv ), gmap%mapCviy( icv ) )
                end if
            else if( present( b2FaceData ) ) then
                !! Face data case
                !! check that it is a face
                call xertst( all( curObj%cls ==             &
                    &   IDS_CLASS_POLOIDALRADIAL_FACE ),    &
                    &   "Assert error 2 in b2_IMAS_Transform_Data_B2_To_IDS_General!" )
                !! get the subobject index for the face in the 2d poloidal
                !! plane space
                ifc = curObj%ind( SPACE_POLOIDALPLANE )
                !! copy data
                idsdata( iobj ) = b2FaceData( gmap%mapFcix( ifc ),  &
                    &   gmap%mapFciy( ifc ), gmap%mapFcIFace( ifc ) )
            else if( present( b2VertexData )) then
                !! Vertex/Node data case
                !! check that it is a vertex
                call xertst( all( curObj%cls == IDS_CLASS_NODE ),   &
                    &   "Assert error 3 in b2_IMAS_Transform_Data_B2_To_IDS_General!" )
                !! get the subobject index for the face in the 2d poloidal
                !! plane space
                ivx = curObj%ind( SPACE_POLOIDALPLANE )
                !! copy data
                idsdata( iobj ) =   &
                    &   b2VertexData( gmap%mapVxix( ivx ), gmap%mapVxiy( ivx ) )
            end if
        end do
    end function b2_IMAS_Transform_Data_B2_To_IDS_General

#endif
#else
#ifdef ITM_ENVIRONMENT_LOADED

    use itm_types , ITM_R8 => R8, ITM_R4 => R4 ! IGNORE
    use euITM_schemas ! IGNORE
    use itm_constants , pi => itm_pi ! IGNORE

  implicit none

  interface b2ITMTransformDataB2ToCPO
      module procedure b2ITMTransformDataB2ToCPOCell, b2ITMTransformDataB2ToCPOFace
  end interface

contains


  !! Below here: service routines to transform data from B2 to CPO

  function b2ITMTransformDataB2ToCPOCell( grid, subgridId, gmap, b2CellData ) result( cpodata )
    real(ITM_R8), dimension(:), pointer :: cpodata

    type(type_complexgrid), intent(in) :: grid
    integer, intent(in) :: subgridId
    type(B2GridMap), intent(in) :: gmap
    real(ITM_R8), intent(in) :: b2CellData(-1:gmap%b2nx, -1:gmap%b2ny)

    cpodata => b2ITMTransformDataB2ToCPOGeneral( grid, subgridId, gmap, b2CellData = b2CellData )
  end function b2ITMTransformDataB2ToCPOCell

  function b2ITMTransformDataB2ToCPOFace( grid, subgridId, gmap, b2FaceData ) result( cpodata )
    real(ITM_R8), dimension(:), pointer :: cpodata

    type(type_complexgrid), intent(in) :: grid
    integer, intent(in) :: subgridId
    type(B2GridMap), intent(in) :: gmap
    real(ITM_R8), intent(in) :: b2FaceData(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)

    cpodata => b2ITMTransformDataB2ToCPOGeneral( grid, subgridId, gmap, b2FaceData = b2FaceData )
  end function b2ITMTransformDataB2ToCPOFace

!! TODO: find a way to include this subroutine in the b2ITMTransformDataB2ToCPO interface

  function b2ITMTransformDataB2ToCPOVertex( grid, subgridId, gmap, b2VertexData ) result( cpodata )
    real(ITM_R8), dimension(:), pointer :: cpodata

    type(type_complexgrid), intent(in) :: grid
    integer, intent(in) :: subgridId
    type(B2GridMap), intent(in) :: gmap
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
    type(B2GridMap), intent(in) :: gmap
    real(ITM_R8), intent(in), optional :: b2CellData(-1:gmap%b2nx, -1:gmap%b2ny)
    real(ITM_R8), intent(in), optional :: b2FaceData(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)
    real(ITM_R8), intent(in), optional :: b2VertexData(-1:gmap%b2nx, -1:gmap%b2ny)

    !! internal
    integer :: nobjs, iobj, ifc, icv, ivx
    type(GridObject) :: curObj

    !! .neqv. is xor (exclusive or)
    !! (http://de.wikibooks.org/wiki/Fortran:_Fortran_95:_Logische_Ausdr%C3%BCcke)
    !! TODO: FIX
    !call assert( present(b2CellData) .neqv. present(b2FaceData) )

    !! Allocate result vector according to subgrid size
    nobjs = gridSubGridSize( grid%subgrids(subgridId) )
    allocate( cpodata(nobjs) )

    !! Collect all data items for the subgrid objects
    do iobj = 1, nobjs
        !! Get the object descriptor
        curObj = subGridGetObject( grid%subgrids(subgridId), iobj )

        if (present(b2CellData)) then
            !! Cell data case
            !! check that it is a cell
            call xertst( all( curObj%cls == CLASS_CELL(1:SPACE_COUNT) ), &
                "Assert error 1 (cell test) in b2ITMTransformDataB2ToCPOGeneral" )
            !! get the subobject index for the face in the 2d poloidal plane space
            icv = curObj%ind(SPACE_POLOIDALPLANE)
            !! copy data
            cpodata(iobj) = b2CellData( gmap%mapCvix(icv), gmap%mapCviy(icv) )
        else if (present(b2FaceData)) then
            !! Face data case
            !! check that it is a face
            call xertst( all( curObj%cls == CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT) ), &
                "Assert error 2 (face test) in b2ITMTransformDataB2ToCPOGeneral" )
            !! get the subobject index for the face in the 2d poloidal plane space
            ifc = curObj%ind(SPACE_POLOIDALPLANE)
            !! copy data
            cpodata(iobj) = b2FaceData( gmap%mapFcix(ifc), gmap%mapFciy(ifc), gmap%mapFcIFace(ifc) )
        else if (present(b2VertexData)) then
            !! Vertex/Node data case
            !! check that it is a vertex
            call xertst( all( curObj%cls == CLASS_NODE(1:SPACE_COUNT) ), &
                "Assert error 3 (vertex test) in b2ITMTransformDataB2ToCPOGeneral" )
            !! get the subobject index for the face in the 2d poloidal plane space
            ivx = curObj%ind(SPACE_POLOIDALPLANE)
            !! copy data
            cpodata(iobj) = b2VertexData( gmap%mapVxix(ivx), gmap%mapVxiy(ivx) )
        end if

    end do

  end function b2ITMTransformDataB2ToCPOGeneral

#endif
#endif

end module b2mod_ual_io_data

!!!Local Variables:
!!! mode: f90
!!! End:
