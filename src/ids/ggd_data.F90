module ggd_data

  use b2mod_types
#ifdef IMAS
  use ids_schemas ! IGNORE
#else
#ifdef ITM
  use euITM_schemas ! IGNORE
#endif
#endif

  implicit none

#ifdef ITM

  interface gridWriteData
      module procedure gridWriteDataScalar, gridWriteDataVector, gridWriteDataMatrix
  end interface

contains

  !> Write a scalar data field given as a scalar data representation to a generic CPO data field. 
  !>
  !> Note: the routine will make sure the required storage is allocated, and will deallocate
  !> and re-allocate fields as necessary.

  subroutine gridWriteDataScalar( cpoField, subgrid, data ) 
    type(type_complexgrid_scalar), intent(inout) :: cpofield
    integer, intent(in) :: subgrid
    real(ITM_R8), intent(in) :: data(:)

    ! set subgrid index
    cpoField % subgrid = subgrid
    
    ! Make sure the data field is properly allocated
    if ( associated(cpoField % scalar) ) then
        if ( .not. all( shape(cpoField % scalar) == shape(data) )) then
            deallocate( cpoField % scalar )
        end if
    end if
    ! If required, allocate storage
    if ( .not. associated(cpoField % scalar) ) then
        allocate(cpoField % scalar( size(data, 1) ))
    end if

    ! copy data
    cpoField % scalar = data

    ! clear unused fields
    if (associated( cpoField % vector )) deallocate(cpoField % vector)
    cpoField % vector => null()
    if (associated( cpoField % matrix )) deallocate(cpoField % matrix)
    cpoField % matrix => null()

  end subroutine gridWriteDataScalar


  !> Write a scalar data field given in a vector data representation to a generic CPO data field. 
  !>
  !> Note: the routine will make sure the required storage is allocated, and will deallocate
  !> and re-allocate fields as necessary.

  subroutine gridWriteDataVector( cpoField, subgrid, data ) 
    type(type_complexgrid_scalar), intent(inout) :: cpofield
    integer, intent(in) :: subgrid
    real(ITM_R8), intent(in) :: data(:,:)

    ! set subgrid index
    cpoField % subgrid = subgrid
    
    ! Make sure the data field is properly allocated
    if ( associated(cpoField % vector) ) then
        if ( .not. all( shape(cpoField % vector) == shape(data) )) then
            deallocate( cpoField % vector )
        end if
    end if
    ! If required, allocate storage
    if ( .not. associated(cpoField % vector) ) then
        allocate(cpoField % vector( size(data, 1), size(data, 2) ))
    end if

    ! copy data
    cpoField % vector = data

    ! clear unused fields
    if (associated( cpoField % scalar )) deallocate(cpoField % scalar)
    cpoField % scalar => null()
    if (associated( cpoField % matrix )) deallocate(cpoField % matrix)
    cpoField % matrix => null()

  end subroutine gridWriteDataVector
  

  !> Write a scalar data field given in a matrix data representation to a generic CPO data field. 
  !>
  !> Note: the routine will make sure the required storage is allocated, and will deallocate
  !> and re-allocate fields as necessary.

  subroutine gridWriteDataMatrix( cpoField, subgrid, data ) 
    type(type_complexgrid_scalar), intent(inout) :: cpofield
    integer, intent(in) :: subgrid
    real(ITM_R8), intent(in) :: data(:,:,:)

    ! set subgrid index
    cpoField % subgrid = subgrid
    
    ! Make sure the data field is properly allocated
    if ( associated(cpoField % matrix) ) then
        if ( .not. all( shape(cpoField % matrix) == shape(data) )) then
            deallocate( cpoField % matrix )
        end if
    end if
    ! If required, allocate storage
    if ( .not. associated(cpoField % matrix) ) then
        allocate(cpoField % matrix( size(data, 1), size(data, 2), size(data, 3) ))
    end if

    ! copy data
    cpoField % matrix = data

    ! clear unused fields
    if (associated( cpoField % scalar )) deallocate(cpoField % scalar)
    cpoField % scalar => null()
    if (associated( cpoField % vector )) deallocate(cpoField % vector)
    cpoField % vector => null()

  end subroutine gridWriteDataMatrix

  !> Write a scalar complex data field given as a scalar data representation to a generic CPO data field. 
  !>
  !> Note: the routine will make sure the required storage is allocated, and will deallocate
  !> and re-allocate fields as necessary.

  subroutine gridWriteDataScalarComplex( cpoField, subgrid, data ) 
    type(type_complexgrid_scalar_cplx), intent(inout) :: cpofield
    integer, intent(in) :: subgrid
    complex(ITM_R8), intent(in) :: data(:)

    ! set subgrid index
    cpoField % subgrid = subgrid
    
    ! Make sure the data field is properly allocated
    if ( associated(cpoField % scalar) ) then
        if ( .not. all( shape(cpoField % scalar) == shape(data) )) then
            deallocate( cpoField % scalar )
        end if
    end if
!    if ( associated(cpoField % scalar ) ) then
!        if ( .not. all( shape(cpoField % scalar) == shape(data) )) then
!            deallocate( cpoField % scalar )
!        end if
!    end if

    ! If required, allocate storage
    if ( .not. associated(cpoField % scalar) ) then
        allocate(cpoField % scalar( size(data, 1) ))
    end if

!    if ( .not. associated(cpoField % scalar % im) ) then
!        allocate(cpoField % scalar % im( size(data, 1) ))
!    end if

    ! copy data
    cpoField % scalar = data
!    cpoField % scalar % im = aimag(data)

    ! clear unused fields
    if (associated( cpoField % vector )) deallocate(cpoField % vector)
    cpoField % vector => null()
    if (associated( cpoField % matrix )) deallocate(cpoField % matrix )
    cpoField % matrix => null()

  end subroutine gridWriteDataScalarComplex

#endif

end module ggd_data

!!!Local Variables:
!!! mode: f90
!!! End:
