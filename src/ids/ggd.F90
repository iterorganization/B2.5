module ggd

  use ggd_assert
  use ggd_common
#ifdef IMAS
  use ids_schemas ! IGNORE
#else
#ifdef ITM
  use euitm_schemas ! IGNORE
  use ggd_access
  use ggd_structured
  use ggd_object
  use ggd_objectlist
  use ggd_subgrid
#endif
#endif
  use string

  implicit none      

#ifdef ITM

contains

  !> Check a grid description for possible errors

  subroutine gridSanityCheck( grid, doStop )
    type(type_complexgrid), intent(in) :: grid
    logical, intent(in), optional :: doStop
    
    ! integer
    integer :: is

    ! disable immediate stopping in assertion module
    call assertSetStopMode( doStop = .false. )   

    do is = 1, size( grid % spaces )
       call assertSetMsgPrefix( "Space "// Int2str( is ) )
       call gridSanityCheckSpace( grid % spaces( is ) ) 
    end do

    ! Output error message if any errors found, and stop (if requested)
    call assertStopOnFailed( "Grid contained errors", doStop )


    ! reset assertion module
    call assertReset()

    ! re-enable immediate stop for assertions
    call assertSetStopMode()   
    
    ! clear message prefix
    call assertSetMsgPrefix()

  end subroutine gridSanityCheck
  

  !> Check a grid space description for possible errors

  subroutine gridSanityCheckSpace( space )
    type(type_complexgrid_space), intent(in) :: space
    
    call assert( ( size( space % coordtype ) <= size( space % objects ) ), &
         & "coordtype, objects,have inconsistent size" )
        

  end subroutine gridSanityCheckSpace

#endif

end module ggd

!!!Local Variables:
!!! mode: f90
!!! End:
