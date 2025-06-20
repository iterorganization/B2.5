      module b2mod_transport_disruption
      use b2mod_types
      use b2mod_version

      integer, save :: quenching
      real (kind=R8), save :: time_quenching
      real (kind=R8), allocatable, save ::
     *  rho_(:), f_quench(:)

      contains

      subroutine alloc_b2mod_transport_disruption(mpg)
      use b2us_map
      implicit none
      type(mapping), intent(in) :: mpg

      allocate(rho_(mpg%nFs))
      allocate(f_quench(mpg%nCv))
      
      end subroutine alloc_b2mod_transport_disruption

      subroutine dealloc_b2mod_transport_disruption

      if(.not.allocated(rho_)) return
      deallocate(rho_)
      deallocate(f_quench)
      
      end subroutine dealloc_b2mod_transport_disruption

      end module b2mod_transport_disruption
