      module b2mod_transport_disruption
      use b2mod_types
      use b2mod_version

      integer, save :: quenching
      real, save :: time_quenching
      real (kind=R8), allocatable, save ::
     *  rho_(:), f_quench(:,:)

      contains

      subroutine alloc_b2mod_transport_disruption(nxd,nyd)
      integer nxd,nyd

      allocate(rho_(-1:nyd))
      allocate(f_quench(-1:nxd,-1:nyd))
      
      end subroutine alloc_b2mod_transport_disruption

      subroutine dealloc_b2mod_transport_disruption

      if(.not.allocated(rho_)) return
      deallocate(rho_)
      deallocate(f_quench)
      
      end subroutine dealloc_b2mod_transport_disruption

      end module b2mod_transport_disruption


