      module b2mod_transport_nspecies
      use b2mod_types
      use b2mod_version
      implicit none
      real (kind=R8), allocatable, save ::
     *  chci_al(:,:,:,:),
     *  cvsahz_hAdp_albe(:,:,:,:,:),
     *  cvsahz_hBdp_al(:,:,:,:)                                       !som 08.10.21

      contains

      subroutine alloc_b2mod_transport_nspecies(nxd,nyd,nspecies)
      implicit none
      integer nxd,nyd,nspecies

      print *, 'nspecies:', nspecies
      allocate(chci_al(-1:nxd,-1:nyd,0:1,0:nspecies-1))
      allocate(cvsahz_hAdp_albe(-1:nxd,-1:nyd,0:1,
     & 0:nspecies-1,0:nspecies-1))
      allocate(cvsahz_hBdp_al(-1:nxd,-1:nyd,0:1,0:nspecies-1))        !som 08.10.21

      return
      end subroutine alloc_b2mod_transport_nspecies

      subroutine dealloc_b2mod_transport_nspecies
      implicit none

      if (.not.allocated(chci_al)) return
      deallocate(chci_al)
      deallocate(cvsahz_hAdp_albe)
      deallocate(cvsahz_hBdp_al)                                 !som 08.10.21

      return
      end subroutine dealloc_b2mod_transport_nspecies

      end module b2mod_transport_nspecies
