      module b2mod_transport_nspecies
      use b2mod_types
      use b2mod_version
      implicit none
      real (kind=R8), allocatable, save ::
     *  chci_al_ast(:,:,:,:,:),
     *  cvsa_hAdp_albe(:,:,:,:,:,:),
     *  cvsa_hBdp_al(:,:,:,:,:),                                      !som 20.10.21
     *  cvsa_uAdp_albe(:,:,:,:,:,:),
     *  cvsa_uBdp_al(:,:,:,:,:),                                      !som 26.04.24
     *  cvsahz_hAdp_albe(:,:,:,:,:,:),
     *  cvsahz_hBdp_al(:,:,:,:,:),                                    !som 08.10.21
     *  cvsahz_uAdp_albe(:,:,:,:,:,:),
     *  cvsahz_uBdp_al(:,:,:,:,:)                                     !som 22.12.23

      contains

      subroutine alloc_b2mod_transport_nspecies(nxd,nyd,nspecies)
      implicit none
      integer nxd,nyd,nspecies

!      print *, 'nspecies:', nspecies
      allocate(chci_al_ast(-1:nxd,-1:nyd,0:1,0:1,0:nspecies-1))
      allocate(cvsa_hAdp_albe(-1:nxd,-1:nyd,0:1,0:1,
     &                         0:nspecies-1,0:nspecies-1))
      allocate(cvsa_hBdp_al(-1:nxd,-1:nyd,0:1,0:1,0:nspecies-1))      !som 20.10.21
      allocate(cvsa_uAdp_albe(-1:nxd,-1:nyd,0:1,0:1,
     &                         0:nspecies-1,0:nspecies-1))
      allocate(cvsa_uBdp_al(-1:nxd,-1:nyd,0:1,0:1,0:nspecies-1))      !som 26.04.24
      allocate(cvsahz_hAdp_albe(-1:nxd,-1:nyd,0:1,0:1,
     &                           0:nspecies-1,0:nspecies-1))
      allocate(cvsahz_hBdp_al(-1:nxd,-1:nyd,0:1,0:1,0:nspecies-1))    !som 08.10.21
      allocate(cvsahz_uAdp_albe(-1:nxd,-1:nyd,0:1,0:1,
     &                           0:nspecies-1,0:nspecies-1))
      allocate(cvsahz_uBdp_al(-1:nxd,-1:nyd,0:1,0:1,0:nspecies-1))    !som 22.12.23

      return
      end subroutine alloc_b2mod_transport_nspecies

      subroutine dealloc_b2mod_transport_nspecies
      implicit none

      if (.not.allocated(chci_al_ast)) return
      deallocate(chci_al_ast)
      deallocate(cvsa_hAdp_albe)
      deallocate(cvsa_hBdp_al)                                        !som 08.10.21
      deallocate(cvsahz_hAdp_albe)
      deallocate(cvsahz_hBdp_al)                                      !som 08.10.21
      deallocate(cvsahz_uAdp_albe)
      deallocate(cvsahz_uBdp_al)                                      !som 22.12.23

      return
      end subroutine dealloc_b2mod_transport_nspecies

      end module b2mod_transport_nspecies
