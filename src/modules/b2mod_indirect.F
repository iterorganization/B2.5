      module b2mod_indirect
      use b2mod_types
      implicit none

      integer, parameter :: NNREGMAX=8

      integer, allocatable, save ::
     *  leftix(:,:),
     *  leftiy(:,:),
     *  rightix(:,:),
     *  rightiy(:,:),
     *  topix(:,:),
     *  topiy(:,:),
     *  bottomix(:,:),
     *  bottomiy(:,:),
     *  region(:,:,:),
     *  resignore(:,:,:),
     *  leftcut(:),
     *  rightcut(:),
     *  topcut(:),
     *  bottomcut(:)
      integer, save ::
     *  nncut, nnreg(0:2), periodic_bc, nncutmax

      contains

      subroutine alloc_b2mod_indirect(nxd,nyd,nncutmax)
      implicit none
      integer nxd,nyd,nncutmax

      allocate(leftix(-1:nxd,-1:nyd))
      allocate(leftiy(-1:nxd,-1:nyd))
      allocate(rightix(-1:nxd,-1:nyd))
      allocate(rightiy(-1:nxd,-1:nyd))
      allocate(topix(-1:nxd,-1:nyd))
      allocate(topiy(-1:nxd,-1:nyd))
      allocate(bottomix(-1:nxd,-1:nyd))
      allocate(bottomiy(-1:nxd,-1:nyd))
      allocate(region(-1:nxd,-1:nyd,0:2))
      nnreg=0
      allocate(resignore(-1:nxd,-1:nyd,2))
      nncutmax=10
      allocate(leftcut(1:nncutmax))
      allocate(rightcut(1:nncutmax))
      allocate(topcut(1:nncutmax))
      allocate(bottomcut(1:nncutmax))
      bottomcut=nyd+1
      topcut   =-2
      rightcut =-2
      leftcut  =nxd+1

      return
      end subroutine alloc_b2mod_indirect

      subroutine dealloc_b2mod_indirect
      implicit none

      deallocate(leftix)
      deallocate(leftiy)
      deallocate(rightix)
      deallocate(rightiy)
      deallocate(topix)
      deallocate(topiy)
      deallocate(bottomix)
      deallocate(bottomiy)
      deallocate(region)
      deallocate(resignore)
      deallocate(leftcut)
      deallocate(rightcut)
      deallocate(topcut)
      deallocate(bottomcut)

      return
      end subroutine dealloc_b2mod_indirect

      end module b2mod_indirect
