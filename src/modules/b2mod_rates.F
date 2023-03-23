      module b2mod_rates
      use b2mod_types
      integer nscxmax
      parameter (nscxmax=3)
      real (kind=R8), allocatable, save ::
     *  rlcx(:,:,:,:,:),
     *  rlqa(:,:,:,:),
     *  rlqr(:,:,:,:),
     *  rlrd(:,:,:,:),
     *  rlbr(:,:,:,:),
     *  rlra(:,:,:,:),
     *  rlsa(:,:,:,:),
     *  rlza(:,:,:,:),
     *  rlz2(:,:,:,:),
     *  rlpt(:,:,:,:),
     *  rlpi(:,:,:,:)
      real (kind=R8), allocatable, save ::
     *  rza(:,:,:),
     *  rz2(:,:,:),
     *  rpt(:,:,:),
     *  rpi(:,:,:)

      contains

      subroutine alloc_b2mod_rates(nxd,nyd,nsd)
      implicit none
      integer nxd,nyd,nsd

      allocate(rlcx(-1:nxd,-1:nyd,0:1,0:nsd-1,0:nscxmax-1))
      allocate(rlqa(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(rlqr(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(rlrd(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(rlbr(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(rlra(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(rlsa(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(rlza(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(rlz2(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(rlpt(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(rlpi(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(rza(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rz2(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rpt(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rpi(-1:nxd,-1:nyd,0:nsd-1))

      return
      end subroutine alloc_b2mod_rates

      subroutine dealloc_b2mod_rates
      implicit none

      deallocate(rlcx)
      deallocate(rlqa)
      deallocate(rlqr)
      deallocate(rlrd)
      deallocate(rlbr)
      deallocate(rlra)
      deallocate(rlsa)
      deallocate(rlza)
      deallocate(rlz2)
      deallocate(rlpt)
      deallocate(rlpi)
      deallocate(rza)
      deallocate(rz2)
      deallocate(rpt)
      deallocate(rpi)

      return
      end subroutine dealloc_b2mod_rates

      end module b2mod_rates
