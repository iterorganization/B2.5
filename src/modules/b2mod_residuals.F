      module b2mod_residuals
      use b2mod_types
      real (kind=R8), allocatable, save ::
     *  resco(:,:,:),
     *  reshe(:,:),
     *  reshi(:,:),
     *  reshn(:,:),
     *  resmo(:,:,:),
     *  resmt(:,:),
     *  respo(:,:),
     *  aresco(:),
     *  aresmo(:)
       real (kind=R8), save ::
     *  aresmt,
     *  areshe,
     *  areshi,
     *  arespo
! csc UNCOMMENT WITH AD
!  REAL(kind=r8), ALLOCATABLE, SAVE :: arescod(:), aresmod(:)
!  REAL(kind=r8), SAVE :: aresmtd, areshed, areshid, arespod

      contains

      subroutine alloc_b2mod_residuals(nxd,nyd,nsd)
      implicit none
      integer nxd,nyd,nsd

      allocate(resco(-1:nxd,-1:nyd,0:nsd-1))
      allocate(reshe(-1:nxd,-1:nyd))
      allocate(reshi(-1:nxd,-1:nyd))
      allocate(reshn(-1:nxd,-1:nyd))
      allocate(resmo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(resmt(-1:nxd,-1:nyd))
      allocate(respo(-1:nxd,-1:nyd))
      allocate(aresco(0:nsd-1))
      allocate(aresmo(0:nsd-1))

! csc UNCOMMENT WITH AD
!    ALLOCATE(arescod(0:nsd-1))
!    ALLOCATE(aresmod(0:nsd-1))
!    arescod = 0.0_R8
!    aresmod = 0.0_R8

      resco=0.0_R8
      reshe=0.0_R8
      reshi=0.0_R8
      reshn=0.0_R8
      resmo=0.0_R8
      resmt=0.0_R8
      respo=0.0_R8
      aresco=0.0_R8
      aresmo=0.0_R8

      return
      end subroutine alloc_b2mod_residuals

      subroutine dealloc_b2mod_residuals
      implicit none

      deallocate(resco)
      deallocate(reshe)
      deallocate(reshi)
      deallocate(reshn)
      deallocate(resmo)
      deallocate(resmt)
      deallocate(respo)
      deallocate(aresco)
      deallocate(aresmo)

! csc UNCOMMENT WITH AD
!      deallocate(arescod)
!      deallocate(aresmod)

      return
      end subroutine dealloc_b2mod_residuals

      subroutine read_b2mod_residuals(nread)
      implicit none
      integer nread

      call cfrure (nread, size(resco), resco, 'resco')
      call cfrure (nread, size(reshe), reshe, 'reshe')
      call cfrure (nread, size(reshi), reshi, 'reshi')
      call cfrure (nread, size(resmo), resmo, 'resmo')
      call cfrure (nread, size(resmt), resmt, 'resmt')
      call cfrure (nread, size(respo), respo, 'respo')

      return
      end subroutine read_b2mod_residuals

      subroutine write_b2mod_residuals(nwrite)
      implicit none
      integer nwrite

      call cfwure (nwrite, size(resco), resco, 'resco')
      call cfwure (nwrite, size(reshe), reshe, 'reshe')
      call cfwure (nwrite, size(reshi), reshi, 'reshi')
      call cfwure (nwrite, size(reshn), reshn, 'reshn')
      call cfwure (nwrite, size(resmo), resmo, 'resmo')
      call cfwure (nwrite, size(resmt), resmt, 'resmt')
      call cfwure (nwrite, size(respo), respo, 'respo')

      return
      end subroutine write_b2mod_residuals

      end module b2mod_residuals
