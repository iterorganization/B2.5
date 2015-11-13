      module b2mod_geo2
      use b2mod_types

      real (kind=R8), allocatable, save ::
     *  bb1(:,:,:),
     *  crx1(:,:,:),
     *  cry1(:,:,:),
     *  ffbz1(:,:,:),
     *  fpsi1(:,:,:),
     *  gs1(:,:,:),
     *  hx1(:,:),
     *  hy1(:,:),
     *  qz1(:,:,:),
     *  qc1(:,:),
     *  pbs1(:,:,:),
     *  wbbl1(:,:,:),
     *  wbbr1(:,:,:),
     *  wbbv1(:,:,:),
     *  vol1(:,:)

      integer, allocatable, save ::
     *  leftix1(:,:),
     *  leftiy1(:,:),
     *  rightix1(:,:),
     *  rightiy1(:,:),
     *  topix1(:,:),
     *  topiy1(:,:),
     *  bottomix1(:,:),
     *  bottomiy1(:,:),
     *  leftcut1(:),
     *  rightcut1(:),
     *  bottomcut1(:),
     *  topcut1(:),
     *  region1(:,:,:),
     *  resignore1(:,:,:)

      contains

      subroutine alloc_b2mod_geo2(nxd,nyd,nncutmax)
      integer nxd, nyd, nncutmax
      allocate(bb1(-1:nxd,-1:nyd,0:3))
      allocate(crx1(-1:nxd,-1:nyd,0:3))
      allocate(cry1(-1:nxd,-1:nyd,0:3))
      allocate(ffbz1(-1:nxd,-1:nyd,0:3))
      allocate(fpsi1(-1:nxd,-1:nyd,0:3))
      allocate(gs1(-1:nxd,-1:nyd,0:2))
      allocate(hx1(-1:nxd,-1:nyd))
      allocate(hy1(-1:nxd,-1:nyd))
      allocate(qz1(-1:nxd,-1:nyd,0:1))
      allocate(qc1(-1:nxd,-1:nyd))
      allocate(pbs1(-1:nxd,-1:nyd,0:1))
      allocate(wbbl1(-1:nxd,-1:nyd,0:3))
      allocate(wbbr1(-1:nxd,-1:nyd,0:3))
      allocate(wbbv1(-1:nxd,-1:nyd,0:3))
      allocate(vol1(-1:nxd,-1:nyd))
      allocate(leftix1(-1:nxd,-1:nyd))
      allocate(leftiy1(-1:nxd,-1:nyd))
      allocate(rightix1(-1:nxd,-1:nyd))
      allocate(rightiy1(-1:nxd,-1:nyd))
      allocate(topix1(-1:nxd,-1:nyd))
      allocate(topiy1(-1:nxd,-1:nyd))
      allocate(bottomix1(-1:nxd,-1:nyd))
      allocate(bottomiy1(-1:nxd,-1:nyd))
      allocate(leftcut1(1:nncutmax))
      allocate(rightcut1(1:nncutmax))
      allocate(topcut1(1:nncutmax))
      allocate(bottomcut1(1:nncutmax))
      allocate(region1(-1:nxd,-1:nyd,0:2))
      allocate(resignore1(-1:nxd,-1:nyd,1:2))
      end subroutine alloc_b2mod_geo2

      subroutine dealloc_b2mod_geo2
      deallocate(bb1)
      deallocate(crx1)
      deallocate(cry1)
      deallocate(ffbz1)
      deallocate(fpsi1)
      deallocate(gs1)
      deallocate(hx1)
      deallocate(hy1)
      deallocate(qz1)
      deallocate(qc1)
      deallocate(pbs1)
      deallocate(wbbl1)
      deallocate(wbbr1)
      deallocate(wbbv1)
      deallocate(vol1)
      deallocate(leftix1)
      deallocate(leftiy1)
      deallocate(rightix1)
      deallocate(rightiy1)
      deallocate(topix1)
      deallocate(topiy1)
      deallocate(bottomix1)
      deallocate(bottomiy1)
      deallocate(leftcut1)
      deallocate(rightcut1)
      deallocate(topcut1)
      deallocate(bottomcut1)
      deallocate(region1)
      deallocate(resignore1)
      return
      end subroutine dealloc_b2mod_geo2

      end module b2mod_geo2
