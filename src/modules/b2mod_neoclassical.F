      module b2mod_neoclassical
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_version

c     max. number of ionisation stages
      integer, save :: nzdecl

      logical, save :: newneo
      integer, save :: nsl
      integer, allocatable, save ::
     *  ncl(:), isort(:,:)
      real (kind=R8), allocatable, save ::
     *  m(:), zsp(:,:),
     *  theta(:,:),
     *  bav_(:), b2av_(:), bi2a_(:), rbt_(:),
     *  bgradp_(:), dpsidr_(:), rnq_(:),
     *  fc_(:), gclass_(:), fm_(:,:),
     *  r2i_(:), rho_(:)
      real (kind=R8), allocatable, save ::
     *  dd(:,:,:), vvn(:,:,:),
     *  chi(:,:,:), vwt(:,:,:)

      real (kind=R8), save ::
     *   bav, b2av, bi2a, rbt, bgradp, dpsidr, rnq, fc, gclass,
     &   r2i
      real (kind=R8), allocatable, save ::
     &   fm(:), cosbl(:), cosb2(:), sinbl(:), sinb2(:)

      contains

      subroutine alloc_b2mod_neoclassical(nxd,nyd,nsd,mmxd)
      integer nxd,nyd,nsd,nzd,mmxd

      nzdecl=nint(maxval(zn(0:nsd-1)))+1
      nzd=nzdecl
      allocate(ncl(0:nsd))
      allocate(isort(0:nsd,0:nzd-1))
      allocate(m(0:nsd))
      allocate(zsp(0:nsd,0:nzd-1))
      allocate(theta(-1:nxd,-1:nyd))
      allocate(bav_(-1:nyd))
      allocate(b2av_(-1:nyd))
      allocate(bi2a_(-1:nyd))
      allocate(rbt_(-1:nyd))
      allocate(bgradp_(-1:nyd))
      allocate(dpsidr_(-1:nyd))
      allocate(rnq_(-1:nyd))
      allocate(fc_(-1:nyd))
      allocate(gclass_(-1:nyd))
      allocate(fm_(1:mmxd,-1:nyd))
      allocate(r2i_(-1:nyd))
      allocate(rho_(-1:nyd))
      allocate(dd(0:nsd,0:nzd-1,-1:nyd))
      allocate(vvn(0:nsd,0:nzd-1,-1:nyd))
      allocate(chi(0:nsd,0:nzd-1,-1:nyd))
      allocate(vwt(0:nsd,0:nzd-1,-1:nyd))

      allocate(fm(0:nmaxgr-1))
      allocate(cosbl(0:nmaxgr-1))
      allocate(cosb2(0:nmaxgr-1))
      allocate(sinbl(0:nmaxgr-1))
      allocate(sinb2(0:nmaxgr-1))

      end subroutine alloc_b2mod_neoclassical

      subroutine dealloc_b2mod_neoclassical

      if (.not.allocated(ncl)) return
      deallocate(ncl)
      deallocate(isort)
      deallocate(m)
      deallocate(zsp)
      deallocate(theta)
      deallocate(bav_)
      deallocate(b2av_)
      deallocate(bi2a_)
      deallocate(rbt_)
      deallocate(bgradp_)
      deallocate(dpsidr_)
      deallocate(rnq_)
      deallocate(fc_)
      deallocate(gclass_)
      deallocate(fm_)
      deallocate(r2i_)
      deallocate(rho_)
      deallocate(dd)
      deallocate(vvn)
      deallocate(chi)
      deallocate(vwt)

      deallocate(fm)
      deallocate(cosbl)
      deallocate(cosb2)
      deallocate(sinbl)
      deallocate(sinb2)

      return
      end subroutine dealloc_b2mod_neoclassical

      end module b2mod_neoclassical


