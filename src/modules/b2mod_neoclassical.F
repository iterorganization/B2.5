      module b2mod_neoclassical
      use b2mod_types
      use b2mod_version

      integer nzdecl
c     max. number of ionization stages
      parameter (nzdecl=92)

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
     *  totchi(:,:,:), totkappa(:,:,:)

      contains

      subroutine alloc_b2mod_neoclassical(nxd,nyd,nsd,nzd,mmxd)
      integer nxd,nyd,nsd,nzd,mmxd

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
      allocate(totchi(0:nsd,0:nzd-1,-1:nyd))
      allocate(totkappa(0:nsd,0:nzd-1,-1:nyd))

      end subroutine alloc_b2mod_neoclassical

      end module b2mod_neoclassical


