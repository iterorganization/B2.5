      module b2mod_nclass
      use b2mod_types
      use b2mod_version

      integer, save :: nsl
      integer, allocatable, save ::
     *  ncl(:), isort(:,:)
      real (kind=R8), allocatable, save ::
     *  m(:), zsp(:,:),
     *  theta(:,:),
     *  p_b2_(:), p_bm2_(:), p_eb_(:), p_fhat_(:),
     *  p_ft_(:), p_grbm2_(:), p_grphi_(:),
     *  p_gr2phi_(:), p_ngrth_(:), p_fm_(:,:), rho_(:)
      real (kind=R8), allocatable, save ::
     *  dd(:,:,:), vvn(:,:,:),
     *  totchi(:,:,:), totkappa(:,:,:)

      contains

      subroutine alloc_b2mod_nclass(nxd,nyd,nsd,nzd)
      integer nxd,nyd,nsd,nzd

      allocate(ncl(0:nsd))
      allocate(isort(0:nsd,0:nzd-1))
      allocate(m(0:nsd))
      allocate(zsp(0:nsd,0:nzd-1))
      allocate(theta(-1:nxd,-1:nyd))
      allocate(p_b2_(-1:nyd))
      allocate(p_bm2_(-1:nyd))
      allocate(p_eb_(-1:nyd))
      allocate(p_fhat_(-1:nyd))
      allocate(p_ft_(-1:nyd))
      allocate(p_grbm2_(-1:nyd))
      allocate(p_grphi_(-1:nyd))
      allocate(p_gr2phi_(-1:nyd))
      allocate(p_ngrth_(-1:nyd))
      allocate(p_fm_(1:mmxd,3))
      allocate(rho_(-1:nyd))
      allocate(dd(0:nsd,0:nzd-1,-1:nyd))
      allocate(vvn(0:nsd,0:nzd-1,-1:nyd))
      allocate(totchi(0:nsd,0:nzd-1,-1:nyd))
      allocate(totkappa(0:nsd,0:nzd-1,-1:nyd))
      


      end subroutine alloc_b2mod_nclass

      end module b2mod_nclass


