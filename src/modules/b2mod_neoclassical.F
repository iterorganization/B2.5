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
     *  theta(:),
     *  bav_(:), b2av_(:), bi2a_(:), rbt_(:),
     *  bgradp_(:), dpsidr_(:), rnq_(:),
     *  fc_(:), gclass_(:), fm_(:,:),
     *  r2i_(:), rho_(:)
      real (kind=R8), allocatable, save ::
     *  dd(:,:,:), vvn(:,:,:),
     *  chi(:,:,:), vwt(:,:,:)

      contains

      subroutine alloc_b2mod_neoclassical(mpg,nsd,mmxd)
      use b2us_map
      implicit none
      type(mapping), intent(in) :: mpg
      integer, intent(in) :: nsd, mmxd
      integer nzd

      nzdecl=nint(maxval(zn(0:nsd-1)))+1
      nzd=nzdecl
      allocate(ncl(0:nsd))
      allocate(isort(0:nsd,0:nzd-1))
      allocate(m(0:nsd))
      allocate(zsp(0:nsd,0:nzd-1))
      allocate(theta(mpg%nCv))
      allocate(bav_(mpg%nFt))
      allocate(b2av_(mpg%nFt))
      allocate(bi2a_(mpg%nFt))
      allocate(rbt_(mpg%nFt))
      allocate(bgradp_(mpg%nFt))
      allocate(dpsidr_(mpg%nFt))
      allocate(rnq_(mpg%nFt))
      allocate(fc_(mpg%nFt))
      allocate(gclass_(mpg%nFt))
      allocate(fm_(1:mmxd,mpg%nFt))
      allocate(r2i_(mpg%nFt))
      allocate(rho_(mpg%nFt))
      allocate(dd(0:nsd,0:nzd-1,mpg%nFt))
      allocate(vvn(0:nsd,0:nzd-1,mpg%nFt))
      allocate(chi(0:nsd,0:nzd-1,mpg%nFt))
      allocate(vwt(0:nsd,0:nzd-1,mpg%nFt))

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

      return
      end subroutine dealloc_b2mod_neoclassical

      end module b2mod_neoclassical


