      module b2us_work
      use b2mod_types

      real (kind=R8), allocatable, save ::
     *  delne(:), delni(:,:), delna(:,:),  delpo(:),
     *  delte(:),  delti(:),  delua(:,:)

      contains

      subroutine alloc_b2us_work(nCv,nsd)
      integer nCv,nsd

      allocate(delne(1:nCv))
      allocate(delni(1:nCv,0:1))
      allocate(delna(1:nCv,0:nsd-1))
      allocate(delpo(1:nCv))
      allocate(delte(1:nCv))
      allocate(delti(1:nCv))
      allocate(delua(1:nCv,0:nsd-1))

      delne=0.0_R8
      delni=0.0_R8
      delna=0.0_R8
      delpo=0.0_R8
      delte=0.0_R8
      delti=0.0_R8
      delua=0.0_R8

      return
      end subroutine alloc_b2us_work

      subroutine dealloc_b2us_work

      if (.not.allocated(delne)) return
      deallocate(delne)
      deallocate(delni)
      deallocate(delna)
      deallocate(delpo)
      deallocate(delte)
      deallocate(delti)
      deallocate(delua)
      return
      end subroutine dealloc_b2us_work

      end module b2us_work
