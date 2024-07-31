      module b2mod_ysmp_sdrv                                    !sv 29.05.02
      use b2mod_types
      implicit none

      real (kind=R8), allocatable, save ::
     * r1_ysmp(:),
     * s1_ysmp(:),
     * a(:),
     * rsp(:)

#ifdef _OPENMP
!$OMP THREADPRIVATE(r1_ysmp, s1_ysmp, a, rsp)
#endif

      integer, allocatable, save ::
     * ia(:),
     * ja(:),
     * p(:),
     * ip(:),
     * isp(:),
     * indexes(:,:),
     * iasave(:),
     * jasave(:)
#ifdef _OPENMP
!$OMP THREADPRIVATE(ia, ja, p, ip, isp, indexes, iasave, jasave)
#endif

      integer, save :: nnd_sdrv, nsp, n_nz_sdrv, esp
      logical, save :: lfirst_sdrv=.true.
#ifdef _OPENMP
!$OMP THREADPRIVATE(nnd_sdrv, nsp, n_nz_sdrv, esp)
!$OMP THREADPRIVATE(lfirst_sdrv)
#endif

      contains

      subroutine alloc_b2mod_ysmp_sdrv
      implicit none

      if (.not.allocated(r1_ysmp)) then

        write(*,*) 'Allocating space for b2mod_ysmp_sdrv'
        write(*,*) 'nnd_sdrv   = ', nnd_sdrv
        write(*,*) 'nsp        = ', nsp
        write(*,*) 'n_nz_sdrv  = ', n_nz_sdrv

        allocate(r1_ysmp(1:nnd_sdrv))
        allocate(s1_ysmp(1:nnd_sdrv))
        allocate(a(1:n_nz_sdrv))
        allocate(rsp(1:nsp))
        allocate(ia(1:nnd_sdrv+1))
        allocate(ja(1:n_nz_sdrv))
        allocate(p(1:nnd_sdrv))
        allocate(ip(1:nnd_sdrv))
        allocate(isp(1:nsp))
        allocate(indexes(1:3,1:n_nz_sdrv))
        allocate(iasave(1:nnd_sdrv+1))
        allocate(jasave(1:n_nz_sdrv))

        write(*,*) 'Allocated space for b2mod_ysmp_sdrv'
      endif
      return
      end subroutine alloc_b2mod_ysmp_sdrv

      subroutine dealloc_b2mod_ysmp_sdrv

#ifdef _OPENMP
!$OMP PARALLEL
#endif
      if (allocated(r1_ysmp)) then
        deallocate(r1_ysmp)
        deallocate(s1_ysmp)
        deallocate(a)
        deallocate(rsp)
        deallocate(ia)
        deallocate(ja)
        deallocate(p)
        deallocate(ip)
        deallocate(isp)
        deallocate(indexes)
        deallocate(iasave)
        deallocate(jasave)
      endif
#ifdef _OPENMP
!$OMP END PARALLEL
#endif
      return
      end subroutine dealloc_b2mod_ysmp_sdrv

      subroutine restart_sdrv_for_7
      lfirst_sdrv=.true.
      write(*,*) ' Restarting sdrv'
      return
      end subroutine restart_sdrv_for_7

      end module b2mod_ysmp_sdrv                                 !sv 29.05.02
