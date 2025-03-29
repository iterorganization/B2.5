      module b2mod_ysmp_sdrv                                    !sv 29.05.02

      use b2mod_types
      implicit none
      real (kind=R8), allocatable, save ::
     * rr1(:),
     * s1(:),
     * a(:),
     * rsp(:)

      integer, allocatable, save ::
     * ia(:),
     * ja(:),
     * p(:),
     * ip(:),
     * isp(:),
     * indexes(:,:),
     * iasave(:),
     * jasave(:)

      integer nnd_sdrv, nsp, n_nz_sdrv, esp
      logical, save :: lfirst_ill7 = .true.


#ifdef _OPENMP
!$omp threadprivate(rr1, s1, a, rsp)
!$omp threadprivate(ia, ja, p, ip, isp, indexes, iasave, jasave)
!$omp threadprivate(nnd_sdrv, nsp, n_nz_sdrv, esp)
!$omp threadprivate(lfirst_ill7)
#endif

      contains

      subroutine alloc_b2mod_ysmp_sdrv
      implicit none

      write(*,*) 'Allocating space for b2mod_ysmp_sdrv'
      write(*,*) 'nnd_sdrv   = ', nnd_sdrv
      write(*,*) 'nsp        = ', nsp
      write(*,*) 'n_nz_sdrv  = ', n_nz_sdrv

      allocate(rr1(1:nnd_sdrv))
      allocate(s1(1:nnd_sdrv))
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

      end subroutine alloc_b2mod_ysmp_sdrv

      subroutine dealloc_b2mod_ysmp_sdrv

      if (.not.allocated(rr1)) return
      deallocate(rr1)
      deallocate(s1)
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

      end subroutine dealloc_b2mod_ysmp_sdrv

      subroutine restart_sdrv_for_7
      implicit none                            
      lfirst_ill7=.true.
      write(*,*) ' Restarting sdrv'
      return
      end subroutine restart_sdrv_for_7

      end module b2mod_ysmp_sdrv                                 !sv 29.05.02
