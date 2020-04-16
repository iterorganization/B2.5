      module b2mod_openmp
!$      use omp_lib       ! IGNORE
      implicit none
      contains
      logical function in_parallel()
#ifdef _OPENMP
      in_parallel = OMP_IN_PARALLEL()
#else
      in_parallel = .false.
#endif
      end function

      integer function get_num_threads()
      get_num_threads = 1
!$OMP parallel
!$OMP single
!$    get_num_threads = omp_get_num_threads()
!$OMP end single
!$OMP end parallel
      end function

      subroutine print_openmp_info()
#ifdef _OPENMP
      write(*,fmt='(A,I3,A)') 'B2.5 running in OpenMP mode using ',
     &  get_num_threads(),' threads.'
#else
      write(*,*) 'B2.5 running in serial mode, without ',
     & 'OpenMP threads.'
#endif
      end subroutine

      end module
