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
      end module