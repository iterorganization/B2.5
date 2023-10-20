      module b2mod_openmp
#if !defined CRAY
!$      use omp_lib       ! IGNORE
#endif
      implicit none
      contains
      logical function in_parallel()
      in_parallel = .false.
!$      in_parallel = OMP_IN_PARALLEL()
      end function

      integer function get_num_threads()
      get_num_threads = 1
!$OMP parallel
!$OMP single
!$    get_num_threads = omp_get_num_threads()
!$OMP end single
!$OMP end parallel
      end function

      integer function get_thread_num()
      get_thread_num = 0
!$    get_thread_num = omp_get_thread_num()
      end function

#if !defined GFORTRAN && !defined NAGFOR && !defined PGF90
#if !defined CRAY
      integer function get_stacksize_s()
      get_stacksize_s = 0
!$    get_stacksize_s = kmp_get_stacksize_s()
#else
      function get_stacksize_s()
      character*12 get_stacksize_s
      character*12 stacksize
      logical streql
      external streql
      call getenv ('OMP_STACKSIZE', stacksize)
      if (streql(stacksize,' ')) stacksize = '128 MB'
      get_stacksize_s = trim(stacksize)
#endif
      end function
#endif

      subroutine print_openmp_info()
#ifdef _OPENMP
      if (get_num_threads().gt.1) then
        write(*,'(A,I3,A)') 'B2.5 running in OpenMP mode using ',
     &   get_num_threads(),' threads.'
      else
        write(*,'(A,I3,A)') 'B2.5 running in OpenMP mode using ',
     &   get_num_threads(),' thread.'
      end if
#if !defined GFORTRAN && !defined NAGFOR && !defined PGF90
      if (get_num_threads().gt.1)
     & write(*,*) ' Stack size allotted per thread is : ',
     & get_stacksize_s()
#endif
#else
      write(*,*) 'B2.5 running in serial mode, without ',
     & 'OpenMP threads.'
#endif
      end subroutine

      end module
