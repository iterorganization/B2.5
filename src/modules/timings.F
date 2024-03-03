      module timings
#ifdef WITH_FTIMINGS
      use ftimings ! IGNORE
      type(timer_t) :: timer
#else
      integer :: x__module_timings_dummy_value
#endif
      end module timings
