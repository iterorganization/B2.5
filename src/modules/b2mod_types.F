      module b2mod_types
#ifdef NAGFOR
      use F90_KIND  ! IGNORE
      INTEGER, PARAMETER :: R16 = real128
#else
#ifdef NOQUAD
! For when quadruple precision is not available
! such as Portland Group compilers and some gfortran versions
      INTEGER, PARAMETER :: R16 = SELECTED_REAL_KIND (14)
#else
      INTEGER, PARAMETER :: R16 = SELECTED_REAL_KIND (33)     !srv 16.01.18
#endif
#endif
      INTEGER, PARAMETER :: R8 = SELECTED_REAL_KIND (14)
      INTEGER, PARAMETER :: R4 = SELECTED_REAL_KIND (6)
      INTEGER, PARAMETER :: I8 = SELECTED_INT_KIND (10)
      end module b2mod_types
