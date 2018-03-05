      module b2mod_types
#ifdef NAGFOR
      use F90_KIND  ! IGNORE
      INTEGER, PARAMETER :: R16 = real128
#else
#ifdef PGF90
! Quadruple precision not available for Portland Group compilers
      INTEGER, PARAMETER :: R16 = SELECTED_REAL_KIND (14)
#else
      INTEGER, PARAMETER :: R16 = SELECTED_REAL_KIND (33)     !srv 16.01.18
#endif
#endif
      INTEGER, PARAMETER :: R8 = SELECTED_REAL_KIND (14)
      INTEGER, PARAMETER :: R4 = SELECTED_REAL_KIND (6)
      INTEGER, PARAMETER :: I8 = SELECTED_INT_KIND (10)
      end module b2mod_types
