      module b2mod_types
#ifdef PGF90
! Quadruple precision not available for Portland Group compilers
      INTEGER, PARAMETER :: R16 = SELECTED_REAL_KIND (14)
#else
      INTEGER, PARAMETER :: R16 = SELECTED_REAL_KIND (33)     !srv 16.01.18
#endif
      INTEGER, PARAMETER :: R8 = SELECTED_REAL_KIND (14)
      INTEGER, PARAMETER :: R4 = SELECTED_REAL_KIND (6)
      end module b2mod_types
