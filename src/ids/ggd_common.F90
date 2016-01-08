module ggd_common

  use b2mod_types
  use ggd_assert

  implicit none      

  integer, parameter :: GRID_UNDEFINED = 0

  ! Data representation definitions
  integer, parameter :: GEO_TYPE_STANDARD = 1
  character(*), parameter :: GEO_TYPE_ID_STANDARD = "Standard"
  integer, parameter :: GEO_TYPE_FOURIER = 2
  character(*), parameter :: GEO_TYPE_ID_FOURIER = "Fourier"

  ! Coordinate type definitions

  ! Cartesian coordinates
  integer, parameter :: COORDTYPE_X = 1 ! (m)
  integer, parameter :: COORDTYPE_Y = 2 ! (m)

  ! The following are part of the ITM convention machine coordinate system
  integer, parameter :: COORDTYPE_R = 4 ! Major radius (m)
  integer, parameter :: COORDTYPE_Z = 5 ! Vertical height (m) 
  integer, parameter :: COORDTYPE_PHI = 6 ! Toroidal angle (rad)

  integer, parameter :: COORDTYPE_PSI = 7   ! Radial flux coordinate
  integer, parameter :: COORDTYPE_THETA = 8 ! Poloidal angle

  ! From 4.09a distribution/distri_vec/dist_func/markers/var_id
  integer, parameter :: COORDTYPE_RHOTOR = 107 ! Square root of the toroidal flux
  integer, parameter :: COORDTYPE_THETAB = 109 ! Boozer poloidal angle     (rad)
  integer, parameter :: COORDTYPE_VX     = 110 ! velocity in X-direction   (m/s)
  integer, parameter :: COORDTYPE_VY     = 111 ! velocity in Y-direction   (m/s)
  integer, parameter :: COORDTYPE_VZ     = 112 ! velocity in Z-direction   (m/s)
  integer, parameter :: COORDTYPE_VEL    = 113 ! total velocity            (m/s)
  integer, parameter :: COORDTYPE_VPHI   = 114 ! velocity in PHI-direction (m/s)
  integer, parameter :: COORDTYPE_VPAR   = 115 ! velocity in parallel-direction (m/s)
  integer, parameter :: COORDTYPE_VPERP  = 116 ! velocity in perpendicular-direction (m/s)
  integer, parameter :: COORDTYPE_E      = 117 ! Hamiltonian energy         (J)
  integer, parameter :: COORDTYPE_Pphi   = 118 ! Canonical toroidal angular momentum (kg m^2/s)
  integer, parameter :: COORDTYPE_mu     = 119 ! magnetic moment           (J/T)
  integer, parameter :: COORDTYPE_lambda = 120 ! mu/E                      (1/T)
  integer, parameter :: COORDTYPE_pitch  = 121 ! vpar/v                     ( )
  integer, parameter :: COORDTYPE_OMNIG  = 122 ! position of the omnigenous plane (generalised equatorial plane)
  integer, parameter :: COORDTYPE_SPIN   = 123 ! particle spin



 ! Field aligned vector definitions
  integer, parameter :: VEC_ALIGN_DEFAULT = 1
  character(len=132), parameter :: VEC_ALIGN_DEFAULT_ID = "DEFAULT"

  integer, parameter :: VEC_ALIGN_POLOIDAL = 1001
  character(len=132), parameter :: VEC_ALIGN_POLOIDAL_ID = "Poloidal"
  integer, parameter :: VEC_ALIGN_RADIAL = 1002
  character(len=132), parameter :: VEC_ALIGN_RADIAL_ID = "Radial"
  integer, parameter :: VEC_ALIGN_PARALLEL = 1003
  character(len=132), parameter :: VEC_ALIGN_PARALLEL_ID = "Parallel"

  integer, parameter :: VEC_ALIGN_TOROIDAL = 1004
  character(len=132), parameter :: VEC_ALIGN_TOROIDAL_ID = "Toroidal"


!contains

  


end module ggd_common

!!!Local Variables:
!!! mode: f90
!!! End:
