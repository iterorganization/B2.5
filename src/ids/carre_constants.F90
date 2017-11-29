module carre_constants

  implicit none

  private

  integer, parameter, public :: GRID_UNDEFINED = 0

  ! Constants used to categorize grid objects
  ! These are used for points, faces and cells
  integer, parameter, public :: GRID_INTERNAL = 1
  integer, parameter, public :: GRID_EXTERNAL = 2
  integer, parameter, public :: GRID_BOUNDARY = 3
  integer, parameter, public :: GRID_BOUNDARY_REFINE = 4
  integer, parameter, public :: GRID_BOUNDARY_REFINE_FIX = 5
  integer, parameter, public :: GRID_INTERNAL_COARSEN = 6
  integer, parameter, public :: GRID_BOUNDARY_COARSEN = 7
  integer, parameter, public :: GRID_REFINE = 8
  integer, parameter, public :: GRID_GUARD = 9
  integer, parameter, public :: GRID_DEAD = 10

  ! Number of flags written out for every cell in cellflags blog of carre.out
  integer, parameter, public :: CARREOUT_NCELLFLAGS = 5
  integer, parameter, public :: CELLFLAG_TYPE = 1
  integer, parameter, public :: CELLFLAG_LEFTFACE = 2
  integer, parameter, public :: CELLFLAG_BOTTOMFACE = 3
  integer, parameter, public :: CELLFLAG_RIGHTFACE = 4
  integer, parameter, public :: CELLFLAG_TOPFACE = 5

  ! Face number marking face on a generic boundary (i.e. not on a specific
  ! structure). This is used to mark boundary faces in standard grids that
  ! do not touch the wall.
  integer, parameter, public :: BOUNDARY_NOSTRUCTURE = -1

end module carre_constants

!!!Local Variables:
!!! mode: f90
!!! End:
