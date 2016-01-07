module tradui_constants
  
  use carre_constants
  implicit none

  private

  integer, public, parameter :: nisomx=1

  ! currently, cell flags computed in Carre are just passed through in traduit
  integer, public, parameter :: GRID_N_CELLFLAGS = CARREOUT_NCELLFLAGS
  
end module tradui_constants

