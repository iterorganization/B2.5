      SUBROUTINE EIRENE_EIRSRT(LSTOP_in,LTIME_in,DELTAT_in,FLUXES_in,
     .                  B2BRM_in,B2RD_in,B2Q_in,B2VP_in,STEP_CPU_in)
      use eirmod_precision
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: FLUXES_in
      REAL(DP), INTENT(IN) :: DELTAT_in, B2BRM_in, B2RD_in, B2Q_in,
     .                        B2VP_in,STEP_CPU_in
      LOGICAL, INTENT(IN) :: LSTOP_in, LTIME_in
      END
