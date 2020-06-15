SUBROUTINE DIM_D(x, xd, y, yd, res, resd)
  USE B2MOD_TYPES
  IMPLICIT NONE
!     ------------------------------------------------------------------
  REAL(kind=r8) :: res, x, y
  REAL(kind=r8) :: resd, xd, yd
!     ------------------------------------------------------------------
!      Differentiated form of intrinsic function dim in forward mode 
!     ------------------------------------------------------------------

  res = dim(x,y)
  IF (res .GE. 0.0_R8) THEN
    resd = xd - yd
  ELSE
    resd = 0.0_R8
  END IF

 RETURN
END SUBROUTINE DIM_D
