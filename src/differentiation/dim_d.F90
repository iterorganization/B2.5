FUNCTION DIM_D(x, xd, y, yd, res)
  USE B2MOD_TYPES
  IMPLICIT NONE
!     ------------------------------------------------------------------
  REAL(kind=r8) :: res, x, y
  REAL(kind=r8) :: dim_d, xd, yd
!     ------------------------------------------------------------------
!      Differentiated form of intrinsic function dim in forward mode 
!     ------------------------------------------------------------------

  res = dim(x,y)
  IF (res .GE. 0.0_R8) THEN
    dim_d = xd - yd
  ELSE
    dim_d = 0.0_R8
  END IF

 RETURN
END FUNCTION DIM_D
