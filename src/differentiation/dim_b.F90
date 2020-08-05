SUBROUTINE DIM_B(x, xb, y, yb, resb)
  USE B2MOD_TYPES
  IMPLICIT NONE
!     ------------------------------------------------------------------
  REAL(kind=r8) :: x, y
  REAL(kind=r8) :: resb, xb, yb
!     ------------------------------------------------------------------
!      B2.5 representation of intrinsic function dim which is 
!      needed for differentiation with TAPENADE       
!     ------------------------------------------------------------------
  REAL(kind=r8) :: dummy
  REAL(kind=r8) :: dummyb
!
  dummy = x - y
  IF (dummy .GT. 0.0_R8) THEN
    dummyb = resb
  ELSE
    dummyb = 0.0
  END IF
  xb = xb + dummyb
  yb = yb - dummyb
END SUBROUTINE DIM_B

