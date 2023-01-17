
      SUBROUTINE TRILINEAR_INT (D, RX, RY, RZ, Z)
      USE B2MOD_TYPES
      IMPLICIT NONE
 
!---- input variables
      REAL(kind=R8), DIMENSION(2,2,2), INTENT(IN):: D
      REAL(kind=R8), INTENT(IN)  :: RX, RY, RZ
!---- output variables
      REAL(kind=R8), INTENT(OUT)  :: Z
 
!---- local variables
      REAL(kind=R8)  :: Z11, Z12, Z21, Z22, Z1, Z2
 
      Z11 = D(1,1,1)*RX + D(2,1,1)*(1-RX)
      Z21 = D(1,2,1)*RX + D(2,2,1)*(1-RX)

      Z12 = D(1,1,2)*RX + D(2,1,2)*(1-RX)
      Z22 = D(1,2,2)*RX + D(2,2,2)*(1-RX)

      Z1  = Z11*RY + Z21*(1-RY)
      Z2  = Z12*RY + Z22*(1-RY)

      Z  = Z1*RZ + Z2*(1-RZ)
 
      RETURN
      END SUBROUTINE TRILINEAR_INT
