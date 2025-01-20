CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ****s* INTERPOLATION/BILINEAR_INT
C NAME
C     BILINEAR_INT (D, RX, RY, Z)
C DESCRIPTION
C     Bilinear interpolation of data
C INPUTS
C     REAL*8, DIMENSION(2,2):: D  ! data to interpolate
C     REAL*8    :: RX, RY ! relative position at which interpolations should be done
C OUTPUT
C     REAL*8    :: Z  ! interpolated data
C     ******
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE BILINEAR_INT (D, RX, RY, Z)
      USE B2MOD_TYPES
      IMPLICIT NONE
 
!---- input variables
      REAL(kind=R8), DIMENSION(2,2), INTENT(IN):: D
      REAL(kind=R8), INTENT(IN)  :: RX, RY
!---- output variables
      REAL(kind=R8), INTENT(OUT)  :: Z
 
!---- local variables
      REAL(kind=R8)  :: Z1, Z2
 
      Z1 = D(1,1)*RX + D(2,1)*(1-RX)
      Z2 = D(1,2)*RX + D(2,2)*(1-RX)
      Z  = Z1*RY + Z2*(1-RY)
 
      RETURN
      END SUBROUTINE BILINEAR_INT
