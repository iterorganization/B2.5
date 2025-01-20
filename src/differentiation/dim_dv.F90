SUBROUTINE DIM_DV(x, xd, y, yd, res, resd, nbdirs)
  USE B2MOD_TYPES
  IMPLICIT NONE
!     ------------------------------------------------------------------
  REAL(kind=r8) :: res, x, y
  REAL(kind=r8) :: resd(nbdirsmax), xd(nbdirsmax), yd(nbdirsmax)
  INTEGER :: nbdirs, nd
!     ------------------------------------------------------------------
!      Differentiated form of intrinsic function dim in forward vector mode 
!     ------------------------------------------------------------------

  res = dim(x,y)
  DO nd = 1, nbdirs
    IF (res .GT. 0.0_R8) THEN
      resd(nd) = xd(nd) - yd(nd)
    ELSE
      resd(nd) = 0.0_R8
    END IF
  END DO
 RETURN
END SUBROUTINE DIM_DV
