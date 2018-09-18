      FUNCTION DSEVAL(N,U,X,Y,B,C,D)
      use b2mod_types
      implicit none
      INTEGER N
      REAL (KIND=R8) :: U,X(N),Y(N),B(N),C(N),D(N), DSEVAL
C
C  this subprogramme calculates the value of cubic
C  spline
C
C  DSEVAL=Y(I)+B(I)*(U-X(I))+C(I)*(U-X(I))**2+D(I)*(U-X(I))**3
C
C  where  X(I).LT.U.LT.X(I+1). The Gorner scheme is used.
C
C  if  U.LT.X(1), then I=1.
C  if  U.GE.X(N), then I=N.
C
C  Input:
C
C     N     -number of given points
C     U     -abscissa, for which value of spline is calculated
C     X,Y   -arrays of given abscissa and ordinates
C     B,C,D -arrays of given spline coefficients, calculated
C            with subprogramme 'spline'
C
C  if in comparison with previous call U is not
C  in the same interval, then for determination
C  proper interval double (binary) search should be used.
C
      INTEGER I,J,K
      REAL (KIND=R8) :: DX
      save i
      DATA I/1/
      IF(I.GE.N) I=1
      IF(U.LT.X(I)) GO TO 10
      IF(U.LE.X(I+1)) GO TO 30
C
C  double (binary) search
C
   10 I=1
      J=N+1
   20 K=(I+J)/2
      IF(U.LT.X(K))J=K
      IF(U.GE.X(K))I=K
      IF(J.GT.I+1)GO TO 20
C
C  calculate spline
C
   30 DX=U-X(I)
      DSEVAL=Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
      RETURN
      END
