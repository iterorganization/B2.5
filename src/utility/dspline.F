        SUBROUTINE DSPLINE(N, X, Y, B, C, D)
	  use b2mod_types
	  implicit none
        INTEGER N
        REAL (kind=R8) :: X(N),Y(N),B(N),C(N),D(N)
C
C     calculates coefficients B(I),C(I), and D(I),
C     I=1,2,...,N, for cubic interpolation spline
C
C       S(X)=Y(I)+B(I)*(X-X(I))+C(I)*(X-X(I))**2+D(I)*(X-X(I))**3
C       for X(I).LE.X.LE.X(I+1)
C
C     INPUT:
C
C       N     = number of input nodes (N.GE.2)
C       X     = abscissae of the nodes in incremental order
C       Y     = ordinates of the nodes
C
C     OUTPUT:
C
C       B,C,D = arrays of spline coefficients
C
C     In the following lines, derivatives are indicated by a prime symbol
C
C       Y(I)=S(X(I))
C       B(I)=S'(X(I))
C       C(I)=S''(X(I))/2
C       D(I)=S'''(X(I))/6 (right-side derivative)
C
C     Subroutine SEVAL calculates the spline values
C
        INTEGER NM1,IB,I
        REAL (kind=R8) :: T
C
        NM1=N-1
        IF(N.LT.2) RETURN
        IF(N.LT.3) GO TO 50
C
C     Build the three-diagonal system
C
C     B=diagonal, D=upper diagonal, C=RHS
C
        D(1)=X(2)-X(1)
        C(2)=(Y(2)-Y(1))/D(1)
        DO 10 I=2,NM1
           D(I)=X(I+1)-X(I)
           B(I)=2.0e0_R8*(D(I-1)+D(I))
           C(I+1)=(Y(I+1)-Y(I))/D(I)
           C(I)=C(I+1)-C(I)
   10   CONTINUE
C
C     At the boundary points X(1) and X(N), the third derivatives
C     are calculated by separated differences
C
        B(1)=-D(1)
        B(N)=-D(N-1)
        C(1)=0.e0_R8
        C(N)=0.e0_R8
        IF(N.EQ.3) GO TO 15
        C(1)=C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
        C(N)=C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
        C(1)=C(1)*D(1)**2/(X(4)-X(1))
        C(N)=-C(N)*D(N-1)**2/(X(N)-X(N-3))
C
C     direct sweeping
C
   15   DO 20 I=2,N
           T=D(I-1)/B(I-1)
           B(I)=B(I)-T*D(I-1)
           C(I)=C(I)-T*C(I-1)
   20   CONTINUE
C
C     reversed substitution
C
        C(N)=C(N)/B(N)
        DO 30 IB=1,NM1
           I=N-IB
           C(I)=(C(I)-D(I)*C(I+1))/B(I)
   30   CONTINUE
C
C     now C(I) contains SIGMA(I) defined in Eq.(4.4) from
C     [G.E. Forsythe,et al., Computer methods for Math. Computations, Prentice-Hall 1977]
C
C     calculation of polynomial coefficients
C
        B(N)=(Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.0e0_R8*C(N))
        DO 40 I=1,NM1
           B(I)=(Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.0e0_R8*C(I))
           D(I)=(C(I+1)-C(I))/D(I)
           C(I)=3.0e0_R8*C(I)
   40   CONTINUE
        C(N)=3.0e0_R8*C(N)
        D(N)=D(N-1)
        RETURN
C
   50   B(1)=(Y(2)-Y(1))/(X(2)-X(1))
        C(1)=0.0e0_R8
        D(1)=0.0e0_R8
        B(2)=B(1)
        C(2)=0.0e0_R8
        D(2)=0.0e0_R8
        RETURN
        END
