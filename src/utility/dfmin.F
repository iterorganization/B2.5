      FUNCTION DFMIN(AX,BX,F,TOL)
      use b2mod_types
      REAL (KIND=R8) :: DFMIN,AX,BX,F,TOL
C
C     calculates an approximattion X to the point where F
C     attains a minimum on the interval (AX,BX)
C
C     input.
C
C     AX   left endpoint of initial interval
C     BX   right endpoint of initial interval
C     F    function subprogramme which evaluates F(X)
C          for any X in the interval (AX,BX)
C    TOL   desired length of the interval of uncertainty
C          of the final result (.GE.0.0)
C
C    output.
C
C    DFMIN  abcissa approximating the point where
C          F attains a minimum
C
C      the method used is a combination of golden section search and
C     successive parabolic interpolation. convergence is never much slower
C     than that for a fibonacci search. if f has a continuous second
C     derivative which is positive at the minimum (which is not at ax or
C     bx), then convergence is superlinear, and usually of the order of
C     about 1.324....
C      the function f is never evaluated at two points closer together
C     than eps*abs(fmin) + (tol/3), where eps is approximately the square
C     root of the relative machine precision. if f is a unimodal
C     function and the computed values of f are always unimodal when
C     separated by at least eps*abs(x) + (tol/3), then fmin approximates
C     the abcissa of the global minimum of f on the interval [ax,bx] with
C     an error less than 3*eps*abs(fmin) + tol. if f is not unimodal,
C     then fmin may approximate a local, but perhaps non-global, minimum to
C     the same accuracy.
C      this function subprogram is a slightly modified version of the
C     algol 60 procedure localmin given in richard brent, algorithms for
C     minimization without derivatives, prentice - hall, inc. (1973).
C
C
      REAL (KIND=R8) :: A,B,C,D,E,EPS,XM,P,Q,R,TOL1,TOL2,U,V,W
      REAL (KIND=R8) :: FU,FV,FW,FX,X
C
C     C is squared value of (1/golden ratio)
C
      C=0.5_R8*(3.0_R8-SQRT(5.0_R8))
C
C     EPS is approximately square root of
C     relative machine precision 
C
      EPS=1.0_R8
   10 EPS=EPS/2.0_R8
      TOL1=1.0_R8+EPS
      IF(TOL1.GT.1.0_R8) GO TO 10
      EPS=SQRT(EPS)
C
C     initialization
C
      A=AX
      B=BX
      V=A+C*(B-A)
      W=V
      X=V
      E=0.0_R8
      FX=F(X)
      FV=FX
      FW=FX
C
C     main loop starts here
C
   20 XM=0.5_R8*(A+B)
      TOL1=EPS*ABS(X)+TOL/3.0_R8
      TOL2=2.0_R8*TOL1
C
C     check stopping criterion
C
      IF(ABS(X-XM).LE.(TOL2-0.5_R8*(B-A))) GO TO 90
C
C     check if golden ratio is necessary
C
      IF(ABS(E).LE.TOL1) GO TO 40
C
C     fit parabola
C
      R=(X-W)*(FX-FV)
      Q=(X-V)*(FX-FW)
      P=(X-V)*Q-(X-W)*R
      Q=2.0_R8*(Q-R)
      IF(Q.GT.0.0) P=-P
      Q=ABS(Q)
      R=E
      E=D
C
C     check if parabola acceptable
C
      IF(ABS(P).GE.ABS(0.5_R8*Q*R)) GO TO 40
      IF(P.LE.Q*(A-X)) GO TO 40
      IF(P.GE.Q*(B-X)) GO TO 40
C
C     step of a parabolic interpolation
C
      D=P/Q
      U=X+D
C
C     F must not be evaluated too close to AX or BX
C
      IF((U-A).LT.TOL2) D=SIGN(TOL1,XM-X)
      IF((B-U).LT.TOL2) D=SIGN(TOL1,XM-X)
      GO TO 50
C
C     golden ration step
C
   40 IF(X.GE.XM) E=A-X
      IF(X.LT.XM)E=B-X
      D=C*E
C
C     F must not be evaluated too close to X
C
   50 IF(ABS(D).GE.TOL1) U=X+D
      IF(ABS(D).LT.TOL1) U=X+SIGN(TOL1,D)
      FU=F(U)
C
C     update A,B,V,W and X
C
      IF(FU.GT.FX) GO TO 60
      IF(U.GE.X) A=X
      IF(U.LT.X) B=X
      V=W
      FV=FW
      W=X
      FW=FX
      X=U
      FX=FU
      GO TO 20
C
   60 IF(U.LT.X) A=U
      IF(U.GE.X) B=U
      IF(FU.LE.FW) GO TO 70
      IF(W.EQ.X) GO TO 70
      IF(FU.LE.FV) GO TO 80
      IF(V.EQ.X) GO TO 80
      IF(V.EQ.W) GO TO 80
      GO TO 20
C
   70 V=W
      FV=FW
      W=U
      FW=FU
      GO TO 20
C
   80 V=U
      FV=FU
      GO TO 20
C
C     end of main loop
C
   90 DFMIN=X
      RETURN
      END
