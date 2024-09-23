c    Subroutine MSTEP solves initial value problems in a system
c    of ODEs. It is a fourth order multistep method with
c    adaptive step size control (Ref: Numerical methods for Scientists
c    and Engineers - H. M. Antia 1991, Tata Mc-Graw Hill publishers).

cRevision 1.1  2000/07/04 07:55:37  mkundu
cInitial revision
c
c$Id: mstep.F 1032 2007-05-21 09:07:53Z dpc $
c_______________________________________________________________________
      SUBROUTINE MSTEP(IBND,NS,N,Y,YT,DY,H,T0,TN,YF,REPS,NSTP,NMAX,
     1  IER,IFLG,iopt,WK,IWK)
      use b2mod_types

c  Added iopt as a input variable instead of passing it thru common
c  as done by Mrityunjay Kundu in his IPR report.

      implicit none
      integer ibnd,n,ns,nstp,nmax,ier,iflg,iopt,nmx
      real (kind=R8) :: YT(N),Y(N,7),DY(N,7),WK(N,3+2*N),
     1  YF(N),DEL(4),H,T0,TN,REPS,EPS
      PARAMETER(NMX=100000,EPS=1.0e-30_R8)
      integer IWK(N+1)

      integer i, ij, ijm1, ijm2, ijm3, ijm4, iflag, j1, j2, i2, nst
      real (kind=R8) :: TSTEP, T, REPS1, R1, R2, RERR, Y1, Y2, TX
      SAVE

      IF(N.LE.0) THEN
        IER=111
        RETURN
      ENDIF
      IF(NMAX.LE.0) NMAX=NMX
      IER=0
      TSTEP=TN-T0

      IF(IFLG.LE.1) THEN
        IF(TSTEP.EQ.0.0_R8) then
!Modification starts.
!Returing Initial values.
          DO I = 1,N
            YF(I) = YT(I)
          ENDDO
          RETURN
        ENDIF
!Modification ends.
        NSTP=0
        IF(H.EQ.0.0_R8) H=TSTEP/8.0_R8
        IF((H.LT.0.0_R8).EQV.(TN.GT.T0)) H=-H

        T=T0
        CALL STRT4(IBND,NS,N,Y,DY,H,T,REPS,IFLG,TSTEP,NSTP,IER,WK)
        IF(IER.GT.0) RETURN
        IJ=4
        IJM1=3
        IJM2=2
        IJM3=1
        T0=T
        IFLAG=0
        IFLG=IFLG+2
        NST=4
      ENDIF

      IF((TN.LE.T0).EQV.(H.GT.0.0d0)) GO TO 6000
      IF(TSTEP.EQ.0.0_R8) GO TO 6000
 3400 IJM4=IJM3
      IJM3=IJM2
      IJM2=IJM1
      IJM1=IJ
      IJ=IJ+1
      IF(IJ.GT.7) IJ=1
      T=T0+H
      REPS1=ABS(REPS*H/TSTEP)

      if(iopt.eq.1) then
        CALL GEAR(IBND,NS,N,Y,DY,H,T,REPS1,NSTP,IJ,IJM1,IJM2,
     1          IJM3,IJM4,IFLAG,IER,WK,WK(1,4),IWK)
      else
        CALL ADAMS(IBND,NS,N,Y,DY,H,T,REPS1,NSTP,IJ,IJM1,IJM2,
     1           IJM3,IJM4,IER,WK)
      endif

      RERR=0
      DO 4500 I=1,N
        R2=ABS(Y(I,IJ))+ABS(Y(I,IJ)-Y(I,IJM1))+EPS
        R1=ABS(WK(I,1)-Y(I,IJ))/R2
        IF(R1.GT.RERR) RERR=R1
 4500 CONTINUE
      RERR=ABS(0.25_R8*RERR*TSTEP/H)

      IF(T.EQ.T0) THEN
        IER=123
        DO 4600 I=1,N
          YF(I)=Y(I,IJ)
 4600   CONTINUE
        RETURN
      ENDIF

      IF(IFLG.LE.2) THEN
        IF(RERR.LT.REPS.AND.IER.EQ.0) THEN
          T0=T0+H
          DO 4700 I=1,N
            WK(I,2)=WK(I,1)
 4700     CONTINUE
          NST=NST+1
          IF((TN.LE.T0).EQV.(H.GT.0.0_R8)) GO TO 6000
          IF(ABS((TN-T0)/TSTEP).LE.REPS) GO TO 6000

          IF(RERR.LT.REPS/32.0_R8.AND.NST.GT.7) THEN
            H=2.0_R8*H
            I2=IJM4-2
            IF(I2.LE.0) I2=I2+7
            DO 4800 I=1,N
              Y(I,IJM1)=Y(I,IJM2)
              Y(I,IJM2)=Y(I,IJM4)
              Y(I,IJM3)=Y(I,I2)
              DY(I,IJM1)=DY(I,IJM2)
              DY(I,IJM2)=DY(I,IJM4)
              DY(I,IJM3)=DY(I,I2)
 4800       CONTINUE
            NST=4
          ENDIF

        ELSE
          IF(NST.GT.-4.0_R8.AND.RERR.LT.16.0_R8*REPS) THEN
            DO 5000 I=1,N

              Y1=(45.0_R8*Y(I,IJM1)+72.0_R8*Y(I,IJM2) +
     &         11.0_R8*Y(I,IJM3))/128.0_R8 + H*(-9.0_R8*DY(I,IJM1) +
     &         36.0_R8*DY(I,IJM2)+3.0_R8*DY(I,IJM3))/128.0_R8

              Y2=(11.0_R8*Y(I,IJM1)+72.0_R8*Y(I,IJM2) +
     &         45.0_R8*Y(I,IJM3))/128.0_R8 + H*(-3.0_R8*DY(I,IJM1) -
     &         36.0_R8*DY(I,IJM2)+9.0_R8*DY(I,IJM3))/128.0_R8

              Y(I,IJ)=Y(I,IJM1)
              Y(I,IJM1)=Y1
              Y(I,IJM3)=Y2
              DY(I,IJ)=DY(I,IJM1)
 5000       CONTINUE
            T=T0-H/2.0_R8
            CALL lines2c(IBND,NS,T,N,Y(1,IJM1),DY(1,IJM1))
            T=T-H
            CALL lines2c(IBND,NS,T,N,Y(1,IJM3),DY(1,IJM3))
            NSTP=NSTP+2
            H=H/2.0_R8
          ELSE

            IFLG=IFLG-2
            T=T0
            H=H/8.0_R8
            DO 5200 I=1,N
              Y(I,1)=Y(I,IJM1)
 5200       CONTINUE
            CALL STRT4(IBND,NS,N,Y,DY,H,T,REPS,IFLG,TSTEP,NSTP,IER,
     .                 WK)
            IF(IER.GT.0) THEN
              DO 5300 I=1,N
                YF(I)=Y(I,IJM1)
 5300         CONTINUE
              RETURN
            ENDIF
            IJ=4
            IJM1=3
            IJM2=2
            IJM3=1
            T0=T
            IFLG=IFLG+2
            IFLAG=0
          ENDIF

          NST=-4
          IF(ABS(H/TSTEP).LT.REPS) THEN
            print*,'came here',h,' ',tstep,' ',reps
            IER=122
            DO 5400 I=1,N
              YF(I)=Y(I,IJM1)
 5400       CONTINUE
            RETURN
          ENDIF

        ENDIF
      ELSE
        IF(IER.GT.0) THEN
          DO 5600 I=1,N
            YF(I)=Y(I,IJM1)
 5600     CONTINUE
          RETURN
        ENDIF
        T0=T0+H
        DO 5700 I=1,N
          WK(I,2)=WK(I,1)
 5700   CONTINUE
        IF((TN.LE.T0).EQV.(H.GT.0.0_R8)) GO TO 6000
        IF(ABS((TN-T0)/TSTEP).LE.REPS) GO TO 6000
      ENDIF

      IF(NSTP.LT.NMAX) GO TO 3400
      IER=133
      DO 5800 I=1,N
        YF(I)=Y(I,IJ)
 5800 CONTINUE
      RETURN

 6000 TX=(TN-T0)/H
      DO 7000 I=1,N
        DEL(1)=DY(I,IJ)*H
        DEL(2)=Y(I,IJ)-Y(I,IJM1)
        DEL(3)=DY(I,IJM1)*H
        DEL(4)=Y(I,IJM1)-Y(I,IJM2)
        DO 6501 J1=2,3
          DO 6500 J2=4,J1,-1
            DEL(J2)=DEL(J2-1)-DEL(J2)
 6500     CONTINUE
 6501   CONTINUE
        DEL(4)=DEL(3)-0.50_R8*DEL(4)
        YF(I)=Y(I,IJ)+TX*(DEL(1)+TX*(DEL(2)+(TX+1)*(DEL(3)+
     1        (TX+1)*DEL(4)/2.0_R8)))
 7000 CONTINUE
      RETURN
      END SUBROUTINE MSTEP
c_______________________________________________________________________
      SUBROUTINE STRT4(IBND,NS,N,Y,DY,H,T,REPS,IFLG,TSTEP,NSTP,IER,
     1 WK)
      use b2mod_types
      implicit none
      integer n,nit,iflg,nstp,ier,ibnd,ns
      real (kind=R8) :: eps, sfac, h, t, reps, tstep
      PARAMETER(EPS=1.0e-30_R8,SFAC=0.90_R8,NIT=20)
      real (kind=R8) :: Y(N,4),DY(N,4),WK(N,2)

      integer i, it
      real (kind=R8) :: T0, H2, RERR, R1, R2

      T0=T
      IER=0
      DO 4000 IT=1,NIT
        CALL lines2c(IBND,NS,T0,N,Y(1,1),DY(1,1))
        CALL RK4(IBND,NS,N,T0,Y(1,1),DY(1,1),H,Y(1,2),WK)
        T=T0+H
        CALL lines2c(IBND,NS,T,N,Y(1,2),DY(1,2))
        CALL RK4(IBND,NS,N,T,Y(1,2),DY(1,2),H,Y(1,3),WK)
        H2=2.0_R8*H
        CALL RK4(IBND,NS,N,T0,Y(1,1),DY(1,1),H2,Y(1,4),WK)
        NSTP=NSTP+11

        RERR=0
        DO 2000 I=1,N
          R2=ABS(Y(I,1))+ABS(Y(I,2))+ABS(Y(I,3))+EPS
          R1=ABS(Y(I,3)-Y(I,4))/R2
          IF(R1.GT.RERR) RERR=R1
 2000   CONTINUE
        RERR=RERR*TSTEP/H

        IF(RERR.LE.REPS.OR.IFLG.GT.0) THEN
          T=T+H
          CALL lines2c(IBND,NS,T,N,Y(1,3),DY(1,3))
          CALL RK4(IBND,NS,N,T,Y(1,3),DY(1,3),H,Y(1,4),WK)
          T=T+H
          CALL lines2c(IBND,NS,T,N,Y(1,4),DY(1,4))
          NSTP=NSTP+5
          DO 3000 I=1,N
            WK(I,2)=Y(I,4)
 3000     CONTINUE
          RETURN
        ELSE
          H=SFAC*H*(REPS/RERR)**0.250_R8
          IF(ABS(H/TSTEP).LT.REPS.OR.T.EQ.T0) THEN
            IER=113
            RETURN
          ENDIF
        ENDIF
 4000 CONTINUE
      IER=115
      RETURN
      END SUBROUTINE STRT4

c_______________________________________________________________________
      SUBROUTINE RK4(IBND,NS,N,T,Y0,DY0,H,Y1,WK)
      use b2mod_types
      implicit none
      integer n, ns, ibnd
      external lines2c
      real (kind=R8) :: T,H
      real (kind=R8) :: Y0(N),DY0(N),Y1(N),WK(N,2)

      integer i
      real (kind=R8) :: H2, T1

      H2=H/2.0_R8
      DO 1000 I=1,N
        Y1(I)=H2*DY0(I)
        WK(I,1)=Y0(I)+Y1(I)
 1000 CONTINUE

      T1=T+H2
      CALL lines2c(IBND,NS,T1,N,WK,WK(1,2))
      DO 1200 I=1,N
        Y1(I)=Y1(I)+H*WK(I,2)
        WK(I,1)=Y0(I)+H2*WK(I,2)
 1200 CONTINUE

      CALL lines2c(IBND,NS,T1,N,WK,WK(1,2))
      DO 1400 I=1,N
        Y1(I)=Y1(I)+H*WK(I,2)
        WK(I,1)=Y0(I)+H*WK(I,2)
 1400 CONTINUE

      T1=T+H
      CALL lines2c(IBND,NS,T1,N,WK,WK(1,2))
      DO 1600 I=1,N
        Y1(I)=Y0(I)+(Y1(I)+H2*WK(I,2))/3.0_R8
 1600 CONTINUE
      RETURN
      END SUBROUTINE RK4
c_______________________________________________________________________
      SUBROUTINE GEAR(IBND,NS,N,Y,DY,H,T,REPS,NSTP,IJ,IJM1,IJM2,
     1                IJM3,IJM4,IFLAG,IER,WK1,WK,IWK)
      use b2mod_types
      implicit none
      integer ns, n, nit, nstp, ij, ijm1, ijm2, ijm3, ijm4
      integer iflag, ier, ibnd
      real (kind=R8) :: eps, cfac, H, T, reps
      PARAMETER(EPS=1.0e-30_R8,NIT=20,CFAC=0.01_R8)
      real (kind=R8) :: Y(N,7),DY(N,7),WK1(N,3),WK(N,3+2*N)
      integer IWK(N)

      integer i, j, k, lj, ipas, ip, num, iflg, ier1
      real (kind=R8) :: T1, X1, DK, S1, S2, SS, R2, RERR, DET

      IER=0
      DO 1000 J=1,N

        T1=Y(J,IJM1)+H*(55.0_R8*DY(J,IJM1)-59.0_R8*DY(J,IJM2) +
     &  37.0_R8*DY(J,IJM3) - 9.0_R8*DY(J,IJM4))/24.0_R8

        Y(J,IJ)=T1+6275.0_R8*(Y(J,IJM1)-WK1(J,2))/8003.0_R8

        WK1(J,1)=T1
 1000 CONTINUE

      LJ=N
      IPAS=0
      IER=0
      NSTP=NSTP+1
      CALL lines2c(IBND,NS,T,N,Y(1,IJ),DY(1,IJ))
      DO 1200 J=1,N
        WK1(J,3)=Y(J,IJ)-(48.0_R8*Y(J,IJM1)-36.0_R8*Y(J,IJM2) +
     &  16.0_R8*Y(J,IJM3) - 3.0_R8*Y(J,IJM4)+12.0_R8*H*DY(J,IJ))/25.0_R8
 1200 CONTINUE

      IF(IFLAG.EQ.0) THEN
        DO 2000 IP=1,N
          X1=Y(IP,IJ)
          DK=100.0_R8*REPS*ABS(X1)
          IF(DK.EQ.0.0_R8) DK=100.0_R8*REPS
          Y(IP,IJ)=X1+DK
          NSTP=NSTP+1
          CALL lines2c(IBND,NS,T,N,Y(1,IJ),DY(1,IJ))
          DO 1800 J=1,N
            WK(J,1)=Y(J,IJ)-(48.0_R8*Y(J,IJM1)-36.0_R8*Y(J,IJM2) +
     &       16.0_R8*Y(J,IJM3) - 3.0_R8*Y(J,IJM4) +
     &       12.0_R8*H*DY(J,IJ))/25.0_R8
           WK(J,N+IP)=(WK(J,1)-WK1(J,3))/DK
 1800     CONTINUE
          Y(IP,IJ)=X1
 2000   CONTINUE

        NUM=N
        IFLG=0
        DO 2500 I=1,N
          DO 2400 J=1,N
            WK(J,I)=0.0_R8
 2400     CONTINUE
          WK(I,I)=1.0_R8
 2500   CONTINUE

        CALL GAUELM(N,NUM,WK(1,N+1),WK,DET,IWK,LJ,IER1,IFLG)
        IF(IER1.GT.0) THEN
          IER=125
          RETURN
        ENDIF
        IFLAG=1
      ENDIF

      DO 2800 I=1,N
        WK(I,N+1)=WK1(I,3)
        WK(I,N+2)=Y(I,IJ)
 2800 CONTINUE

      DO 5000 IPAS=1,NIT

        RERR=0.0_R8
        DO 3500 J=1,N
          S1=0.0_R8
          DO 3200 K=1,N
            S1=S1+WK(J,K)*WK(K,N+1)
 3200     CONTINUE
          R2=ABS(Y(J,IJ))+ABS(Y(J,IJ)-Y(J,IJM1))+EPS
          RERR=MAX(RERR,ABS(S1/R2))
          Y(J,IJ)=WK(J,N+2)-S1
          WK(J,N+3)=-S1
 3500   CONTINUE

        NSTP=NSTP+1
        CALL lines2c(IBND,NS,T,N,Y(1,IJ),DY(1,IJ))
        IF(RERR.LT.CFAC*REPS) RETURN
        DO 3800 J=1,N
          WK1(J,3)=Y(J,IJ)-(48.0_R8*Y(J,IJM1)-36.0_R8*Y(J,IJM2)+
     &     16.0_R8*Y(J,IJM3)- 3.0_R8*Y(J,IJM4)+
     &     12.0_R8*H*DY(J,IJ))/25.0_R8
 3800   CONTINUE

        DO 4000 J=1,N
          WK(J,N+4)=WK1(J,3)-WK(J,N+1)
          WK(J,N+1)=WK1(J,3)
          WK(J,N+2)=Y(J,IJ)
 4000   CONTINUE
        SS=0.0_R8
        DO 4300 J=1,N
          S1=0.0_R8
          S2=0.0_R8
          DO 4100 K=1,N
            S1=S1+WK(K,N+3)*WK(K,J)
            S2=S2+WK(J,K)*WK(K,N+4)
 4100     CONTINUE
          WK1(J,3)=S1
          Y(J,IJ)=S2-WK(J,N+3)
          SS=SS+S1*WK(J,N+4)
 4300   CONTINUE

        IF(SS.EQ.0.0_R8) THEN
          IER=125
          RETURN
        ENDIF
        DO 4400 J=1,N
          DO 4440 K=1,N
            WK(K,J)=WK(K,J)-Y(K,IJ)*WK1(J,3)/SS
 4440     CONTINUE
 4400   CONTINUE
 5000 CONTINUE
      IER=125
      RETURN
      END SUBROUTINE GEAR
c_______________________________________________________________________
      SUBROUTINE GAUELM(N,NUM,A,X,DET,IDINT,LJ,IER,IFLG)
      use b2mod_types
      implicit none
      integer n, num, lj, ier, iflg
      real (kind=R8) :: det
      real (kind=R8) :: A(LJ,N),X(LJ,NUM)
      integer IDINT(N)

      integer j, k, l, l1, km
      real (kind=R8) :: R1, T1

      IF(N.LE.0.OR.N.GT.LJ) THEN
        IER=111
        RETURN
      ENDIF

      IER=122
      IF(IFLG.LE.1) THEN
        DET=1.0_R8
        DO 2600 K=1,N-1
          R1=0.0_R8
          DO 2200 L=K,N
            IF(ABS(A(L,K)).GT.R1) THEN
              R1=ABS(A(L,K))
              KM=L
            ENDIF
 2200     CONTINUE

          IDINT(K)=KM
          IF(KM.NE.K) THEN
            DO 2300 L=K,N
              T1=A(K,L)
              A(K,L)=A(KM,L)
              A(KM,L)=T1
 2300       CONTINUE
            DET=-DET
          ENDIF

          DET=DET*A(K,K)
          IF(A(K,K).EQ.0.0_R8) RETURN
C         IF(ABS(A(K,K)).LT.REPS) RETURN
          DO 2501 L=K+1,N
            A(L,K)=A(L,K)/A(K,K)
            DO 2500 L1=K+1,N
              A(L,L1)=A(L,L1)-A(L,K)*A(K,L1)
 2500       CONTINUE
 2501     CONTINUE

 2600   CONTINUE
        DET=DET*A(N,N)
        IDINT(N)=N
        IF(A(N,N).EQ.0.0_R8) RETURN
C       IF(ABS(A(N,N)).LT.REPS) RETURN

        IER=0
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF

      IER=0
      DO 5000 J=1,NUM
        DO 3000 K=1,N-1
          IF(K.NE.IDINT(K)) THEN
            T1=X(K,J)
            X(K,J)=X(IDINT(K),J)
            X(IDINT(K),J)=T1
          ENDIF
          DO 3100 L=K+1,N
            X(L,J)=X(L,J)-A(L,K)*X(K,J)
 3100     CONTINUE
 3000   CONTINUE

        X(N,J)=X(N,J)/A(N,N)
        DO 3300 K=N-1,1,-1
          DO 3200 L=N,K+1,-1
            X(K,J)=X(K,J)-X(L,J)*A(K,L)
 3200     CONTINUE
          X(K,J)=X(K,J)/A(K,K)
 3300   CONTINUE
 5000 CONTINUE
      RETURN
      END SUBROUTINE GAUELM
c_______________________________________________________________________
      SUBROUTINE ADAMS(IBND,NS,N,Y,DY,H,T,REPS,NSTP,IJ,IJM1,IJM2,
     1                 IJM3,IJM4,IER,WK)
      use b2mod_types
      IMPLICIT none
      integer n, nit, nstp, ij, ijm1, ijm2, ijm3, ijm4, ier, ibnd, ns
      external lines2c
      real (kind=R8) :: eps, cfac, reps, h, t
      PARAMETER(EPS=1.e-30_R8,NIT=10,CFAC=0.20_R8)
      real (kind=R8) :: Y(N,7),DY(N,7),WK(N,2)

      integer j, k
      real (kind=R8) :: T1, T2, R1, RERR

      IER=0
      DO 2000 J=1,N
        T1=Y(J,IJM1) + H*(55.0_R8*DY(J,IJM1) - 59.0_R8*DY(J,IJM2) +
     &                37.0_R8*DY(J,IJM3) - 9.0_R8*DY(J,IJM4))/24.0_R8
        Y(J,IJ)=T1 + 251.0_R8*(Y(J,IJM1) - WK(J,2))/270.0_R8
      WK(J,1)=T1
 2000 CONTINUE

      DO 3000 J=1,NIT
        NSTP=NSTP+1
        CALL lines2c(IBND,NS,T,N,Y(1,IJ),DY(1,IJ))
        RERR=0.0_R8

        DO 2600 K=1,N
          T2=Y(K,IJM1) + H*(9.0_R8*DY(K,IJ) + 19.0_R8*DY(K,IJM1) -
     &        5.0_R8*DY(K,IJM2) + DY(K,IJM3))/24.0_R8
          R1=ABS(T2)+ABS(T2-Y(K,IJM1))+EPS
          T1=ABS((Y(K,IJ)-T2)/R1)
          IF(T1.GT.RERR) RERR=T1
          Y(K,IJ)=T2
 2600   CONTINUE
        IF(RERR.LT.CFAC*REPS) RETURN
 3000 CONTINUE
      IER=125
      RETURN
      END SUBROUTINE ADAMS
c_______________________________________________________________________
      SUBROUTINE RKM(IBND,NS,N,Y,YT,DY,H,T0,TN,REPS,NSTEP,NMAX,IER,
     . WK)
      use b2mod_types
      implicit none
      integer n, nmx, nstep, nmax, ier, ibnd, ns
      real (kind=R8) :: sfac, eps, e1, e2, ce
      real (kind=R8) :: H, T0, TN, REPS
      PARAMETER(NMX=10000,SFAC=0.9_R8,EPS=1.0e-30_R8)
      PARAMETER(E1=0.20_R8,E2=0.25_R8,CE=15.0_R8)
      EXTERNAL lines2c
      real (kind=R8) :: Y(N),YT(N),DY(N),WK(N,5)

      integer I
      real (kind=R8) :: TSTEP, H2, T1, RERR, R1, R2

      IF(N.LE.0) THEN
        IER=111
        RETURN
      ENDIF
      IF(NMAX.LE.0) NMAX=NMX
      NSTEP=0
      TSTEP=TN-T0
      IER=0
      IF(TSTEP.EQ.0.0_R8) then
        DO I = 1,N
          YT(I) = Y(I)
        enddo
        RETURN
      endif
      IF(H.EQ.0.0.OR.ABS(H).GT.ABS(TSTEP)) H=TSTEP
      IF((H.LT.0.0_R8).EQV.(TN.GT.T0)) H=-H
      CALL lines2c(IBND,NS,T0,N,Y,DY)

 1000 NSTEP=NSTEP+1
      H2=H/2.0_R8
      CALL RK4(IBND,NS,N,T0,Y,DY,H2,WK(1,4),WK)
      T1=T0+H2
      CALL lines2c(IBND,NS,T1,N,WK(1,4),WK(1,5))
      CALL RK4(IBND,NS,N,T1,WK(1,4),WK(1,5),H2,WK(1,3),WK)
      CALL RK4(IBND,NS,N,T0,Y,DY,H,WK(1,4),WK)

      RERR=0
      DO 2000 I=1,N
        R2=ABS(Y(I))+ABS(WK(I,3)-Y(I))+EPS
        R1=ABS(WK(I,4)-WK(I,3))/R2
        IF(R1.GT.RERR) RERR=R1
 2000 CONTINUE
      RERR=ABS(RERR*TSTEP/(CE*H))

      IF(T1.EQ.T0) THEN
        IER=123
        RETURN
      ENDIF

      IF(RERR.LT.REPS) THEN
        T0=T0+H
        DO 3000 I=1,N
          Y(I)=WK(I,3)-(WK(I,4)-WK(I,3))/CE
 3000   CONTINUE
        CALL lines2c(IBND,NS,T0,N,Y,DY)
        IF(ABS((TN-T0)/TSTEP).LT.REPS) RETURN

        IF(RERR.EQ.0.0_R8) THEN
          H=2.0_R8*H
        ELSE
          H=SFAC*H*(REPS/RERR)**E1
        ENDIF
        IF((T0+H.GT.TN).EQV.(H.GT.0.0_R8)) H=TN-T0
      ELSE
        H=SFAC*H*(REPS/RERR)**E2
      ENDIF

      IF(ABS(H/TSTEP).LT.REPS) THEN
        IER=122
        RETURN
      ENDIF
      IF(NSTEP.LT.NMAX) GO TO 1000
      IER=133
      RETURN
      END SUBROUTINE RKM
