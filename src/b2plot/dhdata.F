C BREIUP.DHDATA     KOPIE VON RHW:BREIASPL.DHDATA, ANGEPASST AN "DISPEK"
C
       SUBROUTINE DHDATA( YEAR   , TITLF  , IFAIL
     &                  , IZ0    , IZ1    , ICLASS , ITMAX  , IEVCUT
     &                  , ITDIMD , ITMAXD , IDMAXD , IZMAXD
     &                  , DTEV   , DDENS
     &                  , DTEVD  , DDENSD , DRCOFD , ZDATA
     &                  , DRCOFI
     &                  )
       use b2mod_types
       use b2mod_elements
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
       IMPLICIT NONE
C
C-----------------------------------------------------------------------
C
C PURPOSE : TO EXTRACT 'SANC0' COLLISIONAL DIELECTRONIC DATA
C
C NOTE    : THE SOURCE DATA IS CONTAINED AS SEQUENTIAL DATASETS
C           AS FOLLOWS:
C
C         DS-NAME                   ICLASS     RATENKOEFF. FUER
C
C  (1)  ACD<YR>#<IEL).DATA             1       REKOMBINATION
C  (2)  SCD<YR>#<IEL>.DATA             2       IONISATION
C  (3)  CCD<YR>#<IEL>.DATA             3       CHARGE EXCHANGE
C  (4)  PRB<YR>#<IEL>.EV<CUT>.DATA     4       BREMSSTRAHLUNG
C  (5)  PLT<YR>#<IEL>.EV<CUT>.DATA     5       LINIENSTRAHLUNG
C  (6)  PRC<YR>#<IEL>.EV<CUT>.DATA     6       CX-STRAHLUNG?
C  (7)  PLS<YR>#<IEL>.DATA             7       DIAGNOSTIK LINIE
C
C       WHERE, <YR>  = TWO INTEGERS  FOR THE YEAR SELECTED
C              <IEL> = ELEMENT NAME
C              <CUT> = ENERGY CUT-OFF (EV)
C
C              IF <CUT> = 0 THEN .EV<CUT> IS DELETED FROM ABOVE FILES.
C
C INPUT  : (C*2)  YEAR      = YEAR OF DATA
C          (I*4)  IZ0       = NUCLEAR CHARGE
C          (I*4)  IZ1       = MINIMUM ION CHARGE + 1
C          (I*4)  ICLASS    = CLASS OF DATA (1 - 6)
C          (I*4)  ITMAX     = NUMBER OF ( DTEV() , DDENS() ) PAIRS
C          (I*4)  IEVCUT    = ENERGY CUT-OFF (EV)
C          (R*8)  DTEV()    = DLOG10(ELECTRON TEMPERATURES (EV))
C          (R*8)  DDENS()   = DLOG10(ELECTRON DENSITIES (CM-3))
C
C OUTPUT : (C*80) TITLF     = INFORMATION STRING
C          (I*4)  ITDIMD    = MAXIMUM NUMBER OF DATA TEMP & DENS
C          (I*4)  ITMAXD    = NUMBER OF DATA DTEVD()
C          (I*4)  IDMAXD    = NUMBER OF DATA DDENS()
C          (I*4)  IZMAXD    = NUMBER OF DATA ZDATA()
C          (I*4)  ITDIMD    = MAXIMUM NUMBER OF DATA TEMP & DENS
C          (I*4)  ZDATA()   = Z1 CHARGES IN DATASET
C          (I*4)  IFAIL     = 0    IF ROUTINE SUCCESSFUL
C                           = 1    IF ROUTINE OPEN STATEMENT FAILED
C          (R*8)  DTEVD()   = DLOG10(DATA ELECTRON TEMPERATURES (EV))
C          (R*8)  DDENSD()  = DLOG10(DATA ELECTRON DENSITIES (CM-3))
C          (R*8)  DRCOFD()  = DLOG10(DATA RATE COEFFICIENTS (CM-3/S))
C          (R*8)  DRCOFI()  = INTERPOLATION OF DRCOFD(,,) FOR
C                             DTEV() & DDENS()
C
C PROGRAM: (C*2)  SEQUA()   = ION NAMES FOR A PARTICULAR IZ0
C          (C*30) DSNAME    = FILE NAME ( SEE ABOVE TYPES )
C          (C*80) STRING    = GENERAL VARIABLE
C          (C*40) BLANK     = BLANK STRING
C          (C*2)  YEARSV    = LAST YEAR USED IN THIS ROUTINE
C          (I*4)  IREAD     = INPUT STREAM FOR OPEN STATEMENT
C          (I*4)  IZ0SV     = LAST IZ0 USED IN THIS ROUTINE
C          (I*4)  ICLSV     = LAST ICLASS USED IN THIS ROUTINE
C          (I*4)  INDXZ1    = LOCATION OF IZ1 IN ZDATA()
C          (I*4)  LCK       = MUST BE GREATER THAN 'ITMAXD' & 'IDMAXD'
C                             & 'ITMAX' - ARRAY SIZE FOR SPLINE CALCS.
C          (R*8)  A()       = GENERAL ARRAY
C          (R*8)  DRCOF0(,) = INTERPOLATION OF DRCOFD(,,) W.R.T DTEV()
C          (L*8)  LEXIST    = TRUE --- FILE TO OPEN EXISTS ELSE NOT
C
C PE BRIDEN = ADDED VARIABLES (14/01/91)
C
C          (I*4)  L1      = PARAMETER = 1
C          (I*4)  IOPT    = DEFINES THE BOUNDARY DERIVATIVES FOR THE
C                             SPLINE ROUTINE 'XXSPLN', SEE 'XXSPLN'.
C
C          (L*4)  LSETX   = .TRUE.  => SET UP SPLINE PARAMETERS RELATING
C                                      TO X-AXIS.
C                           .FALSE. => DO NOT SET UP SPLINE PARAMETERS
C                                      RELATING TO X-AXIS.
C                                      (I.E. THEY WERE SET IN A PREVIOUS
C                                            CALL )
C                           (VALUE SET TO .FALSE. BY 'XXSPLN')
C
C
C          (R*8)  DY()    = SPLINE INTERPOLATED DERIVATIVES
C
C          (R*8 ADAS FUNCTION - 'R8FUN1' ( X -> X) )
C
C AUTHOR : JAMES SPENCE (TESSELLA SUPPORT SERVICES PLC)
C          K1/0/80
C          JET  EXT. 4866
C
C DATE   : 22/02/90
C
C DATE   : 25/03/06 XP BONNIN - REVISION: SEQUA DEDUCED FROM ELEMENTS ARRAY
C
C DATE   : 24/01/03 XP BONNIN - REVISION: SEQUA(51-54) ADDED ('Sb','Te','I ','Xe')
C
C DATE   : 21/08/90 PE BRIDEN - REVISION: SEQUA(43) CHANGED ('TE'->'TC')
C
C DATE   : 08/10/90 PE BRIDEN - REVISION: RENAMED SUBROUTINE
C
C DATE   : 12/11/90 PE BRIDEN - CORRECTION: MOVE THE SETTING OF 'INDXZ1'
C                                           TO AFTER THE  '20 CONTINUE'
C                                           STATEMENT.   ALSO SAVE  THE
C                                           VALUE OF 'IZ1MIN'.
C
C DATE   : 14/01/91 PE BRIDEN - ADAS91:     CALLS TO NAG SPLINE ROUTINES
C                                           'E01BAF' & 'E02BBF' REPLACED
C                                           BY  CALLS   TO  ADAS  SPLINE
C                                           ROUTINE 'XXSPLN'.
C
C DATE   : 25/06/91 PE BRIDEN - CORRECTION: CHANGED FOLLOWING DIMENSION:
C                                            'DIMENSION DRCOFI(ITDIMD)'
C                                           TO
C                                            'DIMENSION DRCOFI(ITMAX)'
C
C
C DATE   : 25/06/91 HP SUMMERS - REVISION: RENAMED FROM D2DATA TO DHDATA
C                                          SET IOPT=4 FOR IONBAL
C-----------------------------------------------------------------------
C
       INTEGER   L1, IREAD, LCK, IDMAXD_LOC
C
       PARAMETER ( L1    =  1 )
       PARAMETER ( IREAD = 12 , LCK = 5000 , IDMAXD_LOC = 40 )
C
       INTEGER   IOPT
C
       integer IZMAXD, IDMAXD, ITMAXD, IEVCUT, ICLASS, IZ1, IZ0, IFAIL
       integer ITMAX, ITDIMD, IZ1MIN, IZ0SV, ICLSV, IEVSV, LEN
       integer lnblnk, LEVMIN, IZ1MAX, IZMAX, I, ID, IT, IZ, INDXZ1
c
       real (kind=R8) :: DTEV(ITMAX)   , DDENS(ITMAX)
       real (kind=R8) :: DTEVD(ITDIMD) , DDENSD(ITDIMD) , ZDATA(ITDIMD)
       real (kind=R8) :: DRCOFD(ITDIMD,ITDIMD,ITDIMD)
       real (kind=R8) :: DRCOFI(ITMAX)
C
       real (kind=R8) :: A(IDMAXD_LOC)
       real (kind=R8) :: DY(IDMAXD_LOC)
       real (kind=R8) :: DRCOF0(LCK,IDMAXD_LOC)
C
       CHARACTER YEAR*2  , YEARSV*2  , TITLF*80  , EVCUT*6
       CHARACTER SEQUA*2 , DSNAME*256, STRING*80 , BLANKS*40
C
       LOGICAL   LEXIST  , LSETX
C
       character*256 solpstop
       real (kind=R8) :: R8FUN1
       EXTERNAL  R8FUN1
#ifdef USE_PXFGETENV
       integer lenval, ierror
#endif
C
       SAVE      IZ1MIN
C
C-----------------------------------------------------------------------
C
       DATA BLANKS/'                                        '/
C
       DATA YEARSV/'  '/
       DATA IZ0SV /0   /
       DATA ICLSV /0   /
       DATA IEVSV /0   /
C
C-----------------------------------------------------------------------
C
         IFAIL = 0
         if(iz0.le.0) return
C
C-----PARAMETRIC DIMENSIONS CHECKED-------------------------------------
C

C
c-dpc-1997.09.29      IF(YEAR.EQ.YEARSV.AND.IZ0.EQ.IZ0SV.AND.ICLASS.EQ.ICLSV.AND.
c-dpc-1997.09.29     &    IEVCUT.EQ.IEVSV) GOTO 20
       YEARSV = YEAR
       IZ0SV  = IZ0
       ICLSV  = ICLASS
       IEVSV  = IEVCUT
C
C------ENERGY CUTOFF----------------------------------------------------
C
       IF( IEVCUT.GT.0 ) THEN
         WRITE(EVCUT,1050) IEVCUT
         DO 50 I = 1 , 6
           IF( EVCUT(I:I).EQ.' ' ) LEVMIN = I
  50     CONTINUE
         LEVMIN = LEVMIN + 1
       END IF
C
C------FILE NAME--------------------------------------------------------
C
       SEQUA = elements(IZ0)
       call lower_case(SEQUA)
   10  IF (SEQUA(2:2).EQ.' ') THEN
         LEN=1
       ELSE
         LEN=2
       ENDIF
C
#ifdef NO_GETENV
       solpstop=' '
#else
#ifdef USE_PXFGETENV
       CALL PXFGETENV ('SOLPSTOP', 0, solpstop, lenval, ierror)
#else
       call getenv ('SOLPSTOP', solpstop)
#endif
#endif
       if(solpstop.eq.' ') then
         solpstop='/afs/ipp-garching.mpg.de/u/dpc/solps4.0'
         write(*,*) 'DHDATA: no SOLPSTOP set, using: ',
     1    solpstop(1:lnblnk(solpstop))
       endif
       IF     (ICLASS.EQ.1) THEN
           DSNAME='/JETSHP.ACD'//YEAR//'#'//SEQUA(1:LEN)//'.DATA'
       ELSEIF (ICLASS.EQ.2) THEN
           DSNAME='/JETSHP.SCD'//YEAR//'#'//SEQUA(1:LEN)//'.DATA'
       ELSEIF (ICLASS.EQ.3) THEN
           DSNAME='/JETSHP.CCD'//YEAR//'#'//SEQUA(1:LEN)//'.DATA'
       ELSEIF (ICLASS.EQ.4) THEN
         IF( IEVCUT.EQ.0 ) THEN
           DSNAME='/JETSHP.PRB'//YEAR//'#'//SEQUA(1:LEN)//'.DATA'
         ELSE
           DSNAME='/JETSHP.PRB'//YEAR//'#'//SEQUA(1:LEN)//
     &             '.EV'//EVCUT(LEVMIN:6)//'.DATA'
         END IF
       ELSEIF (ICLASS.EQ.5) THEN
         IF( IEVCUT.EQ.0 ) THEN
           DSNAME='/JETSHP.PLT'//YEAR//'#'//SEQUA(1:LEN)//'.DATA'
         ELSE
           DSNAME='/JETSHP.PLT'//YEAR//'#'//SEQUA(1:LEN)//
     &             '.EV'//EVCUT(LEVMIN:6)//'.DATA'
         END IF
       ELSEIF (ICLASS.EQ.6) THEN
         IF( IEVCUT.EQ.0 ) THEN
           DSNAME='/JETSHP.PRC'//YEAR//'#'//SEQUA(1:LEN)//'.DATA'
         ELSE
           DSNAME='/JETSHP.PRC'//YEAR//'#'//SEQUA(1:LEN)//
     &             '.EV'//EVCUT(LEVMIN:6)//'.DATA'
         END IF
       ELSEIF (ICLASS.EQ.7) THEN
C.....     DSNAME='/JETSHP.PLS'//YEAR//'#'//SEQUA(1:LEN)//'.DATA'
           DSNAME='pls'//YEAR//SEQUA(1:LEN)//'.data'
       ENDIF
C
       dsname=solpstop(1:lnblnk(solpstop)) // '/data/summers/' // dsname

       TITLF=BLANKS//BLANKS
CCC    WRITE(TITLF,1000) DSNAME
C
C------DOES FILE TO BE OPEN EXIST OR NOT--------------------------------
C
       INQUIRE(FILE=DSNAME,EXIST=LEXIST)
       IF( .NOT.LEXIST ) THEN
         IFAIL  = 1
         YEARSV = '   '
         IZ0SV  = 0
         ICLSV  = 0
         GOTO 9999
       END IF
C
C------READ FILE # IREAD------------------------------------------------
C
c       OPEN(UNIT=IREAD,FILE=DSNAME,ACTION='READ')
       OPEN(UNIT=IREAD,FILE=DSNAME)
C
       READ(IREAD,1010) IZMAX , IDMAXD , ITMAXD , IZ1MIN , IZ1MAX
C
       IF (IDMAXD_LOC.LT.MAX(IDMAXD,ITMAXD)) THEN
         WRITE(6,1020) 'LOCALS ARRAYS TOO SMALL'
         WRITE(6,1020)
     .    'IDMAXD_LOC MUST BE LARGER THAN ITMAXD AND IDMAXD'
         WRITE(6,1100) IDMAXD_LOC,IDMAXD,ITMAXD
 1100    FORMAT('IDMAXD_LOC = ',I4,' ITMAXD = ',I4,' IDMAXD = ',I4)
         IFAIL=2
         GOTO 9999
       ENDIF
      IF (LCK.LT.MAX(IDMAXD,ITMAXD,ITMAX)) THEN
        WRITE(6,1020) 'LOCALS ARRAYS TOO SMALL'
        WRITE(6,1020)
     .   'LCK MUST BE LARGER THAN ITMAXD, IDMAXD AND ITMAX'
        WRITE(6,1110) LCK,ITMAXD,IDMAXD,ITMAX
 1110   FORMAT('LCK = ',I4,' ITMAXD = ',I4,' IDMAXD = ',I4,
     .   ' ITMAX = ',I4)
        IFAIL=2
        GOTO 9999
      ENDIF
      IF (IDMAXD.GT.ITDIMD .OR. ITMAXD.GT.ITDIMD)   THEN
        if(IDMAXD.GT.ITDIMD)
     1   write(*,*) 'IDMAXD,ITDIMD ',IDMAXD,ITDIMD
        if(ITMAXD.GT.ITDIMD)
     1   write(*,*) 'ITMAXD,ITDIMD ',ITMAXD,ITDIMD
         IFAIL = 2
        GO TO 9999
      ENDIF
      IF (ITDIMD.LT.(IZ1MAX-IZ1MIN+1)) THEN
         WRITE(6,1020) 'LOCALS ARRAYS TOO SMALL'
         WRITE(6,1020)
     .      'ITDIMD MUST BE LARGER THAN IZ1MAX-IZ1MIN+1'
         WRITE(6,1120) ITDIMD,IZ1MAX,IZ1MIN
1120     FORMAT('ITDIMD = ',I4,' IZ1MAX = ',I4,' IZ1MIN = ',I4)
         IFAIL=2
         GOTO 9999
      ENDIF
C
       READ(IREAD,1020) STRING
       READ(IREAD,1040) ( DDENSD(ID) , ID = 1 , IDMAXD )
       READ(IREAD,1040) ( DTEVD(IT)  , IT = 1 , ITMAXD )
C
       IZMAXD = 0
       DO 150 IZ = IZ1MIN , IZ1MAX
          IZMAXD = IZMAXD + 1
          ZDATA(IZMAXD) = IZ
CCCCCCC   IF( IZ .EQ. IZ1 ) INDXZ1 = IZ
          READ(IREAD,1020)STRING
          DO 100 IT = 1 , ITMAXD
             READ(IREAD,1040) ( DRCOFD(IZMAXD,IT,ID) , ID = 1 , IDMAXD )
  100     CONTINUE
  150  CONTINUE
C
       CLOSE(IREAD)
C
C------INTERPOLATE USING SPLINES (NAG ALGORITHM)------------------------
C
   20  CONTINUE
C
C------PE BRIDEN - CORRECTION 12/11/90 - SET INDXZ1 AFTER '20 CONTINUE'-
C
       INDXZ1 = IZ1 - IZ1MIN + 1
       if(indxz1.le.0.0) then
         write(*,*) 'DHDATA: indxz1 = ',indxz1,' <= 0'
         ifail=1
         return
       elseif(indxz1.gt.IZMAXD) then
         write(*,*) 'DHDATA: indxz1 = ',indxz1,' > izmaxd = ',izmaxd
         ifail=1
         return
       endif
C
C-----------------------------------------------------------------------
C
C
C>>>>>>INTERPOLATE DRCOFD(,,,) W.R.T TEMPERATURE
C
      LSETX = .TRUE.
      IOPT  = 4
C
      DO ID = 1 , IDMAXD
C
        DO IT = 1 , ITMAXD
          A(IT) = DRCOFD(INDXZ1,IT,ID)
        END DO
C
       CALL XXSPLN( LSETX  , IOPT  , R8FUN1       ,
     &              ITMAXD , DTEVD , A            ,
     &              ITMAX  , DTEV  , DRCOF0(1,ID) ,
     &              DY )
C
      END DO
C
C>>>>>>INTERPOLATE ABOVE RESULT W.R.T DENSITY
C
      LSETX = .TRUE.
      IOPT  = 4
C
      DO IT = 1 , ITMAX
C
        DO ID = 1 , IDMAXD
          A(ID) = DRCOF0(IT,ID)
        END DO
C
        CALL XXSPLN( LSETX  , IOPT      , R8FUN1     ,
     &               IDMAXD , DDENSD    , A          ,
     &               L1     , DDENS(IT) , DRCOFI(IT) ,
     &               DY )
C
      END DO
C
C-----------------------------------------------------------------------
C
 1000  FORMAT('FILE = ',1A30)
 1010  FORMAT(5I5)
 1020  FORMAT(1A80)
 1040  FORMAT(8F10.5)
 1050  FORMAT(I6)
C
C-----------------------------------------------------------------------
C
 9999  RETURN
       END SUBROUTINE DHDATA
*COPY R8FUN1
      FUNCTION R8FUN1 ( Z )
      use b2mod_types
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  **************** FORTRAN77 REAL FUNCTION: R8FUN1 ******************
C
C  FUNCTION: R8FUN1 =  Z
C
C  CALLING PROGRAM: GENERAL USE
C
C  FUNCTION:
C          (R*8)  R8FUN1  = FUNCTION NAME
C          (R*8)  Z       = INPUT VALUE
C
C AUTHOR:   PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C           K1/0/81
C           JET EXT. 4569
C
C DATE:     13/08/90
C
C-----------------------------------------------------------------------
      real (kind=R8) :: R8FUN1 , Z
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      R8FUN1=Z
C-----------------------------------------------------------------------
      RETURN
      END FUNCTION R8FUN1
*COPY XXSPLN
      SUBROUTINE XXSPLN( LSETX , IOPT  , FINTX ,
     &                   NIN   , XIN   , YIN   ,
     &                   NOUT  , XOUT  , YOUT  ,
     &                   DY
     &                 )
      use b2mod_types
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXSPLN *********************
C
C  PURPOSE:           TO INTERPOLATE/EXTRAPOLATE USING CUBIC SPLINES
C                     (IF IOPT < 0 NO EXTRAPOLATION TAKES PLACE = VALUES
C                      SET TO ZERO).
C
C  CALLING PROGRAMS:  GENERAL USE
C
C  SUBROUTINE:
C
C  I/O   : (L*4)  LSETX   = .TRUE.  => SET UP SPLINE PARAMETERS RELATING
C                                      TO 'XIN' AXIS.
C                           .FALSE. => DO NOT SET UP SPLINE PARAMETERS
C                                      RELATING TO 'XIN' AXIS.
C                                      (I.E. THEY WERE SET IN A PREVIOUS
C                                            CALL )
C                           ( 'LSETX' IS ALWAYS RETURN AS '.FALSE.'  ON
C                             RETURN FROM THE SUBROUTINE ).
C                           ** IMPORTANT: SEE NOTES BELOW ON 'LSETX' **
C  INPUT : (I*4)  IOPT    = SPLINE END CONDITIONS/EXTRAPOLATION CONTROL
C                           SWITCH - SEE NOTES BELOW
C                           I.E. DEFINES THE BOUNDARY DERIVATIVES.
C                                (VALID VALUES = 0, 1, 2, 3, 4)
C                           IF IOPT < 0 THEN NO EXTRAPOLATION TAKES
C                           - ANY VALUES REQUIRING EXTRAPOLATION WILL BE
C                             SET TO ZERO (END CONDITIONS AS FOR IOPT=0)
C  INPUT : (R*8)  FINTX   = INTERPOLATING X-COORDINATE TRANSFORMATION.
C                           EXTERNAL FUNCTION (SEE ROUTINES BELOW)
C
C  INPUT : (I*4)  NIN     = NUMBER OF KNOTS
C  INPUT : (R*8)  XIN()   = X-VALUES OF KNOTS
C  INPUT : (R*8)  YIN()   = Y-VALUES OF KNOTS
C
C  INPUT : (I*4)  NOUT    = NUMBER OF OUTPUT VALUES TO BE INTERPOLATED
C                           EXTRAPOLATED.
C  INPUT : (R*8)  XOUT()  = X-VALUES AT WHICH INTERPOLATION/EXTRAPOLA-
C                           TION REQUIRED
C  OUTPUT: (R*8)  YOUT()  = INTERPOLATED/EXTRAPOLATED Y-VALUES FOR
C                           REQUESTED 'XOUT()' VALUES.
C
C  OUTPUT: (R*8)  DY()    = INTERPOLATED DERIVATIVES
C
C          (I*4)  NKNOTS  = PARAMETER = MAXIMUM  NUMBER OF KNOTS ALLOWED
C          (I*4)  NIOPT   = PARAMETER = MAXIMUM  VALUE OF IOPT ALLOWED
C
C          (I*4)  I       = GENERAL ARRAY USE
C          (I*4)  K       = INDEX OF 'XOUT()' VALUE FOR INTERPOLATION/
C                           EXTRAPOLATION.
C          (I*4)  NIN0    = 'NIN' - 1
C          (I*4)  INTER   = INDEX OF CLOSEST/NEXT HIGHEST VALUE OF
C                           'XIN()' TO THE VALUE OF 'XOUT()' BEING
C                           INTERPOLATED/EXTRAPOLATED.
C          (I*4)  NOPT    = VALUE OF  'IOPT' USED IN  CALCULATING  END-
C                           CONDITIONS   FOR  STORED  'X-VALUE'  SPLINE
C                           PARAMETERS.   (NOTE:  IF  'IOPT < 0',  THEN
C                           'NOPT = 0'.) - I.E. 'NOPT = MAX( 0, IOPT )'.
C
C          (R*8)  XK      = VALUE OF 'XOUT(K)' BEING INTERPOLATED/
C                           EXTRAPOLATED
C          (R*8)  XKK     = TRANSFORMED VALUE OF 'XOUT(K)' BEING
C                           INTERPOLATED/EXTRAPOLATED.
C          (R*8)  T1      = INVERSE OF SEPARATION OF KNOTS EITHER
C                           SIDE OF CURRENT KNOT.
C          (R*8)  T2      = (CURRENT KNOT POSITION TO NEXT HIGHEST KNOT
C                            POSITION) DIVIDED BY 'T1'
C          (R*8)  T3      = (CURRENT KNOT POSITION TO NEXT LOWEST  KNOT
C                            POSITION) DIVIDED BY 'T1'
C          (R*8)  T4      = INTERPOLATION FACTOR FOR CURRENT KNOT
C          (R*8)  DL1     = (REQUESTED 'XOUT()' VALUE TO NEXT HIGHEST
C                            KNOT POSITION) DIVIDED BY SEPERATION OF
C                            KNOTS EITHER SIDE OF 'XOUT(K)'.
C          (R*8)  DL2     = (REQUESTED 'XOUT()' VALUE TO NEXT LOWEST
C                            KNOT POSITION) DIVIDED BY SEPERATION OF
C                            KNOTS EITHER SIDE OF 'XOUT(K)'.
C          (R*8)  DL2     = (REQUESTED 'XOUT()' VALUE TO NEXT LOWEST
C          (R*8)  DL3     =  SEPERATION OF KNOTS EITHER SIDE OF
C                            'XOUT(K)' * 'DL1' * 'DL2'.
C
C          (L*4)  LEXTRP  = .TRUE.  => 'EXTRAPOLATION SWITCHED ON'.
C                           .FALSE. => 'EXTRAPOLATION SWITCHED OFF'.
C
C          (R*8)  QVAL()  = VALUE OF 'Q(1)'   : FUNCTION OF 'NOPT'
C          (R*8)  D2VAL() = VALUE OF 'D2(1)'  : FUNCTION OF 'NOPT'
C          (R*8)  D3VAL() = VALUE OF 'D3(1)'  : FUNCTION OF 'NOPT'
C          (R*8)  UVAL() =  VALUE OF 'U(NIN)' : FUNCTION OF 'NOPT'
C          (R*8)  AGRL() =  POLYNOMIAL CONSTANTS FOR CUBIC SPLINE FOR
C                           GIVEN 'XOUT(K)' VALUE.
C          (R*8)  X()    =  TRANSFORMED VALUES OF 'XIN()'
C          (R*8)  H()    =  SEPERATION, ALONG X-AXIS, OF KNOT FROM NEXT
C                           HIGHEST KNOT.
C          (R*8)  Q()    =  SECOND DERIVATIVE FOR KNOT
C          (R*8)  U()    =  TEMPORARY STORAGE OF DECOMPOSED FACTORS
C          (R*8)  DELY() =  SEPERATION, ALONG Y-AXIS, OF KNOT FROM NEXT
C                           HIGHEST KNOT.
C          (R*8)  D1()   =  MULTIPLICATION FACTOR USED IN CALCULATING
C                           'U()'.
C          (R*8)  D2()   =  MULTIPLICATION FACTOR USED IN CALCULATING
C                           'U()'.
C          (R*8)  D3()   =  MULTIPLICATION FACTOR USED IN CALCULATING
C                           'U()'.
C
C          (L*4)  LUVAL()=  .TRUE. => VALUE OF 'UVAL()' REFERS TO RATE
C                                     OF CHANGE OF SLOPE AT FINAL POINT.
C                           .FALSE.=> VALUE OF 'UVAL()' REFERS TO FINAL
C                                     SLOPE
C                            FUNCTION OF 'NOPT'
C
C NOTES: 'LSETX': SET TO .TRUE. ON ENTRY IF A NEW 'XIN' ARRAY IS BEING
C                 USED.  IF THE 'XIN' AXIS IS THE SAME FOR A NUMBER OF
C                 CALLS THEN DO NOT RESET 'LSETX'  -  THIS  SUBROUTINE
C                 SETS IT TO .FALSE. FOR YOU.   IF THE VALUE OF 'NOPT'
C                 IS CHANGED BETWEEN CALLS THEN THE VALUE  OF  'LSETX'
C                 ON  ENTRY IS TAKEN AS BEING EQUAL TO .TRUE. .
C
C                 THEREFORE 'LSETX' NEED ONLY BE SET TO .TRUE. ON ENTRY
C                 IF EITHER IT IS ITS FIRST CALL OR IF ANY ONE  OF  THE
C                 FOLLOWING VALUES HAS CHANGED:
C
C                 'NIN' , 'FINTX' , 'XIN(I), I=1,NIN'
C
C                 CARE: A VARIABLE MUST BE USED FOR 'LSETX', A CONSTANT,
C                       I.E.  .TRUE. ,  CANNOT BE DIRECTLY TYPED AS  AN
C                       ARGUMENT BECAUSE IT WILL BE CHANGED TO  .FALSE.
C                       ON RETURN.
C
C         SPLINE  END CONDITIONS AND EXTRAPOLATION DEPEND ON 'IOPT' AS
C         FOLLOWS:
C
C         --------------------------------------------------------------
C         | IOPT  | NOPT |  DY(1)  DDY(1)  |  DY(N)   DDY(N)  |EXTRAPLN|
C         |-------|------|-----------------|------------------|--------|
C         | < 0   |   0  |    -     0.0    |    -      0.0    |  NO    |
C         |   0   |   0  |    -     0.0    |    -      0.0    |  YES   |
C         |   1   |   1  |    -     0.0    |  -1.5      -     |  YES   |
C         |   2   |   2  |   0.0     -     |   1.0      -     |  YES   |
C         |   3   |   3  |  -0.5     -     |  -1.5      -     |  YES   |
C         |   4   |   4  |   0.0     -     |    -      0.0    |  YES   |
C         --------------------------------------------------------------
C
C            NB. OPTIONS TO BE EXTENDED FOR POWER AND CX APPLICATION
C
C         -------------------------------------------------------------
C          IF ( IOPT.LT.0 ) - NO EXTRAPOLATION TAKES PLACE VALUES SET
C                             TO ZERO (CARE IF LOG OF OUTPUT IS NEEDED).
C          IF ( IOPT.GT.4 ) PROGRAM STOPS
C         -------------------------------------------------------------
C
C          THIS SUBROUTINE IS AN AMENDED  AND STRUCTURED VERSION OF  THE
C          SUBROUTINE  'ESPLINE'  WRITTEN BY  H.P. SUMMERS,   JET   26TH
C          OCTOBER 1989.   IT REMOVES THE COMMON BLOCK  /IONSPL/ ,   THE
C          SWITCHES 'ISW & ISW2' AND ALSO THE CASE FOR THE INTERPOLATION
C          OF CHARGE STATE VALUES.   IT INTRODUCES THE FEATURE  THAT  AN
C          ARRAY OF INPUT  'X-VALUES'  CAN BE  INTERPOLATED/EXTRAPOLATED
C          IN ONE CALL.
C
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C          FINTX      ------    EXTERNAL  REAL  FUNCTION,   USED  TO
C                               TRANSFORM X-COORDINATES.
C
C
C AUTHOR:   PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C           K1/0/81
C           JET EXT. 4569
C
C DATE:     10/08/90 (30/08/90: IOPT = 4 ADDED & 'LUVAL' PARAMETER)
C
C UPDATE:   17/01/91 - PE BRIDEN: ADAS91 - IOPT < 0 ADDED - NO EXTRAPLN.
C                                        - NOPT DEFINITION CHANGED.
C                                        - INTRODUCED 'LEXTRP'.
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      INTEGER    NKNOTS       , NIOPT
C-----------------------------------------------------------------------
      PARAMETER( NKNOTS = 101 , NIOPT = 4 )
C-----------------------------------------------------------------------
      INTEGER    IOPT         , NIN       , NOUT      , NOPT
      INTEGER    I            , NIN0      , K         , INTER
C-----------------------------------------------------------------------
      real (kind=R8) :: FINTX
      real (kind=R8) :: XK           , XKK
      real (kind=R8) :: T1           , T2         , T3       , T4   ,
     &           DL1          , DL2        , DL3
C-----------------------------------------------------------------------
      LOGICAL    LSETX        , LEXTRP
C-----------------------------------------------------------------------
      real (kind=R8) :: XIN(NIN)       , YIN(NIN)       ,
     &           XOUT(NOUT)     , YOUT(NOUT)     ,
     &           DY(NIN)
      real (kind=R8) :: QVAL(0:NIOPT)  , D2VAL(0:NIOPT) ,
     &           D3VAL(0:NIOPT) , UVAL(0:NIOPT)
      REAL (kind=R8) :: AGRL(4)
      REAL (kind=R8) :: X(NKNOTS)      , DELY(NKNOTS)   ,
     &           H(NKNOTS)      , Q(NKNOTS)      , U(NKNOTS)  ,
     &           D1(NKNOTS)     , D2(NKNOTS)     , D3(NKNOTS)
C-----------------------------------------------------------------------
      LOGICAL    LUVAL(0:NIOPT)
C-----------------------------------------------------------------------
      SAVE QVAL , D2VAL , D3VAL , UVAL , LUVAL
      SAVE NOPT , NIN0
      SAVE X    , H     , Q     , D1   , D2   , D3
C-----------------------------------------------------------------------
      DATA  QVAL / -0.5_R8, -0.5_R8,  0.0_R8,  0.0_R8,  0.0_R8 / ,
     &      D2VAL/  1.5_R8,  1.5_R8,  0.0_R8,  0.0_R8,  0.0_R8 / ,
     &      D3VAL/  0.0_R8,  0.0_R8,  0.0_R8, -0.5_R8,  0.0_R8 / ,
     &      UVAL /  0.0_R8, -1.5_R8,  1.0_R8, -1.5_R8,  0.0_R8 /
      DATA  LUVAL/ .TRUE. , .FALSE. , .FALSE. , .FALSE. , .TRUE. /
      DATA  NOPT / -1     /
C-----------------------------------------------------------------------
C***********************************************************************
      IF (NKNOTS.LT.NIN)
     &                 STOP ' XXSPLN ERROR: TOO MANY KNOTS REQUESTED'
      IF (IOPT.GT.NIOPT)
     &                 STOP ' XXSPLN ERROR: INVALID IOPT VALUE ENTERED'
C-----------------------------------------------------------------------
      IF (IOPT.LT.0) THEN
         LEXTRP = .FALSE.
         IF (NOPT.NE.0)    LSETX=.TRUE.
      ELSE
         LEXTRP = .TRUE.
         IF (IOPT.NE.NOPT) LSETX=.TRUE.
      ENDIF
C-----------------------------------------------------------------------
C
         IF (LSETX) THEN
C
C***********************************************************************
C SET UP PARAMETERS RELATING TO 'XIN'-AXIS
C***********************************************************************
C
            NOPT = MAX0( 0 , IOPT )
C-----------------------------------------------------------------------
               DO 1 I=1,NIN
                  X(I)=FINTX(XIN(I))
    1          CONTINUE
C-----------------------------------------------------------------------
            H(2)  = X(2)-X(1)
            Q(1)  = QVAL(NOPT)
            D2(1) = D2VAL(NOPT)/H(2)
            D3(1) = D3VAL(NOPT)
            NIN0=NIN-1
C-----------------------------------------------------------------------
               DO 2 I=2,NIN0
                  H(I+1) = X(I+1) - X(I)
                  T1     = 1.0 / ( H(I+1) + H(I) )
                  T2     = H(I+1) * T1
                  T3     = 1.0 - T2
                  T4     = 1.0 / ( ( T2 * Q(I-1) ) + 2.0 )
                  Q(I)   = -T3 * T4
                  D1(I)  = ( 3.0 * T4 * T2 ) / H(I)
                  D2(I)  = ( 3.0 * T4 * T3 ) / H(I+1)
                  D3(I)  =  T2 * T4
    2          CONTINUE
C-----------------------------------------------------------------------
            T4      = 1.0 / ( Q(NIN0) + 2.0 )
            D1(NIN) = ( 3.0 * T4 ) / H(NIN)
            D3(NIN) = T4
C-----------------------------------------------------------------------
            LSETX=.FALSE.
C***********************************************************************
         ENDIF
C***********************************************************************
C SET UP CUBIC SPLINE DERIVATIVES
C***********************************************************************
      DELY(2) = YIN(2) - YIN(1)
      U(1)    = ( D2(1) * DELY(2) ) + D3(1)
C-----------------------------------------------------------------------
         DO 3 I=2,NIN0
            DELY(I+1) = YIN(I+1) - YIN(I)
            U(I) = (D1(I)*DELY(I)) + (D2(I)*DELY(I+1)) - (D3(I)*U(I-1))
    3    CONTINUE
C-----------------------------------------------------------------------
         IF (LUVAL(NOPT)) THEN
            U(NIN) = ( D1(NIN)*DELY(NIN) ) - ( D3(NIN)*U(NIN0) )
         ELSE
            U(NIN) = UVAL(NOPT)
         ENDIF
C-----------------------------------------------------------------------
      DY(NIN) = U(NIN)
C-----------------------------------------------------------------------
         DO 4 I=NIN0,1,-1
            DY(I)  = ( Q(I) * DY(I+1) ) + U(I)
    4    CONTINUE
C***********************************************************************
C SET UP PARAMETERS RELATING TO THE REQUESTED 'XOUT' ARRAY VALUES
C***********************************************************************
         DO 5 K=1,NOUT
            XK  = XOUT(K)
            XKK = FINTX(XK)
C-----------------------------------------------------------------------
C EXTRAPOLATE: HIGH 'XOUT' VALUE - IF EXTRAPOLATION SWITCHED ON
C-----------------------------------------------------------------------
               IF     (XK.GE.XIN(NIN)) THEN
                  INTER   = NIN
                  AGRL(1) = 0.0
                  AGRL(3) = 0.0
                     IF (XK.EQ.XIN(NIN)) THEN
                        AGRL(2) = 0.0
                        AGRL(4) = 1.0
                     ELSEIF (LEXTRP) THEN
                        AGRL(2) = XKK-X(NIN)
                        AGRL(4) = 1.0
                     ELSE
                        AGRL(2) = 0.0
                        AGRL(4) = 0.0
                     ENDIF
C-----------------------------------------------------------------------
C EXTRAPOLATE: LOW 'XOUT' VALUE  - IF EXTRAPOLATION SWITCHED ON
C-----------------------------------------------------------------------
               ELSEIF (XK.LE.XIN(1))   THEN
                  INTER   = 2
                  AGRL(2) = 0.0
                  AGRL(4) = 0.0
                     IF (XK.EQ.XIN(1)) THEN
                        AGRL(1) = 0.0
                        AGRL(3) = 1.0
                     ELSEIF (LEXTRP) THEN
                        AGRL(1) = XKK-X(1)
                        AGRL(3) = 1.0
                     ELSE
                        AGRL(1) = 0.0
                        AGRL(3) = 0.0
                     ENDIF
C-----------------------------------------------------------------------
C INTERPOLATE:
C-----------------------------------------------------------------------
               ELSE
                     DO 6 I=NIN,1,-1
                        IF (XK.LT.XIN(I)) INTER=I
    6                CONTINUE
                  DL1     =  ( X(INTER) - XKK    ) / H(INTER)
                  DL2     =  1.0 - DL1
                  DL3     =  H(INTER) * DL1 * DL2
                  AGRL(1) =  DL1 * DL3
                  AGRL(2) = -DL2 * DL3
                  AGRL(3) =  DL1 * DL1 * ( 1.0 + DL2 + DL2 )
                  AGRL(4) =  1.0 - AGRL(3)
               ENDIF
C-----------------------------------------------------------------------
C EVALUATE 'YOUT'-VALUE FOR REQUESTED 'XOUT'-VALUE
C-----------------------------------------------------------------------
            YOUT(K) = ( AGRL(1)* DY(INTER-1) ) + ( AGRL(2)* DY(INTER) )
     &              + ( AGRL(3)*YIN(INTER-1) ) + ( AGRL(4)*YIN(INTER) )
C-----------------------------------------------------------------------
    5    CONTINUE
C***********************************************************************
      RETURN
      END SUBROUTINE XXSPLN
c
