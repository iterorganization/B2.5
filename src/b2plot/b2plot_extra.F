      subroutine get_adas_data(reac,ns,nCv,ne,te,mach)
      use b2mod_types
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_constants
      use b2mod_subsys
      IMPLICIT NONE
*..   input parameters (unchanged on exit)
      integer ns,nCv
      real (kind=R8) :: ne(nCv),te(nCv)
      CHARACTER REAC*9
*..   output data
      real (kind=R8) :: mach(nCv,0:ns-1)
c     helpers for reading in data
      integer IZMAX,IZ1MIN,IZ1MAX
      integer isdimd,iddimd,itdimd
      parameter (isdimd=92+1,iddimd=100,itdimd=100)
      integer nptnl,ndptnl,ndptn,ndptnc,ndcnct,ncnct
      parameter (ndptnl=4,ndptn=128,ndptnc=256,ndcnct=100)
      integer nptn(ndptnl),nptnc(ndptnl,ndptn)
      integer iptnla(ndptnl),iptna(ndptnl,ndptn)
      integer iptnca(ndptnl,ndptn,ndptnc)
      integer icnctv(ndcnct)
      integer iblmx,ismax
      integer isppr(isdimd),ispbr(isdimd),isstgr(isdimd)
      integer idmaxd,itmaxd
      logical lres,lstan,lptn
      real (kind=R8) :: dnr_ams
      character dnr_ele*12
      integer id
      logical first_stage, last_stage, same_element
c     log10 of densities (in 1/cm3), temperatures (in eV)
      real (kind=R8) ::
     *   DDENSD(iddimd),DTEVD(itdimd),logT,logD
c     log10 of rate coefficient (in cm3/s)
      real (kind=R8) ::
     *   DREAC(isdimd,itdimd,iddimd),DTMP(itdimd,iddimd),ans
*..   local variables
      integer is, ic, iCv
      integer, save :: extrap_mode = 1, ncall = 0
      CHARACTER*2 ELNAME
      CHARACTER*200 DSN
*..   procedures
      logical, external :: open_file, streql
      external lower_case, xxdata_11, xerrab

      call subini ('get_adas_data')
      if (ncall.eq.0) then
        call ipgeti ('adas_extrap', extrap_mode)
      end if
      is = 0
      DREAC = 6.0_R8
      do while (is.le.ns-1)
        ELNAME = is_codes(is)
        call lower_case(ELNAME)
        if (trim(ELNAME).eq.'d'.or.trim(ELNAME).eq.'t') ELNAME = 'h'
        DSN = trim(REAC)//'/'//trim(REAC)//'_'//trim(ELNAME)//'.dat'
        if(open_file(29,DSN,
     .   'data.local/adas/adf11|modules/adas/adf11')) then
c     reading in data
cxpb The number 2 is dependent on the type of data being read
          if(index(REAC,'acd').NE.0) then  ! recombination
            id = 1
            first_stage = .false.
            last_stage = .true.
          elseif(index(REAC,'scd').NE.0) then  ! ionisation
            id = 2
            first_stage = .true.
            last_stage = .false.
          elseif(index(REAC,'ccd').NE.0) then  ! charge exchange
            id = 3
            first_stage = .false.
            last_stage = .true.
          elseif(index(REAC,'prb').NE.0) then  ! bremsstrahlung + radiative recombination contribution
            id = 4
            first_stage = .false.
            last_stage = .true.
          elseif(index(REAC,'plt').NE.0) then  ! line radiation
            id = 8
            first_stage = .true.
            last_stage = .false.
          elseif(index(REAC,'zcd').NE.0) then  ! effective charge
            id = 10
            first_stage = .true.
            last_stage = .true.
          elseif(index(REAC,'ycd').NE.0) then  ! effective square charge
            id = 11
            first_stage = .true.
            last_stage = .true.
          elseif(index(REAC,'ecd').NE.0) then  ! effective ionisation potential
            id = 12
            first_stage = .true.
            last_stage = .true.
          else
            write(6,*) 'ADAS data type ',trim(REAC),' not recognized!'
            call subend ()
            return
          endif
          call xxdata_11(29,id,
     *           isdimd,iddimd,itdimd,ndptnl,ndptn,ndptnc,ndcnct,
     *           IZMAX,IZ1MIN,IZ1MAX,
     *           nptnl,nptn,nptnc,iptnla,iptna,iptnca,
     *           ncnct,icnctv,iblmx,ismax,dnr_ele,dnr_ams,
     *           isppr,ispbr,isstgr,
     *           IDMAXD,ITMAXD,DDENSD,DTEVD,DREAC,
     *           lres,lstan,lptn)
          if (IZMAX.ne.nint(zn(is))) then
            write(*,*) 'ADAS tables: wrong element.'
            write(*,*) 'reading file: ',DSN
            call xerrab ('Wrong element in ADAS file')
          endif
          if ((IDMAXD.gt.iddimd).or.(ITMAXD.gt.itdimd)) then
            write(*,*) 'ADAS tables: dimensions too big.'
            write(*,*) 'Change iddimd and itdimd in source code.'
            write(*,*) 'reading file: ',DSN
            call xerrab ('Too big dimensions in ADAS tables')
          endif
c     filling array and computing rate
          do ic=0,ncnct-1
            if(is+ic.ge.ns) cycle
            if(ic.gt.0.and..not.lnext(is,is+ic)) cycle
            if((ic.eq.0.and.first_stage).or.
     &         (ic.eq.ncnct-1.and.last_stage).or.
     &         (ic.gt.0.and.ic.lt.ncnct-1)) then
              if(first_stage) then
                DTMP(1:ITMAXD,1:IDMAXD)=DREAC(ic+1,1:ITMAXD,1:IDMAXD)
              else
                DTMP(1:ITMAXD,1:IDMAXD)=DREAC(ic,1:ITMAXD,1:IDMAXD)
              endif
            elseif (id.lt.10) then
              DTMP(1:ITMAXD,1:IDMAXD)=6.0_R8
            else
              DTMP(1:ITMAXD,1:IDMAXD)=0.0_R8
            endif
            do iCv=1,nCv
              logT=log10(te(iCv)/qe)
              logD=log10(ne(iCv)/1.0e6_R8)
              call int2d(DTMP(1:itdimd,1:iddimd),
     &         extrap_mode,itdimd,iddimd,
     1         DTEVD,ITMAXD,DDENSD,IDMAXD,
     2         logT,logD,ans)
              if(id.lt.10) then
                if((ic.eq.0.and..not.first_stage) .or.
     &             (ic.eq.ncnct-1.and..not.last_stage)) then
                  mach(iCv,is+ic)=0.0_R8
                else
                  mach(iCv,is+ic)=10.0_R8**(ans-6.0_R8)
                end if
              else
                mach(iCv,is+ic)=ans
              endif
            enddo
          enddo
          is = is + ncnct
          close(29)
        else
          write(*,*) 'No ADAS file found : setting data to 0'
          same_element = .true.
          do while (same_element)
            mach(:,is) = 0.0_R8
            ic = is
            is = is + 1
            if (is.eq.ns) then
              same_element = .false.
            else
              same_element = streql(is_codes(is),is_codes(ic))
            end if
          end do
        endif
      enddo
      ncall = ncall + 1
      call subend ()
      return
      end subroutine get_adas_data

      subroutine get_amjuel_data(filnam,h123,reac,crc,
     1 nCv,ne,te,mach)
      use b2mod_types
      use b2mod_constants
      IMPLICIT NONE
#     include <eir_reac.F>
      integer nCv
      real (kind=R8) :: ne(nCv),te(nCv),mach(nCv),D11(0:8,0:8)
      CHARACTER FILNAM*8,H123*4,REAC*9,CRC*3
      integer ir,i,j,ii,iCv,ifex
      real (kind=R8) ::
     1 rmin,rmax,tef,def,dnnn,tei,dej,fp1,fp2,fp3,
     2 s01,s02,ds12,expo1,expo2,ccxm1,ccxm2,coumin,coumax,
     3 demin,demax
      REAL (KIND=R8) :: SLREAC_EXTRAP
      EXTERNAL SLREAC_EXTRAP

      IR=1
      CALL SLREAC_b2plot(IR,FILNAM,H123,REAC,CRC)
      DO J=1,9
        DO I=1,9
          if(h123.eq.'H.-1') then
            D11(J-1,I-1)=CREAC(J,0,IR)
          elseif(
     1     h123.eq.'H.1'.or.
     1     h123.eq.'H.2'.or.
     1     h123.eq.'H.5'.or.
     1     h123.eq.'H.8'.or.
     1     h123.eq.'H.11') then
            D11(J-1,I-1)=CREAC(J,1,IR)
          elseif(
     1     h123.eq.'H.3'.or.
     1     h123.eq.'H.4'.or.
     1     h123.eq.'H.6'.or.
     1     h123.eq.'H.7'.or.
     1     h123.eq.'H.9'.or.
     1     h123.eq.'H.10'.or.
     1     h123.eq.'H.12') then
            D11(J-1,I-1)=CREAC(J,I,IR)
          else
            write(*,*) 'Do not yet know about ',h123
            D11(J-1,I-1)=CREAC(J,I,IR)
          endif
        enddo
      enddo

      rmin=1.0e30_R8
      rmax=1.0e-30_R8
      if (RC2MN(IR).ne.-HUGE(1.0_R8)) then
        demin=RC2MN(IR)+LOG(1.0e-8_R8)
      else
        demin=-huge(1.0_R8)
      endif
      if (RC2MX(IR).ne.HUGE(1.0_R8)) then
        demax=RC2MX(IR)+LOG(1.0e-8_R8)
      else
        demax= huge(1.0_R8)
      endif
      do iCv=1,nCv
        TEF=LOG(TE(iCv)/qe)
        DEF=LOG(ne(iCv)*1.0e-8_R8*1.0e-6_R8)
        dnnn=0.0_R8
        if(
     1   h123.eq.'H.1'.or.
     1   h123.eq.'H.2'.or.
     1   h123.eq.'H.5'.or.
     1   h123.eq.'H.8'.or.
     1   h123.eq.'H.11') then ! single polynomial fit
          IF (TEF.LT.RC1MN(IR)) THEN
            IF (IF1MN(IR).LT.0) THEN ! Linear extrapolation of fit
              S01 = RC1MN(IR)
              S02 = LOG(1.25_R8)+RC1MN(IR)
              DS12 = S02-S01
              EXPO1 = D11(8,0)
              EXPO2 = D11(8,0)
              DO II = 7, 0, -1
                EXPO1 = EXPO1*S01+D11(II,0)
                EXPO2 = EXPO2*S02+D11(II,0)
              END DO
              CCXM1=EXPO1
              CCXM2=EXPO2
              FP1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
              FP2=      (CCXM2-CCXM1)/DS12
              FP3=0.0_R8
              IFEX=3
              COUMIN=EXP(EXPO1)
              DNNN=SLREAC_EXTRAP(TEF,COUMIN,IFEX,FP1,FP2,FP3)
            ELSEIF (IF1MN(IR).GT.0) THEN ! use IF1MN asymptotic expression
              COUMIN=D11(8,0)
              DO II = 7, 0, -1
                COUMIN = COUMIN * RC1MN(IR) + D11(II,0)
              END DO
              IFEX=IF1MN(IR)
              COUMIN=EXP(COUMIN)
              DNNN=SLREAC_EXTRAP(TEF,COUMIN,IFEX,
     .                           FP1L(IR,1),FP1L(IR,2),FP1L(IR,3))
            ELSE
              DO I=0,8
                TEI=RC1MN(IR)**I
                DNNN=DNNN+D11(I,0)*TEI
              ENDDO
            ENDIF
          ELSEIF (TEF.GT.RC1MX(IR)) THEN
            IF (IF1MX(IR).LT.0) THEN
              S01 = RC1MX(IR)
              S02 = LOG(0.75_R8)+RC1MX(IR)
              DS12 = S02-S01
              EXPO1 = D11(8,0)
              EXPO2 = D11(8,0)
              DO II = 7, 0, -1
                EXPO1 = EXPO1*S01+D11(II,0)
                EXPO2 = EXPO2*S02+D11(II,0)
              END DO
              CCXM1=EXPO1
              CCXM2=EXPO2
              FP1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
              FP2=      (CCXM2-CCXM1)/DS12
              FP3=0.0_R8
              IFEX=3
              COUMAX=EXP(EXPO1)
              DNNN=SLREAC_EXTRAP(TEF,COUMAX,IFEX,FP1,FP2,FP3)
            ELSEIF (IF1MX(IR).GT.0) THEN
              COUMAX=D11(8,0)
              DO II = 7, 0, -1
                COUMAX = COUMAX * RC1MX(IR) + D11(II,0)
              END DO
              IFEX=IF1MX(IR)
              COUMAX=EXP(COUMAX)
              DNNN=SLREAC_EXTRAP(TEF,COUMAX,IFEX,
     .                           FP1R(IR,1),FP1R(IR,2),FP1R(IR,3))
            ELSE
              DO I=0,8
                TEI=RC1MX(IR)**I
                DNNN=DNNN+D11(I,0)*TEI
              ENDDO
            ENDIF
          ELSE
            DO I=0,8
              TEI=TEF**I
              DNNN=DNNN+D11(I,0)*TEI
            ENDDO
          ENDIF
          if (earrh(ir).gt.0.0) dnnn=dnnn-earrh(ir)/(te(iCv)/qe)
        elseif(
     1     h123.eq.'H.3'.or.
     1     h123.eq.'H.4'.or.
     1     h123.eq.'H.6'.or.
     1     h123.eq.'H.7'.or.
     1     h123.eq.'H.9'.or.
     1     h123.eq.'H.10'.or.
     1     h123.eq.'H.12') then
          IF (TEF.LT.RC1MN(IR)) THEN
            IF (IF1MN(IR).EQ.0) THEN
              DO I = 0, 8
                TEI=RC1MN(IR)**I
                DO J = 0, 8
                  DEJ=DEF**J
                  DNNN=DNNN+D11(I,J)*TEI*DEJ
                ENDDO
              ENDDO
            ELSEIF (IF1MN(IR).LT.0) THEN
              S01=RC1MN(IR)
              S02=LOG(1.25_R8)+RC1MN(IR)
              DS12 = S02-S01
              EXPO1 = 0.0_R8
              EXPO2 = 0.0_R8
              DO I = 0, 8
                DO J = 0, 8
                  EXPO1=EXPO1+S01**I*DEF**J*D11(I,J)
                  EXPO2=EXPO2+S02**I*DEF**J*D11(I,J)
                ENDDO
              ENDDO
              CCXM1=EXPO1
              CCXM2=EXPO2
              FP1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
              FP2=      (CCXM2-CCXM1)/DS12
              FP3=0.0_R8
              IFEX=3
              COUMIN=EXP(EXPO1)
              DNNN=SLREAC_EXTRAP(TEF,COUMIN,IFEX,FP1,FP2,FP3)
            ELSEIF (IF1MN(IR).GT.0) THEN
              COUMIN=0.0_R8
              DO I = 0, 8
                TEI=RC1MN(IR)**I
                DO J = 0, 8
                  DEJ=DEF**J
                  COUMIN=COUMIN+D11(I,J)*TEI*DEJ
                END DO
              END DO
              IFEX=IF1MN(IR)
              COUMIN=EXP(COUMIN)
              DNNN=SLREAC_EXTRAP(TEF,COUMIN,IFEX,
     .                           FP1L(IR,1),FP1L(IR,2),FP1L(IR,3))
            ENDIF
          ELSEIF (TEF.GT.RC1MX(IR)) THEN
            IF (IF1MX(IR).EQ.0) THEN
              DO I = 0, 8
                TEI=RC1MX(IR)**I
                DO J = 0, 8
                  DEJ=DEF**J
                  DNNN=DNNN+D11(I,J)*TEI*DEJ
                ENDDO
              ENDDO
            ELSEIF (IF1MX(IR).LT.0) THEN
              S01=RC1MX(IR)
              S02=LOG(0.75_R8)+RC1MX(IR)
              DS12 = S02-S01
              EXPO1 = 0.0_R8
              EXPO2 = 0.0_R8
              DO I = 0, 8
                DO J = 0, 8
                  EXPO1=EXPO1+S01**I*DEF**J*D11(I,J)
                  EXPO2=EXPO2+S02**I*DEF**J*D11(I,J)
                ENDDO
              ENDDO
              CCXM1=EXPO1
              CCXM2=EXPO2
              FP1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
              FP2=      (CCXM2-CCXM1)/DS12
              FP3=0.0_R8
              IFEX=3
              COUMAX=EXP(EXPO1)
              DNNN=SLREAC_EXTRAP(TEF,COUMAX,IFEX,FP1,FP2,FP3)
            ELSEIF (IF1MX(IR).GT.0) THEN
              COUMAX=0.0_R8
              DO I = 0, 8
                TEI=RC1MX(IR)**I
                DO J = 0, 8
                  DEJ=DEF**J
                  COUMAX=COUMAX+D11(I,J)*TEI*DEJ
                END DO
              END DO
              IFEX=IF1MX(IR)
              COUMAX=EXP(COUMAX)
              DNNN=SLREAC_EXTRAP(TEF,COUMAX,IFEX,
     .                           FP1R(IR,1),FP1R(IR,2),FP1R(IR,3))
            ENDIF
          ELSEIF (DEF.LT.DEMIN) THEN
            IF (IF2MN(IR).LT.0) THEN
              DO I = 0, 8
                TEI=TEF**I
                DO J = 0, 8
                  DEJ=DEMIN**J
                  DNNN=DNNN+D11(I,J)*TEI*DEJ
                ENDDO
              ENDDO
            ELSEIF (IF2MN(IR).LT.0) THEN
              S01=DEMIN
              S02=LOG(1.25_R8)+DEMIN
              DS12 = S02-S01
              EXPO1 = 0.0_R8
              EXPO2 = 0.0_R8
              DO I = 0, 8
                DO J = 0, 8
                  EXPO1=EXPO1+TEF**I*S01**J*D11(I,J)
                  EXPO2=EXPO2+TEF**I*S02**J*D11(I,J)
                ENDDO
              ENDDO
              CCXM1=EXPO1
              CCXM2=EXPO2
              FP1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
              FP2=      (CCXM2-CCXM1)/DS12
              FP3=0.0_R8
              IFEX=3
              COUMIN=EXP(EXPO1)
              DNNN=SLREAC_EXTRAP(DEF,COUMIN,IFEX,FP1,FP2,FP3)
            ELSEIF (IF2MN(IR).GT.0) THEN
              COUMIN=0.0_R8
              DO I = 0, 8
                TEI=TEF**I
                DO J = 0, 8
                  DEJ=DEMIN**J
                  COUMIN=COUMIN+D11(I,J)*TEI*DEJ
                ENDDO
              ENDDO
              IFEX=IF2MN(IR)
              COUMIN=EXP(COUMIN)
              DNNN=SLREAC_EXTRAP(DEF,COUMIN,IFEX,
     .                           FP2L(IR,1),FP2L(IR,2),FP2L(IR,3))
            ENDIF
          ELSEIF (DEF.GT.DEMAX) THEN
            IF (IF2MX(IR).LT.0) THEN
              DO I = 0, 8
                TEI=TEF**I
                DO J = 0, 8
                  DEJ=DEMAX**J
                  DNNN=DNNN+D11(I,J)*TEI*DEJ
                ENDDO
              ENDDO
            ELSEIF (IF2MX(IR).LT.0) THEN
              S01=DEMAX
              S02=LOG(0.75_R8)+DEMAX
              DS12 = S02-S01
              EXPO1 = 0.0_R8
              EXPO2 = 0.0_R8
              DO I = 0, 8
                DO J = 0, 8
                  EXPO1=EXPO1+TEF**I*S01**J*D11(I,J)
                  EXPO2=EXPO2+TEF**I*S02**J*D11(I,J)
                ENDDO
              ENDDO
              CCXM1=EXPO1
              CCXM2=EXPO2
              FP1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
              FP2=      (CCXM2-CCXM1)/DS12
              FP3=0.0_R8
              IFEX=3
              COUMAX=EXP(EXPO1)
              DNNN=SLREAC_EXTRAP(DEF,COUMAX,IFEX,FP1,FP2,FP3)
            ELSEIF (IF2MX(IR).GT.0) THEN
              COUMAX=0.0_R8
              DO I = 0, 8
                TEI=TEF**I
                DO J = 0, 8
                  DEJ=DEMAX**J
                  COUMAX=COUMAX+D11(I,J)*TEI*DEJ
                ENDDO
              ENDDO
              IFEX=IF2MX(IR)
              COUMAX=EXP(COUMAX)
              DNNN=SLREAC_EXTRAP(DEF,COUMAX,IFEX,
     .                           FP2R(IR,1),FP2R(IR,2),FP2R(IR,3))
            ENDIF
          ELSE
            DO I = 0, 8
              TEI=TEF**I
              DO J = 0, 8
                DEJ=DEF**J
                DNNN=DNNN+D11(I,J)*TEI*DEJ
              ENDDO
            ENDDO
          ENDIF
        ELSE
c         write(*,*) 'Do not yet know about ',h123
        ENDIF
        dnnn=dnnn+log(1.0e-6_R8)
        if(dnnn.gt.log(1.0e30_R8)) then
          dnnn=1.0e30_R8
        elseif(dnnn.lt.log(1.0e-30_R8)) then
          dnnn=1.0e-30_R8
        else
          DNNN=EXP(DNNN)
        endif
        rmin=min(rmin,dnnn)
        rmax=max(rmax,dnnn)
        mach(iCv)=dnnn
      enddo
      write(*,*) 'rmin,rmax= ',rmin,rmax
      return
      end subroutine get_amjuel_data

      SUBROUTINE SLREAC_b2plot(IR,FILNAM,H123,REAC,CRC)
      use b2mod_types
      use b2mod_subsys
c
C  input
C    filnam: read a&m data from file filnam, e.g. AMJUEL
c    ir    : store data on Eirene array creac(...,...,ir)
c    H123  : identifier for reaction in filnam, e.g. H2
c    REAC  : number of reaction in filnam, e.g. 2.2.5
c    crc   : type of data process e.g. EI, CX, OT, etc
C  output
c    iswr  : Eirene flag for type of process  (1,2,...7)
c    creac : Eirene storage array for a&m data CREAC(9,0:9,IR)
c    modcol: see below
C    IFTFLG=IFTFLG(IH): flag for type of fitting expression ("fit-flag=...")
C    IH  internally derived from ISW, for different types of data:
C          0 for interaction potential, or differential cross-sections (ISW=0)
C          1 for cross-section, (ISW=1)
C          2 for rate coeff, (ISW=2,3,4)
C          3 for mom-weighted rate coeff. (ISW=5,6,7)
C          4 for energy-weighted rate coeff. (ISW=8,9,10)
C          5 for other quantities, population densities, etc.. (ISW=11,12)
C
C       IFTFLG(IH) DEFAULTS:
C
C           CASE IH=0  (H.0):
C       IFTFLG(0)   =2,  FOR INTERACTION POTENTIAL (GEN. MORSE)
C
C           CASE  IH=1  (H.1):
C       IFTFLG(1)   =0,  CROSS-SECTION (9-POLYNOMIAL)
C                   =3,  cross-section (ionisation/excitation cross-section
C                           formula (METHANE,...)
C           CASE  IH=2,3....,10 (H.2, H.3,....H.10)
C       IFTFLG(...  =0,  FOR RATE COEFFICIENTS (9-POLYNOMIAL, 9X9-DOUBLE POLYNOMIAL)
C                   =10, FOR RATE COEFFICIENTS (CONSTANT)
C                   =100 FOR RATE, not rate coefficient,
C                   =110 FOR RATE, not rate coefficient, (CONSTANT)
C
      IMPLICIT NONE
C  READ A&M DATA FROM THE FILES INTO EIRENE ARRAY CREAC
C
C  INPUT (PARAMETER LIST):
C    IR: REACTION NUMBER IN CREAC
C    FILNAM: HYDHEL, METHANE, AMJUEL OR CONST.
C
C  OUTPUT (IN COMMON COMXS):
C    READ DATA FROM "FILNAM" INTO ARRAY "CREAC"
C    DEFINE PARAMETER MODCLF(IR) (4 DIGITS NMLK)
C    FIRST DECIMAL  K           =1  CROSS-SECTION AVAILABLE
C                                   (ON CREAC(..,0,IR))
C                   K           =0  ELSE
C    SECOND DECIMAL L           =1  <SIGMA V> FOR ONE
C                                   PARAMETER E (E.G.
C                                   PROJECTILE ENERGY OR ELECTRON
C                                   DENSITY) AVAILABLE
C                                   (ON CREAC(..,1,IR))
C                               =2  <SIGMA V> FOR
C                                   9 PROJECTILE ENERGIES AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                               =3  <SIGMA V> FOR
C                                   9 ELECTRON DENSITIES  AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                   L           =0  ELSE
C    THIRD  DECIMAL M               DATA FOR MOMENTUM EXCHANGE
C                                   TO BE WRITTEN
C    FOURTH DECIMAL N           =1  DELTA-E FOR ONE (E.G. E=0)
C                                   PROJECTILE ENERGY AVAILABLE
C                                   (ON CREAC(..,1,IR))
C                               =2  DELTA-E FOR
C                                   9 PROJECTILE ENERGIES AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                               =3  DELTA-E FOR
C                                   9 ELECTRON DENSITIES  AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                   N           =0  ELSE
C
#     include <eir_reac.F>
      real(kind=R8) :: CONST, FPARM(3)
      REAL(kind=R8) :: R1MN, R1MX, R2MN, R2MX
      integer IR,I0,ISW,J,IND,INDADD,I,K,IH,IANF,IREAC,INDFF,IFLG
      integer IFTFLG(0:5)
      CHARACTER*80 ZEILE, LAST_TEX, ULINE
      CHARACTER FILNAM*8,H123*4,REAC*9,CRC*3
      CHARACTER*1 BACK
      CHARACTER*8 SECTION
      CHARACTER*9 FITFLAG
      CHARACTER*11 REACSTR
      CHARACTER*4 CHR, CHRA, BEND
      CHARACTER*3 CH1L, CH1R, CH2L, CH2R
      CHARACTER*7 C1L, C1R, C2L, C2R
      LOGICAL LCONST
      LOGICAL :: LGC1MIN,LGC1MAX,LGC2MIN,LGC2MAX
      logical open_file
      external open_file, xerrab, xertst
C
      call subini('slreac_b2plot')
      call xertst(ir.le.nreac, 'SLREAC_b2plot: IR > NREAC')

! defining backslash character
      BACK=ACHAR(92)
      SECTION=BACK//'section'
      BEND=BACK//'end'
! set some defaults
      ISWR(IR)=0
      CONST=0.
      CHR=' l0 '
      I0=0
! header of dataset (optional)
      FITFLAG='fit-flag '
      IFTFLG=0
      EARRH(IR) = 0.0_R8
      RC1MN(IR) =-HUGE(1.0_R8)
      RC1MX(IR) = HUGE(1.0_R8)
      RC2MN(IR) =-HUGE(1.0_R8)
      RC2MX(IR) = HUGE(1.0_R8)
C     Defaults: no asymptotics
      CH1L='ll0'
      CH1R='lr0'
      CH2L='ll0'
      CH2R='lr0'
      C1L = 'XXMIN'
      C1R = 'XXMAX'
      C2L = 'YYMIN'
      C2R = 'YYMAX'
C
      IF (INDEX(CRC,'EI').NE.0.OR.
     .    INDEX(CRC,'DS').NE.0) ISWR(IR)=1
      IF (INDEX(CRC,'CX').NE.0) ISWR(IR)=3
      IF (INDEX(CRC,'II').NE.0.OR.
     .    INDEX(CRC,'PI').NE.0) ISWR(IR)=4
      IF (INDEX(CRC,'EL').NE.0) ISWR(IR)=5
      IF (INDEX(CRC,'RC').NE.0) ISWR(IR)=6
      IF (INDEX(CRC,'OT').NE.0) ISWR(IR)=7
C
      IF (INDEX(FILNAM,'AMJUEL').NE.0) THEN
        LCONST=.FALSE.
        call xertst(open_file(29,'AMJUEL',
     1   'data.local|modules/Eirene/Database'), 'No AMJUEL found')
      ELSEIF (INDEX(FILNAM,'METHAN').NE.0) THEN
        LCONST=.FALSE.
        call xertst(open_file(29,'METHANE',
     1   'data.local|modules/Eirene/Database'), 'No METHANE found')
      ELSEIF (INDEX(FILNAM,'HYDHEL').NE.0) THEN
        LCONST=.FALSE.
        call xertst(open_file(29,'HYDHEL',
     1   'data.local|modules/Eirene/Database'), 'No HYDHEL found')
      ELSEIF (INDEX(FILNAM,'AMMONX').NE.0) THEN
        LCONST=.FALSE.
        call xertst(open_file(29,'AMMONX',
     1   'data.local|modules/Eirene/Database'), 'No AMMONX found')
      ELSE
        LCONST=.TRUE.
        READ (REAC(1:9),'(E9.2)') CONST
      ENDIF
C
      IF (H123(4:4).EQ.' ') THEN
        READ (H123(3:3),'(I1)') ISW
      ELSE
        READ (H123(3:4),'(I2)') ISW
      ENDIF
C
      REACSTR=REPEAT(' ',11)
      IANF=VERIFY(REAC,' ')
      IF (IANF > 0) THEN
        IREAC=INDEX(REAC(IANF:),' ')-1
        IF (IREAC.LT.0) IREAC=LEN(REAC(IANF:))
        REACSTR(2:IREAC+1)=REAC(IANF:IREAC+IANF-1)
C  ADD ONE MORE BLANK, IF POSSIBLE
        IREAC=IREAC+2
      ELSE
        IF (.NOT.LCONST) THEN
          call xerrab ('No reaction name specified')
        END IF
      END IF
C
      IF (ISW.EQ.0) THEN
        CHR=' p0 '
        CHRA=' p-1'
        I0=-1
        MODCLF(IR)=MODCLF(IR)+1
        IFLG=0
        IFTFLG(IFLG)=2
      ELSEIF (ISW.EQ.1) THEN
        CHR=' a0 '
        CHRA=' a-1'
        CH1L='al0'
        CH1R='ar0'
        I0=0
        MODCLF(IR)=MODCLF(IR)+10
        IFLG=1
        IFTFLG(IFLG)=0
        C1L = 'ELABMIN'
        C1R = 'ELABMAX'
      ELSEIF (ISW.EQ.2) THEN
        CHR=' b0 '
        CHRA=' b-1'
        CH1L='bl0'
        CH1R='br0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+100
        IFLG=2
        IFTFLG(IFLG)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'
      ELSEIF (ISW.EQ.3) THEN
        CHR=' c0 '
        CHRA=' c-1'
        CH1L='cl0'
        CH1R='cr0'
        CH2L='cb0'
        CH2R='ct0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+200
        IFLG=2
        IFTFLG(IFLG)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'E2MIN'
        C2R = 'E2MAX'
      ELSEIF (ISW.EQ.4) THEN
        CHR=' d0 '
        CHRA=' d-1'
        CH1L='dl0'
        CH1R='dr0'
        CH2L='db0'
        CH2R='dt0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+300
        IFLG=2
        IFTFLG(IFLG)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'N2MIN'
        C2R = 'N2MAX'
      ELSEIF (ISW.EQ.5) THEN
        CHR=' e0 '
        CHRA=' e-1'
        CH1L='el0'
        CH1R='er0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+1000
        IFLG=3
        IFTFLG(IFLG)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'
      ELSEIF (ISW.EQ.6) THEN
        CHR=' f0 '
        CHRA=' f-1'
        CH1L='fl0'
        CH1R='fr0'
        CH2L='fb0'
        CH2R='ft0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+2000
        IFLG=3
        IFTFLG(IFLG)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'E2MIN'
        C2R = 'E2MAX'
      ELSEIF (ISW.EQ.7) THEN
        CHR=' g0 '
        CHRA=' g-1'
        CH1L='gl0'
        CH1R='gr0'
        CH2L='gb0'
        CH2R='gt0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+3000
        IFLG=3
        IFTFLG(IFLG)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'N2MIN'
        C2R = 'N2MAX'
      ELSEIF (ISW.EQ.8) THEN
        CHR=' h0 '
        CHRA=' h-1'
        CH1L='hl0'
        CH1R='hr0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+10000
        IFLG=4
        IFTFLG(IFLG)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'
      ELSEIF (ISW.EQ.9) THEN
        CHR=' i0 '
        CHRA=' i-1'
        CH1L='il0'
        CH1R='ir0'
        CH2L='ib0'
        CH2R='it0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+20000
        IFLG=4
        IFTFLG(IFLG)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'E2MIN'
        C2R = 'E2MAX'
      ELSEIF (ISW.EQ.10) THEN
        CHR=' j0 '
        CHRA=' j-1'
        CH1L='jl0'
        CH1R='jr0'
        CH2L='jb0'
        CH2R='jt0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+30000
        IFLG=4
        IFTFLG(IFLG)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'N2MIN'
        C2R = 'N2MAX'
      ELSEIF (ISW.EQ.11) THEN
        CHR=' k0 '
        CHRA=' k-1'
        CH1L='kl0'
        CH1R='kr0'
        I0=1
        IFLG=5
        IFTFLG(IFLG)=0
        C1L = 'P1MIN'
        C1R = 'P1MAX'
      ELSEIF (ISW.EQ.12) THEN
        CHR=' l0 '
        CHRA=' l-1'
        CH1L='ll0'
        CH1R='lr0'
        CH2L='lb0'
        CH2R='lt0'
        I0=1
        IFLG=5
        IFTFLG(IFLG)=0
        C1L = 'P1MIN'
        C1R = 'P1MAX'
        C2L = 'P2MIN'
        C2R = 'P2MAX'
      ENDIF
C
      IF (.NOT.LCONST) THEN
  100   READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,'##BEGIN DATA HERE##').EQ.0) GOTO 100

        LAST_TEX=REPEAT(' ',80)
    1   READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,H123).EQ.0) THEN
          IF (INDEX(ZEILE,BACK).NE.0) LAST_TEX=ZEILE
          GOTO 1
        ELSE
          IF (INDEX(LAST_TEX,SECTION).EQ.0 .AND.
     .        INDEX(ZEILE,SECTION).EQ.0) GOTO 1
        END IF
C
    2   READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,'H.').NE.0 .and.
     .      INDEX(ZEILE,'section').NE.0) GOTO 990
        IF (INDEX(ZEILE,'Reaction ').EQ.0.or.
     .      INDEX(ZEILE,REACSTR(1:ireac)).EQ.0) GOTO 2
      ENDIF
C
C  SINGLE PARAM. FIT, ISW=1,2,5,8,11
      IF (ISW.EQ.0.OR.ISW.EQ.1.OR.ISW.EQ.2.OR.ISW.EQ.5.OR.ISW.EQ.8.OR.
     .    ISW.EQ.11) THEN
        IF (.NOT.LCONST) THEN
    3     READ (29,'(A80)',END=990) ZEILE
          INDFF=INDEX(ZEILE,FITFLAG)
c  skip empty lines in header
          IF (INDEX(ZEILE,CHR)+INDEX(ZEILE,CHRA)+INDFF.EQ.0) GOTO 3
c  input line found which either contains fit-flag FITFLAG,
c       or the reaction identifier a0,b0,...k0
c       or the Arrhenius prefactor a-1,b-1,...k-1
          IF (INDFF > 0) THEN
            READ (ZEILE(INDFF+8:80),*) IFTFLG(IFLG)
            GOTO 3
          ENDIF
          IF (INDEX(ZEILE,CHRA) > 0) THEN
c  read Arrhenius parameter
            READ (ZEILE((INDFF+8):80),*) EARRH(IR)
            GOTO 3
          ENDIF

C  read only one constant: (FIT-FLAG = 10, 110, ....)
          IF (MOD(IFTFLG(IFLG),100) == 10) THEN
            IND=INDEX(ZEILE,CHR(2:2)) ! position of a,b,...k
            READ (ZEILE((IND+2):80),'(E20.12)') CREAC(1,1,IR)
C  READ 9 FIT COEFFICIENTS, SEPARATED BY 'CHR' FIXED FORMAT E20.12
C  THREE LINES WITH THREE DATA PER LINE
c  READ ONLY THE NON-ZERO PARAMETERS
          ELSE
            DO J=0,2
              IND=0
              DO I=1,3
                INDADD=INDEX(ZEILE((IND+1):80),CHR(2:2))
                IF (INDADD.EQ.0) GOTO 9
                IND=IND+INDEX(ZEILE((IND+1):80),CHR(2:2))
                READ (ZEILE((IND+2):80),'(E20.12)') CREAC(J*3+I,I0,IR)
              ENDDO
              READ (29,'(A80)',END=990) ZEILE
            ENDDO
    9       CONTINUE
          ENDIF
        ELSE
          CREAC(1,I0,IR)=CONST
          DO 6 J=2,9
            CREAC(J,I0,IR)=0.0_R8
    6     CONTINUE
        ENDIF
C
      ELSEIF (ISW.EQ.3.OR.ISW.EQ.4.OR.ISW.EQ.6.OR.ISW.EQ.7.OR.
     .        ISW.EQ.9.OR.ISW.EQ.10.OR.ISW.EQ.12) THEN
C  TWO PARAM. FIT, ISW=3,4,6,7,9,10,12

        DO 11 J=0,2
   16     READ (29,'(A80)',END=990) ZEILE
          INDFF=INDEX(ZEILE,FITFLAG)
          IF (INDEX(ZEILE,'Index')+INDFF.EQ.0) GOTO 16
          IF (INDFF > 0) THEN
c  read parameter for type of fitting expression from data file
            READ (ZEILE((INDFF+8):80),*) IFTFLG(IFLG)
            GOTO 16
          ENDIF
          READ (29,'(1X)')
          IF (MOD(IFTFLG(IFLG),100) == 10) THEN
C  IFTFLG = 10, 110,  210,....ETC: READ ONLY ONE CONSTANT PARAMETER
            READ (29,*) IH,CREAC(1,1,IR)
            EXIT
          ELSE
            DO 17 I=1,9
              READ (29,*) IH,(CREAC(I,K,IR),K=J*3+1,J*3+3)
cdr           write (6,*) IH,(CREAC(I,K,IR),K=J*3+1,J*3+3)
   17       CONTINUE
          END IF
   11   CONTINUE
        READ (29,'(A80)',END=990) ZEILE
C
      ENDIF

C  NEXT: READ ASYMPTOTICS INFORMATION FROM ATOMIC DATA FILE
C        HYDHEL, AMJUEL, AMMONX, H2VIBR, METHANE.

C FOR 1D OR 2D DATASETS. 4 BOUNDARIES, LEFT1, RIGHT1, LEFT2, RIGHT2.
C FOR 1D: ONLY "LEFT1" AND "RIGHT1" ARE USED

      IF (ISW.EQ.0) GOTO 2000    ! NO ASYMPTOTICS FOR POTENTIALS

c  INDICATE, IF EXTRAPOLATION COEFFICIENTS (fp1l,fp1r,fp2l,fp2r) ARE FOUND ON DATA FILE
c  DEFAULT: NO DATA FOUND
      LGC1MIN=.FALSE.
      LGC1MAX=.FALSE.
      LGC2MIN=.FALSE.
      LGC2MAX=.FALSE.

C  FLAG FOR CHOICE OF EXTRAPOLATION OPTION:
C  DEFAULT ASYMPTOTIC EXPRESSION  (...=0): NO ASYMPTOTICS
C  DEFAULT ASYMPTOTIC EXPRESSION  (...=1): SET TO ZERO BEYOND LAST VALID POINT
C  DEFAULT ASYMPTOTIC EXPRESSION  (...=4): TAKE LAST VALID POINT AT r1mn,r1mx,....
C  DEFAULT ASYMPTOTIC EXPRESSION  (...=5): 2ND-ORDER POLYNOMIAL BEYOND LAST VALID POINT
C                                          (OLD DEFAULT FOR CROSS-SECTIONS)
c  AND EXTRAPOLATE CONSTANT FROM THERE
      IF1MN(IR) = 0
      IF1MX(IR) = 0
      IF2MN(IR) = 0
      IF2MX(IR) = 0

cdr Try to read asymptotics, unless already read: jfeximn,jfeximx /= 0
cdr Even if jfeximn,jfeximx /= 0, read anyway, but later: do not use

c  data for one "reaction IR" are located between "\begin" and "\end"
!     BEND="\end" stops looking for asymptotics
      DO WHILE(INDEX(ZEILE,BEND) == 0)

        ULINE = ZEILE
        CALL CHCASE(1,ULINE)
c
c  at this point we have found a card ZEILE which contains
c  one of the extrapolation parameter identifiers al0,ar0,....,k0l,k0r
c  that corresponds to the H.1, ....H.12 type of data IR.
c  next: read up to three fit coefficients FP.L OR FP.R.
c  CHR(2:2) is set to either character a,b,c,....,or k
c
        IF (INDEX(ZEILE,CH1L) /= 0) THEN
          CALL SLREAC_READ_COEFFS (ZEILE,CHR(2:2),FPARM)
          FP1L(IR,1:3) = FPARM(1:3)
          LGC1MIN = .TRUE.
        ENDIF
        IF (INDEX(ZEILE,CH1R) /= 0) THEN
          CALL SLREAC_READ_COEFFS (ZEILE,CHR(2:2),FPARM)
          FP1R(IR,1:3) = FPARM(1:3)
          LGC1MAX = .TRUE.
        END IF

        IF (INDEX(ZEILE,CH2L) /= 0) THEN
          CALL SLREAC_READ_COEFFS (ZEILE,CHR(2:2),FPARM)
          FP2L(IR,1:3) = FPARM(1:3)
          LGC2MIN = .TRUE.
        END IF
        IF (INDEX(ZEILE,CH2R) /= 0) THEN
          CALL SLREAC_READ_COEFFS (ZEILE,CHR(2:2),FPARM)
          FP2R(IR,1:3) = FPARM(1:3)
          LGC2MAX = .TRUE.
        END IF
c
c  currently foreseen asymptotic data identifiers in data files:
c  c1l,c2l,c1r,c2r: ELABMIN, ELABMAX,
c                   T1MIN,T1MAX,E2MIN,E2MAX,N2MIN,N2MAX,
c                   P1MIN,P1MAX,P2MIN,P2MAX
        IF (INDEX(ULINE,TRIM(C1L)) /= 0) THEN
          CALL SLREAC_READ_RANGE(ULINE,C1L,'EXT-FLG',R1MN,IF1MN(IR))
          RC1MN(IR) = LOG(R1MN)
c  old default: 2nd-order polynomial beyond valid range, with coefs. FPL1
          IF (IF1MN(IR) == 0 .AND. LGC1MIN) IF1MN(IR) = 5
c  default extrapolation from r1mn (by constant continuation) will be: jfex1mn=4
          IF (IF1MN(IR) == 0) IF1MN(IR) = 4
        END IF
        IF (INDEX(ULINE,TRIM(C1R)) /= 0) THEN
          CALL SLREAC_READ_RANGE(ULINE,C1R,'EXT-FLG',R1MX,IF1MX(IR))
          RC1MX(IR) = LOG(R1MX)
c  old default: 2nd-order polynomial beyond valid range, with coefs. FPR1
          IF (IF1MX(IR) == 0 .AND. LGC1MAX) IF1MX(IR) = 5
c  default extrapolation from r1mx (by constant continuation) will be: jfex1mx=4
          IF (IF1MX(IR) == 0) IF1MX(IR) = 4
        END IF
        IF (INDEX(ULINE,TRIM(C2L)) /= 0) THEN
          CALL SLREAC_READ_RANGE(ULINE,C2L,'EXT-FLG',R2MN,IF2MN(IR))
          RC2MN(IR) = LOG(R2MN)
c  old default: 2nd-order polynomial beyond valid range, with coefs. FPL2
          IF (IF2MN(IR) == 0 .AND. LGC2MIN) IF2MN(IR) = 5
c  default extrapolation from r2mn (by constant continuation) will be: jfex2mn=4
          IF (IF2MN(IR) == 0) IF2MN(IR) = 4
        END IF
        IF (INDEX(ULINE,TRIM(C2R)) /= 0) THEN
          CALL SLREAC_READ_RANGE(ULINE,C2R,'EXT-FLG',R2MX,IF2MX(IR))
          RC2MX(IR) = LOG(R2MX)
c  old default: 2nd-order polynomial beyond valid range, with coefs. FPr2
          IF (IF2MX(IR) == 0 .AND. LGC2MAX) IF2MX(IR) = 5
c  default extrapolation from r2mx (by constant continuation) will be: jfex2mx=4
          IF (IF2MX(IR) == 0) IF2MX(IR) = 4
        END IF

        READ (29,'(A80)',END=990) ZEILE
      END DO
 2000 CONTINUE
C
      CLOSE (UNIT=29)
C
      call subend ()
      RETURN
C
  990 WRITE (6,*) ' NO DATA FOUND FOR REACTION ',H123,' ',REAC,
     .            ' IN DATASET ',FILNAM
      WRITE (6,*) ' IR,MODCLF(IR) ',IR,MODCLF(IR)
      CLOSE (UNIT=29)
      call subend ()
      RETURN
      END SUBROUTINE SLREAC_b2plot

      SUBROUTINE SLREAC_READ_COEFFS (ZEILE,CH,FP)
      use b2mod_types
      implicit none
c  ZEILE is in upper case, and CH is found.
c  CH is: a,b,c,.....k,l, (depending on H.1, H.2, ...H.12)
c  read up to three parameters FP(i), i=1,3 formatted: '(E20.12)'
c  to be used for extrapolation
      CHARACTER(80), INTENT(IN) :: ZEILE
      CHARACTER(1), INTENT(IN) :: CH
      REAL(kind=R8), INTENT(OUT) :: FP(3)
      INTEGER :: I, IND, INC

      IND=0
      DO I=1,3
        INC=INDEX(ZEILE((IND+1):80),CH)
        IF (INC.GT.0) THEN
           IND=IND+INDEX(ZEILE((IND+1):80),CH)
           READ (ZEILE((IND+3):80),'(E20.12)') FP(I)
        ENDIF
      END DO
      RETURN

      END SUBROUTINE SLREAC_READ_COEFFS

      SUBROUTINE SLREAC_READ_RANGE(ZEILE,KEY1,KEY2,RNG,IFX)
      use b2mod_types
      implicit none
c  called from slreac, after the original fit coefficients for reaction IR
c  are read.
c  At this point an extrapolation card ZEILE belonging to this reaction IR
c  has already been found.
c
c  This routine:
c  Reads validity range from atomic data file:
c  Search in ZEILE for key1, key2 and return: RNG, IFX

c  key1: 'ELABMIN','ELABMAX','T1MIN','T1MAX','E2MIN','E2MAX',
C        'N2MIN','N2MAX','P1MIN','P1MAX','P1MIN','P1MAX'=,
C                    read RNG (unformatted, real)
c  key2: 'EXT-FLG'= ,read IFX (unformatted, integer)

      CHARACTER(80), INTENT(IN) :: ZEILE
      CHARACTER(7), INTENT(IN) :: KEY1, KEY2
      REAL(kind=R8), INTENT(OUT) :: RNG
      INTEGER, INTENT(OUT) :: IFX
      INTEGER :: IND1, IND2, INDE, INDP, INDX, INDA, INDG
      CHARACTER(20) :: FORM

      IND1 = INDEX(ZEILE,TRIM(KEY1))
      IND2 = INDEX(ZEILE,TRIM(KEY2))

      RNG = 0.0_R8
      IFX = 0

      IF (IND1 > 0) THEN
        INDG = INDEX(ZEILE,'=')
        INDE = INDG + VERIFY(ZEILE(INDG+1:),'+-0123456789DEed. ') - 1
        INDP = SCAN(ZEILE(INDG+1:),'.')
        INDX = SCAN(ZEILE(INDG+1:),'EDed')
        INDA = SCAN(ZEILE(INDG+1:),'+-0123456789.')
        FORM=REPEAT(' ',20)
        WRITE (FORM,'(A2,I0,A1,I0,A1)')
     .         '(E',INDX+3-INDA+1,'.',INDX-INDP-1,')'
        READ (ZEILE(INDG+INDA:INDE),FORM) RNG
      END IF

      IF (IND2 > 0) READ (ZEILE(IND2+7:),*) IFX
      RETURN

      END SUBROUTINE SLREAC_READ_RANGE

      FUNCTION SLREAC_EXTRAP(PARM,FPARLIM,IFLAG,FP1,FP2,FP3)
      use b2mod_types
      IMPLICIT NONE

      REAL(kind=R8), INTENT(IN) :: PARM, FPARLIM, FP1, FP2, FP3
      INTEGER, INTENT(IN)  :: IFLAG
      REAL(kind=R8)             :: X, EL, SLREAC_EXTRAP

      IF (IFLAG.EQ.1) THEN
C JUST RETURN ZERO BELOW THRESHOLD
        SLREAC_EXTRAP = 1.0e-30_R8
      ELSEIF (IFLAG.EQ.2) THEN
C  JANEV TYPE 1, HIGH ENERGY EXTRAPOLATION OF CROSS-SECTIONS
C  EXTRAPOLATION AT HIGH ENERGY END FOR REACTIONS WITH A NONZERO THRESHOLD
C  FP1 SHOULD BE = E_THRESHOLD (EV)
C  FP2 SCALING FACTOR, SHOULD BE COMPUTED FROM FPARLIM
C  FP3=-N
        EL = EXP(PARM)
        X = EL/FP1
        SLREAC_EXTRAP=FP2*X**FP3*LOG(X)
      ELSEIF (IFLAG.EQ.3) THEN
C  LINEAR EXTRAPOLATION OF LN(FF) IN PARM= LN(P)
        SLREAC_EXTRAP=EXP(FP1+FP2*PARM)
      ELSEIF (IFLAG.EQ.4) THEN
C   EXTRAPOLATION WITH A CONSTANT
        SLREAC_EXTRAP=FPARLIM
      ELSEIF (IFLAG.EQ.5) THEN
C  QUADRATIC EXTRAPOLATION OF LN(FF) IN PARM=LN(PP)
        SLREAC_EXTRAP=EXP(FP1+FP2*PARM+FP3*PARM**2)
      ELSE
        SLREAC_EXTRAP=0.0_R8
      ENDIF
      RETURN

      END FUNCTION SLREAC_EXTRAP

      subroutine old_chang_hinton(mpg,geo,switch,ns,ismain,
     1 safetyfactor,ne,te,ti,na,rt,transport)
      use b2mod_types
      use b2us_map
      use b2us_geo
      use b2us_plasma
      use b2mod_switches
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo
      type(switches), intent(in) :: switch
      type(B2Rates), intent(in) :: rt
      integer ns,ismain
      real (kind=R8) :: lnlam(mpg%nCv),
     2 ne(mpg%nCv), te(mpg%nCv), ti(mpg%nCv),
     3 na(mpg%nCv,0:ns-1), safetyfactor(mpg%nFt),
     2 transport(mpg%nCv,0:ns-1)

      integer iCv,iFc,iFt,i,is,icmin,icmax
      real (kind=R8) :: crmin,crmax,czmin,czmax,
     1 cr2max,cr2min,cr1max,cr1min,
     2 bave,epsilon,vthi,omegaci,
     1 rhoi,alpha,deltap,f1,f2,taui,nui,nustari,mustari,chi
      external xerrab

      call subini('old_chang_hinton')
      call b2tlnl(mpg%nCv,switch,1,te,ti,ne,lnlam)
      transport=0.0_R8
      do iFt=1,mpg%nFt
        write(*,*) 'iFt: ',iFt
        iCv=mpg%ftCv(mpg%ftCvP(iFt,1))
        if(.not.mpg%cvOnClosedSurface(iCv)) cycle
        icmin=iCv
        icmax=iCv
        crmin=geo%cvX(iCv)
        crmax=crmin
        czmin=geo%cvY(iCv)
        czmax=czmin
        do i=2,mpg%ftCvP(iFt,2)
          iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
          if(geo%cvX(iCv).lt.crmin) icmin=iCv
          if(geo%cvY(iCv).gt.crmax) icmax=iCv
          crmin=min(crmin,geo%cvX(iCv))
          crmax=max(crmax,geo%cvX(iCv))
          czmin=min(czmin,geo%cvY(iCv))
          czmax=max(czmax,geo%cvY(iCv))
        enddo
        write(*,*) 'crmin,crmax,czmin,czmax',
     1   crmin,crmax,czmin,czmax
        iCv=mpg%ftCv(mpg%ftCvP(iFt,1))
        bave=geo%cvX(iCv)*2.0_R8/(crmax+crmin)*geo%cvBb(iCv,3)
        write(*,*) 'ne,te,ti,b,q: ',
     1   ne(iCv),te(iCv)/qe,ti(iCv)/qe,bave,safetyfactor(iFt)
        epsilon=(crmax-crmin)/(crmax+crmin)
        write(*,*) 'epsilon: ',epsilon
        write(*,*) 'lnlam: ',lnlam(iCv)
        cr1max = 0.0_R8
        cr2max = 0.0_R8
        do i=1,mpg%cvFcP(icmax,2)
          iFc = mpg%cvFc(mpg%cvFcP(icmax,1)+i-1)
          if (mpg%fcFs(iFc).ne.0) cycle
          if (cr1max.eq.0.0_R8) then
            cr1max = 0.5_R8*
     .       (geo%cvX(mpg%fcVx(iFc,1))+geo%cvX(mpg%fcVx(iFc,2)))
          else if (cr2max.eq.0.0_R8) then
            cr2max = 0.5_R8*
     .       (geo%cvX(mpg%fcVx(iFc,1))+geo%cvX(mpg%fcVx(iFc,2)))
          end if
        end do
        cr1min = 0.0_R8
        cr2min = 0.0_R8
        do i=1,mpg%cvFcP(icmin,2)
          iFc = mpg%cvFc(mpg%cvFcP(icmin,1)+i-1)
          if (mpg%fcFs(iFc).ne.0) cycle
          if (cr1min.eq.0.0_R8) then
            cr1min = 0.5_R8*
     .       (geo%cvX(mpg%fcVx(iFc,1))+geo%cvX(mpg%fcVx(iFc,2)))
          else if (cr2min.eq.0.0_R8) then
            cr2min = 0.5_R8*
     .       (geo%cvX(mpg%fcVx(iFc,1))+geo%cvX(mpg%fcVx(iFc,2)))
          end if
        end do
        write(*,*) 'cr2max,cr2min,cr1max,cr1min: ',
     1   cr2max,cr2min,cr1max,cr1min
        deltap=((cr2max+cr2min)-(cr1max+cr1min))/
     1         ((cr2max-cr2min)-(cr1max-cr1min))
        write(*,*) 'deltap: ',deltap
        f1=(1_R8+1.5_R8*(epsilon**2+epsilon*deltap)+
     1         0.375_R8*epsilon**3*deltap)/(1_R8+0.5_R8*epsilon*deltap)
        write(*,*) 'f1: ',f1
        f2=sqrt(1_R8-epsilon**2)*(1.0_R8+0.5_R8*epsilon*deltap)/
     1         (1_R8+deltap/epsilon*(sqrt(1_R8-epsilon**2)-1.0_R8))
        write(*,*) 'f2: ',f2
        do is=0,ns-1
          if(.not.is_neutral(is)) then
            write(*,*) 'za,am,na: ', rt%rza(iCv,is),am(is),na(iCv,is)
            vthi=sqrt(ti(iCv)/mp/am(is))
            write(*,*) 'vthi: ',vthi
            omegaci=rt%rza(iCv,is)*qe*bave/mp/am(is)
            write(*,*) 'omegaci: ',omegaci
            rhoi=sqrt(2.0_R8)*vthi/omegaci
            write(*,*) 'rhoi: ',rhoi
            alpha=na(iCv,is)*rt%rz2(iCv,is)/
     1            na(iCv,ismain)/rt%rz2(iCv,ismain)  !!! think about this !!!
            write(*,*) 'alpha: ',alpha
            taui=12.0_R8*pi**1.5_R8*eps0**2*
     1             (am(is)*mp)**0.5_R8*ti(iCv)**1.5_R8/
     2             (na(iCv,is)*rt%rz2(iCv,is)**2*qe**4*lnlam(iCv))
            write(*,*) 'taui: ',taui
            nui=1.0_R8/taui
            nustari=nui*(crmax+crmin)/2.0_R8*safetyfactor(iFt)/
     1             epsilon**(1.5_R8)/vthi
            write(*,*) 'nustari: ',nustari
            mustari=nustari*(1.0_R8+1.54_R8*alpha)
            write(*,*) 'mustari: ',mustari
            chi=epsilon**(-1.5_R8)*
     1             safetyfactor(iFt)**2*rhoi**2/taui*
     2             ((0.66_R8*(1.0_R8+1.54_R8*alpha)+
     3              (1.88_R8*epsilon**0.5_R8-1.54_R8*epsilon)*
     4              (1_R8+3.75_R8*alpha))/
     5              (1.0_R8+1.03_R8*mustari**0.5_R8+0.31*mustari)*f1
     5              +0.59_R8*mustari*epsilon/
     6               (1_R8+0.74_R8*mustari*epsilon**1.5_R8)*
     7               (1_R8+
     8                 1.33_R8*alpha*(1.0_R8+0.60_R8*alpha)/
     9                 (1.0_R8+1.79_R8*alpha))*(f1-f2))
            write(*,*) 'chi: ',chi
            transport(iCv,is)=chi
            do i=2,mpg%ftCvP(iFt,2)
              transport(mpg%ftCv(mpg%ftCvP(iFt,1)+i-1),is)=chi
            enddo
          endif
        enddo
      enddo
      call subend ()
      return
      end subroutine old_chang_hinton
