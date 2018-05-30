      subroutine get_adas_data(reac,ns,nxd,nyd,nx,ny,ne,te,mach)
      use b2mod_types
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_constants
      IMPLICIT NONE
*..   input parameters (unchanged on exit)
      integer ns,nxd,nyd,nx,ny
      real (kind=R8) :: ne(0:nxd+1,0:nyd+1),te(0:nxd+1,0:nyd+1)
      CHARACTER REAC*9
*..   output data
      real (kind=R8) :: mach(0:nxd+1,0:nyd+1,0:ns-1)
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
      logical first_stage, last_stage
c     log10 of densities (in 1/cm3), temperatures (in eV)
      real (kind=R8) ::
     *   DDENSD(iddimd),DTEVD(itdimd),logT,logD
c     log10 of rate coefficient (in cm3/s)
      real (kind=R8) ::
     *   DREAC(isdimd,itdimd,iddimd),DTMP(itdimd,iddimd),ans
*..   local variables
      integer is, ic, ix, iy
      CHARACTER*2 ELNAME
      CHARACTER*200 DSN
*..   procedures
      logical, external :: open_file, streql
      external lower_case, xxdata_11, xerrab, subini, subend

      call subini ('get_adas_data')
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
            do ix=0,nx+1
              do iy=0,ny+1
                logT=log10(te(ix,iy)/qe)
                logD=log10(ne(ix,iy)/1.0e6_R8)
                call int2d(DTMP(1:itdimd,1:iddimd),itdimd,iddimd,
     1           DTEVD,ITMAXD,DDENSD,IDMAXD,
     2           logT,logD,ans)
                if(id.lt.10) then
                  if((ic.eq.0.and..not.first_stage) .or.
     &               (ic.eq.ncnct-1.and..not.last_stage)) then
                    mach(ix,iy,is+ic)=0.0_R8
                  else
                    mach(ix,iy,is+ic)=10.0_R8**(ans-6.0_R8)
                  end if
                else
                  mach(ix,iy,is+ic)=ans
                endif
              enddo
            enddo
          enddo
          is = is + ncnct
          close(29)
        else
          write(*,*) 'No ADAS file found : setting data to 0'
          do while (streql(ELNAME,is_codes(is)))
            mach(0:nx+1,0:ny+1,is) = 0.0_R8
            is = is + 1
          end do
        endif
      enddo
      call subend ()
      return
      end subroutine get_adas_data

      subroutine get_amjuel_data(filnam,h123,reac,crc,
     1 nxd,nyd,nx,ny,ne,te,mach)
      use b2mod_types
      use b2mod_constants
      IMPLICIT NONE
#     include <eir_reac.F>
      integer nxd,nyd,nx,ny
      real (kind=R8) ::
     1 ne(0:nxd+1,0:nyd+1),te(0:nxd+1,0:nyd+1),
     2 mach(0:nxd+1,0:nyd+1),D11(0:8,0:8)
      CHARACTER FILNAM*8,H123*4,REAC*9,CRC*3
      integer ir,i,j,ix,iy
      real (kind=R8) ::
     1 rmin,rmax,tef,def,dnnn,tei,dej

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

      rmin=1.0e30
      rmax=1.0e-30
      do iy=0,ny+1
        do ix=0,nx+1
          TEF=LOG(TE(ix,iy)/qe)
          DEF=LOG(ne(ix,iy)*1.0d-8*1.0d-6)
          dnnn=0.0
          DO I=0,8
            TEI=TEF**I
            if(
     1       h123.eq.'H.1'.or.
     1       h123.eq.'H.2'.or.
     1       h123.eq.'H.5'.or.
     1       h123.eq.'H.8'.or.
     1       h123.eq.'H.11') then
              DNNN=DNNN+D11(I,0)*TEI
            elseif(
     1       h123.eq.'H.3'.or.
     1       h123.eq.'H.4'.or.
     1       h123.eq.'H.6'.or.
     1       h123.eq.'H.7'.or.
     1       h123.eq.'H.9'.or.
     1       h123.eq.'H.10'.or.
     1       h123.eq.'H.12') then
              DO J=0,8
                DEJ=DEF**J
                DNNN=DNNN+D11(I,J)*TEI*DEJ
              enddo
            else
c             write(*,*) 'Do not yet know about ',h123
            endif
          enddo
          if(dnnn.gt.log(1.0e30)) then
            dnnn=1.0e30
          elseif(dnnn.lt.log(1.0e-30)) then
            dnnn=1.0e-30
          else
            DNNN=EXP(DNNN)
          endif
          rmin=min(rmin,dnnn)
          rmax=max(rmax,dnnn)
          mach(ix,iy)=dnnn
        enddo
      enddo
      write(*,*) 'rmin,rmax= ',rmin,rmax
      return
      end subroutine get_amjuel_data

      SUBROUTINE SLREAC_b2plot(IR,FILNAM,H123,REAC,CRC)
      use b2mod_types
c
C  input
C    filnam: read a&m data from file filnam, e.g. AMJUEL
c    ir    : store data on eirene array creac(...,...,ir)
c    H123  : identifier for reaction in filnam, e.g. H2
c    REAC  : number of reaction in filnam, e.g. 2.2.5
c    crc   : type of data process e.g. EI, CX, OT, etc
C  output
c    iswr  : eirene flag for type of process  (1,2,...7)
c    creac : eirene storage array for a&m data CREAC(9,0:9,IR)
c    modcol: see below
c
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
      real (kind=R8) ::
     1 CONST
      integer IR,I0,ISW,J,IND,I,K,IH,ireac
      CHARACTER ZEILE*80
      CHARACTER FILNAM*8,H123*4,REAC*9,CRC*3
      CHARACTER*4 CHR
      CHARACTER*3 CHRL,CHRR
      LOGICAL LCONST,LGEMIN,LGEMAX
      logical open_file
      external open_file, subini, subend, xerrab, xertst
C
      call subini('slreac_b2plot')
      call xertst(ir.le.nreac, 'SLREAC_b2plot: IR > NREAC')
      LGEMIN=.FALSE.
      LGEMAX=.FALSE.
      ISWR(IR)=0
      CONST=0.
      CHR=' l0 '
      I0=0
C
      IF (INDEX(CRC,'EI').NE.0) ISWR(IR)=1
      IF (INDEX(CRC,'DS').NE.0) ISWR(IR)=2
      IF (INDEX(CRC,'CX').NE.0) ISWR(IR)=3
      IF (INDEX(CRC,'II').NE.0) ISWR(IR)=4
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
      IREAC=INDEX(REAC,' ')-1
      IF (IREAC.LT.0) IREAC=LEN(REAC)
C
      IF (ISW.EQ.1) THEN
        CHR=' a0 '
        CHRL='al0'
        CHRR='ar0'
        I0=0
        MODCLF(IR)=MODCLF(IR)+1
      ELSEIF (ISW.EQ.2) THEN
        CHR=' b0 '
        CHRL='bl0'
        CHRR='br0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+10
      ELSEIF (ISW.EQ.3) THEN
        MODCLF(IR)=MODCLF(IR)+20
      ELSEIF (ISW.EQ.4) THEN
        MODCLF(IR)=MODCLF(IR)+30
      ELSEIF (ISW.EQ.5) THEN
        CHR=' e0 '
        CHRL='el0'
        CHRR='er0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+100
      ELSEIF (ISW.EQ.6) THEN
        MODCLF(IR)=MODCLF(IR)+200
      ELSEIF (ISW.EQ.7) THEN
        MODCLF(IR)=MODCLF(IR)+300
      ELSEIF (ISW.EQ.8) THEN
        CHR=' h0 '
        CHRL='hl0'
        CHRR='hr0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+1000
      ELSEIF (ISW.EQ.9) THEN
        MODCLF(IR)=MODCLF(IR)+2000
      ELSEIF (ISW.EQ.10) THEN
        MODCLF(IR)=MODCLF(IR)+3000
      ELSEIF (ISW.EQ.11) THEN
        CHR=' k0 '
        CHRL='kl0'
        CHRR='kr0'
        I0=1
      ELSEIF (ISW.EQ.12) THEN
      ENDIF
C
      IF (.NOT.LCONST) THEN
  100   READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,'##BEGIN DATA HERE##').EQ.0) GOTO 100

    1   READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,H123).EQ.0) GOTO 1
C
    2   READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,'H.').NE.0 .and.
     .      INDEX(ZEILE,'section').NE.0) GOTO 990
        IF (INDEX(ZEILE,'Reaction ').EQ.0.or.
     .      INDEX(ZEILE,REAC(1:ireac)).EQ.0) GOTO 2
      ENDIF
C
C  SINGLE PARAM. FIT, ISW=1,2,5,8,11
      IF (ISW.EQ.1.OR.ISW.EQ.2.OR.ISW.EQ.5.OR.ISW.EQ.8.OR.ISW.EQ.11)
     .THEN
        IF (.NOT.LCONST) THEN
    3     READ (29,'(A80)',END=990) ZEILE
          IF (INDEX(ZEILE,CHR).EQ.0) GOTO 3
          DO 9 J=0,2
            IND=0
            DO 4 I=1,3
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(2:2))
              READ (ZEILE((IND+2):80),'(E20.12)') CREAC(J*3+I,I0,IR)
    4       CONTINUE
            READ (29,'(A80)',END=990) ZEILE
    9     CONTINUE
          IF (INDEX(ZEILE,CHRL).NE.0) THEN
            IND=0
            DO 5 I=1,3
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(2:2))
              READ (ZEILE((IND+3):80),'(E20.12)') fparm(IR,i)
    5       CONTINUE
            LGEMIN=.true.
            READ (29,'(A80)',END=990) ZEILE
          ENDIF
          IF (INDEX(ZEILE,CHRR).NE.0) THEN
            IND=0
            DO 7 I=4,6
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(2:2))
              READ (ZEILE((IND+3):80),'(E20.12)') fparm(IR,i)
    7       CONTINUE
            LGEMAX=.true.
            READ (29,'(A80)',END=990) ZEILE
          ENDIF
c
          rcmn(ir)=-20.
          rcmx(ir)= 20.
          ifexmn(ir)=0
          ifexmx(ir)=0
c
          if (lgemin) then
            IND=INDEX(ZEILE,'=')
            READ (ZEILE((IND+2):80),'(E12.5)') rcmn(IR)
            rcmn(ir)=log(rcmn(ir))
            ifexmn(ir)=5
            READ (29,'(A80)',END=990) ZEILE
          endif
          if (lgemax) then
            IND=INDEX(ZEILE,'=')
            READ (ZEILE((IND+2):80),'(E12.5)') rcmx(IR)
            rcmx(ir)=log(rcmx(ir))
            ifexmx(ir)=5
            READ (29,'(A80)',END=990) ZEILE
          endif
C
        ELSE
          CREAC(1,I0,IR)=CONST
          DO 6 J=2,9
            CREAC(J,I0,IR)=0.
    6     CONTINUE
        ENDIF
C
      ELSE
C  TWO PARAM. FIT, ISW=3,4,6,7,9,10,12

        DO 11 J=0,2
   16     READ (29,'(A80)',END=990) ZEILE
          IF (INDEX(ZEILE,'Index').EQ.0) GOTO 16
          READ (29,'(1X)')
          DO 17 I=1,9
            READ (29,*) IH,(CREAC(I,K,IR),K=J*3+1,J*3+3)
cdr         write (6,*) IH,(CREAC(I,K,IR),K=J*3+1,J*3+3)
   17     CONTINUE
   11   CONTINUE
C
      ENDIF
C
      CLOSE (UNIT=29)
C
      call subend ()
      RETURN
C
  990 WRITE (6,*) ' NO DATA FOUND FOR REACTION ',H123,' ',REAC,
     .            ' IN DATA SET ',FILNAM
      WRITE (6,*) ' IR,MODCLF(IR) ',IR,MODCLF(IR)
      CLOSE (UNIT=29)
      call subend ()
      RETURN
  991 WRITE (6,*) ' INVALID CONSTANT IN SLREAC. CONST= ',CONST
      WRITE (6,*) ' CHECK "REACTION CARDS" FOR REACTION NO. ',IR
      CALL XERRAB ('Invalid constant in reaction card')
      END SUBROUTINE SLREAC_b2plot

      subroutine old_chang_hinton(nx,ny,ns,ismain,jsep,jxa,lambda,
     1 safetyfactor,cr,cz,crx,bb,ne,te,ti,na,transport)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer nx,ny,ns,ismain,jsep,jxa
      real (kind=R8) :: lambda,
     1 cr(-1:nx,-1:ny), cz(-1:nx,-1:ny), crx(-1:nx,-1:ny,0:3),
     2 bb(-1:nx,-1:ny,0:3),
     2 ne(-1:nx,-1:ny), te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     3 na(-1:nx,-1:ny,0:ns-1), safetyfactor(-1:ny),
     2 transport(-1:nx,-1:ny,0:ns-1)

      integer ix,iy,is,irmin,irmax
      real (kind=R8) :: crmin,crmax,czmin,czmax,
     1 cr2max,cr2min,cr1max,cr1min,
     2 bave,epsilon,lnlam,vthi,omegaci,
     1 rhoi,alpha,deltap,f1,f2,taui,nui,nustari,mustari,chi
      external subini,subend,xerrab

      call subini('old_chang_hinton')
      transport=0
      do iy=-1,jsep
        write(*,*) 'iy: ',iy
        ix=jxa
        if(.not.on_closed_surface(ix,iy)) then
          write(*,*) 'We have a problem: not in core!'
          call xerrab ('Chang-Hinton: Not in core!')
        endif
        irmin=jxa
        irmax=jxa
        crmin=cr(ix,iy)
        crmax=crmin
        czmin=cz(ix,iy)
        czmax=czmin
        ix=rightix(ix,iy)
        do while(ix.ne.jxa.and.ix.ge.-1.and.ix.le.nx)
          if(cr(ix,iy).lt.crmin) irmin=ix
          if(cr(ix,iy).gt.crmax) irmax=ix
          crmin=min(crmin,cr(ix,iy))
          crmax=max(crmax,cr(ix,iy))
          czmin=min(czmin,cz(ix,iy))
          czmax=max(czmax,cz(ix,iy))
          ix=rightix(ix,iy)
        enddo
        if(ix.ne.jxa) then
          if(nnreg(0).ge.2) then
            write(*,*) 'We have a problem: core not continuous!'
            call xerrab ('Chang-Hinton: core not continuous!')
          else
            write(*,*) 
     . 'Chang-Hinton calculation only available for closed field lines!'
            call subend ()
            return
          endif
        endif
        write(*,*) 'crmin,crmax,czmin,czmax',
     1   crmin,crmax,czmin,czmax
        bave=cr(ix,iy)*2.0_R8/(crmax+crmin)*bb(ix,iy,3)
        write(*,*) 'ne,te,ti,b,q: ',
     1   ne(ix,iy),te(ix,iy)/qe,ti(ix,iy)/qe,bave,safetyfactor(iy)
        epsilon=(crmax-crmin)/(crmax+crmin)
        write(*,*) 'epsilon: ',epsilon
cxpb    Coulomb logarithm formula from Wesson for ion-ion collisions
        if (lambda.lt.0.0_R8) then
          lnlam=max(-lambda,17.3_R8-
     1     log(ne(ix,iy)/1.0e20_R8)/2.0_R8+
     2     1.5_R8*log(ti(ix,iy)/ev/1000.0_R8))
        else
          lnlam = lambda
        endif
        write(*,*) 'lnlam: ',lnlam
        cr2max=0.5_R8*(crx(irmax,iy,1)+crx(irmax,iy,3))
        cr2min=0.5_R8*(crx(irmin,iy,1)+crx(irmin,iy,3))
        cr1max=0.5_R8*(crx(irmax,iy,0)+crx(irmax,iy,2))
        cr1min=0.5_R8*(crx(irmin,iy,0)+crx(irmin,iy,2))
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
            write(*,*) 'za,am,na: ',rza(ix,iy,is),am(is),na(ix,iy,is)
            vthi=sqrt(ti(ix,iy)/mp/am(is))
            write(*,*) 'vthi: ',vthi
            omegaci=rza(ix,iy,is)*qe*bave/mp/am(is)
            write(*,*) 'omegaci: ',omegaci
            rhoi=sqrt(2.0_R8)*vthi/omegaci
            write(*,*) 'rhoi: ',rhoi
            alpha=na(ix,iy,is)*rz2(ix,iy,is)/
     1            na(ix,iy,ismain)/rz2(ix,iy,ismain)  !!! think about this !!!
            write(*,*) 'alpha: ',alpha
            taui=12.0_R8*pi**1.5_R8*eps0**2*
     1             (am(is)*mp)**0.5_R8*ti(ix,iy)**1.5_R8/
     2             (na(ix,iy,is)*rz2(ix,iy,is)**2*qe**4*lnlam)
            write(*,*) 'taui: ',taui
            nui=1.0_R8/taui
            nustari=nui*(crmax+crmin)/2.0_R8*safetyfactor(iy)/
     1             epsilon**(1.5_R8)/vthi
            write(*,*) 'nustari: ',nustari
            mustari=nustari*(1.0_R8+1.54_R8*alpha)
            write(*,*) 'mustari: ',mustari
            chi=epsilon**(-1.5_R8)*
     1             safetyfactor(iy)**2*rhoi**2/taui*
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
            transport(ix,iy,is)=chi
            ix=rightix(ix,iy)
            do while(ix.ne.jxa.and.ix.ge.-1.and.ix.le.nx)
              transport(ix,iy,is)=chi
              ix=rightix(ix,iy)
            enddo
          endif
        enddo
      enddo
      call subend ()
      return
      end subroutine old_chang_hinton
