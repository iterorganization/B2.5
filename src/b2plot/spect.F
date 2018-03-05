c
      subroutine spect
      use b2mod_types
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      real (kind=R8) ::
     *   dtev((nx+2)*(ny+2)),
     1 ddens((nx+2)*(ny+2)),
     2 drcofi((nx+2)*(ny+2)),dspec(-1:nx,-1:ny,0:ns-1)
      integer ifail,iz0,iz1,iclass,itmax,ievcut,
     1 itdimd,itmaxd,idmaxd,izmaxd,
     2 ix,iy,is
      parameter(itdimd=40)
      real (kind=R8) ::
     *   dtevd(itdimd),ddensd(itdimd),zdata(itdimd),
     1 drcofd(itdimd,itdimd,itdimd)
      real (kind=R8) ::
     *   zmin,zmax
      character year*2, TITLF*80
      data    year   /'03'/

      itmax=0
      do iy=-1,ny
        do ix=-1,nx
          itmax=itmax+1
          dtev(itmax)=log10(te(ix,iy)/ev)
          ddens(itmax)=log10(ne(ix,iy))-6.0_R8
        enddo
      enddo

      do is=0,ns-1
        iz0=nint(zn(is))
        iz1=partition_number(is)+1  ! was nint(za(is))+1
        iclass=7
        ievcut=0

        CALL  DHDATA( YEAR   , TITLF  , IFAIL
     &   , IZ0    , IZ1    , ICLASS , ITMAX  , IEVCUT
     &   , ITDIMD , ITMAXD , IDMAXD , IZMAXD
     &   , DTEV   , DDENS
     &   , DTEVD  , DDENSD , DRCOFD , ZDATA
     &   , DRCOFI )

C
        if (ifail .gt. 0)   then
          write (*,3003) ifail,is,iz0,iz1
 3003     format (
     >     ' >>> error in dhdata: ifail =',i3,/,
     >     '  ifail = 1: no datafile for dhdata available',/,
     >     '        = 2: dimension itdimd too small',/,
     >     '  for species, nuclear charge, ionic charge =',3i3,/)
          do iy=-1,ny
            do ix=-1,nx
              dspec(ix,iy,is)=0.0_R8
            enddo
          enddo
        else
          itmax=0
          zmin=drcofi(1)
          zmax=drcofi(1)
          do iy=-1,ny
            do ix=-1,nx
              itmax=itmax+1
              zmin=min(zmin,drcofi(itmax))
              zmax=max(zmax,drcofi(itmax))
              dspec(ix,iy,is)=ne(ix,iy)*na(ix,iy,is)*
     1         1.0e-6_R8*10.0_R8**drcofi(itmax)/4.0_R8/pi
            enddo
          enddo
          write(*,*) 'is,zmin,zmax ',is,zmin,zmax
        endif

      enddo
      call copy(dspec,1.0_R8,ns,textpl)
      call customise_label('Diagnostic line emission from ADAS'//
     1 ' (photons.m^{-3}.sr^{-1}.s^{-1})',graphlabel(iwrk))
      return
      end subroutine spect


c     -------------------------------------------------------------------------
c     EMISSIVITY
c     -------------------------------------------------------------------------

cc    -------------------------------------------------------------------------
cc    SUBROUTINE InitADAS(IC,NC,ex,fn,ERR)
cc    INPUT:     INTEGER IC:        ion charge
cc               INTEGER NC:        nuclear charge
cc               CHARACTER*3 ex:    adf extension
cc               CHARACTER*80 fn:   file containing adf directories
cc    OUTPUT:    INTEGER ERR: 0=no error, 1=IC/NC-combination not supported
cc    -------------------------------------------------------------------------
cc    InitADAS selects the input file for the IC/NC-combination.
cc    -------------------------------------------------------------------------
cc    03/24/97 bhs
cc    -------------------------------------------------------------------------

      subroutine InitADAS(IC,NC,ex,fn,mie,mir,mic,ERR)

      use b2mod_types
      use b2mod_solpstop
      implicit none

#     include "emissivity.dimensions"
#     include "adas_common.inc"

      INTEGER IC,NC,ERR,counter,MAXNC
      CHARACTER ex*3,fn*80,MIE*2,MIR*2,MIC*2
      logical open_file

      if(.not.open_file(99,fn,
     1 'data.local/adas|modules/adas')) return
      read(99,*) UIDIN
      if(uidin(1:8).eq.'SOLPSTOP') then
        if (uidin(9:9).ne.'/') then
          uidin=trim(solpstop()) // '/' // trim(uidin(9:))
        else
          uidin=trim(solpstop()) // '/' // trim(uidin(10:))
        endif
        write(*,*) 'ADAS UIDIN changed to ',trim(UIDIN)
      endif
      read(99,*) MAXNC
      if (MAXNC .lt. NC) then
         write(*,*) 'InitADAS: Specified data not supported'
         ERR=1
         return
      endif
      TYPIN=' '
      do counter=1,NC
         read(99,*) GRPIN,ex,mie,mir,mic
      enddo
      close(99)
      EXTIN=ex

      if (GRPIN  .eq.  'none    ') then
         write(*,*) 'InitADAS: Specified data not supported'
         ERR=1
         return
      endif
      call XXUID(UIDIN)
      call XXSPEC(GRPIN,TYPIN,EXTIN)
      IZIN=IC
      IZ0IN=NC
      ERR=0

      return
      end subroutine InitADAS

cc    -------------------------------------------------------------------------
cc    SUBROUTINE GetIBSEL(WL,DWL,MIE,MIR,MIC,id,IBSELAE,IBSELAR,IBSELAC)
cc    INPUT:     REAL WL:             wavelength in A
cc               REAL DWL:            half width of wavelength window in A
cc               CHARACTER MIE*2:     metastable index excitation
cc               CHARACTER MIR*2:     metastable index recombination
cc               CHARACTER MIC*2:     metastable index charge exchange
cc               INTEGER id:          dimension for array of indices
cc    OUTPUT:    INTEGER IBSELAE(id): array of indices, excit
cc               INTEGER IBSELAR(id): array of indices, recom
cc               INTEGER IBSELAC(id): array of indices, cexch
cc    -------------------------------------------------------------------------
cc    GetIBSEL selects all data blocks in the wavelength window [WL-DWL,WL+DWL]
cc    InitADAS MUST be called before any call to this subroutine!
cc    -------------------------------------------------------------------------
cc    03/24/97 bhs
cc    -------------------------------------------------------------------------

      subroutine GetIBSEL(WL,DWL,MIE,MIR,MIC,id,
     &                    IBSELAE,IBSELAR,IBSELAC)

      use b2mod_types
      implicit none

#     include "emissivity.dimensions"
#     include "adas_common.inc"

      real (kind=R8) ::
     *   WL,DWL
      character MIE*2,MIR*2,MIC*2
      integer id
      integer IBSELAE(id)
      integer IBSELAR(id)
      integer IBSELAC(id)
      integer count,countEX,countRE,countCE
      integer NBSEL
      character CWAVEL(NSTORE)*10,CTYPE(NSTORE)*8
      character CINDM(NSTORE)*2

      call GetInfo(NBSEL,CWAVEL,CTYPE,CINDM,LOT)
      do count=1,id
         IBSELAE(count)=0
         IBSELAR(count)=0
         IBSELAC(count)=0
      enddo
      countEX=1
      countRE=1
      countCE=1
      do count=1,nbsel
         read(CWAVEL(count),'(f8.1)') wla(count)
         write(*,*) count,wla(count),CWAVEL(count)
         if (abs(wl-wla(count)) .le. dwl) then
            if (CTYPE(count) .eq. 'EXCIT') then
               if (MIE .eq. CINDM(count)) then
                  IBSELAE(countEX)=count
                  countEX=countEX+1
                  write(*,*)
     &            'GetIBSEL: selected transition excitation'
                  write(*,*) 'wavelength:',CWAVEL(count)
                  write(*,*) 'metastable index: ',CINDM(count)
               endif
            endif
            if (CTYPE(count) .eq. 'RECOM') then
               if (MIR .eq. CINDM(count)) then
                  IBSELAR(countRE)=count
                  countRE=countRE+1
                  write(*,*)
     &            'GetIBSEL: selected transition recombination'
                  write(*,*) 'wavelength:',CWAVEL(count)
                  write(*,*) 'metastable index: ',CINDM(count)
               endif
            endif
            if (CTYPE(count) .eq. 'CHEXC') then
               if (MIC .eq. CINDM(count)) then
                  IBSELAC(countCE)=count
                  countCE=countCE+1
                  write(*,*)
     &            'GetIBSEL: selected transition charge exchange'
                  write(*,*) 'wavelength:',CWAVEL(count)
                  write(*,*) 'metastable index: ',CINDM(count)
               endif
            endif
         endif
      enddo
      if ((countEX .eq. 1) .and. (countRE .eq. 1) .and.
     1 (countCE .eq. 1)) then
         write(*,*) 'GetIBSEL: none selected for'
         write(*,*) 'ion charge: ',IZIN,' nuclear charge: ',IZ0IN
      endif
      return
      end subroutine GetIBSEL

cc    -------------------------------------------------------------------------
cc    SUBROUTINE InitIBSEL(id,IBSELAE,IBSELAR,IBSELAC)
cc    INPUT:     INTEGER id:          dimension for array of indices
cc    OUTPUT:    INTEGER IBSELAE(id): array of indices, excit
cc               INTEGER IBSELAR(id): array of indices, recom
cc               INTEGER IBSELAC(id): array of indices, cexch
cc    -------------------------------------------------------------------------
cc    InitIBSEL selects no data block indices, but sets the wavelengths in
cc    adas_common.inc.
cc    -------------------------------------------------------------------------
cc    03/25/97 bhs
cc    -------------------------------------------------------------------------

      subroutine InitIBSEL(id,IBSELAE,IBSELAR,IBSELAC)

      use b2mod_types
      implicit none

#     include "emissivity.dimensions"
#     include "adas_common.inc"

      integer id
      integer IBSELAE(id)
      integer IBSELAR(id)
      integer IBSELAC(id)
      integer count
      integer NBSEL
      character CWAVEL(NSTORE)*10,CTYPE(NSTORE)*8
      character CINDM(NSTORE)*2

      call GetInfo (NBSEL,CWAVEL,CTYPE,CINDM,LOT)
      do count=1,id
         IBSELAE(count)=0
         IBSELAR(count)=0
         IBSELAC(count)=0
      enddo
      do count=1,nbsel
         read(CWAVEL(count),*) wla(count)
      enddo

      return
      end subroutine InitIBSEL


cc    -------------------------------------------------------------------------
cc    SUBROUTINE GetEmissivity(id,IBSELAE,IBSELAR,IBSELAC,NE,NI,NI1,NH,TE,dim,
cc                             E,P)
cc    INPUT:     INTEGER id:          dimension of array of indices
cc               INTEGER IBSELAE(id): array of the selected data indices, excit
cc               INTEGER IBSELAR(id): array of the selected data indices, recom
cc               INTEGER IBSELAC(id): array of the selected data indices, cexch
cc               REAL NE(dim),NI(dim),densities in 1/cm3
cc                    NI1(dim):       electrons,ions,ions higher charge state
cc               REAL NH(dim):        density of neutral H in 1/cm3
cc               TE(dim):             electron temperature in eV
cc               INTEGER dim:         dimension for arrays
cc    OUTPUT:    REAL E(dim):         emissivity in photons/cm3
cc               REAL P(dim):         emissivity in W/cm3
cc    -------------------------------------------------------------------------
cc    GetEmissivity calculates the emissivity of the transitions selected in
cc    IBSELAE, IBSELAR and IBSELAC (see subroutine GetIBSEL)
cc    -------------------------------------------------------------------------
cc    03/24/97 bhs
cc    -------------------------------------------------------------------------

      subroutine GetEmissivity(id,IBSELAE,IBSELAR,IBSELAC,
     &                         cNE,cNI,cNI1,cNH,cTE,dim,cE,cPP)

      use b2mod_types
      use b2mod_constants
      implicit none

#     include "emissivity.constants"
#     include "emissivity.dimensions"
#     include "adas_common.inc"

      integer dim
      real (kind=R8) ::
     *    NE(maxdim),NI(maxdim),NI1(maxdim),
     &       NH(maxdim),TE(maxdim),E(maxdim),
     &       locTE(maxdim)
      real (kind=R8) ::
     *   cNE(dim),cNI(dim),cNI1(dim),cNH(dim),cTE(dim)
      real (kind=R8) ::
     *   cE(dim),cPP(dim)
      integer id
      integer IBSELAE(id),IBSELAR(id),IBSELAC(id)
      integer IRCODE
      character TITLX*296
      real (kind=R8) ::
     *    WL,PECA(maxdim)
      logical LTRNG(maxdim),LDRNG(maxdim)
      real (kind=R8) ::
     *    PECE(maxdim),PECR(maxdim),PECC(maxdim)
      real (kind=R8) ::
     *    PE(maxdim),PR(maxdim),PC(maxdim)
      real (kind=R8) ::
     *   f
      integer index
      integer count,cc
      external subini, subend

      call subini ('GetEmissivity')
      do count=1,dim
         NE(count)=cNE(count)
         NI(count)=cNI(count)
         NI1(count)=cNI1(count)
         NH(count)=cNH(count)
         TE(count)=cTE(count)
         PECE(count)=0.0_R8
         PECR(count)=0.0_R8
         PECC(count)=0.0_R8
         PE(count)=0.0_R8
         PR(count)=0.0_R8
         PC(count)=0.0_R8
      enddo
      f=const_h*c/1.0e-10_R8

      count=1
      do while (IBSELAE(count) .ne. 0)
         index=IBSELAE(count)
         do cc=1,dim
            if (TE(cc) .lt. LOT(index)) then
               locTE(cc)=LOT(index)
            else
               locTE(cc)=TE(cc)
            endif
         enddo
         call SPEC(index,IZIN,IZ0IN,dim,locTE,NE,
     &             WL,PECA(1:dim),LTRNG,LDRNG,TITLX,IRCODE)
         do cc=1,dim
            PECE(cc)=PECE(cc)+PECA(cc)
            PE(cc)=PE(cc)+PECA(cc)*NE(cc)*NI(cc)*f/wla(index)
         enddo
         count=count+1
      enddo
      count=1
      do while (IBSELAR(count) .ne. 0)
         index=IBSELAR(count)
c         do cc=1,dim
c            if (TE(cc) .lt. LOT(index)) then
c               locTE(cc)=LOT(index)
c           else
c               locTE(cc)=TE(cc)
c           endif
c         enddo
         call SPEC(index,IZIN,IZ0IN,dim,locTE,NE,
     &             WL,PECA(1:dim),LTRNG,LDRNG,TITLX,IRCODE)
         do cc=1,dim
            PECR(cc)=PECR(cc)+PECA(cc)
            PR(cc)=PR(cc)+PECA(cc)*NE(cc)*NI1(cc)*f/wla(index)
         enddo
         count=count+1
      enddo
      count=1
      do while (IBSELAC(count) .ne. 0)
         index=IBSELAC(count)
c         do cc=1,dim
c            if (TE(cc) .lt. LOT(index)) then
c               locTE(cc)=LOT(cc)
c            else
c               locTE(cc)=TE(cc)
c            endif
c         enddo
         call SPEC(index,IZIN,IZ0IN,dim,locTE,NE,
     &             WL,PECA(1:dim),LTRNG,LDRNG,TITLX,IRCODE)
         do cc=1,dim
            PECC(cc)=PECC(cc)+PECA(cc)
            PC(cc)=PC(cc)+PECA(cc)*NH(cc)*NI1(cc)*f/wla(index)
         enddo
         count=count+1
      enddo

      do cc=1,dim
         E(cc)=PECE(cc)*NE(cc)*NI(cc)+PECR(cc)*NE(cc)*NI1(cc)+
     &        PECC(cc)*NH(cc)*NI1(cc)
         cE(cc)=E(cc)
         cPP(cc)=PE(cc)+PR(cc)+PC(cc)
      enddo

      call subend ()
      return
      end subroutine GetEmissivity
c     -------------------------------------------------------------------------
c     subroutine GetInfo(NBSEL,CWAVEL,CTYPE,CINDM,T)
c     -------------------------------------------------------------------------
c     OUTPUT: integer NBSEL
c             character CWAVEL(NSTORE)*10
c             character CTYPE (NSTORE)*8
c             character CINDM (NSTORE)*2
c             real T(NSTORE)
c     -------------------------------------------------------------------------
c     InitADAS must have been called before any call to this subroutine.
c     GetInfo provides the following information about the selected input
c     file: number of data blocks NBSEL, wavelengths CWAVEL, data types
c     CTYPE, metastable indices CINDM, lowest temperature in eV T
c     -------------------------------------------------------------------------
c     04/02/97 bhs
c     -------------------------------------------------------------------------

      subroutine GetInfo(NBSEL,CWAVEL,CTYPE,CINDM,T)

      use b2mod_types
      implicit none

#     include "emissivity.dimensions"
#     include "adas_common.inc"

      real (kind=R8) ::
     *    TETA(NTDIM,NSTORE),TEDA(NDDIM,NSTORE)
      real (kind=R8) ::
     *    PEC(NTDIM,NDDIM,NSTORE)
      real (kind=R8) ::
     *   T(NSTORE)
      integer IUNIT,IRCODE
      integer IZ0,IZ,IZ1
      integer NBSEL
      integer ISELA(NSTORE),ITA(NSTORE),IDA(NSTORE)
      integer c
      character CWAVEL(NSTORE)*10
      character CFILE(NSTORE)*8,CTYPE(NSTORE)*8
      character CINDM(NSTORE)*2
      character DSNREQ*256
      character ESYM*2
      external subini, subend

c     initialising
      call subini ('GetInfo')
      IUNIT=4950

c     reading data
      call xxuid(UIDIN)
      call xxspec(GRPIN,TYPIN,EXTIN)
      call e3file(IUNIT,IZ0IN,IZIN,IRCODE,DSNREQ)
      if (IRCODE .ne. 0) then
         write(*,*) 'GetInfo: ERROR'
         return
      endif
      call e3data(IUNIT,DSNREQ,NSTORE,NTDIM,NDDIM,
     &            IZ0,IZ,IZ1,ESYM,NBSEL,ISELA,
     &            CWAVEL,CFILE,CTYPE,CINDM,
     &            ITA,IDA,TETA,TEDA,PEC)
      close(IUNIT)
      do c=1,NBSEL
         T(c)=TETA(1,c)
      enddo

      call subend ()
      return
      end subroutine GetInfo
C  UNIX PORT - SCCS Info : Module @(#)e3spln.for        1.3  Date 4/22/96
C
      SUBROUTINE E3SPLN( NTDIM  , NDDIM  ,
     &                   ITA    , IDA    , ITVAL   ,
     &                   TETA   , TEDA   , TEVA    , DIN    ,
     &                   PEC    ,          PECA    ,
     &                                     LTRNG   , LDRNG
     &                 )
C
      use b2mod_types
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: E3SPLN *********************
C
C  PURPOSE:
C          PERFORMS CUBIC SPLINE ON LOG(TEMPERATURE AND DENSITY)
C          VERSUS LOG(IONISATIONS PER PHOTON)
C          INPUT DATA FOR A GIVEN WAVELENGTH DATA-BLOCK.
C
C          USING  TWO-WAY SPLINES IT CALCULATES  THE  PHOTON EMISSIVITY
C          FOR  'ITVAL'  PAIRS OF  ELECTRON TEMPERATURES  AND  DENSITIES
C          FROM THE TWO-DIMENSIONAL TABLE OF TEMPERATURES/DENSITIES READ
C          IN FROM THE INPUT FILE. IF A  VALUE  CANNOT  BE  INTERPOLATED
C          USING SPLINES IT IS EXTRAPOLATED VIA 'XXSPLE'.
C
C  CALLING PROGRAM: ADAS503/SPEC
C
C
C  SUBROUTINE:
C
C  INPUT : (I*4)  NTDIM   = MAX NUMBER OF ELECTRON TEMPERATURES ALLOWED
C  INPUT : (I*4)  NDDIM   = MAX NUMBER OF ELECTRON DENSITIES    ALLOWED
C
C  INPUT : (I*4)  ITA     = INPUT DATA FILE: NUMBER OF ELECTRON TEMPERA-
C                           TURES READ FOR THE DATA-BLOCK BEING ASSESSED
C  INPUT : (I*4)  IDA     = INPUT DATA FILE: NUMBER OF ELECTRON DENSIT-
C                           IES   READ FOR THE DATA-BLOCK BEING ASSESSED
C  INPUT : (I*4)  ITVAL   = NUMBER OF ISPF ENTERED TEMPERATURE/DENSITY
C                           PAIRS FOR WHICH IONISATIONS PER PHOTON
C                           ARE REQUIRED FOR TABULAR/GRAPHICAL OUTPUT.
C
C  INPUT : (R*8)  TETA()  = INPUT DATA FILE: ELECTRON TEMPERATURES (EV)
C                           FOR THE DATA-BLOCK BEING ASSESSED.
C                           DIMENSION: ELECTRON TEMPERATURE INDEX
C  INPUT : (R*8)  TEDA()  = INPUT DATA FILE: ELECTRON DENSITIES (CM-3)
C                           FOR THE DATA-BLOCK BEING ASSESSED.
C                           DIMENSION: ELECTRON DENSITY INDEX
C  INPUT : (R*8)  TEVA()  = USER ENTERED: ELECTRON TEMPERATURES (EV)
C                           DIMENSION: TEMPERATURE/DENSITY PAIR INDEX
C  INPUT : (R*8)  DIN()   = USER ENTERED: ELECTRON DENSITIES (CM-3)
C                           DIMENSION: TEMPERATURE/DENSITY PAIR INDEX
C
C
C  INPUT : (R*8)  PEC(,)   =INPUT DATA FILE: FULL SET OF IONISATIONS PER
C                           PHOTON VALUES FOR THE DATA-BLOCK BEING
C                           ANALYSED.
C                           1ST DIMENSION: ELECTRON TEMPERATURE INDEX
C                           2ND DIMENSION: ELECTRON DENSITY     INDEX
C  OUTPUT: (R*8)  PECA()  = SPLINE INTERPOLATED OR  EXTRAPOLATED  IONIZ-
C                           ATIONS/PHOTON FOR THE USER ENTERED ELECTRON
C                           TEMPERATURE/DENSITY PAIRS.
C                           DIMENSION: TEMPERATURE/DENSITY PAIR INDEX
C
C  OUTPUT: (L*4)  LTRNG()=  .TRUE.  => OUTPUT 'PECA()' VALUE WAS INTER-
C                                      POLATED  FOR  THE  USER  ENTERED
C                                      ELECTRON TEMPERATURE 'TEVA()'.
C                           .FALSE. => OUTPUT 'PECA()' VALUE WAS EXTRA-
C                                      POLATED  FOR  THE  USER  ENTERED
C                                      ELECTRON TEMPERATURE 'TEVA()'.
C                           DIMENSION: TEMPERATURE/DENSITY PAIR INDEX
C
C  OUTPUT: (L*4)  LDRNG()=  .TRUE.  => OUTPUT 'PECA()' VALUE WAS INTER-
C                                      POLATED  FOR  THE  USER  ENTERED
C                                      ELECTRON DENSITY 'DIN()'.
C                           .FALSE. => OUTPUT 'PECA()' VALUE WAS EXTRA-
C                                      POLATED  FOR  THE  USER  ENTERED
C                                      ELECTRON DENSITY 'DIN()'.
C                           DIMENSION: TEMPERATURE/DENSITY PAIR INDEX
C
C          (I*4)  NIN     = PARAMETER = MAX. NO. OF INPUT  TEMP/DENSITY
C                                       VALUES. MUST BE >= 'ITA'&'IDA'
C          (I*4)  NOUT    = PARAMETER = MAX. NO. OF OUTPUT TEMP/DENSITY
C                                       PAIRS.  MUST BE >= 'ITVAL'
C          (I*4)  L1      = PARAMETER = 1
C
C          (I*4)  IED     = ARRAY SUBSCRIPT USED INPUT  FILE  ELECTRON
C                           DENSITIES.
C          (I*4)  IET     = ARRAY SUBSCRIPT USED INPUT  FILE  ELECTRON
C                           TEMPERATURES.
C          (I*4)  IT      = ARRAY  SUBSCRIPT  USED  FOR  USER  ENTERED
C                           TEMPERATURE/DENSITY PAIRS .
C          (I*4)  IOPT    = DEFINES THE BOUNDARY DERIVATIVES FOR THE
C                           SPLINE ROUTINE 'XXSPLE', SEE 'XXSPLE'.
C                           (VALID VALUES = <0, 0, 1, 2, 3, 4)
C
C          (L*4)  LSETX   = .TRUE.  => SET UP SPLINE PARAMETERS RELATING
C                                      TO 'XIN' AXIS.
C                           .FALSE. => DO NOT SET UP SPLINE PARAMETERS
C                                      RELATING TO 'XIN' AXIS.
C                                      (I.E. THEY WERE SET IN A PREVIOUS
C                                            CALL )
C                           (VALUE SET TO .FALSE. BY 'XXSPLE')
C
C          (R*8)  R8FUN1  = FUNCTION - (SEE ROUTINES SECTION BELOW)
C
C          (R*8)  XIN()   = 1) LOG( DATA FILE ELECTRON DENSITIES    )
C                           2) LOG( DATA FILE ELECTRON TEMPERATURES )
C          (R*8)  YIN()   = LOG( DATA FILE IONISATIONS/PHOTON )
C          (R*8)  XOUT()  = 1) LOG( SCALED USER ENTERED ELECTRON DENS. )
C                           2) LOG( SCALED USER ENTERED ELECTRON TEMPS.)
C          (R*8)  YOUT()  = LOG( OUTPUT GENERATED IONISATIONS/PHOTON )
C          (R*8)  YPASS(,)= LOG( IONISATIONS/PHOTON) INTERMEDIATE ARRAY
C                           WHICH   STORES   INTERPOLATED/EXTRAPOLATED
C                           VALUES  BETWEEN  THE  TWO  SPLINE SECTIONS.
C                           SECTIONS.
C          (R*8)  DF()    = SPLINE INTERPOLATED DERIVATIVES
C
C
C NOTE:
C
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C          XXSPLE     ADAS      SPLINE SUBROUTINE (EXTENDED DIAGNOSTICS)
C          R8FUN1     ADAS      REAL*8 FUNCTION: ( X -> X )
C
C AUTHOR:  H. P. SUMMERS
C          K1/1/57
C          JET EXT. 4941
C
C DATE:    11/10/91
C UNIX-IDL PORT:
C
C VERSION: 1.1                          DATE: 28-02-95
C MODIFIED: LALIT JALOTA (TESSELLA SUPPORT SERVICES PLC)
C               - FIRST VERSION
C
C VERSION: 1.2                          DATE: 15-04-96
C MODIFIED: TIM HAMMOND (TESSELLA SUPPORT SERVICES PLC)
C               - INCREASED PARAMETER NIN 24 -> 26
C                 (CURRENT JET VALUE 35, BUT KEPT TO 26 FOR PRESENT TO
C                  PREVENT CHANGES TO IDL ROUTINES).
C
C VERSION: 1.3                          DATE: 22-04-96
C MODIFIED: TIM HAMMOND
C               - INCREASED NIN TO 35 IN LINE WITH JET.
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      INTEGER    NIN                   , NOUT            , L1
C-----------------------------------------------------------------------
      PARAMETER( NIN = 35              , NOUT = 20       , L1 = 1     )
C-----------------------------------------------------------------------
      INTEGER    NTDIM                 , NDDIM           ,
     &           ITA                   , IDA             , ITVAL
      INTEGER    IET                   , IED             , IT         ,
     &           IOPT
C-----------------------------------------------------------------------
      real (kind=R8) ::
     *        R8FUN1
C-----------------------------------------------------------------------
      LOGICAL    LSETX
C-----------------------------------------------------------------------
      real (kind=R8) ::
     *        TETA(ITA)             , TEDA(IDA)       ,
     &           TEVA(ITVAL)           , DIN(ITVAL)      ,
     &           PECA(ITVAL)           ,
     &           PEC(NTDIM,NDDIM)
      real (kind=R8) ::
     *        DF(NIN)               ,
     &           XIN(NIN)              , YIN(NIN)        ,
     &           XOUT(NOUT)            , YOUT(NOUT)      ,
     &                                   YPASS(NOUT,NIN)
C-----------------------------------------------------------------------
      LOGICAL    LTRNG(ITVAL)          , LDRNG(ITVAL)
C-----------------------------------------------------------------------
      EXTERNAL   R8FUN1
      external   subini, subend, xertst
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      call subini ('e3spln')
      call xertst (NIN .ge. IDA,
     &             ' E1SPLN ERROR: NIN < IDA - INCREASE NIN')
      call xertst (NIN .ge. ITA,
     &             ' E1SPLN ERROR: NIN < ITA - INCREASE NIN')
      call xertst (NOUT .ge. ITVAL,
     &             ' E1SPLN ERROR: NOUT < ITVAL - INCREASE NTOUT')
C
C-----------------------------------------------------------------------
C SET UP FIRST SPLINE BOUNDARY CONDITIONS - SWITCH ON EXTRAPOLATION.
C-----------------------------------------------------------------------
C
      LSETX = .TRUE.
      IOPT  = 0
C
C-----------------------------------------------------------------------
C SET UP ARRAYS CONTAINING LOG VALUES OF INPUT/OUTPUT ELECTRON TEMPS.
C-----------------------------------------------------------------------
C
         DO 1 IET=1,ITA
            XIN(IET) = LOG( TETA(IET) )
    1    CONTINUE
C
         DO 2 IT=1,ITVAL
            XOUT(IT) = LOG( TEVA(IT)  )
    2    CONTINUE
C
C-----------------------------------------------------------------------
C SPLINE OVER ALL DATASET ELECTRON DENSITIES FOR EACH USER ELECTRON TEMP
C-----------------------------------------------------------------------
C
         DO 3 IED=1,IDA
C
               DO 4 IET=1,ITA
                  YIN(IET) = LOG( PEC(IET,IED) )
    4          CONTINUE
C
            CALL XXSPLE( LSETX , IOPT    , R8FUN1 ,
     &                   ITA   , XIN     , YIN    ,
     &                   ITVAL , XOUT    , YOUT   ,
     &                   DF    , LTRNG
     &                 )
C
               DO 5 IT=1,ITVAL
                  YPASS(IT,IED) = YOUT(IT)
    5          CONTINUE
C
    3    CONTINUE
C
C-----------------------------------------------------------------------
C SET UP SECOND SPLINE BOUNDARY CONDITIONS - SWITCH ON EXTRAPOLATION.
C-----------------------------------------------------------------------
C
      LSETX = .TRUE.
      IOPT  = 0
C
C-----------------------------------------------------------------------
C SET UP ARRAYS CONTAINING LOG VALUES OF INPUT/OUTPUT ELECTRON DENSITIES
C-----------------------------------------------------------------------
C
         DO 6 IED=1,IDA
            XIN(IED) = LOG( TEDA(IED) )
    6    CONTINUE
C
         DO 7 IT=1,ITVAL
            XOUT(IT) = LOG( DIN(IT)   )
    7    CONTINUE
C
C-----------------------------------------------------------------------
C SPLINE OVER ALL USER ELECTRON TEMPERATURES FOR EACH USER ELECTRON DENS
C-----------------------------------------------------------------------
C
         DO 8 IT=1,ITVAL
C
               DO 9 IED=1,IDA
                  YIN(IED) = YPASS(IT,IED)
    9          CONTINUE
C
            CALL XXSPLE( LSETX , IOPT       , R8FUN1   ,
     &                   IDA   , XIN        , YIN      ,
     &                   L1    , XOUT(IT)   , YOUT(IT) ,
     &                   DF    , LDRNG(IT)
     &                 )
C
    8    CONTINUE
C
C-----------------------------------------------------------------------
C SET UP OUTPUT IONISATIONS PER PHOTON ARRAY.
C-----------------------------------------------------------------------
C
         DO 10 IT=1,ITVAL
            PECA(IT) = EXP( YOUT(IT) )
   10    CONTINUE
C
C-----------------------------------------------------------------------
C
      call subend ()
      RETURN
      END SUBROUTINE E3SPLN
CX UNIX PORT - SCCS info: Module @(#)spec.for   1.2 Date 4/22/96
CX
      SUBROUTINE SPEC( IBSEL  , IZIN   , IZ0IN  ,
     &                 ITVAL  , TVAL   , DVAL   ,
     &                 WLNGTH ,
     &                 PECA   , LTRNG  , LDRNG  ,
     &                 TITLX  , IRCODE
     &               )
      use b2mod_types
      IMPLICIT NONE

C-----------------------------------------------------------------------
C
C  ******************* FORTRAN77 SUBROUTINE: SPEC **********************
C
C  PURPOSE: TO EXTRACT AND INTERPOLATE  PHOTON EMISSIVITIES FOR
C           EMITTING IONS.
C           USES THE SAME ROUTINES USED BY ADAS503, EXCEPT FOR:
C
C           'E3FILE' - WHICH OPENS THE REQUESTED FILE.
C           'E3CHKB' - WHICH CHECKS INPUT  VALUES  ARE  CONSISTENT  WITH
C                      THE SELECTED DATA-BLOCK 'IBSEL' AND   'IBSEL'  IS
C                      IN RANGE.
C
C            THE FIRST OF THESE FUNCTIONS IS CARRIED  OUT  IN  'ADAS503'
C            VIA ISPF PANELS USING THE ROUTINE 'E3SPF0'  -  ADAS503 DOES
C            NOT REQUIRE THE ROUTINE 'E3CHKB' AS THE USER CANNOT  SELECT
C            AN INVALID VALUE FOR 'IBSEL' OR 'IBSEL'/EMITTER COMBINATION
C
C
C  CALLING PROGRAM: GENERAL USE
C
C  SUBROUTINE:
C
C  INPUT : (I*4)   IBSEL   = INDEX OF DATA-BLOCK SELECTED FOR ANALYSIS
C  INPUT : (I*4)   IZIN    = ION CHARGE OF EMITTING ION
C  INPUT : (I*4)   IZ0IN   = NUCLEAR CHARGE OF EMITTING ION
C
C  INPUT : (I*4)   ITVAL   = NO. OF ELECTRON TEMPERATURE/DENSITY PAIRS
C  INPUT : (R*8)   TVAL()  = ELECTRON TEMPERATURES (UNITS: EV)
C                            DIMENSION: TEMPERATURE/DENSITY PAIR INDEX
C  INPUT : (R*8)   DVAL()  = ELECTRON DENSITIES (UNITS: CM-3)
C                            DIMENSION: TEMPERATURE/DENSITY PAIR INDEX
C
C  OUTPUT: (R*8)   WLNGTH  = SELECTED BLOCK WAVELENGTH (ANGSTROMS)
C
C  OUTPUT: (R*8)   PECA()  = PHOTON EMISSIVITIES.
C                            DIMENSION: TEMPERATURE/DENSITY PAIR INDEX
C  OUTPUT: (L*4)   LTRNG() =.TRUE.  => OUTPUT 'PECA()'  VALUE WAS INTER-
C                                      POLATED  FOR  THE  USER  ENTERED
C                                      ELECTRON TEMPERATURE 'TVAL()'.
C                           .FALSE. => OUTPUT 'PECA()'  VALUE WAS EXTRA-
C                                      POLATED  FOR  THE  USER  ENTERED
C                                      ELECTRON TEMPERATURE 'TVAL()'.
C                            DIMENSION: TEMPERATURE/DENSITY PAIR INDEX
C  OUTPUT: (L*4)   LDRNG() =.TRUE.  => OUTPUT 'PECA()'  VALUE WAS INTER-
C                                      POLATED  FOR  THE  USER  ENTERED
C                                      ELECTRON DENSITY 'DVAL()'.
C                           .FALSE. => OUTPUT 'PECA()'  VALUE WAS EXTRA-
C                                      POLATED  FOR  THE  USER  ENTERED
C                                      ELECTRON DENSITY 'DVAL()'.
C                            DIMENSION: TEMPERATURE/DENSITY PAIR INDEX
C
C  OUTPUT: (C*296) TITLX   = INFORMATION STRING (DSN ETC.)
C  OUTPUT: (I*4)   IRCODE  = RETURN CODE FROM SUBROUTINE:
C                            0 => NORMAL COMPLETION - NO ERROR DETECTED
C                            1 => DATA SET MEMBER FOR EMITTING ION WITH
C                                 CHARGE 'IZIN' & ION CHARGE 'IZ0IN' CAN
C                                 NOT BE FOUND/DOES NOT EXIST.
C                            2 => DISCREPANCY BETWEEN REQUESTED CHARGES
C                                 AND THOSE IN INPUT FILE.
C                            3 => THE SELECTED DATA-BLOCK 'IBSEL' IS OUT
C                                 OF RANGE OR DOES NOT EXIST.
C                            4 => INVALID VALUE FOR 'IZ0IN' ENTERED.
C                                 ('IZ0MIN' <= 'IZ0IN' <= 'IZ0MAX')
C                            5 => INVALID VALUE FOR 'IZIN' ENTERED.
C                                 ( 0  <= 'IZIN' <= 99 )
C                            9 => ERROR ENCOUNTERED WHEN TRYING TO OPEN
C                                 INPUT DATA-SET.
C
C          (I*4)   NSTORE  = PARAMETER= MAXIMUM NUMBER  OF  DATA-BLOCKS
C                                      WHICH CAN BE READ FROM THE INPUT
C                                      DATA-SET.
C          (I*4)   NTDIM   = PARAMETER= MAXIMUM NUMBER OF ELECTRON TEMP-
C                                      ERATURES THAT CAN BE READ  FROM
C                                      AN INPUT DATA-SET DATA-BLOCK.
C          (I*4)   NDDIM   = PARAMETER= MAXIMUM NUMBER OF ELECTRON DENS-
C                                      ITIES  THAT  CAN  BE  READ  FROM
C                                      AN INPUT DATA-SET DATA-BLOCK.
C          (I*4)   IZ0MIN  = PARAMETER: MIN. ALLOWED VALUE FOR 'IZ0IN'
C          (I*4)   IZ0MAX  = PARAMETER: MAX. ALLOWED VALUE FOR 'IZ0IN'
C
C          (I*4)   IZ0LST  = LAST VALUE OF 'IZ0IN' FOR  WHICH  INPUT
C                            DATA WAS READ.
C          (I*4)   IZLAST  = LAST VALUE OF 'IZIN' FOR  WHICH  INPUT
C                            DATA WAS READ.
C          (I*4)   IUNIT   = UNIT TO WHICH INPUT DATA SET IS ALLOCATED
C          (I*4)   NBSEL   = TOTAL NUMBER OF DATA-BLOCKS READ FROM INPUT
C                            DATA SET.
C          (I*4)   IZ0     = INPUT FILE - EMITTING ION - NUCLEAR CHARGE
C          (I*4)   IZ      = INPUT FILE - EMITTING ION - CHARGE
C          (I*4)   IZ1     = INPUT FILE - EMITTING ION - CHARGE + 1
C          (I*4)   IPOS    = USED IN CONVERTING CWAVEL -> WLNGTH
C
C          (L*4)   LOPEN   = .TRUE.  => INPUT DATA SET OPEN.
C                            .FALSE. => INPUT DATA SET CLOSED.
C
C          (C*2)   ESYM    = INPUT FILE - EMITTING ION - ELEMENT SYMBOL
C          (C*3)   EXTIN   = CURRENT ADAS SOURCE DATA FILE EXTENSION
C          (C*3)   EXTLST  = ADAS SOURCE DATA FILE EXT. USED LAST TIME
C                            DATA WAS READ.
CA         (C*80)  UIDIN   = CURRENT ADAS SOURCE DATA USER ID.
CA         (C*80)  UIDLST  = ADAS SOURCE DATA USER ID USED LAST TIME
C                            DATA WAS READ.
C          (C*8)   GRPIN   = CURRENT ADAS SOURCE DATA GROUPNAME
C          (C*8)   GRPLST  = ADAS SOURCE DATA GROUPNAME USED LAST TIME
C                            DATA WAS READ.
CA         (C*80)  TYPIN   = ADAS DATA FILE SUBDIRECTORY (OPTIONAL)
CA         (C*80)  TYPLST  = ADAS DATA FILE SUBDIRECTORY USED LAST TIME
C                            DATA WAS READ.
C          (C*11)  C11     = USED IN CONVERTING CWAVEL -> WLNGTH
CA         (C*80)  DSNREQ  = NAME OF DATA SET REQUESTED
C                            (MAY OR MAY NOT EXIST)
CA         (C*80)  DSNAME  = NAME OF DATA SET INTERROGATED
C
C          (I*4)   ISELA() = INPUT DATA FILE: DATA-BLOCK ENTRY INDICES.
C                            DIMENSION: DATA-BLOCK INDEX
C          (I*4)   ITA()   = INPUT DATA SET-NUMBER OF ELECTRON TEMPERA-
C                            TURES.
C                            DIMENSION: DATA-BLOCK INDEX
C          (I*4)   IDA()   = INPUT DATA SET-NUMBER OF ELECTRON DENSITIES
C                            DIMENSION: DATA-BLOCK INDEX
C
C          (R*8)   TETA(,) = INPUT DATA SET -
C                            ELECTRON TEMPERATURES (UNITS: eV)
C                            1st DIMENSION: ELECTRON TEMPERATURE INDEX
C                            2nd DIMENSION: DATA-BLOCK INDEX
C          (R*8)   TEDA(,) = INPUT DATA SET -
C                            ELECTRON DENSITIES    (UNITS: cm-3)
C                            1st DIMENSION: ELECTRON DENSITY     INDEX
C                            2nd DIMENSION: DATA-BLOCK INDEX
C          (R*8)   PEC(,,)   =INPUT DATA SET -
C                             FULL SET OF IONISATIONS PER PHOTON
C                             1st DIMENSION: ELECTRON TEMPERATURE INDEX
C                             2nd DIMENSION: ELECTRON DENSITY     INDEX
C                             3rd DIMENSION: DATA-BLOCK INDEX
C
C          (C*10)  CWAVEL() = INPUT FILE - WAVELENGTH (ANGSTROMS)
C                             DIMENSION: DATA-BLOCK INDEX
C          (C*8)   CFILE()  = INPUT FILE - SPECIFIC ION FILE SOURCE
C                             DIMENSION: DATA-BLOCK INDEX
C          (C*8)   CTYPE()  = INPUT FILE - TYPE OF DATA (IE EXCIT., ETC)
C                             DIMENSION: DATA-BLOCK INDEX
C          (C*2)   CINDM()  = INPUT FILE - METASTABLE INDEX
C                             DIMENSION: DATA-BLOCK INDEX
C
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C          E3FILE     ADAS      OPEN DATA SET FOR SELECTED EMITTER
C          E3DATA     ADAS      FETCH INPUT DATA FROM SELECTED DATA SET
C          E3CHKB     ADAS      CHECK VALIDITY OF ION AND 'IBSEL'
C          E3SPLN     ADAS      INTERPOLATE DATA WITH TWO WAY SPLINES
C          E3TITL     ADAS      CREATE DESCRIPTIVE TITLE FOR OUTPUT
C          XXUID      ADAS      FETCHES/SETS ADAS SOURCE DATA USER ID
C          XXSPEC     ADAS      FETCHES/SETS ADAS SOURCE DATA FILE NAME+
C
C AUTHOR:  H. P. SUMMERS
C          K1/1/57
C          JET EXT. 4941
C
C DATE:    11/10/91
C
C UPDATE:  05/12/91 - PE BRIDEN: 'NSTORE' INCREASED FROM 10 TO 100
C
C UPDATE:  28/02/92 - PE BRIDEN: 'NSTORE' INCREASED FROM 100 TO 150
C
C UPDATE:  10/03/93 - PE BRIDEN: INTRODUCED CALL TO XXUID TO ESTABLISH
C                                IF USERID OF INPUT DATASET CHANGES
C                                BETWEEN CALLS.
C                                SAVE NAME OF LAST READ DATASET.
C                                (ADDED VARIABLES UIDIN,UIDLST,DSNREQ)
C
C UPDATE:   2/09/93 - HPS      : INTRODUCED CALL TO XXSPEC TO ESTABLISH
C                                IF USRGRP, USRTYP AND USREXT OF INPUT
C                                DATASET CHANGES BETWEEN CALLS.
C                                SAVE NAME OF LAST READ DATASET.
C                                (ADDED VARIABLES GRPIN,GRPLST,TYPIN,
C                                 TYPLST, EXTIN, EXTLST)
C
C UPDATE:   6/05/94 - PEB      : INCREASED PARAMETER NSTORE 150 -> 350
C
C UPDATE:   3/11/94 - L.JALOTA : CHANGED DSNAME, UIDIN SIZE TO 80 CHARS.
C UPDATE:  23/11/94 - L.JALOTA : TIDIED UP STRING LENGTH DEFINITIONS
C
C VERSION: 1.1                  DATE: 25-05-95
C MODIFIED: TIM HAMMOND
C           - PUT UNDER S.C.C.S. CONTROL
C
C VERSION: 1.2                  DATE: 22-04-96
C MODIFIED: TIM HAMMOND/PAUL BRIDEN
C           - INCREASED PARAMETERS: NSTORE: 350 -> 500
C                                   NDDIM:   20 -> 26
C                                   NTDIM:   20 -> 35
C           - MODIFIED CONVERSION OF CWAVEL->WLNGTH
C             (ADDED VARIABLES C11 AND IPOS)
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      INTEGER     NSTORE         , NTDIM             , NDDIM
      INTEGER     IZ0MIN         , IZ0MAX
C-----------------------------------------------------------------------
      PARAMETER(  NSTORE = 500   , NTDIM  = 35       , NDDIM =  26  )
      PARAMETER(  IZ0MIN =   1   , IZ0MAX = 92       )
C-----------------------------------------------------------------------
      INTEGER     IBSEL          ,
     &            IZ0IN          , IZIN              ,
     &            ITVAL          , IRCODE
      INTEGER     IZ0LST         , IZLAST            ,
     &            IUNIT          , NBSEL             ,
     &            IZ0            , IZ                , IZ1          ,
     &            IPOS
C-----------------------------------------------------------------------
      real (kind=R8) ::
     *         WLNGTH
C-----------------------------------------------------------------------
      LOGICAL     LOPEN
C-----------------------------------------------------------------------
      CHARACTER   ESYM*2         , EXTIN*3           , EXTLST*3     ,
     &                             UIDIN*256         , UIDLST*256   ,
     &                             GRPIN*8           , GRPLST*8     ,
     &                             TYPIN*80          , TYPLST*80    ,
     &            C11*11         ,
     &            DSNREQ*256      , DSNAME*256         , TITLX*296
C-----------------------------------------------------------------------
      INTEGER     ISELA(NSTORE)             ,
     &            ITA(NSTORE)               , IDA(NSTORE)
C-----------------------------------------------------------------------
      real (kind=R8) ::
     *     TVAL(ITVAL)               , DVAL(ITVAL)          ,
     &            PECA(ITVAL)
C-----------------------------------------------------------------------
      LOGICAL     LTRNG(ITVAL)              , LDRNG(ITVAL)
C-----------------------------------------------------------------------
      CHARACTER   CINDM(NSTORE)*2           , CFILE(NSTORE)*8      ,
     &            CTYPE(NSTORE)*8           , CWAVEL(NSTORE)*10
C-----------------------------------------------------------------------
      real (kind=R8) ::
     *         TETA(NTDIM,NSTORE)        , TEDA(NDDIM,NSTORE)
      real (kind=R8) ::
     *         PEC(NTDIM,NDDIM,NSTORE)
C-----------------------------------------------------------------------
      SAVE        DSNAME      , UIDLST     , GRPLST     , TYPLST   ,
     &            EXTLST      , IZ0LST     , IZLAST     ,
     &            IUNIT
      SAVE        IZ0,IZ,IZ1,ESYM,NBSEL,ISELA,CWAVEL,CFILE,CTYPE,
     &            CINDM,ITA,IDA,TETA,TEDA,PEC
C-----------------------------------------------------------------------
      DATA        DSNAME /' '/   , UIDLST /' '/ , GRPLST /' '/     ,
     &            TYPLST /' '/   ,
     &            EXTLST /' '/   , IZ0LST /0/   , IZLAST /0/       ,
     &            IUNIT  /15/
C-----------------------------------------------------------------------
      external subini, subend
C-----------------------------------------------------------------------
C
      call subini ('spec')
      IRCODE = 0
      LOPEN  = .FALSE.
C
C-----------------------------------------------------------------------
C ESTABLISH CURRENT ADAS SOURCE DATA UID, FILE AND FILE EXTENSION
C-----------------------------------------------------------------------
C
      UIDIN  = '?'
      CALL XXUID(UIDIN)
C
      GRPIN  = '?'
      TYPIN  = '?'
      EXTIN  = '?'
      CALL XXSPEC(GRPIN , TYPIN , EXTIN)
C
C-----------------------------------------------------------------------
C
         IF ( (IZ0IN .GE. IZ0MIN)  .AND.  (IZ0IN .LE. IZ0MAX) ) THEN
C
C-----------------------------------------------------------------------
C IF NEW EMITTING ION ENTERED OR SOURCE DATA USERID HAS CHANGED:
C - OPEN THE REQUESTED DATA SET & READ IN COMPLETE SET OF RATE DATA.
C-----------------------------------------------------------------------
C
               IF ( (IZLAST .NE. IZIN)   .OR.
     &              (IZ0LST .NE. IZ0IN)  .OR.
     &              (UIDLST .NE. UIDIN)  .OR.
     &              (GRPLST .NE. GRPIN)  .OR.
     &              (TYPLST .NE. TYPIN)  .OR.
     &              (EXTLST .NE. EXTIN)      )  THEN
C
                  CALL E3FILE( IUNIT , IZ0IN , IZIN , IRCODE , DSNREQ )
                     IF (IRCODE .EQ. 0) THEN
                        LOPEN  = .TRUE.
                        IZLAST = IZIN
                        IZ0LST = IZ0IN
                        UIDLST = UIDIN
                        GRPLST = GRPIN
                        TYPLST = TYPIN
                        EXTLST = EXTIN
                        DSNAME = DSNREQ
                        CALL E3DATA( IUNIT  , DSNAME ,
     &                               NSTORE , NTDIM  , NDDIM  ,
     &                               IZ0    , IZ     , IZ1    , ESYM  ,
     &                               NBSEL  , ISELA  ,
     &                               CWAVEL , CFILE  , CTYPE  , CINDM ,
     &                               ITA    , IDA    ,
     &                               TETA   , TEDA   ,
     &                               PEC
     &                             )
                     ENDIF
C
               ENDIF
C
C-----------------------------------------------------------------------
C CHECK VALIDITY OF 'IBSEL' AND ENTERED EMITTING ION CHARGE VALUE.
C-----------------------------------------------------------------------
C
               IF ( IRCODE .EQ. 0 ) THEN
                  CALL E3CHKB( IUNIT , NBSEL  , IBSEL ,
     &                         IZ0IN , IZIN   ,
     &                         IZ0   , IZ     ,
     &                         LOPEN , IRCODE
     &                       )
               ENDIF
C
C-----------------------------------------------------------------------
C 1) INTERPOLATE WITH TWO WAY SPLINES.
C 2) CREATE TITLE FOR OUTPUT.
C-----------------------------------------------------------------------
C
               IF ( IRCODE .EQ. 0 ) THEN
                  CALL E3SPLN( NTDIM             , NDDIM         ,
     &                         ITA(IBSEL)        , IDA(IBSEL)    ,
     &                         ITVAL             ,
     &                         TETA(1,IBSEL)     , TEDA(1,IBSEL) ,
     &                         TVAL              , DVAL          ,
     &                         PEC(1,1,IBSEL)    , PECA          ,
     &                         LTRNG             , LDRNG
     &                       )
                  CALL E3TITL( IBSEL         , DSNAME       ,
     &                         ESYM          , IZ           ,
     &                         CWAVEL(IBSEL) , CINDM(IBSEL) ,
     &                         TITLX
     &                       )
C
C MODIFIED - PEB - 10/04/96
C
                  C11  = CWAVEL(IBSEL)//'A'
                  IPOS = INDEX(C11,'*')
C
                     IF (IPOS .EQ. 0) THEN
                       IPOS = INDEX(C11,'A') - 1
                       READ(C11(1:IPOS),*) WLNGTH
                     ELSE
                       WLNGTH = 0
                     ENDIF
C
C END OF MODIFICATION
C
               ENDIF
C
         ELSE
            IRCODE = 4
         ENDIF
C
C-----------------------------------------------------------------------
C
      call subend ()
      RETURN
      END SUBROUTINE SPEC
CX UNIX PORT - SCCS info: Module @(#)e3file.for 1.2 Date 5/25/95
CX
      SUBROUTINE E3FILE( IUNIT , IZ0 , IZ , IRCODE , DSNAME )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: E3FILE *********************
C
C  PURPOSE: TO OPEN A PHOTON EMISSIVITY 'IONELEC' DATA SET
C           BY DEFAULT, OR AN ALTERNATIVE DATA SET IF REQUIRED, FOR
C           EMITTING ION WITH NUCLEAR CHARGE 'IZ0' AND CHARGE 'IZ'.
C           THIS WILL BE CONNECTED TO UNIT 'IUNIT'.
C
C   DATA SET OPENED: $ADASUSER/<DEFADF>/<GROUP>(OPTIONAL)/<TYPE>/
C                        <GROUP_EXT>#<ELEMENT SYMBOL><CHARGE>
C
C  CALLING PROGRAM: SPEC
C
C  SUBROUTINE:
C
C  INPUT : (I*4)   IUNIT   = UNIT TO WHICH DATA SET WILL BE CONNECTED
C  INPUT : (I*4)   IZ0     = NUCLEAR CHARGE OF EMITTING ION REQUESTED
C  INPUT : (I*4)   IZ      = ION CHARGE OF EMITTING ION REQUESTED
C
C  OUTPUT: (I*4)   IRCODE  = RETURN CODE FROM SUBROUTINE:
C                            0 => DATA SET SUCCESSFULLY CONNECTED
C                            1 => REQUESTED DATA  SET  MEMBER  DOES  NOT
C                                 EXISTS - DATA SET NOT CONNECTED.
C                            6 => INVALID VALUE FOR 'IZ' ENTERED.
C                                 ( 0 <= 'IZ' <= 99 )
C                            9 => REQUESTED DATA  SET  EXISTS BUT CANNOT
C                                 BE OPENED.
C  OUTPUT: (C*256)  DSNAME  = NAME OF OPENED DATA SET UNDER UNIX
C
C          (I*4)   IDLEN   = LENGTH, IN BYTES, OF FIXED 'DSNAME' PREFIX
C          (I*4)   LENF1   = FIRST NON-BLANK CHR OF 'DSNAME' GROUP PART
C          (I*4)   LENF2   = LAST  NON-BLANK CHR OF 'DSNAME' GROUP PART
C          (I*4)   LENF3   = FIRST NON-BLANK CHR OF 'DSNAME' TYPE PART
C          (I*4)   LENF4   = LAST  NON-BLANK CHR OF 'DSNAME' TYPE PART
C          (I*4)   LENF5   = LAST  NON-BLANK CHR OF 'DSNAME' USERID PART
C          (I*4)   LENF6   = LAST  NON-BLANK CHR OF 'DSNAME' USERID PART
C          (I*4)   LENF7   = LAST  NON-BLANK CHR OF 'DSNAME' EXTENSION PART
C          (I*4)   LENF8   = LAST  NON-BLANK CHR OF 'DSNAME' EXTENSION PART
C          (I*4)   IZEND   = LAST BYTE WRITTEN TO IN 'CZ'. (= 1 OR 2)
C          (C*1)   HASH    = '#' IF NON-BLANK EXT, ELSE ' '.
C          (C*2)   CZ      = 'IZ' (NO LEADING BLANKS)
C          (C*2)   XFESYM  = FUNCTION - (SEE ROUTINES SECTION BELOW)
C          (C*2)   ESYM    = ELEMENT SYMBOL FOR NUCLEAR CHARGE 'IZ0'
CA         (C*3)   USREXT  = ADAS SOURCE DATA FILE EXTENSION
CA         (C*80)  USERID  = ADAS SOURCE DATA USER ID
CA         (C*8)   USRGRP  = ADAS SOURCE DATA GROUPNAME
CA         (C*80)  USRTYP  = ADAS SOURCE DATA TYPENAME
C          (C*6)   DEFADF  = DEFAULT DATA DIRECTORY, I.E. ADF13
C
C          (L*4)   LEXIST  = .TRUE.  => REQUESTED  DATA  SET  EXISTS.
C                            .FALSE. => REQUESTED  DATA  SET  DOES  NOT
C                                       EXIST.
C
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C          XXUID      ADAS      FETCHES/SETS ADAS SOURCE DATA USER ID
C          XXSPEC     ADAS      FETCHES/SETS ADAS SOURCE DATA FILENAME
C                               AND FILE EXTENSION
C          XFESYM     ADAS      CHARACTER*2 FUNCTION -
C                               GATHERS ELEMENT SYMBOL FOR NUC. CHARGE
C          XXSLEN     ADAS      FINDS FIRST AND LAST NON-BLANK
C                               CHARACTERS IN STRING.
C
C AUTHOR:  H.P. SUMMERS
C          K1/1/67
C          JET EXT. 4941
C DATE:    11/10/91
C
C UPDATE:  10/03/93 - PE BRIDEN: ADDED CALL TO XXUID AND USERID VARIABLE
C                                - NOW ALLOWS ANY INPUT DATASET USER ID.
C UPDATE:   2/09/93 - HPS      : ADDED CALL TO XXSPEC AND USRGRP, USRTYP
C                                  AND USREXT NAMES
C                                - NOW ALLOWS ANY INPUT DATASET FILENAME
C                                  AND EXTENSION
C UPDATE:  23/11/93 - PEB      : CORRECT ERROR - A '.' HAD MISTAKENLY
C                                BEEN PLACED BEFORE THE MEMBER NAME IN
C                                VARIABLE DSNAME.
C
C UPDATE:   1/11/94 - L. JALOTA: MODIFIED CODE FOR RUNNING UNDER UNIX
C                                USING NEW FILE NAMING CONVENTION.
C                                "ACTION" KEYWORD IN OPEN COMMAND IS IBM
C                                SO REMOVED HERE.
C UPDATE:   8/11/94 - L. JALOTA: ADDED DEFADF
C
C UPDATE:  22/11/94 - L. JALOTA: TIDIED UP CHARACTER STRING LENGTHS
C
C UPDATE:  22/03/95 _ HPS      : INTRODUCED HASH TO ELIMINATE # IN FILE IF
C                                THERE IS NO EXTENSION PART OF THE FILE NAME
C                                ALTER LOGIC TO ALLOW USRTYP, USREXT TO BE A
C                                SINGLE CHARACTER.
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      INTEGER     IUNIT         , IZ0            , IZ         ,
     &            IRCODE        , IDLEN          , IZEND      ,
     &            LENF1         , LENF2          ,
     &            LENF3         , LENF4          ,
     &            LENF5         , LENF6          ,
     &            LENF7         , LENF8
C-----------------------------------------------------------------------
      CHARACTER   XFESYM*2      , ESYM*2         , CZ*2       ,
     &            USERID*80     , DSNAME*256      , DEFADF*5   ,
     &            USRGRP*8      , USRTYP*80      , USREXT*3   ,
     &            HASH*1
C-----------------------------------------------------------------------
      LOGICAL     LEXIST
C-----------------------------------------------------------------------
C SET DEFAULT DIRECTORY.
      PARAMETER (DEFADF = 'adf15')
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C  IDENTIFY ADAS SOURCE DATA USER ID, GROUPNAME, TYPENAME
C  AND FILE EXTENSION
C  NOTE THAT UNDER UNIX TYPENAME IS AN OPTIONAL SUBDIRECTORY.
C-----------------------------------------------------------------------
C
      USERID = '?'
      CALL XXUID(USERID)
C
      USRGRP = '?'
      USRTYP = '?'
      USREXT = '?'
      CALL XXSPEC(USRGRP , USRTYP , USREXT)
C
C-----------------------------------------------------------------------
C FIND THE FIRST AND LAST NON-BLANK CHARACTER STRINGS IN USRGRP & USRTYP
C-----------------------------------------------------------------------
C
      CALL XXSLEN( USRGRP , LENF1 , LENF2 )
      CALL XXSLEN( USRTYP , LENF3 , LENF4 )
      CALL XXSLEN( USERID , LENF5 , LENF6 )
      CALL XXSLEN( USREXT , LENF7 , LENF8 )
C
C-----------------------------------------------------------------------
C
      DSNAME = USERID(LENF5:LENF6)//'/'//DEFADF//'/'//
     &            USRGRP(LENF1:LENF2)//'/'
      IDLEN  = INDEX(DSNAME,' ') - 1
C
      IF (LENF7 .EQ. 0) THEN
         HASH =' '
      ELSE
         HASH = '#'
      ENDIF
C
      IF (LENF3 .EQ. 0) THEN
         DSNAME = DSNAME(1:IDLEN)//USRGRP(LENF1:LENF2)//'_'
     &      //USREXT(LENF7:LENF8)//HASH
      ELSE
         DSNAME = DSNAME(1:IDLEN)//USRTYP(LENF3:LENF4)//'/'
     &      //USRGRP(LENF1:LENF2)//'_'//USREXT(LENF7:LENF8)//HASH
      ENDIF
C
      IDLEN  = INDEX(DSNAME,' ') - 1
      ESYM   = XFESYM( IZ0 )
C
         IF (IZ .LT. 0) THEN
            IRCODE = 6
         ELSEIF (IZ .LT. 10) THEN
            WRITE(CZ,1000) IZ
            IZEND  = 1
            IRCODE = 0
         ELSEIF (IZ .LT. 100) THEN
            WRITE(CZ,1001) IZ
            IZEND  = 2
            IRCODE = 0
         ELSE
            IRCODE = 6
         ENDIF
C
         IF (IRCODE .EQ. 0) THEN
C
            IF ( ESYM(2:2)  .EQ.  ' ' ) THEN
               DSNAME = DSNAME(1:IDLEN)//ESYM(1:1)//CZ(1:IZEND)//'.dat'
            ELSE
               DSNAME = DSNAME(1:IDLEN)//ESYM(1:2)//CZ(1:IZEND)//'.dat'
            ENDIF
            write(*,*) 'Checking on the existence of ',trim(dsname)
C
            call find_file(DSNAME,LEXIST)
C
            IF (.NOT.LEXIST) THEN
               IRCODE = 1
            ELSE
               OPEN(UNIT=IUNIT,FILE=DSNAME,STATUS='OLD',ERR=9999)
               IRCODE = 0
            ENDIF
C
         ENDIF
      RETURN
C
C-----------------------------------------------------------------------
C OPEN STATEMENT ERROR HANDLING
C-----------------------------------------------------------------------
C
 9999 IRCODE = 9
      RETURN
C
C-----------------------------------------------------------------------
C
 1000 FORMAT(I1,' ')
 1001 FORMAT(I2)
C
C-----------------------------------------------------------------------
C
      END SUBROUTINE E3FILE
CX UNIX PORT - SCCS Info : Module @(#)e3titl.for        1.9 Date 11/8/95
CX
      SUBROUTINE E3TITL( IBSEL  , DSNAME ,
     &                   ESYM   , IZ     ,
     &                   CWAVEL , CINDM  ,
     &                   TITLX
     &                 )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: E3TITL *********************
C
C  PURPOSE:  TO CREATE THE DESCRIPTIVE TITLE FOR SELECTED DATA-BLOCK.
C
C  CALLING PROGRAM: ADAS503/SPEC
C
C  SUBROUTINE:
C
C  INPUT : (I*4)  IBSEL    = SELECTED DATA-BLOCK: INDEX
C  INPUT : (C*(*)) DSNAME   = FULL INPUT DATA SET NAME
C
C  INPUT : (C*2)  ESYM     = INPUT DATA SET: EMITTING ION ELEMENT SYMBOL
C  INPUT : (I*4)  IZ       = INPUT DATA SET: EMITTING ION CHARGE-STATE
C
C  INPUT : (C*10) CWAVEL   = SELECTED DATA-BLOCK: WAVELENGTH (ANGS.)
C  INPUT : (C*2)  CINDM    = SELECTED DATA-BLOCK: METASTABLE INDEX
C
C  OUTPUT: (C*296)TITLX    = SELECTED DATA-BLOCK: DESCRIPTIVE TITLE
C
C          (C*2)  C2       = GENERAL USE 2 BYTE  CHARACTER  STRING
C          (C*3)  C3       = GENERAL USE 3 BYTE  CHARACTER  STRING
C          (C*9)  PATH     = UNIX ENVIRONMENT OF CURRENT PATH=$ADASCENT
C          (I*4)  POS_NOW  = CURRENT POSITION IN TITLE STRING
C          (I*4)  LEN_NAME = LENGTH OF FILENAME
C          (I*4)  IFIRST   = POSITION OF FIRST CHARACTER OF DSNAME
C          (I*4)  ILAST    = POSITION OF LAST CHARACTER OF DSNAME
C
C ROUTINES:
C
C AUTHOR:  PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C          K1/0/37
C          JET EXT. 2520
C
C DATE:    01/05/91
C
C UPDATE:  L. JALOTA  13/3/95   MODIFIED FOR USE UNDER UNIX
C
C VERSION 1.7   TIM HAMMOND     13/10/95
C               ADDED VARIABLE C3 TO ALLOW 3 DIGIT DATA BLOCKS
C               TO BE READ FROM THE VARIABLE IBSEL
C
C VERSION 1.8   TIM HAMMOND     08/11/95
C                - REMOVED SUPERFLUOUS VARIABLE TEMP_NAME
C
C VERSION 1.9   TIM HAMMOND     08/11/95
C                - REMOVED SUPERFLUOUS VARIABLE I4UNIT
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER    IBSEL           , IZ
C     INTEGER    XXFCHR
      INTEGER    POS_NOW         , LEN_NAME
      INTEGER    IFIRST          , ILAST
C-----------------------------------------------------------------------
      CHARACTER  DSNAME*(*)      , ESYM*2         ,
     &           CINDM*2         , CWAVEL*10      , TITLX*296  ,
     &           C2*2            , C3*3
C-----------------------------------------------------------------------
      CALL XXSLEN(DSNAME,IFIRST,ILAST)
      LEN_NAME = ILAST-IFIRST+1
      TITLX(1:LEN_NAME+1) = DSNAME(IFIRST:ILAST)
      POS_NOW = LEN_NAME+1
      WRITE(C3,1001) IBSEL
      TITLX(POS_NOW:POS_NOW+8) =  ' BLK=' // C3 // ';'
      POS_NOW = POS_NOW+9
      WRITE(C2,1000) IZ
      TITLX(POS_NOW:POS_NOW+6) = ESYM // '+' // C2
      POS_NOW = POS_NOW + 5
      TITLX(POS_NOW:POS_NOW+29)=' WAVELENGTH= '//CWAVEL//
     &                          ' MET='//CINDM
C
C-----------------------------------------------------------------------
 1000 FORMAT(I2)
 1001 FORMAT(I3)
C-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE E3TITL
CX UNIX PORT - SCCS info: Module @(#)e3chkb.for 1.1 Date 5/25/95
CX
      SUBROUTINE E3CHKB( IUNIT  , NBSEL  , IBSEL  ,
     &                   IZ0IN  , IZIN   ,
     &                   IZ0    , IZ     ,
     &                   LOPEN  , IRCODE
     &                 )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: E3CHKB *********************
C
C  PURPOSE: TO CHECK THE SELECTED BLOCK (IBSEL) OF DATA  EXISTS  IN  THE
C           INPUT DATA SET AND IF SO IT REPRESENTS THE ENTERED VALUES OF
C           'IZ0IN' (NUCLEAR CHARGE OF EMITTING ION)          &
C           'IZIN'  (CHARGE OF EMITTING ION)
C
C           IT ALSO CLOSES THE INPUT DATA SET ALLOCATION IF OPEN.
C
C  CALLING PROGRAM: SPEC
C
C  SUBROUTINE:
C
C  INPUT : (I*4)   IUNIT   = UNIT TO WHICH INPUT DATA SET IS ALLOCATED
C  INPUT : (I*4)   NBSEL   = TOTAL NUMBER OF DATA-BLOCKS READ FROM INPUT
C                            DATA SET.
C  INPUT : (I*4)   IBSEL   = INDEX OF DATA-BLOCK SELECTED FOR ANALYSIS
C
C  INPUT : (I*4)   IZ0IN   = REQUESTED: NUCLEAR CHARGE OF EMITTING ION
C  INPUT : (I*4)   IZIN    = REQUESTED: CHARGE OF EMITTING ION
C
C  INPUT : (I*4)   IZ0     = INPUT FILE: NUCLEAR CHARGE OF EMITTING ION
C  INPUT : (I*4)   IZ      = INPUT FILE: CHARGE OF EMITTING ION
C
C  I/O   : (L*4)   LOPEN   = INPUT : .TRUE.  => INPUT DATA SET OPEN.
C                                    .FALSE. => INPUT DATA SET CLOSED.
C                            OUTPUT: ALWAYS RETURNED AS .FALSE.
C  OUTPUT: (I*4)   IRCODE  = RETURN CODE FROM SUBROUTINE:
C                            0 => NO ERROR DETECTED.
C                            2 => DISCREPANCY BETWEEN REQUESTED CHARGES
C                                 AND THOSE IN INPUT DATA FILE.
C                            3 => SELECTED DATA-BLOCK  OUT OF RANGE  OR
C                                 DOES NOT EXIST.
C
C          (I*4)   I4UNIT  = FUNCTION (SEE ROUTINE SECTION BELOW)
C
C          (C*256)  DSNAME  = UNIX NAME OF DATA SET OPENED
C
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C          E3FILE     ADAS      OPEN DATA SET FOR SELECTED EMITTER
C          I4UNIT     ADAS      FETCH UNIT NUMBER FOR OUTPUT OF MESSAGES
C
C AUTHOR:  H. P. SUMMERS
C          K1/1/57
C          JET EXT. 4941
C
C DATE:    11/10/91
C
C UPDATE:  23/04/93 - PE BRIDEN - ADAS91: ADDED I4UNIT FUNCTION TO WRITE
C                                         STATEMENTS FOR SCREEN MESSAGES
C
C UPDATE:  24/05/93 - PE BRIDEN - ADAS91: CHANGED I4UNIT(0)-> I4UNIT(-1)
C
C UPDATE:   1/11/94 - L. JALOTA - UPDATED TO RUN UNDER UNIX ON DEC-ALPHA
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      INTEGER     I4UNIT
      INTEGER     IUNIT          , IRCODE            ,
     &            NBSEL          , IBSEL             ,
     &            IZ0IN          , IZIN              ,
     &            IZ0            , IZ
C-----------------------------------------------------------------------
      LOGICAL     LOPEN
C-----------------------------------------------------------------------
      CHARACTER   DSNAME*256
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
         IF ( (IBSEL .LE. 0)  .OR.  (IBSEL .GT. NBSEL) ) THEN
C
               IF (.NOT.LOPEN) THEN
                  CALL E3FILE( IUNIT , IZ0IN , IZIN , IRCODE , DSNAME )
                  IF (IRCODE .EQ. 0) LOPEN = .TRUE.
               ENDIF
C
               IF (LOPEN) THEN
                  IRCODE = 3
                  WRITE(I4UNIT(-1),1000) IRCODE , NBSEL , IBSEL
               ENDIF
C
         ELSE
C
               IF ( ( IZ0IN .NE. IZ0 )  .OR.
     &              ( IZIN  .NE. IZ  )      ) THEN
                  IRCODE = 2
               ELSE
                  IRCODE = 0
               ENDIF
C
         ENDIF
C
      IF (LOPEN) CLOSE( IUNIT )
C
C-----------------------------------------------------------------------
C
 1000 FORMAT( 1X,12('*'),' SPEC  RETURN CODE: ',I1,1X,12('*')/
     &        1X,'E3CHKB ERROR: SELECTED DATA BLOCK OUT OF RANGE'/
     &       15X,'MAXIMUM DATA-BLOCK INDEX  = ',I4/
     &       15X,'DATA-BLOCK INDEX SELECTED = ',I4/
     &        1X,46('*'))
C
C-----------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE E3CHKB
CX UNIX PORT - SCCS info: Module @(#)e3data.for 1.2 Date  1/23/96
CX
      SUBROUTINE E3DATA( IUNIT  , DSNAME ,
     &                   NSTORE , NTDIM  , NDDIM  ,
     &                   IZ0    , IZ     , IZ1    , ESYM  ,
     &                   NBSEL  , ISELA  ,
     &                   CWAVEL , CFILE  , CTYPE  , CINDM ,
     &                   ITA    , IDA    ,
     &                   TETA   , TEDA   ,
     &                   PEC
     &                 )
      use b2mod_types
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: E3DATA *********************
C
C  PURPOSE:  TO  FETCH  DATA  FROM  INPUT PHOTON EMISSIVITY FILE
C            FOR A GIVEN EMITTING ION (ELEMENT AND CHARGE).
C            (MEMBER STORED IN IONELEC.DATA - MEMBER PREFIX 'PEC#').
C
C  CALLING PROGRAM: ADAS503/SPEC
C
C  DATA:
C
C           UP TO 'NSTORE' SETS (DATA-BLOCKS) OF DATA MAY BE  READ FROM
C           THE FILE - EACH BLOCK FORMING A COMPLETE SET OF IONISATIONS
C           PER PHOTON VALUES FOR GIVEN TEMP./DENSITY COMBINATION. EACH
C           DATA-BLOCK  IS  ANALYSED INDEPENDENTLY OF ANY  OTHER  DATA-
C           BLOCK.
C
C           THE UNITS USED IN THE DATA FILE ARE TAKEN AS FOLLOWS:
C
C           TEMPERATURES        : EV
C           DENSITIES           : CM-3
C
C  SUBROUTINE:
C
C  INPUT : (I*4)  IUNIT    = UNIT TO WHICH INPUT FILE IS ALLOCATED.
CX INPUT : (C*44) DSNAME   = MVS DATA SET NAME OF DATA SET BEING READ
CA INPUT : (C*256) DSNAME   = NAME OF DATA FILE BEING READ
C
C  INPUT : (I*4)  NSTORE   = MAXIMUM NUMBER  OF  INPUT DATA-BLOCKS  THAT
C                            CAN BE STORED.
C  INPUT : (I*4)  NTDIM    = MAX NUMBER OF ELECTRON TEMPERATURES ALLOWED
C  INPUT : (I*4)  NDDIM    = MAX NUMBER OF ELECTRON DENSITIES    ALLOWED
C
C  OUTPUT: (I*4)  IZ0      = READ - EMITTING ION - NUCLEAR CHARGE
C  OUTPUT: (I*4)  IZ       = READ - EMITTING ION - CHARGE
C  OUTPUT: (I*4)  IZ1      = READ - EMITTING ION - CHARGE + 1
C  OUTPUT: (C*2)  ESYM     = READ - EMITTING ION - ELEMENT SYMBOL
C
C  OUTPUT: (I*4)  NBSEL    = NUMBER OF DATA-BLOCKS ACCEPTED & READ IN.
C  OUTPUT: (I*4)  ISELA()  = READ - DATA-SET DATA-BLOCK ENTRY INDICES
C                            DIMENSION: DATA-BLOCK INDEX
C
C  OUTPUT: (C*10) CWAVEL() = READ - WAVELENGTH (ANGSTROMS)
C                            DIMENSION: DATA-BLOCK INDEX
C  OUTPUT: (C*8)  CFILE()  = READ - SPECIFIC ION FILE SOURCE
C                            DIMENSION: DATA-BLOCK INDEX
C  OUTPUT: (C*8)  CTYPE()  = READ - DATA TYPE
C                            DIMENSION: DATA-BLOCK INDEX
C  OUTPUT: (C*2)  CINDM()  = READ - METASTABLE INDEX
C                            DIMENSION: DATA-BLOCK INDEX
C
C  OUTPUT: (I*4)  ITA()    = READ - NUMBER OF ELECTRON TEMPERATURES
C                            DIMENSION: DATA-BLOCK INDEX
C  OUTPUT: (I*4)  IDA()    = READ - NUMBER OF ELECTRON DENSITIES
C                            DIMENSION: DATA-BLOCK INDEX
C
C  OUTPUT: (R*8)  TETA(,)  = READ - ELECTRON TEMPERATURES (UNITS: eV)
C                            1st DIMENSION: ELECTRON TEMPERATURE INDEX
C                            2nd DIMENSION: DATA-BLOCK INDEX
C  OUTPUT: (R*8)  TEDA(,)  = READ - ELECTRON DENSITIES (UNITS: CM-3)
C                            1st DIMENSION: ELECTRON DENSITY INDEX
C                            2nd DIMENSION: DATA-BLOCK INDEX
C
C  OUTPUT: (R*8)  PEC(,,)   =READ - PHOTON EMISSIVITY VALUES
C                            1st DIMENSION: ELECTRON TEMPERATURE INDEX
C                            2nd DIMENSION: ELECTRON DENSITY INDEX
C                            3rd DIMENSION: DATA-BLOCK INDEX
C
C          (I*4)  I4EIZ0   = FUNCTION - (SEE ROUTINES SECTION BELOW)
C          (I*4)  I4FCTN   = FUNCTION - (SEE ROUTINES SECTION BELOW)
C          (I*4)  I4UNIT   = FUNCTION - (SEE ROUTINES SECTION BELOW)
C          (I*4)  IBLK     = ARRAY INDEX: DATA-BLOCK INDEX
C          (I*4)  ITT      = ARRAY INDEX: ELECTRON TEMPERATURE INDEX
C          (I*4)  ITD      = ARRAY INDEX: ELECTRON DENSITY     INDEX
C          (I*4)  NTNUM    = NUMBER OF ELECTRON TEMPERATURES FOR CURRENT
C                            DATA-BLOCK
C          (I*4)  NDNUM    = NUMBER OF ELECTRON DENSITIES    FOR CURRENT
C                            DATA-BLOCK
C          (I*4)  IABT     = RETURN CODE FROM 'I4FCTN'
C          (I*4)  IPOS1    = GENERAL USE STRING INDEX VARIABLE
C          (I*4)  IPOS2    = GENERAL USE STRING INDEX VARIABLE
C
C          (L*4)  LBEND    = IDENTIFIES WHETHER THE LAST OF THE  INPUT
C                            DATA SUB-BLOCKS HAS BEEN LOCATED.
C                            (.TRUE. => END OF SUB-BLOCKS REACHED)
C
C          (C*1)  CSLASH   = '/' - DELIMITER FOR 'XXHKEY'
C          (C*2)  C2       = GENERAL USE TWO BYTE CHARACTER STRING
C          (C*5)  IONNAM   = EMITTING ION READ FROM DATASET
C          (C*6)  CKEY1    = 'FILMEM' - INPUT BLOCK HEADER KEY
C          (C*4)  CKEY2    = 'TYPE  ' - INPUT BLOCK HEADER KEY
C          (C*4)  CKEY3    = 'INDM  ' - INPUT BLOCK HEADER KEY
C          (C*4)  CKEY4    = 'ISEL  ' - INPUT BLOCK HEADER KEY
C          (C*80) C80      = GENERAL USE 80 BYTE  CHARACTER  STRING  FOR
C                            THE INPUT OF DATA-SET RECORDS.
C
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C          XXHKEY     ADAS      OBTAIN KEY/RESPONSE STRINGS FROM TEXT
C          I4EIZ0     ADAS      INTEGER*4 FUNCTION    -
C                               RETURNS Z0 FOR GIVEN ELEMENT SYMBOL
C          I4FCTN     ADAS      INTEGER*4 FUNCTION    -
C                               CONVERT CHARACTER STRING TO INTEGER
C          I4UNIT     ADAS      INTEGER*4 FUNCTION    -
C                               FETCH UNIT NUMBER FOR OUTPUT OF MESSAGES
C
C AUTHOR:  H. P. SUMMERS
C          K1/1/57
C          JET EXT. 4941
C
C DATE:    11/10/91
C
C UPDATE:  05/12/91 - PE BRIDEN: IONNAM NOW ALLOWED TO OCCUPY EITHER
C                                4 OR 5 SPACES IN THE HEADER.
C
C UPDATE:  23/04/93 - PE BRIDEN - ADAS91: ADDED I4UNIT FUNCTION TO WRITE
C                                         STATEMENTS FOR SCREEN MESSAGES
C
C UPDATE:  24/05/93 - PE BRIDEN - ADAS91: CHANGED I4UNIT(0)-> I4UNIT(-1)
C
C UPDATE:  27/2/95  - L. JALOTA - IDL_ADAS : INCREASED SIZE DSNAME FOR
C                                            USE UNDER UNIX SYSTEMS
C
C UNIX-IDL PORT:
C
C VERSION: 1.2                          DATE: 23-1-96
C MODIFIED: TIM HAMMOND (TESSELLA SUPPORT SERVICES PLC)
C               - CORRECTED FORMAT STATEMENTS FOR DSNAME LENGTH
C
C-----------------------------------------------------------------------
      INTEGER    I4EIZ0                , I4FCTN              , I4UNIT
      INTEGER    IUNIT                 , NSTORE              ,
     &           NTDIM                 , NDDIM               ,
     &           IZ0                   , IZ                  ,
     &           IZ1                   , NBSEL
      INTEGER    IBLK                  ,
     &           ITT                   , ITD                 ,
     &           NTNUM                 , NDNUM               ,
     &           IPOS1                 , IPOS2               , IABT
C-----------------------------------------------------------------------
      LOGICAL    LBEND
C-----------------------------------------------------------------------
CX    CHARACTER  DSNAME*44             , ESYM*2
      CHARACTER  DSNAME*256            , ESYM*2
      CHARACTER  CSLASH*1              , C2*2                ,
     &           CKEY1*6               , CKEY2*4             ,
     &           CKEY3*4               , CKEY4*4             ,
     &           IONNAM*5              , C80*80
C-----------------------------------------------------------------------
      INTEGER    ISELA(NSTORE)         ,
     &           ITA(NSTORE)           , IDA(NSTORE)
C-----------------------------------------------------------------------
      CHARACTER  CINDM(NSTORE)*2       , CFILE(NSTORE)*8      ,
     &           CTYPE(NSTORE)*8       , CWAVEL(NSTORE)*10
C-----------------------------------------------------------------------
      real (kind=R8) ::
     *        TETA(NTDIM,NSTORE)    , TEDA(NDDIM,NSTORE)
      real (kind=R8) ::
     *        PEC(NTDIM,NDDIM,NSTORE)
C-----------------------------------------------------------------------
      SAVE       CSLASH                ,
     &           CKEY1                 , CKEY2                ,
     &           CKEY3                 , CKEY4
      EXTERNAL   SUBINI, SUBEND, XERRAB
C
C-----------------------------------------------------------------------
C
      DATA       CSLASH / '/' /
      DATA       CKEY1  / 'FILMEM' /   , CKEY2  / 'TYPE'   /  ,
     &           CKEY3  / 'INDM'   /   , CKEY4  / 'ISEL'   /
C
C **********************************************************************
C
      CALL SUBINI ('E3DATA')
      LBEND = .FALSE.
C
C-----------------------------------------------------------------------
C READ IN NUMBER OF DATA-BLOCKS PRESENT IN INPUT FILE.
C-----------------------------------------------------------------------
C
      READ(IUNIT,1000) C80
      READ(C80  ,*   ) NBSEL
C
      IPOS1 = INDEX(C80,'/') + 1
      IPOS2 = IPOS1          + 4
      READ(C80(IPOS1:IPOS2),1001) IONNAM
C
         IF ( IONNAM(2:2) .EQ. '+' ) THEN
            ESYM  = IONNAM(1:1) // ' '
         ELSE
            ESYM  = IONNAM(1:2)
         ENDIF
C
      IZ0   = I4EIZ0( ESYM )
      IZ    = I4FCTN( IONNAM(3:5) , IABT )
C
C PEB 05DEC91 - LINE BELOW ADDED
C
      IF (IABT .NE. 0) IZ = I4FCTN( IONNAM(3:4) , IABT )
C
      IZ1   = IZ + 1
C
         IF (NBSEL .GT. NSTORE) THEN
            WRITE(I4UNIT(-1),2000) DSNAME , NSTORE , NBSEL , NSTORE
            NBSEL = NSTORE
         ENDIF
C
C***********************************************************************
C INPUT DATA FOR EACH OF THE DATA-BLOCKS
C***********************************************************************
C
         DO 1 IBLK=1,NBSEL
C
C-----------------------------------------------------------------------
C INPUT TITLE, WAVELENGTH AND OTHER INFORMATION (CHECK BLOCK EXISTS)
C-----------------------------------------------------------------------
C
            IF (.NOT.LBEND) THEN
               READ(IUNIT,1000)   C80
C
                  IF ( C80(1:1) .NE. 'C') THEN
C
                     IPOS1 =  INDEX(C80,'A') + 1
                     IPOS2 =  INDEX(C80,'/') - 1
C
                        IF (IPOS1 .EQ. 1) THEN
                           IPOS1        = 11
                           CWAVEL(IBLK) = C80(1:10)
                        ELSEIF (IPOS1 .GT. 10) THEN
                           CWAVEL(IBLK) = C80(1:10)
                        ELSE
                           CWAVEL(IBLK) = C80(1:IPOS1-1)
                        ENDIF
C
                     READ(C80(IPOS1:IPOS2),*) IDA(IBLK) , ITA(IBLK)
C
                     CALL XXHKEY( C80 , CKEY1 , CSLASH , CFILE(IBLK)  )
                     CALL XXHKEY( C80 , CKEY2 , CSLASH , CTYPE(IBLK)  )
                     CALL XXHKEY( C80 , CKEY3 , CSLASH , CINDM(IBLK)  )
                     CALL XXHKEY( C80 , CKEY4 , CSLASH , C2           )
                     ISELA(IBLK) = I4FCTN( C2 , IABT )
C
                     NDNUM = IDA(IBLK)
                     NTNUM = ITA(IBLK)
C
                         IF (NTNUM .GT. NTDIM) THEN
                            WRITE(I4UNIT(-1),2001) DSNAME , IBLK  ,
     &                                            'TEMPERATURES'  ,
     &                                            NTDIM   , NTNUM
                            CALL XERRAB ('NTNUM .GT. NTDIM')
                         ENDIF
C
                         IF (NDNUM .GT. NDDIM) THEN
                            WRITE(I4UNIT(-1),2001) DSNAME , IBLK   ,
     &                                            'DENSITIES'      ,
     &                                            NDDIM   , NDNUM
                            CALL XERRAB ('NDNUM .GT. NDDIM')
                         ENDIF
C
C-----------------------------------------------------------------------
C INPUT TEMPERATURE, DENSITY AND PHOTON EMISSIVITIES FOR BLOCK
C-----------------------------------------------------------------------
C
                     READ(IUNIT,1002) ( TEDA(ITD,IBLK) , ITD=1,NDNUM )
                     READ(IUNIT,1002) ( TETA(ITT,IBLK) , ITT=1,NTNUM )
C
                        DO 3 ITD=1,NDNUM
                           READ(IUNIT,1002) (PEC(ITT,ITD,IBLK),
     &                                       ITT=1,NTNUM        )
    3                   CONTINUE
C
                  ELSE
                     WRITE(I4UNIT(-1),2002) DSNAME  , NBSEL    ,
     &                                     IBLK - 1 , IBLK - 1
                     LBEND = .TRUE.
                     NBSEL = IBLK - 1
                  ENDIF
C
            ENDIF
C
    1    CONTINUE
C
C***********************************************************************
C
 1000 FORMAT(1A80)
 1001 FORMAT(A5)
 1002 FORMAT(8D9.2)
C
 2000 FORMAT(/1X,31('*'),' E3DATA MESSAGE ',31('*')/
     &        2X,'INPUT PHOT. EMISS. DATA SET NAME: ',A80/
     &        2X,'NUMBER OF WAVELENGTH DATA-BLOCKS',
     &           ' IN INPUT DATA SET TOO GREAT.'/
     &        2X,'THE MAXIMUM NUMBER ALLOWED     = ',I3/
     &        2X,'THE NUMBER PRESENT IN DATA SET = ',I3/
     &        2X,'THEREFORE ONLY THE FIRST ',I3,' HAVE BEEN ACCEPTED'/
     &        1X,31('*'),' END OF MESSAGE ',31('*'))
 2001 FORMAT(/1X,32('*'),' E3DATA ERROR ',32('*')/
     &        2X,'INPUT PHOT. EMISS. DATA SET NAME: ',A80/
     &        2X,'DATA-BLOCK INDEX: ',I3//
     &        2X,'THE NUMBER OF INPUT ELECTRON ',A,' TOO LARGE'/
     &        2X,'THE MAXIMUM NUMBER ALLOWED       = ',I3/
     &        2X,'THE NUMBER PRESENT IN DATA-BLOCK = ',I3/
     &        1X,29('*'),' PROGRAM TERMINATED ',29('*'))
 2002 FORMAT(/1X,31('*'),' E3DATA MESSAGE ',31('*')/
     &        2X,'INPUT IONS./PHOTON DATA SET NAME: ',A80/
     &        2X,'INCONSISTENCY IN THE NUMBER OF WAVELENGTHS',
     &           ' EXPECTED AND READ.'/
     &        2X,'THE NUMBER EXPECTED = ',I3/
     &        2X,'THE NUMBER READ IN  = ',I3/
     &        2X,'THEREFORE ONLY ',I3,' HAVE BEEN ACCEPTED'/
     &        1X,31('*'),' END OF MESSAGE ',31('*'))
C
C-----------------------------------------------------------------------
C
      CALL SUBEND ()
      RETURN
      END SUBROUTINE E3DATA
CX UNIX PORT - SCCS info: Module @(#)xxspec.for 1.1 Date 5/25/95
CX
      SUBROUTINE XXSPEC( USRGRP , USRTYP , USREXT )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXSPEC *********************
C
C  PURPOSE: ADAS ROUTINE - SETS UP THE DEFAULT USEGRP, USRTYP AND USREXT
C           WHICH IDENTIFY THE FILENAME AND EXTENSION TO BE READ IN
C           SUBROUTINE SPEC.  IT WORKS IN THE SAME MANNER AS XXUID WHICH
C           WHICH ALLOWS THE DEFAULT USER SPACE TO BE SET
C
C
C           USRGRP: VALUE ON INPUT  =>  USRGRP: VALUE ON OUTPUT
C
C                     ?                 CURRENT ADAS DATA GROUPNAME
C                     *                 DEFAULT ADAS DATA GROUPNAME
C                  <BLANK>              *** USRGRP VALUE NOT CHANGED ***
C                  <OTHER>              *** USRGRP VALUE NOT CHANGED ***
C
C
C           USRTYP: VALUE ON INPUT  =>  USRTYP: VALUE ON OUTPUT
C
C                     ?                 CURRENT ADAS DATA TYPENAME
C                     *                 DEFAULT ADAS DATA TYPENAME
C                  <BLANK>              *** USRTYP VALUE NOT CHANGED ***
C                  <OTHER>              *** USRTYP VALUE NOT CHANGED ***
C
C
C           USREXT: VALUE ON INPUT  =>  USREXT: VALUE ON OUTPUT
C
C                     ?                 CURRENT ADAS DATA MEMBER EXTENS.
C                     *                 DEFAULT ADAS DATA MEMBER EXTENS
C                  <BLANK>              *** USREXT VALUE NOT CHANGED ***
C                  <OTHER>              *** USREXT VALUE NOT CHANGED ***
C
C         ? => QUERIES CURRENT ADAS DATA USRGRP, USRTYP OR USREXT
C              SETTING.
C         * => SETS ADAS DATA USEGRP, USRTYP OR USREXT SETTING
C              TO DEFAULT VALUE.
C   <BLANK> => SETS ADAS DATA USRGRP, USRTYP OR USREXT SETTING
C              TO DEFAULT VALUE.
C   <OTHER> => SETS ADAS DATA USRGRP, USRTYP OR USREXT SETTING
C              TO INPUT   VALUE.
C
C  CALLING PROGRAM: SPEC AND MAIN PROGRAMS USING SPEC, ADAS503
C
C  SUBROUTINE:
C
C  I/O   : (C*8)  USRGRP   = USRFIL UNDER WHICH ADAS DATA IS STORED
C                            (IF BLANK DEFAULTS TO DEFGRP)
C
CA I/O   : (C*80) USRTYP   = SUBDIRECTORY (OPTIONAL) WHERE ADAS DATA
CA                           FILE IS LOCATED. (IF BLANK DEFAULTS TO
CA                           DEFTYP)
C
C  I/O   : (C*3)  USREXT   = USREXT UNDER WHICH ADAS DATA IS STORED
C                            (IF BLANK DEFAULTS TO DEFEXT)
C
C          (C*8)  DEFGRP   = PARAMETER = DEFAULT USER GROUP FOR ADAS
C                                        DATA SOURCE
C
C          (C*80) DEFTYP   = PARAMETER = DEFAULT SUBDIRECTORY OF ADAS
C                                        DATA SOURCE
C
C          (C*3)  DEFEXT   = PARAMETER = DEFAULT USER EXTENSION FOR ADAS
C                                        DATA SOURCE
C
C          (C*8)  ADASGR   = CURRENT ADAS DATA SOURCE GROUP
CA         (C*80) ADASTY   = CURRENT ADAS DATA SOURCE TYPE
C          (C*3)  ADASEX   = CURRENT ADAS DATA SOURCE EXTENSION
C
C
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C
C NOTE:
C          TO CHECK CURRENT ADAS SOURCE USRGRP, USRTYP AND USREXT
C          CALL XXSPEC WITH ?`S AS INPUTS.
C
C AUTHOR:  HUGH P. SUMMERS, JET
C          K1/1/57
C          JET EXT. 4941
C
C DATE:     2/09/93
C
C UPDATE:  L. JALOTA - 1/11/94  (TESSELLA SUPPORT SERVICES PLC)
C                      CHANGED VALUES OF DEFGRP,DEFTYP, DEFEXT SUITABLE
C                      FOR DEC ALPHA DIRECTORY STRUCTURE.
C
C UPDATE:  L.JALOTA - 23/11/94 : TIDIED UP STRING LENGTH DEFINITIONS.
C-----------------------------------------------------------------------
       CHARACTER  DEFGRP*8             , DEFTYP*80     , DEFEXT*3
C-----------------------------------------------------------------------
      PARAMETER  (DEFGRP = 'ionelec' , DEFTYP = ' ' ,
     &            DEFEXT = 'pec')
C-----------------------------------------------------------------------
      CHARACTER  USRGRP*8             , ADASGR*8      ,
     &           USRTYP*80            , ADASTY*80     ,
     &           USREXT*3             , ADASEX*3
C-----------------------------------------------------------------------
      SAVE       ADASGR , ADASTY , ADASEX
C-----------------------------------------------------------------------
      DATA       ADASGR , ADASTY , ADASEX / DEFGRP , DEFTYP , DEFEXT /
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C SET ADAS DATA SOURCE GROUPNAME
C-----------------------------------------------------------------------
C
      IF (USRGRP(1:1) .EQ. '?') THEN
        USRGRP = ADASGR
      ELSE IF (USRGRP(1:1) .EQ. '*') THEN
        USRGRP = DEFGRP
        ADASGR = DEFGRP
      ELSE IF (USRGRP(1:1) .EQ. ' ') THEN
        ADASGR = DEFGRP
      ELSE
        ADASGR = USRGRP
      ENDIF
C
C-----------------------------------------------------------------------
C SET ADAS DATA SOURCE TYPENAME
C-----------------------------------------------------------------------
C
      IF (USRTYP(1:1) .EQ. '?') THEN
        USRTYP = ADASTY
      ELSE IF (USRTYP(1:1) .EQ. '*') THEN
        USRTYP = DEFTYP
        ADASTY = DEFTYP
      ELSE IF (USRTYP(1:1) .EQ. ' ') THEN
        ADASTY = DEFTYP
      ELSE
        ADASTY = USRTYP
      ENDIF
C
C-----------------------------------------------------------------------
C SET ADAS DATA SOURCE EXTENSION
C-----------------------------------------------------------------------
C
      IF (USREXT(1:1) .EQ. '?') THEN
        USREXT = ADASEX
      ELSE IF (USREXT(1:1) .EQ. '*') THEN
        USREXT = DEFEXT
        ADASEX = DEFEXT
      ELSE IF (USREXT(1:1) .EQ. ' ') THEN
        ADASEX = DEFEXT
      ELSE
        ADASEX = USREXT
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE XXSPEC
CX UNIX PORT - SCCS info: Module @(#)xxuid.for  1.2 Date 8/9/96
CX
      SUBROUTINE XXUID( USERID )
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXUID  *********************
C
C  PURPOSE: ADAS ROUTINE - SETS UP THE DEFAULT USERID WHICH STORES THE
C           DATA TO BE READ USING STANDARD ADAS DATA READING ROUTINES.
CA          UNDER UNIX PORT THIS NOW GETS THE ENVIRONMENT VARIABLE
CA          "ADASCENT".
C
C           USERID: VALUE ON INPUT  =>  USERID: VALUE ON OUTPUT
C
C                     ?                 CURRENT ADAS DATA SOURCE USER ID
C                     *                 DEFAULT ADAS DATA SOURCE USER ID
C                  <BLANK>              *** USERID VALUE NOT CHANGED ***
C                  <OTHER>              *** USERID VALUE NOT CHANGED ***
C
C           ? => QUERIES CURRENT ADAS SOURCE USERID SETTING.
C           * => SETS ADAS SOURCE USERID SETTING TO DEFAULT VALUE.
C     <BLANK> => SETS ADAS SOURCE USERID SETTING TO DEFAULT VALUE.
C     <OTHER> => SETS ADAS SOURCE USERID SETTING TO INPUT   VALUE.
C
C  CALLING PROGRAM: GENERAL USE
C
C  SUBROUTINE:
C
CX  I/O   : (C*6)  USERID   = USERID UNDER WHICH ADAS DATA IS STORED
CX                           (IF BLANK DEFAULTS TO DEFUID)
CX
CX         (C*6)  DEFUID   = PARAMETER = DEFAULT USER ID FOR ADAS DATA
CX                                      SOURCE
CX
CX         (C*6)  ADASID   = CURRENT ADAS DATA SOURCE USER ID
C
C  I/O   : (C*80)  USERID   = USERID UNDER WHICH ADAS DATA IS STORED
C                           (IF BLANK DEFAULTS TO DEFUID)
C
C          (C*80)  DEFUID   = PARAMETER = DEFAULT USER ID FOR ADAS DATA
C                                      SOURCE
C
C          (C*80)  ADASID   = CURRENT ADAS DATA SOURCE USER ID
C
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C
C NOTE:
C          TO CHECK CURRENT ADAS SOURCE USERID CALL XXUID WITH
C          ? AS INPUT.
C
C AUTHOR:  PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C          K1/0/37
C          JET EXT. 5023
C
C DATE:    10/03/93
C
C UPDATE:  L. JALOTA  - 7/11/94 : ADDED CALL TO GETENV TO FETCH UNIX
C                                 ENVIRONMENT VARIABLE ADASUSER
C
C VERSION: 1.2                          DATE: 08-08-96
C MODIFIED: WILLIAM OSBORN (TESSELLA SUPPORT SERVICES PLC)
C               - REMOVED SUPERFLUOUS 'EXTERNAL GETENV' STATEMENT
C
C-----------------------------------------------------------------------
      CHARACTER  DEFUID*80
C-----------------------------------------------------------------------
      CHARACTER  USERID*80          , ADASID*80              , ENV*80
C-----------------------------------------------------------------------
      SAVE       ADASID
C-----------------------------------------------------------------------
#ifdef USE_PXFGETENV
      integer lenval, ierror
#else
#ifdef NAGFOR
      integer lenval, ierror
#endif
#endif
*!!   (We avoid to declare, as external or intrinsic, the subroutine
*     getenv; it is external on some systems (SunOS, Unicos, ...) and
*     intrinsic on others (AIX, HP-UX, ...).  If getenv is declared
*     external on AIX or HP-UX it returns some non-printing garbage,
*     without error message (AIX/XLF Version 2.2, HP-UX Release 9.0))
CA GET THE VALUE OF ADASCENT WHICH SETS UP TOPLEVEL DATA DIRECTORY.
#ifdef NO_GETENV
      env=' '
#else
#ifdef NAGFOR
        call get_environment_variable('ADASCENT',
     .   status=ierror, length=lenval)
        if (ierror.eq.0) then
          call get_environment_variable('ADASCENT',value=env)
        elseif (ierror.eq.1) then
          env=' '
        endif
#else
#ifdef USE_PXFGETENV
      CALL PXFGETENV ('ADASCENT', 0, ENV, lenval, ierror)
#else
      call getenv ('ADASCENT', ENV)
#endif
#endif
#endif
      DEFUID = ENV
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C SET ADAS DATA SOURCE USER ID
C-----------------------------------------------------------------------
C
      IF (USERID(1:1) .EQ. '?') THEN
        USERID = ADASID
      ELSE IF (USERID(1:1) .EQ. '*') THEN
        USERID = DEFUID
        ADASID = DEFUID
      ELSE IF (USERID(1:1) .EQ. ' ') THEN
        ADASID = DEFUID
      ELSE
        ADASID = USERID
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE XXUID
CX UNIX PORT - SCCS info: Module @(#)xxsple.for 1.1 Date 5/1/95
CX
      SUBROUTINE XXSPLE( LSETX , IOPT   , FINTX ,
     &                   NIN   , XIN    , YIN   ,
     &                   NOUT  , XOUT   , YOUT  ,
     &                   DY    , LINTRP
     &                 )
      use b2mod_types
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXSPLE *********************
C
C  PURPOSE:           TO INTERPOLATE/EXTRAPOLATE USING CUBIC SPLINES
C                     (IF IOPT < 0 NO EXTRAPOLATION TAKES PLACE = VALUES
C                      SET TO ZERO).- LOGICAL ARRAY 'LINTRP()' SPECIFIES
C                      WHETHER OUTPUT SPLINE IS INTERPOLATED '.TRUE.' OR
C                      EXTRAPOLATED '.FALSE.'.
C
C                      (AS FOR 'XXSPLN' EXCEPT 'LINTRP' ARGUMENT ADDED).
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
C  OUTPUT: (R*8)  DY()    = DERIVATIVES AT INPUT KNOTS (ARRAY SIZE: NIN)
C  OUTPUT: (L*4)  LINTRP()= .TRUE.  => 'YOUT()' VALUE INTERPOLATED.
C                           .FALSE. => 'YOUT()' VALUE EXTRAPOLATED.
C                           (ARRAY SIZE: NOUT)
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
C          IF ( IOPT .LT. 0 ) - NO EXTRAPOLATION TAKES PLACE VALUES SET
C                             TO ZERO (CARE IF LOG OF OUTPUT IS NEEDED).
C          IF ( IOPT .GT. 4 ) PROGRAM STOPS
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
C          FINTX      ------    EXTERNAL  REAL*8  FUNCTION,   USED  TO
C                               TRANSFORM X-COORDINATES.
C
C
C AUTHOR:   PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C           K1/0/81
C           JET EXT. 4569
C
C DATE:     14/01/91 - ADAS91: AS FOR 'XXSPLN' BUT WITH 'LINTRP()' ADDED
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
      real (kind=R8) ::
     *        FINTX
      real (kind=R8) ::
     *        XK           , XKK
      real (kind=R8) ::
     *        T1           , T2         , T3       , T4   ,
     &           DL1          , DL2        , DL3
C-----------------------------------------------------------------------
      LOGICAL    LSETX        , LEXTRP
C-----------------------------------------------------------------------
      real (kind=R8) ::
     *        XIN(NIN)       , YIN(NIN)       ,
     &           XOUT(NOUT)     , YOUT(NOUT)     ,
     &           DY(NIN)
      real (kind=R8) ::
     *        QVAL(0:NIOPT)  , D2VAL(0:NIOPT) ,
     &           D3VAL(0:NIOPT) , UVAL(0:NIOPT)
      real (kind=R8) ::
     *        AGRL(4)
      real (kind=R8) ::
     *        X(NKNOTS)      , DELY(NKNOTS)   ,
     &           H(NKNOTS)      , Q(NKNOTS)      , U(NKNOTS)  ,
     &           D1(NKNOTS)     , D2(NKNOTS)     , D3(NKNOTS)
C-----------------------------------------------------------------------
      LOGICAL    LINTRP(NOUT)   , LUVAL(0:NIOPT)
C-----------------------------------------------------------------------
      DATA  QVAL / -0.5_R8, -0.5_R8,  0.0_R8,  0.0_R8,  0.0_R8/ ,
     &      D2VAL/  1.5_R8,  1.5_R8,  0.0_R8,  0.0_R8,  0.0_R8/ ,
     &      D3VAL/  0.0_R8,  0.0_R8,  0.0_R8, -0.5_R8,  0.0_R8/ ,
     &      UVAL /  0.0_R8, -1.5_R8,  1.0_R8, -1.5_R8,  0.0_R8/
      DATA  LUVAL/ .TRUE. , .FALSE. , .FALSE. , .FALSE. , .TRUE. /
      DATA  NOPT / -1     /
C-----------------------------------------------------------------------
      SAVE QVAL , D2VAL , D3VAL , UVAL , LUVAL
      SAVE NOPT
      SAVE X    , H     , Q     , D1   , D2   , D3
C-----------------------------------------------------------------------
      EXTERNAL SUBINI, SUBEND, XERTST
C***********************************************************************
      CALL SUBINI ('xxsple')
      CALL XERTST (NKNOTS .GE. NIN,
     &                 ' XXSPLE ERROR: TOO MANY KNOTS REQUESTED')
      CALL XERTST (IOPT .LE. NIOPT,
     &                 ' XXSPLE ERROR: INVALID IOPT VALUE ENTERED')
C-----------------------------------------------------------------------
      IF (IOPT .LT. 0) THEN
         LEXTRP = .FALSE.
         IF (NOPT .NE. 0)    LSETX=.TRUE.
      ELSE
         LEXTRP = .TRUE.
         IF (IOPT .NE. NOPT) LSETX=.TRUE.
      ENDIF
C-----------------------------------------------------------------------
C
C***********************************************************************
C SET UP PARAMETERS RELATING TO 'XIN'-AXIS
C***********************************************************************
C
         NOPT = MAX0( 0 , IOPT )
         NIN0 = NIN-1
C-----------------------------------------------------------------------
         IF (LSETX) THEN
               DO 1 I=1,NIN
                  X(I)=FINTX(XIN(I))
    1          CONTINUE
C-----------------------------------------------------------------------
            H(2)  = X(2)-X(1)
            Q(1)  = QVAL(NOPT)
            D2(1) = D2VAL(NOPT)/H(2)
            D3(1) = D3VAL(NOPT)
C-----------------------------------------------------------------------
               DO 2 I=2,NIN0
                  H(I+1) = X(I+1) - X(I)
                  T1     = 1.0_R8 / ( H(I+1) + H(I) )
                  T2     = H(I+1) * T1
                  T3     = 1.0_R8 - T2
                  T4     = 1.0_R8 / ( ( T2 * Q(I-1) ) + 2.0_R8 )
                  Q(I)   = -T3 * T4
                  D1(I)  = ( 3.0_R8 * T4 * T2 ) / H(I)
                  D2(I)  = ( 3.0_R8 * T4 * T3 ) / H(I+1)
                  D3(I)  =  T2 * T4
    2          CONTINUE
C-----------------------------------------------------------------------
            T4      = 1.0_R8 / ( Q(NIN0) + 2.0_R8 )
            D1(NIN) = ( 3.0_R8 * T4 ) / H(NIN)
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
               IF     (XK .GE. XIN(NIN)) THEN
                  INTER     = NIN
                  LINTRP(K) = .FALSE.
                  AGRL(1)   = 0.0_R8
                  AGRL(3)   = 0.0_R8
                     IF (XK .EQ. XIN(NIN)) THEN
                        LINTRP(K) = .TRUE.
                        AGRL(2)   = 0.0_R8
                        AGRL(4)   = 1.0_R8
                     ELSEIF (LEXTRP) THEN
                        AGRL(2)   = XKK-X(NIN)
                        AGRL(4)   = 1.0_R8
                     ELSE
                        AGRL(2)   = 0.0_R8
                        AGRL(4)   = 0.0_R8
                     ENDIF
C-----------------------------------------------------------------------
C EXTRAPOLATE: LOW 'XOUT' VALUE  - IF EXTRAPOLATION SWITCHED ON
C-----------------------------------------------------------------------
               ELSEIF (XK .LE. XIN(1))   THEN
                  INTER     = 2
                  LINTRP(K) = .FALSE.
                  AGRL(2)   = 0.0_R8
                  AGRL(4)   = 0.0_R8
                     IF (XK .EQ. XIN(1)) THEN
                        LINTRP(K) = .TRUE.
                        AGRL(1)   = 0.0_R8
                        AGRL(3)   = 1.0_R8
                     ELSEIF (LEXTRP) THEN
                        AGRL(1)   = XKK-X(1)
                        AGRL(3)   = 1.0_R8
                     ELSE
                        AGRL(1)   = 0.0_R8
                        AGRL(3)   = 0.0_R8
                     ENDIF
C-----------------------------------------------------------------------
C INTERPOLATE:
C-----------------------------------------------------------------------
               ELSE
                     DO 6 I=NIN,1,-1
                        IF (XK .LT. XIN(I)) INTER=I
    6                CONTINUE
                  LINTRP(K) = .TRUE.
                  DL1       =  ( X(INTER) - XKK    ) / H(INTER)
                  DL2       =  1.0_R8 - DL1
                  DL3       =  H(INTER) * DL1 * DL2
                  AGRL(1)   =  DL1 * DL3
                  AGRL(2)   = -DL2 * DL3
                  AGRL(3)   =  DL1 * DL1 * ( 1.0_R8 + DL2 + DL2 )
                  AGRL(4)   =  1.0_R8 - AGRL(3)
               ENDIF
C-----------------------------------------------------------------------
C EVALUATE 'YOUT'-VALUE FOR REQUESTED 'XOUT'-VALUE
C-----------------------------------------------------------------------
            YOUT(K) = ( AGRL(1)* DY(INTER-1) ) + ( AGRL(2)* DY(INTER) )
     &              + ( AGRL(3)*YIN(INTER-1) ) + ( AGRL(4)*YIN(INTER) )
C-----------------------------------------------------------------------
    5    CONTINUE
C***********************************************************************
      CALL SUBEND ()
      RETURN
      END SUBROUTINE XXSPLE
CX ULTRIX PORT - SCCS info: Module @(#)xxfchr.for       1.2 Date 4/24/95
CX
      SUBROUTINE XXFCHR( CSTRNG , SSTRNG, IFIRST , ILAST )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXFCHR *********************
C
C  PURPOSE: TO IDENTIFY THE FIRST AND LAST OCCURRENCE OF SSTRNG IN
C           CSTRNG, THE VALUES OF WHICH ARE IFIRST , ILAST.
C           - IF    NO    OCCURRENCE OF SSTRNG THEN IFIRST=ILAST=0
C           - IF ONLY ONE OCCURRENCE OF SSTRNG THEN IFIRST=ILAST>0
C           - IF OVER ONE OCCURRENCE OF SSTRNG THEN IFIRST>ILAST>0
C
C  NOTE   : ANY TRAILING BLANKS IN THE SEARCH STRING (SSTRNG) ARE
C           IGNORED.
C
C  CALLING PROGRAM: GENERAL USE
C
C  SUBROUTINE:
C
C  INPUT : (C*(*)) CSTRNG   = INPUT STRING FOR INTERROGATION
C  INPUT : (C*(*)) SSTRNG   = INPUT SEARCH STRING
C
C  OUTPUT: (I*4)   IFIRST   = BYTE POSITION OF FIRST OCCURRENCE OF SSTRNG
C                             IN CSTRNG.
C  OUTPUT: (I*4)   ILAST    = BYTE POSITION OF LAST  OCCURRENCE OF SSTRNG
C                             IN CSTRNG.
C
C          (I*4)   ILENC    = LENGTH OF 'CSTRNG' STRING IN BYTES
C          (I*4)   ILENS    = POSITION OF LAST  NON-BLANK BYTE IN SSTRNG
C          (I*4)   I        = GENERAL USE - INCLUDING DUMP FOR UNWANTED
C                             OUTPUT FROM XXSLEN SUBROUTINE.
C
C NOTE  :  ANY TRAILING BLANKS IN THE SEARCH STRING (SSTRNG) ARE
C          IGNORED.
C
C
C AUTHOR:  LALIT JALOTA (TESSELLA SUPPORT SERVICES PLC)
C
C DATE  :  27/10/94
C
C UPDATE:  06/03/95 - REVISED BY PAUL BRIDEN (TESSELLA SUPPORT SERVICES)
C                     1) ANALYSE ALL OF CSTRNG (NOT JUST NON-BLANK
C                        PART).
C                     2) ONLY IGNORE TRAILING BLANKS FOR SSTRNG
C                        (KEEP LEADING BLANKS).
C                     3) MODIFY DO LOOP INDEX RANGE TO ENSURE THAT
C                        YOU DO NOT GO BEYOND THE END OF CSTRNG.
C                     4) VERIFY LENGTH OF SSTRNG IS NONZERO.
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      INTEGER     IFIRST     , ILAST
      INTEGER     ILENC      , ILENS      , I
C-----------------------------------------------------------------------
      CHARACTER   CSTRNG*(*) , SSTRNG*(*)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
      IFIRST = 0
      ILAST  = 0
C
      ILENC  = LEN(CSTRNG)
      CALL XXSLEN(SSTRNG,I,ILENS)
C
         IF ( (ILENC .GE. ILENS)  .AND.  (ILENS .GT. 0) ) THEN
C
           DO 1 I=1,1+ILENC-ILENS
C
              IF (CSTRNG(I:I-1+ILENS) .EQ. SSTRNG(1:ILENS)) THEN
                 IF (IFIRST .EQ. 0) IFIRST = I
                 ILAST = I
              ENDIF
C
    1      CONTINUE
C
         ENDIF
C
C-----------------------------------------------------------------------
C
       RETURN
       END SUBROUTINE XXFCHR
CX UNIX PORT - SCCS info: Module @(#)xfesym.for 1.1 Date 5/25/95
CX
      FUNCTION XFESYM ( IZ0 )
      use b2mod_elements
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C ************** FORTRAN77 CHARACTER*2 FUNCTION: XFESYM ****************
C
C PURPOSE: TO RETURN THE SYMBOL FOR THE ELEMENT WITH NUCLEAR CHARGE IZ0
C          (CHARACTER*2 FUNCTION VERSION OF 'XXESYM')
C
C CALLING PROGRAM: GENERAL USE
C
C FUNCTION:
C
C          (C*2)  XFESYM  = FUNCTION NAME -
C                           SYMBOL OF ELEMENT WITH NUCLEAR CHARGE 'IZ0'
C          (I*4)  IZ0     = ELEMENT NUCLEAR CHARGE
C
C NOTES:    IF NUCLEAR CHARGE IS OUT OF RANGE, I.E.NOT BETWEEN 1 & 92,
C           THEN THE CHARACTER STRING 'XFESYM' IS RETURNED BLANK.
C
C ROUTINES: NONE
C
C
C AUTHOR:   PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C           K1/0/81
C           JET EXT. 4569
C
C DATE:     12/02/91
C
C UPDATES:  25/10/94 L. JALOTA (TESSELLA SUPPORT SERVICES PLC)
C                    CHANGED CASE OF SYMBOL TO LOWER CASE FOR UNIX
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER      IZ0
C-----------------------------------------------------------------------
      CHARACTER*2  XFESYM
C-----------------------------------------------------------------------
         IF ( (IZ0 .GT. size(elements)) .OR. (IZ0 .LT. 0) ) THEN
            XFESYM = ' '
         ELSE
            XFESYM = elements(IZ0)
         ENDIF
         call lower_case(XFESYM)                  !srv 14.10.13
C-----------------------------------------------------------------------
      RETURN
      END FUNCTION XFESYM
CX UNIX PORT - SCCS Info : Module @(#)xxhkey.for        1.1 Date 5/2/95
CX
      SUBROUTINE XXHKEY( CTEXT  , CKEY , CBREAK , CANS )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXHKEY *********************
C
C  PURPOSE: TO EXTRACT FROM A LINE OF TEXT 'CTEXT' A RESPONSE TO A KEY
C           IN THER FORM OF '<CKEY> = <CANS>'.
C
C  CALLING PROGRAM: GENERAL USE
C
C  SUBROUTINE:
C
C  INPUT : (C*(*)) CTEXT   = INPUT TEXT LINE CONTAINING KEY & RESPONSES
C  INPUT : (C*(*)) CKEY    = KEY TEXT
C  INPUT : (C*1  ) CBREAK  = KEY/RESPONSE PAIR SEPERATOR SYMBOL
C
C  OUTPUT: (C*(*)) CANS    = RERSPONSE FOR GIVEN KEY: BLANK IF NOT FOUND
C
C          (I*4)   LENTXT  = LENGTH IN BYTES OF 'CTEXT' STRING
C          (I*4)   LENKEY  = LENGTH IN BYTES OF 'CKEY' STRING
C          (I*4)   LENANS  = LENGTH IN BYTES OF 'CANS' STRING
C          (I*4)   IKEY    = LENGTH IN BYTES OF 'CKEY' IGNORING TRAILING
C                            BLANKS
C          (I*4)   IPOS1   = USED IN IDENTIFYING RELEVANT BYTES IN CTEXT
C          (I*4)   IPOS2   = USED IN IDENTIFYING RELEVANT BYTES IN CTEXT
C          (I*4)   IPOS3   = USED IN IDENTIFYING RELEVANT BYTES IN CTEXT
C          (I*4)   I       = GENERAL USE INDEX
C
C ROUTINES: NONE
C
C NOTES:    THIS ROUTINE EXTRACTS FROM 'CTEXT' A RESPONSE TO A GIVEN KEY
C           IN THER FORM OF '<CKEY> = <CANS>'. E.G. 'FILE = DSN001'
C           WOULD REQUIRE AS INPUT CKEY='FILE' AND WOULD GIVE AS OUTPUT
C           CANS='DSN001'. ALL KEY/RESPONSE PAIRS MUST BE SEPARATED BY
C           THE CHARACTER GIVEN BY 'CBREAK' E.G. A SLASH, AND EACH KEY
C           MUST BE FOLLOWED BY AN EQUALS SIGN. THE NUMBER OF SPACES
C           BETWEEN THE KEY AND THE EQUAL SIGN AND BETWEEN THE RESPONSE
C           AND THE EQUAL SIGN IS NOT IMPORTANT.
C
C           THE BYTE PRECEEDING THE KEY MUST BE A BLANK OR 'CBREAK'
C           CHARACTER UNLESS IT STARTS AT BYTE ONE IN 'CTEXT'.
C
C           IF A KEY DOES NOT EXIST IN 'CTEXT' THEN 'CANS' IS RETURNED
C           BLANK.
C
C           THE KEY IS TAKEN AS 'CKEY' REMOVING ANY TRAILING BLANKS.
C           LEADING BLANKS ARE LEFT IN PLACE AND WILL USED WHEN THE
C           THE SEARCH FOR THE KEY IS MADE:
C
C           I.E.  'DATA   ' AND 'DATA' ARE THE SAME KEY BUT
C                 ' DATA '  AND 'DATA ' ARE DIFFERENT KEYS ALTHOUGH
C                 BOTH WILL GIVE THE SAME RESULTS IF A SPACE EXISTS
C                 BEFORE 'DATA' IN THE INPUT TEXT LINE.
C
C           AN EXAMPLE OF AN INPUT TEXT LINE IS:
C
C           8524.0 A    5 7 /FILMEM = FBBH91BE/   CODE=   V2B DLN1   /
C
C           THIS WOULD GIVE THE FOLLOWING:
C
C           CKEY='FILMEM'  =>  CANS='FBBH91BE'
C           CKEY=' FILMEM' =>  CANS=' '
C           CKEY='CODE'    =>  CANS='V2B DLN1'
C           CKEY=' CODE'   =>  CANS='V2B DLN1'
C           CKEY='OTHER'   =>  CANS=' '
C
C           (IF THE CHARACTER STRING IS SHORTER THAN THE RESPONSE THEN
C            THE RESPONSE IS TRUNCATED ACCORDINGLY.)
C
C           SPACES CAN EXIST IN THE KEY. I.E. CKEY='PLOT A'. BUT CARE
C           SHOULD BE TAKEN WHEN USING PREFIXES ON A COMMON KEY BASE,
C           I.E. 'A PLOT', 'B PLOT'. THIS IS BECAUSE IF A SUBSEQUENT
C           KEY TO BE FOUND IS 'PLOT' THEN  EITHER OF THESE SATISFY
C           THIS CRITERION AS WELL AS 'PLOT' ITSELF.
C
C           AN EXAMPLE OF AN INPUT TEXT LINE IS:
C
C           A FILE=TEST0/B FILE = TEST1/FILE=TEST2/FILE 1=TEST3/FILE 2=/
C
C           THIS WOULD GIVE THE FOLLOWING:
C
C           CKEY='A FILE'  =>  CANS='TEST0'
C           CKEY='B FILE'  =>  CANS='TEST1'
C           CKEY='FILE'    =>  CANS='TEST0' (WRONG RESPONSE PICKED UP)
C           CKEY='FILE 1'  =>  CANS='TEST3'
C           CKEY='FILE 2'  =>  CANS='TEST4'
C
C           IT IS ALSO POSSIBLE TO IMBED RESPONSES
C
C           AN EXAMPLE OF AN INPUT TEXT LINE IS:
C
C           FILE 1 =  Z1 = 23 / FILE = FILE 1 = 6 /
C
C           THIS WOULD GIVE THE FOLLOWING:
C
C           CKEY='FILE 1'  =>  CANS='Z1 = 23'
C           CKEY=' FILE 1' =>  CANS='6'
C           CKEY='Z1'      =>  CANS='23'
C           CKEY='FILE'    =>  CANS='FILE 1 = 6'
C
C AUTHOR:  PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C          K1/0/37
C          JET EXT. 2520
C
C DATE:    26/04/91
C
C-----------------------------------------------------------------------
      INTEGER    LENTXT      , LENKEY      , LENANS      ,
     &           IPOS1       , IPOS2       , IPOS3       ,
     &           IKEY        , I
C-----------------------------------------------------------------------
      CHARACTER  CTEXT*(*)   , CKEY*(*)    , CBREAK*1    , CANS*(*)
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
      LENTXT = LEN(CTEXT)
      LENKEY = LEN(CKEY)
      LENANS = LEN(CANS)
C
C-----------------------------------------------------------------------
C FIND THE LENGTH OF THE KEY IGNORING ANY TRAILING BLANKS
C-----------------------------------------------------------------------
C
      IKEY   = 0
C
         DO 1 I = 1,LENKEY
            IF (CKEY(I:I) .NE. ' ') IKEY = I
    1    CONTINUE
C
      IF (IKEY .EQ. 0) GOTO 999
C
C-----------------------------------------------------------------------
C ESTABLISH IF A VALID KEY CAN BE FOUND IN 'CTEXT'
C-----------------------------------------------------------------------
C
      IPOS1 = 1
  100 IPOS1 = INDEX( CTEXT(IPOS1:LENTXT) , CKEY(1:IKEY) ) + IPOS1 - 1
      IPOS2 = IPOS1 + IKEY
C
         IF (IPOS2 .GT. LENTXT) THEN
            GOTO 999
         ELSEIF (IPOS1 .EQ. 1) THEN
            IPOS1 = IPOS2
         ELSEIF (IPOS1 .NE. 0) THEN
            IPOS1 = IPOS1 - 1
               IF ( (CTEXT(IPOS1:IPOS1) .EQ. ' '   )
     &                        .OR.
     &              (CTEXT(IPOS1:IPOS1) .EQ. CBREAK) ) THEN
                  IPOS1 = IPOS2
               ELSE
                  GOTO 999
               ENDIF
         ELSE
            GOTO 999
         ENDIF
C
C-----------------------------------------------------------------------
C ESTABLISH IF EQUAL SIGN EXISTS AFTER THE KEY (SEPERATED BY BLANKS).
C-----------------------------------------------------------------------
C
      IPOS2 = INDEX( CTEXT(IPOS1:LENTXT) , '=' )
      IPOS3 = IPOS2 + IPOS1
C
         IF (IPOS3 .GT. LENTXT) THEN
            GOTO 999
         ELSEIF (IPOS2 .EQ. 1) THEN
            IPOS2 = IPOS3
         ELSEIF (IPOS2 .NE. 0) THEN
            IPOS2 = IPOS3 - 2
               IF (CTEXT(IPOS1:IPOS2) .EQ. ' ') THEN
                  IPOS2 = IPOS3
               ELSE
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C SEE IF VALID KEY EXISTS FURTHER DOWN THE 'CTEXT' STRING.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
                  GOTO 100
               ENDIF
         ELSE
            GOTO 999
         ENDIF
C
C-----------------------------------------------------------------------
C FIND SEPERATOR CHARACTER AND IDENTIFY RESPONSE (IF PRESENT).
C-----------------------------------------------------------------------
C
      IPOS1 = 0
      IPOS3 = INDEX( CTEXT(IPOS2:LENTXT) , CBREAK )
C
         IF     (IPOS3 .EQ. 0) THEN
            IPOS3 = LENTXT
         ELSEIF (IPOS3 .EQ. 1) THEN
            GOTO 999
         ELSE
            IPOS3 = IPOS3 + IPOS2 - 2
         ENDIF
C
         DO 3 I = IPOS2,IPOS3
           IF ( ( IPOS1 .EQ. 0 )  .AND.
     1      ( CTEXT(I:I) .NE. ' ' ) ) IPOS1 = I
    3    CONTINUE
C
         IF (IPOS1 .EQ. 0) THEN
            GOTO 999
         ELSEIF ( (IPOS3-IPOS1+1) .GT. LENANS ) THEN
            IPOS3 = IPOS1 + LENANS - 1
         ENDIF
C
C-----------------------------------------------------------------------
C VALID RESPONSE FOUND - SET UP 'CANS'
C-----------------------------------------------------------------------
C
      CANS = CTEXT(IPOS1:IPOS3)
      RETURN
C
C-----------------------------------------------------------------------
C INVALID OR NO RESPONSE/KEY FOUND - RETURN 'CANS' AS BLANK
C-----------------------------------------------------------------------
C
  999 CANS = ' '
      RETURN
C
C-----------------------------------------------------------------------
C
      END SUBROUTINE XXHKEY
CX UNIX PORT - SCCS Info : Module @(#)i4eiz0.for        1.1 Date 5/2/95
CX
      FUNCTION I4EIZ0 ( ESYM )
      use b2mod_elements
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C *************** FORTRAN77 INTEGER*4 FUNCTION: I4EIZ0 *****************
C
C PURPOSE: TO RETURN THE NUCLEAR CHARGE FOR THE ELEMENT SYMBOL ESYM
C          (INTEGER*4 FUNCTION VERSION OF 'XXEIZ0')
C
C CALLING PROGRAM: GENERAL USE
C
C FUNCTION:
C
C          (I*4)  I4EIZ0  = FUNCTION NAME -
C                           ELEMENT NUCLEAR CHARGE
C          (C*2)  ESYM    = SYMBOL OF ELEMENT WITH NUCLEAR CHARGE I4EIZ0
C
C          (I*4)  I       = GENERAL ARRAY USE
C
C          (C*2)  SYMBOL()= SYMBOLS OF FIRST 'NSYM' ELEMENTS.
C                           ARRAY DIMENSION => NUCLEAR CHARGE
C
C NOTES:    IF SYMBOL IS NOT RECOGNISED, I.E. NOT IN Z0 RANGE 1 & 'NSYM',
C           THEN THE INTEGER 'I4EIZ0' IS RETURNED AS ZERO.
C
C ROUTINES: NONE
C
C
C AUTHOR:   PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C           K1/0/81
C           JET EXT. 4569
C
C DATE:     13/02/91
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER    I4EIZ0
C-----------------------------------------------------------------------
      INTEGER    I
C-----------------------------------------------------------------------
      CHARACTER  ESYM*2, SYMBOL*2, SYMBLL*2
C-----------------------------------------------------------------------
      I4EIZ0 = 0
      SYMBOL = ESYM
      call lower_case(SYMBOL)                  !srv 14.10.13
      DO I=1, size(elements)
        SYMBLL = elements(i)
        call lower_case(SYMBLL)                  !srv 14.10.13
        IF ( SYMBOL .EQ. SYMBLL ) I4EIZ0 = I
      END DO
C-----------------------------------------------------------------------
      RETURN
      END  FUNCTION I4EIZ0

      subroutine CalculateEmissivity(is,iwl,idwl)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_elements
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none

      real (kind=R8) ::
     *   wl,dwl
      integer ix,iy,is
      integer ic,nc,nz,iwl,idwl
      integer, allocatable, save :: iae(:),iar(:),iac(:)
      real (kind=R8) ::
     *   lNe(1),lNi(1),lNi1(1),lNH(1),lTe(1),e(1),p(1)
      real (kind=R8) ::
     *   ea(-1:nx,-1:ny)
      integer err
      character extension*3,filename*256,MIE*2,MIR*2,MIC*2,stmp*256
c      logical open_file

c     initialising
      ic=partition_number(is)  ! was nint(za(is))
      nc=nint(zn(is))
      wl=iwl*1.0_R8
      dwl=idwl*1.0_R8
      nz=size(elements)
      allocate(iae(nz),iar(nz),iac(nz))
c      if(.not.open_file(99,'source.location','data.local/adas|modules/adas')) return
c      read(99,*) extension,mie,mir,mic
c      read(99,*) filename
c      close(99)
c     calculating the emissivity
      write(*,*) 'Calculate the emissivity for species ',is
      filename='filelist.adas.adf15'
      call InitADAS(ic,nc,
     &              extension,filename,mie,mir,mic,err)
      if (err .eq. 0) then
         call GetIBSEL(wl,dwl,mie,mir,mic,nz,iae,iar,iac)
         do ix=-1,nx
            do iy=-1,ny
               lNe(1)=NE(ix,iy)/1.0e6_R8
               lNi(1)=NA(ix,iy,is)/1.0e6_R8
               lNi1(1)=NA(ix,iy,is+1)/1.0e6_R8
               lNH(1)=NH(ix,iy)/1.0e6_R8
               lTe(1)=TE(ix,iy)/eV
               call GetEmissivity(nz,iae,iar,iac,
     &          lNe(1),lNi(1),lNi1(1),lNH(1),lTe(1),1,e(1),p(1))
               ea(ix,iy)=e(1)/(4.0_R8*PI)*1.0e6_R8
            enddo
         enddo
      endif
      call copy(ea,1.0_R8,1,textbl)
      if(is.lt.ns) then
        stmp='Emissivity in photons/(m^{3} sr s)'//textpl(is)
      else
        stmp='Emissivity in photons/(m^{3} sr s)'//'  '
      endif
      call customise_label(stmp,graphlabel(iwrk))
      deallocate(iae,iar,iac)
      return

      end subroutine CalculateEmissivity
