      
C =====================================================================
      SUBROUTINE PLOTCONTOUR(NX,NY,NR,NZ,IS,
     . CR,CZ,FIELD,ZLOD,ZHID,PHYSICAL,LINLOG)
C =====================================================================
      use b2mod_types
      use b2mod_b2plot_values
      use b2mod_b2plot_debug
      use b2mod_b2plot_xy
      use b2mod_indirect
      use b2mod_boundary_namelist
      implicit none
#include "graphics_vars.h"
      INTEGER NX,NY,NR,NZ,IS,LINLOG
cxpb  INTEGER NNCUT_L
      REAL (KIND=R8) :: CR(-1:NR,-1:*),CZ(-1:NZ,-1:*),
     *     FIELD(-1:NX,-1:NY),ZLOD,ZHID
      CHARACTER*20 FMTZ
      CHARACTER*256 STRING
      CHARACTER*20 DEBUGSTRING
      LOGICAL PHYSICAL,DEBUG,LBLBR
      
      INTEGER NXSTRT(99),NYSTRT(99),NXEND(99),NYEND(99),NPERD(99)
      INTEGER LRWK,LIWK,LMAP,NWRK,NOGRPS
      PARAMETER (LRWK=50000,LIWK=50000,LMAP=5000000,NWRK=50000,NOGRPS=6)
      REAL (KIND=R4) :: RWRK(LRWK), XWRK(NWRK), YWRK(NWRK)
      INTEGER IWRK(LIWK)
      INTEGER MAP(LMAP),IAREA(NOGRPS),IGRP(NOGRPS)
      INTEGER LFIN(100),choose_color
      CHARACTER*20 LBLS(100)
      
      integer irold,nerrf,nerro,nfmtz,noz,nozm,
     .        i,ic,iiso,ix,ixl,iy,iyt,k,kreg,l,
     .        ixstrt,ixend,iystrt,iyend
      real rlo, rhi, zlo, zhi, zinc, zval, clv, special
      logical iniso
      integer lnblnk, no_available_colors
      real ceil
      external lnblnk, no_available_colors, ceil, xertst, subini, subend
      intrinsic mod
      
      EXTERNAL FILL,CPDRPL

      logical first
      save first
      data first /.true./
      
      INTEGER LNGDEC
      PARAMETER(LNGDEC=200*200)
      REAL (KIND=R4) :: Z(LNGDEC)

      call subini ('plotcontour')
      if(first) then
        call alloc_b2mod_b2plot_xy(LNGDEC)
        first=.false.
      endif
      rwrk=0   ! bug in NCAR: The variable (RWRK(209)) has an undefined value.
      z=0
      
C ========== GLOBAL OPTIONS ==========================================
#ifdef DBG
      DEBUG=.true.
#else
      DEBUG=.false.
#endif
      CALL CPSETI('MAP - MAPPING FLAG',4)
      IF(LABELS) THEN
        CALL CPSETI('LLP - LINE LABEL POSITIONING',3)
        IF(LALONG) CALL CPSETI('LLO - LINE LABEL ORIENTATION',1)
      ELSE
        CALL CPSETI('LLP - LINE LABEL POSITIONING',0)
      ENDIF
      IF(INTPOL) THEN
        CALL CPSETR('T2D - TENSION OF THE 2D SPLINES',1.0)
        CALL CPSETI('PIC - POINT INTERPOLATION FLAG',-10)
        CALL IPOLATE(.FALSE.) ! for further plots
        WRITE (*,*) 'using spline interpolation'
      ELSE
        CALL CPSETR('T2D - TENSION OF THE 2D SPLINES',0.0)
        CALL CPSETI('PIC - POINT INTERPOLATION FLAG',-10)
      ENDIF
      CALL CPSETC('CFT - CONSTANT FIELD LABEL',' ')
      LBLBR=.true.
C ====================================================================
      
      ZLO=ZLOD
      ZHI=ZHID
      rlo=field(-1,-1)
      rhi=field(-1,-1)
      DO IY=-1,NY
        DO IX=-1,NX
          RLO=MIN(RLO,REAL(FIELD(IX,IY)))
          RHI=MAX(RHI,REAL(FIELD(IX,IY)))
        ENDDO
      ENDDO
      IF(zlo.eq.zhi) then
        zlo=rlo
        zhi=rhi
      ENDIF
      
      if(zlo.eq.zhi) then
        zlo=min(zlo*0.9,zlo*1.1)
        zhi=max(zhi*0.9,zhi*1.1)
        if(zlo.eq.zhi) then
          zlo=-1.0
          zhi=+1.0
        endif
      endif
      
      WRITE(*,*) 'zlo =',ZLO,' zhi =',ZHI
      
C =========== Setting up Contour Levels and Labeling =================
      IF(mod(LINLOG,2).EQ.1) THEN   ! ==> do linear plot <====================
        IF(NOFCLS.EQ.0) THEN ! standard contouring
          write (*,*) 'STANDARD CONTOURING'
          CALL NICE(ZLO,ZHI,NOZ,NOZM,FMTZ,NFMTZ,.true.)
          NCONT=NOZ*NOZM
          if(b2plot_debug.ge.10)
     1     write (*,*) 'ncont =',ncont,' noz =',noz,' nozm =',nozm
          do while(NCONT+3.GT.no_available_colors(0).AND.NOZM.GT.1)
            NOZM=NOZM/2
            NCONT=NOZ*NOZM
            if(b2plot_debug.ge.10)
     1       write (*,*) 'ncont =',ncont,' noz =',noz,' nozm =',nozm
          ENDDO
          ZINC=(ZHI-ZLO)/(NCONT)
          if(b2plot_debug.ge.10)
     1       write (*,*) 'zlo =',zlo,' zhi =',zhi,' zinc =',zinc
          CALL CPSETI('CLS - CONTOUR LEVEL SELECTION FLAG',0)
          CALL CPSETI('NCL - NUMBER OF CONTOUR LEVELS',NCONT+1)
          DO I=1,NCONT+1
            CALL CPSETI('PAI - PARAMETER ARRAY INDEX',I)
            CALL CPSETR('CLV - CONTOUR LEVEL VALUE',ZLO+(I-1)*ZINC)
            write (*,*) 'contour level',i,
     1             '=',zlo+(i-1)*zinc,'in color',i
            CALL CPSETR('CLL - CONTOUR LINE LINE WIDTH',LWIDTH)
            CALL CPSETI('AIA - AREA ABOVE IDENTIFIER',I+1)
            CALL CPSETI('AIB - AREA BELOW IDENTIFIER',I)
c              IF(NOZM.NE.0) THEN
c               IF(MOD(I,NOZM).EQ.0)
c     1      CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
c             ENDIF
            CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
          ENDDO
        ELSE ! customised contouring
          write (*,*) 'CUSTOMISED CONTOURING'
          NCONT=NOFCLS-1
          if(ncont+3.gt.no_available_colors(0)) then
C do something to avoid this
          endif
          NOZM=0.2*NCONT
          if(b2plot_debug.ge.3)
     1     write (*,*) 'ncont =',ncont
          do i=1,ncont+1
            write (*,*) 'clval(',i,') =',clval(i)
          enddo
          CALL NICE(ZLO,ZHI,NOZ,NOZM,FMTZ,NFMTZ,.true.)
          ZINC=(ZHI-ZLO)/(NCONT)
          if(b2plot_debug.ge.10)
     1       write (*,*) 'zlo =',zlo,' zhi =',zhi,' zinc =',zinc
          CALL CPSETI('CLS - CONTOUR LEVEL SELECTION FLAG',0)
          CALL CPSETI('NCL - NUMBER OF CONTOUR LEVELS',NCONT+1)
          DO I=1,NCONT+1
            CALL CPSETI('PAI - PARAMETER ARRAY INDEX',I)
            CALL CPSETR('CLV - CONTOUR LEVEL VALUE',CLVAL(I))
            write (*,*) 'contour level',i,'=',clval(i),'in color',i
            CALL CPSETR('CLL - CONTOUR LINE LINE WIDTH',LWIDTH)
            CALL CPSETI('AIA - AREA ABOVE IDENTIFIER',I+1)
            CALL CPSETI('AIB - AREA BELOW IDENTIFIER',I)
c              IF(MOD(I,NOZM).EQ.0)
c     1          CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
            CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
          ENDDO
        ENDIF
        WRITE(STRING,'(a,3(1pg12.5,a))')
     1   'contoured from',ZLO,'by',ZINC,'to',ZHI
      ELSE  ! ==> do log plot <=======================================
        IF(NOFCLS.EQ.0) THEN ! standard contouring
          write (*,*) 'STANDARD CONTOURING'
          if(zhi.le.0.0) then
            write(*,*) 'zhi <= 0 and logarithmic axes chosen'
            zhi=1.0
          endif
          if(zlo.le.0.0) zlo=10.0**(-ndec)*zhi
          ZHI=10.0**CEIL(LOG10(ZHI))
          ZLO=max(10.0**FLOOR(LOG10(ZLO)),10.0**(-ndec)*zhi)
          write (*,*) 'zlo =',zlo,' zhi =',zhi
          NCONT=CEIL(LOG10(ZHI))-FLOOR(LOG10(ZLO))
          FMTZ='(1PG7.1)'
          if(ncont*3+3.lt.no_available_colors(0)) then
            ncont=ncont*3
            call cpseti('CLS - CONTOUR LEVEL SELECTION FLAG',0)
            call cpseti('NCL - NUMBER OF CONTOUR LEVELS',ncont+1)
            do i=1,ncont+1
              call cpseti('PAI - PARAMETER ARRAY INDEX',i)
              zval=10.0**((i-1)/3)*zlo
              if(mod(i-1,3).eq.1) zval=zval*2
              if(mod(i-1,3).eq.2) zval=zval*5
              call cpsetr('CLV - CONTOUR LEVEL VALUE',zval)
              CALL CPSETR('CLL - CONTOUR LINE LINE WIDTH',LWIDTH)
              CALL CPSETI('AIA - AREA ABOVE IDENTIFIER',I+1)
              CALL CPSETI('AIB - AREA BELOW IDENTIFIER',I)
              CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
            enddo
          else
            call cpseti('CLS - CONTOUR LEVEL SELECTION FLAG',0)
            call cpseti('NCL - NUMBER OF CONTOUR LEVELS',ncont+1)
            do i=1,ncont+1
              call cpseti('PAI - PARAMETER ARRAY INDEX',i)
              zval=10.0**(i-1)*zlo
              call cpsetr('CLV - CONTOUR LEVEL VALUE',zval)
              CALL CPSETR('CLL - CONTOUR LINE LINE WIDTH',LWIDTH)
              CALL CPSETI('AIA - AREA ABOVE IDENTIFIER',I+1)
              CALL CPSETI('AIB - AREA BELOW IDENTIFIER',I)
              CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
            enddo
          endif
        ELSE ! customised contouring
          write (*,*) 'CUSTOMISED CONTOURING'
          ZLO=10.0**FLOOR(LOG10(ZLO))
          ZHI=10.0**CEIL(LOG10(ZHI))
          NCONT=NOFCLS-1
          FMTZ='(1PG9.3)'
          CALL CPSETI('CLS - CONTOUR LEVEL SELECTION FLAG',0)
          CALL CPSETI('NCL - NUMBER OF CONTOUR LEVELS',NCONT+1)
          DO I=1,NCONT+1
            CALL CPSETI('PAI - PARAMETER ARRAY INDEX',I)
            CALL CPSETR('CLV - CONTOUR LEVEL VALUE',CLVAL(I))
            write (*,*) 'contour level',i,'=',clval(i),' in color',i
            CALL CPSETR('CLL - CONTOUR LINE LINE WIDTH',LWIDTH)
            CALL CPSETI('AIA - AREA ABOVE IDENTIFIER',I+1)
            CALL CPSETI('AIB - AREA BELOW IDENTIFIER',I)
c              IF(MOD(I-1,3).EQ.0)
c     1       CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
            CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
          ENDDO
        ENDIF
        WRITE(STRING,'(a,2(1pg12.5,a))')
     1   'contoured from',ZLO,'to',ZHI
      ENDIF ! ==> end of lin/log plot <===============================
      
      CALL CPSETC('ILT - INFORMATION LABEL TEXT STRING',
     1 STRING(1:LNBLNK(STRING)))
      CALL CPSETR('ILX - INFORMATION LABEL X POSITION',0.5)
      CALL CPSETR('ILY - INFORMATION LABEL Y POSITION',0.95)
      CALL CPSETI('ILP - INFORMATION LABEL POSITION',-3)
      
      SPECIAL=-1E30
      CALL CPSETR('SPV - SPECIAL VALUE',SPECIAL)

C Initialise Areas
      CALL ARINAM(MAP, LMAP)
      CALL GSFAIS(1)

      IF(PHYSICAL) THEN  ! ========== plot physical domain
        IC=1
        X(IC)=XLO
        Y(IC)=YLO
        IC=IC+1
        X(IC)=XLO
        Y(IC)=YHI
        IC=IC+1
        X(IC)=XHI
        Y(IC)=YHI
        IC=IC+1
        X(IC)=XHI
        Y(IC)=YLO
        IC=IC+1
        X(IC)=XLO
        Y(IC)=YLO
        if (b2plot_debug.ge.8) then
          write(*,*) 'Physical domain bounding box : '
          write(*,*) 'X : ', X(1:5)
          write(*,*) 'Y : ', Y(1:5)
        end if
        CALL AREDAM(MAP,X,Y,IC,2,-1,0)
cxpb Obtaining the outer contour outline of the plasma region to be plotted
        IC=1
        IX=-1
        IY=-1
        X(IC)=CR(IX,IY)
        Y(IC)=CZ(IX,IY)
        do while(iy.lt.ny .and. topiy(ix,iy).eq.iy+1)
          IC=IC+1
          iy=iy+1
          X(IC)=CR(IX,IY)
          Y(IC)=CZ(IX,IY)
        end do
        iyt=iy
        if(iy.eq.ny) then
          iniso=.false.
          do while(ix.lt.nx .and. .not.iniso)
            ix=rightix(ix,iy)
            iniso=.false.
            do iiso = 1, nniso
              iniso = iniso .or.
     &         (ix.ge.nxiso1(iiso).and.ix.le.nxiso2(iiso).and.
     &          iy.ge.nyiso1(iiso).and.iy.le.nyiso2(iiso))
            end do
            iniso = iniso .or.
     &       (rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) .or.
     &       (leftix(ix,iy).eq.-2 .and. ix.ne.-1)
            if(.not.iniso) then
              ic=ic+1
              x(ic)=cr(ix,iy)
              y(ic)=cz(ix,iy)
            else
              ix=ix-1
            end if
          end do
        else
          ix=leftix(0,iyt+1)
          ic=ic+1
          x(ic)=cr(ix,iy)
          y(ic)=cz(ix,iy)
        end if
        do iy = iyt-1, -1, -1
          ic=ic+1
          x(ic)=cr(ix,iy)
          y(ic)=cz(ix,iy)
        end do
        iniso = .false.
        iy=0
        do while(ix.gt.-1 .and. .not.iniso)
          iniso=.false.
          if(iyt.eq.ny) then
            ix=leftix(ix,iy)
            do iiso = 1, nniso
              iniso = iniso .or.
     &         (ix.ge.nxiso1(iiso).and.ix.le.nxiso2(iiso).and.
     &          iy.ge.nyiso1(iiso).and.iy.le.nyiso2(iiso))
            end do
            iniso = iniso .or.
     &       (rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) .or.
     &       (leftix(ix,iy).eq.-2 .and. ix.ne.-1)
          else
            if(leftix(ix,iy).eq.ix-1) then
              ix=leftix(ix,iy)
              do iiso = 1, nniso
                iniso = iniso .or.
     &           (ix.ge.nxiso1(iiso).and.ix.le.nxiso2(iiso).and.
     &            iy.ge.nyiso1(iiso).and.iy.le.nyiso2(iiso))
              end do
              iniso = iniso .or.
     &         (rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) .or.
     &         (leftix(ix,iy).eq.-2 .and. ix.ne.-1)
            else
              ix=ix-1
              ic=ic+1
              x(ic)=cr(ix,iy)
              y(ic)=cz(ix,iy)
              do while(leftix(ix,iy).ne.ix-1)
                ic=ic+1
                x(ic)=cr(ix,iy)
                y(ic)=cz(ix,iy)
                iy=iy+1
              end do
              iy=iy-1
              ix=ix-1
              do while(iy.ge.-1)
                ic=ic+1
                x(ic)=cr(ix,iy)
                y(ic)=cz(ix,iy)
                iy=iy-1
              end do
            end if
            iy=-1
          end if
          if(.not.iniso) then
            ic=ic+1
            x(ic)=cr(ix,iy)
            y(ic)=cz(ix,iy)
          else
            ix=ix+1
          end if
        end do
        do while (iniso)
          DO IY=-1,NY
            IC=IC+1
            X(IC)=CR(IX,IY)
            Y(IC)=CZ(IX,IY)
          ENDDO
          iniso=.false.
          iy=ny
          do while(ix.lt.nx .and. .not.iniso)
            ix=rightix(ix,iy)
            iniso=.false.
            do iiso = 1, nniso
              iniso = iniso .or.
     &               (ix.ge.nxiso1(iiso).and.ix.le.nxiso2(iiso).and.
     &                iy.ge.nyiso1(iiso).and.iy.le.nyiso2(iiso))
            end do
            iniso = iniso .or.
     &       (rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) .or.
     &       (leftix(ix,iy).eq.-2 .and. ix.ne.-1)
            if(.not.iniso) then
              ic=ic+1
              x(ic)=cr(ix,iy)
              y(ic)=cz(ix,iy)
            else
              ix=ix-1
            end if
          end do
          do iy = ny-1, -1, -1
            ic=ic+1
            x(ic)=cr(ix,iy)
            y(ic)=cz(ix,iy)
          end do
          iniso = .false.
          iy=-1
          do while(ix.gt.-1 .and. .not.iniso)
            iniso=.false.
            if(iyt.eq.ny) then
              ix=leftix(ix,iy)
              do iiso = 1, nniso
                iniso = iniso .or.
     &                 (ix.ge.nxiso1(iiso).and.ix.le.nxiso2(iiso).and.
     &                  iy.ge.nyiso1(iiso).and.iy.le.nyiso2(iiso))
              end do
              iniso = iniso .or.
     &         (rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) .or.
     &         (leftix(ix,iy).eq.-2 .and. ix.ne.-1)
            else
              if(leftix(ix,iy).eq.ix-1) then
                ix=leftix(ix,iy)
                do iiso = 1, nniso
                  iniso = iniso .or.
     &             (ix.ge.nxiso1(iiso).and.ix.le.nxiso2(iiso).and.
     &              iy.ge.nyiso1(iiso).and.iy.le.nyiso2(iiso))
                end do
                iniso = iniso .or.
     &           (rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) .or.
     &           (leftix(ix,iy).eq.-2 .and. ix.ne.-1)
              else
                ix=ix-1
                ic=ic+1
                x(ic)=cr(ix,iy)
                y(ic)=cz(ix,iy)
                do while(leftix(ix,iy).ne.ix-1)
                  ic=ic+1
                  x(ic)=cr(ix,iy)
                  y(ic)=cz(ix,iy)
                  iy=iy+1
                end do
                iy=iy-1
                ix=ix-1
                do while(iy.ge.-1)
                  ic=ic+1
                  x(ic)=cr(ix,iy)
                  y(ic)=cz(ix,iy)
                  iy=iy-1
                end do
              end if
              iy=-1
            end if
            do iiso = 1, nniso
              iniso = iniso .or.
     &               (ix.ge.nxiso1(iiso).and.ix.le.nxiso2(iiso).and.
     &                iy.ge.nyiso1(iiso).and.iy.le.nyiso2(iiso))
            end do
            iniso = iniso .or.
     &       (rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) .or.
     &       (leftix(ix,iy).eq.-2 .and. ix.ne.-1)
            if(.not.iniso) then
              ic=ic+1
              x(ic)=cr(ix,iy)
              y(ic)=cz(ix,iy)
            else
              ix=ix+1
            end if
          end do
        end do
        CALL AREDAM(MAP,X,Y,IC,2,0,1)
        IF(.NOT.COLPRT.OR..NOT.CONTPLOT) CALL CURVED(X,Y,IC)
cxpb Obtaining the inner contour outline of the plasma region to be plotted
        if(iyt.eq.ny) then
          iy=-1
          ix=-1
        else
          iy=ny
          ix=0
        end if
        do while(rightix(ix,iy).eq.ix+1 .and. iy.eq.-1)
          ix=ix+1
        end do
        if(iy.eq.-1) ix=ix+1
        ixl=ix
        ic=1
        x(ic)=cr(ixl,iy)
        y(ic)=cz(ixl,iy)
        do while(ix.ne.leftix(ixl,iy))
          ix=rightix(ix,iy)
          ic=ic+1
          x(ic)=cr(ix,iy)
          y(ic)=cz(ix,iy)
        end do
        ic=ic+1
        x(ic)=cr(ixl,iy)
        y(ic)=cz(ixl,iy)
        CALL AREDAM(MAP,X,Y,IC,2,0,1)
        IF(.NOT.COLPRT.OR..NOT.CONTPLOT) CALL CURVED(X,Y,IC)
        
        IF(DEBUG) THEN
          do iiso = 1, nniso
            WRITE (*,*) 'nxiso1(',iiso,') =',nxiso1(iiso),
     1                 ' nxiso2(',iiso,') =',nxiso2(iiso)
            WRITE (*,*) 'nyiso1(',iiso,') =',nyiso1(iiso),
     1                 ' nyiso2(',iiso,') =',nyiso2(iiso)
          end do
        ENDIF
        
        CALL CPSETI('GIC - GROUP IDENTIFIER FOR CONTOURLINES',1)
        CALL CPSETI('GIL - GROUP IDENTIFIER FOR LABELS',2)
        CALL CPSETI('GIS - GROUP IDENTIFIER FOR VSTRIPS',6)
        CALL CPSETI('NVS - NUMBER OF VSTRIPS',0)
        
        IF(.NOT.CONTPLOT) GOTO 5150  ! skip contourplot
        
C ========= Initialise Data Array ===========================
cxpb Region 1 is either the whole plasma (nncut=0 case) or
cxpb the westernmost private region. 
        kreg=1
        iy = -1 
        ix = -1
        nxstrt(kreg)=-1
        nystrt(kreg)=-1
        do while(ix.lt.nx .and. rightix(ix,iy).eq.ix+1)
          ix = ix+1
        end do
        nxend(kreg)=ix
        nperd(kreg)=-1
        do while(iy.lt.ny .and. rightix(ix,iy).ne.ix+1)
          iy = iy+1
        end do
        if(iy.ne.ny) iy=iy-1
        nyend(kreg)=iy
        if(iy.eq.ny .and. ix.eq.nx) goto 5250   ! whole plasma case done
cxpb Region 2 is the core plasma or its western half for double-null cases.
        kreg=kreg+1
        iy=-1
        ix=ix+1
        nxstrt(kreg)=ix
        nystrt(kreg)=iy
        do while(rightix(ix,iy).eq.ix+1)
          ix = ix+1
        end do
        nxend(kreg)=ix
        if(iyt.eq.ny) then
          nperd(kreg)=rightix(ix,iy)
        else
          nperd(kreg)=-1
        end if
        do while(rightix(ix,iy).ne.ix+1)
          iy = iy+1
        end do
        nyend(kreg)=iy-1
cxpb Region 3 is the eastern private region for single nulls or the private
cxpb region to the west of the isolating wall in double nulls.
        kreg=kreg+1                                           
        iy=-1
        ix=ix+1
        nxstrt(kreg)=ix
        nystrt(kreg)=iy
        iniso=.false.
        do while(ix.lt.nx .and. .not.iniso)
          ix=ix+1
          iniso=.false.
          do iiso = 1, nniso
            iniso = iniso .or.
     &             (ix.ge.nxiso1(iiso).and.ix.le.nxiso2(iiso).and.
     &              iy.ge.nyiso1(iiso).and.iy.le.nyiso2(iiso))
          end do                                        
          iniso = iniso .or.
     &     (rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) .or.
     &     (leftix(ix,iy).eq.-2 .and. ix.ne.-1)
        end do
        if(iniso) ix=ix-1
        nxend(kreg)=ix
        if(iyt.eq.ny) then
          nperd(kreg)=leftix(nxstrt(kreg),iy)
        else
          nperd(kreg)=-1
        end if
        ix=nxstrt(kreg)
        do while(leftix(ix,iy).ne.ix-1)
          iy = iy+1
        end do
        nyend(kreg)=iy-1
cxpb Region 4 is the region above region 1, i.e. the first SOL in contact with
cxpb the westernmost plate.
        kreg=kreg+1
        ix=-1
        iy=nyend(1)+1
        nxstrt(kreg)=-1
        nystrt(kreg)=nyend(1)+1
        iniso=.false.
        do while(ix.lt.nx .and. .not.iniso .and.
     &           rightix(ix,iy).eq.ix+1)
          ix = ix+1
          iniso=.false.
          do iiso=1, nniso
            iniso = iniso .or.
     &             (ix.ge.nxiso1(iiso).and.ix.le.nxiso2(iiso).and.
     &              iy.ge.nyiso1(iiso).and.iy.le.nyiso2(iiso))
          end do 
          iniso = iniso .or.
     &     (rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) .or.
     &     (leftix(ix,iy).eq.-2 .and. ix.ne.-1)
        end do
        if(iniso) ix=ix-1
        nxend(kreg)=ix
        nperd(kreg)=-1
        if(iniso .or. ix.eq.nx) then
          iy=iyt
        else
          do while(iy.lt.ny .and. rightix(ix,iy).ne.ix+1)
            iy = iy+1
          end do
          iy=iy-1
        end if
        nyend(kreg)=iy
        if(iy.eq.ny .and. ix.eq.nx .and. iyt.eq.ny)
     >   goto 5250  ! single null case done
cxpb Region 5 in the easternmost private region for double-nulls or 
cxpb the island region for island geometries (iyt<ny case)
        kreg=kreg+1
        if(iyt.eq.ny) then
          iy=-1
          ix=nx
          nxend(kreg)=ix
          nystrt(kreg)=iy
          do while(leftix(ix,iy).eq.ix-1)
            ix=ix-1
          end do
          nxstrt(kreg)=ix
          nperd(kreg)=leftix(ix,iy)
          do while(leftix(ix,iy).ne.ix-1)
            iy=iy+1
          end do
          nyend(kreg)=iy-1
        else
          nxstrt(kreg)=0
          nxend(kreg)=nx
          nystrt(kreg)=iyt+1
          nyend(kreg)=ny-1
          nperd(kreg)=nx
        end if
        if(iyt.ne.ny-1) goto 5250   ! island case done
cxpb Region 6 is the eastern half of the private region
        kreg=kreg+1
        ix=ix-1
        iy=-1
        nxend(kreg)=ix
        nystrt(kreg)=iy
        do while(leftix(ix,iy).eq.ix-1)
          ix=ix-1
        end do
        nxstrt(kreg)=ix
        nperd(kreg)=leftix(ix,iy)
        do while(leftix(ix,iy).ne.ix-1)
          iy=iy+1
        end do
        nyend(kreg)=iy-1
cxpb Region 7 is the private region on the east side of the isolating wall
        kreg=kreg+1
        ix=ix-1
        iy=-1
        nxend(kreg)=ix
        nystrt(kreg)=iy
        iniso=.false.
        do while(.not.iniso)
          ix=ix-1
          iniso=.false.
          do iiso = 1, nniso
            iniso = iniso .or.
     &       (ix.ge.nxiso1(iiso).and.ix.le.nxiso2(iiso).and.
     &        iy.ge.nyiso1(iiso).and.iy.le.nyiso2(iiso))
          end do
          iniso = iniso .or.
     &     (rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) .or.
     &     (leftix(ix,iy).eq.-2 .and. ix.ne.-1)
        end do
        ix=ix+1
        nxstrt(kreg)=ix
        nperd(kreg)=-1
        ix=nxend(kreg)
        do while(rightix(ix,iy).ne.ix+1)
          iy=iy+1
        end do
        nyend(kreg)=iy-1
cxpb Region 8 is the SOL section above region 5
        kreg=kreg+1
        ix=nx
        iy=nyend(5)+1
        iniso=.false.
        do while(.not.iniso .and. leftix(ix,iy).eq.ix-1)
          ix=ix-1
          iniso=.false.
          do iiso = 1, nniso
            iniso = iniso .or.
     &             (ix.ge.nxiso1(iiso).and.ix.le.nxiso2(iiso).and.
     &              iy.ge.nyiso1(iiso).and.iy.le.nyiso2(iiso))
          end do
          iniso = iniso .or.
     &     (rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) .or.
     &     (leftix(ix,iy).eq.-2 .and. ix.ne.-1)
        end do
        if(iniso) ix=ix+1
        nxstrt(kreg)=ix
        if(iniso) then
          nperd(kreg)=-1
          nyend(kreg)=ny
        else
          nperd(kreg)=leftix(ix,iy)
          do while(leftix(ix,iy).ne.ix-1)
            iy=iy+1
          end do
          nyend(kreg)=iy-1
        endif
        if(nyend(5).eq.nyend(7)) goto 5250  ! Connected double null case done
cxpb Region 9 is the missing SOL section on the western half
        kreg=kreg+1
        nxend(kreg)=nxend(3)
        if(nxend(1).gt.nxend(3)) then
          nxstrt(kreg)=nxstrt(2)
          nystrt(kreg)=nystrt(2)+1
          nyend(kreg)=nyend(1)
          nperd(kreg)=leftix(nxstrt(kreg),nystrt(kreg))
        else
          nxstrt(kreg)=-1
          nystrt(kreg)=nyend(4)+1
          nyend(kreg)=ny
          nperd(kreg)=-1
        end if
cxpb Region 10 is the missing SOL section in the eastern half
        kreg=kreg+1
        nxstrt(kreg)=nxstrt(7)
        if(nxend(7).lt.nxend(5)) then
          nxend(kreg)=nxend(6)
          nystrt(kreg)=nystrt(6)+1
          nyend(kreg)=nyend(5)
          nperd(kreg)=leftix(nxstrt(kreg),nystrt(kreg))
        else
          nxend(kreg)=nx
          nystrt(kreg)=nyend(8)+1
          nyend(kreg)=ny
          nperd(kreg)=-1
        end if
C =============================== Set up individual zones
 5250   continue
        do l=1,3
          do k=1,kreg
            ic=0
            do iy=nystrt(k),nyend(k)
              if (nperd(k).ne.-1) then
                ic=ic+1
                x(ic)=  cr(nperd(k),iy)
                y(ic)=  cz(nperd(k),iy)
                z(ic)= field(nperd(k),iy)
              endif
              do ix=nxstrt(k),nxend(k)
                ic=ic+1
                x(ic)= cr(ix,iy)
                y(ic)= cz(ix,iy)
                z(ic)=field(ix,iy)
              enddo
            enddo
            nnx=nxend(k)-nxstrt(k)+1
            if (nperd(k).ne.-1) nnx=nnx+1
            nny=nyend(k)-nystrt(k)+1
            if (b2plot_debug.eq.8) write(*,*) 'k, nnx, nny : ',
     .                                         k, nnx, nny
            call xertst (ic.eq.nnx*nny, 'ic.ne.nnx*nny')
            CALL CPSETI('SET - DO-SET-CALL FLAG',0)
            CALL CPRECT(z, NNX, NNX, NNY, RWRK, LRWK, IWRK, LIWK)
            IF(L.EQ.1) THEN
C Add contours to area map
              CALL CPCLAM(z, RWRK, IWRK, MAP)
C Add label boxes to area map
              IF(LABELS) THEN
                CALL CPLBAM(z, RWRK, IWRK, MAP)
              ENDIF
c            ELSEIF(L.EQ.2) THEN
C Add label boxes to area map
c             IF(LABELS) THEN
c               CALL CPLBAM(z, RWRK, IWRK, MAP)
c             ENDIF
            ELSEIF(L.EQ.3) THEN
C Draw Labels
              IF(LABELS) THEN
                CALL CPLBDR(z,RWRK,IWRK)
              ENDIF
C Draw Contours
              IF(CONTOURS) THEN
C               CALL CPCLDR(z,RWRK,IWRK)
                CALL CPCLDM(z,RWRK,IWRK,MAP,CPDRPL)
              ENDIF
            ENDIF
          ENDDO
          IF(L.EQ.2) THEN
C Set fill style to solid, and fill contours if COLPRT==TRUE
            IF(COLPRT) THEN
              call entsr(irold,1)
              CALL ARPRAM(MAP,1,0,0)
              if(nerro(nerrf).ne.0) then
                call eprin
                call errof
                call retsr(irold)
                call subend ()
                return
              endif
              call retsr(irold)
              CALL ARSCAM(MAP,XWRK,YWRK,NWRK,IAREA,IGRP,NOGRPS,FILL)
            ENDIF
          ENDIF
        ENDDO
        
 5150   CONTINUE
        
      ELSE      ! ======== plot computational domain
        IC=1
        X(IC)=XLO
        Y(IC)=YLO
        IC=IC+1
        X(IC)=XLO
        Y(IC)=YHI
        IC=IC+1
        X(IC)=XHI
        Y(IC)=YHI
        IC=IC+1
        X(IC)=XHI
        Y(IC)=YLO
        IC=IC+1
        X(IC)=XLO
        Y(IC)=YLO
        if (b2plot_debug.ge.8) then
          write(*,*) 'Computational domain bounding box : '
          write(*,*) 'X : ', X(1:5)
          write(*,*) 'Y : ', Y(1:5)
        end if
        CALL AREDAM(MAP,X,Y,IC,2,-1,0)
        IC=1
        X(IC)=max(CR(-1,-1),values(1))
        Y(IC)=max(CZ(-1,-1),values(3))
        IC=IC+1
        X(IC)=max(CR(-1,-1),values(1))
        Y(IC)=min(CZ(NZ,0)+0.5_R4,values(4))
        IC=IC+1
        X(IC)=min(CR(NR,0)+0.5_R4,values(2))
        Y(IC)=min(CZ(NZ,0)+0.5_R4,values(4))
        IC=IC+1
        X(IC)=min(CR(NR,0)+0.5_R4,values(2))
        Y(IC)=min(CZ(-1,-1),values(3))
        IC=IC+1
        X(IC)=min(CR(-1,-1),values(1))
        Y(IC)=min(CZ(-1,-1),values(3))
        if (b2plot_debug.ge.8) then
          write(*,*) 'Computational domain data bounds : '
          write(*,*) 'X : ', X(1:5)
          write(*,*) 'Y : ', Y(1:5)
        end if
        CALL AREDAM(MAP,X,Y,IC,2,0,1)
        IF(.NOT.COLPRT.OR..NOT.CONTPLOT) CALL CURVED(X,Y,IC)
        
        CALL CPSETI('GIC - GROUP IDENTIFIER FOR CONTOURLINES',1)
        CALL CPSETI('GIL - GROUP IDENTIFIER FOR LABELS',2)
        CALL CPSETI('GIS - GROUP IDENTIFIER FOR VSTRIPS',6)
        CALL CPSETI('NVS - NUMBER OF VSTRIPS',0)
        
        IF(.NOT.CONTPLOT) GOTO 5152
        
        IC=0
        ixstrt=max(-1,min(nx,nint(values(1))))
        iystrt=max(-1,min(ny,nint(values(3))))
        ixend =max(-1,min(nx,nint(values(2))))
        iyend =max(-1,min(ny,nint(values(4))))
        DO IY=iystrt,iyend
          DO IX=ixstrt,ixend
            IC=IC+1
            X(IC)=ix    ! was CR(IX,0)
            Y(IC)=iy    ! was CZ(IY,0)
            Z(IC)=FIELD(IX,IY)
          ENDDO
        ENDDO
        NNX=ixend-ixstrt+1
        NNY=iyend-iystrt+1
        IF(IC.NE.NNX*NNY) WRITE(*,*) 'ic <> nnx*nny',IC,NNX,NNY
        IF (b2plot_debug.ge.8) WRITE(*,*) 'NX,NY,NNX,NNY',NX,NY,NNX,NNY
C         CALL CPCNRC(z,NNX,NNX,NNY,ZLO,ZHI,0.0,1,0,0)
        CALL CPSETI('SET - DO-SET-CALL FLAG',0)
        CALL CPSETI('MAP - MAPPING FLAG',1)
        CALL CPSETR('XC1 - X COORDINATE AT I=1', real(ixstrt))
        CALL CPSETR('XCM - X COORDINATE AT I=M', real(ixend))
        CALL CPSETR('YC1 - Y COORDINATE AT J=1', real(iystrt))
        CALL CPSETR('YCN - Y COORDINATE AT J=N', real(iyend))
        CALL CPRECT(z, nnx, nnx, nny, RWRK, LRWK, IWRK, LIWK)
C Add contours to area map
        CALL CPCLAM(z, RWRK, IWRK, MAP)
C Add label boxes to area map
        IF(LABELS) THEN
          CALL CPLBAM(z, RWRK, IWRK, MAP)
        ENDIF
C Set fill style to solid, and fill contours if COLPRT==TRUE
        IF(COLPRT) THEN
          call entsr(irold,1)
          CALL ARPRAM(MAP,1,0,0)
          if(nerro(nerrf).ne.0) then
            call eprin
            call errof
            call retsr(irold)
            call subend ()
            return
          endif
          call retsr(irold)
          CALL ARSCAM(MAP,XWRK,YWRK,NWRK,IAREA,IGRP,NOGRPS,FILL)
        ENDIF
C Draw Labels
        IF(LABELS) THEN
          CALL CPLBDR(z,RWRK,IWRK)
        ENDIF
C Draw Contours
        IF(CONTOURS) THEN
C            CALL CPCLDR(z,RWRK,IWRK)
          CALL CPCLDM(z,RWRK,IWRK,MAP,CPDRPL)
        ENDIF
        
 5152   CONTINUE

      ENDIF
      
      IF(DEBUG) THEN
        DO I=1,2
          WRITE(DEBUGSTRING,'(I5)') I
          CALL ARDBPA(MAP,I,DEBUGSTRING)
        ENDDO
      ENDIF
      
      IF(LBLBR.AND.CONTPLOT.and.colprt) THEN
        DO I=1,NCONT+1
          LFIN(I)=choose_color(i-1,0)
          CALL CPSETI('PAI - PARAMETER ARRAY INDEX',I)
          CALL CPGETR('CLV - CONTOUR LEVEL VALUES',CLV)
          if(nofcls.eq.0) then ! standard contouring
            WRITE (LBLS(I), FMTZ) CLV
          else
            if(abs(clv).ge.1e3.or.abs(clv).lt.1e-2) then
              write(LBLS(I),'(1p,e10.3)') CLV
            elseif(abs(clv).ge.1e2) then
              write(LBLS(I),'(f8.3)') CLV
            elseif(abs(clv).ge.1e1) then
              write(LBLS(I),'(f8.4)') CLV
            else
              write(LBLS(I),'(f8.5)') CLV
            endif
          endif
          CALL SHRINK_LABEL(LBLS(I))
        ENDDO
        i=ncont+2
        LFIN(I)=choose_color(i-1,0)
        LBLS(I:100)=' '
        CALL GSLWSC(1.0)
        if(noax.gt.0) then
          if(rlo.gt.clval(1).and.rhi.le.clval(ncont+1)) then
            CALL LBLBAR(1,
     1       colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2       0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3       NCONT,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),1,LBLS,
     +       NCONT+1,1)
          elseif(rlo.le.clval(1).and.rhi.le.clval(ncont+1)) then
            do i=ncont+2,2,-1
              lbls(i)=lbls(i-1)
            enddo
            lbls(1)=' '
            CALL LBLBAR(1,
     1       colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2       0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3       NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     +       NCONT+2,1)
          elseif(rlo.gt.clval(1).and.rhi.gt.clval(ncont+1)) then
            CALL LBLBAR(1,
     1       colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2       0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3       NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),0,LBLS,
     +       NCONT+2,1)
          else
            CALL LBLBAR(1,
     1       colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2       0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3       NCONT+2,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     +       NCONT+1,1)
          endif
        endif
      ENDIF
      
C Draw Mesh/Separatrix if wanted
      if(wall) call data_copy(is)
      call additional_elements(nx,ny,physical)
      
      call subend ()
      RETURN
      END SUBROUTINE PLOTCONTOUR
