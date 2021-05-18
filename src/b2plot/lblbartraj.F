
C =====================================================================
      subroutine lblbartraj(zlod,zhid,linlog)
C =====================================================================
      use b2mod_types
      use b2mod_b2plot_debug
      implicit none
#include "graphics_vars.h"
      real (kind=R8) :: ZLOD,ZHID,ZVALIN
      integer linlog
      real zlo, zhi
      REAL ZINC, ZVAL, CLV
      INTEGER I, ICOL, NOZ, NOZM, NFMTZ
      REAL CEIL, FLOOR
      EXTERNAL CEIL, FLOOR
      INTRINSIC LOG10

      INTEGER LFIN(100)
      CHARACTER*20 LBLS(100)
      CHARACTER*20 FMTZ

      save zlo,zhi

      zlo=real(zlod,R4)
      zhi=real(zhid,R4)
      
      if(zlo.eq.zhi) then
        zlo=min(zlo*0.9_R4,zlo*1.1_R4)
        zhi=max(zhi*0.9_R4,zhi*1.1_R4)
        if(zlo.eq.zhi) then
          zlo=-1.0_R4
          zhi=+1.0_R4
        endif
      endif

      if(b2plot_debug.ge.10) write(*,*) 'zlo, zhi = ',zlo,zhi

      IF(NOFCLS.EQ.0) THEN    ! standard contouring
        IF(LINLOG.EQ.1) THEN  ! ==> do linear plot <====================
          CALL NICE(ZLO,ZHI,NOZ,NOZM,FMTZ,NFMTZ,.true.)
          NCONT=NOZ*NOZM
          if(b2plot_debug.ge.10)
     1     write (*,*) 'ncont =',ncont,' noz =',noz,' nozm =',nozm
          do while(NCONT+3.GT.NCOL.AND.NOZM.GT.1)
            NOZM=NOZM/2
            NCONT=NOZ*NOZM
            if(b2plot_debug.ge.10)
     1       write (*,*) 'ncont =',ncont,' noz =',noz,' nozm =',nozm
          ENDDO
          ZINC=(ZHI-ZLO)/(NCONT)
          if(b2plot_debug.ge.10)
     1     write (*,*) 'zlo =',zlo,' zhi =',zhi,'zinc =',zinc
          DO I=1,NCONT+1
            clval(i)=ZLO+(I-1)*ZINC
          ENDDO
        ELSE  ! ==> do log plot <=======================================
          if(b2plot_debug.ge.9) write (*,*) 'STANDARD CONTOURING'
          if(zhi.le.0.0_R4) then
            write(*,*) 'zhi <= 0 and logarithmic axes chosen'
            zhi=1.0_R4
          endif
          if(zlo.le.0.0_R4) zlo=10.0_R4**(-ndec)*zhi
          ZHI=10.0_R4**CEIL(LOG10(ZHI))
          ZLO=max(10.0_R4**FLOOR(LOG10(ZLO)),10.0_R4**(-ndec)*zhi)
          if(b2plot_debug.ge.10) write (*,*) 'zlo =',zlo,' zhi =',zhi
          NCONT=INT(CEIL(LOG10(ZHI))-FLOOR(LOG10(ZLO)))
          FMTZ='(1PG7.1)'
          if(ncont*3+3.lt.ncol) then
            ncont=ncont*3
            call cpseti('CLS - CONTOUR LEVEL SELECTION FLAG',0)
            call cpseti('NCL - NUMBER OF CONTOUR LEVELS',ncont+1)
            do i=1,ncont+1
              call cpseti('PAI - PARAMETER ARRAY INDEX',i)
              zval=10.0_R4**((i-1)/3)*zlo
              if(mod(i-1,3).eq.1) zval=zval*2
              if(mod(i-1,3).eq.2) zval=zval*5
              clval(i)=zval
            enddo
          else
            do i=1,ncont+1
              zval=10.0_R4**(i-1)*zlo
              clval(i)=zval
            enddo
          endif
        ENDIF
      ELSE                             ! customised contouring
        NCONT=NOFCLS-1
        DO I=1,NCONT+1
          if(b2plot_debug.ge.4)
     1     write (*,*) 'contour level',i,'=',clval(i),' in color',i
        ENDDO
      ENDIF

      DO I=1,NCONT+1
        if(fill_option.eq.0) then
          LFIN(I)=(I-1)*NCOL/(ncont+3)+2
        else
          LFIN(I)=MOD(I-1,NCOL)+2
        endif
        clv=clval(i)
        if(nofcls.eq.0) then ! standard contouring
          WRITE (LBLS(I), FMTZ) CLV
        else
          if(abs(clv).ge.1.0e3_R4.or.abs(clv).lt.1.0e-2_R4) then
            write(LBLS(I),'(1p,e10.3)') CLV
          elseif(abs(clv).ge.1.0e2_R4) then
            write(LBLS(I),'(f8.3)') CLV
          elseif(abs(clv).ge.1.0e1_R4) then
            write(LBLS(I),'(f8.4)') CLV
          else
            write(LBLS(I),'(f8.5)') CLV
          endif
        endif
        CALL SHRINK_LABEL(LBLS(I))
      ENDDO
      i=ncont+2
      LBLS(I:100)=' '
      if(fill_option.eq.0) then
        LFIN(I)=(I-1)*NCOL/(ncont+3)+2
      else
        LFIN(I)=MOD(I-1,NCOL)+2
      endif
      CALL GSLWSC(1.0)
      call gsfais(1)
      if(noax.gt.0) CALL LBLBAR(1,
     1 colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT+2,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     + NCONT+1,1)

      return

      entry lblbartrajcol(ZVALIN,icol)

      if(zvalin.gt.clval(ncont+1)) then
        icol=ncont+1
      else
        icol=0
        do while(zvalin.gt.clval(icol+1))
          icol=icol+1
        enddo
      endif
      if(fill_option.eq.0) then
        icol=(icol)*NCOL/(ncont+3)+2
      else
        icol=MOD(icol,NCOL)+2
      endif

      return

      end subroutine lblbartraj
