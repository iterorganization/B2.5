C =====================================================================
      subroutine plottripatch(nCv,ntrii,nrknot,is,
     . xtrian,ytrian,necke,work,ZLOD,ZHID,LINLOG)
C =====================================================================
      use b2mod_types
      use b2mod_b2plot_debug
      implicit none
#include "graphics_vars.h"
      integer nCv,ntrii,nrknot,is,linlog
      integer necke(3,ntrii)
      REAL (KIND=R8) :: XTRIAN(NRKNOT),YTRIAN(NRKNOT),
     1 work(NTRII),ZLOD,ZHID

      INTEGER LFIN(100),choose_color
      CHARACTER*20 LBLS(100)
      CHARACTER*20 FMTZ

#ifndef PLPLOT      
      integer nst,nnd
      parameter(nst=16,nnd=24)
      integer ind(nnd)
      real dst(nst)
#endif
      real xbox(8),ybox(8),
     1 zlo,zhi,zval,clv,zinc,rlo,rhi,small
      integer icol,ici,n,i,
     2 no_available_colors,noz,nozm,nfmtz
      real floor,ceil
      external floor,ceil

      zlo=real(zlod,R4)
      zhi=real(zhid,R4)
      RLO=real(work(1),R4)
      RHI=real(work(1),R4)
      do I=1,NTRII
        RLO=MIN(RLO,REAL(work(I),R4))
        RHI=MAX(RHI,REAL(work(I),R4))
      end do
      if(ZLO.EQ.ZHI) then
        ZLO=RLO
        ZHI=RHI
      end if

      small = tiny(1.0_R4)
      if (abs(zlo).lt.small) zlo = 0.0_R4
      if (abs(zhi).lt.small) zhi = 0.0_R4

      if(zlo.eq.zhi) then
        zlo=min(zlo*0.9_R4,zlo*1.1_R4)
        zhi=max(zhi*0.9_R4,zhi*1.1_R4)
        if(zlo.eq.zhi) then
          zlo=-1.0_R4
          zhi=+1.0_R4
        end if
      end if

      if(b2plot_debug.ge.10) write(*,*) 'zlo, zhi = ',zlo,zhi

      IF(NOFECLS.EQ.0) THEN                ! standard contouring
        IF(mod(LINLOG,2).EQ.1) THEN  ! ==> do linear plot <====================
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
     1       write (*,*) 'zlo =',zlo,' zhi =',zhi,'zinc =',zinc
          DO I=1,NCONT+1
            eclval(i)=ZLO+(I-1)*ZINC
          ENDDO
        ELSE  ! ==> do log plot <=======================================
          if(b2plot_debug.ge.10) write (*,*) 'STANDARD CONTOURING'
          if(zhi.le.0.0_R4) then
            write(*,*) 'zhi <= 0 and logarithmic axes chosen'
            zhi=1.0_R4
          endif
          if(zlo.le.0.0_R4) zlo=10.0_R4**(-ndec)*zhi
          ZHI=10.0_R4**CEIL(LOG10(ZHI))
          ZLO=max(10.0_R4**FLOOR(LOG10(ZLO)),10.0_R4**(-ndec)*zhi)
          if(b2plot_debug.ge.10) write (*,*) 'zlo =',zlo,' zhi =',zhi
          NCONT=int(CEIL(LOG10(ZHI))-FLOOR(LOG10(ZLO)))
          FMTZ='(1PG7.1)'
          if(ncont*3+3.lt.no_available_colors(0)) then
            ncont=ncont*3
#ifndef PLPLOT
            call cpseti('CLS - CONTOUR LEVEL SELECTION FLAG',0)
            call cpseti('NCL - NUMBER OF CONTOUR LEVELS',ncont+1)
#endif
            do i=1,ncont+1
#ifndef PLPLOT
              call cpseti('PAI - PARAMETER ARRAY INDEX',i)
#endif
              zval=10.0_R4**((i-1)/3)*zlo
              if(mod(i-1,3).eq.1) zval=zval*2.0_R4
              if(mod(i-1,3).eq.2) zval=zval*5.0_R4
              eclval(i)=zval
            enddo
          else
            do i=1,ncont+1
              zval=10.0_R4**(i-1)*zlo
              eclval(i)=zval
            enddo
          endif
        ENDIF
      ELSE                                ! customised contouring
        NCONT=NOFECLS-1
        DO I=1,NCONT+1
          write (*,*) 'contour level',i,'=',eclval(i),' in color',i
        ENDDO
      ENDIF

      DO I=1,NCONT+1
        LFIN(I)=choose_color(i-1,0)
        clv=eclval(i)
        if(nofecls.eq.0) then ! standard contouring
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
      I=NCONT+2
      LFIN(I)=choose_color(i-1,0)
      LBLS(I:100)=' '
#ifndef PLPLOT
      CALL GSLWSC(1.0)
      call gsfais(1)
      if(noax.gt.0) then
        if(rlo.le.eclval(1).and.rhi.gt.eclval(ncont+1)) then
          CALL LBLBAR(1,
     1     colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2     0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3     NCONT+2,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     +     NCONT+1,1)
        elseif(rlo.gt.eclval(1).and.rhi.le.eclval(ncont+1)) then
          CALL LBLBAR(1,
     1     colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2     0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3     NCONT,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),0,LBLS,
     +     NCONT+1,1)
        elseif(rlo.le.eclval(1).and.rhi.le.eclval(ncont+1)) then
          do i=ncont+2,2,-1
            lbls(i)=lbls(i-1)
          enddo
          lbls(1)=' '
          CALL LBLBAR(1,
     1     colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2     0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3     NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     +     NCONT+2,1)
        elseif(rlo.gt.eclval(1).and.rhi.gt.eclval(ncont+1)) then
          CALL LBLBAR(1,
     1     colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2     0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3     NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),0,LBLS,
     +     NCONT+2,1)
        endif
      endif
#endif

      do i=1,ntrii
        xbox(1)=real(xtrian(necke(1,i)),R4)
        ybox(1)=real(ytrian(necke(1,i)),R4)
        xbox(2)=real(xtrian(necke(2,i)),R4)
        ybox(2)=real(ytrian(necke(2,i)),R4)
        xbox(3)=real(xtrian(necke(3,i)),R4)
        ybox(3)=real(ytrian(necke(3,i)),R4)
        n=3
        call clip_box(xbox,ybox,xlo,ylo,xhi,yhi,n)
        if(n.gt.2) then
          if(work(i).gt.eclval(ncont+1)) then
            icol=ncont+1
          else
            icol=0
            do while(work(i).gt.eclval(icol+1))
              icol=icol+1
            enddo
          endif
          ici=choose_color(icol,0)
#ifdef PLPLOT
          call plcol(ici)
          call plfill(n,xbox,ybox)
#else
          call sfsgfa(xbox,ybox,n,dst,nst,ind,nnd,ici)
#endif
        endif
      enddo
#ifdef PLPLOT
      call plcol(1)
#endif
C Draw Mesh/Separatrix if wanted
      if(wall) call data_copy(is)
      call additional_elements(nCv,.true.)
      return
      end subroutine plottripatch
