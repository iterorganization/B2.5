c
C =====================================================================
      subroutine plttarget
C =====================================================================
      use b2mod_types
      use b2mod_b2plot_debug
      use b2mod_wall
      use b2mod_b2plot
     * ,only : target_stack, itwrk
      implicit none
#include "graphics_vars.h"

#ifndef PLPLOT      
      integer nst,nnd
      parameter(nst=16,nnd=24)
      integer ind(nnd)
      real dst(nst)
#endif
      real (kind=R8) :: zlod,zhid
      real (kind=R4) :: zlo,zhi,rlo,rhi,zinc,zval,clv
      real (kind=R4) :: xbox(8),ybox(8)
      integer i, ibnd, ic, icol, ici, n, noz, nozm, nfmtz
      INTEGER LFIN(100),choose_color
      CHARACTER*20 LBLS(100)
      CHARACTER*20 FMTZ
      real (kind=R4) :: ceil
      integer  no_available_colors
      external no_available_colors, ceil

      if(b2plot_debug.ge.10) 
     > write(*,*) 'plttarget: nbnd, ndepth_nml = ',nbnd, ndepth_nml
      zlod=trgmin
      zhid=trgmax
      if(b2plot_debug.ge.10) 
     > write(*,*) 'plttarget: zlo, zhi = ',zlod,zhid
      if(zlod.eq.zhid.and.zlod.eq.0.0) then
        zlod=target_stack(1,1,itwrk)
        zhid=target_stack(1,1,itwrk)
        do ibnd=1,nbnd
          do i = 1,ndepth_nml
            zlod=min(zlod,target_stack(ibnd,i,itwrk))
            zhid=max(zhid,target_stack(ibnd,i,itwrk))
          enddo
        enddo
      endif
      zlo=zlod
      rlo=zlod
      zhi=zhid
      rhi=zhid
      if(zlo.eq.zhi) then
        zlo=min(zlo*0.9,zlo*1.1)
        zhi=max(zhi*0.9,zhi*1.1)
        if(zlo.eq.zhi) then
          zlo=-1.0
          zhi=+1.0
        end if
      end if
      if(b2plot_debug.ge.10) write(*,*) 'plttarget: zlo, zhi = ',zlo,zhi
      if(b2plot_debug.ge.12) write(*,*) 'trglinlog = ',trglinlog

      IF(NOFTCLS.EQ.0) THEN                ! standard target contouring
        IF(trglinlog.EQ.1) THEN  ! ==> do linear plot <====================
          CALL NICE(ZLO,ZHI,NOZ,NOZM,FMTZ,NFMTZ)
          NCONT=NOZ*NOZM
          if(b2plot_debug.ge.10)
     1     write (*,*) 'ncont =',ncont,' noz =',noz,' nozm =',nozm
          do while(NCONT+3.GT.no_available_colors(1).AND.NOZM.GT.1)
            NOZM=NOZM/2
            NCONT=NOZ*NOZM
            if(b2plot_debug.ge.10)
     1       write (*,*) 'ncont =',ncont,' noz =',noz,' nozm =',nozm
          ENDDO
          ZINC=(ZHI-ZLO)/(NCONT)
          if(b2plot_debug.ge.10)
     1     write (*,*) 'zlo =',zlo,' zhi =',zhi,'zinc =',zinc
          DO I=1,NCONT+1
            tclval(i)=ZLO+(I-1)*ZINC
          ENDDO
        ELSE  ! ==> do log plot <=======================================
          write (*,*) 'STANDARD TARGET CONTOURING'
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
          if(ncont*3+3.lt.no_available_colors(1)) then
            ncont=ncont*3
#ifndef PLPLOT
            call cpseti('CLS - CONTOUR LEVEL SELECTION FLAG',0)
            call cpseti('NCL - NUMBER OF CONTOUR LEVELS',ncont+1)
#endif
            do i=1,ncont+1
#ifndef PLPLOT
              call cpseti('PAI - PARAMETER ARRAY INDEX',i)
#endif
              zval=10.0**((i-1)/3)*zlo
              if(mod(i-1,3).eq.1) zval=zval*2
              if(mod(i-1,3).eq.2) zval=zval*5
              tclval(i)=zval
            enddo
          else
            do i=1,ncont+1
              zval=10.0**(i-1)*zlo
              tclval(i)=zval
            enddo
          endif
        ENDIF
      ELSE                                ! customised target contouring
        NCONT=NOFTCLS-1
        DO I=1,NCONT+1
          write (*,*) 
     .     'target contour level',i,'=',tclval(i),' in color',i
        ENDDO
      ENDIF

      DO I=1,NCONT+1
        LFIN(I)=choose_color(i-1,1)
        clv=tclval(i)
        if(noftcls.eq.0) then ! standard target contouring
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
      LFIN(I)=choose_color(i-1,1)
      LBLS(I)=' '
#ifndef PLPLOT
      CALL GSLWSC(1.0)
      call gsfais(1)
      if (b2plot_debug.ge.12) then
        write(*,*) 'noax = ',noax
        write(*,*) 'targetbar: left, right = ',
     .              targetbar_left,targetbar_right
        write(*,*) 'XSIZE,YSIZE = ',XSIZE,YSIZE
        write(*,*) 'XPOS ,YPOS  = ',XPOS ,YPOS
        write(*,*) 'NCONT,LBSZ  = ',NCONT,LBSZ
        write(*,*) 'LFIN, LBLS  = ',LFIN(1:NCONT+2),LBLS(1:NCONT+1)
      endif
      if(noax.gt.0) then
        if(rlo.le.tclval(1).and.rhi.gt.tclval(ncont+1)) then
          CALL LBLBAR(1,
     1 targetbar_left*XSIZE+XPOS,targetbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT+2,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     + NCONT+1,1)
        elseif(rlo.gt.tclval(1).and.rhi.le.tclval(ncont+1)) then
          CALL LBLBAR(1,
     1 targetbar_left*XSIZE+XPOS,targetbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),0,LBLS,
     + NCONT+1,1)
        elseif(rlo.le.tclval(1).and.rhi.le.tclval(ncont+1)) then
          do i=ncont+2,2,-1
            lbls(i)=lbls(i-1)
          enddo
          lbls(1)=' '
          CALL LBLBAR(1,
     1 targetbar_left*XSIZE+XPOS,targetbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     + NCONT+2,1)
        elseif(rlo.gt.tclval(1).and.rhi.gt.tclval(ncont+1)) then
          CALL LBLBAR(1,
     1 targetbar_left*XSIZE+XPOS,targetbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),0,LBLS,
     + NCONT+2,1)
        endif
      endif
#endif

      do i=1,ndepth_nml
        do ibnd=1,nbnd
          do ic=1,2
            xbox(ic)=targetcr(ibnd,i,ic)
            ybox(ic)=targetcz(ibnd,i,ic)
          enddo
          do ic=1,2
            xbox(2+ic)=targetcr(ibnd,min(i+1,ndepth_nml),ic)
            ybox(2+ic)=targetcz(ibnd,min(i+1,ndepth_nml),ic)
          enddo
          n=4
          call clip_box(xbox,ybox,xlo,ylo,xhi,yhi,n)
          if(n.gt.2) then
            if(target_stack(ibnd,i,itwrk).gt.tclval(ncont+1)) then
              icol=ncont+1
            else
              icol=0
              do while(target_stack(ibnd,i,itwrk).gt.tclval(icol+1))
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
      enddo
#ifdef PLPLOT
      call plcol(1)
#endif
C Draw element boundaries
      if(sptrix) then
        call targetelements(.true.)
      endif
C Draw Target elements outline
      if(outline) then
        call targetoutline(.true.)
      endif
C Draw Mesh/Vessel if wanted
      if(grid) then
        call targetmesh(.true.)
      endif
      return
      end subroutine plttarget
