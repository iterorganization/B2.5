c
C =====================================================================
      subroutine pltwall(nCv,physical,from_target)
C =====================================================================
      use b2mod_types
      use b2mod_indirect
      use b2mod_b2plot_ank
      use b2mod_b2plot_debug
      use b2mod_wall
      implicit none
      integer nCv
      logical physical, from_target
#include "graphics_vars.h"

      real (kind=R8) :: zlod,zhid
      real (kind=R4) :: zlo,zhi,rlo,rhi,zinc,zval,clv
      real (kind=R8) :: wallx1,wallx2,wally1,wally2
      integer i, icol, ici, noz, nozm, nfmtz
      INTEGER LFIN(100),choose_color
      CHARACTER*20 LBLS(100)
      CHARACTER*20 FMTZ
      real (kind=R4) :: ceil, sfloor
      integer  no_available_colors
      external no_available_colors, ceil, sfloor

      if(b2plot_debug.ge.10) write(*,*) 'pltwall: ndata = ',ndata
      zlod=wmin
      zhid=wmax
      if(b2plot_debug.ge.10) write(*,*) 'pltwall: zlo, zhi = ',zlod,zhid
      if(zlod.eq.zhid.and.zlod.eq.0.0) then
        zlod=data(1)
        zhid=data(1)
        do i=2,ndata
          zlod=min(zlod,data(i))
          zhid=max(zhid,data(i))
        enddo
      endif
      zlo=real(zlod,R4)
      rlo=real(zlod,R4)
      zhi=real(zhid,R4)
      rhi=real(zhid,R4)
      if(zlo.eq.zhi) then
        zlo=min(zlo*0.9_R4,zlo*1.1_R4)
        zhi=max(zhi*0.9_R4,zhi*1.1_R4)
        if(zlo.eq.zhi) then
          zlo=-1.0_R4
          zhi=+1.0_R4
        end if
      end if
      if(b2plot_debug.ge.10) write(*,*) 'pltwall: zlo, zhi = ',zlo,zhi
      if(b2plot_debug.ge.12) write(*,*) 'wlinlog = ',wlinlog

      IF(NOFWCLS.EQ.0) THEN                ! standard wall contouring
        IF(wlinlog.EQ.1) THEN  ! ==> do linear plot <====================
          CALL NICE(ZLO,ZHI,NOZ,NOZM,FMTZ,NFMTZ,.true.)
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
            wclval(i)=ZLO+(I-1)*ZINC
          ENDDO
        ELSE  ! ==> do log plot <=======================================
          write (*,*) 'STANDARD WALL CONTOURING'
          if(zhi.le.0.0_R4) then
            write(*,*) 'zhi <= 0 and logarithmic axes chosen'
            zhi=1.0_R4
          endif
          if(zlo.le.0.0_R4) zlo=10.0_R4**(-ndec)*zhi
          ZHI=10.0_R4**CEIL(LOG10(ZHI))
          ZLO=max(10.0_R4**SFLOOR(LOG10(ZLO)),10.0_R4**(-ndec)*zhi)
          write (*,*) 'zlo =',zlo,' zhi =',zhi
          NCONT=INT(CEIL(LOG10(ZHI))-SFLOOR(LOG10(ZLO)))
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
              zval=10.0_R4**((i-1)/3)*zlo
              if(mod(i-1,3).eq.1) zval=zval*2.0_R4
              if(mod(i-1,3).eq.2) zval=zval*5.0_R4
              wclval(i)=zval
            enddo
          else
            do i=1,ncont+1
              zval=10.0_R4**(i-1)*zlo
              wclval(i)=zval
            enddo
          endif
        ENDIF
      ELSE                                ! customised wall contouring
        NCONT=NOFWCLS-1
        DO I=1,NCONT+1
          write (*,*) 'wall contour level',i,'=',wclval(i),' in color',i
        ENDDO
      ENDIF

      DO I=1,NCONT+1
        LFIN(I)=choose_color(i-1,1)
        clv=wclval(i)
        if(nofwcls.eq.0) then ! standard wall contouring
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
      LFIN(I)=choose_color(i-1,1)
      LBLS(I:100)=' '
#ifndef PLPLOT
      CALL GSLWSC(1.0)
      call gsfais(1)
      if (b2plot_debug.ge.12) then
        write(*,*) 'noax = ',noax
        write(*,*) 'wallbar: left, right = ',wallbar_left,wallbar_right
        write(*,*) 'XSIZE,YSIZE = ',XSIZE,YSIZE
        write(*,*) 'XPOS ,YPOS  = ',XPOS ,YPOS
        write(*,*) 'NCONT,LBSZ  = ',NCONT,LBSZ
        write(*,*) 'LFIN, LBLS  = ',LFIN(1:NCONT+2),LBLS(1:NCONT+1)
      endif
      if(noax.gt.0) then
        if(rlo.le.wclval(1).and.rhi.gt.wclval(ncont+1)) then
          CALL LBLBAR(1,
     1 wallbar_left*XSIZE+XPOS,wallbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT+2,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     + NCONT+1,1)
        elseif(rlo.gt.wclval(1).and.rhi.le.wclval(ncont+1)) then
          CALL LBLBAR(1,
     1 wallbar_left*XSIZE+XPOS,wallbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),0,LBLS,
     + NCONT+1,1)
        elseif(rlo.le.wclval(1).and.rhi.le.wclval(ncont+1)) then
          do i=ncont+2,2,-1
            lbls(i)=lbls(i-1)
          enddo
          lbls(1)=' '
          CALL LBLBAR(1,
     1 wallbar_left*XSIZE+XPOS,wallbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     + NCONT+2,1)
        elseif(rlo.gt.wclval(1).and.rhi.gt.wclval(ncont+1)) then
          CALL LBLBAR(1,
     1 wallbar_left*XSIZE+XPOS,wallbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),0,LBLS,
     + NCONT+2,1)
        endif
      endif

      call lastd
      call gslwsc(wwidth)
      do i=1,ndata
        if(data(i).gt.wclval(ncont+1)) then
          icol=ncont+1
        else
          icol=0
          do while(data(i).gt.wclval(icol+1))
            icol=icol+1
          enddo
        endif
        ici=choose_color(icol,1)
        call lastd
        call gsplci(ici)
        if (physical) then
          call lined_d(r1(i),z1(i),r2(i),z2(i))
#ifdef WG_TODO
        else if (.not.from_target) then
          wallx1 = real(imapx(i),R4)
          wallx2 = real(imapx(i),R4)
          wally1 = real(imapy(i),R4)
          wally2 = real(imapy(i),R4)
          if (bottomiy(imapx(i),imapy(i)).eq.-2) then
            wallx1 = wallx1 - 0.5_R4
            wallx2 = wallx2 + 0.5_R4
            wally1 = wally1 - 0.5_R4
            wally2 = wally2 - 0.5_R4
          else if (topiy(imapx(i),imapy(i)).eq.ny+1) then
            wallx1 = wallx1 - 0.5_R4
            wallx2 = wallx2 + 0.5_R4
            wally1 = wally1 + 0.5_R4
            wally2 = wally2 + 0.5_R4
          endif
          if (leftix(imapx(i),imapy(i)).eq.-2) then
            wallx1 = wallx1 - 0.5_R4
            wallx2 = wallx2 - 0.5_R4
            wally1 = wally1 - 0.5_R4
            wally2 = wally2 + 0.5_R4
          elseif (rightix(imapx(i),imapy(i)).eq.nx+1) then
            wallx1 = wallx1 + 0.5_R4
            wallx2 = wallx2 + 0.5_R4
            wally1 = wally1 - 0.5_R4
            wally2 = wally2 + 0.5_R4
          endif
          call lined_d(wallx1,wally1,wallx2,wally2)
#endif
        else if (from_target) then
          wallx1 = i - 0.5_R4
          wallx2 = i + 0.5_R4
          wally1 = 1 - 0.5_R4
          wally2 = 1 - 0.5_R4
          call lined_d(wallx1,wally1,wallx2,wally2)
        endif
      enddo
      call lastd
      call gsplci(1)
      call gslwsc(1.0)
#endif
      return
      end subroutine pltwall
