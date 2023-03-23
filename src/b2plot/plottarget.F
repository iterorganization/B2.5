C =====================================================================
      subroutine plottarget(nx,ny,nbnd,ndepth,
     . targetcr,targetcz,target_stack,
     . ZLOD,ZHID,PHYSICAL,linlog,igraphs)
C =====================================================================
      use b2mod_types
      use b2mod_b2plot_debug
      implicit none
#include "graphics_vars.h"
      integer nx,ny,nbnd,ndepth,igraphs,linlog
      REAL (KIND=R8) :: 
     * targetcr(1:nbnd,1:ndepth,1:2),targetcz(1:nbnd,1:ndepth,1:2),
     1 target_stack(1:nbnd,1:ndepth),ZLOD,ZHID
      LOGICAL PHYSICAL

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
     1 zlo,zhi,zval,clv,zinc,rlo,rhi
      integer ibnd,i,icol,ici,n,ic,
     2 no_available_colors,noz,nozm,nfmtz
      real floor,ceil
      external floor,ceil

      zlo=real(zlod,R4)
      zhi=real(zhid,R4)
      RLO=real(target_stack(1,1),R4)
      RHI=real(target_stack(1,1),R4)
      do I=1,ndepth
        do ibnd=1,nbnd
          if(.not.physical)then
            RLO=MIN(RLO,REAL(target_stack(ibnd,i),R4))
            RHI=MAX(RHI,REAL(target_stack(ibnd,i),R4))
          endif
        end do
      end do
      if(ZLO.EQ.ZHI) then
        ZLO=RLO
        ZHI=RHI
      end if

      if(zlo.eq.zhi) then
        zlo=min(zlo*0.9_R4,zlo*1.1_R4)
        zhi=max(zhi*0.9_R4,zhi*1.1_R4)
        if(zlo.eq.zhi) then
          zlo=-1.0_R4
          zhi=+1.0_R4
        end if
      end if

      if(b2plot_debug.ge.10) write(*,*) 'zlo, zhi = ',zlo,zhi

      IF(NOFTCLS.EQ.0) THEN                ! standard contouring
        IF(mod(LINLOG,2).EQ.1) THEN  ! ==> do linear plot <====================
          CALL NICE(ZLO,ZHI,NOZ,NOZM,FMTZ,NFMTZ)
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
            tclval(i)=ZLO+(I-1)*ZINC
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
          NCONT=INT(CEIL(LOG10(ZHI))-FLOOR(LOG10(ZLO)))
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
              if(mod(i-1,3).eq.1) zval=zval*2
              if(mod(i-1,3).eq.2) zval=zval*5
              tclval(i)=zval
            enddo
          else
            do i=1,ncont+1
              zval=10.0_R4**(i-1)*zlo
              tclval(i)=zval
            enddo
          endif
        ENDIF
      ELSE                                ! customised contouring
        NCONT=NOFTCLS-1
        DO I=1,NCONT+1
          write (*,*) 'target contour level',i,'=',
     .     tclval(i),' in color',i
        ENDDO
      ENDIF

      DO I=1,NCONT+1
        LFIN(I)=choose_color(i-1,0)
        clv=tclval(i)
        if(noftcls.eq.0) then ! standard contouring
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
        if((rlo.le.tclval(1).and.rhi.gt.tclval(ncont+1)).or.
     &   igraphs.gt.1) then
          CALL LBLBAR(1,
     1 colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT+2,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     + NCONT+1,1)
        elseif(rlo.gt.tclval(1).and.rhi.le.tclval(ncont+1)) then
          CALL LBLBAR(1,
     1 colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),0,LBLS,
     + NCONT+1,1)
        elseif(rlo.le.tclval(1).and.rhi.le.tclval(ncont+1)) then
          do i=ncont+2,2,-1
            lbls(i)=lbls(i-1)
          enddo
          lbls(1)=' '
          CALL LBLBAR(1,
     1 colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     + NCONT+2,1)
        elseif(rlo.gt.tclval(1).and.rhi.gt.tclval(ncont+1)) then
          CALL LBLBAR(1,
     1 colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),0,LBLS,
     + NCONT+2,1)
        endif
      endif
#endif

      do i=1,ndepth
        do ibnd=1,nbnd
          if(physical) then
            do ic=1,2
              xbox(ic)=real(targetcr(ibnd,i,ic),R4)
              ybox(ic)=real(targetcz(ibnd,i,ic),R4)
            enddo
            do ic=1,2
              xbox(5-ic)=real(targetcr(ibnd,min(i+1,ndepth),ic),R4)
              ybox(5-ic)=real(targetcz(ibnd,min(i+1,ndepth),ic),R4)
            enddo
          endif
          if(.not.physical) then
            xbox(1)=ibnd-0.5_R4
            xbox(2)=ibnd+1-0.5_R4
            xbox(3)=ibnd+1-0.5_R4
            xbox(4)=ibnd-0.5_R4
            ybox(1)=i-0.5_R4
            ybox(2)=i-0.5_R4
            ybox(3)=i+1-0.5_R4
            ybox(4)=i+1-0.5_R4
          endif
          n=4
          call clip_box(xbox,ybox,xlo,ylo,xhi,yhi,n)
          if(n.gt.2) then
            if(target_stack(ibnd,i).gt.tclval(ncont+1)) then
              icol=ncont+1
            else
              icol=0
              do while(target_stack(ibnd,i).gt.tclval(icol+1))
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
        call targetelements(physical)
      endif
C Draw Target elements outline
      if(outline) then
        call targetoutline(physical)
      endif
C Draw Mesh/Vessel if wanted
      if(grid) then
        call targetmesh(physical)
      endif
      if(vessel.and.physical) then
        call pltvessel
      endif
C Draw Wall elements if wanted
      if(wall) then
        call data_copy(0)
        call pltwall(nx,ny,physical,.true.)
      endif
      return
      end subroutine plottarget
c
      subroutine targetmesh(physical)
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      implicit none
      logical physical
      real (kind=R8) ::
     *   x1,x2,y1,y2
      integer ibnd,i

      if(mcol .gt. 1) then
#ifndef PLPLOT
        call lastd
        call gsplci(mcol)
#endif
      endif
      if(.not.physical) then ! mesh computational domain
        do ibnd=1,nbnd+1
          x1=ibnd-0.5_R8
          x2=ibnd-0.5_R8
          y1=1-0.5_R8
          y2=ndepth+1-0.5_R8
          call lined_d(x1,y1,x2,y2)
        enddo
        do i=1,ndepth+1
          x1=1-0.5_R8
          x2=nbnd+1-0.5_R8
          y1=i-0.5_R8
          y2=i-0.5_R8
          call lined_d(x1,y1,x2,y2)
        enddo
      endif
      if(physical) then ! mesh physical domain
        do ibnd=1,nbnd
          do i=1,ndepth
            call lined_d(targetcr(ibnd,i,1),targetcz(ibnd,i,1),
     1       targetcr(ibnd,i,2),targetcz(ibnd,i,2))  
          enddo
          do i=1,ndepth-1
            call lined_d(targetcr(ibnd,i,1),targetcz(ibnd,i,1),
     1       targetcr(ibnd,i+1,1),targetcz(ibnd,i+1,1))  
            call lined_d(targetcr(ibnd,i,2),targetcz(ibnd,i,2),
     1       targetcr(ibnd,i+1,2),targetcz(ibnd,i+1,2))  
          enddo
        enddo
      endif
      if(mcol .gt. 1) then
#ifndef PLPLOT
        call lastd
        call gsplci(1)
#endif
      endif
      return
      end subroutine targetmesh
c
      subroutine targetoutline(physical)
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      implicit none
      logical physical
      real (kind=R8) ::
     *   x1,x2,y1,y2
      integer ibnd,i,iseg

#ifndef PLPLOT
      call lastd
      if(ocol .gt. 1) call gsplci(ocol)
      call gslwsc_d(owidth)
#endif
      if(.not.physical) then ! outline computational domain
        x1=1-0.5_R8
        x2=1-0.5_R8
        y1=1-0.5_R8
        y2=ndepth+1-0.5_R8
        call lined_d(x1,y1,x2,y2)
        x1=nbnd+1-0.5_R8
        x2=nbnd+1-0.5_R8
        y1=1-0.5_R8
        y2=ndepth+1-0.5_R8
        call lined_d(x1,y1,x2,y2)
        x1=1-0.5_R8
        x2=nbnd+1-0.5_R8
        y1=1-0.5_R8
        y2=1-0.5_R8
        call lined_d(x1,y1,x2,y2)
        x1=1-0.5_R8
        x2=nbnd+1-0.5_R8
        y1=ndepth+1-0.5_R8
        y2=ndepth+1-0.5_R8
        call lined_d(x1,y1,x2,y2)
      endif
      if(physical) then ! outline physical domain
        do ibnd=1,nbnd
          call lined_d(targetcr(ibnd,1,1),targetcz(ibnd,1,1),
     1       targetcr(ibnd,1,2),targetcz(ibnd,1,2))  
          call lined_d(targetcr(ibnd,ndepth,1),targetcz(ibnd,ndepth,1),
     1       targetcr(ibnd,ndepth,2),targetcz(ibnd,ndepth,2))
        enddo  
        do iseg=1,nsegments
          do i=1,ndepth-1
            call lined_d(
     1       targetcr(wall_begin(iseg),i,1),
     2       targetcz(wall_begin(iseg),i,1),
     1       targetcr(wall_begin(iseg),i+1,1),
     2       targetcz(wall_begin(iseg),i+1,1))  
            call lined_d(
     1       targetcr(wall_end(iseg),i,2),
     2       targetcz(wall_end(iseg),i,2),
     1       targetcr(wall_end(iseg),i+1,2),
     2       targetcz(wall_end(iseg),i+1,2))  
          enddo
        enddo
      endif
#ifndef PLPLOT
      call lastd
      if(ocol .gt. 1) call gsplci(1)
      call gslwsc_d(lwidth)
#endif
      return
      end subroutine targetoutline
c
      subroutine targetelements(physical)
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      implicit none
      logical physical
      real (kind=R8) ::
     *   x1,x2,y1,y2
      integer i,iseg

#ifndef PLPLOT
      call lastd
      if(scol .gt. 1) call gsplci(scol)
      call gslwsc_d(sepwidth)
#endif
      if(.not.physical) then ! computational domain
        do iseg=2,nsegments
          x1=wall_begin(iseg)-0.5_R8
          x2=wall_begin(iseg)-0.5_R8
          y1=1-0.5_R8
          y2=ndepth+1-0.5_R8
          call lined_d(x1,y1,x2,y2)
        enddo
      endif
      if(physical) then ! mesh physical domain
        do iseg=1,nsegments
          do i=1,ndepth-1
            call lined_d(
     1       targetcr(wall_begin(iseg),i,1),
     2       targetcz(wall_begin(iseg),i,1),
     1       targetcr(wall_begin(iseg),i+1,1),
     2       targetcz(wall_begin(iseg),i+1,1))  
            call lined_d(
     1       targetcr(wall_end(iseg),i,2),
     2       targetcz(wall_end(iseg),i,2),
     1       targetcr(wall_end(iseg),i+1,2),
     2       targetcz(wall_end(iseg),i+1,2))  
          enddo
        enddo
      endif
#ifndef PLPLOT
      call lastd
      if(scol .gt. 1) call gsplci(1)
      call gslwsc_d(lwidth)
#endif
      return
      end subroutine targetelements
