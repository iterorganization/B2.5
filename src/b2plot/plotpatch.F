C =====================================================================
      subroutine plotpatch(MPG,GEO,IS,
     . work,ZLOD,ZHID,PHYSICAL,LINLOG,igraphs)
C =====================================================================
      use b2mod_types
      use b2us_map
      use b2us_geo
      use b2mod_b2plot_debug
      implicit none
#include "graphics_vars.h"
      type(mapping) :: mpg
      type(geometry) :: geo
      integer is,linlog,igraphs
      REAL (KIND=R8) :: work(mpg%nCv),ZLOD,ZHID
      LOGICAL PHYSICAL

      INTEGER LFIN(100),choose_color
      CHARACTER*20 LBLS(100)
      CHARACTER*20 FMTZ

#ifndef PLPLOT      
      integer nst,nnd
      parameter(nst=16,nnd=24)
      integer ind(nnd)
      real (kind=R4) :: dst(nst)
#endif
      real (kind=R4) :: xbox(8),ybox(8),
     1 zlo,zhi,zval,clv,zinc,rlo,rhi,small
      integer iCv,icol,ici,n,i,ic,
     2 no_available_colors,noz,nozm,nfmtz
      real floor,ceil
      external floor,ceil

      zlo=real(zlod,R4)
      zhi=real(zhid,R4)
      RLO=real(work(1),R4)
      RHI=real(work(1),R4)
      do iCv=1,mpg%nCv
        if(.not.physical.or.mpg%cvReg(iCv).ne.0) then
          RLO=MIN(RLO,REAL(work(iCv),R4))
          RHI=MAX(RHI,REAL(work(iCv),R4))
        endif
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

      IF(NOFCLS.EQ.0) THEN                ! standard contouring
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
            clval(i)=ZLO+(I-1)*ZINC
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
              if(mod(i-1,3).eq.1) zval=zval*2.0_R4
              if(mod(i-1,3).eq.2) zval=zval*5.0_R4
              clval(i)=zval
            enddo
          else
            do i=1,ncont+1
              zval=10.0_R4**(i-1)*zlo
              clval(i)=zval
            enddo
          endif
        ENDIF
      ELSE                                ! customised contouring
        NCONT=NOFCLS-1
        DO I=1,NCONT+1
          write (*,*) 'contour level',i,'=',clval(i),' in color',i
        ENDDO
      ENDIF

      DO I=1,NCONT+1
        LFIN(I)=choose_color(i-1,0)
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
      I=NCONT+2
      LFIN(I)=choose_color(i-1,0)
      LBLS(I:100)=' '
#ifndef PLPLOT
      CALL GSLWSC(1.0)
      call gsfais(1)
      if(noax.gt.0) then
        if((rlo.le.clval(1).and.rhi.gt.clval(ncont+1)).or.
     &   igraphs.gt.1) then
          CALL LBLBAR(1,
     1 colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT+2,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     + NCONT+1,1)
        elseif(rlo.gt.clval(1).and.rhi.le.clval(ncont+1)) then
          CALL LBLBAR(1,
     1 colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),0,LBLS,
     + NCONT+1,1)
        elseif(rlo.le.clval(1).and.rhi.le.clval(ncont+1)) then
          do i=ncont+2,2,-1
            lbls(i)=lbls(i-1)
          enddo
          lbls(1)=' '
          CALL LBLBAR(1,
     1 colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     + NCONT+2,1)
        elseif(rlo.gt.clval(1).and.rhi.gt.clval(ncont+1)) then
          CALL LBLBAR(1,
     1 colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),0,LBLS,
     + NCONT+2,1)
        endif
      endif
#endif

      do iCv=1,mpg%nCv
        if(.not.physical.or.mpg%cvReg(iCv).ne.0) then
cxpb Do not plot dead cells in physical domain
#ifdef WG_TODO
          if(physical) then
            do ic=0,1
              xbox(ic+1)=real(crx(ix,iy,ic),R4)
              ybox(ic+1)=real(cry(ix,iy,ic),R4)
            enddo
            do ic=2,3
              xbox(6-ic)=real(crx(ix,iy,ic),R4)
              ybox(6-ic)=real(cry(ix,iy,ic),R4)
            enddo
          else
            xbox(1)=ix-0.5_R4
            xbox(2)=ix+1-0.5_R4
            xbox(3)=ix+1-0.5_R4
            xbox(4)=ix-0.5_R4
            ybox(1)=iy-0.5_R4
            ybox(2)=iy-0.5_R4
            ybox(3)=iy+1-0.5_R4
            ybox(4)=iy+1-0.5_R4
          endif
          n=4
#endif
          call clip_box(xbox,ybox,xlo,ylo,xhi,yhi,n)
          if(n.gt.2) then
            if(work(iCv).gt.clval(ncont+1)) then
              icol=ncont+1
            else
              icol=0
              do while(work(iCv).gt.clval(icol+1))
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
        endif
      enddo
#ifdef PLPLOT
      call plcol(1)
#endif
C Draw Mesh/Separatrix if wanted
      if(wall) call data_copy(is)
      call additional_elements(mpg%nCv,physical)
      return
      end subroutine plotpatch
c
      subroutine plot_gca(colia,icolia,dimx,dimy,
     1 dpx,dpy,dqx,dqy,dzmin,dzmax,linlog)
      use b2mod_types
      use b2mod_b2plot_debug
      implicit none
#include <graphics_vars.h>

      integer dimx,dimy,linlog
      real (kind=R8) :: colia(dimx*dimy)
      integer icolia(dimx*dimy)
      real (kind=R8) :: dpx,dpy,dqx,dqy,dzmin,dzmax

#ifndef PLPLOT
      real (kind=R4) :: px,py,qx,qy,zmin,zmax
      real (kind=R4) :: clvl(100),sclvl(100)
      real (kind=R4) :: clv,zval,zinc,ceil,floor
      external ceil,floor
      integer noz,nozm,nfmtz,ncnt,i
      integer isclvl(100)
      CHARACTER*20 LBLS(100)
      character*20 fmtz

      zmin=real(dzmin,R4)
      zmax=real(dzmax,R4)
      if(zmin.eq.zmax) then
        zmin=real(colia(1),R4)
        zmax=real(colia(1),R4)
        do i=2,dimx*dimy
          zmin=min(zmin,real(colia(i),R4))
          zmax=max(zmax,real(colia(i),R4))
        enddo
      endif
      write(*,*) 'zmin,zmax = ',zmin,zmax
      if(mod(linlog,2).eq.0) then              ! logarithmic
        if(zmax.le.0.0_R4) then
          write(*,*)
     1     'Non-positive max value for log plot in process_camera'
          return
        endif
        zmax=10.0_R4**ceil(log10(zmax))
        if(zmin.gt.0.0_R4) then
          zmin=max(10.0_R4**floor(log10(zmin)),10.0_R4**(-ndec)*zmax)
        else
          zmin=10.0_R4**(-ndec)*zmax
        endif
        write(*,*) 'zmin,zmax = ',zmin,zmax
        NCNT=INT(CEIL(LOG10(ZMAX))-FLOOR(LOG10(ZMIN)))
        write(*,*) 'ncnt = ',ncnt
        if(ncnt.lt.5) then
          ncnt=ncnt*3
          do i=1,ncnt+1
            zval=10.0_R4**((i-1)/3)*zmin
            if(mod(i-1,3).eq.1) zval=zval*2
            if(mod(i-1,3).eq.2) zval=zval*5
            clvl(i)=zval
          enddo
        else
          do i=1,ncnt+1
            zval=10.0_R4**(i-1)*zmin
            clvl(i)=zval
          enddo
        endif
        fmtz='(1pg8.2)'
        do i=1,dimx*dimy
          colia(i)=(log10(max(real(colia(i),R4),zmin))-log10(zmin))/
     1     (log10(zmax)-log10(zmin))
        enddo      
        do i=1,ncnt+1
          sclvl(i)=(log10(max(clvl(i),zmin))-log10(zmin))/
     1     (log10(zmax)-log10(zmin))
        enddo
      else
        call nice(zmin,zmax,noz,nozm,fmtz,nfmtz)
        NCNT=NOZ*NOZM
        write (*,*) 'ncnt =',ncnt,' noz =',noz,' nozm =',nozm
        do while(NCNT.lt.5.AND.NOZM.GT.1)
          NOZM=NOZM/2
          NCNT=NOZ*NOZM
          write (*,*) 'ncnt =',ncnt,' noz =',noz,' nozm =',nozm
        ENDDO
        ZINC=(ZMAX-ZMIN)/(NCNT)
        if(b2plot_debug.ge.10)
     1   write (*,*) 'zmin =',zmin,' zmax =',zmax,'zinc =',zinc
        DO I=1,NCNT+1
          clvl(i)=ZMIN+(I-1)*ZINC
        ENDDO
        do i=1,dimx*dimy
          colia(i)=(colia(i)-zmin)/(zmax-zmin)
        enddo
        do i=1,ncnt+1
          sclvl(i)=(clvl(i)-zmin)/(zmax-zmin)
        enddo
      endif
      DO I=1,NCNT+1
        clv=clvl(i)
        WRITE (LBLS(I), FMTZ) CLV
        CALL SHRINK_LABEL(LBLS(I))
      ENDDO
      LBLS(NCNT+2:100)=' '
      write(*,*) 'zmin,zmax = ',zmin,zmax

      px=real(dpx,R4)
      py=real(dpy,R4)
      qx=real(dqx,R4)
      qy=real(dqy,R4)
      do i=1,dimx*dimy
        icolia(i)=int(min(1.0_R8,max(0.0_R8,colia(i)))*(ncol-1))+2
      enddo
      call gca(px,py,qx,qy,dimx,dimy,
     1   1,1,dimx,dimy,icolia)

      do i=1,ncnt+1
        isclvl(i)=int(min(1.0_R4,max(0.0_R4,sclvl(i)))*(ncol-1))+2
      enddo
      if(noax.gt.0) CALL LBLBAR(1,
     1 colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 NCNT+1,0.2/(0.2+0.8*lbsz),1.0,isclvl,0,LBLS,
     + NCNT+1,1)
#endif
      return
      end subroutine plot_gca
c
      subroutine specieslblbar(ns,isclvl,lbls)
      use b2mod_types
      implicit none
#include <graphics_vars.h>

      integer ns,isclvl(ns)
      character*(*) lbls(ns)
      character*256 tmplbl(ns)

      integer i

      do i=1,ns
        isclvl(i)=i
        call customise_label(lbls(i),tmplbl(i))
        lbls(i)=trim(tmplbl(i))
      enddo
#ifndef PLPLOT
      if(noax.gt.0) call lblbar(1,
     1 colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2 0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3 ns,0.2/(0.2+0.8*lbsz),1.0,isclvl,0,lbls,
     + ns,1)
#endif
      return
      end subroutine specieslblbar
