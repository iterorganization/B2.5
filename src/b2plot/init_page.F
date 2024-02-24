
      subroutine init_page(title,xlod,xhid,ylod,yhid,linlog,
     1 aspect,colorbar,physical)
      use b2mod_types
      use b2mod_b2plot_debug
      use b2mod_wall
      implicit none
#include <graphics_vars.h>

c     the zoom common
      integer winid
      integer wsnshwin,wsnX11
      common /zoom_com/ winid,wsnshwin,wsnX11

c     the line format common
      character*10 linefmt
      common /line_format/ linefmt

      character*(*) title
      character*20 fmtx,fmty
      real (kind=R8) :: xlod,xhid,ylod,yhid
      logical aspect,colorbar, physical
      integer linlog
      real (kind=R4) :: range,xmid,ymid,ceil,txtwdth,yframe,
     1 reserve_bottom,reserve_top,reserve_left,reserve_right,
     2 test,test1
      integer ix,iy,n_page,itxtwdth,offset,
     1 nox,noxm,nfmtx,nnox,
     2 noy,noym,nfmty,nnoy
      real (kind=R8) :: rratio
      real (kind=R4) :: sfloor
      logical streql
      external streql, rratio, sfloor

      if(b2plot_debug.ge.8) write(*,*) xwidth,yheight

      if(.not.initialised_plotting) call initialise_plotting
#ifndef PLPLOT
      if(wsnshwin.ne.0) call ngpict(wsnshwin,1)
#endif
      if(same_page.eq.0) then
        ngraph=0
        linefmt='$$$$$$$$'
#ifndef PLPLOT
        call dashdc('$',1,0)
#endif
      else
        if(same_page.eq.2) ngraph=0
        ngraph=ngraph+1
        if(b2plot_debug.ge.8) write(*,*) 'ngraph = ',ngraph
        write(linefmt,'(a,i2,a)') '$$$''',ngraph,'''$$$'
#ifndef PLPLOT
        if(label_same_lines) call dashdc(linefmt,8,lfsz)
#endif
      endif
      if(same_page.eq.1) then
        if(b2plot_debug.ge.6) write(*,*) 'same_page set to ',same_page
        return
      endif
      if(dirty_graph.or.(new_page.and.dirty_page)) call end_page
      xlo=real(min(1.0e30_R8,max(-1.0e30_R8,xlod)),R4)
      xhi=real(min(1.0e30_R8,max(-1.0e30_R8,xhid)),R4)
      ylo=real(min(1.0e30_R8,max(-1.0e30_R8,ylod)),R4)
      yhi=real(min(1.0e30_R8,max(-1.0e30_R8,yhid)),R4)
      if(xlo.eq.xhi) then
        xlo=min(xlo*0.9_R4,xlo*1.1_R4)
        xhi=max(xhi*0.9_R4,xhi*1.1_R4)
        if(xlo.eq.xhi) then
          xlo=-1.0_R4
          xhi=+1.0_R4
        endif
      endif
      if(ylo.eq.yhi) then
        ylo=min(ylo*0.9_R4,ylo*1.1_R4)
        yhi=max(yhi*0.9_R4,yhi*1.1_R4)
        if(ylo.eq.yhi) then
          ylo=-1.0_R4
          yhi=+1.0_R4
        endif
      endif
      if(b2plot_debug.ge.8) write(*,*) xlo,xhi,ylo,yhi
#ifdef PLPLOT
      if(aspect.and.want_aspect) then
         call plenv(xlo,xhi,ylo,yhi,1,0)
      else
         call plenv(xlo,xhi,ylo,yhi,0,0)
      endif
#else
      if(aspect.and.want_aspect) then
        range=max((xhi-xlo)*nx_page/xwidth,(yhi-ylo)*ny_page/yheight)
        xmid=(xhi+xlo)/2.0_R4
        xlo=xmid-range/2.0_R4/nx_page*xwidth
        xhi=xmid+range/2.0_R4/nx_page*xwidth
        ymid=(yhi+ylo)/2.0_R4
        ylo=ymid-range/2.0_R4/ny_page*yheight
        yhi=ymid+range/2.0_R4/ny_page*yheight
      endif
      yframe=yheight
      if(.not.streql(globalheader,' ')) then
        yframe=yframe*(1.00_R4-0.02_R4*glsz)
      endif
      xsize=xwidth/real(nx_page,R4)
      ysize=yframe/real(ny_page,R4)
      xpos=(ix_page-1)*xsize
      ypos=(ny_page-iy_page)*ysize

      ix=0
      iy=0

      if(noax.gt.0) then
        if(linlog.eq.1.or.linlog.eq.2) then
          logx=.false.
          call nice(xlo,xhi,nox,noxm,fmtx,nfmtx,nice_axes)
        else
          logx=.true.
          if(xhi.le.0.0_R4) then
            write(*,*) 'xhi <= 0 and logarithmic axes chosen'
            xhi=1.0_R4
          endif
          xlo=max(xlo,10.0_R4**(-ndec)*xhi)
          xlo=10.0_R4**sfloor(log10(xlo))
          xhi=10.0_R4**ceil(log10(xhi))
          nox=1
          fmtx='(1pg8.2)'
          nfmtx=8
          noxm=1
        endif
        if(lbsz.gt.1.0_R4.and.nox.gt.5) then
          nnox=nint(ceil(real(nox,R4)/2.0_R4))
          noxm=noxm*2
          xhi=xlo+(xhi-xlo)*real(nnox,R4)*2.0_R4/real(nox,R4)
          nox=nnox
          ix=ix+1
        endif

        if(linlog.eq.1.or.linlog.eq.3) then
          logy=.false.
          call nice(ylo,yhi,noy,noym,fmty,nfmty,nice_axes)
        else
          logy=.true.
          if(yhi.le.0.0_R4) then
            write(*,*) 'yhi <= 0 and logarithmic axes chosen'
            yhi=1.0_R4
          endif
          yhi=max(yhi,-ylo)
          ylo=max(ylo,10.0_R4**(-ndec)*yhi)
          ylo=10.0_R4**sfloor(log10(ylo))
          yhi=10.0_R4**ceil(log10(yhi))
          noy=1
          fmty='(1pg8.2)'
          nfmty=8
          noym=1
        endif
        if(lbsz.gt.1.0_R4.and.noy.gt.5) then
          nnoy=nint(ceil(real(noy,R4)/2.0_R4))
          noym=noym*2
          yhi=ylo+(yhi-ylo)*real(nnoy,R4)*2.0_R4/real(noy,R4)
          noy=nnoy
          iy=iy+1
        endif
      endif
      if(b2plot_debug.ge.8) write(*,*) xlo,xhi,ylo,yhi

      n_page=max(nx_page,ny_page)
      itxtwdth=max(nint(12*lbsz/N_PAGE),4)
      txtwdth=real(rratio(itxtwdth,1024),R4)
      offset=0
      if(nx_page.gt.ny_page) offset=1
      if(nx_page/ny_page.gt.3) offset=offset+((nx_page/ny_page)-1)/2
      if(ny_page/nx_page.gt.3) offset=offset+1

      if(noax.gt.0) then
        reserve_bottom=txtwdth*5/ysize
        reserve_top=txtwdth*(4+2*offset)/ysize
        reserve_left=txtwdth*(nfmty+3)/xsize
        reserve_right=0.0_R4
        if(target.and.physical) then
          targetbar_right=1.0_R4-reserve_right
          reserve_right=reserve_right+0.08_R4*0.2_R4+0.08_R4*0.8_R4*lbsz
          targetbar_left=1.0_R4-reserve_right
          reserve_right=reserve_right+0.02_R4
        endif
        if(wall) then
          wallbar_right=1.0_R4-reserve_right
          reserve_right=reserve_right+0.08_R4*0.2_R4+0.08_R4*0.8_R4*lbsz
          wallbar_left=1.0_R4-reserve_right
          reserve_right=reserve_right+0.02_R4
        endif
        if(colorbar) then
          colorbar_right=1.0_R4-reserve_right
          reserve_right=reserve_right+0.08_R4*0.2_R4+0.08_R4*0.8_R4*lbsz
          colorbar_left=1.0_R4-reserve_right
          reserve_right=reserve_right+0.02_R4
        endif
        if(.not.streql(extralabel,' ')) reserve_top=reserve_top+
     .    txtwdth*2*(1+offset+nx_page/5)/ysize
        if(reserve_right.eq.0.0_R4) reserve_right=0.05_R4*lbsz
        if(.not.streql(xlabel,' ')) reserve_bottom =
     .   reserve_bottom+txtwdth*2*n_page
        if(.not.streql(ylabel,' ')) reserve_left =
     .   reserve_left+txtwdth*2*n_page
      else
        reserve_bottom=0.0_R4
        reserve_top=0.0_R4
        reserve_left=0.0_R4
        reserve_right=0.0_R4
      endif

      if(aspect.and.want_aspect.and.linlog.eq.1) then
   10   test=(yhi-ylo)/(xhi-xlo)*
     1   (xwidth*
     .    (1.0_R4-(reserve_left+reserve_right))/real(nx_page,R4))/
     2   (yheight*
     .    (1.0_R4-(reserve_bottom+reserve_top))/real(ny_page,R4))
        if(b2plot_debug.ge.9) write(*,*) test
        if(noax.gt.0) then
          if(test.gt.1.0_R4) then
            test1=test*real(nox,R4)/real(nox+1,R4)
            if(abs(1.0_R4-test).gt.abs(1.0_R4-test1)) then
              if(mod(ix,2).eq.0) then
                xhi=xlo+(xhi-xlo)*real(nox+1,R4)/real(nox,R4)
              else
                xlo=xhi-(xhi-xlo)*real(nox+1,R4)/real(nox,R4)
              endif
              ix=ix+1
              nox=nox+1
              goto 10
            endif
          else
            test1=test*real(noy+1,R4)/real(noy,R4)
            if(abs(1.0_R4-test).gt.abs(1.0_R4-test1)) then
              if(mod(iy,2).eq.0) then
                yhi=ylo+(yhi-ylo)*real(noy+1,R4)/real(noy,R4)
              else
                ylo=yhi-(yhi-ylo)*real(noy+1,R4)/real(noy,R4)
              endif
              iy=iy+1
              noy=noy+1
              goto 10
            endif
          endif
        endif
        if(b2plot_debug.ge.9) write(*,*) test
      endif

      if(new_page) then
        call reset
        new_page=.false.
        if(.not.streql(globalheader,' ')) then
          if(hcol.gt.1) call gsplci(hcol)
          CALL SET(0.0_R4*XWIDTH,1.0_R4*XWIDTH,
     1     0.0_R4*YHEIGHT,1.0_R4*YHEIGHT,
     2     0.0_R4,1.0_R4,0.0_R4,1.0_R4,1)
          if(noax.gt.0)
     1     CALL PLCHHQ(0.5_R4,1.0_R4-0.01_R4*glsz,
     1     trim(GLOBALHEADER),
     1     1.0_R4/real(max(50,len_trim(globalheader)),R4)*
     .      xwidth/1.25_R4*glsz,
     1     0.0_R4,0.0_R4)
          call gsplci(1)
        endif
      endif
      dirty_page=.true.
      if (fcol.gt.1) then
        call gaseti('CAX',fcol)
        call gaseti('CLB',fcol)
        call gaseti('CMJ',fcol)
        call gaseti('CMN',fcol)
      end if
      if (bold) then
        call gasetr('WAX',3.)
        call gasetr('WLB',3.)
        call gasetr('WMJ',3.)
        call gasetr('WMN',3.)
        call gstxfp(-30,0)
        call pcseti('FN',30)
      end if

      CALL SET(0.0_R4*XSIZE+XPOS,1.0_R4*XSIZE+XPOS,
     1 0.0_R4*YSIZE+YPOS,1.0_R4*YSIZE+YPOS,
     2 0.0_R4,1.0_R4,0.0_R4,1.0_R4,1)
      dirty_graph=.true.

      if(b2plot_debug.ge.6) write(*,*) 'title ',trim(title)
      if(noax.gt.0) then
        if(gcol.gt.1) call gsplci(gcol)
        CALL PLCHHQ(0.5_R4,(1.0_R4-txtwdth*2*n_page),
     1   trim(TITLE),
     1   1.0_R4/100.0_R4/N_PAGE*lbsz,0.0_R4,0.0_R4)
        call gsplci(1)
        if(bcol.gt.1) call gsplci(bcol)
        if(wall)
     &   CALL PLCHHQ(wallbar_right,1.0_R4-txtwdth*2*n_page,
     1   trim(WALLBARLABEL),
     1   1.0_R4/100.0_R4/N_PAGE*lbsz,0.0_R4,1.0_R4)
        if(target.and.physical)
     &   CALL PLCHHQ(targetbar_right,1.0_R4-txtwdth*2*n_page,
     1   trim(TARGETBARLABEL),
     1   1.0_R4/100.0_R4/N_PAGE*lbsz,0.0_R4,1.0_R4)
        call gsplci(1)
      endif
      if(.not.streql(extralabel,' ')) then
        if(xcol.gt.1) call gsplci(xcol)
        if(noax.gt.0) then
          CALL PLCHHQ(0.5_R4,
     1     1.0_R4-txtwdth*(5-(offset-(nx_page/ny_page)/5))*n_page,
     1     trim(EXTRALABEL),
     1     1.0_R4/100.0_R4/N_PAGE*lbsz,0.0_R4,0.0_R4)
        endif
        call gsplci(1)
      endif
      if(.not.streql(xlabel,' ')) then
        if(xcol.gt.1) call gsplci(xcol)
        if(noax.gt.0) then
         CALL PLCHHQ(
     1        reserve_left+(1.0_R4-reserve_right-reserve_left)/2.0_R4,
     1        txtwdth*1*n_page,
     1        trim(XLABEL),
     1        1.0_R4/100.0_R4/N_PAGE*lbsz,0.0_R4,0.0_R4)
        endif
        call gsplci(1)
      endif
      if(.not.streql(ylabel,' ')) then
        if(xcol.gt.1) call gsplci(xcol)
        if(noax.gt.0) then
         CALL PLCHHQ(txtwdth*1*n_page,
     1        reserve_bottom+(1.0_R4-reserve_top-reserve_bottom)/2.0_R4,
     1        trim(YLABEL),
     1        1.0_R4/100.0_R4/N_PAGE*lbsz,90.0_R4,0.0_R4)
        endif
        call gsplci(1)
      endif

      if(b2plot_debug.ge.6) write(*,*) 'SET: ',
     1   reserve_left*XSIZE+XPOS,
     2   (1.0_R4-reserve_right)*XSIZE+XPOS,
     3   reserve_bottom*YSIZE+YPOS,
     4   (1.0_R4-reserve_top)*YSIZE+YPOS
      CALL SET(reserve_left*XSIZE+XPOS,
     1   (1.0_R4-reserve_right)*XSIZE+XPOS,
     1   reserve_bottom*YSIZE+YPOS,
     2   (1.0_R4-reserve_top)*YSIZE+YPOS,
     2   XLO,XHI,YLO,YHI,LINLOG)
      if(noax.gt.0) then
        CALL GASETI('LTY',1) ! ====> use PLCHHQ for drawing labels
        CALL LABMOD(FMTX,FMTY,NFMTX,NFMTY,
     1   max(nint(12*lbsz/N_PAGE),4),max(nint(12*lbsz/N_PAGE),4),
     2   nint(10*lbsz/n_page),nint(20*lbsz/n_page),0)
        if(b2plot_debug.ge.8) then
          write(*,*) 'nox,noxm,xlo ',nox,noxm,xlo
          write(*,*) 'noy,noym,ylo ',noy,noym,ylo
          write(*,*) 'noax ',noax
        endif
        CALL GRIDAL(NOX,NOXM,NOY,NOYM,1*noax,1*noax,5,XLO,YLO)
      endif
      if(lzero) then
        call lastd
        if (zcol.gt.1) call gsplci(zcol)
        if (zwidth.ne.1.0) call gslwsc(zwidth)
        call gsln(3)
        call dashdc('$$$$''''''''',8,0)
        if(b2plot_debug.ge.9) write(*,*) 'zwidth = ',zwidth
        if(xlo.lt.0.0_R4.and.xhi.gt.0.0_R4) then
          call lined(0.0_R4,ylo,0.0_R4,yhi)
        endif
        if(ylo.lt.0.0_R4.and.yhi.gt.0.0_R4) then
          call lined(xlo,0.0_R4,xhi,0.0_R4)
        endif
        if (label_same_lines.and.same_page.ne.0) then
          call dashdc(linefmt,8,lfsz)
        else
          call dashdc('$',1,0)
        endif
        call gsln(1)
        call lastd
        if (zcol.gt.1) call gsplci(1)
        if (zwidth.ne.1.0) call gslwsc(1.0)
      endif
      CALL CPSETI('MAP - MAPPING FLAG',3)
      call gaseti('CAX',1)
      call gaseti('CLB',1)
      call gaseti('CMJ',1)
      call gaseti('CMN',1)
      call gasetr('WAX',1.)
      call gasetr('WLB',1.)
      call gasetr('WMJ',1.)
      call gasetr('WMN',1.)
      call gstxfp(1,0)
!      call pcseti('FN',0)
#endif
      return
      end subroutine init_page
