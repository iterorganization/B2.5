
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
      real range,xmid,ymid,ceil,txtwdth,
     1 reserve_bottom,reserve_top,reserve_left,reserve_right,
     2 test,test1
      integer ix,iy,n_page,itxtwdth,lnblnk,
     1 nox,noxm,nfmtx,nnox,
     2 noy,noym,nfmty,nnoy
      external lnblnk

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
      xlo=min(1.0e30_R8,max(-1.0e30_R8,xlod))
      xhi=min(1.0e30_R8,max(-1.0e30_R8,xhid))
      ylo=min(1.0e30_R8,max(-1.0e30_R8,ylod))
      yhi=min(1.0e30_R8,max(-1.0e30_R8,yhid))
      if(xlo.eq.xhi) then
        xlo=min(xlo*0.9,xlo*1.1)
        xhi=max(xhi*0.9,xhi*1.1)
        if(xlo.eq.xhi) then
          xlo=-1.0
          xhi=+1.0
        endif
      endif
      if(ylo.eq.yhi) then
        ylo=min(ylo*0.9,ylo*1.1)
        yhi=max(yhi*0.9,yhi*1.1)
        if(ylo.eq.yhi) then
          ylo=-1.0
          yhi=+1.0
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
        xmid=(xhi+xlo)/2.0
        xlo=xmid-range/2.0/nx_page*xwidth
        xhi=xmid+range/2.0/nx_page*xwidth
        ymid=(yhi+ylo)/2.0
        ylo=ymid-range/2.0/ny_page*yheight
        yhi=ymid+range/2.0/ny_page*yheight
      endif
      xsize=xwidth/real(nx_page)
      ysize=yheight/real(ny_page)
      if(globalheader.ne.' ') then
        ysize=ysize*(1.00-0.02*lbsz)
      endif
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
          if(xhi.le.0.0) then
            write(*,*) 'xhi <= 0 and logarithmic axes chosen'
            xhi=1.0
          endif
          xlo=max(xlo,10.0**(-ndec)*xhi)
          xlo=10.0**floor(log10(xlo))
          xhi=10.0**ceil(log10(xhi))
          nox=1
          fmtx='(1pg8.2)'
          nfmtx=8
          noxm=1
        endif
        if(lbsz.gt.1.0.and.nox.gt.5) then
          nnox=nint(ceil(nox/2.0))
          noxm=noxm*2
          xhi=xlo+(xhi-xlo)*real(nnox)*2.0/real(nox)
          nox=nnox
          ix=ix+1
        endif

        if(linlog.eq.1.or.linlog.eq.3) then
          logy=.false.
          call nice(ylo,yhi,noy,noym,fmty,nfmty,nice_axes)
        else
          logy=.true.
          if(yhi.le.0.0) then
            write(*,*) 'yhi <= 0 and logarithmic axes chosen'
            yhi=1.0
          endif
          yhi=max(yhi,-ylo)
          ylo=max(ylo,10.0**(-ndec)*yhi)
          ylo=10.0**floor(log10(ylo))
          yhi=10.0**ceil(log10(yhi))
          noy=1
          fmty='(1pg8.2)'
          nfmty=8
          noym=1
        endif
        if(lbsz.gt.1.0.and.noy.gt.5) then
          nnoy=nint(ceil(noy/2.0))
          noym=noym*2
          yhi=ylo+(yhi-ylo)*real(nnoy)*2.0/real(noy)
          noy=nnoy
          iy=iy+1
        endif
      endif
      if(b2plot_debug.ge.8) write(*,*) xlo,xhi,ylo,yhi

      n_page=max(nx_page,ny_page)
      itxtwdth=max(nint(12*lbsz/N_PAGE),4)
      txtwdth=real(itxtwdth)/real(1024)

      if(noax.gt.0) then
        reserve_bottom=txtwdth*5/ysize
        reserve_top=txtwdth*6/ysize
        reserve_left=txtwdth*(nfmty+3)/xsize
        reserve_right=0.0
        if(target.and.physical) then
          targetbar_right=1.0-reserve_right
          reserve_right=reserve_right+0.08*0.2+0.08*0.8*lbsz
          targetbar_left=1.0-reserve_right
          reserve_right=reserve_right+0.02
        endif
        if(wall) then
          wallbar_right=1.0-reserve_right
          reserve_right=reserve_right+0.08*0.2+0.08*0.8*lbsz
          wallbar_left=1.0-reserve_right
          reserve_right=reserve_right+0.02
        endif
        if(colorbar) then
          colorbar_right=1.0-reserve_right
          reserve_right=reserve_right+0.08*0.2+0.08*0.8*lbsz
          colorbar_left=1.0-reserve_right
          reserve_right=reserve_right+0.02
        endif
        if(extralabel.ne.' ') reserve_top=reserve_top+txtwdth*2/ysize
        if(reserve_right.eq.0.0) reserve_right=0.05*lbsz
      else
        reserve_bottom=0.0
        reserve_top=0.0
        reserve_left=0.0
        reserve_right=0.0
      endif

      if(aspect.and.want_aspect.and.linlog.eq.1) then
   10   test=(yhi-ylo)/(xhi-xlo)*
     1   (xwidth*(1.0-(reserve_left+reserve_right))/real(nx_page))/
     2   (yheight*(1.0-(reserve_bottom+reserve_top))/real(ny_page))
        if(b2plot_debug.ge.9) write(*,*) test
        if(noax.gt.0) then
          if(test.gt.1.0) then
            test1=test*real(nox)/real(nox+1)
            if(abs(1.0-test).gt.abs(1.0-test1)) then
              if(mod(ix,2).eq.0) then
                xhi=xlo+(xhi-xlo)*real(nox+1)/real(nox)
              else
                xlo=xhi-(xhi-xlo)*real(nox+1)/real(nox)
              endif
              ix=ix+1
              nox=nox+1
              goto 10
            endif
          else
            test1=test*real(noy+1)/real(noy)
            if(abs(1.0-test).gt.abs(1.0-test1)) then
              if(mod(iy,2).eq.0) then
                yhi=ylo+(yhi-ylo)*real(noy+1)/real(noy)
              else
                ylo=yhi-(yhi-ylo)*real(noy+1)/real(noy)
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
        if(globalheader.ne.' ') then
          if(hcol.gt.1) call gsplci(hcol)
          CALL SET(0.0*XWIDTH,1.0*XWIDTH,
     1     0.0*YHEIGHT,1.0*YHEIGHT,
     2     0.0,1.0,0.0,1.0,1)
          if(noax.gt.0)
     1     CALL PLCHHQ(0.5,1.0-0.01*lbsz,
     1     GLOBALHEADER(1:LNBLNK(GLOBALHEADER)),
     1     1.0/real(max(50,lnblnk(globalheader)))*xwidth/1.25,0.0,0.0)
          call gsplci(1)
        endif
      endif
      dirty_page=.true.

      CALL SET(0.0*XSIZE+XPOS,1.0*XSIZE+XPOS,
     1 0.0*YSIZE+YPOS,1.0*YSIZE+YPOS,
     2 0.0,1.0,0.0,1.0,1)
      dirty_graph=.true.

      if(b2plot_debug.ge.6) write(*,*) 'title ',title(1:lnblnk(title))
      if(noax.gt.0) then
        if(gcol.gt.1) call gsplci(gcol)
        CALL PLCHHQ(0.5,(1.0-txtwdth*2*n_page),
     1   TITLE(1:LNBLNK(TITLE)),
     1   1.0/100.0/N_PAGE*lbsz,0.0,0.0)
        call gsplci(1)
        if(bcol.gt.1) call gsplci(bcol)
        if(wall)
     &   CALL PLCHHQ(wallbar_right,(1.0-txtwdth*2*n_page),
     1   WALLBARLABEL(1:LNBLNK(WALLBARLABEL)),
     1   1.0/100.0/N_PAGE*lbsz,0.0,1.0)
        if(target.and.physical)
     &   CALL PLCHHQ(targetbar_right,(1.0-txtwdth*2*n_page),
     1   TARGETBARLABEL(1:LNBLNK(TARGETBARLABEL)),
     1   1.0/100.0/N_PAGE*lbsz,0.0,1.0)
        call gsplci(1)
      endif
      IF(EXTRALABEL.NE.' ') THEN
        if(xcol.gt.1) call gsplci(xcol)
        if(noax.gt.0)
     1   CALL PLCHHQ(0.5,(1.0-txtwdth*4*n_page),
     1   EXTRALABEL(1:LNBLNK(EXTRALABEL)),
     1   1.0/100.0/N_PAGE*lbsz,0.0,0.0)
        call gsplci(1)
      ENDIF

      if(b2plot_debug.ge.6) write(*,*) 'SET: ',
     1   reserve_left*XSIZE+XPOS,
     2   (1.0-reserve_right)*XSIZE+XPOS,
     3   reserve_bottom*YSIZE+YPOS,
     4   (1.0-reserve_top)*YSIZE+YPOS
      CALL SET(reserve_left*XSIZE+XPOS,
     1   (1.0-reserve_right)*XSIZE+XPOS,
     1   reserve_bottom*YSIZE+YPOS,
     2   (1.0-reserve_top)*YSIZE+YPOS,
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
        if(xlo.lt.0.0.and.xhi.gt.0.0) then
          call lined(0.0,ylo,0.0,yhi)
        endif
        if(ylo.lt.0.0.and.yhi.gt.0.0) then
          call lined(xlo,0.0,xhi,0.0)
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
#endif
      return
      end subroutine init_page
