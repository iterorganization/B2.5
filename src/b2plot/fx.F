c
      subroutine fx(iy)
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_b2plot_values
      implicit none
      character*32 tmpstring1,tmpstring2
      real (kind=R8) ::
     *   lmin,lmax,dsep,xmin,xmax,xmult
      integer ix,iy,is,lnblnk,ixstart,ixend

c     the line format common
      character*10 linefmt
      common /line_format/ linefmt

      intrinsic min,max,mod,sqrt
      external xertst,lnblnk

      ixstart=-1
      do while (region(ixstart,iy,0).eq.0 .and. ixstart.le.nx)
        ixstart=ixstart+1
      enddo
      call xertst(ixstart.le.nx, 'faulty parameter ixstart')
      ixend=nx
      do while (region(ixend,iy,0).eq.0 .and. ixend.gt.-1)
        ixend=ixend-1
      enddo
      call xertst(ixstart.le.ixend, 'faulty parameter ixend')
      do is=0,ndim(3,iwrk)-1
        tmpstring1=textwk(is,iwrk)
        lmin=work(ixstart,iy,is,iwrk)
        lmax=work(ixstart,iy,is,iwrk)
        xmult=1.0_R8
        if(xparallel) xmult=bb(ixstart,iy,3)/
     1   sqrt(bb(ixstart,iy,0)**2+bb(ixstart,iy,1)**2)
        wrk(ixstart)=0.0_R8
        do ix=ixstart+1,ixend
          wrk(ix)=wrk(ix-1)+hx(ix-1,iy)/2.0_R8*xmult
          if(xparallel) xmult=bb(ix,iy,3)/
     1     sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)
          wrk(ix)=wrk(ix)+hx(ix,iy)/2.0_R8*xmult
          lmin=min(lmin,work(ix,iy,is,iwrk))
          lmax=max(lmax,work(ix,iy,is,iwrk))
        enddo
        dsep=(1.0_R8-values(13))*wrk(ixstart)+(values(13))*wrk(ixend)
        do ix=ixstart,ixend
          wrk(ix)=wrk(ix)-dsep
        enddo
        write(tmpstring2,'('' <iy = '',i4,''>'')') iy
        if(values(9) .ne. values(10)) then
          lmin=values(9)
          lmax=values(10)
        endif
        if(lmin.eq.lmax) then
          lmin=min(lmin*0.9_R8,lmin*1.1_R8)
          lmax=max(lmax*0.9_R8,lmax*1.1_R8)
          if(lmin.eq.lmax) then
            lmin=-1.0_R8
            lmax=+1.0_R8
          end if
        end if
        if(mod(linlog,2).eq.0.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
        if(mod(linlog,2).eq.0.and.lmax.le.0.0_R8) lmax=1.0_R8
        if(mod(typeofplot,2) .eq. 1) then
          if(values(11) .ne. values(12)) then
            xmin=values(11)
            xmax=values(12)
          else
            xmin=wrk(ixstart)
            xmax=wrk(ixend)
          endif
          if(xmin.eq.xmax) then
            xmin=min(xmin*0.9_R8,xmin*1.1_R8)
            xmax=max(xmax*0.9_R8,xmax*1.1_R8)
            if(xmin.eq.xmax) then
              xmin=-1.0_R8
              xmax=+1.0_R8
            end if
          end if
          if(linlog.gt.2.and.xmin.le.0.0_R8) xmin=1.0e-30_R8
          if(linlog.gt.2.and.xmax.le.0.0_R8) xmax=1.0_R8
          call init_page(
     1     graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//
     2     tmpstring2(1:lnblnk(tmpstring2))//
     3     tmpstring1(1:lnblnk(tmpstring1)),
     1     xmin,xmax,lmin,lmax,
     4     linlog,.false.,.false.,.true.)
          if(writing_to_file) then
            write(96,'(''# '',a)')                          !srv  12.05.08
     1       graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) //
     2       tmpstring2(1:lnblnk(tmpstring2)) //
     2       tmpstring1(1:lnblnk(tmpstring1))
            do ix=ixstart,ixend
              write(96,*) wrk(ix),work(ix,iy,is,iwrk)        !srv  12.05.08
            enddo
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
          if (dcol.gt.1) call gsplci(dcol)
#endif
          do ix=ixstart,ixend-1
            call lined_d(
     1       wrk(ix),
     2       min(lmax,max(lmin,work(ix,iy,is,iwrk))),
     3       wrk(ix+1),
     4       min(lmax,max(lmin,work(ix+1,iy,is,iwrk))))
            if(rightix(ix,iy).ne.ix+1) then
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. dcol.gt.1) call gsplci(scol)
              call gslwsc_d(sepwidth)
              call dashdc('$',1,0)
#endif
              if ((wrk(ix)+wrk(ix+1)).gt.0.0_R8.or.linlog.le.2)
     >          call lined_d((wrk(ix)+wrk(ix+1))/2.0_R8,lmin,
     1                       (wrk(ix)+wrk(ix+1))/2.0_R8,lmax)
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. dcol.gt.1) call gsplci(dcol)
              call gslwsc_d(lwidth)
              call dashdc(linefmt,8,0)
#endif
            endif
          enddo
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (dcol.gt.1) call gsplci(1)
#endif
          call end_page
        endif
        if(typeofplot/2 .eq. 1) then
          do ix=-1,nx
            wrk(ix)=xx(ix,2)
            lmin=min(lmin,work(ix,iy,is,iwrk))
            lmax=max(lmax,work(ix,iy,is,iwrk))
          enddo
          if(values(9) .ne. values(10)) then
            lmin=values(9)
            lmax=values(10)
          endif
          if(lmin.eq.lmax) then
            lmin=min(lmin*0.9_R8,lmin*1.1_R8)
            lmax=max(lmax*0.9_R8,lmax*1.1_R8)
            if(lmin.eq.lmax) then
              lmin=-1.0_R8
              lmax=+1.0_R8
            end if
          end if
          if(mod(linlog,2).eq.0.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
          if(mod(linlog,2).eq.0.and.lmax.le.0.0_R8) lmax=1.0_R8
          if(values(1).ne.values(2)) then
            xmin=values(1)
            xmax=values(2)
          else
            xmin=wrk(ixstart)
            xmax=wrk(ixend)
          endif
          if(xmin.eq.xmax) then
            xmin=min(xmin*0.9_R8,xmin*1.1_R8)
            xmax=max(xmax*0.9_R8,xmax*1.1_R8)
            if(xmin.eq.xmax) then
              xmin=-1.0_R8
              xmax=+1.0_R8
            end if
          end if
          if(linlog.gt.2.and.xmin.le.0.0_R8) xmin=1.0e-30_R8
          if(linlog.gt.2.and.xmax.le.0.0_R8) xmax=1.0_R8
          call init_page(graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//
     2     tmpstring2(1:lnblnk(tmpstring2))//
     3     tmpstring1(1:lnblnk(tmpstring1)),
     1     xmin,xmax,lmin,lmax,
     4     linlog,.false.,.false.,.false.)
          if(writing_to_file) then
            write(96,'(''# '',a)')                           !srv  12.05.08
     1       graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) //
     2       tmpstring2(1:lnblnk(tmpstring2)) //
     2       tmpstring1(1:lnblnk(tmpstring1))
            do ix=-1,nx
              write(96,*) wrk(ix),work(ix,iy,is,iwrk)        !srv  12.05.08
            enddo
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
          if (dcol.gt.1) call gsplci(dcol)
#endif
          do ix=-1,nx-1
            call lined_d(
     1       min(values(2),max(values(1),wrk(ix))),
     2       min(lmax,max(lmin,work(ix,iy,is,iwrk))),
     3       min(values(2),max(values(1),wrk(ix+1))),
     4       min(lmax,max(lmin,work(ix+1,iy,is,iwrk))))
            if(rightix(ix,iy).ne.ix+1) then
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. dcol.gt.1) call gsplci(scol)
              call gslwsc_d(sepwidth)
              call dashdc('$',1,0)
#endif
              if ((wrk(ix)+wrk(ix+1)).gt.0.0_R8.or.linlog.le.2)
     >          call lined_d((wrk(ix)+wrk(ix+1))/2.0_R8,lmin,
     1                       (wrk(ix)+wrk(ix+1))/2.0_R8,lmax)
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. dcol.gt.1) call gsplci(dcol)
              call gslwsc_d(lwidth)
              call dashdc(linefmt,8,0)
#endif
            endif
          enddo
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (dcol.gt.1) call gsplci(1)
#endif
          call end_page
        endif
      enddo
      return
      end
