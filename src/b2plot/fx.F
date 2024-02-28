c
      subroutine fx(iy)
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_b2plot_values
      implicit none
      integer, intent(in) :: iy
      character*12 hlp_frm
      character*32 tmpstring1,tmpstring2
      character*256 tmpstring,titlestring
      real (kind=R8) ::
     *   lmin,lmax,dsep,xmin,xmax,xmult
      integer ix,is,ixstart,ixend,len_of_digits

c     the line format common
      character*10 linefmt
      common /line_format/ linefmt

      intrinsic min,max,mod,sqrt
      external xertst,len_of_digits

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
      write(hlp_frm,'(a,i1,a)') '(a,i',len_of_digits(iy),',a)'
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
        if (values(13).ne.init_value) then
          dsep=(1.0_R8-values(13))*wrk(ixstart)+(values(13))*wrk(ixend)
        else if (ixref.ge.ixstart .and. ixref.le.ixend) then
cxpb : new default behaviour : zero of distance set at ixref
          dsep=wrk(ixref)
        else ! old default behaviour : values(13)=0.5
          dsep=0.5_R8*(wrk(ixstart)+wrk(ixend))
        endif
        do ix=ixstart,ixend
          wrk(ix)=wrk(ix)-dsep
        enddo
        write(tmpstring2,hlp_frm) ' <iy = ',iy,'>'
        tmpstring =
     1    trim(graphlabel(iwrk))//trim(tmpstring2)//trim(tmpstring1)
        call customise_label(tmpstring,titlestring)
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
          if(values(11).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) then
            xmin=values(11)
          else
            xmin=min(wrk(ixstart),wrk(ixend))
          endif
          if(values(12).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) then
            xmax=values(12)
          else
            xmax=max(wrk(ixstart),wrk(ixend))
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
          call init_page(titlestring,
     1     xmin,xmax,lmin,lmax,
     4     linlog,.false.,.false.,.true.)
          if(writing_to_file) then
            write(nwtf,'(''# '',a)') trim(titlestring)
            do ix=ixstart,ixend
              write(nwtf,*) wrk(ix),work(ix,iy,is,iwrk)
            enddo
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
          if (dcol.gt.1) then
            call gsplci(dcol)
          elseif(dcol.lt.0) then
            call gsplci(-dcol)
            dcol=dcol-1
            if (dcol.lt.-17) then
              dcol=-mod(abs(dcol),17)
              if (dcol.eq.0) dcol=-17
            endif
          endif
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
              if (scol.gt.1 .or. dcol.ne.1) call gsplci(scol)
              call gslwsc_d(sepwidth)
              call dashdc('$',1,0)
#endif
              if ((wrk(ix)+wrk(ix+1)).gt.0.0_R8.or.linlog.le.2)
     >          call lined_d((wrk(ix)+wrk(ix+1))/2.0_R8,lmin,
     1                       (wrk(ix)+wrk(ix+1))/2.0_R8,lmax)
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. dcol.ne.1) call gsplci(abs(dcol))
              call gslwsc_d(lwidth)
              call dashdc(linefmt,8,lfsz)
#endif
            endif
          enddo
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (dcol.ne.1) call gsplci(1)
#endif
          call end_page
        endif
        if(typeofplot/2 .eq. 1) then
          do ix=ixstart,ixend
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
          call init_page(titlestring,
     1     xmin,xmax,lmin,lmax,
     4     linlog,.false.,.false.,.false.)
          if(writing_to_file) then
            write(nwtf,'(''# '',a)') trim(titlestring)
            do ix=ixstart,ixend
              write(nwtf,*) wrk(ix),work(ix,iy,is,iwrk)
            enddo
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
          if (dcol.gt.1) then
            call gsplci(dcol)
          elseif(dcol.lt.0) then
            if(typeofplot.eq.3) dcol=dcol+1
            call gsplci(-dcol)
            dcol=dcol-1
            if (dcol.lt.-17) then
              dcol=-mod(abs(dcol),17)
              if (dcol.eq.0) dcol=-17
            endif
          endif
#endif
          do ix=ixstart,ixend-1
            if(wrk(ix+1).lt.xmin) cycle
            if(wrk(ix).gt.xmax) cycle
            call lined_d(
     1       wrk(ix),
     2       min(lmax,max(lmin,work(ix,iy,is,iwrk))),
     3       wrk(ix+1),
     4       min(lmax,max(lmin,work(ix+1,iy,is,iwrk))))
            if(rightix(ix,iy).ne.ix+1) then
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. dcol.ne.1) call gsplci(scol)
              call gslwsc_d(sepwidth)
              call dashdc('$',1,0)
#endif
              if ((wrk(ix)+wrk(ix+1)).gt.0.0_R8.or.linlog.le.2)
     >          call lined_d((wrk(ix)+wrk(ix+1))/2.0_R8,lmin,
     1                       (wrk(ix)+wrk(ix+1))/2.0_R8,lmax)
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. dcol.ne.1) call gsplci(abs(dcol))
              call gslwsc_d(lwidth)
              call dashdc(linefmt,8,lfsz)
#endif
            endif
          enddo
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (dcol.ne.1) call gsplci(1)
#endif
          call end_page
        endif
      enddo

      return
      end subroutine fx
