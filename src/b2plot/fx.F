c
      subroutine fx(iy)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      use b2us_data
      implicit none
      integer, intent(in) :: iy
      character*12 hlp_frm
      character*32 tmpstring1,tmpstring2
      character*256 tmpstring,titlestring
      real (kind=R8) ::
     *   lmin,lmax,dsep,xmin,xmax,xmult
      real (kind=R8) :: wrk(-1:mpg%nx)
      integer ix,is,iCv,iCn,ixstart,ixend,len_of_digits

c     the line format common
      character*10 linefmt
      common /line_format/ linefmt

      intrinsic min,max,mod,sqrt
      external xertst,len_of_digits

      if(mpg%isClassicalGrid.eq.0) return
      ixstart=-1
      iCv=mpg%imapCv(ixstart,iy)
      do while (iCv.eq.0 .and. ixstart.le.mpg%nx)
        ixstart=ixstart+1
        iCv=mpg%imapCv(ixstart,iy)
      enddo
      call xertst(ixstart.le.mpg%nx, 'faulty parameter ixstart')
      ixend=mpg%nx
      iCv=mpg%imapCv(ixend,iy)
      do while (iCv.eq.0 .and. ixend.gt.-1)
        ixend=ixend-1
        iCv=mpg%imapCv(ixend,iy)
      enddo
      call xertst(ixstart.le.ixend, 'faulty parameter ixend')
      write(hlp_frm,'(a,i1,a)') '(a,i',len_of_digits(iy),',a)'
      do is=0,ndim(2,iwrk)-1
        iCv=mpg%imapCv(ixstart,iy)
        tmpstring1=textwk(is,iwrk)
        lmin=work(iCv,is,iwrk)
        lmax=work(iCv,is,iwrk)
        xmult=1.0_R8
        if(xparallel) xmult=geo%cvBb(iCv,3)/
     1   sqrt(geo%cvBb(iCv,0)**2+geo%cvBb(iCv,1)**2)
        wrk(ixstart)=0.0_R8
        do ix=ixstart+1,ixend
          wrk(ix)=wrk(ix-1)+geo%cvHx(iCv)/2.0_R8*xmult
          iCv=mpg%imapCv(ix,iy)
          if(xparallel) xmult=geo%cvBb(iCv,3)/
     1     sqrt(geo%cvBb(iCv,0)**2+geo%cvBb(iCv,1)**2)
          wrk(ix)=wrk(ix)+geo%cvHx(iCv)/2.0_R8*xmult
          lmin=min(lmin,work(iCv,is,iwrk))
          lmax=max(lmax,work(iCv,is,iwrk))
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
              iCv=mpg%imapCv(ix,iy)
              write(nwtf,*) wrk(ix),work(iCv,is,iwrk)
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
            iCv=mpg%imapCv(ix,iy)
            iCn=mpg%imapCv(ix+1,iy)
            call lined_d(
     1       wrk(ix),
     2       min(lmax,max(lmin,work(iCv,is,iwrk))),
     3       wrk(ix+1),
     4       min(lmax,max(lmin,work(iCn,is,iwrk))))
            if(mpg%cvOnClosedSurface(iCv).neqv.
     .         mpg%cvOnClosedSurface(iCn)) then
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
        if(typeofplot/2 .eq. 1 .and. mpg%isClassicalGrid .eq. 1) then
          do ix=ixstart,ixend
            iCv=mpg%imapCv(ix,iy)
            wrk(ix)=xx(iCv,2)
            lmin=min(lmin,work(iCv,is,iwrk))
            lmax=max(lmax,work(iCv,is,iwrk))
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
              iCv=mpg%imapCv(ix,iy)
              write(nwtf,*) wrk(ix),work(iCv,is,iwrk)
            enddo
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
          if (dcol.gt.1) then
            call gsplci(dcol)
          elseif(dcol.lt.0) then
            if(typeofplot.eq.3.and.mpg%isClassicalGrid.eq.1)
     >        dcol=dcol+1
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
            iCv=mpg%imapCv(ix,iy)
            iCn=mpg%imapCv(ix+1,iy)
            call lined_d(
     1       wrk(ix),
     2       min(lmax,max(lmin,work(iCv,is,iwrk))),
     3       wrk(ix+1),
     4       min(lmax,max(lmin,work(iCn,is,iwrk))))
            if(mpg%cvOnClosedSurface(iCv).neqv.
     .         mpg%cvOnClosedSurface(iCn)) then
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
c
      subroutine ft(iFt)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      use b2us_data
      implicit none
      integer, intent(in) :: iFt
      character*12 hlp_frm
      character*32 tmpstring1,tmpstring2
      character*256 tmpstring,titlestring
      real (kind=R8) :: dsep,lmin,lmax,xmin,xmax,xmult
      real (kind=R8), allocatable :: wrk(:)
      integer i,is,iCv,iCn
      integer len_of_digits
      external len_of_digits

c     the line format common
      character*10 linefmt
      common /line_format/ linefmt

      write(hlp_frm,'(a,i1,a)') '(a,i',len_of_digits(iFt),',a)'
      allocate(wrk(1:mpg%ftCvP(iFt,2)))
      do is=0,ndim(2,iwrk)-1
        iCv=mpg%ftCv(mpg%ftCvP(iFt,1))
        tmpstring1=textwk(is,iwrk)
        lmin=work(iCv,is,iwrk)
        lmax=work(iCv,is,iwrk)
        xmult=1.0_R8
        if(xparallel) xmult=geo%cvBb(iCv,3)/
     1   sqrt(geo%cvBb(iCv,0)**2+geo%cvBb(iCv,1)**2)
        wrk(1)=0.0_R8
        do i=2,mpg%ftCvP(iFt,2)
          wrk(i)=wrk(i-1)+geo%cvHx(iCv)/2.0_R8*xmult
          iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
          if(xparallel) xmult=geo%cvBb(iCv,3)/
     1     sqrt(geo%cvBb(iCv,0)**2+geo%cvBb(iCv,1)**2)
          wrk(i)=wrk(i)+geo%cvHx(iCv)/2.0_R8*xmult
          lmin=min(lmin,work(iCv,is,iwrk))
          lmax=max(lmax,work(iCv,is,iwrk))
        end do
        if (values(13).ne.init_value) then
          dsep=(1.0_R8-values(13))*wrk(1)+
     .                (values(13))*wrk(mpg%ftCvP(iFt,2))
        else ! old default behaviour : values(13)=0.5
          dsep=0.5_R8*(wrk(1)+wrk(mpg%ftCvP(iFt,2)))
        end if
        do i=1,mpg%ftCvP(iFt,2)
          wrk(i)=wrk(i)-dsep
        enddo
        write(tmpstring2,hlp_frm) ' <iFt = ',iFt,'>'
        tmpstring =
     1   trim(graphlabel(iwrk))//trim(tmpstring2)//trim(tmpstring1)
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
            xmin=min(wrk(1),wrk(mpg%ftCvP(iFt,2)))
          endif
          if(values(12).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) then
            xmax=values(12)
          else
            xmax=max(wrk(1),wrk(mpg%ftCvP(iFt,2)))
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
            do i=1,mpg%ftCvP(iFt,2)
              iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
              write(nwtf,*) wrk(i),work(iCv,is,iwrk)
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
          do i=1,mpg%ftCvP(iFt,2)-1
            iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
            iCn=mpg%ftCv(mpg%ftCvP(iFt,1)+i)
            call lined_d(
     1       wrk(i),
     2       min(lmax,max(lmin,work(iCv,is,iwrk))),
     3       wrk(i+1),
     4       min(lmax,max(lmin,work(iCn,is,iwrk))))
          end do
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (dcol.ne.1) call gsplci(1)
#endif
          call end_page
        end if
        if(typeofplot/2 .eq. 1) then
          do i=1,mpg%ftCvP(iFt,2)
            iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
            wrk(i)=real(i,R8)
            lmin=min(lmin,work(iCv,is,iwrk))
            lmax=max(lmax,work(iCv,is,iwrk))
          end do
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
            xmin=wrk(1)
            xmax=wrk(mpg%ftCvP(iFt,2))
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
            do i=1,mpg%ftCvP(iFt,2)
              iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
              write(nwtf,*) wrk(i),work(iCv,is,iwrk)
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
          do i=1,mpg%ftCvP(iFt,2)-1
            if(wrk(i+1).lt.xmin) cycle
            if(wrk(i).gt.xmax) cycle
            iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
            iCn=mpg%ftCv(mpg%ftCvP(iFt,1)+i)
            call lined_d(
     1       wrk(i),
     2       min(lmax,max(lmin,work(iCv,is,iwrk))),
     3       wrk(i+1),
     4       min(lmax,max(lmin,work(iCn,is,iwrk))))
          enddo
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (dcol.ne.1) call gsplci(1)
#endif
          call end_page
        end if
      end do
      deallocate(wrk)

      return
      end subroutine ft
