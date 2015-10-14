c
      subroutine fy(ix)
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2mod_boundary_namelist
      implicit none
      character*44 tmpstring1,tmpstring2
      character*256 tmpstring,titlestring
      real (kind=R8) ::
     *   lmin,lmax,x1,x2,y1,y2,dsep,dist,xmin,xmax,
     1   y_transform
      integer ix,ixl,ixr,iy,is,iiso,lnblnk,ysepl,ysepr,iystart,iyend
      integer ixx
      logical iniso
      intrinsic min,max,sqrt
      external xertst, lnblnk
      dist(x1,y1,x2,y2)=sqrt((x1-x2)**2+(y1-y2)**2)
c     the line format common
      character*10 linefmt
      common /line_format/ linefmt

c      if(nidim .ne. 2) call line(2)
      iystart=-1
      do while (region(ix,iystart,0).eq.0 .and. iystart.le.ny)
        iystart=iystart+1
      enddo
      call xertst(iystart.le.ny, 'faulty parameter iystart')
      iyend=ny
      do while (region(ix,iyend,0).eq.0 .and. iyend.gt.-1)
        iyend=iyend-1
      enddo
      call xertst(iystart.le.iyend, 'faulty parameter iyend')
      if (b2plot_debug.ge.6) 
     > write(6,*) ' iystart, iyend : ', iystart, iyend
      if (.not.map_to_midplane) then
        ixx = ix
      else
        if (b2plot_debug.ge.5) 
     >   write(*,*) 'Mapping to midplane : jxa = ',jxa
        ixx = jxa
      endif
      do is=0,ndim(3,iwrk)-1
        tmpstring1=textwk(is,iwrk)
        lmin=work(ix,iystart,is,iwrk)
        lmax=work(ix,iystart,is,iwrk)
        if(y_transform_l .ne. 0) then
          do iy=iystart,iyend
            wrk(iy)=y_transform(cr(ixx,iy),cz(ixx,iy))
            write(*,*) wrk(iy),cz(ixx,iy)
            lmin=min(lmin,work(ix,iy,is,iwrk))
            lmax=max(lmax,work(ix,iy,is,iwrk))
          enddo
        else
          wrk(iystart)=dist(cr(ixx,iystart),cz(ixx,iystart),
     1          0.5_R8*(crx(ixx,iystart,0)+crx(ixx,iystart,1)),     
     2          0.5_R8*(cry(ixx,iystart,0)+cry(ixx,iystart,1)))
          do iy=iystart+1,iyend
            wrk(iy)=wrk(iy-1)+dist(cr(ixx,iy-1),cz(ixx,iy-1),
     1                             cr(ixx,iy),cz(ixx,iy))
            lmin=min(lmin,work(ix,iy,is,iwrk))
            lmax=max(lmax,work(ix,iy,is,iwrk))
          enddo
          if (iyref.eq.-2) then
            ixl=ixx
            ysepl=iystart
            iy=iystart
            iniso=.false.
            do while (ixl.gt.-1 .and. .not.iniso .and.
     &                leftix(ixl,iy).eq.ixl-1)
              ixl=ixl-1
              iniso=.false.
              do iiso = 1,nniso
                iniso = iniso .or.
     &           (ixl.ge.nxiso1(iiso).and.ixl.le.nxiso2(iiso).and.
     &            iy .ge.nyiso1(iiso).and.iy .le.nyiso2(iiso))
              end do
            end do
            if(.not.iniso .and. ixl.ne.-1) then
              do while (ysepl.lt.iyend .and.
     &                  leftix(ixl,ysepl).ne.ixl-1)
                ysepl=ysepl+1
              end do
            end if
            if(ysepl.eq.iystart) ysepl=iyend
            ixr=ixx
            ysepr=iystart
            iniso=.false.
            do while (ixr.lt.nx .and. .not.iniso .and.
     &                rightix(ixr,iy).eq.ixr+1)
              ixr=ixr+1
              iniso=.false.
              do iiso = 1,nniso
                iniso = iniso .or.
     &           (ixr.ge.nxiso1(iiso).and.ixr.le.nxiso2(iiso).and.
     &            iy .ge.nyiso1(iiso).and.iy .le.nyiso2(iiso))
              end do
            end do
            if(.not.iniso .and. ixr.ne.nx) then
              do while(ysepr.lt.iyend .and.
     &                 rightix(ixr,ysepr).ne.ixr+1)
                ysepr=ysepr+1
              end do
            end if
            if(ysepr.eq.iystart) ysepr=iyend
            if((ysepl.ne.iyend .or. ysepr.ne.iyend) .and.
     &         (iystart.le.min(ysepr,ysepl).and.
     &          iyend  .ge.min(ysepr,ysepl))) then
!            dsep=0.5_R8*(wrk(min(ysepr,ysepl))+wrk(min(ysepr,ysepl)-1))
              dsep=wrk(min(ysepr,ysepl))-
     1         dist(
     1          cr(ixx,min(ysepr,ysepl)),cz(ixx,min(ysepr,ysepl)),
     1            0.5_R8*(crx(ixx,min(ysepr,ysepl),0)+
     1                    crx(ixx,min(ysepr,ysepl),1)),     
     2            0.5_R8*(cry(ixx,min(ysepr,ysepl),0)+
     2                    cry(ixx,min(ysepr,ysepl),1)))
              if (b2plot_debug.ge.6) 
     >          write(6,*) ' ysepr, ysepl : ', ysepr, ysepl
            else
              dsep=0.0_R8
            endif
          else
            if (iystart.le.iyref .and. iyend.ge.iyref) then
              dsep=wrk(iyref)-
     1         dist(cr(ixx,iyref),cz(ixx,iyref),
     1            0.5_R8*(crx(ixx,iyref,0)+crx(ixx,iyref,1)),     
     2            0.5_R8*(cry(ixx,iyref,0)+cry(ixx,iyref,1)))
            else
              dsep=0.0_R8
            endif
          endif
          do iy=iystart,iyend
            wrk(iy)=wrk(iy)-dsep
          enddo
        endif
        if (.not.map_to_midplane) then
          write(tmpstring2,'('' <ix = '',i4,''>'')') ix
        else
          write(tmpstring2,
     .     '('' <ix = '',i4,''> (mapped to midplane)'')') ix
        endif
        tmpstring = 
     1     graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//
     2     tmpstring2(1:lnblnk(tmpstring2))//
     3     tmpstring1(1:lnblnk(tmpstring1))
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
            xmin=min(wrk(iystart),wrk(iyend))
          endif
          if(values(12).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) then
            xmax=values(12)
          else
            xmax=max(wrk(iystart),wrk(iyend))
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
            write(nwtf,'(''# '',a)')
     1       graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) //
     2       tmpstring2(1:lnblnk(tmpstring2)) //
     2       tmpstring1(1:lnblnk(tmpstring1))
            do iy=iystart,iyend
              write(nwtf,*) wrk(iy),work(ix,iy,is,iwrk)
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
          endif
#endif
          do iy=iystart,iyend-1
            call lined_d(
     1       wrk(iy),
     2       min(lmax,max(lmin,work(ix,iy,is,iwrk))),
     3       wrk(iy+1),
     4       min(lmax,max(lmin,work(ix,iy+1,is,iwrk))))
            if(topiy(ix,iy).ne.iy+1) then
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. dcol.ne.1) call gsplci(scol)
              call gslwsc_d(sepwidth)
              call dashdc('$',1,0)
#endif
              if ((wrk(iy)+wrk(iy+1)).gt.0.0_R8.or.linlog.le.2)
     >          call lined_d((wrk(iy)+wrk(iy+1))/2.0_R8,lmin,
     1                       (wrk(iy)+wrk(iy+1))/2.0_R8,lmax)
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
          do iy=iystart,iyend
            wrk(iy)=yy(iy,2)
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
          if(values(11).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) then
            xmin=values(11)
          else
            xmin=min(wrk(iystart),wrk(iyend))
          endif
          if(values(12).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) then
            xmax=values(12)
          else
            xmax=max(wrk(iystart),wrk(iyend))
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
          write(tmpstring2,'('' <ix = '',i4,''>'')') ix
          call init_page(titlestring,
     1     xmin,xmax,lmin,lmax,
     4     linlog,.false.,.false.,.false.)
          if(writing_to_file) then
            write(nwtf,'(''# '',a)')
     1       graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) //
     2       tmpstring2(1:lnblnk(tmpstring2)) //
     2       tmpstring1(1:lnblnk(tmpstring1))
            do iy=iystart,iyend
              write(nwtf,*) wrk(iy),work(ix,iy,is,iwrk)
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
          endif
#endif
          do iy=iystart,iyend-1
            call lined_d(
     1       wrk(iy),
     2       min(lmax,max(lmin,work(ix,iy,is,iwrk))),
     3       wrk(iy+1),
     4       min(lmax,max(lmin,work(ix,iy+1,is,iwrk))))
            if(topiy(ix,iy).ne.iy+1) then
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. dcol.ne.1) call gsplci(scol)
              call gslwsc_d(sepwidth)
              call dashdc('$',1,0)
#endif
              if ((wrk(iy)+wrk(iy+1)).gt.0.0_R8.or.linlog.le.2)
     >          call lined_d((wrk(iy)+wrk(iy+1))/2.0_R8,lmin,
     1                       (wrk(iy)+wrk(iy+1))/2.0_R8,lmax)
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
      end subroutine fy
c
      function y_transform(r,z)
      use b2mod_types
      real (kind=R8) ::
     *   y_transform
      real (kind=R8) ::
     *   r,z
      real (kind=R8) ::
     *   z1(100),s1(100),z2(100),s2(100)
      logical first,open_file
      integer n1,n2,i1,i2
      real (kind=R8) ::
     *   x, r0r
      save z1,s1,z2,s2,n1,n2,first
      data first /.true./

      write(*,*) 'y_transform'
      if(first) then
        if(.not.open_file(99,'n.s.r.z.a_inner',
     1   'data.local|data')) then
          stop 'Couldn''t open n.s.r.z.a_inner'
        endif
        n1=0
        read(99,*)
   10   read(99,*,end=20) x,s1(n1+1),x,z1(n1+1)
        n1=n1+1
        goto 10
   20   continue
        write(*,*) 'Read n.s.r.z.a_inner',n1
        close(99)
        if(.not.open_file(99,'n.s.r.z.a_outer',
     1   'data.local|data')) then
          stop 'Couldn''t open n.s.r.z.a_outer'
        endif
        n2=0
        read(99,*)
   30   read(99,*,end=40) x,s2(n2+1),x,z2(n2+1)
        n2=n2+1
        goto 30
   40   continue
        write(*,*) 'Read n.s.r.z.a_outer',n2
        close(99)
      endif
      first=.false.

c     if(r0r .eq. 0.0_R8) r0r=1.5_R8
      r0r=1.5_R8
      if(r .lt. r0r) then
        write(*,*) 'r <'
        i1=(z-z1(1))/(z1(n1)-z1(1))*n1+1
        i1=min(n1-1,max(1,i1))
        do while(z1(i1) .lt. z .and. i1 .gt. 1)
          i1=i1-1
          write(*,*) '-: i1',i1
        enddo
        do while(z1(i1+1) .gt. z .and. i1 .lt. n1-1)
          i1=i1+1
          write(*,*) '+: i1',i1
        enddo
        y_transform=(z-z1(i1))/(z1(i1+1)-z1(i1))*
     1   (s1(i1+1)-s1(i1))+s1(i1)
      else
        write(*,*) 'r >'
        i2=(z-z2(1))/(z2(n2)-z2(1))*n2+1
        i2=min(n2-1,max(1,i2))
        do while(z2(i2) .gt. z .and. i2 .gt. 1)
          i2=i2-1
          write(*,*) '-: i1',i1
        enddo
        do while(z2(i2+1) .lt. z .and. i2 .lt. n2-1)
          i2=i2+1
          write(*,*) '+: i1',i1
        enddo
        y_transform=(z-z2(i2))/(z2(i2+1)-z2(i2))*
     1   (s2(i2+1)-s2(i2))+s2(i2)
      endif

      return
      end function y_transform

      subroutine plot_expfile
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      implicit none
      character*256 line, title
      real (kind=R8) ::
     *   xmin,xmax,ymin,ymax
      logical open_file, is_comment
      external is_comment
      integer lngarr,lnblnk,i,inc
      parameter (lngarr=100 000)
      real (kind=R8) ::
     *   xvalues(lngarr),yvalues(lngarr)

      write(*,*) 'Called plot_expfile'

      if(.not.open_file(99,expfile,
     1 'data.local/exp|data/exp')) then
        write(*,*) 'Couldn''t open ',trim(expfile)
      else
        title = ''
   10   read(99,'(a)',end=40) line
        if(is_comment(line)) then
          title = line(1:lnblnk(line))
          write(*,*) line(1:lnblnk(line))
          goto 10
        endif
        backspace(99)
        inc=0
   20   read(99,*,end=30) xvalues(inc+1),yvalues(inc+1)
        inc=inc+1
        goto 20
   30   continue
        close(99)
        xmin=xvalues(1)
        xmax=xvalues(1)
        ymin=yvalues(1)
        ymax=yvalues(1)
        do i=2,inc
          xmin=min(xmin,xvalues(i))
          xmax=max(xmax,xvalues(i))
          ymin=min(ymin,yvalues(i))
          ymax=max(ymax,yvalues(i))
        enddo
        if(values(11).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) 
     >   xmin=values(11)
        if(values(12).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) 
     >   xmax=values(12)
        if(xmin.eq.xmax) then
          xmin=min(xmin*0.9_R8,xmin*1.1_R8)
          xmax=max(xmax*0.9_R8,xmax*1.1_R8)
          if(xmin.eq.xmax) then
            xmin=-1.0_R8
            xmax=+1.0_R8
          end if
        end if
        if(values(14) .ne. values(15)) then
          ymin=values(14)
          ymax=values(15)
        endif
        if(ymin.eq.ymax) then
          ymin=min(ymin*0.9_R8,ymin*1.1_R8)
          ymax=max(ymax*0.9_R8,ymax*1.1_R8)
          if(ymin.eq.ymax) then
            ymin=-1.0_R8
            ymax=+1.0_R8
          end if
        end if
        if(mod(linlog,2).eq.0.and.ymin.le.0.0_R8) ymin=1.0e-30_R8
        if(mod(linlog,2).eq.0.and.ymax.le.0.0_R8) ymax=1.0_R8
        call init_page(
     1   title(1:lnblnk(title)),
     1   xmin,xmax,ymin,ymax,linlog,.false.,.false.,.true.)
         if(b2plot_debug.ge.4)
     >   write(*,*) 'xmin,xmax,ymin,ymax ',xmin,xmax,ymin,ymax
        if(writing_to_file) then
          write(nwtf,'(''# '',a)')
     1     title(1:lnblnk(title))
          do i=1,inc
            write(nwtf,*) xvalues(i),yvalues(i)
          enddo
        endif
#ifndef PLPLOT
        call lastd
        if(ecol .gt. 1) then
          if (expline) then
            call gsplci(ecol)
          else
            call gspmci(ecol)
          endif
        elseif(ecol .lt. 0) then
          if (expline) then
            call gsplci(-ecol)
          else
            call gspmci(-ecol)
          endif
          ecol=ecol-1
        endif
        if (ewidth.ne.1.0_R8) then
          call lastd
          if (expline) then
            call gslwsc_d(ewidth)
          else
            call gsmksc_d(ewidth)
          endif
        endif
#endif
        if (expline) then
          do i=1,inc-1
            call lined_d(xvalues(i),yvalues(i),
     .                   xvalues(i+1),yvalues(i+1))
          enddo
        else
          do i=1,inc
            call grpts_d(xvalues(i),yvalues(i),expsymbol)
          enddo
        endif
#ifndef PLPLOT
        if (ewidth.ne.1.0_R8) then
          call lastd
          if (expline) then
            call gslwsc_d(1.0_R8)
          else
            call gsmksc_d(1.0_R8)
          endif
        endif
#endif
        call end_page
   40   continue
        if(ecol .ne. 1) then
#ifndef PLPLOT
          call lastd
          if (expline) then
            call gsplci(1)
          else
            call gspmci(1)
          endif
#endif
        endif
      endif
      return
      end subroutine plot_expfile

      subroutine plot_wallfile
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      implicit none
      character*256 line, title
      real (kind=R8) ::
     *   xmin,xmax,ymin,ymax
      logical open_file, is_comment
      external is_comment
      integer lngarr,lnblnk,i,inc
      parameter (lngarr=100 000)
      real (kind=R8) ::
     *   xvalues(lngarr),yvalues(lngarr)

      write(*,*) 'Called plot_wallfile'

      if(.not.open_file(99,wallfile,
     1 'data.local/exp|data/exp')) then
        write(*,*) 'Couldn''t open ',trim(wallfile)
      else
        title = ''
        if(ecol .gt. 1) then
#ifndef PLPLOT
          call lastd
          call gsplci(ecol)
#endif
        elseif(ecol.lt.0) then
#ifndef PLPLOT
          call lastd
          call gsplci(-ecol)
#endif
          ecol=ecol-1
        endif
   10   read(99,'(a)',end=40) line
        if(is_comment(line)) then
          title = line(1:lnblnk(line))
          write(*,*) line(1:lnblnk(line))
          goto 10
        endif
        backspace(99)
        inc=0
   20   read(99,*,end=30) xvalues(inc+1),yvalues(inc+1)
        inc=inc+1
        goto 20
   30   continue
        close(99)
        xmin=xvalues(1)
        xmax=xvalues(1)
        ymin=yvalues(1)
        ymax=yvalues(1)
        do i=2,inc
          xmin=min(xmin,xvalues(i))
          xmax=max(xmax,xvalues(i))
          ymin=min(ymin,yvalues(i))
          ymax=max(ymax,yvalues(i))
        enddo
        if(values(11).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) 
     >   xmin=values(11)
        if(values(12).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) 
     >   xmax=values(12)
        if(xmin.eq.xmax) then
          xmin=min(xmin*0.9_R8,xmin*1.1_R8)
          xmax=max(xmax*0.9_R8,xmax*1.1_R8)
          if(xmin.eq.xmax) then
            xmin=-1.0_R8
            xmax=+1.0_R8
          end if
        end if
        if(values(14) .ne. values(15)) then
          ymin=values(14)
          ymax=values(15)
        endif
        if(ymin.eq.ymax) then
          ymin=min(ymin*0.9_R8,ymin*1.1_R8)
          ymax=max(ymax*0.9_R8,ymax*1.1_R8)
          if(ymin.eq.ymax) then
            ymin=-1.0_R8
            ymax=+1.0_R8
          end if
        end if
        if(mod(wlinlog,2).eq.0.and.ymin.le.0.0_R8) ymin=1.0e-30_R8
        if(mod(wlinlog,2).eq.0.and.ymax.le.0.0_R8) ymax=1.0_R8
        call init_page(
     1   title(1:lnblnk(title)),
     1   xmin,xmax,ymin,ymax,wlinlog,.false.,.false.,.true.)
        if(b2plot_debug.ge.4)
     >   write(*,*) 'xmin,xmax,ymin,ymax ',xmin,xmax,ymin,ymax
        if(writing_to_file) then
          write(nwtf,'(''# '',a)')
     1     title(1:lnblnk(title))
          do i=1,inc
            write(nwtf,*) xvalues(i),yvalues(i)
          enddo
        endif
#ifndef PLPLOT
        if (wwidth.ne.1.0_R8) then
          call lastd
          call gslwsc_d(wwidth)
        endif
#endif
        do i=1,inc-1
          call lined_d(xvalues(i),yvalues(i),xvalues(i+1),yvalues(i+1))
        enddo
#ifndef PLPLOT
        if (wwidth.ne.1.0_R8) then
          call lastd
          call gslwsc_d(1.0_R8)
        endif
#endif
        call end_page
   40   continue
        if(ecol .ne. 1) then
#ifndef PLPLOT
          call lastd
          call gsplci(1)
#endif
        endif
      endif
      return
      end subroutine plot_wallfile
