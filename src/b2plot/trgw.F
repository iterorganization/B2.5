c
      subroutine trgw(i)
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer i
      character*32 tmpstring
      real (kind=R8) ::
     *   lmin,lmax,x1,x2,y1,y2,dist,xmin,xmax
      integer ibnd,iseg,lnblnk

c     the line format common
      character*10 linefmt
      common /line_format/ linefmt

      intrinsic min,max,mod,sqrt
      external xertst
      dist(x1,y1,x2,y2)=sqrt((x1-x2)**2+(y1-y2)**2)
c      if(nidim .ne. 1) call line(1)
      wwrk(1)=dist(targetr(1,i),targetz(1,i),
     1          targetcr(1,i,1),targetcz(1,i,1))     
      lmin=target_stack(1,i,itwrk)
      lmax=target_stack(1,i,itwrk)
      do ibnd=2,nbnd
        wwrk(ibnd)=wwrk(ibnd-1)+
     1   dist(targetr(ibnd-1,i),targetz(ibnd-1,i),
     2        targetcr(ibnd-1,i,2),targetcz(ibnd-1,i,2))+
     1   dist(targetcr(ibnd-1,i,2),targetcz(ibnd-1,i,2),
     2        targetcr(ibnd  ,i,1),targetcz(ibnd  ,i,1))+
     1   dist(targetr(ibnd,i),targetz(ibnd,i),
     2        targetcr(ibnd,i,1),targetcz(ibnd,i,1))
        lmin=min(lmin,target_stack(ibnd,i,itwrk))
        lmax=max(lmax,target_stack(ibnd,i,itwrk))
      enddo
      write(tmpstring,'('' <i = '',i4,''>'')') i
      if(values(32) .ne. values(33)) then
        lmin=values(32)
        lmax=values(33)
      endif
      if(lmin.eq.lmax) then
        lmin=min(lmin*0.9,lmin*1.1)
        lmax=max(lmax*0.9,lmax*1.1)
        if(lmin.eq.lmax) then
          lmin=-1.0
          lmax=+1.0
        end if
      end if
      if(mod(trglinlog,2).eq.0.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
      if(mod(trglinlog,2).eq.0.and.lmax.le.0.0_R8) lmax=1.0e-30_R8
      if(mod(typeofplot,2) .eq. 1) then
        if(values(34) .ne. values(35)) then
          xmin=values(34)
          xmax=values(35)
        else
          xmin=wwrk(1)
          xmax=wwrk(nbnd)
        endif
        if(xmin.eq.xmax) then
          xmin=min(xmin*0.9,xmin*1.1)
          xmax=max(xmax*0.9,xmax*1.1)
          if(xmin.eq.xmax) then
            xmin=-1.0
            xmax=+1.0
          end if
        end if
        if(trglinlog.gt.2.and.xmin.le.0.0_R8) xmin=1.0e-30_R8
        if(trglinlog.gt.2.and.xmax.le.0.0_R8) xmax=1.0_R8
        call init_page(
     1     target_label(itwrk)(1:lnblnk(target_label(itwrk)))//
     2     tmpstring(1:lnblnk(tmpstring)),
     1     values(34),values(35),values(32),values(33),
     2     trglinlog,.false.,.false.,.true.)
        if(writing_to_file) then
          write(96,'(''# '',a)')                                  !srv  12.05.08
     1	   target_label(itwrk)(1:lnblnk(target_label(itwrk)))//
     2     tmpstring(1:lnblnk(tmpstring))
          do ibnd=1,nbnd
            write(96,*) wwrk(ibnd),target_stack(ibnd,i,itwrk)     !srv  12.05.08
          enddo
        endif
#ifndef PLPLOT
        call lastd
        if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
        if (tcol.gt.1) call gsplci(tcol)
#endif
        do ibnd=1,nbnd-1
          call lined_d(
     1       min(xmax,max(xmin,wwrk(ibnd))),
     2       min(lmax,max(lmin,target_stack(ibnd  ,i,itwrk))),
     3       min(xmax,max(xmin,wwrk(ibnd+1))),
     4       min(lmax,max(lmin,target_stack(ibnd+1,i,itwrk))))
        end do
        do iseg = 1, nsegments-1
#ifndef PLPLOT
          call lastd
          if (scol.gt.1 .or. tcol.gt.1) call gsplci(scol)
          call gslwsc_d(sepwidth)
          call dashdc('$',1,0)
#endif
          if ((wwrk(wall_end(iseg))+wwrk(wall_begin(iseg+1))).gt.0.0_R8
     &      .or.trglinlog.le.2) call lined_d(
     .       (wwrk(wall_end(iseg))+wwrk(wall_begin(iseg+1)))/2.0_R8,
     1       lmin,
     1       (wwrk(wall_end(iseg))+wwrk(wall_begin(iseg+1)))/2.0_R8,
     2       lmax)
#ifndef PLPLOT
          call lastd
          if (scol.gt.1 .or. tcol.gt.1) call gsplci(tcol)
          call gslwsc_d(lwidth)
          call dashdc(linefmt,8,0)
#endif
        enddo
#ifndef PLPLOT
        call lastd
        if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
        if (tcol.gt.1) call gsplci(1)
#endif
        call end_page
      endif
      if(typeofplot/2 .eq. 1) then
        do ibnd=1,nbnd
          wwrk(ibnd)=targetx(ibnd,2)
          lmin=min(lmin,target_stack(ibnd,i,itwrk))
          lmax=max(lmax,target_stack(ibnd,i,itwrk))
        enddo
        if(values(32) .ne. values(33)) then
          lmin=values(32)
          lmax=values(33)
        endif
        if(lmin.eq.lmax) then
          lmin=min(lmin*0.9,lmin*1.1)
          lmax=max(lmax*0.9,lmax*1.1)
          if(lmin.eq.lmax) then
            lmin=-1.0
            lmax=+1.0
          end if
        end if
        if(mod(trglinlog,2).eq.0.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
        if(mod(trglinlog,2).eq.0.and.lmax.le.0.0_R8) lmax=1.0_R8
        if(values(24).eq.values(25)) then
          values(24)=min(values(24)*0.9,values(24)*1.1)
          values(25)=max(values(25)*0.9,values(25)*1.1)
          if(values(24).eq.values(25)) then
            values(24)=-1.0
            values(25)=+1.0
          end if
        end if
        if(mod(trglinlog,2).eq.0.and.values(24).le.0.0_R8) 
     >   values(24)=1.0e-30_R8
        if(mod(trglinlog,2).eq.0.and.values(25).le.0.0_R8) 
     >   values(25)=1.0_R8
        call init_page(
     1     target_label(itwrk)(1:lnblnk(target_label(itwrk)))//
     2     tmpstring(1:lnblnk(tmpstring)),
     1     values(24),values(25),values(32),values(33),
     3     trglinlog,.false.,.false.,.false.)
        if(writing_to_file) then
          write(96,'(''# '',a)')                                 !srv  12.05.08
     1	   target_label(itwrk)(1:lnblnk(target_label(itwrk)))//
     2     tmpstring(1:lnblnk(tmpstring))
          do ibnd=1,nbnd
            write(96,*) wwrk(ibnd),target_stack(ibnd,i,itwrk)    !srv  12.05.08
          enddo
        endif
#ifndef PLPLOT
        call lastd
        if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
        if (tcol.gt.1) call gsplci(tcol) 
#endif
        do ibnd=1,nbnd-1
          call lined_d(
     1       min(values(25),max(values(24),wwrk(ibnd))),
     2       min(lmax,max(lmin,target_stack(ibnd  ,i,itwrk))),
     3       min(values(25),max(values(24),wwrk(ibnd+1))),
     4       min(lmax,max(lmin,target_stack(ibnd+1,i,itwrk))))
        enddo
        do iseg = 1, nsegments-1
#ifndef PLPLOT
          call lastd
          if (scol.gt.1 .or. tcol.gt.1) call gsplci(scol)
          call gslwsc_d(sepwidth)
          call dashdc('$',1,0)
#endif
          if ((wwrk(wall_end(iseg))+wwrk(wall_begin(iseg+1))).gt.0.0_R8
     &      .or.trglinlog.le.2) call lined_d(
     .       (wwrk(wall_end(iseg))+wwrk(wall_begin(iseg+1)))/2.0_R8,
     1       lmin,
     1       (wwrk(wall_end(iseg))+wwrk(wall_begin(iseg+1)))/2.0_R8,
     2       lmax)
#ifndef PLPLOT
          call lastd
          if (scol.gt.1 .or. tcol.gt.1) call gsplci(tcol)
          call gslwsc_d(lwidth)
          call dashdc(linefmt,8,0)
#endif
        enddo
#ifndef PLPLOT
        call lastd
        if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
        if (tcol.gt.1) call gsplci(1)
#endif
        call end_page
      endif
      return
      end subroutine trgw
