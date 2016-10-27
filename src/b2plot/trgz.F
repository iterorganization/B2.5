c
      subroutine trgz(ibnd)
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer ibnd
      character*14 ibstring
      character*256 tmpstring, titlestring
      real (kind=R8) ::
     *   lmin,lmax,x1,x2,y1,y2,dist,xmin,xmax
      integer i,lnblnk
      external lnblnk

c     the line format common
      character*10 linefmt
      common /line_format/ linefmt

      intrinsic min,max,mod,sqrt
      external xertst
      dist(x1,y1,x2,y2)=sqrt((x1-x2)**2+(y1-y2)**2)
c      if(nidim .ne. 1) call line(1)
      if(ndepth_nml.le.1) then
        write(*,*) 'This requires a 2-D target array!'
        return
      endif
      wwrk(1)=0.0_R8     
      lmin=target_stack(ibnd,1,itwrk)
      lmax=target_stack(ibnd,1,itwrk)
      do i=2,ndepth_nml
        wwrk(i)=wwrk(i-1)+
     1   dist(targetr(ibnd,i-1),targetz(ibnd,i-1),
     2        targetr(ibnd,i  ),targetz(ibnd,i  ))
        lmin=min(lmin,target_stack(ibnd,i,itwrk))
        lmax=max(lmax,target_stack(ibnd,i,itwrk))
      enddo
      write(ibstring,'('' <ibnd = '',i4,''>'')') ibnd
      if(values(32) .ne. values(33)) then
        lmin=values(32)
        lmax=values(33)
      endif
      if(lmin.eq.lmax) then
        lmin=min(lmin*0.9_R8,lmin*1.1_R8)
        lmax=max(lmax*0.9_R8,lmax*1.1_R8)
        if(lmin.eq.lmax) then
          lmin=-1.0_R8
          lmax=+1.0_R8
        end if
      end if
      if(mod(trglinlog,2).eq.0.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
      if(mod(trglinlog,2).eq.0.and.lmax.le.0.0_R8) lmax=1.0_R8
      if(mod(typeofplot,2) .eq. 1) then
        if(values(34).ne.init_value .and.
     &     .not.(values(34).eq.0.0_R8 .and. values(35).eq.0.0_R8)) then
          xmin=values(34)
        else
          xmin=wwrk(1)
        endif
        if(values(35).ne.init_value .and.
     &     .not.(values(34).eq.0.0_R8 .and. values(35).eq.0.0_R8)) then
          xmax=values(35)
        else
          xmax=wwrk(ndepth_nml)
        endif
        if(xmin.eq.xmax) then
          xmin=min(xmin*0.9_R8,xmin*1.1_R8)
          xmax=max(xmax*0.9_R8,xmax*1.1_R8)
          if(xmin.eq.xmax) then
            xmin=-1.0_R8
            xmax=+1.0_R8
          end if
        end if
        if(trglinlog.gt.2.and.xmin.le.0.0_R8) xmin=1.0e-30_R8
        if(trglinlog.gt.2.and.xmax.le.0.0_R8) xmax=1.0_R8
        tmpstring=
     1     target_label(itwrk)(1:lnblnk(target_label(itwrk)))//
     2     ibstring
        call customise_label(tmpstring,titlestring)
        call init_page(titlestring,
     1     xmin,xmax,lmin,lmax,
     2     trglinlog,.false.,.false.,.true.)
        if(writing_to_file) then
          write(nwtf,'(''# '',a)') titlestring(1:lnblnk(titlestring))
          do i=1,ndepth_nml
            write(nwtf,*) wwrk(i),target_stack(ibnd,i,itwrk)
          enddo
        endif
#ifndef PLPLOT
        call lastd
        if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
        if (tcol.gt.1) then
          call gsplci(tcol)
        elseif (tcol.lt.0) then
          call gsplci(-tcol)
          tcol=tcol-1
        endif
#endif
        do i=1,ndepth_nml-1
          call lined_d(
     1       min(xmax,max(xmin,wwrk(i))),
     2       min(lmax,max(lmin,target_stack(ibnd,i  ,itwrk))),
     3       min(xmax,max(xmin,wwrk(i+1))),
     4       min(lmax,max(lmin,target_stack(ibnd,i+1,itwrk))))
        end do
#ifndef PLPLOT
        call lastd
        if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
        if (tcol.ne.1) call gsplci(1)
#endif
        call end_page
      endif
      if(typeofplot/2 .eq. 1) then
        do i=1,ndepth_nml
          wwrk(i)=targety(i,2)
          lmin=min(lmin,target_stack(ibnd,i,itwrk))
          lmax=max(lmax,target_stack(ibnd,i,itwrk))
        enddo
        if(values(32) .ne. values(33)) then
          lmin=values(32)
          lmax=values(33)
        endif
        if(lmin.eq.lmax) then
          lmin=min(lmin*0.9_R8,lmin*1.1_R8)
          lmax=max(lmax*0.9_R8,lmax*1.1_R8)
          if(lmin.eq.lmax) then
            lmin=-1.0_R8
            lmax=+1.0_R8
          end if
        end if
        if(mod(trglinlog,2).eq.0.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
        if(mod(trglinlog,2).eq.0.and.lmax.le.0.0_R8) lmax=1.0_R8
        if(values(26).eq.values(27)) then
          values(26)=min(values(26)*0.9_R8,values(26)*1.1_R8)
          values(27)=max(values(27)*0.9_R8,values(27)*1.1_R8)
          if(values(26).eq.values(27)) then
            values(26)=-1.0_R8
            values(27)=+1.0_R8
          end if
        end if
        if(mod(trglinlog,2).eq.0.and.values(26).le.0.0_R8) 
     >   values(26)=1.0e-30_R8
        if(mod(trglinlog,2).eq.0.and.values(27).le.0.0_R8) 
     >   values(27)=1.0_R8
        call init_page(titlestring,
     1     values(26),values(27),lmin,lmax,
     3     trglinlog,.false.,.false.,.false.)
        if(writing_to_file) then
          write(nwtf,'(''# '',a)') titlestring(1:lnblnk(titlestring))
          do i=1,ndepth_nml
            write(nwtf,*) wwrk(i),target_stack(ibnd,i,itwrk)
          enddo
        endif
#ifndef PLPLOT
        call lastd
        if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
        if (tcol.gt.1) then
          call gsplci(tcol)
        elseif (tcol.lt.0) then
          if(typeofplot.eq.3) tcol=tcol+1
          call gsplci(-tcol)
          tcol=tcol-1
        endif
#endif
        do i=1,ndepth_nml-1
          call lined_d(
     1       min(values(27),max(values(26),wwrk(i))),
     2       min(lmax,max(lmin,target_stack(ibnd,i  ,itwrk))),
     3       min(values(27),max(values(26),wwrk(i+1))),
     4       min(lmax,max(lmin,target_stack(ibnd,i+1,itwrk))))
        enddo
#ifndef PLPLOT
        call lastd
        if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
        if (tcol.ne.1) call gsplci(1)
#endif
        call end_page
      endif
      return
      end subroutine trgz
