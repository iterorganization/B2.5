c
      subroutine plot_wall_elements

      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      implicit none
#include "ank.h"
c     the line format common
      character*10 linefmt
      common /line_format/ linefmt

      character*256 extralabel
      integer i, is, iseg
      real (kind=R8) ::
     2   xvalues(nwall),
     1   xmin,xmax,ymin,ymax,lmin,lmax
      integer lnblnk
      external lnblnk

      call xertst(nwall.ge.ndata, 'increase nwall')

      do is = 0, nwdim(2,iwwrk)-1
        if (typeofplot .ge. 2) then
          if(b2plot_debug.ge.7) write(*,*) 'plot along wall index'
          do i = 1, ndata
            xvalues(i) = i   
          end do
          xmin=xvalues(1)
          xmax=xvalues(1)
          ymin=wall_stack(1,is,iwwrk)
          ymax=wall_stack(1,is,iwwrk)
          do i=2,ndata
            xmin=min(xmin,xvalues(i))
            xmax=max(xmax,xvalues(i))
            ymin=min(ymin,wall_stack(i,is,iwwrk))
            ymax=max(ymax,wall_stack(i,is,iwwrk))
          enddo
          lmin=ymin
          lmax=ymax
          if(lmin.eq.lmax) then
            lmin=min(lmin*0.9,lmin*1.1)
            lmax=max(lmax*0.9,lmax*1.1)
            if(lmin.eq.lmax) then
              lmin=-1.0
              lmax=+1.0
            end if
          end if
          if(wlinlog.eq.2.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
          if(wlinlog.eq.2.and.lmax.le.0.0_R8) lmax=1.0_R8
          if(values(11) .ne. values(12)) then
            xmin=values(11)
            xmax=values(12)
          endif
          if(xmin.eq.xmax) then
            xmin=min(xmin*0.9,xmin*1.1)
            xmax=max(xmax*0.9,xmax*1.1)
            if(xmin.eq.xmax) then
              xmin=-1.0
              xmax=+1.0
            end if
          end if
          if(wlinlog.gt.2.and.xmin.le.0.0_R8) xmin=1.0e-30_R8
          if(wlinlog.gt.2.and.xmax.le.0.0_R8) xmax=1.0_R8
          if(values(14) .ne. values(15)) then
            ymin=values(14)
            ymax=values(15)
          endif
          call init_page(
     2     wall_label(iwwrk)(1:lnblnk(wall_label(iwwrk)))
     .      //' '//textwwk(is,iwwrk),
     1     xmin,xmax,ymin,ymax,
     2     wlinlog,.false.,.false.,.false.)
          if (b2plot_debug.ge.6) 
     >     write(*,*) 'xmin,xmax,ymin,ymax ',xmin,xmax,ymin,ymax
          if(writing_to_file) then
            write(96,'(''# '',a)')           !srv  12.05.08
     2       wall_label(iwwrk)(1:lnblnk(wall_label(iwwrk)))
     .        //' '//textwwk(is,iwwrk)
            call get_extralabel(extralabel)
            if(len_trim(extralabel).ne.0) then
              write(96,'(''# '',a)') trim(extralabel)        !srv  12.05.08
            endif
            do i=1,ndata
              write(96,*) xvalues(i),wall_stack(i,is,iwwrk)        !srv  12.05.08
            enddo
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
          if (wcol.gt.1) call gsplci(wcol)
#endif
          do iseg = 1, nsegments
            do i = wall_begin(iseg), wall_end(iseg)-1
              call lined_d(
     .         min(xmax,max(xmin,xvalues(i))),
     .         min(ymax,max(ymin,wall_stack(i,is,iwwrk))),
     .         min(xmax,max(xmin,xvalues(i+1))),
     .         min(ymax,max(ymin,wall_stack(i+1,is,iwwrk))))
            enddo
            if (iseg.lt.nsegments) then
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. wcol.gt.1) call gsplci(scol)
              call gslwsc_d(sepwidth)
              call dashdc('$',1,0)
#endif
              if (wlinlog.le.2.or.
     .(xvalues(wall_end(iseg))+xvalues(wall_begin(iseg+1))).gt.0.0_R8)
     >          call lined_d(
     .           (xvalues(wall_end(iseg))+xvalues(wall_begin(iseg+1)))
     .            /2.0_R8,lmin,
     .           (xvalues(wall_end(iseg))+xvalues(wall_begin(iseg+1)))
     .            /2.0_R8,lmax)
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. wcol.gt.1) call gsplci(wcol)
              call gslwsc_d(lwidth)
              call dashdc(linefmt,8,0)
#endif
            endif
          enddo
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (wcol.gt.1) call gsplci(1)
#endif
          call end_page
        endif
        if(mod(typeofplot,2) .eq. 1) then
          if(b2plot_debug.ge.7) write(*,*) 'plot along wall distance'
          do i = 1, ndata
            if (i.eq.1) then
              xvalues(i) = 0.5_R8*
     .         sqrt((r2(i)-r1(i))**2+(z2(i)-z1(i))**2)
            else
              xvalues(i) = xvalues(i-1) + 0.5_R8*(
     .         sqrt((r2(i-1)-r1(i-1))**2+(z2(i-1)-z1(i-1))**2)+
     .         sqrt((r2(i)-r1(i))**2+(z2(i)-z1(i))**2))+
     .         sqrt((r1(i)-r2(i-1))**2+(z1(i)-z2(i-1))**2)
            endif
          enddo
          xmin=xvalues(1)
          xmax=xvalues(1)
          ymin=wall_stack(1,is,iwwrk)
          ymax=wall_stack(1,is,iwwrk)
          do i=2,ndata
            xmin=min(xmin,xvalues(i))
            xmax=max(xmax,xvalues(i))
            ymin=min(ymin,wall_stack(i,is,iwwrk))
            ymax=max(ymax,wall_stack(i,is,iwwrk))
          enddo
          lmin=ymin
          lmax=ymax
          if(lmin.eq.lmax) then
            lmin=min(lmin*0.9,lmin*1.1)
            lmax=max(lmax*0.9,lmax*1.1)
            if(lmin.eq.lmax) then
              lmin=-1.0
              lmax=+1.0
            end if
          end if
          if(wlinlog.eq.2.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
          if(wlinlog.eq.2.and.lmax.le.0.0_R8) lmax=1.0_R8
          if(values(11) .ne. values(12)) then
            xmin=values(11)
            xmax=values(12)
          endif
          if(xmin.eq.xmax) then
            xmin=min(xmin*0.9,xmin*1.1)
            xmax=max(xmax*0.9,xmax*1.1)
            if(xmin.eq.xmax) then
              xmin=-1.0
              xmax=+1.0
            end if
          end if
          if(wlinlog.gt.2.and.xmin.le.0.0_R8) xmin=1.0e-30_R8
          if(wlinlog.gt.2.and.xmax.le.0.0_R8) xmax=1.0_R8
          if(values(14) .ne. values(15)) then
            ymin=values(14)
            ymax=values(15)
          endif
          call init_page(
     2     wall_label(iwwrk)(1:lnblnk(wall_label(iwwrk)))
     .      //' '//textwwk(is,iwwrk),
     1     xmin,xmax,ymin,ymax,
     2     wlinlog,.false.,.false.,.true.)
          if (b2plot_debug.ge.6) 
     >     write(*,*) 'xmin,xmax,ymin,ymax ',xmin,xmax,ymin,ymax
          if(writing_to_file) then
            write(96,'(''# '',a)')                           !srv  12.05.08
     2       wall_label(iwwrk)(1:lnblnk(wall_label(iwwrk)))
     .        //' '//textwwk(is,iwwrk)
            call get_extralabel(extralabel)
            if(len_trim(extralabel).ne.0) then
              write(96,'(''# '',a)') trim(extralabel)        !srv  12.05.08
            endif
            do i=1,ndata
              write(96,*) xvalues(i),wall_stack(i,is,iwwrk)  !srv  12.05.08
            enddo
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
          if (wcol.gt.1) call gsplci(wcol)
#endif
          do iseg = 1, nsegments
            do i = wall_begin(iseg), wall_end(iseg)-1
              call lined_d(
     .         min(xmax,max(xmin,xvalues(i))),
     .         min(ymax,max(ymin,wall_stack(i,is,iwwrk))),
     .         min(xmax,max(xmin,xvalues(i+1))),
     .         min(ymax,max(ymin,wall_stack(i+1,is,iwwrk))))
            enddo
            if (iseg.lt.nsegments) then
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. wcol.gt.1) call gsplci(scol)
              call gslwsc_d(sepwidth)
              call dashdc('$',1,0)
#endif
              if (wlinlog.le.2.or.
     .(xvalues(wall_end(iseg))+xvalues(wall_begin(iseg+1))).gt.0.0_R8)
     >          call lined_d(
     .          (xvalues(wall_end(iseg))+xvalues(wall_begin(iseg+1)))
     .          /2.0_R8,lmin,
     .          (xvalues(wall_end(iseg))+xvalues(wall_begin(iseg+1)))
     .          /2.0_R8,lmax)
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. wcol.gt.1) call gsplci(wcol)
              call gslwsc_d(lwidth)
              call dashdc(linefmt,8,0)
#endif
            endif
          enddo
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (wcol.gt.1) call gsplci(1)
#endif
          call end_page
        endif
      enddo
      close(96)                                     !srv 12.05.08
      return
      end subroutine plot_wall_elements
