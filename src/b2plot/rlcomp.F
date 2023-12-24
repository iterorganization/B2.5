      subroutine rlcomp
      use b2mod_types
      use b2mod_constants
      use b2mod_diag
      use b2mod_subsys
      use b2mod_dimensions
      implicit none
c
c  version : 18.02.98 19:00
c
c======================================================================
c*** Computation of the radiation loads.
c*** Returns the total radiation power coming to each of the specified
c*** wall segments with possible shadowing taken into account.
c*** The radiation loads are calculated using a cylindrical approximation.
c***
c*** Input:
c***  xysrc   R and Z coordinates of the point sources
c***  ssrcm   their intensities
c***  xp,yp   coordinates of the boundary polygon P
c***  nbndry  number of points in shadowing structure S
c***  lbnd    indices of the S segments (assumed closed) in P
c***  ntrg    number of segments representing the radiation target T
c***  ltrg    their indices in P (must belong to S!)
c***
c*** Output
c***  radld   radiation power coming to the segments of T
c***  ttsrc   the total power radiated from the sources
c***  tttrg   the total power coming to T
c======================================================================
#include "rdload.inc"
      real (kind=R8) :: upi, udpi, ueps
      parameter (upi=pi)
      parameter (udpi=2.0e0_R8*upi, ueps=1.0e-20_R8)
      integer i, j, k
      real (kind=R8) :: x1, x2, y1, y2, angle, anglb, uda, utx, uxx
      intrinsic mod, atan2
      external xerrab
c<<<
      integer lun_s, lun_t, lun_l
      parameter (lun_s=97, lun_t=98, lun_l=99)
      real (kind=R8) :: uu_xe, uu_ye
      logical uu_dbg
      character tab
c>>>
c======================================================================
c
c<<<
       call subini ('rlcomp')
       tab=char(9)
#ifdef DBG
       uu_dbg=.true.
#else
       uu_dbg=.false.
#endif
c>>>
c
c*** Loop over the sources
c
c<<<
      if (uu_dbg) then
        print *,'rlcomp: nsrc,ntrg = ',nsrc,ntrg
        open(lun_s,file='sources.dat')
        write(lun_s,'(2a)') '# Radiation sources : ',trim(label)
        write(lun_s,'(5a)') '# R',tab,'Z',tab,'S'
      end if
c>>>
      ttsrc=0.e0_R8
      radld=0.e0_R8
      do i=1,nsrc
        ssrc=ssrcm(i)
        ttsrc=ttsrc+ssrc
        xsrc=xysrc(1,i)
        ysrc=xysrc(2,i)
c<<<
c        print *
c        print *,'Source number, ssrc, xsrc, ysrc'
c        print '(1p,2x,i4,3e10.2)',i, ssrc, xsrc, ysrc
        if (uu_dbg)
     >   write(lun_s,*) xsrc,tab,ysrc,tab,ssrc
c>>>
c
c*** Loop over the target segments
c
c<<<
        if(uu_dbg .and. i.eq.1) then
          open(lun_t,file='targ-rlc.dat')
          write(lun_t,'(2a)') '# Target geometry : ',trim(label)
          write(lun_t,
     ,     '(a1,4x,a2,12x,a2,12x,a2,12x,a2,11x,'//
     ,      'a3,11x,a3,11x,a3,10x,a5,9x,a5,9x,a5)')
     ,      '#','Rb','Zb','Re','Ze',
     ,      'ddx','ddy','uda','anglb','angle','radld'
          open(lun_l,file='layout-t.dat')
          write(lun_l,'(2a)') '# Target layout : ',trim(label)
          write(lun_l,'(3a)') '# R',tab,'Z'
        end if
c>>>
        do j=1,ntrg
          jb=ltrg(j)
          do k=1,nbndry
            if(lbnd(k) .eq. jb) then
              jcb=k
              go to 10
            end if
          end do
          write(*,*) 'Error in rlcomp: target not in structure',j,jb
          call xerrab ('Target not in structure')
c
   10     continue

          jcb=mod(jcb,nbndry)+1
          jc=lbnd(jcb)
          xb=xp(jb)
          yb=yp(jb)
          xc=xp(jc)
          yc=yp(jc)
          ddx=xc-xb
          ddy=yc-yb
          jm=1
          jsc=1
          sb=0.0e0_R8
          se=1.0e0_R8
c<<<
          uu_xe=xc
          uu_ye=yc
c          print *
c          print *,'Target number',j
c          print '(2x,4a4,6(3x,a4,4x))',' jb ',' jc ',' jm ',' jsc',
c     ,                      ' xb ',' yb ',' xc ',' yc ',' ddx ',' ddy '
c          print '(1p,1x,4i4,6e11.3)', jb, jc, jm, jsc,
c     ,                                        xb, yb, xc, yc, ddx, ddy
c>>>
c
c*** Calculate the shadowing...
c
          call rlshdw
c<<<
c      print *
c      print *,'After rlshdw.  j, se, sb:',j,se,sb
c>>>
          if(se .gt. sb) then
c
c*** At least something is seen. Calculate the load
c
            if(sb .le. 0.0e0_R8) then
              anglb=atan2(yb-ysrc,xb-xsrc)
            else
              x1=xb+ddx*sb
              y1=yb+ddy*sb
              anglb=atan2(y1-ysrc,x1-xsrc)
            end if
            x2=xb+ddx*se
            y2=yb+ddy*se
            angle=atan2(y2-ysrc,x2-xsrc)
            uda=angle-anglb
            if(abs(ddy) .gt. ueps) then
              utx=(ysrc-yb)/ddy
              uxx=xb+utx*ddx-xsrc
              if(uxx .lt. 0.0e0_R8) then
                if(utx .ge. 0.0e0_R8 .and.  utx .le. 1.0e0_R8) then
                  if(angle .lt. 0.0e0_R8) uda=uda+udpi
                  if(anglb .lt. 0.0e0_R8) uda=uda-udpi
                end if
              end if
            end if
            radld(j)=radld(j)+abs(uda)/udpi*ssrc
c<<<
c            print *
c            print *,'End of do 100. anglb, angle, uda, radld(j)'
c            print *, anglb, angle, uda, radld(j)
c            print *, 'xb, yb, ddx, ddy, uxx, utx'
c            print *, xb, yb, ddx, ddy, uxx, utx
c>>>
          else
            uda=0.0e0_R8
            anglb=0.0e0_R8
            angle=0.0e0_R8
          end if
c<<<
          if(uu_dbg) then
            write(lun_t,'(1p,10(e12.5,2x))')
     ,         xb,yb,uu_xe,uu_ye,ddx,ddy,uda/udpi,anglb,angle,radld(j)
            if (i.eq.1) then
              write(lun_l,*) xb,tab,yb
              write(lun_l,*) uu_xe,tab,uu_ye
              write(lun_l,*)
            end if
          end if
c>>>
        end do
      end do
c
c*** Calculate the total power coming to the target segments
c
      tttrg=0.0e0_R8
      do j=1,ntrg
        tttrg=tttrg+radld(j)
      end do
c======================================================================
c<<<
      if (uu_dbg) then
        close(lun_s)
        close(lun_t)
        close(lun_l)
      end if
c>>>
      call subend ()
      return
      end subroutine rlcomp
