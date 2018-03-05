      subroutine rlcomp
      use b2mod_types
      use b2mod_constants
      use b2mod_diag
      implicit none
c
c  version : 18.02.98 19:00
c
c======================================================================
c*** Computation of the radiation loads.
c*** Returns the total radiation power coming to each of the specified
c*** wall segments with possible shadowing taken into account.
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
#include "DIMENSIONS.F"
#include "rdload.inc"
      real (kind=R8) :: upi, udpi, ueps
      parameter (upi=pi)
      parameter (udpi=2.0e0_R8*upi, ueps=1.0e-20_R8)
      integer i, j, k
      real (kind=R8) :: x1, x2, y1, y2, angle, anglb, uda, utx, uxx
      intrinsic mod, atan2
      external subini, subend, xerrab
c<<<
c      integer lun_s, lun_t, lun_l
c      parameter (lun_s=76, lun_t=75, lun_l=74)
c      real (kind=R8) :: uu_xe, uu_ye
      logical uu_dbg
c      character tab, uu_file*48
c>>>
c======================================================================
c
c<<<
       call subini ('rlcomp')
c      tab=char(9)
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
c      print *,'rlcomp: nsrc,ntrg = ',nsrc,ntrg
c      open(lun_s,file='sources.dat')
c      write(lun_s,'(6a)') '#Radiation sources : ',label
c      write(lun_s,'(6a)') '#R',tab,'Z',tab,'S'
c>>>
      ttsrc=0.e0_R8
      radld=0.e0_R8
      do i=1,nsrc
        ssrc=ssrcm(i)
        ttsrc=ttsrc+ssrc
        xsrc=xysrc(1,i)
        ysrc=xysrc(2,i)
c<<<
c        uu_dbg= i .eq. 1
c        print *
c        print *,'Source number, ssrc, xsrc, ysrc'
c        print '(1p,2x,i4,3e10.2)',i, ssrc, xsrc, ysrc
c        write(lun_s,*) xsrc,tab,ysrc,tab,ssrc
c>>>
c
c*** Loop over the target segments
c
c<<<
c        if(uu_dbg) then
c          open(lun_t,file='targ-rlc.dat')
c          write(lun_t,'(6a)') '#Target geometry : ',label
c          write(lun_t,'(20a)') '#Rb',tab,'Zb',tab,'Re',tab,'Ze',tab,
c     ,    'ddx',tab,'ddy',tab,'uda',tab,'anglb',tab,'angle',tab,'radld'
c          open(lun_l,file='layout-t.dat')
c          write(lun_l,'(6a)') '#Target layout : ',label
c          write(lun_l,'(20a)') '#R',tab,'Z'
c        end if
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
c          uu_xe=xc
c          uu_ye=yc
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
c          if(uu_dbg) then
c            write(lun_t,*) xb,tab,yb,tab,uu_xe,tab,uu_ye,tab,
c     ,        ddx,tab,ddy,tab,uda/udpi,tab,anglb,tab,angle,tab,radld(j)
c            write(lun_l,*) xb,tab,yb
c            write(lun_l,*) uu_xe,tab,uu_ye
c            write(lun_l,*)
c           end if
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
c      close(lun_s)
c      close(lun_t)
c      close(lun_l)
c>>>
      call subend ()
      return
      end subroutine rlcomp
