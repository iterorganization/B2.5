      subroutine rlshdw
      use b2mod_types
      use b2mod_subsys
      implicit none
c
c  version : 07.11.97 13:19
c
c======================================================================
c*** Computation of the shadowing of the base segment by the structure.
c*** This is a 2-D calculation that uses a cylindrical approximation.
c*** The structure is supposed to be a closed polygon without self-
c*** intersections. The base belongs to the same polygon.
c*** Orientation of the polygon is important: the radiation sources are
c*** supposed to be located "inside", and this is the left-hand side
c*** when moving along the polygon.
c***
c*** All the data are transferred via the common blocks.
c***
c*** Input
c***  nplgs   number of the segments in the structure
c***  xp,yp   coordinates of the corners of the structure
c***  jcb     index of the "base" segment
c***  [xy]src coordinates of the source
c***
c***  xb,yb   coordinates of the base starting point
c***  sb,se   starting and ending points of the irradiated part of the
c***          base segment normalised to the base length.
c***          Should be set to 0 and 1 before the call.
c***  ddx,ddy   base dimensions in x and y directions
c***
c*** Output
c***  sb,se   starting and ending points of the irradiated part of the
c***          base segment normalised to the base length.
c***          sb >= se means total shadowing (no irradiation)
c======================================================================
#include "DIMENSIONS.F"
#include "rdload.inc"
      real (kind=R8) :: ueps
      parameter (ueps=1.0e-20_R8)
      integer i, jj, js, k
      real (kind=R8) :: d, ux, uy, utc, utb, udx, udy, uxe, uye
c======================================================================
      intrinsic min, max, mod, sign
      external xertst
c======================================================================
c*** Initialise the code. The base segment is at the position jb
c*** in the structure list, we need its endpoints (xb,yb,uxe,uye).
c*** jm and jsc are set to 1 because we start from the base endpoint.
c*** jcb is the segment counter for the structure, we set it thus to
c*** mod(jb,nbndry)+1 to specify the segment next to the base.
c
      call subini('rlshdw')
c      jm=1
c      jsc=1
c      jcb=mod(jb,nbndry)+1
c      k=lbnd(jb)
c      xb=xp(k)
c      yb=yp(k)
c      k=lbnd(jcb)
c      uxe=xp(k)
c      uye=yp(k)
      uxe=xc
      uye=yc
c
c*** Check the orientation of the base with respect to the source.
c*** We assume that all the radiation comes from the "inner" side of
c*** the wall segments; therefore, they must be oriented properly.
c
      udx=xsrc-xb
      udy=ysrc-yb
      d=ddx*udy-ddy*udx
c<<<
c      print *
c      print *,'rlshdw called.'
c      print '(2x,3a4,9(3x,a4,4x))','jsc','jm','jcb','xsrc','ysrc',
c     ,               ' xb ',' yb ',' ddx ',' ddy ','udx ','udy ',' d  '
c      print '(1x,3i4,1p,9e11.3)',
c     ,                     jsc,jm,jcb,xsrc,ysrc,xb,yb,ddx,ddy,udx,udy,d
c>>>
      if(d .le. ueps) then
c
c*** The base has its back turned to the source - no irradiation!
c
        se=0.0_R8
        sb=1.0_R8
        call subend ()
        return
      end if
c
c*** Now check the shadowing by all the segments.
c*** Here we assume that the "base" (i.e. the segment which is being
c*** examined) belongs to the same structure which may shadow it.
c*** Therefore, we start the procedure from the end of the base segment
c*** and go around the structure until the beginning of the base.
c
      do i=2,nbndry
c
c***  Take the next point (js, xs, ys) and determine its sector (jj)
c
        js=mod(jcb,nbndry)+1
        k=lbnd(js)
        xs=xp(k)
        ys=yp(k)
        ux=xs-xsrc
        uy=ys-ysrc
        d=ux*ddy-uy*ddx
c<<<
c      print *
c      print '(2x,3a4,5(3x,a4,4x))','  i ',' jcb',' js ',
c     ,                               ' xs ',' ys ',' ux ',' uy ',' det'
c      print '(1x,3i4,1p,5e11.3)',i,jcb,js,xs,ys,ux,uy,d
c>>>
        if(abs(d) .gt. ueps) then
          d=1.0_R8/d
c
c*** utb and utc should be the parameters dertermining the point of
c*** intersection of the base segment with the straight line passing
c*** through the source S(xsrc,ysrc) and the current point s(xs,ys)
c***
c*** utb is the distance along the base from the base beginning (xb,yb)
c***     normalised to the base length. Thus 0 =< utb =< 1 means that
c***     the intersection point belongs to the base segment.
c*** utc is the distance from S to the intersection point expressed in
c***     the units of [S,s]. Thus 0 < utc < 1 means that the point is
c***     hidden by the base, and utc > 1 means that the point
c***     creates a shadow on the base
c
          utb=(ux*udy-uy*udx)*d
          utc=(ddx*udy-ddy*udx)*d
c          utc=(ddx*uy-ddy*ux)*d
c
c*** Determine the sector JJ to which the current point belongs.
c*** Only the segments protruding the sector 1 can cause shadowing.
c***
c***     The sectors layout         The point modes
c***            \ 3 /
c***             \ /              jm=0  no interference with the base
c***           4  S  2               1  shadowing the end point
c***             / \                 2  shadowing the starting point
c***            / 1 \
c***           /=====>
c***          /  base \
c
          if(utb .gt. 1.0_R8) then
              if(utc .gt. 0.0_R8) then
                  jj=2
              else
                  jj=4
              end if
          else if(utb .lt. 0.0_R8) then
              if(utc .lt. 0.0_R8) then
                  jj=2
              else
                  jj=4
              end if
          else if(utc .lt. 0.0_R8) then
              jj=3
          else
              jj=1
          end if
c<<<
c      print *,'abs(d) > ueps. jj, utb, utc: ',jj,utb,utc
c>>>
        else
c
c*** No intersection: the ray is parallel to the base
c
          utc=10.0_R8
          if(ux*ddx+uy*ddy .gt. 0.0_R8) then
              jj=2
              utb=2.0_R8
          else
              jj=4
              utb=-2.0_R8
          end if
c<<<
c      print *,'abs(d) < ueps. jj, utb, utc: ',jj,utb,utc
c>>>
        end if
c
c*** Analyse the transition between sectors and the current mode
c
c<<<
c      print *,'Before analysis. jsc, jj, jm: ', jsc, jj, jm
c>>>
c
c*** Here jsc is the sector to which the other end of the current
c*** segment - that is, the previous current point - belongs, and
c*** jm is the "mode" of the previous point.
c***
c*** The idea is to determine the mode of the current point on the
c*** basis of the analysis of the transition between the sectors
c*** when following the segments of the structue, and to calculate
c*** the strongest shadowing from both ends of the base.
c***
c*** The following table shows the mode of the point after transition
c***
c***----------------------------------------------------------\
c***  To       (jj) :     1           2         3         4   |
c***  From (jsc,jm) /---------------------------------------- |
c***         1   0  |     0           0         0         0   |
c***         1   1  |     1           0        0-*        *   |
c***         1   2  |     2           *        0-*        0   |
c***         2      |    0-1          0         0        0-*  |
c***         3      |   0-1-2         0         0         0   |
c***         4      |    0-2         0-*        0         0   |
c***----------------------------------------------------------/
c
        if(jsc .eq. 1) then
c
c*** The segment starts in the sector 1. Where does it finish?
c
c<<<
c      print *
c      print *,'jsc = 1. jj, jm:',jj,jm
c      print '(8(4x,a4,3x))',' xc ',' yc ',' xs ',' ys ',
c     ,                                      ' xb ',' yb ','xsrc','ysrc'
c      print '(1p,1x,8e11.3)',xc, yc, xs, ys, xb, yb, xsrc, ysrc
c>>>
          if(jj .eq. 1) then
c*** The same sector. jm remains the same - unless for the first segment
            if(i .eq. 2) then
              if(utc .gt. 1.0_R8) then
                jm=1
              else
                jm=0
              end if
            end if
          else if(jj .eq. 2) then
c*** ... the sector 2. Check the history
            if(jm .eq. 2) then
c***        the structure crosses sight completely: total shading
              se=0.0_R8
              sb=1.0_R8
              call subend ()
              return
            end if
c***        no shading in sector 2
            jm=0
          else if(jj .eq. 3) then
c*** ... the sector 3. Check the history
            if(jm .eq. 1) then
c<<<
c      print *,'(1,1) -> 3.  f(b),f(S):',flin(xb,yb),flin(xsrc,ysrc)
c>>>
              if(sign(1.0_R8,flin(xb,yb))  .ne. 
     .         sign(1.0_R8,flin(xsrc,ysrc))) then
c***        the structure crosses sight completely: total shading
                se=0.0_R8
                sb=1.0_R8
                call subend ()
                return
              end if
c***        no shading
              jm=0
            else if(jm .eq. 2) then
c<<<
c      print *,'(1,2) -> 3.  f(e),f(S):',flin(uxe,uye),flin(xsrc,ysrc)
c>>>
              if(sign(1.0_R8,flin(uxe,uye))  .ne. 
     .         sign(1.0_R8,flin(xsrc,ysrc))) then
c***        the structure crosses sight completely: total shading
                se=0.0_R8
                sb=1.0_R8
                call subend ()
                return
              end if
c***        no shading
              jm=0
            end if
          else if(jj .eq. 4) then
c*** ... the sector 4. Check the history
            if(jm .eq. 1) then
c***     the structure crosses sight completely: total shading
              se=0.0_R8
              sb=1.0_R8
              call subend ()
              return
            end if
c***        no shading
            jm=0
          end if
        else if(jsc .eq. 2) then
c
c*** The segment starts in the sector 2. Where does it finish?
c
          if(jj .eq. 1) then
c*** ... the sector 1. Check the shading
            if(utc .lt. 1) then
                jm=0
            else
                jm=1
            end if
          else if(jj .eq. 4) then
c*** ... the sector 4. Check the shading
c<<<
c      print *,'2 -> 4.  f(b),f(S):',flin(xb,yb),flin(xsrc,ysrc)
c>>>
            if(sign(1.0_R8,flin(xb,yb))  .ne. 
     .       sign(1.0_R8,flin(xsrc,ysrc))) then
c***        the structure crosses sight completely: total shading
              se=0.0_R8
              sb=1.0_R8
              call subend ()
              return
            end if
            jm=0
          end if
        else if(jsc .eq. 3) then
c
c*** The segment starts in the sector 3. Where does it finish?
c
          if(jj .eq. 1) then
c*** ... the sector 1. Check the shading
c<<<
c      print *,'3 -> 1.  utc,f(b),f(S):',utc,flin(xb,yb),flin(xsrc,ysrc)
c>>>
            if(utc .lt. 1.0_R8) then
              jm=0
            else
              if(sign(1.0_R8,flin(xb,yb))  .ne. 
     .         sign(1.0_R8,flin(xsrc,ysrc))) then
                jm=2
              else
                jm=1
              end if
            end if
          end if
        else if(jsc .eq. 4) then
c
c*** The segment starts in the sector 4. Where does it finish?
c
c<<<
c      print *
c      print *,'jsc = 4. jj,utc:',jj,utc
c      print '(8(4x,a4,3x))',' xc ',' yc ',' xs ',' ys ',
c     ,                                      ' xb ',' yb ','xsrc','ysrc'
c      print '(1p,1x,8e11.3)',xc, yc, xs, ys, xb, yb, xsrc, ysrc
c>>>
          if(jj .eq. 1) then
c*** ... the sector 1. Check the shading
            if(utc .gt. 1.0_R8) then
              jm=2
            else
              jm=0
            end if
          else if(jj .eq. 2) then
c*** ... the sector 1. Check the shading
c<<<
c      print *,'4 -> 1.  f(e),f(S):',flin(uxe,uye),flin(xsrc,ysrc)
c>>>
            if(sign(1.0_R8,flin(xb,yb))  .ne. 
     .       sign(1.0_R8,flin(xsrc,ysrc))) then
c***        the structure crosses sight completely: total shading
              se=0.0_R8
              sb=1.0_R8
              call subend ()
              return
            end if
          end if
        else
          write(*,*) 'Unexpected error in rlshdw: jsc =',jsc
          call xerrab ('Unexpected error in rlshdw')
        end if
c<<<
c      print *,'Before correction. jm, jj, sb, se, utb:',jm,jj,sb,se,utb
c>>>
c
c*** Correct the shadowing data.
c*** Only points in sector 1, mode 1 or 2, can produce a partial shadow
c
        if(jj .eq. 1) then
          if(jm .eq. 1) then
            se=min(se,utb)
            if(se .le. sb) then
              call subend ()
              return
            end if
          else if(jm .eq. 2) then
            sb=max(sb,utb)
            if(sb .ge. se) then
              call subend ()
              return
            end if
          end if
        end if
        jsc=jj
        jcb=js
        xc=xs
        yc=ys
      end do
      call subend ()
c======================================================================
      return

      contains

      function flin(x,y)
      implicit none
      real (kind=R8), intent(in) :: x, y
      real (kind=R8) :: flin

      flin = (x-xc)*(ys-yc)-(y-yc)*(xs-xc)

      return
      end function flin

      end subroutine rlshdw
