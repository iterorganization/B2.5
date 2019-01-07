      subroutine pagini (style, xmin, xmax, ymin, ymax, xlab, ylab,
     &  title)
      use b2mod_types
      implicit none
      integer style
      real (kind=R8) ::
     *   xmin, xmax, ymin, ymax
      character xlab*(*), ylab*(*), title*(*)
*     -----------------------------------------------------------------
*     PAGINI starts a new plot.
*     The argument (style) specifies the intended style of the plot.
*     - style.eq.0: lin(x)-lin(y) curves
*     - style.eq.1: lin(x)-log(y) curves
*     - style.eq.2: log(x)-lin(y) curves
*     - style.eq.3: log(x)-log(y) curves
*     - style.eq.4: lin(x)-lin(y) contours
*     - style.eq.5: lin(x)-log(y) contours
*     - style.eq.6: log(x)-lin(y) contours
*     - style.eq.7: log(x)-log(y) contours
*     This version is for NCARGKS.
*     -----------------------------------------------------------------
      integer styl4
      character xlab0*16, ylab0*16, titl0*60, titl1*60
      real (kind=R8) ::
     *   aspmin, t0, xa, xb, ya, yb
      parameter (aspmin=1.0_R8)
      intrinsic max, abs, mod, len
#ifndef PLPLOT
      real (kind=R4) xmin4, ymin4
      external xertst, pagdef, gsln, gschup_d, gtx_d
      external set_d, labmod, tick4, gridal
#endif
*     -----------------------------------------------------------------
*!!   The parameter aspmin may be changed within the range 0.lt.aspmin
*     .le.1.0 It defines the minimum aspect ratio on contour plots.
*     aspmin=1.0 means the domain will always be expanded to fill a
*     square. Reasonable values for aspmin are in the range 0.6-1.0.
*     -----------------------------------------------------------------
*   ..test domain parameters
      call xertst (xmin.lt.xmax.and.ymin.lt.ymax,
     &  'pagini--faulty domain specification for '//title)
*   ..set styl4
      if (0.le.style.and.style.lt.8) then
       styl4 = mod(style,4)
      else
       styl4 = 0
      endif
*   ..set xlab0, ylab0, titl0, titl1
      xlab0 = xlab
      ylab0 = ylab
      titl0 = title
      if (len(titl0).lt.len(title)) then
       titl1 = title(len(titl0)+1:len(title))
      else
       titl1 = ' '
      endif
#ifndef PLPLOT
*   ..set GKS and NCAR parameters to default
      call pagdef ()
#endif
*   ..compute box dimensions
      xa = 0.14
      ya = 0.10
      if (style.eq.4) then
       if (abs(ymax-ymin).le.abs(xmax-xmin)) then
        t0 = abs(ymax-ymin)/abs(xmax-xmin)
        xb = xa+0.80
        yb = ya+max(aspmin,t0)*0.80
       else
        t0 = abs(xmax-xmin)/abs(ymax-ymin)
        xb = xa+max(aspmin,t0)*0.80
        yb = ya+0.80
       endif
      else
       xb = xa+0.80
       yb = ya+0.80
      endif
#ifndef PLPLOT
*   ..write title
      if (titl1.eq.' ') then
*      (one line)
       call gtx_d (xa-2.0e-2_R8, yb+4.0e-2_R8, titl0)
      else
*      (two lines)
       call gtx_d (xa-2.0e-2_R8, yb+6.0e-2_R8, titl0)
       call gtx_d (xa-2.0e-2_R8, yb+3.0e-2_R8, titl1)
      endif
*   ..write x label
      call gtx_d (xb-2.0e-1_R8, ya-5.4e-2_R8, xlab0)
*   ..write y label
      call gschup_d (-1.0_R8, 0.0_R8)
      call gtx_d (xa-4.0e-2_R8, yb-2.0e-1_R8, ylab0)
      call gschup_d (0.0_R8, 1.0_R8)
*   ..specify window and viewport
      call set_d (xa, xb, ya, yb, xmin, xmax, ymin, ymax, styl4+1)
*   ..unobtrusive ticks
      call tick4 (-14, -9, -14, -9)
*   ..format of tick labels
      call labmod ('(1pe9.2)', '(1pe9.2)', 0, 0, 14, 14, 18, 20, 0)
*   ..draw background
*!!   (The result is not satisfactory. The labels should be placed
*     at 'nice' values of the abscissae, rather than at points of a
*     uniform subdivision of the axes. In the case of log-plots I
*     would like to have three minor divisions per major division,
*     with ticks at values 2 and 5. I would like to see the numeric
*     y-labels rotated by 90 degrees.)
      xmin4=xmin
      ymin4=ymin
      if (style.eq.0) then
       call gridal (5,  2, 5,  2, 1, 1, 5, xmin4, ymin4)
      else if (style.eq.1) then
       call gridal (5,  2, 1, 10, 1, 1, 5, xmin4, ymin4)
       call gsln (2)
       call gridal (5,  2, 1, 10,-1, 0, 4, xmin4, ymin4)
       call gsln (1)
      else if (style.eq.2) then
       call gridal (1, 11, 5,  2, 1, 1, 5, xmin4, ymin4)
       call gsln (2)
       call gridal (1, 11, 5,  2, 0,-1, 1, xmin4, ymin4)
       call gsln (1)
      else if (style.eq.3) then
       call gridal (1, 11, 1, 10, 1, 1, 5, xmin4, ymin4)
       call gsln (2)
       call gridal (1, 11, 1, 10, 0, 0, 0, xmin4, ymin4)
       call gsln (1)
      else if (style.eq.4) then
       call gridal (5,  2, 5,  2, 1, 1, 5, xmin4, ymin4)
      else if (style.eq.5) then
       call gridal (5,  2, 1, 11, 1, 1, 5, xmin4, ymin4)
      else if (style.eq.6) then
       call gridal (1, 11, 5,  2, 1, 1, 5, xmin4, ymin4)
      else if (style.eq.7) then
       call gridal (1, 11, 1, 11, 1, 1, 5, xmin4, ymin4)
      else
       call gridal (1, 11, 1, 11, 1, 1, 5, xmin4, ymin4)
      endif
*   ..initialise clipping?
*!!   (clipping should be done by GKS or at the device level?)
*$$$  call gsclip (1)
#endif
      return
*     -----------------------------------------------------------------
      end subroutine pagini
