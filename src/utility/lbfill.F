      subroutine lbfill (iftp, acrx, acry, nacra, indx)
      use b2mod_types
      implicit none
      integer iftp, iftp1, nacra, indx
      real (kind=R4) acrx(0:nacra-1), acry(0:nacra-1)
*     ------------------------------------------------------------------
*     LBFILL is a substitute for a routine of the same name in the
*     LABELBAR utility of NCAR graphics. It fills a specified area.
*     LBFILX and LBFILQ are alternate entries in LBFILL. LBFILX is
*     called to define the filling style and LBFILQ is called to query
*     the filling style.
*     ------------------------------------------------------------------
*     The parameter iftp in the call to lbfill is ignored.
*     The parameter iftp in the call to lbfilx is copied into a local
*     'saved' variable iftp0. iftp0 defines the style of filling:
*      iftp0.le.0: No action.
*      iftp0.eq.1: Solid color fill.
*      iftp0.eq.2: Pattern fill, no color.
*      iftp0.eq.3: Hatch fill, no color.
*      iftp0.eq.4: Pattern fill using color.
*      iftp0.eq.5: Hash fill using color.
*      iftp0.gt.5: No action.
*     (Some styles may require special initialisation action when
*     GKS is opened; some styles may not be available on some output
*     devices.)
*     The parameters acrx, acry and nacra specify the fill region.
*     The parameter indx defines the color and pattern to be used.
*     (Color and pattern may depend upon GKS environment and on the
*     specific initialisation of GKS parameters.)
*     The parameter iftp1 in the call to lbfilq returns the saved
*     value of iftp0.
*     ------------------------------------------------------------------
      integer iftp0
      save iftp0
#ifndef PLPLOT
      external gsfais, gsfaci, gsfasi, gfa
#endif
      data iftp0/0/
*     ------------------------------------------------------------------
#ifndef PLPLOT
*   ..branch on style
      if (indx.lt.0.or.iftp0.lt.1) then
*!!    (or maybe indx.lt.1.or.iftp0.lt.1?)
*    ..no action
      else if (iftp0.eq.1) then
*    ..solid color fill
*      (It is assumed here that we have 30 genuine colors, each one
*      selected by a color index in the range (2:31); color indices 0
*      and 1 correspond to (black and white) or to (white and black).)
       call gsfais (1)
       if (indx.lt.2) then
        call gsfaci (indx)
       else
        call gsfaci (2+mod(indx+28,30))
       endif
       call gfa (nacra, acrx, acry)
      else if (iftp0.eq.2) then
*    ..pattern fill, no color
*      (assume we have 6 patterns)
       call gsfais (2)
       call gsfaci (1)
       call gsfasi (1+mod(indx+5,6))
       call gfa (nacra, acrx, acry)
      else if (iftp0.eq.3) then
*    ..hatch fill, no color
*      (assume we have 6 hatches)
       call gsfais (3)
       call gsfaci (1)
       call gsfasi (1+mod(indx+5,6))
       call gfa (nacra, acrx, acry)
      else if (iftp0.eq.4) then
*    ..pattern color fill
*      (assume we have 6 patterns and 30 colors)
       call gsfais (2)
       if (indx.lt.2) then
        call gsfaci (indx)
       else
        call gsfaci (2+mod(indx+28,30))
       endif
       call gsfasi (1+mod(indx+5,6))
       call gfa (nacra, acrx, acry)
      else if (iftp0.eq.5) then
*    ..hatch color fill
*      (assume we have 6 hatches and 30 colors)
       call gsfais (3)
       if (indx.lt.2) then
        call gsfaci (indx)
       else
        call gsfaci (2+mod(indx+28,30))
       endif
       call gsfasi (1+mod(indx+5,6))
       call gfa (nacra, acrx, acry)
      endif
      call gsfais (1)
      call gsfaci (1)
      call gsfasi (1)
#endif
      return
*     ------------------------------------------------------------------
      entry lbfilx (iftp)
      iftp0 = iftp
      return
*     ------------------------------------------------------------------
      entry lbfilq (iftp1)
      iftp1 = iftp0
      return
*     ------------------------------------------------------------------
      end
