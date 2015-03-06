c  VERSION : 31.08.97 11:43
c
c=====================================================
c*** The common blocks for rdload code.
c***
c*** 1. General geometry
c***  nplgm   the maximum number of segments in the "structure" polygon
c***  nplg    the actual number of segments in the "structure" polygon
c***  xp,yp   coordinates of the polygon corners (m)
c***  nbnd    the actual number of segments in the boundary polygon
c***  lbnd    list of indices of those segments
c***
c*** 2. Sources
c***  nsrcm   the maximum number of the source points
c***  nsrc    the actual number of the source points
c***  nsrcs   the number of the user-specified source points
c***  xysrc   coordinates of the source points (m)
c***  ssrcm   corresponding source intensities (MW)
c***  ttsrc   total source intensity (MW)
c***
c*** 3. Target specification
c***  ntrg    number of the segments belonging to the target
c***  ltrg    list of indices of those segments
c***  chords  additional lines to appear on the geometry layout
c***  nchrd   the actual number of the chords
c***
c*** 4. Output
c***  ztrg    coordinate - along the wall - of the target segments
c***  radld   radiation load (MW/m**2)
c***  tttrg   total power coming to the target (MW)
c***
c*** 5. Internal data
c***  xsrc,ysrc   coordinates of the current source point
c***  ssrc    corresponding source intensity
c***  jb      index of the base segment
c***  jc      index of the current segment
c***  jcb     index of the current segment in the boundary list
c***  jm      mode of the corner point (see rlshdw.F)
c***  jsc     sector of the corner point
c***  xb,yb   coordinates of the beginning of the base segment
c***  xc,yc   coordinates of the corner point
c***  xs,ys   coordinates of the end of the current segment
c***  ddx,ddy base segment dimensions in x and y directions
c***  sb,se   beginning and end of the visible part of the base
c***          segment normalized to its length
c=====================================================
      integer nplgm, nsrcm, nchrdm
      parameter (nplgm=400, nsrcm=4000, nchrdm=200)
      integer jb,jc,jm,jsc,jcb,nplg,nbndry,nsrc,nsrcs,ntrg,nchrd,
     ,  lbnd(nplgm),ltrg(nplgm)
      real (kind=R8) ::
     *   xsrc,ysrc,ssrc,sb,se,xb,yb,ddx,ddy,xc,yc,xs,ys,
     ,  ztrg(nplgm),radld(nplgm),xp(nplgm),yp(nplgm),
     ,  xysrc(2,nsrcm),ssrcm(nsrcm),ttsrc,tttrg,chords(2,2,nchrdm)
      common/shdcom/xsrc,ysrc,ssrc,sb,se,xb,yb,ddx,ddy,xc,yc,xs,ys,
     ,   jb,jc,jm,jsc,jcb
      common/rldcom/ztrg,radld,xp,yp,xysrc,ssrcm,ttsrc,tttrg,chords,
     i  nplg,nbndry,nsrc,nsrcs,ntrg,nchrd,lbnd,ltrg
