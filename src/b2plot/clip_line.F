
      subroutine clip_line(xas,yas,xbs,ybs,xlo,ylo,xhi,yhi)
      use b2mod_types
      implicit none
      real (kind=R4) :: xas,yas,xbs,ybs,xlo,ylo,xhi,yhi
      real (kind=R4) :: dy,dx,xt,yt

      dy=ybs-yas
      dx=xbs-xas
      if(xas.lt.xlo.or.xbs.lt.xlo) then
        yt=ybs-dy/dx*(xbs-xlo)
        if(yt.ge.ylo.and.yt.le.yhi) then
          if(xas.lt.xbs) then
            xas=xlo
            yas=yt
          else
            xbs=xlo
            ybs=yt
          endif
        endif
      endif
      if(xas.gt.xhi.or.xbs.gt.xhi) then
        yt=ybs-dy/dx*(xbs-xhi)
        if(yt.ge.ylo.and.yt.le.yhi) then
          if(xas.gt.xbs) then
            xas=xhi
            yas=yt
          else
            xbs=xhi
            ybs=yt
          endif
        endif
      endif
      if(yas.lt.ylo.or.ybs.lt.ylo) then
        xt=xbs-dx/dy*(ybs-ylo)
        if(xt.ge.xlo.and.xt.le.xhi) then
          if(yas.lt.ybs) then
            xas=xt
            yas=ylo
          else
            xbs=xt
            ybs=ylo
          endif
        endif
      endif
      if(yas.gt.yhi.or.ybs.gt.yhi) then
        xt=xbs-dx/dy*(ybs-yhi)
        if(xt.ge.xlo.and.xt.le.xhi) then
          if(yas.gt.ybs) then
            xas=xt
            yas=yhi
          else
            xbs=xt
            ybs=yhi
          endif
        endif
      endif
      return
      end
