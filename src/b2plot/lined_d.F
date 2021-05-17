
      subroutine lined_d(xa,ya,xb,yb)
      use b2mod_types
#include <graphics_vars.h>
      real (kind=R8) :: xa,xb,ya,yb
      real (kind=R4) :: xas,xbs,yas,ybs

      xas=min(1e30_R8,max(-1e30_R8,xa))
      xbs=min(1e30_R8,max(-1e30_R8,xb))
      yas=min(1e30_R8,max(-1e30_R8,ya))
      ybs=min(1e30_R8,max(-1e30_R8,yb))
#ifdef PLPLOT
      call pljoin(xas,yas,xbs,ybs)
#else
c      write(*,*) 'lined: ',xas,yas,xbs,ybs
      if(
     1 (xas.lt.xlo.and.xbs.lt.xlo).or.
     2 (xas.gt.xhi.and.xbs.gt.xhi).or.
     3 (yas.lt.ylo.and.ybs.lt.ylo).or.
     4 (yas.gt.yhi.and.ybs.gt.yhi)) return
      if(.not.logy.or.(yas.gt.0.and.ybs.gt.0)) then
        if(logy) then
          call lined(xas,yas,xbs,ybs)
        else
          call clip_line(xas,yas,xbs,ybs,xlo,ylo,xhi,yhi)
          if(xlo.le.xas.and.xas.le.xhi.and.
     1       ylo.le.yas.and.yas.le.yhi.and.
     1       xlo.le.xbs.and.xbs.le.xhi.and.
     1       ylo.le.ybs.and.ybs.le.yhi) call lined(xas,yas,xbs,ybs)
        endif
      elseif(yas.lt.0.and.ybs.lt.0) then
        call lastd
        call gsln(3)
        call lined(xas,-yas,xbs,-ybs)
        call gsln(1)
        call lastd
      elseif(abs(yas).gt.0.and.abs(ybs).gt.0) then
        call lastd
        call gsln(4)
        call lined(xas,abs(yas),xbs,abs(ybs))
        call gsln(1)
        call lastd
      else
      endif
#endif
      return
      end subroutine lined_d

      subroutine grpts_d(xa,ya,is)
      use b2mod_types
#include <graphics_vars.h>
      real (kind=R8) :: xa,ya
      real (kind=R4) :: xas,yas
      integer is

      xas=min(1e30_R8,max(-1e30_R8,xa))
      yas=min(1e30_R8,max(-1e30_R8,ya))
      if((xas.lt.xlo).or.(xas.gt.xhi).or.
     3   (yas.lt.ylo).or.(yas.gt.yhi)) return
      if(.not.logy.or.(yas.gt.0)) then
        call gsmk(is)
        call gpm(1,xas,yas)
      endif
      return
      end subroutine grpts_d
          
      

