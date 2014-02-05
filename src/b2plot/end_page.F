
      subroutine end_page
      use b2mod_types
#include <graphics_vars.h>

c     the zoom common
      integer winid
      integer wsnshwin,wsnX11
      common /zoom_com/ winid,wsnshwin,wsnX11

      integer wkid

      if(same_page.ne.0) then
        same_page=1
        write(*,*) 'same_page set to ',same_page
        return
      endif
      if(ix_page.ne.nx_page.and..not.new_page) then
        ix_page=ix_page+1
      elseif(iy_page.ne.ny_page.and..not.new_page) then
        iy_page=iy_page+1
        ix_page=1
      else
        ix_page=1
        iy_page=1
#ifndef PLPLOT
cxpb          write(*,*)'b2psh: vx   ',cfux(0.0)
cxpb          write(*,*)'b2psh: vy   ',cfuy(0.0)
cxpb          write(*,*)'b2psh: vxx  ',cfux(0.9999)
cxpb          write(*,*)'b2psh: vyy  ',cfuy(0.9999)
cxpb          write(*,*)'b2psh: read values!'
#ifdef NCAR4
        do wkid=1,nwork
          if (wkid.eq.wsnX11) then
            call ngpict(wkid,3)
          else
            if (wkid.eq.wsnshwin) then
              call ngpict(wkid,0)
            else
              call ngpict(wkid,1)
            endif
          endif
        enddo
#else
        call frame
#endif
#endif
        new_page=.true.
        dirty_page=.false.
      endif
      dirty_graph=.false.
      return
      end
