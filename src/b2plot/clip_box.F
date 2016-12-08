
      subroutine clip_box(xbox,ybox,xlo,ylo,xhi,yhi,n)
      use b2mod_types
      implicit none
      real xbox(8),ybox(8),xlo,ylo,xhi,yhi
      integer n
      real xtmp(8),ytmp(8)
      integer i,j

      j=0
      do i=1,n
        j=j+1
        xtmp(j)=xbox(i)
        ytmp(j)=ybox(i)
        xtmp(j+1)=xbox(mod(i,n)+1)
        ytmp(j+1)=ybox(mod(i,n)+1)
        if(xtmp(j).lt.xlo.and.xtmp(j+1).lt.xlo) then
          xtmp(j)=xlo
          xtmp(j+1)=xlo
        elseif(xtmp(j).gt.xhi.and.xtmp(j+1).gt.xhi) then
          xtmp(j)=xhi
          xtmp(j+1)=xhi
        endif
        if(ytmp(j).lt.ylo.and.ytmp(j+1).lt.ylo) then
          ytmp(j)=ylo
          ytmp(j+1)=ylo
        elseif(ytmp(j).gt.yhi.and.ytmp(j+1).gt.yhi) then
          ytmp(j)=yhi
          ytmp(j+1)=yhi
        endif
        call clip_line(xtmp(j),ytmp(j),xtmp(j+1),ytmp(j+1),
     1	 xlo,ylo,xhi,yhi)
        if(xtmp(j).lt.xlo.or.xtmp(j).gt.xhi) then
          if(xtmp(j).lt.xlo) then 
            xtmp(j)=xlo
          else
            xtmp(j)=xhi
          endif
          if(ytmp(j+1).lt.ylo) then
            ytmp(j)=ylo
          else
            ytmp(j)=yhi
          endif
        elseif(ytmp(j).lt.ylo.or.ytmp(j).gt.yhi) then
          if(ytmp(j).lt.ylo) then
            ytmp(j)=ylo
          else
            ytmp(j)=yhi
          endif
          if(xtmp(j+1).lt.xlo) then
            xtmp(j)=xlo
          else
            xtmp(j)=xhi
          endif
        else
          j=j+1
        endif
      enddo
      if(xtmp(j).eq.xtmp(1).and.ytmp(j).eq.ytmp(1)) j=j-1
      n=1
      i=1
      xbox(n)=xtmp(i)
      ybox(n)=ytmp(i)
      do i=2,j
        if(xtmp(i).ne.xtmp(i-1).or.ytmp(i).ne.ytmp(i-1)) then
          n=n+1
          xbox(n)=xtmp(i)
          ybox(n)=ytmp(i)
        endif
      enddo
      if(xbox(n).eq.xbox(1).and.ybox(n).eq.ybox(1)) n=n-1
      return
      end
