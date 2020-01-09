      subroutine binsearch(xx,n,x,i)
c     ***********************************************************
c     * ermittlung des Feldindex i mit vorgegebener zahl x,     *
c     * so dass x zwischen xx(i) und xx(i+1)                    *
c     * ferner: falls x<=xx(1) i=1, und falls x>=xx(n) i=n      *
c     * *********************************************************
      use B2MOD_TYPES
      implicit none
 
      integer, intent(in) :: n
      integer, intent(out) :: i
      real(kind=R8), intent(in) :: xx(n), x
      integer :: bl, bm, bu
 
c  savest version:
c  all tests included
      entry binsearch_0(xx,n,x,i)
 
      bl=0
      bu=n+1
 
      if (xx(n).ge.xx(1)) then
! monoton increasing
        if(x.le.xx(1))then
          i=1
        else if(x.ge.xx(n))then
          i=n
        else
c  binary search
          do while (bu-bl.gt.1)
            bm=(bu+bl)*0.5
            if(x.ge.xx(bm)) then
              bl=bm
            else
              bu=bm
            endif
          enddo
          i=bl
        end if
      else
! monoton decreasing
        if(x.ge.xx(1))then
          i=1
        else if(x.le.xx(n))then
          i=n
        else
c  binary search
          do while (bu-bl.gt.1)
            bm=(bu+bl)*0.5
            if(x.le.xx(bm)) then
              bl=bm
            else
              bu=bm
            endif
          enddo
          i=bl
        end if
      end if
c
      return
 
c  fast version:
c  we already know: a)  xx is monotonically increasing (not decreasing)
c                   b)  x  lies between xx(1) and xx(n)
      entry binsearch_2(xx,n,x,i)
 
      bl=0
      bu=n+1
 
c  binary search
      do while (bu-bl.gt.1)
        bm=(bu+bl)*0.5
        if(x.ge.xx(bm))then
          bl=bm
        else
          bu=bm
        endif
      end do
c
      i=bl
 
      return
      end
