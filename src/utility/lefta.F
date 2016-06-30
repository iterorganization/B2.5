      subroutine lefta(c,n,m)
c
c  version : 05.07.89
c
c=====================================================
      implicit none
      character*1 c(*),b
      integer :: i,n,k,m
c=====================================================
c
      b=char(0)
      k=0
      do 100 i=1,n
         if(c(i).eq.' ' .or. c(i).eq.b) go to 100
         k=k+1
         c(k)=c(i)
 100  continue
      m=k
      if(m.lt.n) call bfill(c(k+1),n-k,' ')
      return
c=====================================================
      entry rghta(c,n,m)
c=====================================================
      b=char(0)
      k=n+1
      do 110 i=n,1,-1
         if(c(i).eq.' ' .or. c(i).eq.b) go to 110
         k=k-1
         c(k)=c(i)
 110  continue
      m=n+1-k
      if(m.lt.n) call bfill(c,n-k,' ')
c=====================================================
      end
