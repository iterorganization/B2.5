
      subroutine lubksb(a,n,np,indx,b)

      use b2mod_types

      implicit none
      
      integer ::
     *     n,np,indx(0:np)
      real (kind=R8) ::
     *     a(0:np,0:np),b(0:np)
      integer ::
     *     i,ii,j,ll
      real (kind=R8) ::  sum
      ii=-1
      do i=0,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.-1) then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
      enddo
      do i=n,0,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
      enddo
      return
      end
