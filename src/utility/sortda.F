      subroutine sortda(a,ind,n)
      use b2mod_types
c version : 08.06.2001 00:35
c=======================================================================
c*** sorting real*8 array in ascending order - "quick and dirty"
c
      implicit none
      integer n,ind(n),i,j,k,l,m,ll
      real (kind=R8) :: a(n),r,u
c=======================================================================

      do i=1,n !{
        ind(i)=i
      end do !}
      if(n.le.1) return

      do k=2,n !{
        m=k-1
        ll=ind(k-1)
        r=a(ll)
        do j=k,n !{
          l=ind(j)
          u=a(l)
          if(u.lt.r) then !{
            ll=l
            r=u
            m=j
          end if !}
        end do !}
        if(m.ne.k-1) then !{
          i=ind(k-1)
          ind(k-1)=ll
          ind(m)=i
        end if !}
      end do !}
c=======================================================================
      end subroutine sortda
