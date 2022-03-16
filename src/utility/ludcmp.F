
      subroutine ludcmp(a,n,np,indx,d)
      
      use b2mod_types
      use b2mod_subsys
      implicit none
      
      integer ::
     *     n,np,indx(0:np),nmax
      real (kind=r8) ::
     *     d,a(0:np,0:np),tiny
      parameter (nmax=500,tiny=1.0e-20)
      integer ::
     *     i,imax,j,k
      real (kind=r8) ::
     *      aamax,dum,sum,vv(0:nmax-1)

      call subini('ludcmp')
      d=1.
      do i=0,n
        aamax=0.
        do j=0,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo
        call xertst (aamax.ne.0.0_R8, 'singular matrix in ludcmp')
        vv(i)=1./aamax
      enddo
      do j=0,n
        do i=0,j-1
          sum=a(i,j)
          do k=0,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo
        aamax=0.
        do i=j,n
          sum=a(i,j)
          do k=0,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
        enddo
        if (j.ne.imax)then
          do k=0,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=tiny
        if(j.ne.n)then
          dum=1./a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo
      call subend ()
      return
      end
