c
      subroutine sumy()
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      use b2us_data
      implicit none
      integer iCv,iC1,is,ix,iy
c
c apply operator
c
      if(mpg%isClassicalGrid.eq.0) return
      do is=0,ndim(2,iwrk)-1
        do ix=-1,mpg%nx
          iC1=0
          do iy=-1,mpg%ny
            iCv=mpg%imapCv(ix,iy)
            if (iCv.eq.0) cycle
            if (iC1.eq.0) then
              iC1 = iCv
            else
              work(iC1,is,iwrk)=work(iC1,is,iwrk)+work(iCv,is,iwrk)
            endif
          enddo
          if (iC1.eq.0) cycle
          do iy=-1,mpg%ny
            iCv=mpg%imapCv(ix,iy)
            if (iCv.eq.0) cycle
            work(iCv,is,iwrk)=work(iC1,is,iwrk)
          enddo
        enddo
      enddo
      values(9)=minval(work(:,0:ndim(2,iwrk)-1,iwrk))
      values(10)=maxval(work(:,0:ndim(2,iwrk)-1,iwrk))

      return
      end subroutine sumy
c
      subroutine inty()
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      use b2us_data
      implicit none
      integer ix,iy,is,iCv,iCn
c
c apply operator
c
      if(mpg%isClassicalGrid.eq.0) return
      do is=0,ndim(2,iwrk)-1
        do ix=-1,mpg%nx
          do iy=0,mpg%ny
            iCv=mpg%imapCv(ix,iy)
            iCn=mpg%imapCv(ix,iy-1)
            if(iCv.eq.0.or.iCn.eq.0) cycle
            work(iCv,is,iwrk)=work(iCn,is,iwrk)+work(iCv,is,iwrk)
          enddo
        enddo
      enddo
      values(9)=minval(work(:,0:ndim(2,iwrk)-1,iwrk))
      values(10)=maxval(work(:,0:ndim(2,iwrk)-1,iwrk))

      return
      end subroutine inty
