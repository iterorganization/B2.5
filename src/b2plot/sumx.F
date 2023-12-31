c
      subroutine sumx()
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
        do iy=-1,mpg%ny
          iC1=0
          do ix=-1,mpg%nx
            iCv=mpg%imapCv(ix,iy)
            if (iCv.eq.0) cycle
            if (iC1.eq.0) then
              iC1 = iCv
            else
              work(iC1,is,iwrk)=work(iC1,is,iwrk)+work(iCv,is,iwrk)
            endif
          enddo
          if (iC1.eq.0) cycle
          do ix=-1,mpg%nx
            iCv=mpg%imapCv(ix,iy)
            if (iCv.eq.0) cycle
            work(iCv,is,iwrk)=work(iC1,is,iwrk)
          enddo
        enddo
      enddo
      values(9)=minval(work(:,0:ndim(2,iwrk)-1,iwrk))
      values(10)=maxval(work(:,0:ndim(2,iwrk)-1,iwrk))

      return
      end subroutine sumx
c
      subroutine sumt()
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer i,is
      if(nedim(1,iewrk) .lt. 0) then
        write(*,*) 'To apply this operator, an array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
      do is=0,nedim(2,iewrk)-1
        do i= 2,nedim(1,iewrk)
          eirene_stack(1,is,iewrk)=
     .     eirene_stack(1,is,iewrk)+eirene_stack(i,is,iewrk)
        enddo
      enddo
      values(38)=eirene_stack(1,0,iewrk)
      values(39)=eirene_stack(1,0,iewrk)
      do is=0,nedim(2,iewrk)-1
        do i=2,nedim(1,iewrk)
          eirene_stack(i,is,iewrk)=eirene_stack(1,is,iewrk)
        enddo
        values(38)=min(values(38),eirene_stack(1,is,iewrk))
        values(39)=max(values(39),eirene_stack(1,is,iewrk))
      enddo

      return
      end subroutine sumt
c
      subroutine sumw()
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer i,is
      if(nwdim(1,iwwrk) .lt. 0) then
        write(*,*) 'To apply this operator, an array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
      do is=0,nwdim(2,iwwrk)-1
        do i= 2,nwdim(1,iwwrk)
          wall_stack(1,is,iwwrk)=
     .     wall_stack(1,is,iwwrk)+wall_stack(i,is,iwwrk)
        enddo
      enddo
      values(14)=wall_stack(1,0,iwwrk)
      values(15)=wall_stack(1,0,iwwrk)
      do is=0,nwdim(2,iwwrk)-1
        do i=2,nwdim(1,iwwrk)
          wall_stack(i,is,iwwrk)=wall_stack(1,is,iwwrk)
        enddo
        values(14)=min(values(14),wall_stack(1,is,iwwrk))
        values(15)=max(values(15),wall_stack(1,is,iwwrk))
      enddo

      return
      end subroutine sumw
c
      subroutine intx()
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      use b2us_data
      implicit none
      integer i,iCv,iCp,iFt,is
      intrinsic min, max
c
c apply operator
c
      do is=0,ndim(2,iwrk)-1
        do iFt=1,mpg%nFt
          do i = 2, mpg%ftCvP(iFt,2)
            iCv = mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
            iCp = mpg%ftCv(mpg%ftCvP(iFt,1)+i-2)
            work(iCv,is,iwrk)=work(iCv,is,iwrk)+work(iCp,is,iwrk)
          enddo
        enddo
      enddo
      values(9)=minval(work(:,0:ndim(2,iwrk)-1,iwrk))
      values(10)=maxval(work(:,0:ndim(2,iwrk)-1,iwrk))

      return
      end subroutine intx
c
      subroutine intw()
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer iw,iseg,is
      intrinsic min, max
      if(nwdim(1,iwwrk) .le. 1) then
        write(*,*) 'To apply this operator, an array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
      do is=0,nwdim(2,iwwrk)-1
        do iseg = 1, nsegments
          do iw=wall_begin(iseg)+1,wall_end(iseg)
            wall_stack(iw,is,iwwrk)=
     .       wall_stack(iw-1,is,iwwrk)+
     .       wall_stack(iw,is,iwwrk)
          enddo
        enddo
      enddo
      values(14)=wall_stack(1,0,iwwrk)
      values(15)=wall_stack(1,0,iwwrk)
      do is=0,nwdim(2,iwwrk)-1
        do iw=1,nbnd
          values(14)=min(values(14),wall_stack(iw,is,iwwrk))
          values(15)=max(values(15),wall_stack(iw,is,iwwrk))
        enddo
      enddo

      return
      end subroutine intw
c
      subroutine intwstar()
c same as intw but no account of segments
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer iw,is
      intrinsic min, max
      if(nwdim(1,iwwrk) .le. 1) then
        write(*,*) 'To apply this operator, an array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
      do is=0,nwdim(2,iwwrk)-1
        do iw = 2, nwdim(1,iwwrk)
          wall_stack(iw,is,iwwrk)=
     1     wall_stack(iw-1,is,iwwrk)+wall_stack(iw,is,iwwrk)
        enddo
      enddo
      values(14)=wall_stack(1,0,iwwrk)
      values(15)=wall_stack(1,0,iwwrk)
      do is=0,nwdim(2,iwwrk)-1
        do iw=1,nbnd
          values(14)=min(values(14),wall_stack(iw,is,iwwrk))
          values(15)=max(values(15),wall_stack(iw,is,iwwrk))
        enddo
      enddo

      return
      end subroutine intwstar
