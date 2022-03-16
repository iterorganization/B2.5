      subroutine interp2d(nx,ny,nx1,ny1,xx,yy,xx1,yy1,fa,fb)
      use b2mod_types
      use b2mod_subsys
      implicit none
*..   input arguments
      integer nx, ny, nx1, ny1
      real (kind=R8) :: fa(nx,ny)
      real (kind=R8) :: xx(nx), yy(ny)
      real (kind=R8) :: xx1(nx1), yy1(ny1)
*..   output arguments
      real (kind=R8) :: fb(nx1,ny1)
*..   local variables
      integer i, j
      integer ib, ih, it, jb, jh, jt
      real (kind=R8) :: fx, fy
*..   procedures
      external xertst

cxpb
c     This routine is intended to do a 2D bi-linear interpolation
c     of the fa(nx,ny) array onto the fb(nx1,ny1) array, with their
c     nodes located on two rectangular grids provided by the
c     xx and yy, and xx1 and yy1 arrays respectively.
cxpb

      call subini('interp2d')
      call xertst(1.le.nx, 'faulty input nx')
      call xertst(1.le.ny, 'faulty input ny')
      call xertst(1.le.nx1, 'faulty input nx1')
      call xertst(1.le.ny1, 'faulty input ny1')
      do i = 1, nx1
        call xertst(xx1(i).ge.xx(1).and.xx1(i).le.xx(nx),
     .   'interp2d not meant for X extrapolation!')
        ib=1
        it=nx
        do while (it-ib.gt.1)
          ih=(ib+it)/2
	  if(xx1(i).eq.xx(it)) then
	    ib=it-1
          elseif(xx1(i).lt.xx(ih)) then
            it=ih
          elseif(xx1(i).gt.xx(ih)) then
            ib=ih
          elseif(xx1(i).eq.xx(ih)) then
	    ib=ih
	    it=ih+1
          endif
        end do
        call xertst(xx(it).gt.xx(ib), 'XX array non-monotonic!')
        fx=(xx1(i)-xx(ib))/(xx(it)-xx(ib))
        do j = 1, ny1
          call xertst(yy1(j).ge.yy(1).and.yy1(j).le.yy(ny),
     .     'interp2d not meant for Y extrapolation!')
          jb=1
          jt=ny
          do while (jt-jb.gt.1)
            jh=(jb+jt)/2
            if(yy1(j).eq.yy(jt)) then
              jb=jt-1
            elseif(yy1(j).lt.yy(jh)) then
              jt=jh
            elseif(yy1(j).gt.yy(jh)) then
              jb=jh
            elseif(yy1(j).eq.yy(jh)) then
              jb=jh
              jt=jh+1
            endif
          end do
          call xertst(yy(jt).gt.yy(jb), 'YY array non-monotonic!')
          fy=(yy1(j)-yy(jb))/(yy(jt)-yy(jb))
          fb(i,j)=fx*fy*fa(it,jt)+fx*(1.0_R8-fy)*fa(it,jb)+
     &     (1.0_R8-fx)*fy*fa(ib,jt)+(1.0_R8-fx)*(1.0_R8-fy)*fa(ib,jb)
        end do
      end do

      call subend()
      return
      end subroutine interp2d
