      subroutine interp1d(nn,nn1,xx,xx1,fa,fb)
      use b2mod_types
      use b2mod_subsys
      implicit none
*..   input arguments
      integer nn, nn1
      real (kind=R8) :: fa(nn)
      real (kind=R8) :: xx(nn)
      real (kind=R8) :: xx1(nn1)
*..   output arguments
      real (kind=R8) :: fb(nn1)
*..   local variables
      integer i
      integer ib, ih, it
*..   procedures
      external xertst

csc
c     This routine performs simple 1D interpolation
csc

      call subini('interp1d')
      call xertst(1.le.nn, 'faulty input nn')
      call xertst(1.le.nn1, 'faulty input nn1')
      do i = 1, nn1
        call xertst(xx1(i).ge.xx(1).and.xx1(i).le.xx(nn),
     .   'interp1d not meant for X extrapolation!')
        ib=1
        it=nn
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
        fb(i)=(xx1(i)-xx(ib))/(xx(it)-xx(ib))*(fa(it)-fa(ib))+fa(ib)
      end do

      call subend()
      return
      end subroutine interp1d
