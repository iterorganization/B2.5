      subroutine b2pwrldw
      use b2mod_b2plot_wall_loading
c
c  version : 06.09.97 21:35
c
c======================================================================
c*** Calculates radiation power load onto the specified wall segments
c*** Power is in W/m**2
c======================================================================
      implicit none
#include <COUPLE/KOPPLDIM.F>
      integer i,j,k,l,m
      logical ltg
c======================================================================
c*** This routine works only once
      if(l2rdldw) return
c======================================================================
c
      print *,'b2pwrldw starting...'
c*** Input the DG data and prepare the geometry
c
      call b2pwinp
c<<<
c      print '(a/5a10/1p,l10,3i10,e10.2)','After b2pwinp',
c     ,          'l2pwinp','lppp','nwldprt','nadsrc','add_rad',
c     ,                              l2pwinp,lppp,nwldprt,nadsrc,add_rad
c      print '(a/(16i5))','Shadowing structure',
c     ,                                     (lwldlis(i,1),i=1,nwldlis(1)
c      print '(a/(5(2f8.3,1x)))','Inner target',
c     ,                             (pp1(1,lppp+i),pp1(2,lppp+i),i=1,ny)
c      print '(a/(5(2f8.3,1x)))','Outer target',
c     ,                       (pp1(1,lppp+ny+i),pp1(2,lppp+ny+i),i=1,ny)
c>>>
      if(.not.l2pwinp) then
        write(*,*) '**** b2pwrldw: error detected in the input data ',
     ,                            '==> no radiation load calculated!!!'
        do i=1,nsgmx
          wldrdld(i)=0.0_R8
        end do
        return
      end if
c----------------------------------------------------------------------
c*** Fill in the input arrays for b2pwrdld
c
      ltg=.false.
      k=0
      do j=2,nwldprt
        do 100 i=1,nwldlis(j)
          l=iabs(lwldlis(i,j))
          if(l.gt.lppp) then
            ltg=.true.
          else
            do m=1,k
              if(l.eq.lwldrad(m)) go to 100
            end do
            k=k+1
            lwldrad(k)=l
          end if
 100    continue
      end do
      nwldrad=k
c
c*** Calculate the load
c
      call b2pwrdld
      l2rdldw=.true.
c
c*** If requested, calculate the target load
c
      if(ltg) then
        irdparm=1
        call b2pwrldt
      end if
c======================================================================
      end
