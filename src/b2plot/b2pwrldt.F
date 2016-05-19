      subroutine b2pwrldt
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_wall_loading
c
c  version : 06.09.97 20:48
c
c======================================================================
c*** Calculates the radiation power load onto the targets
c*** Power is in W/m**2
c======================================================================
      implicit none
      integer i
c======================================================================
c*** This routine works only once
      if(l2rdldt) return
c======================================================================
c
      print *,'b2pwrldt starting...'
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
        write(*,*) '**** b2pwrldt: error detected in the input data ',
     ,                            '==> no radiation load calculated!!!'
        do i=1,ny
          wldrdld(lppp+ny+1-i)=0.0_R8
          wldrdld(lppp+ny+i)=0.0_R8
        end do
        return
      end if
c----------------------------------------------------------------------
c*** Fill in the input arrays for b2pwrdld
c
      nwldrad=2*ny
      do i=1,ny
        lwldrad(i)=lppp+i
        lwldrad(i+ny)=lppp+ny+i
      end do
c
c*** Calculate the load
c
      call b2pwrdld
      l2rdldt=.true.
c<<<
      print '(a/(1p,5e17.4))','Inner target: load',
     ,                                      (wldrdld(lppp+i),i=1,ny)
      print '(a/(1p,5e17.4))','Outer target: load',
     ,                                      (wldrdld(lppp+ny+i),i=1,ny)
c>>>
c======================================================================
      end
