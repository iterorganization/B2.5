      subroutine b2pwld(ig,BoRiS) 
      use b2mod_types
      use b2mod_b2plot_wall_loading
c
c  version : 26.08.97 22:43
c
c======================================================================
      implicit none
      integer ig,ii,iret, use_eirene ! iyv
      real (kind=R8) :: BoRiS
      logical from_b2plot
      data from_b2plot / .true. /
      data use_eirene /0/
c======================================================================
c
      print *,'b2pwld: ',ig
      ildparm=mod(ig,100)
      ii=ig
      if(mod(ii,10) .gt. 0) call b2prdl
      ii=ii/10
      if(mod(ii,10) .gt. 0) call b2pnwl(iret)
      ii=ii/10
      irdparm=1
      if(mod(ii,10) .gt. 0) call b2ptrgl(BoRiS, from_b2plot, use_eirene)
      ii=ii/10
      irdparm=2
      if(mod(ii,10) .gt. 0) call b2pwldl(BoRiS)
c======================================================================
      return
      end subroutine b2pwld
