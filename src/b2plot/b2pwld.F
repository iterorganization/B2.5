      subroutine b2pwld(ig,BoRiS)
      use b2mod_types
c
c  version : 26.08.97 22:43
c
c======================================================================
      implicit none
      integer ig,ii
      real (kind=R8) :: BoRiS
c======================================================================
c
      print *,'b2pwld: ',ig
      ii=ig
      if(mod(ii,10) .gt. 0) call b2prdl
      ii=ii/10
      if(mod(ii,10) .gt. 0) call b2pnwl(*50)
 50   ii=ii/10
      if(mod(ii,10) .gt. 0) call b2ptrgl(BoRiS)
c======================================================================
      end
