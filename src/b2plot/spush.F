c
      subroutine spush(sstack,isstack,sval)
      use b2mod_types
      implicit none
      character*256 sstack(1000),sval
      integer isstack
      if(isstack .eq. 1000) then
        write(*,*) 'sstack overflow'
      else
        isstack=isstack+1
        sstack(isstack)=sval
      endif
      return
      end
