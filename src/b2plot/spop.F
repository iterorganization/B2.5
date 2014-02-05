c
      function spop(sstack,isstack)
      use b2mod_types
      implicit none
      character*256 sstack(1000)
      character*256 spop
      integer isstack
      if(isstack .eq. 0) then
        write(*,*) 'sstack underflow'
        spop=' '
      else
        spop=sstack(isstack)
        isstack=isstack-1
      endif
      return
      end
