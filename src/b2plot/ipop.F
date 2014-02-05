c
      integer function ipop(istack,iistack)
      use b2mod_types
      implicit none
      integer istack(1000),iistack
      if(iistack .eq. 0) then
        write(*,*) 'istack underflow'
        ipop=0
      else
        ipop=istack(iistack)
        iistack=iistack-1
      endif
      return
      end
