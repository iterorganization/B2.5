      subroutine ipush(istack,iistack,ival)
      use b2mod_types
      implicit none
      integer istack(1000),ival,iistack
      if(iistack .eq. 1000) then
        write(*,*) 'istack overflow'
      else
        iistack=iistack+1
        istack(iistack)=ival
      endif
      return
      end
