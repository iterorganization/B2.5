c
      subroutine rpush(rstack,irstack,rval)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   rstack(1000),rval
      integer irstack
      if(irstack .eq. 1000) then
        write(*,*) 'rstack overflow'
      else
        irstack=irstack+1
        rstack(irstack)=rval
      endif
      return
      end
