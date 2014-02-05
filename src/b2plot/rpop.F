c
      function rpop(rstack,irstack)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   rstack(1000),rpop
      integer irstack
      if(irstack .eq. 0) then
        write(*,*) 'rstack underflow'
        rpop=0.0_R8
      else
        rpop=rstack(irstack)
        irstack=irstack-1
      endif
      return
      end
