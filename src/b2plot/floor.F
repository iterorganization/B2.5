      real function floor(x)      
      use b2mod_types
      implicit none
      real :: x
      real, parameter :: onemeps=1.0-1.0e-5
      if(x.ge.0.) then
        floor=int(x)
      else
        floor=int(x-onemeps)
      endif
      return
      end
