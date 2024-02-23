      function ceil(x)
      use b2mod_types
      implicit none
      real (kind=R4) :: ceil
      real (kind=R4) :: x, onemeps
      parameter (onemeps=1.0_R4-1.0e-5_R4)
      if(x.ge.0.0_R4) then
        ceil=int(x+onemeps)
      else
        ceil=int(x)
      endif
      return
      end function ceil
