      function sfloor(x)
      use b2mod_types
      implicit none
      real (kind=R4) :: sfloor
      real (kind=R4) :: x
      real, parameter :: onemeps=1.0_R4-1.0e-5_R4
      if(x.ge.0._R4) then
        sfloor=int(x)
      else
        sfloor=int(x-onemeps)
      endif
      return
      end function sfloor
