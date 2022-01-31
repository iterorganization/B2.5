! $AD NOCHECKPOINT
      function DIM(x, y)
      use b2mod_types
      implicit none
      REAL(kind=r8) :: dim, x, y, dum

      dim=max(x-y,0.0_R8)

      return
      end function dim

      
    
