      real function ceil(x)
      use b2mod_types
      real x, onemeps
      parameter (onemeps=1.0-1.0e-5)
      if(x.ge.0.0) then
        ceil=int(x+onemeps)
      else
        ceil=int(x)
      endif
      return
      end
