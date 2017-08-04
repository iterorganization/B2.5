
      function b2ye_tprf (val,nn) result(log)

      use b2mod_types
      implicit none

      integer, intent(in) :: nn
      real (kind=R8), intent(in) :: val(nn)
      logical :: log

      real (kind=R8) :: vmin, vmax

      vmin = minval(val(1:nn))
      vmax = maxval(val(1:nn))

      log = abs(vmax - vmin)/abs(max(vmin,vmax)+1.0e-30_R8) > 1.0e-5_R8

      return
      end function b2ye_tprf
