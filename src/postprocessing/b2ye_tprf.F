
      function b2ye_tprf (val,nn) result(log)

      implicit none

      integer, intent(in) :: nn
      real*8, intent(in) :: val(nn)
      logical :: log

      real*8 :: vmin, vmax

      vmin = minval(val(1:nn))
      vmax = maxval(val(1:nn))

      log = abs(vmax - vmin)/abs(max(vmin,vmax)+1.d-30) > 1.d-5

      return
      end function b2ye_tprf
