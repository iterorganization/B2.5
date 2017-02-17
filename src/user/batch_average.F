
      subroutine batch_average (n,var,aver,it,batch_size)

      use b2mod_types

      implicit none

      integer :: n, it, batch_size
      real (kind=R8) :: var(1:n), aver(1:n), gnew, gold
      integer :: fac

      fac = mod(it,batch_size)
      if (fac.eq.0) fac = batch_size
      if (fac.eq.1) then
        aver = var
      else
        gnew = 1.0_R8 / fac
        gold = (fac - 1) * gnew
        aver = gold*aver + gnew*var
      endif

      return
      end subroutine batch_average
