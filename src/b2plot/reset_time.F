c
      subroutine reset_time(time)
      use b2mod_types
      use b2mod_graphics
      implicit none
      real (kind=R8) :: time

      runtime=real(time,R4)
      return
      end subroutine reset_time
