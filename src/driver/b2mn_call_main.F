      subroutine b2mn_call_main(J)
      use b2mod_main
      use b2mod_types
      use b2mod_b2cmpt
      use b2mod_boundary_namelist
      implicit none
      real(kind=R8) :: J


      call b2mn_step
      call b2mn_pp(J)


      end subroutine b2mn_call_main

