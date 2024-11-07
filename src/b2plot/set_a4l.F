
      subroutine set_a4l
      use b2mod_types
      use b2mod_graphics
      implicit none

      xwidth=1.0_R4
      yheight=1.0_R4/sqrt(2.0_R4)
#ifdef PLPLOT
      call plsori(0)
#endif
      return
      end subroutine set_a4l
