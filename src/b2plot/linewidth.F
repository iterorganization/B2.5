c
      subroutine linewidth(width)
      use b2mod_types
      use b2mod_graphics
      implicit none
      real (kind=R8) :: width

      lwidth=real(width,R4)
      write (*,*) 'Line Width set to',width
      return
      end subroutine linewidth
c
      subroutine axislinewidth(width)
      use b2mod_types
      use b2mod_graphics
      implicit none
      real (kind=R8) :: width

      zwidth=real(width,R4)
      write (*,*) 'Axis Line Width set to',width
      return
      end subroutine axislinewidth
c
      subroutine walllinewidth(width)
      use b2mod_types
      use b2mod_graphics
      implicit none
      real (kind=R8) :: width

      wwidth=real(width,R4)
      write (*,*) 'Wall Line Width set to',width
      return
      end subroutine walllinewidth
