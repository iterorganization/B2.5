c     ======================================================================
c     SUBROUTINES FOR INTERACTION OF B2PLOT AND B2PSH
c     ======================================================================

c     ----------------------------------------------------------------------
c     GETTING THE PLOTAREA
c     ----------------------------------------------------------------------

      subroutine getxyrz
      use b2mod_b2plot_values
      implicit none

      write(*,*) 'b2psh: xmin ',values(1)
      write(*,*) 'b2psh: xmax ',values(2)
      write(*,*) 'b2psh: ymin ',values(3)
      write(*,*) 'b2psh: ymax ',values(4)
      write(*,*) 'b2psh: rmin ',values(5)
      write(*,*) 'b2psh: rmax ',values(6)
      write(*,*) 'b2psh: zmin ',values(7)
      write(*,*) 'b2psh: zmax ',values(8)
      write(*,*) 'b2psh: read values!'

      return
      end subroutine getxyrz
      
c     ----------------------------------------------------------------------
c     SETTING THE WINDOW ID
c     ----------------------------------------------------------------------

      subroutine setwinid(n)
      use b2mod_graphics
      implicit none

c     INPUT
      integer n

      winid=n
      write(*,*) 'Setting WindowID to: ',winid

      return
      end subroutine setwinid

c     ----------------------------------------------------------------------
c     RESETTING OUTPUT FILE
c     ----------------------------------------------------------------------

      subroutine resetoutput
      use b2mod_types
      implicit none
      call terminate_plotting
      call initialise_plotting_flags(0.0_R8)
      write(*,*) 'OUTPUT FILE RESET'

      return
      end subroutine resetoutput
