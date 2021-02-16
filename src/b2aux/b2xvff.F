      subroutine b2xvff (nx, ny, funf, vrname)
      use b2mod_types
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *   funf(-1:nx,-1:ny,0:1)
      character vrname*(*)
*     -----------------------------------------------------------------
*     B2XVFF verifies the zero entries in the edge locations of some
*     quantity defined on faces of the B2 mesh. This routine is an
*     interface to b2xvfx and b2xvfy.
*     funf(,,0) is defined on x-faces and funf(,,1) on y-faces. The
*     entries funf(-1,-1:ny,0) and funf(-1:nx,-1,1) must vanish or the
*     code will halt via a call to xertst/xerrab.
*     vrname must be the name of the program variable in the calling
*     routine that is associated with funf.
*     -----------------------------------------------------------------
      external b2xvfx, b2xvfy
*     -----------------------------------------------------------------
*$$$  call subini ('b2xvff')
      call b2xvfx (nx, ny, funf(-1,-1,0), vrname)
      call b2xvfy (nx, ny, funf(-1,-1,1), vrname)
*$$$  call subend ()
      return
*     -----------------------------------------------------------------
      end subroutine b2xvff
