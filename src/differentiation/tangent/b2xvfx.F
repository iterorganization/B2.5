      subroutine b2xvfx (nx, ny, funx, vrname)
      use b2mod_types
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *   funx(-1:nx,-1:ny)
      character vrname*(*)
*     -----------------------------------------------------------------
*     B2XVFX verifies the zero entries in the edge locations of some
*     quantity defined on x-faces of the B2 mesh.
*     The entries funx(-1,-1:ny) must vanish or the code will halt via
*     a call to xertst/xerrab.
*     vrname must be the name of the program variable in the calling
*     routine that is associated with funx.
*     -----------------------------------------------------------------
      character vrnam0*8
      real (kind=R8) ::
     *   damax
      external xertst, damax
*     -----------------------------------------------------------------
*$$$  call subini ('b2xvfx')
      vrnam0 = vrname
      call xertst (damax(ny+2,funx(-1,-1),nx+2).eq.0,
     &  'b2xvfx--faulty argument x-edge '//vrnam0)
*$$$  call subend ()
      return
*     -----------------------------------------------------------------
      end subroutine b2xvfx
