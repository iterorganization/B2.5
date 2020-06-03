      subroutine b2xvfv (nx, ny, funv, vrname)
      use b2mod_types
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *   funv(-1:nx,-1:ny,0:1)
      character vrname*(*)
*     -----------------------------------------------------------------
*     B2XVFV verifies the zero entries in the edge locations of some
*     quantity that is defined in the manner of a viscosity on the
*     cells and corners of the B2 mesh.
*     The entries funv(-1,-1:ny,0), funv(nx,-1:ny,0), funv(-1:nx,-1,1)
*     and funv(-1,-1:ny,1) must vanish or the code will halt via a
*     call to xertst/xerrab.
*     vrname must be the name of the program variable in the calling
*     routine that is associated with funv.
*     -----------------------------------------------------------------
      character vrnam0*8
      real (kind=R8) ::
     *   samax
      external xertst, samax
*     -----------------------------------------------------------------
*$$$  call subini ('b2xvfv')
      vrnam0 = vrname
      call xertst (samax(ny+2,funv(-1,-1,0),nx+2).eq.0,
     &  'b2xvfv--faulty argument edge (-1,,0) '//vrnam0)
      call xertst (samax(ny+2,funv(nx,-1,0),nx+2).eq.0,
     &  'b2xvfv--faulty argument edge (nx,,0) '//vrnam0)
      call xertst (samax(ny+2,funv(-1,-1,1),nx+2).eq.0,
     &  'b2xvfv--faulty argument edge (-1,,1) '//vrnam0)
      call xertst (samax(nx+2,funv(-1,-1,1),1).eq.0,
     &  'b2xvfv--faulty argument edge (,-1,1) '//vrnam0)
*$$$  call subend ()
      return
*     -----------------------------------------------------------------
      end
