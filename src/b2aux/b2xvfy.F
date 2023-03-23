      subroutine b2xvfy (nx, ny, funy, vrname)
      use b2mod_types
      use b2mod_indirect
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *   funy(-1:nx,-1:ny)
      character vrname*(*)
*     -----------------------------------------------------------------
*     B2XVFY verifies the zero entries in the edge locations of some
*     quantity defined on y-faces of the B2 mesh.
*     The entries funy(-1:nx,-1) must vanish or the code will halt via
*     a call to xertst/xerrab.
*     vrname must be the name of the program variable in the calling
*     routine that is associated with funy.
*     -----------------------------------------------------------------
      character vrnam0*8
      integer ix
      external xertst
*     -----------------------------------------------------------------
*$$$  call subini ('b2xvfy')
      vrnam0 = vrname
      if(allocated(bottomiy)) then
        do ix = -1, nx
          call xertst(bottomiy(ix,-1).ne.-2 .or. funy(ix,-1).eq.0,
     &     'b2xvfy--faulty argument y-edge '//vrnam0)
        enddo
      else
        do ix = -1, nx
          call xertst(funy(ix,-1).eq.0,
     &     'b2xvfy--faulty argument y-edge '//vrnam0)
        enddo
      endif
*$$$  call subend ()
      return
*     -----------------------------------------------------------------
      end subroutine b2xvfy

*     -----------------------------------------------------------------
      subroutine b2xvfy1(nx, ny, funy, vrname, bottomiy1)
      use b2mod_types
      implicit none
      integer nx, ny, bottomiy1(-1:nx,-1:ny)
      real (kind=R8) ::
     *   funy(-1:nx,-1:ny)
      character vrname*(*)
*     -----------------------------------------------------------------
*     B2XVFY1 verifies the zero entries in the edge locations of some
*     quantity defined on y-faces of the B2 coarse mesh.
*     The entries funy(-1:nx,-1) must vanish or the code will halt via
*     a call to xertst/xerrab.
*     vrname must be the name of the program variable in the calling
*     routine that is associated with funy.
*     -----------------------------------------------------------------
      character vrnam0*8
      integer ix
      external xertst
*     -----------------------------------------------------------------
*$$$  call subini ('b2xvfy1')
      vrnam0 = vrname
      do ix = -1, nx
        call xertst(bottomiy1(ix,-1).ne.-2 .or. funy(ix,-1).eq.0,
     &  'b2xvfy1--faulty argument y-edge '//vrnam0)
      enddo
*$$$  call subend ()
      return
*     -----------------------------------------------------------------
      end subroutine b2xvfy1
