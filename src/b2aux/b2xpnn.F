      subroutine b2xpnn (nCv, ns, ismain0, na, nn)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer nCv, ns, ismain0
      real (kind=R8) ::
     *   na(nCv,0:ns-1), nn(nCv)
*   ..common blocks
*     ------------------------------------------------------------------
*     B2XPNN computes the total density of main hydrogenic neutrals, nn
*     ------------------------------------------------------------------
      integer is
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpnn')
      nn = 0.0_R8
      if (is_neutral(ismain0).and.nint(zn(ismain0)).eq.1)
     > nn = na(:,ismain0)
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpnn
