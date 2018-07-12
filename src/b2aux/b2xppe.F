      subroutine b2xppe (nCv, ne, te, pe)
      use b2mod_types
      implicit none
      integer nCv
      real (kind=R8) ::
     *  ne(nCv), te(nCv), pe(nCv)
*     ------------------------------------------------------------------
*     B2XPPE computes the electron pressure, pe = ne*te.
*     ------------------------------------------------------------------

*     ------------------------------------------------------------------
*$$$  call subini ('b2xppe')
      pe = ne*te
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xppe
