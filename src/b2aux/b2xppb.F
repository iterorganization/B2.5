      subroutine b2xppb (nCv, rzb, nb, te, ti, pb)
      use b2mod_types
      implicit none
      integer nCv
      real (kind=R8) ::
     *  rzb(nCv), nb(nCv), te(nCv), ti(nCv), pb(nCv)
*     ------------------------------------------------------------------
*     B2XPPB computes the partial pressure, pb = nb*(rzb*te+ti).
*     ------------------------------------------------------------------
      
*     ------------------------------------------------------------------
*$$$  call subini ('b2xppb')
      pb = nb*(rzb*te+ti)
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end
