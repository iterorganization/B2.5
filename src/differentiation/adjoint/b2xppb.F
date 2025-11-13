! $AD NOCHECKPOINT
      subroutine b2xppb (nCv, rzb, nb, te, ti, tn, is, pb)
      use b2mod_types
      use b2mod_b2cmpa_diff
      implicit none
      integer nCv, is
      real (kind=R8) ::
     *  rzb(nCv), nb(nCv), te(nCv), ti(nCv), tn(nCv), pb(nCv)
*     ------------------------------------------------------------------
*     B2XPPB computes the partial pressure, pb = nb*(rzb*te+ti).
*     ------------------------------------------------------------------
      
*     ------------------------------------------------------------------
*$$$  call subini ('b2xppb')
      if (.not.is_neutral(is)) then
        pb = nb*(rzb*te+ti)
      else
        pb = nb*tn
      endif
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xppb
