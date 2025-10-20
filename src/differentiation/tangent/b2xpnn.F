      subroutine b2xpnn (nCv, ns, na, nn)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa_diffv
      implicit none
      integer nCv, ns
      real (kind=R8) ::
     *   na(nCv,0:ns-1), nn(nCv)
*   ..common blocks
*     ------------------------------------------------------------------
*     B2XPNN computes the total density of hydrogenic neutrals, nn:
*       nn = (sum is, is_neutral(is).and.zn(is).eq.1 :: na(,is)) .
*     ------------------------------------------------------------------
      integer is
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpnn')
      nn = 0.0_R8
      do is = 0, ns-1
       if(is_neutral(is).and.nint(zn(is)).eq.1) nn = nn + na(:,is)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpnn
