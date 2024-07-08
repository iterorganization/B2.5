      subroutine b2xzf2 (nCv, ns, ismain, rz2, na, am1, zf2)            !som 13.11.20 {
      use b2mod_types
      implicit none
      integer nCv, ns, ismain
      real (kind=R8) ::
     *  rz2(nCv,0:ns-1), na(nCv,0:ns-1),
     *  zf2(nCv), am1(0:ns-1)
*     ------------------------------------------------------------------
*     B2XZF2 computes the zf2 quantity needed for the Grad-Zhdanov
*     closure scheme:
*       zf2() = (sum over j /= i :: (mu(i,j)/m(i))^2*sqrt(mu(i,j)/m(i))*
*                                    Z_2_n_i(j)/Z_2_n_i(i))
*     ------------------------------------------------------------------
      integer is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xzf2')
      call sfill (nCv, 0.0_R8, zf2, 1)
      do is = 0, ns-1
       if(is.ne.ismain) then
         zf2 = zf2 +
     .       rz2(:,is)*na(:,is)*
     .       sqrt(am1(is)/(am1(ismain)+am1(is)))*
     .       (am1(is)/(am1(ismain)+am1(is)))**2
        endif
      enddo
      zf2=zf2/(rz2(:,ismain)*na(:,ismain))
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xzf2                                             !som 13.11.20 }
