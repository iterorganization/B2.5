      subroutine b2xzf1 (nCv, ns, ismain, rz2, na, am1, zf1)            !som 13.11.20 {
      use b2mod_types
      implicit none
      integer nCv, ns, ismain
      real (kind=R8) ::
     *  rz2(nCv,0:ns-1), na(nCv,0:ns-1),
     *  zf1(nCv), am1(0:ns-1)
*     ------------------------------------------------------------------
*     B2XZF1 computes the zf1 quantity needed for the Grad-Zhdanov
*     closure scheme:
*       zf1() = (sum over j /= i :: mu(i,j)/m(i)*sqrt(mu(i,j)/m(i))*
*                                   Z_2_n_i(j)/Z_2_n_i(i))
*     ------------------------------------------------------------------
      integer is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xzf1')
      call sfill (nCv, 0.0_R8, zf1, 1)
      do is = 0, ns-1
       if(is.ne.ismain) then
         zf1 = zf1 +
     .       rz2(:,is)*na(:,is)*
     .       sqrt(am1(is)/(am1(ismain)+am1(is)))*
     .       (am1(is)/(am1(ismain)+am1(is)))
        endif
      enddo
      zf1=zf1/(rz2(:,ismain)*na(:,ismain))
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xzf1                                             !som 13.11.20 }
