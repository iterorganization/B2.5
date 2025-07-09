      subroutine b2xzmc (nCv, ns, ismain, rz2, na, am1, zmc)            !som 13.11.20 {
      use b2mod_types
      implicit none
      integer nCv, ns, ismain
      real (kind=R8) ::
     *  rz2(nCv,0:ns-1), na(nCv,0:ns-1),
     *  zmc(nCv), am1(0:ns-1)
*     ------------------------------------------------------------------
*     B2XZMC computes the zmc quantity needed for the Grad-Zhdanov
*     closure scheme:
*       zmc() = (sum over j /= i :: sqrt(mu(i,j)/m(i))*
*                                    Z_2_n_i(j)/Z_2_n_i(i))
*     ------------------------------------------------------------------
      integer is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xzmc')
      call sfill (nCv, 0.0_R8, zmc, 1)
      do is = 0, ns-1
       if(is.ne.ismain) then
         zmc = zmc +
     .       rz2(:,is)*na(:,is)*
     .       sqrt(am1(is)/(am1(ismain)+am1(is)))
        endif
      enddo
      zmc=zmc/(rz2(:,ismain)*na(:,ismain))
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xzmc                                             !som 13.11.20 }
