      subroutine b2xz12 (nCv, ns, ismain, rz2, na, am1, z12)            !som 13.11.20 {
      use b2mod_types
      implicit none
      integer nCv, ns, ismain
      real (kind=R8) ::
     *  rz2(nCv,0:ns-1), na(nCv,0:ns-1),
     *  z12(nCv), am1(0:ns-1)
*     ------------------------------------------------------------------
*     B2XZ12 computes the z12 quantity needed for the Grad-Zhdanov
*     closure scheme:
*       z12() = (sum over j /= i ::
*                mu(i,j)/m(i)*(m(j)/m(i)*kappa(i,j))*sqrt(mu(i,j)/m(i))*
*                Z_2_n_i(j)/Z_2_n_i(i)*(1+32/23*m(i)/m(j)+
*                                         84/23*(m(i)/m(j))^2))
*     ------------------------------------------------------------------
      integer is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xz12')
      call sfill (nCv, 0.0_R8, z12, 1)
      do is = 0, ns-1
       if(is.ne.ismain) then
         z12 = z12 +
     .       rz2(:,is)*na(:,is)*
     .       sqrt(am1(is)/(am1(ismain)+am1(is)))*
     .       (am1(is)/(am1(ismain)+am1(is)))**3*
     .       (1.0_R8+1.3913_R8*am1(ismain)/am1(is)+
     .       3.6522_R8*(am1(ismain)/am1(is))**2)
        endif
      enddo
      z12=z12/(rz2(:,ismain)*na(:,ismain))
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xz12                                             !som 13.11.20 }
