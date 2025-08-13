      subroutine b2xz11 (nCv, ns, ismain, rz2, na, am1, z11)            !som 13.11.20 {
      use b2mod_types
      implicit none
      integer nCv, ns, ismain
      real (kind=R8) ::
     *  rz2(nCv,0:ns-1), na(nCv,0:ns-1),
     *  z11(nCv), am1(0:ns-1)
*     ------------------------------------------------------------------
*     B2XZ11 computes the z11 quantity needed for the Grad-Zhdanov
*     closure scheme:
*       z11() = (sum over j /= i ::
*               (m(j)/m(i)*kappa(i,j))*sqrt(mu(i,j)/m(i))*
*               Z_2_n_i(j)/Z_2_n_i(i)*(1+16/13*m(i)/m(j)+
*                                        30/13*(m(i)/m(j))^2))
*     ------------------------------------------------------------------
      integer is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xz11')
      call sfill (nCv, 0.0_R8, z11, 1)
      do is = 0, ns-1
       if(is.ne.ismain) then
         z11 = z11 +
     .       rz2(:,is)*na(:,is)*
     .       sqrt(am1(is)/(am1(ismain)+am1(is)))*
     .       (am1(is)/(am1(ismain)+am1(is)))**2*
     .       (1.0_R8+1.2308_R8*am1(ismain)/am1(is)+
     .       2.3077_R8*(am1(ismain)/am1(is))**2)
        endif
      enddo
      z11=z11/(rz2(:,ismain)*na(:,ismain))
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xz11                                             !som 13.11.20 }
