      subroutine b2xz22 (nx, ny, ns, ismain, rz2, na, am1, z22)             !som 13.11.20 {
      use b2mod_types
      implicit none
      integer nx, ny, ns, ismain
      real (kind=R8) ::
     *  rz2(-1:nx,-1:ny,0:ns-1), na(-1:nx,-1:ny,0:ns-1),
     *  z22(-1:nx,-1:ny), am1(0:ns-1)
*     ------------------------------------------------------------------
*     B2XZ22 computes the z22 quantity needed for the Grad-Zhdanov
*     closure scheme:
*       z22(,) = (sum over j /= i ::
*                (m(j)/m(i)*kappa(i,j))^2*sqrt(mu(i,j)/m(i))*
*                 Z_2_n_i(j)/Z_2_n_i(i)*(1+1088/433*m(i)/m(j)+
*                 3672/433*(m(i)/m(j))^2+1792/433*(m(i)/m(j))^3+
*                 1400/433*(m(i)/m(j))^4))
*     ------------------------------------------------------------------
      integer n2, is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xz22')
      n2 = (nx+2)*(ny+2)
      call sfill (n2, 0.0_R8, z22, 1)
      do is = 0, ns-1
       if(is.ne.ismain) then
         z22(-1:nx,-1:ny) = z22(-1:nx,-1:ny) +
     .       rz2(-1:nx,-1:ny,is)*na(-1:nx,-1:ny,is)*
     .       sqrt(am1(is)/(am1(ismain)+am1(is)))*
     .       (am1(is)/(am1(ismain)+am1(is)))**4*
     .       (1.0_R8+2.5127_R8*am1(ismain)/am1(is)+
     .       8.4804_R8*(am1(ismain)/am1(is))**2+
     .       4.1386_R8*(am1(ismain)/am1(is))**3+
     .       3.2333_R8*(am1(ismain)/am1(is))**4)
        endif
      enddo
      z22=z22/(rz2(-1:nx,-1:ny,ismain)*na(-1:nx,-1:ny,ismain))
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xz22                                             !som 13.11.20 }
