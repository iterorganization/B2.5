      subroutine b2xzmc (nx, ny, ns, ismain, rz2, na, am1, zmc)             !som 13.11.20 {
      use b2mod_types
      implicit none
      integer nx, ny, ns, ismain
      real (kind=R8) ::
     *  rz2(-1:nx,-1:ny,0:ns-1), na(-1:nx,-1:ny,0:ns-1),
     *  zmc(-1:nx,-1:ny), am1(0:ns-1)
*     ------------------------------------------------------------------
*     B2XZMC computes the zmc quantity needed for the Grad-Zhdanov
*     closure scheme:
*       zmc(,) = (sum over j /= i :: sqrt(mu(i,j)/m(i))*
*                                    Z_2_n_i(j)/Z_2_n_i(i))
*     ------------------------------------------------------------------
      integer n2, is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xzmc')
      n2 = (nx+2)*(ny+2)
      call sfill (n2, 0.0_R8, zmc, 1)
      do is = 0, ns-1
       if(is.ne.ismain) then
         zmc(-1:nx,-1:ny) = zmc(-1:nx,-1:ny) +
     .       rz2(-1:nx,-1:ny,is)*na(-1:nx,-1:ny,is)*
     .       sqrt(am1(is)/(am1(ismain)+am1(is)))
        endif
      enddo
      zmc=zmc/(rz2(-1:nx,-1:ny,ismain)*na(-1:nx,-1:ny,ismain))
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xzmc                                             !som 13.11.20 }
