      subroutine b2xpic (nx, ny, ns, qe, rza, bb, na, ua, zeff,       !srv 13.01.17
     &                   fch_pi_c)
      use b2mod_types
      use b2mod_external
      implicit none
      integer nx, ny, ns
      real (kind=R8) ::
     *  qe, rza(-1:nx,-1:ny,0:ns-1), 
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     *  bb(-1:nx,-1:ny,0:3), zeff(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fch_pi_c(-1:nx,-1:ny,0:1)
*     ------------------------------------------------------------------
*     B2XPIC computes the poloidal projection of the ion parallel
*     current on cell centers.
*     Attention! There is no multiplication by vol/hx. 
*     It will be done later in B2XPIF.
*     ------------------------------------------------------------------
      integer n2, is
      integer, save :: ncall = 0
      external sfill
*     ------------------------------------------------------------------
*   ..set n2
      n2 = (nx+2)*(ny+2)
      fch_pi_c(-1:nx,-1:ny,0) = 0.0_R8
      do is = 0, ns-1
       fch_pi_c(-1:nx,-1:ny,0) = fch_pi_c(-1:nx,-1:ny,0) +
     &  rza(-1:nx,-1:ny,is)*(rza(-1:nx,-1:ny,is)-zeff(-1:nx,-1:ny))*
     &  na(-1:nx,-1:ny,is)*ua(-1:nx,-1:ny,is)
      enddo
      do is = 0, ns_ext-1
       fch_pi_c(-1:nx,-1:ny,0) = fch_pi_c(-1:nx,-1:ny,0) +
     &  za_ext(-1:nx,-1:ny,is)*
     & (za_ext(-1:nx,-1:ny,is)-zeff(-1:nx,-1:ny))*
     &  na_ext(-1:nx,-1:ny,is)*ua_ext(-1:nx,-1:ny,is)
      enddo
      fch_pi_c(-1:nx,-1:ny,0) = fch_pi_c(-1:nx,-1:ny,0)*qe*
     & (bb(-1:nx,-1:ny,0)/bb(-1:nx,-1:ny,3))/zeff(-1:nx,-1:ny)

      fch_pi_c(-1:nx,-1:ny,1) = 0.0_R8

      ncall = ncall+1
      return
*     ------------------------------------------------------------------
      end subroutine b2xpic
