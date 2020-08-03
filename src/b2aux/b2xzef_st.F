      subroutine b2xzef_st (nx, ny, ns, rz2, na, ne, zeff)             !srv 13.01.17 {
      use b2mod_types
      use b2mod_external
      implicit none
      integer nx, ny, ns
      real (kind=R8) ::
     *  rz2(-1:nx,-1:ny,0:ns-1), na(-1:nx,-1:ny,0:ns-1), 
     *  ne(-1:nx,-1:ny), zeff(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XZEF_ST computes zeff:
*       zeff(,) = (sum is :: rz2(,,is)*na(,,is))/ne .
*     The routines only works on structured grids.
*     ------------------------------------------------------------------
      integer n2, is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xzef_st')
      n2 = (nx+2)*(ny+2)
      call sfill (n2, 0.0_R8, zeff, 1)
      do is = 0, ns-1
       zeff(-1:nx,-1:ny) = zeff(-1:nx,-1:ny) +
     .  rz2(-1:nx,-1:ny,is)*na(-1:nx,-1:ny,is)
      enddo
      do is = 0, ns_ext-1
       zeff(-1:nx,-1:ny) = zeff(-1:nx,-1:ny) +
     .  za2_ext(-1:nx,-1:ny,is)*na_ext(-1:nx,-1:ny,is)
      enddo
      zeff=zeff/ne
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xzef_st                                         !srv 13.01.17 }
