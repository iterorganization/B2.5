      subroutine b2xppz_st (nx, ny, ns, ne, na, te, ti, pz)     !sv 18.01.02 {
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_external
      implicit none
      integer nx, ny, ns
      real (kind=R8) ::
     *  ne(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1), 
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), pz(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XPPZ_ST computes pressure of charged particles, pz:
*       pz(,) = ne(,)*te(,) + (sum is :: na(,,is)*ti(,) for za(is).ne.0)
*     The routines only works on structured grids.
*     ------------------------------------------------------------------
      integer n2, ix, iy, is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xppz_st')
      n2 = (nx+2)*(ny+2)
      call sfill (n2, 0.0_R8, pz, 1)
      pz = ne*te
      do is = 0, ns-1
        if (.not.is_neutral(is)) then
          do iy = -1, ny
            do ix = -1, nx
              pz(ix,iy) = pz(ix,iy)+na(ix,iy,is)*ti(ix,iy)
            enddo
          enddo
        endif
      enddo
      do is = 0, ns_ext-1
        if (.not.is_ext_neutral(is)) then
          do iy = -1, ny
            do ix = -1, nx
              pz(ix,iy) = pz(ix,iy)+na_ext(ix,iy,is)*ta_ext(ix,iy,is)
            enddo
          enddo
        endif
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end                                                        !sv 18.01.02 }
