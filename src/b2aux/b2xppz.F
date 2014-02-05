      subroutine b2xppz (nx, ny, ns, rza, na, te, ti, pz)         !sv 18.01.02 {
      use b2mod_types
      use b2mod_b2cmpa
      implicit none
      integer nx, ny, ns
      real (kind=R8) ::
     *  rza(-1:nx,-1:ny,0:ns-1), na(-1:nx,-1:ny,0:ns-1), 
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), pz(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XPPZ computes pressure of charged particles, pz:
*       pz(,) = (sum is :: na(,,is)*(za(is)*te(,)+ti(,))) for za(is) .ne. 0.
*     ------------------------------------------------------------------
      integer n2, ix, iy, is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xppz')
      n2 = (nx+2)*(ny+2)
      call sfill (n2, 0.0_R8, pz, 1)
      do is = 0, ns-1
	 if (.not.is_neutral(is)) then
        do iy = -1, ny
         do ix = -1, nx
          pz(ix,iy) = pz(ix,iy)+
     &     na(ix,iy,is)*(rza(ix,iy,is)*te(ix,iy)+ti(ix,iy))
         enddo
        enddo
	 endif
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end                                                        !sv 18.01.02 }
