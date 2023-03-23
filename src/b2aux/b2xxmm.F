      subroutine b2xxmm (nx, ny, fun, iloc, ix0, iy0, ix1, iy1)
      use b2mod_types
      implicit none
      integer nx, ny, iloc, ix0, iy0, ix1, iy1
      real (kind=R8) ::
     *   fun(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XXMM computes the locations of the minimum and maximum entries
*     of a quantity defined on the B2 mesh.
*     ------------------------------------------------------------------
      integer kxb, kyb, kx0, ky0, kx1, ky1, kx, ky
      external xertst
*     ------------------------------------------------------------------
*$$$  call subini ('b2xxmm')
      call xertst (0.le.nx.and.0.le.ny,
     &  'b2xxmm--faulty argument nx, ny')
      call xertst (0.le.iloc.and.iloc.lt.4,
     &  'b2xxmm--faulty argument iloc')
*     If iloc.eq.0 then fun is defined on (-1:nx,-1:ny) (centers).
*     If iloc.eq.1 then fun is defined on ( 0:nx,-1:ny) (x-faces).
*     If iloc.eq.2 then fun is defined on (-1:nx, 0:ny) (y-faces).
*     If iloc.eq.3 then fun is defined on ( 0:nx, 0:ny) (corners).
      kxb = mod(iloc,2)-1
      kyb = iloc/2-1
      kx0 = kxb
      ky0 = kyb
      kx1 = kxb
      ky1 = kyb
      do ky = kyb, ny
       do kx = kxb, nx
        if (fun(kx,ky).lt.fun(kx0,ky0)) then
         kx0 = kx
         ky0 = ky
        else if (fun(kx1,ky1).le.fun(kx,ky)) then
         kx1 = kx
         ky1 = ky
        endif
       enddo
      enddo
      ix0 = kx0
      iy0 = ky0
      ix1 = kx1
      iy1 = ky1
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xxmm
