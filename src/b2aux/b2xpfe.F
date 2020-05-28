      subroutine b2xpfe (nx, ny, ns, ns_ext, qe, 
     &  rza, za_ext, fna, fa_ext, fch, fne)
      use b2mod_types
      use b2mod_geo
      use b2mod_indirect
      implicit none
      integer nx, ny, ns, ns_ext
      real (kind=R8) ::
     *  qe, rza(-1:nx,-1:ny,0:ns-1), fna(-1:nx,-1:ny,0:1,0:ns-1),
     *  za_ext(-1:nx,-1:ny,0:max(ns_ext,1)-1), 
     *  fa_ext(-1:nx,-1:ny,0:1,0:max(ns_ext,1)-1),
     &  fch(-1:nx,-1:ny,0:1), fne(-1:nx,-1:ny,0:1)
*     ------------------------------------------------------------------
*     B2XPFE computes the electron particle flux, fne:
*       fne(,,) = (sum is :: rza(,,is)*fna(,,,is))-(1/qe)*fch(,,) .
*     A contribution from external species can also be included.
*     ------------------------------------------------------------------
      integer n2, ix, iy, is
      external sfill, b2saxpy
*     ------------------------------------------------------------------
      call subini ('b2xpfe')
      n2 = (nx+2)*(ny+2)
      call sfill (n2*2, 0.0_R8, fne, 1)
      do is = 0, ns-1
        do iy = -1, ny
          do ix = -1, nx
            if (leftix(ix,iy).ne.-2) 
     >        fne(ix,iy,0) = fne(ix,iy,0) + fna(ix,iy,0,is)*
     .         (rza(ix,iy,is)*vol(leftix(ix,iy),leftiy(ix,iy))+
     .          rza(leftix(ix,iy),leftiy(ix,iy),is)*vol(ix,iy))/
     .         (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))
            if (bottomiy(ix,iy).ne.-2) 
     >        fne(ix,iy,1) = fne(ix,iy,1) + fna(ix,iy,1,is)*
     .         (rza(ix,iy,is)*vol(bottomix(ix,iy),bottomiy(ix,iy))+
     .          rza(bottomix(ix,iy),bottomiy(ix,iy),is)*vol(ix,iy))/
     .         (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))
          enddo
        enddo
      enddo
      do is = 0, ns_ext-1
        do iy = -1, ny
          do ix = -1, nx
            if (leftix(ix,iy).ne.-2) 
     >        fne(ix,iy,0) = fne(ix,iy,0) + fa_ext(ix,iy,0,is)*
     .         (za_ext(ix,iy,is)*vol(leftix(ix,iy),leftiy(ix,iy))+
     .          za_ext(leftix(ix,iy),leftiy(ix,iy),is)*vol(ix,iy))/
     .         (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))
            if (bottomiy(ix,iy).ne.-2) 
     >        fne(ix,iy,1) = fne(ix,iy,1) + fa_ext(ix,iy,1,is)*
     .         (za_ext(ix,iy,is)*vol(bottomix(ix,iy),bottomiy(ix,iy))+
     .          za_ext(bottomix(ix,iy),bottomiy(ix,iy),is)*vol(ix,iy))/
     .         (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))
          enddo
        enddo
      enddo
      call b2saxpy (n2*2, -1.0_R8/qe, fch, 1, fne, 1)
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpfe
