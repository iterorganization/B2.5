      subroutine b2xpif (nx, ny, hx, vol, fch_pi_c, fch_pi_f)      !srv 13.01.17
      use b2mod_types
      use b2mod_indirect
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *  hx(-1:nx,-1:ny),
     *  vol(-1:nx,-1:ny), fch_pi_c(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fch_pi_f(-1:nx,-1:ny,0:1)
*     ------------------------------------------------------------------
*     B2XPIF computes the ion parallel current on cell faces.
*     ------------------------------------------------------------------
      integer ix, iy
      integer, save :: ncall = 0
      real (kind=R8), save  :: prl_cur = 1.0_R8
*     ------------------------------------------------------------------
      if (ncall.eq.0) then
       call ipgetr ('b2tfhe_prl_cur', prl_cur)
      endif

      call intfaceh(nx,ny,leftix,leftiy,vol,fch_pi_c(-1,-1,0),
     .                                      fch_pi_f(-1,-1,0))

      do iy = -1, ny
        do ix = -1, nx
          if(leftix(ix,iy).ne.-2) then
            fch_pi_f(ix,iy,0) = fch_pi_f(ix,iy,0)*
     &       (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))/
     &       (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*prl_cur
          endif
        enddo
      enddo

      fch_pi_f(-1:nx,-1:ny,1) = 0.0_R8

      ncall = ncall+1
      return
*     ------------------------------------------------------------------
      end subroutine b2xpif
