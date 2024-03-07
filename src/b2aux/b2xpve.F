      subroutine b2xpve (nx, ny, ns, qe, rza, rz2, alfx_c, sigx_c,      !srv 17.11.16 {
     &                   po, ne, te, na, ua, ue_ext, ne_ext, ve)
      use b2mod_types
      use b2mod_geo
      use b2mod_indirect
      use b2mod_plasma
     &    , only : fch_pi_c                                             !srv 13.01.17
      implicit none
      integer nx, ny, ns
      real (kind=R8) ::
     &   qe, rza(-1:nx,-1:ny,0:ns-1), rz2(-1:nx,-1:ny,0:ns-1),
     &   alfx_c(-1:nx,-1:ny), sigx_c(-1:nx,-1:ny),
     &   te(-1:nx,-1:ny), po(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     &   na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &   ue_ext(-1:nx,-1:ny), ne_ext(-1:nx,-1:ny)
      real (kind=R8) ::
     *   ve(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XPVE computes the parallel electron velocity on cell centers.
*     ------------------------------------------------------------------
      integer ix, iy, is
      integer, save :: ncall = 0                                        !srv 03.04.15
      real (kind=R8), save :: prl_cur = 1.0_R8, phm2 = 1.0_R8           !srv 03.04.15
      integer, save :: style = 2                                        !srv 13.01.17
      real (kind=R8) ::
     &   wrksig(-1:nx,-1:ny), wrkalf(-1:nx,-1:ny), zeff(-1:nx,-1:ny)
      external sfill
*     ------------------------------------------------------------------
      if (ncall.eq.0) then
       call ipgetr ('b2tfhe_prl_cur', prl_cur)
       call ipgetr ('b2sifr_phm2', phm2)
       call ipgeti ('b2sigp_style', style)                              !srv 13.01.17
      endif

      if(style.eq.2) then                                               !srv 13.01.17 {
       call b2xzef (nx, ny, ns, rz2, na, ne, zeff)
       call b2xpic (nx, ny, ns, qe, rza, bb, na, ua, zeff, fch_pi_c)
      else
       fch_pi_c = 0.0_R8
      endif                                                             !srv 13.01.17 }
*   ..compute ve, electron parallel velocity                            !iyv 19.11.13 {
      do iy=-1,ny
        do ix=-1,nx
          if(region(ix,iy,0).eq.0.or.leftix(ix,iy).eq.-2) then
            wrksig(ix,iy) = 0.0_R8
            wrkalf(ix,iy) = 0.0_R8
          else
            wrksig(ix,iy) = (
     &          -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))+
     &           (ne(ix,iy)*te(ix,iy)
     &           -ne(leftix(ix,iy),leftiy(ix,iy))*
     &            te(leftix(ix,iy),leftiy(ix,iy)))*
     &    (2.0_R8/qe)/(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy)) )
     &    /(hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8

            wrkalf(ix,iy) =
     &         (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))
     &        /(hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
          endif
        enddo
      enddo
cxpb Special treatment to remove spurious values in top right corners
      do iy=-1,ny
        do ix=-1,nx
          if(topiy(ix,iy).eq.ny+1.and.rightix(ix,iy).eq.nx+1) then
            wrksig(ix,iy) = wrksig(leftix(ix,iy),leftiy(ix,iy))
            wrkalf(ix,iy) = wrkalf(leftix(ix,iy),leftiy(ix,iy))
          endif
        enddo
      enddo

      do iy=-1,ny
        do ix=-1,nx
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            ve(ix,iy)=-( sigx_c(ix,iy)*
     &         (wrksig(ix,iy)+wrksig(rightix(ix,iy),rightiy(ix,iy)))-
     &                       phm2*alfx_c(ix,iy)*
     &         (wrkalf(ix,iy)+wrkalf(rightix(ix,iy),rightiy(ix,iy))))
     &          /2.0_R8
          else if(leftix(ix,iy).eq.-2) then
            ve(ix,iy)=-( sigx_c(ix,iy)*
     &         (wrksig(rightix(ix,iy),rightiy(ix,iy)))-
     &                       phm2*alfx_c(ix,iy)*
     &         (wrkalf(rightix(ix,iy),rightiy(ix,iy))))
          else if(rightix(ix,iy).eq.nx+1) then
            ve(ix,iy)=-( sigx_c(ix,iy)*wrksig(ix,iy)
     &          -phm2*alfx_c(ix,iy)*wrkalf(ix,iy))
          endif
        enddo
      enddo
      ve(:,:) = (ve(:,:)+fch_pi_c(:,:,0,0))/(bb(:,:,0)/bb(:,:,3))/qe    !srv 13.01.17
     &                                                    *prl_cur      !srv 21.12.15
      do is=0,ns-1
        ve(:,:)=ve(:,:)+rza(:,:,is)*ua(:,:,is)*na(:,:,is)
      enddo
      ve(:,:)=ve(:,:)+ue_ext(:,:)*ne_ext(:,:)
      ve(:,:)=ve(:,:)/ne(:,:)

      ncall = ncall+1
      return
*     ------------------------------------------------------------------
      end subroutine b2xpve                                             !srv 17.11.16 }
