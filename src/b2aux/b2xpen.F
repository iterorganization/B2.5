      subroutine b2xpen (nx, ny, ns, fna,
     1 fna_53, fne_53, fch, fhe, fhi, kinrgy, ua, cvsa, rpt, te, ti, po,
     2 BoRiS, fhm, fhp, fhj, fht)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      implicit none
      integer nx, ny, ns
! inputs
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     *  fna_53(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     1  fne_53(-1:nx,-1:ny,0:1,0:1),
     2  fch(-1:nx,-1:ny,0:1,0:1),
     2  fhe(-1:nx,-1:ny,0:1,0:1), fhi(-1:nx,-1:ny,0:1,0:1),
     3  kinrgy(-1:nx,-1:ny,0:ns-1), rpt(-1:nx,-1:ny,0:ns-1),
     4  ua(-1:nx,-1:ny,0:ns-1), cvsa(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     5  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny), BoRiS
! outputs
      real (kind=R8) ::
     &  fhm(-1:nx,-1:ny,0:1,0:1,0:ns-1),fhp(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     1  fhj(-1:nx,-1:ny,0:1,0:1),       fht(-1:nx,-1:ny,0:1,0:1)
*     ------------------------------------------------------------------
*     B2XPEN computes contributions to the total energy fluxes,
*     as well as the total energy fluxes
*     The contributions from external species are not included and must
*     be added explicitly in the calling routine when needed.
*     ------------------------------------------------------------------
      integer ix, iy, is
      real (kind=R8) ::
     &  fnt(-1:nx,-1:ny,0:1,0:1)
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpen')
      fhm = 0.0_R8
      fhp = 0.0_R8
      fhj = 0.0_R8
      fnt = 0.0_R8
      fht = 0.0_R8
c calculate additional fluxes (fhj, fhm, fhp)
      do iy=-1,ny
        do ix=-1,nx
          if(leftix(ix,iy).ge.-1) then
            fhj(ix,iy,0,0)=-fch(ix,iy,0,0)*
     1       (po(ix,iy)+po(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8
            fhj(ix,iy,0,1)=-fch(ix,iy,0,1)*
     1       (po(ix,iy)+po(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8
            fnt(ix,iy,0,0)=fne_53(ix,iy,0,0)*
     1       (te(ix,iy)+te(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8
            fnt(ix,iy,0,1)=fne_53(ix,iy,0,1)*
     1       (te(ix,iy)+te(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8
          endif
          if(bottomiy(ix,iy).ge.-1) then
            fhj(ix,iy,1,0)=-fch(ix,iy,1,0)*
     1       (po(ix,iy)+po(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8
            fhj(ix,iy,1,1)=-fch(ix,iy,1,1)*
     1       (po(ix,iy)+po(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8
            fnt(ix,iy,1,0)=fne_53(ix,iy,1,0)*
     1       (te(ix,iy)+te(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8
            fnt(ix,iy,1,1)=fne_53(ix,iy,1,1)*
     1       (te(ix,iy)+te(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8
          endif
        enddo
      enddo
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            if(leftix(ix,iy).ge.-1) then
              fhm(ix,iy,0,0,is)=0.5_R8*fna(ix,iy,0,0,is)*
     2         (kinrgy(leftix(ix,iy),leftiy(ix,iy),is)+kinrgy(ix,iy,is))
     3        -0.5_R8*cvsa(ix,iy,0,0,is)*
     4         (ua(ix,iy,is)+ua(leftix(ix,iy),leftiy(ix,iy),is))*
     5         (ua(ix,iy,is)-ua(leftix(ix,iy),leftiy(ix,iy),is))
              fhm(ix,iy,0,1,is)=0.5_R8*fna(ix,iy,0,1,is)*
     2         (kinrgy(leftix(ix,iy),leftiy(ix,iy),is)+kinrgy(ix,iy,is))
     3        -0.5_R8*cvsa(ix,iy,0,1,is)*
     4         (ua(ix,iy,is)+ua(leftix(ix,iy),leftiy(ix,iy),is))*
     5         (ua(ix,iy,is)-ua(leftix(ix,iy),leftiy(ix,iy),is))
              fhp(ix,iy,0,0,is)=0.5_R8*fna(ix,iy,0,0,is)*
     2         (rpt(leftix(ix,iy),leftiy(ix,iy),is)+rpt(ix,iy,is))*ev
              fhp(ix,iy,0,1,is)=0.5_R8*fna(ix,iy,0,1,is)*
     2         (rpt(leftix(ix,iy),leftiy(ix,iy),is)+rpt(ix,iy,is))*ev
              fnt(ix,iy,0,0)=fnt(ix,iy,0,0)+0.5_R8*fna_53(ix,iy,0,0,is)*
     2         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))
              fnt(ix,iy,0,1)=fnt(ix,iy,0,1)+0.5_R8*fna_53(ix,iy,0,1,is)*
     2         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))
            endif
            if(bottomiy(ix,iy).ge.-1) then
              fhm(ix,iy,1,0,is)=0.5_R8*fna(ix,iy,1,0,is)*
     2         (kinrgy(bottomix(ix,iy),bottomiy(ix,iy),is)+
     3          kinrgy(ix,iy,is)) - 0.5_R8*cvsa(ix,iy,1,0,is)*
     4         (ua(ix,iy,is)+ua(bottomix(ix,iy),bottomiy(ix,iy),is))*
     5         (ua(ix,iy,is)-ua(bottomix(ix,iy),bottomiy(ix,iy),is))
              fhm(ix,iy,1,1,is)=0.5_R8*fna(ix,iy,1,1,is)*
     2         (kinrgy(bottomix(ix,iy),bottomiy(ix,iy),is)+
     3          kinrgy(ix,iy,is)) - 0.5_R8*cvsa(ix,iy,1,1,is)*
     4         (ua(ix,iy,is)+ua(bottomix(ix,iy),bottomiy(ix,iy),is))*
     5         (ua(ix,iy,is)-ua(bottomix(ix,iy),bottomiy(ix,iy),is))
              fhp(ix,iy,1,0,is)=0.5_R8*fna(ix,iy,1,0,is)*
     2         (rpt(bottomix(ix,iy),bottomiy(ix,iy),is)+
     1          rpt(ix,iy,is))*ev
              fhp(ix,iy,1,1,is)=0.5_R8*fna(ix,iy,1,1,is)*
     2         (rpt(bottomix(ix,iy),bottomiy(ix,iy),is)+
     1          rpt(ix,iy,is))*ev
              fnt(ix,iy,1,0)=fnt(ix,iy,1,0)+0.5_R8*fna_53(ix,iy,1,0,is)*
     2         (ti(bottomix(ix,iy),bottomiy(ix,iy))+
     1          ti(ix,iy))
              fnt(ix,iy,1,1)=fnt(ix,iy,1,1)+0.5_R8*fna_53(ix,iy,1,1,is)*
     2         (ti(bottomix(ix,iy),bottomiy(ix,iy))+
     1          ti(ix,iy))
            endif
          enddo
        enddo
      enddo
c calculate total fluxes (fht)
      fht = fhe + fhi + fhj + fnt*(1.0_R8-BoRiS)
      do is=0,ns-1
        fht = fht + fhp(:,:,:,:,is) + fhm(:,:,:,:,is)*(1.0_R8-BoRiS)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpen
