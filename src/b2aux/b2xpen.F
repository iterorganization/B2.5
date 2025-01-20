      subroutine b2xpen (nCv, nFc, ns, mpg, switch,
     1 fna, fna_53, fne_53, fch, fhe, fhi, fhn, rpt, te, ti, tn, po,
     2 fhm, fhp, fhj, fht)
      use b2mod_types
      use b2mod_constants
      use b2us_map
      use b2mod_b2cmpa
      use b2mod_switches
      implicit none
      integer nCv, nFc, ns
! inputs
      type (mapping), intent (in) :: mpg
      type (switches), intent (inout) :: switch
      real (kind=R8) ::
     *  fna(nFc,0:1,0:ns-1), fna_53(nFc,0:1,0:ns-1), fne_53(nFc,0:1),
     *  fch(nFc,0:1), fhe(nFc,0:1), fhi(1:nFc,0:1), fhn(1:nFc,0:1),
     *  rpt(nCv,0:ns-1), te(nCv), ti(nCv), tn(nCv), po(nCv)
! outputs
      real (kind=R8) ::
     *  fhm(nFc,0:1,0:ns-1),fhp(nFc,0:1,0:ns-1),
     *  fhj(nFc,0:1),       fht(nFc,0:1)
*     ------------------------------------------------------------------
*     B2XPEN computes contributions to the total energy fluxes,
*     as well as the total energy fluxes
*     The contributions from external species are not included and must
*     be added explicitly in the calling routine when needed.
*     ------------------------------------------------------------------
      integer iFc, is
      real (kind=R8) ::
     &  fnt(nFc,0:1)
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpen')

c     calculate additional fluxes (fhj, fhp)
cwdk  note that fhm is now computed in b2urmo
      do iFc=1,nFc
        fhj(iFc,0:1)=-0.5_R8*fch(iFc,0:1)*
     1    (po(mpg%fcCv(iFc,1))+po(mpg%fcCv(iFc,2)))
        fnt(iFc,0:1)= 0.5_R8*fne_53(iFc,0:1)*
     1    (te(mpg%fcCv(iFc,1))+te(mpg%fcCv(iFc,2)))
      enddo
      do is=0,ns-1
        do iFc=1,nFc
          fhp(iFc,0:1,is)=0.5_R8*fna(iFc,0:1,is)*
     2      (rpt(mpg%fcCv(iFc,1),is)+rpt(mpg%fcCv(iFc,2),is))*ev
          if (.not.is_neutral(is)) then
            fnt(iFc,0:1)=fnt(iFc,0:1)+0.5_R8*fna_53(iFc,0:1,is)*
     2        (ti(mpg%fcCv(iFc,1))+ti(mpg%fcCv(iFc,2)))
          else
            fnt(iFc,0:1)=fnt(iFc,0:1)+0.5_R8*fna_53(iFc,0:1,is)*
     2        (tn(mpg%fcCv(iFc,1))+tn(mpg%fcCv(iFc,2)))
          endif
        enddo
      enddo
c calculate total fluxes (fht)
      fht = fhe + fhi + fhn + fhj + fnt*(1.0_R8-switch%BoRiS)
      do is=0,ns-1
        fht = fht + fhp(:,:,is) + fhm(:,:,is)*(1.0_R8-switch%BoRiS)
      enddo

      if (switch%iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,fht(1,0),'b2xpen_fht_th')
        call my_out_us(70,nFc,1,fht(1,1),'b2xpen_fht_r')
      endif
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpen
