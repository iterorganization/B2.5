      subroutine b2xpen (nCv, nFc, ns, mpg,
     1 fna, fne, fch, fhe, fhi, kinrgy, rpt, te, ti, po, BoRiS,
     2 fhm, fhp, fhj, fht)
      use b2mod_types
      use b2mod_constants
      use b2us_map
      implicit none
      integer nCv, nFc, ns
! inputs
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  fna(nFc,0:1,0:ns-1), fne(nFc,0:1), fch(nFc,0:1),
     *  fhe(nFc,0:1), fhi(1:nFc,0:1), kinrgy(nCv,0:ns-1), 
     *  rpt(nCv,0:ns-1), te(nCv), ti(nCv), po(nCv), BoRiS
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
c calculate additional fluxes (fhj, fhm, fhp)
      do iFc=1,nFc
        fhj(iFc,0:1)=-0.5_R8*fch(iFc,0:1)*
     1    (po(mpg%fcCv(iFc,1))+po(mpg%fcCv(iFc,2)))
        fnt(iFc,0:1)= 0.5_R8*fne(iFc,0:1)*
     1    (te(mpg%fcCv(iFc,1))+te(mpg%fcCv(iFc,2)))
      enddo
      do is=0,ns-1
        do iFc=1,nFc
          fhm(iFc,0:1,is)=0.5_R8*fna(iFc,0:1,is)*
     2      (kinrgy(mpg%fcCv(iFc,1),is)+kinrgy(mpg%fcCv(iFc,2),is))
          fhp(iFc,0:1,is)=0.5_R8*fna(iFc,0:1,is)*
     2      (rpt(mpg%fcCv(iFc,1),is)+rpt(mpg%fcCv(iFc,2),is))*ev
          fnt(iFc,0:1)=fnt(iFc,0:1)+0.5_R8*fna(iFc,0:1,is)*
     2      (ti(mpg%fcCv(iFc,1))+ti(mpg%fcCv(iFc,2)))
        enddo
      enddo
c calculate total fluxes (fht)
      fht = fhe + fhi + fhj + fnt*(1.0_R8-BoRiS)
      do is=0,ns-1
	fht = fht + fhp(:,:,is) + fhm(:,:,is)*(1.0_R8-BoRiS)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpen
