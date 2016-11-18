*-----------------------------------------------------------------------
*.specification

      subroutine b2tfhe (nx, ny, ns, guard_flows, facdrift,
     1  BoRiS, vol, hx, hy, qz, te, po, ne, vedia, veecrb,
     2  fne_32, fne_52, fch_p, csig, chce, chve, calf, chvemx,           !xpb
     &  cddi, fhe, fhe_mdf, fhePSch, fhet, floe, cone, fllime)           !srv 13.10.06
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_anomalous_transport                             	 !xpb
      use b2mod_geo_corner                                      	 !xpb
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, guard_flows                                    !xpb
      real (kind=R8) ::
     * te(-1:nx,-1:ny), po(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     & vol(-1:nx,-1:ny), facdrift(-1:nx,-1:ny),
     & fne_32(-1:nx,-1:ny,0:1), fne_52(-1:nx,-1:ny,0:1),
     & fch_p(-1:nx,-1:ny,0:1),
     1 csig(-1:nx,-1:ny,0:1), 
     & chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),
     & calf(-1:nx,-1:ny,0:1), chvemx(-1:nx,-1:ny),
     & cddi(-1:nx,-1:ny,0:1,0:ns-1),                                     !srv 13.10.06
     1 hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1),
     * vedia(-1:nx,-1:ny,0:1), veecrb(-1:nx,-1:ny,0:1),
     2 BoRiS
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  fhe(-1:nx,-1:ny,0:1),
     &  fhet(-1:nx,-1:ny,0:1), 
     &  fhe_mdf(-1:nx,-1:ny,0:1), fhePSch(-1:nx,-1:ny,0:1),              !srv 13.10.06
     1  floe(-1:nx,-1:ny,0:1), cone(-1:nx,-1:ny,0:1,0:2),            !srv 02.01.07
     2  fllime(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFHE computes the electron heat flux between neighboring cells,
*     and the coefficients that specify the linearised response to
*     changes in electron temperature.
*
*
*  8. further comments
*
*     In the evaluation of the thermo-electric term in the electron
*     heat flux, the expression harm(te(),te()) cannot be replaced
*     by the arithmetic average (this is an empirical result).
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, ib, is, ixl, iyl, ixb, iyb,
     *   no_current, conduction_only, no_hybr, pot_eq
      integer, save :: iout = 0                                         !srv 18.06.08
      integer, save :: lim_flux = 0, mdf_fhe = 0                        !srv 04.02.03 13.10.06
      integer, save :: drift_style = 1
      integer, save :: conductive_limit = 1
      real (kind=R8), save :: alfTeEh = 0.0_R8                          !srv 28.04.08 30.05.08
      real (kind=R8), save :: fch_pTe = 1.0_R8                          !srv 30.12.09
      real (kind=R8) ::
     *   ehx, ehy, flo, dif, floe0, flomx, cone0, alfe0, t0, t1, 
     &   tex, tex0, cutlo,
     1   flo_vedia, flo_veecrb, prefactor,
     &   floe0_mdf, cone0_mdf, alfe0_mdf, t0_mdf, tex0_mdf              !srv 22.11.99 13.10.06
      real (kind=R8), save :: PSch=1.0e0_R8                             !srv 17.01.07
      parameter (cutlo = 1.0e-60_R8)
#ifdef CLASSICAL_CURRENT
      real (kind=R8) ::
     *   gampvd(-1:nx,-1:ny), gamrvd(-1:nx,-1:ny),
     *   gampvdc(-1:nx,-1:ny)
#endif
      save ncall, no_current, conduction_only, no_hybr, pot_eq
*   ..procedures
      intrinsic dim, min, max, abs, sqrt, log, mod
      real (kind=R8) ::
     *   harm, hybr, upwind, hy1
      external subini, subend, xertst, ipgeti
      external b2xvsg, b2xvff, b2xvfx, intcor
      harm(t0,t1) = 2*t0*t1/(t0+t1)
      hybr(flo,dif) = max(dif,abs(flo/2),sqrt(dif**2+flo**2/4)) !!!
      upwind(flo,t0,t1) = max(flo,0.0_R8)*t0+min(flo,0.0_R8)*t1
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
*   ..initialisation
      data ncall/0/, no_current/0/, conduction_only/0/, no_hybr/0/, 
     *  pot_eq/1/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfhe')
      if (ncall.eq.0) then                                               !xpb
       call ipgeti ('b2tfhe_no_current', no_current)                     !xpb
       call ipgeti ('b2news_poteq', pot_eq)                         !srv 27.02.13
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
       call ipgeti ('b2tfhe_conduction_only', conduction_only)           !xpb
       call ipgeti ('b2tfhe_lim_flux', lim_flux)
       call ipgeti ('b2trcl_conductive_limit', conductive_limit) 
       call ipgeti ('b2tfhe_no_hybr', no_hybr)                           !xpb
       call ipgeti ('b2tfnb_drift_style', drift_style)
       call ipgeti ('b2tfhe_mdf', mdf_fhe)                               !srv 13.10.06
       call ipgetr ('b2tfhe_PSch', PSch)                                 !srv 17.01.07
       call ipgeti ('b2tfhe_iout', iout)                                 !srv 18.06.08
       call ipgetr ('b2tfhe_alfTeEh', alfTeEh)                           !srv 28.04.08
       call ipgetr ('b2tfhe_fch_pTe', fch_pTe)                           !srv 30.12.09
       if(fch_pTe .eq. 1.0_R8) alfTeEh = 0.0_R8                          !srv 16.04.10
       call xertst (.not.(0.0_R8.ne.alfTeEh.and.0.0_R8.ne.fch_pTe),
     &        ' Do not set 0.0_R8.ne.alfTeEh.and.0.0_R8.ne.fch_pTe')
      endif                                                              !xpb
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of fne
       call b2xvff (nx, ny, fne_32, 'fne_32')
       call b2xvff (nx, ny, fne_52, 'fne_52')
*    ..test sign of ne, te, csig, chce
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2*2, csig, 1, 'csig', '.ge.')
       call b2xvsg (n2*2, chce, 1, 'chce', '.ge.')
*    ..test edge values of csig, chce, chve, chvemx, calf, calf_an
       call b2xvff (nx, ny, csig, 'csig')
       call b2xvff (nx, ny, chce, 'chce')
       call b2xvff (nx, ny, chve, 'chve')
       call b2xvfx (nx, ny, chvemx, 'chvemx')
       call b2xvff (nx, ny, calf, 'calf')
       call b2xvff (nx, ny, calf_an, 'calf_an')
      endif

      fhePSch = 0.0e0_R8                                                 !sv 13.10.06
      cone    = 0.0e0_R8                                                 !sv 02.01.07
      if(maxval(facdrift).ne.0.0e0_R8) then
*   ..compute P.Sch heat x-flux                                          !sv 13.10.06 {
       do is = 0,ns 
	if (nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then
!        iy = -1
!         do ix = -1, nx
!	  if(leftix(ix,iy).eq.-2) then
!           fhePSch(ix,iy,0)=0.0e0_R8
!	  else
!           fhePSch(ix,iy,0)=-0.5e0_R8*(
!       &    cddi(ix,topiy(ix,iy),1,is)*
!       &    (ne(ix,topiy(ix,iy))*te(ix,topiy(ix,iy))**2-
!       &                        ne(ix,iy)*te(ix,iy)**2)+
!       &    cddi(leftix(ix,iy),topiy(ix,iy),1,is)*
!       &     (ne(leftix(ix,iy),topiy(ix,iy))*
!       &      te(leftix(ix,iy),topiy(ix,iy))**2
!       &                    -ne(leftix(ix,iy),iy)*te(leftix(ix,iy),iy)**2)
!       &                               )
!	  endif
!	 enddo

         do iy = 0, ny-1
          do ix = -1, nx
	   if(leftix(ix,iy).eq.-2) then
            fhePSch(ix,iy,0)=0.0e0_R8
	   else
	    ixl = leftix(ix,iy)
	    iyl = leftiy(ix,iy)
            fhePSch(ix,iy,0)=-0.25e0_R8*PSch*(
     &       cddi(topix(ix,iy),topiy(ix,iy),1,is)*
     &        (ne(topix(ix,iy),topiy(ix,iy))*
     &         te(topix(ix,iy),topiy(ix,iy))**2-ne(ix,iy)*te(ix,iy)**2)+
     &       cddi(topix(ixl,iyl),topiy(ixl,iyl),1,is)*
     &        (ne(topix(ixl,iyl),topiy(ixl,iyl))*
     &         te(topix(ixl,iyl),topiy(ixl,iyl))**2
     &                     -ne(ixl,iyl)*te(ixl,iyl)**2)+
     &       cddi(ixl,iyl,1,is)*
     &       (ne(ixl,iyl)*te(ixl,iyl)**2
     &                     -ne(bottomix(ixl,iyl),bottomiy(ixl,iyl))*
     &                      te(bottomix(ixl,iyl),bottomiy(ixl,iyl))**2)+
     &       cddi(ix,iy,1,is)*
     &       (ne(ix,iy)*te(ix,iy)**2-ne(bottomix(ix,iy),bottomiy(ix,iy))
     &                          *te(bottomix(ix,iy),bottomiy(ix,iy))**2)
     &                                  )
            cone(ix,iy,0,1) = -0.5e0_R8*PSch*cddi(ix,iy,1,is)*  ! 0.25*2  !sv 02.01.07
     *                                               ne(ix,iy)*te(ix,iy)
            cone(ix,iy,0,2) = -0.5e0_R8*PSch*cddi(ix,iy,1,is)*            ! 0.25*2  !sv 02.01.07
     &                         ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &                         te(bottomix(ix,iy),bottomiy(ix,iy))
	   endif
	  enddo
         enddo

!         iy = ny
!          do ix = -1, nx
!	   if(leftix(ix,iy).eq.-2) then
!            fhePSch(ix,iy,0)=0.0e0_R8
!	   else
!            fhePSch(ix,iy,0)=-0.5e0_R8*(
!        &    cddi(leftix(ix,iy),iy,1,is)*
!        &    (ne(leftix(ix,iy),iy)*te(leftix(ix,iy),iy)**2
!        &                         -ne(leftix(ix,iy),bottomiy(ix,iy))*
!        &                            te(leftix(ix,iy),bottomiy(ix,iy))**2)+
!        &    cddi(ix,iy,1,is)*
!        &    (ne(ix,iy)*te(ix,iy)**2
!        &                -ne(ix,bottomiy(ix,iy))*te(ix,bottomiy(ix,iy))**2)
!        &                               )
!	   endif
!	  enddo

*   .....compute P.Sch y-flux
         iy = -1
         fhePSch(:,iy,1)=0.0e0_R8
         do iy = 0, ny
          do ix = -1, nx
!	   if(leftix(ix,iy).eq.-2) then
!            fhePSch(ix,iy,1) = 0.5e0_R8*(
!        &    cddi(rightix(ix,iy),iy,0,is)*
!        &    (ne(rightix(ix,iy),iy)*te(rightix(ix,iy),iy)**2
!        &                                        -ne(ix,iy)*te(ix,iy)**2)+
!        &    cddi(rightix(ix,iy),bottomiy(ix,iy),0,is)*
!        &    (ne(rightix(ix,iy),bottomiy(ix,iy))*
!        &             te(rightix(ix,iy),bottomiy(ix,iy))**2
!        &                -ne(ix,bottomiy(ix,iy))*te(ix,bottomiy(ix,iy))**2)
!        &                               )
!	   elseif(rightix(ix,iy).eq.nx+1) then
!            fhePSch(ix,iy,1) = 0.5e0_R8*(
!        &    cddi(ix,iy,0,is)*
!        &    (ne(ix,iy)*te(ix,iy)**2
!        &                   -ne(leftix(ix,iy),iy)*te(leftix(ix,iy),iy)**2)+
!        &    cddi(ix,bottomiy(ix,iy),0,is)*
!        &    (ne(ix,bottomiy(ix,iy))*te((ix),bottomiy(ix,iy))**2
!        &                         -ne(leftix(ix,iy),bottomiy(ix,iy))*
!        &                            te(leftix(ix,iy),bottomiy(ix,iy))**2)
!        &                               )
!	   else
	   if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            fhePSch(ix,iy,1) = 0.0e0_R8
	   else
	    ixb = bottomix(ix,iy)
	    iyb = bottomiy(ix,iy)
            fhePSch(ix,iy,1) = 0.25e0_R8*PSch*(
     &       cddi(rightix(ix,iy),rightiy(ix,iy),0,is)*
     &        (ne(rightix(ix,iy),rightiy(ix,iy))*
     &         te(rightix(ix,iy),rightiy(ix,iy))**2
     &                                         -ne(ix,iy)*te(ix,iy)**2)+
     &       cddi(ix,iy,0,is)*
     &       (ne(ix,iy)*te(ix,iy)**2
     &                      -ne(leftix(ix,iy),leftiy(ix,iy))*
     &                       te(leftix(ix,iy),leftiy(ix,iy))**2)+
     &       cddi(bottomix(ix,iy),bottomiy(ix,iy),0,is)*
     &        (ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &         te(bottomix(ix,iy),bottomiy(ix,iy))**2
     &                         -ne(leftix(ixb,iyb),leftiy(ixb,iyb))*
     &                          te(leftix(ixb,iyb),leftiy(ixb,iyb))**2)+
     &       cddi(rightix(ixb,iyb),rightiy(ixb,iyb),0,is)*
     &        (ne(rightix(ixb,iyb),rightiy(ixb,iyb))*
     &                          te(rightix(ixb,iyb),rightiy(ixb,iyb))**2
     &                         -ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &                          te(bottomix(ix,iy),bottomiy(ix,iy))**2)
     &                                  )
            cone(ix,iy,1,1) = -0.5e0_R8*PSch*cddi(ix,iy,0,is)*   ! 0.25*2 !sv 02.01.07 {
     &                                              ne(ix,iy)*te(ix,iy)
            cone(ix,iy,1,2) = -0.5e0_R8*PSch*cddi(ix,iy,0,is)*             ! 0.25*2
     &                           ne(leftix(ix,iy),leftiy(ix,iy))*
     &                           te(leftix(ix,iy),leftiy(ix,iy))      !sv 02.01.07 }
	    endif
	   enddo
          enddo
	 endif
        enddo
        do iy = -1, ny                                                  !srv 02.01.07 {
        do ix = -1, nx
         fhePSch(ix,iy,:)=2.5e0_R8*fhePSch(ix,iy,:)*facdrift(ix,iy)
         cone(ix,iy,:,1) = cone(ix,iy,:,1)*facdrift(ix,iy)
         cone(ix,iy,:,2) = cone(ix,iy,:,2)*facdrift(ix,iy)
	enddo
       enddo                                                           !srv 02.01.07 }

c  ...set null flux in guard cells  ?
!       fhePSch(-1, :,:)=0.0e0_R8    ! MODIFY for Double Case 
!       fhePSch(nx, :,1)=0.0e0_R8
!       fhePSch( :,-1,:)=0.0e0_R8
!       fhePSch( :,ny,0)=0.0e0_R8
      endif

      fhe_mdf = 0.0e0_R8                                              !sv 13.10.06

* ..compute fluxes and linearisation
*  ..do not forget equivalent code in B2TFCH !!!!!!!!!!!!!!!!!!!!!
*   ..compute fhex, floex, conex
      do iy = -1, ny
        do ix = -1, nx
          if(leftix(ix,iy).eq.-2) then
            floe(ix,iy,0) = 0.0_R8
            cone(ix,iy,0,0) = 0.0_R8                                 !srv 02.01.07
            fhe_mdf(ix,iy,0) = 0.0e0_R8                              !srv 13.10.06 
            fhe(ix,iy,0) = 0.0_R8
            fhet(ix,iy,0) = 0.0_R8
            fllime(ix,iy) = 1.0_R8
          else
!xpb        ehx = -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))+
!xpb     &         (ne(ix,iy)*te(ix,iy)
!xpb     &         -ne(leftix(ix,iy),leftiy(ix,iy))*
!xpb     &         te(leftix(ix,iy),leftiy(ix,iy)))*
!xpb     &         (2/qe)/(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))
            if (no_current.eq.1) then
              ehx = 0.0_R8
            else
              ehx = -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))+ 	 !xpb
     &         (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))+      	 !xpb
     &         (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))*0.5_R8* 	 !xpb
     &         (log(ne(ix,iy))-log(ne(leftix(ix,iy),leftiy(ix,iy)))) 	 !xpb 
     &         )/qe                                             	 !xpb 
            endif
            tex = te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))
cxpb Changed the solution to (3/2+BoRiS) for benchmarking !!!!!!!!!!!!!!!!!!!
            if (conduction_only.eq.1) then
              floe0 = 0.0_R8
	      floe0_mdf = 0.0_R8                                         !srv 18.08.08
              flo_vedia=0.0_R8
              flo_veecrb=0.0_R8
            else
              floe0 = (1.5_R8+BoRiS)*fne_32(ix,iy,0)+
     1         2.5_R8*fne_52(ix,iy,0)+
     &         chve(ix,iy,0)*
     1           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))/2.0_R8
     &         -0.71e0_R8*fch_p(ix,iy,0)/qe*fch_pTe                      !srv 30.12.09
              prefactor=
     1         (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))    	 !xpb
     &         /(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))     	 !xpb
     &         *(ne(leftix(ix,iy),leftiy(ix,iy))                	 !xpb
     &          +ne(ix,iy))/2.0_R8
              if (drift_style.eq.0) then
                flo_vedia=0.5_R8*prefactor*(vedia(ix,iy,0) + 
     1            vedia(leftix(ix,iy),leftiy(ix,iy),0))
                flo_veecrb=0.5_R8*prefactor*(veecrb(ix,iy,0) + 
     1            veecrb(leftix(ix,iy),leftiy(ix,iy),0))
              else
                flo_vedia=prefactor*vedia(ix,iy,0)
                flo_veecrb=prefactor*veecrb(ix,iy,0)
              endif
            endif
            floe0_mdf=floe0+1.5_R8*flo_veecrb                            !srv 18.08.08
            floe0=floe0+2.5_R8*flo_vedia+1.5_R8*flo_veecrb
            if (no_hybr.eq.0) then
              cone0 = hybr(floe0,chce(ix,iy,0))-abs(floe0/2.0_R8)
              cone0_mdf = hybr(floe0_mdf,chce(ix,iy,0))-
     &                    abs(floe0_mdf/2.0_R8)                          !srv 18.08.08
            else
              cone0 = chce(ix,iy,0)
            endif
            if (chce(ix,iy,0).eq.0.0_R8) then
               alfe0 = 0.0_R8
               alfe0_mdf = 0.0_R8                                        !srv 18.08.08
            else
              alfe0 = calf(ix,iy,0)*cone0/chce(ix,iy,0)
              alfe0_mdf = calf(ix,iy,0)*cone0_mdf/chce(ix,iy,0)
            endif
cwdk 10.03.16 unified treatment of chvemx in b2tlhe.F. Test for b2tqce_model no longer required
            flomx = chvemx(ix,iy)*(ne(leftix(ix,iy),leftiy(ix,iy))+
     &       ne(ix,iy))/2.0_R8
            if (flomx.gt.0.0_R8 .and. 
     &        lim_flux.eq.1 .and. conductive_limit.eq.0) then
              t0 = cone0*abs(tex)/(flomx*max(
     &         te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
              t0_mdf = cone0_mdf*abs(tex)/(flomx*max(
     &         te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))          !srv 18.08.08
            else
              t0 = 0.0_R8
              t0_mdf = 0.0_R8                                       !srv 18.08.08
            endif

            tex0 = tex/(1.0_R8+t0)
            tex0_mdf = tex/(1.0_R8+t0_mdf)                           !srv 18.08.08
	    if(mdf_fhe.eq.0) then
              fllime(ix,iy) = 1.0_R8 + t0
            else
              fllime(ix,iy) = 1.0_R8 + t0_mdf
	    endif
            if (no_hybr.eq.0) then
              fhe(ix,iy,0) =
     &         floe0*(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
     &         -abs(floe0/2.0_R8)*tex-
     &          cone0*tex0 + alfTeEh*alfe0*
     &          harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              fhe_mdf(ix,iy,0) =floe0_mdf*                              !srv 18.08.08
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8-
     &         abs(floe0_mdf/2.0_R8)*tex-cone0_mdf*tex0_mdf+
     &         alfTeEh*alfe0_mdf*
     &         harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
     &         -fhePSch(ix,iy,0)                                       !srv 13.06.06 }
            else
              t1 = dim ( cone0/fllime(ix,iy),
     .         abs(floe0) *
     .         max(0.0_R8, 
     .             min(te(ix,iy),te(leftix(ix,iy),leftiy(ix,iy))))/
     .         max(cutlo, te(ix,iy)+te(leftix(ix,iy),leftiy(ix,iy))) )
              fhe(ix,iy,0) = 
     .         upwind(floe0, te(leftix(ix,iy),leftiy(ix,iy)), te(ix,iy))
     .         - t1 * (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))) 
     &         +alfTeEh*alfe0*
     &          harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            endif
            fhet(ix,iy,0) = fhe(ix,iy,0) +
     1       ((1.0_R8-BoRiS)*fne_32(ix,iy,0)+flo_veecrb)*
     &       (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
! modify fhet for mdf case
            if (no_hybr.eq.0) then
	      if(mdf_fhe.eq.0) then                                  !srv 13.10.06 {
                floe(ix,iy,0) = floe0
                cone(ix,iy,0,0) = abs(floe0/2)+cone0/(1.0_R8+t0)     !srv 02.01.07
	      else
                floe(ix,iy,0) = floe0_mdf
                cone(ix,iy,0,0)=abs(floe0_mdf/2)+
     .           cone0_mdf/(1.0_R8+t0_mdf)                           !srv 02.01.07
	      endif
            else
              cone(ix,iy,0,0) = t1                                   !srv 02.01.07
            endif
          endif
        enddo
      enddo
*   ..compute fhey, floey, coney
      do iy = -1, ny
        do ix = -1, nx
          if (bottomiy(ix,iy).eq.-2) then
            floe(ix,iy,1) = 0.0_R8
            cone(ix,iy,1,0) = 0.0_R8                                 !srv 02.01.07
            fhe_mdf(ix,iy,1) = 0.0e0_R8                              !srv 13.10.06
            fhe(ix,iy,1) = 0.0_R8
            fhet(ix,iy,1) = 0.0_R8
          else
            if (no_current.eq.1) then
              ehy = 0.0_R8
            else
              ehy = -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))+
     &         (ne(ix,iy)*te(ix,iy)-ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &         te(bottomix(ix,iy),bottomiy(ix,iy)))*
     &         (2.0_R8/qe)/
     &         (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))
            endif
            if (conduction_only.eq.1) then
              floe0 = 0.0_R8
              flo_vedia=0.0_R8
              flo_veecrb=0.0_R8
            else
              floe0 = (1.5_R8+BoRiS)*fne_32(ix,iy,1)+ 
     1         2.5_R8*fne_52(ix,iy,1)+ 
     &         chve(ix,iy,1)*
     1         (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))/2.0_R8
              prefactor=
     &         (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))	 !xpb
     &         /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy)) 	 !xpb
     &         *(ne(bottomix(ix,iy),bottomiy(ix,iy))            	 !xpb
     &          +ne(ix,iy))/2.0_R8                              	 !xpb
              if (drift_style.eq.0) then
                flo_vedia=0.5_R8*prefactor*(vedia(ix,iy,1) + 
     1            vedia(bottomix(ix,iy),bottomiy(ix,iy),1))
                flo_veecrb=0.5_R8*prefactor*(veecrb(ix,iy,1) + 
     1            veecrb(bottomix(ix,iy),bottomiy(ix,iy),1))
              else
                flo_vedia=prefactor*vedia(ix,iy,1)
                flo_veecrb=prefactor*veecrb(ix,iy,1)
              endif
            endif
	    if(mdf_fhe.eq.0) then                                        !srv 13.10.06 {
              floe(ix,iy,1) = floe0+2.5_R8*flo_vedia+1.5_R8*flo_veecrb
	    else
              floe(ix,iy,1) = floe0+1.5_R8*flo_veecrb
	    endif                                                        !srv 18.08.08 }
            if (no_hybr.eq.0) then
              cone0 = hybr(floe(ix,iy,1),chce(ix,iy,1))
              fhe(ix,iy,1) = floe(ix,iy,1)*
     &          (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8-
     &           cone0*(te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
     &          +alfTeEh*calf(ix,iy,1)*
     &           harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
              fhe_mdf(ix,iy,1) = floe(ix,iy,1)*
     &          (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8-
     &           cone0*(te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
     &          +alfTeEh*calf(ix,iy,1)*
     &           harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
     &            -fhePSch(ix,iy,1)                                     !srv 13.10.06
            else
              cone0 = chce(ix,iy,1)
              t1 = dim ( cone0,
     .         abs(floe(ix,iy,1)) *
     .         max(0.0_R8, 
     .             min(te(ix,iy),te(bottomix(ix,iy),bottomiy(ix,iy))))/
     .         max(cutlo, 
     .                 te(ix,iy)+te(bottomix(ix,iy),bottomiy(ix,iy))))
              fhe(ix,iy,1) = 
     .         upwind(floe(ix,iy,1), 
     .                te(bottomix(ix,iy),bottomiy(ix,iy)), te(ix,iy))
     .          - t1 * (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
     &          +alfTeEh*calf(ix,iy,1)*
     &           harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
            endif
            fhet(ix,iy,1) = fhe(ix,iy,1)+
     1    ((1.0_R8-BoRiS)*fne_32(ix,iy,1)+flo_veecrb)*
     &    (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
            cone(ix,iy,1,0) = cone0                                      !xpb !srv 02.01.07
!xpb     &         2.5_R8*((te(bottomix(ix,iy),bottomiy(ix,iy))+
!xpb     &     te(ix,iy))/(2*qe))*csig(ix,iy,1)/qe
          endif                                                 	 !xpb
        enddo
      enddo

      if(guard_flows.eq.0) then                                 	 !xpb
cxpb turn off flows between boundary cells
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
              fhe(ix,iy,1)=0.0_R8
              fhe_mdf(ix,iy,1)=0.0_R8
              fhet(ix,iy,1)=0.0_R8
              floe(ix,iy,1)=0.0_R8
              cone(ix,iy,1,0)=0.0_R8
            endif
          enddo        
        enddo

cxpb turn off poloidal flows in core boundary
        do ix = -1, nx
          if(mod(region(ix,-1,0),4).eq.1) then
            fhe(ix,-1,0) = 0.0_R8
            fhe_mdf(ix,-1,0) = 0.0_R8
            fhet(ix,-1,0)=0.0_R8
            floe(ix,-1,0) = 0.0_R8
            cone(ix,-1,0,0) = 0.0_R8
            fllime(ix,iy) = 1.0_R8
          endif
        enddo
      endif                                                     	 !xpb

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fhe(ix,iy,0)=0.0_R8
            fhe_mdf(ix,iy,0)=0.0_R8
            fhet(ix,iy,0)=0.0_R8
            floe(ix,iy,0)=0.0_R8
            cone(ix,iy,0,0)=0.0_R8
            fllime(ix,iy)=1.0_R8
            fhe(ix,iy,1)=0.0_R8
            fhe_mdf(ix,iy,1)=0.0_R8
            fhet(ix,iy,1)=0.0_R8
            floe(ix,iy,1)=0.0_R8
            cone(ix,iy,1,0)=0.0_R8
          end do
        end do
      end do

      if(iout.ne.0) then                                             !srv 18.06.08 {
        call my_out(70,nx,ny,fhe(-1,-1,0),'b2tfhe_fhex')
        call my_out(70,nx,ny,fhe(-1,-1,1),'b2tfhe_fhey')
        call my_out(70,nx,ny,fhePSch(-1,-1,0),'b2tfhe_fhePSchx')
        call my_out(70,nx,ny,fhePSch(-1,-1,1),'b2tfhe_fhePSchy')
        call my_out(70,nx,ny,fhe_mdf(-1,-1,0),'b2tfhe_fhe_mdfx')
        call my_out(70,nx,ny,fhe_mdf(-1,-1,1),'b2tfhe_fhe_mdfy')
        call my_out(70,nx,ny,fhet(-1,-1,0),'b2tfhe_fhetx')
        call my_out(70,nx,ny,fhet(-1,-1,1),'b2tfhe_fhety')
	call my_out(70,nx,ny,cone(-1,-1,0,0),'b2tfhe_conex0')
	call my_out(70,nx,ny,cone(-1,-1,0,1),'b2tfhe_conex1')
	call my_out(70,nx,ny,cone(-1,-1,0,2),'b2tfhe_conex2')
	call my_out(70,nx,ny,cone(-1,-1,1,0),'b2tfhe_coney0')
	call my_out(70,nx,ny,cone(-1,-1,1,1),'b2tfhe_coney1')
	call my_out(70,nx,ny,cone(-1,-1,1,2),'b2tfhe_coney2')
	call my_out(70,nx,ny,floe(-1,-1,0),'b2tfhe_floex')
	call my_out(70,nx,ny,floe(-1,-1,1),'b2tfhe_floey')
	call my_out(70,nx,ny,fllime,'b2tfhe_fllime')
      endif                                                          !srv 18.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

cxpb  This version of the routine uses a (1.5+BoRiS) factor for the
cxpb  poloidal convective heat flux instead of the usual 3/2 or 5/2 
cxpb  so as to solve the same equation as in BoRiS for benchmark purposes. 
cxpb  The default is BoRiS = 0.0.
cxpb  For comparison, set BoRiS = 1.0.
cxpb  Routine rewritten to explicitly separate 3/2 and 5/2 terms.
cxpb  The contributions from fne_32 and fne_52 are NO DRIFTS
cxpb  contributions, i.e. diamagnetic and ExB terms are treated
cxpb  separately.

*-----------------------------------------------------------------------
*.end b2tfhe

      end subroutine b2tfhe
