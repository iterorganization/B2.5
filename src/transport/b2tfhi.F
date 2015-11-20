*-----------------------------------------------------------------------
*.specification

      subroutine b2tfhi (nx, ny, ns, guard_flows, facdrift, BoRiS,
     1  vol, hx, hy, gs, qz, qc, bb, ti, ni, na,
     2  fna, fni_32, fni_52, kinrgy,
     3  chci, cddi, chvi, chvimx, vsal, vadia, vaecrb,
     &  fhi, fhi_mdf, fhiPSch, fhit, floi, coni, fllimi)             !srv 13.10.06
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, guard_flows
      real (kind=R8) ::
     *  ti(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  fna(-1:nx,-1:ny,0:1,0:ns-1), cddi(-1:nx,-1:ny,0:1,0:ns-1),    !srv 13.10.06
     1  fni_32(-1:nx,-1:ny,0:1), fni_52(-1:nx,-1:ny,0:1),
     &  vol(-1:nx,-1:ny), chci(-1:nx,-1:ny,0:1), chvi(-1:nx,-1:ny,0:1),
     &  chvimx(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny), 
     1  qz(-1:nx,-1:ny,0:1), bb(-1:nx,-1:ny,0:3), 
     .  na(-1:nx,-1:ny,0:ns-1), vsal(-1:nx,-1:ny,0:1,0:ns-1),
     2  vadia(-1:nx,-1:ny,0:1,0:ns-1), vaecrb(-1:nx,-1:ny,0:1,0:ns-1),
     3  facdrift(-1:nx,-1:ny), kinrgy(-1:nx,-1:ny,0:ns-1), BoRiS
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  fhi(-1:nx,-1:ny,0:1), fhit(-1:nx,-1:ny,0:1), 
     &  fhi_mdf(-1:nx,-1:ny,0:1), fhiPSch(-1:nx,-1:ny,0:1),          !srv 13.10.06
     1  floi(-1:nx,-1:ny,0:1), coni(-1:nx,-1:ny,0:1,0:2),            !srv 02.01.07
     2  fllimi(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFHI computes the all atom heat flux between neighboring cells,
*     and the coefficients that specify its linearisation.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the atom heat flux, fhi, and the
*     coefficients floi and coni that linearise the response of the
*     atom heat flux to changes in the atom temperature.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, ib, conduction_only, no_hybr, 
     &        ixl, iyl, ixb, iyb
      integer, save :: iout = 0                                        !sv 16.09.02
      integer, save :: lim_flux = 0, mdf_fhi = 0                       !srv 04.02.03 13.10.06
      real (kind=R8), save :: PSch=1.0e0_R8                            !srv 17.01.07
      real (kind=R8) ::
     *  flo, dif, floi0, flomx, coni0, t0, t1, dti, dti0, w1,
     1  cum_flo_vadia, cum_flo_vaecrb, prefactor, cutlo,
     &   floi0_mdf, coni0_mdf, t0_mdf,dti0_mdf, t1_mdf                  !srv 13.10.06
      parameter (cutlo=1.0e-60_R8)
      save ncall, conduction_only, no_hybr
*   ..procedures
      intrinsic min, max, dim, abs, sqrt, sign
      real (kind=R8) ::
     *  hybr, upwind, hy1, ave
      external subini, subend, xertst, ipgeti
      external b2xvsg, b2xvff
      hybr(flo,dif) = max(dif,abs(flo/2),sqrt(dif**2+flo**2/4)) !!!
      upwind(flo,t0,t1) = max(flo,0.0_R8)*t0+min(flo,0.0_R8)*t1
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
      ave(t0,t1) = 2.0_R8*t0*t1 / (cutlo+t0+t1)
*   ..initialisation
      data ncall/0/, conduction_only/0/, no_hybr/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfhi')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test facdrift
      call xertst (0.0_R8.le.minval(facdrift) .and.
     &       1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
*   ..get input parameters
      if (ncall.eq.0) then
        call ipgeti ('b2tfhe_conduction_only', conduction_only)
        call ipgeti ('b2tfhi_no_hybr', no_hybr)
        call ipgeti ('b2tfhi_mdf', mdf_fhi)                               !srv 13.10.06
        call ipgetr ('b2tfhi_PSch', PSch)                                 !srv 17.01.07
        call ipgeti ('b2tfhi_iout', iout)                                 !sv 16.09.02
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of fni
       call b2xvff (nx, ny, fni_32, 'fni_32')
       call b2xvff (nx, ny, fni_52, 'fni_52')
*    ..test sign of na, ni, ti, chci
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2*2, chci, 1, 'chci', '.ge.')
*    ..test edge values of chci, chvi
       call b2xvff (nx, ny, chci, 'chci')
       call b2xvff (nx, ny, chvi, 'chvi')
      endif

      fhiPSch = 0.0e0_R8                                              !sv 13.10.06
      coni    = 0.0e0_R8                                              !sv 02.01.07
      if(maxval(facdrift).ne.0.0e0_R8) then
*   ..compute P.Sch heat x-flux                                 !sv 13.10.06 {
!      iy = -1
!       do ix = -1, nx
!	  if(leftix(ix,iy).eq.-2) then
!         fhiPSch(ix,iy,0)=0.0e0_R8
!	  else
!	   do is = 0, ns-1
!	    if(za(is) .ne. 0.0e0_R8) then
!          fhiPSch(ix,iy,0) = fhiPSch(ix,iy,0) + 0.5e0_R8*(
!     &    cddi(ix,topiy(ix,iy),1,is)*
!     &    (na(ix,topiy(ix,iy),is)*ti(ix,topiy(ix,iy))**2-
!     &                        na(ix,iy,is)*ti(ix,iy)**2)+
!     &    cddi(leftix(ix,iy),topiy(ix,iy),1,is)*
!     &     (na(leftix(ix,iy),topiy(ix,iy),is)*
!     &      ti(leftix(ix,iy),topiy(ix,iy))**2
!     &                 -na(leftix(ix,iy),iy,is)*ti(leftix(ix,iy),iy)**2)
!     &                               )
!	    endif
!	   enddo
!	  endif
!	 enddo

      do iy = 0, ny-1
       do ix = -1, nx
	  if(leftix(ix,iy).eq.-2) then
         fhiPSch(ix,iy,0) = 0.0e0_R8
	   coni(ix,iy,0,:) = 0.0e0_R8
	  else
	   do is = 0, ns-1
	    if(.not.is_neutral(is)) then
	     ixl = leftix(ix,iy)
	     iyl = leftiy(ix,iy)
           fhiPSch(ix,iy,0) = fhiPSch(ix,iy,0) + 0.25e0_R8*PSch*(
     &     cddi(topix(ix,iy),topiy(ix,iy),1,is)*
     &     (na(topix(ix,iy),topiy(ix,iy),is)*
     &      ti(topix(ix,iy),topiy(ix,iy))**2-na(ix,iy,is)*ti(ix,iy)**2)+
     &     cddi(topix(ixl,iyl),topiy(ixl,iyl),1,is)*
     &     (na(topix(ixl,iyl),topiy(ixl,iyl),is)*
     &      ti(topix(ixl,iyl),topiy(ixl,iyl))**2
     &                -na(ixl,iyl,is)*ti(ixl,iyl)**2)+
     &     cddi(ixl,iyl,1,is)*
     &     (na(ixl,iyl,is)*ti(ixl,iyl)**2
     &                     -na(bottomix(ixl,iyl),bottomiy(ixl,iyl),is)*
     &                      ti(bottomix(ixl,iyl),bottomiy(ixl,iyl))**2)+
     &     cddi(ix,iy,1,is)*
     &     (na(ix,iy,is)*ti(ix,iy)**2
     &             -na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &              ti(bottomix(ix,iy),bottomiy(ix,iy))**2)
     &                               )
           coni(ix,iy,0,1) = coni(ix,iy,0,1) + 0.5e0_R8*PSch*       ! 0.25*2  !sv 02.01.07
     *               cddi(ix,iy,1,is)*na(ix,iy,is)*ti(ix,iy)
           coni(ix,iy,0,2) = coni(ix,iy,0,2) + 0.5e0_R8*PSch*       ! 0.25*2  !sv 02.01.07
     &     cddi(ix,iy,1,is)*na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &                      ti(bottomix(ix,iy),bottomiy(ix,iy))
	    endif
	   enddo
	  endif
	 enddo
      enddo

!      iy = ny
!       do ix = -1, nx
!	  if(leftix(ix,iy).eq.-2) then
!         fhiPSch(ix,iy,0)=0.0e0_R8
!	  else
!	   do is = 0, ns-1
!          fhiPSch(ix,iy,0) = fhiPSch(ix,iy,0) + 0.5e0_R8*(
!     &     cddi(leftix(ix,iy),iy,1,is)*
!     &    (na(leftix(ix,iy),iy,is)*ti(leftix(ix,iy),iy)**2
!     &                         -na(leftix(ix,iy),bottomiy(ix,iy),is)*
!     &                            ti(leftix(ix,iy),bottomiy(ix,iy))**2)+
!     &     cddi(ix,iy,1,is)*
!     &    (na(ix,iy,is)*ti(ix,iy)**2
!     &             -na(ix,bottomiy(ix,iy),is)*ti(ix,bottomiy(ix,iy))**2)
!     &                               )
!	   enddo
!	  endif
!	 enddo

*   ..compute P.Sch y-flux
      iy = -1
      fhiPSch(:,iy,1)=0.0e0_R8
      do iy = 0, ny
       do ix = -1, nx
	  do is = 0, ns-1
!	   if(leftix(ix,iy).eq.-2) then
!          fhiPSch(ix,iy,1) = fhiPSch(ix,iy,1) - 0.5e0_R8*(
!     &     cddi(rightix(ix,iy),iy,0,is)*
!     &     (na(rightix(ix,iy),iy,is)*ti(rightix(ix,iy),iy)**2
!     &                                     -na(ix,iy,is)*ti(ix,iy)**2)+
!     &     cddi(rightix(ix,iy),bottomiy(ix,iy),0,is)*
!     &     (na(rightix(ix,iy),bottomiy(ix,iy),is)*
!     &             ti(rightix(ix,iy),bottomiy(ix,iy))**2
!     &           -na(ix,bottomiy(ix,iy),is)*ti(ix,bottomiy(ix,iy))**2)
!     &                               )
!	   elseif(rightix(ix,iy).eq.nx+1) then
!          fhiPSch(ix,iy,1) = fhiPSch(ix,iy,1) - 0.5e0_R8*(
!     &     cddi(ix,iy,0,is)*(na(ix,iy,is)*ti(ix,iy)**2
!     &                -na(leftix(ix,iy),iy,is)*ti(leftix(ix,iy),iy)**2)+
!     &     cddi(ix,bottomiy(ix,iy),0,is)*
!     &     (na(ix,bottomiy(ix,iy),is)*ti((ix),bottomiy(ix,iy))**2
!     &                         -na(leftix(ix,iy),bottomiy(ix,iy),is)*
!     &                            ti(leftix(ix,iy),bottomiy(ix,iy))**2)
!     &                               )
!	   else
	  if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
	    fhiPSch(ix,iy,1) = 0.0e0_R8
	    coni(ix,iy,1,:)  = 0.0e0_R8
	  else
	    ixb = bottomix(ix,iy)
	    iyb = bottomiy(ix,iy)
          fhiPSch(ix,iy,1) = fhiPSch(ix,iy,1) - 0.25e0_R8*PSch*(
     &     cddi(rightix(ix,iy),rightiy(ix,iy),0,is)*
     &      (na(rightix(ix,iy),rightiy(ix,iy),is)*
     &       ti(rightix(ix,iy),rightiy(ix,iy))**2
     &                                     -na(ix,iy,is)*ti(ix,iy)**2)+
     &     cddi(ix,iy,0,is)*(na(ix,iy,is)*ti(ix,iy)**2
     &                -na(leftix(ix,iy),leftiy(ix,iy),is)*
     &                 ti(leftix(ix,iy),leftiy(ix,iy))**2)+
     &     cddi(bottomix(ix,iy),bottomiy(ix,iy),0,is)*
     &      (na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &       ti(bottomix(ix,iy),bottomiy(ix,iy))**2
     &                         -na(leftix(ixb,iyb),leftiy(ixb,iyb),is)*
     &                          ti(leftix(ixb,iyb),leftiy(ixb,iyb))**2)+
     &     cddi(rightix(ixb,iyb),rightiy(ixb,iyb),0,is)*
     &      (na(rightix(ixb,iyb),rightiy(ixb,iyb),is)*
     &       ti(rightix(ixb,iyb),rightiy(ixb,iyb))**2
     &             -na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &              ti(bottomix(ix,iy),bottomiy(ix,iy))**2)
     &                               )
          coni(ix,iy,1,1) = coni(ix,iy,1,1)+0.5e0_R8*PSch*          ! 0.25*2 !sv 02.01.07 {
     &                      cddi(ix,iy,0,is)*na(ix,iy,is)*ti(ix,iy)
          coni(ix,iy,1,2) = coni(ix,iy,1,2)+0.5e0_R8*PSch*          ! 0.25*2
     &       cddi(ix,iy,0,is)*na(leftix(ix,iy),leftiy(ix,iy),is)*
     &                        ti(leftix(ix,iy),leftiy(ix,iy))       !sv 02.01.07 }
	   endif
	  enddo
	 enddo
      enddo

	do iy = -1, ny
       do ix = 0, nx
        fhiPSch(ix,iy,:) = 2.5e0_R8*fhiPSch(ix,iy,:)*facdrift(ix,iy)
        coni(ix,iy,:,1:2) = coni(ix,iy,:,1:2)*facdrift(ix,iy)          !sv 02.01.07
	 enddo
	enddo

c  ..set null flux in guard cells  ?
!      fhiPSch(-1, :,:)=0.0e0_R8    ! MODIFY for Double Case 
!      fhiPSch(nx, :,1)=0.0e0_R8
!      fhiPSch( :,-1,:)=0.0e0_R8
!      fhiPSch( :,ny,0)=0.0e0_R8
	endif

      fhi_mdf = 0.0e0_R8                                              !srv 13.10.06

* ..compute fluxes and linearisation
*   ..compute fhix, floix, conix
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         floi(ix,iy,0) = 0.0_R8
         coni(ix,iy,0,0) = 0.0_R8                                    !srv 02.01.07
         fhi_mdf(ix,iy,0) = 0.0_R8
         fhi(ix,iy,0) = 0.0_R8
         fhit(ix,iy,0) = 0.0_R8
         fllimi(ix,iy) = 1.0_R8
        else
         dti = ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy))
         w1 = 0.0_R8
cxpb Changed the solution to (3/2+BoRiS) for benchmarking !!!!!!!!!!!!!!!!!!!
         if (conduction_only.eq.1) then
           floi0 = 0.0_R8
           cum_flo_vadia=0.0_R8
           cum_flo_vaecrb=0.0_R8
         else
           if (BoRiS.eq.1.0_R8) then   ! Add the kinetic and viscous terms
             do is = 0, ns-1
               floi0 = kinrgy(ix,iy,is)*fna(ix,iy,0,is)
               coni0 = gs(ix,iy,0)*qc(ix,iy)*
     .          ave(vsal(ix,iy,0,is)/hx(ix,iy),
     .              vsal(leftix(ix,iy),leftiy(ix,iy),0,is)/
     .              hx(leftix(ix,iy),leftiy(ix,iy)))
               t1 = dim ( coni0,
     .          abs(floi0) *
     .          max(0.0_R8, 
     .             min(ti(ix,iy),ti(leftix(ix,iy),leftiy(ix,iy))))/
     .          max(cutlo, ti(ix,iy)+ti(leftix(ix,iy),leftiy(ix,iy))) )
               w1 = w1 +
     .          upwind(floi0,ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy))
     .           - t1 * (ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))
             enddo
           endif
           floi0 = (1.5_R8+BoRiS)*fni_32(ix,iy,0)+
     1      2.5_R8*fni_52(ix,iy,0)+
     &      chvi(ix,iy,0)*(ni(leftix(ix,iy),leftiy(ix,iy),0)+
     1                     ni(ix,iy,0))/2.0_R8
           cum_flo_vadia=0.0_R8
           cum_flo_vaecrb=0.0_R8
           do is=0,ns-1
             prefactor=
     1        (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))     !sv 17.05.99
     &        /(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))      !sv 17.05.99
     &        *(na(leftix(ix,iy),leftiy(ix,iy),is)              !sv 17.05.99
     &        *bb(leftix(ix,iy),leftiy(ix,iy),2)                !sv 17.05.99
     &        /bb(leftix(ix,iy),leftiy(ix,iy),3)                !sv 17.05.99
     &        +na(ix,iy,is)                                     !sv 17.05.99
     &        *bb(ix,iy,2)/bb(ix,iy,3))/2.0_R8
             cum_flo_vadia=cum_flo_vadia+
     1        prefactor*vadia(ix,iy,0,is)
             cum_flo_vaecrb=cum_flo_vaecrb+
     1        prefactor*vaecrb(ix,iy,0,is)
           enddo
         endif
         floi0_mdf = floi0+1.5_R8*cum_flo_vaecrb                !srv 18.08.08
         floi0 = floi0+2.5_R8*cum_flo_vadia+1.5_R8*cum_flo_vaecrb
         flomx = chvimx(ix,iy)*(ni(leftix(ix,iy),leftiy(ix,iy),0)
     &                         +ni(ix,iy,0))/2.0_R8
         if (no_hybr.eq.0) then
           coni0 = hybr(floi0,chci(ix,iy,0))-abs(floi0/2.0_R8)
           coni0_mdf = hybr(floi0_mdf,chci(ix,iy,0))-
     &                  abs(floi0_mdf/2.0_R8)                   !srv 18.08.08
         else
           coni0 = chci(ix,iy,0)
         endif
         if (flomx.gt.0.0_R8) then
           t0 = coni0*abs(dti)/
     &      (flomx*max(ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy)))
           t0_mdf = coni0*abs(dti)/                             !srv 18.08.08
     &      (flomx*max(ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy)))
         else
           t0 = 0.0_R8
           t0_mdf = 0.0_R8                                      !srv 18.08.08
         endif
         dti0 = dti/(1.0_R8+t0)
         t1 = t0/(1.0_R8+t0)
         dti0_mdf = dti/(1.0_R8+t0_mdf)             !srv 18.08.08 {
         t1_mdf = t0_mdf/(1.0_R8+t0_mdf)
	 if(mdf_fhi.eq.0) then
          fllimi(ix,iy) = 1.0_R8+t0
	 else
          fllimi(ix,iy) = 1.0_R8+t0_mdf
	 endif                                      !srv 18.08.08 }
         if (no_hybr.eq.0) then
           fhi(ix,iy,0) = w1 +
     .                    floi0*(ti(leftix(ix,iy),leftiy(ix,iy))
     &                          +ti(ix,iy))/2.0_R8-
     &                  abs(floi0/2.0_R8)*dti-coni0*dti0
           fhi_mdf(ix,iy,0) = w1 +
     .                    floi0_mdf*(ti(leftix(ix,iy),leftiy(ix,iy))
     &                          +ti(ix,iy))/2.0_R8-
     &                  abs(floi0_mdf/2.0_R8)*dti-coni0_mdf*dti0_mdf
     &                       - fhiPSch(ix,iy,0)                      !srv 13.10.06 }
         else
           t1 = dim ( coni0/fllimi(ix,iy),
     .         abs(floi0) *
     .         max(0.0_R8, 
     .             min(ti(ix,iy),ti(leftix(ix,iy),leftiy(ix,iy))))/
     .         max(cutlo, ti(ix,iy)+ti(leftix(ix,iy),leftiy(ix,iy))) )
           fhi(ix,iy,0) = w1 +
     .         upwind(floi0, ti(leftix(ix,iy),leftiy(ix,iy)), ti(ix,iy))
     .         - t1 * (ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))
         endif
         fhit(ix,iy,0) = fhi(ix,iy,0) +
     1    ((1.0_R8-BoRiS)*fni_32(ix,iy,0)+cum_flo_vaecrb)*
     &    (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
	 if(mdf_fhi.eq.0) then
          floi(ix,iy,0) = floi0-
     &    sign(flomx,ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))*
     &    min(0.0_R8,t1)**2
	 else
          floi(ix,iy,0) = floi0_mdf-
     &    sign(flomx,ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))*
     &    min(0.0_R8,t1_mdf)**2
	 endif
         if (no_hybr.eq.0) then
	  if(mdf_fhi.eq.0) then
           coni(ix,iy,0,0) = abs(floi0/2.0_R8)+coni0/(1.0_R8+t0)     !srv 02.01.07
	  else
           coni(ix,iy,0,0) = abs(floi0_mdf/2.0_R8)+
     &                       coni0_mdf/(1.0_R8+t0_mdf)     !srv 02.01.07
	  endif
         else
           coni(ix,iy,0,0) = coni0/(1.0_R8+t0)                       !srv 02.01.07
         endif
        endif
       enddo
      enddo
*   ..compute fhiy, floiy, coniy
      do iy = -1, ny
       do ix = -1, nx
        if(bottomiy(ix,iy).eq.-2) then
         floi(ix,iy,1) = 0.0_R8
         coni(ix,iy,1,0) = 0.0_R8                                    !srv 02.01.07
         fhi_mdf(ix,iy,1) = 0.0e0_R8                                 !srv 13.10.06
         fhi(ix,iy,1) = 0.0_R8
         fhit(ix,iy,1) = 0.0_R8
        else
         w1 = 0.0_R8
         if (conduction_only.eq.1) then
           floi0 = 0.0_R8
           floi0_mdf = 0.0_R8
           cum_flo_vadia=0.0_R8
           cum_flo_vaecrb=0.0_R8
         else
           if (BoRiS.eq.1.0_R8) then   ! Add the kinetic and viscous terms
             do is = 0, ns-1
               floi0 = kinrgy(ix,iy,is)*fna(ix,iy,1,is)
               coni0 = gs(ix,iy,1)*
     .          ave(vsal(ix,iy,1,is)/hy1(ix,iy),
     .              vsal(bottomix(ix,iy),bottomiy(ix,iy),1,is)/
     .              hy1(bottomix(ix,iy),bottomiy(ix,iy)))
               t1 = dim ( coni0,
     .                   abs(floi0) *
     .          max(0.0_R8, 
     .             min(ti(ix,iy),ti(bottomix(ix,iy),bottomiy(ix,iy))))/
     .          max(cutlo, 
     .                 ti(ix,iy)+ti(bottomix(ix,iy),bottomiy(ix,iy))))
               w1 = w1 +
     .          upwind(floi0, 
     .                ti(bottomix(ix,iy),bottomiy(ix,iy)), ti(ix,iy))
     .           - t1 * (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
             enddo
           endif
           floi0 = (1.5_R8+BoRiS)*fni_32(ix,iy,1)+
     1      2.5_R8*fni_52(ix,iy,1)+
     2      chvi(ix,iy,1)*(ni(bottomix(ix,iy),bottomiy(ix,iy),0)
     &                    +ni(ix,iy,0))/2.0_R8
           cum_flo_vadia=0.0_R8
           cum_flo_vaecrb=0.0_R8
           do is=0,ns-1
             prefactor=
     &        (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))    !sv 17.05.99
     &        /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))   !sv 17.05.99
     &        *(na(bottomix(ix,iy),bottomiy(ix,iy),is)             !sv 17.05.99
     &         +na(ix,iy,is))/2.0_R8                               !sv 17.05.99
             cum_flo_vadia=cum_flo_vadia+
     1        prefactor*vadia(ix,iy,1,is)
             cum_flo_vaecrb=cum_flo_vaecrb+
     1        prefactor*vaecrb(ix,iy,1,is)
           enddo
         endif
	 if(mdf_fhi.eq.0) then
          floi(ix,iy,1) = floi0+
     1     2.5_R8*cum_flo_vadia+1.5_R8*cum_flo_vaecrb
	 else
          floi(ix,iy,1) = floi0+1.5_R8*cum_flo_vaecrb
	 endif
         if (no_hybr.eq.0) then
           coni(ix,iy,1,0) = hybr(floi(ix,iy,1),chci(ix,iy,1))    !srv 02.01.07
           fhi(ix,iy,1) = w1 + floi(ix,iy,1)*
     1      (ti(bottomix(ix,iy),bottomiy(ix,iy))
     &                                +ti(ix,iy))/2.0_R8-
     &       coni(ix,iy,1,0)*                                     !srv 02.01.07
     1      (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
           fhi_mdf(ix,iy,1) = w1 + floi(ix,iy,1)*
     1      (ti(bottomix(ix,iy),bottomiy(ix,iy))
     &                                +ti(ix,iy))/2.0_R8-
     &       coni(ix,iy,1,0)*                                     !srv 02.01.07
     1      (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
     &      -fhiPSch(ix,iy,1)                                     !sv 13.10.06 }
         else
           coni(ix,iy,1,0) = chci(ix,iy,1)                        !srv 02.01.07
           t1 = dim ( coni(ix,iy,1,0),                            !srv 02.01.07
     .         abs(floi(ix,iy,1)) *
     .         max(0.0_R8, 
     .             min(ti(ix,iy),ti(bottomix(ix,iy),bottomiy(ix,iy))))/
     .         max(cutlo, 
     .                 ti(ix,iy)+ti(bottomix(ix,iy),bottomiy(ix,iy))))
           fhi(ix,iy,1) = w1 + 
     .         upwind(floi(ix,iy,1), 
     .                ti(bottomix(ix,iy),bottomiy(ix,iy)), ti(ix,iy))
     .         - t1 * (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
         endif
         fhit(ix,iy,1) = fhi(ix,iy,1)+
     1    ((1.0_R8-BoRiS)*fni_32(ix,iy,1)+cum_flo_vaecrb)*
     &    (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
        endif
       enddo
      enddo

      if(guard_flows.eq.0) then                                    !xpb
cxpb turn off flows between boundary cells
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
              fhi(ix,iy,1)=0.0_R8
              floi(ix,iy,1)=0.0_R8
            endif
          enddo
        enddo

cxpb turn off poloidal flows in core boundary
        do ix = -1, nx
          if(mod(region(ix,-1,0),4).eq.1) then
            fhi(ix,-1,0) = 0.0_R8
            floi(ix,-1,0) = 0.0_R8
          endif
        enddo
      endif                                                        !xpb

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fhi(ix,iy,0)=0.0_R8
            floi(ix,iy,0)=0.0_R8
            fhi(ix,iy,1)=0.0_R8
            floi(ix,iy,1)=0.0_R8
          end do
        end do
      end do

      if(iout.ne.0) then                                             !srv 18.06.08 {
        call my_out(70,nx,ny,floi(-1,-1,0),'b2tfhi_floix')
        call my_out(70,nx,ny,floi(-1,-1,1),'b2tfhi_floiy')
        call my_out(70,nx,ny,coni(-1,-1,0,0),'b2tfhi_conix0')
        call my_out(70,nx,ny,coni(-1,-1,0,1),'b2tfhi_conix1')
        call my_out(70,nx,ny,coni(-1,-1,0,2),'b2tfhi_conix2')
        call my_out(70,nx,ny,coni(-1,-1,1,0),'b2tfhi_coniy0')
        call my_out(70,nx,ny,coni(-1,-1,1,1),'b2tfhi_coniy1')
        call my_out(70,nx,ny,coni(-1,-1,1,2),'b2tfhi_coniy2')
        call my_out(70,nx,ny,fhi(-1,-1,0),'b2tfhi_fhix')
        call my_out(70,nx,ny,fhi(-1,-1,1),'b2tfhi_fhiy')
        call my_out(70,nx,ny,fhi_mdf(-1,-1,0),'b2tfhi_fhi_mdfx')
        call my_out(70,nx,ny,fhi_mdf(-1,-1,1),'b2tfhi_fhi_mdfy')
        call my_out(70,nx,ny,fhiPSch(-1,-1,0),'b2tfhi_fhiPSchx')
        call my_out(70,nx,ny,fhiPSch(-1,-1,1),'b2tfhi_fhiPSchy')
      endif                                                          !srv 18.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

*!!   Reconsider factor 2.5_R8 versus 1.5_R8 in transport and linearisation.
cxpb  Factor changed to 1.5_R8.
cxpb  This version of the routine uses a (1.5+BoRiS) factor for the
cxpb  poloidal convective heat flux instead of the usual 3/2 or 5/2 
cxpb  so as to solve the same equation as in BoRiS for benchmark purposes. 
cxpb  The default is BoRiS = 0.0.
cxpb  For comparison, set BoRiS = 1.0.
cxpb  Routine rewritten to explicitly separate 3/2 and 5/2 terms.
cxpb  The contributions from fni_32 and fni_52 are NO DRIFTS
cxpb  contributions, i.e. diamagnetic and ExB terms are treated
cxpb  separately.
cxpb  The BoRiS treatment must be completed with addition of the 
cxpb  viscous and kinetic terms, as per the SOLPS4.0 treatment.

*-----------------------------------------------------------------------
*.end b2tfhi

      end subroutine b2tfhi
