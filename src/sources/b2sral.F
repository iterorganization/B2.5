*-----------------------------------------------------------------------
*.specification

      subroutine b2sral (nx, ny, ns, nscx, nscxmax, iscx, ismain,        !xpb
     &  dtim, BoRiS, facdrift, fac_ExB,                                  !xpb
     &  vol, crx, hx, hy, hz, qz, qc, gs, pbs, pbshz, bb,                !srv 09.01.01
     &  na, ua, uadia, vedia, vadia, wadia, veecrb, vaecrb, ve, wedia,
     &  te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,                  !srv 28.10.02
     &  fna, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52, 
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,        !srv 13.10.06 25.12.06 23.10.09
     &  fhe, fhi, fnaPSch,fhePSch, fhiPSch, fch, 
     &  fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,   !srv 21.10.02
     &  fchanml, fna_eir, fne_eir, fhe_eir, fhi_eir, 
     &  cdna, cvsa, cvla, chce, chve, chci, chvi, calf,   !xpb
     &  rlsa, rlra, rlqa, rlcx, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &  rza, rz2, rpt, rpi,
     &  sna, smo, smq, she, shi, sch,
     &  wrong_flow, main_call)						 !xpb, dpc
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_tallies
      use b2mod_sources 
     &     , only: rcxna, rcxhi, rcxmo
      use b2mod_input_profile
      use b2mod_anomalous_transport
      use b2mod_b2cmpa
      use b2mod_b2cmfs
#ifdef ASTRA
      use b2mod_b2_to_astra
#endif	                       
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nscx, nscxmax, iscx(0:nscxmax-1), ismain       !xpb
      real (kind=R8) ::                                                  !xpb
     *  dtim, BoRiS, facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny)         !xpb
*     (geometry)
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), crx(-1:nx,-1:ny,0:3), hx(-1:nx,-1:ny), 
     &  hy(-1:nx,-1:ny), hz(-1:nx,-1:ny), pbshz(-1:nx,-1:ny,0:1),        !srv 09.01.01
     &  qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny), pbs(-1:nx,-1:ny,0:1),
     1  bb(-1:nx,-1:ny,0:3), gs(-1:nx,-1:ny,0:2)
*     (state)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1), ve(-1:nx,-1:ny),
     1  uadia(-1:nx,-1:ny,0:1,0:ns-1), vadia(-1:nx,-1:ny,0:1,0:ns-1),
     2  vedia(-1:nx,-1:ny,0:1), veecrb(-1:nx,-1:ny,0:1), 
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), 
     &  kinrgy(-1:nx,-1:ny,0:ns-1), wedia(-1:nx,-1:ny,0:1),
     &  wadia(-1:nx,-1:ny,0:1,0:ns-1), vaecrb(-1:nx,-1:ny,0:1,0:ns-1)
*     (fluxes)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1), fna_32(-1:nx,-1:ny,0:1,0:ns-1), 
     &  fna_52(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     & fch_p(-1:nx,-1:ny,0:1), fchvispar(-1:nx,-1:ny,0:1),               !srv 11.10.02 21.10.02 {
     & fchvisper(-1:nx,-1:ny,0:1), fchin(-1:nx,-1:ny,0:1),
     & fchvisq(-1:nx,-1:ny,0:1), fchdia(-1:nx,-1:ny,0:1),
     & fchinert(-1:nx,-1:ny,0:1), fchanml(-1:nx,-1:ny,0:1),              !srv 21.10.02 }
     *  fni_32(-1:nx,-1:ny,0:1), fni_52(-1:nx,-1:ny,0:1),
     *  fne_32(-1:nx,-1:ny,0:1), fne_52(-1:nx,-1:ny,0:1),
     &  floe_noc(-1:nx,-1:ny,0:1), floi_noc(-1:nx,-1:ny,0:1),            !srv 28.10.02
     & fna_fcor(-1:nx,-1:ny,0:1,0:ns-1), fna_he(-1:nx,-1:ny,0:1,0:ns-1), !srv 25.02.00 24.10.02
     & fna_nodrift(-1:nx,-1:ny,0:1,0:ns-1),
     & fnaPSch(-1:nx,-1:ny,0:1,0:ns-1), fna_mdf(-1:nx,-1:ny,0:1,0:ns-1), !srv 13.10.06 25.12.06
     & fhePSch(-1:nx,-1:ny,0:1), fhiPSch(-1:nx,-1:ny,0:1),               !srv 13.10.06
     & fhe_mdf(-1:nx,-1:ny,0:1), fhi_mdf(-1:nx,-1:ny,0:1),                !srv 23.10.09
     & fna_eir(-1:nx,-1:ny,0:1, 0:ns-1), fne_eir(-1:nx,-1:ny,0:1),
     & fhe_eir(-1:nx,-1:ny,0:1), fhi_eir(-1:nx,-1:ny,0:1)
*     (transport coefficients)						 !xpb
      real (kind=R8) ::							 !xpb
     1  chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),			 !xpb
     2  chci(-1:nx,-1:ny,0:1), chvi(-1:nx,-1:ny,0:1),			 !xpb
     3  calf(-1:nx,-1:ny,0:1), cdna(-1:nx,-1:ny,0:1,0:ns-1),         !xpb
     4  cvsa(-1:nx,-1:ny,0:1,0:ns-1), cvla(-1:nx,-1:ny,0:1,0:ns-1),
     5  cdpa(-1:nx,-1:ny,0:1,0:ns-1)
*     (log-log rate coefficients)
      real (kind=R8) ::
     *  rlsa(-1:nx,-1:ny,0:1,0:ns-1), rlra(-1:nx,-1:ny,0:1,0:ns-1),
     &  rlqa(-1:nx,-1:ny,0:1,0:ns-1), rlrd(-1:nx,-1:ny,0:1,0:ns-1),
     &  rlbr(-1:nx,-1:ny,0:1,0:ns-1),
     &  rlcx(-1:nx,-1:ny,0:1,0:ns-1,0:nscxmax-1),
     &  rlza(-1:nx,-1:ny,0:1,0:ns-1), rlz2(-1:nx,-1:ny,0:1,0:ns-1),
     &  rlpt(-1:nx,-1:ny,0:1,0:ns-1), rlpi(-1:nx,-1:ny,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  rza(-1:nx,-1:ny,0:ns-1), rz2(-1:nx,-1:ny,0:ns-1), 
     &  rpt(-1:nx,-1:ny,0:ns-1), rpi(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *  sna(-1:nx,-1:ny,0:1,0:ns-1), smo(-1:nx,-1:ny,0:3,0:ns-1),
     &  smq(-1:nx,-1:ny,0:3,0:ns-1), she(-1:nx,-1:ny,0:3),
     &  shi(-1:nx,-1:ny,0:3), sch(-1:nx,-1:ny,0:3)
      logical wrong_flow						 !xpb
      logical main_call                                                  !dpc
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SRAL provides linearised expressions for the sources of
*     density, momentum, electron heat, atom heat, and charge.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, k, is, ix, iy, sources_inputfile
      integer guard_flows
      integer, save :: style = 0                                 !srv 28.06.08
	integer, save :: iout = 0                                  !srv 30.06.08
	integer, save :: use_eirene =0, save_f30 = 0, save_f31 = 0 !srv 14.07.10
      integer, save :: anomalous_core_only = 0                   !srv 19.04.07
      integer, save :: no_fluid_neutral_sources = 0
 	character chns*2, chk*1                                    !srv 30.06.08
      real (kind=R8) ::
!     *  snatot(0:ns-1),                                         !srv 18.05.09
     1  rsa(-1:nx,-1:ny,0:ns-1), wrk1(-1:nx,-1:ny),
     2  wrk0(-1:nx,-1:ny), rra(-1:nx,-1:ny,0:ns-1),
     3  rcx(-1:nx,-1:ny,0:ns-1),
     4  smq0(-1:nx,-1:ny,0:3,0:ns-1), smq0_ion(-1:nx,-1:ny,0:3,0:ns-1), !srv 11.09.09
     4  smq0_rec(-1:nx,-1:ny,0:3,0:ns-1),
     5  rqa(-1:nx,-1:ny,0:ns-1), rrd(-1:nx,-1:ny,0:ns-1), 
     6  rbr(-1:nx,-1:ny,0:ns-1), fni(-1:nx,-1:ny,0:1),
     7  fne(-1:nx,-1:ny,0:1),
     8  sna0(-1:nx,-1:ny,0:1,0:ns-1), sna0_ion(-1:nx,-1:ny,0:1,0:ns-1), !srv 11.09.09
     8  sna0_rec(-1:nx,-1:ny,0:1,0:ns-1),
     9  shi0(-1:nx,-1:ny,0:3), shi0_ion(-1:nx,-1:ny,0:3),
     &  shi0_rec(-1:nx,-1:ny,0:3), she0(-1:nx,-1:ny,0:3),
     &  smo0(-1:nx,-1:ny,0:3,0:ns-1), sch0(-1:nx,-1:ny,0:3),
     1  wpo(-1:nx,-1:ny), wpb(-1:nx,-1:ny), pb(-1:nx,-1:ny,0:ns-1)
      save ncall, sources_inputfile, guard_flows
*   ..procedures
      real (kind=R8) ::
     *  smin, smax, hy1
      external subini, subend, xertst, sfill, saxpy, smin, smax
      external b2xvsg, b2xvff, b2xvps, b2xpfe, b2xpfi, b2sqel, b2sqcx,
     &  b2stbc, b2stbm, b2stbr, b2stel, b2stcx, ipgeti, b2xzdd
      hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)
*   ..initialisation
      data ncall/0/, sources_inputfile/0/, guard_flows/1/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sral')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns, 
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
*   ..test facdrift (not on first call)
      if (ncall.gt.0) then
        call xertst (0.0_R8.le.minval(facdrift) .and.
     &       1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
        call xertst (0.0_R8.le.minval(fac_ExB) .and.
     &       1.0_R8.ge.maxval(fac_ExB), 'faulty argument fac_ExB')
       endif
*   n.inputfile for profile sources
      if (ncall.eq.0) then
       call ipgeti('b2sral_inputfile', sources_inputfile)
       call ipgeti('b2sral_style', style)                               !srv 28.07.08
       call ipgeti('b2sral_iout', iout)                                 !srv 30.06.08
       call ipgeti ('b2tfnb_anomalous_core_only', anomalous_core_only)  !srv 19.04.07
       call ipgeti ('b2mndr_eirene',use_eirene)                         !srv 14.07.10 {
       call ipgeti ('eirene_savef30', save_f30)
       call ipgeti ('eirene_savef31', save_f31)                         !srv 14.07.10 }
       call ipgeti ('b2sral_no_fluid_neutral_sources', 
     .                                       no_fluid_neutral_sources)
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol, hx, hy
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
*    ..test range of qz
       call xertst (-1.0_R8.lt.smin(n2,qz(-1,-1,0),1).and.
     &   smax(n2,qz(-1,-1,0),1).lt.1.0_R8, 
     &   'faulty argument range qz(,,0)')
       call xertst (0.0_R8.lt.smin(n2,qz(-1,-1,1),1).and.
     &   smax(n2,qz(-1,-1,1),1).le.1.0_R8, 
     &   'faulty argument range qz(,,1)')
*    ..test state
       call b2xvps (nx, ny, ns, 
     &   ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
       call b2xvsg (n2*ns, kinrgy, 1, 'kinrgy', '.ge.')
      endif
*   ..compute fne, fni
      if(style.eq.0) then                                               !srv 28.07.08
       call b2xpfe (nx, ny, ns, qe, rza, fna, fch, fne)
	else
       call b2xpfe (nx, ny, ns, qe, rza, fna_nodrift, fch_p, fne)       !srv 11.10.02
	endif
      call b2xpfi (nx, ny, ns, fna, fni)

*   ..test facdrift and fac_ExB
      if (ncall.eq.0) then                                               !xpb
       call ipgeti ('b2news_guard_flows', guard_flows)                   !xpb
      endif                                                              !xpb
      call xertst (0.0_R8.le.minval(facdrift) .and.
     &     1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
      call xertst (0.0_R8.le.minval(fac_ExB) .and.
     &     1.0_R8.ge.maxval(fac_ExB), 'faulty argument fac_ExB')
	
      if(style.eq.0) then                                               !srv 28.07.08
*   ..compute temporary drift velocities on first call
       if (ncall.eq.0) then
         call intfacev(nx,ny,bottomix,bottomiy,vol,po,wpo)
         do is = 0, ns-1
           pb(-1:nx,-1:ny,is) = na(-1:nx,-1:ny,is)*ti(-1:nx,-1:ny)
           call intfacev(nx,ny,bottomix,bottomiy,vol,pb(-1,-1,is),wpb)
           do ix = -1, nx
             do iy = -1, ny
               if(topiy(topix(ix,iy),min(topiy(ix,iy),ny)).ne.iy .and.
     &            rightix(ix,iy).ne.nx+1 .and. leftix(ix,iy).ne.-2 .and.
     &            topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &            region(ix,iy,0).ne.0 .and. .not.is_neutral(is)) then
                vaecrb(ix,iy,0,is) = - 1.0_R8/hy1(ix,iy)*  ! signmf !srv 11.01.13
     1            1.0_R8/bb(ix,iy,3)*  
     2            (wpo(topix(ix,iy),topiy(ix,iy))-wpo(ix,iy))*
     &            fac_ExB(ix,iy)
                wadia(ix,iy,0,is) = - bb(ix,iy,2)/bb(ix,iy,3)/hy1(ix,iy)
     1            /bb(ix,iy,3)/qe/rza(ix,iy,is)/na(ix,iy,is)*
     2            (wpb(topix(ix,iy),topiy(ix,iy))-wpb(ix,iy))*
     3            facdrift(ix,iy)
               else
                vaecrb(ix,iy,0,is) = 0.0_R8
                wadia(ix,iy,0,is) = 0.0_R8
               endif
             enddo
           enddo
         enddo
         call intfaceh(nx,ny,leftix,leftiy,vol,po,wpo)
         do is = 0, ns-1
           call intfaceh(nx,ny,leftix,leftiy,vol,pb(-1,-1,is),wpb)
           do ix = -1, nx
             do iy = -1, ny                                               !xpb
               if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &            bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &            region(ix,iy,0).ne.0 .and. .not.is_neutral(is)) then
                 vaecrb(ix,iy,1,is) = 1.0_R8/hx(ix,iy)*                   !xpb
     1             bb(ix,iy,2)/bb(ix,iy,3)**2*                            !xpb
     2            (wpo(rightix(ix,iy),rightiy(ix,iy))-wpo(ix,iy))*
     &            fac_ExB(ix,iy)
                 wadia(ix,iy,1,is) = bb(ix,iy,2)/
     &            (qe*rza(ix,iy,is)*hx(ix,iy)*na(ix,iy,is))/
     1             bb(ix,iy,3)**2*
     2            (wpb(rightix(ix,iy),rightiy(ix,iy))-wpb(ix,iy))*
     3             facdrift(ix,iy)
               else
                 vaecrb(ix,iy,1,is) = 0.0_R8
                 wadia(ix,iy,1,is) = 0.0_R8
               endif
             enddo
           enddo
         enddo
* drift velocities at walls are extrapolated from those at neighbouring cells
* drift velocities at iy=ny are extrapolated from those at iy=ny-1       !xpb
* drift velocities at iy=-1 are extrapolated from those at iy=0          !xpb
         do ix = -1, nx                                                   !xpb
           do iy = -1, ny                                                 !xpb
             if(leftix(ix,iy).eq.-2) then
               vaecrb(ix,iy,:,:) = 
     .          vaecrb(rightix(ix,iy),rightiy(ix,iy),:,:)                 !xpb
               wadia(ix,iy,:,:) = 
     .          wadia(rightix(ix,iy),rightiy(ix,iy),:,:)                  !xpb
             elseif(rightix(ix,iy).eq.nx+1) then                          !xpb
               vaecrb(ix,iy,:,:) = 
     .          vaecrb(leftix(ix,iy),leftiy(ix,iy),:,:)                   !xpb
               wadia(ix,iy,:,:) = wadia(leftix(ix,iy),leftiy(ix,iy),:,:)  !xpb
             endif
             if(topiy(ix,iy).eq.ny+1) then                                !xpb
               vaecrb(ix,iy,:,:) = 
     .          vaecrb(bottomix(ix,iy),bottomiy(ix,iy),:,:)               !xpb
               wadia(ix,iy,:,:) = 
     .          wadia(bottomix(ix,iy),bottomiy(ix,iy),:,:)                !xpb
             elseif(bottomiy(ix,iy).eq.-2) then
               vaecrb(ix,iy,:,:) = vaecrb(topix(ix,iy),topiy(ix,iy),:,:)  !xpb
               wadia(ix,iy,:,:) = wadia(topix(ix,iy),topiy(ix,iy),:,:)    !xpb
             endif
           enddo
         enddo
* turn off side drift flows in guard cells                               !xpb
         do iy = -1, ny                                                   !xpb
          do ix = -1, nx
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            vaecrb(ix,iy,1,:) = 0.0_R8                                    !xpb
            wadia(ix,iy,1,:) = 0.0_R8                                     !xpb
           endif
           if(bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1) then
            vaecrb(ix,iy,0,:) = 0.0_R8                                    !xpb
            wadia(ix,iy,0,:) = 0.0_R8                                     !xpb
           endif
          enddo
         enddo                                                            !xpb 
         if(guard_flows.eq.0) then                                        !xpb
cxpb turn off flows between boundary cells
           do iy = -1, ny
            do ix = -1, nx
             if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
              vaecrb(ix,iy,1,:)=0.0_R8
              wadia(ix,iy,1,:)=0.0_R8
             endif
            enddo
           enddo
cxpb turn off poloidal flows in core boundary
           do ix = -1, nx
            if(mod(region(ix,-1,0),4).eq.1) then
             vaecrb(ix,-1,0,:) = 0.0_R8
             wadia(ix,-1,0,:) = 0.0_R8
            endif
           enddo
         endif                                                            !xpb
       endif
	endif                                                               !srv 28.07.08
      if(use_eirene.ne.0 .or. save_f31.ne.0) then                       !srv 14.07.10 {
*   ..recompute currents
       if((maxval(facdrift).ne.0.0_R8 .or.
     &    maxval(fac_ExB).ne.0.0_R8)) then
        call b2tvspa (nx, ny, ns, vol, hx, hy, qz, bb, ua, cvsa,        !srv 20.09.06  {
     &                       facdrift, fchvispar)
	  call b2tiner (nx, ny, ns, vol,
     &                       gs, hx, hy, bb, na, ua, facdrift, fchinert)
        call b2tanml (nx, ny, ns, vol, hx, hy, csig_an, po,
     &                       fchanml)                                   !srv 20.09.06  }
	  if(anomalous_core_only.ne.0) then                               !srv 19.04.07  {
	     do iy = 0, ny                                                !lk 19.09.07
             do ix = -1, nx
     			 if(mod(region(ix,bottomiy(ix,iy),0),4).ne.1) then    !lk 19.09.07
                 fchanml(ix,iy,1) = 0.0e0_R8
               endif 
		   enddo
	     enddo
	  endif                                                           !srv 19.04.07  }
	 else
        fchvispar = 0.0_R8
        fchinert  = 0.0_R8
        fchanml   = 0.0_R8
       endif

csw 30sep2013, added cdpa in call, uadia
      call b2tfrn (nx, ny, ns, guard_flows,
     &  facdrift, fac_ExB, rza,
     &  vol, hx, hy, hz, qz, gs, qc, bb, pbs, ua, na, ne, te, ti, po,   
     &  fch_p, fchvispar, fchinert, fchanml,
     &  fne, cvla, cdna, cdpa, chce, chci, vaecrb, wadia, ve, wedia,
     &  uadia,
     &  fna_eir, fne_eir, fhe_eir, fhi_eir)
      endif                                                         !srv 14.07.10 }

* ..compute source coefficients
*   ..compute standard contributions due to boundaries
      wrong_flow = .false.						 !xpb
      if(style.eq.0) then                                               !srv 28.07.08
       call b2stbc (nx, ny, ns, ismain, BoRiS, facdrift, fac_ExB,         !xpb
     &  vol, crx, hx, hy, hz, qz, qc, gs, pbs, pbshz, bb,                 !xpb !srv 09.01.01
     &  na, ua, uadia, vedia, vadia, veecrb, vaecrb, 
     .  te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,
     1  fna, fna_mdf, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52, 
     &  fhe, fhi, fhe_mdf, fhi_mdf, fch, fne, fni, fna_fcor,            !srv 16.02.12
     &  fnaPSch, fhePSch, fhiPSch,                                      !srv 25.02.00 13.10.06
     &  fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,   !srv 21.10.02
     2  cdna, chce, chve, chci, chvi, calf,				               !xpb
     &  sna, smo, she, shi, sch, wrong_flow, main_call)			 !xpb
	elseif(style.eq.1) then                                           !srv 28.07.08 {
       call b2stbc_spb (nx, ny, ns, facdrift, vol, gs, hx, hy, hz, qz,  !srv 09.01.01
     &  pbs, pbshz, bb, na, ua, uadia, vaecrb, vadia, te, ti, po, ne,   !srv 14.10.99 12.05.00 09.01.01 28.10.02
     1  ni, floe_noc, floi_noc, fna, fna_fcor, fnaPSch, fhe, fhi,       !srv 13.10.06
     &  fhePSch, fhiPSch, fne, fni, fch,                                !srv 25.02.00 21.10.02
     &  fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,  !srv 21.10.02
     &  sna, smo, she, shi, sch, wrong_flow)                            !srv 09.07.99 01.07.05
       call b2xpfe (nx, ny, ns, qe, rza, fna, fch_p, fne)               !srv 25.12.06
      else
       if(ncall.eq.0) write(*,*) ' style=',style,' b2stbc called'
       call b2stbc (nx, ny, ns, ismain, BoRiS, facdrift, fac_ExB,         !xpb
     &  vol, crx, hx, hy, hz, qz, qc, gs, pbs, pbshz, bb,                 !xpb !srv 09.01.01
     &  na, ua, uadia, vedia, vadia, veecrb, vaecrb, 
     .  te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,
     1  fna, fna_mdf, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52, 
     &  fhe, fhi, fhe_mdf, fhi_mdf, fch, fne, fni, fna_fcor,             !srv 16.02.12
     &  fnaPSch, fhePSch, fhiPSch,                                       !srv 25.02.00 13.10.06
     &  fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,   !srv 21.10.02
     2  cdna, chce, chve, chci, chvi, calf,					             !xpb
     &  sna, smo, she, shi, sch, wrong_flow, main_call)			 !xpb
       call b2xpfe (nx, ny, ns, qe, rza, fna, fch_p, fne)               !srv 25.12.06
	endif                                              !srv 28.07.08 }
#ifdef ASTRA
      se  =0.0_R8                                         !srv 28.07.11
      seh =0.0_R8                                         !srv 11.06.13
      sih =0.0_R8                                         !srv 11.06.13
	sheNB_astra=0.0_R8                                  !srv 21.06.12 {
	shiNB_astra=0.0_R8
      do iy = -1,ny
	 do ix = -1,nx
        if(mod(region(ix,iy,0),4).eq.1) then                            !srv 23.10.12
	   sheNb_astra(ix,iy,0)=vol(ix,iy)*PeNB_astra(iy)*1.0e6_R8
	   shiNb_astra(ix,iy,0)=vol(ix,iy)*PiNB_astra(iy)*1.0e6_R8
	  endif
	 enddo
	enddo
      call saxpy (n2*4, 1.0_R8, sheNb_astra, 1, she, 1)
      call saxpy (n2*4, 1.0_R8, shiNb_astra, 1, shi, 1)  !srv 21.06.12 }
#endif	                       
*   ..initialise smq to 0
      call sfill (n2*4*ns, 0.0_R8, smq, 1)
*   ..compute special contributions (boundaries and recycling)
      call b2stbr (nx, ny, ns, nscx, nscxmax, iscx, ismain,              !srv 25.10.08
     &  dtim, BoRiS, vol, hx, hy, gs, qz, pbs, na, ua, kinrgy, te, ti,
     &  po, ne, ni, fna, fhe, fhi, fch, fne, fni,                        !srv 29.07.08
     &  sna0, smo0, she0, shi0, sch0, wrk1, main_call)                   !dpc
*   ..zero out contributions in dead cells
      do is = 0, ns-1
        call b2xzdd (nx, ny, 1, sna0(-1,-1,0,is))
        call b2xzdd (nx, ny, 3, smo0(-1,-1,0,is))
      enddo
      call b2xzdd (nx, ny, 3, she0)
      call b2xzdd (nx, ny, 3, shi0)
      call b2xzdd (nx, ny, 3, sch0)
*   ..add to previous contributions
      call saxpy (n2*2*ns, 1.0_R8, sna0, 1, sna, 1)
      call saxpy (n2*4*ns, 1.0_R8, smo0, 1, smo, 1)
      call saxpy (n2*4, 1.0_R8, she0, 1, she, 1)
      call saxpy (n2*4, 1.0_R8, shi0, 1, shi, 1)
      call saxpy (n2*4, 1.0_R8, sch0, 1, sch, 1)
!      do is=0,ns-1                                                          !srv 18.05.09 {
!        snatot(is)=0.0_R8
!        do ix=-1,nx
!          do iy=-1,ny
!            snatot(is)=snatot(is)+sna(ix,iy,0,is)+
!     &	     sna(ix,iy,1,is)*na(ix,iy,is)
!          enddo
!        enddo
!      enddo
!      if(main_call.or.wrong_flow)                                       !xpb !srv 18.05.09 }
!     > write(*,*) 'snatot',(snatot(is),is=0,ns-1)
*   ..compute more special contributions
      call b2stbm (nx, ny, ns, vol, hx, hy, qz, pbs, na, ua, te, ti,
     &  po, ne, ni, fna, fhe, fhi, fch, fne, fni,
     &  sna0, smo0, she0, shi0, sch0,
     &  wrk0, wrk1)
*   ..zero out contributions in dead cells and core boundary
      do is = 0, ns-1
        call b2xzdd (nx, ny, 1, sna0(-1,-1,0,is))
        call b2xzdd (nx, ny, 3, smo0(-1,-1,0,is))
      enddo
      call b2xzdd (nx, ny, 3, she0)
      call b2xzdd (nx, ny, 3, shi0)
      call b2xzdd (nx, ny, 3, sch0)
*   ..add to previous contributions
      call saxpy (n2*2*ns, 1.0_R8, sna0, 1, sna, 1)
      call saxpy (n2*4*ns, 1.0_R8, smo0, 1, smo, 1)
      call saxpy (n2*4, 1.0_R8, she0, 1, she, 1)
      call saxpy (n2*4, 1.0_R8, shi0, 1, shi, 1)
      call saxpy (n2*4, 1.0_R8, sch0, 1, sch, 1)
*   ..compute sources by inputfiles (heating,...)
      if (sources_inputfile.eq.1) then
        call source_input(nx,ny,ns,sna0,smo0,she0,shi0,sch0,vol,ncall)
*   ..zero out contributions in dead cells and core boundary
        do is = 0, ns-1
          call b2xzdd (nx, ny, 1, sna0(-1,-1,0,is))
          call b2xzdd (nx, ny, 3, smo0(-1,-1,0,is))
        enddo
        call b2xzdd (nx, ny, 3, she0)
        call b2xzdd (nx, ny, 3, shi0)
        call b2xzdd (nx, ny, 3, sch0)
*   ..add to previous contributions
        call saxpy (n2*2*ns, 1.0_R8, sna0, 1, sna, 1)
        call saxpy (n2*4*ns, 1.0_R8, smo0, 1, smo, 1)
        call saxpy (n2*4, 1.0_R8, she0, 1, she, 1)
        call saxpy (n2*4, 1.0_R8, shi0, 1, shi, 1)
        call saxpy (n2*4, 1.0_R8, sch0, 1, sch, 1)
      endif

*   ..compute electron rate coefficients
      call b2sqel (nx, ny, ns, ismain, ev, te, 
     &  rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &  rsa, rra, rqa, rrd, rbr, rza, rz2, rpt, rpi, wrk0)
*   ..compute sources due to electron-atom processes
*     (momentum source must be pure friction)
      call b2stel (nx, ny, ns, ismain, BoRiS, bb, vol, hx, hy, hz, na,     !srv 09.01.01
     & ua, wadia, vadia, vaecrb, te, ti, po, ne, ni, facdrift, fac_ExB,   !srv 07.07.09
     & rsa, rra, rqa, rrd, rbr, rza, rz2, rpt, rpi, rlsa, rlra, rlqa, 
     & sna0_ion, sna0_rec, smq0_ion, smq0_rec, she0, shi0_ion, shi0_rec) !srv 11.09.09 {
*   ..zero out contributions in dead cells and core boundary
      do is = 0, ns-1
        call b2xzdd (nx, ny, 1, sna0_ion(-1,-1,0,is))
        call b2xzdd (nx, ny, 1, sna0_rec(-1,-1,0,is))
        call b2xzdd (nx, ny, 3, smq0_ion(-1,-1,0,is))
        call b2xzdd (nx, ny, 3, smq0_rec(-1,-1,0,is))
      enddo
      call b2xzdd (nx, ny, 3, she0)
      call b2xzdd (nx, ny, 3, shi0_ion)
      call b2xzdd (nx, ny, 3, shi0_rec)
*   ..add to previous contributions
      call saxpy (n2*2*ns, 1.0_R8, sna0_ion, 1, sna, 1)
      call saxpy (n2*2*ns, 1.0_R8, sna0_rec, 1, sna, 1)
      call saxpy (n2*4*ns, 1.0_R8, smq0_ion, 1, smq, 1)
      call saxpy (n2*4*ns, 1.0_R8, smq0_rec, 1, smq, 1)
      call saxpy (n2*4, 1.0_R8, she0, 1, she, 1)
      call saxpy (n2*4, 1.0_R8, shi0_ion, 1, shi, 1)
      call saxpy (n2*4, 1.0_R8, shi0_rec, 1, shi, 1)                     !srv 11.09.09 }
*     (no contribution to sch)
*   ..compute sources due to charge exchange
      rcxnareg=0.0_R8
      rcxhireg=0.0_R8
      rcxmoreg=0.0_R8
      rcxna=0.0_R8
      rcxhi=0.0_R8
      rcxmo=0.0_R8
      do k = 0, nscx-1
*    ..compute charge-exchange rate coefficients
        call b2sqcx (nx, ny, ns, ev, am(iscx(k)), ti, rlcx(-1,-1,0,0,k),
     &   rcx, wrk0)
*    ..compute sources due to charge exchange
*      (momentum source must be pure friction)
        call b2stcx (nx, ny, ns, iscx(k), ismain, hz, BoRiS, vol, na,     !srv 09.01.01
     &    ua, wadia, vadia, vaecrb, ti, ni, rcx, sna0, smq0, shi0)       !srv 07.07.09
*   ..zero out contributions in dead cells and core boundary
        do is = 0, ns-1
         call b2xzdd (nx, ny, 1, sna0(-1,-1,0,is))
         call b2xzdd (nx, ny, 3, smq0(-1,-1,0,is))
        enddo
        call b2xzdd (nx, ny, 3, shi0)
*    ..add to previous contributions
        call saxpy (n2*2*ns, 1.0_R8, sna0, 1, sna, 1)
        call saxpy (n2*4*ns, 1.0_R8, smq0, 1, smq, 1)
        call saxpy (n2*4, 1.0_R8, shi0, 1, shi, 1)
*      (no contribution to she or to sch)
      enddo

      if(iout.ne.0) then                                            !srv 30.06.08 {
       do is = 0,ns-1                              
	   write (chns,'(i2.2)') is                          
	   do k=0,1
	    write (chk,'(i1)') k 
          call my_out(70,nx,ny,sna(-1,-1,k,is),'b2sral_sna'//chk//chns)
         enddo
 	   do k=0,3
	    write (chk,'(i1)') k 
          call my_out(70,nx,ny,smo(-1,-1,k,is),'b2sral_smo'//chk//chns)
          call my_out(70,nx,ny,smq(-1,-1,k,is),'b2sral_smq'//chk//chns)
	   enddo
	 enddo 
  	 do k=0,3
	   write (chk,'(i1)') k 
         call my_out(70,nx,ny,she(-1,-1,k),'b2sral_she'//chk)
         call my_out(70,nx,ny,shi(-1,-1,k),'b2sral_shi'//chk)
         call my_out(70,nx,ny,sch(-1,-1,k),'b2sral_sch'//chk)
#ifdef ASTRA
         call my_out(70,nx,ny,sheNb_astra(-1,-1,k),
     &        'b2sral_sheNb_astra'//chk)
         call my_out(70,nx,ny,shiNb_astra(-1,-1,k),
     &        'b2sral_shiNb_astra'//chk)
#endif
	 enddo
	endif                                                         !srv 30.06.08 }

csw 17oct2012, erase sna, smo, smq, (she, shi, sch) for fluid neutrals
      if(no_fluid_neutral_sources .eq. 1) then
       do is=0,ns-1
        if(is_neutral(is)) then
          sna(:,:,:,is) =0.d0
          smo(:,:,:,is) =0.d0
          smq(:,:,:,is) =0.d0
        endif
       enddo
      endif

* ..return
      if(ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sral

      end subroutine b2sral
