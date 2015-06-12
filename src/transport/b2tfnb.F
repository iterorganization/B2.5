*-----------------------------------------------------------------------
*.specification

      subroutine b2tfnb (nx, ny, isb, ismain,                           !srv 08.09.01 !lk 02.06.11
     &  guard_flows, facdrift, fac_ExB, mb, vol, zb,                    !xpb
     &  pbs, pbshz, nb, pb, ub, ubdia, vbecrb, vbdia, wbdia, ne,        !xpb !srv 09.01.01
     &  te, ti, bb, hx, hy, hz, qz, gs, qc, po, cvlb, cdnb, cdpb, cddb, !srv 09.01.01 13.10.06
     &  cvlbhz, cvmbhz, cdnbhz, cdpbhz,                                 !srv 09.01.01
     &  fchvispar, fchinert, fchanml,                                   !srv 20.09.06
     &  fnb, fnb_nodrift, fnb_fcor, fnb_he, fnbPSch, fnb_mdf,           !srv 13.10.06
     &  flob, conb, fnb_32, fnb_52, fnb_32nd,  
     &  fnb_52nd,kbnrgy, fllim0fnb)
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_work
     , , only : iffboundx, iffboundy
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb, guard_flows, ismain                           !xpb !srv 08.09.01 !lk 02.06.11
      real (kind=R8) ::
     *  pbs(-1:nx,-1:ny,0:1), pbshz(-1:nx,-1:ny,0:1), nb(-1:nx,-1:ny),   !srv 09.01.01
     &  pb(-1:nx,-1:ny),
     &  ub(-1:nx,-1:ny), cvlb(-1:nx,-1:ny,0:1), cdnb(-1:nx,-1:ny,0:1),
     &  cdpb(-1:nx,-1:ny,0:1), ne(-1:nx,-1:ny),
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), hz(-1:nx,-1:ny),               !srv 09.01.01
     &  po(-1:nx,-1:ny), bb(-1:nx,-1:ny,0:3),
     2  qz(-1:nx,-1:ny,0:1), gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny),
     2  facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny), zb(-1:nx,-1:ny),
     3  vol(-1:nx,-1:ny), te(-1:nx,-1:ny), ti(-1:nx,-1:ny), mb,
     &  cvlbhz(-1:nx,-1:ny,0:1), cvmbhz(-1:nx,-1:ny,0:1),                !srv 09.01.01
     &  cdnbhz(-1:nx,-1:ny,0:1), cdpbhz(-1:nx,-1:ny,0:1),                !srv 09.01.01
     &  cddb(-1:nx,-1:ny,0:1),                                           !srv 13.10.06
     &  fchvispar(-1:nx,-1:ny,0:1), fchinert(-1:nx,-1:ny,0:1),           !srv 20.09.06
     &  fchanml(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  fnb(-1:nx,-1:ny,0:1), fnb_32(-1:nx,-1:ny,0:1),
     1  fnb_52(-1:nx,-1:ny,0:1), fnb_nodrift(-1:nx,-1:ny,0:1),
     1  fnb_32nd(-1:nx,-1:ny,0:1), fnb_52nd(-1:nx,-1:ny,0:1), 
     1  flob(-1:nx,-1:ny,0:1),  fnb_fcor(-1:nx,-1:ny,0:1),           !srv 01.10.99
     2  fnb_he(-1:nx,-1:ny,0:1),                                     !srv 21.01.00
     &  fnbPSch(-1:nx,-1:ny,0:1), fnb_mdf(-1:nx,-1:ny,0:1),          !srv 13.10.06
     &  conb(-1:nx,-1:ny,0:1,0:4), ubdia(-1:nx,-1:ny,0:1),           !srv 13.10.06          
     &  vbecrb(-1:nx,-1:ny,0:1), vbdia(-1:nx,-1:ny,0:1),                 !xpb
     &  wbdia(-1:nx,-1:ny,0:1), kbnrgy(-1:nx,-1:ny),
     &  fllim0fnb(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFNB computes the particle flux of one atomic species between
*     neighbouring cells, and the coefficients that specify its
*     linearisation.
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the particle flux, fnb, of one atomic
*     species, and the coefficients flob and conb that linearise the
*     response of the particle flux to changes in the pressure.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, ib, area_fix, no_hybr,
     &        ixl, iyl, ixb, iyb                                         !srv 13.10.06
      integer, save :: iout = 0,                                         !srv 16.09.02
     &                 iout_b2npc9 = 0, iout_b2npco = 0, iout_b2npmo =0
      integer, save :: mdf_fnb = 1                                       !srv 13.10.06
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0              !sv 06.07.06
      integer, save :: PSch_style = 1                                    !srv 12.11.08
      character chns*3                                                   !srv 16.09.02
      integer, save :: fnb_nodrift_style = 1                             !srv 20.06.08
      integer, save :: drift_style = 1                                   !srv 21.06.08
      real (kind=R8), save ::  cfvma = 0.0_R8                            !srv 15.12.05
      real (kind=R8), save ::  alpha_ = 0.0_R8                           !srv 09.07.03
      real (kind=R8), save ::  ycur = 0.0_R8                             !srv 07.07.08
      real (kind=R8), save :: xvecrb = 1.0_R8, xwdia = 1.0_R8            !srv 27.01.12
      real (kind=R8), save ::  
     &   poleldr = 1.0_R8, radeldr = 1.0_R8,                             !srv 25.04.00
     &   poldidr = 1.0_R8, raddidr = 1.0_R8,                             !srv 25.04.00
     &   vis_par = 1.0_R8, fch_inert = 1.0_R8,                           !srv 20.09.06
     &   fch_anomalous = 1.0_R8, vis_per = 0.0_R8, PSch = 1.0_R8         !srv 17.01.07
      real (kind=R8), save :: xfrhie_no_drift = 0.0_R8                   !srv 20.11.09
      real (kind=R8), save :: xfrhiehz = 1.0_R8                          !srv 05.11.09 30.01.10
      real (kind=R8) ::
     *   con(-1:nx,-1:ny), dif, flo, con_nodrift, flo_nodrift, xfrhie, 
     *   cutlo, alpha, gamma, vbar, t0, t1, t2, t3, maxfluxlimit(0:1),
     *   flux_limit_min_ti,
     *   flo2dia, con_2dia(-1:nx,-1:ny), flo_mdf, con_mdf(-1:nx,-1:ny),  !srv 01.10.99 21.01.00 22.05.03 04.02.05 13.10.06
     *   term(-1:nx,-1:ny), flo_he, con_he(-1:nx,-1:ny),                 !srv 25.04.00 13.09.02 11.09.09
     &   termhz(-1:nx,-1:ny)                                             !srv 05.11.09
      real (kind=R8) :: 
     &   vbecrb_(-1:nx,-1:ny,0:1), wbdia_(-1:nx,-1:ny,0:1),              !srv 14.07.10 {
     &   wedia_(-1:nx,-1:ny,0:1), wedia(-1:nx,-1:ny,0:1),                !srv 14.07.10 }
     &   kmprti(-1:nx,-1:ny,0:1), dkmprti(-1:nx,-1:ny,0:1)               !srv 21.10.10
      parameter (cutlo=1.0e-60_R8)
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny), wrk1(-1:nx,-1:ny), wrk2(-1:nx,-1:ny),
     *   wrk3(-1:nx,-1:ny), wrk4(-1:nx,-1:ny), wrk5(-1:nx,-1:ny),        !srv 11.09.09
     *   wrk6(-1:nx,-1:ny), wrk7(-1:nx,-1:ny), wrk8(-1:nx,-1:ny),        !srv 11.09.09 !iys 17.04.12
     &   wpo(-1:nx,-1:ny), wpb(-1:nx,-1:ny)
      real (kind=R8) ::
     *   poc(-1:nx,-1:ny), ub2dia(-1:nx,-1:ny,0:1), nbcor(-1:nx,-1:ny),  !srv 01.10.99 !lk 09.06.11
     *   nbti(-1:nx,-1:ny), nbtic(-1:nx,-1:ny)                           !srv 08.10.04
      real (kind=R8) :: weight(-1:nx,-1:ny)                              !lk 
      save ncall, alpha, gamma, flux_limit_min_ti, area_fix, no_hybr
*   ..procedures
      real (kind=R8) ::
     *   kmpr, hybr, kap, cvmgp, hy1
      external subini, subend, xertst, ipgeti, ipgetr                    !xpb
      external b2xvsg, b2xvff, intfaceh, intfacev, cvmgp
      intrinsic abs, max, min, mod, sqrt
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
      kmpr(ix,iy) = nb(ix,iy)/pb(ix,iy)
      kap(ix,iy) = 1.0_R8 / (zb(ix,iy)*te(ix,iy) + ti(ix,iy))
      hybr(flo,dif) = max(dif,abs(flo/2.),sqrt(dif**2+flo**2/4.)) !!!
*   ..initialisation
      data ncall/0/, xfrhie/1.0_R8/, alpha/0.0_R8/, gamma/2.0_R8/, 
     *     flux_limit_min_ti/0.0_R8/, area_fix/3/, no_hybr/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfnb')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgetr ('b2tfnb_xfrhie', xfrhie)
        call xertst (0.0_R8.le.xfrhie,
     &   'faulty internal parameter xfrhie')
        call ipgetr ('b2tfnb_xfrhie_no_drift', xfrhie_no_drift)       !srv 20.11.09 } 
        call xertst (0.0_R8.le.xfrhie_no_drift,
     &   'faulty internal parameter xfrhie_no_drift')                 !srv 20.11.09 } 
        call ipgetr ('b2tfnb_xfrhiehz', xfrhiehz)                     !srv 05.11.09 {
        call xertst (0.0_R8.le.xfrhiehz,
     &   'faulty internal parameter xfrhiehz')                        !srv 05.11.09 }
        call ipgetr ('b2tfnb_alpha', alpha)
        call xertst (0.0_R8.le.alpha,
     &   'faulty internal parameter alpha')
        call ipgetr ('b2tfnb_gamma', gamma)
        call xertst (0.0_R8.lt.gamma,
     &   'faulty internal parameter gamma')
        call ipgetr ('b2tfnb_flux_limit_min_ti', flux_limit_min_ti)
        flux_limit_min_ti=flux_limit_min_ti*ev
        call xertst (0.0_R8.le.flux_limit_min_ti,
     &   'faulty internal parameter flux_limit_min_ti')
        call ipgeti ('b2news_area_fix', area_fix)
        call ipgeti ('b2tfnb_no_hybr', no_hybr)
        call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)                  !sv 06.07.06
        call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)                  !sv 06.07.06
        call ipgetr ('b2tfnb_poleldr', poleldr)                         !srv 25.04.00 {
        call ipgetr ('b2tfnb_radeldr', radeldr)
        call ipgetr ('b2tfnb_poldidr', poldidr)
        call ipgetr ('b2tfnb_raddidr', raddidr)                         !srv 25.04.00 }
        call ipgeti ('b2tfnb_iout', iout)                               !srv 16.09.02
        call ipgetr ('b2tlc0_alpha', alpha_)                            !srv 09.07.03 {
        call xertst (0.0_R8.le.alpha_,
     &   'faulty internal parameter alpha')
        call xertst (.not.(0.0_R8.ne.alpha .and. 0.0_R8.ne.alpha_),     !srv 14.10.09
     &   'b2tfnb_alpha and b2tlc0_alpha cannot be differed from zero ') !srv 09.07.03 }
        call ipgetr ('b2tqna_cfvma', cfvma)                             !srv 15.12.05
        call ipgeti ('b2tfnb_mdf', mdf_fnb)                             !srv 13.10.06
        call ipgetr ('b2tfnb_PSch', PSch)                               !srv 17.01.07
        call ipgeti ('b2tfnb_PSch_style', PSch_style)                   !srv 12.11.08
        call ipgeti ('b2tfnb_fnb_nodrift_style', fnb_nodrift_style)     !srv 20.06.08
        call ipgeti ('b2tfnb_drift_style', drift_style)                 !srv 21.06.08
        call ipgetr ('b2tfnb_ycur', ycur)                               !srv 17.07.08
        call ipgeti ('b2npc9_iout', iout_b2npc9)                        !srv 11.09.09
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgetr ('b2tfrn_xvecrb', xvecrb)                             !srv 27.01.12
        call ipgetr ('b2tfrn_xwdia', xwdia)                               !srv 27.01.12
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of nb, ne, pb, cdnb, cdpb
       call b2xvsg (n2, nb, 1, 'nb', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, pb, 1, 'pb', '.gt.')
       call b2xvsg (n2*2, cdnb, 1, 'cdnb', '.ge.')
       call b2xvsg (n2*2, cdpb, 1, 'cdpb', '.ge.')
*    ..test edge values of cvlb, cdnb, cdpb
       call b2xvff (nx, ny, cvlb, 'cvlb')
       call b2xvff (nx, ny, cdnb, 'cdnb')
       call b2xvff (nx, ny, cdpb, 'cdpb')
      endif

      if(mdf_fnb.eq.0) then                                     !srv 21.10.10 {
       kmprti  = 0.0_R8
       dkmprti = 0.0_R8
      else
       if((maxval(facdrift).ne.0.0_R8 .or. maxval(fac_ExB).ne.0.0_R8)
     &    .and. minval(zb).gt.0.0_R8) then
*   ..compute kmpr*ti on the left face cells and the bottom cells
        do iy = -1, ny
         do ix = -1, nx
          wrk3(ix,iy) = kmpr(ix,iy)*ti(ix,iy)
         enddo
        enddo
        call intfaceh(nx,ny,leftix,leftiy,vol,wrk3,kmprti(-1,-1,0))
        call intfacev(nx,ny,bottomix,bottomiy,vol,wrk3,kmprti(-1,-1,1))
*   ..compute d(kmpr*ti)/dy at the left face cells
        do iy = -1, ny
         do ix = -1, nx
          if(bottomiy(ix,iy).eq.-2) then
           wrk4(ix,iy) = 0.0_R8
          else
           wrk4(ix,iy)=
     &      (wrk3(bottomix(ix,iy),bottomiy(ix,iy))-wrk3(ix,iy))
!     &    /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))*2.0_R8 because cddi contains factor vol/hx/hy
          endif
         enddo
        enddo
        do iy=-1,ny
         do ix=-1,nx
          if(bottomix(ix,iy).eq.-2) then
           wrk5(ix,iy)=wrk4(topix(ix,iy),topiy(ix,iy))
          elseif(topiy(ix,iy).eq.ny+1) then
           wrk5(ix,iy)=wrk4(ix,iy)
          else
           wrk5(ix,iy)=0.5_R8*
     1         (wrk4(topix(ix,iy),topiy(ix,iy))+wrk4(ix,iy))
          endif
         enddo
        enddo
        call intfaceh(nx,ny,leftix,leftiy,vol,wrk5,dkmprti(-1,-1,0))
*   ..compute d(kmpr*ti)/dx at the bottom face cells
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2) then
           wrk4(ix,iy) = 0.0_R8
          else
           wrk4(ix,iy) = (wrk3(leftix(ix,iy),leftiy(ix,iy))-wrk3(ix,iy))
!     &    /(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*2.0_R8  because cddi contains factor vol/hx/hy
          endif
         enddo
        enddo
        do iy=-1,ny
         do ix=-1,nx
          if(leftix(ix,iy).eq.-2) then
           wrk5(ix,iy)=wrk4(rightix(ix,iy),rightiy(ix,iy))
          elseif(rightix(ix,iy).eq.nx+1) then
           wrk5(ix,iy)=wrk4(ix,iy)
          else
           wrk5(ix,iy)=0.5_R8*
     1         (wrk4(rightix(ix,iy),rightiy(ix,iy))+wrk4(ix,iy))
          endif
         enddo
        enddo
        call intfacev(nx,ny,bottomix,bottomiy,vol,wrk5,dkmprti(-1,-1,1))
       else
        kmprti  = 0.0_R8
        dkmprti = 0.0_R8
       endif
      endif                                                        !srv 21.10.10 }

c Drift terms
      maxfluxlimit=0.0_R8
      ubdia =0.0_R8
      ub2dia=0.0_R8                                                 !srv 01.10.99
      vbdia =0.0_R8
      wbdia =0.0_R8
      wedia =0.0_R8
      wedia_=0.0_R8
      vbecrb=0.0_R8
      fnbPSch   = 0.0_R8                                              !srv 13.10.06
      fnb_mdf = 0.0_R8                                                !srv 19.11.09
      conb  = 0.0e0_R8                                                !srv 13.10.06
      if((maxval(facdrift).ne.0.0_R8 .or. maxval(fac_ExB).ne.0.0_R8)
     &    .and. minval(zb).gt.0.0_R8) then
       if(drift_style.eq.0) then                                      !srv 21.06.08
        call intfacev(nx,ny,bottomix,bottomiy,vol,po,wpo)
        call intfacev(nx,ny,bottomix,bottomiy,vol,nb*ti,wpb)
        do iy = -1, ny
          do ix = -1, nx                                                 !xpb
cxpb  Must check that we are not at the island X-point or at the mesh edge !xpb
            if(topiy(topix(ix,iy),min(topiy(ix,iy),ny)).ne.iy .and.
     &         topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &         rightix(ix,iy).ne.nx+1 .and. leftiy(ix,iy).ne.-2 .and.
     &         region(ix,iy,0).ne.0) then                                !xpb
cxpb  Standard case
              vbecrb(ix,iy,0) = - signmf/hy1(ix,iy)*                     !xpb
     1         1.0_R8/bb(ix,iy,3)*                                       !xpb
     2         (wpo(topix(ix,iy),topiy(ix,iy))-wpo(ix,iy))*
     &         fac_ExB(ix,iy)                                            !xpb
              vbdia(ix,iy,0) = - bb(ix,iy,2)/bb(ix,iy,3)/hy1(ix,iy)*     !xpb
     1         1.0_R8/bb(ix,iy,3)*                                       !xpb
     2         2.0_R8*ti(ix,iy)/bb(ix,iy,3)/qe/zb(ix,iy)*                !xpb
     3        (wbbv(topix(ix,iy),topiy(ix,iy),3)-wbbv(ix,iy,3))*
     &        facdrift(ix,iy)                                            !xpb
              wbdia(ix,iy,0) = - bb(ix,iy,2)/bb(ix,iy,3)/hy1(ix,iy)
     1         /bb(ix,iy,3)/qe/zb(ix,iy)/nb(ix,iy)*
     2         (wpb(topix(ix,iy),topiy(ix,iy))-wpb(ix,iy))*
     &         facdrift(ix,iy)
            else                                                         !xpb
cxpb  Case of the island X-point cells
cxpb          vbecrb(ix,iy,0) = - signmf/hy1(ix,iy)*fac_ExB(ix,iy)*      !xpb
cxpb     1             1.0_R8/bb(ix,iy,3)*                               !xpb
cxpb     2             ((po(topix(ix,iy),topiy(ix,iy))*vol(ix,iy)+       !xpb
cxpb     2               po(ix,iy)*vol(topix(ix,iy),topiy(ix,iy)))/      !xpb
cxpb     2          (vol(ix,iy)+vol(topix(ix,iy),topiy(ix,iy)))-wpo(ix,iy)) !xpb
cxpb          vbdia(ix,iy,0) = - signmf/hy1(ix,iy)*facdrift(ix,iy)*      !xpb
cxpb     1             1.0_R8/bb(ix,iy,3)*                               !xpb
cxpb     2             2.0_R8*ti(ix,iy)/bb(ix,iy,3)/qe/zb(ix,iy)*        !xpb
cxpb     3         ((vol(ix,iy)*bb(topix(ix,iy),topiy(ix,iy),3)+         !xpb
cxpb     3           vol(topix(ix,iy),topiy(ix,iy))*bb(ix,iy,3))/        !xpb
cxpb     3              (vol(topix(ix,iy),topiy(ix,iy))+vol(ix,iy))-     !xpb
cxpb     3              wbbv(ix,iy,3))                                   !xpb
              vbecrb(ix,iy,0) = 0.0_R8                                   !xpb
              vbdia(ix,iy,0) = 0.0_R8                                    !xpb
            endif                                                        !xpb
            ubdia(ix,iy,0) = vbecrb(ix,iy,0) + vbdia(ix,iy,0)    !srv 17.05.99
            ub2dia(ix,iy,0) = vbecrb(ix,iy,0) + 2.0_R8*vbdia(ix,iy,0)   !srv 01.10.99
          enddo
        enddo
        
        call intfaceh(nx,ny,leftix,leftiy,vol,po,wpo)
        call intfaceh(nx,ny,leftix,leftiy,vol,nb*ti,wpb)
        do iy = -1, ny                                                   !xpb
          do ix = -1, nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
             vbecrb(ix,iy,1) = 1.0_R8/hx(ix,iy)*                        !xpb
     1          bb(ix,iy,2)/bb(ix,iy,3)**2*                              !xpb
     2         (wpo(rightix(ix,iy),rightiy(ix,iy))-wpo(ix,iy))*
     &         fac_ExB(ix,iy)                                            !xpb
             vbdia(ix,iy,1) = ti(ix,iy)*bb(ix,iy,2)
     &         /(qe*zb(ix,iy)*hx(ix,iy))                                 !xpb
     &         *(2.0_R8/bb(ix,iy,3)**3)                                  !xpb
     &         *(wbbr(ix,iy,3)-wbbl(ix,iy,3))*facdrift(ix,iy)            !xpb
             wbdia(ix,iy,1) = bb(ix,iy,2)
     &         /(qe*zb(ix,iy)*hx(ix,iy)*nb(ix,iy))                       !xpb
     &         /bb(ix,iy,3)**2*                                          !xpb
     2         (wpb(rightix(ix,iy),rightiy(ix,iy))-wpb(ix,iy))*
     &         facdrift(ix,iy)
             ubdia(ix,iy,1) = vbecrb(ix,iy,1) + vbdia(ix,iy,1)         !srv 17.05.99
             ub2dia(ix,iy,1) = vbecrb(ix,iy,1) + 2.0_R8*vbdia(ix,iy,1) !srv 01.10.99
            endif
          enddo
        enddo

* drift velocities at walls are extrapolated from those at neighbouring cells
* drift velocities at iy=ny are extrapolated from those at iy=ny-1       !xpb
* drift velocities at iy=-1 are extrapolated from those at iy=0          !xpb
        do ix = -1, nx                                                   !xpb
          do iy = -1, ny                                                 !xpb
            if(leftix(ix,iy).eq.-2) then
              vbecrb(ix,iy,:) = vbecrb(rightix(ix,iy),rightiy(ix,iy),:)  !xpb
              vbdia(ix,iy,:) = vbdia(rightix(ix,iy),rightiy(ix,iy),:)    !xpb
              wbdia(ix,iy,:) = wbdia(rightix(ix,iy),rightiy(ix,iy),:)    !xpb
              ubdia(ix,iy,:) = vbecrb(ix,iy,:) + vbdia(ix,iy,:) 
              ub2dia(ix,iy,:) = vbecrb(ix,iy,:) + 2.0_R8*vbdia(ix,iy,:)  !srv 17.07.08
            elseif(rightix(ix,iy).eq.nx+1) then                          !xpb
              vbecrb(ix,iy,:) = vbecrb(leftix(ix,iy),leftiy(ix,iy),:)    !xpb
              vbdia(ix,iy,:) = vbdia(leftix(ix,iy),leftiy(ix,iy),:)      !xpb
              wbdia(ix,iy,:) = wbdia(leftix(ix,iy),leftiy(ix,iy),:)      !xpb
              ubdia(ix,iy,:) = vbecrb(ix,iy,:) + vbdia(ix,iy,:) 
              ub2dia(ix,iy,:) = vbecrb(ix,iy,:) + 2.0_R8*vbdia(ix,iy,:)  !srv 17.07.08 
            endif
            if(topiy(ix,iy).eq.ny+1) then                                !xpb
              vbecrb(ix,iy,:) = 
     .         vbecrb(bottomix(ix,iy),bottomiy(ix,iy),:)                 !xpb
              vbdia(ix,iy,:) = vbdia(bottomix(ix,iy),bottomiy(ix,iy),:)  !xpb
              wbdia(ix,iy,:) = wbdia(bottomix(ix,iy),bottomiy(ix,iy),:)  !xpb
              ubdia(ix,iy,:) = vbecrb(ix,iy,:) + vbdia(ix,iy,:) 
              ub2dia(ix,iy,:) = vbecrb(ix,iy,:) + 2.0_R8*vbdia(ix,iy,:)  !srv 17.07.08 
            elseif(bottomiy(ix,iy).eq.-2) then
              vbecrb(ix,iy,:) = vbecrb(topix(ix,iy),topiy(ix,iy),:)      !xpb
              vbdia(ix,iy,:) = vbdia(topix(ix,iy),topiy(ix,iy),:)        !xpb
              wbdia(ix,iy,:) = wbdia(topix(ix,iy),topiy(ix,iy),:)        !xpb
              ubdia(ix,iy,:) = vbecrb(ix,iy,:) + vbdia(ix,iy,:) 
              ub2dia(ix,iy,:) = vbecrb(ix,iy,:) + 2.0_R8*vbdia(ix,iy,:)  !srv 17.07.08 
            endif
          enddo
        enddo

* turn off side drift flows in guard cells                               !xpb
        do iy = -1, ny                                                   !xpb
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
           vbecrb(ix,iy,1) = 0.0_R8                                      !xpb
           vbdia(ix,iy,1) = 0.0_R8                                       !xpb
           wbdia(ix,iy,1) = 0.0_R8                                       !xpb
           ubdia(ix,iy,1) = 0.0_R8                                       !xpb
           ub2dia(ix,iy,1) = 0.0_R8                                      !srv 17.07.08
          endif
          if(bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1) then
           vbecrb(ix,iy,0) = 0.0_R8                                      !xpb
           vbdia(ix,iy,0) = 0.0_R8                                       !xpb
           wbdia(ix,iy,0) = 0.0_R8                                       !xpb
           ubdia(ix,iy,0) = 0.0_R8                                       !xpb
           ub2dia(ix,iy,0) = 0.0_R8                                      !srv 17.07.08
          endif
         enddo
        enddo                                                            !xpb 
       else                                                              !srv 21.06.08 {
        weight = 1.0e0_R8                                                !lk   
!lk       call intcor(nx,ny,vol,po,poc)
        call intcor(nx,ny,weight,po,poc)                                 !lk 12.12.07
        nbti=nb*ti
!lk       call intcor(nx,ny,vol,nbti,nbtic)                                !srv 08.10.04
        call intcor(nx,ny,weight,nbti,nbtic)                             !lk 12.12.07
        do iy = 0, ny-1
         do ix = 0, nx
          if(leftix(ix,iy).ne.-2) then                                   !srv 26.11.02
!            vbecrb(ix,iy,0) = - poleldr*fac_ExB(ix,iy)*signmf/            !srv 17.05.99 25.08.99 07.06.00 05.02.05
!     1      (hy(ix,iy)+hy(leftix(ix,iy),leftiy(ix,iy)))*             !srv
!     2      (bb(ix,iy,2)/bb(ix,iy,3)**2                              !srv 05.02.05
!     2         +bb(leftix(ix,iy),leftiy(ix,iy),2)/
!     2          bb(leftix(ix,iy),leftiy(ix,iy),3)**2)*                  !srv 05.02.05
!     2       (poc(topix(ix,iy),topiy(ix,iy))-poc(ix,iy))             !srv 
           vbecrb(ix,iy,0) = - poleldr*fac_ExB(ix,iy)* ! signmf*           !lk 12.12.07 } !srv 11.01.13
     1       (hz(ix,iy)*bb(ix,iy,2)
     2         +hz(leftix(ix,iy),leftiy(ix,iy))*
     3         bb(leftix(ix,iy),leftiy(ix,iy),2))/2.0e0_R8*              
     4       (vol(leftix(ix,iy),leftiy(ix,iy))/bb(ix,iy,3)**2      
     5         +vol(ix,iy)/bb(leftix(ix,iy),leftiy(ix,iy),3)**2)/    
     6         (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))*
     7         (poc(topix(ix,iy),topiy(ix,iy))-poc(ix,iy))/
     8         (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))*
     9         (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))               !lk 12.12.07)
           vbdia(ix,iy,0) =    poldidr*facdrift(ix,iy)* ! signmf*      !srv 07.06.00 17.05.99 25.04.00 05.02.05 11.01.13
     &          (ti(leftix(ix,iy),leftiy(ix,iy))                       !srv 05.03.99
     &           *bb(leftix(ix,iy),leftiy(ix,iy),2)                    !srv 05.02.05
     &           /(qe*zb(ix,iy)*hy1(leftix(ix,iy),leftiy(ix,iy)))              !srv
     &          +ti(ix,iy)*bb(ix,iy,2)/(qe*zb(ix,iy)*hy1(ix,iy)))*0.5_R8    !srv 05.02.05
     &         *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))  !srv 05.03.99
           wbdia(ix,iy,0) = - bb(ix,iy,2)/bb(ix,iy,3)/hy1(ix,iy)         !srv 08.10.04 {
     1         /bb(ix,iy,3)/qe/zb(ix,iy)/nb(ix,iy)*
     2         (nbtic(topix(ix,iy),topiy(ix,iy))-nbtic(ix,iy))*         !srv 08.10.04 }
     &          facdrift(ix,iy)
           ubdia(ix,iy,0) = vbecrb(ix,iy,0) + vbdia(ix,iy,0)            !srv 17.05.99 
           ub2dia(ix,iy,0) = vbecrb(ix,iy,0) + 2.0e0_R8*vbdia(ix,iy,0)   !srv 01.10.99
          endif                                                         !srv 26.11.02
         enddo
        enddo

        wrk0=0.0_R8  ! compute .... on the cell faces                   !srv 14.07.10 {
        do iy=0,ny
         do ix=-1,nx
          wrk0(ix,iy)=(ne(ix,iy)*te(ix,iy)-
     &                 ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &                 te(bottomix(ix,iy),bottomiy(ix,iy)))/
     &          (hy(ix,iy)+hy(bottomix(ix,iy),bottomiy(ix,iy)))*2.0_R8
         enddo
        enddo
        do iy=-1,ny  ! compute ... on the cell centers
         do ix=-1,nx
          if(bottomiy(ix,iy).eq.-2) then
           wedia_(ix,iy,0) = wrk0(ix,topiy(ix,iy))/
     &           (qe*ne(ix,iy)*bb(ix,iy,2)*bb(ix,iy,3)**2)
          elseif(topiy(ix,iy).eq.ny+1) then
           wedia_(ix,iy,0) = wrk0(ix,iy)/
     &           (qe*ne(ix,iy)*bb(ix,iy,2)*bb(ix,iy,3)**2)
          else
           wedia_(ix,iy,0) = 0.5_R8*
     &        (wrk0(topix(ix,iy),topiy(ix,iy))+wrk0(ix,iy))/
     &           (qe*ne(ix,iy)*bb(ix,iy,2)*bb(ix,iy,3)**2)
          endif
         enddo
        enddo
        call intfaceh(nx,ny,leftix,leftiy,vol,wedia_(-1,-1,0),
     &                                        wedia(-1,-1,0))

        wrk0=0.0_R8  ! compute ... on the cell faces
        do iy=-1,ny
         do ix=-1,nx
          if(leftix(ix,iy).ne.-2) then
           wrk0(ix,iy)=(ne(ix,iy)*te(ix,iy)-
     &         ne(leftix(ix,iy),leftiy(ix,iy))*
     &         te(leftix(ix,iy),leftiy(ix,iy)))/
     &        (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0_R8
          endif
         end do
        end do
        do iy=-1,ny  ! compute .... on the cell centers
         do ix=-1,nx
          if(leftix(ix,iy).eq.-2) then
           wedia_(ix,iy,1)=-wrk0(rightix(ix,iy),iy)/
     &         (qe*ne(ix,iy)*bb(ix,iy,2)*bb(ix,iy,3)**2)
          elseif(rightix(ix,iy).eq.nx+1) then
           wedia_(ix,iy,1)=-wrk0(ix,iy)/
     &         (qe*ne(ix,iy)*bb(ix,iy,2)*bb(ix,iy,3)**2)
          else
           wedia_(ix,iy,1)=-0.5_R8*
     &        (wrk0(rightix(ix,iy),rightiy(ix,iy))+
     &         wrk0(ix,iy))/
     &        (qe*ne(ix,iy)*bb(ix,iy,2)*bb(ix,iy,3)**2)
          endif
         enddo
        enddo
        call intfacev(nx,ny,bottomix,bottomiy,vol,wedia_(-1,-1,1),
     &                                             wedia(-1,-1,1))       !srv 14.07.10 }


        do iy = 0, ny
         do ix = 0, nx-1
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then       !srv 26.11.02
           vbecrb(ix,iy,1) = radeldr*fac_ExB(ix,iy)* ! signmf*            !srv 11.01.13 07.06.00 21.04.00 25.08.99 17.05.99
     1       (hz(ix,iy)*bb(ix,iy,2)                                    !lk ( 12.12.07
     2          +hz(bottomix(ix,iy),bottomiy(ix,iy))*
     3          bb(bottomix(ix,iy),bottomiy(ix,iy),2))/2.0e0_R8*              
     4       (vol(bottomix(ix,iy),bottomiy(ix,iy))/bb(ix,iy,3)**2      
     5          +vol(ix,iy)/bb(bottomix(ix,iy),bottomiy(ix,iy),3)**2)/    
     6          (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))*
     7        (poc(rightix(ix,iy),rightiy(ix,iy))-poc(ix,iy))/          !lk 12.12.07)   
     8          (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))*
     9          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
           vbdia(ix,iy,1) =  -raddidr*facdrift(ix,iy)* ! signmf*          !srv 11.01.13 07.06.00 17.05.99 25.04.00
     &           (ti(bottomix(ix,iy),bottomiy(ix,iy))                     !srv 05.03.99
     &           *bb(bottomix(ix,iy),bottomiy(ix,iy),2)                  !srv
     &           /(qe*zb(ix,iy)*hx(bottomix(ix,iy),bottomiy(ix,iy)))            !srv
     &           +ti(ix,iy)*bb(ix,iy,2)/(qe*zb(ix,iy)*hx(ix,iy)))*0.5_R8         !srv
     &        *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy)) !srv 05.03.99
           wbdia(ix,iy,1) = bb(ix,iy,2)/                                 !srv 08.10.04 {
     &          (qe*zb(ix,iy)*hx(ix,iy)*nb(ix,iy))
     &          /bb(ix,iy,3)**2*
     2          (nbtic(rightix(ix,iy),rightiy(ix,iy))-nbtic(ix,iy))*       !srv 08.10.04 }
     &          facdrift(ix,iy)
           ubdia(ix,iy,1) = vbecrb(ix,iy,1) + vbdia(ix,iy,1)              !srv 17.05.99
           ub2dia(ix,iy,1) = vbecrb(ix,iy,1) + 2.0e0_R8*vbdia(ix,iy,1)    !srv 01.10.99
          endif                                                           !srv 26.11.02
         enddo
        enddo

       endif                                                              !srv 21.06.08 {

*   ..compute P.Sch x-flux                                 !srv 13.10.06 {
       if(PSch_style.eq.0) then                            !srv 12.11.08
        call intcor(nx,ny,vol,po,poc)
        nbti=nb*ti
        call intcor(nx,ny,vol,nbti,nbtic)                                !srv 08.10.04
        do iy = 0, ny-1
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. 
     &     bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 ) then     !srv 26.05.08
           fnbPSch(ix,iy,0)=0.0e0_R8
           conb(ix,iy,0,3) = 0.0e0_R8                                  !srv 26.05.08
           conb(ix,iy,0,4) = 0.0e0_R8
          else
           ixl = leftix(ix,iy)
           iyl = leftiy(ix,iy)
           fnbPSch(ix,iy,0) = 0.25e0_R8*PSch*( 
     &      cddb(topix(ix,iy),topiy(ix,iy),1)*
     &       (nb(topix(ix,iy),topiy(ix,iy))*
     &        ti(topix(ix,iy),topiy(ix,iy))- nb(ix,iy)*ti(ix,iy))+
     &      cddb(topix(ixl,iyl),topiy(ixl,iyl),1)*
     &       (nb(topix(ixl,iyl),topiy(ixl,iyl))*
     &        ti(topix(ixl,iyl),topiy(ixl,iyl))-
     &        nb(ixl,iyl)*ti(ixl,iyl))+
     &      cddb(ixl,iyl,1)*
     &       (nb(ixl,iyl)*ti(ixl,iyl)
     &                   -nb(bottomix(ixl,iyl),bottomiy(ixl,iyl))*
     &                    ti(bottomix(ixl,iyl),bottomiy(ixl,iyl)))+
     &      cddb(ix,iy,1)*
     &       (nb(ix,iy)*ti(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy))*
     &                         ti(bottomix(ix,iy),bottomiy(ix,iy)))
     &                                       )
!           conb(ix,iy,0,3) = 0.25e0_R8*PSch*cddb(ix,iy,1)*ti(ix,iy)  !srv 21.10.10 {
!           conb(ix,iy,0,4) = 0.25e0_R8*PSch*cddb(ix,iy,1)*
!     &                      ti(bottomix(ix,iy),bottomiy(ix,iy))
           conb(ix,iy,0,3) = 0.25_R8*PSch*cddb(ix,iy,1)*kmprti(ix,iy,1)
           conb(ix,iy,0,4) = conb(ix,iy,0,3)                          !srv 21.10.10 }
          endif
         enddo
        enddo
       else
        weight = 1.0e0_R8                                                !lk   
        call intcor(nx,ny,weight,po,poc)                                 !lk 12.12.07
        nbti=nb*ti
        call intcor(nx,ny,weight,nbti,nbtic)                             !lk 12.12.07
        call intcor(nx,ny,weight,nb,nbcor)                               !lk 09.06.11  
        do iy = 0, ny-1
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. 
     &     bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 ) then     !srv 26.05.08
           fnbPSch(ix,iy,0)=0.0e0_R8
           conb(ix,iy,0,3) = 0.0e0_R8                                  !srv 26.05.08
           conb(ix,iy,0,4) = 0.0e0_R8
          else
           ixl = leftix(ix,iy)
           iyl = leftiy(ix,iy)
!           fnbPSch(ix,iy,0) = PSch*                                           !lk {
!     &      (nbtic(topix(ix,iy),topiy(ix,iy))-nbtic(ix,iy))*
!     &      (bzb(ix,iy)*hz(ix,iy)*vol(ixl,iyl)
!     &      +bzb(ixl,iyl)*hz(ixl,iyl)*vol(ix,iy))
!     &      /(vol(ixl,iyl)+vol(ix,iy))/(qe*zb(ix,iy))                          !lk }
           fnbPSch(ix,iy,0) = PSch*                                           !lk {
     &      (nbtic(topix(ix,iy),topiy(ix,iy))-nbtic(ix,iy))*
     &      (bzb(ix,iy)*(vol(ix,iy)/gs(ix,iy,2))*vol(ixl,iyl)
     &      +bzb(ixl,iyl)*(vol(ixl,iyl)/gs(ixl,iyl,2))*vol(ix,iy))
     &      /(vol(ixl,iyl)+vol(ix,iy))/(qe*zb(ix,iy))                          !lk }
     &      *(nb(ix,iy)+nb(ixl,iyl))/(nbcor(ix,iy)+                       !lk 09.06.11
     &      nbcor(topix(ix,iy),topiy(ix,iy)))                             !lk 09.06.11
!            conb(ix,iy,0,3) = 0.25e0_R8*PSch*cddb(ix,iy,1)*ti(ix,iy)   !srv 21.10.10 {
!            conb(ix,iy,0,4) = 0.25e0_R8*PSch*cddb(ix,iy,1)*
!     &                      ti(bottomix(ix,iy),bottomiy(ix,iy))
           conb(ix,iy,0,3) = 0.25_R8*PSch*cddb(ix,iy,1)*kmprti(ix,iy,1)
           conb(ix,iy,0,4) = conb(ix,iy,0,3)                            !srv 21.10.10 }
          endif
         enddo
        enddo
       endif

*   ..compute P.Sch y-flux
       iy = -1
       fnbPSch(:,iy,1)=0.0e0_R8
       if(PSch_style.eq.0) then                            !srv 12.11.08
        do iy = 0, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
           fnbPSch(ix,iy,1)= 0.0e0_R8
           conb(ix,iy,1,3) = 0.0e0_R8                           !srv 26.05.08
           conb(ix,iy,1,4) = 0.0e0_R8
          else
           ixb = bottomix(ix,iy)
           iyb = bottomiy(ix,iy)
           fnbPSch(ix,iy,1) = -0.25e0_R8*PSch*(
     &    cddb(rightix(ix,iy),rightiy(ix,iy),0)*
     &     (nb(rightix(ix,iy),rightiy(ix,iy))*
     &      ti(rightix(ix,iy),rightiy(ix,iy))-nb(ix,iy)*ti(ix,iy))+
     &    cddb(ix,iy,0)*
     &    (nb(ix,iy)*ti(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy))*
     &                         ti(leftix(ix,iy),leftiy(ix,iy)))+
     &    cddb(bottomix(ix,iy),bottomiy(ix,iy),0)*
     &     (nb(bottomix(ix,iy),bottomiy(ix,iy))*
     &      ti(bottomix(ix,iy),bottomiy(ix,iy))
     &                         -nb(leftix(ixb,iyb),leftiy(ixb,iyb))*
     &                          ti(leftix(ixb,iyb),leftiy(ixb,iyb)))+
     &    cddb(rightix(ixb,iyb),rightiy(ixb,iyb),0)*
     &     (nb(rightix(ixb,iyb),rightiy(ixb,iyb))*
     &      ti(rightix(ixb,iyb),rightiy(ixb,iyb))
     &                  -nb(bottomix(ix,iy),bottomiy(ix,iy))*
     &                   ti(bottomix(ix,iy),bottomiy(ix,iy)))
     &                                         )

!           conb(ix,iy,1,3) = 0.25e0_R8*PSch*cddb(ix,iy,0)*ti(ix,iy)  !srv 21.10.10 {
!           conb(ix,iy,1,4) = 0.25e0_R8*PSch*cddb(ix,iy,0)*
!     &                        ti(leftix(ix,iy),leftiy(ix,iy))
           conb(ix,iy,1,3) = 0.25_R8*PSch*cddb(ix,iy,0)*kmprti(ix,iy,0)
           conb(ix,iy,1,4) = conb(ix,iy,1,3)                          !srv 21.10.10 }
          endif
         enddo
        enddo
       else
        do iy = 0, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
           fnbPSch(ix,iy,1)= 0.0e0_R8
           conb(ix,iy,1,3) = 0.0e0_R8                           !srv 26.05.08
           conb(ix,iy,1,4) = 0.0e0_R8
          else
           ixb = bottomix(ix,iy)
           iyb = bottomiy(ix,iy)
!           fnbPSch(ix,iy,1) = -PSch*                                          !lk {
!     &      (nbtic(rightix(ix,iy),rightiy(ix,iy))-nbtic(ix,iy))*
!     &      (bzb(ix,iy)*hz(ix,iy)*vol(ixb,iyb)
!     &      +bzb(ixb,iyb)*hz(ixb,iyb)*vol(ix,iy))
!     &      /(vol(ixb,iyb)+vol(ix,iy))/(qe*zb(ix,iy))                          !lk }
           fnbPSch(ix,iy,1) = -PSch*                                          !lk {
     &      (nbtic(rightix(ix,iy),rightiy(ix,iy))-nbtic(ix,iy))*
     &      (bzb(ix,iy)*(vol(ix,iy)/gs(ix,iy,2))*vol(ixb,iyb)
     &      +bzb(ixb,iyb)*(vol(ixb,iyb)/gs(ixb,iyb,2))*vol(ix,iy))
     &      /(vol(ixb,iyb)+vol(ix,iy))/(qe*zb(ix,iy))                          !lk }
     &      *(nb(ix,iy)+nb(ixb,iyb))/(nbcor(ix,iy)+                       !lk 09.06.11
     &      nbcor(rightix(ix,iy),rightiy(ix,iy)))                         !lk 09.06.11

!           conb(ix,iy,1,3) = 0.25e0_R8*PSch*cddb(ix,iy,0)*ti(ix,iy)    !srv 21.10.10 {
!           conb(ix,iy,1,4) = 0.25e0_R8*PSch*cddb(ix,iy,0)*
!     &                        ti(leftix(ix,iy),leftiy(ix,iy))
           conb(ix,iy,1,3) = 0.25_R8*PSch*cddb(ix,iy,0)*kmprti(ix,iy,0)
           conb(ix,iy,1,4) = conb(ix,iy,1,3)                            !srv 21.10.10 }
          endif
         enddo
        enddo
       endif

       do iy = -1, ny                                                  !srv 02.01.07 {
        do ix = -1, nx
         fnbPSch(ix,iy,0:1) = fnbPSch(ix,iy,0:1)*facdrift(ix,iy)
         conb(ix,iy,0:1,3) = conb(ix,iy,0:1,3)*facdrift(ix,iy)
         conb(ix,iy,0:1,4) = conb(ix,iy,0:1,4)*facdrift(ix,iy)
        enddo
       enddo                                                           !srv 02.01.07 }

      endif

      if (isb.ne.ismain .and. .not.is_neutral(isb)) then             !lk 02.06.11 {
       do iy = 0, ny                                                  
        do ix = -1, nx
          if (leftix(ix,iy).ne.-2) then
           t0 = max(ti(ix,iy), ti(leftix(ix,iy),leftiy(ix,iy)))
           t1 = (nb(ix,iy) + nb(leftix(ix,iy),leftiy(ix,iy)))
     &          /2.0_R8
           vbar = abs(bb(ix,iy,0) + bb(leftix(ix,iy),leftiy(ix,iy),0))   !srv 11.01.13
     &          /(bb(ix,iy,3) + bb(leftix(ix,iy),leftiy(ix,iy),3))
     &          *sqrt(8.0_R8 * t0 / (pi * am(isb)*mp))
           t2 = 1.0_R8 + abs(fnbPSch(ix,iy,0))/(
     &         vbar*t1*(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     &         (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))
     &          )
           fnbPSch(ix,iy,0)= fnbPSch(ix,iy,0)/t2
             vbdia(ix,iy,0)= vbdia(ix,iy,0)/
     &                              (1.0_R8 + abs(vbdia(ix,iy,0))/vbar) 
             wbdia(ix,iy,0)= wbdia(ix,iy,0)/
     &                              (1.0_R8 + abs(wbdia(ix,iy,0))/vbar) 
             ubdia(ix,iy,0)= ubdia(ix,iy,0)/
     &                              (1.0_R8 + abs(ubdia(ix,iy,0))/vbar) 
             ub2dia(ix,iy,0)= ub2dia(ix,iy,0)/
     &                              (1.0_R8 + abs(ub2dia(ix,iy,0))/vbar) 
          endif
          t0 = max(ti(ix,iy), ti(bottomix(ix,iy),bottomiy(ix,iy)))
          t1 = (nb(ix,iy) + nb(bottomix(ix,iy),bottomiy(ix,iy)))
     &          /2.0_R8
           vbar=abs(bb(ix,iy,0) + bb(bottomix(ix,iy),bottomiy(ix,iy),0)) !srv 11.01.13
     &          /(bb(ix,iy,3) + bb(bottomix(ix,iy),bottomiy(ix,iy),3))
     &          *sqrt(8.0_R8 * t0 / (pi * am(isb)*mp))
          t2 = 1.0_R8 + abs(fnbPSch(ix,iy,1))/(
     &        vbar*t1*(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     &         (hy(ix,iy)+hy(bottomix(ix,iy),bottomiy(ix,iy)))
     &          )
          fnbPSch(ix,iy,1)= fnbPSch(ix,iy,1)/t2
            vbdia(ix,iy,1)= vbdia(ix,iy,1)/
     &                               (1.0_R8 + abs(vbdia(ix,iy,1))/vbar) 
            wbdia(ix,iy,1)= wbdia(ix,iy,1)/
     &                               (1.0_R8 + abs(wbdia(ix,iy,1))/vbar) 
            ubdia(ix,iy,1)= ubdia(ix,iy,1)/
     &                               (1.0_R8 + abs(ubdia(ix,iy,1))/vbar) 
            ub2dia(ix,iy,1)= ub2dia(ix,iy,1)/
     &                              (1.0_R8 + abs(ub2dia(ix,iy,1))/vbar) 
           enddo
       enddo
      endif                                                            !lk 02.06.11 }

* ..compute fluxes
*   ..compute wrk0
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
          wrk0(ix,iy) = 0.0_R8
        else
          wrk0(ix,iy) = !pbs(ix,iy,0)*                                !srv 20.06.08
     &     (vol(ix,iy)*ub(leftix(ix,iy),leftiy(ix,iy))+
     &      vol(leftix(ix,iy),leftiy(ix,iy))*ub(ix,iy))/
     &     (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
        endif
       enddo
      enddo
      call b2xxgs (nx, ny, 0.0_R8, wrk0, 1)
*   ..compute wrk1
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         wrk1(ix,iy) = 0.0_R8
        else
         wrk1(ix,iy) = (pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))/
     &    (2.0_R8*(pb(leftix(ix,iy),leftiy(ix,iy))+pb(ix,iy)))
        endif
       enddo
      enddo
      call b2xxgs (nx, ny, 0.0_R8, wrk1, 1)
*   ..compute wrk2
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         wrk2(ix,iy) = 0.0_R8
        elseif(leftix(leftix(ix,iy),leftiy(ix,iy)).eq.-2) then
         wrk2(ix,iy) = 0.0_R8
        elseif(rightix(ix,iy).eq.nx+1) then
         wrk2(ix,iy) = 0.0_R8
        else
         wrk2(ix,iy) = wrk1(leftix(ix,iy),leftiy(ix,iy))-
     &    2.0_R8*wrk1(ix,iy)+wrk1(rightix(ix,iy),rightiy(ix,iy))
        endif
       enddo
      enddo
      call b2xxgs (nx, ny, 0.0_R8, wrk2, 1)

*   ..compute fnbx, flobx, conbx
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         fnb(ix,iy,0)         = 0.0_R8
         fnb_32(ix,iy,0)      = 0.0_R8
         fnb_52(ix,iy,0)      = 0.0_R8
         fnb_nodrift(ix,iy,0) = 0.0_R8
         fnb_32nd(ix,iy,0)    = 0.0_R8
         fnb_52nd(ix,iy,0)    = 0.0_R8
         flob(ix,iy,0)        = 0.0_R8
         fnb_fcor(ix,iy,0)    = 0.0_R8                              !srv 01.10.99
         fnb_he(ix,iy,0)      = 0.0_R8                              !srv 21.01.00
         fnb_mdf(ix,iy,0)     = 0.0_R8                              !srv 13.10.06 
         conb(ix,iy,0,0:2)    = 0.0_R8                              !srv 13.10.06
         term(ix,iy)          = 0.0_R8
         termhz(ix,iy)        = 0.0_R8                              !srv 05.11.09
         con(ix,iy)           = 0.0_R8
         con_2dia(ix,iy)      = 0.0_R8
         con_he(ix,iy)        = 0.0_R8
         con_mdf(ix,iy)       = 0.0_R8
        else
         term(ix,iy) = xfrhie*wrk2(ix,iy)*0.5_R8*abs(pbs(ix,iy,0))*    !iyv 19.04.07 11.01.13
     &    sqrt((pb(leftix(ix,iy),leftiy(ix,iy))+pb(ix,iy))/
     &    (mb*(nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))))*
     &    min(vol(leftix(ix,iy),leftiy(ix,iy)),vol(ix,iy))/
     &    (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
         termhz(ix,iy) = xfrhiehz*wrk2(ix,iy)*0.5_R8*abs(pbshz(ix,iy,0))  !srv 05.11.09 11.01.13
     &    *sqrt((pb(leftix(ix,iy),leftiy(ix,iy))+pb(ix,iy))/
     &    (mb*(nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))))*
     &    min(vol(leftix(ix,iy),leftiy(ix,iy)),vol(ix,iy))/
     &    (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
         if(drift_style.eq.0)then                                     !srv 21.06.08
          flo = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+term(ix,iy)+   !srv 20.06.08
     1     (bb(ix,iy,2)/bb(ix,iy,3)+                                  !srv 11.01.13
     2     bb(leftix(ix,iy),leftiy(ix,iy),2)/
     3     bb(leftix(ix,iy),leftiy(ix,iy),3))/2.0_R8
     4     *(ubdia(ix,iy,0)+ubdia(leftix(ix,iy),leftiy(ix,iy),0))/
     5     2.0_R8*0.5_R8*(vol(ix,iy)/hx(ix,iy)
     6     +vol(leftix(ix,iy),leftiy(ix,iy))/
     7     hx(leftix(ix,iy),leftiy(ix,iy)))
          flo_mdf = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+           !srv 17.07.08 {
     &      term(ix,iy) +                                                    !iyv 19.04.07
     1     (bb(ix,iy,2)/bb(ix,iy,3)+                                   !srv 11.0113
     2     bb(leftix(ix,iy),leftiy(ix,iy),2)/
     3     bb(leftix(ix,iy),leftiy(ix,iy),3))/2.0_R8
     4     *(vbecrb(ix,iy,0)+vbecrb(leftix(ix,iy),leftiy(ix,iy),0))/
     5     2.0_R8*0.5_R8*(vol(ix,iy)/hx(ix,iy)
     6     +vol(leftix(ix,iy),leftiy(ix,iy))/
     7     hx(leftix(ix,iy),leftiy(ix,iy)))                           !srv 17.07.08 }
          flo_he = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+   !srv 13.09.02 { 01.07.05 20.01.11
     &     term(ix,iy) +                                                     !iyv 30.05.08 !srv 03.10.08 
     1     (bb(ix,iy,2)/bb(ix,iy,3)+                                  !srv 11.01.13
     2     bb(leftix(ix,iy),leftiy(ix,iy),2)/
     3     bb(leftix(ix,iy),leftiy(ix,iy),3))/2.0_R8
     4     *(vbecrb(ix,iy,0)+vbecrb(leftix(ix,iy),leftiy(ix,iy),0))/
     5     2.0_R8*0.5_R8*(vol(ix,iy)/hx(ix,iy)
     6     +vol(leftix(ix,iy),leftiy(ix,iy))/
     7     hx(leftix(ix,iy),leftiy(ix,iy)))                           !srv 17.07.08 }
          flo2dia = cvlbhz(ix,iy,0)+pbshz(ix,iy,0)*wrk0(ix,iy)+          !srv 01.10.99 01.07.05
     &     termhz(ix,iy)+                                                !srv 05.11.09
     1     (bb(ix,iy,2)/bb(ix,iy,3)+                                  !srv 11.01.13
     2     bb(leftix(ix,iy),leftiy(ix,iy),2)/
     3     bb(leftix(ix,iy),leftiy(ix,iy),3))/2.0_R8
     4     *(ub2dia(ix,iy,0)+ub2dia(leftix(ix,iy),leftiy(ix,iy),0))/
     5     2.0_R8*0.5_R8*(hz(ix,iy)*vol(ix,iy)/hx(ix,iy)
     6     +hz(leftix(ix,iy),leftiy(ix,iy))*
     7     vol(leftix(ix,iy),leftiy(ix,iy))/
     8     hx(leftix(ix,iy),leftiy(ix,iy)))
          flo_nodrift = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+
     &                  term(ix,iy)   !srv 20.06.08
         else                                                         !srv 21.06.08 {
          flo = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+               !srv 01.07.05
     &      term(ix,iy) +                                                    !iyv 19.04.07
     4      ubdia(ix,iy,0)*
     5      (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/            !srv 01.03.00
     6      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))               !srv 01.03.00
          flo_mdf = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+            !srv 13.10.06 {
     &      term(ix,iy) +                                                    !iyv 19.04.07
     &      vbecrb(ix,iy,0)*
     &      (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     &      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))
          flo_he = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+ !srv 13.09.02 { 01.07.05 20.01.11
     &     term(ix,iy) +                                                     !iyv 30.05.08 !srv 03.10.08 
     4     vbecrb(ix,iy,0)*
     5     (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/               
     6      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))
          flo2dia = cvlbhz(ix,iy,0)+pbshz(ix,iy,0)*wrk0(ix,iy)+          !srv 01.10.99 01.07.05
     &      termhz(ix,iy) +                                                !srv 20.10.09 05.11.09
     4      ub2dia(ix,iy,0)*
     5      (hz(ix,iy)*vol(ix,iy)+hz(leftix(ix,iy),leftiy(ix,iy))
     5                          *vol(leftix(ix,iy),leftiy(ix,iy)))/     !srv 01.03.00 09.01.01
     6      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))                 !srv 01.03.00
          flo_nodrift = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy) +
     &                  xfrhie_no_drift*term(ix,iy)                     !srv 20.06.08 20.10.09
         endif                                                        !srv 21.06.08 }
         if (no_hybr.eq.0) then
           if(no_Ptncr_x_co.eq.0) then                                   !sv 13.10.06 {
            con_mdf(ix,iy)  = hybr(flo_mdf,cdnb(ix,iy,0))                  !srv 13.10.06 }
            con_he(ix,iy)   = hybr(flo_he,cdnb(ix,iy,0))                    !iyv 30.05.08 !srv 03.10.08
            con(ix,iy)      = hybr(flo,cdnb(ix,iy,0))
            con_nodrift     = hybr(flo_nodrift,cdnb(ix,iy,0))
            con_2dia(ix,iy) = hybr(flo2dia,cdnbhz(ix,iy,0))                 !srv 09.01.01
           else
            con_mdf(ix,iy)  = cdnb(ix,iy,0)
            con_he(ix,iy)   = cdnb(ix,iy,0)                  !iyv 30.05.08 !srv 03.10.08
            con(ix,iy)      = cdnb(ix,iy,0)
            con_nodrift     = cdnb(ix,iy,0)
            con_2dia(ix,iy) = cdnbhz(ix,iy,0)                 !srv 09.01.01
           endif                                                         !sv 13.10.06 }
           if(fnb_nodrift_style.eq.0) then                               !srv 20.06.08
            fnb_nodrift(ix,iy,0) = flo_nodrift*
     1      (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &      con_nodrift*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
           else                                                          !srv 20.06.08 {
            fnb_nodrift(ix,iy,0) = flo_nodrift*
     1      (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
           endif                                                         !srv 20.06.08 }
           fnb(ix,iy,0) =
     1      flo*(nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &      con(ix,iy)*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
           fnb_fcor(ix,iy,0) = flo2dia*(nb(leftix(ix,iy),leftiy(ix,iy))  !srv 01.10.99
     &      +nb(ix,iy))/2.0_R8-
     &      con_2dia(ix,iy)*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpbhz(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy))) !srv 01.10.99 09.01.01
           fnb_mdf(ix,iy,0) = flo_mdf*(nb(leftix(ix,iy),leftiy(ix,iy))   !srv 13.10.06 {
     &      +nb(ix,iy))/2.0_R8-
     &      con_mdf(ix,iy)*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))-
     &      fnbPSch(ix,iy,0)                                            !srv 13.10.06 }
*      fnb_he will be used only for heat equations                      !srv 21.01.00 {      
           fnb_he(ix,iy,0) = flo_he*(nb(leftix(ix,iy),leftiy(ix,iy))     !srv 13.09.02
     &      +nb(ix,iy))/2.0_R8-
csw 30sep2014 (is this 5/3 factor correct here? I think not...)
     &      (5.0e0_R8/3.0e0_R8)*
     &      (con_he(ix,iy)*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))+         !srv 13.09.02
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy))) ) !srv 21.01.00 }
           fnb_32(ix,iy,0) =
     1      flo*(nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &      abs(flo/2.0_R8)*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))
           fnb_52(ix,iy,0) = -
     &      (con(ix,iy)-abs(flo/2.0_R8))*
     1      (nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
           fnb_32nd(ix,iy,0) =
     1      flo_nodrift*
     2      (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &      abs(flo_nodrift/2.0_R8)*
     2      (nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))
           fnb_52nd(ix,iy,0) = -
     &      (con_nodrift-abs(flo_nodrift/2.0_R8))*
     1      (nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
           if(mdf_fnb.eq.0) then                                     !srv 13.10.06  {
             flob(ix,iy,0) = flo*(kmpr(leftix(ix,iy),leftiy(ix,iy))+
     &        kmpr(ix,iy))/2.0_R8-
     &        con(ix,iy)*(kmpr(ix,iy)-kmpr(leftix(ix,iy),leftiy(ix,iy)))
             conb(ix,iy,0,0) = con(ix,iy)*                              !srv 13.10.06
     &        (kmpr(leftix(ix,iy),leftiy(ix,iy))+kmpr(ix,iy))/2.0_R8-   !srv 30.01.10
     &        flo*(kmpr(ix,iy)-kmpr(leftix(ix,iy),leftiy(ix,iy)))/4.0_R8
     &        +cdpb(ix,iy,0)
           else
            flob(ix,iy,0) = flo_mdf*(kmpr(leftix(ix,iy),leftiy(ix,iy))
     &         +kmpr(ix,iy))/2.0_R8 - con_mdf(ix,iy)*
     &         (kmpr(ix,iy)-kmpr(leftix(ix,iy),leftiy(ix,iy)))
     &         -cddb(ix,iy,1)*dkmprti(ix,iy,1)                         !srv 21.10.10
            conb(ix,iy,0,0) = con_mdf(ix,iy)*                          !srv 13.10.06
     &         (kmpr(leftix(ix,iy),leftiy(ix,iy))+kmpr(ix,iy))/2.0_R8-
     &          flo_mdf*
     &         (kmpr(ix,iy)-kmpr(leftix(ix,iy),leftiy(ix,iy)))/4.0_R8+
     &          cdpb(ix,iy,0)
            conb(ix,iy,0,1) = cdpb(ix,iy,0)*kmprti(ix,iy,0)
            conb(ix,iy,0,2) = conb(ix,iy,0,1)
           endif
         else  !  SOLPS4.0 basic treatment
           con(ix,iy) = cdnb(ix,iy,0)
           con_nodrift = cdnb(ix,iy,0)
csw 14mar2013 missing con_he added --> needed in b2tfhi_ / b2tfhe_ !!!
           con_he(ix,iy)   = cdnb(ix,iy,0) 
!         T0 = (NI(IX,IY) + NI(IX+1,IY)) / 2
!         T1 = (2 * NI(IX,IY) * NI(IX+1,IY)) /
!     .    (CUTLO + NI(IX,IY) + NI(IX+1,IY))
!         T2 = CVMGP(T0, T1, (UU(IX,IY)+VP(IX,IY))*
!     1    (NI(IX,IY)-NI(IX+1,IY)))
!         FNIX(IX,IY) = (UU(IX,IY)+VP(IX,IY)) * SX(IX,IY) * T2
!         FLOX(IX,IY) = MIN(KAP(IX,IY), KAP(IX+1,IY)) *
!     .    (UU(IX,IY) + VP(IX,IY)) *
!     1    SX(IX,IY) * T2 / (T0 + CUTLO)
!         CONX(IX,IY) = (CCX(IX,IY)+CCX1(IX,IY)) * SX(IX,IY)
           t0 = (nb(ix,iy) + nb(leftix(ix,iy),leftiy(ix,iy))) / 2.0_R8
           t1 = (2.0_R8 * nb(ix,iy) * nb(leftix(ix,iy),leftiy(ix,iy))) /
     .      (cutlo + nb(ix,iy) + nb(leftix(ix,iy),leftiy(ix,iy)))
           t2 = cvmgp(t0, t1, flo*
     1      (nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy))))
           t3 = cvmgp(t0, t1, flo_nodrift*
     1      (nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy))))
           fnb(ix,iy,0) = flo * t2-
     &      con(ix,iy)*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
           fnb_nodrift(ix,iy,0) = flo_nodrift * t3-
     &      con_nodrift*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
           flob(ix,iy,0) = 
     .      min(kap(ix,iy), kap(leftix(ix,iy),leftiy(ix,iy))) *
     .      flo * t2 / (t0 + cutlo)
csw 14mar2013 missing fnb_he added --> needed in b2tfhi_ / b2tfhe_ !!!
           fnb_he(ix,iy,0) = flo_he * t2 -
csw 30sep2014 (is this 5/3 factor correct here? I think not...)
     &      (5.0e0_R8/3.0e0_R8)*
     &      (con_he(ix,iy)*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))+         !srv 13.09.02
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy))) ) !srv 21.01.00 }
!         conb(ix,iy,0,0) = (CCX(IX,IY)+CCX1(IX,IY)) * SX(IX,IY)
           fnb_32(ix,iy,0) = flo * t2
           fnb_52(ix,iy,0) = -
     &      con(ix,iy)*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
           fnb_32nd(ix,iy,0) = flo_nodrift * t3
           fnb_52nd(ix,iy,0) = -
     &      con_nodrift*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
           flob(ix,iy,0) = flo*(kmpr(leftix(ix,iy),leftiy(ix,iy))
     &       +kmpr(ix,iy))/2.0_R8
           conb(ix,iy,0,0) =                                         !srv 13.10.06
     &      con(ix,iy)*(kmpr(leftix(ix,iy),leftiy(ix,iy))+
     &                 kmpr(ix,iy))/2.0_R8 + cdpb(ix,iy,0)
         endif
         if ((alpha.gt.0.0_R8 .and. zb(ix,iy).eq.0.0_R8) .and.       !srv 20.11.08
     *        alpha_.eq.0.0_R8 ) then                                ! neutral flux limit
           t1 = max(flux_limit_min_ti,
     *           ti(ix,iy), ti(leftix(ix,iy),leftiy(ix,iy)))
           vbar = sqrt(8.0_R8 * t1 / (pi * mb))
           if (area_fix.ge.1) then
             t0 = (1.0_R8 + (abs(fnb(ix,iy,0)) / 
     &        (((nb(ix,iy)+nb(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8)*
     &         alpha*vbar/4.0_R8*gs(ix,iy,0)*qc(ix,iy)))**gamma)**
     &        (1.0_R8/gamma)
           else
             t0 = (1.0_R8 + (abs(fnb(ix,iy,0)) / 
     &        (((nb(ix,iy)+nb(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8)*
     &         alpha*vbar/4.0_R8*gs(ix,iy,0)))**gamma)**
     &        (1.0_R8/gamma)
           endif
           maxfluxlimit(0)=max(maxfluxlimit(0),t0)
           fnb(ix,iy,0) = fnb(ix,iy,0) / t0
           fnb_32(ix,iy,0) = fnb_32(ix,iy,0) / t0
           fnb_52(ix,iy,0) = fnb_52(ix,iy,0) / t0
           fnb_nodrift(ix,iy,0) = fnb_nodrift(ix,iy,0) / t0
           fnb_32nd(ix,iy,0) = fnb_32nd(ix,iy,0) / t0
           fnb_52nd(ix,iy,0) = fnb_52nd(ix,iy,0) / t0
           flob(ix,iy,0) = flob(ix,iy,0) / t0
           conb(ix,iy,0,0:4) = conb(ix,iy,0,0:4) / t0                 !srv 13.10.06 {
           fnb_he(ix,iy,0)   = fnb_he(ix,iy,0)/t0
           fnb_fcor(ix,iy,0) = fnb_fcor(ix,iy,0)/t0
           fnb_mdf(ix,iy,0) = fnb_mdf(ix,iy,0)/t0
           fnbPSch(ix,iy,0)  = fnbPSch(ix,iy,0)/t0                    !srv 13.10.06 }
           fllim0fnb(ix,iy,0) = t0
         endif
        endif
       enddo
      enddo

      if(iout_b2npc9.ne.0 .or. iout_b2npco.ne.0) then                 !srv 11.09.09 {
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          wrk1(ix,iy) = 0.0_R8
          wrk2(ix,iy) = 0.0_R8
          wrk3(ix,iy) = 0.0_R8
          wrk4(ix,iy) = 0.0_R8
          wrk5(ix,iy) = 0.0_R8
          wrk7(ix,iy) = 0.0_R8
          wrk8(ix,iy) = 0.0_R8                                        !iys 17.04.2012
         else
          wrk1(ix,iy) = pbs(ix,iy,0)*wrk0(ix,iy)*
     *     (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
          wrk2(ix,iy) = term(ix,iy)*
     *     (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
          wrk3(ix,iy) = vbecrb(ix,iy,0)*
     &     (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     &     (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*
     &     (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
          wrk4(ix,iy) = -con_2dia(ix,iy)*(nb(ix,iy)-
     &     nb(leftix(ix,iy),leftiy(ix,iy)))
          wrk5(ix,iy) = -cdpb(ix,iy,0)*(pb(ix,iy)-
     &     pb(leftix(ix,iy),leftiy(ix,iy)))
          wrk6(ix,iy) = termhz(ix,iy)*
     *     (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
          wrk7(ix,iy) = -con_mdf(ix,iy)*(nb(ix,iy)-
     &     nb(leftix(ix,iy),leftiy(ix,iy)))
          wrk8(ix,iy) = cvlb(ix,iy,0)*                                !iys 17.04.2012
     &     (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
         endif
        enddo
       enddo
       write(chns,'(i3.3)') isb
       if(.not. is_neutral(isb)) then
        call my_out(70,nx,ny,wrk3,'b2tfnb_vaecrbnax'//chns)
        call my_out(70,nx,ny,fnbPSch(-1,-1,0),'b2tfnb_fnbPSchx'//chns)
       endif
       call my_out(70,nx,ny,wrk1,'b2tfnb_bxuanax'//chns)
       call my_out(70,nx,ny,wrk2,'b2tfnb_dpccornax'//chns)
       call my_out(70,nx,ny,wrk4,'b2tfnb_dPat_2diagradnax'//chns)
       call my_out(70,nx,ny,wrk5,'b2tfnb_dgradpbx'//chns)
       call my_out(70,nx,ny,wrk6,'b2tfnb_dpccorhznax'//chns)         !srv 05.11.09
       call my_out(70,nx,ny,wrk7,'b2tfnb_dPat_mdf_gradnax'//chns)
       call my_out(70,nx,ny,wrk8,'b2tfnb_cvlbnax'//chns)             !iys 17.04.2012
!       call my_out(70,nx,ny,con_2dia,'b2tfnb_dPat_2diax'//chns)
!       call my_out(70,nx,ny,con_mdf,'b2tfnb_dPat_mdfx'//chns)
      endif
      if(iout_b2npmo.ne.0) then
       write(chns,'(i3.3)') isb
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          wrk6(ix,iy) = 0.0_R8
         else
          wrk6(ix,iy) = vbdia(ix,iy,0)*
     5     (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     6     (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*
     &     (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8    
         endif
        enddo
       enddo
       if(.not. is_neutral(isb)) 
     &  call my_out(70,nx,ny,wrk6,'b2tfnb_vadianax'//chns)
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2) then
           wrk6(ix,iy) = 0.0_R8
          else
           wrk6(ix,iy) = vbdia(ix,iy,0)*
     5      (hz(ix,iy)*vol(ix,iy)+hz(leftix(ix,iy),leftiy(ix,iy))
     5      *vol(leftix(ix,iy),leftiy(ix,iy)))/
     6      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*
     &      (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8    
          endif
         enddo
        enddo
        if(.not. is_neutral(isb)) 
     &   call my_out(70,nx,ny,wrk6,'b2tfnb_vadianahzx'//chns)
        call my_out(70,nx,ny,fnb_fcor(-1,-1,0),'b2tfnb_fnb_fcorx'//chns)
      endif                                                     !srv 11.09.09 }

*   ..compute fnby, floby, conby
      do iy = -1, ny
       do ix = -1, nx
        if(bottomiy(ix,iy).eq.-2) then
         fnb(ix,iy,1)         = 0.0_R8
         fnb_32(ix,iy,1)      = 0.0_R8
         fnb_52(ix,iy,1)      = 0.0_R8
         fnb_nodrift(ix,iy,1) = 0.0_R8
         fnb_32nd(ix,iy,1)    = 0.0_R8
         fnb_52nd(ix,iy,1)    = 0.0_R8
         flob(ix,iy,1)        = 0.0_R8
         fnb_fcor(ix,iy,1)    = 0.0_R8                                   !srv 01.10.99
         fnb_he(ix,iy,1)      = 0.0_R8                                   !srv 21.01.00
         fnb_mdf(ix,iy,1)     = 0.0_R8                                   !srv 13.10.06
         conb(ix,iy,1,0:2)    = 0.0_R8                                   !srv 13.10.06
         con(ix,iy)           = 0.0_R8
         con_2dia(ix,iy)      = 0.0_R8
         con_he(ix,iy)        = 0.0_R8
         con_mdf(ix,iy)       = 0.0_R8
        else
         if(drift_style.eq.0)then                                     !srv 21.06.08
          if(nint(zn(isb)).eq.1.and.nint(zamax(isb)).eq.1) then          !srv 20.09.06
           flo = cvlb(ix,iy,1)+
     &      0.5_R8*(ubdia(ix,iy,1)*vol(ix,iy)/hy1(ix,iy)
     &      +ubdia(bottomix(ix,iy),bottomiy(ix,iy),1)*
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
     &       /hy1(bottomix(ix,iy),bottomiy(ix,iy)))+
     &        (fchvispar(ix,iy,1)+fchinert(ix,iy,1)+ ! fchvisper(ix,iy,1)+  !srv 17.07.08 {
     &         fchanml(ix,iy,1))*ycur/
     &        (qe*(nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2)        !srv 17.07.08 }
           flo_mdf = cvlb(ix,iy,1)+
     &      0.5_R8*(vbecrb(ix,iy,1)*vol(ix,iy)/hy1(ix,iy)
     &      +vbecrb(bottomix(ix,iy),bottomiy(ix,iy),1)*
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
     &       /hy1(bottomix(ix,iy),bottomiy(ix,iy)))+
     &        (fchvispar(ix,iy,1)+fchinert(ix,iy,1)+ ! fchvisper(ix,iy,1)+
     &         fchanml(ix,iy,1))*ycur/
     &        (qe*(nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2)  !srv 20.09.06
           if(cfvma .eq. 0.0_R8) then                                     !srv 15.12.05
            flo2dia = cvlbhz(ix,iy,1)+
     &       (ub2dia(ix,iy,1)+ub2dia(bottomix(ix,iy),bottomiy(ix,iy),1))
     &       /2.0_R8*0.5_R8*(hz(ix,iy)*vol(ix,iy)/hy1(ix,iy)+          !srv 01.10.99 09.01.01
     &        hz(bottomix(ix,iy),bottomiy(ix,iy))
     5                        *vol(bottomix(ix,iy),bottomiy(ix,iy))/
     &        hy1(bottomix(ix,iy),bottomiy(ix,iy)))+    !srv 01.03.00 09.01.01
     &        (hz(ix,iy)+hz(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8* !srv 18.09.08
     &        (fchvispar(ix,iy,1)+fchinert(ix,iy,1)+ ! fchvisper(ix,iy,1)+
     &         fchanml(ix,iy,1))*ycur/
     &        (qe*(nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2)   !srv 20.09.06
           else
            flo2dia = cvmbhz(ix,iy,1)+
     &       (ub2dia(ix,iy,1)+ub2dia(bottomix(ix,iy),bottomiy(ix,iy),1))
     &       /2.0_R8*0.5_R8*(hz(ix,iy)*vol(ix,iy)/hy1(ix,iy)+          !srv 01.10.99 09.01.01
     &        hz(bottomix(ix,iy),bottomiy(ix,iy))
     5                        *vol(bottomix(ix,iy),bottomiy(ix,iy))/
     &        hy1(bottomix(ix,iy),bottomiy(ix,iy)))    !srv 01.03.00 09.01.01
           endif
          else
           flo = cvlb(ix,iy,1)+
     &      0.5_R8*(ubdia(ix,iy,1)*vol(ix,iy)/hy1(ix,iy)
     &      +ubdia(bottomix(ix,iy),bottomiy(ix,iy),1)*
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
     &       /hy1(bottomix(ix,iy),bottomiy(ix,iy)))
           flo_mdf = cvlb(ix,iy,1)+
     &      0.5_R8*(vbecrb(ix,iy,1)*vol(ix,iy)/hy1(ix,iy)
     &      +vbecrb(bottomix(ix,iy),bottomiy(ix,iy),1)*
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
     &       /hy1(bottomix(ix,iy),bottomiy(ix,iy)))
           if(cfvma .eq. 0.0_R8) then                                     !srv 15.12.05
            flo2dia = cvlbhz(ix,iy,1)+
     &       (ub2dia(ix,iy,1)+ub2dia(bottomix(ix,iy),bottomiy(ix,iy),1))
     &       /2.0_R8*0.5_R8*(hz(ix,iy)*vol(ix,iy)/hy1(ix,iy)+          !srv 01.10.99 09.01.01
     &       hz(bottomix(ix,iy),bottomiy(ix,iy))
     5       *vol(bottomix(ix,iy),bottomiy(ix,iy))/
     &       hy1(bottomix(ix,iy),bottomiy(ix,iy)))    !srv 01.03.00 09.01.01
           else
            flo2dia = cvmbhz(ix,iy,1)+
     &       (ub2dia(ix,iy,1)+ub2dia(bottomix(ix,iy),bottomiy(ix,iy),1))
     &       /2.0_R8*0.5_R8*(hz(ix,iy)*vol(ix,iy)/hy1(ix,iy)+          !srv 01.10.99 09.01.01
     &        hz(bottomix(ix,iy),bottomiy(ix,iy))
     5                        *vol(bottomix(ix,iy),bottomiy(ix,iy))/
     &        hy1(bottomix(ix,iy),bottomiy(ix,iy)))    !srv 01.03.00 09.01.01
           endif
          endif
          flo_he = cvlb(ix,iy,1)+                      !srv 13.09.02 ( 20.01.11
     &      0.5_R8*(vbecrb(ix,iy,1)*vol(ix,iy)/hy1(ix,iy)
     &      +vbecrb(bottomix(ix,iy),bottomiy(ix,iy),1)*
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
     &       /hy1(bottomix(ix,iy),bottomiy(ix,iy)))
         else                                                        !srv 21.06.08 {
          if(nint(zn(isb)).eq.1.and.nint(zamax(isb)).eq.1) then          !srv 20.09.06
           flo = cvlb(ix,iy,1)+ubdia(ix,iy,1)*
     5      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/          !srv 01.03.00
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))+            !srv 01.03.00
     &      (fchvispar(ix,iy,1)+fchinert(ix,iy,1)+ ! fchvisper(ix,iy,1)+
     &      fchanml(ix,iy,1))/
     &      (qe*(nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2) !srv 20.09.06
           flo_mdf = cvlb(ix,iy,1)+vbecrb(ix,iy,1)*                     !srv 01.06.06 {
     5      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))+
     &      (fchvispar(ix,iy,1)+fchinert(ix,iy,1)+ ! fchvisper(ix,iy,1)+
     &      fchanml(ix,iy,1))/
     &      (qe*(nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2)  !srv 20.09.06
           if(cfvma .eq. 0.0_R8) then                                     !srv 15.12.05
            flo2dia = cvlbhz(ix,iy,1)+ub2dia(ix,iy,1)*                   !srv 01.10.99 09.01.01
     5       (hz(ix,iy)*vol(ix,iy)+hz(bottomix(ix,iy),bottomiy(ix,iy))
     5       *vol(bottomix(ix,iy),bottomiy(ix,iy)))/   !srv 01.03.00 09.01.01
     6       (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))+            !srv 01.03.00
     &       (hz(ix,iy)+hz(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8*   !srv 18.09.08
     &       (fchvispar(ix,iy,1)+fchinert(ix,iy,1)+ ! fchvisper(ix,iy,1)+
     &       fchanml(ix,iy,1))/
     &       (qe*(nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2)   !srv 20.09.06
           else
            flo2dia = cvmbhz(ix,iy,1)+ub2dia(ix,iy,1)*                    !srv 01.10.99 09.01.01
     5       (hz(ix,iy)*vol(ix,iy)+hz(bottomix(ix,iy),bottomiy(ix,iy))
     5       *vol(bottomix(ix,iy),bottomiy(ix,iy)))/   !srv 01.03.00 09.01.01
     6       (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))            !srv 01.03.00
           endif
          else
           flo = cvlb(ix,iy,1)+ubdia(ix,iy,1)*
     5      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/          !srv 01.03.00
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))            !srv 01.03.00
           flo_mdf = cvlb(ix,iy,1)+vbecrb(ix,iy,1)*                     !srv 01.06.06 {
     5      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
           if(cfvma .eq. 0.0_R8) then                                     !srv 15.12.05
            flo2dia = cvlbhz(ix,iy,1)+ub2dia(ix,iy,1)*                   !srv 01.10.99 09.01.01
     5       (hz(ix,iy)*vol(ix,iy)+hz(bottomix(ix,iy),bottomiy(ix,iy))
     5       *vol(bottomix(ix,iy),bottomiy(ix,iy)))/   !srv 01.03.00 09.01.01
     6       (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))            !srv 01.03.00
           else
            flo2dia = cvmbhz(ix,iy,1)+ub2dia(ix,iy,1)*                    !srv 01.10.99 09.01.01
     5       (hz(ix,iy)*vol(ix,iy)+hz(bottomix(ix,iy),bottomiy(ix,iy))
     5       *vol(bottomix(ix,iy),bottomiy(ix,iy)))/   !srv 01.03.00 09.01.01
     6       (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))            !srv 01.03.00
           endif
          endif
          flo_he = cvlb(ix,iy,1)+vbecrb(ix,iy,1)*         !srv 13.09.02 ( 20.01.11
     5     (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6     (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))

         endif                                                       !srv 21.06.08 }
         flo_nodrift = cvlb(ix,iy,1)
         if (no_hybr.eq.0) then
          if(no_Ptncr_y_co.eq.0) then                                  !srv 06.07.06 {
           con_he(ix,iy)   = hybr(flo_he,cdnb(ix,iy,1))                   !iyv 30.05.08 !srv 03.10.08
           con_2dia(ix,iy) = hybr(flo2dia,cdnbhz(ix,iy,1))                !srv 01.10.99 09.01.01
           con_mdf(ix,iy)  = hybr(flo_mdf,cdnb(ix,iy,1))                  !srv 13.10.06
           con(ix,iy)      = hybr(flo,cdnb(ix,iy,1))
           con_nodrift     = hybr(flo_nodrift,cdnb(ix,iy,1))
          else
           con_he(ix,iy)   = cdnb(ix,iy,1)                            !iyv 30.05.08 !srv 03.10.08
           con_2dia(ix,iy) = cdnbhz(ix,iy,1)                       !srv 01.10.99 09.01.01
           con_mdf(ix,iy)  = cdnb(ix,iy,1)                         !srv 13.10.06
           con(ix,iy)      = cdnb(ix,iy,1)
           con_nodrift     = cdnb(ix,iy,1)
          endif                                                        !sv 06.07.06 }
          fnb(ix,iy,1) =
     1      flo*(nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &      con(ix,iy)*(nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &      cdpb(ix,iy,1)*(pb(ix,iy)
     &                    -pb(bottomix(ix,iy),bottomiy(ix,iy)))
*      fnb_he will be used only for heat equations                       !srv 21.01.00 {     
          fnb_he(ix,iy,1) = flo_he*(nb(bottomix(ix,iy),bottomiy(ix,iy)) !srv 13.09.02 {
     &          +nb(ix,iy))/2.0_R8-
csw 30sep2014 (is this 5/3 factor correct here? I think not...)
     &    (5.0_R8/3.0_R8)*
     &    (con_he(ix,iy)*(nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))        !srv 13.09.02 }
     &    +cdpb(ix,iy,1)*(pb(ix,iy)
     &                   -pb(bottomix(ix,iy),bottomiy(ix,iy))) )         !srv 21.01.00 }
          if(cfvma .eq. 0.0_R8) then                                    !srv 15.12.05
            fnb_fcor(ix,iy,1) = flo2dia*                                 !srv 01.10.99
     &            (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2-
     &       con_2dia(ix,iy)*(nb(ix,iy)-
     &                        nb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &       cdpbhz(ix,iy,1)*(pb(ix,iy)                                  !srv 09.01.01
     &                     -pb(bottomix(ix,iy),bottomiy(ix,iy)))         !srv 01.10.99
          else                                                          !srv 15.12.05
            fnb_fcor(ix,iy,1) = flo2dia*                                 !srv 01.10.99
     &            (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2
          endif
          if(fnb_nodrift_style.eq.0) then                               !srv 20.06.08
           fnb_nodrift(ix,iy,1) = flo_nodrift*
     1      (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &       con_nodrift*(nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))
     &      -cdpb(ix,iy,1)*(pb(ix,iy)
     &                     -pb(bottomix(ix,iy),bottomiy(ix,iy)))
          else                                                          !srv 20.06.08 {
           fnb_nodrift(ix,iy,1) = flo_nodrift*
     1      (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
          endif                                                         !srv 20.06.08 }
          fnb_mdf(ix,iy,1) = flo_mdf*                                   !srv 13.10.06 {
     &     (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &     con_mdf(ix,iy)*(nb(ix,iy)-
     &                     nb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &     cdpb(ix,iy,1)*(pb(ix,iy)
     &                   -pb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &                    fnbPSch(ix,iy,1)                              !srv 13.10.06 }
          fnb_32(ix,iy,1) =
     1      flo*(nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &      abs(flo/2.0_R8)*
     1      (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))
          fnb_52(ix,iy,1) = -
     1      (con(ix,iy)-abs(flo/2.0_R8))*
     2      (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &      cdpb(ix,iy,1)*(pb(ix,iy)
     &                    -pb(bottomix(ix,iy),bottomiy(ix,iy)))
          fnb_32nd(ix,iy,1) =
     1      flo_nodrift*
     2      (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &      abs(flo_nodrift/2.0_R8)*
     1      (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))
          fnb_52nd(ix,iy,1) = -
     1      (con_nodrift-abs(flo_nodrift/2.0_R8))*
     2      (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &      cdpb(ix,iy,1)*(pb(ix,iy)
     &                    -pb(bottomix(ix,iy),bottomiy(ix,iy)))
          if(mdf_fnb.eq.0) then                                     !srv 13.10.06  {
           flob(ix,iy,1) = flo*(kmpr(bottomix(ix,iy),bottomiy(ix,iy))
     &        +kmpr(ix,iy))/2.0_R8-
     &        con(ix,iy)*(kmpr(ix,iy)-
     &                    kmpr(bottomix(ix,iy),bottomiy(ix,iy)))
           conb(ix,iy,1,0) =                                         !srv 13.10.06
     &        con(ix,iy)*(kmpr(bottomix(ix,iy),bottomiy(ix,iy))
     &        +kmpr(ix,iy))/2.0_R8-
     &        flo*(kmpr(ix,iy)
     &        -kmpr(bottomix(ix,iy),bottomiy(ix,iy)))/4.0_R8+
     &        cdpb(ix,iy,1)
          else
           flob(ix,iy,1) = flo_mdf*
     &       (kmpr(bottomix(ix,iy),bottomiy(ix,iy))+kmpr(ix,iy))/2.0_R8-
     &       con_mdf(ix,iy)*(kmpr(ix,iy)-
     &                       kmpr(bottomix(ix,iy),bottomiy(ix,iy)))
     &         +cddb(ix,iy,0)*dkmprti(ix,iy,0)                         !srv 21.10.10
           conb(ix,iy,1,0) =                                         !srv 13.10.06
     &        con_mdf(ix,iy)*(kmpr(bottomix(ix,iy),bottomiy(ix,iy))
     &        +kmpr(ix,iy))/2.0_R8-
     &        flo_mdf*(kmpr(ix,iy)
     &        -kmpr(bottomix(ix,iy),bottomiy(ix,iy)))/4.0_R8+
     &        cdpb(ix,iy,1)
           conb(ix,iy,1,1) = cdpb(ix,iy,1)*kmprti(ix,iy,1)
           conb(ix,iy,1,2) = conb(ix,iy,1,1)
          endif                                                      !srv 23.06.08 }
         else ! SOLPS4.0 basic treatment
           con(ix,iy) = cdnb(ix,iy,1)
           con_nodrift = cdnb(ix,iy,1)
csw 14mar2013 missing con_he added --> needed in b2tfhi_ / b2tfhe_ !!!
           con_he(ix,iy)   = cdnb(ix,iy,1)     
*  -- COMPUTE FNIY, FLOY, CONY --
           t0 = (nb(ix,iy) + nb(bottomix(ix,iy),bottomiy(ix,iy))) /
     1      2.0_R8
           t1 = 
     .      (2.0_R8 * nb(ix,iy) * nb(bottomix(ix,iy),bottomiy(ix,iy))) /
     .      (cutlo + nb(ix,iy) + nb(bottomix(ix,iy),bottomiy(ix,iy)))
           t2 = cvmgp(t0, t1, 
     .      flo*(nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy))))
           t3 = cvmgp(t0, t1, 
     .      flo_nodrift*(nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy))))
           fnb(ix,iy,1) = flo * t2-
     &      con(ix,iy)*(nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &      cdpb(ix,iy,1)*(pb(ix,iy)
     &                    -pb(bottomix(ix,iy),bottomiy(ix,iy)))
           fnb_nodrift(ix,iy,1) = flo_nodrift*t3-
     &       con_nodrift*(nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))
     &      -cdpb(ix,iy,1)*(pb(ix,iy)
     &                     -pb(bottomix(ix,iy),bottomiy(ix,iy)))
           flob(ix,iy,1) = min(kap(ix,iy), 
     .                         kap(bottomix(ix,iy),bottomiy(ix,iy))) *
     .      flo * t2 / (t0 + cutlo)
csw 14mar2013 missing fnb_he added --> needed in b2tfhi_ / b2tfhe_ !!!
           fnb_he(ix,iy,1) = flo_he * t2 -
csw 30sep2014 (is this 5/3 factor correct here? I think not...)
     &      (5.0_R8/3.0_R8)*
     &      (con_he(ix,iy)*(nb(ix,iy)-
     &      nb(bottomix(ix,iy),bottomiy(ix,iy)))+        !srv 13.09.02 }
     &      cdpb(ix,iy,1)*(pb(ix,iy)
     &      -pb(bottomix(ix,iy),bottomiy(ix,iy))) )         !srv 21.01.00 }
!          CONY(IX,IY) = CCY(IX,IY) * SY(IX,IY)
           fnb_32(ix,iy,1) = flo * t2
           fnb_52(ix,iy,1) = -
     &      con(ix,iy)*(nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &      cdpb(ix,iy,1)*(pb(ix,iy)
     &                    -pb(bottomix(ix,iy),bottomiy(ix,iy)))
           fnb_32nd(ix,iy,1) = flo_nodrift*t3
           fnb_52nd(ix,iy,1) = -
     &       con_nodrift*(nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))
     &      -cdpb(ix,iy,1)*(pb(ix,iy)
     &                     -pb(bottomix(ix,iy),bottomiy(ix,iy)))
           conb(ix,iy,1,0) =                                         !srv 13.10.06
     &      con(ix,iy)*(kmpr(bottomix(ix,iy),bottomiy(ix,iy))+
     &      kmpr(ix,iy))/2.0_R8+cdpb(ix,iy,1)
         endif
         if ((alpha.gt.0.0_R8 .and. zb(ix,iy).eq.0.0_R8) .and.       !srv 20.11.08
     *       alpha_.eq.0.0_R8 ) then                                 ! neutral flux limit
          t1 = max(flux_limit_min_ti,
     *           ti(ix,iy), ti(bottomix(ix,iy),bottomiy(ix,iy)))
          vbar = sqrt(8.0_R8 * t1 / (pi * mb))
          t0 = (1.0_R8 + (abs(fnb(ix,iy,1)) / 
     &      (((nb(ix,iy)+nb(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8)*
     &       alpha*vbar/4.0_R8*gs(ix,iy,1)))**gamma)**(1.0_R8/gamma)
          maxfluxlimit(1)=max(maxfluxlimit(1),t0)
          fnb(ix,iy,1) = fnb(ix,iy,1) / t0
          fnb_32(ix,iy,1) = fnb_32(ix,iy,1) / t0
          fnb_52(ix,iy,1) = fnb_52(ix,iy,1) / t0
          fnb_nodrift(ix,iy,1) = fnb_nodrift(ix,iy,1) / t0
          fnb_32nd(ix,iy,1) = fnb_32nd(ix,iy,1) / t0
          fnb_52nd(ix,iy,1) = fnb_52nd(ix,iy,1) / t0
          flob(ix,iy,1) = flob(ix,iy,1) / t0
          conb(ix,iy,1,0:4) = conb(ix,iy,1,0:4) / t0                 !srv 13.10.06 {
          fnb_he(ix,iy,1)   = fnb_he(ix,iy,1)/t0
          fnb_fcor(ix,iy,1) = fnb_fcor(ix,iy,1)/t0
          fnb_mdf(ix,iy,1) = fnb_mdf(ix,iy,1)/t0
          fnbPSch(ix,iy,1)  = fnbPSch(ix,iy,1)/t0                    !srv 13.10.06 }
          fllim0fnb(ix,iy,1) = t0
         endif
        endif
       enddo
      enddo

      if (alpha.gt.0.0_R8 .and. maxval(zb).eq.0.0_R8 .and. ncall.lt.100) 
     > then
        write(*,*) 'b2tfnb: maxfluxlimit for zb = ',maxval(zb),
     1   ' mb = ',mb/mp,' are ',maxfluxlimit
      endif

      if(iout_b2npc9.ne.0 .or. iout_b2npco.ne.0) then                 !srv 11.09.09 {
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).eq.-2) then
          wrk3(ix,iy) = 0.0_R8
          wrk4(ix,iy) = 0.0_R8
          wrk5(ix,iy) = 0.0_R8
          wrk7(ix,iy) = 0.0_R8
          wrk8(ix,iy) = 0.0_R8                                        !iys 17.04.2012
         else
          wrk3(ix,iy) = vbecrb(ix,iy,1)*
     5     (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6     (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &     (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
          wrk4(ix,iy) = -con_2dia(ix,iy)*
     &                   (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))
          wrk5(ix,iy) = -cdpb(ix,iy,1)*(pb(ix,iy)
     &                    -pb(bottomix(ix,iy),bottomiy(ix,iy)))
          wrk7(ix,iy) = -con_mdf(ix,iy)*
     &                   (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))
          wrk8(ix,iy) = cvlb(ix,iy,1)*                             !iys 17.04.2012
     &     (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
         endif
        enddo
       enddo
       write(chns,'(i3.3)') isb
       if(.not. is_neutral(isb)) then
        call my_out(70,nx,ny,fnbPSch(-1,-1,1),'b2tfnb_fnbPSchy'//chns)
        call my_out(70,nx,ny,wrk3,'b2tfnb_vaecrbnay'//chns)
       endif
        call my_out(70,nx,ny,wrk4,'b2tfnb_dPat_2diagradnay'//chns)
       call my_out(70,nx,ny,wrk5,'b2tfnb_dgradpby'//chns)
        call my_out(70,nx,ny,wrk7,'b2tfnb_dPat_mdf_gradnay'//chns)
        call my_out(70,nx,ny,wrk8,'b2tfnb_cvlbnay'//chns)               !iys 17.04.2012
!        call my_out(70,nx,ny,con_2dia,'b2tfnb_dPat_2diay'//chns)
!         call my_out(70,nx,ny,con_mdf,'b2tfnb_dPat_mdfy'//chns)
      endif
      if(iout_b2npmo.ne.0) then
       write(chns,'(i3.3)') isb
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).eq.-2) then
           wrk6(ix,iy) = 0.0_R8
         else
           wrk6(ix,iy) = vbdia(ix,iy,1)*
     5      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &      (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
         endif
        enddo
       enddo
       if(.not. is_neutral(isb)) 
     &  call my_out(70,nx,ny,wrk6,'b2tfnb_vadianay'//chns)
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).eq.-2) then
          wrk6(ix,iy) = 0.0_R8
         else
          wrk6(ix,iy) = vbdia(ix,iy,1)*
     5     (hz(ix,iy)*vol(ix,iy)+hz(bottomix(ix,iy),bottomiy(ix,iy))
     5     *vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6     (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &     (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
         endif
        enddo
       enddo
       if(.not. is_neutral(isb)) 
     &  call my_out(70,nx,ny,wrk6,'b2tfnb_vadianahzy'//chns)
       call my_out(70,nx,ny,fnb_fcor(-1,-1,1),'b2tfnb_fnb_fcory'//chns)
      endif                                                     !srv 11.09.09 }

      if(guard_flows.eq.0) then                                        !xpb
cxpb turn off flows between boundary cells
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            fnb(ix,iy,1)=0.0_R8
            fnb_32(ix,iy,1)=0.0_R8
            fnb_52(ix,iy,1)=0.0_R8
            fnb_nodrift(ix,iy,1)=0.0_R8
            fnb_32nd(ix,iy,1)=0.0_R8
            fnb_52nd(ix,iy,1)=0.0_R8
            flob(ix,iy,1)=0.0_R8
            vbecrb(ix,iy,1)=0.0_R8
            vbdia(ix,iy,1)=0.0_R8
            wbdia(ix,iy,1)=0.0_R8
            ubdia(ix,iy,1)=0.0_R8
          endif
         enddo
        enddo

cxpb turn off poloidal flows in core boundary
        do ix = -1, nx
          if(mod(region(ix,-1,0),4).eq.1) then
            fnb(ix,-1,0) = 0.0_R8
            fnb_32(ix,-1,0) = 0.0_R8
            fnb_52(ix,-1,0) = 0.0_R8
            fnb_nodrift(ix,-1,0) = 0.0_R8
            fnb_32nd(ix,-1,0) = 0.0_R8
            fnb_52nd(ix,-1,0) = 0.0_R8
            flob(ix,-1,0) = 0.0_R8
            vbecrb(ix,-1,0) = 0.0_R8
            vbdia(ix,-1,0) = 0.0_R8
            wbdia(ix,-1,0) = 0.0_R8
            ubdia(ix,-1,0) = 0.0_R8
          endif
        enddo
       elseif(guard_flows.eq.1) then                                  !srv 20.06.08 {
*   .. do nothing
       elseif(guard_flows.eq.2) then
*     turn off flows between boundary cells in west and east
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            fnb(ix,iy,1)        = 0.0_R8
!            fnb_32(ix,iy,1)=0.0_R8
!            fnb_52(ix,iy,1)=0.0_R8
!            fnb_nodrift(ix,iy,1)=0.0_R8
!            fnb_32nd(ix,iy,1)=0.0_R8
!            fnb_52nd(ix,iy,1)=0.0_R8
!            flob(ix,iy,1)=0.0_R8
!            vbecrb(ix,iy,1)=0.0_R8
!            vbdia(ix,iy,1)=0.0_R8
!            wbdia(ix,iy,1)=0.0_R8
!            ubdia(ix,iy,1)=0.0_R8
            fnb_nodrift(ix,iy,1) = 0.0_R8
            fnb_fcor(ix,iy,1)    = 0.0_R8
            fnb_he(ix,iy,1)      = 0.0_R8
            fnb_mdf(ix,iy,1)     = 0.0_R8
            flob(ix,iy,1)        = 0.0_R8
            conb(ix,iy,1,0:2)    = 0.0_R8
          endif
         enddo
        enddo
*     turn off flows between boundary cells in south and in north
         do iy = -1, ny, ny+1
          fnb(:,iy,0)         = 0.0_R8
!          fnb_32(:,iy,0) = 0.0_R8
!          fnb_52(:,iy,0) = 0.0_R8
!          fnb_nodrift(:,iy,0) = 0.0_R8
!          fnb_32nd(:,iy,0) = 0.0_R8
!          fnb_52nd(:,iy,0) = 0.0_R8
!          flob(:,iy,0) = 0.0_R8
!          vbecrb(:,iy,0) = 0.0_R8
!          vbdia(:,iy,0) = 0.0_R8
!          wbdia(:,iy,0) = 0.0_R8
!          ubdia(:,iy,0) = 0.0_R8
          fnb_nodrift(:,iy,0) = 0.0_R8
          fnb_fcor(:,iy,0)    = 0.0_R8
          fnb_he(:,iy,0)      = 0.0_R8
          fnb_mdf(:,iy,0)     = 0.0_R8
          flob(:,iy,0)        = 0.0_R8
          conb(:,iy,0,0:2)    = 0.0_R8
         enddo                                                        !srv 20.06.08 }
      endif                                                           !xpb

csw 18feb2014 remove fluxes between guard cells 
        do ix=-1,nx
         do iy=-1,ny
          if(iffboundx(ix,iy).ne.0) then
            fnb(ix,iy,0)=0.0_R8
            fnb_32(ix,iy,0)=0.0_R8
            fnb_52(ix,iy,0)=0.0_R8
            fnb_nodrift(ix,iy,0)=0.0_R8
            fnb_32nd(ix,iy,0)=0.0_R8
            fnb_52nd(ix,iy,0)=0.0_R8
            fnb_fcor(ix,iy,0)=0.0_R8
            fnb_he(ix,iy,0)=0.0_R8
            fnb_mdf(ix,iy,0)=0.0_R8
            flob(ix,iy,0)=0.0_R8
            vbecrb(ix,iy,0)=0.0_R8
            vbdia(ix,iy,0)=0.0_R8
            wbdia(ix,iy,0)=0.0_R8
            ubdia(ix,iy,0)=0.0_R8
            conb(ix,iy,0,:)=0.0_R8
          endif
          if(iffboundy(ix,iy).ne.0) then
            fnb(ix,iy,1)=0.0_R8
            fnb_32(ix,iy,1)=0.0_R8
            fnb_52(ix,iy,1)=0.0_R8
            fnb_nodrift(ix,iy,1)=0.0_R8
            fnb_32nd(ix,iy,1)=0.0_R8
            fnb_52nd(ix,iy,1)=0.0_R8
            fnb_fcor(ix,iy,1)=0.0_R8
            fnb_he(ix,iy,1)=0.0_R8
            fnb_mdf(ix,iy,1)=0.0_R8
            flob(ix,iy,1)=0.0_R8
            vbecrb(ix,iy,1)=0.0_R8
            vbdia(ix,iy,1)=0.0_R8
            wbdia(ix,iy,1)=0.0_R8
            ubdia(ix,iy,1)=0.0_R8
            conb(ix,iy,1,:)=0.0_R8
          endif
         enddo
        enddo

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            vbecrb(ix,iy,0)=0.0_R8
            vbdia(ix,iy,0)=0.0_R8
            wbdia(ix,iy,0)=0.0_R8
            ubdia(ix,iy,0)=0.0_R8
            fnb(ix,iy,0)=0.0_R8
            fnb_32(ix,iy,0)=0.0_R8
            fnb_52(ix,iy,0)=0.0_R8
            fnb_nodrift(ix,iy,0)=0.0_R8
            fnb_32nd(ix,iy,0)=0.0_R8
            fnb_52nd(ix,iy,0)=0.0_R8
            flob(ix,iy,0)=0.0_R8
            fnb(ix,iy,1)=0.0_R8
            fnb_32(ix,iy,1)=0.0_R8
            fnb_52(ix,iy,1)=0.0_R8
            fnb_nodrift(ix,iy,1)=0.0_R8
            fnb_32nd(ix,iy,1)=0.0_R8
            fnb_52nd(ix,iy,1)=0.0_R8
            flob(ix,iy,1)=0.0_R8
            vbecrb(ix,iy,1)=0.0_R8
            vbdia(ix,iy,1)=0.0_R8
            wbdia(ix,iy,1)=0.0_R8
            ubdia(ix,iy,1)=0.0_R8
          end do
        end do
      end do

* ..compute kinetic energy of species
      if(drift_style.eq.0) then                                      !srv 14.07.10 {
       do ix=-1,nx
        do iy=-1,ny
          kbnrgy(ix,iy) = 0.5_R8*mb*(ub(ix,iy)**2+  !all values are calculated in the cell center
     .     (wbdia(ix,iy,0)+vbecrb(ix,iy,0))**2+
     .     (wbdia(ix,iy,1)+vbecrb(ix,iy,1))**2)
        end do
       end do
      else
       if((maxval(facdrift).ne.0.0_R8 .or. maxval(fac_ExB).ne.0.0_R8)
     &    .and. minval(zb).gt.0.0_R8) then
        wrk0=0.0_R8  ! compute elecrical field on the cell faces
        do iy=0,ny
         do ix=-1,nx
          wrk0(ix,iy)=(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))/
     &          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*2.0_R8
         end do
        end do
        do iy=-1,ny  ! compute perpendicular projection V_ExB on the cell centers
         do ix=-1,nx
          if(bottomiy(ix,iy).eq.-2) then
           vbecrb_(ix,iy,0) = wrk0(ix,topiy(ix,iy))
          elseif(topiy(ix,iy).eq.ny+1) then
           vbecrb_(ix,iy,0) = wrk0(ix,iy)
          else
           vbecrb_(ix,iy,0) = 0.5_R8*
     &        (wrk0(topix(ix,iy),topiy(ix,iy))+wrk0(ix,iy))
          endif
            vbecrb_(ix,iy,0)=-vbecrb_(ix,iy,0)*
     &       bb(ix,iy,2)/bb(ix,iy,3)**2
         enddo
        enddo
        wrk0=0.0_R8  ! compute .... on the cell faces
        do iy=0,ny
         do ix=-1,nx
          wrk0(ix,iy)=(nb(ix,iy)*ti(ix,iy)-
     &                 nb(bottomix(ix,iy),bottomiy(ix,iy))*
     &                 ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*2.0_R8
         end do
        end do
        do iy=-1,ny  ! compute perpendicular projection V_ExB on the cell centers
         do ix=-1,nx
          if(bottomiy(ix,iy).eq.-2) then
           wbdia_(ix,iy,0)=wrk0(ix,topiy(ix,iy))
          elseif(topiy(ix,iy).eq.ny+1) then
           wbdia_(ix,iy,0)=wrk0(ix,iy)
          else
           wbdia_(ix,iy,0)=0.5_R8*
     &        (wrk0(topix(ix,iy),topiy(ix,iy))+
     &         wrk0(ix,iy))
          endif
            wbdia_(ix,iy,0)=-wbdia_(ix,iy,0)/
     &         (zb(ix,iy)*qe*nb(ix,iy))*bb(ix,iy,3)/bb(ix,iy,3)**2
         enddo
        enddo

        wrk0=0.0_R8  ! compute y-component of elecrical field on the cell faces
        do iy=-1,ny
         do ix=-1,nx
          if(leftix(ix,iy).ne.-2) then
           wrk0(ix,iy)=(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))/
     &          (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0_R8
          endif
         end do
        end do
        do iy=-1,ny  ! compute .... on the cell centers
         do ix=-1,nx
          if(leftix(ix,iy).eq.-2) then
           vbecrb_(ix,iy,1)=wrk0(rightix(ix,iy),iy)
          elseif(rightix(ix,iy).eq.nx+1) then
           vbecrb_(ix,iy,1)=wrk0(ix,iy)
          else
           vbecrb_(ix,iy,1)=0.5_R8*
     &        (wrk0(rightix(ix,iy),rightiy(ix,iy))+wrk0(ix,iy))
          endif
            vbecrb_(ix,iy,1)= vbecrb_(ix,iy,1)*
     &       bb(ix,iy,2)/bb(ix,iy,3)**2    !srv 15.02.12
         enddo
        enddo

        wrk0=0.0_R8  ! compute ... on the cell faces
        do iy=-1,ny
         do ix=-1,nx
          if(leftix(ix,iy).ne.-2) then
           wrk0(ix,iy)=(nb(ix,iy)*ti(ix,iy)-
     &         nb(leftix(ix,iy),leftiy(ix,iy))*
     &         ti(leftix(ix,iy),leftiy(ix,iy)))/
     &        (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0_R8
          endif
         enddo
        enddo
        do iy=-1,ny  ! compute .... on the cell centers
         do ix=-1,nx
          if(leftix(ix,iy).eq.-2) then
           wbdia_(ix,iy,1)=wrk0(rightix(ix,iy),iy)
          elseif(rightix(ix,iy).eq.nx+1) then
           wbdia_(ix,iy,1)=wrk0(ix,iy)
          else
           wbdia_(ix,iy,1)=0.5_R8*
     &        (wrk0(rightix(ix,iy),rightiy(ix,iy))+wrk0(ix,iy))
          endif
          wbdia_(ix,iy,1)=wbdia_(ix,iy,1)/
     &        (zb(ix,iy)*qe*nb(ix,iy))*bb(ix,iy,2)/bb(ix,iy,3)**2       !srv 15.02.12
         enddo
        enddo
       else
        wbdia_ = 0.0_R8
        vbecrb_= 0.0_R8
       endif

       do ix=-1,nx
        do iy=-1,ny
          kbnrgy(ix,iy) = 0.5_R8*mb*(ub(ix,iy)**2+   ! ua is a value in the cell center
     .      (xwdia*wbdia_(ix,iy,0)+xvecrb*vbecrb_(ix,iy,0))**2+   ! wbdia_ is a value in the cell center
     .      (xwdia*wbdia_(ix,iy,1)+xvecrb*vbecrb_(ix,iy,1))**2)   ! vbecrb_ is a value in the cell center
        end do
       end do

      endif                                                          !srv 14.07.10 }

      if(iout.ne.0) then                                             !srv 18.06.08 {
       write (chns,'(i3.3)') isb
       call my_out(70,nx,ny,fnb(-1,-1,0),'b2tfnb_fnbx'//chns)
       call my_out(70,nx,ny,fnb(-1,-1,1),'b2tfnb_fnby'//chns)
       call my_out(70,nx,ny,fnb_nodrift(-1,-1,0),
     &             'b2tfnb_fnb_nodriftx'//chns)
       call my_out(70,nx,ny,fnb_nodrift(-1,-1,1),
     &             'b2tfnb_fnb_nodrifty'//chns)
       call my_out(70,nx,ny,fnb_fcor(-1,-1,0),'b2tfnb_fnb_fcorx'//chns)
       call my_out(70,nx,ny,fnb_fcor(-1,-1,1),'b2tfnb_fnb_fcory'//chns)
       call my_out(70,nx,ny,fnb_he(-1,-1,0),'b2tfnb_fnb_hex'//chns)
       call my_out(70,nx,ny,fnb_he(-1,-1,1),'b2tfnb_fnb_hey'//chns)
       call my_out(70,nx,ny,fnbPSch(-1,-1,0),'b2tfnb_fnbPSchx'//chns)
       call my_out(70,nx,ny,fnbPSch(-1,-1,1),'b2tfnb_fnbPSchy'//chns)
       call my_out(70,nx,ny,fnb_mdf(-1,-1,0),'b2tfnb_fnb_mdfx'//chns)
       call my_out(70,nx,ny,fnb_mdf(-1,-1,1),'b2tfnb_fnb_mdfy'//chns)
       call my_out(70,nx,ny,fnb_32(-1,-1,0),'b2tfnb_fnb_32x'//chns)
       call my_out(70,nx,ny,fnb_32(-1,-1,1),'b2tfnb_fnb_32y'//chns)
       call my_out(70,nx,ny,fnb_52(-1,-1,0),'b2tfnb_fnb_52x'//chns)
       call my_out(70,nx,ny,fnb_52(-1,-1,1),'b2tfnb_fnb_52y'//chns)
       call my_out(70,nx,ny,fnb_32nd(-1,-1,0),
     &             'b2tfnb_fnb_32ndx'//chns)
       call my_out(70,nx,ny,fnb_32nd(-1,-1,1),
     &             'b2tfnb_fnb_32ndy'//chns)
       call my_out(70,nx,ny,fnb_52nd(-1,-1,0),
     &             'b2tfnb_fnb_52ndx'//chns)
       call my_out(70,nx,ny,fnb_52nd(-1,-1,1),
     &             'b2tfnb_fnb_52ndy'//chns)
       call my_out(70,nx,ny,fllim0fnb(-1,-1,0),
     &             'b2tfnb_fllim0fnbx'//chns)
       call my_out(70,nx,ny,fllim0fnb(-1,-1,1),
     &             'b2tfnb_fllim0fnby'//chns)
       call my_out(70,nx,ny,vbecrb(-1,-1,0),'b2tfnb_vbecrbx'//chns)
       call my_out(70,nx,ny,vbecrb(-1,-1,1),'b2tfnb_vbecrby'//chns)
       call my_out(70,nx,ny,ubdia(-1,-1,0),'b2tfnb_ubdiax'//chns)
       call my_out(70,nx,ny,ubdia(-1,-1,1),'b2tfnb_ubdiay'//chns)
       call my_out(70,nx,ny,vbdia(-1,-1,0),'b2tfnb_vbdiax'//chns)
       call my_out(70,nx,ny,vbdia(-1,-1,1),'b2tfnb_vbdiay'//chns)
       call my_out(70,nx,ny,wbdia(-1,-1,0),'b2tfnb_wbdiax'//chns)
       call my_out(70,nx,ny,wbdia(-1,-1,1),'b2tfnb_wbdiay'//chns)
       call my_out(70,nx,ny,kbnrgy,'b2tfnb_kbnrgy'//chns)
       call my_out(70,nx,ny,cvlb(-1,-1,0),'b2tfnb_cvlbx'//chns)
       call my_out(70,nx,ny,cvlb(-1,-1,1),'b2tfnb_cvlby'//chns)
       call my_out(70,nx,ny,cvlbhz(-1,-1,0),'b2tfnb_cvlbhzx'//chns)
       call my_out(70,nx,ny,cvlbhz(-1,-1,1),'b2tfnb_cvlbhzy'//chns)
       call my_out(70,nx,ny,cdpbhz(-1,-1,0),'b2tfnb_cdpbhzx'//chns)
       call my_out(70,nx,ny,cdpbhz(-1,-1,1),'b2tfnb_cdpbhzy'//chns)
       call my_out(70,nx,ny,cdpb(-1,-1,0),'b2tfnb_cdpbx'//chns)
       call my_out(70,nx,ny,cdpb(-1,-1,1),'b2tfnb_cdpby'//chns)
       call my_out(70,nx,ny,cdnb(-1,-1,0),'b2tfnb_cdnbx'//chns)
       call my_out(70,nx,ny,cdnb(-1,-1,1),'b2tfnb_cdnby'//chns)
       call my_out(70,nx,ny,cddb(-1,-1,0),'b2tfnb_cddbx'//chns)
       call my_out(70,nx,ny,cddb(-1,-1,1),'b2tfnb_cddby'//chns)
       call my_out(70,nx,ny,pbshz(-1,-1,0),'b2tfnb_pbshzx'//chns)
       call my_out(70,nx,ny,pbshz(-1,-1,1),'b2tfnb_pbshzy'//chns)
       call my_out(70,nx,ny,ub2dia(-1,-1,0),'b2tfnb_ub2diax'//chns)
       call my_out(70,nx,ny,ub2dia(-1,-1,1),'b2tfnb_ub2diay'//chns)
       call my_out(70,nx,ny,conb(-1,-1,0,0),'b2tfnb_conbx0'//chns)
       call my_out(70,nx,ny,conb(-1,-1,0,1),'b2tfnb_conbx1'//chns)
       call my_out(70,nx,ny,conb(-1,-1,0,2),'b2tfnb_conbx2'//chns)
       call my_out(70,nx,ny,conb(-1,-1,0,3),'b2tfnb_conbx3'//chns)
       call my_out(70,nx,ny,conb(-1,-1,0,4),'b2tfnb_conbx4'//chns)
       call my_out(70,nx,ny,conb(-1,-1,1,0),'b2tfnb_conby0'//chns)
       call my_out(70,nx,ny,conb(-1,-1,1,1),'b2tfnb_conby1'//chns)
       call my_out(70,nx,ny,conb(-1,-1,1,2),'b2tfnb_conby2'//chns)
       call my_out(70,nx,ny,conb(-1,-1,1,3),'b2tfnb_conby3'//chns)
       call my_out(70,nx,ny,conb(-1,-1,1,4),'b2tfnb_conby4'//chns)
       call my_out(70,nx,ny,flob(-1,-1,0),'b2tfnb_flobx'//chns)
       call my_out(70,nx,ny,flob(-1,-1,1),'b2tfnb_floby'//chns)
       call my_out(70,nx,ny,kmprti(-1,-1,0),'b2tfnb_kmprtix'//chns)
       call my_out(70,nx,ny,kmprti(-1,-1,1),'b2tfnb_kmprtiy'//chns)
       call my_out(70,nx,ny,dkmprti(-1,-1,0),'b2tfnb_dkmprtix'//chns)
       call my_out(70,nx,ny,dkmprti(-1,-1,1),'b2tfnb_dkmprtiy'//chns)
       call my_out(70,nx,ny,ub,'b2tfnb_ua'//chns)
       call my_out(70,nx,ny,nb,'b2tfnb_na'//chns)
       call my_out(70,nx,ny,te,'b2tfnb_te'//chns)
       call my_out(70,nx,ny,ti,'b2tfnb_ti'//chns)
       call my_out(70,nx,ny,po,'b2tfnb_po'//chns)
       call my_out(70,nx,ny,zb,'b2tfnb_zb'//chns)
       call my_out(70,nx,ny,facdrift,'b2tfnb_facdrift'//chns)
       call my_out(70,nx,ny,fac_ExB,'b2tfnb_fac_ExB'//chns)
       call my_out(70,nx,ny,fchvispar(-1,-1,0),
     &             'b2tfnb_fchvisparx'//chns)
       call my_out(70,nx,ny,fchvispar(-1,-1,1),
     &             'b2tfnb_fchvispary'//chns)
       call my_out(70,nx,ny,fchinert(-1,-1,0),
     &             'b2tfnb_fchinertx'//chns)
       call my_out(70,nx,ny,fchinert(-1,-1,1),
     &             'b2tfnb_fchinerty'//chns)
       call my_out(70,nx,ny,fchanml(-1,-1,0),'b2tfnb_fchanmlx'//chns)
       call my_out(70,nx,ny,fchanml(-1,-1,1),'b2tfnb_fchanmly'//chns)
      endif                                                          !srv 18.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfnb

      end subroutine b2tfnb


      function cvmgp (x1,x2,x3)
      use b2mod_types
      implicit none
      real (kind=R8) :: cvmgp
      real (kind=R8) :: x1,x2,x3
*     ------------------------------------------------------------------
*     CVMGP( CRAY - routinen - ersatz).(RFS 27.5.1992)
*     x1 returned if x3 >= 0
*     x2 returned if x3 < 0
*     ------------------------------------------------------------------
      real (kind=R8) :: t
*     ------------------------------------------------------------------
      if (x3.ge.0.0_R8) then
       t = x1
      else
       t = x2
      endif
      cvmgp = t
      return
*     ------------------------------------------------------------------
      end function cvmgp
