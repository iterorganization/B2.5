*-----------------------------------------------------------------------
*.specification

      subroutine b2tfnb (nx, ny, isb, ismain,                           !srv 08.09.01 !lk 02.06.11
     &  guard_flows, facdrift, fac_ExB, mb, vol, cssb, zb,              !xpb !iyv 06.02.14
     &  pbs, pbshz, nb, pb, ub, ubdia, vbecrb, vbdia, wbdia, ne,        !xpb !srv 09.01.01
     &  te, ti, bb, hx, hy, hz, gs, qc, po,                             !srv 09.01.01 13.10.06
     &  cvlb, cdnb, cdpb, cddb, cvlbhz, cdnbhz, cdpbhz,                 !srv 09.01.01
     &  fchvispar_b, fchvisper_b, fchinert_b, fchanml_b, fchvisq_b,     !srv 20.09.06 08.09.21
     &  fchin_b, fnb, fnb_nodrift, fnb_fcor, fnb_he, fnbPSch, fnb_mdf,  !srv 13.10.06, IYS 17.10.2023
     &  flob, conb, fnb_32, fnb_52, fnb_32nd,
     &  fnb_52nd, fnb_53, kbnrgy, fllim0fnb,
     &  balance)                                                        !djm Jan2017
      use b2mod_types
      use b2mod_math
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_geo
     &   , only : hy1, hyl, hxb, hzl, hzb
      use b2mod_balance                                                 !djm Jan2017
     & , only : fna_pinch, fna_pll, fna_drift, fna_nanom, fna_panom,
     &          fna_ch, fna_pschused, balance_netcdf
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb, guard_flows, ismain                          !xpb !srv 08.09.01 !lk 02.06.11
      real (kind=R8) ::
     *  pbs(-1:nx,-1:ny,0:1), pbshz(-1:nx,-1:ny,0:1), nb(-1:nx,-1:ny),  !srv 09.01.01
     &  pb(-1:nx,-1:ny),
     &  ub(-1:nx,-1:ny), cvlb(-1:nx,-1:ny,0:1),
     &  cdnb(-1:nx,-1:ny,0:1), cdpb(-1:nx,-1:ny,0:1),
     &  ne(-1:nx,-1:ny),
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), hz(-1:nx,-1:ny),              !srv 09.01.01
     &  po(-1:nx,-1:ny), bb(-1:nx,-1:ny,0:3),
     2  gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny),
     2  facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny), zb(-1:nx,-1:ny),
     3  vol(-1:nx,-1:ny), te(-1:nx,-1:ny), ti(-1:nx,-1:ny), mb,
     &  cvlbhz(-1:nx,-1:ny,0:1),                                        !srv 09.01.01
     &  cdnbhz(-1:nx,-1:ny,0:1), cdpbhz(-1:nx,-1:ny,0:1),               !srv 09.01.01
     &  cddb(-1:nx,-1:ny,0:1), cssb(-1:nx,-1:ny,0:1),                   !srv 13.10.06 !iyv 06.02.14
     &  fchvispar_b(-1:nx,-1:ny,0:1),
     &  fchinert_b(-1:nx,-1:ny,0:1),                                    !srv 20.09.06
     &  fchvisper_b(-1:nx,-1:ny,0:1),
     &  fchvisq_b(-1:nx,-1:ny,0:1),
     &  fchanml_b(-1:nx,-1:ny,0:1),                                     !srv 20.09.06
     &  fchin_b(-1:nx,-1:ny,0:1)                                        !IYS 17.10.2023
      logical balance                                                   !djm Jan2017
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  fnb(-1:nx,-1:ny,0:1), fnb_32(-1:nx,-1:ny,0:1),
     1  fnb_52(-1:nx,-1:ny,0:1), fnb_nodrift(-1:nx,-1:ny,0:1),
     1  fnb_32nd(-1:nx,-1:ny,0:1), fnb_52nd(-1:nx,-1:ny,0:1),
     1  fnb_53(-1:nx,-1:ny,0:1),                                        !ajs 06.04.20
     1  flob(-1:nx,-1:ny,0:1), fnb_fcor(-1:nx,-1:ny,0:1),               !srv 01.10.99
     2  fnb_he(-1:nx,-1:ny,0:1),                                        !srv 21.01.00
     &  fnbPSch(-1:nx,-1:ny,0:1), fnb_mdf(-1:nx,-1:ny,0:1),             !srv 13.10.06
     &  conb(-1:nx,-1:ny,0:1,0:4), ubdia(-1:nx,-1:ny,0:1),              !srv 13.10.06
     &  vbecrb(-1:nx,-1:ny,0:1), vbdia(-1:nx,-1:ny,0:1),                !xpb
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
      integer ncall, n2, ix, iy, ib, area_fix, no_hybr
      integer, save :: iout = 0, iout_b2wdat = 0,                       !srv 16.09.02
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2npc7 = 0, iout_b2npc11 = 0,               !srv 02.04.18
     &                 iout_b2npmo = 0
      integer, save :: mdf_fnb = 0                                      !srv 13.10.06
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0             !srv 06.07.06
      integer, save :: pflux_cor = 1                                    !srv 05.03.16
      integer, save :: ion_vlct_restrict = 0                            !iyv 06.02.14
      integer, save :: vlct_diagno = 1
      character chns*3                                                  !srv 16.09.02
      integer, save :: drift_style                                      !srv 21.06.08
      integer, save :: dia_lim_all_ions = 1                             !som 04.11.21
      real (kind=R8), save :: alpha_ = 0.0_R8                           !srv 09.07.03
      real (kind=R8), save :: xcur = 0.0_R8                             !xpb 27.09.17
      real (kind=R8), save :: ycur = 1.0_R8                             !srv 07.07.08
      real (kind=R8), save :: vis_per = 1.0_R8, vis_q = 1.0_R8          !srv 02.10.17
      real (kind=R8), save :: xvecrb = 1.0_R8, xwdia = 1.0_R8           !srv 27.01.12
      real (kind=R8), save ::
     &   poleldr = 1.0_R8, radeldr = 1.0_R8,                            !srv 25.04.00
     &   poldidr = 1.0_R8, raddidr = 1.0_R8                             !srv 25.04.00
      real (kind=R8) ::
     *   con(-1:nx,-1:ny), flo, con_nodrift, flo_nodrift, xfrhie,
     *   alpha, gamma, vbar, t0, t1, t2, t3, maxfluxlimit(0:1),
     *   flux_limit_min_ti,
     *   flo2dia, con_2dia(-1:nx,-1:ny), flo_mdf, con_mdf(-1:nx,-1:ny), !srv 01.10.99 21.01.00 22.05.03 04.02.05 13.10.06
     *   term(-1:nx,-1:ny), flo_he, flo_53, con_he(-1:nx,-1:ny),        !srv 25.04.00 13.09.02 11.09.09
     *   termhz(-1:nx,-1:ny)                                            !srv 05.11.09
      real (kind=R8) ::
     &   kmprti(-1:nx,-1:ny,0:1), dkmprti(-1:nx,-1:ny,0:1),             !srv 21.10.10
     &   tbfl_fnb(1:nbc), div_fnb, tbfl_fnb_mdf(1:nbc), div_fnb_mdf,    !srv 05.03.16
     &   rdiv, scur, vsprcur                                            !srv 05.03.16 02.04.18
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny), wrk1(-1:nx,-1:ny), wrk2(-1:nx,-1:ny),
     *   wrk3(-1:nx,-1:ny), wrk4(-1:nx,-1:ny), wrk5(-1:nx,-1:ny),       !srv 11.09.09
     *   wrk6(-1:nx,-1:ny), wrk7(-1:nx,-1:ny), wrk8(-1:nx,-1:ny),       !srv 11.09.09 !iys 17.04.12
     *   wpo(-1:nx,-1:ny), wpb(-1:nx,-1:ny)
      real (kind=R8) ::
     *   poc(-1:nx,-1:ny), ub2dia(-1:nx,-1:ny,0:1),                     !srv 01.10.99 !lk 09.06.11
     *   nbti(-1:nx,-1:ny), nbtic(-1:nx,-1:ny),                         !srv 08.10.04
     *   nbzb(-1:nx,-1:ny), nbzbc(-1:nx,-1:ny)                          !srv 13.05.16
      real (kind=R8) :: weight(-1:nx,-1:ny)                             !lk
      integer :: xcount, ycount
      save ncall, alpha, gamma, flux_limit_min_ti, area_fix, no_hybr
!!FR Local save variables that get updated need to be THREADPRIVATE
!!FR Local save variables that are read only are okay - no change needed.
!!FR Because alpha, gamma, flux_limit_min_ti, area_fix & others could
!!FR potentially be written to by multiple threads on their ncall==0
!!FR iteration we need to make them THREADPRIVATE too.
!!TF Note that variables that appear in DATA statement have implied SAVE attribute.
      !$OMP THREADPRIVATE(ncall)
      !$OMP THREADPRIVATE(alpha, gamma)
      !$OMP THREADPRIVATE(flux_limit_min_ti, area_fix, no_hybr)
      !$OMP THREADPRIVATE(iout, iout_b2wdat, iout_b2npc9, iout_b2npco)
      !$OMP THREADPRIVATE(iout_b2npmo, mdf_fnb)
      !$OMP THREADPRIVATE(no_Ptncr_x_co, no_Ptncr_y_co)
      !$OMP THREADPRIVATE(pflux_cor, ion_vlct_restrict)
      !$OMP THREADPRIVATE(vlct_diagno, drift_style, alpha_, xcur, ycur)
      !$OMP THREADPRIVATE(vis_per, vis_q, xvecrb, xwdia)
      !$OMP THREADPRIVATE(poldidr, raddidr, poleldr, radeldr)
*   ..procedures
      real (kind=R8) ::
     *   cvmgp
      external xertst, ipgeti, ipgetr
      external b2xvsg, b2xvff, intfaceh, intfacev, cvmgp
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
        call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)                 !srv 06.07.06
        call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)                 !srv 06.07.06
        call ipgetr ('b2tfnb_poleldr', poleldr)                          !srv 25.04.00 {
        call ipgetr ('b2tfnb_radeldr', radeldr)
        call ipgetr ('b2tfnb_poldidr', poldidr)
        call ipgetr ('b2tfnb_raddidr', raddidr)                          !srv 25.04.00 }
        call ipgeti ('b2tfnb_iout', iout)                                !srv 16.09.02
        call ipgetr ('b2tlc0_alpha', alpha_)                             !srv 09.07.03 {
        call xertst (0.0_R8.le.alpha_,
     &   'faulty internal parameter alpha')
        call xertst (.not.(0.0_R8.ne.alpha .and. 0.0_R8.ne.alpha_),      !srv 14.10.09
     &   'b2tfnb_alpha and b2tlc0_alpha cannot be both different '//
     &   'from zero ')                                                   !srv 09.07.03 }
        call ipgeti ('b2tfnb_mdf', mdf_fnb)                              !srv 13.10.06
        call ipgeti ('b2tfnb_pflux_cor', pflux_cor)                      !srv 05.03.16
        if (redef_gmtry.eq.0) then
          drift_style = 1
        else
          drift_style = 2
        end if
        call ipgeti ('b2tfnb_drift_style', drift_style)                  !srv 21.06.08
        call ipgetr ('b2tfnb_xcur', xcur)                                !xpb 27.09.17
        call ipgetr ('b2tfnb_ycur', ycur)                                !srv 17.07.08
        call ipgetr ('b2tfnb_vis_per', vis_per)                          !srv 02.10.17
        call ipgetr ('b2tfnb_vis_q', vis_q)                              !srv 02.10.17
        call ipgeti ('b2wdat_iout', iout_b2wdat)                         !iys 17.04.17
        call ipgeti ('b2npc9_iout', iout_b2npc9)                         !srv 11.09.09
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npc7_iout', iout_b2npc7)                         !srv 02.04.18
        call ipgeti ('b2npc11_iout', iout_b2npc11)                       !srv 02.04.18
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgetr ('b2tfrn_xvecrb', xvecrb)                            !srv 27.01.12
        call ipgetr ('b2tfrn_xwdia', xwdia)                              !srv 27.01.12
        call ipgeti ('b2npmo_ion_vlct_restrict', ion_vlct_restrict)      !iyv 06.02.14
        call ipgeti ('b2npmo_vlct_diagno', vlct_diagno)
        call ipgeti ('b2tfnb_dia_lim_all_ions', dia_lim_all_ions)        !som 04.11.21
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

      write(chns,'(i3.3)') isb
      if(mdf_fnb.eq.0) then                                              !srv 21.10.10 {
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
     &      (wrk3(ix,iy)-wrk3(bottomix(ix,iy),bottomiy(ix,iy)))
!     &    /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))*2.0_R8 because cddi contains factor vol/hx/hy1
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
           wrk4(ix,iy) = (wrk3(ix,iy)-wrk3(leftix(ix,iy),leftiy(ix,iy)))
!     &    /(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*2.0_R8  because cddi contains factor vol/hx/hy1
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
      endif                                                              !srv 21.10.10 }

c Drift terms
      maxfluxlimit=0.0_R8
      ubdia =0.0_R8
      ub2dia=0.0_R8                                                      !srv 01.10.99
      vbdia =0.0_R8
      wbdia =0.0_R8
      vbecrb=0.0_R8
      fnbPSch = 0.0_R8                                                   !srv 13.10.06
      fnb_mdf = 0.0_R8                                                   !srv 19.11.09
      conb  = 0.0e0_R8                                                   !srv 13.10.06
      nbti=nb*ti
      if((maxval(facdrift).ne.0.0_R8 .or. maxval(fac_ExB).ne.0.0_R8)
     &    .and. minval(zb).gt.0.0_R8) then

       if(drift_style.eq.0) then                                         !srv 21.06.08

        call intfacev(nx,ny,bottomix,bottomiy,vol,po,wpo)
        call intfacev(nx,ny,bottomix,bottomiy,vol,nbti,wpb)
        do iy = -1, ny
          do ix = -1, nx                                                 !xpb
cxpb  Must check that we are not at the island X-point or at the mesh edge
            if(topiy(topix(ix,iy),min(topiy(ix,iy),ny)).ne.iy .and.
     &         topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &         rightix(ix,iy).ne.nx+1 .and. leftiy(ix,iy).ne.-2 .and.
     &         region(ix,iy,0).ne.0) then                                !xpb
cxpb  Standard case
cwdk correct computation of drift components to always represent the poloidal
cwdk projection of the drift
              vbecrb(ix,iy,0) = - fac_ExB(ix,iy)*poleldr/hy1(ix,iy)*     !xpb
     1         bb(ix,iy,2)/bb(ix,iy,3)**2*                               !xpb
     2         (wpo(topix(ix,iy),topiy(ix,iy))-wpo(ix,iy))
              vbdia(ix,iy,0) = - bb(ix,iy,2)/bb(ix,iy,3)/hy1(ix,iy)/     !xpb
     1         bb(ix,iy,3)*2.0_R8*ti(ix,iy)/bb(ix,iy,3)/qe/zb(ix,iy)*    !xpb
     3         (wbbv(topix(ix,iy),topiy(ix,iy),3)-wbbv(ix,iy,3))*
     &         facdrift(ix,iy)*poldidr                                   !xpb
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
              wbdia(ix,iy,0) = 0.0_R8
            endif                                                        !xpb
            ubdia(ix,iy,0) = vbecrb(ix,iy,0) + vbdia(ix,iy,0)            !srv 17.05.99
            ub2dia(ix,iy,0) = vbecrb(ix,iy,0) +
     &        2.0_R8*vbdia(ix,iy,0)                                      !srv 01.10.99
          enddo
        enddo

        call intfaceh(nx,ny,leftix,leftiy,vol,po,wpo)
        call intfaceh(nx,ny,leftix,leftiy,vol,nbti,wpb)
        do iy = -1, ny                                                   !xpb
          do ix = -1, nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
             vbecrb(ix,iy,1) = fac_ExB(ix,iy)*radeldr/hx(ix,iy)*         !xpb
     1          bb(ix,iy,2)/bb(ix,iy,3)**2*                              !xpb
     2         (wpo(rightix(ix,iy),rightiy(ix,iy))-wpo(ix,iy))
             vbdia(ix,iy,1) = facdrift(ix,iy)*raddidr*
     &         ti(ix,iy)*bb(ix,iy,2)/(qe*zb(ix,iy)*hx(ix,iy))            !xpb
     &         *(2.0_R8/bb(ix,iy,3)**3)                                  !xpb
     &         *(wbbr(ix,iy,3)-wbbl(ix,iy,3))                            !xpb
             wbdia(ix,iy,1) = bb(ix,iy,2)
     &         /(qe*zb(ix,iy)*hx(ix,iy)*nb(ix,iy))                       !xpb
     &         /bb(ix,iy,3)**2*                                          !xpb
     2         (wpb(rightix(ix,iy),rightiy(ix,iy))-wpb(ix,iy))*
     &         facdrift(ix,iy)
             ubdia(ix,iy,1) = vbecrb(ix,iy,1) + vbdia(ix,iy,1)           !srv 17.05.99
             ub2dia(ix,iy,1) = vbecrb(ix,iy,1) +
     &         2.0_R8*vbdia(ix,iy,1)                                     !srv 01.10.99
            endif
          enddo
        enddo

* drift velocities at walls are extrapolated from those at neighbouring cells
* drift velocities at iy=ny are extrapolated from those at iy=ny-1       !xpb
* drift velocities at iy=-1 are extrapolated from those at iy=0          !xpb
        do ix = -1, nx                                                   !xpb
          do iy = -1, ny                                                 !xpb
            if(leftix(ix,iy).eq.-2) then
              vbecrb(ix,iy,:) =
     &          vbecrb(rightix(ix,iy),rightiy(ix,iy),:)                  !xpb
              vbdia(ix,iy,:) =
     &          vbdia(rightix(ix,iy),rightiy(ix,iy),:)                   !xpb
              wbdia(ix,iy,:) =
     &          wbdia(rightix(ix,iy),rightiy(ix,iy),:)                   !xpb
              ubdia(ix,iy,:) = vbecrb(ix,iy,:) + vbdia(ix,iy,:)
              ub2dia(ix,iy,:) =
     &          vbecrb(ix,iy,:) + 2.0_R8*vbdia(ix,iy,:)                  !srv 17.07.08
            elseif(rightix(ix,iy).eq.nx+1) then                          !xpb
              vbecrb(ix,iy,:) =
     &          vbecrb(leftix(ix,iy),leftiy(ix,iy),:)                    !xpb
              vbdia(ix,iy,:) = vbdia(leftix(ix,iy),leftiy(ix,iy),:)      !xpb
              wbdia(ix,iy,:) = wbdia(leftix(ix,iy),leftiy(ix,iy),:)      !xpb
              ubdia(ix,iy,:) = vbecrb(ix,iy,:) + vbdia(ix,iy,:)
              ub2dia(ix,iy,:) =
     &          vbecrb(ix,iy,:) + 2.0_R8*vbdia(ix,iy,:)                  !srv 17.07.08
            endif
            if(topiy(ix,iy).eq.ny+1) then                                !xpb
              vbecrb(ix,iy,:) =
     .          vbecrb(bottomix(ix,iy),bottomiy(ix,iy),:)                !xpb
              vbdia(ix,iy,:) =
     &          vbdia(bottomix(ix,iy),bottomiy(ix,iy),:)                 !xpb
              wbdia(ix,iy,:) =
     &          wbdia(bottomix(ix,iy),bottomiy(ix,iy),:)                 !xpb
              ubdia(ix,iy,:) = vbecrb(ix,iy,:) + vbdia(ix,iy,:)
              ub2dia(ix,iy,:) =
     &          vbecrb(ix,iy,:) + 2.0_R8*vbdia(ix,iy,:)                  !srv 17.07.08
            elseif(bottomiy(ix,iy).eq.-2) then
              vbecrb(ix,iy,:) = vbecrb(topix(ix,iy),topiy(ix,iy),:)      !xpb
              vbdia(ix,iy,:) = vbdia(topix(ix,iy),topiy(ix,iy),:)        !xpb
              wbdia(ix,iy,:) = wbdia(topix(ix,iy),topiy(ix,iy),:)        !xpb
              ubdia(ix,iy,:) = vbecrb(ix,iy,:) + vbdia(ix,iy,:)
              ub2dia(ix,iy,:) =
     &          vbecrb(ix,iy,:) + 2.0_R8*vbdia(ix,iy,:)                  !srv 17.07.08
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

       else if (drift_style .eq. 1 .or. redef_gmtry .eq. 0) then        !srv 21.06.08 {

        weight = 1.0e0_R8                                                !lk
!lk       call intcor(nx,ny,vol,po,poc)
        call intcor(nx,ny,weight,po,poc)                                 !lk 12.12.07
!lk       call intcor(nx,ny,vol,nbti,nbtic)                              !srv 08.10.04
        call intcor(nx,ny,weight,nbti,nbtic)                             !lk 12.12.07
        nbzb = nb*zb                                                     !srv 13.05.16
        call intcor(nx,ny,weight,nbzb,nbzbc)                             !srv 13.05.16
        xcount = 0
        do iy = 0, ny-1
         do ix = 0, nx
          if(leftix(ix,iy).ne.-2) then                                   !srv 26.11.02
             vbecrb(ix,iy,0) = - poleldr*fac_ExB(ix,iy)/                 !srv 17.05.99 25.08.99 07.06.00 05.02.05 13.05.16 {
     1      (hy1(ix,iy)+hy1(leftix(ix,iy),leftiy(ix,iy)))* ! signmf*     !srv
     2      (bb(ix,iy,2)/bb(ix,iy,3)**2                                  !srv 05.02.05
     2         +bb(leftix(ix,iy),leftiy(ix,iy),2)/
     2          bb(leftix(ix,iy),leftiy(ix,iy),3)**2)*                   !srv 05.02.05
     2       (poc(topix(ix,iy),topiy(ix,iy))-poc(ix,iy))                 !srv 13.05.16 }
!           vbecrb(ix,iy,0) = - poleldr*fac_ExB(ix,iy)* ! signmf*         !lk 12.12.07 } !srv 11.01.13
!     1       (hz(ix,iy)*bb(ix,iy,2)
!     2         +hz(leftix(ix,iy),leftiy(ix,iy))*
!     3         bb(leftix(ix,iy),leftiy(ix,iy),2))/2.0e0_R8*
!     4       (vol(leftix(ix,iy),leftiy(ix,iy))/bb(ix,iy,3)**2
!     5         +vol(ix,iy)/bb(leftix(ix,iy),leftiy(ix,iy),3)**2)/
!     6         (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))*
!     7         (poc(topix(ix,iy),topiy(ix,iy))-poc(ix,iy))/
!     8         (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))*
!     9         (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))               !lk 12.12.07)
           if(ion_vlct_restrict.eq.2) then                               !iyv 06.02.14 {
             if( abs(vbecrb(ix,iy,0)) .gt. 2.0_R8*cssb(ix,iy,0) ) then
               xcount = xcount+1
               if (vlct_diagno.ge.2) then
                 write(6,'(a,a,i3,a,i4,a,i4,a)')
     .            'Poloidal ExB velocity restriction applied ',
     .            'to species ',isb,' in cell (',ix,',',iy,'):'
                 write(6,'(a,1p,1e12.5,a,1e12.5)')
     .            'Mach 2: ',2.0_R8*cssb(ix,iy,0),
     .            ' Vx_ExB:',vbecrb(ix,iy,0)
               end if
               vbecrb(ix,iy,0)=
     .          sign(2.0_R8*cssb(ix,iy,0),vbecrb(ix,iy,0))
             endif
           endif                                                         !iyv 06.02.14  }
           vbdia(ix,iy,0) =    poldidr*facdrift(ix,iy)* ! signmf*        !srv 07.06.00 17.05.99 25.04.00 05.02.05 11.01.13
     &          (ti(leftix(ix,iy),leftiy(ix,iy))                         !srv 05.03.99
     &           *bb(leftix(ix,iy),leftiy(ix,iy),2)                      !srv 05.02.05
     &           /(qe*zb(ix,iy)*hy1(leftix(ix,iy),leftiy(ix,iy)))        !srv
     &          +ti(ix,iy)*bb(ix,iy,2)/(qe*zb(ix,iy)*hy1(ix,iy)))*0.5_R8 !srv 05.02.05
     &         *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))    !srv 05.03.99
           wbdia(ix,iy,0) = - facdrift(ix,iy)/                           !srv 13.05.16 {
     1        (hy1(ix,iy)+hy1(leftix(ix,iy),leftiy(ix,iy)))*
     2        (bb(ix,iy,2)/bb(ix,iy,3)**2
     2        +bb(leftix(ix,iy),leftiy(ix,iy),2)/
     2         bb(leftix(ix,iy),leftiy(ix,iy),3)**2)*
     2         (nbtic(topix(ix,iy),topiy(ix,iy))-nbtic(ix,iy))/
     &        (nbzbc(topix(ix,iy),topiy(ix,iy))+nbzbc(ix,iy))*2.0_R8/qe  !srv 13.05.16 }
           ubdia(ix,iy,0) = vbecrb(ix,iy,0) + vbdia(ix,iy,0)             !srv 17.05.99
           ub2dia(ix,iy,0) = vbecrb(ix,iy,0) + 2.0e0_R8*vbdia(ix,iy,0)   !srv 01.10.99
          endif                                                          !srv 26.11.02
         enddo
        enddo
        if (vlct_diagno.ge.1) then
          if (xcount.gt.1) then
            write(6,'(a,a,i6,a,i3)') 'In b2tfnb:',
     .       ' Poloidal ExB velocity limit of Mach 2 was applied ',
     .         xcount,' times for species ',isb
          else if (xcount.eq.1) then
            write(6,'(a,a,i6,a,i3)') 'In b2tfnb:',
     .       ' Poloidal ExB velocity limit of Mach 2 was applied ',
     .         xcount,' time for species ',isb
          end if
        end if

        ycount = 0
        do iy = 0, ny
         do ix = 0, nx-1
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then      !srv 26.11.02
             vbecrb(ix,iy,1) = radeldr*fac_ExB(ix,iy)/                   !srv 17.05.99 25.08.99 07.06.00 05.02.05 13.05.16 {
     1       (hx(ix,iy)+hx(bottomix(ix,iy),bottomiy(ix,iy)))* !signmf*   !srv
     2       (bb(ix,iy,2)/bb(ix,iy,3)**2                                 !srv 05.02.05
     2          +bb(bottomix(ix,iy),bottomiy(ix,iy),2)/
     2           bb(bottomix(ix,iy),bottomiy(ix,iy),3)**2)*              !srv 05.02.05
     2      (poc(rightix(ix,iy),rightiy(ix,iy))-poc(ix,iy))              !srv 13.05.16 }
!           vbecrb(ix,iy,1) = radeldr*fac_ExB(ix,iy)* ! signmf*          !srv 11.01.13 07.06.00 21.04.00 25.08.99 17.05.99
!     1       (hz(ix,iy)*bb(ix,iy,2)                                     !lk( 12.12.07
!     2          +hz(bottomix(ix,iy),bottomiy(ix,iy))*
!     3          bb(bottomix(ix,iy),bottomiy(ix,iy),2))/2.0e0_R8*
!     4       (vol(bottomix(ix,iy),bottomiy(ix,iy))/bb(ix,iy,3)**2
!     5          +vol(ix,iy)/bb(bottomix(ix,iy),bottomiy(ix,iy),3)**2)/
!     6          (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))*
!     7        (poc(rightix(ix,iy),rightiy(ix,iy))-poc(ix,iy))/          !lk 12.12.07)
!     8          (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))*
!     9          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
           if(ion_vlct_restrict.eq.2) then                               !iyv 06.02.14 {
             if( abs(vbecrb(ix,iy,1)) .gt. 2.0_R8*cssb(ix,iy,1) ) then
               ycount = ycount + 1
               if (vlct_diagno.ge.2) then
                 write(6,'(a,a,i3,a,i4,a,i4,a)')
     .            'Radial ExB velocity restriction applied ',
     .            'to species ',isb,' in cell (',ix,',',iy,'):'
                 write(6,'(a,1p,1e12.5,a,1e12.5)')
     .            'Mach 2: ',2.0_R8*cssb(ix,iy,1),
     .            ' Vy_ExB:',vbecrb(ix,iy,1)
               end if
               vbecrb(ix,iy,1)=
     .          sign(2.0_R8*cssb(ix,iy,1),vbecrb(ix,iy,1))
             endif
           endif                                                         !iyv 06.02.14  }
           vbdia(ix,iy,1) =  -raddidr*facdrift(ix,iy)* ! signmf*         !srv 11.01.13 07.06.00 17.05.99 25.04.00
     &           (ti(bottomix(ix,iy),bottomiy(ix,iy))                    !srv 05.03.99
     &           *bb(bottomix(ix,iy),bottomiy(ix,iy),2)                  !srv
     &           /(qe*zb(ix,iy)*hx(bottomix(ix,iy),bottomiy(ix,iy)))     !srv
     &           +ti(ix,iy)*bb(ix,iy,2)/(qe*zb(ix,iy)*hx(ix,iy)))*0.5_R8 !srv
     &        *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy)) !srv 05.03.99
           wbdia(ix,iy,1) = facdrift(ix,iy)/                             !srv 13.05.16 {
     1      (hx(ix,iy)+hx(bottomix(ix,iy),bottomiy(ix,iy)))*
     2      (bb(ix,iy,2)/bb(ix,iy,3)**2
     2       +bb(bottomix(ix,iy),bottomiy(ix,iy),2)/
     2       bb(bottomix(ix,iy),bottomiy(ix,iy),3)**2)*
     2     (nbtic(rightix(ix,iy),rightiy(ix,iy))-nbtic(ix,iy))/
     &     (nbzbc(rightix(ix,iy),rightiy(ix,iy))+nbzbc(ix,iy))*2.0_R8/qe !srv 13.05.16 }

           ubdia(ix,iy,1) = vbecrb(ix,iy,1) + vbdia(ix,iy,1)             !srv 17.05.99
           ub2dia(ix,iy,1) =
     &       vbecrb(ix,iy,1) + 2.0e0_R8*vbdia(ix,iy,1)                   !srv 01.10.99
          endif                                                          !srv 26.11.02
         enddo
        enddo

        if (vlct_diagno.ge.1) then
          if (ycount.gt.1) then
            write(6,'(a,a,i6,a,i3)') 'In b2tfnb:',
     .       ' Radial ExB velocity limit of Mach 2 was applied ',
     .         ycount,' times for species ',isb
          else if (ycount.eq.1) then
            write(6,'(a,a,i6,a,i3)') 'In b2tfnb:',
     .       ' Radial ExB velocity limit of Mach 2 was applied ',
     .         ycount,' time for species ',isb
          end if
        end if

       else if (drift_style .eq. 2 .and. redef_gmtry .eq. 1) then        ! IYS 04.03.2019 - 11.02.2020

        weight = 1.0e0_R8                                                !lk

        call intcor(nx,ny,weight,po,poc)                                 !lk 12.12.07
        call intcor(nx,ny,weight,nbti,nbtic)                             !lk 12.12.07
        nbzb = nb*zb                                                     !srv 13.05.16
        call intcor(nx,ny,weight,nbzb,nbzbc)                             !srv 13.05.16
        xcount = 0
        do iy = 0, ny-1
         do ix = 0, nx
          if(leftix(ix,iy).ne.-2) then                                   !srv 26.11.02
             vbecrb(ix,iy,0) = - poleldr*fac_ExB(ix,iy)                  !srv 17.05.99 25.08.99 07.06.00 05.02.05 13.05.16 {
     2      *(poc(topix(ix,iy),topiy(ix,iy))-poc(ix,iy))                 !srv 13.05.16 }
     &      *wbbl(ix,iy,2)/(wbbl(ix,iy,3))**2/hyl(ix,iy)/qc(ix,iy)       ! IYS 04.03.2019

           if(ion_vlct_restrict.eq.2) then                               !iyv 06.02.14 {
             if( abs(vbecrb(ix,iy,0)) .gt. 2.0_R8*cssb(ix,iy,0) ) then
               xcount = xcount+1
               if (vlct_diagno.ge.2) then
                 write(6,'(a,a,i3,a,i4,a,i4,a)')
     .            'Poloidal ExB velocity restriction applied ',
     .            'to species ',isb,' in cell (',ix,',',iy,'):'
                 write(6,'(a,1p,1e12.5,a,1e12.5)')
     .            'Mach 2: ',2.0_R8*cssb(ix,iy,0),
     .            ' Vx_ExB:',vbecrb(ix,iy,0)
               end if
               vbecrb(ix,iy,0)=
     .          sign(2.0_R8*cssb(ix,iy,0),vbecrb(ix,iy,0))
             endif
           endif                                                         !iyv 06.02.14  }

           vbdia(ix,iy,0) =    poldidr*facdrift(ix,iy)* ! signmf*        !srv 07.06.00 17.05.99 25.04.00 05.02.05 11.01.13
     &       0.5_R8 * (ti(leftix(ix,iy),leftiy(ix,iy)) + ti(ix,iy)) *
     &         wbbl(ix,iy,2) / (qe*zb(ix,iy) * hyl(ix,iy) * qc(ix,iy))   ! IYS 04.03.2019
     &         *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))    !srv 05.03.99
           wbdia(ix,iy,0) = - facdrift(ix,iy)/                           ! IYS 04.03.2019 {
     1        (qe*zb(ix,iy) * hyl(ix,iy) * qc(ix,iy))  !(hy1(ix,iy)+hy1(leftix(ix,iy),leftiy(ix,iy)))*
!     2        (bb(ix,iy,2)/bb(ix,iy,3)**2
!     2        +bb(leftix(ix,iy),leftiy(ix,iy),2)/
!    2         bb(leftix(ix,iy),leftiy(ix,iy),3)**2)*
     &         * wbbl(ix,iy,2) / (wbbl(ix,iy,3))**2 *
     2         (nbtic(topix(ix,iy),topiy(ix,iy)) - nbtic(ix,iy)) /
     &         (nb(ix,iy) + nb(leftix(ix,iy),leftiy(ix,iy))) * 2.0_R8
!     &        (nbzbc(topix(ix,iy),topiy(ix,iy))+nbzbc(ix,iy))*2.0_R8/qe !IYS 04.03.2019 }
           ubdia(ix,iy,0) = vbecrb(ix,iy,0) + vbdia(ix,iy,0)             !srv 17.05.99
           ub2dia(ix,iy,0) = vbecrb(ix,iy,0) + 2.0e0_R8*vbdia(ix,iy,0)   !srv 01.10.99
          endif                                                          !srv 26.11.02
         enddo
        enddo

        if (vlct_diagno.ge.1) then
          if (xcount.gt.1) then
            write(6,'(a,a,i6,a,i3)') 'In b2tfnb:',
     .       ' Poloidal ExB velocity limit of Mach 2 was applied ',
     .         xcount,' times for species ',isb
          else if (xcount.eq.1) then
            write(6,'(a,a,i6,a,i3)') 'In b2tfnb:',
     .       ' Poloidal ExB velocity limit of Mach 2 was applied ',
     .         xcount,' time for species ',isb
          end if
        end if

        ycount = 0
        do iy = 0, ny
         do ix = 0, nx-1
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then      !srv 26.11.02
             vbecrb(ix,iy,1) = radeldr*fac_ExB(ix,iy)                    !srv 17.05.99 25.08.99 07.06.00 05.02.05 13.05.16 {
     2      *(poc(rightix(ix,iy),rightiy(ix,iy))-poc(ix,iy))             !srv 13.05.16 }
     &      * wbbv(ix,iy,2) / wbbv(ix,iy,3)**2 / hxb(ix,iy)              ! IYS 04.03.2019

           if(ion_vlct_restrict.eq.2) then                               !iyv 06.02.14 {
             if( abs(vbecrb(ix,iy,1)) .gt. 2.0_R8*cssb(ix,iy,1) ) then
               ycount = ycount + 1
               if (vlct_diagno.ge.2) then
                 write(6,'(a,a,i3,a,i4,a,i4,a)')
     .            'Radial ExB velocity restriction applied ',
     .            'to species ',isb,' in cell (',ix,',',iy,'):'
                 write(6,'(a,1p,1e12.5,a,1e12.5)')
     .            'Mach 2: ',2.0_R8*cssb(ix,iy,1),
     .            ' Vy_ExB:',vbecrb(ix,iy,1)
               end if
               vbecrb(ix,iy,1)=
     .          sign(2.0_R8*cssb(ix,iy,1),vbecrb(ix,iy,1))
             endif
           endif                                                         !iyv 06.02.14  }

           vbdia(ix,iy,1) =  -raddidr*facdrift(ix,iy)* ! signmf*         !srv 11.01.13 07.06.00 17.05.99 25.04.00
     &     0.5_R8*(ti(bottomix(ix,iy),bottomiy(ix,iy)) + ti(ix,iy))*
     &       wbbv(ix,iy,2) / (qe*zb(ix,iy)*hxb(ix,iy))
     &        *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy)) !srv 05.03.99
           wbdia(ix,iy,1) = facdrift(ix,iy)/                             !IYS 04.03.2019 {
     1      (qe*zb(ix,iy)*hxb(ix,iy))   !(hx(ix,iy)+hx(bottomix(ix,iy),bottomiy(ix,iy)))*
     2      * wbbv(ix,iy,2) / wbbv(ix,iy,3)**2 *
!     2       +bb(bottomix(ix,iy),bottomiy(ix,iy),2)/
!     2       bb(bottomix(ix,iy),bottomiy(ix,iy),3)**2)*
     2     (nbtic(rightix(ix,iy),rightiy(ix,iy)) - nbtic(ix,iy)) /
     &     (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))*2.0_R8 !/qe   !IYS 04.03.2019 }

           ubdia(ix,iy,1) = vbecrb(ix,iy,1) + vbdia(ix,iy,1)             !srv 17.05.99
           ub2dia(ix,iy,1) = vbecrb(ix,iy,1) + 2.0e0_R8*vbdia(ix,iy,1)   !srv 01.10.99
          endif                                                          !srv 26.11.02
         enddo
        enddo

        if (vlct_diagno.ge.1) then
          if (ycount.gt.1) then
            write(6,'(a,a,i6,a,i3)') 'In b2tfnb:',
     .       ' Radial ExB velocity limit of Mach 2 was applied ',
     .         ycount,' times for species ',isb
          else if (ycount.eq.1) then
            write(6,'(a,a,i6,a,i3)') 'In b2tfnb:',
     .       ' Radial ExB velocity limit of Mach 2 was applied ',
     .         ycount,' time for species ',isb
          end if
        end if

       else           ! drift_style                                      ! IYS 11.02.2020 {

         write(*,'(2(a20,i4))') 'drift_style = ',drift_style,
     &      'redef_gmtry = ',redef_gmtry
         call xerrab (
     &   'faulty definition of drift_style and redef_gmtry, ' //
     &   'or logical mistake occurred')

       endif          ! drift_style                                      ! IYS 11.02.2020 }

       call b2txnpsch (nx, ny, isb, guard_flows, facdrift,               !srv 02.04.18 {
     &        kmprti, cddb, conb)

*   ..compute Pfirsch-Schlueter particle fluxes                          !srv 13.10.06
       call b2tnpsch (nx, ny, isb, guard_flows, facdrift, zb,            !srv 02.04.18
     &   nb, ti, cddb, fnbPSch)

      endif

      if (dia_lim_all_ions .eq. 1) then                                  !som 04.11.21
       if (.not.is_neutral(isb)) then                                    !som 04.11.21 {
        do iy = 0, ny
          do ix = -1, nx
            if (leftix(ix,iy).ne.-2) then
              t0 = max(ti(ix,iy), ti(leftix(ix,iy),leftiy(ix,iy)))
              t1 = (nb(ix,iy) + nb(leftix(ix,iy),leftiy(ix,iy)))
     &          /2.0_R8
              vbar = sqrt(8.0_R8 * t0 / (pi * am(isb)*mp))               !som 04.11.21
              t2 = 1.0_R8 + abs(fnbPSch(ix,iy,0))/(
     &         vbar*t1*(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     &         (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))
     &          )
              fnbPSch(ix,iy,0)= fnbPSch(ix,iy,0)/t2
              vbdia(ix,iy,0)= vbdia(ix,iy,0)/
     &                          (1.0_R8 + abs(vbdia(ix,iy,0))/vbar)
              wbdia(ix,iy,0)= wbdia(ix,iy,0)/
     &                          (1.0_R8 + abs(wbdia(ix,iy,0))/vbar)
              ubdia(ix,iy,0)= ubdia(ix,iy,0)/
     &                          (1.0_R8 + abs(ubdia(ix,iy,0))/vbar)
              ub2dia(ix,iy,0)= ub2dia(ix,iy,0)/
     &                          (1.0_R8 + abs(ub2dia(ix,iy,0))/vbar)
            endif
            if (bottomiy(ix,iy).ne.-2) then
              t0 = max(ti(ix,iy), ti(bottomix(ix,iy),bottomiy(ix,iy)))
              t1 = (nb(ix,iy) + nb(bottomix(ix,iy),bottomiy(ix,iy)))
     &         /2.0_R8
              vbar = sqrt(8.0_R8 * t0 / (pi * am(isb)*mp))               !som 04.11.21
              t2 = 1.0_R8 + abs(fnbPSch(ix,iy,1))/(
     &         vbar*t1*(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))
     &         /(hy(ix,iy)+hy(bottomix(ix,iy),bottomiy(ix,iy)))
     &         )
              fnbPSch(ix,iy,1)= fnbPSch(ix,iy,1)/t2
              vbdia(ix,iy,1)= vbdia(ix,iy,1)/
     &         (1.0_R8 + abs(vbdia(ix,iy,1))/vbar)
              wbdia(ix,iy,1)= wbdia(ix,iy,1)/
     &         (1.0_R8 + abs(wbdia(ix,iy,1))/vbar)
              ubdia(ix,iy,1)= ubdia(ix,iy,1)/
     &         (1.0_R8 + abs(ubdia(ix,iy,1))/vbar)
              ub2dia(ix,iy,1)= ub2dia(ix,iy,1)/
     &         (1.0_R8 + abs(ub2dia(ix,iy,1))/vbar)
            endif
          enddo
        enddo
       endif                                                             !som 04.11.21 }
      else
       if (isb.ne.ismain .and. .not.is_neutral(isb)) then                !lk 02.06.11 {
        do iy = 0, ny
          do ix = -1, nx
            if (leftix(ix,iy).ne.-2) then
              t0 = max(ti(ix,iy), ti(leftix(ix,iy),leftiy(ix,iy)))
              t1 = (nb(ix,iy) + nb(leftix(ix,iy),leftiy(ix,iy)))
     &          /2.0_R8
              vbar =
     .         abs(bb(ix,iy,0) + bb(leftix(ix,iy),leftiy(ix,iy),0))      !srv 11.01.13
     &           /(bb(ix,iy,3) + bb(leftix(ix,iy),leftiy(ix,iy),3))
     &          *sqrt(8.0_R8 * t0 / (pi * am(isb)*mp))
              t2 = 1.0_R8 + abs(fnbPSch(ix,iy,0))/(
     &         vbar*t1*(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     &         (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))
     &          )
              fnbPSch(ix,iy,0)= fnbPSch(ix,iy,0)/t2
              vbdia(ix,iy,0)= vbdia(ix,iy,0)/
     &                            (1.0_R8 + abs(vbdia(ix,iy,0))/vbar)
              wbdia(ix,iy,0)= wbdia(ix,iy,0)/
     &                            (1.0_R8 + abs(wbdia(ix,iy,0))/vbar)
              ubdia(ix,iy,0)= ubdia(ix,iy,0)/
     &                            (1.0_R8 + abs(ubdia(ix,iy,0))/vbar)
              ub2dia(ix,iy,0)= ub2dia(ix,iy,0)/
     &                            (1.0_R8 + abs(ub2dia(ix,iy,0))/vbar)
            endif
            if (bottomiy(ix,iy).ne.-2) then
              t0 = max(ti(ix,iy), ti(bottomix(ix,iy),bottomiy(ix,iy)))
              t1 = (nb(ix,iy) + nb(bottomix(ix,iy),bottomiy(ix,iy)))
     &         /2.0_R8
              vbar=
     .         abs(bb(ix,iy,0) + bb(bottomix(ix,iy),bottomiy(ix,iy),0))  !srv 11.01.13
     &         /(bb(ix,iy,3) + bb(bottomix(ix,iy),bottomiy(ix,iy),3))
     &         *sqrt(8.0_R8 * t0 / (pi * am(isb)*mp))
              t2 = 1.0_R8 + abs(fnbPSch(ix,iy,1))/(
     &         vbar*t1*(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))
     &         /(hy(ix,iy)+hy(bottomix(ix,iy),bottomiy(ix,iy)))
     &         )
              fnbPSch(ix,iy,1)= fnbPSch(ix,iy,1)/t2
              vbdia(ix,iy,1)= vbdia(ix,iy,1)/
     &         (1.0_R8 + abs(vbdia(ix,iy,1))/vbar)
              wbdia(ix,iy,1)= wbdia(ix,iy,1)/
     &         (1.0_R8 + abs(wbdia(ix,iy,1))/vbar)
              ubdia(ix,iy,1)= ubdia(ix,iy,1)/
     &         (1.0_R8 + abs(ubdia(ix,iy,1))/vbar)
              ub2dia(ix,iy,1)= ub2dia(ix,iy,1)/
     &         (1.0_R8 + abs(ub2dia(ix,iy,1))/vbar)
            endif
          enddo
        enddo
       endif                                                            !lk 02.06.11 }
      endif

* ..compute fluxes
*   ..compute wrk0
      call intfaceh(nx,ny,leftix,leftiy,vol,ub,wrk0)
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
         fnb_53(ix,iy,0)      = 0.0_R8                                  !ajs 06.04.20
         flob(ix,iy,0)        = 0.0_R8
         fnb_fcor(ix,iy,0)    = 0.0_R8                                  !srv 01.10.99
         fnb_he(ix,iy,0)      = 0.0_R8                                  !srv 21.01.00
         fnb_mdf(ix,iy,0)     = 0.0_R8                                  !srv 13.10.06
         conb(ix,iy,0,0:2)    = 0.0_R8                                  !srv 13.10.06
         term(ix,iy)          = 0.0_R8
         termhz(ix,iy)        = 0.0_R8                                  !srv 05.11.09
         con(ix,iy)           = 0.0_R8
         con_2dia(ix,iy)      = 0.0_R8
         con_he(ix,iy)        = 0.0_R8
         con_mdf(ix,iy)       = 0.0_R8
         fllim0fnb(ix,iy,0)   = 1.0_R8
         if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
          fna_pinch(ix,iy,0,isb)    = 0.0_R8
          fna_pll(ix,iy,0,isb)      = 0.0_R8
          fna_drift(ix,iy,0,isb)    = 0.0_R8
          fna_nanom(ix,iy,0,isb)    = 0.0_R8
          fna_panom(ix,iy,0,isb)    = 0.0_R8
          fna_ch(ix,iy,0,isb)       = 0.0_R8
          fna_pschused(ix,iy,0,isb) = 0.0_R8
         endif
        else
         term(ix,iy) = xfrhie*wrk2(ix,iy)*0.5_R8*abs(pbs(ix,iy,0))*     !iyv 19.04.07 11.01.13
     &    sqrt((pb(leftix(ix,iy),leftiy(ix,iy))+pb(ix,iy))/
     &    (mb*(nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))))*
     &    min(vol(leftix(ix,iy),leftiy(ix,iy)),vol(ix,iy))/
     &    (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
         termhz(ix,iy) = xfrhie*wrk2(ix,iy)*0.5_R8*abs(pbshz(ix,iy,0))  !srv 05.11.09 11.01.13
     &    *sqrt((pb(leftix(ix,iy),leftiy(ix,iy))+pb(ix,iy))/
     &    (mb*(nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))))*
     &    min(vol(leftix(ix,iy),leftiy(ix,iy)),vol(ix,iy))/
     &    (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
!         if (nint(zn(isb)).eq.1.and.nint(zamax(isb)).eq.1) then
!         if(isb.eq.ismain) then                                        !srv 19.02.21
         if(.not.is_neutral(isb)) then                                  !srv 08.09.21
          scur =
     &     (fchinert_b(ix,iy,0)+fchvispar_b(ix,iy,0)
     &     +fchvisper_b(ix,iy,0)+fchin_b(ix,iy,0)
     &     +fchanml_b(ix,iy,0)+fchvisq_b(ix,iy,0))*xcur/
     &     (zb(ix,iy)*qe*
     &     (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8)          !srv 08.09.21
         else
          scur = 0.0_R8
         endif

         if(drift_style.eq.0)then                                       !srv 21.06.08
          t0 = (hz(ix,iy)*vol(leftix(ix,iy),leftiy(ix,iy))+
     &         hz(leftix(ix,iy),leftiy(ix,iy))*vol(ix,iy))/
     &        (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))
          flo = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+term(ix,iy)+     !srv 20.06.08
     1     (ubdia(ix,iy,0)+ubdia(leftix(ix,iy),leftiy(ix,iy),0))/
     5     2.0_R8*0.5_R8*(vol(ix,iy)/hx(ix,iy)
     6     +vol(leftix(ix,iy),leftiy(ix,iy))/
     7     hx(leftix(ix,iy),leftiy(ix,iy))) + scur
          flo_mdf = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+             !srv 17.07.08 {
     &      term(ix,iy) +                                               !iyv 19.04.07
     1     (vbecrb(ix,iy,0)+vbecrb(leftix(ix,iy),leftiy(ix,iy),0))/
     5     2.0_R8*0.5_R8*(vol(ix,iy)/hx(ix,iy)
     6     +vol(leftix(ix,iy),leftiy(ix,iy))/
     7     hx(leftix(ix,iy),leftiy(ix,iy))) + scur                      !srv 17.07.08 }
          flo_he = 5.0_R8/3.0_R8*cvlb(ix,iy,0)+                         !srv 13.09.02 { 01.07.05 20.01.11
     &     pbs(ix,iy,0)*wrk0(ix,iy)+term(ix,iy) +                       !iyv 30.05.08 !srv 03.10.08
     1     (vbecrb(ix,iy,0)+vbecrb(leftix(ix,iy),leftiy(ix,iy),0))/
     5     2.0_R8*0.5_R8*(vol(ix,iy)/hx(ix,iy)
     6     +vol(leftix(ix,iy),leftiy(ix,iy))/
     7     hx(leftix(ix,iy),leftiy(ix,iy)))                             !srv 17.07.08 }
          flo_53 = pbs(ix,iy,0)*wrk0(ix,iy)+term(ix,iy) +               !ajs 06.04.20 Part of the convective flow added in total energy flux
     1     (vbecrb(ix,iy,0)+vbecrb(leftix(ix,iy),leftiy(ix,iy),0))/
     2     2.0_R8*0.5_R8*(vol(ix,iy)/hx(ix,iy)
     3     +vol(leftix(ix,iy),leftiy(ix,iy))/
     7     hx(leftix(ix,iy),leftiy(ix,iy)))
          flo2dia = cvlbhz(ix,iy,0)+pbshz(ix,iy,0)*wrk0(ix,iy)+         !srv 01.10.99 01.07.05
     &     termhz(ix,iy)+                                               !srv 05.11.09
     1     (ub2dia(ix,iy,0)+ub2dia(leftix(ix,iy),leftiy(ix,iy),0))/
     5     2.0_R8*0.5_R8*(hz(ix,iy)*vol(ix,iy)/hx(ix,iy)
     6     +hz(leftix(ix,iy),leftiy(ix,iy))*
     7     vol(leftix(ix,iy),leftiy(ix,iy))/
     8     hx(leftix(ix,iy),leftiy(ix,iy))) + scur*t0
          flo_nodrift = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+
     &                  term(ix,iy)                                     !srv 20.06.08
         elseif (drift_style .eq. 1 .or. redef_gmtry .eq. 0) then       !srv 21.06.08 {
          t0 = (hz(ix,iy)*vol(leftix(ix,iy),leftiy(ix,iy))+
     &         hz(leftix(ix,iy),leftiy(ix,iy))*vol(ix,iy))/
     &        (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))
          flo = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+                 !srv 01.07.05
     &      term(ix,iy) +                                               !iyv 19.04.07
     4      ubdia(ix,iy,0)*
     5      (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/              !srv 01.03.00
     6      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy))) + scur          !srv 01.03.00
          flo_mdf = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+             !srv 13.10.06 {
     &      term(ix,iy) +                                               !iyv 19.04.07
     &      vbecrb(ix,iy,0)*
     &      (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     &      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy))) + scur
          flo_he = 5.0_R8/3.0_R8*cvlb(ix,iy,0)+                         !srv 13.09.02 { 01.07.05 20.01.11
     &     pbs(ix,iy,0)*wrk0(ix,iy)+term(ix,iy) +                       !iyv 30.05.08 !srv 03.10.08
     4     vbecrb(ix,iy,0)*
     5     (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     6      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))
          flo_53 = pbs(ix,iy,0)*wrk0(ix,iy)+term(ix,iy) +               !ajs 06.04.20 Part of the convective flow added in total energy flux
     1     vbecrb(ix,iy,0)*
     2     (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     3      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))
          flo2dia = cvlbhz(ix,iy,0)+pbshz(ix,iy,0)*wrk0(ix,iy)+
     &      termhz(ix,iy) +
     &   (ub2dia(ix,iy,0)/(hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy))) !srv this is correct. Keep in mind hz*vol/hx and interpolation of hz
     &      + scur / (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy))))
     &      *(hz(ix,iy)*vol(leftix(ix,iy),leftiy(ix,iy))+
     &      hz(leftix(ix,iy),leftiy(ix,iy))*vol(ix,iy))
          flo_nodrift = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy) +
     &                  term(ix,iy)                                     !srv 20.06.08 20.10.09
         elseif (drift_style .eq. 2 .and. redef_gmtry .eq. 1) then      !IYS 04.03.2019 - 11.02.2020
          t0 = hzl(ix,iy)
          flo = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+                 !srv 01.07.05
     &      term(ix,iy) +                                               !iyv 19.04.07
     &      ubdia(ix,iy,0)*gs(ix,iy,0)*qc(ix,iy) + scur                 !IYS 04.03.2019
          flo_mdf = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy)+             !srv 13.10.06 {
     &      term(ix,iy) +                                               !iyv 19.04.07
     &      vbecrb(ix,iy,0)*gs(ix,iy,0)*qc(ix,iy) + scur                !IYS 04.03.2019
          flo_he = 5.0_R8/3.0_R8*cvlb(ix,iy,0)+                         !srv 13.09.02 { 01.07.05 20.01.11
     &     pbs(ix,iy,0)*wrk0(ix,iy)+term(ix,iy) +                       !iyv 30.05.08 !srv 03.10.08
     &     vbecrb(ix,iy,0)*gs(ix,iy,0)*qc(ix,iy)                        !IYS 04.03.2019
          flo_53 = pbs(ix,iy,0)*wrk0(ix,iy)+term(ix,iy) +               !ajs 06.04.20 Part of the convective flow added in total energy flux
     1     vbecrb(ix,iy,0)*gs(ix,iy,0)*qc(ix,iy)
          flo2dia = cvlbhz(ix,iy,0)+pbshz(ix,iy,0)*wrk0(ix,iy)+
     &     termhz(ix,iy) +
     &     (ub2dia(ix,iy,0)*gs(ix,iy,0)*qc(ix,iy) + scur)*hzl(ix,iy)    ! IYS 04.03.2019 - 11.02.2020
!     &      *(hz(ix,iy) + hz(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8      !IYS 04.03.2019 - 11.02.2020
          flo_nodrift = cvlb(ix,iy,0)+pbs(ix,iy,0)*wrk0(ix,iy) +
     &                  term(ix,iy)                                     !srv 20.06.08 20.10.09
         endif                                                          !srv 21.06.08 }
!         flo = flo + scur
!         flo_mdf = flo_mdf + scur
!         flo2dia = flo2dia + scur*t0
         if (no_hybr.eq.0) then
           if(no_Ptncr_x_co.eq.0) then                                  !srv 13.10.06 {
             con_mdf(ix,iy)  = hybr(flo_mdf,cdnb(ix,iy,0))
             con(ix,iy)      = con_mdf(ix,iy)
             con_he(ix,iy)   = con_mdf(ix,iy)
             con_nodrift     = hybr(flo_nodrift,cdnb(ix,iy,0))          !did not change
             con_2dia(ix,iy) = con_mdf(ix,iy)*t0
           else
             con_mdf(ix,iy)  = cdnb(ix,iy,0)
             con_he(ix,iy)   = cdnb(ix,iy,0)                            !iyv 30.05.08 !srv 03.10.08
             con(ix,iy)      = cdnb(ix,iy,0)
             con_nodrift     = cdnb(ix,iy,0)
             con_2dia(ix,iy) = cdnbhz(ix,iy,0)                          !srv 09.01.01
           endif                                                        !srv 13.10.06 }
           fnb_nodrift(ix,iy,0) = flo_nodrift*
     1      (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &      con_nodrift*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
           fnb(ix,iy,0) =
     1      flo*(nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &      con(ix,iy)*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
           fnb_fcor(ix,iy,0) =flo2dia*(nb(leftix(ix,iy),leftiy(ix,iy))  !srv 01.10.99
     &      +nb(ix,iy))/2.0_R8-
     &      con_2dia(ix,iy)*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpbhz(ix,iy,0)*
     &     (pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))                  !srv 01.10.99 09.01.01
           fnb_mdf(ix,iy,0) = flo_mdf*(nb(leftix(ix,iy),leftiy(ix,iy))  !srv 13.10.06 {
     &      +nb(ix,iy))/2.0_R8-
     &      con_mdf(ix,iy)*
     &           (nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))-
     &      fnbPSch(ix,iy,0)                                            !srv 13.10.06 }
           if (balance.and.(balance_netcdf.ne.0)) then
cdjm Jan2017 store components of fnb for balance
             fna_pinch(ix,iy,0,isb) = cvlb(ix,iy,0)*
     &         (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
             fna_pll(ix,iy,0,isb) = (pbs(ix,iy,0)*wrk0(ix,iy)+
     &                               term(ix,iy))*
     &         (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
             if (mdf_fnb.eq.0) then
               if (drift_style.eq.0) then
                fna_drift(ix,iy,0,isb) =
     &           (ubdia(ix,iy,0)+
     &            ubdia(leftix(ix,iy),leftiy(ix,iy),0))/
     &           2.0_R8*0.5_R8*(vol(ix,iy)/hx(ix,iy)+
     &                          vol(leftix(ix,iy),leftiy(ix,iy))/
     &                          hx(leftix(ix,iy),leftiy(ix,iy)))*
     &           (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
               elseif (drift_style.eq.1 .or. redef_gmtry.eq.0) then
                fna_drift(ix,iy,0,isb) = ubdia(ix,iy,0)*
     &           (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     &           (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*
     &           (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
               elseif (drift_style.eq.2 .and. redef_gmtry.eq.1) then
                fna_drift(ix,iy,0,isb) = ubdia(ix,iy,0)*
     &           gs(ix,iy,0)*qc(ix,iy)*
     &           (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
               endif
               fna_nanom(ix,iy,0,isb) = -con(ix,iy)*
     &          (nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))
               fna_pschused(ix,iy,0,isb) = 0.0_R8
             else
               if (drift_style.eq.0) then
                fna_drift(ix,iy,0,isb) =
     &           (vbecrb(ix,iy,0)+
     &            vbecrb(leftix(ix,iy),leftiy(ix,iy),0))/
     &           2.0_R8*0.5_R8*(vol(ix,iy)/hx(ix,iy)+
     &                          vol(leftix(ix,iy),leftiy(ix,iy))/
     &                          hx(leftix(ix,iy),leftiy(ix,iy)))*
     &           (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
               elseif (drift_style.eq.1 .or. redef_gmtry.eq.0) then
                fna_drift(ix,iy,0,isb) = vbecrb(ix,iy,0)*
     &           (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     &           (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*
     &           (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
               elseif (drift_style.eq.2 .and. redef_gmtry.eq.1) then
                fna_drift(ix,iy,0,isb) = vbecrb(ix,iy,0)*
     &                    gs(ix,iy,0)*qc(ix,iy)*                       !IYS 04.03.2019
     &           (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
               endif
               fna_nanom(ix,iy,0,isb) = -con_mdf(ix,iy)*
     &          (nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))
               fna_pschused(ix,iy,0,isb) = -fnbPSch(ix,iy,0)
             endif
             fna_panom(ix,iy,0,isb) = -cdpb(ix,iy,0)*
     &        (pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
             fna_ch(ix,iy,0,isb) = 0.0_R8
           endif
*      fnb_he will be used only for heat equations                      !srv 21.01.00 {
           fnb_he(ix,iy,0) = flo_he*(nb(leftix(ix,iy),leftiy(ix,iy))    !srv 13.09.02
     &      +nb(ix,iy))/2.0_R8-
csw 30sep2014 (is this 5/3 factor correct here? I think not...)
     &      (5.0e0_R8/3.0e0_R8)*
     &      (con_he(ix,iy)*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))+ !srv 13.09.02
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy))))  !srv 21.01.00 }
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
           fnb_53(ix,iy,0) = flo_53*(nb(leftix(ix,iy),leftiy(ix,iy))    !ajs 06.04.20
     &      +nb(ix,iy))/2.0_R8
           if(mdf_fnb.eq.0) then                                        !srv 13.10.06 {
             flob(ix,iy,0) = flo*
     &        (kmpr(leftix(ix,iy),leftiy(ix,iy))+kmpr(ix,iy))/2.0_R8-
     &        con(ix,iy)*
     &        (kmpr(ix,iy)-kmpr(leftix(ix,iy),leftiy(ix,iy)))
             conb(ix,iy,0,0) = con(ix,iy)*                              !srv 13.10.06
     &        (kmpr(leftix(ix,iy),leftiy(ix,iy))+kmpr(ix,iy))/2.0_R8-   !srv 30.01.10
     &        flo*
     &        (kmpr(ix,iy)-kmpr(leftix(ix,iy),leftiy(ix,iy)))/4.0_R8
     &        +cdpb(ix,iy,0)
           else
             flob(ix,iy,0) = flo_mdf*
     &         (kmpr(leftix(ix,iy),leftiy(ix,iy))+kmpr(ix,iy))/2.0_R8 -
     &          con_mdf(ix,iy)*
     &         (kmpr(ix,iy)-kmpr(leftix(ix,iy),leftiy(ix,iy)))
     &         -cddb(ix,iy,1)*dkmprti(ix,iy,1)                          !srv 21.10.10
             conb(ix,iy,0,0) =
     &          con_mdf(ix,iy)*                          !srv 13.10.06
     &         (kmpr(leftix(ix,iy),leftiy(ix,iy))+kmpr(ix,iy))/2.0_R8-
     &          flo_mdf*
     &         (kmpr(ix,iy)-kmpr(leftix(ix,iy),leftiy(ix,iy)))/4.0_R8
     &         +cdpb(ix,iy,0)
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
           if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
             fna_pinch(ix,iy,0,isb) = cvlb(ix,iy,0)*t2
             fna_pll(ix,iy,0,isb) = (pbs(ix,iy,0)*wrk0(ix,iy)+
     &                               term(ix,iy))*t2
             if (drift_style.eq.0) then
              fna_drift(ix,iy,0,isb) =
     &           (ubdia(ix,iy,0)+
     &            ubdia(leftix(ix,iy),leftiy(ix,iy),0))/
     &           2.0_R8*0.5_R8*(vol(ix,iy)/hx(ix,iy)+
     &                          vol(leftix(ix,iy),leftiy(ix,iy))/
     &                          hx(leftix(ix,iy),leftiy(ix,iy)))*t2
             elseif (drift_style .eq. 1 .or. redef_gmtry .eq. 0) then
              fna_drift(ix,iy,0,isb) = ubdia(ix,iy,0)*
     &           (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     &           (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*t2
             elseif (drift_style .eq. 2 .and. redef_gmtry .eq. 1) then
              fna_drift(ix,iy,0,isb) = ubdia(ix,iy,0)*
     &           gs(ix,iy,0)*qc(ix,iy)*t2
             endif
             fna_nanom(ix,iy,0,isb) = -con(ix,iy)*
     &          (nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))
             fna_pschused(ix,iy,0,isb) = 0.0_R8
             fna_panom(ix,iy,0,isb) = -cdpb(ix,iy,0)*
     &        (pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
             fna_ch(ix,iy,0,isb) = scur*(nb(leftix(ix,iy),leftiy(ix,iy)) 
     &         + nb(ix,iy))/2.0_R8
           endif
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
     &      (con_he(ix,iy)*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))+  !srv 13.09.02
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy))) )  !srv 21.01.00 }
!         conb(ix,iy,0,0) = (CCX(IX,IY)+CCX1(IX,IY)) * SX(IX,IY)
           fnb_32(ix,iy,0) = flo * t2
           fnb_52(ix,iy,0) = -
     &      con(ix,iy)*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
           fnb_32nd(ix,iy,0) = flo_nodrift * t3
           fnb_52nd(ix,iy,0) = -
     &      con_nodrift*(nb(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy)))-
     &      cdpb(ix,iy,0)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))
           fnb_53(ix,iy,0) = flo_53*(nb(leftix(ix,iy),leftiy(ix,iy))    !ajs 06.04.20
     &      +nb(ix,iy))/2.0_R8
           flob(ix,iy,0) = flo*(kmpr(leftix(ix,iy),leftiy(ix,iy))
     &       +kmpr(ix,iy))/2.0_R8
           conb(ix,iy,0,0) =                                            !srv 13.10.06
     &      con(ix,iy)*(kmpr(leftix(ix,iy),leftiy(ix,iy))+
     &                 kmpr(ix,iy))/2.0_R8 + cdpb(ix,iy,0)
         endif
         if ((alpha.gt.0.0_R8 .and. zb(ix,iy).eq.0.0_R8) .and.          !srv 20.11.08
     *        alpha_.eq.0.0_R8 ) then                            ! neutral flux limit
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
           conb(ix,iy,0,0:4) = conb(ix,iy,0,0:4) / t0                   !srv 13.10.06 {
           fnb_he(ix,iy,0)   = fnb_he(ix,iy,0)/t0
           fnb_53(ix,iy,0)   = fnb_53(ix,iy,0)/t0                       !ajs 06.04.20
           fnb_fcor(ix,iy,0) = fnb_fcor(ix,iy,0)/t0
           fnb_mdf(ix,iy,0) = fnb_mdf(ix,iy,0)/t0
           fnbPSch(ix,iy,0)  = fnbPSch(ix,iy,0)/t0                      !srv 13.10.06 }
           fllim0fnb(ix,iy,0) = t0
         endif
        endif
       enddo
      enddo

      if((iout_b2npc9.ne.0 .and. balance) .or.
     &   (iout_b2npco.ne.0 .and. balance) .or.
     &   (iout_b2npc7.ne.0 .and. balance) .or.
     &   (iout_b2npc11.ne.0 .and. balance) .or.                         !srv 02.04.18
     &    iout .eq. 1 .or. (iout_b2wdat.eq.4 .and. balance)) then       !srv 11.09.09 {
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          wrk1(ix,iy) = 0.0_R8
          wrk2(ix,iy) = 0.0_R8
          wrk3(ix,iy) = 0.0_R8
          wrk4(ix,iy) = 0.0_R8
          wrk5(ix,iy) = 0.0_R8
          wrk6(ix,iy) = 0.0_R8
          wrk7(ix,iy) = 0.0_R8
          wrk8(ix,iy) = 0.0_R8                                          !iys 17.04.2012
         else
          wrk1(ix,iy) = pbs(ix,iy,0)*wrk0(ix,iy)*
     *     (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
          wrk2(ix,iy) = term(ix,iy)*
     *     (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
          if (drift_style.eq.0) then
            wrk3(ix,iy) =
     1         (vbecrb(ix,iy,0)+
     2          vbecrb(leftix(ix,iy),leftiy(ix,iy),0))/
     5         2.0_R8*0.5_R8*(vol(ix,iy)/hx(ix,iy)
     6         +vol(leftix(ix,iy),leftiy(ix,iy))/
     7         hx(leftix(ix,iy),leftiy(ix,iy)))*
     &         (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
          elseif (drift_style .eq. 1 .or. redef_gmtry .eq. 0) then
            wrk3(ix,iy) = vbecrb(ix,iy,0) *
     &         (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     &         (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy))) *            !IYS 04.03.2019
     &         (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
          elseif (drift_style .eq. 2 .and. redef_gmtry .eq. 1) then
            wrk3(ix,iy) = vbecrb(ix,iy,0) *
     &          gs(ix,iy,0) * qc(ix,iy) *                               !IYS 04.03.2019
     &         (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
          endif
          wrk4(ix,iy) = -con_2dia(ix,iy)*(nb(ix,iy)-
     &     nb(leftix(ix,iy),leftiy(ix,iy)))
          wrk5(ix,iy) = -cdpb(ix,iy,0)*(pb(ix,iy)-
     &     pb(leftix(ix,iy),leftiy(ix,iy)))
          wrk6(ix,iy) = termhz(ix,iy)*
     *     (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
          wrk7(ix,iy) = -con_mdf(ix,iy)*(nb(ix,iy)-
     &     nb(leftix(ix,iy),leftiy(ix,iy)))
          wrk8(ix,iy) = cvlb(ix,iy,0)*                                  !iys 17.04.2012
     &     (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
         endif
        enddo
       enddo
       if(.not. is_neutral(isb)) then
        call my_out(70,nx,ny,wrk3,'b2tfnb_vaecrbnax'//chns)
        call my_out(70,nx,ny,fnbPSch(-1,-1,0),'b2tfnb_fnbPSchx'//chns)
        call my_out(70,nx,ny,wrk7,'b2tfnb_dPat_mdf_gradnax'//chns)
        call my_out(70,nx,ny,wrk4,'b2tfnb_dPat_2diagradnax'//chns)
       else
        call my_out(70,nx,ny,wrk5,'b2tfnb_dgradpbx'//chns)
       endif
       call my_out(70,nx,ny,wrk1,'b2tfnb_bxuanax'//chns)
       call my_out(70,nx,ny,wrk2,'b2tfnb_dpccornax'//chns)
       if (iout_b2wdat.eq.4 .and. balance)
     &           call my_out(70,nx,ny,wrk6,'b2tfnb_dpccorhznax'//chns)  !srv 05.11.09
       if (iout.eq.1) call my_out(70,nx,ny,wrk8,'b2tfnb_cvlbnax'//chns) !iys 17.04.2012
!       call my_out(70,nx,ny,con_2dia,'b2tfnb_dPat_2diax'//chns)
!       call my_out(70,nx,ny,con_mdf,'b2tfnb_dPat_mdfx'//chns)
      endif
      if((iout_b2npmo.ne.0 .and. balance)
     &   .or. (iout_b2wdat.eq.4 .and. balance) .or. iout.eq.1) then
       write(chns,'(i3.3)') isb
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          wrk6(ix,iy) = 0.0_R8
          wrk8(ix,iy) = 0.0_R8
         else
          if (drift_style.eq.0) then
           wrk6(ix,iy) =
     &      (vbdia(ix,iy,0)+vbdia(leftix(ix,iy),leftiy(ix,iy),0))/
     5      2.0_R8*0.5_R8*(vol(ix,iy)/hx(ix,iy)
     6      +vol(leftix(ix,iy),leftiy(ix,iy))/
     7      hx(leftix(ix,iy),leftiy(ix,iy)))*
     &      (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
           wrk8(ix,iy) =
     &      (vbdia(ix,iy,0)+vbdia(leftix(ix,iy),leftiy(ix,iy),0))/
     5      2.0_R8*0.5_R8*(hz(ix,iy)*vol(ix,iy)/hx(ix,iy)
     6      +hz(leftix(ix,iy),leftiy(ix,iy))*
     7      vol(leftix(ix,iy),leftiy(ix,iy))/
     8      hx(leftix(ix,iy),leftiy(ix,iy)))*
     &      (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
          elseif (drift_style .eq. 1 .or. redef_gmtry .eq. 0) then
           wrk6(ix,iy) = vbdia(ix,iy,0)*
     5      (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     6      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*
     &      (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
           wrk8(ix,iy) = vbdia(ix,iy,0)*
     5      (hz(ix,iy)*vol(ix,iy)+hz(leftix(ix,iy),leftiy(ix,iy))
     5      *vol(leftix(ix,iy),leftiy(ix,iy)))/
     6      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*
     &      (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
          elseif (drift_style .eq. 2 .and. redef_gmtry .eq. 1) then
           wrk6(ix,iy) = vbdia(ix,iy,0)*gs(ix,iy,0)*qc(ix,iy)*          ! IYS 07.05.2020
     &      (nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
           wrk8(ix,iy)= vbdia(ix,iy,0)*
     &      gs(ix,iy,0)*qc(ix,iy)*hzl(ix,iy)                            ! IYS 07.05.2020
     &      *(nb(leftix(ix,iy),leftiy(ix,iy))+nb(ix,iy))/2.0_R8
          endif
         endif
        enddo
       enddo
       if(.not. is_neutral(isb)) then
        call my_out(70,nx,ny,wrk6,'b2tfnb_vadianax'//chns)
        call my_out(70,nx,ny,wrk8,'b2tfnb_vadianahzx'//chns)
       endif
       call my_out(70,nx,ny,fnb_fcor(-1,-1,0),
     &              'b2tfnb_fnb_fcorx'//chns)
      endif                                                             !srv 11.09.09 }

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
         fnb_53(ix,iy,1)      = 0.0_R8                                  !ajs 06.04.20
         flob(ix,iy,1)        = 0.0_R8
         fnb_fcor(ix,iy,1)    = 0.0_R8                                  !srv 01.10.99
         fnb_he(ix,iy,1)      = 0.0_R8                                  !srv 21.01.00
         fnb_mdf(ix,iy,1)     = 0.0_R8                                  !srv 13.10.06
         conb(ix,iy,1,0:2)    = 0.0_R8                                  !srv 13.10.06
         con(ix,iy)           = 0.0_R8
         con_2dia(ix,iy)      = 0.0_R8
         con_he(ix,iy)        = 0.0_R8
         con_mdf(ix,iy)       = 0.0_R8
         fllim0fnb(ix,iy,1)   = 1.0_R8
         if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
          fna_pinch(ix,iy,1,isb)    = 0.0_R8
          fna_pll(ix,iy,1,isb)      = 0.0_R8
          fna_drift(ix,iy,1,isb)    = 0.0_R8
          fna_nanom(ix,iy,1,isb)    = 0.0_R8
          fna_panom(ix,iy,1,isb)    = 0.0_R8
          fna_ch(ix,iy,1,isb)       = 0.0_R8
          fna_pschused(ix,iy,1,isb) = 0.0_R8
         endif
        else
!         if(nint(zn(isb)).eq.1.and.nint(zamax(isb)).eq.1) then         !srv 20.09.06
!         if(isb.eq.ismain) then                                        !srv 19.02.21
         if(.not.is_neutral(isb)) then                                  !srv 08.09.21
          scur =
     &    (fchinert_b(ix,iy,1)+fchvispar_b(ix,iy,1)+fchin_b(ix,iy,1)+
     &     fchanml_b(ix,iy,1)+vis_q*fchvisq_b(ix,iy,1))*ycur/           !srv 17.07.08 {  02.10.17
     &    (zb(ix,iy)*qe*
     &    (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8)       !srv 17.07.08 }
          vsprcur = vis_per*fchvisper_b(ix,iy,1)/(zb(ix,iy)*qe)*ycur    !srv 02.04.18 08.09.21
         else
           scur = 0.0_R8
           vsprcur = 0.0_R8
         endif                                                          !srv 08.09.21 }
         if(drift_style.eq.0)then                                       !srv 21.06.08
          t0 = (hz(ix,iy)*vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &          hz(bottomix(ix,iy),bottomiy(ix,iy))*vol(ix,iy))/
     &         (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))
          flo = cvlb(ix,iy,1)+
     &      0.5_R8*(ubdia(ix,iy,1)*vol(ix,iy)/hy1(ix,iy)
     &      +ubdia(bottomix(ix,iy),bottomiy(ix,iy),1)*
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
     &       /hy1(bottomix(ix,iy),bottomiy(ix,iy)))
     &      +scur
          flo_mdf = cvlb(ix,iy,1)+
     &      0.5_R8*(vbecrb(ix,iy,1)*vol(ix,iy)/hy1(ix,iy)
     &      +vbecrb(bottomix(ix,iy),bottomiy(ix,iy),1)*
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
     &       /hy1(bottomix(ix,iy),bottomiy(ix,iy)))
     &      +scur
          flo2dia = cvlbhz(ix,iy,1)+
     &       (ub2dia(ix,iy,1)+
     &        ub2dia(bottomix(ix,iy),bottomiy(ix,iy),1))
     &       /2.0_R8*0.5_R8*(hz(ix,iy)*vol(ix,iy)/hy1(ix,iy)+           !srv 01.10.99 09.01.01
     &        hz(bottomix(ix,iy),bottomiy(ix,iy))
     5                        *vol(bottomix(ix,iy),bottomiy(ix,iy))/
     &        hy1(bottomix(ix,iy),bottomiy(ix,iy)))+                    !srv 01.03.00 09.01.01
     &        (hz(ix,iy)+hz(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8*   !srv 18.09.08
     &         scur
          flo_he = 5.0_R8/3.0_R8*cvlb(ix,iy,1)+                         !srv 13.09.02 ( 20.01.11
     &      0.5_R8*(vbecrb(ix,iy,1)*vol(ix,iy)/hy1(ix,iy)
     &      +vbecrb(bottomix(ix,iy),bottomiy(ix,iy),1)*
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
     &       /hy1(bottomix(ix,iy),bottomiy(ix,iy)))
          flo_53 = 0.5_R8*(vbecrb(ix,iy,1)*vol(ix,iy)/hy1(ix,iy)        !ajs 06.04.20
     &      +vbecrb(bottomix(ix,iy),bottomiy(ix,iy),1)*
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
     &       /hy1(bottomix(ix,iy),bottomiy(ix,iy)))
         elseif (drift_style.eq.1 .or. redef_gmtry.eq.0) then           !srv 21.06.08 {
           t0 = (hz(ix,iy)*vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &           hz(bottomix(ix,iy),bottomiy(ix,iy))*vol(ix,iy))/
     &          (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))
           flo = cvlb(ix,iy,1)+ubdia(ix,iy,1)*
     5      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/          !srv 01.03.00
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))           !srv 01.03.00
     &      +scur
           flo_mdf = cvlb(ix,iy,1)+vbecrb(ix,iy,1)*                     !srv 01.06.06 {
     5      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
     &      +scur
           flo2dia = cvlbhz(ix,iy,1)+ub2dia(ix,iy,1)*
     5      (hz(ix,iy)*vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &       hz(bottomix(ix,iy),bottomiy(ix,iy))*vol(ix,iy))/
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))+          !srv it is correct. Keep in mind hz*vol/hy and interpolation of hz
     &      (hz(ix,iy)*vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &       hz(bottomix(ix,iy),bottomiy(ix,iy))*vol(ix,iy))/
     &      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))*
     &       scur
           flo_he = 5.0_R8/3.0_R8*cvlb(ix,iy,1)+vbecrb(ix,iy,1)*        !srv 13.09.02 ( 20.01.11
     5      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
           flo_53 = vbecrb(ix,iy,1)*                                    !ajs 06.04.20
     5      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
         elseif (drift_style.eq.2 .and. redef_gmtry.eq.1) then          !IYS 04.03.2019 - 11.02.2020
           t0 = hzb(ix,iy)
           flo = cvlb(ix,iy,1)+ubdia(ix,iy,1)*gs(ix,iy,1)               !IYS 04.03.2019
     &      +scur
           flo_mdf = cvlb(ix,iy,1)+vbecrb(ix,iy,1)*gs(ix,iy,1)          !IYS 04.03.2019, srv 01.06.06 {
     &      +scur
           flo2dia = cvlbhz(ix,iy,1)+                                   !IYS 04.03.2019
     &       (ub2dia(ix,iy,1)*gs(ix,iy,1) + scur) * hzb(ix,iy)          !IYS 04.03.2019
!     &      (hz(ix,iy)*vol(bottomix(ix,iy),bottomiy(ix,iy))+
!     &       hz(bottomix(ix,iy),bottomiy(ix,iy))*vol(ix,iy))/
!     &      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))*
!     &       scur
           flo_he = 5.0_R8/3.0_R8*cvlb(ix,iy,1)+
     &      vbecrb(ix,iy,1)*gs(ix,iy,1)          !IYS 04.03.2019        !srv 13.09.02 ( 20.01.11
           flo_53 = vbecrb(ix,iy,1)*gs(ix,iy,1)
         endif                                                          !srv 21.06.08 }
         flo_nodrift = cvlb(ix,iy,1)
         if (no_hybr.eq.0) then
          if(no_Ptncr_y_co.eq.0) then                                   !srv 06.07.06 {
            con_mdf(ix,iy)  = hybr(flo_mdf,cdnb(ix,iy,1))
            con(ix,iy)      = con_mdf(ix,iy)
            con_he(ix,iy)   = con_mdf(ix,iy)
            con_2dia(ix,iy) = con_mdf(ix,iy)*t0
            con_nodrift     = hybr(flo_nodrift,cdnb(ix,iy,1))           ! did not change
          else
           con_he(ix,iy)   = cdnb(ix,iy,1)                !iyv 30.05.08 !srv 03.10.08
           con_2dia(ix,iy) = cdnbhz(ix,iy,1)                            !srv 01.10.99 09.01.01
           con_mdf(ix,iy)  = cdnb(ix,iy,1)                              !srv 13.10.06
           con(ix,iy)      = cdnb(ix,iy,1)
           con_nodrift     = cdnb(ix,iy,1)
          endif                                                         !srv 06.07.06 }
          fnb(ix,iy,1) =
     1      flo*(nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &      con(ix,iy)*(nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &      cdpb(ix,iy,1)*(pb(ix,iy)
     &                    -pb(bottomix(ix,iy),bottomiy(ix,iy)))+vsprcur !srv 02.04.18
*      fnb_he will be used only for heat equations                      !srv 21.01.00 {
          fnb_he(ix,iy,1) = flo_he*
     &      (nb(bottomix(ix,iy),bottomiy(ix,iy)) + nb(ix,iy))/2.0_R8 -
csw 30sep2014 (is this 5/3 factor correct here? I think not...)
     &      (5.0_R8/3.0_R8)*
     &      (con_he(ix,iy)*
     &      (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy))) +
     &      cdpb(ix,iy,1)*(pb(ix,iy)
     &                   -pb(bottomix(ix,iy),bottomiy(ix,iy))) )        !srv 21.01.00 }
          fnb_fcor(ix,iy,1) = flo2dia*                                  !srv 01.10.99
     &       (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8 -
     &       con_2dia(ix,iy)*(nb(ix,iy)-
     &                        nb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &       cdpbhz(ix,iy,1)*(pb(ix,iy)                                 !srv 09.01.01
     &                     -pb(bottomix(ix,iy),bottomiy(ix,iy)))+       !srv 01.10.99
     &      t0*vsprcur                                                  !srv 02.04.18
          fnb_nodrift(ix,iy,1) = flo_nodrift*
     1      (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &       con_nodrift*
     &       (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &       cdpb(ix,iy,1)*(pb(ix,iy)
     &                     -pb(bottomix(ix,iy),bottomiy(ix,iy)))
          fnb_mdf(ix,iy,1) = flo_mdf*                                   !srv 13.10.06 {
     &     (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &     con_mdf(ix,iy)*(nb(ix,iy)-
     &                     nb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &     cdpb(ix,iy,1)*(pb(ix,iy)
     &                   -pb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &                    fnbPSch(ix,iy,1) + vsprcur                    !srv 13.10.06 } 02.04.18
          fnb_53(ix,iy,1) = flo_53*                                     !ajs 06.04.20
     &      (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
          if (balance.and.(balance_netcdf.ne.0)) then
cdjm Jan2017 store components of fnb for balance
            fna_pinch(ix,iy,1,isb) = cvlb(ix,iy,1)*
     &       (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
            fna_pll(ix,iy,1,isb) = 0.0_R8
            if (mdf_fnb.eq.0) then
              if (drift_style.eq.0) then
                fna_drift(ix,iy,1,isb) =
     &              0.5_R8*(ubdia(ix,iy,1)*vol(ix,iy)/hy1(ix,iy)+
     &                      ubdia(bottomix(ix,iy),bottomiy(ix,iy),1)*
     &                      vol(bottomix(ix,iy),bottomiy(ix,iy))/
     &                      hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &            (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
              elseif (drift_style.eq.1 .or. redef_gmtry.eq.0) then
                fna_drift(ix,iy,1,isb) = ubdia(ix,iy,1)*
     &            (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     &            (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &            (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
              elseif (drift_style.eq.2 .and. redef_gmtry.eq.1) then
                fna_drift(ix,iy,1,isb) = ubdia(ix,iy,1)*gs(ix,iy,1)*
     &            (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
              endif
            else
              if (drift_style.eq.0) then
                fna_drift(ix,iy,1,isb) =
     &              0.5_R8*(vbecrb(ix,iy,1)*vol(ix,iy)/hy1(ix,iy)+
     &                      vbecrb(bottomix(ix,iy),bottomiy(ix,iy),1)*
     &                      vol(bottomix(ix,iy),bottomiy(ix,iy))
     &                      /hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &            (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
              elseif (drift_style.eq.1 .or. redef_gmtry.eq.0) then
                fna_drift(ix,iy,1,isb) = vbecrb(ix,iy,1)*
     &            (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     &            (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &            (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
              elseif (drift_style.eq.2 .and. redef_gmtry.eq.1) then
                fna_drift(ix,iy,1,isb) = vbecrb(ix,iy,1)*gs(ix,iy,1)*    ! IYS 11.02.2020
     &            (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
              endif
            endif
!            if(nint(zn(isb)).eq.1.and.nint(zamax(isb)).eq.1) then
!            if(isb.eq.ismain) then                                      !srv 19.02.21
              fna_ch(ix,iy,1,isb) = scur*
     &            (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
!            else
!              fna_ch(ix,iy,1,isb) = 0.0_R8
!            endif
            if (mdf_fnb.eq.0) then
              fna_nanom(ix,iy,1,isb) = -con(ix,iy)*
     &          (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))
              fna_pschused(ix,iy,1,isb) = 0.0_R8
            else
              fna_nanom(ix,iy,1,isb) = -con_mdf(ix,iy)*
     &          (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))
              fna_pschused(ix,iy,1,isb) = -fnbPSch(ix,iy,1)
            endif
            fna_panom(ix,iy,1,isb) = -cdpb(ix,iy,1)*
     &          (pb(ix,iy)-pb(bottomix(ix,iy),bottomiy(ix,iy)))
          endif
          fnb_32(ix,iy,1) =
     &      flo*(nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &      abs(flo/2.0_R8)*
     &      (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))
          fnb_52(ix,iy,1) = -
     &      (con(ix,iy)-abs(flo/2.0_R8))*
     &      (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &      cdpb(ix,iy,1)*(pb(ix,iy)
     &                    -pb(bottomix(ix,iy),bottomiy(ix,iy)))
          fnb_32nd(ix,iy,1) =
     1      flo_nodrift*
     2      (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8-
     &      abs(flo_nodrift/2.0_R8)*
     1      (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))
          fnb_52nd(ix,iy,1) = -
     1      (con_nodrift-abs(flo_nodrift/2.0_R8))*
     &      (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))-
     &      cdpb(ix,iy,1)*(pb(ix,iy)
     &                    -pb(bottomix(ix,iy),bottomiy(ix,iy)))
          if(mdf_fnb.eq.0) then                                         !srv 13.10.06 {
           flob(ix,iy,1) = flo*(kmpr(bottomix(ix,iy),bottomiy(ix,iy))
     &        +kmpr(ix,iy))/2.0_R8
     &       -con(ix,iy)*(kmpr(ix,iy)-
     &                    kmpr(bottomix(ix,iy),bottomiy(ix,iy)))
           conb(ix,iy,1,0) =                                            !srv 13.10.06
     &        con(ix,iy)*(kmpr(bottomix(ix,iy),bottomiy(ix,iy))
     &        +kmpr(ix,iy))/2.0_R8-
     &        flo*(kmpr(ix,iy)
     &            -kmpr(bottomix(ix,iy),bottomiy(ix,iy)))/4.0_R8+
     &        cdpb(ix,iy,1)
          else
           flob(ix,iy,1) = flo_mdf*
     &       (kmpr(bottomix(ix,iy),bottomiy(ix,iy))+kmpr(ix,iy))/2.0_R8
     &      -con_mdf(ix,iy)*(kmpr(ix,iy)-
     &                       kmpr(bottomix(ix,iy),bottomiy(ix,iy)))
     &         +cddb(ix,iy,0)*dkmprti(ix,iy,0)                          !srv 21.10.10
           conb(ix,iy,1,0) =                                            !srv 13.10.06
     &        con_mdf(ix,iy)*(kmpr(bottomix(ix,iy),bottomiy(ix,iy))
     &        +kmpr(ix,iy))/2.0_R8-
     &        flo_mdf*(kmpr(ix,iy)
     &        -kmpr(bottomix(ix,iy),bottomiy(ix,iy)))/4.0_R8
     &        +cdpb(ix,iy,1)
           conb(ix,iy,1,1) = cdpb(ix,iy,1)*kmprti(ix,iy,1)
           conb(ix,iy,1,2) = conb(ix,iy,1,1)
          endif                                                         !srv 23.06.08 }
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
     &                    -pb(bottomix(ix,iy),bottomiy(ix,iy)))+vsprcur !srv 02.04.18
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
     &      nb(bottomix(ix,iy),bottomiy(ix,iy)))+                       !srv 13.09.02 }
     &      cdpb(ix,iy,1)*(pb(ix,iy)
     &      -pb(bottomix(ix,iy),bottomiy(ix,iy))) )                     !srv 21.01.00 }
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
           fnb_53(ix,iy,1) = flo_53*t0                                  !ajs 06.04.20
           conb(ix,iy,1,0) =                                            !srv 13.10.06
     &      con(ix,iy)*(kmpr(bottomix(ix,iy),bottomiy(ix,iy))+
     &      kmpr(ix,iy))/2.0_R8+cdpb(ix,iy,1)
         endif
         if ((alpha.gt.0.0_R8 .and. zb(ix,iy).eq.0.0_R8) .and.          !srv 20.11.08
     *       alpha_.eq.0.0_R8 ) then                             ! neutral flux limit
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
          conb(ix,iy,1,0:4) = conb(ix,iy,1,0:4) / t0                    !srv 13.10.06 {
          fnb_he(ix,iy,1)   = fnb_he(ix,iy,1)/t0
          fnb_53(ix,iy,1)   = fnb_53(ix,iy,1)/t0                        !ajs 06.04.20
          fnb_fcor(ix,iy,1) = fnb_fcor(ix,iy,1)/t0
          fnb_mdf(ix,iy,1) = fnb_mdf(ix,iy,1)/t0
          fnbPSch(ix,iy,1)  = fnbPSch(ix,iy,1)/t0                       !srv 13.10.06 }
          fllim0fnb(ix,iy,1) = t0
         endif
        endif
       enddo
      enddo

      if (alpha.gt.0.0_R8 .and. maxval(zb).eq.0.0_R8 .and. ncall.lt.100)
     > then
        write(*,'(a,1p,e12.5,a,e12.5,a,2e12.5)')
     .   'b2tfnb: maxfluxlimit for zb = ',maxval(zb),
     1   ' mb = ',mb/mp,' are ',maxfluxlimit
      endif

      if((iout_b2npc9.ne.0 .and. balance) .or.
     &   (iout_b2npco.ne.0 .and. balance) .or.
     &   (iout_b2npc7.ne.0 .and. balance) .or.
     &   (iout_b2npc11.ne.0 .and. balance) .or.                         !srv 02.04.18
     &    iout .eq. 1 .or. (iout_b2wdat.eq.4 .and. balance)) then       !srv 11.09.09 {
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).eq.-2) then
          wrk3(ix,iy) = 0.0_R8
          wrk4(ix,iy) = 0.0_R8
          wrk5(ix,iy) = 0.0_R8
          wrk7(ix,iy) = 0.0_R8
          wrk8(ix,iy) = 0.0_R8                                          !iys 17.04.2012
         else
           if (drift_style.eq.0) then
             wrk3(ix,iy) =
     &         0.5_R8*(vbecrb(ix,iy,1)*vol(ix,iy)/hy1(ix,iy)
     &         +vbecrb(bottomix(ix,iy),bottomiy(ix,iy),1)*
     &         vol(bottomix(ix,iy),bottomiy(ix,iy))
     &         /hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &         (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
           elseif (drift_style .eq. 1 .or. redef_gmtry .eq. 0) then
             wrk3(ix,iy) = vbecrb(ix,iy,1)*
     5         (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6         (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &         (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
           elseif (drift_style .eq. 2 .and. redef_gmtry .eq. 1) then
             wrk3(ix,iy) = vbecrb(ix,iy,1)*gs(ix,iy,1)*
     &         (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
           endif
           wrk4(ix,iy) = -con_2dia(ix,iy)*
     &                   (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))
           wrk5(ix,iy) = -cdpb(ix,iy,1)*(pb(ix,iy)
     &                    -pb(bottomix(ix,iy),bottomiy(ix,iy)))
           wrk7(ix,iy) = -con_mdf(ix,iy)*
     &                   (nb(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy)))
           wrk8(ix,iy) = cvlb(ix,iy,1)*                                !iys 17.04.2012
     &     (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
         endif
        enddo
       enddo
       write(chns,'(i3.3)') isb
       if(.not. is_neutral(isb)) then
        call my_out(70,nx,ny,fnbPSch(-1,-1,1),'b2tfnb_fnbPSchy'//chns)
        call my_out(70,nx,ny,wrk3,'b2tfnb_vaecrbnay'//chns)
        call my_out(70,nx,ny,wrk7,'b2tfnb_dPat_mdf_gradnay'//chns)
        call my_out(70,nx,ny,wrk4,'b2tfnb_dPat_2diagradnay'//chns)
       else
        call my_out(70,nx,ny,wrk5,'b2tfnb_dgradpby'//chns)
       endif
       call my_out(70,nx,ny,wrk8,'b2tfnb_cvlbnay'//chns)               !iys 17.04.2012
!       call my_out(70,nx,ny,con_2dia,'b2tfnb_dPat_2diay'//chns)
!       call my_out(70,nx,ny,con_mdf,'b2tfnb_dPat_mdfy'//chns)
      endif
      if((iout_b2npmo.ne.0 .and. balance) .or.
     &     (iout_b2wdat.eq.4 .and. balance) .or. iout.eq.1) then
       write(chns,'(i3.3)') isb
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).eq.-2) then
           wrk6(ix,iy) = 0.0_R8
           wrk8(ix,iy) = 0.0_R8
         else
          if (drift_style.eq.0) then
           wrk6(ix,iy) = 0.5_R8*(vbdia(ix,iy,1)*vol(ix,iy)/hy1(ix,iy)
     &      +vbdia(bottomix(ix,iy),bottomiy(ix,iy),1)*
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
     &       /hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &       (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
           wrk8(ix,iy) =
     &       (vbdia(ix,iy,1)
     &       +vbdia(bottomix(ix,iy),bottomiy(ix,iy),1))
     &       /2.0_R8*0.5_R8*(hz(ix,iy)*vol(ix,iy)/hy1(ix,iy)+           !srv 01.10.99 09.01.01
     &        hz(bottomix(ix,iy),bottomiy(ix,iy))
     5                        *vol(bottomix(ix,iy),bottomiy(ix,iy))/
     &        hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &       (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
          elseif (drift_style .eq. 1 .or. redef_gmtry .eq. 0) then
           wrk6(ix,iy) = vbdia(ix,iy,1)*
     5       (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6       (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &       (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
           wrk8(ix,iy) = vbdia(ix,iy,1)*
     5      (hz(ix,iy)*vol(ix,iy)+hz(bottomix(ix,iy),bottomiy(ix,iy))
     5      *vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &      (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
          elseif (drift_style .eq. 2 .and. redef_gmtry .eq. 1) then
            wrk6(ix,iy) = vbdia(ix,iy,1)*gs(ix,iy,1)*
     &       (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
            wrk8(ix,iy) = vbdia(ix,iy,1)*gs(ix,iy,1)*hzb(ix,iy)*
     &       (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
          endif
         endif
        enddo
       enddo
       if(.not. is_neutral(isb)) then
        call my_out(70,nx,ny,wrk6,'b2tfnb_vadianay'//chns)
        call my_out(70,nx,ny,wrk8,'b2tfnb_vadianahzy'//chns)
       endif
!       do iy = -1, ny
!        do ix = -1, nx
!         if(bottomiy(ix,iy).eq.-2) then
!          wrk6(ix,iy) = 0.0_R8
!         else
!          wrk6(ix,iy) = vbdia(ix,iy,1)*
!     5     (hz(ix,iy)*vol(ix,iy)+hz(bottomix(ix,iy),bottomiy(ix,iy))
!     5     *vol(bottomix(ix,iy),bottomiy(ix,iy)))/
!     6     (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
!     &     (nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2.0_R8
!         endif
!        enddo
!       enddo
!       if(.not. is_neutral(isb) .and. (iout_b2wdat.eq.4 .and. balance))
!     &  call my_out(70,nx,ny,wrk6,'b2tfnb_vadianahzy'//chns)
        call my_out(70,nx,ny,fnb_fcor(-1,-1,1),
     &                                        'b2tfnb_fnb_fcory'//chns)
      endif                                                             !srv 11.09.09 }

      if(guard_flows.eq.0) then                                         !xpb
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
            conb(ix,iy,1,:)=0.0_R8
            fllim0fnb(ix,iy,1) = 1.0_R8
            vbecrb(ix,iy,1)=0.0_R8
            vbdia(ix,iy,1)=0.0_R8
            wbdia(ix,iy,1)=0.0_R8
            ubdia(ix,iy,1)=0.0_R8
            if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
               fna_pinch(ix,iy,1,isb)=0.0_R8
               fna_pll(ix,iy,1,isb)=0.0_R8
               fna_drift(ix,iy,1,isb)=0.0_R8
               fna_nanom(ix,iy,1,isb)=0.0_R8
               fna_panom(ix,iy,1,isb)=0.0_R8
               fna_ch(ix,iy,1,isb)=0.0_R8
               fna_pschused(ix,iy,1,isb)=0.0_R8
            endif
          endif
         enddo
        enddo

cxpb turn off poloidal flows in core boundary
        if (boundary_namelist.eq.1) then
          do ib = 1, nbc
            if(core_boundary(ib)) then
              do ix = 1, bc_list_size(ib)
                fnb(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                fnb_32(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                fnb_52(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                fnb_nodrift(bc_list_x(ix,ib),bc_list_y(ix,ib),0) =
     .                                                          0.0_R8
                fnb_32nd(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                fnb_52nd(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                flob(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                conb(bc_list_x(ix,ib),bc_list_y(ix,ib),0,:) = 0.0_R8
                fllim0fnb(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 1.0_R8
                vbecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                vbdia(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                wbdia(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                ubdia(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
                  fna_pinch(bc_list_x(ix,ib),bc_list_y(ix,ib),
     .                                               0,isb)=0.0_R8
                  fna_pll(bc_list_x(ix,ib),bc_list_y(ix,ib),
     .                                               0,isb)=0.0_R8
                  fna_drift(bc_list_x(ix,ib),bc_list_y(ix,ib),
     .                                               0,isb)=0.0_R8
                  fna_nanom(bc_list_x(ix,ib),bc_list_y(ix,ib),
     .                                               0,isb)=0.0_R8
                  fna_panom(bc_list_x(ix,ib),bc_list_y(ix,ib),
     .                                               0,isb)=0.0_R8
                  fna_ch(bc_list_x(ix,ib),bc_list_y(ix,ib),
     .                                               0,isb)=0.0_R8
                  fna_pschused(bc_list_x(ix,ib),bc_list_y(ix,ib),
     .                                               0,isb)=0.0_R8
                endif
              enddo
            endif
          enddo
        else
          do ix = -1, nx
            if (on_closed_surface(ix,-1)) then
              fnb(ix,-1,0) = 0.0_R8
              fnb_32(ix,-1,0) = 0.0_R8
              fnb_52(ix,-1,0) = 0.0_R8
              fnb_nodrift(ix,-1,0) = 0.0_R8
              fnb_32nd(ix,-1,0) = 0.0_R8
              fnb_52nd(ix,-1,0) = 0.0_R8
              flob(ix,-1,0) = 0.0_R8
              conb(ix,-1,0,:)=0.0_R8
              fllim0fnb(ix,-1,0) = 1.0_R8
              vbecrb(ix,-1,0) = 0.0_R8
              vbdia(ix,-1,0) = 0.0_R8
              wbdia(ix,-1,0) = 0.0_R8
              ubdia(ix,-1,0) = 0.0_R8
              if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
                fna_pinch(ix,-1,0,isb)=0.0_R8
                fna_pll(ix,-1,0,isb)=0.0_R8
                fna_drift(ix,-1,0,isb)=0.0_R8
                fna_nanom(ix,-1,0,isb)=0.0_R8
                fna_panom(ix,-1,0,isb)=0.0_R8
                fna_ch(ix,-1,0,isb)=0.0_R8
                fna_pschused(ix,-1,0,isb)=0.0_R8
              endif
            endif
          enddo
        endif

      elseif(guard_flows.eq.1) then                                     !srv 20.06.08 {
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
            fnb_53(ix,iy,1)      = 0.0_R8                               !ajs 06.04.20
            fnb_mdf(ix,iy,1)     = 0.0_R8
            flob(ix,iy,1)        = 0.0_R8
            conb(ix,iy,1,0:2)    = 0.0_R8
            fllim0fnb(ix,iy,0)   = 1.0_R8
            if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
               fna_pinch(ix,iy,1,isb)=0.0_R8
               fna_pll(ix,iy,1,isb)=0.0_R8
               fna_drift(ix,iy,1,isb)=0.0_R8
               fna_nanom(ix,iy,1,isb)=0.0_R8
               fna_panom(ix,iy,1,isb)=0.0_R8
               fna_ch(ix,iy,1,isb)=0.0_R8
               fna_pschused(ix,iy,1,isb)=0.0_R8
            endif
          endif
         enddo
        enddo
*     turn off flows between boundary cells in south and in north
        do iy = -1, ny, ny+1
          fnb(:,iy,0)       = 0.0_R8
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
          fnb_53(:,iy,0)      = 0.0_R8                                  !ajs 06.04.20
          fnb_mdf(:,iy,0)     = 0.0_R8
          flob(:,iy,0)        = 0.0_R8
          conb(:,iy,0,0:2)    = 0.0_R8
          fllim0fnb(:,iy,0)   = 1.0_R8
          if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
             fna_pinch(:,iy,0,isb)=0.0_R8
             fna_pll(:,iy,0,isb)=0.0_R8
             fna_drift(:,iy,0,isb)=0.0_R8
             fna_nanom(:,iy,0,isb)=0.0_R8
             fna_panom(:,iy,0,isb)=0.0_R8
             fna_ch(:,iy,0,isb)=0.0_R8
             fna_pschused(:,iy,0,isb)=0.0_R8
          endif
        enddo                                                           !srv 20.06.08 }
      endif                                                             !xpb

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
            fnb_53(ix,iy,0)=0.0_R8                                      !ajs 06.04.20
            fnb_mdf(ix,iy,0)=0.0_R8
            flob(ix,iy,0)=0.0_R8
            vbecrb(ix,iy,0)=0.0_R8
            vbdia(ix,iy,0)=0.0_R8
            wbdia(ix,iy,0)=0.0_R8
            ubdia(ix,iy,0)=0.0_R8
            conb(ix,iy,0,:)=0.0_R8
            fllim0fnb(ix,iy,0)=1.0_R8
            if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
               fna_pinch(ix,iy,0,isb)=0.0_R8
               fna_pll(ix,iy,0,isb)=0.0_R8
               fna_drift(ix,iy,0,isb)=0.0_R8
               fna_nanom(ix,iy,0,isb)=0.0_R8
               fna_panom(ix,iy,0,isb)=0.0_R8
               fna_ch(ix,iy,0,isb)=0.0_R8
               fna_pschused(ix,iy,0,isb)=0.0_R8
            endif
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
            fnb_53(ix,iy,1)=0.0_R8                                      !ajs 06.04.20
            fnb_mdf(ix,iy,1)=0.0_R8
            flob(ix,iy,1)=0.0_R8
            vbecrb(ix,iy,1)=0.0_R8
            vbdia(ix,iy,1)=0.0_R8
            wbdia(ix,iy,1)=0.0_R8
            ubdia(ix,iy,1)=0.0_R8
            conb(ix,iy,1,:)=0.0_R8
            fllim0fnb(ix,iy,1)=1.0_R8
            if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
               fna_pinch(ix,iy,1,isb)=0.0_R8
               fna_pll(ix,iy,1,isb)=0.0_R8
               fna_drift(ix,iy,1,isb)=0.0_R8
               fna_nanom(ix,iy,1,isb)=0.0_R8
               fna_panom(ix,iy,1,isb)=0.0_R8
               fna_ch(ix,iy,1,isb)=0.0_R8
               fna_pschused(ix,iy,1,isb)=0.0_R8
            endif
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
            fnb_fcor(ix,iy,0)=0.0_R8
            fnb_he(ix,iy,0)=0.0_R8
            fnb_53(ix,iy,0)=0.0_R8
            fnb_mdf(ix,iy,0)=0.0_R8
            flob(ix,iy,0)=0.0_R8
            conb(ix,iy,0,:)=0.0_R8
            fllim0fnb(ix,iy,0)=1.0_R8
            fnb(ix,iy,1)=0.0_R8
            fnb_32(ix,iy,1)=0.0_R8
            fnb_52(ix,iy,1)=0.0_R8
            fnb_nodrift(ix,iy,1)=0.0_R8
            fnb_32nd(ix,iy,1)=0.0_R8
            fnb_52nd(ix,iy,1)=0.0_R8
            fnb_fcor(ix,iy,1)=0.0_R8
            fnb_he(ix,iy,1)=0.0_R8
            fnb_53(ix,iy,1)=0.0_R8
            fnb_mdf(ix,iy,1)=0.0_R8
            flob(ix,iy,1)=0.0_R8
            conb(ix,iy,1,:)=0.0_R8
            fllim0fnb(ix,iy,1) = 1.0_R8
            vbecrb(ix,iy,1)=0.0_R8
            vbdia(ix,iy,1)=0.0_R8
            wbdia(ix,iy,1)=0.0_R8
            ubdia(ix,iy,1)=0.0_R8
            if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
               fna_pinch(ix,iy,0:1,isb)=0.0_R8
               fna_pll(ix,iy,0:1,isb)=0.0_R8
               fna_drift(ix,iy,0:1,isb)=0.0_R8
               fna_nanom(ix,iy,0:1,isb)=0.0_R8
               fna_panom(ix,iy,0:1,isb)=0.0_R8
               fna_ch(ix,iy,0:1,isb)=0.0_R8
               fna_pschused(ix,iy,0:1,isb)=0.0_R8
            endif
          end do
        end do
      end do

* ..compute kinetic energy of species
      do ix=-1,nx
       do iy=-1,ny
         kbnrgy(ix,iy) = 0.5_R8*mb*ub(ix,iy)**2
       end do
      end do

*   ..flux correction
*   ..we rescale the fnb fluxes to have the same integral over the boundaries
*   ..as the fnb_mdf fluxes
      if(pflux_cor.ne.0 .and. mdf_fnb.ne.0) then                        !srv 05.03.16 {
        if(.not. is_neutral(isb)) then
          call b2tbfl(nx, ny,.false.,fnb, tbfl_fnb, div_fnb)
          call b2tbfl(nx, ny,.false.,fnb_mdf, tbfl_fnb_mdf, div_fnb_mdf)
          if(div_fnb .ne. 0.0_R8) then
            rdiv= div_fnb_mdf/div_fnb
            fnb = rdiv*fnb
          else
            if (ncall.eq.0) write(*,"(a,i4)")
     >       ' Boundary flux correction not attempted for is = ',isb
          endif
          if(iout.ne.0) then
            write(*,*)
            write(*,*) 'Before flux correction in b2tfnb, is = ', isb
            write(*,'(a,i4,a,es18.11,a,es18.11)')
     $       ' isb = ',isb,' fnb_tot = ',div_fnb,
     $       ' fnb_mdf_tot = ',div_fnb_mdf

            call b2tbfl(nx, ny,.false., fnb, tbfl_fnb, div_fnb)
            call b2tbfl(nx, ny,.false., fnb_mdf, tbfl_fnb_mdf,
     .                                                 div_fnb_mdf)
            write(*,*) 'After flux correction in b2tfnb, is = ', isb
            write(*,'(a,i4,a,es18.11,a,es18.11)')
     $       ' isb = ',isb,' fnb_tot = ',div_fnb,
     $       ' fnb_mdf_tot = ',div_fnb_mdf
            write(*,*)
!            write(*,'(a,i3,a,es15.7,a,es15.7,a,es15.7)')
!     &       ' is =',isb,' rdiv =',rdiv,' div_fnb =',div_fnb,
!     &       ' div_fnb_mdf =',div_fnb_mdf
          endif
        endif
      endif                                                             !srv 05.03.16 }

      if (iout_b2wdat.eq.4 .and. balance) then
       call my_out(70,nx,ny,fnb(-1,-1,0),'b2tfnb_fnbx'//chns)
       call my_out(70,nx,ny,fnb(-1,-1,1),'b2tfnb_fnby'//chns)
       call my_out(70,nx,ny,fnb_fcor(-1,-1,0),
     &            'b2tfnb_fnb_fcorx'//chns)
       call my_out(70,nx,ny,fnb_fcor(-1,-1,1),
     &            'b2tfnb_fnb_fcory'//chns)
       call my_out(70,nx,ny,fnb_he(-1,-1,0),'b2tfnb_fnb_hex'//chns)
       call my_out(70,nx,ny,fnb_he(-1,-1,1),'b2tfnb_fnb_hey'//chns)
       call my_out(70,nx,ny,fnb_53(-1,-1,0),'b2tfnb_fnb_53x'//chns)
       call my_out(70,nx,ny,fnb_53(-1,-1,1),'b2tfnb_fnb_53y'//chns)
       call my_out(70,nx,ny,fnbPSch(-1,-1,0),'b2tfnb_fnbPSchx'//chns)
       call my_out(70,nx,ny,fnbPSch(-1,-1,1),'b2tfnb_fnbPSchy'//chns)
!      if (isb .eq. ismain) then
          call my_out(70,nx,ny,vbecrb(-1,-1,0),'b2tfnb_vbecrbx'//chns)
          call my_out(70,nx,ny,vbecrb(-1,-1,1),'b2tfnb_vbecrby'//chns)
!      endif
!       call my_out(70,nx,ny,ubdia(-1,-1,0),'b2tfnb_ubdiax'//chns)
!       call my_out(70,nx,ny,ubdia(-1,-1,1),'b2tfnb_ubdiay'//chns)
       call my_out(70,nx,ny,vbdia(-1,-1,0),'b2tfnb_vbdiax'//chns)
       call my_out(70,nx,ny,vbdia(-1,-1,1),'b2tfnb_vbdiay'//chns)
       call my_out(70,nx,ny,wbdia(-1,-1,0),'b2tfnb_wbdiax'//chns)
       call my_out(70,nx,ny,wbdia(-1,-1,1),'b2tfnb_wbdiay'//chns)
       call my_out(70,nx,ny,kbnrgy,'b2tfnb_kbnrgy'//chns)
!       call my_out(70,nx,ny,fchvispar_b(-1,-1,0),
!     &                                 'b2tfnb_fchvispar_bx'//chns)
!       call my_out(70,nx,ny,fchvispar_b(-1,-1,1),
!     &                                 'b2tfnb_fchvispar_by'//chns)
!       call my_out(70,nx,ny,fchvisq_b(-1,-1,0),
!     &                                  'b2tfnb_fchvisq_bx'//chns)
!       call my_out(70,nx,ny,fchvisq_b(-1,-1,1),
!     &                                  'b2tfnb_fchvisq_by'//chns)
!       call my_out(70,nx,ny,fchvisper_b(-1,-1,0),
!     &                                  'b2tfnb_fchvisper_bx'//chns)
!       call my_out(70,nx,ny,fchvisper_b(-1,-1,1),
!     &                                  'b2tfnb_fchvisper_by'//chns)
!       call my_out(70,nx,ny,fchinert_b(-1,-1,0),
!     &                                 'b2tfnb_fchinert_bx'//chns)
!       call my_out(70,nx,ny,fchinert_b(-1,-1,1),
!     &                                 'b2tfnb_fchinert_by'//chns)
!       call my_out(70,nx,ny,fchin_b(-1,-1,0),
!     &                                 'b2tfnb_fchin_bx'//chns)
!       call my_out(70,nx,ny,fchin_b(-1,-1,1),
!     &                                 'b2tfnb_fchin_by'//chns)
!       call my_out(70,nx,ny,fchanml_b(-1,-1,0),
!     &                                 'b2tfnb_fchanml_bx'//chns)
!       call my_out(70,nx,ny,fchanml_b(-1,-1,1),
!     &                                 'b2tfnb_fchanml_by'//chns)
      endif

      if(iout.ne.0) then                                                !srv 18.06.08 {
       write (chns,'(i3.3)') isb
       call my_out(70,nx,ny,fnb(-1,-1,0),'b2tfnb_fnbx'//chns)
       call my_out(70,nx,ny,fnb(-1,-1,1),'b2tfnb_fnby'//chns)
       call my_out(70,nx,ny,fnb_nodrift(-1,-1,0),
     &             'b2tfnb_fnb_nodriftx'//chns)
       call my_out(70,nx,ny,fnb_nodrift(-1,-1,1),
     &             'b2tfnb_fnb_nodrifty'//chns)
       call my_out(70,nx,ny,fnb_fcor(-1,-1,0),
     &             'b2tfnb_fnb_fcorx'//chns)
       call my_out(70,nx,ny,fnb_fcor(-1,-1,1),
     &             'b2tfnb_fnb_fcory'//chns)
       call my_out(70,nx,ny,fnb_he(-1,-1,0),'b2tfnb_fnb_hex'//chns)
       call my_out(70,nx,ny,fnb_he(-1,-1,1),'b2tfnb_fnb_hey'//chns)
       call my_out(70,nx,ny,fnb_53(-1,-1,0),'b2tfnb_fnb_53x'//chns)
       call my_out(70,nx,ny,fnb_53(-1,-1,1),'b2tfnb_fnb_53y'//chns)
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
       call my_out(70,nx,ny,term,'b2tfnb_term'//chns)
       call my_out(70,nx,ny,termhz,'b2tfnb_termhz'//chns)
       call my_out(70,nx,ny,pbs(-1,-1,0),'b2tfnb_pbsx')
       call my_out(70,nx,ny,pbshz(-1,-1,0),'b2tfnb_pbshzx')
       call my_out(70,nx,ny,pbshz(-1,-1,1),'b2tfnb_pbshzy')
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
       call my_out(70,nx,ny,pb,'b2tfnb_pa'//chns)
       call my_out(70,nx,ny,te,'b2tfnb_te')
       call my_out(70,nx,ny,ti,'b2tfnb_ti')
       call my_out(70,nx,ny,po,'b2tfnb_po')
       call my_out(70,nx,ny,zb,'b2tfnb_zb'//chns)
       call my_out(70,nx,ny,facdrift,'b2tfnb_facdrift'//chns)
       call my_out(70,nx,ny,fac_ExB,'b2tfnb_fac_ExB'//chns)
       call my_out(70,nx,ny,fchvispar_b(-1,-1,0),
     &                           'b2tfnb_fchvispar_bx'//chns)
       call my_out(70,nx,ny,fchvispar_b(-1,-1,1),
     &                           'b2tfnb_fchvispar_by'//chns)
       call my_out(70,nx,ny,fchvisq_b(-1,-1,0),
     &                           'b2tfnb_fchvisq_bx'//chns)
       call my_out(70,nx,ny,fchvisq_b(-1,-1,1),
     &                           'b2tfnb_fchvisq_by'//chns)
       call my_out(70,nx,ny,fchvisper_b(:,:,0),
     &                           'b2tfnb_fchvisper_bx'//chns)
       call my_out(70,nx,ny,fchvisper_b(:,:,1),
     &                           'b2tfnb_fchvisper_by'//chns)
       call my_out(70,nx,ny,fchinert_b(-1,-1,0),
     &                           'b2tfnb_fchinert_bx'//chns)
       call my_out(70,nx,ny,fchinert_b(-1,-1,1),
     &                           'b2tfnb_fchinert_by'//chns)
       call my_out(70,nx,ny,fchin_b(-1,-1,0),
     &                           'b2tfnb_fchin_bx'//chns)
       call my_out(70,nx,ny,fchin_b(-1,-1,1),
     &                           'b2tfnb_fchin_by'//chns)
       call my_out(70,nx,ny,fchanml_b(-1,-1,0),
     &                           'b2tfnb_fchanml_bx'//chns)
       call my_out(70,nx,ny,fchanml_b(-1,-1,1),
     &                           'b2tfnb_fchanml_by'//chns)
      endif                                                             !srv 18.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfnb

      contains

      function kap(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: kap
      kap = 1.0_R8 / (zb(ix,iy)*te(ix,iy) + ti(ix,iy))
      return
      end function kap

      function kmpr(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: kmpr
      kmpr = nb(ix,iy)/pb(ix,iy)
      return
      end function kmpr

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
