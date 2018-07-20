*-----------------------------------------------------------------------
*.specification

      subroutine b2tfnb (nCv, nFc, nVx, isb, ismain, geo, mpg, pl, dv,
     &  co, rt, balance)
      use b2mod_types
!      use b2mod_boundary_namelist
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2us_geo
      use b2us_map
      use b2us_plasma
!      use b2mod_balance                                                  !djm Jan2017
!     & , only : fna_pinch, fna_pll, fna_drift, fna_nanom, fna_panom,
!     &          fna_ch, fna_pschused, balance_netcdf
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, isb, ismain
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type(B2Plasma), intent (in) :: pl 
      type(B2Rates), intent (in) :: rt
      logical balance                                                   !djm Jan2017
*   ..input/output arguments
      type(B2Coeff), intent (inout) :: co
      type(B2Derivatives), intent (inout) :: dv
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
      integer ncall, iCv, iFc, ib, no_hybr,
     &        ixl, iyl, ixb, iyb                                         !srv 13.10.06
      integer, save :: iout = 0, iout_b2wdat = 0,                        !srv 16.09.02
     &                 iout_b2npc9 = 0, iout_b2npco = 0, iout_b2npmo =0
      integer, save :: mdf_fnb = 0                                       !srv 13.10.06
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0              !srv 06.07.06
      integer, save :: pflux_cor = 1                                     !srv 05.03.16
      integer, save :: PSch_style = 1                                    !srv 12.11.08
      integer, save :: ion_vlct_restrict = 0                             !iyv 06.02.14
      integer, save :: vlct_diagno = 1
      integer, save :: meth = 2
      character chns*3                                                   !srv 16.09.02
      integer, save :: drift_style = 1                                   !srv 21.06.08
      real (kind=R8), save :: alpha_ = 0.0_R8                            !srv 09.07.03
      real (kind=R8), save :: xcur = 0.0_R8                              !xpb 27.09.17
      real (kind=R8), save :: ycur = 1.0_R8                              !srv 07.07.08
      real (kind=R8), save :: vis_per = 0.0_R8, vis_q = 1.0_R8           !srv 02.10.17
      real (kind=R8), save :: xvecrb = 1.0_R8, xwdia = 1.0_R8            !srv 27.01.12
      real (kind=R8), save ::  
     &   poleldr = 1.0_R8, radeldr = 1.0_R8,                             !srv 25.04.00
     &   poldidr = 1.0_R8, raddidr = 1.0_R8,                             !srv 25.04.00
     &   PSch = 1.0_R8                                                   !srv 17.01.07
      real (kind=R8) ::
     *   flo(nFc,0:1), flo_mdf(nFc,0:1), flo_nodrift(nFc,0:1), 
     *   flo2dia(nFc,0:1), flo_he(nFc,0:1), 
     *   cutlo, alpha, gamma, vbar, t0, t1, t2, t3, maxfluxlimit(0:1),
     *   flux_limit_min_ti, xfrhie, term(nFc), termhz(nFc),
     *   flob0(nFc,0:1), conb0(nFc,0:1), scur(nFc,0:1), fl
      real (kind=R8) :: 
     &   vbecrb_(nCv,0:1), wbdia_(nCv,0:1),
     &   kmprti(nFc), dkmprti(nFc,0:1),              !srv 21.10.10
     &   kmpr(nCv),kmprf(nFc),dkmpr(nFc,0:1)!,
!     &   tbfl_fnb(1:nbc), div_fnb, tbfl_fnb_mdf(1:nbc), div_fnb_mdf,     !srv 05.03.16
!     &   rdiv                                                            !srv 05.03.16
      parameter (cutlo=1.0e-60_R8)
      real (kind=R8) ::
     *   wrk0(nFc), wrk1(nFc), wrk2(nFc), wrk3(nFc), wrk4(nFc), 
     *   wrk5(nFc),wrk6(nFc), wrk7(nFc), wrk8(nFc),
     *   wrk(nCv), wrkVx(nVx), wpo(nCv), wpb(nCv),
     *   dumm0(nFc,0:1), dumm1(nFc,0:1)
      real (kind=R8) ::
     *   poc(nVx), ub2dia(nFc,0:1), 
     *   nbcor(nVx),                                             !srv 01.10.99 !lk 09.06.11
     *   nbti(nCv), nbtic(nVx),                          !srv 08.10.04
     *   nbzb(nCv), nbzbc(nVx),                          !srv 13.05.16
     *   dnb(nFc,0:1), dpb(nFc,0:1),
     *   dpo(nFc,0:1), dnbti(nFc,0:1),
     *   nbf(nFc), pbf(nFc),
     *   dOnedBsq(nFc,0:1),
     *   fna_mdf0(nFc,0:1)
      real (kind=R8) :: weight(nFc,2)                              !lk
      integer :: xcount, ycount
      save ncall, alpha, gamma, flux_limit_min_ti, no_hybr
*   ..procedures
      external subini, subend, xertst, ipgeti, ipgetr                    !xpb
      external b2xvsg, b2xvff, intface, diff
      intrinsic abs, max, min, mod, sqrt, sign
*   ..initialisation
      data ncall/0/, xfrhie/1.0_R8/, alpha/0.0_R8/, gamma/2.0_R8/, 
     *     flux_limit_min_ti/0.0_R8/, no_hybr/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfnb')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
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
        call ipgeti ('b2tfnb_no_hybr', no_hybr)
          if (no_hybr .eq.1) meth = 0
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
        call ipgetr ('b2tfnb_PSch', PSch)                                !srv 17.01.07
        call ipgeti ('b2tfnb_pflux_cor', pflux_cor)                      !srv 05.03.16
        call ipgeti ('b2tfnb_PSch_style', PSch_style)                    !srv 12.11.08
        call ipgetr ('b2tfnb_xcur', xcur)                                !xpb 27.09.17
        call ipgetr ('b2tfnb_ycur', ycur)                                !srv 17.07.08
        call ipgetr ('b2tfnb_vis_per', vis_per)                          !srv 02.10.17
        call ipgetr ('b2tfnb_vis_q', vis_q)                              !srv 02.10.17
        call ipgeti ('b2wdat_iout', iout_b2wdat)                         !iys 17.04.17
        call ipgeti ('b2npc9_iout', iout_b2npc9)                         !srv 11.09.09
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgetr ('b2tfrn_xvecrb', xvecrb)                            !srv 27.01.12
        call ipgetr ('b2tfrn_xwdia', xwdia)                              !srv 27.01.12
        call ipgeti ('b2npmo_ion_vlct_restrict', ion_vlct_restrict)      !iyv 06.02.14
        call ipgeti ('b2npmo_vlct_diagno', vlct_diagno)
        call ipgeti ('b2tfnb_discr_meth', meth)
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of nb, ne, pb, cdnb, cdpb
       call b2xvsg (nCv, pl % na(:,isb), 1, 'nb', '.gt.')
       call b2xvsg (nCv, dv % ne       , 1, 'ne', '.gt.')
       call b2xvsg (nCv, dv % pa(:,isb), 1, 'pb', '.gt.')
       call b2xvsg (nCv*2,co%cdna(1,0,isb)*geo%fcQalf, 1, 'cdna','.ge.')
       call b2xvsg (nCv*2,co%cdpa(1,0,isb)*geo%fcQalf, 1, 'cdpa','.ge.')
      endif

*   ..precompute kmpr and its differences
      kmpr = pl % na(:,isb) / dv % pa(:,isb)
      call diff(nCv, nFc, nVx, 0, geo, mpg, kmpr, wrkVx, dkmpr)

      write(chns,'(i3.3)') isb
      if(mdf_fnb .ne. 0 .and.
     &  (maxval(dv % facdrift)   .ne. 0.0_R8  .or. 
     &   maxval(dv % fac_ExB )   .ne. 0.0_R8) .and.
     &   minval(rt % rza(:,isb)) .gt. 0.0_R8) then
*   ..compute kmpr*ti on the left face cells and the bottom cells
        wrk = kmpr*pl%ti
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,wrk,kmprti)
*   ..compute d(kmpr*ti) at the left face cells
        call diff(nCv, nFc, nVx, 0, geo, mpg, wrk, wrkVx, dkmprti)
      else
        kmprti  = 0.0_R8
        dkmprti = 0.0_R8
      endif                                                              !srv 21.10.10 }

c Drift terms
      maxfluxlimit        = 0.0_R8
      dv%uadia(:,:,isb)   = 0.0_R8
      dv%vadia(:,:,isb)   = 0.0_R8
      dv%wadia(:,:,isb)   = 0.0_R8
      dv%vaecrb(:,:,isb)  = 0.0_R8
      dv%fnaPSch(:,:,isb) = 0.0_R8                                                 !srv 13.10.06
      dv%fna_mdf(:,:,isb) = 0.0_R8                                                 !srv 19.11.09
      dv%flob  = 0.0_R8
      dv%conb  = 0.0_R8                                                   !srv 13.10.06
      ub2dia=0.0_R8                                                      !srv 01.10.99
#ifdef WG_DRIFTS
      if((maxval(facdrift).ne.0.0_R8 .or. maxval(fac_ExB).ne.0.0_R8)
     &    .and. minval(zb).gt.0.0_R8) then
       if(drift_style.eq.0) then                                         !srv 21.06.08
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
cwdk correct computation of drift components to always represent the poloidal
cwdk projection of the drift
              vbecrb(ix,iy,0,0) = - fac_ExB(ix,iy)*poleldr/hy1(ix,iy)*   !xpb
     1         bb(ix,iy,2)/bb(ix,iy,3)**2*                               !xpb
     2         (wpo(topix(ix,iy),topiy(ix,iy))-wpo(ix,iy))
              vbdia(ix,iy,0,0) = - bb(ix,iy,2)/bb(ix,iy,3)/hy1(ix,iy)/   !xpb
     1         bb(ix,iy,3)*2.0_R8*ti(ix,iy)/bb(ix,iy,3)/qe/zb(ix,iy)*    !xpb
     3         (wbbv(topix(ix,iy),topiy(ix,iy),3)-wbbv(ix,iy,3))*
     &         facdrift(ix,iy)*poldidr                                   !xpb
              wbdia(ix,iy,0,0) = - bb(ix,iy,2)/bb(ix,iy,3)/hy1(ix,iy)
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
              vbecrb(ix,iy,0,0) = 0.0_R8                                 !xpb
              vbdia(ix,iy,0,0) = 0.0_R8                                  !xpb
              wbdia(ix,iy,0,0) = 0.0_R8
            endif                                                        !xpb
            ubdia(ix,iy,0,0) = vbecrb(ix,iy,0,0) + vbdia(ix,iy,0,0)      !srv 17.05.99
            ub2dia(ix,iy,0,0) = vbecrb(ix,iy,0,0) + 
     &        2.0_R8*vbdia(ix,iy,0,0)                                    !srv 01.10.99
          enddo
        enddo
        
        call intfaceh(nx,ny,leftix,leftiy,vol,po,wpo)
        call intfaceh(nx,ny,leftix,leftiy,vol,nb*ti,wpb)
        do iy = -1, ny                                                   !xpb
          do ix = -1, nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
             vbecrb(ix,iy,1,1) = fac_ExB(ix,iy)*radeldr/hx(ix,iy)*       !xpb
     1          bb(ix,iy,2)/bb(ix,iy,3)**2*                              !xpb
     2         (wpo(rightix(ix,iy),rightiy(ix,iy))-wpo(ix,iy))
             vbdia(ix,iy,1,1) = facdrift(ix,iy)*raddidr* 
     &         ti(ix,iy)*bb(ix,iy,2)/(qe*zb(ix,iy)*hx(ix,iy))            !xpb
     &         *(2.0_R8/bb(ix,iy,3)**3)                                  !xpb
     &         *(wbbr(ix,iy,3)-wbbl(ix,iy,3))                            !xpb
             wbdia(ix,iy,1,1) = bb(ix,iy,2)
     &         /(qe*zb(ix,iy)*hx(ix,iy)*nb(ix,iy))                       !xpb
     &         /bb(ix,iy,3)**2*                                          !xpb
     2         (wpb(rightix(ix,iy),rightiy(ix,iy))-wpb(ix,iy))*
     &         facdrift(ix,iy)
             ubdia(ix,iy,1,1) = vbecrb(ix,iy,1,1) + vbdia(ix,iy,1,1)     !srv 17.05.99
             ub2dia(ix,iy,1,1) = vbecrb(ix,iy,1,1) + 
     &         2.0_R8*vbdia(ix,iy,1,1)                                   !srv 01.10.99
            endif
          enddo
        enddo

* drift velocities at walls are extrapolated from those at neighbouring cells
* drift velocities at iy=ny are extrapolated from those at iy=ny-1       !xpb
* drift velocities at iy=-1 are extrapolated from those at iy=0          !xpb
        do ix = -1, nx                                                   !xpb
          do iy = -1, ny                                                 !xpb
            if(leftix(ix,iy).eq.-2) then
              vbecrb(ix,iy,:,:) = 
     &          vbecrb(rightix(ix,iy),rightiy(ix,iy),:,:)                !xpb
              vbdia(ix,iy,:,:) = 
     &          vbdia(rightix(ix,iy),rightiy(ix,iy),:,:)                 !xpb
              wbdia(ix,iy,:,:) = 
     &          wbdia(rightix(ix,iy),rightiy(ix,iy),:,:)                 !xpb
              ubdia(ix,iy,:,:) = vbecrb(ix,iy,:,:) + vbdia(ix,iy,:,:) 
              ub2dia(ix,iy,:,:) = 
     &          vbecrb(ix,iy,:,:) + 2.0_R8*vbdia(ix,iy,:,:)  !srv 17.07.08
            elseif(rightix(ix,iy).eq.nx+1) then                          !xpb
              vbecrb(ix,iy,:,:) = 
     &          vbecrb(leftix(ix,iy),leftiy(ix,iy),:,:)                  !xpb
              vbdia(ix,iy,:,:) = vbdia(leftix(ix,iy),leftiy(ix,iy),:,:)  !xpb
              wbdia(ix,iy,:,:) = wbdia(leftix(ix,iy),leftiy(ix,iy),:,:)  !xpb
              ubdia(ix,iy,:,:) = vbecrb(ix,iy,:,:) + vbdia(ix,iy,:,:) 
              ub2dia(ix,iy,:,:) = 
     &          vbecrb(ix,iy,:,:) + 2.0_R8*vbdia(ix,iy,:,:)              !srv 17.07.08
            endif
            if(topiy(ix,iy).eq.ny+1) then                                !xpb
              vbecrb(ix,iy,:,:) = 
     .          vbecrb(bottomix(ix,iy),bottomiy(ix,iy),:,:)              !xpb
              vbdia(ix,iy,:,:) = 
     &          vbdia(bottomix(ix,iy),bottomiy(ix,iy),:,:)               !xpb
              wbdia(ix,iy,:,:) = 
     &          wbdia(bottomix(ix,iy),bottomiy(ix,iy),:,:)               !xpb
              ubdia(ix,iy,:,:) = vbecrb(ix,iy,:,:) + vbdia(ix,iy,:,:) 
              ub2dia(ix,iy,:,:) = 
     &          vbecrb(ix,iy,:,:) + 2.0_R8*vbdia(ix,iy,:,:)              !srv 17.07.08
            elseif(bottomiy(ix,iy).eq.-2) then
              vbecrb(ix,iy,:,:) = vbecrb(topix(ix,iy),topiy(ix,iy),:,:)  !xpb
              vbdia(ix,iy,:,:) = vbdia(topix(ix,iy),topiy(ix,iy),:,:)    !xpb
              wbdia(ix,iy,:,:) = wbdia(topix(ix,iy),topiy(ix,iy),:,:)    !xpb
              ubdia(ix,iy,:,:) = vbecrb(ix,iy,:,:) + vbdia(ix,iy,:,:) 
              ub2dia(ix,iy,:,:) = 
     &          vbecrb(ix,iy,:,:) + 2.0_R8*vbdia(ix,iy,:,:)              !srv 17.07.08
            endif
          enddo
        enddo

* turn off side drift flows in guard cells                               !xpb
        do iy = -1, ny                                                   !xpb
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
           vbecrb(ix,iy,1,:) = 0.0_R8                                    !xpb
           vbdia(ix,iy,1,:) = 0.0_R8                                     !xpb
           wbdia(ix,iy,1,:) = 0.0_R8                                     !xpb
           ubdia(ix,iy,1,:) = 0.0_R8                                     !xpb
           ub2dia(ix,iy,1,:) = 0.0_R8                                    !srv 17.07.08
          endif
          if(bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1) then
           vbecrb(ix,iy,0,:) = 0.0_R8                                    !xpb
           vbdia(ix,iy,0,:) = 0.0_R8                                     !xpb
           wbdia(ix,iy,0,:) = 0.0_R8                                     !xpb
           ubdia(ix,iy,0,:) = 0.0_R8                                     !xpb
           ub2dia(ix,iy,0,:) = 0.0_R8                                    !srv 17.07.08
          endif
         enddo
        enddo                                                            !xpb
       else                                                              !srv 21.06.08 {
        weight = 1.0e0_R8                                                !lk
cwdk add computation of 9-point corrections for drifts
cwdk precomputation of differences on faces
!lk       call intcor(nx,ny,vol,po,poc)
!        call intcor(nx,ny,weight,po,poc)                                 !lk 12.12.07
        nbti=nb*ti
!lk       call intcor(nx,ny,vol,nbti,nbtic)                              !srv 08.10.04
!        call intcor(nx,ny,weight,nbti,nbtic)                             !lk 12.12.07
        nbzb = nb*zb                                                     !srv 13.05.16
        call intcor(nx,ny,weight,nbzb,nbzbc)                             !srv 13.05.16
        call diff(nx, ny, hx, hy, weight, qs, po,   dpo)
        call diff(nx, ny, hx, hy, weight, qs, nbti, dnbti)
        call diff(nx, ny, hx, hy, vol, qs, OnedBsq, dOnedBsq)
!        call diff_p(nx,ny, 1, hx, hy, weight, qs, OnedBsq, OnedBsqc, 
!     .    dOnedBsq(-1,-1,0,0))
!        call diff_r(nx,ny, 1, hx, hy, weight, qs, OnedBsq, OnedBsqc, 
!     .    dOnedBsq(-1,-1,0,1))
        xcount = 0
        do iy = 0, ny-1
         do ix = 0, nx
          if(leftix(ix,iy).ne.-2) then                                   !srv 26.11.02
             vbecrb(ix,iy,0,0) = - poleldr*fac_ExB(ix,iy)/               !srv 17.05.99 25.08.99 07.06.00 05.02.05 13.05.16 {
     1      (hy1(ix,iy)+hy1(leftix(ix,iy),leftiy(ix,iy)))* ! signmf*     !srv
     2      (bb(ix,iy,2)/bb(ix,iy,3)**2                                  !srv 05.02.05
     2         +bb(leftix(ix,iy),leftiy(ix,iy),2)/
     2          bb(leftix(ix,iy),leftiy(ix,iy),3)**2)*                   !srv 05.02.05
     2       dpo(ix,iy,0,1)                                              !srv 13.05.16 }
cwdk computation of radial drifts on poloidal faces
cwdk somewhat messy interpolation of geometric factors, for exact
cwdk cancelation of some terms in poloidal and radial drift components
             vbecrb(ix,iy,0,1) =   radeldr*fac_ExB(ix,iy)/
     1      (hy1(ix,iy)+hy1(leftix(ix,iy),leftiy(ix,iy)))*
     &      (hy(ix,iy) + hy(leftix(ix,iy),leftiy(ix,iy)))/
     &      (hx(ix,iy) + hx(leftix(ix,iy),leftiy(ix,iy)))*
     &       qc(ix,iy,0)*
     2      (bb(ix,iy,2)/bb(ix,iy,3)**2
     2         +bb(leftix(ix,iy),leftiy(ix,iy),2)/
     2          bb(leftix(ix,iy),leftiy(ix,iy),3)**2)*
     2       dpo(ix,iy,0,0)
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
             if( abs(vbecrb(ix,iy,0,0)) .gt. 2.0_R8*cssb(ix,iy,0) ) then
               xcount = xcount+1
               if (vlct_diagno.ge.2) then
                 write(6,'(a,a,i3,a,i4,a,i4,a)')
     .            'Poloidal ExB velocity restriction applied ',
     .            'to species ',isb,' in cell (',ix,',',iy,'):'
                 write(6,'(a,1p,1e12.5,a,1e12.5)')
     .            'Mach 2: ',2.0_R8*cssb(ix,iy,0),
     .            ' Vx_ExB:',vbecrb(ix,iy,0,0)
               end if
               vbecrb(ix,iy,0,0)=
     .          sign(2.0_R8*cssb(ix,iy,0),vbecrb(ix,iy,0,0))
             endif
             if( abs(vbecrb(ix,iy,0,1)) .gt. 2.0_R8*cssb(ix,iy,0) ) then
               if (vlct_diagno.ge.2) then
                 write(6,'(a,a,i3,a,i4,a,i4,a)')
     .            'Radial ExB velocity restriction applied ',
     .            'to species ',isb,', cell (',ix,',',iy,'), left face:'
                 write(6,'(a,1p,1e12.5,a,1e12.5)')
     .            'Mach 2: ',2.0_R8*cssb(ix,iy,0),
     .            ' Vy_ExB:',vbecrb(ix,iy,0,1)
               end if
               vbecrb(ix,iy,0,1)=
     .          sign(2.0_R8*cssb(ix,iy,0),vbecrb(ix,iy,0,1))
             endif
           endif                                                         !iyv 06.02.14  }
           vbdia(ix,iy,0,0) =    poldidr*facdrift(ix,iy)* ! signmf*      !srv 07.06.00 17.05.99 25.04.00 05.02.05 11.01.13
     &          (ti(leftix(ix,iy),leftiy(ix,iy))                         !srv 05.03.99
     &           *bb(leftix(ix,iy),leftiy(ix,iy),2)                      !srv 05.02.05
     &           /(qe*zb(ix,iy)*hy1(leftix(ix,iy),leftiy(ix,iy)))        !srv
     &          +ti(ix,iy)*bb(ix,iy,2)/(qe*zb(ix,iy)*hy1(ix,iy)))*0.5_R8 !srv 05.02.05
     &         *dOnedBsq(ix,iy,0,1)                                      !srv 05.03.99
cwdk computation of radial drifts on poloidal faces
cwdk somewhat messy interpolation of geometric factors, for exact
cwdk cancelation of some terms in poloidal and radial drift components
           vbdia(ix,iy,0,1) = -  raddidr*facdrift(ix,iy)* ! signmf*
     &          (ti(leftix(ix,iy),leftiy(ix,iy))
     &           *bb(leftix(ix,iy),leftiy(ix,iy),2)
     &           /(qe*zb(ix,iy)*hy1(leftix(ix,iy),leftiy(ix,iy))) 
     &         +ti(ix,iy)*bb(ix,iy,2)/(qe*zb(ix,iy)*hy1(ix,iy)))*0.5_R8*
     &      (hy(ix,iy) + hy(leftix(ix,iy),leftiy(ix,iy)))/
     &      (hx(ix,iy) + hx(leftix(ix,iy),leftiy(ix,iy)))*
     &       qc(ix,iy,0)
     &         *dOnedBsq(ix,iy,0,0) 
           wbdia(ix,iy,0,0) = - facdrift(ix,iy)/                         !srv 13.05.16 {
     1	    (hy1(ix,iy)+hy1(leftix(ix,iy),leftiy(ix,iy)))*
     2	    (bb(ix,iy,2)/bb(ix,iy,3)**2
     2        +bb(leftix(ix,iy),leftiy(ix,iy),2)/
     2         bb(leftix(ix,iy),leftiy(ix,iy),3)**2)*
     2         dnbti(ix,iy,0,1)/
     &        (nbzbc(topix(ix,iy),topiy(ix,iy))+nbzbc(ix,iy))*2.0_R8/qe  !srv 13.05.16 }
cwdk computation of radial drifts on poloidal faces
cwdk somewhat messy interpolation of geometric factors, for exact
cwdk cancelation of some terms in poloidal and radial drift components
           wbdia(ix,iy,0,1) =   facdrift(ix,iy)/
     1	    (hy1(ix,iy)+hy1(leftix(ix,iy),leftiy(ix,iy)))*
     &      (hy(ix,iy) + hy(leftix(ix,iy),leftiy(ix,iy)))/
     &      (hx(ix,iy) + hx(leftix(ix,iy),leftiy(ix,iy)))*
     &       qc(ix,iy,0)*
     2	    (bb(ix,iy,2)/bb(ix,iy,3)**2
     2        +bb(leftix(ix,iy),leftiy(ix,iy),2)/
     2         bb(leftix(ix,iy),leftiy(ix,iy),3)**2)*
     2         dnbti(ix,iy,0,0)/
     &        (nbzbc(topix(ix,iy),topiy(ix,iy))+nbzbc(ix,iy))*2.0_R8/qe
           ubdia(ix,iy,0,0) = vbecrb(ix,iy,0,0) + vbdia(ix,iy,0,0)       !srv 17.05.99
           ubdia(ix,iy,0,1) = vbecrb(ix,iy,0,1) + vbdia(ix,iy,0,1)
           ub2dia(ix,iy,0,0) = 
     &       vbecrb(ix,iy,0,0) + 2.0e0_R8*vbdia(ix,iy,0,0)
           ub2dia(ix,iy,0,1) = 
     &       vbecrb(ix,iy,0,1) + 2.0e0_R8*vbdia(ix,iy,0,1)
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
cwdk compuation based on dpo, dOnedBsq, ...
cwdk as long as y-faces are flux-surface aligned, no actual change in computation
             vbecrb(ix,iy,1,1) = radeldr*fac_ExB(ix,iy)/                 !srv 17.05.99 25.08.99 07.06.00 05.02.05 13.05.16 {
     1       (hx(ix,iy)+hx(bottomix(ix,iy),bottomiy(ix,iy)))* !signmf*   !srv
     2       (bb(ix,iy,2)/bb(ix,iy,3)**2                                 !srv 05.02.05
     2          +bb(bottomix(ix,iy),bottomiy(ix,iy),2)/
     2           bb(bottomix(ix,iy),bottomiy(ix,iy),3)**2)*              !srv 05.02.05
     2        dpo(ix,iy,1,0)                                             !srv 13.05.16 }
!           vbecrb(ix,iy,1) = radeldr*fac_ExB(ix,iy)* ! signmf*           !srv 11.01.13 07.06.00 21.04.00 25.08.99 17.05.99
!     1       (hz(ix,iy)*bb(ix,iy,2)                                      !lk( 12.12.07
!     2          +hz(bottomix(ix,iy),bottomiy(ix,iy))*
!     3          bb(bottomix(ix,iy),bottomiy(ix,iy),2))/2.0e0_R8*              
!     4       (vol(bottomix(ix,iy),bottomiy(ix,iy))/bb(ix,iy,3)**2      
!     5          +vol(ix,iy)/bb(bottomix(ix,iy),bottomiy(ix,iy),3)**2)/    
!     6          (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))*
!     7        (poc(rightix(ix,iy),rightiy(ix,iy))-poc(ix,iy))/           !lk 12.12.07)
!     8          (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))*
!     9          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
           if(ion_vlct_restrict.eq.2) then                               !iyv 06.02.14 {
             if( abs(vbecrb(ix,iy,1,1)) .gt. 2.0_R8*cssb(ix,iy,1) ) then
               ycount = ycount + 1
               if (vlct_diagno.ge.2) then
                 write(6,'(a,a,i3,a,i4,a,i4,a)')
     .            'Radial ExB velocity restriction applied ',
     .            'to species ',isb,' in cell (',ix,',',iy,'):'
                 write(6,'(a,1p,1e12.5,a,1e12.5)')
     .            'Mach 2: ',2.0_R8*cssb(ix,iy,1),
     .            ' Vy_ExB:',vbecrb(ix,iy,1,1)
               end if
               vbecrb(ix,iy,1,1)=
     .          sign(2.0_R8*cssb(ix,iy,1),vbecrb(ix,iy,1,1))
             endif
           endif                                                         !iyv 06.02.14  }
           vbdia(ix,iy,1,1) =  -raddidr*facdrift(ix,iy)* ! signmf*       !srv 11.01.13 07.06.00 17.05.99 25.04.00
     &           (ti(bottomix(ix,iy),bottomiy(ix,iy))                    !srv 05.03.99
     &           *bb(bottomix(ix,iy),bottomiy(ix,iy),2)                  !srv
     &           /(qe*zb(ix,iy)*hx(bottomix(ix,iy),bottomiy(ix,iy)))     !srv
     &           +ti(ix,iy)*bb(ix,iy,2)/(qe*zb(ix,iy)*hx(ix,iy)))*0.5_R8 !srv
     &        *dOnedBsq(ix,iy,1,0)                                       !srv 05.03.99
           wbdia(ix,iy,1,1) = facdrift(ix,iy)/                           !srv 13.05.16 {
     1      (hx(ix,iy)+hx(bottomix(ix,iy),bottomiy(ix,iy)))*
     2      (bb(ix,iy,2)/bb(ix,iy,3)**2
     2       +bb(bottomix(ix,iy),bottomiy(ix,iy),2)/
     2       bb(bottomix(ix,iy),bottomiy(ix,iy),3)**2)*
     2      dnbti(ix,iy,1,0)/
     &     (nbzbc(rightix(ix,iy),rightiy(ix,iy))+nbzbc(ix,iy))*2.0_R8/qe !srv 13.05.16 }

           ubdia(ix,iy,1,0) = vbecrb(ix,iy,1,0) + vbdia(ix,iy,1,0)
           ubdia(ix,iy,1,1) = vbecrb(ix,iy,1,1) + vbdia(ix,iy,1,1)       !srv 17.05.99
           ub2dia(ix,iy,1,0) = 
     &       vbecrb(ix,iy,1,0) + 2.0e0_R8*vbdia(ix,iy,1,0)
           ub2dia(ix,iy,1,1) = 
     &       vbecrb(ix,iy,1,1) + 2.0e0_R8*vbdia(ix,iy,1,1)               !srv 01.10.99
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
       endif                                                             !srv 21.06.08 {

*   ..compute P.Sch x-flux                                               !srv 13.10.06 {
       if(PSch_style.eq.0) then                                          !srv 12.11.08
        nbti=nb*ti
        call diff(nx, ny, hx, hy, vol, qs, nbti, dnbti)
        do iy = 0, ny-1
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. 
     &     bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 ) then        !srv 26.05.08
           fnbPSch(ix,iy,0,0:1)=0.0e0_R8
           conb(ix,iy,0,0:1,3) = 0.0e0_R8                                !srv 26.05.08
           conb(ix,iy,0,0:1,4) = 0.0e0_R8
          else
           ixl = leftix(ix,iy)
           iyl = leftiy(ix,iy)
           fnbPSch(ix,iy,0,0) = 0.25e0_R8*PSch*( 
     &      cddb(topix(ix,iy),topiy(ix,iy),1,1)*
     &       dnbti(topix(ix,iy),topiy(ix,iy),1,1)+
     &      cddb(topix(ixl,iyl),topiy(ixl,iyl),1,1)*
     &       dnbti(topix(ixl,iyl),topiy(ixl,iyl),1,1)+
     &      cddb(ixl,iyl,1,1)*
     &       dnbti(ixl,iyl,1,1)+
     &      cddb(ix,iy,1,1)*
     &       dnbti(ix,iy,1,1)
     &                                       )
cwdk computation of radial PSch flux on x-faces
cwdk somewhat messy interpolation of geometric factors, for exact
cwdk cancelation of some terms in poloidal and radial flux components
           fnbPSch(ix,iy,0,1) = 0.25e0_R8*PSch*( 
     &      cddb(topix(ix,iy),topiy(ix,iy),1,1)*
     &       dnbti(topix(ix,iy),topiy(ix,iy),1,0)*
     &       qs(topix(ix,iy),topiy(ix,iy),1)*
     &      (hy(ix,iy) + hy(topix(ix,iy),topiy(ix,iy)))/
     &      (hx(ix,iy) + hx(topix(ix,iy),topiy(ix,iy)))+
     &      cddb(topix(ixl,iyl),topiy(ixl,iyl),1,1)*
     &       dnbti(topix(ixl,iyl),topiy(ixl,iyl),1,0)*
     &       qs(topix(ixl,iyl),topiy(ixl,iyl),1)*
     &      (hy(ixl,iyl) + hy(topix(ixl,iyl),topiy(ixl,iyl)))/
     &      (hx(ixl,iyl) + hx(topix(ixl,iyl),topiy(ixl,iyl)))+
     &      cddb(ixl,iyl,1,1)*
     &       dnbti(ixl,iyl,1,0)*qs(ixl,iyl,1)*
     &      (hy(ixl,iyl) + hy(bottomix(ixl,iyl),bottomiy(ixl,iyl)))/
     &      (hx(ixl,iyl) + hx(bottomix(ixl,iyl),bottomiy(ixl,iyl)))+
     &      cddb(ix,iy,1,1)*
     &       dnbti(ix,iy,1,0)*qs(ix,iy,1)*
     &      (hy(ix,iy) + hy(bottomix(ix,iy),bottomiy(ix,iy)))/
     &      (hx(ix,iy) + hx(bottomix(ix,iy),bottomiy(ix,iy)))
     &                                       )
!           conb(ix,iy,0,0,3) = 0.25e0_R8*PSch*cddb(ix,iy,1,1)*ti(ix,iy)  !srv 21.10.10 {
!           conb(ix,iy,0,0,4) = 0.25e0_R8*PSch*cddb(ix,iy,1,1)*
!     &                      ti(bottomix(ix,iy),bottomiy(ix,iy))
           conb(ix,iy,0,0,3) = 0.25_R8*PSch*cddb(ix,iy,1,1)*
     &                         kmprti(ix,iy,1)
           conb(ix,iy,0,1,3) = -conb(ix,iy,0,0,3)
           conb(ix,iy,0,0,4) = conb(ix,iy,0,0,3)                          !srv 21.10.10 }
           conb(ix,iy,0,1,4) = conb(ix,iy,0,1,3)
          endif
         enddo
        enddo
       else
        weight = 1.0e0_R8                                               !lk
        nbti=nb*ti
        call diff(nx, ny, hx, hy, weight, qs, nbti, dnbti)
        call intcor(nx,ny,weight,nb,nbcor)                              !lk 09.06.11
        do iy = 0, ny-1
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. 
     &     bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 ) then     !srv 26.05.08
           fnbPSch(ix,iy,0,0:1)=0.0e0_R8
           conb(ix,iy,0,0:1,3) = 0.0e0_R8                                  !srv 26.05.08
           conb(ix,iy,0,0:1,4) = 0.0e0_R8
          else
           ixl = leftix(ix,iy)
           iyl = leftiy(ix,iy)
!           fnbPSch(ix,iy,0,0) = PSch*                                           !lk {
!     &      (nbtic(topix(ix,iy),topiy(ix,iy))-nbtic(ix,iy))*
!     &      (bzb(ix,iy)*hz(ix,iy)*vol(ixl,iyl)
!     &      +bzb(ixl,iyl)*hz(ixl,iyl)*vol(ix,iy))
!     &      /(vol(ixl,iyl)+vol(ix,iy))/(qe*zb(ix,iy))                          !lk }
           fnbPSch(ix,iy,0,0) = PSch*dnbti(ix,iy,0,1)*                          !lk {
     &      (bzb(ix,iy)*(vol(ix,iy)/gs(ix,iy,2))*vol(ixl,iyl)
     &      +bzb(ixl,iyl)*(vol(ixl,iyl)/gs(ixl,iyl,2))*vol(ix,iy))
     &      /(vol(ixl,iyl)+vol(ix,iy))/(qe*zb(ix,iy))                          !lk }
     &      *(nb(ix,iy)+nb(ixl,iyl))/(nbcor(ix,iy)+                       !lk 09.06.11
     &      nbcor(topix(ix,iy),topiy(ix,iy)))                             !lk 09.06.11
cwdk computation of radial PSch flux on x-faces
cwdk somewhat messy interpolation of geometric factors, for exact
cwdk cancelation of some terms in poloidal and radial flux components
           fnbPSch(ix,iy,0,1) = PSch*dnbti(ix,iy,0,0)*
     &      (hy(ix,iy) + hy(leftix(ix,iy),leftiy(ix,iy)))/
     &      (hx(ix,iy) + hx(leftix(ix,iy),leftiy(ix,iy)))*
     &       qs(ix,iy,0)*
     &      (bzb(ix,iy)*(vol(ix,iy)/gs(ix,iy,2))*vol(ixl,iyl)
     &      +bzb(ixl,iyl)*(vol(ixl,iyl)/gs(ixl,iyl,2))*vol(ix,iy))
     &      /(vol(ixl,iyl)+vol(ix,iy))/(qe*zb(ix,iy))
     &      *(nb(ix,iy)+nb(ixl,iyl))/(nbcor(ix,iy)+
     &      nbcor(topix(ix,iy),topiy(ix,iy)))

!            conb(ix,iy,0,0,3) = 0.25e0_R8*PSch*cddb(ix,iy,1,1)*ti(ix,iy)   !srv 21.10.10 {
!            conb(ix,iy,0,0,4) = 0.25e0_R8*PSch*cddb(ix,iy,1,1)*
!     &                      ti(bottomix(ix,iy),bottomiy(ix,iy))
           conb(ix,iy,0,0,3) = 0.25_R8*PSch*cddb(ix,iy,1,1)*
     &                         kmprti(ix,iy,1)
           conb(ix,iy,0,1,3) = -conb(ix,iy,0,0,3)
           conb(ix,iy,0,0,4) = conb(ix,iy,0,0,3)                            !srv 21.10.10 }
           conb(ix,iy,0,1,4) = conb(ix,iy,0,1,3)
          endif
         enddo
        enddo
       endif

*   ..compute P.Sch y-flux
       iy = -1
       fnbPSch(:,iy,1,1)=0.0e0_R8
       if(PSch_style.eq.0) then                                           !srv 12.11.08
        do iy = 0, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
           fnbPSch(ix,iy,1,0:1)= 0.0e0_R8
           conb(ix,iy,1,0:1,3) = 0.0e0_R8                         !srv 26.05.08
           conb(ix,iy,1,0:1,4) = 0.0e0_R8
          else
           ixb = bottomix(ix,iy)
           iyb = bottomiy(ix,iy)
           fnbPSch(ix,iy,1,1) = -0.25e0_R8*PSch*(
     &    cddb(rightix(ix,iy),rightiy(ix,iy),0,0)*
     &     dnbti(rightix(ix,iy),rightiy(ix,iy),0,0) +
     &    cddb(ix,iy,0,0)*
     &     dnbti(ix,iy,0,0)+
     &    cddb(bottomix(ix,iy),bottomiy(ix,iy),0,0)*
     &     dnbti(bottomix(ix,iy),bottomiy(ix,iy),0,0)+
     &    cddb(rightix(ixb,iyb),rightiy(ixb,iyb),0,0)*
     &     dnbti(rightix(ixb,iyb),rightiy(ixb,iyb),0,0)
     &                                         )

!           conb(ix,iy,1,1,3) = 0.25e0_R8*PSch*cddb(ix,iy,0,0)*ti(ix,iy)  !srv 21.10.10 {
!           conb(ix,iy,1,1,4) = 0.25e0_R8*PSch*cddb(ix,iy,0,0)*
!     &                        ti(leftix(ix,iy),leftiy(ix,iy))
           conb(ix,iy,1,1,3) = 0.25_R8*PSch*cddb(ix,iy,0,0)*
     &                         kmprti(ix,iy,0)
           conb(ix,iy,1,1,4) = conb(ix,iy,1,1,3)                          !srv 21.10.10 }
          endif
         enddo
        enddo
       else
        do iy = 0, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
           fnbPSch(ix,iy,1,0:1)= 0.0e0_R8
           conb(ix,iy,1,0:1,3) = 0.0e0_R8                           !srv 26.05.08
           conb(ix,iy,1,0:1,4) = 0.0e0_R8
          else
           ixb = bottomix(ix,iy)
           iyb = bottomiy(ix,iy)
!           fnbPSch(ix,iy,1,1) = -PSch*                                          !lk {
!     &      (nbtic(rightix(ix,iy),rightiy(ix,iy))-nbtic(ix,iy))*
!     &      (bzb(ix,iy)*hz(ix,iy)*vol(ixb,iyb)
!     &      +bzb(ixb,iyb)*hz(ixb,iyb)*vol(ix,iy))
!     &      /(vol(ixb,iyb)+vol(ix,iy))/(qe*zb(ix,iy))                          !lk }
           fnbPSch(ix,iy,1,1) = -PSch*                                          !lk {
     &      dnbti(ix,iy,1,0)*
     &      (bzb(ix,iy)*(vol(ix,iy)/gs(ix,iy,2))*vol(ixb,iyb)
     &      +bzb(ixb,iyb)*(vol(ixb,iyb)/gs(ixb,iyb,2))*vol(ix,iy))
     &      /(vol(ixb,iyb)+vol(ix,iy))/(qe*zb(ix,iy))                   !lk }
     &      *(nb(ix,iy)+nb(ixb,iyb))/(nbcor(ix,iy)+                     !lk 09.06.11
     &      nbcor(rightix(ix,iy),rightiy(ix,iy)))                       !lk 09.06.11

!           conb(ix,iy,1,1,3) = 0.25e0_R8*PSch*cddb(ix,iy,0,0)*ti(ix,iy)    !srv 21.10.10 {
!           conb(ix,iy,1,1,4) = 0.25e0_R8*PSch*cddb(ix,iy,0,0)*
!     &                        ti(leftix(ix,iy),leftiy(ix,iy))
           conb(ix,iy,1,1,3) = 0.25_R8*PSch*cddb(ix,iy,0,0)*
     &                         kmprti(ix,iy,0)
           conb(ix,iy,1,1,4) = conb(ix,iy,1,1,3)                          !srv 21.10.10 }
          endif
         enddo
        enddo
       endif

        do iy = -1, ny                                                  !srv 02.01.07 {
          do ix = -1, nx
            fnbPSch(ix,iy,0:1,0:1) = 
     .        fnbPSch(ix,iy,0:1,0:1)*facdrift(ix,iy)
            conb(ix,iy,:,:,3) = conb(ix,iy,:,:,3)*facdrift(ix,iy)
            conb(ix,iy,:,:,4) = conb(ix,iy,:,:,4)*facdrift(ix,iy)
          enddo
        enddo                                                           !srv 02.01.07 }

      endif

      if (isb.ne.ismain .and. .not.is_neutral(isb)) then                !lk 02.06.11 {
        do iy = 0, ny
          do ix = -1, nx
            if (leftix(ix,iy).ne.-2) then
              t0 = max(ti(ix,iy), ti(leftix(ix,iy),leftiy(ix,iy)))
              t1 = (nb(ix,iy) + nb(leftix(ix,iy),leftiy(ix,iy)))
     &          /2.0_R8
              vbar =
     .         abs(bb(ix,iy,0) + bb(leftix(ix,iy),leftiy(ix,iy),0))     !srv 11.01.13
     &           /(bb(ix,iy,3) + bb(leftix(ix,iy),leftiy(ix,iy),3))
     &          *sqrt(8.0_R8 * t0 / (pi * am(isb)*mp))
              t2 = 1.0_R8 + abs(fnbPSch(ix,iy,0,0))/(
     &         vbar*t1*(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     &         (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))
     &          )
              fnbPSch(ix,iy,0,0)= fnbPSch(ix,iy,0,0)/t2
              fnbPSch(ix,iy,0,1)= fnbPSch(ix,iy,0,1)/t2
              vbdia(ix,iy,0,0)= vbdia(ix,iy,0,0)/
     &                            (1.0_R8 + abs(vbdia(ix,iy,0,0))/vbar) 
              wbdia(ix,iy,0,0)= wbdia(ix,iy,0,0)/
     &                            (1.0_R8 + abs(wbdia(ix,iy,0,0))/vbar) 
              ubdia(ix,iy,0,0)= ubdia(ix,iy,0,0)/
     &                            (1.0_R8 + abs(ubdia(ix,iy,0,0))/vbar) 
              ub2dia(ix,iy,0,0)= ub2dia(ix,iy,0,0)/
     &                            (1.0_R8 + abs(ub2dia(ix,iy,0,0))/vbar) 
              vbdia(ix,iy,0,1)= vbdia(ix,iy,0,1)/
     &                            (1.0_R8 + abs(vbdia(ix,iy,0,1))/vbar) 
              wbdia(ix,iy,0,1)= wbdia(ix,iy,0,1)/
     &                            (1.0_R8 + abs(wbdia(ix,iy,0,1))/vbar) 
              ubdia(ix,iy,0,1)= ubdia(ix,iy,0,1)/
     &                            (1.0_R8 + abs(ubdia(ix,iy,0,1))/vbar) 
              ub2dia(ix,iy,0,1)= ub2dia(ix,iy,0,1)/
     &                            (1.0_R8 + abs(ub2dia(ix,iy,0,1))/vbar) 
            endif
	    if (bottomiy(ix,iy).ne.-2) then
	      t0 = max(ti(ix,iy), ti(bottomix(ix,iy),bottomiy(ix,iy)))
	      t1 = (nb(ix,iy) + nb(bottomix(ix,iy),bottomiy(ix,iy)))
     &	       /2.0_R8
	      vbar=
     .	       abs(bb(ix,iy,0) + bb(bottomix(ix,iy),bottomiy(ix,iy),0))  !srv 11.01.13
     &	       /(bb(ix,iy,3) + bb(bottomix(ix,iy),bottomiy(ix,iy),3))
     &	       *sqrt(8.0_R8 * t0 / (pi * am(isb)*mp))
	      t2 = 1.0_R8 + abs(fnbPSch(ix,iy,1,1))/(
     &	       vbar*t1*(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))
     &	       /(hy(ix,iy)+hy(bottomix(ix,iy),bottomiy(ix,iy)))
     &	       )
	      fnbPSch(ix,iy,1,1)= fnbPSch(ix,iy,1,1)/t2
	      vbdia(ix,iy,1,1)= vbdia(ix,iy,1,1)/
     &	       (1.0_R8 + abs(vbdia(ix,iy,1,1))/vbar) 
	      wbdia(ix,iy,1,1)= wbdia(ix,iy,1,1)/
     &	       (1.0_R8 + abs(wbdia(ix,iy,1,1))/vbar) 
	      ubdia(ix,iy,1,1)= ubdia(ix,iy,1,1)/
     &	       (1.0_R8 + abs(ubdia(ix,iy,1,1))/vbar) 
	      ub2dia(ix,iy,1,1)= ub2dia(ix,iy,1,1)/
     &	       (1.0_R8 + abs(ub2dia(ix,iy,1,1))/vbar) 
	    endif
	  enddo
        enddo
      endif                                                             !lk 02.06.11 }
#endif

*   ..interpolate na and pa to cell faces
      weight = 1.0_R8
      call intface(nCv,nFc,mpg%fcCv,weight,pl%na(:,isb),nbf)
      call intface(nCv,nFc,mpg%fcCv,weight,dv%pa(:,isb),pbf)

*   ..compute differences of na and pa at cell faces
!      call diff(nCv, nFc, nVx, 0, geo, mpg, pl%na(:,isb), wrkVx, dnb)
      call diff(nCv, nFc, nVx, 0, geo, mpg, dv%pa(:,isb), wrkVx, dpb)

*   ..interpolate cvH ('hz') to cell faces
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,geo%cvH,wrk3)

* ..poloidal fluxes

*   ..interpolate ua to cell faces
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,pl%ua(:,isb),wrk0)

*   ..compute wrk1
      do iFc = 1, nFc
        wrk1(iFc) = 0.5_R8*
     &    (dv%pa(mpg%fcCv(iFc,2),isb)-dv%pa(mpg%fcCv(iFc,1),isb))/
     &    (dv%pa(mpg%fcCv(iFc,2),isb)+dv%pa(mpg%fcCv(iFc,1),isb))
      enddo
!      call b2xxgs (nx, ny, 0.0_R8, wrk1, 1)
*   ..compute wrk2
!      do iy = -1, ny
!       do ix = -1, nx
!        if(leftix(ix,iy).eq.-2) then
!         wrk2(ix,iy) = 0.0_R8
!        elseif(leftix(leftix(ix,iy),leftiy(ix,iy)).eq.-2) then
!         wrk2(ix,iy) = 0.0_R8
!        elseif(rightix(ix,iy).eq.nx+1) then
!         wrk2(ix,iy) = 0.0_R8
!        else
!         wrk2(ix,iy) = wrk1(leftix(ix,iy),leftiy(ix,iy))-
!     &    2.0_R8*wrk1(ix,iy)+wrk1(rightix(ix,iy),rightiy(ix,iy))
!        endif
!       enddo
!      enddo
!      call b2xxgs (nx, ny, 0.0_R8, wrk2, 1)
!cwdk detailed computation of V_corr_dpc to be checked... tmp: 0
!cwdk still need a volume ratio at the end of term and termhz, to make
!cwdk it zero on faces with guard cells.
      wrk2 = 0.0_R8
      term   = 
     &  xfrhie*wrk2*0.5_R8*abs(geo%fcPbs)*sqrt(pbf/(am(isb)*mp*nbf))
      termhz = 
     &  xfrhie*wrk2*0.5_R8*abs(geo%fcPbshz)*sqrt(pbf/(am(isb)*mp*nbf))

*   ..compute convective coefficients (poloidal and radial)
      flo(:,0) = geo%fcPbs*wrk0 + term +
     &  co%cvla(:,0,isb) + dv%uadia(:,0,isb)*geo%fcS*geo%fcQalf(:,0)
      flo(:,1) = 
     &  co%cvla(:,1,isb) + dv%uadia(:,1,isb)*geo%fcS*geo%fcQalf(:,1)
      flo_mdf(:,0) = geo%fcPbs*wrk0 + term +
     &  co%cvla(:,0,isb) + dv%vaecrb(:,0,isb)*geo%fcS*geo%fcQalf(:,0)
      flo_mdf(:,1) = 
     &  co%cvla(:,1,isb) + dv%vaecrb(:,1,isb)*geo%fcS*geo%fcQalf(:,1)
      flo_he(:,0) = geo%fcPbs*wrk0 + term +
     &  5.0_R8/3.0_R8*co%cvla(:,0,isb) + 
     &  dv%vaecrb(:,0,isb)*geo%fcS*geo%fcQalf(:,0)
      flo_he(:,1)  = 5.0_R8/3.0_R8*co%cvla(:,1,isb) +
     &  dv%vaecrb(:,1,isb)*geo%fcS*geo%fcQalf(:,1)
      flo2dia(:,0) = geo%fcPbshz*wrk0 + termhz +
     &  co%cvlahz(:,0,isb) + ub2dia(:,0)*geo%fcS*geo%fcQalf(:,0)*wrk3
      flo2dia(:,1) = 
     &  co%cvlahz(:,1,isb) + ub2dia(:,1)*geo%fcS*geo%fcQalf(:,1)*wrk3
      flo_nodrift(:,0) = co%cvla(:,0,isb) + geo%fcPbs*wrk0 + term
      flo_nodrift(:,1) = co%cvla(:,1,isb)

*   ..add contribution from currents
      if (nint(zn(isb)).eq.1.and.nint(zamax(isb)).eq.1) then
        scur(:,0) = 
     &     (dv%fchinert(:,0) + dv%fchvispar(:,0) + dv%fchvisper(:,0) +
     &      dv%fchanml(:,0)  + dv%fchvisq(:,0))*xcur/(qe*nbf)
        scur(:,1) = 
     &     (dv%fchinert(:,1) + dv%fchvispar(:,1) + 
     &      vis_per*dv%fchvisper(:,0) + dv%fchanml(:,1) + 
     &      vis_q*dv%fchvisq(:,1))*ycur/(qe*nbf)
        flo     = flo     + scur
        flo_mdf = flo_mdf + scur
        flo2dia(:,0) = flo2dia(:,0) + scur(:,0)*wrk3
        flo2dia(:,1) = flo2dia(:,1) + scur(:,1)*wrk3
      end if

*   ..apply discretization scheme for density-part of fna_mdf
      call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%na(:,isb),flo_mdf,
     &  co%cdna(:,:,isb),fna_mdf0(:,:),dv%fna_32(:,:,isb),
     &  dv%fna_52(:,:,isb))

*   ..compute fluxes
      dv%fna(:,0,isb)      = fna_mdf0(:,0) +  
     &      (flo(:,0) - flo_mdf(:,0))*nbf - co%cdpa(:,0,isb)*dpb(:,0)
      dv%fna(:,1,isb)      = fna_mdf0(:,1) +  
     &      (flo(:,1) - flo_mdf(:,1))*nbf - co%cdpa(:,1,isb)*dpb(:,1)

      dv%fna_mdf(:,0,isb)  = fna_mdf0(:,0) - 
     &      co%cdpa(:,0,isb)*dpb(:,0) - dv%fnaPSch(:,0,isb)
      dv%fna_mdf(:,1,isb)  = fna_mdf0(:,1) - 
     &      co%cdpa(:,1,isb)*dpb(:,1) - dv%fnaPSch(:,1,isb)

      dv%fna_fcor(:,0,isb) = fna_mdf0(:,0)*wrk3 + 
     &      (flo2dia(:,0) - flo_mdf(:,0)*wrk3)*nbf - 
     &       co%cdpahz(:,0,isb)*dpb(:,0)
      dv%fna_fcor(:,1,isb) = fna_mdf0(:,1)*wrk3 + 
     &      (flo2dia(:,1) - flo_mdf(:,1)*wrk3)*nbf - 
     &       co%cdpahz(:,1,isb)*dpb(:,1)

      dv%fna_he(:,0,isb)   = fna_mdf0(:,0)*(5.0_R8/3.0_R8) + 
     &      (flo_he(:,0) - (5.0_R8/3.0_R8)*flo_mdf(:,0))*nbf - 
     &      (5.0_R8/3.0_R8)*co%cdpa(:,0,isb)*dpb(:,0)
      dv%fna_he(:,1,isb)   = fna_mdf0(:,1)*(5.0_R8/3.0_R8) + 
     &      (flo_he(:,1) - (5.0_R8/3.0_R8)*flo_mdf(:,1))*nbf - 
     &      (5.0_R8/3.0_R8)*co%cdpa(:,1,isb)*dpb(:,1)

      dv%fna_32(:,0,isb)   = dv%fna_32(:,0,isb) + 
     &     (flo(:,0) - flo_mdf(:,0))*nbf
      dv%fna_32(:,1,isb)   = dv%fna_32(:,1,isb) + 
     &     (flo(:,1) - flo_mdf(:,1))*nbf

      dv%fna_52(:,0,isb)   = dv%fna_52(:,0,isb) - 
     &     co%cdpa(:,0,isb)*dpb(:,0)
      dv%fna_52(:,1,isb)   = dv%fna_52(:,1,isb) - 
     &     co%cdpa(:,1,isb)*dpb(:,1)

*   ..apply discretization scheme for density-part of fna_nodrift
      call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%na(:,isb),flo_nodrift,
     &  co%cdna(:,:,isb),dv%fna_nodrift(:,:,isb),dv%fna_32nd(:,0,isb),
     &  dv%fna_52nd(:,0,isb))

*   ..compute _nodrift fluxes
      dv%fna_nodrift(:,0,isb) = dv%fna_nodrift(:,0,isb) - 
     &      co%cdpa(:,0,isb)*dpb(:,0)
      dv%fna_nodrift(:,1,isb) = dv%fna_nodrift(:,1,isb) - 
     &      co%cdpa(:,1,isb)*dpb(:,1)

      dv%fna_52nd(:,0,isb) = dv%fna_52nd(:,0,isb) - 
     &      co%cdpa(:,0,isb)*dpb(:,0)
      dv%fna_52nd(:,1,isb) = dv%fna_52nd(:,1,isb) - 
     &      co%cdpa(:,1,isb)*dpb(:,1)

*   ..compute flob, conb
      call calccoef(nCv,nFc,nVx,meth,geo,mpg,
     &    flo_mdf,co%cdna(:,:,isb),flob0,conb0)
      if (mdf_fnb .eq. 0) then
        flob0 = flob0 + flo - flo_mdf
      endif
      call intface(nCv,nFc,mpg%fcCv,weight,kmpr,kmprf)
      dv%flob(:,0)   = flob0(:,0)*kmprf 
      dv%flob(:,1)   = flob0(:,1)*kmprf 
      dv%conb(:,0,0) = conb0(:,0)*kmprf 
      dv%conb(:,1,0) = conb0(:,1)*kmprf 

*   ..contribution from pressure diffusion
      dumm0 = 0.0_R8
      call calccoef(nCv,nFc,nVx,2,geo,mpg,
     &  dumm0,co%cdpa(:,:,isb),flob0,conb0)
      dv%flob = dv%flob + flob0                                     !wdk flob0 should be zero here
      dv%conb(:,:,0) = dv%conb(:,:,0) + conb0
      dv%conb(:,0,1) = conb0(:,0)*kmprti
      dv%conb(:,1,1) = conb0(:,1)*kmprti
      dv%conb(:,:,2) = dv%conb(:,:,1)

*   ..flux limit
      if ((alpha.gt.0.0_R8 .and. maxval(rt%rza(:,isb)).eq.0.0_R8) .and.          !srv 20.11.08
     *     alpha_.eq.0.0_R8 ) then                                   ! neutral flux limit
cwdk try to introduce a common flux limit here?
       do iFc = 1,nFc
        t2 = max(flux_limit_min_ti,
     *        pl%ti(mpg%fcCv(iFc,1)), pl%ti(mpg%fcCv(iFc,2)))
        vbar = sqrt(8.0_R8 * t2 / (pi * am(isb)*mp))
        t0 = (1.0_R8 + (abs(dv%fna(iFc,0,isb) / (0.5_R8*
     &   (pl%na(mpg%fcCv(iFc,1),isb) + pl%na(mpg%fcCv(iFc,2),isb))* 
     &    alpha*vbar/4.0_R8*geo%fcS(iFc)*geo%fcQalf(iFc,0))))**gamma)**
     &   (1.0_R8/gamma)
        t1 = (1.0_R8 + (abs(dv%fna(iFc,1,isb) / (0.5_R8*
     &   (pl%na(mpg%fcCv(iFc,1),isb) + pl%na(mpg%fcCv(iFc,2),isb))* 
     &    alpha*vbar/4.0_R8*geo%fcS(iFc)*geo%fcQalf(iFc,1))))**gamma)**
     &   (1.0_R8/gamma)
*     ..poloidal direction
        maxfluxlimit(0)=max(maxfluxlimit(0),t0)
        dv%fna(iFc,0,isb) = dv%fna(iFc,0,isb) / t0
        dv%fna_32(iFc,0,isb) = dv%fna_32(iFc,0,isb) / t0
        dv%fna_52(iFc,0,isb) = dv%fna_52(iFc,0,isb) / t0
        dv%fna_nodrift(iFc,0,isb) = dv%fna_nodrift(iFc,0,isb) / t0
        dv%fna_32nd(iFc,0,isb) = dv%fna_32nd(iFc,0,isb) / t0
        dv%fna_52nd(iFc,0,isb) = dv%fna_52nd(iFc,0,isb) / t0
        dv%flob(iFc,0) = dv%flob(iFc,0) / t0                         ! to check: flob(:,0), conb(:,0,:) refers to connector, not poloidal... do flux limit before computation of flob/conb?
        dv%conb(iFc,0,0:4) = dv%conb(iFc,0,0:4) / t0                 !srv 13.10.06 {
        dv%fna_he(iFc,0,isb)   = dv%fna_he(iFc,0,isb)/t0
        dv%fna_fcor(iFc,0,isb) = dv%fna_fcor(iFc,0,isb)/t0
        dv%fna_mdf(iFc,0,isb) = dv%fna_mdf(iFc,0,isb)/t0
        dv%fnaPSch(iFc,0,isb)  = dv%fnaPSch(iFc,0,isb)/t0            !srv 13.10.06 }
        co%fllim0fna(iFc,0,isb) = t0
*     ..radial direction
        maxfluxlimit(1)=max(maxfluxlimit(1),t1)
        dv%fna(iFc,1,isb) = dv%fna(iFc,1,isb) / t1
        dv%fna_32(iFc,1,isb) = dv%fna_32(iFc,1,isb) / t1
        dv%fna_52(iFc,1,isb) = dv%fna_52(iFc,1,isb) / t1
        dv%fna_nodrift(iFc,1,isb) = dv%fna_nodrift(iFc,1,isb) / t1
        dv%fna_32nd(iFc,1,isb) = dv%fna_32nd(iFc,1,isb) / t1
        dv%fna_52nd(iFc,1,isb) = dv%fna_52nd(iFc,1,isb) / t1
        dv%flob(iFc,1) = dv%flob(iFc,1) / t1                         ! to check: flob(:,1), conb(:,1,:) refers to tangential, not radial...
        dv%conb(iFc,1,0:4) = dv%conb(iFc,1,0:4) / t1                 !srv 13.10.06 {
        dv%fna_he(iFc,1,isb)   = dv%fna_he(iFc,1,isb)/t1
        dv%fna_fcor(iFc,1,isb) = dv%fna_fcor(iFc,1,isb)/t1
        dv%fna_mdf(iFc,1,isb) = dv%fna_mdf(iFc,1,isb)/t1
        dv%fnaPSch(iFc,1,isb)  = dv%fnaPSch(iFc,1,isb)/t1            !srv 13.10.06 }
        co%fllim0fna(iFc,1,isb) = t1
       enddo
      endif



      if(iout_b2npc9.ne.0 .or. iout_b2npco.ne.0
     &   .or. iout .eq. 1 .or. iout_b2wdat.eq.4) then                   !srv 11.09.09 {
cwdk to be considered: con_mdf, con_2dia not easily recoverable? 
cwdk should this be conb in general? (but without the rescaling for pressure)
       call diff(nCv, nFc, nVx, 0, geo, mpg, pl%na(:,isb), wrkVx, dnb) 
       wrk1 = geo%fcPbs*wrk0*nbf
       wrk2 = term*nbf
       wrk3 = dv%vaecrb(:,0,isb)*nbf*geo%fcS*geo%fcQalf(:,0)
!       wrk4 = -con_2dia*dnb(:,0)
       wrk5 = -co%cdpa(:,0,isb)*dpb(:,0)
       wrk6 = termhz*nbf
!       wrk7 = -con_mdf*dnb(:,0)
       wrk8 = co%cvla(:,0,isb)*nbf
       if(.not. is_neutral(isb)) then
        call my_out_us(70,nFc,1,wrk3,'b2tfnb_vaecrbnax'//chns)
        call my_out_us(70,nFc,1,dv%fnaPSch(1,0,isb),
     &   'b2tfnb_fnbPSchx'//chns)
        call my_out_us(70,nFc,1,wrk7,'b2tfnb_dPat_mdf_gradnax'//chns)
        call my_out_us(70,nFc,1,wrk4,'b2tfnb_dPat_2diagradnax'//chns)
       else
        call my_out_us(70,nFc,1,wrk5,'b2tfnb_dgradpbx'//chns)
       endif
       call my_out_us(70,nFc,1,wrk1,'b2tfnb_bxuanax'//chns)
       call my_out_us(70,nFc,1,wrk2,'b2tfnb_dpccornax'//chns)
       if (iout_b2wdat .ne. 4) 
     &  call my_out_us(70,nFc,1,wrk6,'b2tfnb_dpccorhznax'//chns)            !srv 05.11.09
       if (iout.eq.1) 
     &  call my_out_us(70,nFc,1,wrk8,'b2tfnb_cvlbnax'//chns)                !iys 17.04.2012
!       call my_out(70,nx,ny,con_2dia,'b2tfnb_dPat_2diax'//chns)
!       call my_out(70,nx,ny,con_mdf,'b2tfnb_dPat_mdfx'//chns)
      endif
      if(iout_b2npmo.ne.0 .or. iout_b2wdat.eq.4 .or. iout.eq.1) then
       write(chns,'(i3.3)') isb
       wrk6 = dv%vadia(:,0,isb)*nbf*geo%fcS*geo%fcQalf(:,0)
       if(.not. is_neutral(isb)) 
     &  call my_out_us(70,nFc,1,wrk6,'b2tfnb_vadianax'//chns)
       call intface(nCv,nFc,mpg%fcCv,geo%fcVol,geo%cvH,wrk3)
       wrk6 = dv%vadia(:,0,isb)*nbf*geo%fcS*geo%fcQalf(:,0)*wrk3
       if(.not. is_neutral(isb) .and. iout_b2wdat.ne.4)
     &   call my_out_us(70,nFc,1,wrk6,'b2tfnb_vadianahzx'//chns)
       call my_out_us(70,nFc,1,dv%fna_fcor(1,0,isb),
     &              'b2tfnb_fnb_fcorx'//chns)
      endif                                                     !srv 11.09.09 }

      if (alpha.gt.0.0_R8 .and. maxval(rt%rza(:,isb)).eq.0.0_R8 
     > .and. ncall.lt.100) then 
       write(*,'(a,1p,e12.5,a,e12.5,a,2e12.5)') 
     .   'b2tfnb: maxfluxlimit for zb = ',maxval(rt%rza(:,isb)),
     1   ' mb = ',am(isb),' are ',maxfluxlimit
      endif

      if(iout_b2npc9.ne.0 .or. iout_b2npco.ne.0
     &       .or. iout .eq. 1 .or. iout_b2wdat.eq.4) then               !srv 11.09.09 {
       wrk3 = dv%vaecrb(:,1,isb)*nbf*geo%fcS*geo%fcQalf(:,1)
!       wrk4 = -con_2dia*dnb(:,1)
       wrk5 = -co%cdpa(:,1,isb)*dpb(:,1)
!       wrk7 = -con_mdf*dnb(:,1)
       wrk8 = co%cvla(:,1,isb)*nbf
       write(chns,'(i3.3)') isb
       if(.not. is_neutral(isb)) then
        call my_out_us(70,nFc,1,dv%fnaPSch(:,1,isb),
     &   'b2tfnb_fnbPSchy'//chns)
        call my_out_us(70,nFc,1,wrk3,'b2tfnb_vaecrbnay'//chns)
        call my_out_us(70,nFc,1,wrk7,'b2tfnb_dPat_mdf_gradnay'//chns)
        call my_out_us(70,nFc,1,wrk4,'b2tfnb_dPat_2diagradnay'//chns)
       else
        call my_out_us(70,nFc,1,wrk5,'b2tfnb_dgradpby'//chns)
       endif
       call my_out_us(70,nFc,1,wrk8,'b2tfnb_cvlbnay'//chns)               !iys 17.04.2012
!       call my_out_us(70,nFc,1,con_2dia,'b2tfnb_dPat_2diay'//chns)
!       call my_out_us(70,nFc,1,con_mdf,'b2tfnb_dPat_mdfy'//chns)
      endif
      if(iout_b2npmo.ne.0 .or. iout_b2wdat.eq.4 .or. iout.eq.1) then
       write(chns,'(i3.3)') isb
       wrk6 = dv%vadia(:,1,isb)*nbf*geo%fcS*geo%fcQalf(:,1)
       if(.not. is_neutral(isb)) 
     &  call my_out_us(70,nFc,1,wrk6,'b2tfnb_vadianay'//chns)
       call intface(nCv,nFc,mpg%fcCv,geo%fcVol,geo%cvH,wrk3)
       wrk6 = dv%vadia(:,1,isb)*nbf*geo%fcS*geo%fcQalf(:,1)*wrk3 
       if(.not. is_neutral(isb) .and. iout_b2wdat.ne.4)
     &  call my_out_us(70,nFc,1,wrk6,'b2tfnb_vadianahzy'//chns)
      call my_out_us(70,nFc,1,dv%fna_fcor(1,1,isb),
     &  'b2tfnb_fnb_fcory'//chns)
      endif                                                              !srv 11.09.09 }


* ..compute kinetic energy of species
      vbecrb_ = 0.0_R8
      wbdia_  = 0.0_R8
#ifdef WG_DRIFTS
      if((maxval(dv%facdrift).ne.0.0_R8 .or. maxval(dv%fac_ExB).ne.0.0_R8)
     &    .and. minval(zb).gt.0.0_R8) then
        wrk0=0.0_R8  ! compute electric field on the cell faces
        weight = 1.0e0_R8
        call grad(nx, ny, hx, hy, weight, qz, qs, po, dpo)
        do iy=-1,ny  ! compute perpendicular projection V_ExB on the cell centers
         do ix=-1,nx
          if(bottomiy(ix,iy).eq.-2) then
           vbecrb_(ix,iy,0) = dpo(ix,topiy(ix,iy),1,1)
          elseif(topiy(ix,iy).eq.ny+1) then
           vbecrb_(ix,iy,0) = dpo(ix,iy,1,1)
          else
           vbecrb_(ix,iy,0) = 0.5_R8*
     &        (dpo(topix(ix,iy),topiy(ix,iy),1,1)+dpo(ix,iy,1,1))
          endif
cwdk need full diamagnetic component for energy
           vbecrb_(ix,iy,0)=-vbecrb_(ix,iy,0)/
     &       bb(ix,iy,3)
         enddo
        enddo
        call grad(nx, ny, hx, hy, weight, qz, qs, nbti, dnbti) ! compute .... on the cell faces
        do iy=-1,ny  ! compute perpendicular projection V_dia on the cell centers
         do ix=-1,nx
          if(bottomiy(ix,iy).eq.-2) then
           wbdia_(ix,iy,0)=dnbti(ix,topiy(ix,iy),1,1)
          elseif(topiy(ix,iy).eq.ny+1) then
           wbdia_(ix,iy,0)=dnbti(ix,iy,1,1)
          else
           wbdia_(ix,iy,0)=0.5_R8*
     &        (dnbti(topix(ix,iy),topiy(ix,iy),1,1)+dnbti(ix,iy,1,1))
          endif
cwdk need full diamagnetic component for energy
           wbdia_(ix,iy,0)=-wbdia_(ix,iy,0)/
     &         (zb(ix,iy)*qe*nb(ix,iy)*bb(ix,iy,3))
         enddo
        enddo

        ! compute y-component of electric field on the cell faces
        do iy=-1,ny  ! compute .... on the cell centers
         do ix=-1,nx
          if(leftix(ix,iy).eq.-2) then
           vbecrb_(ix,iy,1)=dpo(rightix(ix,iy),iy,0,0)
          elseif(rightix(ix,iy).eq.nx+1) then
           vbecrb_(ix,iy,1)=dpo(ix,iy,0,0)
          else
           vbecrb_(ix,iy,1)=0.5_R8*
     &        (dpo(rightix(ix,iy),rightiy(ix,iy),0,0)+dpo(ix,iy,0,0))
          endif
            vbecrb_(ix,iy,1)= vbecrb_(ix,iy,1)*
     &       bb(ix,iy,2)/bb(ix,iy,3)**2                                 !srv 15.02.12
         enddo
        enddo

        do iy=-1,ny  ! compute .... on the cell centers
         do ix=-1,nx
          if(leftix(ix,iy).eq.-2) then
           wbdia_(ix,iy,1)=dnbti(rightix(ix,iy),iy,0,0)
          elseif(rightix(ix,iy).eq.nx+1) then
           wbdia_(ix,iy,1)=dnbti(ix,iy,0,0)
          else
           wbdia_(ix,iy,1)=0.5_R8*
     &       (dnbti(rightix(ix,iy),rightiy(ix,iy),0,0)+dnbti(ix,iy,0,0))
          endif
          wbdia_(ix,iy,1)=wbdia_(ix,iy,1)/
     &       (zb(ix,iy)*qe*nb(ix,iy))*bb(ix,iy,2)/bb(ix,iy,3)**2        !srv 15.02.12
         enddo
        enddo
       else
        wbdia_ = 0.0_R8
        vbecrb_= 0.0_R8
       endif

      endif                                                             !srv 14.07.10 }
#endif
      dv%kinrgy(:,isb) = 0.5_R8*am(isb)*mp*pl%ua(:,isb)**2 + 
     .     (xwdia*wbdia_(:,0)+xvecrb*vbecrb_(:,0))**2+                  ! wbdia_ is a value in the cell center
     .     (xwdia*wbdia_(:,1)+xvecrb*vbecrb_(:,1))**2                   ! vbecrb_ is a value in the cell center
#ifdef WG_DRIFTS
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
            write(*,fmt=
     $"('isb = ',i4,' fnb_tot = ',es18.11,' fnb_mdf_tot = ',es18.11)")
     $       isb, div_fnb, div_fnb_mdf

            call b2tbfl(nx, ny,.false., fnb, tbfl_fnb, div_fnb)
            call b2tbfl(nx, ny,.false., fnb_mdf, tbfl_fnb_mdf,
     .                                                 div_fnb_mdf)
            write(*,*) 'After flux correction in b2tfnb, is = ', isb
            write(*,fmt=
     $"('isb = ',i4,' fnb_tot = ',es18.11,' fnb_mdf_tot = ',es18.11)")
     $       isb, div_fnb, div_fnb_mdf
            write(*,*)
!            write(*,"(' is =',i3,' rdiv =',es15.7,' div_fnb =',es15.7,
!     &       ' div_fnb_mdf =',es15.7)") isb, rdiv, div_fnb, div_fnb_mdf
          endif
        endif
      endif                                                             !srv 05.03.16 }
#endif
      if (iout_b2wdat .eq. 4) then
       call my_out_us(70,nFc,1,dv%fna(1,0,isb),'b2tfnb_fnbx'//chns)
       call my_out_us(70,nFc,1,dv%fna(1,1,isb),'b2tfnb_fnby'//chns)
       call my_out_us(70,nFc,1,dv%fna_fcor(1,0,isb),
     &            'b2tfnb_fnb_fcorx'//chns)
       call my_out_us(70,nFc,1,dv%fna_fcor(1,1,isb),
     &            'b2tfnb_fnb_fcory'//chns)
      call my_out_us(70,nFc,1,dv%fna_he(1,0,isb),'b2tfnb_fnb_hex'//chns)
      call my_out_us(70,nFc,1,dv%fna_he(1,1,isb),'b2tfnb_fnb_hey'//chns)
       call my_out_us(70,nFc,1,dv%fnaPSch(1,0,isb),
     &            'b2tfnb_fnbPSchx'//chns)
       call my_out_us(70,nFc,1,dv%fnaPSch(1,1,isb),
     &            'b2tfnb_fnbPSchy'//chns)
       if (isb .eq. ismain) then
          call my_out_us(70,nFc,1,dv%vaecrb(1,0,isb),
     &            'b2tfnb_vbecrbx'//chns)
          call my_out_us(70,nFc,1,dv%vaecrb(1,1,isb),
     &            'b2tfnb_vbecrby'//chns)
       endif
!       call my_out_us(70,nFc,1,dv%uadia(1,0,isb),'b2tfnb_ubdiax'//chns)
!       call my_out_us(70,nFc,1,dv%uadia(1,1,isb),'b2tfnb_ubdiay'//chns)
       call my_out_us(70,nFc,1,dv%vadia(1,0,isb),'b2tfnb_vbdiax'//chns)
       call my_out_us(70,nFc,1,dv%vadia(1,1,isb),'b2tfnb_vbdiay'//chns)
       call my_out_us(70,nFc,1,dv%wadia(1,0,isb),'b2tfnb_wbdiax'//chns)
       call my_out_us(70,nFc,1,dv%wadia(1,1,isb),'b2tfnb_wbdiay'//chns)
       call my_out_us(70,nCv,0,dv%kinrgy,'b2tfnb_kbnrgy'//chns)
      endif

      if(iout.ne.0) then                                                !srv 18.06.08 {
       write (chns,'(i3.3)') isb
       call my_out_us(70,nFc,1,dv%fna(1,0,isb),'b2tfnb_fnbx'//chns)
       call my_out_us(70,nFc,1,dv%fna(1,1,isb),'b2tfnb_fnby'//chns)
       call my_out_us(70,nFc,1,dv%fna_nodrift(1,0,isb),
     &             'b2tfnb_fnb_nodriftx'//chns)
       call my_out_us(70,nFc,1,dv%fna_nodrift(1,1,isb),
     &             'b2tfnb_fnb_nodrifty'//chns)
      call my_out_us(70,nFc,1,dv%fna_fcor(1,0,isb),
     &             'b2tfnb_fnb_fcorx'//chns)
      call my_out_us(70,nFc,1,dv%fna_fcor(1,1,isb),
     &             'b2tfnb_fnb_fcory'//chns)
      call my_out_us(70,nFc,1,dv%fna_he(1,0,isb),'b2tfnb_fnb_hex'//chns)
      call my_out_us(70,nFc,1,dv%fna_he(1,1,isb),'b2tfnb_fnb_hey'//chns)
       call my_out_us(70,nFc,1,dv%fnaPSch(1,0,isb),
     &             'b2tfnb_fnbPSchx'//chns)
       call my_out_us(70,nFc,1,dv%fnaPSch(1,1,isb),
     &             'b2tfnb_fnbPSchy'//chns)
       call my_out_us(70,nFc,1,dv%fna_mdf(1,0,isb),
     &             'b2tfnb_fnb_mdfx'//chns)
       call my_out_us(70,nFc,1,dv%fna_mdf(1,1,isb),
     &             'b2tfnb_fnb_mdfy'//chns)
      call my_out_us(70,nFc,1,dv%fna_32(1,0,isb),'b2tfnb_fnb_32x'//chns)
      call my_out_us(70,nFc,1,dv%fna_32(1,1,isb),'b2tfnb_fnb_32y'//chns)
      call my_out_us(70,nFc,1,dv%fna_52(1,0,isb),'b2tfnb_fnb_52x'//chns)
      call my_out_us(70,nFc,1,dv%fna_52(1,1,isb),'b2tfnb_fnb_52y'//chns)
      call my_out_us(70,nFc,1,dv%fna_32nd(1,0,isb),
     &             'b2tfnb_fnb_32ndx'//chns)
       call my_out_us(70,nFc,1,dv%fna_32nd(1,1,isb),
     &             'b2tfnb_fnb_32ndy'//chns)
       call my_out_us(70,nFc,1,dv%fna_52nd(1,0,isb),
     &             'b2tfnb_fnb_52ndx'//chns)
       call my_out_us(70,nFc,1,dv%fna_52nd(1,1,isb),
     &             'b2tfnb_fnb_52ndy'//chns)
       call my_out_us(70,nFc,1,co%fllim0fna(1,0,isb),
     &             'b2tfnb_fllim0fnbx'//chns)
       call my_out_us(70,nFc,1,co%fllim0fna(1,1,isb),
     &             'b2tfnb_fllim0fnby'//chns)
      call my_out_us(70,nFc,1,dv%vaecrb(1,0,isb),'b2tfnb_vbecrbx'//chns)
      call my_out_us(70,nFc,1,dv%vaecrb(1,1,isb),'b2tfnb_vbecrby'//chns)
       call my_out_us(70,nFc,1,dv%uadia(1,0,isb),'b2tfnb_ubdiax'//chns)
       call my_out_us(70,nFc,1,dv%uadia(1,1,isb),'b2tfnb_ubdiay'//chns)
       call my_out_us(70,nFc,1,dv%vadia(1,0,isb),'b2tfnb_vbdiax'//chns)
       call my_out_us(70,nFc,1,dv%vadia(1,1,isb),'b2tfnb_vbdiay'//chns)
       call my_out_us(70,nFc,1,dv%wadia(1,0,isb),'b2tfnb_wbdiax'//chns)
       call my_out_us(70,nFc,1,dv%wadia(1,1,isb),'b2tfnb_wbdiay'//chns)
       call my_out_us(70,nCv,0,dv%kinrgy,'b2tfnb_kbnrgy'//chns)
       call my_out_us(70,nFc,1,co%cvla(1,0,isb),'b2tfnb_cvlbx'//chns)
       call my_out_us(70,nFc,1,co%cvla(1,1,isb),'b2tfnb_cvlby'//chns)
      call my_out_us(70,nFc,1,co%cvlahz(1,0,isb),'b2tfnb_cvlbhzx'//chns)
      call my_out_us(70,nFc,1,co%cvlahz(1,1,isb),'b2tfnb_cvlbhzy'//chns)
      call my_out_us(70,nFc,1,co%cdpahz(1,0,isb),'b2tfnb_cdpbhzx'//chns)
      call my_out_us(70,nFc,1,co%cdpahz(1,1,isb),'b2tfnb_cdpbhzy'//chns)
       call my_out_us(70,nFc,1,co%cdpa(1,0,isb),'b2tfnb_cdpbx'//chns)
       call my_out_us(70,nFc,1,co%cdpa(1,1,isb),'b2tfnb_cdpby'//chns)
       call my_out_us(70,nFc,1,co%cdna(1,0,isb),'b2tfnb_cdnbx'//chns)
       call my_out_us(70,nFc,1,co%cdna(1,1,isb),'b2tfnb_cdnby'//chns)
       call my_out_us(70,nFc,1,co%cddi(1,0,isb),'b2tfnb_cddbx'//chns)
       call my_out_us(70,nFc,1,co%cddi(1,1,isb),'b2tfnb_cddby'//chns)
       call my_out_us(70,nFc,1,term,'b2tfnb_term'//chns)
       call my_out_us(70,nFc,1,termhz,'b2tfnb_termhz'//chns)
       call my_out_us(70,nFc,1,geo%fcPbs,'b2tfnb_pbsx')
       call my_out_us(70,nFc,1,geo%fcPbshz,'b2tfnb_pbshzx')
       call my_out_us(70,nFc,1,ub2dia(1,0),'b2tfnb_ub2diax'//chns)
       call my_out_us(70,nFc,1,ub2dia(1,1),'b2tfnb_ub2diay'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,0,0),'b2tfnb_conbx0'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,0,1),'b2tfnb_conbx1'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,0,2),'b2tfnb_conbx2'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,0,3),'b2tfnb_conbx3'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,0,4),'b2tfnb_conbx4'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,1,0),'b2tfnb_conby0'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,1,1),'b2tfnb_conby1'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,1,2),'b2tfnb_conby2'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,1,3),'b2tfnb_conby3'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,1,4),'b2tfnb_conby4'//chns)
       call my_out_us(70,nFc,1,dv%flob(1,0),'b2tfnb_flobx'//chns)
       call my_out_us(70,nFc,1,dv%flob(1,1),'b2tfnb_floby'//chns)
       call my_out_us(70,nFc,1,kmprti(1),'b2tfnb_kmprtix'//chns)
       call my_out_us(70,nFc,1,dkmprti(1,0),'b2tfnb_dkmprtix'//chns)
       call my_out_us(70,nFc,1,dkmprti(1,1),'b2tfnb_dkmprtiy'//chns)
       call my_out_us(70,nCv,0,pl%ua(1,isb),'b2tfnb_ua'//chns)
       call my_out_us(70,nCv,0,pl%na(1,isb),'b2tfnb_na'//chns)
       call my_out_us(70,nCv,0,dv%pa(1,isb),'b2tfnb_pa'//chns)
       call my_out_us(70,nCv,0,pl%te,'b2tfnb_te')
       call my_out_us(70,nCv,0,pl%ti,'b2tfnb_ti')
       call my_out_us(70,nCv,0,pl%po,'b2tfnb_po')
       call my_out_us(70,nCv,0,rt%rza(1,isb),'b2tfnb_zb'//chns)
       call my_out_us(70,nFc,1,dv%facdrift,'b2tfnb_facdrift'//chns)
       call my_out_us(70,nFc,1,dv%fac_ExB,'b2tfnb_fac_ExB'//chns)
       call my_out_us(70,nFc,1,dv%fchvispar(1,0),'b2tfnb_fchvisparx')
       call my_out_us(70,nFc,1,dv%fchvispar(1,1),'b2tfnb_fchvispary')
       call my_out_us(70,nFc,1,dv%fchinert(1,0),'b2tfnb_fchinertx')
       call my_out_us(70,nFc,1,dv%fchinert(1,1),'b2tfnb_fchinerty')
       call my_out_us(70,nFc,1,dv%fchanml(1,0),'b2tfnb_fchanmlx')
       call my_out_us(70,nFc,1,dv%fchanml(1,1),'b2tfnb_fchanmly')
      endif                                                             !srv 18.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfnb

      end subroutine b2tfnb

