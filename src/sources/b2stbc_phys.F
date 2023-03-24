*-----------------------------------------------------------------------
*.specification
      subroutine b2stbc_phys (nCv, nFc, nVx, ns, ismain, ismain0,
     &  switch, geo, mpg, pl, dv, co, rt, rtw, st_ext, srw,
     &  wrong_flow, main_call, praverage)
      use b2mod_types
      use b2mod_math
      use b2mod_feedback
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_constants
!      use b2mod_geo_corner
      use b2mod_time
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_external
      use b2mod_boundary_sources                                        !srv 26.02.18
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
#ifdef B25_EIRENE
      use eirmod_extrab25
      use eirmod_braeir
#endif
      use b2mod_astra_to_b2                                          !srv 28.07.11
      use b2mod_ad
     & , only : ncall_b2stbc_phys,
     &          UBIG, VBIG, PBIG, UP00, EPSLON, NE00, US
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : senepar, senipar, sconpar, scurpar,
     &  prev_cur_delta, cur_delta, po_step, po_prev, po_curr, coeff_16,
     *  NUEOUT, FS, VS, NAS, TES, TIS, NIVS, FFAC, DFS, POS, S2
#endif
      use b2mod_subsys
      implicit none
#     include "DIMENSIONS.F"
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain, ismain0
      real (kind=R8) :: praverage(2)
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2RatesWork), intent (in) :: rtw
      type (B2Coeff), intent (inout) :: co
      type (B2StateExt), intent (in) :: st_ext
      logical main_call                                              !xpb
*   ..output arguments (unspecified on entry)
      type (B2SourceWork), intent (inout) :: srw
      logical wrong_flow                                             !xpb
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBC_PHYS provides linearised expressions for the sources of
*     density, parallel momentum, electron heat, atom heat, and charge,
*     due to the boundary conditions.
*     Here the boundary conditions are expressed by use of the old B2
*     namelist format.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iCv1, iCv2, iFc, iVx1, iVx2, is, ib, ibc, ib0,
     &  ism, ibc0, k
#ifdef WG_TODO
     & ,coreregno, coreregn2, ibw, ireg, isb, isa,
     1  bc_type13_ref, bc_type13_ix, bc_type13_iy, bc_type13_i,
     2  bc_type13_ixx,
     1  bc_type16_ref, bc_type16_ix, bc_type16_iy, bc_type16_i,
     2  bc_type16_ixx,
     1  bc_type20_ref, bc_type20_ix, bc_type20_iy, bc_type20_i,
     2  bc_type20_ixx,
     1  bc_type21_ref, bc_type21_ix, bc_type21_iy, bc_type21_i,
     2  bc_type21_ixx,
     1  bc_type16_kinetic_energy
      integer       :: num_inner                                        !srv 26.11.02
      integer, save :: istyle_average_pr = 0, istyle_fchi = 0           !srv 08.08.02 14.02.05
      integer, save :: iform_corr_dia = 0                               !srv 25.07.05
      integer, save :: bcpot_16_step = 50                               !Solovyev 23.05.14 !srv 29.09.15 }
#endif
      integer, save :: bccon_flux_fb_style =0
      integer, save :: bcene_15_style =0                                !srv 19.04.13
      integer, save :: bceni_15_style =0                                !srv 19.04.13
      integer, save :: bcpot_12_style =0                                !srv 19.04.13
      integer, save :: bcenk_17_style =0
      logical once_bccon(0:ns-1), once_bcene, once_bceni, once_bcpot    !srv 16.02.12 12.02.15 29.09.15
      real (kind=R8) ::
     *  vte, s1, s1hz, fche, fchi, t0, t1, t2,                          !srv 15.11.07
     1  t3, t4, totflux, dsh, vt,
     2  vbnd, totpar,
     &  naua2sum(0:ns-1), naua2average(0:ns-1),
     &  pzaverage, pzsum, volsum                                        !srv 24.03.99 18.01.02
#ifdef WG_TODO
     1  ,weight, volave
     &  ,s0, s0hz, ua2sum, ua2average(0:ns-1)                           !srv 07.04.00 15.02.02
     &  ,uasum, uaaverage(0:ns-1), totfluxn, totfluxi                   !srv 15.02.02
     &  ,nesum, neave, niave, tisum, tiave, tiaverage, tesum, teave     !srv 15.12.00 25.01.02
     &  ,nasum, naave, naaverage(0:ns-1), poave, bbsum, bbaverage       !srv 22.01.02 18.10.04
      real (kind=R8) :: bzhzsum, bzhzaverage
#endif
      real (kind=R8) :: bzsum, bzaverage                                !lk 20.11.07
      real (kind=R8) ::
     *  wdia(nCv), wrk0(nCv), wrkf(nFc),
#ifdef WG_TODO
     *  uu(nFc,0:ns+st_ext%ns-1), zarg,
     *  ud(nFc,0:ns-1), !sx(-1:nx,-1:ny),
#endif
     *  pz(nCv), !sxhz(-1:nx,-1:ny),                                    !srv 24.08.08
     *  gOnedBsq(nFc,0:1),
     &  sna0_no_mdf(nCv,0:1,0:ns-1)
      real (kind=R8), save :: cor9 = 0.0_R8                             !srv 23.09.08
      character chns*3, chk*1                                           !srv 11.07.05
      logical decay_length_ok, ishigh
#ifdef WG_TODO
      logical fedback_species
      logical, save :: BCCON14_W(NBCD), BCCON14_E(NBCD)                 !srv 20.09.17
      character*24 hlp_frm
      real (kind=R8) ::
     *  upout, uuout, uutest, zbar
#endif
      real (kind=R8) ::
     *  rz(nCv), cs, csb,                                               !srv 01.02.09
     *  fchy_dia, fchy_inertia, seec                                    !srv 07.04.00
      character*26 boundary_location(nbcd)
#ifdef WG_TODO
      integer isfeedback, fix_fch_in_fhe_sheath, is_kinetic
      integer, save :: fix_sx_in_fhe_sheath = 1                         !srv 08.05.09
      save coreregno, coreregn2,  bc_type13_ref,
     1  bc_type16_ref, bc_type20_ref,
     2  bc_type21_ref,
     *  isfeedback, fix_fch_in_fhe_sheath,
     *  is_kinetic
#endif
*   ..procedures
      intrinsic abs, sqrt
#ifdef WG_TODO
      real (kind=R8) :: bcint
      external bcint
#endif
      real (kind=R8) ::
     *  sy, intvertex_s                                                 !srv 24.08.08
      external xertst, ipgeti, eseec0
      external ipgetr, b2xvsg, b2xvff, b2xvps
*   ..initialisation
#ifdef WG_TODO
      integer ibnd
      integer, save :: conduction_only, no_hybr_e, no_hybr_i
      real (kind=R8) ::
     .  Tguess, fy, ehy, fylimt, prefactor, w1,
     .  fyt, gam0, gam1, tlimit, cone0, coni0, flo,
     .  cum_flo_vadia, cum_flo_vaecrb, flo_vedia, flo_veecrb
      data coreregno/1/, coreregn2/4/, bc_type13_ref/1/,
     & bc_type16_ref/1/, bc_type20_ref/1/,
     & isfeedback/0/, conduction_only/0/, no_hybr_e/0/, no_hybr_i/0/,
     & fix_fch_in_fhe_sheath/2/, bc_type16_kinetic_energy/0/,
     & bc_type21_ref/1/, is_kinetic/0/
      real (kind=R8) :: harm, upwind
      external harm, upwind
#endif
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbc_phys')
*   ..set internal parameters on first call
      if (ncall_b2stbc_phys.eq.0) then
       call ipgeti ('b2stbc_bccon_flux_fb_style', bccon_flux_fb_style)
       call ipgeti ('b2stbc_bcene_15_style', bcene_15_style)
       call ipgeti ('b2stbc_bceni_15_style', bceni_15_style)
       call ipgeti ('b2stbc_bcpot_12_style', bcpot_12_style)
       call ipgeti ('b2stbc_bcenk_17_style', bcenk_17_style)
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c       call ipgeti ('b2tfhe_conduction_only', conduction_only)
c       call ipgeti ('b2stbc_iform_corr_dia', iform_corr_dia)        !srv 25.07.05
c       call ipgeti ('b2stbc_coreregno', coreregno)
c       call xertst (cbirso.le.coreregno.and.
c     1   coreregno.le.cbirso+cbnrso-1,
c     2   'faulty internal parameter coreregno')
c       if (nncut.ge.2) then
c         call ipgeti ('b2stbc_coreregn2', coreregn2)
c         call xertst (cbirso.le.coreregn2.and.
c     1     coreregn2.le.cbirso+cbnrso-1,
c     2     'faulty internal parameter coreregn2')
c       endif
c       call ipgeti ('b2tfhe_no_hybr', no_hybr_e)
c       call ipgeti ('b2tfhi_no_hybr', no_hybr_i)
c       call ipgeti ('b2stbc_type13_ref', bc_type13_ref)
c       call ipgeti ('b2stbc_type16_ref', bc_type16_ref)
c       call ipgeti ('b2stbc_type20_ref', bc_type20_ref)
c       call ipgeti ('b2stbc_type21_ref', bc_type21_ref)
c       call ipgeti ('b2stbc_fix_fch_in_fhe_sheath',
c     1  fix_fch_in_fhe_sheath)
c       call ipgeti ('b2stbc_fix_sx_in_fhe_sheath',
c     1  fix_sx_in_fhe_sheath)                                       !srv 08.05.09
c       call ipgeti ('b2stbc_istyle_fchi', istyle_fchi)              !srv 14.02.05
c       call ipgeti ('b2stbc_bcpot_16_step', bcpot_16_step)          !Solovyev 23.05.14 !srv 29.09.15
c       call ipgeti ('b2stbc_type16_kinetic_energy',
c     1  bc_type16_kinetic_energy)
c       call ipgeti ('b2stbc_average_pr', istyle_average_pr)         !srv 08.08.02
        call ipgetr ('b2stbc_BC2_cor9', cor9)
      endif
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
      call xertst (0.le.ismain0.and.ismain0.lt.ns,
     . 'invalid main neutral species index ismain0')
      call xertst (is_neutral(ismain0).or.ismain.eq.ismain0,
     . 'invalid main neutral species ismain0; must be neutral')
*   ..test facdrift (not on first call)                                  !xpb
      if (ncall_b2stbc_phys.gt.0)                                                    !xpb
     > call xertst (0.0_R8.le.minval(dv%facdrift) .and.
     &       1.0_R8.ge.maxval(dv%facdrift), 'faulty argument facdrift')
*   ..extensive tests on first few calls
      if (ncall_b2stbc_phys.lt.3) then
*    ..test sign of vol
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
       call b2xvsg (nFc*2, geo%fcVol, 1, 'vol', '.gt.')
*    ..test state
       call b2xvps (nCv, nFc, ns, pl, dv)
*    ..test sign of chce, chci
       wrkf(:) = co%chce(:,0)*geo%fcQalf(:,0)
       call b2xvsg (nFc,wrkf,1,'chce0','.ge.')
       wrkf(:) = co%chce(:,1)*geo%fcQalf(:,1)
       call b2xvsg (nFc,wrkf,1,'chce1','.ge.')
       wrkf(:) = co%chci(:,0)*geo%fcQalf(:,0)
       call b2xvsg (nFc,wrkf,1,'chci0','.ge.')
       wrkf(:) = co%chci(:,1)*geo%fcQalf(:,1)
       call b2xvsg (nFc,wrkf,1,'chci1','.ge.')
      endif
*   ..set cutll, cutlo
      if (.not.b2mod_math_initialised) call b2mod_math_init

*  -- INITIALISE CERTAIN CONSTANTS --
      EPSLON=1.e-10_R8
*  -- SET NE00, UP00, UBIG, VBIG --
      NE00=1.0e20_R8
      UP00=SQRT(100*EV/MP)
      UBIG=UP00 / EPSLON
      VBIG=1.0e-2_R8*UBIG
      PBIG=UP00

      sna0_no_mdf = 0.0_R8

      call b2xppz (nCv, ns, dv%ne, pl%na, pl%te, pl%ti, pz, st_ext)!srv 18.01.02

*   ..pre-compute some quantities for core BCs with drifts
      call grad(nCv, nFc, nVx, 1, geo, mpg, geo%cvOnedBsq,
     &  geo%vxOnedBsq, gOnedBsq)
      pzsum     = 0.0_R8
      naua2sum  = 0.0e0_R8
      volsum    = 0.0_R8
      bzsum     = 0.0_R8
      do ib = 1, nbc
       ! identify core boundaries; assume all core boundaries
       ! form single flux surface
       iCv1 = mpg%bcCv(mpg%bcCvP(ib,1),1)
       if (mpg%cvOnClosedSurface(iCv1)) then
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
         iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
         iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
*      ..accumulate the average pressure for inner flux surface
         pzsum    = pzsum  + geo%fcS(iFc)*pz(iCv1)
         volsum   = volsum + geo%fcS(iFc)
*      ..accumulate the average na*ua**2, na for inner flux surface
         do is = 0, ns-1
           naua2sum(is)  = naua2sum(is) + geo%fcS(iFc)*
     .       pl%na(iCv1,is)*pl%ua(iCv1,is)**2
         enddo
*      ..accumulate the average fcBb(:,2) for inner flux surface
         bzsum = bzsum + geo%fcS(iFc)*geo%fcBb(iFc,2)
        enddo
       endif
      enddo
      if (volsum.gt.0.0_R8) then
        pzaverage    = pzsum/volsum
        naua2average = naua2sum/volsum
        bzaverage    = bzsum/volsum
      endif

#ifdef WG_TODO
cwdk code to be removed after full adaption to unstructured solver...

      if (istyle_average_pr .eq.0) then                           !srv 08.08.02
       iy = -1                                                    !srv 24.03.99
       pzsum =0.0e0_R8                                            !srv 24.03.99 18.01.02
       volsum=0.0e0_R8                                            !srv 24.03.99
       do ix = -1, nx                                             !srv 26.11.02
        if (on_closed_surface(ix,iy)) then                        !srv 26.11.02
         pzsum=pzsum + pz(ix,iy)*vol(ix,iy)                       !srv 24.03.99 18.01.02
         volsum=volsum + vol(ix,iy)                               !srv 24.03.99
        endif                                                     !srv 26.11.02
       enddo                                                      !srv 24.03.99
       if (volsum.gt.0.0_R8) pzaverage=pzsum/volsum               !srv 24.03.99 18.01.02
      else                                                        !srv 08.08.02
       iy = 0                                                     !srv 24.03.99  07.06.02 {
       pzsum =0.0e0_R8                                            !srv 24.03.99 18.01.02
       volsum=0.0e0_R8                                            !srv 24.03.99
       num_inner=0                                                !srv 26.11.02
       do ix = -1, nx                                             !srv 26.11.02
        if (on_closed_surface(ix,iy)) then                        !srv 26.11.02
         pzsum=pzsum
     *   + (pz(ix,iy)+pz(bottomix(ix,iy),bottomiy(ix,iy)))        !srv 24.03.99 18.01.02
         volsum=volsum + vol(ix,iy)                               !srv 24.03.99
         num_inner=num_inner+1                                    !srv 26.11.02
        endif                                                     !srv 26.11.02
       enddo                                                      !srv 24.03.99
       if (num_inner.gt.0) pzaverage=pzsum*0.5e0_R8/num_inner     !srv 24.03.99 18.01.02 07.06.02 } 26.11.02
      endif                                                       !srv 08.08.02
*   ..find the average ne, te, ti, bb for inner flux surface      !srv 15.12.00 18.10.04
      bbsum = 0.0e0_R8                                            !srv 18.10.04
      nesum = 0.0e0_R8                                            !srv 15.12.00 25.01.02
      tesum = 0.0e0_R8                                            !srv 15.12.00
      tisum = 0.0e0_R8                                            !srv 15.12.00
      bzhzsum  = 0.0e0_R8                                         !lk 20.11.07
      bbaverage  = 0.0e0_R8                                       !lk 20.11.07
      do ix = -1, nx                                              !srv 26.11.02
       if (on_closed_surface(ix,iy)) then                         !srv 26.11.02
        bbsum=bbsum + bb(ix,iy,3)*vol(ix,iy)                      !srv 18.10.04
        nesum=nesum + ne(ix,iy)*vol(ix,iy)                        !srv 15.12.00 25.01.02
        tisum=tisum + ti(ix,iy)*vol(ix,iy)                        !srv 15.12.00
        tesum=tesum + te(ix,iy)*vol(ix,iy)                        !srv 15.12.00
        bzhzsum = bzhzsum+bb(ix,iy,2)*vol(ix,iy)**2/gs(ix,iy,2)   !lk 20.11.07 !srv 07.05.13
       endif                                                      !srv 26.11.02
      enddo                                                       !srv 15.12.00
      if (volsum.gt.0.0_R8) then
       bbaverage=bbsum/volsum                                     !srv 18.10.04
       tiaverage=tisum/volsum                                     !srv 15.12.00
       bzhzaverage=bzhzsum/volsum                                 !lk 20.11.07
      end if
*   ..find the average ua**2, na for inner flux surface           !srv 22.01.02 {
      do is = 0, ns-1
       uasum  = 0.0e0_R8                                           !srv 15.02.02
       ua2sum = 0.0e0_R8                                           !srv 07.04.00 15.02.02
       nasum  = 0.0e0_R8
       do ix = -1,nx                                               !srv 26.11.02
        if (on_closed_surface(ix,iy)) then                         !srv 26.11.02
         ua2sum = ua2sum + ua(ix,iy,is)**2*vol(ix,iy)              !srv 07.04.00 15.02.02
         uasum  = uasum  + ua(ix,iy,is)*vol(ix,iy)                 !srv 15.02.02
         nasum  = nasum  + na(ix,iy,is)*vol(ix,iy)
        endif                                                      !srv 26.11.02
       enddo
       ua2average(is) = ua2sum/volsum                              !srv 07.04.00 15.02.02
       uaaverage(is) = uasum/volsum                                !srv 07.04.00 15.02.02 16.04.02
       naaverage(is) = nasum/volsum
      enddo                                                        !srv 22.01.02 }
#endif
cwdk  ro and pr based on charged particles only!
cwdk  no neutral contribution to sound speed definition
      call b2xprz (nCv, ns, mp, am, pl%na, rz, st_ext)             !srv 01.02.09
!      call b2xppz (nCv, ns, dv%ne, pl%na, pl%te, pl%ti, pz, st_ext)
      wdia=(dv%ne*pl%te+dv%ni(:,1)*pl%ti)*geo%cvBb(:,2)
      do is = 0, st_ext%ns-1
        if (.not.st_ext%is_neutral(is)) wdia=wdia+
     .   st_ext%na(:,is)*st_ext%ta(:,is)*geo%cvBb(:,2)
      end do
#ifdef WG_TODO
cwdk do we want to keep this kind of back-mapping to SOLPS4???
      do iy = -1, ny
        do ix = -1, nx
          if (rightix(ix,iy).ne.nx+1) then
            sx(ix,iy)=gs(rightix(ix,iy),rightiy(ix,iy),0)
            if (area_fix.ge.3)
     >       sx(ix,iy) = sx(ix,iy)*qc(rightix(ix,iy),rightiy(ix,iy),0)
            sxhz(ix,iy)=gs(rightix(ix,iy),rightiy(ix,iy),0)*       !srv 24.08.08
     &                  hz(ix,iy)
!     &                  hz(rightix(ix,iy),rightiy(ix,iy))
            if (area_fix.ge.3)
     >       sxhz(ix,iy) = sxhz(ix,iy)*
     >         qc(rightix(ix,iy),rightiy(ix,iy),0)
          else
            sx(ix,iy)=0.0_R8
            sxhz(ix,iy)=0.0_R8                                     !srv 24.08.08
          endif
          if (leftix(ix,iy).ne.-2) then
            do is = 0, ns-1
              uu(ix,iy,is)=fna(ix,iy,0,0,is)/gs(ix,iy,0)/
     .         ((na(ix,iy,is)*hx(leftix(ix,iy),leftiy(ix,iy))+
     .           na(leftix(ix,iy),leftiy(ix,iy),is)*hx(ix,iy))/
     .          (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy))))
            end do
            do is = 0, ns_ext-1
              uu(ix,iy,is+ns)=fa_ext(ix,iy,0,0,is)/gs(ix,iy,0)/
     .         ((na_ext(ix,iy,is)*hx(leftix(ix,iy),leftiy(ix,iy))+
     .           na_ext(leftix(ix,iy),leftiy(ix,iy),is)*hx(ix,iy))/
     .          (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy))))
            end do
            if (area_fix.ge.1) uu(ix,iy,:) = uu(ix,iy,:)/qc(ix,iy,0)
          else
            uu(ix,iy,:)=0.0_R8
          endif
        enddo
      enddo
cxpb  Following suggestion from Russians, the drift velocity
cxpb  should only include ExB and not diamagnetic contribution
cwdk  With redefinition of drifts as being the poloidal projections,
cwdk  now need to recompute full diamagnetic component for ud
      if (isymm.eq.2.or.isymm.eq.4) then
        ud(:,:) = 0.0_R8
      else if (switch%sheath_drift_fix.eq.0) then
        do is = 0,ns-1
         ud(:,is) = dv%uadia(:,0,is)*geo%fcBb(:,3)/geo%fcBb(:,2)*signmf
        enddo
      else
        do is=0,ns-1
         ud(:,is) = dv%vaecrb(:,0,is)*geo%fcBb(:,3)/geo%fcBb(:,2)*signmf
        enddo
      endif
      if (main_call .and. switch%b2stbc_diagno.ge.1) then             !xpb
        write(*,'(a,1p,(10000e12.4))') 'dna0',
     1   (dna0(nx/2,ny,is),is=0,ns-1)
        write(*,'(a,1p,(10000e12.4))') 'dpa0',
     1   (dpa0(nx/2,ny,is),is=0,ns-1)
        write(*,'(a,1p,(10000e12.4))') 'vsa0',
     1   (vsa0(nx/2,ny,is),is=0,ns-1)
        write(*,'(a,1p,(10000e12.4))') 'hce0',
     1   hce0(nx/2,ny)
        write(*,'(a,1p,(10000e12.4))') 'hci0',
     1   hci0(nx/2,ny)
        write(*,'(a,1p,(10000e12.4))') 'sig0',
     1   sig0(nx/2,ny)
      endif                                                           !xpb
#endif
      if(ncall_b2stbc_phys.eq.0) write(*,*) 'Doing B2-style bc'

      once_bccon = .false.                                            !srv 16.02.12 {
      once_bcene = .false.
      once_bceni = .false.                                            !srv 16.02.12 }
      once_bcpot = .false.                                            !Solovyev 13.05.14 !srv 29.09.15
* -- Prepare boundary_location array in advance                       !srv 05.08.22 {
      DO IB = 1, NBC
        if (bcpos(ib).eq.-2) then
          write(boundary_location(ib),'(a23,i2)')
     .     ' on boundary cell list ',IB
        else
          write(boundary_location(ib),'(a4,i4,a6,i4,a4,i4)')
     .     ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
        endif
      ENDDO                                                           !srv 05.08.22 }

* -- LOOP OVER ALL BOUNDARIES
      DO IB = 1, NBC

* -- Moving core tally according to B.C.
* -- Must be done only once!

#ifdef WG_TODO
        IF (ncall_b2stbc_phys.eq.0 .and. bc_type16_ref.gt.1) THEN
          IF (BCCHAR(IB).EQ.'S' .AND.
     &       (BCENE(IB).EQ.16.or.BCENE(IB).EQ.17.or.BCENE(IB).EQ.20.or.
     &        BCENI(IB).EQ.16.or.BCENI(IB).EQ.17.or.BCENI(IB).EQ.20))
     &     THEN
            DO IX = 1, bc_list_size(ib)
              bc_type16_ix=bc_list_x(ix,ib)
              bc_type16_iy=bc_list_y(ix,ib)
              do bc_type16_i=1,bc_type16_ref      ! how many steps away from the boundary?
                bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
                bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
                bc_type16_ix=bc_type16_ixx
                call xertst (-2.lt.bc_type16_ix.and.
     &                       bc_type16_ix.lt.nx+1,
     1           'stepped off x-grid in type 16 bc')
                call xertst (-2.lt.bc_type16_iy.and.
     &                       bc_type16_iy.lt.ny+1,
     1           'stepped off y-grid in type 16 bc')
              enddo
              region(bc_type16_ix,bc_type16_iy,2) =
     .         region(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2)
              region(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2) = 0
            ENDDO
            if (bc_type16_ref.gt.2) then
              write(*,'(a,i2,a)')
     .         'Moved core radial flux tally position by ',
     .          bc_type16_ref-1,' radial steps.'
            else if (bc_type16_ref.eq.2) then
              write(*,'(a)')
     .         'Moved core radial flux tally position by '//
     .         '1 radial step.'
            endif
          ENDIF
        ENDIF
#endif
        DO IS=0, NS-1
        if (is.eq.ns-1) then
          ishigh = .true.
        else
          ishigh = .not.lnext(is,is+1)
        end if

* -- BOUNDARY CONDITIONS FOR CONTINUITY EQUATIONS --
        select case (BCCON(IS,IB))

        case (0)
          continue

* -- BCCON=1 -- PRESCRIBE THE VALUE OF THE DENSITY

        case (1)
          if(ncall_b2stbc_phys.eq.0) then
            if (CONPAR(IS,IB,1).eq.0.0_R8) then
              CONPAR(IS,IB,1) = switch%b2mndr_na_min
              write(*,*) 'Zero density B.C. no longer allowed!'
              write(*,*) 'Changed to na_min.'
            endif
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCCON =  1 : constant density ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1       'BCCON = 1, CONPAR(,,1) !> 0')
            icv = mpg%bcCv(mpg%bcCvP(ib,1),1)
            if (mpg%cvOnClosedSurface(icv)) then
             if(is.eq.ismain.and.switch%use_astra.ne.0) then
              CONPAR(IS,IB,1) = neb_astra                               !srv 28.07.11
             endif
             if (ishigh) call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 1')
            ENDIF
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%sna0(iCv1,0,is) =
     .        VBIG*geo%fcS(iFc)*conpar(is,ib,1)
            srw%sna0(iCv1,1,is) =
     .       -VBIG*geo%fcS(iFc)
            sna0_no_mdf(iCv1,0,is) = srw%sna0(iCv1,0,is)
            sna0_no_mdf(iCv1,1,is) = srw%sna0(iCv1,1,is)
          enddo

* -- BCCON=2 -- PRESCRIBE THE GRADIENT OF THE DENSITY

        case (2)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCCON =  2 : density gradient ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            icv = mpg%bcCv(mpg%bcCvP(ib,1),1)
            if (mpg%cvOnClosedSurface(icv).and.ishigh)
     >       call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 2')
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                    ! number of the corresponding domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            iVx1 = mpg%fcVx(iFc,1)
            iVx2 = mpg%fcVx(iFc,2)
            t1 = intvertex_s(iVx1,nCv,nVx,mpg,geo%vxVol,pl%na(:,is))
            t2 = intvertex_s(iVx2,nCv,nVx,mpg,geo%vxVol,pl%na(:,is))
            srw%sna0(iCv1,0,is)=  VBIG*geo%fcS(iFc)*
     .         (conpar(is,ib,1)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2))*
     .          geo%fcQgam(iFc,0) - cor9*(t2 - t1)*geo%fcQgam(iFc,1)*
     .          (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     .          mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1) + pl%na(iCv2,is))
            srw%sna0(iCv1,1,is)= -VBIG*geo%fcS(iFc)
            sna0_no_mdf(iCv1,0,is) = srw%sna0(iCv1,0,is)
            sna0_no_mdf(iCv1,1,is) = srw%sna0(iCv1,1,is)
          enddo

* -- BCCON=3 -- SHEATH CONDITIONS

        case (3)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,a,a,a,i3)')
     .       'BCCON =  3 : sheath density bc on ',BCCHAR(IB),
     1       boundary_location(ib),
     2       ' for species ',IS
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            srw%sna0(iCv1,0,is) =
     .        UBIG*geo%fcS(iFc)*pl%na(iCv2,is)
            srw%sna0(iCv1,1,is) =
     .       -UBIG*geo%fcS(iFc)
            sna0_no_mdf(iCv1,0,is) = srw%sna0(iCv1,0,is)
            sna0_no_mdf(iCv1,1,is) = srw%sna0(iCv1,1,is)
          enddo

* -- BCCON=4 -- PRESCRIBE THE VALUE OF THE DENSITY, WEAKLY
*                 A MIXED BOUNDARY CONDITION

        case (4)
          call xerrab ('b2stbc_phys -- BCCON = 4 not adapted for WG')
          if(ncall_b2stbc_phys.eq.0) then
            if (CONPAR(IS,IB,1).eq.0.0_R8) then
              CONPAR(IS,IB,1) = switch%b2mndr_na_min
              write(*,*) 'Zero density B.C. no longer allowed!'
              write(*,*) 'Changed to na_min.'
            endif
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCCON =  4 : weak density ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1       'BCCON = 4, CONPAR(,,1) !> 0')
            call xertst(CONPAR(IS,IB,2).gt.0.0_R8,
     1       'BCCON = 4, CONPAR(,,2) !> 0')
            icv = mpg%bcCv(mpg%bcCvP(ib,1),1)
            if (mpg%cvOnClosedSurface(icv).and.ishigh)
     .       call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 4')
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%sna0(iCv1,0,is) =
     .        conpar(is,ib,2)*geo%fcS(iFc)*conpar(is,ib,1)
            srw%sna0(iCv1,1,is) =
     .       -conpar(is,ib,2)*geo%fcS(iFc)
            sna0_no_mdf(iCv1,0,is) = srw%sna0(iCv1,0,is)
            sna0_no_mdf(iCv1,1,is) = srw%sna0(iCv1,1,is)
          enddo

* -- BCCON=5 -- PRESCRIBE THE PARTICLE FLUX PER UNIT AREA

        case (5)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCCON =  5 : particle flux per area ',CONPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
cwdk      Possible backwards-compatibility issue:
cwdk      In the unstructured solver, we impose flux per unit area normal
cwdk      to the surface/boundary.
cwdk      In the original code, the flux per unit area in the poloidal
cwdk      direction was imposed at E/W boundaries (and normal to the
cwdk      radial direction at N/S - unchanged for aligned grids)
          do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            srw%sna0(iCv1,0,is) = conpar(is,iB,1)*geo%fcS(iFc)
            srw%sna0(iCv1,1,is) = 0.0_R8
            sna0_no_mdf(iCv1,0,is) = srw%sna0(iCv1,0,is)
            sna0_no_mdf(iCv1,1,is) = srw%sna0(iCv1,1,is)
            if (switch%mdf_fnb.ne.0) then
              !wdk Added correction to make mdf_fnb options equivalent
              t0 = mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &          ((dv%fnaPSch(iFc,0,is) + dv%fnaPSch(iFc,1,is)) +
     &           (geo%fcS(iFc)*geo%fcQalf(iFc,0)*dv%vadia(iFc,0,is) +
     &            geo%fcS(iFc)*geo%fcQalf(iFc,1)*dv%vadia(iFc,1,is))*
     &           0.5_R8*(pl%na(iCv1,is)+pl%na(iCv2,is)))
              srw%sna0(iCv1,0,is) = srw%sna0(iCv1,0,is) +
     &          max(0.0_R8,t0)
              srw%sna0(iCv1,1,is) = srw%sna0(iCv1,1,is) +
     &          min(0.0_R8,t0)/pl%na(iCv1,is)
            endif
          enddo

* -- BCCON=6 -- PRESCRIBE THE TOTAL PARTICLE FLUX FOR A
*                 CONSTANT DENSITY

        case (6)
          call xerrab ('b2stbc_phys -- BCCON = 6 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
           write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .      'BCCON =  6 : total particle flux ',CONPAR(IS,IB,1),
     1      ' with const density on ',BCCHAR(IB),
     2      boundary_location(ib),' for species ',IS
           icv = mpg%bcCv(mpg%bcCvP(ib,1),1)
           if (mpg%cvOnClosedSurface(icv).and.ishigh)
     >      call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 6')
          endif
          IF(BCCHAR(IB).EQ.'N') THEN
            s1=0.0_R8
            s2=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              s1=s1+SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                 BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/
     .         (0.5_R8*
     .         (HY1(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     .          HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
              s2=s2+SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                 BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .          na(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .             BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS)/
     .         (0.5_R8*
     .         (HY1(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     .          HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
            ENDDO
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)= VBIG*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .         (CONPAR(IS,IB,1)/dna0(bc_list_x(ix,ib),
     .                               bc_list_y(ix,ib),is) + s2) / s1
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=-VBIG*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            end do
          ELSEIF(BCCHAR(IB).EQ.'S') THEN
            s1=0.0_R8
            s2=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              s1=s1+SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .         (0.5_R8*(HY1(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                      TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     .                  HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
              s2=s2+SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         na(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS)/
     .         (0.5_R8*(HY1(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                      TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     .                  HY1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))
            ENDDO
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .          VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         (CONPAR(IS,IB,1)/dna0(bc_list_x(ix,ib),
     .                               bc_list_y(ix,ib),is) + s2) / s1
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .          - VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            end do
          ELSE
            CALL XERRAB ('B2STBC -- BCCON=6, NOT AVAILABLE FOR W,E')
          ENDIF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                     !srv 23.09.08 }
#endif

* -- BCCON=7 -- PRESCRIBE THE GIVEN PROFILE OF DENSITY

        case (7)
          if(ncall_b2stbc_phys.eq.0) then                              !srv 26.02.18 {
            write(*,'(a,a,a,a,a,i3)')
     .       'profile of density ',' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            srw%sna0(iCv1,0,is) =
     .        VBIG*geo%fcS(iFc)*bv_na(iCv1,is)
            srw%sna0(iCv1,1,is) =
     .       -VBIG*geo%fcS(iFc)
            sna0_no_mdf(iCv1,0,is) = srw%sna0(iCv1,0,is)
            sna0_no_mdf(iCv1,1,is) = srw%sna0(iCv1,1,is)
          enddo

* -- BCCON=8 -- PRESCRIBE THE TOTAL PARTICLE FLUX WITH
*                 CONSTANT FLUX DENSITY

        case (8)
          if(ncall_b2stbc_phys.eq.0) then
            if(is.eq.ismain.and.BCCHAR(IB).eq.'S'
     &                                   .and.switch%use_astra.ne.0)
     >       CONPAR(IS,IB,1) = fneb_astra                               !srv 04.10.11
            write(*,'(a,1p,g14.7,a4,a1,a,a13,i3)')
     .       'BCCON =  8 : total particle flux ', CONPAR(IS,IB,1),
     1       ' on ', BCCHAR(IB), boundary_location(ib),
     2       ' for species ',IS
          endif
cwdk      Note: take full area for face, whether "poloidal" or "radial" or mixed
cwdk      Backwards compatibility issue for E/W and 5-pt stencil?
cwdk      Total flux will be fine, possibly distributed somewhat differently...
          US = 0.0_R8
          do ibc = 1, mpg%bcCvP(ib,2)
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            US   = US + geo%fcS(iFc)
          enddo
          do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            srw%sna0(iCv1,0,is) = srw%sna0(iCv1,0,is) +
     .        conpar(is,iB,1)/US*geo%fcS(iFc)
            srw%sna0(iCv1,1,is) = 0.0_R8
            sna0_no_mdf(iCv1,0,is) = sna0_no_mdf(iCv1,0,is) +
     .        srw%sna0(iCv1,0,is)
            sna0_no_mdf(iCv1,1,is) = sna0_no_mdf(iCv1,1,is) +
     .        srw%sna0(iCv1,1,is)
            if (switch%mdf_fnb.ne.0) then
              !wdk Correction now involves vadia instead of ExB velocity
              !wdk This makes implementation equivalent to mdf_fnb = 0 option
              t0 = mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &          ((dv%fnaPSch(iFc,0,is) + dv%fnaPSch(iFc,1,is)) +
     &           (geo%fcS(iFc)*geo%fcQalf(iFc,0)*dv%vadia(iFc,0,is) +
     &            geo%fcS(iFc)*geo%fcQalf(iFc,1)*dv%vadia(iFc,1,is))*
     &           0.5_R8*(pl%na(iCv1,is)+pl%na(iCv2,is)))
              srw%sna0(iCv1,0,is) = srw%sna0(iCv1,0,is) +
     &          max(0.0_R8,t0)
              srw%sna0(iCv1,1,is) = srw%sna0(iCv1,1,is) +
     &          min(0.0_R8,t0)/pl%na(iCv1,is)
            endif
          enddo

* -- BCCON=9 -- PRESCRIBE THE DECAY LENGTH FOR THE DENSITY

        case (9)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCCON =  9 : density decay length ',CONPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).ne.0.0_R8,
     1       'B2STBC -- BCCON=9, CONPAR(,,1) = 0')
            icv = mpg%bcCv(mpg%bcCvP(ib,1),1)
            if (mpg%cvOnClosedSurface(icv).and.ishigh)
     >       call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 9')
          endif

          do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            sy = geo%fcS(iFc)

            srw%sna0(iCv1,0,is) = 0.0_R8
            srw%sna0(iCv1,1,is) = -(co%dna0(iCv1,is)+co%dpa0(iCv1,is)*
     &             (rt%rza(iCv1,is)*pl%te(iCv1)+pl%ti(iCv1)))/
     &             conpar(is,ib,1)*sy
            if (switch%mdf_fnb.eq.1) then
              t0 = mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     .          (dv%fnaPSch(iFc,0,is) + dv%fnaPSch(iFc,1,is))
              if(t0.ge.0.0_R8) then
                srw%sna0(iCv1,0,is) = srw%sna0(iCv1,0,is) + t0
              else
                srw%sna0(iCv1,1,is) = srw%sna0(iCv1,1,is) +
     .            t0/pl%na(iCv1,is)
              endif
              t1 = -mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &          (geo%fcS(iFc)*geo%fcQalf(iFc,0)*dv%vaecrb(iFc,0,is) +
     &           geo%fcS(iFc)*geo%fcQalf(iFc,1)*dv%vaecrb(iFc,1,is))*
     &          0.5_R8*(pl%na(iCv1,is) + pl%na(iCv2,is))
              if(t1.ge.0.0_R8) then
                srw%sna0(iCv1,0,is) = srw%sna0(iCv1,0,is) + t1
              else
                srw%sna0(iCv1,1,is) = srw%sna0(iCv1,1,is) +
     &            t1/pl%na(iCv1,is)
              endif
            endif
          enddo
#ifdef WG_TODO

          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=0.0_R8
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     1         -(dna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)+
     2           dpa0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)*
     3           (rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .             te(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .             ti(bc_list_x(ix,ib),bc_list_y(ix,ib))))/
     2         CONPAR(IS,IB,1)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     .         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)            !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)            !srv 23.09.08
              if(switch%mdf_fnb.ne.0) then                                    !srv 08.09.08 {
                t0 =
     .           -fnaPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1,is)
                if (t0.ge.0.0_R8) then
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     .             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) + t0
                else
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) + t0/
     &               na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                end if
                t1 =
     .           vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1,is)
                if (t1.le.0.0_R8) then
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     .             vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .             hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*t1
                else
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     .             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+
     .             vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .             hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*t1*
     .              na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                end if
              endif                                                    !srv 08.09.08 }
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=0.0_R8
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     1         -(dna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)+
     2           dpa0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)*
     3           (rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .             te(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .             ti(bc_list_x(ix,ib),bc_list_y(ix,ib))))/
     2         CONPAR(IS,IB,1)*
     2         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     .         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)            !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     .         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)            !srv 23.09.08
              if(mdf_fnb.ne.0) then                                    !srv 08.09.08 {
                t0 = fnaPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1,is)
                if (t0.ge.0.0_R8) then
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) + t0
                else
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) + t0/
     &               na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                end if
                t1 = -vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1,is)
                if (t1.le.0.0_R8) then
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &             vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &             hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*t1
                else
                  sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &             sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+
     &             vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &             hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*t1*
     &              na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                end if
              endif                                                    !srv 08.09.08 }
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCCON=9, NOT AVAILABLE FOR W,E')
          ENDIF

#endif
* -- BCCON=10 -- leakage option for density, recommended for cases with drifts

        case (10)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCCON = 10 : leakage ',CONPAR(IS,IB,1),' on ',BCCHAR(IB),
     1       boundary_location(ib),
     2       ' for species ',IS
            icv = mpg%bcCv(mpg%bcCvP(ib,1),1)
            if (mpg%cvOnClosedSurface(icv).and.ishigh)
     >       call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 10')
          endif
cwdk      Note: now total face area taken. Need to check backwards
cwdk      compatibility for E/W boundaries
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                    ! number of the plasma cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%sna0(iCv1,0,is) = 0.0_R8
            srw%sna0(iCv1,1,is) = geo%fcS(iFc)*
     &         sqrt((rt%rza(iCv1,is)*pl%te(iCv1)+pl%ti(iCv1))/
     &         (am(is)*mp))*conpar(is,ib,1)
            sna0_no_mdf(iCv1,:,is) =
     &         sna0_no_mdf(iCv1,:,is) + srw%sna0(iCv1,:,is)
            if (switch%mdf_fnb.eq.0) then
              srw%sna0(iCv1,1,is) = srw%sna0(iCv1,1,is) -
     &          max(0.0e0_R8,mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &          (geo%fcS(iFc)*geo%fcQalf(iFc,0)*
     &          (dv%vadia(iFc,0,is) + dv%vaecrb(iFc,0,is)) +
     &           geo%fcS(iFc)*geo%fcQalf(iFc,1)*
     &          (dv%vadia(iFc,1,is) + dv%vaecrb(iFc,1,is)))*
     &          0.5_R8*(pl%na(iCv1,is) + pl%na(iCv2,is))/pl%na(iCv1,is))
            else
              srw%sna0(iCv1,1,is) = srw%sna0(iCv1,1,is) -
     &          max(0.0e0_R8,mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &          (geo%fcS(iFc)*geo%fcQalf(iFc,0)*dv%vaecrb(iFc,0,is) +
     &           geo%fcS(iFc)*geo%fcQalf(iFc,1)*dv%vaecrb(iFc,1,is))*
     &          0.5_R8*(pl%na(iCv1,is) + pl%na(iCv2,is))/pl%na(iCv1,is))
              t0 = mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     .          (dv%fnaPSch(iFc,0,is) + dv%fnaPSch(iFc,1,is))
              if(t0.ge.0.0_R8) then
                srw%sna0(iCv1,0,is) = srw%sna0(iCv1,0,is) + t0
              else
                srw%sna0(iCv1,1,is) = srw%sna0(iCv1,1,is) +
     .            t0/pl%na(iCv1,is)
              endif
            endif
          enddo

* -- BCCON=11 -- particle flux feedback boundary condition

        case (11)
          call xerrab ('b2stbc_phys -- BCCON = 11 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
           write(*,'(a,a,a,a,i3)')
     .      'BCCON = 11 : particle flux feedback boundary condition '//
     .      'on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
           call ipgeti ('b2stbc_isfeedback', isfeedback)
           call xertst (0.le.isfeedback.and.isfeedback.lt.ns,
     1      'faulty internal parameter isfeedback')
           if (switch%use_eirene.ne.0) is_kinetic = 1
           fedback_species = is.eq.isfeedback
           do isa = 0, nspecies-1
             if (na_feedback_ib(isa).eq.ib) then
               fedback_species = fedback_species .or.
     .                          (is.ge.eb2spcr(isa+1)+is_kinetic .and.
     .                           is.le.eb2spcr(isa+1)+nfluids(isa+1))
             end if
           end do
           call xertst(fedback_species,
     1       'BCCON(,IS,1) = 11 for non-feedback species!')
           call xertst(lfeedback, 'LFEEDBACK must be turned on!')
          endif
          fedback_species = .false.
          do isa = 0, nspecies-1
            if (na_feedback_target(isa).gt.0.0_R8 .and.
     &          na_feedback_ib(isa).eq.ib .and.
     &          is.ge.eb2spcr(isa+1)+is_kinetic .and.
     .          is.le.eb2spcr(isa+1)+nfluids(isa+1)) then
              t0 = conpar(is,ib,1)
              fedback_species = .true.
            end if
          end do
! csc Warning! Here nothing prevents feedback being applied on E-W boundaries
!     Needs redefinition also in view of extended grids or hybrid
          if (.not.fedback_species) then
*   ..initialise ireg, not sure is asbolutely correct for WG code
            ireg = mpg%cvReg(mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1))
            t0 = cbsna(0,IS,ireg)
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            srw%sna0(iCv1,0,is) = t0*geo%fcS(iFc)
            srw%sna0(iCv1,1,is) = 0.0_R8
            sna0_no_mdf(iCv1,0,is) = srw%sna0(iCv1,0,is)
            sna0_no_mdf(iCv1,1,is) = srw%sna0(iCv1,1,is)
          enddo
          IF (BCCHAR(IB).EQ.'S') THEN
            if (.not.fedback_species) then
*   ..initialise ireg
              ireg = cbirso
*   ..update ireg
 1101         if (cbrbrk(ireg)*nx.lt.
     .          maxval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8) then
                call xertst (ireg.lt.cbirso+cbnrso-1,
     &           'error in computing region on south boundary')
                ireg = ireg+1
                goto 1101
              endif
              if (ireg.gt.cbirso) then
                call xertst (cbrbrk(ireg-1)*nx.lt.
     &           minval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8,
     &           'feedback bc covers more than one region!')
              endif
              t0 = cbsna(0,IS,ireg)
            end if
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         t0*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            ENDDO
          ELSEIF (BCCHAR(IB).EQ.'N') THEN
            if (.not.fedback_species) then
*   ..initialise ireg
              ireg = cbirno
*   ..update ireg
 1102         if (cbrbrk(ireg)*nx.lt.
     .          maxval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8) then
                call xertst (ireg.lt.cbirno+cbnrno-1,
     &           'error in computing region on north boundary')
                ireg = ireg+1
                goto 1102
              endif
              if (ireg.gt.cbirno) then
                call xertst (cbrbrk(ireg-1)*nx.lt.
     &           minval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8,
     &           'feedback bc covers more than one region!')
              endif
              t0 = cbsna(0,IS,ireg)
            end if
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         t0*SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .               BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            ENDDO
          ELSE
            CALL XERRAB
     .       ('B2STBC -- BCCON=11, NOT AVAILABLE FOR W,E')
          ENDIF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                       !srv 23.09.08 }
#endif
* -- BCCON=12 -- particle density feedback boundary condition

        case (12)
          call xerrab ('b2stbc_phys -- BCCON = 12 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
           write(*,'(a,a,a,a,a,i3)')
     .      'BCCON = 12 : particle density feedback ',
     1      'boundary condition on ',BCCHAR(IB),boundary_location(ib),
     2      ' for species ',IS
           call ipgeti ('b2stbc_isfeedback', isfeedback)
           call xertst (0.le.isfeedback.and.isfeedback.lt.ns,
     1      'faulty internal parameter isfeedback')
           call xertst(is.eq.isfeedback,
     1       'BCCON(,IS,1) = 12 for IS <> ISFEEDBACK!')
           call xertst(lfeedback, 'LFEEDBACK must be turned on!')
           call xertst ((0.0_R8.eq.switch%nesepm.and.0.0_R8.eq.
     &      switch%ndes_sol).or.0.0_R8.eq.switch%nesepm_overshoot,
     1      'BCCON=12 coding incompatible with nesepm_overshoot!')
           icv = mpg%bcCv(mpg%bcCvP(ib,1),1)
           if (mpg%cvOnClosedSurface(icv).and.ishigh)
     >      call xertst (switch%ionising_core.eq.0,
     .      'ionising_core switch usage inconsistent with BCCON = 12')
          endif
          IF (BCCHAR(IB).EQ.'S') THEN
*   ..initialise ireg
            ireg = cbirso
*   ..update ireg
 1201       if (cbrbrk(ireg)*nx.lt.
     .          maxval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8) then
              call xertst (ireg.lt.cbirso+cbnrso-1,
     &         'error in computing region on south boundary')
              ireg = ireg+1
              goto 1201
            endif
            if (ireg.gt.cbirso) then
              call xertst (cbrbrk(ireg-1)*nx.lt.
     .          minval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8,
     &         'feedback bc covers more than one region!')
            endif
            t0=cbsna(0,IS,ireg)
            call xertst (t0.gt.0.0_R8,
     &       'BCCON=12, trying to set a negative density!')
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*t0*
     .         SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            ENDDO
          ELSEIF (BCCHAR(IB).EQ.'N') THEN
*   ..initialise ireg
            ireg = cbirno
*   ..update ireg
 1202       if (cbrbrk(ireg)*nx.lt.
     .          maxval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8) then
              call xertst (ireg.lt.cbirno+cbnrno-1,
     &         'error in computing region on north boundary')
              ireg = ireg+1
              goto 1202
            endif
            if (ireg.gt.cbirno) then
              call xertst (cbrbrk(ireg-1)*nx.lt.
     .          minval(bc_list_x(1:bc_list_size(ib),ib))+0.5_R8,
     &         'feedback bc covers more than one region!')
            endif
            t0=cbsna(0,IS,ireg)
            call xertst (t0.gt.0.0_R8,
     &       'BCCON=12, trying to set a negative density!')
            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*t0*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=-VBIG*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            ENDDO
          ELSE
            CALL XERRAB
     .       ('B2STBC -- BCCON=12, NOT AVAILABLE FOR W,E')
          ENDIF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                       !srv 23.09.08 }
#endif
* -- BCCON=13 -- Feedback constant density for a given particle flux

        case (13)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,a,a,a,i3,a,1p,2g14.7)')
     .       'BCCON = 13 : particle density to achieve '//
     .       'specified total flux on ',
     1       BCCHAR(IB),boundary_location(ib),
     2       ' for species ',IS,' with parameters ',
     3       CONPAR(IS,IB,1),CONPAR(IS,IB,2)
            icv = mpg%bcCv(mpg%bcCvP(ib,1),1)
            if (switch%ionising_core.gt.0) write(*,*)
     .       'imposed flux adapted due to ionising_core switch usage'
          endif
*       ..compute average density on boundary
          t0 = 0.0e0_R8
          t1 = 0.0e0_R8
          totflux = 0.0e0_R8
          !wdk Note on backwards compatibility: feedback now always on
          !wdk flux through the boundary itself (no offset)
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            t0 = t0 + geo%fcS(iFc)
            t1 = t1 + geo%fcS(iFc)*pl%na(iCv1,is)
            if (bccon_flux_fb_style.eq.0) then
              totflux = totflux - mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &          (dv%fna(iFc,0,is) + dv%fna(iFc,1,is))
            else
              totflux = totflux - mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &          (dv%fna_eir(iFc,0,is) + dv%fna_eir(iFc,1,is))
            endif
          enddo
          t1 = t1/t0
cnh conpar(is,ib,3) contains the contribution from the ionized neutrals
c   (ionising_core option)
          totpar = conpar(is,ib,1) + conpar(is,ib,3)
          if (conpar(is,ib,2).gt.0.0_R8) t3 = conpar(is,ib,2)
          if (conpar(is,ib,2).le.0.0_R8) t3 = 1.0e-2_R8
          t2 = (1.0_R8-t3*(totflux-totpar)/
     &         max(abs(totflux),abs(totpar),1.0e10_R8))*t1
          if (switch%b2stbc_diagno.ge.2) then
           write(*,'(a,i3,1p,4(1x,e15.6))')
     &      'bccon=13: is, curr_flux, desired_flux, old/new density, ',
     &      is,totflux,totpar,t1,t2
          endif
*       ..impose density
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%sna0(iCv1,0,is) =  VBIG*geo%fcS(iFc)*t2
            srw%sna0(iCv1,1,is) = -VBIG*geo%fcS(iFc)
            sna0_no_mdf(iCv1,0,is) = srw%sna0(iCv1,0,is)
            sna0_no_mdf(iCv1,1,is) = srw%sna0(iCv1,1,is)
          enddo

* -- BCCON=14 -- Sound speed flux

        case (14)
          if (ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCCON = 14 : sound speed flux ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1       'BCCON = 14, CONPAR(,,1) <= 0 not allowed!')
            call xertst(CONPAR(IS,IB,2).ge.0.0_R8,
     1       'BCCON = 14, CONPAR(,,2) < 0 not allowed!')
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                    ! number of the corresponding domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            s1 = abs(geo%fcPbs(iFc))
            cs = sqrt(pz(iCv1)/rz(iCv1))
            srw%sna0(iCv1,0,is) = 0.0_R8
            srw%sna0(iCv1,1,is) = -s1*conpar(is,ib,1)*cs
            if (conpar(is,ib,2).gt.0.0_R8) then
              srw%sna0(iCv1,1,is) = srw%sna0(iCv1,1,is) -
     .          co%dna0(iCv1,is)/conpar(is,ib,2)*geo%fcS(iFc)
            endif
            sna0_no_mdf(iCv1,:,is) = srw%sna0(iCv1,:,is)
            if (switch%mdf_fnb.ne.0) then
              t0 = mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &          ((dv%fnaPSch(iFc,0,is) + dv%fnaPSch(iFc,1,is)) +
     &           (geo%fcS(iFc)*geo%fcQalf(iFc,0)*dv%vadia(iFc,0,is) +
     &            geo%fcS(iFc)*geo%fcQalf(iFc,1)*dv%vadia(iFc,1,is))*
     &           0.5_R8*(pl%na(iCv1,is)+pl%na(iCv2,is)))
              srw%sna0(iCv1,0,is) = srw%sna0(iCv1,0,is) +
     &          max(0.0_R8,t0)
              srw%sna0(iCv1,1,is) = srw%sna0(iCv1,1,is) +
     &          min(0.0_R8,t0)/pl%na(iCv1,is)
            endif
          enddo

* -- BCCON=15 -- prescribe a radial leakage velocity, as a multiplier to the local thermal velocity

        case (15)                                                   !srv 25.09.09 16.02.12
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(2a,1p,g14.7,a,a,a,a,i3)')
     .       'BCCON = 15 : prescribe the decay length for the density',
     &       ' for drifs switched on ',CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            icv = mpg%bcCv(mpg%bcCvP(ib,1),1)
            if (mpg%cvOnClosedSurface(icv).and.ishigh)
     .       call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 15')
          endif
cwdk      Notes:
cwdk       * need correction for misaligned grids -> t1
cwdk       * need to correct this BC for hybrid scheme
          do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            csb  = sqrt((rt%rza(iCv1,is)*pl%te(iCv1) + pl%ti(iCv1))/
     &                  (am(is)*mp))
            t0   = geo%fcQalf(iFc,1)
            t1   = (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     &              geo%fcQbet(iFc,0)
            srw%sna0(iCv1,0,is) =  VBIG*
     &          co%cdna(iFc,1,is)*(t0*pl%na(iCv2,is) -
     &          mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*t1*
     &         (pl%na(mpg%fcVx(iFc,2),is) - pl%na(mpg%fcVx(iFc,1),is)))
            srw%sna0(iCv1,1,is) = -VBIG*
     &         (co%cdna(iFc,1,is)*t0 - CONPAR(is,iB,1)*geo%fcS(iFc)*csb)
            sna0_no_mdf(iCv1,:,is) = srw%sna0(iCv1,:,is)
          enddo                                                        !srv 25.09.09 }

* -- BCCON=16 -- Feedback boundary condition: given total particle flux with constant density  !srv 16.02.12 {

        case (16)
          call xerrab ('b2stbc_phys -- BCCON = 16 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
#ifdef ASTRA
            ! IYS 08.04.2014 temporary solution while we do not have impurities in ASTRA
            call xertst(is.eq.ismain,
     1         'BCCON = 16 is applied for main ion only for ASTRA') !srv 29.09.15
#endif
            write(*,'(a,a,i3)')
     .       'BCCON = 16 : given total particle flux with constant',
     &       ' density for ion is =', is
!           write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     1           boundary_location(ibw),' for species ',IS
                call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                        BC_LIST_Y(1,IBW)),
     .           'BCCON = 16 only applicable to core boundaries !')
              endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.19.or.
     1            BCCON(IS,IBW).EQ.21.or.BCCON(IS,IBW).EQ.22.or.
     1            BCCON(IS,IBW).EQ.23.or.BCCON(IS,IBW).EQ.25.or.
     1            BCCON(IS,IBW).EQ.26.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab
     .('BCCON types 16, 19, 21, 22, 23, 25, 26 and 27 cannot be mixed!')
              endif
            ENDDO
            icv = mpg%bcCv(mpg%bcCvP(ib,1),1)
            if (mpg%cvOnClosedSurface(icv).and.ishigh)
     .       call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 16')
          endif
          if(.not.once_bccon(is)) then
            IF(BCCHAR(IB).EQ.'S') THEN
              if (switch%use_astra.ne.0) then
#ifdef ASTRA
                SCONPAR(IS) = fneb_astra
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SCONPAR(IS) = 0.0_R8
                DO IBW=IB,NBC
                  if(BCCHAR(IBW).EQ.'S'.and.
     >               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    SCONPAR(IS) = SCONPAR(IS) + CONPAR(IS,IBW,1)
                  endif
                ENDDO
              endif
              write(*,*) 'BCCON=16: Summed value is ', SCONPAR(IS)
              NAS = 0.0_R8
              VS = 0.0_R8
              DO IBW=IB,NBC
                if(BCCHAR(IBW).EQ.'S'.and.
     &             BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    NAS = NAS+NA(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)
     .                      *VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    VS = VS +VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  END DO
                endif
              END DO
              NAS=NAS/VS
              totflux = 0.0_R8
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  if(switch%mdf_fnb.ne.0) then
                    DO IX = 1, BC_LIST_SIZE(IBW)
                      totflux = totflux +
     &        fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,is)
                    ENDDO
                  else
                    DO IX = 1, BC_LIST_SIZE(IBW)
                      totflux = totflux +
     &            fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,is)
                    ENDDO
                  endif
                ENDIF
              ENDDO
              IF (max(abs(totflux),abs(SCONPAR(IS))).GT.0.0_R8)
     >         NAS = NAS *
     &         (1.0_R8-switch%b2stbc_bc_ref*(totflux-SCONPAR(IS))/
     &                 max(abs(totflux),abs(SCONPAR(IS))))
!              T0 = (totflux-CONPAR(IS,IB,1))/
!     &           max(abs(totflux),CONPAR(IS,IB,1))
!              NAS = NAS *
!     &         (1.0_R8-max(-switch%b2stbc_bc_ref,min(switch%b2stbc_bc_ref,T0)))
              write(*,'(a,es14.7,a,es14.7,a,es14.7,a,i3)')
     &         ' NAS=',NAS,' totflux=',totflux,
     &         ' SCONPAR=',SCONPAR(IS),' IS=',IS
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     *                VBIG*NAS*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     *               -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  END DO
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)=
     .               SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  END DO                                              !srv 23.09.08 }
                END IF
              END DO
              once_bccon(is) = .true.
            ELSE                                                      !srv 16.02.12 }
              CALL XERRAB
     1         ('B2STBC -- BCCON=16, NOT AVAILABLE FOR N,W,E')
            ENDIF
          endif                                                       !srv 16.02.12 }
#endif
* -- BCCON=17 -- FILLER

        case (17)
          CALL XERRAB ('B2STBC -- BCCON=17, NOT YET AVAILABLE')

#ifdef ASTRA
* -- BCCON=18 -- PRESCRIBE THE MODIFIED TOTAL PARTICLE FLUX

        case (18)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,a,a,a,i3)')
     .       'BCCON = 18 : modified total particle flux on ',
     1        BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(is.eq.ismain,
     1       'BCCON = 18 is applied for main ion only ')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            t0 = 0.0_R8
            US=0.0_R8
            isb = ismain0
            DO IX = 1, BC_LIST_SIZE(IB)
              US=US+SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              if (is_neutral(isb)) t0 = t0 +
     .             fna(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                 topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1,isb)
            END DO
            CONPAR(IS,IB,1) = (fSeNBI_astra - t0)                     !srv 29.11.11

!            if(ncall_b2stbc_phys.eq.0) then
               write (*,'(a,a,i3,1p,3(a,g11.4))')
     &        ' BCCON=18', ' IS =', IS,
     &        ' CONPAR=',CONPAR(IS,IB,1), ' fSeNBI_astra=',fSeNBI_astra,
     &        ' integrated neutral flux =',t0
!            endif

            DO IX = 1, BC_LIST_SIZE(IB)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     &         SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)+          !srv 01.02.09
     1          CONPAR(IS,IB,1)/US*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &         sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+
     &         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)           !srv 23.09.08
              sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &         sna0_no_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)+
     &         sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)           !srv 23.09.08
              if(switch%mdf_fnb.ne.0) then                                   !srv 12.09.08 {
               sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &          sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) -
     &          fnaPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1,is)
               sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) =
     &          sna0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is) +
     &           vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &           vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1,is)
              endif                                                   !srv 12.09.08 }
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCCON=18, ONLY FOR SOUTH')
          ENDIF
#else
* -- BCCON=18 -- FILLER

        case (18)
          CALL XERRAB ('B2STBC -- BCCON=18, NOT YET AVAILABLE')
#endif

* -- BCCON=19 -- Feedback boundary condition: given sum of integrated neutrals and
*                main ion particle fluxes with constant ion density.

        case (19)
          call xerrab ('b2stbc_phys -- BCCON = 19 not adapted for WG')
          if (switch%spatial_hybrid.ne.0) then
            call xerrab ('BCCON = 19 not compatible with spatial hybr')
          endif
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,*) 'BCCON = 19 : given sum of particle flux '//
     .                 'for isonuclear sequence with constant '//
     .                 'highest ionisation stage density '
!             write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     1           boundary_location(ibw),' for species ',IS
                call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                        BC_LIST_Y(1,IBW)),
     .            'BCCON = 19 only applicable to core boundaries !')
              endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.21.or.BCCON(IS,IBW).EQ.22.or.
     1            BCCON(IS,IBW).EQ.23.or.BCCON(IS,IBW).EQ.25.or.
     1            BCCON(IS,IBW).EQ.26.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab('BCCON types 16, 19, 21, 22, 23, 25, '//
     .                      '26 and 27 cannot be mixed!')
              endif
            ENDDO
            if (is.lt.ns-1) call xertst(.not.lnext(is,is+1),
     1       'BCCON = 19 can only be applied to ions in their '//
     2       'highest ionisation stage')
            if (ishigh) call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 19')
          endif
          if(.not.once_bccon(is)) then
            IF(BCCHAR(IB).EQ.'S') THEN
              if (switch%use_astra.ne.0) then
#ifdef ASTRA
                SCONPAR(IS) = fSeNBI_astra
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SCONPAR(IS) = 0.0_R8
                DO IBW=IB,NBC
                  if(BCCHAR(IBW).EQ.'S'.and.
     >               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    SCONPAR(IS) = SCONPAR(IS) + CONPAR(IS,IBW,1)
                  endif
                ENDDO
              endif
              write(*,*) 'BCCON=19: Summed value is ', SCONPAR(IS)
              NAS = 0.0_R8
              VS = 0.0_R8
              DO IBW=IB,NBC
                if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    NAS = NAS+NA(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)
     .                      *VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    VS = VS +VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  END DO
                endif
              ENDDO
              NAS = NAS/VS
              totfluxn = 0.0_R8
              totfluxi = 0.0_R8
#ifdef B25_EIRENE
cxpb  Add contribution from Eirene core boundary fluxes
              do isa = 1, nnatmi
                if (latmscl(isa).eq.b2espcr(is)) then
                  do k = 1, ncns
                    ix = nlim + lcns(k)
                    totfluxn = totfluxn -
     &               wldna(ix,isa,0) + wldra(ix,isa,0) + wldpa(ix,isa,0)
                    do ism = 1, nnmoli
                      if (mlcmp(isa,ism).gt.0) then
                        totfluxn = totfluxn - mlcmp(isa,ism) *
     &                   (wldnm(ix,ism,0) - wldrm(ix,ism,0)
     &                                    - wldpm(ix,ism,0))
                      end if
                    end do
                  end do
                end if
              end do
#endif
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
! totfluxn here means the sum of fluxes over all ions (including fluid
! neutrals) of this isotope except for the highest ionisation stage.
! For coupled runs, both fna and fna_mdf are rescaled for neutrals
! Note that for neutrals fna = fna_mdf
                    do isb = 0, ns-1
                      if(b2espcr(is).eq.b2espcr(isb).and.is.ne.isb) then
                        if(.not.is_neutral(isb).or.
     &                                      switch%use_eirene.eq.0) then
                          if(switch%mdf_fnb.ne.0) then
                            totfluxn = totfluxn +
     &       fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,isb)
                          else
                            totfluxn = totfluxn +
     &           fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,isb)
                          endif
                        endif
                      endif
                    enddo
! totfluxi here means the flux of the ion in the highest ionisation stage.
                    if(switch%mdf_fnb.ne.0) then
                      totfluxi = totfluxi +
     &       fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,is)
                    else
                      totfluxi = totfluxi +
     &           fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,is)
                    endif
                  ENDDO
                ENDIF
              ENDDO
              IF(max(abs(totfluxn+totfluxi),abs(SCONPAR(IS))).GT.0.0_R8)
     >         NAS = NAS *
     &         (1.0_R8-switch%b2stbc_bc_ref*
     &           (totfluxn+totfluxi-SCONPAR(IS))/
     &           max(abs(totfluxn+totfluxi),SCONPAR(IS)))
              write(*,'(a,es14.7,a,es14.7,a,es14.7,a,es14.7,a,es14.7)')
     *         ' NAS=',NAS,' totfluxn=',totfluxn,
     &         ' totfluxi=',totfluxi,' SCONPAR=',SCONPAR(IS),
     &         ' totfluxn+totfluxi=',totfluxn+totfluxi
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*NAS*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  ENDDO
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)=
     &               SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  ENDDO                                              !srv 23.09.08 }
                ENDIF
              ENDDO
              once_bccon(is) = .true.
            ELSE                                                     !srv 16.02.12 }
              CALL XERRAB
     1         ('B2STBC -- BCCON=19, NOT AVAILABLE FOR N,W,E')
            ENDIF
          endif                                                      !srv 16.02.12 }
#endif
* -- BCCON=20 -- constant density feedback condition

        case (20)
          call xerrab ('b2stbc_phys -- BCCON = 20 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCCON = 20 : constant density (neo) ',
     1       CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1       'BCCON = 1, CONPAR(,,1) !> 0')
            icv = mpg%bcCv(mpg%bcCvP(ib,1),1)
            if (mpg%cvOnClosedSurface(icv).and.ishigh)
     >       call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 20')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            naave = 0.0_R8
            volave = 0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type20_ix=bc_list_x(ix,ib)
              bc_type20_iy=bc_list_y(ix,ib)
              do bc_type20_i=1,bc_type20_ref  ! how many steps away from the boundary?
                bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
                bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
                bc_type20_ix=bc_type20_ixx
                call xertst (-2.lt.bc_type20_ix.and.
     &                             bc_type20_ix.lt.nx+1,
     1           'stepped off x-grid in type 20 bc')
                call xertst (-2.lt.bc_type20_iy.and.
     &                             bc_type20_iy.lt.ny+1,
     1           'stepped off y-grid in type 20 bc')
              enddo
              naave=naave+na(bc_type20_ix,bc_type20_iy,is)*
     1         vol(bc_type20_ix,bc_type20_iy)
              volave=volave+vol(bc_type20_ix,bc_type20_iy)
            end do
            naave=naave/volave
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type20_ix=bc_list_x(ix,ib)
              bc_type20_iy=bc_list_y(ix,ib)
              do bc_type20_i=1,bc_type20_ref  ! how many steps away from the boundary?
                bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
                bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
                bc_type20_ix=bc_type20_ixx
                call xertst (-2.lt.bc_type20_ix.and.
     &                             bc_type20_ix.lt.nx+1,
     1           'stepped off x-grid in type 20 bc')
                call xertst (-2.lt.bc_type20_iy.and.
     &                             bc_type20_iy.lt.ny+1,
     1           'stepped off y-grid in type 20 bc')
              enddo
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*CONPAR(IS,IB,1)*
     .         SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         na(bc_type20_ix,bc_type20_iy,is)/naave
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=20, NOT AVAILABLE FOR N,W,E')
          ENDIF
          DO IX = 1, BC_LIST_SIZE(IB)
            SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)=
     .       SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),:,IS)
          ENDDO                                                     !srv 23.09.08 }
#endif
cxpb Careful here because of BC number clash
* From 5.0: starting at version 01.001.030
* -- BCCON=21 -- constant density scaled by density on the ring bc_type21_ref away
*                has been moved up to next available number (BCCON=24)

*  From 5.2, after version number 01.001.078:
*            -- BCCON=21 -- PRESCRIBE THE VALUE OF THE DENSITY !lk and IYS 01.03.2013
*                           AND ADD A POLOIDAL DENSITY VARIATION TO GET A SOLUTION
*                           WHICH IS AS CLOSE AS POSSIBLE TO NEOCLASSICAL THEORY
*                           DERIVED FOR PURE H/D/T PLASMA ONLY.
*
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH
* CORRESPONDING CONDITION ON ION TEMPERATURE (BCENI=23,24)

        case (21)
          call xerrab ('b2stbc_phys -- BCCON = 21 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            if (boundary_version.ge.'01.001.030' .and.
     &          boundary_version.lt.'02.000.000') then
              write(*,*) 'BCCON=21 type B.C. moved to 24 ',
     &         ' because of clash with 5.2 numbering. '
              write(*,*) 'Please modify b2.boundary.parameters file'
              call xerrab('B2STBC -- BCCON=21, DISPLACED')
            endif
            write(*,'(a,1p,g14.7)')
     .       'BCCON = 21 : constant density',CONPAR(IS,IB,1)
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                        BC_LIST_Y(1,IBW)),
     .           'BCCON = 21 only applicable to core boundaries !')
                call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1           'BCCON = 21, CONPAR(,,1) !> 0')
                call xertst(CONPAR(IS,IB,1).eq.CONPAR(IS,IBW,1),  ! IYS 22.01.2015
     1           'BCCON = 21, Average density on all parts ' //
     1             'of the surface should be the same')
                if (bceni(ibw).ne.23.and.bceni(ibw).ne.24) then
                  write(0,*) 'Should use BCCON(IS,IB).EQ.21 '//
     &                       'together with BCENI=23 or 24'
                  write(*,*) 'Should use BCCON(IS,IB).EQ.21 '//
     &                       'together with BCENI=23 or 24'
                  write(*,'(a)')
     &             'Please ensure that this boundary condition '//
     &             'is consistent with the '//
     &             'ion energy boundary condition'                  !srv 21.09.12
                endif
              endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.22.or.
     1            BCCON(IS,IBW).EQ.23.or.BCCON(IS,IBW).EQ.25.or.
     1            BCCON(IS,IBW).EQ.26.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab('BCCON types 16, 19, 21, 22, 23, 25, '//
     .                      '26 and 27 cannot be mixed!')
              endif
            END DO
            call xertst(is.eq.ismain,
     1       'BCCON = 21 is applied for main ion only ')
            if (.not.(ns.eq.1.or.(ns.eq.2.and.is_neutral(0)))) then
              write(0,*)
     &         'Use of BCCON=21 only meant for single-fluid runs!'
              write(0,*) 'Use BCCON=25 instead!'
              write(6,*)
     &         'Use of BCCON=21 only meant for single-fluid runs!'
              write(6,*) 'Use BCCON=25 instead!'
            endif
            if (ishigh) call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 21')
          endif
          IF(BCCHAR(IB).EQ.'S' .and.
     &       on_closed_surface(BC_LIST_X(1,IB),BC_LIST_Y(1,IB))) THEN
            if (.not.once_bccon(is)) then
              if(is.eq.ismain.and.switch%use_astra.ne.0) then
#ifdef ASTRA
                DO IBW=IB,NBC
                  if(BCCHAR(IBW).EQ.'S'.and.
     &               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    CONPAR(IS,IBW,1) = neb_astra                       !srv 28.07.11
                  endif
                ENDDO
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              endif
              t0=0.0e0_R8   ! Volume between two surfaces               !lk 12.12.12 {
              t1=0.0e0_R8   ! Sum of volume-averaged electron and ion temperatures
              DO IBW=IB,NBC
                IF (BCCHAR(IBW).EQ.'S' .and.
     &              BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    t0=t0 +
     .                vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    t1=t1 +
     .                (ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &                 te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &                vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  ENDDO
                ENDIF
              ENDDO
              if (t0.gt.0.0_R8) t1=t1/t0        ! Sum of volume-averaged electron and ion temperatures
              DO IBW=IB,NBC
                IF (BCCHAR(IBW).EQ.'S' .and.
     &              BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                  do IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &                 VBIG*CONPAR(IS,IBW,1)*t1
     &                 /(ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &                   te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &                   SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &                -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  enddo                                                !lk 12.12.12 }
                  DO IX = 1, BC_LIST_SIZE(IBW)                         !srv 23.09.08
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)=
     &               SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  ENDDO
                ENDIF
              ENDDO
              once_bccon(is) = .true.
            endif
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=21, VALID FOR CORE ONLY')
          ENDIF      !lk and IYS 01.03.2013}
#endif
* -- BCCON=22 -- Feedback boundary condition: given total particle flux with
*    constant average density. The poloidal density variation is added to make
*    the solution as close as possible to the neoclassical value derived for
*    pure H/D/T plasma.
*    It is recommended to use this boundary condition
*    together with corresponding boundary condition on ion temperature
*    (BCENI=23,24) !srv 16.02.12, IYS 01.03.2013 {

        case (22)
          call xerrab ('b2stbc_phys -- BCCON = 22 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,a,i3)')
     &       'BCCON = 22 : given total particle flux with constant',
     &       ' density for ion is =', is
!           write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                        BC_LIST_Y(1,IBW)),
     .           'BCCON = 22 only applicable to core boundaries !')
                if (bceni(ibw).ne.23.and.bceni(ibw).ne.24) then
                  write(0,*) 'Should use BCCON(IS,IB).EQ.22 '//
     &                       'together with BCENI=23 or 24'
                  write(*,*) 'Should use BCCON(IS,IB).EQ.22 '//
     &                       'together with BCENI=23 or 24'
                  write(*,'(a)')
     &             'Please ensure that this boundary condition '//
     &             'is consistent with the '//
     &             'ion energy boundary condition'                    !srv 21.09.12
                endif
              endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.21.or.
     1            BCCON(IS,IBW).EQ.23.or.BCCON(IS,IBW).EQ.25.or.
     1            BCCON(IS,IBW).EQ.26.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab('BCCON types 16, 19, 21, 22, 23, 25, '//
     .                      '26 and 27 cannot be mixed!')
              endif
            ENDDO
            call xertst(is.eq.ismain,
     1       'BCCON = 22 is applied for main ion only ')
            if (.not.(ns.eq.1.or.(ns.eq.2.and.is_neutral(0)))) then
              write(0,*)
     &         'Use of BCCON=22 only meant for single-fluid runs!'
              write(0,*) 'Use BCCON=26 instead!'
              write(6,*)
     &         'Use of BCCON=22 only meant for single-fluid runs!'
              write(6,*) 'Use BCCON=26 instead!'
            endif
            if (ishigh) call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 22')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            if(.not.once_bccon(is)) then
              if (use_astra.ne.0) then
#ifdef ASTRA
                SCONPAR(IS) = fneb_astra
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SCONPAR(IS) = 0.0_R8
                DO IBW=IB,NBC
                  if(BCCHAR(IBW).EQ.'S'.and.
     >               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    SCONPAR(IS) = SCONPAR(IS) + CONPAR(IS,IBW,1)
                  endif
                ENDDO
              endif
              write(*,*) 'BCCON=22: Summed value is ', SCONPAR(IS)
              totflux = 0.0_R8
              t0=0.0e0_R8   ! Volume between two surfaces             !lk 12.12.12 {
              t1=0.0e0_R8   ! Sum of volume-averaged electron and ion temperatures
              NAS=0.0e0_R8
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    t0=t0 +
     &                vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    t1=t1 +
     &                (ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &                 te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     .                vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    NAS=NAS+na(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    if(mdf_fnb.ne.0) then
                      totflux = totflux +
     &         fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,is)
                    else
                      totflux = totflux +
     &         fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &             topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,is)
                    endif
                  END DO
                ENDIF
              ENDDO
              t1=t1/t0        ! Sum of volume-averaged electron and ion temperatures
              NAS=NAS/t0      ! Volume-averaged density
              IF (max(abs(totflux),abs(SCONPAR(IS))).GT.0.0_R8)
     >         NAS = NAS*
     &         (1.0_R8-switch%b2stbc_bc_ref*(totflux-SCONPAR(IS))/
     &           max(abs(totflux),abs(SCONPAR(IS))))
!              T0 = (totflux-CONPAR(IS,IB,1))/
!     &           max(abs(totflux),CONPAR(IS,IB,1))
!              NAS = NA((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB),is) *
!     &         (1.0_R8-max(-switch%b2stbc_bc_ref,min(switch%b2stbc_bc_ref,T0)))
              write(*,'(a,es14.7,a,es14.7,a,es14.7,a,i3)')
     .         ' NAS=',NAS,' totflux=',totflux,
     &         ' SCONPAR=',SCONPAR(IS),' IS=',IS
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*NAS*t1
     &                   /(ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &                     te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
     &                    *SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))   ! lk
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  ENDDO
                  DO IX = 1, BC_LIST_SIZE(IBW)                       !srv 23.09.08 {
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)
     &                =SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  ENDDO                                              !srv 23.09.08 }
                ENDIF
              ENDDO
              once_bccon(is) = .true.
            endif                                                    !srv 16.02.12 }
          ELSE                                                       !srv 16.02.12 }
            CALL XERRAB
     1       ('B2STBC -- BCCON=22 IS AVAILABLE FOR CORE ONLY')
          ENDIF                                                      !iys 01.03.13 }
#endif
* -- BCCON=23 -- Feedback boundary condition: given sum of integrated neutrals and
*    main ion particle fluxes with constant average density.
*    The poloidal density variation is added to make the solution
*    as close as possible to the neoclassical value derived for pure H/D/T plasma.
*    It is recommended to use this boundary condition with corresponding
*    boundary condition on ion temperature (BCENI=23,24) !srv 16.02.12, IYS 04.03.2013 {
*    Use this boundary condition for pure H/D/T plasma only,
*    If a minor impurity is present, please use the BCCON=27 instead. ! IYS 21.01.2016

        case (23)
          call xerrab ('b2stbc_phys -- BCCON = 23 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,a,i3)')
     .       'BCCON = 23 : given total particle flux with constant',
     &       ' density for main ion is =', ismain
!          write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                        BC_LIST_Y(1,IBW)),
     .           'BCCON = 23 only applicable to core boundaries !')
                call xertst(CONPAR(IS,IBW,1).ge.0.0_R8,
     1           'BCCON = 23, CONPAR(,,1) < 0')
                if (bceni(ibw).ne.23.and.bceni(ibw).ne.24) then
                  write(0,*) 'Should use BCCON(IS,IB).EQ.23 '//
     &                       'together with BCENI=23 or 24'
                  write(*,*) 'Should use BCCON(IS,IB).EQ.23 '//
     &                       'together with BCENI=23 or 24'
                  write(*,'(a)')
     &             'Please ensure that this boundary condition '//
     &             'is consistent with the '//
     &             'ion energy boundary condition'                  !srv 21.09.12
                endif
              endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.21.or.
     1            BCCON(IS,IBW).EQ.22.or.BCCON(IS,IBW).EQ.25.or.
     1            BCCON(IS,IBW).EQ.26.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab('BCCON types 16, 19, 21, 22, 23, 25, '//
     .                      '26 and 27 cannot be mixed!')
              endif
            ENDDO
            if (.not.(ns.eq.1.or.(ns.eq.2.and.is_neutral(0)))) then
              write(0,*)
     &         'Use of BCCON=23 only meant for single-fluid runs!'
              write(0,*) 'Use BCCON=27 instead!'
              write(6,*)
     &         'Use of BCCON=23 only meant for single-fluid runs!'
              write(0,*) 'Use BCCON=27 instead!'
            endif
            call xertst(is.eq.ismain,
     1       'BCCON = 23 is applied for main ion only ')
            if (ishigh) call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 23')
          endif

          IF(BCCHAR(IB).EQ.'S') THEN
            if(.not.once_bccon(is)) then
              if (use_astra.ne.0) then
#ifdef ASTRA
                SCONPAR(IS) = fSeNBI_astra
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SCONPAR(IS) = 0.0_R8
                DO IBW=IB,NBC
                  if(BCCHAR(IBW).EQ.'S'.and.
     &               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    SCONPAR(IS) = SCONPAR(IS) + CONPAR(IS,IBW,1)
                  endif
                ENDDO
              endif
              write(*,*) 'BCCON=23: Summed value is ', SCONPAR(IS)
              totfluxn = 0.0_R8
              totfluxi = 0.0_R8
              t0=0.0e0_R8   ! Volume between two surfaces               !lk 12.12.12 {
              t1=0.0e0_R8   ! Sum of volume-averaged electron and ion temperatures
              NAS=0.0e0_R8
#ifdef B25_EIRENE
cxpb  Add contribution from Eirene core boundary fluxes
              do isa = 1, nnatmi
                if (latmscl(isa).eq.b2espcr(ismain0)) then
                  do k = 1, ncns
                    ix = nlim + lcns(k)
                    totfluxn = totfluxn -
     &               wldna(ix,isa,0) + wldra(ix,isa,0) + wldpa(ix,isa,0)
                    do ism = 1, nnmoli
                      if (mlcmp(isa,ism).gt.0) then
                        totfluxn = totfluxn - mlcmp(isa,ism) *
     &                   (wldnm(ix,ism,0) - wldrm(ix,ism,0)
     &                                    - wldpm(ix,ism,0))
                      end if
                    end do
                  end do
                end if
              end do
#endif
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    isb = ismain0
                    if (is_neutral(isb)) totfluxn = totfluxn +
     &           fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,isb)
                    isb = ismain
                    if(mdf_fnb.ne.0) then
                      totfluxi = totfluxi +
     &       fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,isb)
                    else
                      totfluxi = totfluxi +
     &           fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,isb)
                    endif
                    t0=t0+vol(topix(BC_LIST_X(IX,IBW),
     .                              BC_LIST_Y(IX,IBW)),
     &                        topiy(BC_LIST_X(IX,IBW),
     .                              BC_LIST_Y(IX,IBW)))
                    t1=t1 + (ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     .                       te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     .                vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    NAS=NAS+na(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  ENDDO
                ENDIF
              ENDDO
              t1=t1/t0        ! Sum of volume-averaged electron and ion temperatures
              NAS=NAS/t0      ! Volume-averaged density
              IF(abs(totfluxn)+abs(totfluxi)+abs(SCONPAR(IS)).GT.0.0_R8)
     >         NAS = NAS*
     &         (1.0_R8-switch%b2stbc_bc_ref*
     &           (totfluxn+totfluxi-SCONPAR(IS))/
     &           (abs(totfluxn)+abs(totfluxi)+abs(SCONPAR(IS))))
              write(*,'(5(a,es14.7))')
     *          ' NAS=',NAS,' totfluxn=',totfluxn,
     &          ' totfluxi=',totfluxi,' SCONPAR=', SCONPAR(IS),
     &          ' totfluxn+totfluxi=',totfluxn+totfluxi
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*NAS*t1
     &               /(ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &                 te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &                 SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))         ! lk
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  ENDDO
                  DO IX =1, BC_LIST_SIZE(IBW)                          !srv 23.09.08 {
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)=
     &               SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  ENDDO                                                !srv 23.09.08 }
                ENDIF
              ENDDO
              once_bccon(is) = .true.
            endif                                                      !srv 16.02.12 }
          ELSE                                                         !srv 16.02.12 }
            CALL XERRAB
     1        ('B2STBC -- BCCON=23 IS AVAILABLE FOR CORE ONLY')
          ENDIF                                                        !IYS 04.03.2013 }
#endif
* -- BCCON=24 -- constant density scaled by density on the ring bc_type21_ref away

        case (24)
          call xerrab ('b2stbc_phys -- BCCON = 24 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCCON = 24 : constant density (modified by nearby ring) ',
     1       CONPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1       'BCCON = 24, CONPAR(,,1) !> 0')
            call xertst(CONPAR(is,ib,2).gt.0.0_R8,
     &           'faulty parameter conpar2')
            icv = mpg%bcCv(mpg%bcCvP(ib,1),1)
            if (mpg%cvOnClosedSurface(icv).and.ishigh)
     >       call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 24')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            naave = 0.0_R8
            volave = 0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type21_ix=bc_list_x(ix,ib)
              bc_type21_iy=bc_list_y(ix,ib)
              do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
                bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
                bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
                bc_type21_ix=bc_type21_ixx
                call xertst (-2.lt.bc_type21_ix.and.
     &                             bc_type21_ix.lt.nx+1,
     1           'stepped off x-grid in type 24 bc')
                call xertst (-2.lt.bc_type21_iy.and.
     &                             bc_type21_iy.lt.ny+1,
     1           'stepped off y-grid in type 24 bc')
              enddo
              naave=naave+na(bc_type21_ix,bc_type21_iy,is)*
     1                   vol(bc_type21_ix,bc_type21_iy)
              volave=volave+vol(bc_type21_ix,bc_type21_iy)
            end do
            naave=naave/volave
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type21_ix=bc_list_x(ix,ib)
              bc_type21_iy=bc_list_y(ix,ib)
              do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
                bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
                bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
                bc_type21_ix=bc_type21_ixx
                call xertst (-2.lt.bc_type21_ix.and.
     &                             bc_type21_ix.lt.nx+1,
     1           'stepped off x-grid in type 24 bc')
                call xertst (-2.lt.bc_type21_iy.and.
     &                             bc_type21_iy.lt.ny+1,
     1           'stepped off y-grid in type 24 bc')
              enddo
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*CONPAR(IS,IB,1)*
     .         SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         ((1.0_R8-conpar(is,ib,2))+
     .           conpar(is,ib,2)*na(bc_type21_ix,bc_type21_iy,is)/naave)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            end do
          ELSE IF(BCCHAR(IB).EQ.'N') THEN
            naave = 0.0_R8
            volave = 0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type21_ix=bc_list_x(ix,ib)
              bc_type21_iy=bc_list_y(ix,ib)
              do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
                bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
                bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
                bc_type21_ix=bc_type21_ixx
                call xertst (-2.lt.bc_type21_ix.and.
     &                             bc_type21_ix.lt.nx+1,
     1           'stepped off x-grid in type 24 bc')
                call xertst (-2.lt.bc_type21_iy.and.
     &                             bc_type21_iy.lt.ny+1,
     1           'stepped off y-grid in type 24 bc')
              enddo
              naave=naave+na(bc_type21_ix,bc_type21_iy,is)*
     1                   vol(bc_type21_ix,bc_type21_iy)
              volave=volave+vol(bc_type21_ix,bc_type21_iy)
            end do
            naave=naave/volave
            DO IX = 1, BC_LIST_SIZE(IB)
              bc_type21_ix=bc_list_x(ix,ib)
              bc_type21_iy=bc_list_y(ix,ib)
              do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
                bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
                bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
                bc_type21_ix=bc_type21_ixx
                call xertst (-2.lt.bc_type21_ix.and.
     &                             bc_type21_ix.lt.nx+1,
     1           'stepped off x-grid in type 24 bc')
                call xertst (-2.lt.bc_type21_iy.and.
     &                             bc_type21_iy.lt.ny+1,
     1           'stepped off y-grid in type 24 bc')
              enddo
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*CONPAR(IS,IB,1)*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .         ((1.0_R8-conpar(is,ib,2))+
     .           conpar(is,ib,2)*na(bc_type21_ix,bc_type21_iy,is)/naave)
              SNA0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=-VBIG*
     .         SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=24, NOT AVAILABLE FOR W,E')
          ENDIF
#endif
* -- BCCON=25  Feedback boundary condition: PRESCRIBE THE AVERAGE VALUE OF THE DENSITY  !lk and IYS 01.03.2013{
*              AND ADD A POLOIDAL DENSITY VARIATION FROM NEIGHBOURING RADIAL CELL.
*              THIS BOUNDARY CONDITION IS SUITABLE FOR ANY SPECIES OF
*              MULTI-SPECIES PLASMA (I.E. IN THE CASE WHEN THE CONDITION USED
*              IN BCCON=21 FAILS)
*              IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH
*              CORRESPONDING CONDITION ON ION TEMPERATURE (BCENI=26,27)

        case (25)
          call xerrab ('b2stbc_phys -- BCCON = 25 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7)')
     .       'BCCON = 25 : constant average density',CONPAR(IS,IB,1)
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                        BC_LIST_Y(1,IBW)),
     .           'BCCON = 25 only applicable to core boundaries !')
                call xertst(CONPAR(IS,IBW,1).gt.0.0_R8,
     1           'BCCON = 25, CONPAR(,,1) !> 0')
                call xertst(CONPAR(IS,IB,1).eq.CONPAR(IS,IBW,1),  ! IYS 22.01.2015
     1           'BCCON = 25, Average density on all parts ' //
     1             'of the surface should be the same')
                if (bceni(ibw).ne.26.and.bceni(ibw).ne.27) then
                  write(0,*) 'Should use BCCON(IS,IB).EQ.25 '//
     &             'together with BCENI=26 or 27'
                  write(*,*) 'Should use BCCON(IS,IB).EQ.25 '//
     &             'together with BCENI=26 or 27'
                  write(*,'(a)')
     &             'Please ensure that this boundary condition '//
     &             'is consistent with the '//
     &             'ion energy boundary condition'                  !srv 21.09.12
                endif
              endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.21.or.
     1            BCCON(IS,IBW).EQ.22.or.BCCON(IS,IBW).EQ.23.or.
     1            BCCON(IS,IBW).EQ.26.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab('BCCON types 16, 19, 21, 22, 23, 25, '//
     .                      '26 and 27 cannot be mixed!')
              endif
            END DO
            if (ishigh) call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 25')
          endif

          IF(BCCHAR(IB).EQ.'S') THEN
            if (.not.once_bccon(is)) then
              if(is.eq.ismain.and.use_astra.ne.0) then
#ifdef ASTRA
                DO IBW=IB,NBC
                  if(BCCHAR(IBW).EQ.'S'.and.
     &               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    CONPAR(IS,IBW,1) = neb_astra                       !srv 28.07.11
                  endif
                ENDDO
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              endif
              t0=0.0e0_R8   ! Volume between two surfaces              !lk 12.12.12 {
              t1=0.0e0_R8   ! Sum of volume-averaged density
              DO IBW=IB,NBC
                if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  do IX = 1, BC_LIST_SIZE(IBW)
                    t0=t0 +
     .                vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    t1=t1 +
     .                 na(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),is)
     .               *vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  enddo
                endif
              enddo
              t1=t1/t0
              DO IBW=IB,NBC
                do IX = 1, BC_LIST_SIZE(IBW)
                  if(on_closed_surface(BC_LIST_X(IX,IBW),
     &                                 BC_LIST_Y(IX,IBW))) then
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*(CONPAR(IS,IBW,1) +
     &                 na(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),IS)
     .                -t1)*
     &               SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  endif
                enddo                                                  !lk 12.12.12 }
                DO IX = 1, BC_LIST_SIZE(IBW)                           !srv 23.09.08
                  SNA0_NO_MDF(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)=
     &             SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                ENDDO
              ENDDO
              once_bccon(is) = .true.
            endif
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCCON=25, VALID FOR CORE ONLY')
          ENDIF      !lk and IYS 01.03.2013}
#endif
* -- BCCON=26  Feedback boundary condition: PRESCRIBE THE TOTAL ION FLUX AND
*              FIND THE AVERAGE DENSITY.  !lk and IYS 01.03.2013
*              THE POLOIDAL DENSITY VARIATION FROM THE NEIGHBOURING RADIAL
*              RING IS APPLIED ON TOP OF THE AVERAGE DENSITY.
*              THIS BOUNDARY CONDITION IS SUITABLE FOR ANY SPECIES OF
*              MULTI-SPECIES PLASMA (I.E. IN THE CASE WHEN THE CONDITION USED
*              IN BCCON=22 FAILS)
*              IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH
*              A CORRESPONDING CONDITION ON THE ION TEMPERATURE (BCENI=26,27)

        case (26)
          call xerrab ('b2stbc_phys -- BCCON = 26 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,a,i3)')                                     !IYS 21.01.2016
     .       'BCCON = 26 : given total particle flux with ',        !IYS 21.01.2016
     &       'constant average density for ion is =', is            !IYS 21.01.2016
!           write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                        BC_LIST_Y(1,IBW)),
     .           'BCCON = 26 only applicable to core boundaries !')
                if (bceni(ibw).ne.26.and.bceni(ibw).ne.27) then
                  write(0,*) 'Should use BCCON(IS,IB).EQ.26 '//
     &                       'together with BCENI=26 or 27'
                  write(*,*) 'Should use BCCON(IS,IB).EQ.26 '//
     &                        'together with BCENI=26 or 27'
                  write(*,'(a)')
     &           'Please ensure that this boundary condition '//
     &           'is consistent with the '//
     &           'ion energy boundary condition'                      !srv 21.09.12
                endif
              endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.21.or.
     1            BCCON(IS,IBW).EQ.22.or.BCCON(IS,IBW).EQ.23.or.
     1            BCCON(IS,IBW).EQ.25.or.BCCON(IS,IBW).EQ.27)) then  ! IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab('BCCON types 16, 19, 21, 22, 23, 25, '//
     .                      '26 and 27 cannot be mixed!')
              endif
            ENDDO
            if (ishigh) call xertst (switch%ionising_core.eq.0,
     .        'ionising_core switch usage inconsistent with BCCON = 26')
          endif  ! ncall_b2stbc_phys == 0
          IF(BCCHAR(IB).EQ.'S') THEN
            if(.not.once_bccon(is)) then
              if (use_astra.ne.0) then
#ifdef ASTRA
                SCONPAR(IS) = fneb_astra
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SCONPAR(IS) = 0.0_R8
                DO IBW=IB,NBC
                  if(BCCHAR(IBW).EQ.'S'.and.
     >               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    SCONPAR(IS) = SCONPAR(IS) + CONPAR(IS,IBW,1)
                  endif
                ENDDO
              endif
              write(*,*) 'BCCON=26: Summed value is ', SCONPAR(IS)
              totflux = 0.0_R8
              t0=0.0e0_R8   ! Volume between two surfaces               !lk 12.12.12 {
              t1=0.0e0_R8   ! Volume-averaged density on previous surface
              NAS=0.0e0_R8
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    t0=t0 +
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    t1=t1 +
     .                na(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),IS)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    NAS=NAS+na(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    if(mdf_fnb.ne.0) then
                      totflux = totflux +
     &        fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,is)
                    else
                      totflux = totflux +
     &        fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &            topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,is)
                    endif
                  ENDDO
                ENDIF
              ENDDO
              t1=t1/t0        ! Sum of volume-averaged electron and ion temperatures
              NAS=NAS/t0      ! Volume-averaged density
              if (max(abs(totflux),abs(SCONPAR(IS))).GT.0.0_R8)
     >         NAS = NAS*
     &           (1.0_R8-switch%b2stbc_bc_ref*(totflux-SCONPAR(IS))/
     &             max(abs(totflux),abs(SCONPAR(IS))))
              write(*,'(a,es14.7,a,es14.7,a,es14.7,a,i3)')
     .         ' NAS=',NAS,' totflux=',totflux,
     .         ' SCONPAR=',SCONPAR(IS),' IS=',IS
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*(NAS +
     &                 na(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),IS)
     &                 -t1)
     &               *SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))          ! lk
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  ENDDO
                  DO IX = 1, BC_LIST_SIZE(IBW)                         !srv 23.09.08 {
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)
     &               =SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  ENDDO                                                !srv 23.09.08 }
                ENDIF
              ENDDO
              once_bccon(is) = .true.
            endif                                                      !srv 16.02.12 }
          ELSE                                                         !srv 16.02.12 }
            CALL XERRAB
     1        ('B2STBC -- BCCON=26 IS AVAILABLE FOR CORE ONLY')
          ENDIF                                                        !IYS 01.03.2013 }
#endif
* -- BCCON=27  Feedback boundary condition: prescribe the particle flux sum
*              for all neutrals and ions belonging to a given isonuclear sequence
*              and find the average density of the highest ionisation state.
*              THE POLOIDAL DENSITY VARIATION FROM THE NEIGHBOURING RADIAL
*              CELL IS ADDED TO THE AVERAGE DENSITY.
*              IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH
*              A CORRESPONDING CONDITION ON THE ION TEMPERATURE (BCENI=26,27)

        case (27)
          call xerrab ('b2stbc_phys -- BCCON = 27 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,a,i3)')
     &       'BCCON = 27 : given total particle flux with constant',
     &       ' density for highest ionisation stage ion is =', is
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >         THEN
                write(*,'(a,a,a,a,i3)') ' on ', BCCHAR(IBW),
     .           boundary_location(ibw),' for species ',IS
                call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                        BC_LIST_Y(1,IBW)),
     .           'BCCON = 27 only applicable to core boundaries !')
                call xertst(CONPAR(IS,IBW,1).ge.0.0_R8,
     1           'BCCON = 27, CONPAR(,,1) < 0')
                if (bceni(ibw).ne.26.and.bceni(ibw).ne.27) then
                  write(0,*) 'Should use BCCON(IS,IB).EQ.27 '//
     &                       'together with BCENI=26 or 27'
                  write(*,*) 'Should use BCCON(IS,IB).EQ.27 '//
     &                       'together with BCENI=26 or 27'
                  write(*,'(a)')
     &             'Please ensure that this boundary condition '//
     &             'is consistent with the '//
     &             'ion energy boundary condition'                    !srv 21.09.12
              endif
              if(BCCHAR(IBW).EQ.'S'.and.(BCCON(IS,IBW).EQ.16.or.
     1            BCCON(IS,IBW).EQ.19.or.BCCON(IS,IBW).EQ.21.or.
     1            BCCON(IS,IBW).EQ.22.or.BCCON(IS,IBW).EQ.23.or.
     1            BCCON(IS,IBW).EQ.25.or.BCCON(IS,IBW).EQ.26)) then   !IYS 22.01.2016 added BCCON = 25,26 and 27
                call xerrab('BCCON types 16, 19, 21, 22, 23, 25, '//
     .                      '26 and 27 cannot be mixed!')
              endif
            ENDDO
            if (is.lt.ns-1) call xertst(.not.lnext(is,is+1),
     1       'BCCON = 27 is applied to highest ionisation stages only')
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            if(.not.once_bccon(is)) then
              if (use_astra.ne.0) then
#ifdef ASTRA
                SCONPAR(IS) = fSeNBI_astra
#else
                call xerrab
     .           ('Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SCONPAR(IS) = 0.0_R8
                DO IBW=IB,NBC
                  if(BCCHAR(IBW).EQ.'S'.and.
     >               BCCON(IS,IBW).EQ.BCCON(IS,IB)) THEN
                    SCONPAR(IS) = SCONPAR(IS) + CONPAR(IS,IBW,1)
                  endif
                ENDDO
              endif
              write(*,*) 'BCCON=27: Summed value is ', SCONPAR(IS)
              totfluxn = 0.0_R8
              totfluxi = 0.0_R8
              t0=0.0e0_R8   ! Volume between two surfaces              !lk 12.12.12 {
              t1=0.0e0_R8   ! Volume-averaged on neighbouring surface
              NAS=0.0e0_R8
#ifdef B25_EIRENE
cxpb  Add contribution from Eirene core boundary fluxes
              do isa = 1, nnatmi
                if (latmscl(isa).eq.b2espcr(is)) then
                  do k = 1, ncns
                    ix = nlim + lcns(k)
                    totfluxn = totfluxn -
     &                 (wldna(ix,isa,0) - wldra(ix,isa,0)
     &                                  - wldpa(ix,isa,0))
                    do ism = 1, nnmoli
                      if (mlcmp(isa,ism).gt.0) then
                        totfluxn = totfluxn - mlcmp(isa,ism) *
     &                   (wldnm(ix,ism,0) - wldrm(ix,ism,0)
     &                                    - wldpm(ix,ism,0))
                      end if
                    end do
                  end do
                end if
              end do
#endif
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
! totfluxn here means the sum of fluxes over all ions (including fluid
! neutrals) of this isotope except for the highest ionisation stage.
! For coupled runs, both fna and fna_mdf are rescaled for neutrals
! Note that for neutrals fna = fna_mdf
                    do isb = 0, ns-1
                      if(b2espcr(is).eq.b2espcr(isb).and.is.ne.isb) then
                        if(.not.is_neutral(isb).or.use_eirene.eq.0) then
                          if(mdf_fnb.ne.0) then
                            totfluxn = totfluxn +
     &       fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,isb)
                          else
                            totfluxn = totfluxn +
     &           fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,isb)
                          endif
                        endif
                      endif
                    enddo
! totfluxi here means the flux of the ion in the highest ionisation stage.
                    if(mdf_fnb.ne.0) then
                      totfluxi = totfluxi +
     &       fna_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,is)
                    else
                      totfluxi = totfluxi +
     &           fna(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1,is)
                    endif
                    t0=t0 +
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    t1=t1 +
     .                na(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   IS)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                    NAS=NAS+na(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),is)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  ENDDO
                ENDIF
              ENDDO
              t1=t1/t0        !
              NAS = NAS/t0    ! Volume-averaged density
              IF(abs(totfluxn)+abs(totfluxi)+abs(SCONPAR(IS)).GT.0.0_R8)
     >         NAS = NAS*
     &         (1.0_R8-switch%b2stbc_bc_ref*
     &           (totfluxn+totfluxi-SCONPAR(IS))/
     &           (abs(totfluxn)+abs(totfluxi)+abs(SCONPAR(IS))))
              write(*,'(5(a,es14.7))')
     *         ' NAS=',NAS,' totfluxn=',totfluxn,
     &         ' totfluxi=',totfluxi,' SCONPAR=',SCONPAR(IS),
     &         ' totfluxn+totfluxi=',totfluxn+totfluxi
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))
     >           THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0,IS)=
     &               VBIG*(NAS +
     &                 na(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                    topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),IS)
     &                -t1)
     &               *SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))          ! lk
                    SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1,IS)=
     &              -VBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  ENDDO
                  DO IX =1, BC_LIST_SIZE(IBW)                          !srv 23.09.08 {
                    SNA0_NO_MDF(BC_LIST_X(IX,IBW),
     .                          BC_LIST_Y(IX,IBW),:,IS)=
     &               SNA0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),:,IS)
                  ENDDO                                                !srv 23.09.08 }
                ENDIF
              ENDDO
              once_bccon(is) = .true.
            endif                                                      !srv 16.02.12 }
          ELSE                                                         !srv 16.02.12 }
            CALL XERRAB
     1       ('B2STBC -- BCCON=27 IS AVAILABLE FOR CORE ONLY')
          ENDIF                                                        !IYS 04.03.2013 }
#endif

* -- BCCON=28  Prescribe the flux surface averaged density, with poloidal
*              perturbation taken from flux tube just inside domain

          case (28)
          if(ncall_b2stbc_phys.eq.0) then
            do ibc = 1, mpg%bcCvP(ib,2)
             iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
             call xertst(mpg%cvOnClosedSurface(iCv1),
     .        'BCCON=28 intended for core boundaries only!')
            enddo
            if(is.eq.ismain.and.switch%use_astra.ne.0) then
              conpar(is,iB,1) = neb_astra                               !srv 28.07.11
            endif
            if (conpar(is,iB,1).eq.0.0_R8) then
              conpar(is,iB,1) = switch%b2mndr_na_min
              write(*,*) 'Zero density B.C. no longer allowed!'
              write(*,*) 'Changed to na_min.'
            endif
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'constant density + perturbation',conpar(is,iB,1),' on ',
     1       BCCHAR(iB),boundary_location(iB),' for species ',is
            call xertst(conpar(is,iB,1).gt.0.0_R8,
     1       'BCCON = 28, conpar(,,1) !> 0')
          endif
*       ..compute average density on core boundary
          t0 = 0.0e0_R8
          t1 = 0.0e0_R8
          do ib0 = 1, nbc
           ! identify core boundaries; assume all core boundaries
           ! form single flux surface
           ! implicit assumption at the moment: first ring in the
           ! core is a 'structured grid' ring with quad cells
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
           if (mpg%cvOnClosedSurface(iCv1)) then
            do ibc0 = 1, mpg%bcCvP(ib0,2)                               ! loop over number of cells in the boundary
             iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                 ! number of the guard cell
             iCv2 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,2)                 ! number of the domain cell
             iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
             t0 = t0 + geo%cvVol(iCv2)
             t1 = t1 + geo%cvVol(iCv2)*pl%na(iCv2,is)
            enddo
           endif
          enddo
          t1 = t1/t0
*       ..impose density with perturbation
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                 ! number of the domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%sna0(iCv1,0,is) =
     .        VBIG*geo%fcS(iFc)*(pl%na(iCv2,is) + conpar(is,iB,1) - t1)
            srw%sna0(iCv1,1,is) =
     .       -VBIG*geo%fcS(iFc)
            sna0_no_mdf(iCv1,0,is) = srw%sna0(iCv1,0,is)
            sna0_no_mdf(iCv1,1,is) = srw%sna0(iCv1,1,is)
          enddo

* -- BCCON=29  Feedback on the total particle flux, by imposing a flux surface
*              averaged density, with poloidal perturbation taken from flux
*              tube just inside domain

          case (29)
          if(ncall_b2stbc_phys.eq.0) then
            do ibc = 1, mpg%bcCvP(ib,2)
             iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
             call xertst(mpg%cvOnClosedSurface(iCv1),
     .        'BCCON=29 intended for core boundaries only!')
            enddo
          endif
*       ..compute average density on core boundary
          t0 = 0.0e0_R8
          t1 = 0.0e0_R8
          t2 = 0.0e0_R8
          totflux = 0.0e0_R8
          totpar  = 0.0e0_R8
          do ib0 = 1, nbc
           ! identify core boundaries; assume all core boundaries
           ! form single flux surface
           ! implicit assumption at the moment: first ring in the
           ! core is a 'structured grid' ring with quad cells
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
           if (mpg%cvOnClosedSurface(iCv1)) then
            if (bccon(is,ib0).ne.29)
     .        call xerrab ('BCCON=29: incompatible core BCs')
            do ibc0 = 1, mpg%bcCvP(ib0,2)                               ! loop over number of cells in the boundary
             iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                 ! number of the guard cell
             iCv2 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,2)                 ! number of the domain cell
             iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
             t0 = t0 + geo%cvVol(iCv2)
             t1 = t1 + geo%cvVol(iCv2)*pl%na(iCv1,is)
             t2 = t2 + geo%cvVol(iCv2)*pl%na(iCv2,is)
             if (bccon_flux_fb_style.eq.0) then
               totflux = totflux - mpg%bcFcOr(mpg%bcCvP(ib0,1)+ibc0-1)*
     &           (dv%fna(iFc,0,is) + dv%fna(iFc,1,is))
             else
               totflux = totflux - mpg%bcFcOr(mpg%bcCvP(ib0,1)+ibc0-1)*
     &           (dv%fna_eir(iFc,0,is) + dv%fna_eir(iFc,1,is))
             endif
            enddo
cnh conpar(is,ib0,3) contains the contribution from the ionized neutrals
c   (ionising_core option)
            totpar = totpar + conpar(is,ib0,1) + conpar(is,ib0,3)
           endif
          enddo
          t1 = t1/t0
          t2 = t2/t0
          if (conpar(is,ib,2).gt.0.0_R8) t3 = conpar(is,ib,2)
          if (conpar(is,ib,2).le.0.0_R8) t3 = 1.0e-2_R8
          t4 = (1.0_R8-t3*(totflux-totpar)/
     &         max(abs(totflux),abs(totpar),1.0e10_R8))*t1
*       ..impose density with perturbation
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                    ! number of the domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%sna0(iCv1,0,is) =
     .        VBIG*geo%fcS(iFc)*(pl%na(iCv2,is) + t4 - t2)
            srw%sna0(iCv1,1,is) =
     .       -VBIG*geo%fcS(iFc)
            sna0_no_mdf(iCv1,0,is) = srw%sna0(iCv1,0,is)
            sna0_no_mdf(iCv1,1,is) = srw%sna0(iCv1,1,is)
          enddo

        case default
          write(*,*) BCCON(IS,IB)
          call xerrab ('BCCON out of bounds')
        end select

* -- BOUNDARY CONDITIONS FOR PARALLEL MOMENTUM EQUATIONS --

        select case (BCMOM(IS,IB))
        case (0)
          continue

* -- BCMOM=1 -- PRESCRIBE THE VALUE OF THE PARALLEL VELOCITY

        case (1)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCMOM =  1 : specified parallel velocity ',
     1        MOMPAR(IS,IB,1),' on ',BCCHAR(IB),boundary_location(ib),
     2       ' for species ',IS
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%smo0(iCv1,0,is)=
     .         VBIG*NE00*am(is)*mompar(is,ib,1)*geo%fcS(iFc)*
     .         geo%fcHz(iFc)
            srw%smo0(iCv1,1,is)=
     .        -VBIG*NE00*am(is)*geo%fcS(iFc)*geo%fcHz(iFc)
          enddo

* -- BCMOM=2 -- PRESCRIBE THE GRADIENT OF THE PARALLEL VELOCITY

        case (2)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCMOM =  2 : parallel velocity gradient ',MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                    ! number of the corresponding domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            iVx1 = mpg%fcVx(iFc,1)
            iVx2 = mpg%fcVx(iFc,2)
            t1 = intvertex_s(iVx1,nCv,nVx,mpg,geo%vxVol,pl%ua(:,is))
            t2 = intvertex_s(iVx2,nCv,nVx,mpg,geo%vxVol,pl%ua(:,is))
            srw%smo0(iCv1,0,is)=
     .         VBIG*am(is)*pl%na(iCv1,is)*geo%fcS(iFc)*geo%fcHz(iFc)*
     .         (mompar(is,ib,1)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2))*
     .          geo%fcQgam(iFc,0) - cor9*(t2 - t1)*geo%fcQgam(iFc,1)*
     .          (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     .          mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1) + pl%ua(iCv2,is))
            srw%smo0(iCv1,1,is)=
     .        -VBIG*am(is)*pl%na(iCv1,is)*geo%fcS(iFc)*geo%fcHz(iFc)
          enddo

* -- BCMOM=3 -- SHEATH CONDITIONS, MACH NUMBER AS INPUT

        case (3)
          call xerrab ('b2stbc_phys -- BCMOM = 3 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCMOM =  3 : sheath velocity bc ',MOMPAR(IS,IB,2),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              IF(MOMPAR(IS,IB,2).LT.0.5_R8) THEN
                CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))
                UUOUT=MOMPAR(IS,IB,1)*CS*
     .           abs(wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .               wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3))
                if(.not.is_neutral(is)) then                           !xpb
                  UPOUT=max(0.0_R8, UUOUT/
     .             abs(wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .                 wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)) +    !xpb
     .             bcint(ud(-1,-1,is),nx,ny,iy,ib)*
     .                (wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/
     .                 wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)))     !xpb
                else                                                   !xpb
                  UPOUT = UUOUT/
     .             abs(wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .                 wbbr(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3))      !xpb
                endif                                                  !xpb
cxpb  Added multiplicative factor to UD according to the Chodura criterion
cxpb  and max(,,) statement for stability
              ELSE
                ZBAR=bcint(na(-1,-1,is),nx,ny,iy,ib)/
     .           bcint(NE,nx,ny,IY,IB)
                ZARG=(GAMMAI*bcint(ti,nx,ny,iy,ib)
     .           +bcint(rz2(-1,-1,is),nx,ny,iy,ib)*ZBAR
     .           *bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                CS=SQRT(ZARG)
                UUTEST=abs(
     .           ua(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .           (1.0_R8+
     .            HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .    HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .       RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .    ua(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .       RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .            HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .    HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .       RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))
                UUOUT=MAX(MOMPAR(IS,IB,1)*CS,UUTEST)
                if(.not.is_neutral(is)) then                           !xpb
                  UPOUT=max(0.0_R8, UUOUT+                             !xpb
     .             bcint(ud(-1,-1,is),nx,ny,iy,ib)*
     .                (wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/
     .                 wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)))     !xpb
                else                                                   !xpb
                  UPOUT = UUOUT                                        !xpb
                endif                                                  !xpb
cxpb  Added multiplicative factor to UD according to the Chodura criterion
cxpb  and max(,,) statement for stability
              ENDIF
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         (UBIG/EPSLON**2)*AM(IS)*mp*
     1         bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .         SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     2         (-UPOUT)*
     .         SIGN(1.0_R8,BB(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))    !csc
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .         -(UBIG/EPSLON**2)*AM(IS)*mp*
     .         bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .         SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            end do
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)
              IF(MOMPAR(IS,IB,2).LT.0.5_R8) THEN
                CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))
                UUOUT=MOMPAR(IS,IB,1)*CS*
     .           abs(wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .               wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3))
                if(.not.is_neutral(is)) then                           !xpb
                  UPOUT=max(0.0_R8, UUOUT/                             !xpb
     .              abs(wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .                  wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3))     !xpb
     .             -bcint(ud(-1,-1,is),nx,ny,iy,ib)*
     .                 (wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/
     .                  wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)))    !xpb
                else                                                   !xpb
                  UPOUT = UUOUT/
     .             abs(wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)/
     .                 wbbl(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3))      !xpb
                endif                                                  !xpb
cxpb  Added multiplicative factor to UD according to the Chodura criterion
cxpb  and max(,,) statement for stability
              ELSE
                ZBAR=bcint(na(-1,-1,is),nx,ny,iy,ib)/
     .           bcint(NE,nx,ny,IY,IB)
                ZARG=(GAMMAI*bcint(ti,nx,ny,iy,ib)
     .           +bcint(rz2(-1,-1,is),nx,ny,iy,ib)*ZBAR
     .           *bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                CS=SQRT(ZARG)
                UUTEST=abs(
     .           ua(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .           (1.0_R8+
     .            HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .       HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .       ua(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .            HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .       HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))
                UUOUT=MAX(MOMPAR(IS,IB,1)*CS,UUTEST)
                if(.not.is_neutral(is)) then                           !xpb
                  UPOUT=max(0.0_R8, UUOUT -                            !xpb
     .             bcint(ud(-1,-1,is),nx,ny,iy,ib)*
     .                (wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),2)/
     .                 wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)))     !xpb
                else                                                   !xpb
                  UPOUT = UUOUT                                        !xpb
                endif                                                  !xpb
cxpb  Added multiplicative factor to UD according to the Chodura criterion
cxpb  and max(,,) statement for stability
              ENDIF
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         (UBIG/EPSLON**2)*AM(IS)*mp*
     1         bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     2         SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))*UPOUT*
     .         SIGN(1.0_R8,BB(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))    !csc
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=
     .        -(UBIG/EPSLON**2)*AM(IS)*mp*
     .         bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .         SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCMOM=3, NO SHEATH IN RADIAL DIR.')
          END IF
#endif
* -- BCMOM=4 -- PRESCRIBE THE VALUE OF THE VELOCITY, WEAKLY
*                 A MIXED BOUNDARY CONDITION

        case (4)
          call xerrab ('b2stbc_phys -- BCMOM = 4 not adapted for WG')
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCMOM =  4 : weak velocity ',MOMPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(MOMPAR(IS,IB,2).gt.0.0_R8,
     1       'BCMOM = 4, MOMPAR(,,2) !> 0')
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%smo0(iCv1,0,is)=
     .         MOMPAR(IS,IB,2)*mompar(is,ib,1)*geo%fcS(iFc)*
     .         geo%fcHz(iFc)
            srw%smo0(iCv1,1,is)=
     .        -MOMPAR(IS,IB,2)*geo%fcS(iFc)*geo%fcHz(iFc)
          enddo

* -- BCMOM=5 -- PRESCRIBE THE PARALLEL MOMENTUM FLUX PER UNIT AREA

        case (5)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCMOM =  5 : parallel momentum flux per unit area ',
     .        MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
cwdk      Possible backwards-compatibility issue:
cwdk      In the unstructured solver, we impose flux per unit area normal
cwdk      to the surface/boundary.
cwdk      In the original code, the flux per unit area in the poloidal
cwdk      direction was imposed at E/W boundaries (and normal to the
cwdk      radial direction at N/S - unchanged for aligned grids)
          do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            srw%smo0(iCv1,0,is) = mompar(is,iB,1)*geo%fcS(iFc)*
     .        geo%fcHz(iFc)
            srw%smo0(iCv1,1,is) = 0.0_R8
          enddo

* -- BCMOM=6 -- PRESCRIBE THE TOTAL PARALLEL MOMENTUM FLUX FOR A
*                 CONSTANT PARALLEL VELOCITY

        case (6)
          CALL XERRAB ('B2STBC -- BCMOM=6, NOT YET AVAILABLE')

* -- BCMOM=7 -- PRESCRIBE THE GIVEN PROFILE OF VELOSITY

        case (7)
          if(ncall_b2stbc_phys.eq.0) then                               !srv 26.02.18 {
            write(*,'(4a,i3)')
     .      'BCMOM = 7 : specified parallel velocity on ',
     .       BCCHAR(IB), boundary_location(ib), ' for species ',IS
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%smo0(iCv1,0,is)=
     .         VBIG*NE00*am(is)*bv_ua(iCv1,is)*geo%fcS(iFc)*
     .         geo%fcHz(iFc)
            srw%smo0(iCv1,1,is)=
     .        -VBIG*NE00*am(is)*geo%fcS(iFc)*geo%fcHz(iFc)
          enddo

* -- BCMOM=8 -- SPECIAL : LIMITED SHEAR

        case (8)
          call xerrab ('b2stbc_phys -- BCMOM = 8 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,a,a,a,i3)')
     .       'BCMOM =  8 : limited shear velocity on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              A=MIN(1.0_R8,
     .         SQRT(TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .              TI(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                 TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*AM(IS) *
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)*A*
     .         UA(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        - VBIG*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*AM(IS) *
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            end do
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              A=MIN(1.0_R8,
     .         SQRT(TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     .          TI(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .             BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)= VBIG*
     .         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .         AM(IS) * na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)*A*
     .         UA(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .            BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),IS)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=-VBIG*
     .         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .         AM(IS) * na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCMOM=8, NOT AVAILABLE FOR W,E')
          END IF
#endif
* -- BCMOM=9 -- PRESCRIBE THE TOTAL PARALLEL MOMENTUM FLUX WITH
*                 CONSTANT FLUX DENSITY

        case (9)
          call xerrab ('b2stbc_phys -- BCMOM = 9 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCMOM =  9 : total parallel momentum flux ',
     .       MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            US=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              US=US+SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            ENDDO
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS) =
     1         MOMPAR(IS,IB,1)/US*
     1         SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            US=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              US=US+SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .                   BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            ENDDO
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS) =
     1         MOMPAR(IS,IB,1)/US*
     2         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=0.0_R8
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'W') THEN
            US=0.0_R8
            DO IY = 1, BC_LIST_SIZE(IB)
              US=US+SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            ENDDO
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         MOMPAR(IS,IB,1)/US*
     1         SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            US=0.0_R8
            DO IY = 1, BC_LIST_SIZE(IB)
              US=US+SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            ENDDO
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         MOMPAR(IS,IB,1)/US*
     1         SXHZ(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'X') THEN
            US=0.0_R8
            DO IX = 1, BC_LIST_SIZE(IB)
              US=US+SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            ENDDO
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         MOMPAR(IS,IB,1)/US*
     1         SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'Y') THEN
            US=0.0_R8
            DO IY = 1, BC_LIST_SIZE(IB)
              US=US+SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            ENDDO
            DO IY = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,IS)=
     .         MOMPAR(IS,IB,1)/US*
     1         SXHZ(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              SMO0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1,IS)=0.0_R8
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCMOM=9, NOT AVAILABLE')
          ENDIF
#endif
* -- BCMOM=10 -- PRESCRIBE THE DECAY LENGTH FOR THE PARALLEL MOMENTUM

        case (10)
          call xerrab ('b2stbc_phys -- BCMOM = 10 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCMOM = 10 : parallel momentum decay length ',
     .       MOMPAR(IS,IB,1),' on ',
     1       BCCHAR(IB),boundary_location(ib),' for species ',IS
            call xertst(MOMPAR(IS,IB,1).ne.0.0_R8,
     1       'BCMOM = 10, MOMPAR(,,1) = 0')
            if (bccon(is,ib).ne.9) then
              write(0,*)
     &         'Should use BCCON(IS,IB).EQ.9 together with BCMOM=10'
              write(*,*)
     &         'Should use BCCON(IS,IB).EQ.9 together with BCMOM=10'
              write(*,'(1x,a)')
     &                   'Please ensure that this boundary condition '//
     &                   'is consistent with '//
     &                   'the density boundary condition'  !srv 21.09.12
            endif
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)=
     .         AM(IS)*mp*ua(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
     .         *SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -vsa0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)/
     1         MOMPAR(IS,IB,1)*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     .         AM(IS)*mp*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            ENDDO
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         AM(IS)*mp*ua(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
     .         *SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     .        -vsa0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)/
     1         MOMPAR(IS,IB,1)*
     2         SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .              BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))) +
     .         AM(IS)*mp*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            ENDDO
          ELSE
            CALL XERRAB ('B2STBC -- BCMOM=10, NOT AVAILABLE FOR W,E')
          ENDIF
#endif
* -- BCMOM=11 -- ROZHANSKY VISCOSITY CONDITION FOR THE PARALLEL MOMENTUM

        case (11)
          call xerrab ('b2stbc_phys -- BCMOM = 11 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCMOM = 11 : ROZ parallel momentum bc on ',BCCHAR(IB),
     1       boundary_location(ib),
     2       ' for species ',IS
            write(*,*) ' BCMOM=11 no longer supported !'
            write(0,*) ' BCMOM=11 no longer supported !'
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              IF (BC_LIST_X(IX,IB).GE.0) THEN
                smo0(bc_list_x(ix,ib),bc_list_y(ix,ib),0,is)=
     1           vsa0(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .           na(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*am(is)*
     2           sqrt(1.0_R8-pit(bc_list_x(ix,ib),bc_list_y(ix,ib))**2)*
     .           ua(bc_list_x(ix,ib),bc_list_y(ix,ib),is)
                smo0(bc_list_x(ix,ib),bc_list_y(ix,ib),1,is)=
     1           -dna0(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*am(is)*
     2           (na(bc_list_x(ix,ib),bc_list_y(ix,ib),is)-
     2            na(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),is))/
     3           hy1(bc_list_x(ix,ib),bc_list_y(ix,ib))-
     3           vsa0(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .           na(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     4           am(is)*
     .           sqrt(1.0_R8-pit(bc_list_x(ix,ib),bc_list_y(ix,ib))**2)/
     .           hy1(bc_list_x(ix,ib),bc_list_y(ix,ib))
              end if
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCMOM=11, NOT AVAILABLE FOR N,W,E')
          ENDIF
#endif
* -- BCMOM=12 -- CONDITION from b2stbc_spb FOR THE PARALLEL MOMENTUM   !srv 01.02.09 {
!srv added accumulation in order to account twice b.c. for corner cells
        case (12)
          call xerrab ('b2stbc_phys -- BCMOM = 12 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a,a,i3)')
     .       'BCMOM = 12 : from b2stbc_spb parallel momentum bc '//
     .       'with mompar1 =', MOMPAR(IS,IB,1),
     1       ' and mompar2 =', MOMPAR(IS,IB,2),' on ',
     2       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              s0hz = hz(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &              vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &              hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))             !srv 09.01.01
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),2,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3,IS)=
     &         SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3,IS)+
     &             MOMPAR(IS,IB,2)*s0hz +
!     &         MOMPAR(IS,IB,2)*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &         MOMPAR(IS,IB,1)*
     &         max(0.0e0_R8,-fna_fcor(
     &          TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &          TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1,is))/
     &          na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
            END DO
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              s0hz = hz(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &              vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &              hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))             !srv 09.01.01
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),2,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3,IS)=
     &         SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3,IS)+
     &         MOMPAR(IS,IB,2)*s0hz+
!     &           MOMPAR(IS,IB,2)*
!     &           SYHZ(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
!     &                BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))+
     &         MOMPAR(IS,IB,1)*
     &         max(0.0e0_R8,
     &          fna_fcor(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1,is))/
     &                na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
            end do
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCMOM=12, NOT AVAILABLE FOR W,E')
          ENDIF
#endif
* -- BCMOM=13 -- CONDITION from b2stbc_spb FOR THE PARALLEL MOMENTUM
!srv added accumulation in order to account twice b.c. for corner cells
        case (13)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a,a,i3)')
     .       'BCMOM = 13 : from b2stbc_spb parallel momentum bc '//
     .       'with mompar1 =', MOMPAR(IS,IB,1),
     1       ' and mompar2 =', MOMPAR(IS,IB,2),' on ',
     2       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            cs   = sqrt(pz(iCv1)/rz(iCv1))
            s1hz = geo%fcPbshz(iFc)
!            vbnd = cs*s1hz*mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1) -
!     &        dv%vaecrb(iFc,0,is) -
!     &        dv%vaecrb(iFc,1,is)*geo%fcQalf(iFc,1)/geo%fcQalf(iFc,0)
            vbnd = cs*s1hz*mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1) -
     &        (dv%vaecrb(iFc,0,is)*geo%fcQalf(iFc,0) +
     &         dv%vaecrb(iFc,1,is)*geo%fcQalf(iFc,1))*
     &        sign(1.0_R8,geo%fcQalf(iFc,0))*geo%fcHz(iFc)*geo%fcS(iFc)
            srw%smo0(iCv1,0,is) = srw%smo0(iCv1,0,is) +                        !srv 07.05.09 {
     &         (UBIG/EPSLON**2)*am(is)*mp*pl%na(iCv1,is)*
     &         mompar(is,ib,1)*vbnd
            srw%smo0(iCv1,1,IS) = srw%smo0(iCv1,1,IS) -
     &         (UBIG/EPSLON**2)*am(is)*mp*pl%na(iCv1,is)*abs(s1hz)        !srv 07.05.09 }
            ! Avoid setting to sound speed at very shallow incidence
            ! Set mom. flux to leakage instead
            if (abs(geo%fcQalf(iFc,0)*geo%cvBb(iCv1,0)/geo%cvBb(iCv1,3))
     &          .lt.switch%Qalfmin) then
              srw%smo0(iCv1,0,is) = 0.0_R8
              srw%smo0(iCv1,1,is) = 0.0_R8
              if (conpar(is,ib,2).gt.0.0_R8) then
                srw%smo0(iCv1,1,is) = srw%smo0(iCv1,1,is) -
     &            mompar(is,ib,2)/conpar(is,ib,2)*am(is)*mp*
     &            co%dna0(iCv1,is)*pl%na(iCv1,is)*geo%fcS(iFc)*
     &            geo%fcHz(iFc)
              endif
            endif
          enddo                                                        !srv 01.02.09 }

* -- BCMOM=14 -- CONDITION from b2stbc_spb FOR THE PARALLEL MOMENTUM  !srv 01.02.09 {
!srv added accumulation in order to account twice b.c. for corner cells
        case (14)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a,a,i3)')
     .       'BCMOM = 14 : from b2stbc_spb momentum bc with mompar1 =',
     1       MOMPAR(IS,IB,1),' and mompar2 =', MOMPAR(IS,IB,2),' on ',
     2       BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%smo0(iCv1,0,is) = 0.0_R8 !TODO: drift contributions!
            srw%smo0(iCv1,1,IS) = 0.0_R8
            srw%smo0(iCv1,2,IS) = 0.0_R8
            srw%smo0(iCv1,3,IS) = srw%smo0(iCv1,3,IS) +
     &         MOMPAR(IS,IB,2)*geo%fcS(iFc)*geo%fcHz(iFc) +
     &         MOMPAR(IS,IB,1)*
     &         max(0.0e0_R8,mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &         dv%fna_fcor(iFc,1,is))/pl%na(iCv1,is)
          enddo

#ifdef WG_TODO
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              s0hz = hz(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &              vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &              hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))              !srv 09.01.01
              if (.not.is_neutral(is)) then                             !srv 15.02.02 {
                select case (iform_corr_dia)                            !srv 25.07.05
                case (0)
                  smo0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &             smo0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+        !srv 09.01.01 01.07.05
     &             s0hz*(am(is)*mp)*naaverage(is)*                      !srv 15.02.02
     &            (-2.0e0_R8*
     &             facdrift(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*phm2*    !srv 07.06.00 17.05.99 25.04.00 15.02.02
     &             tiaverage*bb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),2)/
     &             (qe*rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &             hx(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*              !srv 15.02.02
     &            (OnedBsqc(rightix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                      rightiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))-
     &              OnedBsqc(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))        !srv 05.03.99 15.02.02
     &             )*uaaverage(is) + corr_flux*s0hz                     !srv 15.02.02 08.09.03 25.07.05
                case (1)
                  smo0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is) =
     &             smo0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,is)+        !srv 09.01.01 01.07.05
     &             s0hz*(am(is)*mp)*naaverage(is)*                      !srv 15.02.02
     &           (-2.0e0_R8*facdrift(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &             phm2*                                                !srv 07.06.00 17.05.99 25.04.00 15.02.02
     &             tiaverage*bb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),2)/
     &            (qe*rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &             hx(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*              !srv 15.02.02
     &            (OnedBsqc(rightix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                      rightiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))-
     &             OnedBsqc(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))         !srv 05.03.99 15.02.02
     &            )*ua(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)+corr_flux
     &             *s0hz                                                !srv 08.09.03 25.07.05
                end select                                              !srv 25.07.05
              endif                                                     !srv 15.02.02
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),2,IS)= 0.0_R8
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3,IS)=
     &         SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3,IS)+
     &         MOMPAR(IS,IB,2)*s0hz +
!     &         MOMPAR(IS,IB,2)*SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &         MOMPAR(IS,IB,1)*
     &         max(0.0e0_R8,
     &      -fna_fcor(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1,is))/
     &                   na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
            END DO
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCMOM=14, NOT AVAILABLE FOR N,W,E')
          ENDIF
#endif

* -- BCMOM=15 -- PRESCRIBE THE VALUE OF THE PARALLEL VELOCITY, SCALED WITH
*                B_average/B_local

        case (15)
          call xerrab ('b2stbc_phys -- BCMOM = 15 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,a,1p,g14.7,a,a,a,a,i3)')
     .       'BCMOM = 15 : specified parallel velocity as a function',
     &       ' of poloidal coordinate ',MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)=
     .         VBIG*NE00*AM(IS)*MOMPAR(IS,IB,1)*
     .         (bbaverage/bb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3))*
     &         SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))                  !srv 24.08.08
              SMO0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)=
     &        -VBIG*NE00*AM(IS)*
     .         SYHZ(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))                  !srv 24.08.08
            END DO
          ELSE
            CALL XERRAB ('B2STBC -- BCMOM=15, ONLY FOR S')
          ENDIF                                                         !srv 02.04.09 }
#endif
* -- BCMOM=16 -- PRESCRIBE THE AVERAGE VALUE OF THE PARALLEL VELOCITY   !lk 12.12.12 {

        case (16)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCMOM = 16 : specified average parallel velocity ',
     .       MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
           do ibc = 1, mpg%bcCvP(ib,2)
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
            call xertst(mpg%cvOnClosedSurface(iCv1),
     .       'BCMOM=16 intended for core boundaries only!')
           enddo
          endif
*       ..compute average parallel velocity on core boundary
          t0 = 0.0e0_R8
          t1 = 0.0e0_R8
          do ib0 = 1, nbc
           ! identify core boundaries; assume all core boundaries
           ! form single flux surface
           ! implicit assumption at the moment: first ring in the
           ! core is a 'structured grid' ring with quad cells
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
           if (mpg%cvOnClosedSurface(iCv1)) then
            do ibc0 = 1, mpg%bcCvP(ib0,2)                               ! loop over number of cells in the boundary
             iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                 ! number of the guard cell
             iCv2 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,2)                 ! number of the domain cell
             iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
             t0 = t0 + geo%cvVol(iCv2)
             t1 = t1 + geo%cvVol(iCv2)*pl%ua(iCv2,is)
            enddo
           endif
          enddo
          t1 = t1/t0
*       ..impose velocity
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                    ! number of the domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%smo0(iCv1,0,is)=
     .         VBIG*NE00*am(is)*geo%fcS(iFc)*geo%fcHz(iFc)*
     .         (mompar(is,ib,1) + pl%ua(iCv2,is) - t1)
            srw%smo0(iCv1,1,is)=
     .        -VBIG*NE00*am(is)*geo%fcS(iFc)*geo%fcHz(iFc)
          enddo

* -- BCMOM=17 -- leakage option for parallel momentum

        case (17)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCMOM = 17 : parallel momentum leakage ',MOMPAR(IS,IB,1),
     1       ' on ',BCCHAR(IB),boundary_location(ib),' for species ',IS
          endif
cwdk      Possible backwards-compatibility issue:
cwdk      In the unstructured solver, we impose flux per unit area normal
cwdk      to the surface/boundary.
cwdk      In the original code, the flux per unit area in the poloidal
cwdk      direction was imposed at E/W boundaries (and normal to the
cwdk      radial direction at N/S - unchanged for aligned grids)
          do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            srw%smo0(iCv1,0,is) = 0.0_R8
            srw%smo0(iCv1,1,is) = geo%fcS(iFc)*geo%fcHz(iFc)*
     .        sqrt((rt%rza(iCv1,is)*pl%te(iCv1)+pl%ti(iCv1))/
     .        (am(is)*mp))*(am(is)*mp)*mompar(is,ib,1)*pl%na(iCv1,is)
          enddo


* -- BCMOM=18 -- PRESCRIBE ZERO SOURCE TERMS (-> no flux)

        case (18)
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a,a,i3)')
     .       'BCMOM = 18 : specified parallel velocity ',
     .        MOMPAR(IS,IB,1),' on ',BCCHAR(IB),boundary_location(ib),
     1       ' for species ',IS
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            srw%smo0(iCv1,0,is)= 0.0_R8
            srw%smo0(iCv1,1,is)= 0.0_R8

          enddo

        case default
          write(*,*) BCMOM(IS,IB)
          call xerrab ('BCMOM out of bounds')

        end select
* -- END LOOP OVER SPECIES
      enddo

* -- BOUNDARY CONDITIONS FOR ELECTRON ENERGY EQUATION --

      select case (BCENE(IB))

      case (0)
        continue

* -- BCENE=1 -- PRESCRIBE THE VALUE OF THE ELECTRON TEMPERATURE

      case (1)
        if(ncall_b2stbc_phys.eq.0) then
          IF(BCCHAR(IB).EQ.'S'.and.switch%use_astra.ne.0) THEN
            ENEPAR(IB,1) = teb_astra*1.0e3_R8                        !srv 14.06.12
          ENDIF
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENE =  1 : specified electron temperature ',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENEPAR(IB,1).gt.0.0_R8,
     1     'BCENE = 1, ENEPAR(,1) !> 0')
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          srw%she0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*enepar(iB,1)*EV
          srw%she0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
        enddo

* -- BCENE=2 -- PRESCRIBE THE GRADIENT OF THE ELECTRON TEMPERATURE

      case (2)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENE =  2 : electron temperature gradient ',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                    ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
          iVx1 = mpg%fcVx(iFc,1)
          iVx2 = mpg%fcVx(iFc,2)
          t1 = intvertex_s(iVx1,nCv,nVx,mpg,geo%vxVol,pl%te)
          t2 = intvertex_s(iVx2,nCv,nVx,mpg,geo%vxVol,pl%te)
          srw%she0(iCv1,0)=  VBIG*NE00*geo%fcS(iFc)*
     .       (EV*enepar(ib,1)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2))*
     .        geo%fcQgam(iFc,0) - cor9*(t2 - t1)*geo%fcQgam(iFc,1)*
     .        (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     .        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1) + pl%te(iCv2))
          srw%she0(iCv1,1)= -VBIG*NE00*geo%fcS(iFc)
        enddo

* -- BCENE=3 -- SHEATH CONDITIONS, ELECTRON ENERGY TRANSMISSION

      case (3)
        call xerrab ('b2stbc_phys -- BCENE = 3 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,g14.7,a,g14.7,a,a,a)')
     .     'BCENE =  3 : sheath electron with MOMPAR(ismain,IB,2) =',
     1     MOMPAR(ismain,IB,2),' and MOMPAR(ismain,IB,1) =',
     2     MOMPAR(ismain,IB,1),' and POTPAR(IB,2) =',POTPAR(IB,2),
     3     ' and ENEPAR(IB,1) =',ENEPAR(IB,1),
     4     ' on ',BCCHAR(IB),boundary_location(ib)
          if(fix_fch_in_fhe_sheath.eq.0 .and.
     &       fix_sx_in_fhe_sheath.eq.0) then !srv 08.05.09
            write(*,*) 'WARNING: Wrong (old) implementation of ',
     2       'the FCH contribution to FHE used'
          else if (fix_fch_in_fhe_sheath.eq.1 .and.
     &       fix_sx_in_fhe_sheath.eq.0) then
            write(*,*) 'WARNING: Wrong (second) implementation of ',
     2       'the FCH contribution to FHE used'
          endif
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
              CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))  !srv 01.02.09
              UUOUT=abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .                  wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3))*
     .         MOMPAR(ismain,IB,1)*CS
              if(fix_fch_in_fhe_sheath.eq.0) then
                NUEOUT=MIN(0.0_R8,-UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,0)/QE)
              elseif(fix_fch_in_fhe_sheath.eq.1 .and.
     .               fix_sx_in_fhe_sheath.eq.0) then
                NUEOUT=MIN(0.0_R8,-UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,0)/QE)
              elseif(fix_fch_in_fhe_sheath.eq.2 .or.
     .               fix_sx_in_fhe_sheath.eq.1) then                     !srv 08.05.09 {
                NUEOUT=MIN(0.0_R8,-UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,0)/QE/
     .            sx(bc_list_x(iy,ib),bc_list_y(iy,ib)))
              endif
            ELSE
              NUEOUT=0.0_R8
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  ZBAR=bcint(na(-1,-1,IS),nx,ny,IY,IB)/
     .             bcint(NE,nx,ny,IY,IB)
                  ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .             bcint(rz2(-1,-1,IS),nx,ny,IY,IB)*ZBAR*
     .             bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                  CS=SQRT(ZARG)
                  UUTEST=abs(
     .             UU(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .             (1.0_R8+
     .              HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .      HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .         RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .      UU(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .         RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .              HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .      HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .         RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))
                  UUOUT=MAX( MOMPAR(IS,IB,1)*CS*
     .             abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .                 wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .             UUTEST)
                  NUEOUT=NUEOUT - UUOUT*bcint(rza(-1,-1,IS),nx,ny,IY,IB)
     .             *bcint(na(-1,-1,IS),nx,ny,IY,IB)
                endif
              end do
              DO IS=0,NS_EXT-1
                if (.not.is_ext_neutral(is)) then
                  NUEOUT=NUEOUT +
     .             UU(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS+NS)
     .             *bcint(za_ext(-1,-1,IS),nx,ny,IY,IB)
     .             *bcint(na_ext(-1,-1,IS),nx,ny,IY,IB)
                endif
              end do
              if(fix_fch_in_fhe_sheath.eq.0) then
                NUEOUT=MIN(0.0_R8, NUEOUT
     .           - fch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,0)/QE )
              elseif(fix_fch_in_fhe_sheath.eq.1 .and.
     .               fix_sx_in_fhe_sheath.eq.0) then
                NUEOUT=MIN(0.0_R8, NUEOUT -
     &           fch(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,0)/QE)
              elseif(fix_fch_in_fhe_sheath.eq.2 .or.
     .               fix_sx_in_fhe_sheath.eq.1) then
                NUEOUT=MIN(0.0_R8, NUEOUT -
     &           fch(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,0)/QE/
     &            sx(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)) )
              endif
            ENDIF
            if(enepar(ib,1).gt.0.0_R8) then                             !xpb
              she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .        she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) +
     .         (ENEPAR(IB,1)+
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT*
     .          SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            else           ! xpb new pressure linearisation as per B2.5-style
              she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)=
     .        she0(bc_list_x(iy,ib),bc_list_y(iy,ib),3) +               !xpb
     .         (-ENEPAR(IB,1)+                                          !xpb
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT*
     .          SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     .          ne(bc_list_x(iy,ib),bc_list_y(iy,ib))
            endif                                                       !xpb
            if(switch%mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &         she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &         fhePSch(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                 rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,0)
            endif                                                !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
              CS=SQRT(bcint(pz,nx,ny,IY,IB)/bcint(rz,nx,ny,IY,IB))  !srv 01.02.09
              UUOUT=MOMPAR(ismain,IB,1)*CS*
     .         abs(wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .             wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),3))
              if(fix_fch_in_fhe_sheath.eq.0) then
                NUEOUT=MAX(0.0_R8,UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,0)/QE)
              elseif(fix_fch_in_fhe_sheath.eq.1 .and.
     .               fix_sx_in_fhe_sheath.eq.0) then
                NUEOUT=MAX(0.0_R8,UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,0)/QE)
              elseif(fix_fch_in_fhe_sheath.eq.2 .or.
     .               fix_sx_in_fhe_sheath.eq.1) then
                NUEOUT=MAX(0.0_R8,UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,0)/QE/
     .            sx(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
              endif
            ELSE
              NUEOUT=0.0_R8
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  ZBAR=bcint(na(-1,-1,IS),nx,ny,IY,IB)/
     .             bcint(NE,nx,ny,IY,IB)
                  ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .             bcint(rz2(-1,-1,IS),nx,ny,IY,IB)*ZBAR*
     .             bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                  CS=SQRT(ZARG)
                  UUTEST=abs(
     .             UU(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .             (1.0_R8+
     .              HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .       HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .       UU(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .              HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .       HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))
                  UUOUT=MAX( MOMPAR(IS,IB,1)*CS*
     .             abs(wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .                 wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .             UUTEST)
                  NUEOUT=NUEOUT + UUOUT*bcint(rza(-1,-1,IS),nx,ny,IY,IB)
     .             *bcint(na(-1,-1,IS),nx,ny,IY,IB)
                endif
              end do
              DO IS=0,NS_EXT-1
                if (.not.is_ext_neutral(is)) then
                  NUEOUT=NUEOUT +
     .              UU(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),IS+NS)
     .             *bcint(za_ext(-1,-1,IS),nx,ny,IY,IB)
     .             *bcint(na_ext(-1,-1,IS),nx,ny,IY,IB)
                endif
              end do
              if(fix_fch_in_fhe_sheath.eq.0) then
                NUEOUT=MAX(0.0_R8, NUEOUT -
     .           fch(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,0)/QE )
              elseif(fix_fch_in_fhe_sheath.eq.1 .and.
     .               fix_sx_in_fhe_sheath.eq.0) then
                NUEOUT=MAX(0.0_R8,
     .          NUEOUT - fch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,0)/QE )
              elseif(fix_fch_in_fhe_sheath.eq.2 .or.
     .               fix_sx_in_fhe_sheath.eq.1) then
                NUEOUT=MAX(0.0_R8,
     .           NUEOUT - fch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,0)/QE/
     .           SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))) )
              endif
            ENDIF
            if(enepar(ib,1).gt.0.0_R8) then                             !xpb
              she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .        she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) -
     .         (ENEPAR(IB,1) +
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT*
     .          SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .             LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            else           ! xpb new pressure linearisation as per B2.5-style
              she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)=
     .        she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3) -              !xpb
     .         (-ENEPAR(IB,1) +                                        !xpb
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT*
     .          SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .             LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
     .         /ne(bc_list_x(iy,ib),bc_list_y(iy,ib))
            endif                                                      !xpb
            if(switch%mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &         she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) +
     &         fhePSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,0)
            endif                                                !srv 12.09.08 }
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=3, NO SHEATH IN RADIAL DIR.')
        END IF
#endif
* -- BCENE=4 -- PRESCRIBE THE VALUE OF THE ELECTRON TEMPERATURE,
*                 WEAKLY A MIXED BOUNDARY CONDITION

      case (4)
        call xerrab ('b2stbc_phys -- BCENE = 4 not adapted for WG')
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENE =  4 : weak electron temperature ',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENEPAR(IB,1).gt.0.0_R8,
     1     'BCENE = 4, ENEPAR(,1) !> 0')
          call xertst(ENEPAR(IB,2).gt.0.0_R8,
     1     'BCENE = 4, ENEPAR(,2) !> 0')
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          srw%she0(iCv1,0) =  ENEPAR(IB,2)*geo%fcS(iFc)*enepar(ib,1)
          srw%she0(iCv1,1) = -ENEPAR(IB,2)/EV*geo%fcS(iFc)
        enddo

* -- BCENE=5 -- PRESCRIBE THE ELECTRON ENERGY FLUX PER UNIT AREA

      case (5)
        call xerrab ('b2stbc_phys -- BCENE = 5 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENE =  5 : electron energy flux per area ',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENEPAR(IB,1)*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENEPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       ENEPAR(IB,1)*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=ENEPAR(IB,1)*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=5, NOT YET AVAILABLE')
        END IF
#endif
* -- BCENE=6 -- PRESCRIBE THE TOTAL ELECTRON ENERGY FLUX FOR A
*                 CONSTANT ELECTRON TEMPERATURE

      case (6)
          call xerrab ('b2stbc_phys -- BCENE = 6 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENE =  6 : total electron energy flux ',ENEPAR(IB,1),
     1     ' with const Te on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'N') THEN
          Tguess=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            Tguess=Tguess+TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
          Tguess=Tguess/BC_LIST_SIZE(IB)
          if(switch%b2stbc_diagno.ge.2)
     >     write(*,*) 'Initial Tguess:',Tguess/ev
          T0=Tguess*0.95_R8
          T1=Tguess*1.01_R8
          Tlimit=Tguess*0.0001_R8
          fylimt=ENEPAR(IB,1)*0.0001_R8
          ibnd=-1
20610     ibnd=ibnd+1
          if (ibnd.eq.0) then
            Tguess=T0
          elseif (ibnd.eq.1) then
            Tguess=T1
          else
            Tguess=T1-GAM1*(T1-T0)/(EPSLON + GAM1 - GAM0)
          endif
          fyt=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
cxpb code segment analog to b2tfhe
            ehy = -(po(bc_list_x(ix,ib),bc_list_y(ix,ib))-
     &     po(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .        bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))+
     &    (ne(bc_list_x(ix,ib),bc_list_y(ix,ib))*Tguess-
     &     ne(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .        bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))*
     &     te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .        bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))*
     &    (2.0_R8/qe)/(ne(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                    bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &                +ne(bc_list_x(ix,ib),bc_list_y(ix,ib)))
            if (conduction_only.eq.1) then
              flo = 0.0_R8
            else
              prefactor=
     &         (vol(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &          vol(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &         /(hy1(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &           hy1(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &         *(ne(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &          +ne(bc_list_x(ix,ib),bc_list_y(ix,ib)))/2.0_R8
              flo_vedia=prefactor*
     .         vedia(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)
              flo_veecrb=prefactor*
     .         veecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)
              flo = (1.5_R8+switch%BoRiS)*
     .         fne_32(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)+
     1         2.5_R8*fne_52(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)+
     &       chve(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)*
     .       (ne(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     &       ne(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .          bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))/2.0_R8
              if(switch%mdf_fhe.eq.0) then
                flo = flo+2.5_R8*flo_vedia+1.5_R8*flo_veecrb
              else
                flo = flo+1.5_R8*flo_veecrb
              endif
            endif
            if (no_hybr_e.eq.0) then
              cone0 =
     .         hybr(flo,chce(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1))
              fy=flo*(te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &         Tguess)/2.0_R8-cone0*(Tguess-
     &         te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .            bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))+
     &         switch%alfTeEh*calf(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)*ehy*
     &         harm(te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     &         Tguess)
            else
              cone0 = chce(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)
              t2 = dim ( cone0,
     .         abs(flo) *
     .         max(0.0_R8,
     .          min(Tguess,
     .              te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))))/
     .         max(cutlo, Tguess+
     .          te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))) )
              fy =
     .         upwind(flo,
     .                te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     .                Tguess)
     .         - t2 * (Tguess-te(bottomix(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)),
     .                           bottomiy(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)))) +
     &         switch%alfTeEh*calf(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)*
     .         harm(te(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     &                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     &                            Tguess)*ehy
            endif
            fyt=fyt + fy
          ENDDO
          if(switch%b2stbc_diagno.ge.2)
     >     write(*,'(a,i4,f15.6,1pe14.6)') 'ibnd,Tg,fyt(e):',
     .     ibnd,Tguess/ev,fyt
          if (ibnd.eq.0) then
            GAM0 = fyt - ENEPAR(IB,1)
            goto 20610
          elseif (ibnd.eq.1) then
            GAM1 = fyt - ENEPAR(IB,1)
            goto 20610
          else
            T0 = T1
            GAM0 = GAM1 - ENEPAR(IB,1)
            T1 = Tguess
            GAM1 = fyt
            if(switch%b2stbc_diagno.ge.2)
     >       write(*,*) 'GAM0, GAM1 = ',GAM0,GAM1
            if(switch%b2stbc_diagno.ge.2)
     >       write(*,*) 'fyt,t1,t0 = ',fyt,t1/ev,t0/ev
            if (ibnd .gt. 100)
     >       call xerrab ('B2STBC -- BCENE=6, NOT CONVERGING')
            if (abs(t1-t0).gt.Tlimit .and.
     &          abs(fyt-ENEPAR(IB,1)).gt.fylimt) goto 20610
          endif
          if(switch%b2stbc_diagno.ge.1)
     >     write(*,'(a,f15.6,i4,1pe14.6)') ' tebound:',T1/ev,ibnd,fyt
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*T1*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'S') THEN
          Tguess=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            Tguess=Tguess+TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
          Tguess=Tguess/BC_LIST_SIZE(IB)
          if(switch%b2stbc_diagno.ge.2)
     >     write(*,*) 'Initial Tguess:',Tguess/ev
          T0=Tguess*0.95_R8
          T1=Tguess*1.01_R8
          Tlimit=Tguess*0.0001_R8
          fylimt=ENEPAR(IB,1)*0.0001_R8
          ibnd=-1
20620     ibnd=ibnd+1
          if (ibnd.eq.0) then
            Tguess=T0
          elseif (ibnd.eq.1) then
            Tguess=T1
          else
            Tguess=T1-GAM1*(T1-T0)/(EPSLON + GAM1 - GAM0)
          endif
          fyt=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
cxpb code segment analog to b2tfhe
            ehy = -(po(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     &                po(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &         (ne(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))*
     &          te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     &          ne(bc_list_x(ix,ib),bc_list_y(ix,ib))*Tguess)*
     &         (2.0_R8/qe)/(ne(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                         topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &                     +ne(bc_list_x(ix,ib),bc_list_y(ix,ib)))
            if (conduction_only.eq.1) then
              flo = 0.0_R8
            else
              prefactor=
     &         (vol(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &          vol(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &         /(hy1(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &           hy1(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &         *(ne(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))                    !xpb
     &          +ne(bc_list_x(ix,ib),bc_list_y(ix,ib)))/2.0_R8                   !xpb
              flo_vedia=prefactor*
     .         vedia(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)
              flo_veecrb=prefactor*
     .         veecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)
              flo = (1.5_R8+switch%BoRiS)*
     .         fne_32(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1)+
     1         2.5_R8*fne_52(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1)+
     &         chve(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1)*
     &        (ne(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     &         ne(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .            topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))/2.0_R8
              if(switch%mdf_fhe.eq.0) then
                flo = flo+2.5_R8*flo_vedia+1.5_R8*flo_veecrb
              else
                flo = flo+1.5_R8*flo_veecrb
              endif
            endif
            if (no_hybr_e.eq.0) then
              cone0 = hybr(flo,
     &         chce(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1))
              fy=flo*(te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &         Tguess)/2.0_R8-cone0*
     &         (te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-Tguess)+
     &         switch%alfTeEh*calf(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1)*ehy*
     &         harm(te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),Tguess)
            else
              cone0 = chce(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                     topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1)
              t2 = dim ( cone0,
     .         abs(flo) *
     .         max(0.0_R8,
     .          min(Tguess,
     .              te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))))/
     .         max(cutlo, Tguess+
     .          te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))) )
              fy =
     .         upwind(flo, Tguess,
     .                te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .         - t2 * (te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                    topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     .                 Tguess) + switch%alfTeEh *
     &         calf(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1)*
     &         harm(te(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     &              Tguess)*ehy
            endif
            fyt=fyt + fy
          ENDDO
          if(switch%b2stbc_diagno.ge.2)
     >     write(*,'(a,i4,f15.6,1pe14.6)') 'ibnd,Tg,fyt(e):',
     .     ibnd,Tguess/ev,fyt
          if (ibnd.eq.0) then
            GAM0 = fyt - ENEPAR(IB,1)
            goto 20620
          elseif (ibnd.eq.1) then
            GAM1 = fyt - ENEPAR(IB,1)
            goto 20620
          else
            T0 = T1
            GAM0 = GAM1 - ENEPAR(IB,1)
            T1 = Tguess
            GAM1 = fyt
            if (ibnd .gt. 100)
     >       call xerrab ('B2STBC -- BCENE=6, NOT CONVERGING')
            if (abs(t1-t0).gt.Tlimit .and.
     &          abs(fyt-ENEPAR(IB,1)).gt.fylimt) goto 20620
          endif
          if(switch%b2stbc_diagno.ge.1)
     >     write(*,'(a,f15.6,i4,1pe14.6)') ' tebound:',T1/ev,ibnd,fyt
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=6, NOT AVAILABLE FOR W,E')
        ENDIF
#endif
* -- BCENE=7 -- PRESCRIBE THE ELECTRON TEMPERATURE PROFILE

      case (7)
        if(ncall_b2stbc_phys.eq.0) then                               !srv 26.02.18 {
          write(*,'(3a)')
     .     'specified electron temperature on ',
     1     BCCHAR(IB), boundary_location(ib)
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
          srw%she0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*bv_te(iCv1)*EV
          srw%she0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
        enddo

* -- BCENE=8 -- PRESCRIBE THE TOTAL ELECTRON HEAT FLUX WITH CONSTANT
*                 FLUX DENSITY
!wdk BC not recommended for drifts!

      case (8)
        if(ncall_b2stbc_phys.eq.0) then
          if (switch%use_astra.ne.0) then
#ifdef ASTRA
            ENEPAR(IB,1) = fheb_astra                                 !srv 28.07.11
#else
            call xerrab
     .       ('Compile with -DASTRA option to couple to ASTRA!')
#endif
          endif
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENE =  8 : total electron energy flux ',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
        endif
        US = 0.0_R8
        do ibc = 1, mpg%bcCvP(ib,2)
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          US   = US + geo%fcS(iFc)
        enddo
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          srw%she0(iCv1,0) = enepar(iB,1)/US*geo%fcS(iFc)
          srw%she0(iCv1,1) = 0.0_R8
          if(switch%mdf_fhe.ne.0) then
            !wdk Correction now involves vedia instead of ExB velocity
            !wdk This makes implementation equivalent to mdf_f = 0 option
            srw%she0(iCv1,0) = srw%she0(iCv1,0) +
     &        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &       (dv%fhePSch(iFc,0) + dv%fhePSch(iFc,1))
            t0 = mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*2.5_R8*
     &        (geo%fcS(iFc)*geo%fcQalf(iFc,0)*dv%vedia(iFc,0) +
     &         geo%fcS(iFc)*geo%fcQalf(iFc,1)*dv%vedia(iFc,1))*
     &         0.5_R8*(dv%ne(iCv1)+dv%ne(iCv2))*
     &         0.5_R8*(pl%te(iCv1)+pl%te(iCv2))
            srw%she0(iCv1,0) = srw%she0(iCv1,0) +
     &        max(0.0_R8,t0)
            srw%she0(iCv1,1) = srw%she0(iCv1,1) +
     &        min(0.0_R8,t0)/pl%te(iCv1)
          endif
        enddo

* -- BCENE=9 -- PRESCRIBE THE DECAY LENGTH FOR THE ELECTRON TEMPERATURE

      case (9)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENE =  9 : electron temperature decay length ',
     1     ENEPAR(IB,1),' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENEPAR(IB,1).ne.0.0_R8,
     1     'BCENE = 9, ENEPAR(,1) = 0')
          ism = 0
          do is = 0, ns-1
            if (bccon(is,ib).ne.9.and..not.is_neutral(is)) then
              ism = ism + 1
              write(0,*)
     &         'Should use BCCON(IS,IB).EQ.9 together with BCENE=9'
              write(*,*)
     &         'Should use BCCON(IS,IB).EQ.9 together with BCENE=9'
            endif
            if (bccon(is,ib).eq.1) then
              write(0,*)
     &         'Should use BCENE=19 in conjunction with BCCON=1'
            end if
            if (ism.gt.0) then
              write(*,'(1x,a)')
     &         'Please ensure that this boundary condition '//
     &         'is consistent with the density boundary conditions'      !srv 21.09.12
            end if
          enddo
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          sy = geo%fcS(iFc)
          srw%she0(iCv1,0) = 0.0_R8
          srw%she0(iCv1,1) = -co%hce0(iCv1)/enepar(ib,1)*
     &                    sy
          do is = 0, ns-1
            srw%she0(iCv1,0) = srw%she0(iCv1,0) + rt%rza(iCv1,is)*
     &                     2.5_R8*pl%te(iCv1)*srw%sna0(iCv1,0,is)
            srw%she0(iCv1,1) = srw%she0(iCv1,1) + rt%rza(iCv1,is)*
     &                     2.5_R8*srw%sna0(iCv1,1,is)*pl%na(iCv1,is)
          enddo !! end species loop
          if(switch%mdf_fhe.ne.0) then
            srw%she0(iCv1,0) = srw%she0(iCv1,0) +
     &        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &       (dv%fhePSch(iFc,0) + dv%fhePSch(iFc,1))
            t0 = -mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*1.5_R8*
     &        (geo%fcS(iFc)*geo%fcQalf(iFc,0)*dv%veecrb(iFc,0) +
     &         geo%fcS(iFc)*geo%fcQalf(iFc,1)*dv%veecrb(iFc,1))*
     &         0.5_R8*(dv%ne(iCv1)+dv%ne(iCv2))*
     &         0.5_R8*(pl%te(iCv1)+pl%te(iCv2))
            srw%she0(iCv1,0) = srw%she0(iCv1,0) +
     &        max(0.0_R8,t0)
            srw%she0(iCv1,1) = srw%she0(iCv1,1) +
     &        min(0.0_R8,t0)/pl%te(iCv1)
          endif
        enddo   !! end loop over cells
#ifdef WG_TODO
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -hce0(bc_list_x(ix,ib),bc_list_y(ix,ib))/ENEPAR(IB,1)*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            DO IS=0,NS-1
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .         rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*2.5_R8*
     .         TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .        she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     .         rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*2.5_R8*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            ENDDO
            if(switch%mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &         fhePSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1)
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     &         1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &         veecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1)*
     &                   ne(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            endif                                                !srv 12.09.08 }
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -hce0(bc_list_x(ix,ib),bc_list_y(ix,ib))/ENEPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            DO IS=0,NS-1
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .         rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*2.5_R8*
     .         TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .        she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     .         rza(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*2.5_R8*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*
     .         na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            ENDDO
            if(switch%mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     .         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &         fhePSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1)
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     .         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) -
     &         1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &         veecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1)*
     &             ne(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            endif                                                !srv 12.09.08 }
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=9, NOT AVAILABLE FOR W,E')
        ENDIF
#endif
* -- BCENE=10 -- Feedback option for core

      case (10)
          call xerrab ('b2stbc_phys -- BCENE = 10 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,a,a)')
     .     'BCENE = 10 : electron heat flux feedback bc on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(lfeedback, 'LFEEDBACK must be turned on!')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if (region(bc_list_x(ix,ib),bc_list_y(ix,ib),0).eq.1) then
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .         cbshe(0,coreregno)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            elseif (region(bc_list_x(ix,ib),bc_list_y(ix,ib),0).eq.5)
     .       then
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .         cbshe(0,coreregn2)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            endif
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=10, NOT AVAILABLE FOR N,W,E')
        ENDIF
#endif
* -- BCENE=11 -- FROM B2STBC_SPB                                       !srv 01.02.09 {
!srv added accumulation in order to account twice b.c. for corner cells
      case (11)
          call xerrab ('b2stbc_phys -- BCENE = 11 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'BCENE = 11 : from b2stbc_spb '//
     1     'with ENEPAR(,1) =',ENEPAR(IB,1),
     1     ' and ENEPAR(,2) =',ENEPAR(IB,2),
     2     ' on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
!            vte = sqrt(bcint(te,nx,ny,ix,ib)/me)
            vte = sqrt(te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/me)
            s0 = vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            if(switch%mdf_fhe.ne.0) then                                    !srv 21.03.09 {
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &          fhePSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1)
            endif                                                    !srv 21.03.09 }
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &       she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)+ENEPAR(IB,1)*
     &       max(0.0e0_R8,-floe_noc(
     &         topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &         topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1))       !srv 28.10.02
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)=
     &       she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)+ENEPAR(IB,2)*
     &       vte*s0
!     &       vte*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
!            vte = sqrt(bcint(te,nx,ny,ix,ib)/me)
            vte = sqrt(te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/me)
            s0 = vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            if(switch%mdf_fhe.ne.0) then                             !srv 21.03.09 {
              she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &         fhePSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1)
          endif                                                      !srv 21.03.09 }
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &       she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)+ENEPAR(IB,1)*
     &       max(0.0e0_R8,
     &        floe_noc(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1))       !srv 28.10.02
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)=
     &       she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),3)+ENEPAR(IB,2)*
     &       vte*s0
!     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
!     2          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=11, NOT AVAILABLE FOR W,E')
        ENDIF                                                         !srv 01.02.09 }
#endif
* -- BCENE=12 -- SHEATH CONDITIONS, EL ENERGY TRANSMISSION COEFFICIENT

      case (12)
        call xerrab ('b2stbc_phys -- BCENE = 12 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,a,a)')
     .     'BCENE = 12 : sheath electron bc on ',BCCHAR(IB),
     1     boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(bc_list_x(ix,ib),bc_list_y(ix,ib),3) =
     .      she0(bc_list_x(ix,ib),bc_list_y(ix,ib),3)-
     1       ENEPAR(IB,1)*max(0.0_R8,
     1        -fne(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1))/
     2                ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(bc_list_x(ix,ib),bc_list_y(ix,ib),3) =
     .      she0(bc_list_x(ix,ib),bc_list_y(ix,ib),3)-
     1             ENEPAR(IB,1)*max(0.0_R8,
     .       fne(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1))/
     2        ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(bc_list_x(iy,ib),bc_list_y(iy,ib),3) =
     .      she0(bc_list_x(iy,ib),bc_list_y(iy,ib),3)-
     1       ENEPAR(IB,1)*max(0.0_R8,
     1        -fne(rightix(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .             rightiy(bc_list_x(iy,ib),bc_list_y(iy,ib)),0,0))/
     2                  ne(bc_list_x(iy,ib),bc_list_y(iy,ib))
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(bc_list_x(iy,ib),bc_list_y(iy,ib),3) =
     .      she0(bc_list_x(iy,ib),bc_list_y(iy,ib),3)-
     1             ENEPAR(IB,1)*max(0.0_R8,
     .       fne(bc_list_x(iy,ib),bc_list_y(iy,ib),0,0))/
     2        ne(bc_list_x(iy,ib),bc_list_y(iy,ib))
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=12, INVALID OPTION')
        END IF
#endif
* -- BCENE=13 -- PRESCRIBE THE ELECTRON ENERGY FLUX PER UNIT AREA
*             -- PROPORTIONAL TO TEMPERATURE

      case (13)
          call xerrab ('b2stbc_phys -- BCENE = 13 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     1     'BCENE = 13 : electron energy flux per area '//
     .     'proportional to temperature ',
     2     ENEPAR(IB,1),' on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       ENEPAR(IB,1)*TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       ENEPAR(IB,1)*TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       ENEPAR(IB,1)*TE(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     1       SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .       ENEPAR(IB,1)*TE(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=13, NOT YET AVAILABLE')
        END IF
#endif
* -- BCENE=14 -- leakage option for electron energy

      case (14)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENE = 14 : electron energy leakage ',ENEPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          cs   = sqrt(pz(iCv1)/rz(iCv1))
          srw%she0(iCv1,0) = geo%fcS(iFc)*cs*dv%ne(iCv1)*
     &      enepar(ib,1)*enepar(ib,2)*pl%po(iCv1)*qe
          srw%she0(iCv1,3) = geo%fcS(iFc)*cs*enepar(ib,1)
          if(switch%mdf_fhe.ne.0) then
            srw%she0(iCv1,0) = srw%she0(iCv1,0) +
     &        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &       (dv%fhePSch(iFc,0) + dv%fhePSch(iFc,1))
            t0 = mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*2.5_R8*
     &        (geo%fcS(iFc)*geo%fcQalf(iFc,0)*dv%vedia(iFc,0) +
     &         geo%fcS(iFc)*geo%fcQalf(iFc,1)*dv%vedia(iFc,1))*
     &         0.5_R8*(dv%ne(iCv1)+dv%ne(iCv2))*
     &         0.5_R8*(pl%te(iCv1)+pl%te(iCv2))
            srw%she0(iCv1,0) = srw%she0(iCv1,0) +
     &        max(0.0_R8,t0)
            srw%she0(iCv1,1) = srw%she0(iCv1,1) +
     &        min(0.0_R8,t0)/pl%te(iCv1)
          endif
        enddo

* -- BCENE=15 -- FROM B2STBC_SPB                                        !srv 01.02.09 {
      case (15)
        if(ncall_b2stbc_phys.eq.0) then                                 !srv 20.09.17 {
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENE = 15 (style = 0) : with GAMMAE =',GAMMAE,
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          write(*,*) 'BCENE=15 - assuming BCCON=14 for all ions'
#ifdef WG_TODO
          IF(BCCHAR(IB).EQ.'W') THEN
            BCCON14_W(IB) = .TRUE.  ! IF BCCON=14 HAS BEEN IMPOSED AT WEST BOUNDARY IB FOR ALL IONS
            DO IS = 0, NS-1
              IF(.NOT.IS_NEUTRAL(IS)) THEN
                IF(BCCON(IS,IB).NE.14) THEN
                  BCCON14_W(IB) = .FALSE.
                ENDIF
              ENDIF
            ENDDO
            if(BCCON14_W(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 HAS ',
     &      'BEEN IMPOSED AT WEST BOUNDARY', IB,' FOR ALL IONS'
            if(.not.BCCON14_W(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 ',
     &      'HAS NOT BEEN IMPOSED AT WEST BOUNDARY', IB,' FOR ALL IONS'
          ENDIF
          IF(BCCHAR(IB).EQ.'E') THEN
            BCCON14_E(IB) = .TRUE.  ! IF BCCON=14 HAS BEEN IMPOSED AT WEST BOUNDARY IB FOR ALL IONS
            DO IS = 0, NS-1
              IF(.NOT.IS_NEUTRAL(IS)) THEN
                IF(BCCON(IS,IB).NE.14) THEN
                  BCCON14_E(IB) = .FALSE.
                ENDIF
              ENDIF
            ENDDO
            if(BCCON14_E(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 HAS ',
     &       'BEEN IMPOSED AT EAST BOUNDARY', IB,' FOR ALL IONS'
            if(.not.BCCON14_E(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 ',
     &       'HAS NOT BEEN IMPOSED AT EAST BOUNDARY', IB,' FOR ALL IONS'
          ENDIF
#endif
        endif                                                           !srv 20.09.17 }
!srv added accumulation in order to account twice b.c. for corner cells
        if(bcene_15_style.eq.0) then
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a)')
     .       'BCENE = 15 (style = 0) : with GAMMAE =',GAMMAE,
     1       ' on ',BCCHAR(IB),boundary_location(ib)
          endif
          do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            cs   = sqrt(pz(iCv1)/rz(iCv1))
            t0 = 0.0_R8
            t1 = 0.0_R8
            do is = 0, ns-1
              t0 = t0 + rt%rza(iCv1,is)*conpar(is,ib,1)*
     *          pl%na(iCv1,is)
              if (.not.is_neutral(is) .and. conpar(is,ib,2).gt.0.0_R8)
     &          t1 = t1 + co%dna0(iCv1,is)/conpar(is,ib,2)*
     &            pl%na(iCv1,is)
            enddo
            t0 = max(0.0e0_R8, (abs(geo%fcPbs(iFc))*t0*cs -           !srv 05.05.09
     &        dv%fch_p(iFc,0)*mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)/qe))  !srv 25.12.06 11.01.13 20.09.17
            t1 = t1*geo%fcS(iFc)
            if(switch%pot_eq.eq.1) then                               !srv 27.02.13 {
              srw%she0(iCv1,1) = srw%she0(iCv1,1) -                   !srv 11.03.09 {
     &          (1.0e0_R8+GAMMAE)/(1.0e0_R8-GAMMAE)*t0 -              !srv 19.04.13
     &           enepar(ib,2)*t1                                      !wdk to check: correct for gammae?
              srw%she0(iCv1,0) = srw%she0(iCv1,0) -
     &           qe*pl%po(iCv1)*(t0 + t1)                               !srv 07.11.12
            else                                                        !srv 27.02.13
              srw%she0(iCv1,1) = srw%she0(iCv1,1) -                     !srv 11.03.09 {
     &          ((1.0e0_R8+GAMMAE)/(1.0e0_R8-GAMMAE)+switch%delpo)*t0 - !srv 19.04.13
     &          (enepar(ib,2) + switch%delpo)*t1                       !wdk to check: correct for gammae?
            endif                                                       !srv 27.02.13 }
            if(switch%mdf_fhe.ne.0) then
              srw%she0(iCv1,0) = srw%she0(iCv1,0) +
     &          mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &         (dv%fhePSch(iFc,0) + dv%fhePSch(iFc,1))
              t0 = mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*2.5_R8*
     &          (geo%fcS(iFc)*geo%fcQalf(iFc,0)*dv%vedia(iFc,0) +
     &           geo%fcS(iFc)*geo%fcQalf(iFc,1)*dv%vedia(iFc,1))*
     &           0.5_R8*(dv%ne(iCv1)+dv%ne(iCv2))*
     &           0.5_R8*(pl%te(iCv1)+pl%te(iCv2))
              srw%she0(iCv1,0) = srw%she0(iCv1,0) +
     &          max(0.0_R8,t0)
              srw%she0(iCv1,1) = srw%she0(iCv1,1) +
     &          min(0.0_R8,t0)/pl%te(iCv1)
            endif
          enddo
cwdk    Todo: corrections for BCCON.ne.14
        elseif(bcene_15_style.eq.1) then
!srv added accumulation in order to account twice b.c. for corner cells
        call xerrab ('BCENE_15_STYLE 1 not adapted for WG')
#ifdef WG_TODO
          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,a,a)')
     &       'BCENE = 15 (style = 1) : with GAMMAE =',GAMMAE,
     1       ' on ',BCCHAR(IB),boundary_location(ib)
          endif
          IF(BCCHAR(IB).EQ.'W') THEN
            DO IY = 1, BC_LIST_SIZE(IB)                                 !srv 20.09.17
!              CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))
              s1=abs(pbs(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                   rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)) !srv 18.04.13 {
              vte=sqrt(te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/me)
              if(bccon14_w(ib)) then                                    !srv 20.09.17 {
*    .. apply bcene=15 according to the form of bccon=14
                t0 = 0.0_R8
                do is = 0, ns-1
                  t0 = t0 + rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*
     *            CONPAR(IS,IB,1)*
     *            na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
                enddo
              else
                t0 = ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              endif                                                     !srv 20.09.17
              t0 =s1*t0*sqrt(1.0_R8/(2.0_R8*pi))*vte*                   !srv 20.09.17
     &          expu2(max(-50.0_R8,
     &         -qe*(po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))-potpar(ib,2))/
     &              te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
              if(switch%pot_eq.eq.1) then                               !srv 27.02.13 {
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-
     &            (1.0e0_R8+GAMMAE)*t0
                 she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &            (1.0e0_R8-GAMMAE)*
     &            qe*po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*t0
              else                                                      !srv 27.02.13
                cs=sqrt(pz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &                  rz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-
     &               s1*((1.0e0_R8+GAMMAE)+(1.0e0_R8-GAMMAE)*delpo)*
     &           ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*cs
              endif                                                     !srv 27.02.13 }
              if(switch%mdf_fhe.ne.0) then                                     !srv 01.06.06 {
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &           fhePSch(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                   rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,0)
            endif                                                       !srv 01.06.06 }
            END DO
          ELSEIF(BCCHAR(IB).EQ.'E') THEN
            DO IY = 1, BC_LIST_SIZE(IB)                                 !srv 20.09.17
!              CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))
              s1=abs(pbs(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0))          !srv 18.04.13 {
              vte=sqrt(te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/me)
              if(bccon14_e(ib)) then                                    !srv 20.09.17 {
*    .. apply bcene=15 according to the form of bccon=14
                t0 = 0.0_R8
                do is = 0, ns-1
                  t0 = t0+rza(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)*
     &            CONPAR(IS,IB,1)*
     &            na(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),is)
                enddo
              else
                t0 = ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              endif                                                     !srv 20.09.17 }
              t0 =s1*t0*sqrt(1.0_R8/(2.0_R8*pi))*vte*                   !srv 20.09.17
     &          expu2(max(-50.0_R8,
     &         -qe*(po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))-potpar(ib,2))/
     &              te(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
              if(switch%pot_eq.eq.1) then                                !srv 27.02.13 {
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-
     &               (1.0e0_R8+GAMMAE)*t0
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &           (1.0e0_R8-GAMMAE)*qe*
     &           po(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*t0
              else                                                      !srv 27.02.13
                cs=sqrt(pz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     &                  rz(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)-
     &               s1*((1.0e0_R8+GAMMAE)+(1.0e0_R8-GAMMAE)*delpo)*
     &           ne(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*cs
              endif                                                     !srv 27.02.13 }
              if(switch%mdf_fhe.ne.0) then                                     !srv 01.06.06 {
                she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &           she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) +
     &           fhePSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,0)         !srv 20.09.17 }
            endif                                                       !srv 01.06.06 }
            end do
          ELSE
            CALL XERRAB
     .       ('B2STBC -- BCENE=15, not available for S and N')
          END IF
#endif
        else
          CALL XERRAB
     .     ('B2STBC -- BCENE=15, incorrect bcene_15_style value')
        endif                                                           !srv 01.02.09 }

* -- BCENE=16 -- Feedback boundary condition with constant temperature

      case (16)

        if(switch%bcene_16_style.eq.0) then

          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .       'BCENE = 16 (style = 0) : '//
     .       'electron temperature feedback by flux ',
     1       ENEPAR(IB,1),' with ',ENEPAR(IB,2),' on ',
     1       BCCHAR(IB),boundary_location(ib)
            call xertst(enepar(ib,1).ne.0, 'faulty parameter enepar1')
            call xertst(enepar(ib,2).ge.0, 'faulty parameter enepar2')
          endif
*       ..compute average temperature on boundary
          t0 = 0.0e0_R8
          t1 = 0.0e0_R8
          totflux = 0.0e0_R8
          !wdk Note on backwards compatibility: feedback now always on
          !wdk flux through the boundary itself (no offset)
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
           iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                     ! number of the guard cell
           iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
           t0 = t0 + geo%fcS(iFc)*dv%ne(iCv1)
           t1 = t1 + geo%fcS(iFc)*dv%ne(iCv1)*pl%te(iCv1)
           totflux = totflux - (dv%fhe(iFc,0) + dv%fhe(iFc,1))*
     &       mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)
          enddo
          t1 = t1/t0
          totpar = enepar(ib,1)
          if (enepar(ib,2).gt.0.0_R8) t3 = enepar(ib,2)
          if (enepar(ib,2).le.0.0_R8) t3 = 1.0e-2_R8
          t2 = (1.0_R8-t3*(totflux-totpar)/
     &         max(abs(totflux),abs(totpar)))*t1
          if (switch%b2stbc_diagno.ge.2) then
           write(*,'(a,1p,4(1x,e15.6))')
     &      'bcene=16: curr_flux, desired_flux, old/new temperature, ',
     &      totflux,totpar,t1/qe,t2/qe
          endif
*       ..impose temperature
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%she0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*t2
            srw%she0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
          enddo
        elseif(switch%bcene_16_style.eq.1)then
          call xerrab ('bcene_16_style.eq.1 not adapted for WG')
#ifdef WG_TODO

          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a)') 'BCENE = 16 (style = 1) : '//
     &       'Given integrated electron heat flux with constant Te '
!          write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCENE(IBW).EQ.BCENE(IB)) THEN
                write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                    boundary_location(ibw)
                call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                        BC_LIST_Y(1,IBW)),
     .           'BCENE = 16 only applicable to core boundaries !')
                call xertst(enepar(ibw,1).ne.0.0_R8,
     .           'faulty parameter enepar1')
              endif
            END DO
          endif
          IF(BCCHAR(IB).EQ.'S') THEN
            if(.not.once_bcene) then
              if (switch%use_astra.ne.0) then
#ifdef ASTRA
                SENEPAR = fheb_astra
#else
                call xerrab(
     .           'Compile with -DASTRA option to couple to ASTRA!')
#endif
              else
                SENEPAR = 0.0_R8
                DO IBW=IB,NBC
                  if(BCCHAR(IBW).EQ.'S'.and.BCENE(IBW).EQ.BCENE(IB))THEN
                    SENEPAR = SENEPAR + ENEPAR(IBW,1)
                  endif
                END DO
              endif
              write(*,*) 'BCENE=16: Summed value is ', SENEPAR
              TES = 0.0_R8
              volsum = 0.0_R8
              DO IBW=IB,NBC
                if(BCCHAR(IBW).EQ.'S'.and.BCENE(IBW).EQ.BCENE(IB))THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    TES = TES +
     &                  TE(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                 VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    VOLSUM = VOLSUM +
     &                 VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  END DO
                endif
              END DO
              TES=TES/VOLSUM
              totflux = 0.0_R8
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCENE(IBW).EQ.BCENE(IB)) THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    if(switch%mdf_fhe.ne.0) then
                      totflux = totflux +
     &           fhe_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1)
     &                +fac_ExB(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*   !srv 11.04.13 {
     &                    ((po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))-
     &               po(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))/
     &                     (hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &               hx(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))+
     &             (po(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))-
     &                      po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))/
     &             (hx(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))+
     &                      hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))*
     &                      ne(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                      te(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                      SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                      bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),2)/
     &                      bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),3)**2 !srv 11.04.13 }
                    else
                      totflux = totflux +
     .             fhe(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1)
                    endif
                  ENDDO
                ENDIF
              ENDDO
              TES = TES *
     &         (1.0_R8-switch%b2stbc_bc_ref_te*(totflux-SENEPAR)/                     !srv 11.04.13
     &           max(abs(totflux),SENEPAR))
              write(*,'(3(a,es14.7))')
     &         ' TES=',TES/ev,' totflux=',totflux,
     &         ' SENEPAR=',SENEPAR
              DO IBW=IB,NBC
                IF(BCCHAR(IBW).EQ.'S'.AND.BCENE(IBW).EQ.BCENE(IB)) THEN
                  DO IX = 1, BC_LIST_SIZE(IBW)
                    she0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)=
     .               VBIG*NE00*TES*
     .               SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    she0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)=
     .              -VBIG*NE00*
     .               SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  ENDDO
                ENDIF
              END DO
              once_bcene = .true.
            endif                                                      !srv 16.02.12 }
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCENE=16, ONLY AVAILABLE FOR CORE!')
          ENDIF
#endif
        else
          CALL XERRAB
     1       (' Incorrect style of boundary condition')
        endif

* -- BCENE=17 -- Feedback boundary condition with constant shared temperature

      case (17)
        call xerrab ('b2stbc_phys -- BCENE = 17 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'BCENE = 17 : shared temperature feedback by flux ',
     1     ENEPAR(IB,1),' with ',ENEPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(enepar(ib,1).ne.0, 'faulty parameter enepar1')
          call xertst(enepar(ib,2).gt.0, 'faulty parameter enepar2')
          call xertst(bcene(ib).eq.bceni(ib),
     1     'both bcene and bceni should be set to 17')
          call xertst(enipar(ib,1).ne.0, 'faulty parameter enipar1')
          call xertst(enipar(ib,2).gt.0, 'faulty parameter enipar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          totflux = 0.0_R8
          t0 = 0.0_R8
          weight = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type16_ix=bc_list_x(ix,ib)
            bc_type16_iy=bc_list_y(ix,ib)
            do bc_type16_i=1,bc_type16_ref    ! how many steps away from the boundary?
              bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
              bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
              bc_type16_ix=bc_type16_ixx
              call xertst (-2.lt.bc_type16_ix.and.bc_type16_ix.lt.nx+1,
     1         'stepped off x-grid in type 16 bc')
              call xertst (-2.lt.bc_type16_iy.and.bc_type16_iy.lt.ny+1,
     1         'stepped off y-grid in type 16 bc')
            enddo
            totflux=totflux+fhe(bc_type16_ix,bc_type16_iy,1,1)+
     .       fhi(bc_type16_ix,bc_type16_iy,1,1)+
     .       fhi_ext(bc_type16_ix,bc_type16_iy,1,1)
            if (bc_type16_kinetic_energy.eq.1) then
              do is = 0, ns-1
                totflux=totflux+(1.0_R8-switch%BoRiS)*
     .           fna(bc_type16_ix,bc_type16_iy,1,1,is)*0.5_R8*
     .           (kinrgy(bc_type16_ix,bc_type16_iy,is)+
     .            kinrgy(bottomix(bc_type16_ix,bc_type16_iy),
     .                   bottomiy(bc_type16_ix,bc_type16_iy),is))
              enddo
              do is = 0, ns_ext-1
                totflux=totflux+(1.0_R8-switch%BoRiS)*
     .           fa_ext(bc_type16_ix,bc_type16_iy,1,1,is)*0.5_R8*
     .           0.5_R8*am_ext(is)*mp*
     .           (ua_ext(bc_type16_ix,bc_type16_iy,is)**2+
     .            ua_ext(bottomix(bc_type16_ix,bc_type16_iy),
     .                   bottomiy(bc_type16_ix,bc_type16_iy),is)**2)
              enddo
            endif
            t0 = t0 +
     .        te(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .        ne(bc_list_x(ix,ib),bc_list_y(ix,ib)) +
     .        ti(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .        ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)
            weight = weight +
     1       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .        ne(bc_list_x(ix,ib),bc_list_y(ix,ib)) +
     2       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .        ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)
          ENDDO
          t0 = t0 / weight
          if(switch%b2stbc_diagno.ge.2) write(*,*) 'Totflux :',totflux
          T1 = (1.0_R8 -
     .     enepar(ib,2)*
     1       (totflux*enepar(ib,1)/(enepar(ib,1)+enipar(ib,1))-
     2         enepar(ib,1))/enepar(ib,1) -
     .     enipar(ib,2)*
     1       (totflux*enipar(ib,1)/(enepar(ib,1)+enipar(ib,1))-
     2         enipar(ib,1))/enipar(ib,1)
     .     ) * t0
          if(switch%b2stbc_diagno.ge.2)
     1     write(*,*) 'Old/New shared temp :',t0/qe, t1/qe
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=17, NOT AVAILABLE FOR N,W,E')
        ENDIF
#endif
* -- BCENE=18 -- FRACTIONAL DROP CONDITION (FOR EDGE2D COMPARISON)

      case (18)
        call xerrab ('b2stbc_phys -- BCENE = 18 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENE = 18 : electron temperature fractional drop ',
     1     ENEPAR(IB,1),' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENEPAR(IB,1).lt.1.0_R8,
     1     'BCENE = 18, ENEPAR(,1) => 1')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENEPAR(IB,1)) *
     .       TE(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENEPAR(IB,1)) *
     .       TE(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*T1*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENEPAR(IB,1)) *
     .       TE(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .        VBIG*NE00*T1*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .       -VBIG*NE00*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENEPAR(IB,1)) *
     .       TE(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= VBIG*NE00*T1*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            she0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-VBIG*NE00*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=18, INVALID OPTION')
        ENDIF
#endif
* -- BCENE=19 -- PRESCRIBE THE DECAY LENGTH FOR THE ELECTRON TEMPERATURE

      case (19)
        call xerrab ('b2stbc_phys -- BCENE = 19 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENE = 19 : electron temperature decay length (new) ',
     .     ENEPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(ENEPAR(IB,1).ne.0.0_R8,
     1     'BCENE = 19, ENEPAR(,1) = 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -hce0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/ENEPAR(IB,1)*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            t0=0.0_R8
            DO IS=0,NS-1
               t0=t0-bcint(rza(-1,-1,is),nx,ny,IX,IB)*(
     .              1.5_R8*fna_32(topix(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),
     .                            topiy(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),1,1,is)+
     .              2.5_R8*fna_52(topix(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),
     .                            topiy(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),1,1,is))
            ENDDO
            DO IS=0,NS_EXT-1
               t0=t0-bcint(za_ext(-1,-1,is),nx,ny,IX,IB)*
     .              2.5_R8*fa_ext(topix(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),
     .                            topiy(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),1,1,is)
            ENDDO
            if(t0.gt.0.0_R8) then
               she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .              t0 * TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            else
               she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .         she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) + t0
            endif
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
             she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
             she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -hce0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/ENEPAR(IB,1)*
     1            SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .               BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
             t0=0.0_R8
             DO IS=0,NS-1
                t0=t0+bcint(rza(-1,-1,is),nx,ny,IX,IB)*(
     .               1.5_R8*fna_32(bc_list_x(ix,ib),
     .                             bc_list_y(ix,ib),1,1,is)+
     .               2.5_R8*fna_52(bc_list_x(ix,ib),
     .                             bc_list_y(ix,ib),1,1,is))
             ENDDO
             DO IS=0,NS_EXT-1
                t0=t0+bcint(za_ext(-1,-1,is),nx,ny,IX,IB)*
     .               2.5_R8*fa_ext(bc_list_x(ix,ib),
     .                             bc_list_y(ix,ib),1,1,is)
             ENDDO
             if(t0.gt.0.0_R8) then
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .          she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .              t0 * TE(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
             else
                she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .          she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) + t0
             endif
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=19, NOT AVAILABLE FOR W,E')
        ENDIF
#endif
* -- BCENE=20 -- Feedback boundary condition with constant temperature

      case (20)
        call xerrab ('b2stbc_phys -- BCENE = 20 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'BCENE = 20 : electron temperature feedback by flux ',
     1     ENEPAR(IB,1),' with ',ENEPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(enepar(ib,1).ne.0, 'faulty parameter enepar1')
          call xertst(enepar(ib,2).gt.0, 'faulty parameter enepar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          totflux = 0.0_R8
          t0 = 0.0_R8
          weight = 0.0_R8
          teave = 0.0_R8
          neave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type16_ix=bc_list_x(ix,ib)
            bc_type16_iy=bc_list_y(ix,ib)
            do bc_type16_i=1,bc_type16_ref    ! how many steps away from the boundary?
              bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
              bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
              bc_type16_ix=bc_type16_ixx
              call xertst (-2.lt.bc_type16_ix.and.bc_type16_ix.lt.nx+1,
     1         'stepped off x-grid in type 16 bc')
              call xertst (-2.lt.bc_type16_iy.and.bc_type16_iy.lt.ny+1,
     1         'stepped off y-grid in type 16 bc')
            enddo
            totflux=totflux+
     .       fhe(bc_type16_ix,bc_type16_iy,1,1)
            t0 = t0 +
     .       te(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
            weight = weight + vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .                        ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
            bc_type20_ix=bc_list_x(ix,ib)
            bc_type20_iy=bc_list_y(ix,ib)
            do bc_type20_i=1,bc_type20_ref    ! how many steps away from the boundary?
              bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
              bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
              bc_type20_ix=bc_type20_ixx
              call xertst (-2.lt.bc_type20_ix.and.bc_type20_ix.lt.nx+1,
     1         'stepped off x-grid in type 20 bc')
              call xertst (-2.lt.bc_type20_iy.and.bc_type20_iy.lt.ny+1,
     1         'stepped off y-grid in type 20 bc')
            enddo
            teave=teave+te(bc_type20_ix,bc_type20_iy)*
     1        ne(bc_type20_ix,bc_type20_iy)*
     2       vol(bc_type20_ix,bc_type20_iy)
            neave=neave+ne(bc_type20_ix,bc_type20_iy)*
     1       vol(bc_type20_ix,bc_type20_iy)
          ENDDO
          t0 = t0 / weight
          teave=teave/neave
          if(switch%b2stbc_diagno.ge.2) write(*,*) 'Totflux :',totflux
          T1 = (1.0_R8 -
     .     enepar(ib,2)*(totflux-enepar(ib,1))/enepar(ib,1)) * t0
          if(switch%b2stbc_diagno.ge.2)
     1     write(*,*) 'Old/New electron temp :',t0/qe, t1/qe
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type20_ix=bc_list_x(ix,ib)
            bc_type20_iy=bc_list_y(ix,ib)
            do bc_type20_i=1,bc_type20_ref    ! how many steps away from the boundary?
              bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
              bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
              bc_type20_ix=bc_type20_ixx
              call xertst (-2.lt.bc_type20_ix.and.bc_type20_ix.lt.nx+1,
     1         'stepped off x-grid in type 20 bc')
              call xertst (-2.lt.bc_type20_iy.and.bc_type20_iy.lt.ny+1,
     1         'stepped off y-grid in type 20 bc')
            enddo
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*T1*
     1       te(bc_type20_ix,bc_type20_iy)/teave
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .       -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=20, NOT AVAILABLE FOR N,W,E')
        ENDIF
#endif
* -- BCENE=21 -- constant temperature scaled by temperature on the ring bc_type21_ref away

      case (21)
        call xerrab ('b2stbc_phys -- BCENE = 21 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'BCENE = 21 : fixed electron temperature '//
     .     '(modified by nearby ring) ',
     1     ENEPAR(IB,1),' with ',ENEPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(enepar(ib,1).ne.0.0_R8,
     1     'faulty parameter enepar1')
          call xertst(enepar(ib,2).gt.0.0_R8,
     1     'faulty parameter enepar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          teave = 0.0_R8
          neave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref    ! how many steps away from the boundary?
              bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            teave=teave+te(bc_type21_ix,bc_type21_iy)*
     1       ne(bc_type21_ix,bc_type21_iy)*
     2       vol(bc_type21_ix,bc_type21_iy)
            neave=neave+ne(bc_type21_ix,bc_type21_iy)*
     1       vol(bc_type21_ix,bc_type21_iy)
          ENDDO
          teave=teave/neave
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref    ! how many steps away from the boundary?
              bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*ENEPAR(IB,1)*EV*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .       ((1.0_R8-enepar(ib,2))+
     .           enepar(ib,2)*te(bc_type21_ix,bc_type21_iy)/teave)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE IF(BCCHAR(IB).EQ.'N') THEN
          teave = 0.0_R8
          neave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref    ! how many steps away from the boundary?
              bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            teave=teave+te(bc_type21_ix,bc_type21_iy)*
     1       ne(bc_type21_ix,bc_type21_iy)*
     2       vol(bc_type21_ix,bc_type21_iy)
            neave=neave+ne(bc_type21_ix,bc_type21_iy)*
     1       vol(bc_type21_ix,bc_type21_iy)
          ENDDO
          teave=teave/neave
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref    ! how many steps away from the boundary?
              bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*ENEPAR(IB,1)*EV*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .       ((1.0_R8-enepar(ib,2))+
     .           enepar(ib,2)*te(bc_type21_ix,bc_type21_iy)/teave)
            she0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=21, NOT AVAILABLE FOR W,E')
        ENDIF
#endif
* -- BCENE=22 -- RADIAL LEAKAGE CONDITION FOR THE ELECTRON TEMPERATURE   !srv 23.07.10 {
*                SPECIFIES A TEMPERATURE GRADIENT SUCH THAT THE DIFFUSIVE
*                FLUX IS SET TO BE EQUIVALENT TO A LEAKAGE VELOCITY, IN
*                UNITS OF THE ELECTRON THERMAL VELOCITY

      case (22)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENE = 22 : with ENEPAR(IB,1)= ',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
        endif
cwdk    Notes:
cwdk     * need correction for misaligned grids -> t1
cwdk     * need to correct this BC for hybrid scheme
        do ibc = 1, mpg%bcCvP(ib,2)                                    ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                     ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                     ! number of the domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
          vte  = sqrt(pl%te(iCv1)/me)
          t0   = geo%fcQalf(iFc,1)
          t1   = (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     &            geo%fcQbet(iFc,0)
          srw%she0(iCv1,0) =  VBIG*NE00*EV*
     &        co%chce(iFc,1)*(t0*pl%te(iCv2) -
     &        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*t1*
     &       (pl%te(mpg%fcVx(iFc,2)) - pl%te(mpg%fcVx(iFc,1))))
          srw%she0(iCv1,1) = -VBIG*NE00*EV*
     &       (co%chce(iFc,1)*t0 -
     &       ENEPAR(iB,1)*geo%fcS(iFc)*vte*dv%ne(iCv1))
        enddo


* -- BCENE=23 -- Prescribe the value of the electron temperature, with poloidal
*                perturbation taken from flux tube just inside domain

      case (23)
        if(ncall_b2stbc_phys.eq.0) then
          do ibc = 1, mpg%bcCvP(ib,2)
           iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
           call xertst(mpg%cvOnClosedSurface(iCv1),
     .      'BCENE=23 intended for core boundaries only!')
          enddo
          IF(switch%use_astra.ne.0) THEN
          ENEPAR(IB,1) = teb_astra*1.0e3_R8
          ENDIF
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified electron temp. + perturbation',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENEPAR(IB,1).gt.0.0_R8,
     1     'BCENE = 23, ENEPAR(,1) !> 0')
        endif
*     ..compute average temperature on core boundary
        t0 = 0.0e0_R8
        t1 = 0.0e0_R8
        do ib0 = 1, nbc
         ! identify core boundaries; assume all core boundaries
         ! form single flux surface
         ! implicit assumption at the moment: first ring in the
         ! core is a 'structured grid' ring with quad cells
         iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
         if (mpg%cvOnClosedSurface(iCv1)) then
          do ibc0 = 1, mpg%bcCvP(ib0,2)                                 ! loop over number of cells in the boundary
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                   ! number of the guard cell
           iCv2 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,2)                   ! number of the domain cell
           iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
           t0 = t0 + geo%cvVol(iCv2)
           t1 = t1 + geo%cvVol(iCv2)*pl%te(iCv2)
          enddo
         endif
        enddo
        t1 = t1/t0
*     ..impose temperature with perturbation
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                      ! number of the domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          srw%she0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*
     .      (pl%te(iCv2) + enepar(iB,1)*EV - t1)
          srw%she0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
        enddo


* -- BCENE=24 -- Feedback on the total electron energy flux, by imposing a flux surface
*                averaged temperature, with poloidal perturbation taken from flux
*                tube just inside domain

      case (24)
        if(ncall_b2stbc_phys.eq.0) then
          do ibc = 1, mpg%bcCvP(ib,2)
           iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
           call xertst(mpg%cvOnClosedSurface(iCv1),
     .      'BCENE=24 intended for core boundaries only!')
          enddo
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified electron heat flux',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENEPAR(IB,1).gt.0.0_R8,
     1     'BCENE = 24, ENEPAR(,1) !> 0')
        endif
*     ..compute average temperature on core boundary
        t0 = 0.0e0_R8
        t1 = 0.0e0_R8
        t2 = 0.0e0_R8
        totflux = 0.0e0_R8
        totpar  = 0.0e0_R8
        do ib0 = 1, nbc
         ! identify core boundaries; assume all core boundaries
         ! form single flux surface
         ! implicit assumption at the moment: first ring in the
         ! core is a 'structured grid' ring with quad cells
         iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
         if (mpg%cvOnClosedSurface(iCv1)) then
          if (bcene(ib0).ne.24)
     .      call xerrab ('BCENE=24: incompatible core BCs')
          do ibc0 = 1, mpg%bcCvP(ib0,2)                                 ! loop over number of cells in the boundary
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                   ! number of the guard cell
           iCv2 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,2)                   ! number of the domain cell
           iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
           t0 = t0 + geo%cvVol(iCv2)
           t1 = t1 + geo%cvVol(iCv2)*pl%te(iCv1)
           t2 = t2 + geo%cvVol(iCv2)*pl%te(iCv2)
           totflux = totflux - (dv%fhe(iFc,0) + dv%fhe(iFc,1))*
     &       mpg%bcFcOr(mpg%bcCvP(ib0,1)+ibc0-1)
          enddo
          totpar = totpar + enepar(ib0,1)
         endif
        enddo
        t1 = t1/t0
        t2 = t2/t0
        if (enepar(ib,2).gt.0.0_R8) t3 = enepar(ib,2)
        if (enepar(ib,2).le.0.0_R8) t3 = 1.0e-2_R8
        t4 = (1.0_R8-t3*(totflux-totpar)/
     &         max(abs(totflux),abs(totpar)))*t1
*     ..impose temperature with perturbation
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                      ! number of the domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          srw%she0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*
     .      (pl%te(iCv2) + t4 - t2)
          srw%she0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
        enddo

      case default
        write(*,*) BCENE(IB)
        call xerrab ('BCENE out of bounds')

      end select

* -- BOUNDARY CONDITIONS FOR ION ENERGY EQUATIONS --

      select case (BCENI(IB))

      case (0)
        continue

* -- BCENI=1 -- PRESCRIBE THE VALUE OF THE ION TEMPERATURE

      case (1)
        if(ncall_b2stbc_phys.eq.0) then
          IF(BCCHAR(IB).EQ.'S'.and.switch%use_astra.ne.0) THEN
            ENIPAR(IB,1) = tib_astra*1.0e3_R8                          !srv 14.06.12
          ENDIF
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENI =  1 : specified ion temperature ',ENIPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENIPAR(IB,1).gt.0.0_R8,
     1     'BCENI = 1, ENIPAR(,1) !> 0')
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          srw%shi0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*enipar(iB,1)*EV
          srw%shi0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
        enddo

* -- BCENI=2 -- PRESCRIBE THE GRADIENT OF THE ION TEMPERATURE

      case (2)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENI =  2 : ion temperature gradient ',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                    ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
          iVx1 = mpg%fcVx(iFc,1)
          iVx2 = mpg%fcVx(iFc,2)
          t1 = intvertex_s(iVx1,nCv,nVx,mpg,geo%vxVol,pl%ti)
          t2 = intvertex_s(iVx2,nCv,nVx,mpg,geo%vxVol,pl%ti)
          srw%shi0(iCv1,0)=  VBIG*NE00*geo%fcS(iFc)*
     .       (EV*enipar(ib,1)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2))*
     .        geo%fcQgam(iFc,0) - cor9*(t2 - t1)*geo%fcQgam(iFc,1)*
     .        (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     .        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1) + pl%ti(iCv2))
          srw%shi0(iCv1,1)= -VBIG*NE00*geo%fcS(iFc)
        enddo

* -- BCENI=3 -- SHEATH CONDITIONS, ION ENERGY TRANSMISSION COEFFICIENT

      case (3)
        call xerrab ('b2stbc_phys -- BCENI = 3 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,g14.7,a,g14.7,a,a,a)')
     .     'BCENI =  3 : sheath ion with MOMPAR(ismain,IB,2) =',
     1     MOMPAR(ismain,IB,2),' and MOMPAR(ismain,IB,1) =',
     2     MOMPAR(ismain,IB,1),
     .     ' and ENIPAR(IB,1) =',enipar(ib,1),               !srv 24.02.10 {
     1     ' and ENIPAR(IB,2) =',enipar(ib,2),
     2     ' on ',BCCHAR(IB),boundary_location(ib)           !srv 24.02.10 }
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
              CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))  !srv 01.02.09
              UUOUT=MOMPAR(ismain,IB,1)*CS*
     .         abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .             wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3))
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  if(enipar(ib,1).gt.0.0_R8) then                        !xpb
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) - UUOUT *
     .               ENIPAR(IB,1)*bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .               SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
                  else     ! xpb new pressure linearisation as per B2.5-style
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)=
     .              SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3) - UUOUT *  !xpb
     .               (-ENIPAR(IB,1))*bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .                SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     .                ni(bc_list_x(iy,ib),bc_list_y(iy,ib),0)
                  endif                                                  !xpb
                endif
              end do
            ELSE
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  ZBAR=bcint(na(-1,-1,is),nx,ny,IY,IB)/
     .             bcint(NE,nx,ny,IY,IB)
                  ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .             bcint(rz2(-1,-1,is),nx,ny,IY,IB)*ZBAR*
     .             bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                  CS=SQRT(ZARG)
                  UUTEST=abs(
     .             UU(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .             (1.0_R8+
     .              HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .      HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .         RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .      UU(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .         RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .              HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .      HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .         RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))
                  UUOUT=MAX(MOMPAR(IS,IB,1)*CS*
     .             abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .                 wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .             UUTEST)
                  if(enipar(ib,1).gt.0.0_R8) then                        !xpb
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) - UUOUT *
     .               ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .               SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
                  else     ! xpb new pressure linearisation as per B2.5-style
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)=
     .              SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3) - UUOUT *  !xpb
     .               (-ENIPAR(IB,1))*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .                SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))/
     .                ni(bc_list_x(iy,ib),bc_list_y(iy,ib),0)
                  endif                                                  !xpb
                endif
              end do
            ENDIF
            if(switch%mdf_fhi.ne.0) then                                        !srv 12.09.08 {
              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &         shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) -
     &         fhiPSch(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                 rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,0)
            endif                                                        !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
              CS=SQRT(bcint(PZ,nx,ny,IY,IB)/bcint(RZ,nx,ny,IY,IB))       !srv 01.02.09
              UUOUT=MOMPAR(ismain,IB,1)*CS*
     .         abs(wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .             wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),3))
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  if(enipar(ib,1).gt.0.0_R8) then                        !xpb
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) - UUOUT *
     .               ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .               SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                  LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
                  else     ! xpb new pressure linearisation as per B2.5-style
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)=
     .              SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3) - UUOUT *  !xpb
     .               (-ENIPAR(IB,1))*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .                SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .                ni(bc_list_x(iy,ib),bc_list_y(iy,ib),0)
                  endif                                                  !xpb
                endif
              end do
            ELSE
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  ZBAR=bcint(na(-1,-1,IS),nx,ny,IY,IB)/
     .             bcint(NE,nx,ny,IY,IB)
                  ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .             bcint(rz2(-1,-1,IS),nx,ny,IY,IB)*ZBAR*
     .             bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                  CS=SQRT(ZARG)
                  UUTEST=abs(
     .             UU(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     1             (1.0_R8+
     2              HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .       HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .       UU(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .              HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .       HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .          LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                 LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))
                  UUOUT=MAX(MOMPAR(IS,IB,1)*CS*
     .             abs(wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .                 wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .             UUTEST)
                  if(enipar(ib,1).gt.0.0_R8) then                        !xpb
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) - UUOUT *
     .               ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .               SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
                  else     ! xpb new pressure linearisation as per B2.5-style
                    shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3)=
     .              SHI0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),3) - UUOUT *  !xpb
     .               (-ENIPAR(IB,1))*
     .                 bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .                 SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                    LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .                ni(bc_list_x(iy,ib),bc_list_y(iy,ib),0)
                  endif                                                  !xpb
                endif
              end do
            ENDIF
            if(switch%mdf_fhi.ne.0) then                                 !srv 12.09.08 {
              shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) =
     &         shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0) +
     &         fhiPSch(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,0)
            endif                                                        !srv 12.09.08 }
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=3, NO SHEATH IN RADIAL DIR.')
        END IF
#endif
* -- BCENI=4 -- PRESCRIBE THE VALUE OF THE ION TEMPERATURE, WEAKLY
*                 A MIXED BOUNDARY CONDITION

      case (4)
        call xerrab ('b2stbc_phys -- BCENI = 4 not adapted for WG')
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENI =  4 : weak ion temperature ',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(ENIPAR(IB,1).gt.0.0_R8,
     1     'BCENI = 4, ENIPAR(,1) !> 0')
          call xertst(ENIPAR(IB,2).gt.0.0_R8,                !srv 02.02.11
     1     'BCENI = 4, ENIPAR(,2) !> 0')
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          srw%shi0(iCv1,0) =  ENIPAR(IB,2)*geo%fcS(iFc)*
     .                        ENIPAR(IB,1)
          srw%shi0(iCv1,1) = -ENIPAR(IB,2)/EV*geo%fcS(iFc)
        enddo

* -- BCENI=5 -- PRESCRIBE THE ION ENERGY FLUX PER UNIT AREA

      case (5)
        call xerrab ('b2stbc_phys -- BCENI = 5 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENI =  5 : ion energy flux per area ',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENIPAR(IB,1)*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENIPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=ENIPAR(IB,1)*
     .        SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=ENIPAR(IB,1)*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=5, NOT YET AVAILABLE')
        END IF
#endif
* -- BCENI=6 -- PRESCRIBE THE TOTAL ION ENERGY FLUX FOR A
*                 CONSTANT ION TEMPERATURE

      case (6)
        call xerrab ('b2stbc_phys -- BCENI = 6 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENI =  6 : total ion energy flux ',ENIPAR(IB,1),
     1     ' with const Ti on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'N') THEN
          Tguess=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            Tguess=Tguess+TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
          Tguess=Tguess/BC_LIST_SIZE(IB)
          if(switch%b2stbc_diagno.ge.2)
     >     write(*,*) 'Initial Tguess:',Tguess/ev
          T0=Tguess*0.95_R8
          T1=Tguess*1.01_R8
          Tlimit=Tguess*0.0001_R8
          fylimt=ENIPAR(IB,1)*0.0001_R8
          ibnd=-1
30610     ibnd=ibnd+1
          if (ibnd.eq.0) then
            Tguess=T0
          elseif (ibnd.eq.1) then
            Tguess=T1
          else
            Tguess=T1-GAM1*(T1-T0)/(EPSLON + GAM1 - GAM0)
          endif
          fyt=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
cxpb code segment analog to b2tfhi
            if (conduction_only.eq.1) then
              flo = 0.0_R8
            else
              cum_flo_vadia=0.0_R8
              cum_flo_vaecrb=0.0_R8
              do is=0,ns-1
                prefactor=
     &           (vol(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &            vol(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &           /(hy1(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &             hy1(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &           *(na(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),is)
     &            +na(bc_list_x(ix,ib),bc_list_y(ix,ib),is))/2.0_R8
                cum_flo_vadia=cum_flo_vadia+
     1         prefactor*vadia(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1,is)
                cum_flo_vaecrb=cum_flo_vaecrb+prefactor*
     1           vaecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1,is)
              enddo
              flo = (1.5_R8+switch%BoRiS)*
     .         fni_32(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)+
     1         2.5_R8*fni_52(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)+
     &         chvi(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)*
     &         (ni(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),0)+
     &          ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0))/2.0_R8
              if (switch%mdf_fhi.eq.0) then
                flo = flo + 2.5_R8*cum_flo_vadia+1.5_R8*cum_flo_vaecrb
              else
                flo = flo + 1.5_R8*cum_flo_vaecrb
              endif
            endif
            w1 = 0.0_R8
            if (switch%BoRiS.eq.1.0_R8) then
              do is = 0, ns-1
                t2 = dim ( coni0,
     .           abs(kinrgy(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .                  fna(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1,is))*
     .           max(0.0_R8,
     .            min(ti(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))))/
     .           max(cutlo,ti(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .                ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))))
                w1 = w1 +
     .            upwind(kinrgy(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .                   fna(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1,is),
     .             ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     .             ti(bc_list_x(ix,ib),bc_list_y(ix,ib))) - t2 *
     .            (ti(bc_list_x(ix,ib),bc_list_y(ix,ib)) -
     .             ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
              end do
            end if
            if (no_hybr_i.eq.0) then
              coni0 = hybr(flo,
     .         chci(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1))
              fy = w1 + flo*
     .         (ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &         +Tguess)/2.0_R8-coni0*
     1         (Tguess-ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     2                    bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
            else
              coni0 = chci(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1)
              t2 = dim ( coni0,
     .         abs(flo) *
     .         max(0.0_R8,
     .             min(Tguess,ti(bottomix(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)),
     .                           bottomiy(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)))))/
     .         max(cutlo,
     .                 Tguess+ti(bottomix(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)),
     .                           bottomiy(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)))) )
              fy = w1 +
     .         upwind(flo,
     .                ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     .                Tguess)
     .         - t2 * (Tguess-ti(bottomix(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib)),
     .                           bottomiy(bc_list_x(ix,ib),
     .                                    bc_list_y(ix,ib))))
            endif
            fyt=fyt + fy
          ENDDO
          if(switch%b2stbc_diagno.ge.2)
     >     write(*,'(a,i4,1pe15.6,1pe14.6)') 'ibnd,Tg,fyt(i):',
     .     ibnd,Tguess/ev,fyt
          if (ibnd.eq.0) then
            GAM0 = fyt - ENIPAR(IB,1)
            goto 30610
          elseif (ibnd.eq.1) then
            GAM1 = fyt - ENIPAR(IB,1)
            goto 30610
          else
            T0 = T1
            GAM0 = GAM1 - ENIPAR(IB,1)
            T1 = Tguess
            GAM1 = fyt
            if (ibnd .gt. 100)
     >       call xerrab ('B2STBC -- BCENI=6, NOT CONVERGING')
            if (abs(t1-t0).gt.Tlimit .and.
     &          abs(fyt-ENIPAR(IB,1)).gt.fylimt) goto 30610
          endif
          if(switch%b2stbc_diagno.ge.1)
     >     write(*,'(a,f15.6,i4,1pe14.6)') ' tibound:',T1/ev,ibnd,fyt
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*T1*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'S') THEN
          Tguess=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            Tguess=Tguess+TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
          Tguess=Tguess/BC_LIST_SIZE(IB)
          if(switch%b2stbc_diagno.ge.2)
     >     write(*,*) 'Initial Tguess:',Tguess/ev
          T0=Tguess*0.95_R8
          T1=Tguess*1.01_R8
          Tlimit=Tguess*0.0001_R8
          fylimt=ENIPAR(IB,1)*0.0001_R8
          ibnd=-1
30620     ibnd=ibnd+1
          if (ibnd.eq.0) then
            Tguess=T0
          elseif (ibnd.eq.1) then
            Tguess=T1
          else
            Tguess=T1-GAM1*(T1-T0)/(EPSLON + GAM1 - GAM0)
          endif
          fyt=0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
cxpb code segment analog to b2tfhi
            if (conduction_only.eq.1) then
              flo = 0.0_R8
            else
             cum_flo_vadia = 0.0_R8
             cum_flo_vaecrb = 0.0_R8
             do is=0,ns-1
               prefactor=
     &          (vol(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &           vol(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &          /(hy1(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     &            hy1(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &          *(na(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),is)
     &           +na(bc_list_x(ix,ib),bc_list_y(ix,ib),is))/2.0_R8
               cum_flo_vadia=cum_flo_vadia+prefactor*
     1          vadia(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     2                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1,is)
               cum_flo_vaecrb=cum_flo_vaecrb+
     1        prefactor*vaecrb(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1,is)
             enddo
             flo = (1.5_R8+switch%BoRiS)*
     1         fni_32(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1)+
     1         2.5_R8*fni_52(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1)+
     &         chvi(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1)*
     &         (ni(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),0)+
     &          ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0))/2.0_R8
              if (switch%mdf_fhi.eq.0) then
                flo = flo + 2.5_R8*cum_flo_vadia+1.5_R8*cum_flo_vaecrb
              else
                flo = flo + 1.5_R8*cum_flo_vaecrb
              endif
            endif
            w1 = 0.0_R8
            if (switch%BoRiS.eq.1.0_R8) then
              do is = 0, ns-1
                t2 = dim ( coni0,
     .           abs(kinrgy(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .            fna(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1,is))
     .           *max(0.0_R8,
     .            min(ti(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))))/
     .           max(cutlo,ti(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .                ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))))
                w1 = w1 +
     .            upwind(kinrgy(bc_list_x(ix,ib),bc_list_y(ix,ib),is)*
     .                fna(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                    topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                    1,1,is),
     .             ti(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .             - t2 *
     .            (ti(bc_list_x(ix,ib),bc_list_y(ix,ib)) -
     .             ti(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))))
              end do
            end if
            if (no_hybr_i.eq.0) then
              coni0 = hybr(flo,
     &         chci(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1))
              fy = w1 +
     .          flo*(ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                     topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     &                 +Tguess)/2.0_R8-coni0*
     1         (ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-Tguess)
            else
              coni0 = chci(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                     topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1)
              t2 = dim ( coni0,
     .         abs(flo) *
     .         max(0.0_R8,
     .             min(Tguess,ti(topix(bc_list_x(ix,ib),
     .                                 bc_list_y(ix,ib)),
     .                           topiy(bc_list_x(ix,ib),
     .                                 bc_list_y(ix,ib)))))/
     .         max(cutlo,
     .                 Tguess+ti(topix(bc_list_x(ix,ib),
     .                                 bc_list_y(ix,ib)),
     .                           topiy(bc_list_x(ix,ib),
     .                                 bc_list_y(ix,ib)))) )
              fy = w1 +
     .         upwind(flo, Tguess,
     .                ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .         - t2 * (ti(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                    topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     .                 Tguess)
            endif
            fyt=fyt + fy
          ENDDO
          if(switch%b2stbc_diagno.ge.2)
     >     write(*,'(a,i4,f15.6,1pe14.6)') 'ibnd,Tg,fyt(i):',
     .     ibnd,Tguess/ev,fyt
          if (ibnd.eq.0) then
            GAM0 = fyt - ENIPAR(IB,1)
            goto 30620
          elseif (ibnd.eq.1) then
            GAM1 = fyt - ENIPAR(IB,1)
            goto 30620
          else
            T0 = T1
            GAM0 = GAM1 - ENIPAR(IB,1)
            T1 = Tguess
            GAM1 = fyt
            if (ibnd .gt. 100)
     >       call xerrab ('B2STBC -- BCENI=6, NOT CONVERGING')
            if (abs(t1-t0).gt.Tlimit .and.
     &          abs(fyt-ENIPAR(ib,1)).gt.fylimt) goto 30620
          endif
          if(switch%b2stbc_diagno.ge.1)
     >     write(*,'(a,f15.6,i4,1pe14.6)') ' tibound:',T1/ev,ibnd,fyt

          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .          VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .         -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=6, NOT AVAILABLE FOR W,E')
        ENDIF
#endif

* -- BCENI=7 -- PRESCRIBE THE ION TEMPERATURE PROFILE

      case (7)
        if(ncall_b2stbc_phys.eq.0) then                               !srv 26.02.18 {
          write(*,'(3a)') 'specified ion temperature on ',
     1     BCCHAR(IB), boundary_location(ib)
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          srw%shi0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*bv_ti(iCv1)*EV
          srw%shi0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
        enddo

* -- BCENI=8 -- PRESCRIBE THE TOTAL ION HEAT FLUX WITH CONSTANT
*                 FLUX DENSITY

      case (8)
        if(ncall_b2stbc_phys.eq.0) then
          if (switch%use_astra.ne.0) ENIPAR(IB,1) = fhib_astra           !srv 28.07.11
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENI =  8 : total ion energy flux',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        US = 0.0_R8
        do ibc = 1, mpg%bcCvP(ib,2)
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          US   = US + geo%fcS(iFc)
        enddo
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          srw%shi0(iCv1,0) = enipar(iB,1)/US*geo%fcS(iFc)
          srw%shi0(iCv1,1) = 0.0_R8
          if(switch%mdf_fhi.ne.0) then
            !wdk Correction now involves vedia instead of ExB velocity
            !wdk This makes implementation equivalent to mdf_fhi = 0 option
            srw%shi0(iCv1,0) = srw%shi0(iCv1,0) +
     &        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &       (dv%fhiPSch(iFc,0) + dv%fhiPSch(iFc,1))
            t0 = 0.0_R8
            do is=0,ns-1
              t0 = t0 + mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*2.5_R8*
     &          (geo%fcS(iFc)*geo%fcQalf(iFc,0)*dv%vadia(iFc,0,is) +
     &           geo%fcS(iFc)*geo%fcQalf(iFc,1)*dv%vadia(iFc,1,is))*
     &           0.5_R8*(pl%na(iCv1,is)+pl%na(iCv2,is))*
     &           0.5_R8*(pl%ti(iCv1)+pl%ti(iCv2))
            end do
            srw%shi0(iCv1,0) = srw%shi0(iCv1,0) +
     &        max(0.0_R8,t0)
            srw%shi0(iCv1,1) = srw%shi0(iCv1,1) +
     &        min(0.0_R8,t0)/pl%ti(iCv1)
          endif
        enddo

* -- BCENI=9 -- PRESCRIBE THE DECAY LENGTH FOR THE ION TEMPERATURE

      case (9)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENI =  9 : ion temperature decay length ',ENIPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENIPAR(IB,1).ne.0.0_R8,
     1     'BCENI = 9, ENIPAR(,1) = 0')
          decay_length_ok = .true.
          do is = 0, ns-1
            decay_length_ok = decay_length_ok .and.
     1       (bccon(is,ib).eq.9 .or. bccon(is,ib).eq.0 .or.
     2       (is_neutral(is) .and. switch%use_eirene.eq.1))
            if (bccon(is,ib).eq.1) then
              write(0,*)
     &         'Should use BCENI=19 in conjunction with BCCON=1'
            end if
          end do
          if (.not.decay_length_ok) then
            write(0,*)
     &       'Should use BCENI.EQ.9 together with BCCON=9'
            write(*,*)
     &       'Should use BCENI.EQ.9 together with BCCON=9'
            write(*,'(1x,a)')
     &                 'Please ensure that this boundary condition '//
     &                 'is consistent with '//
     &                 'the density boundary condition'  !srv 21.09.12
          endif
        endif
        do ibc = 1, mpg%bcCvP(ib,2)
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          sy = geo%fcS(iFc)
          srw%shi0(iCv1,0) = 0.0_R8
          if (switch%recycle_afn.eq.0 .or. switch%tn_style.ne.0) then
            !! Note that
            !! if tn_style=0: hci0 includes contributions from all ions and all fluid neutrals
            !! if tn_style=1: hci0 includes contributions from ions only
            !! if tn_style=2: hci0 includes contributions from all ions and non-hydrogenic neutrals
            srw%shi0(iCv1,1) = -co%hci0(iCv1)/enipar(ib,1)*sy
          else                                                       !! tn_style=0.and.recycle_afn=1 => take out hydrogenic atoms
            srw%shi0(iCv1,1) = -(co%hci0(iCv1)-co%hcn0(iCv1))/
     &        enipar(ib,1)*sy
          endif
          do is = 0, ns-1
            if(.not.(is_neutral(is)) .or. switch%tn_style.eq.0 .or.
     &         (switch%tn_style.eq.2.and.zn(is).ne.1)) then
              srw%shi0(iCv1,0) = srw%shi0(iCv1,0) +
     &                     2.5_R8*pl%ti(iCv1)*srw%sna0(iCv1,0,is)
              srw%shi0(iCv1,1) = srw%shi0(iCv1,1) +
     &                     2.5_R8*pl%na(iCv1,is)*srw%sna0(iCv1,1,is)
            endif
          enddo !! end species loop
          if(switch%mdf_fhi.ne.0) then
            srw%shi0(iCv1,0) = srw%shi0(iCv1,0) +
     &        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &       (dv%fhiPSch(iFc,0) + dv%fhiPSch(iFc,1))
            do is=0,ns-1
              t0 = -mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*1.5_R8*
     &          (geo%fcS(iFc)*geo%fcQalf(iFc,0)*dv%vaecrb(iFc,0,is) +
     &           geo%fcS(iFc)*geo%fcQalf(iFc,1)*dv%vaecrb(iFc,1,is))*
     &           0.5_R8*(pl%na(iCv1,is)+pl%na(iCv2,is))*
     &           0.5_R8*(pl%ti(iCv1)+pl%ti(iCv2))
              srw%shi0(iCv1,0) = srw%shi0(iCv1,0) +
     &          max(0.0_R8,t0)
              srw%shi0(iCv1,1) = srw%shi0(iCv1,1) +
     &          min(0.0_R8,t0)/pl%ti(iCv1)
            enddo
          endif
        enddo   !! end boundary cell loop
#ifdef WG_TODO
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .         -hci0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/ENIPAR(IB,1)*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            DO IS=0,NS-1
cxpb          if (.not.is_neutral(is)) then   ! dpc --- need to check
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .           2.5_R8*TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .           SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)    !srv 23.09.08
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     .           2.5_R8*
     .           SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*   !srv 23.09.08
     &                    na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
cxpb          endif
            ENDDO
            if(switch%mdf_fhi.ne.0) then                                       !srv 12.09.08 {
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &        shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &        fhiPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1)
              do is=0,ns-1
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     &           1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                  hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &          vaecrb(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1,is)
     &                   *na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
              end do
            endif                                                       !srv 12.09.08 }
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .         -hci0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/ENIPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            DO IS=0,NS-1
cxpb          if (.not.is_neutral(is)) then   ! dpc --- need to check
cxpb            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
cxpb     .       shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
cxpb     .           rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*2.5_R8*
cxpb     .            TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
cxpb     .           SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)
cxpb            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
cxpb     .       shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
cxpb     .           rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*2.5_R8*
cxpb     .       SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*
cxpb     .           na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
cxpb          endif
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .           2.5_R8*TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0,IS)      !srv 23.09.08
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) +
     .           2.5_R8*
     .         SNA0_NO_MDF(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,IS)*     !srv 23.09.08
     &                  na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),IS)
            ENDDO
            if(switch%mdf_fhi.ne.0) then                                       !srv 12.09.08 {
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     .         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &         fhiPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1)
              do is=0,ns-1
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) =
     &           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) -
     &           1.5_R8*vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &                  hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     &           vaecrb(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1,is)*
     &               na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
             end do
            endif                                                       !srv 12.09.08 }
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=9, NOT AVAILABLE FOR W,E')
        ENDIF

#endif
* -- BCENI=10 -- Feedback option for core

      case (10)
        call xerrab ('b2stbc_phys -- BCENI = 10 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENI = 10 : ion heat flux feedback bc on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(lfeedback, 'LFEEDBACK must be turned on!')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if (region(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0).eq.1) then
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     1         cbshi(0,ismain,coreregno)*
     .         SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            elseif (region(bc_list_x(ix,ib),bc_list_y(ix,ib),0).eq.5)
     .       then
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     1         cbshi(0,ismain,coreregn2)*
     .         SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            endif
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=10, NOT AVAILABLE FOR N,W,E')
        ENDIF
#endif
* -- BCENI=11 -- SHEATH CONDITIONS, ION ENERGY TRANSMISSION COEFFICIENT

      case (11)
        call xerrab ('b2stbc_phys -- BCENI = 11 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,a,a)')
     .     'BCENI = 11 : sheath ion bc on ',BCCHAR(IB),
     1     boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            DO IS=0,NS-1
              if (.not.is_neutral(is)) then
                ZBAR=bcint(na(-1,-1,IS),nx,ny,IY,IB)/
     .           bcint(NE,nx,ny,IY,IB)
                ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .           bcint(rz2(-1,-1,IS),nx,ny,IY,IB)*ZBAR*
     .           bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                CS=SQRT(ZARG)
                UUTEST=
     .           UU(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     .           (1.0_R8+
     .            HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .    HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .       RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .    UU(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .       RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .            HX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .    HX(RIGHTIX(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .       RIGHTIY(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .               RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
                UUOUT=MAX(MOMPAR(IS,IB,1)*CS*
     .           abs(wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .               wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .           -UUTEST)
                shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .           shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) - UUOUT *
     .           ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .           SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
              endif
            end do
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            DO IS=0,NS-1
              if (.not.is_neutral(is)) then
                ZBAR=bcint(na(-1,-1,IS),nx,ny,IY,IB)/
     .           bcint(NE,nx,ny,IY,IB)
                ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .           bcint(rz2(-1,-1,IS),nx,ny,IY,IB)*ZBAR*
     .           bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                CS=SQRT(ZARG)
                UUTEST=UU(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                    LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),IS)*
     1        (1.0_R8+HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .         HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .            LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .            LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))))-
     .         UU(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .            LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),IS)*
     .                HX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))/
     .         HX(LEFTIX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))),
     .            LEFTIY(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .                   LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))))
                UUOUT=MAX(MOMPAR(IS,IB,1)*CS*
     .           abs(wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0)/
     .               wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),3)),
     .           UUTEST)
                shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .           shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1) - UUOUT *
     .           ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .           SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .              LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
              endif
            end do
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=11, NO SHEATH IN RADIAL DIR.')
        END IF
#endif
* -- BCENI=12 -- SHEATH CONDITIONS, ION ENERGY TRANSMISSION COEFFICIENT

      case (12)
        call xerrab ('b2stbc_phys -- BCENI = 12 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,a,a)')
     .     'BCENI = 12 : sheath ion bc on ',BCCHAR(IB),
     1     boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            DO IS=0,NS-1
              shi0(bc_list_x(ix,ib),bc_list_y(ix,ib),3) =
     .        shi0(bc_list_x(ix,ib),bc_list_y(ix,ib),3)-
     1         ENIPAR(IB,1)*max(0.0_R8,
     1          -fna(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1,1,is))/
     2         ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)
            end do
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            DO IS=0,NS-1
              shi0(bc_list_x(ix,ib),bc_list_y(ix,ib),3) =
     .        shi0(bc_list_x(ix,ib),bc_list_y(ix,ib),3)-
     1         ENIPAR(IB,1)*max(0.0_R8,
     .         fna(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1,is))/
     2         ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)
            end do
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            DO IS=0,NS-1
              shi0(bc_list_x(iy,ib),bc_list_y(iy,ib),3) =
     .        shi0(bc_list_x(iy,ib),bc_list_y(iy,ib),3)-
     1         ENIPAR(IB,1)*max(0.0_R8,
     1         -fna(rightix(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .              rightiy(bc_list_x(iy,ib),bc_list_y(iy,ib)),0,0,is))
     2         /ni(bc_list_x(iy,ib),bc_list_y(iy,ib),0)
            end do
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            DO IS=0,NS-1
              shi0(bc_list_x(iy,ib),bc_list_y(iy,ib),3) =
     .        shi0(bc_list_x(iy,ib),bc_list_y(iy,ib),3)-
     1         ENIPAR(IB,1)*max(0.0_R8,
     .         fna(bc_list_x(iy,ib),bc_list_y(iy,ib),0,0,is))/
     2         ni(bc_list_x(iy,ib),bc_list_y(iy,ib),0)
            end do
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=12, INVALID OPTION')
        END IF
#endif
* -- BCENI=13 -- PRESCRIBE THE ION ENERGY FLUX PER UNIT AREA
*             -- PROPORTIONAL TO TEMPERATURE

      case (13)
        call xerrab ('b2stbc_phys -- BCENI = 13 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     1     'BCENI = 13 : ion energy flux per area '//
     1     'proportional to temperature ',
     2     ENIPAR(IB,1),' on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENIPAR(IB,1)*
     .       TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=ENIPAR(IB,1)*
     .       TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=ENIPAR(IB,1)*
     .       TI(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     1       SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=ENIPAR(IB,1)*
     .       TI(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=13, INVALID OPTION')
        END IF
#endif
* -- BCENI=14 -- leakage option for ion energy

      case (14)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCCON = 14 : ion energy leakage ',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          cs   = sqrt(pz(iCv1)/rz(iCv1))
          srw%shi0(iCv1,0) = 0.0_R8
          srw%shi0(iCv1,3) = geo%fcS(iFc)*cs*enipar(ib,1)
          if(switch%mdf_fhi.ne.0) then
            srw%shi0(iCv1,0) = srw%shi0(iCv1,0) +
     &        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &       (dv%fhiPSch(iFc,0) + dv%fhiPSch(iFc,1))
            t0 = 0.0_R8
            do is=0,ns-1
              t0 = t0 + mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*2.5_R8*
     &          (geo%fcS(iFc)*geo%fcQalf(iFc,0)*dv%vadia(iFc,0,is) +
     &           geo%fcS(iFc)*geo%fcQalf(iFc,1)*dv%vadia(iFc,1,is))*
     &           0.5_R8*(pl%na(iCv1,is)+pl%na(iCv2,is))*
     &           0.5_R8*(pl%ti(iCv1)+pl%ti(iCv2))
            end do
            srw%shi0(iCv1,0) = srw%shi0(iCv1,0) +
     &        max(0.0_R8,t0)
            srw%shi0(iCv1,1) = srw%shi0(iCv1,1) +
     &        min(0.0_R8,t0)/pl%ti(iCv1)
          endif
        enddo

* -- BCENI=15 -- from b2stbc_spb                                      !srv 01.02.09 {

      case (15)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     &     'BCENI = 15 : from b2stbc_spb for ion temperature '//
     &     'with ENIPAR(IB,1) =',
     &     ENIPAR(IB,1),' on ', BCCHAR(IB),boundary_location(ib)      !srv 20.09.17
          write(*,*) 'BCENI=15 - assuming BCCON=14 for all ions'
          call xertst(ENIPAR(IB,1).gt.0.0_R8,
     1     'BCENI = 15, ENIPAR(IB,1) <= 0 not allowed!')              !srv 20.09.17
#ifdef WG_TODO
          IF(BCCHAR(IB).EQ.'W') THEN
            BCCON14_W(IB) = .TRUE.  ! IF BCCON=14 HAS BEEN IMPOSED AT WEST BOUNDARY IB FOR ALL IONS
            DO IS = 0, NS-1
              IF(.NOT.IS_NEUTRAL(IS)) THEN
                IF(BCCON(IS,IB).NE.14) THEN
                  BCCON14_W(IB) = .FALSE.
                ENDIF
              ENDIF
            ENDDO
            if(BCCON14_W(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 HAS ',
     &      'BEEN IMPOSED AT WEST BOUNDARY', IB,' FOR ALL IONS'
            if(.not.BCCON14_W(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 ',
     &      'HAS NOT BEEN IMPOSED AT WEST BOUNDARY', IB,' FOR ALL IONS'
          ENDIF
          IF(BCCHAR(IB).EQ.'E') THEN
            BCCON14_E(IB) = .TRUE.  ! IF BCCON=14 HAS BEEN IMPOSED AT WEST BOUNDARY IB FOR ALL IONS
            DO IS = 0, NS-1
              IF(.NOT.IS_NEUTRAL(IS)) THEN
                IF(BCCON(IS,IB).NE.14) THEN
                  BCCON14_E(IB) = .FALSE.
                ENDIF
              ENDIF
            ENDDO
            if(BCCON14_E(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 HAS ',
     &       'BEEN IMPOSED AT EAST BOUNDARY', IB,' FOR ALL IONS'
            if(.not.BCCON14_E(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 ',
     &       'HAS NOT BEEN IMPOSED AT EAST BOUNDARY', IB,' FOR ALL IONS'
          ENDIF
#endif
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          cs   = sqrt(pz(iCv1)/rz(iCv1))
          t0 = 0.0_R8
          t1 = 0.0_R8
          t2 = 0.0_R8
          if (bceni_15_style == 0) then
           do is = 0, ns-1
            if(.not.is_neutral(is)) then
              t0 = t0 + conpar(is,ib,1)*pl%na(iCv1,is)
              if (conpar(is,ib,2).gt.0.0_R8)
     &         t1 = t1 + co%dna0(iCv1,is)/conpar(is,ib,2)*pl%na(iCv1,is)
            endif
           enddo
           srw%shi0(iCv1,1) = srw%shi0(iCv1,1) -                        !srv 11.03.09 {
     &      enipar(ib,1)*t0*cs*abs(geo%fcPbs(iFc)) -
     &      enipar(ib,2)*t1*geo%fcS(iFc)                                !srv 19.04.13
          else
           !wdk sheath transmission factor consistent with truncated
           !    Maxwellian distribution on Eirene side (NEMODS=7)
           !    See eirene manual for details.
           do is = 0, ns-1
            if(.not.is_neutral(is)) then
             vt = sqrt(2.0_R8*pl%ti(iCv1)/(am(is)*mp))
             t3 = pit(iCv1)*abs(geo%fcQalf(iFc,0))*cs/vt
             t4 = sqrt((pit(iCv1)*abs(geo%fcQalf(iFc,1))*cs/vt)**2 +
     &                 (geo%cvBb(iCv1,2)/geo%cvBb(iCv1,3)*cs/vt)**2)    !WG_TODO correct for drift contributions!
             dsh = 1.0_R8 + t3**2 + t4**2 +                             !corrected for factor 3/2 in parallel convection
     &         0.5_R8*(sheath_g(t3) - 1.0_R8)/sheath_g(t3)
             t0 = t0 + dsh*conpar(is,ib,1)*pl%na(iCv1,is)
             if (conpar(is,ib,2).gt.0.0_R8)
     &        t1 = t1 + enipar(ib,2)*co%dna0(iCv1,is)/conpar(is,ib,2)*
     &          pl%na(iCv1,is)
             t2 = t2 + dv%fhm(iFc,0,is) + dv%fhm(iFc,1,is)
            endif
           enddo
           srw%shi0(iCv1,1) = srw%shi0(iCv1,1) -
     &      t0*cs*abs(geo%fcPbs(iFc)) - t1*geo%fcS(iFc) +
     &      t2*mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)
          endif
          if(switch%mdf_fhi.ne.0) then
            srw%shi0(iCv1,0) = srw%shi0(iCv1,0) +
     &        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &       (dv%fhiPSch(iFc,0) + dv%fhiPSch(iFc,1))
            t0 = 0.0_R8
            do is=0,ns-1
              t0 = t0 + mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*2.5_R8*
     &          (geo%fcS(iFc)*geo%fcQalf(iFc,0)*dv%vadia(iFc,0,is) +
     &           geo%fcS(iFc)*geo%fcQalf(iFc,1)*dv%vadia(iFc,1,is))*
     &           0.5_R8*(pl%na(iCv1,is)+pl%na(iCv2,is))*
     &           0.5_R8*(pl%ti(iCv1)+pl%ti(iCv2))
            end do
            srw%shi0(iCv1,0) = srw%shi0(iCv1,0) +
     &        max(0.0_R8,t0)
            srw%shi0(iCv1,1) = srw%shi0(iCv1,1) +
     &        min(0.0_R8,t0)/pl%ti(iCv1)
          endif
cwdk    Todo: corrections for BCCON.ne.14
        enddo                                                           !srv 01.02.09 }

* -- BCENI=16 -- Feedback boundary condition with constant temperature

      case (16)

        if(switch%bceni_16_style.eq.0) then

          if(ncall_b2stbc_phys.eq.0) then
            write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .       'ion temperature feedback by flux ',
     1       ENIPAR(IB,1),' with ',ENIPAR(IB,2),' on ',
     1       BCCHAR(IB),boundary_location(ib)
            call xertst(enipar(ib,1).ne.0, 'faulty parameter enipar1')
            call xertst(enipar(ib,2).ge.0, 'faulty parameter enipar2')
          endif
*       ..compute average temperature on boundary
          t0 = 0.0e0_R8
          t1 = 0.0e0_R8
          totflux = 0.0e0_R8
          !wdk Note on backwards compatibility: feedback now always on
          !wdk flux through the boundary itself (no offset)
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
           iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                     ! number of the guard cell
           iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
           t0 = t0 + geo%fcS(iFc)*dv%ni(iCv1,0)
           t1 = t1 + geo%fcS(iFc)*dv%ni(iCv1,0)*pl%ti(iCv1)
           !wdk at the moment, totflux only contains ion/neutral heat flux;
           !wdk bc_type16_kinetic_energy.eq.1 option not yet included
           totflux = totflux - (dv%fhi(iFc,0) + dv%fhi(iFc,1))*
     &       mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)
          enddo
          t1 = t1/t0
          totpar = enipar(ib,1)
          if (enipar(ib,2).gt.0.0_R8) t3 = enipar(ib,2)
          if (enipar(ib,2).le.0.0_R8) t3 = 1.0e-2_R8
          t2 = (1.0_R8-t3*(totflux-totpar)/
     &         max(abs(totflux),abs(totpar)))*t1
          if (switch%b2stbc_diagno.ge.2) then
           write(*,'(a,1p,4(1x,e15.6))')
     &      'bceni=16: curr_flux, desired_flux, old/new temperature, ',
     &      totflux,totpar,t1/qe,t2/qe
          endif
*       ..impose temperature
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%shi0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*t2
            srw%shi0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
          enddo

        elseif(switch%bceni_16_style.eq.1)then                          !srv 16.02.12 {

         call xerrab ('bceni_16_style.eq.1 not adapted for WG')
#ifdef WG_TODO
         if(ncall_b2stbc_phys.eq.0) then
           write(*,*) 'BCENI = 16 (style = 1) : '//
     &      'Given integrated ion heat flux with constant Ti '
!          write all boundaries with this type of boundary condition
           DO IBW=IB,NBC
             if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
               write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                   boundary_location(ibw)
               call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                       BC_LIST_Y(1,IBW)),
     .           'BCENI = 16 only applicable to core boundaries !')
               call xertst(enipar(ibw,1).ne.0,
     &           'faulty parameter enipar1')
             endif
             if(BCCHAR(IBW).EQ.'S'.and.
     &         (BCENI(IBW).EQ.23.or.BCENI(IBW).EQ.24.or.
     &          BCENI(IBW).EQ.26.or.BCENI(IBW).EQ.27)) then  ! IYS added BCENI= 26 and 27
               call xerrab
     .          ('BCENI types 16, 23, 24, 26 and 27 cannot be mixed!')
             endif
           ENDDO
         endif
         IF(BCCHAR(IB).EQ.'S') THEN
           if(.not.once_bceni) then
             if (switch%use_astra.ne.0) then
#ifdef ASTRA
               SENIPAR = fhib_astra
#else
               call xerrab
     &          ('Compile with -DASTRA option to couple to ASTRA!')
#endif
             else
               SENIPAR = 0.0_R8
               DO IBW=IB,NBC
                 if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
                   SENIPAR = SENIPAR + ENIPAR(IBW,1)
                 endif
               ENDDO
             endif
             write(*,*) 'BCENI=16: Summed value is ', SENIPAR
             TIS = 0.0_R8
             NIVS = 0.0_R8
             DO IBW=IB,NBC
               if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
                 DO IX = 1, BC_LIST_SIZE(IBW)
                   TIS = TIS +
     .                   TI(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     .                   NI(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     .                  VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                   NIVS = NIVS +
     .                   NI(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     .                  VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                 END DO
               endif
             END DO
             TIS = TIS/NIVS
             totflux = 0.0_R8
             DO IBW=IB,NBC
               IF(BCCHAR(IBW).EQ.'S'.AND.BCENI(IBW).EQ.BCENI(IB)) THEN
                 DO IX = 1, BC_LIST_SIZE(IBW)
                   if(switch%mdf_fhi.ne.0) then
                      totflux = totflux +
     &          fhi_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1)
     &               +fac_ExB(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*    !srv 11.04.13 {
     &                   ((po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))-
     &              po(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))/
     &                    (hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &              hx(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))+
     &            (po(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))-
     &                     po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))/
     &            (hx(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))+
     &                     hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))*
     &                     ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)*
     &                     ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                     SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                     bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),2)/
     &                     bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),3)**2 !srv 11.04.13 }
                   else
                     totflux = totflux +
     &             fhi(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1)
                   endif
                 ENDDO
               ENDIF
             ENDDO
             TIS = TIS *
     &         (1.0_R8-switch%b2stbc_bc_ref_ti*(totflux-SENIPAR)/                     !srv 11.04.13
     &           max(abs(totflux),SENIPAR))
             write(*,'(3(a,es14.7))')
     &         ' TIS=',TIS/ev,' totflux=',totflux,
     &         ' SENIPAR=',SENIPAR
             DO IBW=IB,NBC
               IF(BCCHAR(IBW).EQ.'S'.AND.BCENE(IBW).EQ.BCENE(IB)) THEN
                 DO IX = 1, BC_LIST_SIZE(IBW)
                   shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)=
     .               VBIG*NE00*TIS*
     .               SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                   shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)=
     .              -VBIG*NE00*
     .               SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                 ENDDO
               ENDIF
             ENDDO
             once_bceni = .true.
           endif                                                       !srv 16.02.12 }
         ELSE
           CALL XERRAB
     1         ('B2STBC -- BCENI=16, ONLY AVAILABLE FOR CORE!')
         ENDIF
#endif
        else
         CALL XERRAB (' Incorrect style of boundary condition')
        endif

* -- BCENI=17 -- Feedback boundary condition with constant temperature

      case (17)
        call xerrab ('b2stbc_phys -- BCENI = 17 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'BCENI = 17 : ion temperature feedback by flux ',
     1     ENIPAR(IB,1),' with ',ENIPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(enipar(ib,1).ne.0, 'faulty parameter enipar1')
          call xertst(enipar(ib,2).gt.0, 'faulty parameter enipar2')
          call xertst(bcene(ib).eq.bceni(ib),
     1     'both bcene and bceni should be set to 17')
          call xertst(enepar(ib,1).ne.0, 'faulty parameter enepar1')
          call xertst(enepar(ib,2).gt.0, 'faulty parameter enepar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
c already handled when the electron heat bc was handled
c see label 2980
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=17, NOT AVAILABLE FOR N,W,E')
        ENDIF
#endif
* -- BCENI=18 -- FRACTIONAL DROP CONDITION (FOR EDGE2D COMPARISON)

      case (18)
        call xerrab ('b2stbc_phys -- BCENI = 18 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENI = 18 : Ion temperature fractional drop',
     1     ENIPAR(IB,1),' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENIPAR(IB,1).lt.1.0_R8,
     1     'BCENI = 18, ENIPAR(,1) => 1')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENIPAR(IB,1)) *
     .       TI(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .           VBIG*NE00*T1*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .          -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENIPAR(IB,1)) *
     .       TI(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*T1*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENIPAR(IB,1)) *
     .       TI(RIGHTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          RIGHTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=
     .           VBIG*NE00*T1*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=
     .          -VBIG*NE00*SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            T1 = (1.0_R8 - ENIPAR(IB,1)) *
     .       TI(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)= VBIG*NE00*T1*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            shi0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=-VBIG*NE00*
     .       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=18, INVALID OPTION')
        ENDIF
#endif
* -- BCENI=19 -- PRESCRIBE THE DECAY LENGTH FOR THE ION TEMPERATURE

      case (19)
        call xerrab ('b2stbc_phys -- BCENI = 19 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCENI = 19 : ion temperature decay length (new) ',
     1     ENIPAR(IB,1),' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENIPAR(IB,1).ne.0.0_R8,
     1     'BCENI = 19, ENIPAR(,1) = 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .          -hci0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/ENIPAR(IB,1)*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            t0=0.0_R8
            DO IS=0,NS-1
               t0=t0-
     .              1.5_R8*fna_32(topix(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),
     .                            topiy(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),1,1,is)-
     .              2.5_R8*fna_52(topix(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),
     .                            topiy(bc_list_x(ix,ib),
     .                                  bc_list_y(ix,ib)),1,1,is)
            ENDDO
            if(t0.gt.0.0_R8) then
               shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .          shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .              t0 * TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            else
               shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .          shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) + t0
            endif
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
             shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
             shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .          -hci0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/ENIPAR(IB,1)*
     1            SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .               BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
             t0=0.0_R8
             DO IS=0,NS-1
                t0=t0+
     .          1.5_R8*fna_32(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1,is)+
     .          2.5_R8*fna_52(bc_list_x(ix,ib),bc_list_y(ix,ib),1,1,is)
             ENDDO
            if(t0.gt.0.0_R8) then
               shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .          shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     .              t0 * TI(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
             else
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .          shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) + t0
             endif
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=19, NOT AVAILABLE FOR W,E')
        ENDIF
#endif
* -- BCENI=20 -- Feedback boundary condition with constant temperature scaled by temperature on the ring bc_type20_ref away

      case (20)
        call xerrab ('b2stbc_phys -- BCENI = 20 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'BCENI = 20 : ion temperature feedback by flux (neo) ',
     1     ENIPAR(IB,1),' with ',ENIPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(enipar(ib,1).ne.0, 'faulty parameter enipar1')
          call xertst(enipar(ib,2).gt.0, 'faulty parameter enipar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          totflux = 0.0_R8
          t0 = 0.0_R8
          weight = 0.0_R8
          tiave = 0.0_R8
          niave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type16_ix=bc_list_x(ix,ib)
            bc_type16_iy=bc_list_y(ix,ib)
            do bc_type16_i=1,bc_type16_ref      ! how many steps away from the boundary?
              bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
              bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
              bc_type16_ix=bc_type16_ixx
              call xertst (-2.lt.bc_type16_ix.and.bc_type16_ix.lt.nx+1,
     1         'stepped off x-grid in type 16 bc')
              call xertst (-2.lt.bc_type16_iy.and.bc_type16_iy.lt.ny+1,
     1         'stepped off y-grid in type 16 bc')
            enddo
            totflux = totflux +
     .       fhi(bc_type16_ix,bc_type16_iy,1,1)
            t0 = t0 +
     .       ti(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)
            weight = weight + vol(bc_list_x(ix,ib),bc_list_y(ix,ib))*
     .       ni(bc_list_x(ix,ib),bc_list_y(ix,ib),0)
            bc_type20_ix=bc_list_x(ix,ib)
            bc_type20_iy=bc_list_y(ix,ib)
            do bc_type20_i=1,bc_type20_ref      ! how many steps away from the boundary?
              bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
              bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
              bc_type20_ix=bc_type20_ixx
              call xertst (-2.lt.bc_type20_ix.and.bc_type20_ix.lt.nx+1,
     1         'stepped off x-grid in type 20 bc')
              call xertst (-2.lt.bc_type20_iy.and.bc_type20_iy.lt.ny+1,
     1         'stepped off y-grid in type 20 bc')
            enddo
            tiave=tiave+ti(bc_type20_ix,bc_type20_iy)*
     1       ni(bc_type20_ix,bc_type20_iy,0)*
     2       vol(bc_type20_ix,bc_type20_iy)
            niave=niave+ni(bc_type20_ix,bc_type20_iy,0)*
     1       vol(bc_type20_ix,bc_type20_iy)
          ENDDO
          t0 = t0 / weight
          tiave=tiave/niave
          if(switch%b2stbc_diagno.ge.2) write(*,*) 'Totflux :',totflux
          T1 = (1.0_R8 -
     .     enipar(ib,2)*(totflux-enipar(ib,1))/enipar(ib,1)) * t0
          if(switch%b2stbc_diagno.ge.2)
     1     write(*,*) 'Old/New ion temp :',t0/qe, t1/qe
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type20_ix=bc_list_x(ix,ib)
            bc_type20_iy=bc_list_y(ix,ib)
            do bc_type20_i=1,bc_type20_ref      ! how many steps away from the boundary?
              bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
              bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
              bc_type20_ix=bc_type20_ixx
              call xertst (-2.lt.bc_type20_ix.and.bc_type20_ix.lt.nx+1,
     1         'stepped off x-grid in type 20 bc')
              call xertst (-2.lt.bc_type20_iy.and.bc_type20_iy.lt.ny+1,
     1         'stepped off y-grid in type 20 bc')
            enddo
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*T1*
     1       ti(bc_type20_ix,bc_type20_iy)/tiave
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=20, NOT AVAILABLE FOR N,W,E')
        ENDIF
#endif
!! FIXME: Clash here with 5.0/5.1 BCENI=21 type
cxpb Careful here because of BC number clash
* From 5.0: starting at version 01.001.030
* -- BCCON=21 -- constant ion temperature scaled by ion temperature on the ring bc_type21_ref away
* has been moved up to next available number (BCENI=25)

* From 5.2, after version number 01.001.078:
* -- BCENI=21 -- from b2stbc_spb                                     !srv 01.02.09 {

      case (21)
        call xerrab ('b2stbc_phys -- BCENI = 21 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          if (boundary_version.ge.'01.001.030' .and.
     &        boundary_version.lt.'02.000.000') then
            write(*,*) 'BCENI=21 type B.C. moved to 25 ',
     &       ' because of clash with 5.2 numbering. '
            write(*,*) 'Please modify b2.boundary.parameters file'
            call xerrab('B2STBC -- BCENI=25, DISPLACED')
          endif
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     &     'BCENI = 21 : from b2stbc_spb for BCENI with ENIPAR(IB,1) =',
     &     ENIPAR(IB,1),' and ENIPAR(IB,2) =', ENIPAR(IB,2),' on ',
     &     BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            s0 = vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            if(switch%mdf_fhi.ne.0) then                                    !srv 21.03.09 {
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &         fhiPSch(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1)
          endif                                                      !srv 21.03.09 }
            do is = 0, ns-1
              if(.not.is_neutral(is)) then                           !srv 19.03.09 {
                t0 = na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                csb =
     .           sqrt((rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &                  te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &                  ti(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/
     &               (am(is)*mp))
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)+ENIPAR(IB,2)*
     &                            t0*csb*s0
!     &                            t0*csb*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
              endif
            enddo
            SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &       SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1) -
     &       max(0.0e0_R8,-floi_noc(
     &        topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &        topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1))         !srv 28.10.02
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            s0 = vol(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/
     &           hy1(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            if(switch%mdf_fhi.ne.0) then                             !srv 21.03.09 {
              shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &         shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &         fhiPSch(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1)
            endif                                                    !srv 21.03.09 }
            do is = 0, ns-1
              if(.not.is_neutral(is)) then                           !srv 19.03.09 }
                t0 = na(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)
                csb = sqrt((rza(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),is)*
     &                       te(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))+
     &                       ti(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/
     &                     (am(is)*mp))
                shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &           shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)+ENIPAR(IB,2)*
     &          t0*csb*s0
!     &          t0*csb*SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
!     &                    BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
              endif                                                  !srv 19.03.09 {
            enddo
            SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     &       SHI0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)-
     &       max(0.0e0_R8,
     &        floi_noc(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1))       !srv 28.10.02
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=21, NOT AVAILABLE FOR W,E')
        ENDIF                                                        !srv 01.02.09 }
#endif
* -- BCENI=22 -- RADIAL LEAKAGE CONDITION FOR THE ION TEMPERATURE    !srv 23.07.10 {
*                SPECIFIES A TEMPERATURE GRADIENT SUCH THAT THE DIFFUSIVE
*                FLUX IS SET TO BE EQUIVALENT TO A LEAKAGE VELOCITY, IN
*                UNITS OF THE COLLECTIVE ION THERMAL VELOCITY

      case (22)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     &     'BCENI = 22 : with ENIPAR(IB,1)= ',ENIPAR(IB,1),
     &     ' on ', BCCHAR(IB),boundary_location(ib)
        endif
cwdk    Notes:
cwdk     * need correction for misaligned grids -> t1
cwdk     * need to correct this BC for hybrid scheme
        do ibc = 1, mpg%bcCvP(ib,2)                                    ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                     ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                     ! number of the domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
          t0   = geo%fcQalf(iFc,1)
          t1   = (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     &            geo%fcQbet(iFc,0)
          t2   = 0.0_R8
          do is = 0, ns-1
            if(.not.is_neutral(is)) then
              csb = sqrt((rt%rza(iCv1,is)*pl%te(iCv1) + pl%ti(iCv1))/
     &                  (am(is)*mp))
              t2  = t2 + csb*pl%na(iCv1,is)
            endif
          enddo
          srw%shi0(iCv1,0) =  VBIG*NE00*EV*
     &        co%chci(iFc,1)*(t0*pl%ti(iCv2) -
     &        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*t1*
     &       (pl%ti(mpg%fcVx(iFc,2)) - pl%ti(mpg%fcVx(iFc,1))))
          srw%shi0(iCv1,1) = -VBIG*NE00*EV*
     &       (co%chci(iFc,1)*t0 -
     &       ENIPAR(iB,1)*geo%fcS(iFc)*t2)
        enddo                                                        !srv 23.07.10 }

* -- BCENI=23 -- PRESCRIBE THE POLOIDALLY AVERAGED VALUE OF THE ION TEMPERATURE
* AND INTRODUCE A POLOIDAL VARIATION AS CLOSE AS POSSIBLE TO NEOCLASSICAL SOLUTION.
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH CORRESPONDING
* CONDITION ON ION DENSITY BCCON = 21,22,23

      case (23)
        call xerrab ('b2stbc_phys -- BCENI = 23 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7)')
     .     'BCENI = 23 : specified ion temperature ',ENIPAR(IB,1)
!          write all boundaries with this type of boundary condition
          DO IBW=IB,NBC
            IF (BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB)) THEN
              write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                  boundary_location(ibw)
              call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                      BC_LIST_Y(1,IBW)),
     .         'BCENI = 23 only applicable to core boundaries !')
              call xertst(switch%use_astra.ne.0.or.
     .         ENIPAR(IB,1).eq.ENIPAR(IBW,1),
     1         'BCENI = 23, inconsistent ENIPAR(,1) values!')
              if (.not.(bccon(ismain,ib).eq.21 .or.
     .                  bccon(ismain,ib).eq.22 .or.
     .                  bccon(ismain,ib).eq.23)) then
                write(0,'(a,i3,a,i2,a)') 'Should use BCCON(',ismain,
     &             ',',IB,').EQ.21,22,23 together with BCENI=23'
                write(*,'(a,i3,a,i2,a)') 'Should use BCCON(',ismain,
     &             ',',IB,').EQ.21,22,23 together with BCENI=23'
                write(*,'(1x,a)')
     &           'Please ensure that this boundary condition '//
     &           'is consistent with '//
     &           'the ion density boundary condition'  !srv 21.09.12
              end if
              if (.not.(ns.eq.1.or.(ns.eq.2.and.is_neutral(0)))) then
                write(0,*)
     &           'Use of BCENI=23 only meant for single-fluid runs!'
                write(0,*) 'Use BCENI=26 instead!'
                write(6,*)
     &           'Use of BCENI=23 only meant for single-fluid runs!'
                write(6,*) 'Use BCENI=26 instead!'
              endif
            endif
            if(BCCHAR(IBW).EQ.'S'.and.
     &         (BCENI(IBW).EQ.16.or.BCENI(IBW).EQ.24.or.
     &          BCENI(IBW).EQ.26.or.BCENI(IBW).EQ.27)) then  ! IYS added BCENI= 26 and 27
                call xerrab
     .           ('BCENI types 16, 23, 24, 26 and 27 cannot be mixed!')
            endif
          ENDDO
          call xertst(switch%use_astra.ne.0.or.ENIPAR(IBW,1).gt.0.0_R8,
     1     'BCENI = 23, ENIPAR(,1) !> 0')
        endif

        if(.not.once_bceni) then
          if (switch%use_astra.ne.0) THEN
#ifdef ASTRA
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
                ENIPAR(IBW,1) = tib_astra*1.0e3_R8                     !srv 14.06.12
              endif
            ENDDO
#else
            call xerrab
     .       ('Compile with -DASTRA option to couple with ASTRA!')
#endif
          ENDIF

          IF(BCCHAR(IB).EQ.'S' .and.
     &     on_closed_surface(BC_LIST_X(1,IB),BC_LIST_Y(1,IB))) THEN
!          DO IX = 1, BC_LIST_SIZE(IB)
!            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)= VBIG*NE00*ENIPAR(IB,1)
!     .       *EV*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
!            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
!     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
!          end do
            t0=0.0e0_R8   ! Volume                                     !lk 12.12.12 {
            t1=0.0e0_R8   ! Poloidally volume-averaged ion temperature in the cell next to boundary
            t2=0.0e0_R8   ! Poloidally volume-averaged electron density in the cell next to boundary
            t3=0.0e0_R8   ! Poloidally volume-averaged electron density in the boundary cell
            DO IBW=IB,NBC
              IF(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  t0=t0 +
     &              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t1=t1 +
     &               ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t2=t2 +
     &               ne(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t3=t3 + ne(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &             vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                END DO
              ENDIF
            ENDDO
            t1=t1/t0
            t2=t2/t0
            t3=t3/t0
            DO IBW=IB,NBC
              IF(BCCHAR(IBW).EQ.'S'.AND.BCENI(IBW).EQ.BCENI(IB)) THEN
                do IX = 1, BC_LIST_SIZE(IBW)
                  if(on_closed_surface(BC_LIST_X(IX,IBW),
     .                                 BC_LIST_Y(IX,IBW))) then
                    shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)=
     &               shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)
     &                +VBIG*NE00*(ENIPAR(IBW,1)*EV
     &                +(ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                     topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
     &                  -t1)*t1*t2/t3/ENIPAR(IBW,1)/EV)
     &                *SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)=
     &               -VBIG*NE00*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  endif
                enddo
              ENDIF                                                 !lk 12.12.12 }
            ENDDO
            once_bceni = .true.
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCENI=23 IS AVAILABLE FOR CORE BOUNDARY ONLY.')
          END IF
        endif
#endif
* -- BCENI=24 -- Feedback boundary condition with prescribed total ion heat flux,
*                constant poloidally averaged ion temperature and a
*                poloidal variation as close as possible to neoclassical solution.
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH CORRESPONDING
* CONDITION ON ION DENSITY BCCON = 21,22,23

      case (24)
        call xerrab ('b2stbc_phys -- BCENI = 24 not adapted for WG')
#ifdef WG_TODO
!IYS 04.03.2013 {
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a)') 'BCENI = 24 : Given integrated ion heat flux',
     &                   ' with constant Ti '
!          write all boundaries with this type of boundary condition
          DO IBW=IB,NBC
            if (BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB)) THEN
              write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                  boundary_location(ibw)
              call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                      BC_LIST_Y(1,IBW)),
     .         'BCENI = 24 only applicable to core boundaries !')
              call xertst(enipar(ibw,1).ne.0,'faulty parameter enipar1')
              if (.not.(bccon(ismain,ib).eq.21 .or.
     .                  bccon(ismain,ib).eq.22 .or.
     .                  bccon(ismain,ib).eq.23)) then
                  write(0,'(a,i3,a,i2,a)') 'Should use BCCON(',ismain,
     &             ',',IB,').EQ.21,22,23 together with BCENI=24'
                  write(*,'(a,i3,a,i2,a)') 'Should use BCCON(',ismain,
     &             ',',IB,').EQ.21,22,23 together with BCENI=24'
                  write(*,'(1x,a)')
     &             'Please ensure that this boundary condition '//
     &             'is consistent with '//
     &             'the ion density boundary condition'  !srv 21.09.12
              end if
              if (.not.(ns.eq.1.or.(ns.eq.2.and.is_neutral(0)))) then
               write(0,*)
     &          'Use of BCENI=24 only meant for single-fluid runs!'
               write(0,*) 'Use BCENI=27 instead!'
               write(6,*)
     &          'Use of BCENI=24 only meant for single-fluid runs!'
               write(6,*) 'Use BCENI=27 instead!'
              endif
            endif
            if(BCCHAR(IBW).EQ.'S'.and.
     &         (BCENI(IBW).EQ.16.or.BCENI(IBW).EQ.23.or.
     &          BCENI(IBW).EQ.26.or.BCENI(IBW).EQ.27)) then  ! IYS added BCENI= 26 and 27
                call xerrab
     .           ('BCENI types 16, 23, 24, 26 and 27 cannot be mixed!')
            endif
          ENDDO
        endif

        if(.not.once_bceni) then
          if (switch%use_astra.ne.0) then
#ifdef ASTRA
            SENIPAR = fhib_astra
#else
            call xerrab
     .       ('Compile with -DASTRA option to couple with ASTRA!')
#endif
          else
            SENIPAR = 0.0_R8
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
                SENIPAR = SENIPAR + ENIPAR(IBW,1)
              endif
            ENDDO
          endif
          write(*,*) 'BCENI=24: Summed value is ', SENIPAR

          IF(BCCHAR(IB).EQ.'S' .and.
     $     on_closed_surface(BC_LIST_X(1,IB),BC_LIST_Y(1,IB))) THEN
            TIS = 0.0_R8  ! Poloidally density-averaged ion temperature in the boundary cell     !IYS 21.01.2016
            NIVS = 0.0_R8
            totflux = 0.0_R8
            t0=0.0e0_R8   ! Volume                                     !lk 12.12.12 {
            t1=0.0e0_R8   ! Poloidally volume-averaged ion temperature in the cell next to boundary
            t2=0.0e0_R8   ! Poloidally volume-averaged electron density in the cell next to boundary
            t3=0.0e0_R8   ! Poloidally volume-averaged electron density in the boundary cell
            DO IBW=IB,NBC
              IF(BCCHAR(IBW).EQ.'S'.AND.BCENI(IBW).EQ.BCENI(IB)) THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  t0=t0 +
     &               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t1=t1 +
     &                ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  TIS=TIS + ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                      ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     &               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  NIVS=NIVS +
     &                     ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     &              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t2=t2 +
     &                ne(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t3=t3 + ne(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
     &              *vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  if(switch%mdf_fhi.ne.0) then
                    totflux = totflux +
     &        fhi_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1)
     &             +fac_ExB(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*      !srv 11.04.13 {
     &                 ((po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))-
     &            po(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))/
     &                  (hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &            hx(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))+
     &          (po(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &              rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))-
     &                   po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))/
     &          (hx(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &              rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))+
     &                   hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))*
     &                   ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)*
     &                   ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                   SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                   bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),2)/
     &                   bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),3)**2  !srv 11.04.13 }
                  else
                    totflux = totflux +
     &               fhi(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1)
                  endif
                ENDDO
              ENDIF
            ENDDO
            TIS=TIS/NIVS
            t1=t1/t0
            t2=t2/t0
            t3=t3/t0
            TIS=TIS*
     &          (1.0_R8-switch%b2stbc_bc_ref_ti*(totflux-SENIPAR)/
     &            max(abs(totflux),SENIPAR))
            write(*,'(3(a,es14.7))')
     &        ' TIS=',TIS/ev,' totflux=',totflux,
     &        ' SENIPAR=',SENIPAR
            DO IBW=IB,NBC
              IF(BCCHAR(IBW).EQ.'S'.AND.BCENE(IBW).EQ.BCENE(IB)) THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)=
     .             VBIG * NE00 * (TIS +
     .             (ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                 topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))-t1)*
     .              t1*t2/t3/TIS)*
     .              SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)=
     .             -VBIG*NE00*
     .              SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                ENDDO
              ENDIF
            ENDDO
            once_bceni = .true.
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCENI=24 IS AVAILABLE FOR CORE ONLY')
          ENDIF
        endif                                                          !srv 16.02.12 }
#endif
* -- BCENI=25 -- constant ion temperature scaled by temperature on the ring bc_type21_ref away
* was BCENI=21 in 5.0/5.1

      case (25)
        call xerrab ('b2stbc_phys -- BCENI = 25 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'BCENI = 25 : fixed ion temperature '//
     .     '(modified by nearby ring) ',
     1     ENIPAR(IB,1),' with ',ENIPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(enipar(ib,1).ne.0.0_R8,
     1     'faulty parameter enipar1')
          call xertst(enipar(ib,2).gt.0.0_R8,
     1     'faulty parameter enipar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          tiave = 0.0_R8
          niave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            tiave=tiave+ti(bc_type21_ix,bc_type21_iy)*
     1       ni(bc_type21_ix,bc_type21_iy,0)*
     2       vol(bc_type21_ix,bc_type21_iy)
            niave=niave+ni(bc_type21_ix,bc_type21_iy,0)*
     1       vol(bc_type21_ix,bc_type21_iy)
          ENDDO
          tiave=tiave/niave
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*ENIPAR(IB,1)*EV*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .       ((1.0_R8-enipar(ib,2))+
     .           enipar(ib,2)*ti(bc_type21_ix,bc_type21_iy)/tiave)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -VBIG*NE00*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE IF(BCCHAR(IB).EQ.'N') THEN
          tiave = 0.0_R8
          niave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            tiave=tiave+ti(bc_type21_ix,bc_type21_iy)*
     1       ni(bc_type21_ix,bc_type21_iy,0)*
     2       vol(bc_type21_ix,bc_type21_iy)
            niave=niave+ni(bc_type21_ix,bc_type21_iy,0)*
     1       vol(bc_type21_ix,bc_type21_iy)
          ENDDO
          tiave=tiave/niave
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*NE00*ENIPAR(IB,1)*EV*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .       ((1.0_R8-enipar(ib,2))+
     .           enipar(ib,2)*ti(bc_type21_ix,bc_type21_iy)/tiave)
            shi0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=25, NOT AVAILABLE FOR W,E')
        ENDIF
#endif
* -- BCENI=26 -- PRESCRIBE THE POLOIDALLY AVERAGED VALUE OF THE ION TEMPERATURE
* AND INTRODUCE A POLOIDAL VARIATION IN A SIMPLIFIED MANNER.
* THIS BOUNDARY CONDITION IS SUITABLE FOR ANY PLASMA COMPOSITION (i.e. when BCENI=23 fails).
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH CORRESPONDING
* CONDITION ON ION DENSITY BCCON = 25,26,27

      case (26)
        call xerrab ('b2stbc_phys -- BCENI = 26 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7)')
     .     'BCENI = 26 : specified ion temperature ',ENIPAR(IB,1)
!         write all boundaries with this type of boundary condition
          DO IBW=IB,NBC
            IF (BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB)) THEN
              write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                  boundary_location(ibw)
              call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                      BC_LIST_Y(1,IBW)),
     .         'BCENI = 26 only applicable to core boundaries !')
              call xertst(switch%use_astra.ne.0.or.
     .                                   ENIPAR(IBW,1).gt.0.0_R8,
     1         'BCENI = 26, ENIPAR(,1) !> 0')
              call xertst(switch%use_astra.ne.0.or.
     .         ENIPAR(IB,1).eq.ENIPAR(IBW,1),
     1         'BCENI = 26, inconsistent ENIPAR(,1) values!')
            endif
            if(BCCHAR(IBW).EQ.'S'.and.
     &         (BCENI(IBW).EQ.16.or.BCENI(IBW).EQ.23.or.
     &          BCENI(IBW).EQ.24.or.BCENI(IBW).EQ.27)) then  ! IYS added BCENI= 26 and 27
                call xerrab
     .           ('BCENI types 16, 23, 24, 26 and 27 cannot be mixed!')
            endif
          ENDDO
        endif

        if(.not.once_bceni) then
          if (switch%use_astra.ne.0) THEN
#ifdef ASTRA
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB)) THEN
                ENIPAR(IBW,1) = tib_astra*1.0e3_R8                     !srv 14.06.12
              endif
            ENDDO
#else
            call xerrab
     .       ('Compile with -DASTRA option to couple with ASTRA!')
#endif
          ENDIF
          IF(BCCHAR(IB).EQ.'S' .and.
     &     on_closed_surface(BC_LIST_X(1,IB),BC_LIST_Y(1,IB))) THEN
            t0=0.0e0_R8   ! Volume                                     !lk 12.12.12 {
            t1=0.0e0_R8   ! Poloidally volume-averaged ion temperature in the cell next to boundary
            DO IBW=IB,NBC
              IF(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
                DO IX= 1, BC_LIST_SIZE(IBW)
                  t0 = t0 +
     .              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t1 = t1 +
     .               ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     .              vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                  topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                END DO
              END IF
            END DO
            t1=t1/t0
            DO IBW=IB,NBC
              IF(BCCHAR(IBW).EQ.'S'.AND.BCENI(IBW).EQ.BCENI(IB)) THEN
                do IX = 1, BC_LIST_SIZE(IBW)
                  if(on_closed_surface(BC_LIST_X(IX,IBW),
     .                                 BC_LIST_Y(IX,IBW))) then
                    shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)=
     &               shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)
     &                +VBIG*NE00*(ENIPAR(IBW,1)*EV
     &                +(ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                     topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
     &                 -t1))
     &                *SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)=
     &               -VBIG*NE00*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  endif
                enddo
              ENDIF                                                 !lk 12.12.12 }
            ENDDO
            once_bceni = .true.
          ELSE
            CALL XERRAB
     1       ('B2STBC -- BCENI=26 IS AVAILABLE FOR CORE BOUNDARY ONLY.')
          END IF
        endif
#endif
* -- BCENI=27 -- Feedback boundary condition with prescribed total ion heat flux,
*                constant poloidally averaged ion temperature and a
*                poloidal variation in a simplified manner.
* THIS BOUNDARY CONDITION IS SUITABLE FOR ANY PLASMA COMPOSITION (i.e. when BCENI=24 fails).
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH CORRESPONDING
* CONDITION ON ION DENSITY BCCON = 25,26,27

      case (27)
        call xerrab ('b2stbc_phys -- BCENI = 27 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,*) 'BCENI = 27 : Given integrated ion heat flux',
     &               ' with constant Ti '
!         write all boundaries with this type of boundary condition
          DO IBW=IB,NBC
            if (BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
              write(*,'(a,a,a)') ' on ',BCCHAR(IBW),
     .                                  boundary_location(ibw)
              call xertst(on_closed_surface(BC_LIST_X(1,IBW),
     .                                      BC_LIST_Y(1,IBW)),
     .         'BCENI = 27 only applicable to core boundaries !')
              call xertst(enipar(ibw,1).ne.0,'faulty parameter enipar1')
            endif
            if(BCCHAR(IBW).EQ.'S'.and.
     &        (BCENI(IBW).EQ.16.or.BCENI(IBW).EQ.23.or.
     &         BCENI(IBW).EQ.24.or.BCENI(IBW).EQ.26)) then  ! IYS added BCENI= 26 and 27
              call xerrab
     .         ('BCENI types 16, 23, 24, 26 and 27 cannot be mixed!')
            endif
          ENDDO
        endif

        IF(BCCHAR(IB).EQ.'S' .and.
     &    on_closed_surface(BC_LIST_X(1,IB),BC_LIST_Y(1,IB))) THEN
          if(.not.once_bceni) then
            if (switch%use_astra.ne.0) then
#ifdef ASTRA
              SENIPAR = fhib_astra
#else
              call xerrab
     .         ('Compile with -DASTRA option to couple with ASTRA!')
#endif
            else
              SENIPAR = 0.0_R8
              DO IBW=IB,NBC
                if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN
                  SENIPAR = SENIPAR + ENIPAR(IBW,1)
                endif
              ENDDO
            endif
            write(*,*) 'BCENI=27: Summed value is ', SENIPAR
            TIS = 0.0_R8 ! Poloidally density-averaged ion temperature in the boundary cell
            NIVS = 0.0_R8
            t0=0.0e0_R8   ! Volume                                      !lk 12.12.12 {
            t1=0.0e0_R8   ! Poloidally volume-averaged ion temperature in the cell next to boundary
            totflux = 0.0_R8   ! IYS 22.03.2017
            DO IBW=IB,NBC
              if(BCCHAR(IBW).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB)) THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  TIS=TIS + ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     .                      ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  NIVS=NIVS +
     .                      ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)*
     .               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     .                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t0=t0 +
     &               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  t1=t1 +
     &                ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))*
     &               vol(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                   topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
                  if(switch%mdf_fhi.ne.0) then
                    totflux = totflux +
     &        fhi_mdf(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1)
     &             +fac_ExB(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*      !srv 11.04.13 {
     &                 ((po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))-
     &            po(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))/
     &                  (hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))+
     &            hx(leftix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               leftiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))+
     &          (po(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &              rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))-
     &                   po(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))/
     &          (hx(rightix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &              rightiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))+
     &                   hx(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))))*
     &                   ni(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)*
     &                   ti(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                   SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     &                   bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),2)/
     &                   bb(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),3)**2  !srv 11.04.13 }
                  else
                    totflux = totflux +
     &           fhi(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &               topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),1,1)
                  endif
                END DO
              endif
            END DO
            TIS=TIS/NIVS
            t1=t1/t0
            TIS=TIS*
     &          (1.0_R8-switch%b2stbc_bc_ref_ti*(totflux-SENIPAR)/
     &            max(abs(totflux),SENIPAR))
            write(*,'(3(a,es14.7))')
     &        ' TIS=',TIS/ev,' totflux=',totflux,
     &        ' SENIPAR=',SENIPAR
            DO IBW=IB,NBC
              IF(BCCHAR(IBW).EQ.'S'.AND.BCENE(IBW).EQ.BCENE(IB)) THEN
                DO IX = 1, BC_LIST_SIZE(IBW)
                  shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0)=
     .             VBIG*NE00*
     $             (TIS+(ti(topix(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)),
     &                      topiy(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)))
     &                  -t1))*
     .              SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                  shi0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1)=
     .             -VBIG*NE00*
     .              SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                ENDDO
              ENDIF
            ENDDO
            once_bceni = .true.
          endif                                                        !srv 16.02.12 }
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCENI=27 IS AVAILABLE FOR CORE ONLY')
        ENDIF
#endif

* -- BCENI=28 -- Prescribe the value of the ion temperature, with poloidal
*                perturbation taken from flux tube just inside domain

      case (28)
        if(ncall_b2stbc_phys.eq.0) then
          do ibc = 1, mpg%bcCvP(ib,2)
           iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
           call xertst(mpg%cvOnClosedSurface(iCv1),
     .      'BCENI=28 intended for core boundaries only!')
          enddo
          IF(switch%use_astra.ne.0) THEN
          ENIPAR(IB,1) = teb_astra*1.0e3_R8
          ENDIF
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified ion temperature + perturbation',ENIPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENIPAR(IB,1).gt.0.0_R8,
     1     'BCENI = 28, ENIPAR(,1) !> 0')
        endif
*     ..compute average temperature on core boundary
        t0 = 0.0e0_R8
        t1 = 0.0e0_R8
        do ib0 = 1, nbc
         ! identify core boundaries; assume all core boundaries
         ! form single flux surface
         ! implicit assumption at the moment: first ring in the
         ! core is a 'structured grid' ring with quad cells
         iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
         if (mpg%cvOnClosedSurface(iCv1)) then
          do ibc0 = 1, mpg%bcCvP(ib0,2)                                 ! loop over number of cells in the boundary
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                   ! number of the guard cell
           iCv2 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,2)                   ! number of the domain cell
           iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
           t0 = t0 + geo%cvVol(iCv2)
           t1 = t1 + geo%cvVol(iCv2)*pl%ti(iCv2)
          enddo
         endif
        enddo
        t1 = t1/t0
*     ..impose temperature with perturbation
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                      ! number of the domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          srw%shi0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*
     .      (pl%ti(iCv2) + enipar(iB,1)*EV - t1)
          srw%shi0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
        enddo

* -- BCENI=29 -- Feedback on the total ion energy flux, by imposing a flux surface
*                averaged temperature, with poloidal perturbation taken from flux
*                tube just inside domain

      case (29)
        if(ncall_b2stbc_phys.eq.0) then
          do ibc = 1, mpg%bcCvP(ib,2)
           iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
           call xertst(mpg%cvOnClosedSurface(iCv1),
     .      'BCENI=29 intended for core boundaries only!')
          enddo
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified ion heat flux',ENIPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENIPAR(IB,1).gt.0.0_R8,
     1     'BCENI = 29, ENIPAR(,1) !> 0')
        endif
*     ..compute average temperature on core boundary
        t0 = 0.0e0_R8
        t1 = 0.0e0_R8
        t2 = 0.0e0_R8
        totflux = 0.0e0_R8
        totpar  = 0.0e0_R8
        do ib0 = 1, nbc
         ! identify core boundaries; assume all core boundaries
         ! form single flux surface
         ! implicit assumption at the moment: first ring in the
         ! core is a 'structured grid' ring with quad cells
         iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
         if (mpg%cvOnClosedSurface(iCv1)) then
          if (bceni(ib0).ne.29)
     .      call xerrab ('BCENI=29: incompatible core BCs')
          do ibc0 = 1, mpg%bcCvP(ib0,2)                                 ! loop over number of cells in the boundary
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                   ! number of the guard cell
           iCv2 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,2)                   ! number of the domain cell
           iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
           t0 = t0 + geo%cvVol(iCv2)
           t1 = t1 + geo%cvVol(iCv2)*pl%ti(iCv1)
           t2 = t2 + geo%cvVol(iCv2)*pl%ti(iCv2)
           totflux = totflux - (dv%fhi(iFc,0) + dv%fhi(iFc,1))*
     &       mpg%bcFcOr(mpg%bcCvP(ib0,1)+ibc0-1)
          enddo
          totpar = totpar + enipar(ib0,1)
         endif
        enddo
        t1 = t1/t0
        t2 = t2/t0
        if (enipar(ib,2).gt.0.0_R8) t3 = enipar(ib,2)
        if (enipar(ib,2).le.0.0_R8) t3 = 1.0e-2_R8
        t4 = (1.0_R8-t3*(totflux-totpar)/
     &         max(abs(totflux),abs(totpar)))*t1
*     ..impose temperature with perturbation
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                      ! number of the domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          srw%shi0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*
     .      (pl%ti(iCv2) + t4 - t2)
          srw%shi0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
        enddo

      case default
        write(*,*) BCENI(IB)
        call xerrab ('BCENI out of bounds')

      end select

* -- BOUNDARY CONDITIONS FOR POTENTIAL EQUATION --

      select case (BCPOT(IB))

      case (0)
        continue

* -- BCPOT=1 -- PRESCRIBE THE VALUE OF THE POTENTIAL

      case (1)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCPOT =  1 : specified potential ',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
          srw%sch0(iCv1,0) =  PBIG*geo%fcS(iFc)*potpar(iB,1)
          srw%sch0(iCv1,1) = -PBIG*geo%fcS(iFc)
        enddo

* -- BCPOT=2 -- PRESCRIBE THE GRADIENT OF THE POTENTIAL

      case (2)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCPOT =  2 : potential gradient ',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                    ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
          iVx1 = mpg%fcVx(iFc,1)
          iVx2 = mpg%fcVx(iFc,2)
          t1 = intvertex_s(iVx1,nCv,nVx,mpg,geo%vxVol,pl%po)
          t2 = intvertex_s(iVx2,nCv,nVx,mpg,geo%vxVol,pl%po)
          srw%sch0(iCv1,0)=  PBIG*geo%fcS(iFc)*
     .       (potpar(ib,1)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2))*
     .        geo%fcQgam(iFc,0) - cor9*(t2 - t1)*geo%fcQgam(iFc,1)*
     .        (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     .        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1) + pl%po(iCv2))
          srw%sch0(iCv1,1) = -PBIG*geo%fcS(iFc)
#ifdef B25_EIRENE
*     ..set sheath potential drop for Eirene
          if (switch%eirene_sheath_pot.eq.1) then
            delta_sheathb(iFc) =
     &          qe*(pl%po(iCv1)-potpar(ib,2))/pl%te(iCv1)
          endif
#endif
        enddo

* -- BCPOT=3 -- SHEATH CONDITIONS, ELECTRON ENERGY TRANSMISSION

      case (3)
        call xerrab ('b2stbc_phys -- BCPOT = 3 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,a,a,1p,g14.7,a)')
     .     'BCPOT =  3 : sheath potential bc on ',BCCHAR(IB),
     1     boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            s1 = abs(pbs(rightix(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .                   rightiy(bc_list_x(iy,ib),bc_list_y(iy,ib)),0))
            vte = sqrt(bcint(te,nx,ny,iy,ib)/me)
            t0 = 0.0_R8
            do is = 0, ns-1
              t0 = t0-bcint(rza(-1,-1,is),nx,ny,iy,ib)*
     1         fna(rightix(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .             rightiy(bc_list_x(iy,ib),bc_list_y(iy,ib)),0,0,is)
            enddo
            do is = 0, ns_ext-1
              t0 = t0-bcint(za_ext(-1,-1,is),nx,ny,iy,ib)*
     1         fa_ext(rightix(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .                rightiy(bc_list_x(iy,ib),bc_list_y(iy,ib)),0,0,is)
            enddo
            fchi = qe*t0
            if(t0.lt.0.0_R8 .and.
     .       bottomiy(bc_list_x(iy,ib),bc_list_y(iy,ib)).ne.-2 .and.
     &       topiy(bc_list_x(iy,ib),bc_list_y(iy,ib)).ne.ny+1) then
              wrong_flow = .true.                                     !xpb
              if(switch%b2stbc_diagno.ge.0) then
                write(*,'(a,2i4,1x,1p,g14.7)') 'b2stbc_west:ix,iy,fchi',
     1           bc_list_x(iy,ib),bc_list_y(iy,ib),fchi
                write(hlp_frm,'(a,i3,a)') '(a,1p,',ns,'(g14.7,1x))'
                write(*,hlp_frm) 'rza : ',
     .           (rza(rightix(bc_list_x(iy,ib),
     .                        bc_list_y(iy,ib)),
     .                rightiy(bc_list_x(iy,ib),
     .                        bc_list_y(iy,ib)),is),is=0,ns-1)
                write(*,hlp_frm) 'fna (A) : ',
     &           (fna(rightix(bc_list_x(iy,ib),
     .                        bc_list_y(iy,ib)),
     .                rightiy(bc_list_x(iy,ib),
     .                        bc_list_y(iy,ib)),0,0,is)*qe,is=0,ns-1)
                write(*,'(a,1p,4g14.7)') 'Bx,Bz,Ud,-Ua =',
     &           wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0),
     .           wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),2),
     &           ud(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain),
     .          -ua(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain)
               endif
c              fchi=0.0_R8
            endif
            fche = sqrt(1/(2*pi))*qe*bcint(ne,nx,ny,iy,ib)*vte*s1*
     &       expu2(-qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
     .       bcint(te,nx,ny,iy,ib))
            if (switch%secmodel.eq.1) then
              if (ne(bc_list_x(iy,ib),bc_list_y(iy,ib))/
     .           (rza(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain)*
     .            na(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain))
     .           .gt.1.1_R8) then
                write(*,*) 'Large ne contribution from impurities!'
                write(*,*) 'Kinetic s.e.c. model may not be valid!'
              endif
              call eseec0 (bcint(te,nx,ny,iy,ib),
     &         bcint(ti,nx,ny,iy,ib),bcint(ne,nx,ny,iy,ib),
     &         bcint(po,nx,ny,iy,ib),
     .         wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),0),
     &         wbbr(bc_list_x(iy,ib),bc_list_y(iy,ib),3),0,
     .         bc_list_x(iy,ib),bc_list_y(iy,ib),
     &         nx,ny,ns,ismain,seec,switch%b2stbc_diagno)               !cfs
            else
              seec = gammae
            endif
            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/bcint(te,nx,ny,iy,ib)
*    ..compute charge source
            sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),0) =
     .       sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),0)+
     &       (1.0_R8-seec)*fche-fchi+
     .       t0*po(bc_list_x(iy,ib),bc_list_y(iy,ib))
            sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),1) =
     .      sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),1)-t0
#ifdef B25_EIRENE
*     ..set sheath potential drop for Eirene
            if (switch%eirene_sheath_pot.eq.1) then
              ! B2 indexing for delta_sheathxb
              delta_sheathxb(BC_LIST_X(IY,IB)+1,BC_LIST_Y(IY,IB)+1) =
     .          qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
     .          bcint(te,nx,ny,iy,ib)
            endif
#endif
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            s1 = abs(pbs(bc_list_x(iy,ib),bc_list_y(iy,ib),0))
            vte = sqrt(bcint(te,nx,ny,iy,ib)/me)
            t0 = 0.0_R8
            do is = 0, ns-1
              t0 = t0+
     &         bcint(rza(-1,-1,is),nx,ny,iy,ib)*
     .         fna(bc_list_x(iy,ib),bc_list_y(iy,ib),0,0,is)
            enddo
            do is = 0, ns_ext-1
              t0 = t0+bcint(za_ext(-1,-1,is),nx,ny,iy,ib)*
     &         fa_ext(bc_list_x(iy,ib),bc_list_y(iy,ib),0,0,is)
            enddo
            fchi = qe*t0
            if(t0.lt.0.0_R8 .and.
     .       bottomiy(bc_list_x(iy,ib),bc_list_y(iy,ib)).ne.-2 .and.
     &       topiy(bc_list_x(iy,ib),bc_list_y(iy,ib)).ne.ny+1) then
              wrong_flow = .true.                                     !xpb
              if(switch%b2stbc_diagno.ge.0) then
                write(*,'(a,2i4,1x,1p,g14.7)') 'b2stbc_east:ix,iy,fchi',
     1           bc_list_x(iy,ib),bc_list_y(iy,ib),fchi
                write(hlp_frm,'(a,i3,a)') '(a,1p,',ns,'(g14.7,1x))'
                write(*,hlp_frm) 'rza : ',
     &           (rza(bc_list_x(iy,ib),bc_list_y(iy,ib),is),is=0,ns-1)
                write(*,hlp_frm) 'fna (A) : ',
     .           (fna(bc_list_x(iy,ib),bc_list_y(iy,ib),0,0,is)*qe,
     .            is=0,ns-1)
                write(*,'(a,1p,4g14.7)') 'Bx,Bz,Ud,Ua =',
     &           wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0),
     .           wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),2),
     &           ud(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain),
     .           ua(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain)
               endif
c              fchi=0.0_R8
            endif
            fche = sqrt(1/(2*pi))*qe*bcint(ne,nx,ny,iy,ib)*vte*s1*
     &       expu2(-qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
     .       bcint(te,nx,ny,iy,ib))
            if (switch%secmodel.eq.1) then
              if (ne(bc_list_x(iy,ib),bc_list_y(iy,ib))/
     .           (rza(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain)*
     .            na(bc_list_x(iy,ib),bc_list_y(iy,ib),ismain))
     .           .gt.1.1_R8) then
                write(*,*) 'Large ne contribution from impurities!'
                write(*,*) 'Kinetic s.e.c. model may not be valid!'
              endif
              call eseec0 (bcint(te,nx,ny,iy,ib),
     &         bcint(ti,nx,ny,iy,ib),bcint(ne,nx,ny,iy,ib),
     &         bcint(po,nx,ny,iy,ib),
     .         wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),0),
     &         wbbl(bc_list_x(iy,ib),bc_list_y(iy,ib),3),0,
     .         bc_list_x(iy,ib),bc_list_y(iy,ib),
     &         nx,ny,ns,ismain,seec,switch%b2stbc_diagno)               !cfs
            else
              seec = gammae
            endif
            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/bcint(te,nx,ny,iy,ib)
*    ..compute charge source
            sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),0) =
     .      sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),0)+
     &       (1.0_R8-seec)*fche-fchi+
     .       t0*po(bc_list_x(iy,ib),bc_list_y(iy,ib))
            sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),1) =
     .      sch0(bc_list_x(iy,ib),bc_list_y(iy,ib),1)-t0
#ifdef B25_EIRENE
*     ..set sheath potential drop for Eirene
            if (switch%eirene_sheath_pot.eq.1) then
              ! B2 indexing for delta_sheathxb
              delta_sheathxb(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)+1) =
     .          qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
     .          bcint(te,nx,ny,iy,ib)
            endif
#endif
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=3, NO SHEATH IN RADIAL DIR.')
        END IF
#endif
* -- BCPOT=4 -- PRESCRIBE THE VALUE OF THE POTENTIAL,
*                 WEAKLY A MIXED BOUNDARY CONDITION

      case (4)
        call xerrab ('b2stbc_phys -- BCPOT = 4 not adapted for WG')
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCPOT =  4 : weak potential ',POTPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(POTPAR(IB,2).gt.0.0_R8,
     1     'BCPOT = 4, POTPAR(,2) !> 0')
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
          srw%sch0(iCv1,0) =  POTPAR(IB,2)*geo%fcS(iFc)*potpar(iB,1)
          srw%sch0(iCv1,1) = -POTPAR(IB,2)*geo%fcS(iFc)
        enddo

* -- BCPOT=5 -- PRESCRIBE THE CURRENT FLUX DENSITY PER UNIT AREA

      case (5)
        call xerrab ('b2stbc_phys -- BCPOT = 5 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCPOT =  5 : current flux density per area ',POTPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=POTPAR(IB,1)*
     .       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=POTPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=POTPAR(IB,1)*
     .       SX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)=POTPAR(IB,1)*
     1       SX(LEFTIX(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     .          LEFTIY(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)))
            sch0(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=5, INVALID OPTION')
        END IF
#endif
* -- BCPOT=6 -- PRESCRIBE THE TOTAL CURRENT FLUX DENSITY FOR A
*                 CONSTANT POTENTIAL

      case (6)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCPOT =  6 : total current flux density ',POTPAR(IB,1),
     1     ' with const potential on ',BCCHAR(IB),boundary_location(ib)
        endif
        CALL XERRAB ('B2STBC -- BCPOT=6, NOT YET AVAILABLE')

* -- BCPOT=7 -- PRESCRIBE THE GIVEN PROFILE OF POTENTIAL

      case (7)
        if(ncall_b2stbc_phys.eq.0) then                              !srv 26.02.18 {
          write(*,'(3a)') 'specified potential on ',
     1     BCCHAR(IB), boundary_location(ib)
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          srw%sch0(iCv1,0) =  PBIG*geo%fcS(iFc)*bv_po(iCv1)
          srw%sch0(iCv1,1) = -PBIG*geo%fcS(iFc)
        enddo

* -- BCPOT=8 -- PRESCRIBE THE TOTAL ELECTRON CURRENT WITH CONSTANT
*                 FLUX DENSITY

      case (8)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCPOT =  8 : total current ',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        US = 0.0_R8
        do ibc = 1, mpg%bcCvP(ib,2)
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          US   = US + geo%fcS(iFc)
        enddo
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          srw%sch0(iCv1,0) = potpar(iB,1)/US*geo%fcS(iFc)
          srw%sch0(iCv1,1) = 0.0_R8
        enddo

* -- BCPOT=9 -- PRESCRIBE THE DECAY LENGTH FOR
*                 THE POTENTIAL

      case (9)
        call xerrab ('b2stbc_phys -- BCPOT = 9 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCPOT =  9 : potential decay length ',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(POTPAR(IB,1).ne.0.0_R8,
     1     'BCPOT = 9, POTPAR(,1) = 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -sig0(bc_list_x(ix,ib),bc_list_y(ix,ib))/POTPAR(IB,1)*
     1       SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))/qe/
     .       ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=0.0_R8
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -sig0(bc_list_x(ix,ib),bc_list_y(ix,ib))/POTPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))/qe/
     2       ne(bc_list_x(ix,ib),bc_list_y(ix,ib))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=9, NOT AVAILABLE FOR W,E')
        ENDIF
#endif
* -- BCPOT=10 -- Feedback option for core
*                NOT YET TESTED!!!!!!!!!

      case (10)
        call xerrab ('b2stbc_phys -- BCPOT = 10 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,a,a)') 'BCPOT = 10 : current feedback bc on ',
     1     BCCHAR(IB),boundary_location(ib)
          write(*,*) 'NOT YET TESTED !!!!!!!!!!'
          call xertst(lfeedback, 'LFEEDBACK must be turned on!')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if (region(bc_list_x(ix,ib),bc_list_y(ix,ib),0).eq.1) then
              sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        cbsch(0,coreregno)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            elseif (region(bc_list_x(ix,ib),bc_list_y(ix,ib),0).eq.5)
     .       then
              sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .        cbsch(0,coreregn2)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            endif
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=0.0_R8
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=10, NOT AVAILABLE FOR N,W,E')
        ENDIF
#endif
* -- BCPOT=11 -- SHEATH CONDITIONS, CURRENT CONTINUITY; FROM b2stbc_spb !srv 01.02.09 {

      case (11)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,a,1p,g14.7,a,g14.7,a,a,a)')
     &     'BCPOT = 11 : sheath potential bc from b2stbc_spb',
     &     ' with POTPAR(IB,2) =',POTPAR(IB,2),' gammae =', gammae,
     &     ' on ',BCCHAR(IB),boundary_location(ib)
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          s1   = abs(geo%fcPbs(iFc))
          cs   = sqrt(pz(iCv1)/rz(iCv1))
          vte  = sqrt(pl%te(iCv1)/me)
          t0 = 0.0_R8
          do is = 0, ns-1
            t0 = t0 + rt%rza(iCv1,is)*conpar(is,ib,1)*cs*s1*
     *        pl%na(iCv1,is)
          enddo
          fchi = qe*t0
          fche = sqrt(1/(2*pi))*qe*dv%ne(iCv1)*vte*s1*
     &      expu2(-qe*(pl%po(iCv1)- potpar(ib,2))/pl%te(iCv1))
          seec = gammae
          t0 = max(fchi,(1.0_R8-seec)*fche)*qe/pl%te(iCv1)
*    ..compute charge source
          srw%sch0(iCv1,0) = srw%sch0(iCv1,0) +
     &      (1.0_R8-seec)*fche-fchi + t0*pl%po(iCv1)
          srw%sch0(iCv1,1) = srw%sch0(iCv1,1) - t0
#ifdef B25_EIRENE
*     ..set sheath potential drop for Eirene
            if (switch%eirene_sheath_pot.eq.1) then
              ! B2 indexing for delta_sheathxb
              delta_sheathb(ifc) =
     &          qe*(pl%po(iCv1)-potpar(ib,2))/pl%te(iCv1)
            endif
#endif
cwdk Todo: seecmodel
          wrong_flow = .false.
        enddo                                                         !srv 01.02.09 }

* -- BCPOT=12 -- Test condition for south core region ! Solovyev 07.04.14   !srv 03.03.15 {
c code here is a part of code from 50000 for istyle_cur_contr_on_S_and_N.eq.1

      case (12)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,a,a,a,a)') 'BCPOT = 12 : Adding contributions ',
     &     'from diamagnetic and inertial currents ',
     &     'on ',BCCHAR(IB),boundary_location(ib)
          if(switch%istyle_cur_contr_on_S_and_N.eq.0 .or.
     &       switch%istyle_cur_contr_on_S_and_N.eq.1) then
            write(*,*) ' BCPOT=12 cannot be used when ',
     &                 'istyle_cur_contr_on_S_and_N = 0 or 1'
            call xerrab
     .       ('faulty parameter istyle_cur_contr_on_S_and_N !')
          else
           do ibc = 1, mpg%bcCvP(ib,2)
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
            call xertst(mpg%cvOnClosedSurface(iCv1),
     .       'BCPOT=12 intended for core boundaries only!')
           enddo
          endif
        endif
        if(maxval(dv%facdrift) .ne. 0.0_R8) then
*       ..compute the source
         if (bcpot_12_style.eq.-1) then
          !wdk with <fcBb(:,2)>
          do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            fchy_dia = -mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &        switch%b2stbc_phm0*
     &        pzaverage*dv%facdrift(iFc)*bzaverage*
     &        gOnedBsq(iFc,0)*geo%fcS(iFc)*geo%fcQalf(iFc,1)
            srw%sch0(iCv1,0) = srw%sch0(iCv1,0) - fchy_dia
            fchy_inertia = 0.0_R8
            do is = 0, ns-1
             if(.not.is_neutral(is))then
              fchy_inertia = fchy_inertia -
     &          mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &          switch%b2stbc_phm1*
     &          0.5e0_R8*am(is)*mp*naua2average(is)*
     &          dv%facdrift(iFc)*bzaverage*gOnedBsq(iFc,0)*
     &          geo%fcS(iFc)*geo%fcQalf(iFc,1)
             endif
            enddo
            srw%sch0(iCv1,0) = srw%sch0(iCv1,0) - fchy_inertia
          enddo
         elseif (bcpot_12_style.eq.0) then
          do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            !wdk now uses fcBb(iFc,2) instead of <fcBb(:,2)>
            fchy_dia = -mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &        switch%b2stbc_phm0*
     &        pzaverage*dv%facdrift(iFc)*geo%fcBb(iFc,2)*
     &        gOnedBsq(iFc,0)*geo%fcS(iFc)*geo%fcQalf(iFc,1)
            srw%sch0(iCv1,0) = srw%sch0(iCv1,0) - fchy_dia
            fchy_inertia = 0.0_R8
            do is = 0, ns-1
             if(.not.is_neutral(is))then
              fchy_inertia = fchy_inertia -
     &          mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &          switch%b2stbc_phm1*
     &          0.5e0_R8*am(is)*mp*naua2average(is)*
     &          dv%facdrift(iFc)*geo%fcBb(iFc,2)*gOnedBsq(iFc,0)*
     &          geo%fcS(iFc)*geo%fcQalf(iFc,1)
             endif
            enddo
            srw%sch0(iCv1,0) = srw%sch0(iCv1,0) - fchy_inertia
          enddo
         elseif (bcpot_12_style.eq.1) then
          !wdk experimental - use fchdia / fchinert as computed
          do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
            !wdk now uses fcBb(iFc,2) instead of <fcBb(:,2)>
            fchy_dia = mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &        (dv%fchdia(iFc,0)  + dv%fchdia(iFc,1))
            srw%sch0(iCv1,0) = srw%sch0(iCv1,0) - fchy_dia
            fchy_inertia = mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &        (dv%fchinert(iFc,0) + dv%fchinert(iFc,1))
            srw%sch0(iCv1,0) = srw%sch0(iCv1,0) - fchy_inertia
          enddo
         elseif (bcpot_12_style.eq.2) then
          !wdk experimental - find constant potential (with
          !    perturbation taken from flux tube just inside
          !    domain) such that total current through core
          !    is sum of dia + inert currents
*     ..compute average potential on core boundary
          t0 = 0.0e0_R8
          t1 = 0.0e0_R8
          t2 = 0.0e0_R8
          totflux = 0.0e0_R8
          totpar  = 0.0e0_R8
          do ib0 = 1, nbc
           ! identify core boundaries; assume all core boundaries
           ! form single flux surface
           ! implicit assumption at the moment: first ring in the
           ! core is a 'structured grid' ring with quad cells
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
           if (mpg%cvOnClosedSurface(iCv1)) then
            if (bcpot(ib0).ne.12)
     .        call xerrab ('BCPOT=12: incompatible core BCs')
            do ibc0 = 1, mpg%bcCvP(ib0,2)                               ! loop over number of cells in the boundary
             iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                 ! number of the guard cell
             iCv2 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,2)                 ! number of the domain cell
             iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
             t0 = t0 + geo%cvVol(iCv2)
             t1 = t1 + geo%cvVol(iCv2)*pl%po(iCv1)
             t2 = t2 + geo%cvVol(iCv2)*pl%po(iCv2)
             totflux = totflux - mpg%bcFcOr(mpg%bcCvP(ib0,1)+ibc0-1)*
     &         (dv%fch(iFc,0) + dv%fch(iFc,1))
             totpar = totpar - mpg%bcFcOr(mpg%bcCvP(ib0,1)+ibc0-1)*
     &         (dv%fchdia(iFc,0) + dv%fchdia(iFc,1) +
     &          dv%fchinert(iFc,0) + dv%fchinert(iFc,1))
            enddo
           endif
          enddo
          t1 = t1/t0
          t2 = t2/t0
          if (potpar(ib,2).gt.0.0_R8) t3 = potpar(ib,2)
          if (potpar(ib,2).le.0.0_R8) t3 = 1.0e-2_R8
          t4 = (1.0_R8-t3*(totflux-totpar)/
     &         max(abs(totflux),totpar,1.0e-10_R8))*t1
*       ..impose potential with perturbation
          do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                      ! number of the domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
            srw%sch0(iCv1,0) =  PBIG*geo%fcS(iFc)*
     &        (pl%po(iCv2) + t4 - t2)
            srw%sch0(iCv1,1) = -PBIG*geo%fcS(iFc)
          enddo
         else
          call xerrab ('wrong value bcpot_12_style')
         endif
        endif

* -- BCPOT=13 -- Test condition for south private & north regions ! Solovyev 07.04.14
c code here is a part of code from 50000 for istyle_cur_contr_on_S_and_N.eq.1

      case (13)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,a,a,a,a)') 'Adding contributions ',
     &     'from all currents ',
     &     'on ',BCCHAR(IB),boundary_location(ib)
          if(switch%istyle_cur_contr_on_S_and_N.eq.0 .or.
     &       switch%istyle_cur_contr_on_S_and_N.eq.1) then
            write(*,*) ' BCPOT=13 cannot be used when ',
     &                 'istyle_cur_contr_on_S_and_N = 0 or 1'
            call xerrab
     .       ('faulty parameter istyle_cur_contr_on_S_and_N !')
          else
           do ibc = 1, mpg%bcCvP(ib,2)
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
            call xertst(.not.mpg%cvOnClosedSurface(iCv1),
     .       'BCPOT=13 not intended for core boundaries!')
           enddo
          endif
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          srw%sch0(iCv1,0) = srw%sch0(iCv1,0) -
     &      mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &     (dv%fchdia(iFc,0)    + dv%fchdia(iFc,1)    +
     &      dv%fchin(iFc,0)     + dv%fchin(iFc,1)     +
     &      dv%fchvispar(iFc,0) + dv%fchvispar(iFc,1) +
     &      dv%fchvisq(iFc,0)   + dv%fchvisq(iFc,1)   +
     &      dv%fchviskt(iFc,0)  + dv%fchviskt(iFc,1)   +
     &      dv%fchvisper(iFc,0) + dv%fchvisper(iFc,1) +
     &      dv%fchinert(iFc,0)  + dv%fchinert(iFc,1))
        enddo                                                           !srv 03.03.15 }

* -- BCPOT=14 -- PRESCRIBE THE TOTAL CURRENT WITHOUT IMPOSING A
*                 SPECIFIC PROFILE
cwdk    Intended use: core boundary, total PF boundary,...
cwdk    For testing purposes, use with care!

      case (14)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'BCPOT = 14 : specified radial sheath potential ',
     1     POTPAR(IB,1),' on ',BCCHAR(IB),boundary_location(ib)
        endif
        US = 0.0_R8
        t0 = 0.0_R8
        do ibc = 1, mpg%bcCvP(ib,2)
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          US   = US + geo%fcS(iFc)
          t0   = t0 + mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &      (dv%fch(iFc,0) + dv%fch(iFc,1))
        enddo
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          srw%sch0(iCv1,0) = - mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)*
     &      (dv%fch(iFc,0) + dv%fch(iFc,1)) +
     &      (t0-POTPAR(IB,1))/US*geo%fcS(iFc)
          srw%sch0(iCv1,1) = 0.0_R8
        enddo

* -- BCPOT=15 -- Feedback boundary condition on total current with constant potential

      case (15)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'zero net anomalous current on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
*     ..compute average potential on boundary
        t0 = 0.0e0_R8
        t1 = 0.0e0_R8
        totflux = 0.0e0_R8
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          t0 = t0 + geo%fcS(iFc)
          t1 = t1 + geo%fcS(iFc)*pl%po(iCv1)
          totflux = totflux - (dv%fch(iFc,0) + dv%fch(iFc,1))*
     &      mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)
        enddo
        t1 = t1/t0
        totpar = potpar(ib,1)
        if (potpar(ib,2).gt.0.0_R8) t3 = potpar(ib,2)
        if (potpar(ib,2).le.0.0_R8) t3 = 1.0e-2_R8
        t2 = (1.0_R8-t3*(totflux-totpar)/
     &       max(abs(totflux),abs(totpar),1.0_R8))*t1
        if (switch%b2stbc_diagno.ge.2) then
          write(*,'(a,1p,4(1x,e15.6))')
     &      'bcpot=15: curr_flux, desired_flux, old/new potential, ',
     &      totflux,totpar,t1,t2
        endif
*     ..impose potential
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          srw%sch0(iCv1,0) =  PBIG*geo%fcS(iFc)*t2
          srw%sch0(iCv1,1) = -PBIG*geo%fcS(iFc)
        enddo

* -- BCPOT=16 -- Feedback boundary condition with constant potential    !Solovyev 22.03.14 !srv 29.09.15

      case (16)
        call xerrab ('b2stbc_phys -- BCPOT = 16 not adapted for WG')
#ifdef WG_TODO
        if (ncall_b2stbc_phys.eq.0) then
          write(*,'(a,a,a)')
     .     'BCPOT = 16 : electric potential feedback by current on ',
     1      BCCHAR(IB),boundary_location(ib)
! add smth about #ifdef ASTRA here
c Init variables for BCPOT=16                                  { Solovyev 01.05.14
          cur_delta = 0.0_R8              ! square of delta between desired and calculated current
          prev_cur_delta = 0.0_R8         ! previous value of upper value
          po_step = 0.0_R8                ! step of approximation (set after 1st iteration)
          coeff_16 = 1.0_R8
c************************************************************* } Solovyev 01.05.14
          scurpar = 0.0_R8
          do ibw=ib,nbc
            if (BCCHAR(ibw).eq.'S'.and.BCPOT(ibw).eq.BCPOT(ib)) then
              scurpar = scurpar + potpar(ibw,1)
            endif
          enddo
          write (*,*) 'SCURPAR = ', scurpar
          POS = 0.0_R8                                                  !srv 29.09.15
          VS  = 0.0_R8
          DO IBW=IB,NBC
            if(BCCHAR(IBW).EQ.'S'.and.BCPOT(IBW).EQ.BCPOT(IB)) THEN
              DO IX=1, BC_LIST_SIZE(IBW)
                POS = POS + PO(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))*
     *                     VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                VS = VS + VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
              ENDDO
            endif
          ENDDO
          POS=POS/VS                                                    !srv 29.09.15 }
          po_prev = POS
          write (*,*) 'POS at 0 = ', POS
        endif
        if (.not.once_bcpot) then
          if (BCCHAR(ib).eq.'S') then
            totflux = 0.0_R8
            do ibw=ib,nbc
              if (BCCHAR(ibw).eq.'S'.and.BCPOT(ibw).eq.BCPOT(ib)) then
                do ix=1, bc_list_size(ibw)
                  totflux = totflux +
     &              fch(topix(bc_list_x(ix,ibw),bc_list_y(ix,ibw)),
     &                  topiy(bc_list_x(ix,ibw),bc_list_y(ix,ibw)),1,1)
                enddo
              endif
            enddo
c The implementation is based on one step of Newton method
            cur_delta = totflux-scurpar                                 ! { Solovyev 03.05.14
            if (ncall_b2stbc_phys.eq.0) then
              write (*, *) '# Newton'
              prev_cur_delta = cur_delta         ! previous value of upper value
            else
              if (ncall_b2stbc_phys.eq.1) then                          ! Solovyev 23.05.14
                po_curr = 0.0_R8                                        !srv 29.09.15
                VS  = 0.0_R8
                DO IBW=IB,NBC
                  if(BCCHAR(IBW).EQ.'S'.and.BCENE(IBW).EQ.BCENE(IB))THEN
                    DO IX=1, BC_LIST_SIZE(IBW)
                      po_curr=po_curr+
     &                   PO(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
     &                 *VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                      VS = VS + VOL(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))
                    ENDDO
                  endif
                ENDDO
                po_curr=po_curr/VS                                      !srv 29.09.15 }
              endif                                                     !Solovyev 23.05.14
              if (mod(ncall_b2stbc_phys,bcpot_16_step).eq.0) then       !Solovyev 23.05.14
                if (abs(cur_delta-prev_cur_delta).gt.1.0e-20_R8) then   !{ Solovyev 24.05.14
                  coeff_16 = (po_prev - po_curr)/
     &              (cur_delta - prev_cur_delta)
                endif
                po_step = coeff_16*cur_delta                            !} Solovyev 24.05.14
                POS = po_curr + po_step
                prev_cur_delta = cur_delta
                po_prev = po_curr
                po_curr = POS
              endif                                                     ! Solovyev 23.05.14
            endif                                                       ! } Solovyev 03.05.14
c end of implementation
            write (*, *) 'ncall, totflux, POS',
     &        ncall_b2stbc_phys,totflux,POS                              !Solovyev 29.04.14 - test output to b2fdebug
            do ibw=ib,nbc
              if (BCCHAR(ibw).eq.'S'.and.BCPOT(ibw).eq.BCPOT(ib)) then
                do ix=1, BC_LIST_SIZE(IBW)
                  sch0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),0) =
     &                 PBIG*POS*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW)) !Solovyev 27.04.14
                  sch0(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW),1) =
     &                -PBIG*SY(BC_LIST_X(IX,IBW),BC_LIST_Y(IX,IBW))     !Solovyev 27.04.14
                enddo
              endif
            enddo
            once_bcpot = .true.
          else
            call xerrab('B2STBC -- BCPOT=16, NOT AVAILABLE FOR N,W,E')
          endif
        endif                                                               !srv 29.09.15 }
#endif
* -- BCPOT=17 -- Feedback boundary condition on total current, imposing potential,
*                with same profile as neighbouring flux tube

      case (17)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'imposed current on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
*     ..compute average potential on boundary
        t0 = 0.0e0_R8
        t1 = 0.0e0_R8
        t2 = 0.0e0_R8
        totflux = 0.0e0_R8
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                      ! number of the plasma cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          t0 = t0 + geo%fcS(iFc)
          t1 = t1 + geo%fcS(iFc)*pl%po(iCv1)
          t1 = t2 + geo%fcS(iFc)*pl%po(iCv2)
          totflux = totflux - (dv%fch(iFc,0) + dv%fch(iFc,1))*
     &      mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)
        enddo
        t1 = t1/t0
        t2 = t2/t0
        totpar = potpar(ib,1)
        if (potpar(ib,2).gt.0.0_R8) t3 = potpar(ib,2)
        if (potpar(ib,2).le.0.0_R8) t3 = 1.0e-2_R8
        t4 = (1.0_R8-t3*(totflux-totpar)/
     &       max(abs(totflux),abs(totpar),1.0_R8))*t1
        if (switch%b2stbc_diagno.ge.2) then
          write(*,'(a,1p,4(1x,e15.6))')
     &      'bcpot=17: curr_flux, desired_flux, old/new potential, ',
     &      totflux,totpar,t1,t4
        endif
*     ..impose potential
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                      ! number of the domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          srw%sch0(iCv1,0) =  PBIG*geo%fcS(iFc)*
     .      (pl%po(iCv2) + t4 - t2)
          srw%sch0(iCv1,1) = -PBIG*geo%fcS(iFc)
        enddo

* -- BCPOT=18 -- Nothing implemented

      case (18)
        CALL XERRAB ('B2STBC -- BCPOT=18, NOT AVAILABLE FOR S,N,W,E')

* -- BCPOT=19 -- Nothing implemented

      case (19)
        CALL XERRAB ('B2STBC -- BCPOT=19, NOT AVAILABLE FOR S,N,W,E')

* -- BCPOT=20 -- Nothing implemented

      case (20)
        CALL XERRAB ('B2STBC -- BCPOT=20, NOT AVAILABLE FOR S,N,W,E')

* -- BCPOT=21 -- constant potential scaled by potential on the ring bc_type21_ref away

      case (21)
        call xerrab ('b2stbc_phys -- BCPOT = 21 not adapted for WG')
#ifdef WG_TODO
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     .     'BCPOT = 21 : fixed potential (modified by nearby ring) ',
     1     POTPAR(IB,1),' with ',POTPAR(IB,2),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(potpar(ib,2).gt.0.0_R8,
     1     'faulty parameter potpar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          poave = 0.0_R8
          neave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            poave=poave+po(bc_type21_ix,bc_type21_iy)*
     1       ne(bc_type21_ix,bc_type21_iy)*
     2       vol(bc_type21_ix,bc_type21_iy)
            neave=neave+ne(bc_type21_ix,bc_type21_iy)*
     1       vol(bc_type21_ix,bc_type21_iy)
          ENDDO
          poave=poave/neave
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=topix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=topiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*POTPAR(IB,1)*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*
     .       ((1.0_R8-potpar(ib,2))+
     .           potpar(ib,2)*po(bc_type21_ix,bc_type21_iy)/poave)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=
     .      -VBIG*SY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
          ENDDO
        ELSE IF(BCCHAR(IB).EQ.'N') THEN
          poave = 0.0_R8
          neave = 0.0_R8
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            poave=poave+po(bc_type21_ix,bc_type21_iy)*
     1       ne(bc_type21_ix,bc_type21_iy)*
     2       vol(bc_type21_ix,bc_type21_iy)
            neave=neave+ne(bc_type21_ix,bc_type21_iy)*
     1       vol(bc_type21_ix,bc_type21_iy)
          ENDDO
          poave=poave/neave
          DO IX = 1, BC_LIST_SIZE(IB)
            bc_type21_ix=bc_list_x(ix,ib)
            bc_type21_iy=bc_list_y(ix,ib)
            do bc_type21_i=1,bc_type21_ref      ! how many steps away from the boundary?
              bc_type21_ixx=bottomix(bc_type21_ix,bc_type21_iy)
              bc_type21_iy=bottomiy(bc_type21_ix,bc_type21_iy)
              bc_type21_ix=bc_type21_ixx
              call xertst (-2.lt.bc_type21_ix.and.bc_type21_ix.lt.nx+1,
     1         'stepped off x-grid in type 21 bc')
              call xertst (-2.lt.bc_type21_iy.and.bc_type21_iy.lt.ny+1,
     1         'stepped off y-grid in type 21 bc')
            enddo
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)=
     .       VBIG*POTPAR(IB,1)*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))*
     .       ((1.0_R8-potpar(ib,2))+
     .           potpar(ib,2)*po(bc_type21_ix,bc_type21_iy)/poave)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)=-VBIG*
     1       SY(BOTTOMIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     .          BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=21, NOT AVAILABLE FOR W,E')
        ENDIF
#endif

* -- BCPOT=22 -- Prescribe the value of the potential, with perturbation taken
*                from first flux tube in the domain
      case(22)
        if(ncall_b2stbc_phys.eq.0) then
          do ibc = 1, mpg%bcCvP(ib,2)
           iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
           call xertst(mpg%cvOnClosedSurface(iCv1),
     .      'BCPOT=22 intended for core boundaries only!')
          enddo
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified potential ',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
*     ..compute average potential on core boundary
        t0 = 0.0e0_R8
        t1 = 0.0e0_R8
        do ib0 = 1, nbc
         ! identify core boundaries; assume all core boundaries
         ! form single flux surface
         ! implicit assumption at the moment: first ring in the
         ! core is a 'structured grid' ring with quad cells
         iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
         if (mpg%cvOnClosedSurface(iCv1)) then
          do ibc0 = 1, mpg%bcCvP(ib0,2)                                 ! loop over number of cells in the boundary
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                   ! number of the guard cell
           iCv2 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,2)                   ! number of the domain cell
           iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
           t0 = t0 + geo%cvVol(iCv2)
           t1 = t1 + geo%cvVol(iCv2)*pl%po(iCv2)
          enddo
         endif
        enddo
        t1 = t1/t0
*     ..impose potential with perturbation
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                      ! number of the domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          srw%sch0(iCv1,0) =  PBIG*geo%fcS(iFc)*
     &      (pl%po(iCv2) + potpar(iB,1) - t1)
          srw%sch0(iCv1,1) = -PBIG*geo%fcS(iFc)
        enddo

* -- BCPOT=23 -- Feedback on the total current, by imposing a flux surface
*                averaged potential, with poloidal perturbation taken from flux
*                tube just inside domain

        case(23)
        if(ncall_b2stbc_phys.eq.0) then
          do ibc = 1, mpg%bcCvP(ib,2)
           iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
           call xertst(mpg%cvOnClosedSurface(iCv1),
     .      'BCPOT=23 intended for core boundaries only!')
          enddo
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified potential ',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
*     ..compute average potential on core boundary
        t0 = 0.0e0_R8
        t1 = 0.0e0_R8
        t2 = 0.0e0_R8
        totflux = 0.0e0_R8
        totpar  = 0.0e0_R8
        do ib0 = 1, nbc
         ! identify core boundaries; assume all core boundaries
         ! form single flux surface
         ! implicit assumption at the moment: first ring in the
         ! core is a 'structured grid' ring with quad cells
         iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
         if (mpg%cvOnClosedSurface(iCv1)) then
          if (bcpot(ib0).ne.23)
     .      call xerrab ('BCPOT=23: incompatible core BCs')
          do ibc0 = 1, mpg%bcCvP(ib0,2)                                 ! loop over number of cells in the boundary
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                   ! number of the guard cell
           iCv2 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,2)                   ! number of the domain cell
           iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
           t0 = t0 + geo%cvVol(iCv2)
           t1 = t1 + geo%cvVol(iCv2)*pl%po(iCv1)
           t2 = t2 + geo%cvVol(iCv2)*pl%po(iCv2)
           totflux = totflux - (dv%fch(iFc,0) + dv%fch(iFc,1))*
     &       mpg%bcFcOr(mpg%bcCvP(ib0,1)+ibc0-1)
          enddo
          totpar = totpar + potpar(ib0,1)
         endif
        enddo
        t1 = t1/t0
        t2 = t2/t0
        if (potpar(ib,2).gt.0.0_R8) t3 = potpar(ib,2)
        if (potpar(ib,2).le.0.0_R8) t3 = 1.0e-2_R8
        t4 = (1.0_R8-t3*(totflux-totpar)/
     &         max(abs(totflux),abs(totpar),1.0e-10_R8))*t1
*     ..impose potential with perturbation
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                      ! number of the domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          srw%sch0(iCv1,0) =  PBIG*geo%fcS(iFc)*
     &      (pl%po(iCv2) + t4 - t2)
          srw%sch0(iCv1,1) = -PBIG*geo%fcS(iFc)
        enddo

      case default
        write(*,*) BCPOT(IB)
        call xerrab ('BCPOT out of bounds')
      end select

#ifdef WG_TODO
      if(istyle_cur_contr_on_S_and_N.eq.0) then                       !srv 01.02.09
!xpb
!xpb  Adding fchy_dia contributions on south and north boundaries
!xpb  following "sv 01.07.99" coding (see B2.5-style implementation)
!xpb  Special attention is paid to the case where there is a discontinuity in B
!xpb
        IF (BCPOT(IB).NE.0.AND.BCPOT(IB).NE.5.AND.BCPOT(IB).NE.8) THEN
          if (ncall_b2stbc_phys.eq.0) write(*,'(1x,a,i2,a,i2)')
     &          'Disabled special treatment of fchy_dia on boundary ',
     &          ib,' because potential is imposed by BCPOT =',BCPOT(IB)
        ELSE IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if(maxval(facdrift) .ne. 0.0_R8) then
              if(on_closed_surface(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .         then
                if(periodic_bc.ne.1 .or.
     .             rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)).eq.
     .             bc_list_x(ix,ib)+1 .or. nncut.eq.0) then
                  fchy_dia = -praverage(1)*
     .             facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))
     &             *switch%b2stbc_fchy_dia                               !sv 01.07.99
     &             *(OnedBsqc(rightix(topix(bc_list_x(ix,ib),
     .                                      bc_list_y(ix,ib)),
     1                                topiy(bc_list_x(ix,ib),
     .                                      bc_list_y(ix,ib))),
     &                        rightiy(topix(bc_list_x(ix,ib),
     .                                      bc_list_y(ix,ib)),
     1                                topiy(bc_list_x(ix,ib),
     .                                      bc_list_y(ix,ib))))
     &              -OnedBsqc(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                        topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))) !sv
                else
                  if(ncall_b2stbc_phys.eq.0  )write(*,*)
     >             'b2stbc: special treatment for cell ix,iy = ',
     .              bc_list_x(ix,ib),bc_list_y(ix,ib)
                  fchy_dia = -praverage(1)
     .             *facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))
     &             *switch%b2stbc_fchy_dia  !sv 01.07.99
     &             *2.0_R8*(1.0_R8/
     1              wbbv(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)**2
     &             -OnedBsqc(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                       topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))  !xpb
                endif
              else
                if(periodic_bc.ne.1 .or. nncut.eq.0 .or.
     &     rightix(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))).eq.
     .             bc_list_x(ix,ib)+1)
     >           then
                  fchy_dia = -bcint(wdia,nx,ny,ix,ib)
     .             *(vol(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     .               vol(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .             /(gs(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                  topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2)+
     .               gs(bc_list_x(ix,ib),bc_list_y(ix,ib),2))
     .      *(OnedBsqc(rightix(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     .                 rightiy(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .             -OnedBsqc(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                       topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .             *facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))
     .             *switch%b2stbc_fchy_dia
                else
                  if (ncall_b2stbc_phys.eq.0) write(*,*)
     .             'b2stbc: special treatment for cell ix,iy = ',
     .              bc_list_x(ix,ib),bc_list_y(ix,ib)
                  fchy_dia = -bcint(wdia,nx,ny,ix,ib)
     .             *(vol(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))+
     .               vol(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .             /(gs(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                  topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2)+
     .               gs(bc_list_x(ix,ib),bc_list_y(ix,ib),2))*2.0_R8
     .     *(1.0_R8/wbbv(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)**2
     .     -OnedBsqc(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .             *facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))
     .             *switch%b2stbc_fchy_dia
                endif
              endif
            else
              if(switch%b2stbc_neoclassical .ne. 0.0_R8 .and.
     &           on_closed_surface(bc_list_x(ix,ib),
     .                             bc_list_y(ix,ib))) then
                if(periodic_bc.ne.1 .or. nncut.eq.0 .or.
     &     rightix(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .             topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))).eq.
     .             bc_list_x(ix,ib)+1)
     >           then
                  fchy_dia = -praverage(1)
     .    *(OnedBsqc(rightix(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                       topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))),
     .               rightiy(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                       topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .              -OnedBsqc(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                        topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .              *switch%b2stbc_neoclassical
                else
                  fchy_dia = -praverage(1)*2.0_R8*
     .     *(1.0_R8/wbbv(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                   topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),3)**2
     .              -OnedBsqc(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                        topiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .              *switch%b2stbc_neoclassical
                endif
              else
                fchy_dia = 0.0_R8
              endif
            endif
            if(.not.on_closed_surface(bc_list_x(ix,ib),
     .                                bc_list_y(ix,ib)) .and.
     &         switch%b2stbc_fchy_dia_coreonly.eq.1) fchy_dia = 0.0_R8
            sch0(bc_list_x(ix,ib),bc_list_y(ix,ib),0) =
     .      sch0(bc_list_x(ix,ib),bc_list_y(ix,ib),0) + fchy_dia
          END DO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if(maxval(facdrift) .ne. 0.0_R8) then
              if(on_closed_surface(bc_list_x(ix,ib),
     .                             bc_list_y(ix,ib))) then
                fchy_dia = -praverage(2)
     .           *facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))
     &           *switch%b2stbc_fchy_dia     !sv 01.07.99
     &           *(OnedBsqc(rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                      rightiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))-
     1             OnedBsqc(bc_list_x(ix,ib),bc_list_y(ix,ib)))   !sv
              else
                fchy_dia = -bcint(wdia,nx,ny,ix,ib)
     .           *(vol(bc_list_x(ix,ib),bc_list_y(ix,ib))+
     .             vol(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                 bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib))))
     .           /(gs(bc_list_x(ix,ib),bc_list_y(ix,ib),2)+
     .             gs(bottomix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                bottomiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),2))
     .           *(OnedBsqc(rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                      rightiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .            -OnedBsqc(bc_list_x(ix,ib),bc_list_y(ix,ib)))*
     .             facdrift(bc_list_x(ix,ib),bc_list_y(ix,ib))
     &            *switch%b2stbc_fchy_dia
              endif
            else
              if(switch%b2stbc_neoclassical .ne. 0.0_R8 .and.
     &           on_closed_surface(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     >         then
                fchy_dia = -praverage(2)
     .           *(OnedBsqc(rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .                      rightiy(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .            -OnedBsqc(bc_list_x(ix,ib),bc_list_y(ix,ib)))
     .           *switch%b2stbc_neoclassical
              else
                fchy_dia = 0.0_R8
              endif
            endif
            if(.not.on_closed_surface(bc_list_x(ix,ib),bc_list_y(ix,ib))
     .                      .and.
     &         switch%b2stbc_fchy_dia_coreonly.eq.1) fchy_dia = 0.0_R8
            sch0(bc_list_x(ix,ib),bc_list_y(ix,ib),0) =
     .      sch0(bc_list_x(ix,ib),bc_list_y(ix,ib),0) - fchy_dia
          END DO
        ENDIF
!xpb
      elseif(istyle_cur_contr_on_S_and_N.eq.1) then                   !srv 01.02.09 { 27.03.14

        IF (BCPOT(IB).NE.0.AND.BCPOT(IB).NE.5.AND.BCPOT(IB).NE.8) THEN
          if (ncall_b2stbc_phys.eq.0) write(*,'(a,i2,a,i2)')
     &          'Disabled special treatment of fchy_dia on boundary ',
     &          ib,' because potential there is imposed by BCPOT=',
     &          BCPOT(IB)
        ELSE IF(BCCHAR(IB).EQ.'S') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            if(on_closed_surface(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)))
     >        then
*        for inner surfaces
              if(maxval(facdrift) .ne. 0.0_R8) then
                fchy_dia =                                             !srv 01.07.99
     &          -pzaverage
     &          *facdrift(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*phm0      ! *signmf        !srv 24.03.99 01.07.99 21.04.00 18.01.02 2
     &          *bzhzaverage                                           !lk 20.11.07 }
     &          *(OnedBsqc(
     &            rightix(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))),
     &            rightiy(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))-
     &            OnedBsqc(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                     topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))   !srv 01.07.99
                sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) + fchy_dia
                fchy_inertia = 0.0_R8                                   !srv 25.06.03
                do is = 0, ns-1                                         !srv 25.06.03
                  if(.not.is_neutral(is))then                           !srv 25.06.03
                    fchy_inertia = fchy_inertia                         !srv 07.04.00 25.06.03
     &               -0.5e0_R8*am(is)*mp*
     &                facdrift(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))*phm1  ! *signmf     !srv 07.04.00 21.04.00 21.10.02
     &               *naaverage(is)                                     !lk 20.11.07  !srv 07.04.00
     &               *ua2average(is)                                    !srv 07.04.00 15.02.02
     &               *bzhzaverage                                       !lk 20.11.07
     &               *(OnedBsqc(
     &            rightix(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))),
     &            rightiy(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))-
     &           OnedBsqc(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                    topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))))    !srv 07.04.00
                  endif
                enddo
                sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &           sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0)+fchy_inertia
              endif
            else
*        for PR region
              sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &          sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) +
     &           fchdia(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1)+
     &           fchin(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1)+
     &        fchvispar(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1)+
     &          fchvisq(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1)+
     &        fchvisper(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1)+
     &         fchinert(topix(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                  topiy(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,1)
            endif
          END DO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) =
     &        sch0(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0) -
     &             fchdia(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1)-
     &             fchin(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1)-
     &             fchvispar(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1)-
     &             fchvisq(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1)-
     &             fchvisper(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1)-
     &             fchinert(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,1)
          END DO
        ENDIF
      else if (switch%pot_eq.eq.1.and.
     &   BCPOT(IB).eq.0.or.BCPOT(IB).eq.5.or.BCPOT(IB).eq.8) then !srv 27.03.14
        if(ncall_b2stbc_phys.eq.0) then
          IF (BCCHAR(IB).EQ.'S') THEN
            IF (on_closed_surface(bc_list_x(1,ib),
     &                            bc_list_y(1,ib))) THEN
              write(0,'(a/a)')
     &        'Core boundary potential B.C. does not '//
     &        'ensure drift current continuity',
     &        'Use BCPOT=12 to impose it.'                           !srv 03.03.15
            ELSE
              write(0,'(a/a)')
     &        'South boundary potential B.C. does not '//
     &        'ensure drift current continuity',
     &        'Use BCPOT=13 to impose it.'                           !srv 03.03.15
            ENDIF
          ELSE IF (BCCHAR(IB).EQ.'N') THEN
            write(0,'(a/a)')
     &        'North boundary potential B.C. does not '//
     &        'ensure drift current continuity',
     &        'Use BCPOT=13 to impose it.'                           !srv 03.03.15
           ENDIF
         endif
        endif                                                        !srv 01.02.09 }
#endif
* -- BOUNDARY CONDITIONS FOR TURBULENT KINETIC ENERGY EQUATION --

        select case (BCENK(IB))

        case (0)

* -- BCENK=1 -- PRESCRIBE THE VALUE OF KT

        case (1)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified kt ',ENKPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(ENKPAR(IB,1).gt.0.0_R8,
     1     'BCENK = 1, ENKPAR(,1) !> 0')
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          srw%skt0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*enkpar(iB,1)*EV
          srw%skt0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
        enddo

* -- BCENK=2 -- PRESCRIBE THE GRADIENT OF THE TURB. KIN. ENERGY
        case (2)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'turb. kin. en. gradient ',ENKPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
cwdk    Note: specifies the gradient in the outward normal direction
cwdk    Tmp: in "connector" direction - ok for orthogonal grids
        do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                    ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
          iVx1 = mpg%fcVx(iFc,1)
          iVx2 = mpg%fcVx(iFc,2)
          t1 = intvertex_s(iVx1,nCv,nVx,mpg,geo%vxVol,pl%kt)
          t2 = intvertex_s(iVx2,nCv,nVx,mpg,geo%vxVol,pl%kt)
          srw%skt0(iCv1,0)=  VBIG*NE00*geo%fcS(iFc)*
     .       (EV*enkpar(ib,1)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2))*
     .        geo%fcQgam(iFc,0) - cor9*(t2 - t1)*geo%fcQgam(iFc,1)*
     .        (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     .        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1) + pl%kt(iCv2))
          srw%skt0(iCv1,1)= -VBIG*NE00*geo%fcS(iFc)
        enddo

* -- BCENK=3 -- not implemented
        case (3)
        CALL XERRAB ('B2STBC -- BCENK=3, NOT AVAILABLE')

* -- BCENK=4 -- not implemented
        case (4)
        CALL XERRAB ('B2STBC -- BCENK=4, NOT AVAILABLE')

* -- BCENK=5 -- not implemented
        case (5)
        CALL XERRAB ('B2STBC -- BCENK=5, NOT AVAILABLE')

* -- BCENK=6 -- Prescribe the total flux of kt, following the
*               profile of fhi + fhe

        case (6)
        if(ncall_b2stbc_phys.eq.0) then
          do ibc = 1, mpg%bcCvP(ib,2)
           iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
           call xertst(mpg%cvOnClosedSurface(iCv1),
     .      'BCENK=6 intended for core boundaries only!')
          enddo
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'turb. kin. energy flux ',ENKPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
*     ..compute total ion+electron heat flux, and total kt flux
        totflux = 0.0e0_R8
        totpar  = 0.0e0_R8
        do ib0 = 1, nbc
         ! identify core boundaries; assume all core boundaries
         ! form single flux surface
         ! implicit assumption at the moment: first ring in the
         ! core is a 'structured grid' ring with quad cells
         iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
         if (mpg%cvOnClosedSurface(iCv1)) then
          if (bcenk(ib0).ne.6)
     .      call xerrab ('BCENK=6: incompatible core BCs')
          do ibc0 = 1, mpg%bcCvP(ib0,2)                                 ! loop over number of cells in the boundary
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                   ! number of the guard cell
           iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
           totflux = totflux -
     &      (dv%fhe(iFc,0) + dv%fhe(iFc,1) +
     &       dv%fhi(iFc,0) + dv%fhi(iFc,1))*
     &       mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)
          enddo
          totpar = totpar + enkpar(ib0,1)
         endif
        enddo
*     ..impose the flux
        write (*,*) 'BCENK=6, totflux = ',totflux,
     &                     ', totpar = ',totpar
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          cs   = sqrt(pz(iCv1)/rz(iCv1))
          srw%skt0(iCv1,0) = -totpar/totflux*
     &      (dv%fhe(iFc,0) + dv%fhe(iFc,1) +
     &       dv%fhi(iFc,0) + dv%fhi(iFc,1))*
     &       mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)
          srw%skt0(iCv1,1) = 0.0_R8
        enddo

* -- BCENK=7 -- Prescribe the total flux of kt, following the
*               profile of fna

        case (7)
        if(ncall_b2stbc_phys.eq.0) then
          do ibc = 1, mpg%bcCvP(ib,2)
           iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
           call xertst(mpg%cvOnClosedSurface(iCv1),
     .      'BCENK=7 intended for core boundaries only!')
          enddo
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'turb. kin. energy flux ',ENKPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
*     ..compute total particle flux, and total kt flux
        totflux = 0.0e0_R8
        totpar  = 0.0e0_R8
        do ib0 = 1, nbc
         ! identify core boundaries; assume all core boundaries
         ! form single flux surface
         ! implicit assumption at the moment: first ring in the
         ! core is a 'structured grid' ring with quad cells
         iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
         if (mpg%cvOnClosedSurface(iCv1)) then
          if (bcenk(ib0).ne.7)
     .      call xerrab ('BCENK=7: incompatible core BCs')
          do ibc0 = 1, mpg%bcCvP(ib0,2)                                 ! loop over number of cells in the boundary
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                   ! number of the guard cell
           iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
           totflux = totflux -
     &      (dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))*
     &       mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)
          enddo
          totpar = totpar + enkpar(ib0,1)
         endif
        enddo
*     ..impose the flux
        write (*,*) 'BCENK=7, totflux = ',totflux,
     &                     ', totpar = ',totpar
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          cs   = sqrt(pz(iCv1)/rz(iCv1))
          srw%skt0(iCv1,0) = -totpar/totflux*
     &      (dv%fna(iFc,0,ismain) + dv%fna(iFc,1,ismain))*
     &       mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)
          srw%skt0(iCv1,1) = 0.0_R8
        enddo

* -- BCENK=8 -- PRESCRIBE THE TOTAL FLUX OF KT WITH CONSTANT
*               FLUX DENSITY

        case (8)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'total turb. kin. energy flux ',ENKPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        US = 0.0_R8
        do ibc = 1, mpg%bcCvP(ib,2)
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          US   = US + geo%fcS(iFc)
        enddo
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          srw%skt0(iCv1,0) = enkpar(iB,1)/US*geo%fcS(iFc)
          srw%skt0(iCv1,1) = 0.0_R8
        enddo

* -- BCENK=9 -- not implemented
        case (9)
        CALL XERRAB ('B2STBC -- BCENK=9, NOT AVAILABLE')

* -- BCENK=10 -- not implemented
        case (10)
        CALL XERRAB ('B2STBC -- BCENK=10, NOT AVAILABLE')

* -- BCENK=11 -- not implemented
        case (11)
        CALL XERRAB ('B2STBC -- BCENK=11, NOT AVAILABLE')

* -- BCENK=12 -- not implemented
        case (12)
        CALL XERRAB ('B2STBC -- BCENK=12, NOT AVAILABLE')

* -- BCENK=13 -- not implemented
        case (13)
        CALL XERRAB ('B2STBC -- BCENK=13, NOT AVAILABLE')

* -- BCENK=14 -- Leakage option for turbulent kinetic energy
        case (14)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'turb. kin. energy leakage ',ENKPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          cs   = sqrt(pz(iCv1)/rz(iCv1))
          srw%skt0(iCv1,0) = 0.0_R8
          srw%skt0(iCv1,3) = geo%fcS(iFc)*cs*enkpar(ib,1)
        enddo

* -- BCENK=15 -- Sheath tranmission factor
        case (15)
        if(ncall_b2stbc_phys.eq.0) then
            write(*,*) 'BCENK=15 - assuming BCCON=14 for all ions'
#ifdef WG_TODO
          IF(BCCHAR(IB).EQ.'W') THEN
            BCCON14_W(IB) = .TRUE.  ! IF BCCON=14 HAS BEEN IMPOSED AT WEST BOUNDARY IB FOR ALL IONS
            DO IS = 0, NS-1
              IF(.NOT.IS_NEUTRAL(IS)) THEN
                IF(BCCON(IS,IB).NE.14) THEN
                  BCCON14_W(IB) = .FALSE.
                ENDIF
              ENDIF
            ENDDO
            if(BCCON14_W(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 HAS ',
     &      'BEEN IMPOSED AT WEST BOUNDARY', IB,' FOR ALL IONS'
            if(.not.BCCON14_W(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 ',
     &      'HAS NOT BEEN IMPOSED AT WEST BOUNDARY', IB,' FOR ALL IONS'
          ENDIF
          IF(BCCHAR(IB).EQ.'E') THEN
            BCCON14_E(IB) = .TRUE.  ! IF BCCON=14 HAS BEEN IMPOSED AT WEST BOUNDARY IB FOR ALL IONS
            DO IS = 0, NS-1
              IF(.NOT.IS_NEUTRAL(IS)) THEN
                IF(BCCON(IS,IB).NE.14) THEN
                  BCCON14_E(IB) = .FALSE.
                ENDIF
              ENDIF
            ENDDO
            if(BCCON14_E(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 HAS ',
     &       'BEEN IMPOSED AT EAST BOUNDARY', IB,' FOR ALL IONS'
            if(.not.BCCON14_E(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 ',
     &       'HAS NOT BEEN IMPOSED AT EAST BOUNDARY', IB,' FOR ALL IONS'
          ENDIF
          call xertst(ENKPAR(IB,1).gt.0.0_R8,
     1       'BCENK = 15, ENKPAR(I,,1) <= 0 not allowed!')              !srv 20.09.17
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     &     'from b2stbc_spb for ion temperature with ENIPAR(IB,1) =',
     &     ENKPAR(IB,1),' on ', BCCHAR(IB),boundary_location(ib)        !srv 20.09.17
#endif
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          cs   = sqrt(pz(iCv1)/rz(iCv1))
          t0 = 0.0_R8
          t1 = 0.0_R8
          do is = 0, ns-1
            if(.not.is_neutral(is)) then
              t0 = t0 + conpar(is,ib,1)*pl%na(iCv1,is)
              if (conpar(is,ib,2).gt.0.0_R8)
     &         t1 = t1 + co%dna0(iCv1,is)/conpar(is,ib,2)*pl%na(iCv1,is)
            endif
          enddo
          srw%skt0(iCv1,1) = srw%skt0(iCv1,1) -
     &      enkpar(ib,1)*t0*cs*abs(geo%fcPbs(iFc))*
     &      (1.0_R8 + switch%b2tfhi_fflokt)  -
     &      enkpar(ib,2)*t1*geo%fcS(iFc)
cwdk    Todo: corrections for BCCON.ne.14
        enddo

* -- BCENK=16 -- Prescribe the value of the turb. kin. energy, with poloidal
*                perturbation taken from flux tube just inside domain

        case (16)
        if(ncall_b2stbc_phys.eq.0) then
          do ibc = 1, mpg%bcCvP(ib,2)
           iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
           call xertst(mpg%cvOnClosedSurface(iCv1),
     .      'BCENK=16 intended for core boundaries only!')
          enddo
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified turb. kin. energy + perturbation',ENKPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENKPAR(IB,1).gt.0.0_R8,
     1     'BCENK = 16, ENKPAR(,1) !> 0')
        endif
*     ..compute average kt on core boundary
        t0 = 0.0e0_R8
        t1 = 0.0e0_R8
        do ib0 = 1, nbc
         ! identify core boundaries; assume all core boundaries
         ! form single flux surface
         ! implicit assumption at the moment: first ring in the
         ! core is a 'structured grid' ring with quad cells
         iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
         if (mpg%cvOnClosedSurface(iCv1)) then
          do ibc0 = 1, mpg%bcCvP(ib0,2)                                 ! loop over number of cells in the boundary
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                   ! number of the guard cell
           iCv2 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,2)                   ! number of the domain cell
           iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
           t0 = t0 + geo%cvVol(iCv2)
           t1 = t1 + geo%cvVol(iCv2)*pl%kt(iCv2)
          enddo
         endif
        enddo
        t1 = t1/t0
*     ..impose kt with perturbation
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                      ! number of the domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          srw%skt0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*
     .      (pl%kt(iCv2) + enkpar(iB,1)*EV - t1)
          srw%skt0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
        enddo


* -- BCENK=17 -- Feedback on turb. kin. energy flux, by imposing a flux surface
*                averaged kt, with poloidal perturbation taken from flux
*                tube just inside domain

        case (17)
        if(ncall_b2stbc_phys.eq.0) then
          do ibc = 1, mpg%bcCvP(ib,2)
           iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
           call xertst(mpg%cvOnClosedSurface(iCv1),
     .      'BCENK=17 intended for core boundaries only!')
          enddo
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified turb. kin. energy flux ',ENKPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENKPAR(IB,1).gt.0.0_R8,
     1     'BCENK = 17, ENKPAR(,1) !> 0')
        endif
*     ..compute average kt on core boundary
        t0 = 0.0e0_R8
        t1 = 0.0e0_R8
        t2 = 0.0e0_R8
        totflux = 0.0e0_R8
        totpar  = 0.0e0_R8
        do ib0 = 1, nbc
         ! identify core boundaries; assume all core boundaries
         ! form single flux surface
         ! implicit assumption at the moment: first ring in the
         ! core is a 'structured grid' ring with quad cells
         iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
         if (mpg%cvOnClosedSurface(iCv1)) then
          if (bcenk(ib0).ne.17)
     .      call xerrab ('BCENK=17: incompatible core BCs')
          do ibc0 = 1, mpg%bcCvP(ib0,2)                                 ! loop over number of cells in the boundary
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                   ! number of the guard cell
           iCv2 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,2)                   ! number of the domain cell
           iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
           t0 = t0 + geo%cvVol(iCv2)*dv%ni(iCv1,0)
           t1 = t1 + geo%cvVol(iCv2)*dv%ni(iCv1,0)*pl%kt(iCv1)
           t2 = t2 + geo%cvVol(iCv2)*dv%ni(iCv2,0)*pl%kt(iCv2)
           totflux = totflux - (dv%fkt(iFc,0) + dv%fkt(iFc,1))*
     &       mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1)
          enddo
          totpar = totpar + enkpar(ib0,1)
         endif
        enddo
        t1 = t1/t0
        t2 = t2/t0
        if (enkpar(ib,2).gt.0.0_R8) t3 = enkpar(ib,2)
        if (enkpar(ib,2).le.0.0_R8) t3 = 1.0e-2_R8
        t4 = (1.0_R8-t3*(totflux-totpar)/
     &         max(abs(totflux),totpar))*t1
        if (bcenk_17_style.eq.0) then
*       ..impose kt with perturbation
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                    ! number of the domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%skt0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*
     .        (pl%kt(iCv2) + t4 - t2)
            srw%skt0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
          enddo
        else
*       ..impose kt w/o perturbation
          do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
            iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
            iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                    ! number of the domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
            srw%skt0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*t2
            srw%skt0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
          enddo
        endif

        case default
          write(*,*) BCENK(IB)
          call xerrab ('BCENK out of bounds')

        end select

* -- BOUNDARY CONDITIONS FOR TURBULENT ENSTROPHY EQUATION --

        select case (BCENZ(IB))

        case (0)

* -- BCENZ=1 -- PRESCRIBE THE VALUE OF ZT

        case (1)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified zt ',ENZPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
          call xertst(ENZPAR(IB,1).ge.0.0_R8,
     1     'BCENZ = 1, ENZPAR(,1) !> 0')
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          srw%szt0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*enzpar(iB,1)*EV
          srw%szt0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
        enddo

* -- BCENZ=2 -- PRESCRIBE THE GRADIENT OF THE TURB. ENSTROPHY
        case (2)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'turb. enstr. gradient ',ENZPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                   ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                    ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                    ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                          ! number of the guard cell face
          iVx1 = mpg%fcVx(iFc,1)
          iVx2 = mpg%fcVx(iFc,2)
          t1 = intvertex_s(iVx1,nCv,nVx,mpg,geo%vxVol,pl%zt)
          t2 = intvertex_s(iVx2,nCv,nVx,mpg,geo%vxVol,pl%zt)
          srw%szt0(iCv1,0)=  VBIG*NE00*geo%fcS(iFc)*
     .       (EV*enzpar(ib,1)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2))*
     .        geo%fcQgam(iFc,0) - cor9*(t2 - t1)*geo%fcQgam(iFc,1)*
     .        (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     .        mpg%bcFcOr(mpg%bcCvP(ib,1)+ibc-1) + pl%zt(iCv2))
          srw%szt0(iCv1,1)= -VBIG*NE00*geo%fcS(iFc)
        enddo

* -- BCENZ=3 -- not implemented
        case (3)
        CALL XERRAB ('B2STBC -- BCENZ=3, NOT AVAILABLE')

* -- BCENZ=4 -- not implemented
        case (4)
        CALL XERRAB ('B2STBC -- BCENZ=4, NOT AVAILABLE')

* -- BCENZ=5 -- not implemented
        case (5)
        CALL XERRAB ('B2STBC -- BCENZ=5, NOT AVAILABLE')

* -- BCENZ=6 -- Prescribe the total flux of kt, following the
*               profile of fhi + fhe
        case (6)
        CALL XERRAB ('B2STBC -- BCENZ=6, NOT AVAILABLE')

* -- BCENZ=7 -- Prescribe the total flux of kt, following the
*               profile of fna
        case (7)
        CALL XERRAB ('B2STBC -- BCENZ=7, NOT AVAILABLE')

* -- BCENZ=8 -- PRESCRIBE THE TOTAL FLUX OF KT WITH CONSTANT
*               FLUX DENSITY
        case (8)
        CALL XERRAB ('B2STBC -- BCENZ=8, NOT AVAILABLE')

* -- BCENZ=9 -- not implemented
        case (9)
        CALL XERRAB ('B2STBC -- BCENZ=9, NOT AVAILABLE')

* -- BCENZ=10 -- not implemented
        case (10)
        CALL XERRAB ('B2STBC -- BCENZ=10, NOT AVAILABLE')

* -- BCENZ=11 -- not implemented
        case (11)
        CALL XERRAB ('B2STBC -- BCENZ=11, NOT AVAILABLE')

* -- BCENZ=12 -- not implemented
        case (12)
        CALL XERRAB ('B2STBC -- BCENZ=12, NOT AVAILABLE')

* -- BCENZ=13 -- not implemented
        case (13)
        CALL XERRAB ('B2STBC -- BCENZ=13, NOT AVAILABLE')

* -- BCENZ=14 -- Leakage option for turbulent enstrophy
        case (14)
        if(ncall_b2stbc_phys.eq.0) then
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'turb. enstrophy leakage ',ENZPAR(IB,1),' on ',
     1     BCCHAR(IB),boundary_location(ib)
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          cs   = sqrt(pz(iCv1)/rz(iCv1))
          srw%szt0(iCv1,0) = 0.0_R8
          srw%szt0(iCv1,3) = geo%fcS(iFc)*cs*enzpar(ib,1)
        enddo

* -- BCENZ=15 -- Sheath tranmission factor
        case (15)
        if(ncall_b2stbc_phys.eq.0) then
            write(*,*) 'BCENZ=15 - assuming BCCON=14 for all ions'
#ifdef WG_TODO
          IF(BCCHAR(IB).EQ.'W') THEN
            BCCON14_W(IB) = .TRUE.  ! IF BCCON=14 HAS BEEN IMPOSED AT WEST BOUNDARY IB FOR ALL IONS
            DO IS = 0, NS-1
              IF(.NOT.IS_NEUTRAL(IS)) THEN
                IF(BCCON(IS,IB).NE.14) THEN
                  BCCON14_W(IB) = .FALSE.
                ENDIF
              ENDIF
            ENDDO
            if(BCCON14_W(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 HAS ',
     &      'BEEN IMPOSED AT WEST BOUNDARY', IB,' FOR ALL IONS'
            if(.not.BCCON14_W(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 ',
     &      'HAS NOT BEEN IMPOSED AT WEST BOUNDARY', IB,' FOR ALL IONS'
          ENDIF
          IF(BCCHAR(IB).EQ.'E') THEN
            BCCON14_E(IB) = .TRUE.  ! IF BCCON=14 HAS BEEN IMPOSED AT WEST BOUNDARY IB FOR ALL IONS
            DO IS = 0, NS-1
              IF(.NOT.IS_NEUTRAL(IS)) THEN
                IF(BCCON(IS,IB).NE.14) THEN
                  BCCON14_E(IB) = .FALSE.
                ENDIF
              ENDIF
            ENDDO
            if(BCCON14_E(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 HAS ',
     &       'BEEN IMPOSED AT EAST BOUNDARY', IB,' FOR ALL IONS'
            if(.not.BCCON14_E(IB)) write(*,"(a,a,i3,a)") 'BCCON=14 ',
     &       'HAS NOT BEEN IMPOSED AT EAST BOUNDARY', IB,' FOR ALL IONS'
          ENDIF
          call xertst(ENKPAR(IB,1).gt.0.0_R8,
     1       'BCENZ = 15, ENZPAR(I,,1) <= 0 not allowed!')              !srv 20.09.17
          write(*,'(a,1p,g14.7,a,g14.7,a,a,a)')
     &     'from b2stbc_spb for ion temperature with ENIPAR(IB,1) =',
     &     ENZPAR(IB,1),' on ', BCCHAR(IB),boundary_location(ib)        !srv 20.09.17
#endif
        endif
        do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
          cs   = sqrt(pz(iCv1)/rz(iCv1))
          t0 = 0.0_R8
          t1 = 0.0_R8
          do is = 0, ns-1
            if(.not.is_neutral(is)) then
              t0 = t0 + conpar(is,ib,1)*pl%na(iCv1,is)
              if (conpar(is,ib,2).gt.0.0_R8)
     &         t1 = t1 + co%dna0(iCv1,is)/conpar(is,ib,2)*pl%na(iCv1,is)
            endif
          enddo
          srw%szt0(iCv1,1) = srw%szt0(iCv1,1) -
     &      enzpar(ib,1)*t0*cs*abs(geo%fcPbs(iFc))*
     &      (1.0_R8 + switch%b2tfhi_fflozt)  -
     &      enzpar(ib,2)*t1*geo%fcS(iFc)
cwdk    Todo: corrections for BCCON.ne.14
        enddo

* -- BCENZ=16 -- Prescribe the value of the turb. enstropy, with poloidal
*                perturbation taken from flux tube just inside domain
        case (16)
        if(ncall_b2stbc_phys.eq.0) then
          do ibc = 1, mpg%bcCvP(ib,2)
           iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)
           call xertst(mpg%cvOnClosedSurface(iCv1),
     .      'BCENZ=16 intended for core boundaries only!')
          enddo
          write(*,'(a,1p,g14.7,a,a,a)')
     .     'specified turb. enstrophy + perturbation',ENZPAR(IB,1),
     1     ' on ',BCCHAR(IB),boundary_location(ib)
          call xertst(ENZPAR(IB,1).gt.0.0_R8,
     1     'BCENZ = 16, ENZPAR(,1) !> 0')
        endif
*     ..compute average zt on core boundary
        t0 = 0.0e0_R8
        t1 = 0.0e0_R8
        do ib0 = 1, nbc
         ! identify core boundaries; assume all core boundaries
         ! form single flux surface
         ! implicit assumption at the moment: first ring in the
         ! core is a 'structured grid' ring with quad cells
         iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1),1)
         if (mpg%cvOnClosedSurface(iCv1)) then
          do ibc0 = 1, mpg%bcCvP(ib0,2)                                 ! loop over number of cells in the boundary
           iCv1 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,1)                   ! number of the guard cell
           iCv2 = mpg%bcCv(mpg%bcCvP(ib0,1)+ibc0-1,2)                   ! number of the domain cell
           iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                           ! number of the guard cell face
           t0 = t0 + geo%cvVol(iCv2)
           t1 = t1 + geo%cvVol(iCv2)*pl%zt(iCv2)
          enddo
         endif
        enddo
        t1 = t1/t0
*     ..impose zt with perturbation
        do ibc = 1, mpg%bcCvP(ib,2)                                     ! loop over number of cells in the boundary
          iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                      ! number of the guard cell
          iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                      ! number of the domain cell
          iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                            ! number of the guard cell face
          srw%szt0(iCv1,0) =  VBIG*NE00*geo%fcS(iFc)*
     .      (pl%zt(iCv2) + enzpar(iB,1)*EV - t1)
          srw%szt0(iCv1,1) = -VBIG*NE00*geo%fcS(iFc)
        enddo

        case default
          write(*,*) BCENZ(IB)
          call xerrab ('BCENZ out of bounds')

        end select


      enddo
#ifdef WG_TODO
* ---------------------------------------------------------------------
*  INTERIOR BOUNDARY AT DEAD REGIONS
*  ---------------------------------------------------------------------
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            DO IS=0, NS-1
              sna0(IX,IY,0,IS)= UBIG*NIISO(IS)
              sna0(IX,IY,1,IS)=-UBIG
              smo0(ix,iy,0,is)= ubig*ne00*am(is)*0.0_R8*
     .         sx(min(nx-1,ix),iy)
              smo0(ix,iy,1,is)=-ubig*ne00*am(is)*
     .         sx(min(nx-1,ix),iy)
            end do
            she0(IX,IY,0)= UBIG*NE(IX,IY)*TEISO*EV
            she0(IX,IY,1)=-UBIG*NE(IX,IY)
            shi0(IX,IY,0)= UBIG*NE(IX,IY)*TIISO*EV
            shi0(IX,IY,1)=-UBIG*NE(IX,IY)
            sch0(IX,IY,0)= UBIG*PHIISO
            sch0(IX,IY,1)=-UBIG
          end do
        end do
      end do
#endif
      if (switch%b2npmo_iout.eq.1 .or. switch%iout_b2wdat.eq.4) then               !srv 11.09.09 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        wrk0=srw%smo0(:,0,is)+
     &       srw%smo0(:,1,is)*pl%ua(:,is)+
     &       srw%smo0(:,2,is)*pl%na(:,is)*mp*am(is)+
     &       srw%smo0(:,3,is)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
        call my_out_us(70,nCv,0,wrk0,'b2stbc_phys_smo'//chns)
       enddo
      endif
      if ((switch%b2npco_iout.eq.1) .or.
     &                           switch%iout_b2wdat.eq.4) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        wrk0=srw%sna0(:,0,is)+srw%sna0(:,1,is)*pl%na(:,is)
        call my_out_us(70,nCv,0,wrk0,'b2stbc_phys_sna'//chns)
       enddo
      endif
      if (switch%b2npht_iout.eq.1 .or.
     &                           switch%iout_b2wdat.eq.4) then
       wrk0=srw%she0(:,0)+srw%she0(:,1)*pl%te+
     &      srw%she0(:,2)*dv%ne + srw%she0(:,3)*pl%te*dv%ne
       call my_out_us(70,nCv,0,wrk0,'b2stbc_phys_she')
       wrk0=srw%shi0(:,0)+srw%shi0(:,1)*pl%ti+
     &      srw%shi0(:,2)*dv%ni(:,0)+srw%shi0(:,3)*pl%ti*dv%ni(:,0)
       call my_out_us(70,nCv,0,wrk0,'b2stbc_phys_shi')
      endif
      if ((switch%b2npp7_iout.eq.1 .or. switch%b2nppo_iout.eq.1) .or.
     &                           switch%iout_b2wdat.eq.4) then
       wrk0=srw%sch0(:,0)+srw%sch0(:,1)*pl%po+
     &      srw%sch0(:,2)*dv%ne+srw%sch0(:,3)*pl%po*dv%ne
       call my_out_us(70,nCv,0,wrk0,'b2stbc_phys_sch')
      endif

      if (switch%b2stbc_iout.eq.1) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        do k=0,3
          write (chk,'(i1)') k
          call my_out_us(70,nCv,0,srw%smo0(1,k,is),
     &               'b2stbc_phys_smo'//chk//chns)
        enddo
        do k=0,1
          write (chk,'(i1)') k
          call my_out_us(70,nCv,0,srw%sna0(1,k,is),
     &               'b2stbc_phys_sna'//chk//chns)
        enddo
       enddo
       do k=0,3
         write (chk,'(i1)') k
         call my_out_us(70,nCv,0,srw%sch0(1,k),'b2stbc_phys_sch'//chk)
         call my_out_us(70,nCv,0,srw%she0(1,k),'b2stbc_phys_she'//chk)
         call my_out_us(70,nCv,0,srw%shi0(1,k),'b2stbc_phys_shi'//chk)
       enddo
      endif                                                           !srv 11.09.09 }

* ..return
      if(ncall_b2stbc_phys.eq.0 .or. main_call)
     &  ncall_b2stbc_phys = ncall_b2stbc_phys+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbc_phys

      contains

      function sheath_g (t0)
      use b2mod_types
      implicit none
      real (kind=R8) :: sheath_g, t0
      intrinsic sqrt, erf, exp

      sheath_g = 1.0_R8 + sqrt(pi)*t0*(1+erf(t0))*exp(t0**2)

      return
      end function sheath_g

      function pit(iCv)
      implicit none
      integer, intent(in) :: iCv
      real (kind=R8) :: pit

      pit = abs(geo%cvBb(iCv,0)/geo%cvBb(iCv,3))

      return
      end function pit

      end subroutine b2stbc_phys
