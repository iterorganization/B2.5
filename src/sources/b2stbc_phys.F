
*-----------------------------------------------------------------------
*.specification
      subroutine b2stbc_phys (nx, ny, ns, ismain, 
     1  BoRiS, facdrift, fac_ExB,
     1  vol, crx, hx, hy, hz, qz, qc, gs, pbs, pbshz, bb, na, ua,     !srv 24.08.08
     1  uadia, vedia, vadia, veecrb, vaecrb,
     &  te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc, fna, fna_mdf, !srv 28.10.02
     1  fnaPSch, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,      !srv 13.10.06
     &  fhe, fhi, fhe_mdf, fhi_mdf, fna_fcor, fhePSch, fhiPSch,
     &  fch, fne, fni,                                                !srv 13.10.06
     &  fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,!srv 21.10.02
     2  cdna, chce, chve, chci, chvi, calf,                           !srv 25.09.09
     &  sna0, smo0, she0, shi0, sch0, wrong_flow, main_call, diagno,
     &  b2stbc_fchy_dia, b2stbc_neoclassical, praverage,
     1  b2stbc_fchy_dia_coreonly, secmodel, sound_model, area_fix)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_anomalous_transport
      use b2mod_boundary_namelist
      use b2mod_constants
      use b2mod_geo_corner                                               !xpb
      use b2mod_time
      use b2mod_b2cmpa
      use b2mod_b2cmpb
#ifdef PAR_OPT
	use b2mod_par_opt
#endif
#ifdef ASTRA
	use b2mod_astra_to_b2                                          !srv 28.07.11
#endif
      implicit none
#ifdef USE_EIRENE
#     include "DIMENSIONS.F"
#     include "COUPLE/KOPPLDIM.F"
#     include "COUPLE/EIRDIAG.F"
#endif
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, diagno
      real (kind=R8) ::
     *  BoRiS, facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny), 
     &  vol(-1:nx,-1:ny), crx(-1:nx,-1:ny,0:3),
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), hz(-1:nx,-1:ny),           !srv 24.08.08
     &  qz(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1),
     &  pbshz(-1:nx,-1:ny,0:1),
     1  bb(-1:nx,-1:ny,0:3), gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     1  uadia(-1:nx,-1:ny,0:1,0:ns-1), 
     &  vedia(-1:nx,-1:ny,0:1), vadia(-1:nx,-1:ny,0:1,0:ns-1),
     2  veecrb(-1:nx,-1:ny,0:1), vaecrb(-1:nx,-1:ny,0:1,0:ns-1), 
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), 
     &  kinrgy(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1), fna_32(-1:nx,-1:ny,0:1,0:ns-1), 
     *  fna_mdf(-1:nx,-1:ny,0:1,0:ns-1), 
     &  fna_52(-1:nx,-1:ny,0:1,0:ns-1), fni_32(-1:nx,-1:ny,0:1), 
     &  fni_52(-1:nx,-1:ny,0:1), fne_32(-1:nx,-1:ny,0:1), 
     &  fne_52(-1:nx,-1:ny,0:1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fhe_mdf(-1:nx,-1:ny,0:1),
     &  fhi_mdf(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     &  fch_p(-1:nx,-1:ny,0:1), fchvispar(-1:nx,-1:ny,0:1),           !srv 11.10.02 21.10.02 {
     &  fchvisper(-1:nx,-1:ny,0:1), fchin(-1:nx,-1:ny,0:1),
     &  fchvisq(-1:nx,-1:ny,0:1), fchdia(-1:nx,-1:ny,0:1),
     &  fchinert(-1:nx,-1:ny,0:1),                                    !srv 21.10.02 }
     &  fne(-1:nx,-1:ny,0:1), fni(-1:nx,-1:ny,0:1),
     &  fna_fcor(-1:nx,-1:ny,0:1,0:ns-1),                            !srv 25.02.00
     &  floe_noc(-1:nx,-1:ny,0:1), floi_noc(-1:nx,-1:ny,0:1),        !srv 28.10.02
     &  fhePSch(-1:nx,-1:ny,0:1), fhiPSch(-1:nx,-1:ny,0:1),          !srv 13.10.06
     &  fnaPSch(-1:nx,-1:ny,0:1,0:ns-1)                              !srv 13.10.06
      real (kind=R8) ::
     1  chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),
     2  chci(-1:nx,-1:ny,0:1), chvi(-1:nx,-1:ny,0:1),
     3  calf(-1:nx,-1:ny,0:1), cdna(-1:nx,-1:ny,0:1,0:ns-1)          !srv 25.09.09
      real (kind=R8) ::
     &  b2stbc_fchy_dia, b2stbc_neoclassical, praverage(2)
      logical main_call                                                  !xpb
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3)
      integer b2stbc_fchy_dia_coreonly, secmodel, sound_model, area_fix
      logical wrong_flow                                                 !xpb
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBC_PHYS provides linearised expressions for the sources of
*     density, momentum, electron heat, atom heat, and charge, due to
*     the boundary conditions.
*     Here the boundary conditions are expressed by use of the old B2
*     namelist format.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, ib, ireg, coreregno, coreregn2, 
     1  bc_type13_ref, bc_type13_ix, bc_type13_iy, 
     1  bc_type13_i, bc_type13_ixx, bc_type16_ref, 
     1  bc_type16_ix, bc_type16_iy, bc_type16_i, 
     1  bc_type16_ixx,bc_type20_ref, bc_type20_ix, 
     1  bc_type20_iy, bc_type20_i, bc_type20_ixx,
     1  bc_type16_kinetic_energy, ibw
     	integer       :: num_inner, k, isb, isa, ism                      !srv 26.11.02
      integer, save :: istyle_average_pr = 0, istyle_fchi = 0           !srv 08.08.02 14.02.05
      integer, save :: iform_corr_dia = 0                               !srv 25.07.05
	integer, save :: mdf_fnb = 0, mdf_fhe = 0, mdf_fhi = 0, iout = 0, !srv 13.10.06
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0, iout_b2npmo =0
	integer, save :: istyle_cur_contr_on_S_and_N = 0                  !srv 01.02.09
	integer, save :: bcene_16_style =0, bceni_16_style = 0            !srv 16.02.12
	integer, save :: bcene_15_style =0                                !srv 19.04.13
	integer, save :: pot_eq = 1                                       !srv 27.02.13
	logical once_bccon_16, once_bccon_19, once_bcene_16, once_bceni_16 !srv 16.02.12
      real (kind=R8) , save :: cbc = 0.0e0_R8, bc_ref = 0.01            !srv 14.10.99
      real (kind=R8) , save :: bc_ref_te = 0.01, bc_ref_ti=0.01         !IYS 07.03.2013 11.04.13
      real (kind=R8) , save :: senepar, senipar, sconpar                !srv 16.02.12
      real (kind=R8) ::
     *  cutll, cutlo, vte, s0, s1, s0hz, s1hz, fche, fchi, t0, t1, t2,      !srv 15.11.07
     1  t3, totflux, weight, naave, volave, teave, neave, tiave, niave,
     2  vwest, veast, totfluxn, totfluxi
     &  ,pzaverage, pzsum, volsum                                     !srv 24.03.99 18.01.02
     &  ,ua2sum, ua2average(0:ns-1)                                   !srv 07.04.00 15.02.02
     &  ,uasum, uaaverage(0:ns-1)                                     !srv 15.02.02
     &  ,nesum, neaverage, tisum, tiaverage, tesum, teaverage         !srv 15.12.00 25.01.02
     &  ,nasum, naaverage(0:ns-1), bbsum, bbaverage                   !srv 22.01.02 18.10.04
      real (kind=R8) :: bzhzsum, bzhzaverage                          !lk 20.11.07 
      real (kind=R8) ::
     *  pr(-1:nx,-1:ny), wdia(-1:nx,-1:ny), uu(-1:nx,-1:ny,0:ns-1),
     *  sx(-1:nx,-1:ny), sxhz(-1:nx,-1:ny), pz(-1:nx,-1:ny),         !srv 24.08.08
     &  sna0_no_mdf(-1:nx,-1:ny,0:1,0:ns-1)                          !srv 23.09.08
      real (kind=R8), save ::                                        !srv 16.02.12
     *  NUEOUT, UBIG, VBIG, UP00, EPSLON, NE00, S2, US, FS, VS,
     *  NAS, TES, TIS, FFAC, DFS
      real (kind=R8) , save :: phm0 = 1.0e0_R8, phm1 = 1.0e0_R8,
     &                         phm2 = 1.0e0_R8,
     &                         corr_flux   = 0.0e0_R8                 !srv 25.07.05
 	character chns*2, chk*1                                        !srv 11.07.05
      logical done
      real (kind=R8) ::
     *  ro(-1:nx,-1:ny), cs, csb, uuout, upout, zbar, zarg, uutest,  !srv 01.02.09
     *  fchy_dia, fchy_inertia, seec                                 !srv 07.04.00
      integer elm_count, isfeedback, fix_fch_in_fhe_sheath
      integer, save :: fix_sx_in_fhe_sheath = 1                       !srv 08.05.09
      real (kind=R8) :: nesepm, nesepm_overshoot,
     1 bc_type13_norm, bc_type13_fac
      logical catch_up
      save ncall, elm_count, coreregno, coreregn2,  bc_type13_ref,
     1  bc_type16_ref, bc_type20_ref, bc_type13_norm, bc_type13_fac,
     *  isfeedback, nesepm, nesepm_overshoot, fix_fch_in_fhe_sheath
*   ..procedures
      intrinsic min, max, abs, sqrt, exp, log, dim, mod
      real (kind=R8) ::
     *  machsfr, expu2, hy1, sy, ud, bcint, syhz                   !srv 24.08.08
      external subini, subend, xertst, ipgeti, machsfr, bcint, eseec0
      external ipgetr, b2xvsg, b2xvff, b2xvps
      expu2(t0) = dim(exp(min(20.0_R8,max(t0,cutll))),cutlo)
*     (expu2 is the exponential function with underflow set to 0 and
*     large values limited at a somewhat arbitrary level.)
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
      sy(ix,iy)=gs(topix(ix,iy),topiy(ix,iy),1)
      syhz(ix,iy)=gs(topix(ix,iy),topiy(ix,iy),1)*                    !srv 24.08.08
     &            hz(ix,iy)
!     &            hz(topix(ix,iy),topiy(ix,iy))
      ud(ix,iy,is)=uadia(ix,iy,0,is)
*   ..initialisation
      data ncall/0/, coreregno/1/, coreregn2/4/, bc_type13_ref/1/, 
     1 bc_type16_ref/1/, bc_type20_ref/1/,
     * isfeedback/0/, nesepm/0.0_R8/, nesepm_overshoot/0.0_R8/,
     5 fix_fch_in_fhe_sheath/1/, bc_type16_kinetic_energy/0/,
     & bc_type13_norm/1.0e15_R8/, bc_type13_fac/1.0_R8/
#ifdef WANT_THIS
      real (kind=R8) ::
     .  Tguess, fy, ehy, fylimt, prefactor,
     .  fyt, gam0, gam1, t2, tlimit, cone0, coni0, flo, a,               !srv 15.02.08
     .  cum_flo_vadia, cum_flo_vaecrb, flo_vedia, flo_veecrb
      real (kind=R8) :: hybr, harm, pit, upwind
      integer ibnd, no_current, conduction_only, no_hybr_e, no_hybr_i
      data no_current/0/, conduction_only/0/, no_hybr_e/0/, no_hybr_i/0/
      save no_current, conduction_only, no_hybr_e, no_hybr_i
      hybr(flo,t1) = max(t1,abs(flo/2.),sqrt(t1**2+flo**2/4.))
      harm(t1,t2)= 2.0_R8*t1*t2/(t1+t2)
      pit(ix,iy)=abs(bb(ix,iy,0)/bb(ix,iy,3))
      upwind(flo,t1,t2) = max(flo,0.0_R8)*t1+min(flo,0.0_R8)*t2
#endif
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbc_phys')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2stbc_phm0', phm0)                                  !srv 17.01.07
       call ipgetr ('b2stbc_phm1', phm1)                                  !srv 17.01.07
       call ipgetr ('b2stbc_phm2', phm2)                                  !srv 17.01.07
        call ipgeti ('b2stbc_coreregno', coreregno)
        call xertst (cbirso.le.coreregno.and.
     1   coreregno.le.cbirso+cbnrso-1,
     2   'faulty internal parameter coreregno')
        if (nncut.ge.2) then
          call ipgeti ('b2stbc_coreregn2', coreregn2)
          call xertst (cbirso.le.coreregn2.and.
     1     coreregn2.le.cbirso+cbnrso-1,
     2     'faulty internal parameter coreregn2')
        endif
#ifdef WANT_THIS
        call ipgeti ('b2tfhe_no_current', no_current)
        call ipgeti ('b2tfhe_conduction_only', conduction_only)
        call ipgeti ('b2tfhe_no_hybr', no_hybr_e)
        call ipgeti ('b2tfhi_no_hybr', no_hybr_i)
#endif
        call ipgeti ('b2stbc_type13_ref', bc_type13_ref)
        call ipgeti ('b2stbc_type16_ref', bc_type16_ref)
        call ipgeti ('b2stbc_type20_ref', bc_type20_ref)
	  call ipgetr ('b2stbc_type13_norm', bc_type13_norm)
        call ipgetr ('b2stbc_type13_fac', bc_type13_fac)
        call ipgeti ('b2stbc_fix_fch_in_fhe_sheath',
     1   fix_fch_in_fhe_sheath)
        call ipgeti ('b2stbc_fix_sx_in_fhe_sheath',
     1   fix_sx_in_fhe_sheath)                                       !srv 08.05.09
        call ipgeti ('b2stbc_type16_kinetic_energy',
     1   bc_type16_kinetic_energy)
        call ipgetr ('b2stbc_cbc', cbc)                              !srv 14.10.99
        call ipgeti ('b2stbc_average_pr', istyle_average_pr)         !srv 08.08.02
        call ipgeti ('b2stbc_istyle_fchi', istyle_fchi)              !srv 14.02.05
        call ipgetr ('b2stbc_corr_flux_in_core', corr_flux)          !srv 25.07.05
        call ipgeti ('b2stbc_iform_corr_dia', iform_corr_dia)        !srv 25.07.05
        call ipgeti ('b2stbc_iout', iout)                            !srv 11.07.05
        call ipgeti ('b2tfnb_mdf', mdf_fnb)                          !srv 13.10.06
        call ipgeti ('b2tfhe_mdf', mdf_fhe)                          !srv 13.10.06
        call ipgeti ('b2tfhi_mdf', mdf_fhi)                          !srv 13.10.06
        call ipgeti ('b2stbc_istyle_cur_contr_on_S_and_N', 
     &                       istyle_cur_contr_on_S_and_N)            !srv 01.02.09
        call ipgeti ('b2npc9_iout', iout_b2npc9)                     !srv 11.09.09 {
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)
        call ipgeti ('b2npp7_iout', iout_b2npp7)
        call ipgeti ('b2nppo_iout', iout_b2nppo)                     !srv 11.09.09 }
        call ipgetr ('b2stbc_bc_ref', bc_ref)
        call ipgetr ('b2stbc_bc_ref_te', bc_ref_te)                  !srv 11.04.13
        call ipgetr ('b2stbc_bc_ref_ti', bc_ref_ti)                  !srv 11.04.13
        call ipgeti ('b2stbc_bcene_15_style', bcene_15_style)
        call ipgeti ('b2stbc_bcene_16_style', bcene_16_style)
        call ipgeti ('b2stbc_bceni_16_style', bceni_16_style)
        call ipgeti ('b2news_poteq', pot_eq)                         !srv 27.02.13
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
*   ..test facdrift (not on first call)                                  !xpb
      if (ncall.gt.0)                                                    !xpb
     > call xertst (0.0_R8.le.minval(facdrift) .and.
     &       1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol, hx, hy
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
*    ..test state
       call b2xvps (nx, ny, ns, 
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
*    ..test edge values of fne, fni, chce, chve, chci, chvi, calf
       call b2xvff (nx, ny, fne, 'fne')
       call b2xvff (nx, ny, fni, 'fni')
       call b2xvff (nx, ny, chce, 'chce')
       call b2xvff (nx, ny, chve, 'chve')
       call b2xvff (nx, ny, chci, 'chci')
       call b2xvff (nx, ny, chvi, 'chvi')
       call b2xvff (nx, ny, calf, 'calf')
*    ..test sign of chce, chci
       call b2xvsg (n2*2, chce, 1, 'chce', '.ge.')
       call b2xvsg (n2*2, chci, 1, 'chci', '.ge.')
      endif
*   ..set cutll, cutlo
      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))

      if(ncall.eq.0) then
        call alloc_b2mod_boundary(ns,20,10)
        call read_b2mod_boundary_namelist(nx, ny, ns, .true.)
        call xertst (nbc.le.20, 'increase size in b2stbc for nbcd')
        call xertst (nniso.le.10, 'increase size in b2stbc for nnisod')
        call write_b2mod_boundary_namelist
        if(boundary_time_mod.gt.0.0_R8) then
          elm_count=tim/boundary_time_mod
        else
          elm_count=0
        endif
        catch_up=boundary_time_mod.gt.0.0_R8
        do while(catch_up)
          if(mod(tim,boundary_time_mod).gt.boundary_time_switch)
     1     then
            write(*,*) 'b2stbc: catching up ...'
            call read_b2mod_boundary_namelist(nx, ny, ns, .true.)
            call xertst (nbc.le.20, 'increase size in b2stbc for nbcd')
            call xertst (nniso.le.10, 
     1       'increase size in b2stbc for nnisod')
            call write_b2mod_boundary_namelist
            catch_up=boundary_time_mod.gt.0.0_R8
          else
            catch_up=.false.
          endif
        enddo
        if(nniso.gt.0) then
          done=.true.
          do ib = 1, nniso
            do ix = nxiso1(ib), nxiso2(ib)
              do iy = nyiso1(ib), nyiso2(ib)
                if(region(ix,iy,0).ne.0) then
                  done=.false.
                  region(ix,iy,0) = 0
                  resignore(ix,iy,1)=0
                  resignore(ix,iy,2)=0
                  write(*,*)
     .             'REGION set to zero for dead cell ',ix,iy
                endif
              enddo
            enddo
          enddo
          if(.not.done) then
            write(*,*) 'Isolated region was unknown from geometry'
            do ib = 1, nniso
              do ix = nxiso1(ib), nxiso2(ib)
                bottomiy(ix,nyiso1(ib))=nyiso2(ib)
                write(*,*) 'ix,iy,iyb',ix,nyiso1(ib),nyiso2(ib)
                topiy(ix,nyiso2(ib))=nyiso1(ib)
                write(*,*) 'ix,iy,iyt',ix,nyiso2(ib),nyiso1(ib)
              enddo
              do iy = nyiso1(ib),nyiso2(ib)
                leftix(nxiso1(ib),iy)=nxiso2(ib)
                write(*,*) 'ix,iy,ixl',nxiso1(ib),iy,nxiso2(ib)
                rightix(nxiso2(ib),iy)=nxiso1(ib)
                write(*,*) 'ix,iy,ixr',nxiso2(ib),iy,nxiso1(ib)
              enddo
            enddo
          endif
        endif
        do ib = 1, nbc
          if(bcchar(ib).eq.'N' .or. bcchar(ib).eq.'S' .or.
     &       bcchar(ib).eq.'X') then
            do is = 0, ns-1
              if (bccon(is,ib).ne.0)
     >         resignore(bcstart(ib):bcend(ib),bcpos(ib),1) = 0
              if (bcmom(is,ib).ne.0)
     >         resignore(bcstart(ib):bcend(ib),bcpos(ib),2) = 0
            end do
            if (bcene(ib).ne.0 .or. bceni(ib).ne.0)
     >       resignore(bcstart(ib):bcend(ib),bcpos(ib),1) = 0
          elseif(bcchar(ib).eq.'E' .or. bcchar(ib).eq.'W' .or.
     &           bcchar(ib).eq.'Y') then
            do is = 0, ns-1
              if (bccon(is,ib).ne.0)
     >         resignore(bcpos(ib),bcstart(ib):bcend(ib),1) = 0
              if (bcmom(is,ib).ne.0)
     >         resignore(bcpos(ib),bcstart(ib):bcend(ib),2) = 0
            end do
            if (bcene(ib).ne.0 .or. bceni(ib).ne.0)
     >       resignore(bcpos(ib),bcstart(ib):bcend(ib),1) = 0
          endif
        enddo
      else
        if(boundary_time_mod.gt.0.0_R8) then
          write(*,*) 'ELM: ',tim,mod(tim,boundary_time_mod),
     1     boundary_time_switch,int(tim/boundary_time_mod),elm_count
          if(mod(tim,boundary_time_mod).gt.boundary_time_switch.or.
     1     int(tim/boundary_time_mod).gt.elm_count) then
            elm_count=tim/boundary_time_mod
            call read_b2mod_boundary_namelist(nx, ny, ns, .true.)
            call xertst (nbc.le.20, 'increase size in b2stbc for nbcd')
            call xertst (nniso.le.10, 
     1       'increase size in b2stbc for nnisod')
            call write_b2mod_boundary_namelist
          endif
        endif
      endif
#ifdef PAR_OPT
      conpar(0,2,1) = par_opt_phys(4)    !  Be carefully!
	conpar(0,3,1) = par_opt_phys(4)    !  It depends on the file b2.boundary.parameters
	conpar(0,5,1) = par_opt_phys(4)
	conpar(0,6,1) = par_opt_phys(4)

      conpar(1,2,1) = par_opt_phys(5)
      conpar(1,3,1) = par_opt_phys(5)
      conpar(1,5,1) = par_opt_phys(5)
      conpar(1,6,1) = par_opt_phys(5)
      conpar(1,7,1) = par_opt_phys(5)
      conpar(1,8,1) = par_opt_phys(5)

	enepar(2,1) = par_opt_phys(6)
	enepar(3,1) = par_opt_phys(6)
	enepar(5,1) = par_opt_phys(6)
	enepar(6,1) = par_opt_phys(6)
	enepar(7,1) = par_opt_phys(6)
	enepar(8,1) = par_opt_phys(6)

	enipar(2,1) = par_opt_phys(7)
	enipar(3,1) = par_opt_phys(7)
	enipar(5,1) = par_opt_phys(7)
	enipar(6,1) = par_opt_phys(7)
	enipar(7,1) = par_opt_phys(7)
	enipar(8,1) = par_opt_phys(7)
#endif

*  -- INITIALISE CERTAIN CONSTANTS --
      EPSLON=1.e-10_R8
*  -- SET NE00, UP00, UBIG, VBIG --
      NE00=1.0e20_R8
      UP00=SQRT(100*EV/MP)
      UBIG=UP00 / EPSLON
      VBIG=1.0e-2_R8*UBIG

      sna0_no_mdf = 0.0_R8

	call b2xppz (nx, ny, ns, rza, na, te, ti, pz)               !srv 18.01.02
*   ..find the average pressure for inner flux surface            !srv 24.03.99
      if (istyle_average_pr .eq.0) then                            !srv 08.08.02
       iy = -1                                                     !srv 24.03.99
	 pzsum =0.0e0_R8                                             !srv 24.03.99 18.01.02
	 volsum=0.0e0_R8                                             !srv 24.03.99
	 do ix = -1,nx                                               !srv 26.11.02
        if (mod(region(ix,iy,0),4).eq.1) then                      !srv 26.11.02
	   pzsum=pzsum + pz(ix,iy)*vol(ix,iy)                        !srv 24.03.99 18.01.02
	   volsum=volsum + vol(ix,iy)                                !srv 24.03.99
        endif                                                      !srv 26.11.02
	 enddo                                                       !srv 24.03.99
       pzaverage=pzsum/volsum                                      !srv 24.03.99 18.01.02
      else                                                         !srv 08.08.02
       iy = 0                                                      !srv 24.03.99  07.06.02 {
	 pzsum =0.0e0_R8                                             !srv 24.03.99 18.01.02
	 volsum=0.0e0_R8                                             !srv 24.03.99
	 num_inner=0                                                 !srv 26.11.02
	 do ix = -1,nx                                               !srv 26.11.02
        if (mod(region(ix,iy,0),4).eq.1) then                      !srv 26.11.02
	   pzsum=pzsum
     *   + (pz(ix,iy)+pz(bottomix(ix,iy),bottomiy(ix,iy)))         !srv 24.03.99 18.01.02
	   volsum=volsum + vol(ix,iy)                                !srv 24.03.99
         num_inner=num_inner+1                                     !srv 26.11.02
        endif                                                      !srv 26.11.02
	 enddo                                                       !srv 24.03.99
	 pzaverage=pzsum*0.5e0_R8/num_inner                         !srv 24.03.99 18.01.02 07.06.02 } 26.11.02
      endif                                                        !srv 08.08.02
*   ..find the average ne, te, ti, bb for inner flux surface      !srv 15.12.00 18.10.04
      bbsum = 0.0e0_R8                                            !srv 18.10.04
	nesum = 0.0e0_R8                                            !srv 15.12.00 25.01.02
	tesum = 0.0e0_R8                                            !srv 15.12.00
	tisum = 0.0e0_R8                                            !srv 15.12.00
      bzhzsum  = 0.0e0_R8                                         !lk 20.11.07
      bbaverage  = 0.0e0_R8                                       !lk 20.11.07 
	do ix = -1,nx                                               !srv 26.11.02
       if (mod(region(ix,iy,0),4).eq.1) then                      !srv 26.11.02
	  bbsum=bbsum + bb(ix,iy,3)*vol(ix,iy)                      !srv 18.10.04
	  nesum=nesum + ne(ix,iy)*vol(ix,iy)                        !srv 15.12.00 25.01.02
	  tisum=tisum + ti(ix,iy)*vol(ix,iy)                        !srv 15.12.00
	  tesum=tesum + te(ix,iy)*vol(ix,iy)                        !srv 15.12.00
        bzhzsum = bzhzsum+abs(bb(ix,iy,2))*vol(ix,iy)**2/gs(ix,iy,2)   !lk 20.11.07 !srv 07.05.13
       endif                                                      !srv 26.11.02
	enddo                                                       !srv 15.12.00
	bbaverage=bbsum/volsum                                      !srv 18.10.04
	neaverage=nesum/volsum                                      !srv 15.12.00 25.01.02
	tiaverage=tisum/volsum                                      !srv 15.12.00
	teaverage=tesum/volsum                                      !srv 15.12.00
	bzhzaverage=bzhzsum/volsum                                  !lk 20.11.07
*   ..find the average ua**2, na for inner flux surface           !srv 22.01.02 {
      do is = 0, ns-1
	 uasum  = 0.0e0_R8                                           !srv 15.02.02
	 ua2sum = 0.0e0_R8                                           !srv 07.04.00 15.02.02
	 nasum  = 0.0e0_R8
	 do ix = -1,nx                                               !srv 26.11.02
        if (mod(region(ix,iy,0),4).eq.1) then                      !srv 26.11.02
	   ua2sum = ua2sum + ua(ix,iy,is)**2*vol(ix,iy)              !srv 07.04.00 15.02.02
	   uasum  = uasum  + ua(ix,iy,is)*vol(ix,iy)                 !srv 15.02.02
	   nasum  = nasum  + na(ix,iy,is)*vol(ix,iy)
        endif                                                      !srv 26.11.02
	 enddo
	 ua2average(is) = ua2sum/volsum                              !srv 07.04.00 15.02.02
	 uaaverage(is) = uasum/volsum                                !srv 07.04.00 15.02.02 16.04.02
	 naaverage(is) = nasum/volsum
	enddo                                                        !srv 22.01.02 }

      call b2xpro (nx, ny, ns, mp, am, na, ro)                       !srv 01.02.09
      call b2xppr (nx, ny, ns, ne, na, te, ti, pr)
      wdia(:,:)=(ne(:,:)*te(:,:)+ni(:,:,1)*ti(:,:))*bb(:,:,2)
      do iy = -1, ny
        do ix = -1, nx
          if (rightix(ix,iy).ne.nx+1) then
            sx(ix,iy)=gs(rightix(ix,iy),rightiy(ix,iy),0)
            sxhz(ix,iy)=gs(rightix(ix,iy),rightiy(ix,iy),0)*         !srv 24.08.08
     &                  hz(ix,iy)
!     &                  hz(rightix(ix,iy),rightiy(ix,iy))
          else
            sx(ix,iy)=0.0_R8
            sxhz(ix,iy)=0.0_R8                                       !srv 24.08.08
          endif
          if (leftix(ix,iy).ne.-2) then
            uu(ix,iy,:)=fna(ix,iy,0,:)/na(ix,iy,:)/gs(ix,iy,0)
          else
            uu(ix,iy,:)=0.0_R8
          endif
        enddo
      enddo
      if (area_fix.ge.3) sx(:,:) = sx(:,:)*qc(:,:)                     !srv 09.03.09
      if (area_fix.ge.3) sxhz(:,:) = sxhz(:,:)*qc(:,:)                 !srv 27.04.09
      if (area_fix.ge.1) then
        do is = 0, ns-1
          uu(:,:,is) = uu(:,:,is)/qc(:,:)
        enddo
      endif
      if (main_call .and. diagno.ge.1) then                            !xpb
        write(*,'(a,1p,(10000e12.4))') 'dna0',
     1   dna0(nx/2,ny,0:ns-1)
        write(*,'(a,1p,(10000e12.4))') 'dpa0',
     1   dpa0(nx/2,ny,0:ns-1)
        write(*,'(a,1p,(10000e12.4))') 'vsa0',
     1   vsa0(nx/2,ny,0:ns-1)
        write(*,'(a,1p,(10000e12.4))') 'hce0',
     1   hce0(nx/2,ny)
        write(*,'(a,1p,(10000e12.4))') 'hci0',
     1   hci0(nx/2,ny)
        write(*,'(a,1p,(10000e12.4))') 'sig0',
     1   sig0(nx/2,ny)
      endif                                                            !xpb
        if(ncall.eq.0) write(*,*) 'Doing B2-style bc'
        once_bccon_16 = .false.                                          !srv 16.02.12 {
        once_bccon_19 = .false.
        once_bcene_16 = .false.
        once_bceni_16 = .false.                                          !srv 16.02.12 }
* -- LOOP OVER ALL BOUNDARIES
        DO IB=1, NBC

          DO IS=0, NS-1

* -- BOUNDARY CONDITIONS FOR CONTINUITY EQUATIONS --
            if(BCCON(IS,IB).eq.0) goto 10000
	    GOTO(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,
     1	     1400,1500,1600,1700,1800,1900,2000,2100,2200,2300),     !srv 19.04.13
     .       BCCON(IS,IB)
            write(*,*) BCCON(IS,IB)
            call xerrab ('BCCON out of bounds')

* -- BCCON=1 -- PRESCRIBE THE VALUE OF THE DENSITY

  100       CONTINUE
            if(ncall.eq.0) then
#ifdef ASTRA
             IF(BCCHAR(IB).EQ.'S') THEN
	        if(is.eq.ismain) then
	         CONPAR(IS,IB,1) = neb_astra                               !srv 28.07.11
	        endif
             ENDIF
#endif
              write(*,*) 'constant density',CONPAR(IS,IB,1),' on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS
              call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1         'BCCON = 1, CONPAR(,,1) !> 0')
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)= VBIG*CONPAR(IS,IB,1)*
     .           SY(IX,BCPOS(IB))
                SNA0(IX,BCPOS(IB),1,IS)=-VBIG*SY(IX,BCPOS(IB))
              end do
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)= VBIG*CONPAR(IS,IB,1)*
     .           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
                SNA0(IX,BCPOS(IB),1,IS)=-VBIG*
     .           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
              end do
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'W') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SNA0(BCPOS(IB),IY,0,IS)= UBIG*CONPAR(IS,IB,1)*
     .           SX(BCPOS(IB),IY)
                SNA0(BCPOS(IB),IY,1,IS)=-UBIG*SX(BCPOS(IB),IY)
              end do
              DO IY=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(BCPOS(IB),IY,:,IS)=SNA0(BCPOS(IB),IY,:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SNA0(BCPOS(IB),IY,0,IS)= UBIG*CONPAR(IS,IB,1)*
     .           SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
                SNA0(BCPOS(IB),IY,1,IS)=-UBIG*
     .           SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
              end do
              DO IY=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(BCPOS(IB),IY,:,IS)=SNA0(BCPOS(IB),IY,:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'X') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)= UBIG*
     .           CONPAR(IS,IB,1)*VOL(IX,BCPOS(IB))
                SNA0(IX,BCPOS(IB),1,IS)=-UBIG*
     .           VOL(IX,BCPOS(IB))
              end do
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'Y') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SNA0(BCPOS(IB),IY,0,IS)= UBIG*
     .           CONPAR(IS,IB,1)*VOL(BCPOS(IB),IY)
                SNA0(BCPOS(IB),IY,1,IS)=-UBIG*
     .           VOL(BCPOS(IB),IY)
              end do
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(BCPOS(IB),IY,:,IS)=SNA0(BCPOS(IB),IY,:,IS)
              ENDDO                                                     !srv 23.09.08 }
           ELSE
              CALL XERRAB
     1         ('B2STBC -- BCCON=1, INVALID OPTION.')
            ENDIF
            GOTO 10000

* -- BCCON=2 -- PRESCRIBE THE GRADIENT OF THE DENSITY

  200       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'density gradient',CONPAR(IS,IB,1),' on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)= VBIG*SY(IX,BCPOS(IB))*
     .           (-CONPAR(IS,IB,1)*0.5_R8*
     .           (HY1(IX,BCPOS(IB))+
     .            HY1(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)))) +
     .           na(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)),IS))
                SNA0(IX,BCPOS(IB),1,IS)=-VBIG*SY(IX,BCPOS(IB))
              end do
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)= VBIG*
     .           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))*
     .           (CONPAR(IS,IB,1)*0.5_R8*
     .           (HY1(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))+
     .            HY1(IX,BCPOS(IB))) +
     .           na(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)),IS))
                SNA0(IX,BCPOS(IB),1,IS)=-VBIG*
     .           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
              end do
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'W') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SNA0(BCPOS(IB),IY,0,IS)= UBIG*SX(BCPOS(IB),IY)*
     .           (-CONPAR(IS,IB,1)*0.5_R8*
     .           (HX(BCPOS(IB),IY)+
     .            HX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))) +
     .           na(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY),IS))
                SNA0(BCPOS(IB),IY,1,IS)=-UBIG*SX(BCPOS(IB),IY)
              end do
              DO IY=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(BCPOS(IB),IY,:,IS)=SNA0(BCPOS(IB),IY,:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SNA0(BCPOS(IB),IY,0,IS)= UBIG*
     .           SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))*
     .           (CONPAR(IS,IB,1)*0.5_R8*
     .           (HX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))+
     .            HX(BCPOS(IB),IY)) +
     .           na(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY),IS))
                SNA0(BCPOS(IB),IY,1,IS)=-UBIG*
     .           SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
              end do
              DO IY=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(BCPOS(IB),IY,:,IS)=SNA0(BCPOS(IB),IY,:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCCON=2, INVALID OPTION.')
            END IF
            GOTO 10000

* -- BCCON=3 -- SHEATH CONDITIONS

  300       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'sheath density bc on ',BCCHAR(IB),
     1         ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB),
     2         ' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'W') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SNA0(BCPOS(IB),IY,0,IS)= UBIG*
     .           na(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY),IS)*
     .           SX(BCPOS(IB),IY)
                SNA0(BCPOS(IB),IY,1,IS)=-UBIG*SX(BCPOS(IB),IY)
              end do
              DO IY=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(BCPOS(IB),IY,:,IS)=SNA0(BCPOS(IB),IY,:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SNA0(BCPOS(IB),IY,0,IS)= UBIG*
     .           na(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY),IS)*
     .           SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
                SNA0(BCPOS(IB),IY,1,IS)=-UBIG*
     .           SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
              end do
              DO IY=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(BCPOS(IB),IY,:,IS)=SNA0(BCPOS(IB),IY,:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCCON=3, NO SHEATH IN RADIAL DIR.')
            END IF
            GOTO 10000

* -- BCCON=4 -- PRESCRIBE THE VALUE OF THE DENSITY, WEAKLY
*                 A MIXED BOUNDARY CONDITION

  400       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'weak density',CONPAR(IS,IB,1),' on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2         ' to ',BCEND(IB),' for species ',IS
              call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1         'BCCON = 4, CONPAR(,,1) !> 0')
              call xertst(CONPAR(IS,IB,2).gt.0.0_R8,
     1         'BCCON = 4, CONPAR(,,2) !> 0')
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)= CONPAR(IS,IB,2)
     .           *CONPAR(IS,IB,1)*SY(IX,BCPOS(IB))
                SNA0(IX,BCPOS(IB),1,IS)=
     1           -CONPAR(IS,IB,2)*SY(IX,BCPOS(IB))
              end do
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)= CONPAR(IS,IB,2)*
     .           CONPAR(IS,IB,1)*
     .           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
                SNA0(IX,BCPOS(IB),1,IS)=-CONPAR(IS,IB,2)*
     .           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
              end do
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'W') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SNA0(BCPOS(IB),IY,0,IS)= CONPAR(IS,IB,2)
     .           *CONPAR(IS,IB,1)*SX(BCPOS(IB),IY)
                SNA0(BCPOS(IB),IY,1,IS)=
     1           -CONPAR(IS,IB,2)*SX(BCPOS(IB),IY)
              end do
              DO IY=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(BCPOS(IB),IY,:,IS)=SNA0(BCPOS(IB),IY,:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SNA0(BCPOS(IB),IY,0,IS)= CONPAR(IS,IB,2)*
     .           CONPAR(IS,IB,1)*
     .           SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
                SNA0(BCPOS(IB),IY,1,IS)=-CONPAR(IS,IB,2)*
     .           SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
              end do
              DO IY=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(BCPOS(IB),IY,:,IS)=SNA0(BCPOS(IB),IY,:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCCON=4, NOT YET AVAILABLE')
            END IF
            GOTO 10000

* -- BCCON=5 -- PRESCRIBE THE PARTICLE FLUX PER UNIT AREA

  500       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'particle flux per area',CONPAR(IS,IB,1),
     1         ' on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',
     2         BCSTART(IB),' to ',BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)=
     1           CONPAR(IS,IB,1)*SY(IX,BCPOS(IB))
                SNA0(IX,BCPOS(IB),1,IS)=0.0_R8
              end do
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)=CONPAR(IS,IB,1)*
     1           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
                SNA0(IX,BCPOS(IB),1,IS)=0.0_R8
              end do
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'W') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SNA0(BCPOS(IB),IY,0,IS)=
     1           CONPAR(IS,IB,1)*SX(BCPOS(IB),IY)
                SNA0(BCPOS(IB),IY,1,IS)=0.0_R8
              end do
              DO IY=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(BCPOS(IB),IY,:,IS)=SNA0(BCPOS(IB),IY,:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SNA0(BCPOS(IB),IY,0,IS)=CONPAR(IS,IB,1)*
     1           SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
                SNA0(BCPOS(IB),IY,1,IS)=0.0_R8
              end do
              DO IY=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(BCPOS(IB),IY,:,IS)=SNA0(BCPOS(IB),IY,:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCCON=5, NOT YET AVAILABLE')
            END IF

            GOTO 10000

* -- BCCON=6 -- PRESCRIBE THE TOTAL PARTICLE FLUX FOR A
*                 CONSTANT DENSITY

cghgh
  600       CONTINUE
            if(ncall.eq.0) then
             write(*,*) 'total particle flux',CONPAR(IS,IB,1),
     1         ' with const density on ',BCCHAR(IB),' at ',BCPOS(IB),
     2         ' from ',BCSTART(IB),' to ',BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'N') THEN
              s1=0.0_R8
              s2=0.0_R8
              DO IX=BCSTART(IB), BCEND(IB)
                s1=s1+SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))/
     .           (0.5_R8*
     .           (HY1(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))+
     .            HY1(IX,BCPOS(IB))))
                s2=s2+SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))*
     .            na(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)),IS)/
     .           (0.5_R8*
     .           (HY1(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))+
     .            HY1(IX,BCPOS(IB))))
              ENDDO
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)= VBIG*
     .           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))*
     .           (CONPAR(IS,IB,1)/dna0(ix,bcpos(ib),is) + s2) / s1
                SNA0(IX,BCPOS(IB),1,IS)=-VBIG*
     .           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
              end do
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'S') THEN
              s1=0.0_R8
              s2=0.0_R8
              DO IX=BCSTART(IB), BCEND(IB)
                s1=s1+SY(IX,BCPOS(IB))/
     .           (0.5_R8*(HY1(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)))+
     .                    HY1(IX,BCPOS(IB))))
                s2=s2+SY(IX,BCPOS(IB))*
     .           na(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)),IS)/
     .           (0.5_R8*(HY1(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)))+
     .                    HY1(IX,BCPOS(IB))))
              ENDDO
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)= VBIG*SY(IX,BCPOS(IB))*
     .           (CONPAR(IS,IB,1)/dna0(ix,bcpos(ib),is) + s2) / s1
                SNA0(IX,BCPOS(IB),1,IS)= - VBIG*SY(IX,BCPOS(IB))
              end do
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSE
              CALL XERRAB ('B2STBC -- BCCON=6, NOT AVAILABLE FOR W,E')
            ENDIF
            GOTO 10000
cghgh

* -- BCCON=7 -- PRESCRIBE THE DENSITY AS A FUNCTION OF OTHER
*                 PLASMA PARAMETERS

  700       CONTINUE
            CALL XERRAB ('B2STBC -- BCCON=7, NOT YET AVAILABLE')
            GOTO 10000

* -- BCCON=8 -- PRESCRIBE THE TOTAL PARTICLE FLUX WITH
*                 CONSTANT FLUX DENSITY

  800       CONTINUE
            if(ncall.eq.0) then
#ifdef ASTRA
       if(is.eq.ismain.and.BCCHAR(IB).eq.'S') then                  !srv 20.06.12
        CONPAR(IS,IB,1) = fneb_astra                                !srv 04.10.11
	 endif
#endif
              write(*,*) 'total particle flux', CONPAR(IS,IB,1),' on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              US=0.0_R8
              DO IX=BCSTART(IB), BCEND(IB)
                US=US+SY(IX,BCPOS(IB))
              ENDDO
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)=SNA0(IX,BCPOS(IB),0,IS)+         !srv 01.02.09
     1           CONPAR(IS,IB,1)/US*SY(IX,BCPOS(IB))
                SNA0(IX,BCPOS(IB),1,IS)=0.0_R8
                sna0_no_mdf(ix,bcpos(ib),0,is) =
     &            sna0_no_mdf(ix,bcpos(ib),0,is)+sna0(ix,bcpos(ib),0,is) !srv 23.09.08
                sna0_no_mdf(ix,bcpos(ib),1,is) =
     &            sna0_no_mdf(ix,bcpos(ib),1,is)+sna0(ix,bcpos(ib),1,is) !srv 23.09.08
	          if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                 sna0(ix,bcpos(ib),0,is) = sna0(ix,bcpos(ib),0,is) -
     &                              fnaPSch(ix,topiy(ix,bcpos(ib)),1,is)
                 sna0(ix,bcpos(ib),1,is) = sna0(ix,bcpos(ib),1,is) +
     &            vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &               vaecrb(ix,topiy(ix,bcpos(ib)),1,is)
	          endif                                                    !srv 12.09.08 }
              ENDDO
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              US=0.0_R8
              DO IX=BCSTART(IB), BCEND(IB)
                US=US+SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
              ENDDO
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)=SNA0(IX,BCPOS(IB),0,IS)+         !srv 01.02.09
     1           CONPAR(IS,IB,1)/US*
     2           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
                SNA0(IX,BCPOS(IB),1,IS)=0.0_R8
                sna0_no_mdf(ix,bcpos(ib),0,is) =
     &            sna0_no_mdf(ix,bcpos(ib),0,is)+sna0(ix,bcpos(ib),0,is) !srv 23.09.08
                sna0_no_mdf(ix,bcpos(ib),1,is) =
     &            sna0_no_mdf(ix,bcpos(ib),1,is)+sna0(ix,bcpos(ib),1,is) !srv 23.09.08
	          if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                 sna0(ix,bcpos(ib),0,is) = sna0(ix,bcpos(ib),0,is) +
     &                                  fnaPSch(ix,bcpos(ib),1,is)
                 sna0(ix,bcpos(ib),1,is) = sna0(ix,bcpos(ib),1,is) -
     &            vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &               vaecrb(ix,bcpos(ib),1,is)
	          endif                                                    !srv 12.09.08 }
              ENDDO
            ELSEIF(BCCHAR(IB).EQ.'W') THEN
              US=0.0_R8
              DO IY=BCSTART(IB), BCEND(IB)
                US=US+SX(BCPOS(IB),IY)
              ENDDO
              DO IY=BCSTART(IB), BCEND(IB)
                SNA0(BCPOS(IB),IY,0,IS)=
     1           CONPAR(IS,IB,1)/US*SX(BCPOS(IB),IY)
                SNA0(BCPOS(IB),IY,1,IS)=0.0_R8
                sna0_no_mdf(bcpos(ib),iy,0,is) = sna0(bcpos(ib),iy,0,is) !srv 23.09.08
                sna0_no_mdf(bcpos(ib),iy,1,is) = sna0(bcpos(ib),iy,1,is) !srv 23.09.08
	          if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                 sna0(bcpos(ib),iy,0,is) = sna0(bcpos(ib),iy,0,is) -
     &                          fnaPSch(rightix(bcpos(ib),iy),iy,0,is)
                 sna0(bcpos(ib),iy,1,is) = sna0(bcpos(ib),iy,1,is) +
     &            vol(bcpos(ib),iy)/hy1(bcpos(ib),iy)*
     &               vaecrb(bcpos(ib),iy,0,is)
	          endif                                                    !srv 12.09.08 }
              ENDDO
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              US=0.0_R8
              DO IY=BCSTART(IB), BCEND(IB)
                US=US+SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
              ENDDO
              DO IY=BCSTART(IB), BCEND(IB)
                SNA0(BCPOS(IB),IY,0,IS)=
     1           CONPAR(IS,IB,1)/US*
     2           SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
                SNA0(BCPOS(IB),IY,1,IS)=0.0_R8
                sna0_no_mdf(bcpos(ib),iy,0,is) = sna0(bcpos(ib),iy,0,is) !srv 23.09.08
                sna0_no_mdf(bcpos(ib),iy,1,is) = sna0(bcpos(ib),iy,1,is) !srv 23.09.08
	          if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                 sna0(bcpos(ib),iy,0,is) = sna0(bcpos(ib),iy,0,is) +
     &                          fnaPSch(bcpos(ib),iy,0,is)
                 sna0(bcpos(ib),iy,1,is) = sna0(bcpos(ib),iy,1,is) -
     &            vol(bcpos(ib),iy)/hy1(bcpos(ib),iy)*
     &               vaecrb(bcpos(ib),iy,0,is)
	          endif                                                    !srv 12.09.08 }
              ENDDO
            ELSE
              CALL XERRAB ('B2STBC -- BCCON=8, NOT AVAILABLE')
            ENDIF
            GOTO 10000

* -- BCCON=9 -- PRESCRIBE THE DECAY LENGTH FOR THE DENSITY

  900       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'density decay length',CONPAR(IS,IB,1),' on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS
              call xertst(CONPAR(IS,IB,1).ne.0.0_R8,
     1         'B2STBC -- BCCON=9, CONPAR(,,1) = 0')
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)=0.0_R8
                SNA0(IX,BCPOS(IB),1,IS)=
     1           -(dna0(IX,BCPOS(IB),IS)+
     2           dpa0(IX,BCPOS(IB),IS)*
     3           (rza(ix,bcpos(ib),is)*te(ix,bcpos(ib))+
     .            ti(ix,bcpos(ib))))/
     2           CONPAR(IS,IB,1)*SY(IX,BCPOS(IB))
                sna0_no_mdf(ix,bcpos(ib),0,is) = sna0(ix,bcpos(ib),0,is) !srv 23.09.08
                sna0_no_mdf(ix,bcpos(ib),1,is) = sna0(ix,bcpos(ib),1,is) !srv 23.09.08
	          if(mdf_fnb.ne.0) then                                    !srv 08.09.08 {
                 sna0(ix,bcpos(ib),0,is) = sna0(ix,bcpos(ib),0,is) -
     &                              fnaPSch(ix,topiy(ix,bcpos(ib)),1,is)
                 sna0(ix,bcpos(ib),1,is) = sna0(ix,bcpos(ib),1,is)+
     &            vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &               vaecrb(ix,topiy(ix,bcpos(ib)),1,is)
	          endif                                                    !srv 08.09.08 }
              ENDDO
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)=0.0_R8
                SNA0(IX,BCPOS(IB),1,IS)=
     1           -(dna0(IX,BCPOS(IB),IS)+
     2           dpa0(IX,BCPOS(IB),IS)*
     3           (rza(ix,bcpos(ib),is)*te(ix,bcpos(ib))+
     .            ti(ix,bcpos(ib))))/
     2           CONPAR(IS,IB,1)*
     2           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
                sna0_no_mdf(ix,bcpos(ib),0,is) = sna0(ix,bcpos(ib),0,is) !srv 23.09.08
                sna0_no_mdf(ix,bcpos(ib),1,is) = sna0(ix,bcpos(ib),1,is) !srv 23.09.08
	          if(mdf_fnb.ne.0) then                                    !srv 08.09.08 {
                 sna0(ix,bcpos(ib),0,is) = sna0(ix,bcpos(ib),0,is) +
     &                                  fnaPSch(ix,bcpos(ib),1,is)
                 sna0(ix,bcpos(ib),1,is) = sna0(ix,bcpos(ib),1,is)-
     &            vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &               vaecrb(ix,bcpos(ib),1,is)
	          endif                                                    !srv 08.09.08 }
              ENDDO
            ELSE
              CALL XERRAB ('B2STBC -- BCCON=9, NOT AVAILABLE FOR W,E')
            ENDIF
            GOTO 10000

* -- BCCON=10 -- leakage option for density

 1000       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'leakage',CONPAR(IS,IB,1),' on ',BCCHAR(IB),
     1         ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB),
     2         ' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                s0 = vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))                 !srv 26.02.09
                SNA0(IX,BCPOS(IB),0,IS)=0.0_R8
                SNA0(IX,BCPOS(IB),1,IS)=SNA0(IX,BCPOS(IB),1,IS)+         !srv 01.02.09
     1           s0*sqrt((rza(ix,bcpos(ib),is)*te(ix,bcpos(ib))+         !srv 26.02.09
     2                 ti(ix,bcpos(ib)))/
     3           (am(is)*mp))*
     4           CONPAR(IS,IB,1)
                sna0_no_mdf(ix,bcpos(ib),:,is) =
     &            sna0_no_mdf(ix,bcpos(ib),:,is)+sna0(ix,bcpos(ib),:,is) !srv 23.09.08
*    ..drifts to the wall                                                !srv 02.03.05 {
	          if(mdf_fnb.eq.0) then                                    !srv 09.11.06 08.09.08 {
                 sna0(ix,bcpos(ib),1,is) = sna0(ix,bcpos(ib),1,is)+
     *            min(0.0e0_R8,s0*(vadia(ix,topiy(ix,bcpos(ib)),1,is)+
     &                vaecrb(ix,topiy(ix,bcpos(ib)),1,is)))
	          else
                 sna0(ix,bcpos(ib),0,is) = sna0(ix,bcpos(ib),0,is) -
     &                              fnaPSch(ix,topiy(ix,bcpos(ib)),1,is)
                 sna0(ix,bcpos(ib),1,is) = sna0(ix,bcpos(ib),1,is)+
     &            vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &               min(0.0e0_R8,vaecrb(ix,topiy(ix,bcpos(ib)),1,is))   !srv 02.03.05 }
	          endif                                                    !srv 09.11.06 08.09.08 }
              ENDDO
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
	          s0=vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))                  !srv 26.02.09
                SNA0(IX,BCPOS(IB),0,IS)=0.0_R8
                SNA0(IX,BCPOS(IB),1,IS)=SNA0(IX,BCPOS(IB),1,IS)+        !srv 01.02.09
     1           s0*sqrt((rza(ix,bcpos(ib),is)*te(ix,bcpos(ib))+        !srv 26.02.09
     2                 ti(ix,bcpos(ib)))/
     3           (am(is)*mp))*
     4           CONPAR(IS,IB,1)
                sna0_no_mdf(ix,bcpos(ib),:,is) =
     &            sna0_no_mdf(ix,bcpos(ib),:,is)+sna0(ix,bcpos(ib),:,is) !srv 23.09.08
*    ..drifts to the wall                                                !srv 02.03.05 {
	          if(mdf_fnb.eq.0) then                                    !srv 09.11.06 08.09.08 {
                 sna0(ix,bcpos(ib),1,is) = sna0(ix,bcpos(ib),1,is)-
     *           max(0.0e0_R8,s0*(vadia(ix,bcpos(ib),1,is)+
     *                           vaecrb(ix,bcpos(ib),1,is)))
	          else
                 sna0(ix,bcpos(ib),0,is) = sna0(ix,bcpos(ib),0,is) +
     &                                  fnaPSch(ix,bcpos(ib),1,is)
                 sna0(ix,bcpos(ib),1,is) = sna0(ix,bcpos(ib),1,is)-
     &            max(0.0e0_R8,(vol(ix,bcpos(ib))+
     &             vol(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib))))/
     &            (hy1(ix,bcpos(ib))+
     &             hy1(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib))))      !srv 07.05.08
     &             *vaecrb(ix,bcpos(ib),1,is))                                   !lk 12.12.07 }
	          endif                                                    !srv 09.11.06 08.09.08 }
              ENDDO
            ELSEIF(BCCHAR(IB).EQ.'W') THEN
              CALL XERRAB ('B2STBC -- BCCON=10, NOT AVAILABLE FOR WEST') !srv 08.09.08 }
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              CALL XERRAB ('B2STBC -- BCCON=10, NOT AVAILABLE FOR EAST') !srv 08.09.08 }
            ELSE
              CALL XERRAB ('B2STBC -- BCCON=10, NOT AVAILABLE')
            ENDIF
            GOTO 10000

* -- BCCON=11 -- particle flux feedback boundary condition

 1100       CONTINUE
            if(ncall.eq.0) then
             write(*,*) 'particle flux feedback boundary condition on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS
             call ipgeti ('b2stbc_isfeedback', isfeedback)
             call xertst (0.le.isfeedback.and.isfeedback.lt.ns,
     1        'faulty internal parameter isfeedback')
             call xertst(is.eq.isfeedback,
     1         'BCCON(,IS,1) = 11 for IS <> ISFEEDBACK!')
             call xertst(lfeedback, 'LFEEDBACK must be turned on!')
            endif
            IF (BCCHAR(IB).EQ.'S') THEN
*   ..initialise ireg
              ireg = cbirso
*   ..update ireg
 1101	        if (cbrbrk(ireg)*nx.lt.bcend(ib)+0.5_R8) then
                call xertst (ireg.lt.cbirso+cbnrso-1,
     &           'error in computing region on south boundary')
                ireg = ireg+1
                goto 1101
              endif
              if (ireg.gt.cbirso) then
                call xertst (cbrbrk(ireg-1)*nx.lt.bcstart(ib)+0.5_R8,
     &           'feedback bc covers more than one region!')
              endif
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)=cbsna(0,IS,ireg)*
     .           SY(IX,BCPOS(IB))
                SNA0(IX,BCPOS(IB),1,IS)=0.0_R8
              ENDDO
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF (BCCHAR(IB).EQ.'N') THEN
*   ..initialise ireg
              ireg = cbirno
*   ..update ireg
 1102	      if (cbrbrk(ireg)*nx.lt.bcend(ib)+0.5_R8) then
                call xertst (ireg.lt.cbirno+cbnrno-1,
     &           'error in computing region on north boundary')
                ireg = ireg+1
                goto 1102
              endif
              if (ireg.gt.cbirno) then
                call xertst (cbrbrk(ireg-1)*nx.lt.bcstart(ib)+0.5_R8,
     &           'feedback bc covers more than one region!')
              endif
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)=cbsna(0,IS,ireg)*
     .           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
                SNA0(IX,BCPOS(IB),1,IS)=0.0_R8
              ENDDO
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSE
              CALL XERRAB 
     .         ('B2STBC -- BCCON=11, NOT AVAILABLE FOR W,E')
            ENDIF
            GOTO 10000

* -- BCCON=12 -- particle density feedback boundary condition

 1200       CONTINUE
            if(ncall.eq.0) then
             write(*,*) 
     .         'particle density feedback boundary condition on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS
             call ipgeti ('b2stbc_isfeedback', isfeedback)
             call xertst (0.le.isfeedback.and.isfeedback.lt.ns,
     1        'faulty internal parameter isfeedback')
             call xertst(is.eq.isfeedback,
     1         'BCCON(,IS,1) = 12 for IS <> ISFEEDBACK!')
             call xertst(lfeedback, 'LFEEDBACK must be turned on!')
             call ipgetr ('b2stbc_nesepm', nesepm)
             call ipgetr ('b2stbc_nesepm_overshoot', nesepm_overshoot)
             call xertst (0.0_R8.eq.nesepm.or.
     1        0.0_R8.eq.nesepm_overshoot,
     1        'BCCON=12 coding incompatible with nesepm_overshoot!')
            endif
            IF (BCCHAR(IB).EQ.'S') THEN
*   ..initialise ireg
              ireg = cbirso
*   ..update ireg
 1201	      if (cbrbrk(ireg)*nx.lt.bcend(ib)+0.5_R8) then
                call xertst (ireg.lt.cbirso+cbnrso-1,
     &           'error in computing region on south boundary')
                ireg = ireg+1
                goto 1201
              endif
              if (ireg.gt.cbirso) then
                call xertst (cbrbrk(ireg-1)*nx.lt.bcstart(ib)+0.5_R8,
     &           'feedback bc covers more than one region!')
              endif
              call xertst (cbsna(0,is,ireg).gt.0.0_R8,
     &         'BCCON=12, trying to set a negative density!')
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)= VBIG*cbsna(0,IS,ireg)*
     .           SY(IX,BCPOS(IB))
                SNA0(IX,BCPOS(IB),1,IS)=-VBIG*SY(IX,BCPOS(IB))
              ENDDO
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSEIF (BCCHAR(IB).EQ.'N') THEN
*   ..initialise ireg
              ireg = cbirno
*   ..update ireg
 1202	      if (cbrbrk(ireg)*nx.lt.bcend(ib)+0.5_R8) then
                call xertst (ireg.lt.cbirno+cbnrno-1,
     &           'error in computing region on north boundary')
                ireg = ireg+1
                goto 1202
              endif
              if (ireg.gt.cbirno) then
                call xertst (cbrbrk(ireg-1)*nx.lt.bcstart(ib)+0.5_R8,
     &           'feedback bc covers more than one region!')
              endif
              call xertst (cbsna(0,is,ireg).gt.0.0_R8,
     &         'BCCON=12, trying to set a negative density!')
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)= VBIG*cbsna(0,IS,ireg)*
     .           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
                SNA0(IX,BCPOS(IB),1,IS)=-VBIG*
     .           SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
              ENDDO
               DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
           ELSE
              CALL XERRAB 
     .         ('B2STBC -- BCCON=12, NOT AVAILABLE FOR W,E')
            ENDIF
            GOTO 10000

* -- BCCON=13 -- IN test --- DPC

 1300       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 
     .         'particle density to achieve specified total flux on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS, ' with parameters ',
     3         CONPAR(IS,IB,1),CONPAR(IS,IB,2)
            endif
c dpc
            IF(BCCHAR(IB).EQ.'S') THEN
              FS=0.0_R8
              VS=0.0_R8
              NAS=0.0_R8
              DO IX=BCSTART(IB), BCEND(IB)
                FS=FS+FNA(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1,is)
                VS=VS+VOL(ix,BCPOS(IB))
                NAS=NAS+NA(ix,BCPOS(IB),is)*vol(ix,BCPOS(IB))
              ENDDO
              NAS=NAS/VS
              FFAC=(1.0_R8+CONPAR(IS,IB,2)*
     1         (CONPAR(IS,IB,1)-FS)/
     2         (abs(CONPAR(IS,IB,1))+abs(FS)+
     3	       abs(bc_type13_norm)+abs(NAS*bc_type13_fac)))
              NAS=NAS*FFAC
              write(*,*)
     1         'b2stbc_conpar=13: is,curr_flux,desired_flux,density ',
     2         IS,FS,CONPAR(IS,IB,1),NAS
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)=VBIG*NAS*SY(IX,BCPOS(IB))
                SNA0(IX,BCPOS(IB),1,IS)=-VBIG*SY(IX,BCPOS(IB))
              ENDDO
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCCON=13, ONLY CODED FOR SOUTH')
            ENDIF
            GOTO 10000

* -- BCCON=14 -- sound speed flux                               !srv 14.04.08 {

 1400       CONTINUE
            if(ncall.eq.0) then
              write(*,*) ' sound speed flux',CONPAR(IS,IB,1),' on ',
     1          BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2          BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'W') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                S1 = ABS(PBS(RIGHTIX(BCPOS(IB),IY),
     *                       RIGHTIY(BCPOS(IB),IY),0))
                CS=SQRT(PR(BCPOS(IB),IY)/RO(BCPOS(IB),IY))
!                CS=SQRT(BCINT(PR,NX,NY,IY,IB)/BCINT(RO,NX,NY,IY,IB))
                SNA0(BCPOS(IB),IY,0,IS)=0.0_R8
                SNA0(BCPOS(IB),IY,1,IS)=SNA0(BCPOS(IB),IY,1,IS)+
     &                                  CONPAR(IS,IB,1)*CS*S1
                sna0_no_mdf(BCPOS(IB),IY,:,IS)=
     &	    sna0_no_mdf(BCPOS(IB),IY,:,IS) + SNA0(BCPOS(IB),IY,:,IS)
	          if(mdf_fnb.ne.0) then                                   !srv 13.10.06 {
                 sna0(BCPOS(IB),IY,0,is) = sna0(BCPOS(IB),IY,0,is)-
     &                          fnaPSch(rightix(BCPOS(IB),iy),iy,0,is)
	          endif                                                   !srv 13.10.06 }
              ENDDO
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                S1 = ABS(PBS(BCPOS(IB),IY,0))
                CS=SQRT(PR(BCPOS(IB),IY)/RO(BCPOS(IB),IY))
!                CS=SQRT(BCINT(PR,NX,NY,IY,IB)/BCINT(RO,NX,NY,IY,IB))
                SNA0(BCPOS(IB),IY,0,IS)=0.0_R8
                SNA0(BCPOS(IB),IY,1,IS)=SNA0(BCPOS(IB),IY,1,IS)+
     &                                  CONPAR(IS,IB,1)*CS*S1
                sna0_no_mdf(BCPOS(IB),IY,:,IS)=
     &	    sna0_no_mdf(BCPOS(IB),IY,:,IS) + SNA0(BCPOS(IB),IY,:,IS)
	          if(mdf_fnb.ne.0) then                                   !srv 13.10.06 {
                 sna0(BCPOS(IB),IY,0,is) = sna0(BCPOS(IB),IY,0,is)+
     &                                  fnaPSch(BCPOS(IB),iy,0,is)
	          endif                                                   !srv 13.10.06 }
              ENDDO
            ELSE
              CALL XERRAB ('B2STBC -- BCCON=14, NO IN RADIAL DIR')
            ENDIF
            GOTO 10000                                                  !srv 14.04.08 }

* -- BCCON=15 -- prescribe the decay length for density for drifts switched on

 1500       CONTINUE                                                  !srv 25.09.09 16.02.12{
            if(ncall.eq.0) then
              write(*,*) 'prescribe the decay length for the density',
     &         ' for drifs switched on',CONPAR(IS,IB,1),' on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'S') THEN                                !srv 23.07.10 {
              DO IX=BCSTART(IB), BCEND(IB)
                csb=sqrt((rza(ix,bcpos(ib),is)*te(ix,bcpos(ib))+
     &                   ti(ix,bcpos(ib)))/(am(is)*mp))
                SNA0(IX,BCPOS(IB),0,IS)= VBIG*
     &           CDNA(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)),1,IS)*
     .           NA(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)),IS)
                SNA0(IX,BCPOS(IB),1,IS)=-VBIG*
     &           (CDNA(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)),1,IS)-
     &            CONPAR(IS,IB,1)*SY(IX,BCPOS(IB))*csb)
              end do
              DO IX=BCSTART(IB), BCEND(IB)
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                csb=sqrt((rza(ix,bcpos(ib),is)*te(ix,bcpos(ib))+
     &                   ti(ix,bcpos(ib)))/(am(is)*mp))
                SNA0(IX,BCPOS(IB),0,IS)= VBIG*CDNA(IX,BCPOS(IB),1,IS)*
     .           NA(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)),IS)
                SNA0(IX,BCPOS(IB),1,IS)= -VBIG*           
     &           (CDNA(IX,BCPOS(IB),1,IS)-CONPAR(IS,IB,1)*
     .            SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))*csb)
              end do
              DO IX=BCSTART(IB), BCEND(IB)                            !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                   !srv 23.09.08 }
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCCON=15, NO IN POLOIDAL DIR.')
            END IF                                                    !srv 23.07.10 }
            GOTO 10000                                                !srv 25.09.09 }

* -- BCCON=16 -- Feedback boundary condition: given total particle flux with constant density  !srv 16.02.12 {

 1600       CONTINUE
            if(ncall.eq.0) then
             call xertst(is.eq.ismain,
     1         'BCCON = 16 is applied for main ion only ')
	       write(*,*) 'given total particle flux with constant',
     &       ' density for main ion is =', ismain
!            write all boundaries with this type of boundary condition
             DO IBW=IB,NBC
		 if(BCCHAR(IB).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN	       
              write(*,*) ' on ', BCCHAR(IBW),' at ',BCPOS(IBW),' from ',
     1         BCSTART(IBW),' to ', BCEND(IBW),' for species ',IS
               call xertst(CONPAR(IS,IBW,1).gt.0.0_R8,
     1         'BCCON = 16, CONPAR(,,1) !> 0')
	       endif
	       ENDDO
#ifdef ASTRA
             SCONPAR = fneb_astra
#else
             SCONPAR = 0.0_R8
             DO IBW=IB,NBC
		 if(BCCHAR(IB).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN	       
                SCONPAR = SCONPAR + CONPAR(IS,IBW,1)
	       endif
	       ENDDO
#endif
		   write(*,*) ' Summary value is ', SCONPAR
	       NAS = NA((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB),is)
            endif
	     if(.not.once_bccon_16) then
            IF(BCCHAR(IB).EQ.'S') THEN
	        totflux = 0.0_R8
            DO IBW=IB,NBC
            IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN
	        if(mdf_fnb.ne.0) then
               DO IX=BCSTART(IBW), BCEND(IBW)
	          totflux = totflux + 
     &          fna_mdf(topix(ix,BCPOS(IBW)),topiy(ix,BCPOS(IBW)),1,is) 
               ENDDO
	        else
               DO IX=BCSTART(IBW), BCEND(IBW)
	          totflux = totflux + 
     &          fna(topix(ix,BCPOS(IBW)),topiy(ix,BCPOS(IBW)),1,is) 
               ENDDO
	        endif
	       ENDIF
	      ENDDO
              NAS = NA((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB),is) * 
     &         (1.0_R8-bc_ref*(totflux-SCONPAR)/
     &           max(abs(totflux),SCONPAR))
!              T0 = (totflux-CONPAR(IS,IB,1))/
!     &           max(abs(totflux),CONPAR(IS,IB,1))
!              NAS = NA((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB),is) * 
!     &         (1.0_R8-max(-bc_ref,min(bc_ref,T0)))
             write(*,"(' BCCON=16: NAS=',es14.7,' totflux=',es14.7,
     *       ' SCONPAR=',es14.7,' IS=',i3)") NAS, totflux, SCONPAR, IS
             DO IBW=IB,NBC
           IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN
               DO IX=BCSTART(IBW), BCEND(IBW)
                SNA0(IX,BCPOS(IBW),0,IS)=VBIG*NAS*SY(IX,BCPOS(IBW))
                SNA0(IX,BCPOS(IBW),1,IS)=-VBIG*SY(IX,BCPOS(IBW))
               ENDDO

               DO IX=BCSTART(IBW), BCEND(IBW)                            !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IBW),:,IS)=SNA0(IX,BCPOS(IBW),:,IS)
               ENDDO                                                     !srv 23.09.08 }
	        ENDIF
	       ENDDO
		   once_bccon_16 = .true.
            ELSE                                                       !srv 16.02.12 }
             CALL XERRAB
     1         ('B2STBC -- BCCON=16, NOT AVAILABLE FOR N,W,E')
            ENDIF
	     endif                                                          !srv 16.02.12 }
            GOTO 10000

* -- BCCON=17 -- FILLER

 1700       CONTINUE
            CALL XERRAB ('B2STBC -- BCCON=17, NOT YET AVAILABLE')
            GOTO 10000

#ifdef ASTRA

* -- BCCON=18 -- PRESCRIBE THE MODIFIED TOTAL PARTICLE FLUX

 1800       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'modyfied total particle flux on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
	        t0 = 0.0_R8
              US=0.0_R8
              DO IX=BCSTART(IB), BCEND(IB)
                US=US+SY(IX,BCPOS(IB))
	          do isb = 0, ns-1
	           if(is_neutral(isb)) then
		      t0 = t0 + fna(ix,topiy(ix,bcpos(ib)),1,isb)
                 endif
                ENDDO
              ENDDO
              CONPAR(IS,IB,1) = (fSeNBI_astra - t0)                     !srv 29.11.11

!              if(ncall.eq.0) then
               write (*,*) ' BCCON=18', ' IS =', IS, 
     &        ' CONPAR=',CONPAR(IS,IB,1), ' fSeNBI_astra=',fSeNBI_astra,
     &        ' integrated neutral flux =',t0
!              endif

              DO IX=BCSTART(IB), BCEND(IB)
                SNA0(IX,BCPOS(IB),0,IS)=SNA0(IX,BCPOS(IB),0,IS)+         !srv 01.02.09
     1           CONPAR(IS,IB,1)/US*SY(IX,BCPOS(IB))
                SNA0(IX,BCPOS(IB),1,IS)=0.0_R8
                sna0_no_mdf(ix,bcpos(ib),0,is) =
     &            sna0_no_mdf(ix,bcpos(ib),0,is)+sna0(ix,bcpos(ib),0,is) !srv 23.09.08
                sna0_no_mdf(ix,bcpos(ib),1,is) =
     &            sna0_no_mdf(ix,bcpos(ib),1,is)+sna0(ix,bcpos(ib),1,is) !srv 23.09.08
	          if(mdf_fnb.ne.0) then                                    !srv 12.09.08 {
                 sna0(ix,bcpos(ib),0,is) = sna0(ix,bcpos(ib),0,is) -
     &                              fnaPSch(ix,topiy(ix,bcpos(ib)),1,is)
                 sna0(ix,bcpos(ib),1,is) = sna0(ix,bcpos(ib),1,is) +
     &            vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &               vaecrb(ix,topiy(ix,bcpos(ib)),1,is)
	          endif                                                    !srv 12.09.08 }
              ENDDO
            ELSE
             CALL XERRAB ('B2STBC -- BCCON=18, ONLY FOR SOUTH')
            ENDIF
            GOTO 10000
#else
* -- BCCON=18 -- FILLER

 1800       CONTINUE
            CALL XERRAB ('B2STBC -- BCCON=18, NOT YET AVAILABLE')
            GOTO 10000
#endif


* -- BCCON=19 -- Feedback boundary condition: given sum of integrated neutrals and main ion particle flux
*                                             with constant ion density

 1900      CONTINUE
           if(ncall.eq.0) then
             call xertst(is.eq.ismain,
     1         'BCCON = 19 is applied for main ion only ')
	       write(*,*) 'given sum of integrated neutrals and main ',
     &                  'ion particle flux with constant ion density '
!            write all boundaries with this type of boundary condition
            DO IBW=IB,NBC
		 if(BCCHAR(IB).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN	       
              write(*,*) ' on ', BCCHAR(IBW),' at ',BCPOS(IBW),' from ',
     1         BCSTART(IBW),' to ', BCEND(IBW),' for species ',IS
               call xertst(CONPAR(IS,IBW,1).gt.0.0_R8,
     1         'BCCON = 19, CONPAR(,,1) !> 0')
	        endif
	      ENDDO
            SCONPAR = 0.0_R8
            DO IBW=IB,NBC
		 if(BCCHAR(IB).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN	       
                SCONPAR = SCONPAR + CONPAR(IS,IBW,1)
	       endif
	      ENDDO
#ifdef ASTRA
            SCONPAR = fSeNBI_astra
#else
            SCONPAR = 0.0_R8
            DO IBW=IB,NBC
		 if(BCCHAR(IB).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN	       
                SCONPAR = SCONPAR + CONPAR(IS,IBW,1)
	       endif
	      ENDDO
#endif
	      write(*,*) ' Summary value is ', SCONPAR
	       NAS = NA((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB),is)
           endif
	    if(.not.once_bccon_19) then
           IF(BCCHAR(IB).EQ.'S') THEN
	         totfluxn = 0.0_R8
	         totfluxi = 0.0_R8
            DO IBW=IB,NBC
            IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN
              DO IX=BCSTART(IBW), BCEND(IBW)
#ifdef USE_EIRENE
                do isa = 1,nnatmi
		     totfluxn = totfluxn + rfluxa(topix(ix,bcpos(ibw))+1,
     &                 topiy(ix,bcpos(ibw))+1,isa,1)*sy(ix,bcpos(ibw))
                enddo
                do ism = 1,nnmoli
                 do isa = 1,nnatmi
		      totfluxn = totfluxn + !  mlcmp(isa,ism)*
     &                 rfluxm(topix(ix,bcpos(ibw))+1,
     &                 topiy(ix,bcpos(ibw))+1,ism,1)*sy(ix,bcpos(ibw))
                 enddo
                enddo
#else
 	          do isb = 0, ns-1
	           if(is_neutral(isb)) then
		      totfluxn = totfluxn + 
     &            fna(topix(ix,BCPOS(IBW)),topiy(ix,bcpos(ibw)),1,isb)
                 endif
                enddo
#endif
	          do isb = 0, ns-1
	           if(isb.eq.ismain) then
	            if(mdf_fnb.ne.0) then
		       totfluxi = totfluxi + fna_mdf(topix(ix,BCPOS(IBW)),
     &                                   topiy(ix,bcpos(ibw)),1,isb)
	            else
		       totfluxi = totfluxi + 
     &             fna(topix(ix,BCPOS(IBW)),topiy(ix,bcpos(ibw)),1,isb)
	            endif
                 endif
                enddo
              ENDDO
	      ENDIF
	      ENDDO
            NAS = NA((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB),is) * 
     &         (1.0_R8-bc_ref*(totfluxn+totfluxi-SCONPAR)/
     &           max(abs(totfluxn+totfluxi),SCONPAR))
             write(*,"(' BCCON=19: NAS=',es14.7,' totfluxn=',es14.7,
     &        ' totfluxi=',es14.7,' SCONPAR=',es14.7,
     &        ' totfluxn+totfluxi=',es14.7)")
     *      NAS, totfluxn, totfluxi, SCONPAR, totfluxn+totfluxi
           DO IBW=IB,NBC
            IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN
               DO IX=BCSTART(IBW), BCEND(IBW)                              !srv 01.11.12
                SNA0(IX,BCPOS(IBW),0,IS)=VBIG*NAS*SY(IX,BCPOS(IBW))
                SNA0(IX,BCPOS(IBW),1,IS)=-VBIG*SY(IX,BCPOS(IBW))
               ENDDO

               DO IX=BCSTART(IBW), BCEND(IBW)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IBW),:,IS)=SNA0(IX,BCPOS(IBW),:,IS)
               ENDDO                                                     !srv 23.09.08 }
	      ENDIF
	     ENDDO
		   once_bccon_19 = .true.
           ELSE                                                       !srv 16.02.12 }
             CALL XERRAB
     1         ('B2STBC -- BCCON=19, NOT AVAILABLE FOR N,W,E')
           ENDIF
	    endif                                                          !srv 16.02.12 }
          GOTO 10000

* -- BCCON=20 -- constant density feedback condition

 2000       CONTINUE
            if(ncall.eq.0) then
	      write(*,*) 'constant density (neo)',
     1	       CONPAR(IS,IB,1),' on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS
              call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1         'BCCON = 1, CONPAR(,,1) !> 0')
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
	      naave = 0.0_R8
	      volave = 0.0_R8
              DO IX=BCSTART(IB), BCEND(IB)
		bc_type20_ix=ix
		bc_type20_iy=BCPOS(IB)
		do bc_type20_i=1,bc_type20_ref      ! how many steps away from the boundary?
		  bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
		  bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
		  bc_type20_ix=bc_type20_ixx
		  call xertst (-2.lt.bc_type20_ix.and.bc_type20_ix.lt.nx+1,
     1		   'stepped off x-grid in type 20 bc')
		  call xertst (-2.lt.bc_type20_iy.and.bc_type20_iy.lt.ny+1,
     1		   'stepped off y-grid in type 20 bc')
		enddo
		naave=naave+na(bc_type20_ix,bc_type20_iy,is)*
     &                 vol(bc_type20_ix,bc_type20_iy)
		volave=volave+vol(bc_type20_ix,bc_type20_iy)
              end do
	      naave=naave/volave
              DO IX=BCSTART(IB), BCEND(IB)
		bc_type20_ix=ix
		bc_type20_iy=BCPOS(IB)
		do bc_type20_i=1,bc_type20_ref      ! how many steps away from the boundary?
		  bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
		  bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
		  bc_type20_ix=bc_type20_ixx
		  call xertst (-2.lt.bc_type20_ix.and.bc_type20_ix.lt.nx+1,
     1		   'stepped off x-grid in type 20 bc')
		  call xertst (-2.lt.bc_type20_iy.and.bc_type20_iy.lt.ny+1,
     1		   'stepped off y-grid in type 20 bc')
		enddo
                SNA0(IX,BCPOS(IB),0,IS)= VBIG*CONPAR(IS,IB,1)*
     .           SY(IX,BCPOS(IB))*na(bc_type20_ix,bc_type20_iy,is)/naave
                SNA0(IX,BCPOS(IB),1,IS)=-VBIG*SY(IX,BCPOS(IB))
              end do
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO                                                     !srv 23.09.08 }
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCCON=20, NOT AVAILABLE FOR N,W,E')
            ENDIF
            GOTO 10000

* -- BCCON=21 -- PRESCRIBE THE VALUE OF THE DENSITY  !lk and IYS 01.03.2013{
* AND ADD A DENSITY PERTURBATION TO GET A SOLUTION
* WHICH IS AS CLOSE AS POSSIBLE TO NEOCLASSICAL THEORY

* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH
* CORRESPONDING CONDITION ON IONS TEMPERATURE

 2100       CONTINUE
            if(ncall.eq.0) then
#ifdef ASTRA
             IF(BCCHAR(IB).EQ.'S') THEN
	        if(is.eq.ismain) then
	         CONPAR(IS,IB,1) = neb_astra                               !srv 28.07.11
	        endif
             ENDIF
#endif
              write(*,*) 'constant density',CONPAR(IS,IB,1),' on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS
              call xertst(CONPAR(IS,IB,1).gt.0.0_R8,
     1         'BCCON = 21, CONPAR(,,1) !> 0')
            endif
            IF(BCCHAR(IB).EQ.'S' .and. 
     $     (mod(region(BCSTART(IB),bcpos(ib),0),4).eq.1)) THEN
!              DO IX=BCSTART(IB), BCEND(IB)
!                SNA0(IX,BCPOS(IB),0,IS)= VBIG*CONPAR(IS,IB,1)*
!     .           SY(IX,BCPOS(IB))
!                SNA0(IX,BCPOS(IB),1,IS)=-VBIG*SY(IX,BCPOS(IB))
!              end do
	        t0=0.0e0_R8   ! Volume between two surfaces                  !lk 12.12.12 {              
	        t1=0.0e0_R8   ! Sum of volume averaged electrons ans ions temperatures
              do ix = -1,nx                                               
                if (mod(region(ix,bcpos(ib),0),4).eq.1) then                      
	           t0=t0 + vol(ix,0)                        
	           t1=t1 + (ti(ix,-1)+te(ix,-1))*vol(ix,0)                        
                endif                                                      
	      enddo
	      t1=t1/t0        ! Sum of volume averaged electrons ans ions temperatures
              do IX=BCSTART(IB), BCEND(IB)                                              	                                                       
		if(mod(region(ix,bcpos(ib),0),4).eq.1) then
                   SNA0(IX,BCPOS(IB),0,IS)= VBIG*CONPAR(IS,IB,1)*t1
     &			   /(ti(ix,-1)+te(ix,-1))*SY(IX,BCPOS(IB))
                   SNA0(IX,BCPOS(IB),1,IS)=-VBIG*SY(IX,BCPOS(IB))
	        endif
	      enddo
                                                     !lk 12.12.12 }
              DO IX=BCSTART(IB), BCEND(IB)                              !srv 23.09.08 
                SNA0_NO_MDF(IX,BCPOS(IB),:,IS)=SNA0(IX,BCPOS(IB),:,IS)
              ENDDO
           ELSE
              CALL XERRAB
     1         ('B2STBC -- BCCON=21, VALID FOR CORE ONLY')
            ENDIF
            GOTO 10000 !lk and IYS 01.03.2013}

* -- BCCON=22 -- Feedback boundary condition: given total particle flux with constant 
*    average density. The poloidal density variation is added to make the solution
*    as close as possible to the neoclassical value. It is recommended to use this boundary
*    condition with corresponding boundary condition on ions temperature  !srv 16.02.12, IYS 01.03.2013 {

 2200       CONTINUE
            if(ncall.eq.0) then
             call xertst(is.eq.ismain,
     1         'BCCON = 23 is applied for main ion only ')
	       write(*,*) 'given total particle flux with constant',
     &       ' density for main ion is =', ismain
!            write all boundaries with this type of boundary condition
             DO IBW=IB,NBC
	   if(BCCHAR(IB).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB).and. 
     $     (mod(region(BCSTART(IB),bcpos(ib),0),4).eq.1)) THEN	       
              write(*,*) ' on ', BCCHAR(IBW),' at ',BCPOS(IBW),' from ',
     1         BCSTART(IBW),' to ', BCEND(IBW),' for species ',IS
               call xertst(CONPAR(IS,IBW,1).gt.0.0_R8,
     1         'BCCON = 22, CONPAR(,,1) !> 0')
	       endif
	       ENDDO
#ifdef ASTRA
             SCONPAR = fneb_astra
#else
             SCONPAR = 0.0_R8
             DO IBW=IB,NBC
		 if(BCCHAR(IB).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN	       
                SCONPAR = SCONPAR + CONPAR(IS,IBW,1)
	       endif
	       ENDDO
#endif
	     write(*,*) ' Summary value is ', SCONPAR

	     NAS = NA((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB),is)

             endif
	     if(.not.once_bccon_16) then
            IF(BCCHAR(IB).EQ.'S') THEN
	        totflux = 0.0_R8
	        t0=0.0e0_R8   ! Volume between two surfaces                  !lk 12.12.12 {              
	        t1=0.0e0_R8   ! Sum of volume averaged electrons ans ions temperatures
                NAS=0.0e0_R8
            DO IBW=IB,NBC
            IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN
	        if(mdf_fnb.ne.0) then
               DO IX=BCSTART(IBW), BCEND(IBW)
	          totflux = totflux + 
     &          fna_mdf(topix(ix,BCPOS(IBW)),topiy(ix,BCPOS(IBW)),1,is) 
	           t0=t0 + vol(ix,0)                        
	           t1=t1 + (ti(ix,-1)+te(ix,-1))*vol(ix,0)
                   NAS=NAS+na(ix,-1,is)*vol(ix,0)                        
               ENDDO
	        else
               DO IX=BCSTART(IBW), BCEND(IBW)
	          totflux = totflux + 
     &          fna(topix(ix,BCPOS(IBW)),topiy(ix,BCPOS(IBW)),1,is) 
	           t0=t0 + vol(ix,0)                        
	           t1=t1 + (ti(ix,-1)+te(ix,-1))*vol(ix,0)                        
                   NAS=NAS+na(ix,-1,is)*vol(ix,0)                        
               ENDDO
	        endif
	       ENDIF
	      ENDDO
	      t1=t1/t0        ! Sum of volume averaged electrons ans ions temperatures
              NAS=NAS/t0      ! Volume averaged density

              NAS = NAS* !NA((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB),is) * 
     &         (1.0_R8-bc_ref*(totflux-SCONPAR)/
     &           max(abs(totflux),SCONPAR))
!              T0 = (totflux-CONPAR(IS,IB,1))/
!     &           max(abs(totflux),CONPAR(IS,IB,1))
!              NAS = NA((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB),is) * 
!     &         (1.0_R8-max(-bc_ref,min(bc_ref,T0)))
             write(*,"(' BCCON=22: NAS=',es14.7,' totflux=',es14.7,
     *       ' SCONPAR=',es14.7,' IS=',i3)") NAS, totflux, SCONPAR, IS
             DO IBW=IB,NBC
           IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN
               DO IX=BCSTART(IBW), BCEND(IBW)
                SNA0(IX,BCPOS(IBW),0,IS)=VBIG*NAS*t1
     &			   /(ti(ix,-1)+te(ix,-1))*SY(IX,BCPOS(IBW))  ! lk
                SNA0(IX,BCPOS(IBW),1,IS)=-VBIG*SY(IX,BCPOS(IBW))
               ENDDO

               DO IX=BCSTART(IBW), BCEND(IBW)                            !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IBW),:,IS)=SNA0(IX,BCPOS(IBW),:,IS)
               ENDDO                                                     !srv 23.09.08 }
	        ENDIF
	       ENDDO
		   once_bccon_16 = .true.
            ELSE                                                       !srv 16.02.12 }
             CALL XERRAB
     1         ('B2STBC -- BCCON=22 IS AVAILABLE FOR CORE ONLY')
            ENDIF
	     endif                                                          !srv 16.02.12 }
            GOTO 10000                            !IYS 01.03.2013 }

* -- BCCON=23 -- Feedback boundary condition: given total particle flux with constant 
*    average density. The poloidal density variation is added to make the solution
*    as close as possible to the neoclassical value. It is recommended to use this boundary
*    condition with corresponding boundary condition on ions temperature  !srv 16.02.12, IYS 04.03.2013 {

 2300       CONTINUE
            if(ncall.eq.0) then
             call xertst(is.eq.ismain,
     1         'BCCON = 23 is applied for main ion only ')
	       write(*,*) 'given total particle flux with constant',
     &       ' density for main ion is =', ismain
!            write all boundaries with this type of boundary condition
             DO IBW=IB,NBC
		 if(BCCHAR(IB).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN	       
              write(*,*) ' on ', BCCHAR(IBW),' at ',BCPOS(IBW),' from ',
     1         BCSTART(IBW),' to ', BCEND(IBW),' for species ',IS
               call xertst(CONPAR(IS,IBW,1).gt.0.0_R8,
     1         'BCCON = 23, CONPAR(,,1) !> 0')
	       endif
	       ENDDO
#ifdef ASTRA
            SCONPAR = fSeNBI_astra
#else
            SCONPAR = 0.0_R8
            DO IBW=IB,NBC
		 if(BCCHAR(IB).EQ.'S'.and.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN	       
                SCONPAR = SCONPAR + CONPAR(IS,IBW,1)
	       endif
	      ENDDO
#endif
	     write(*,*) ' Summary value is ', SCONPAR

	     NAS = NA((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB),is)

             endif
	     if(.not.once_bccon_16) then
            IF(BCCHAR(IB).EQ.'S') THEN
	        totflux = 0.0_R8
	        t0=0.0e0_R8   ! Volume between two surfaces                  !lk 12.12.12 {              
	        t1=0.0e0_R8   ! Sum of volume averaged electrons ans ions temperatures

                NAS=0.0e0_R8
            DO IBW=IB,NBC
            IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN
              DO IX=BCSTART(IBW), BCEND(IBW)
#ifdef USE_EIRENE
                do isa = 1,nnatmi
		     totfluxn = totfluxn + rfluxa(topix(ix,bcpos(ibw))+1,
     &                 topiy(ix,bcpos(ibw))+1,isa,1)*sy(ix,bcpos(ibw))
                enddo
                do ism = 1,nnmoli
                 do isa = 1,nnatmi
		      totfluxn = totfluxn + !  mlcmp(isa,ism)*
     &                 rfluxm(topix(ix,bcpos(ibw))+1,
     &                 topiy(ix,bcpos(ibw))+1,ism,1)*sy(ix,bcpos(ibw))
                 enddo
                enddo
#else
 	          do isb = 0, ns-1  ! incorrect for multispecies plasma; IYS
	           if(is_neutral(isb)) then
		      totfluxn = totfluxn + 
     &            fna(topix(ix,BCPOS(IBW)),topiy(ix,bcpos(ibw)),1,isb)
                 endif
                enddo
#endif
	          do isb = 0, ns-1  ! incorrect for multispecies plasma; IYS
	           if(isb.eq.ismain) then
	            if(mdf_fnb.ne.0) then
		       totfluxi = totfluxi + fna_mdf(topix(ix,BCPOS(IBW)),
     &                                   topiy(ix,bcpos(ibw)),1,isb)
	            else
		       totfluxi = totfluxi + 
     &             fna(topix(ix,BCPOS(IBW)),topiy(ix,bcpos(ibw)),1,isb)
	            endif

                 endif
                enddo
	           t0=t0 + vol(ix,0)                        
	           t1=t1 + (ti(ix,-1)+te(ix,-1))*vol(ix,0)
                   NAS=NAS+na(ix,-1,is)*vol(ix,0)                        
              ENDDO
	      ENDIF
	      ENDDO
	      t1=t1/t0        ! Sum of volume averaged electrons ans ions temperatures
              NAS=NAS/t0      ! Volume averaged density

              NAS = NAS* 
     &         (1.0_R8-bc_ref*(totfluxn+totfluxi-SCONPAR)/
     &           (abs(totfluxn)+abs(totfluxi)+abs(SCONPAR)))
             write(*,"(' BCCON=23: NAS=',es14.7,' totfluxn=',es14.7,
     &        ' totfluxi=',es14.7,' SCONPAR=',es14.7,
     &        ' totfluxn+totfluxi=',es14.7)")
     *      NAS, totfluxn, totfluxi, SCONPAR, totfluxn+totfluxi
             DO IBW=IB,NBC
           IF(BCCHAR(IBW).EQ.'S'.AND.BCCON(IS,IBW).EQ.BCCON(IS,IB))THEN
               DO IX=BCSTART(IBW), BCEND(IBW)
                SNA0(IX,BCPOS(IBW),0,IS)=VBIG*NAS*t1
     &			   /(ti(ix,-1)+te(ix,-1))*SY(IX,BCPOS(IBW))  ! lk
                SNA0(IX,BCPOS(IBW),1,IS)=-VBIG*SY(IX,BCPOS(IBW))
               ENDDO

               DO IX=BCSTART(IBW), BCEND(IBW)                            !srv 23.09.08 {
                SNA0_NO_MDF(IX,BCPOS(IBW),:,IS)=SNA0(IX,BCPOS(IBW),:,IS)
               ENDDO                                                     !srv 23.09.08 }
	        ENDIF
	       ENDDO
		   once_bccon_16 = .true.
            ELSE                                                       !srv 16.02.12 }
             CALL XERRAB
     1         ('B2STBC -- BCCON=23 IS AVAILABLE FOR CORE ONLY')
            ENDIF
	     endif                                                          !srv 16.02.12 }
            GOTO 10000                            !IYS 04.03.2013 }

* -- BOUNDARY CONDITIONS FOR MOMENTUM EQUATIONS --

10000       continue
            if(BCMOM(IS,IB).eq.0) goto 20000
	    GOTO(10100,10200,10300,10400,10500,
     1	     10600,10700,10800,10900,11000,
     1         11100,11200,11300,11400,11500,11600),BCMOM(IS,IB)      !srv 01.02.09 19.04.13
            write(*,*) BCMOM(IS,IB)
            call xerrab ('BCMOM out of bounds')

* -- BCMOM=1 -- PRESCRIBE THE VALUE OF THE PARALLEL VELOCITY

10100	    CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'specified parallel velocity',MOMPAR(IS,IB,1),
     1         ' on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2         ' to ',BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SMO0(IX,BCPOS(IB),0,IS)= VBIG*NE00*AM(IS)*
     .           MOMPAR(IS,IB,1)*SYHZ(IX,BCPOS(IB))                   !srv 24.08.08
                SMO0(IX,BCPOS(IB),1,IS)=-VBIG*NE00*AM(IS)*
     .           SYHZ(IX,BCPOS(IB))                                   !srv 24.08.08
              end do
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SMO0(IX,BCPOS(IB),0,IS)= VBIG*NE00*AM(IS)*
     .           MOMPAR(IS,IB,1)*
     .           SYHZ(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))  !srv 24.08.08
                SMO0(IX,BCPOS(IB),1,IS)=-VBIG*NE00*AM(IS)*
     .           SYHZ(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))  !srv 24.08.08
              end do
            ELSEIF(BCCHAR(IB).EQ.'W') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SMO0(BCPOS(IB),IY,0,IS)= UBIG*NE00*AM(IS)*
     .           MOMPAR(IS,IB,1)*SXHZ(BCPOS(IB),IY)                   !srv 24.08.08
                SMO0(BCPOS(IB),IY,1,IS)=-UBIG*NE00*AM(IS)*
     .           SXHZ(BCPOS(IB),IY)                   !srv 24.08.08
              end do
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SMO0(BCPOS(IB),IY,0,IS)= UBIG*NE00*AM(IS)*
     .           MOMPAR(IS,IB,1)*
     .           SXHZ(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))      !srv 24.08.08
                SMO0(BCPOS(IB),IY,1,IS)=-UBIG*NE00*AM(IS)*
     .           SXHZ(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))      !srv 24.08.08
              end do
            ELSEIF(BCCHAR(IB).EQ.'X') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SMO0(IX,BCPOS(IB),0,IS)= UBIG*NE00*AM(IS)*
     .           MOMPAR(IS,IB,1)*VOL(IX,BCPOS(IB))
                SMO0(IX,BCPOS(IB),1,IS)=-UBIG*NE00*AM(IS)*
     .           VOL(IX,BCPOS(IB))
              end do
            ELSEIF(BCCHAR(IB).EQ.'Y') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SMO0(BCPOS(IB),IY,0,IS)= UBIG*NE00*AM(IS)*
     .           MOMPAR(IS,IB,1)*VOL(BCPOS(IB),IY)
                SMO0(BCPOS(IB),IY,1,IS)=-UBIG*NE00*AM(IS)*
     .           VOL(BCPOS(IB),IY)
              end do
            ELSE
              CALL XERRAB ('B2STBC -- BCMOM=1, INVALID OPTION.')
            ENDIF
            GOTO 20000

* -- BCMOM=2 -- PRESCRIBE THE GRADIENT OF THE PARALLEL VELOCITY

10200       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'parallel velocity gradient',MOMPAR(IS,IB,1),
     1         ' on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2         ' to ',BCEND(IB),' for species ',IS
            endif
!#ifdef WANT_THIS                                                     !srv 21.03.09
            IF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SMO0(IX,BCPOS(IB),0,IS)= VBIG*
     1           SYHZ(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))* !srv 24.08.08
     2           AM(IS)*na(IX,BCPOS(IB),IS)*
     .           (MOMPAR(IS,IB,1)*
     1           (HY1(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))+
     1            HY1(IX,BCPOS(IB)))/2 + 
     2           UA(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)),IS))
                SMO0(IX,BCPOS(IB),1,IS)=-VBIG*
     .           SYHZ(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))* !srv 24.08.08
     .           AM(IS)*na(IX,BCPOS(IB),IS)
              end do
            ELSEIF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SMO0(IX,BCPOS(IB),0,IS)=VBIG*SYHZ(IX,BCPOS(IB))*      !srv 24.08.08
     .           AM(IS)*na(IX,BCPOS(IB),IS) *
     .           (-MOMPAR(IS,IB,1)*
     1           (HY1(IX,BCPOS(IB))+
     2            HY1(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB))))/2 +
     .           UA(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)),IS))
                SMO0(IX,BCPOS(IB),1,IS)=- VBIG*SYHZ(IX,BCPOS(IB)) *   !srv 24.08.08
     .           AM(IS) *
     .           na(IX,BCPOS(IB),IS)
              end do
            ELSEIF(BCCHAR(IB).EQ.'W') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SMO0(BCPOS(IB),IY,0,IS)=VBIG*SXHZ(BCPOS(IB),IY)*      !srv 24.08.08
     .           AM(IS)*na(BCPOS(IB),IY,IS) *
     .           (-MOMPAR(IS,IB,1)*
     1           (HX(BCPOS(IB),IY)+
     2            HX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)))/2 +
     .           UA(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY),IS))
                SMO0(BCPOS(IB),IY,1,IS)=- VBIG*SXHZ(BCPOS(IB),IY) *   !srv 24.08.08
     .           AM(IS) *
     .           na(BCPOS(IB),IY,IS)
              end do
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SMO0(BCPOS(IB),IY,0,IS)= VBIG*
     .           SXHZ(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))*AM(IS)*  !srv 24.08.08
     .           na(BCPOS(IB),IY,IS) *
     .           (MOMPAR(IS,IB,1)*
     1           (HX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))+
     2            HX(BCPOS(IB),IY))/2 +
     .           UA(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY),IS))
                SMO0(BCPOS(IB),IY,1,IS)=-VBIG*
     .           SXHZ(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)) *    !srv 24.08.08
     .           AM(IS) *
     .           na(BCPOS(IB),IY,IS)
              end do
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCMOM=2, INVALID OPTION.')
            END IF
!#else                                                                !srv 21.03.09
!            CALL XERRAB ('B2STBC -- BCMOM=2, NOT IMPLEMENTED.')
!#endif                                                               !srv 21.03.09
            GOTO 20000

* -- BCMOM=3 -- SHEATH CONDITIONS, MACH NUMBER AS INPUT

10300       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'sheath velocity',MOMPAR(IS,IB,2),' on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'W') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                IF(MOMPAR(IS,IB,2).LT.0.5_R8) THEN
                  if(sound_model.eq.0) then
                    CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))  !srv 01.02.09
                  else
                    zarg = 0.0_R8
                    do ix = 0, ns-1
                       zarg = zarg + bcint(na(-1,-1,ix),nx,ny,IY,IB)*
     .                        (bcint(rza(-1,-1,ix),nx,ny,IY,IB)*
     .                         bcint(te,nx,ny,IY,IB)+
     .                         bcint(ti,nx,ny,IY,IB))
     .                        /(am(ix)*
     .                        bcint(ni(-1,-1,0),nx,ny,IY,IB)*mp)
                    enddo
                    CS=SQRT(zarg)
                  endif
                  UUOUT=MOMPAR(IS,IB,1)*CS*
     .             abs(wbbr(BCPOS(IB),IY,0)/wbbr(BCPOS(IB),IY,3))
                  if(.not.is_neutral(is)) then                           !xpb
                    UPOUT=max(0.0_R8, UUOUT/
     .               abs(wbbr(BCPOS(IB),IY,0)/wbbr(BCPOS(IB),IY,3)) +    !xpb
     .               bcint(uadia(-1,-1,0,is),nx,ny,iy,ib)*
     .               abs(wbbr(bcpos(ib),iy,2)/wbbr(bcpos(ib),iy,0)))     !xpb
                  else                                                   !xpb
                    UPOUT = UUOUT/
     .               abs(wbbr(BCPOS(IB),IY,0)/wbbr(BCPOS(IB),IY,3))      !xpb
                  endif                                                  !xpb
cxpb  Added multiplicative factor to UADIA according to the Chodura criterion
cxpb  and max(,,) statement for stability
                ELSE
                  ZBAR=bcint(na(-1,-1,is),nx,ny,iy,ib)/
     .             bcint(NE,nx,ny,IY,IB)
                  ZARG=(GAMMAI*bcint(ti,nx,ny,iy,ib)
     .             +bcint(rz2(-1,-1,is),nx,ny,iy,ib)*ZBAR
     .             *bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                  CS=SQRT(ZARG)
                  UUTEST=
     .             ua(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY),IS)*
     .             (1.0_R8+
     .              HX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))/
     .      HX(RIGHTIX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),
     .         RIGHTIY(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))))-
     .      ua(RIGHTIX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),
     .         RIGHTIY(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),IS)*
     .              HX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))/
     .      HX(RIGHTIX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),
     .         RIGHTIY(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)))
                  UUOUT=MAX(MOMPAR(IS,IB,1)*CS,-UUTEST)
                  if(.not.is_neutral(is)) then                           !xpb
                    UPOUT=max(0.0_R8, UUOUT+                             !xpb
     .               bcint(uadia(-1,-1,0,is),nx,ny,iy,ib)*
     .               abs(wbbr(bcpos(ib),iy,2)/wbbr(bcpos(ib),iy,0)))     !xpb
                  else                                                   !xpb
                    UPOUT = UUOUT                                        !xpb
                  endif                                                  !xpb
cxpb  Added multiplicative factor to UADIA according to the Chodura criterion
cxpb  and max(,,) statement for stability
                ENDIF
                SMO0(BCPOS(IB),IY,0,IS)= (UBIG/EPSLON**2)*AM(IS)*mp*
     1           bcint(na(-1,-1,IS),nx,ny,IY,IB)*SXHZ(BCPOS(IB),IY)*  !srv 24.08.08
     2           (-UPOUT)
                SMO0(BCPOS(IB),IY,1,IS)=-(UBIG/EPSLON**2)*AM(IS)*mp*
     .           bcint(na(-1,-1,IS),nx,ny,IY,IB)*SXHZ(BCPOS(IB),IY)   !srv 24.08.08
              end do
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                IF(MOMPAR(IS,IB,2).LT.0.5_R8) THEN
                  if(sound_model.eq.0) then
                    CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB)) !srv 01.02.09
                  else
                    zarg = 0.0_R8
                    do ix = 0, ns-1
                      zarg = zarg + bcint(na(-1,-1,ix),nx,ny,IY,IB)*
     .                        (bcint(rza(-1,-1,ix),nx,ny,IY,IB)*
     .                         bcint(te,nx,ny,IY,IB)+
     .                         bcint(ti,nx,ny,IY,IB))
     .                        /(am(ix)*
     .                        bcint(ni(-1,-1,0),nx,ny,IY,IB)*mp)
                    enddo
                    CS=SQRT(zarg)
                  endif
                  UUOUT=MOMPAR(IS,IB,1)*CS*
     .             abs(wbbl(BCPOS(IB),IY,0)/wbbl(BCPOS(IB),IY,3))
                  if(.not.is_neutral(is)) then                           !xpb
                    UPOUT=max(0.0_R8, UUOUT/                             !xpb
     .                abs(wbbl(BCPOS(IB),IY,0)/wbbl(BCPOS(IB),IY,3))     !xpb
     .               -bcint(uadia(-1,-1,0,is),nx,ny,iy,ib)*
     .                abs(wbbl(bcpos(ib),iy,2)/wbbl(bcpos(ib),iy,0)))    !xpb
                  else                                                   !xpb
                    UPOUT = UUOUT/
     .               abs(wbbl(BCPOS(IB),IY,0)/wbbl(BCPOS(IB),IY,3))      !xpb
                  endif                                                  !xpb
cxpb  Added multiplicative factor to UADIA according to the Chodura criterion
cxpb  and max(,,) statement for stability
                ELSE
                  ZBAR=bcint(na(-1,-1,is),nx,ny,iy,ib)/
     .             bcint(NE,nx,ny,IY,IB)
                  ZARG=(GAMMAI*bcint(ti,nx,ny,iy,ib)
     .             +bcint(rz2(-1,-1,is),nx,ny,iy,ib)*ZBAR
     .             *bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                  CS=SQRT(ZARG)
                  UUTEST=
     .             ua(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY),IS)*
     .             (1.0_R8+
     .              HX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))/
     .         HX(LEFTIX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),
     .            LEFTIY(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))))-
     .         ua(LEFTIX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),
     .            LEFTIY(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),IS)*
     .              HX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))/
     .         HX(LEFTIX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),
     .            LEFTIY(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)))
                  UUOUT=MAX(MOMPAR(IS,IB,1)*CS,UUTEST)
                  if(.not.is_neutral(is)) then                           !xpb
                    UPOUT=max(0.0_R8, UUOUT -                            !xpb
     .               bcint(uadia(-1,-1,0,is),nx,ny,iy,ib)*
     .               abs(wbbl(bcpos(ib),iy,2)/wbbl(bcpos(ib),iy,0)))     !xpb
                  else                                                   !xpb
                    UPOUT = UUOUT                                        !xpb
                  endif                                                  !xpb
cxpb  Added multiplicative factor to UADIA according to the Chodura criterion
cxpb  and max(,,) statement for stability
                ENDIF
                SMO0(BCPOS(IB),IY,0,IS)= (UBIG/EPSLON**2)*AM(IS)*mp*
     1           bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     2           SXHZ(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))*UPOUT  !srv 24.08.08
                SMO0(BCPOS(IB),IY,1,IS)=-(UBIG/EPSLON**2)*AM(IS)*mp*
     .           bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .           SXHZ(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))      !srv 24.08.08
              end do
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCMOM=3, NO SHEATH IN RADIAL DIR.')
            END IF
            GOTO 20000

* -- BCMOM=4 -- PRESCRIBE THE VALUE OF THE VELOCITY, WEAKLY
*                 A MIXED BOUNDARY CONDITION

10400       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'weak velocity with mompar1 =',MOMPAR(IS,IB,1),
     1         ' and mompar2 =', MOMPAR(IS,IB,2),' on ',
     2         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     3         BCEND(IB),' for species ',IS
              call xertst(MOMPAR(IS,IB,2).gt.0.0_R8,
     1         'BCMOM = 4, MOMPAR(,,2) !> 0')
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SMO0(IX,BCPOS(IB),0,IS)= MOMPAR(IS,IB,2)
     .           *MOMPAR(IS,IB,1)*SYHZ(IX,BCPOS(IB))                  !srv 24.08.08
                SMO0(IX,BCPOS(IB),1,IS)=-MOMPAR(IS,IB,2)
     .           *SYHZ(IX,BCPOS(IB))                                  !srv 24.08.08
              end do
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SMO0(IX,BCPOS(IB),0,IS)= MOMPAR(IS,IB,2)*
     .           MOMPAR(IS,IB,1)*
     .           SYHZ(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))  !srv 24.08.08
                SMO0(IX,BCPOS(IB),1,IS)=-MOMPAR(IS,IB,2)*
     1           SYHZ(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))  !srv 24.08.08
              end do
            ELSEIF(BCCHAR(IB).EQ.'W') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SMO0(BCPOS(IB),IY,0,IS)= MOMPAR(IS,IB,2)
     .           *MOMPAR(IS,IB,1)*SXHZ(BCPOS(IB),IY)                  !srv 24.08.08
                SMO0(BCPOS(IB),IY,1,IS)=
     1           -MOMPAR(IS,IB,2)*SXHZ(BCPOS(IB),IY)                  !srv 24.08.08
              end do
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SMO0(BCPOS(IB),IY,0,IS)= MOMPAR(IS,IB,2)*
     .           MOMPAR(IS,IB,1)*
     .           SXHZ(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))      !srv 24.08.08
                SMO0(BCPOS(IB),IY,1,IS)=-MOMPAR(IS,IB,2)*
     .           SXHZ(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))      !srv 24.08.08
              end do
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCMOM=4, OPTION NOT AVAILABLE')
            ENDIF
            GOTO 20000

* -- BCMOM=5 -- PRESCRIBE THE MOMENTUM FLUX PER UNIT AREA

10500       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'momentum flux per unit area',MOMPAR(IS,IB,1),
     1         ' on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2         ' to ',BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
             DO IX=BCSTART(IB), BCEND(IB)
              SMO0(IX,BCPOS(IB),0,IS)=MOMPAR(IS,IB,1)*SYHZ(IX,BCPOS(IB)) !srv 24.08.08
              SMO0(IX,BCPOS(IB),1,IS)=0.0_R8
             end do
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SMO0(IX,BCPOS(IB),0,IS)=MOMPAR(IS,IB,1)*
     1           SYHZ(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB))) !srv 24.08.08
                SMO0(IX,BCPOS(IB),1,IS)=0.0_R8
              end do
            ELSEIF(BCCHAR(IB).EQ.'W') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SMO0(BCPOS(IB),IY,0,IS)=
     1           MOMPAR(IS,IB,1)*SXHZ(BCPOS(IB),IY)                  !srv 24.08.08
                SMO0(BCPOS(IB),IY,1,IS)=0.0_R8
              end do
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              DO IY=BCSTART(IB), BCEND(IB)
                SMO0(BCPOS(IB),IY,0,IS)=MOMPAR(IS,IB,1)*
     1           SXHZ(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))     !srv 24.08.08
                SMO0(BCPOS(IB),IY,1,IS)=0.0_R8
              end do
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCMOM=5, OPTION NOT AVAILABLE')
            ENDIF
            GOTO 20000

* -- BCMOM=6 -- PRESCRIBE THE TOTAL MOMENTUM FLUX FOR A
*                 CONSTANT PARALLEL VELOCITY

10600       CONTINUE
            CALL XERRAB ('B2STBC -- BCMOM=6, NOT YET AVAILABLE')
            GOTO 20000

* -- BCMOM=7 -- PRESCRIBE THE MOMENTUM AS A FUNCTION OF OTHER
*                 PLASMA PARAMETERS

10700       CONTINUE
            CALL XERRAB ('B2STBC -- BCMOM=7, NOT YET AVAILABLE')
            GOTO 20000

* -- BCMOM=8 -- SPECIAL : LIMITED SHEAR

10800       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'limited shear velocity',MOMPAR(IS,IB,1),
     1         ' on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2         ' to ',BCEND(IB),' for species ',IS
            endif
#ifdef WANT_THIS
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                A=MIN(1.0_R8,
     .           SQRT(TI(IX,BCPOS(IB))/
     .                TI(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)))))
                SMO0(IX,BCPOS(IB),0,IS)=VBIG*SYHZ(IX,BCPOS(IB))*      !srv 24.08.08
     .           AM(IS)*na(IX,BCPOS(IB),IS)*A*
     .           UA(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)),IS)
                SMO0(IX,BCPOS(IB),1,IS)=-VBIG*SYHZ(IX,BCPOS(IB))*     !srv 24.08.08
     .           AM(IS)*na(IX,BCPOS(IB),IS)
              end do
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                A=MIN(1.0_R8,
     .           SQRT(TI(IX,BCPOS(IB))/
     .            TI(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))))
                SMO0(IX,BCPOS(IB),0,IS)= VBIG*
     .           SYHZ(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))* !srv 24.08.08
     .           AM(IS) * na(IX,BCPOS(IB),IS)*A*
     .           UA(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)),IS)
                SMO0(IX,BCPOS(IB),1,IS)=-VBIG*
     .           SYHZ(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB))) !srv 24.08.08
     .           *AM(IS) * na(IX,BCPOS(IB),IS)
              ENDDO
            ELSE
              CALL XERRAB ('B2STBC -- BCMOM=8, NOT AVAILABLE FOR W,E')
            END IF
#else
            CALL XERRAB ('B2STBC -- BCMOM=8, NOT IMPLEMENTED.')
#endif
            GOTO 20000

* -- BCMOM=9 -- PRESCRIBE THE TOTAL MOMENTUM FLUX WITH
*                 CONSTANT FLUX DENSITY

10900       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'total momentum flux',MOMPAR(IS,IB,1),' on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              US=0.0_R8
              DO IX=BCSTART(IB), BCEND(IB)
                US=US+SY(IX,BCPOS(IB))
              ENDDO
              DO IX=BCSTART(IB), BCEND(IB)
                SMO0(IX,BCPOS(IB),0,IS) =
     1           MOMPAR(IS,IB,1)/US*SYHZ(IX,BCPOS(IB))                !srv 24.08.08
                SMO0(IX,BCPOS(IB),1,IS)=0.0_R8
              ENDDO
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              US=0.0_R8
              DO IX=BCSTART(IB), BCEND(IB)
                US=US+SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
              ENDDO
              DO IX=BCSTART(IB), BCEND(IB)
                SMO0(IX,BCPOS(IB),0,IS) =
     1           MOMPAR(IS,IB,1)/US*
     2           SYHZ(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))  !srv 24.08.08
                SMO0(IX,BCPOS(IB),1,IS)=0.0_R8
              ENDDO
            ELSEIF(BCCHAR(IB).EQ.'W') THEN
              US=0.0_R8
              DO IY=BCSTART(IB), BCEND(IB)
                US=US+SX(BCPOS(IB),IY)
              ENDDO
              DO IY=BCSTART(IB), BCEND(IB)
                SMO0(BCPOS(IB),IY,0,IS)=MOMPAR(IS,IB,1)/US*
     1           SXHZ(BCPOS(IB),IY)                   !srv 24.08.08
                SMO0(BCPOS(IB),IY,1,IS)=0.0_R8
              ENDDO
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              US=0.0_R8
              DO IY=BCSTART(IB), BCEND(IB)
                US=US+SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
              ENDDO
              DO IY=BCSTART(IB), BCEND(IB)
                SMO0(BCPOS(IB),IY,0,IS)=MOMPAR(IS,IB,1)/US*
     1           SXHZ(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))      !srv 24.08.08
                SMO0(BCPOS(IB),IY,1,IS)=0.0_R8
              ENDDO
            ELSE
              CALL XERRAB ('B2STBC -- BCMOM=9, NOT AVAILABLE')
            ENDIF
            GOTO 20000

* -- BCMOM=10 -- PRESCRIBE THE DECAY LENGTH FOR THE MOMENTUM

11000       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'momentum decay length',MOMPAR(IS,IB,1),' on ',
     1         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2         BCEND(IB),' for species ',IS
              call xertst(MOMPAR(IS,IB,1).ne.0.0_R8,
     1         'BCMOM = 10, MOMPAR(,,1) = 0')
!              call xertst(BCCON(IS,IB).EQ.9,
!     1         'Should use BCCON(IS,IB).EQ.9 together with BCMOM=10 ')   !srv 16.07.12
          write(*,*) 'Please be assured that this boundary condition',
     &               'is consistent with density boundary condition'      !srv 21.09.12
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
crfs check if additional terms due to tine are necessary for smo0 and SMO0V
                SMO0(IX,BCPOS(IB),0,is)=AM(IS)*mp*ua(IX,BCPOS(IB),IS)
     .           *SNA0_NO_MDF(IX,BCPOS(IB),0,IS)                          !srv 23.09.08
                SMO0(IX,BCPOS(IB),1,IS)=-vsa0(IX,BCPOS(IB),IS)/
     1           MOMPAR(IS,IB,1)*SYHZ(IX,BCPOS(IB))   +                   !srv 24.08.08
     .           AM(IS)*mp*SNA0_NO_MDF(IX,BCPOS(IB),1,IS)*                !srv 23.09.08
     .           na(IX,BCPOS(IB),IS)
              ENDDO
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SMO0(IX,BCPOS(IB),0,IS)=AM(IS)*mp*ua(IX,BCPOS(IB),IS)
     .           *SNA0_NO_MDF(IX,BCPOS(IB),0,IS)                          !srv 23.09.08
                SMO0(IX,BCPOS(IB),1,IS)=-vsa0(IX,BCPOS(IB),IS)/
     1           MOMPAR(IS,IB,1)*
     2           SYHZ(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB))) +    !srv 24.08.08
     .           AM(IS)*mp*SNA0_NO_MDF(IX,BCPOS(IB),1,IS)*                !srv 23.09.08
     .           na(IX,BCPOS(IB),IS)
              ENDDO
            ELSE
              CALL XERRAB ('B2STBC -- BCMOM=10, NOT AVAILABLE FOR W,E')
            ENDIF
            GOTO 20000

* -- BCMOM=11 -- ROZHANSKY VISCOSITY CONDITION FOR THE MOMENTUM

11100       CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'ROZ momentum bc on ',BCCHAR(IB),
     1         ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB),
     2         ' for species ',IS
            endif
#ifdef WANT_THIS
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=max(0,BCSTART(IB)), BCEND(IB)
                smo0(ix,bcpos(ib),0,is)=
     1           vsa0(ix,bcpos(ib),is)*na(ix,bcpos(ib),is)*am(is)*
     2           sqrt(1-pit(ix,bcpos(ib))**2)*ua(ix,bcpos(ib),is)
                smo0(ix,bcpos(ib),1,is)=
     1           -dna0(ix,bcpos(ib),is)*am(is)*
     2           (na(ix,bcpos(ib),is)-
     2            na(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),is))/
     3           hy1(ix,bcpos(ib))-
     3           vsa0(ix,bcpos(ib),is)*na(ix,bcpos(ib),is)*
     4           am(is)*sqrt(1-pit(ix,bcpos(ib))**2)/hy1(ix,bcpos(ib))
              end do
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCMOM=11, NOT AVAILABLE FOR N,W,E')
            ENDIF
#else
            CALL XERRAB ('B2STBC -- BCMOM=11, NOT IMPLEMENTED.')
#endif
            GOTO 20000

* -- BCMOM=12 -- CONDITION from b2stbc_spb FOR THE MOMENTUM         !srv 01.02.09 {

11200       CONTINUE
            if(ncall.eq.0) then !srv added accumulation in order to account twice b.c. for corner cells
              write(*,*) 'from b2stbc_spb momentum bc with mompar1 =',
     1         MOMPAR(IS,IB,1),' and mompar2 =', MOMPAR(IS,IB,2),' on ',
     2         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     3         BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                s0hz = hz(ix,BCPOS(IB))*vol(ix,BCPOS(IB))/
     &                                  hy1(ix,BCPOS(IB))             !srv 09.01.01
                SMO0(IX,BCPOS(IB),0,IS)= 0.0_R8
                SMO0(IX,BCPOS(IB),1,IS)= 0.0_R8
                SMO0(IX,BCPOS(IB),2,IS)= 0.0_R8
                SMO0(IX,BCPOS(IB),3,IS)= SMO0(IX,BCPOS(IB),3,IS)+
     &		   MOMPAR(IS,IB,2)*s0hz +
!			   MOMPAR(IS,IB,2)*SYHZ(IX,BCPOS(IB))+
     &           MOMPAR(IS,IB,1)*
     &           max(0.0e0_R8,-fna_fcor(IX,TOPIY(IX,BCPOS(IB)),1,is))/
     &                          na(IX,BCPOS(IB),is)
              end do
            ELSEIF(BCCHAR(IB).EQ.'N') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                s0hz = hz(ix,BCPOS(IB))*vol(ix,BCPOS(IB))/
     &                                  hy1(ix,BCPOS(IB))             !srv 09.01.01
                SMO0(IX,BCPOS(IB),0,IS)= 0.0_R8
                SMO0(IX,BCPOS(IB),1,IS)= 0.0_R8
                SMO0(IX,BCPOS(IB),2,IS)= 0.0_R8
                SMO0(IX,BCPOS(IB),3,IS)= SMO0(IX,BCPOS(IB),3,IS)+
     &           MOMPAR(IS,IB,2)*s0hz+
!     &           MOMPAR(IS,IB,2)*
!     &           SYHZ(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))+
     &           MOMPAR(IS,IB,1)*
     &           max(0.0e0_R8,fna_fcor(IX,BCPOS(IB),1,is))/
     &                          na(IX,BCPOS(IB),is)
              end do
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCMOM=12, NOT AVAILABLE FOR W,E')
            ENDIF
            GOTO 20000

* -- BCMOM=13 -- CONDITION from b2stbc_spb FOR THE MOMENTUM

11300       CONTINUE
            if(ncall.eq.0) then   !srv added accumulation in order to account twice b.c. for corner cells
              write(*,*) 'from b2stbc_spb momentum bc with mompar1 =',
     1         MOMPAR(IS,IB,1),' and mompar2 =', MOMPAR(IS,IB,2),' on ',
     2         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     3         BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'W') THEN
              DO IY=BCSTART(IB), BCEND(IB)
!                CS=SQRT(BCINT(PR,NX,NY,IY,IB)/BCINT(RO,NX,NY,IY,IB))
               cs=sqrt(pr(bcpos(ib),iy)/ro(bcpos(ib),iy))
               s1hz = pbshz(rightix(bcpos(ib),iy),
     &                           rightiy(bcpos(ib),iy),0)           !srv 11.01.13
               vwest=-cs*s1hz -                                !srv 17.02.04 19.02.08 11.01.13
     &          hz(bcpos(ib),iy)*vol(rightix(bcpos(ib),iy),iy)/
     &          hx(rightix(bcpos(ib),iy),iy)*                       !*signmf*bb(ix,iy,2)/bb(ix,iy,3)   !srv 05.02.05
     &          vaecrb(rightix(bcpos(ib),iy),iy,0,is)*cbc*          !srv 17.02.04 19.02.08
     &          (abs(bb(bcpos(ib),iy,0))/bb(bcpos(ib),iy,0))        !srv 11.01.13
               SMO0(BCPOS(IB),IY,0,IS) = SMO0(BCPOS(IB),IY,0,IS) +  !srv 07.05.09 {
     &            MOMPAR(IS,IB,1)*(UBIG/EPSLON**2)*AM(IS)*mp*
     &            na(bcpos(ib),iy,is)*vwest
               SMO0(BCPOS(IB),IY,1,IS)= SMO0(BCPOS(IB),IY,1,IS)-
     &          (UBIG/EPSLON**2)*AM(IS)*mp*na(bcpos(ib),iy,is)*abs(s1hz) !srv 07.05.09 }
!                SMO0(BCPOS(IB),IY,2,IS)=SMO0(BCPOS(IB),IY,2,IS)+      !srv 07.05.09 { equivalence form
!     &                   MOMPAR(IS,IB,1)*(UBIG/EPSLON**2)*vwest
!                SMO0(BCPOS(IB),IY,3,IS)= SMO0(BCPOS(IB),IY,3,IS)-     !srv 07.05.09 }
!     &                   (UBIG/EPSLON**2)*S1HZ
              end do
            ELSEIF(BCCHAR(IB).EQ.'E') THEN
              DO IY=BCSTART(IB), BCEND(IB)
!               CS=SQRT(BCINT(PR,NX,NY,IY,IB)/BCINT(RO,NX,NY,IY,IB))
               cs=sqrt(pr(bcpos(ib),iy)/ro(bcpos(ib),iy))
               s1hz = pbshz(bcpos(ib),iy,0)                          !srv 11.01.13
               veast = cs*s1hz-                                      !srv 17.02.04 19.02.08 25.04.08 { 11.01.13
     &          hz(bcpos(ib),iy)*vol(leftix(bcpos(ib),iy),iy)/
     &          hx(leftix(bcpos(ib),iy),iy)                          !*signmf*bb(ix,iy,2)/bb(ix,iy,3) !srv 05.02.05
     &          *vaecrb(bcpos(ib),iy,0,is)*cbc*                      !srv 17.02.04 19.02.08
     &          (abs(bb(bcpos(ib),iy,0))/bb(bcpos(ib),iy,0))         !srv 11.01.13
               SMO0(BCPOS(IB),IY,0,IS) = SMO0(BCPOS(IB),IY,0,IS) +   !srv 07.05.09 {
     &            MOMPAR(IS,IB,1)*(UBIG/EPSLON**2)*AM(IS)*mp*
     &            na(bcpos(ib),iy,is)*veast
               SMO0(BCPOS(IB),IY,1,IS)= SMO0(BCPOS(IB),IY,1,IS)-
     &          (UBIG/EPSLON**2)*AM(IS)*mp*na(bcpos(ib),iy,is)*abs(s1hz)  !srv 07.05.09 }
!                SMO0(BCPOS(IB),IY,2,IS)=SMO0(BCPOS(IB),IY,2,IS)+      !srv 07.05.09 { equivalence form
!     &               MOMPAR(IS,IB,1)*(UBIG/EPSLON**2)*veast
!                SMO0(BCPOS(IB),IY,3,IS)= SMO0(BCPOS(IB),IY,3,IS) -    !srv 07.05.09 }
!     &               (UBIG/EPSLON**2)*S1HZ
              END DO
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCMOM=13, NOT AVAILABLE FOR S,N')
            ENDIF
            GOTO 20000                                                 !srv 01.02.09 }

* -- BCMOM=14 -- CONDITION from b2stbc_spb FOR THE MOMENTUM         !srv 01.02.09 {

11400       CONTINUE
            if(ncall.eq.0) then !srv added accumulation in order to account twice b.c. for corner cells
              write(*,*) 'from b2stbc_spb momentum bc with mompar1 =',
     1         MOMPAR(IS,IB,1),' and mompar2 =', MOMPAR(IS,IB,2),' on ',
     2         BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     3         BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                s0hz = hz(ix,BCPOS(IB))*vol(ix,BCPOS(IB))/
     &                                  hy1(ix,BCPOS(IB))              !srv 09.01.01
                if (.not.is_neutral(is)) then                          !srv 15.02.02 {
	           select case (iform_corr_dia)                          !srv 25.07.05
	            case (0)
                   smo0(ix,BCPOS(IB),0,is) = smo0(ix,BCPOS(IB),0,is)+  !srv 09.01.01 01.07.05
     &               s0hz*(am(is)*mp)*naaverage(is)*                   !srv 15.02.02
     &               (-2.0e0_R8*facdrift(ix,BCPOS(IB))*phm2* ! signmf* !srv 07.06.00 17.05.99 25.04.00 15.02.02
     &                 tiaverage*bb(ix,BCPOS(IB),2)/
     &                 (qe*rza(ix,BCPOS(IB),is)*hx(ix,BCPOS(IB)))*     !srv 15.02.02
     &                (OnedBsqc(rightix(ix,BCPOS(IB)),
     &                          rightiy(ix,BCPOS(IB)))-
     &                 OnedBsqc(ix,BCPOS(IB)))                         !srv 05.03.99 15.02.02
     &               )*uaaverage(is) + corr_flux*s0hz                  !srv 15.02.02 08.09.03 25.07.05
                  case (1)
                   smo0(ix,BCPOS(IB),0,is) = smo0(ix,BCPOS(IB),0,is)+  !srv 09.01.01 01.07.05
     &               s0hz*(am(is)*mp)*naaverage(is)*                   !srv 15.02.02
     &               (-2.0e0_R8*facdrift(ix,BCPOS(IB))*phm2* ! signmf* !srv 07.06.00 17.05.99 25.04.00 15.02.02
     &                 tiaverage*bb(ix,BCPOS(IB),2)/
     &               (qe*rza(ix,BCPOS(IB),is)*hx(ix,BCPOS(IB)))*                     !srv 15.02.02
     &               (OnedBsqc(rightix(ix,BCPOS(IB)),
     &                    rightiy(ix,BCPOS(IB)))-OnedBsqc(ix,BCPOS(IB))) !srv 05.03.99 15.02.02
     &                )*ua(ix,BCPOS(IB),is)+corr_flux*s0hz               !srv 08.09.03 25.07.05
	           end select                                           !srv 25.07.05
	          endif                                                 !srv 15.02.02
                SMO0(IX,BCPOS(IB),1,IS)= 0.0_R8
                SMO0(IX,BCPOS(IB),2,IS)= 0.0_R8
                SMO0(IX,BCPOS(IB),3,IS)= SMO0(IX,BCPOS(IB),3,IS)+
     &		   MOMPAR(IS,IB,2)*s0hz +
!			   MOMPAR(IS,IB,2)*SYHZ(IX,BCPOS(IB))+
     &           MOMPAR(IS,IB,1)*
     &           max(0.0e0_R8,-fna_fcor(IX,TOPIY(IX,BCPOS(IB)),1,is))/
     &                          na(IX,BCPOS(IB),is)
              end do
            ELSE
              CALL XERRAB
     1         ('B2STBC -- BCMOM=14, NOT AVAILABLE FOR N,W,E')
            ENDIF
            GOTO 20000

* -- BCMOM=15 -- PRESCRIBE THE VALUE OF THE PARALLEL VELOCITY AS A FUNCTION OF !srv 02.04.09 {
*                POLOIDAL COORDINATE

11500	    CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'specified parallel velocity as a function of',
     &         ' poloidal coordinate ',MOMPAR(IS,IB,1),
     1         ' on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2         ' to ',BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SMO0(IX,BCPOS(IB),0,IS)= VBIG*NE00*AM(IS)*
     .           MOMPAR(IS,IB,1)*(bbaverage/bb(ix,BCPOS(IB),3))*
     &           SYHZ(IX,BCPOS(IB))                                   !srv 24.08.08
                SMO0(IX,BCPOS(IB),1,IS)=-VBIG*NE00*AM(IS)*
     .           SYHZ(IX,BCPOS(IB))                                   !srv 24.08.08
              end do
           ELSE
              CALL XERRAB ('B2STBC -- BCMOM=15, ONLY FOR S')
            ENDIF
            GOTO 20000                                                !srv 02.04.09 }

* -- BCMOM=16 -- PRESCRIBE THE AVERAGE VALUE OF THE PARALLEL VELOCITY !lk 12.12.12 { 

11600	    CONTINUE
            if(ncall.eq.0) then
              write(*,*) 'specified parallel velocity',MOMPAR(IS,IB,1),
     1         ' on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2         ' to ',BCEND(IB),' for species ',IS
            endif
            IF(BCCHAR(IB).EQ.'S') THEN
              DO IX=BCSTART(IB), BCEND(IB)
                SMO0(IX,BCPOS(IB),0,IS)= VBIG*NE00*AM(IS)*
     .           MOMPAR(IS,IB,1)*SYHZ(IX,BCPOS(IB))                   !srv 24.08.08
                SMO0(IX,BCPOS(IB),1,IS)=-VBIG*NE00*AM(IS)*
     .           SYHZ(IX,BCPOS(IB))                                   !srv 24.08.08
              end do
		    t0=0.0e0_R8                                                 
	        t1=0.0e0_R8
              do ix = -1,nx                                               
                if (mod(region(ix,iy,0),4).eq.1) then                      
	           t0=t0 + vol(ix,0)                        
	           t1=t1 + UA(IX,TOPIY(IX,BCPOS(IB)),IS)*vol(ix,0)                        
                endif                                                      
	        enddo
		    t1=t1/t0                                                       
              do ix = -1,nx                                               
		     if(mod(region(ix,bcpos(ib),0),4).eq.1) then
                 SMO0(IX,BCPOS(IB),0,IS)= SMO0(IX,BCPOS(IB),0,IS)+
     .		   VBIG*NE00*AM(IS)*
     .           (UA(IX,TOPIY(IX,BCPOS(IB)),IS)-t1)*SYHZ(IX,BCPOS(IB))                   
	         endif
	        enddo
	                                                     !lk 12.12.12 }
            ELSE
              CALL XERRAB ('B2STBC -- BCMOM=16, IS VALID FOR CORE ONLY')
            ENDIF
            GOTO 20000

* -- END LOOP OVER SPECIES
20000     continue
        enddo

* -- BOUNDARY CONDITIONS FOR ELECTRON ENERGY EQUATION --

        if(BCENE(IB).eq.0) goto 30000
	GOTO(20100,20200,20300,20400,20500,
     1	 20600,20700,20800,20900,21000,
     2	 21100,21200,21300,21400,21500,
     3	 21600,21700,21800,21900,22000,22100,22200), BCENE(IB)   !srv 23.07.10
        write(*,*) BCENE(IB)
        call xerrab ('BCENE out of bounds')

* -- BCENE=1 -- PRESCRIBE THE VALUE OF THE ELECTRON TEMPERATURE

20100   CONTINUE
        if(ncall.eq.0) then
#ifdef ASTRA
          IF(BCCHAR(IB).EQ.'S') THEN
	      ENEPAR(IB,1) = teb_astra*1.0e3_R8                        !srv 14.06.12
          ENDIF
#endif
          write(*,*) 'specified electron temperature',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2     ' to ',BCEND(IB)
          call xertst(ENEPAR(IB,1).gt.0.0_R8,
     1     'BCENE = 1, ENEPAR(,1) !> 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)= VBIG*NE00*ENEPAR(IB,1)*EV
     .       *SY(IX,BCPOS(IB))
            she0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)= VBIG*NE00*ENEPAR(IB,1)*EV*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            she0(IX,BCPOS(IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            she0(BCPOS(IB),IY,0)= UBIG*NE00*ENEPAR(IB,1)*EV
     .       *SX(BCPOS(IB),IY)
            she0(BCPOS(IB),IY,1)=-UBIG*NE00*SX(BCPOS(IB),IY)
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            she0(BCPOS(IB),IY,0)= UBIG*NE00*ENEPAR(IB,1)*EV*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            she0(BCPOS(IB),IY,1)=-UBIG*NE00*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          end do
        ELSEIF(BCCHAR(IB).EQ.'X') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)= UBIG*NE00*
     .       ENEPAR(IB,1)*EV*VOL(IX,BCPOS(IB))
            she0(IX,BCPOS(IB),1)=-UBIG*NE00*
     .       VOL(IX,BCPOS(IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'Y') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            she0(BCPOS(IB),IY,0)= UBIG*NE00*
     .       ENEPAR(IB,1)*EV*VOL(BCPOS(IB),IY)
            she0(BCPOS(IB),IY,1)=-UBIG*NE00*
     .       VOL(BCPOS(IB),IY)
          end do
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCENE=1, INVALID OPTION.')
        ENDIF

        GOTO 30000

* -- BCENE=2 -- PRESCRIBE THE GRADIENT OF THE ELECTRON TEMPERATURE

20200   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'electron temperature gradient',ENEPAR(IB,1),
     1     ' on',BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2     ' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)=VBIG*NE00*EV*SY(IX,BCPOS(IB))*
     .       (-ENEPAR(IB,1)*0.5_R8*
     .       (HY1(IX,BCPOS(IB))+
     .        HY1(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)))) +
     .       TE(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB))))
            she0(IX,BCPOS(IB),1)=- VBIG*NE00*EV*SY(IX,BCPOS(IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)= VBIG*NE00*EV*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))*
     .       (ENEPAR(IB,1)*0.5_R8*
     .       (HY1(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))+
     .        HY1(IX,BCPOS(IB))) +
     .       TE(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB))))
            she0(IX,BCPOS(IB),1)=-VBIG*NE00*EV*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            she0(BCPOS(IB),IY,0)= VBIG*NE00*EV*SX(BCPOS(IB),IY)*
     .       (-ENEPAR(IB,1)*0.5_R8*
     .       (HX(BCPOS(IB),IY)+
     .        HX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))) +
     .       TE(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)))
            she0(BCPOS(IB),IY,1)=-VBIG*NE00*EV*SX(BCPOS(IB),IY)
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            she0(BCPOS(IB),IY,0)= VBIG*NE00*EV*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))*
     .       (ENEPAR(IB,1)*0.5_R8*
     .       (HX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))+
     .        HX(BCPOS(IB),IY)) +
     .       TE(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)))
            she0(BCPOS(IB),IY,1)=-VBIG*NE00*EV*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          end do
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCENE=2, INVALID OPTION.')
        END IF
        GOTO 30000

* -- BCENE=3 -- SHEATH CONDITIONS, ELECTRON ENERGY TRANSMISSION

20300   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'sheath electron with MOMPAR(ismain,IB,2) =',
     1     MOMPAR(ismain,IB,2),' and MOMPAR(ismain,IB,1) =',
     2     MOMPAR(ismain,IB,1),' and POTPAR(IB,2) =',POTPAR(IB,2),
	3     ' and ENEPAR(IB,1) =',ENEPAR(IB,1),
     3     ' on ',BCCHAR(IB),' at ',BCPOS(IB),
     3     ' from ',BCSTART(IB),' to ',BCEND(IB)
          if(fix_fch_in_fhe_sheath.eq.0 .or. fix_sx_in_fhe_sheath.eq.0) !srv 08.05.09
     1     write(*,*) 'WARNING: Wrong (old) implementation of ',
     2     'the FCH contribution to FHE used'
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
               if(sound_model.eq.0) then
                 CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))  !srv 01.02.09
               else
                  zarg = 0.0_R8
                  do is = 0, ns-1
                     zarg = zarg + bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .                    (bcint(rza(-1,-1,is),nx,ny,IY,IB)*
     .                     bcint(te,nx,ny,IY,IB)+bcint(ti,nx,ny,IY,IB))
     .                    /(am(is)*bcint(ni(-1,-1,0),nx,ny,IY,IB)*mp)
                  enddo
                  CS=SQRT(zarg)
               endif
              UUOUT=abs(wbbr(bcpos(ib),iy,0)/wbbr(bcpos(ib),iy,3))*
     .         MOMPAR(ismain,IB,1)*CS
              if(fix_fch_in_fhe_sheath.eq.0) then
                NUEOUT=MIN(0.0_R8,-UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(BCPOS(IB),IY,0)/QE)
              else
			 if(fix_sx_in_fhe_sheath.eq.0) then                        !srv 08.05.09 {
                NUEOUT=MIN(0.0_R8,-UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY),0)/QE)
               else
                NUEOUT=MIN(0.0_R8,
     &		    -UUOUT*bcint(NE,nx,ny,IY,IB)*SX(BCPOS(IB),IY) -
     &           fch(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY),0)/QE)
               endif                                                     !srv 08.05.09 }
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
                  UUTEST=
     .             UU(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY),IS)*
     .             (1.0_R8+
     .              HX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))/
     .      HX(RIGHTIX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),
     .         RIGHTIY(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))))-
     .      UU(RIGHTIX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),
     .         RIGHTIY(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),IS)*
     .              HX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))/
     .      HX(RIGHTIX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),
     .         RIGHTIY(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)))
                  UUOUT=MAX( MOMPAR(IS,IB,1)*CS*
     .             abs(wbbr(bcpos(ib),iy,0)/wbbr(bcpos(ib),iy,3)),
     .             -UUTEST)
                  NUEOUT=NUEOUT - UUOUT*bcint(rza(-1,-1,IS),nx,ny,IY,IB)
     .             *bcint(na(-1,-1,IS),nx,ny,IY,IB)
                endif
              end do
              if(fix_fch_in_fhe_sheath.eq.0) then
                NUEOUT=MIN(0.0_R8, NUEOUT - fch(BCPOS(IB),IY,0)/QE )
              else
			 if(fix_sx_in_fhe_sheath.eq.0) then                        !srv 08.05.09 {
                NUEOUT=MIN(0.0_R8, NUEOUT - 
     &           fch(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY),0)/QE )
	         else
                NUEOUT=MIN(0.0_R8, NUEOUT*SX(BCPOS(IB),IY) - 
     &           fch(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY),0)/QE )
               endif                                                     !srv 08.05.09 }
	        endif
            ENDIF
            if(enepar(ib,1).gt.0.0_R8) then                             !xpb
	       if(fix_sx_in_fhe_sheath.eq.0) then                        !srv 08.05.09 {
              she0(BCPOS(IB),IY,1)=she0(BCPOS(IB),IY,1) +
     .         (ENEPAR(IB,1)+
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT*SX(BCPOS(IB),IY)
	       else
              she0(BCPOS(IB),IY,1)=she0(BCPOS(IB),IY,1) +
     .         (ENEPAR(IB,1)+
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT
	       endif                                                      !srv 08.05.09 }
            else           ! xpb new pressure linearisation as per B2.5-style
	       if(fix_sx_in_fhe_sheath.eq.0) then                         !srv 08.05.09 {
              she0(BCPOS(IB),IY,3)= she0(bcpos(ib),iy,3) +              !xpb
     .         (-ENEPAR(IB,1)+                                          !xpb
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT*SX(BCPOS(IB),IY)/
     .          ne(bcpos(ib),iy)
	       else
              she0(BCPOS(IB),IY,3)= she0(bcpos(ib),iy,3) +              !xpb
     .         (-ENEPAR(IB,1)+                                          !xpb
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT/ne(bcpos(ib),iy)
	       endif                                                      !srv 08.05.09 }
            endif                                                       !xpb
	      if(mdf_fhe.ne.0) then                                !srv 12.09.08 {
             she0(BCPOS(IB),IY,0) = she0(BCPOS(IB),IY,0) -
     &                   fhePSch(rightix(BCPOS(IB),iy),iy,0)
	      endif                                                !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
               if(sound_model.eq.0) then
                 CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))  !srv 01.02.09
               else
                  zarg = 0.0_R8
                  do is = 0, ns-1
                     zarg = zarg + bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .                    (bcint(rza(-1,-1,is),nx,ny,IY,IB)*
     .                     bcint(te,nx,ny,IY,IB)+bcint(ti,nx,ny,IY,IB))
     .                    /(am(is)*
     .                    bcint(ni(-1,-1,0),nx,ny,IY,IB)*mp)
                  enddo
                  CS=SQRT(zarg)
               endif
              UUOUT=MOMPAR(ismain,IB,1)*CS*
     .         abs(wbbl(bcpos(ib),iy,0)/wbbl(bcpos(ib),iy,3))
              if(fix_fch_in_fhe_sheath.eq.0) then
                NUEOUT=MAX(0.0_R8,UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY),0)/QE)
              else
               if(fix_sx_in_fhe_sheath.eq.0) then                    !srv 08.05.09 {
                NUEOUT=MAX(0.0_R8,UUOUT*bcint(NE,nx,ny,IY,IB) -
     .           fch(BCPOS(IB),IY,0)/QE)
	         else
                NUEOUT=MAX(0.0_R8,UUOUT*bcint(NE,nx,ny,IY,IB)*
     &		         SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)) -
     &           fch(BCPOS(IB),IY,0)/QE)
	         endif                                                 !srv 08.05.09 }
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
                  UUTEST=
     .             UU(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY),IS)*
     .             (1.0_R8+
     .              HX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))/
     .       HX(LEFTIX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),
     .          LEFTIY(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))))-
     .       UU(LEFTIX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),
     .          LEFTIY(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),IS)*
     .              HX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))/
     .       HX(LEFTIX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),
     .          LEFTIY(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)))
                  UUOUT=MAX( MOMPAR(IS,IB,1)*CS*
     .             abs(wbbl(bcpos(ib),iy,0)/wbbl(bcpos(ib),iy,3)),
     .             +UUTEST)
                  NUEOUT=NUEOUT + UUOUT*bcint(rza(-1,-1,IS),nx,ny,IY,IB)
     .             *bcint(na(-1,-1,IS),nx,ny,IY,IB)
                endif
              end do
              if(fix_fch_in_fhe_sheath.eq.0) then
                NUEOUT=MAX(0.0_R8, NUEOUT -
     .           fch(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY),0)/QE )
              else
               if(fix_sx_in_fhe_sheath.eq.0) then                    !srv 08.05.09 {
                NUEOUT=MAX(0.0_R8, NUEOUT - fch(BCPOS(IB),IY,0)/QE )
	         else
                NUEOUT=MAX(0.0_R8, NUEOUT*
     &		   SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)) -
     &		   fch(BCPOS(IB),IY,0)/QE )
	         endif                                                 !srv 08.05.09 }
              endif
            ENDIF
            if(enepar(ib,1).gt.0.0_R8) then                             !xpb
             if(fix_sx_in_fhe_sheath.eq.0) then                      !srv 08.05.09 {
              she0(BCPOS(IB),IY,1)=she0(BCPOS(IB),IY,1) -
     .         (ENEPAR(IB,1) +
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT*
     .          SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
	       else
              she0(BCPOS(IB),IY,1)=she0(BCPOS(IB),IY,1) -
     .         (ENEPAR(IB,1) +
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT
     	       endif                                                   !srv 08.05.09 }
            else           ! xpb new pressure linearisation as per B2.5-style
             if(fix_sx_in_fhe_sheath.eq.0) then                      !srv 08.05.09 {
              she0(BCPOS(IB),IY,3)=she0(BCPOS(IB),IY,3) -              !xpb
     .         (-ENEPAR(IB,1) +                                        !xpb
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT*
     .          SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
     .         /ne(bcpos(ib),iy)
	       else
              she0(BCPOS(IB),IY,3)=she0(BCPOS(IB),IY,3) -              !xpb
     .         (-ENEPAR(IB,1) +                                        !xpb
     .          QE*(bcint(PO,nx,ny,IY,IB)-POTPAR(IB,2))/
     .          bcint(TE,nx,ny,IY,IB))*NUEOUT/ne(bcpos(ib),iy)
	       endif                                                   !srv 08.05.09 }
            endif                                                      !xpb
	      if(mdf_fhe.ne.0) then                                !srv 12.09.08 {
             she0(BCPOS(IB),IY,0) = she0(BCPOS(IB),IY,0) +
     &                           fhePSch(BCPOS(IB),iy,0)
	      endif                                                !srv 12.09.08 }
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=3, NO SHEATH IN RADIAL DIR.')
        END IF
        GOTO 30000

* -- BCENE=4 -- PRESCRIBE THE VALUE OF THE ELECTRON TEMPERATURE,
*                 WEAKLY A MIXED BOUNDARY CONDITION

20400   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'weak electron temperature',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',
     2     BCSTART(IB),' to ',BCEND(IB)
          call xertst(ENEPAR(IB,1).gt.0.0_R8,
     1     'BCENE = 4, ENEPAR(,1) !> 0')
          call xertst(ENEPAR(IB,2).gt.0.0_R8,
     1     'BCENE = 4, ENEPAR(,2) !> 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)= ENEPAR(IB,2)*ENEPAR(IB,1)*
     .       SY(IX,BCPOS(IB))
            she0(IX,BCPOS(IB),1)=-ENEPAR(IB,2)/EV*SY(IX,BCPOS(IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)= ENEPAR(IB,2)*ENEPAR(IB,1)*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            she0(IX,BCPOS(IB),1)=-ENEPAR(IB,2)/EV*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            she0(BCPOS(IB),IY,0)= ENEPAR(IB,2)*ENEPAR(IB,1)*
     1       SX(BCPOS(IB),IY)
            she0(BCPOS(IB),IY,1)=-ENEPAR(IB,2)/EV*SX(BCPOS(IB),IY)
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
           DO IY=BCSTART(IB), BCEND(IB)
            she0(BCPOS(IB),IY,0)= ENEPAR(IB,2)*ENEPAR(IB,1)*
     1       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            she0(BCPOS(IB),IY,1)=-ENEPAR(IB,2)/EV*
     2       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=4, NOT YET AVAILABLE')
        END IF
        GOTO 30000

* -- BCENE=5 -- PRESCRIBE THE ELECTRON ENERGY FLUX PER UNIT AREA

20500   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'electron energy flux per area',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2     ' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)=ENEPAR(IB,1)*SY(IX,BCPOS(IB))
            she0(IX,BCPOS(IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)=ENEPAR(IB,1)*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            she0(IX,BCPOS(IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            she0(BCPOS(IB),IY,0)=ENEPAR(IB,1)*SX(BCPOS(IB),IY)
            she0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            she0(BCPOS(IB),IY,0)=ENEPAR(IB,1)*
     1       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            she0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=5, NOT YET AVAILABLE')
        END IF
        GOTO 30000

* -- BCENE=6 -- PRESCRIBE THE TOTAL ELECTRON ENERGY FLUX FOR A
*                 CONSTANT ELECTRON TEMPERATURE

cghgh
20600   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'total electron energy flux',ENEPAR(IB,1),
     1     ' with const Te on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',
     2     BCSTART(IB),' to ',BCEND(IB)
        endif
#ifdef WANT_THIS
        IF(BCCHAR(IB).EQ.'N') THEN
          Tguess=0.0_R8
          DO IX=BCSTART(IB), BCEND(IB)
            Tguess=Tguess+
     +       TE(bottomix(IX,BCPOS(IB)),bottomiy(ix,bcpos(ib)))
          ENDDO
          Tguess=Tguess/(BCEND(IB)-BCSTART(IB)+1)
          if(diagno.ge.2) write(*,*) 'Initial Tguess:',Tguess/ev
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
          DO IX=BCSTART(IB), BCEND(IB)
cxpb code segment analog to b2tfhe
            if(ncall.eq.0) then
             vedia(ix,bcpos(ib),1) = 
     .        ti(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)))*
     &        bb(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)),2)       !xpb
     &        /(qe*hx(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib))))   !xpb
     &        *(2.0_R8/
     &          bb(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)),3)**3) !xpb
     &        *(wbbr(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)),3)-
     &          wbbl(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)),3))
     &        *facdrift(ix,bcpos(ib))                                   !xpb
             veecrb(ix,bcpos(ib),1) = 1.0_R8/hx(ix,bcpos(ib))*          !xpb
     1          bb(ix,bcpos(ib),2)/bb(ix,bcpos(ib),3)**2*               !xpb
     2         (po(rightix(ix,bcpos(ib)),rightiy(ix,bcpos(ib)))+
     3          po(leftix(ix,bcpos(ib)),leftiy(ix,bcpos(ib))))*0.5_R8*
     &         fac_ExB(ix,bcpos(ib))                                    !xpb
            endif
!xpb            if(no_current.eq.1) then
!xpb              ehy = 0.0_R8
!xpb            else
            ehy = -(po(ix,bcpos(ib))-
     &     po(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib))))+
     &    (ne(ix,bcpos(ib))*Tguess-
     &     ne(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)))*
     &     te(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib))))*
     &    (2.0_R8/qe)/(ne(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)))
     &                +ne(ix,bcpos(ib)))
!xpb            endif
            if (conduction_only.eq.1) then
              flo = 0.0_R8
            else
              prefactor=
     &         (vol(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)))+
     &          vol(ix,bcpos(ib)))
     &         /(hy1(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)))+
     &           hy1(ix,bcpos(ib)))
     &         *(ne(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)))
     &          +ne(ix,bcpos(ib)))/2.0_R8
              flo_vedia=prefactor*vedia(ix,bcpos(ib),1)
              flo_veecrb=prefactor*veecrb(ix,bcpos(ib),1)
              flo = (1.5_R8+BoRiS)*fne_32(ix,bcpos(ib),1)+ 
     1         2.5_R8*fne_52(ix,bcpos(ib),1)+
     &       chve(ix,bcpos(ib),1)*(ne(ix,bcpos(ib))+
     &       ne(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib))))/2.0_R8+
     &         2.5_R8*flo_vedia+1.5_R8*flo_veecrb
            endif
            if (no_hybr_e.eq.0) then
              cone0 = hybr(flo,chce(ix,bcpos(ib),1))
              fy=flo*(te(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)))+
     &         Tguess)/2.0_R8-cone0*(Tguess-
     &         te(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib))))+
     &         calf(ix,bcpos(ib),1)*ehy*
     &         harm(te(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib))),
     &         Tguess)
            else
              cone0 = chce(ix,bcpos(ib),1)
              t2 = dim ( cone0,
     .         abs(flo) *
     .         max(0.0_R8, 
     .          min(Tguess,
     .              te(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)))))/
     .         max(cutlo, Tguess+
     .          te(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)))) )
              fy = 
     .         upwind(flo, 
     .                te(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib))), 
     .                Tguess)
     .         - t2 * (Tguess-te(bottomix(ix,bcpos(ib)),
     .                           bottomiy(ix,bcpos(ib)))) +
     &         calf(ix,bcpos(ib),1)*harm(te(bottomix(ix,bcpos(ib)),
     &                                      bottomiy(ix,bcpos(ib))),
     &                            Tguess)*ehy
            endif
            fyt=fyt + fy
          ENDDO
          if(diagno.ge.2)
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
            if(diagno.ge.2) write(*,*) 'GAM0, GAM1 = ',GAM0,GAM1
            if(diagno.ge.2) write(*,*) 'fyt,t1,t0 = ',fyt,t1/ev,t0/ev
            if (ibnd .gt. 100)
     >       call xerrab ('B2STBC -- BCENE=6, NOT CONVERGING')
            if (abs(t1-t0).gt.Tlimit .and.
     &          abs(fyt-ENEPAR(IB,1)).gt.fylimt) goto 20610
          endif
          if(diagno.ge.1)
     >     write(*,'(a,f15.6,i4,1pe14.6)') ' tebound:',T1/ev,ibnd,fyt
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)= VBIG*NE00*T1*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            she0(IX,BCPOS(IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'S') THEN
          Tguess=0.0_R8
          DO IX=BCSTART(IB), BCEND(IB)
            Tguess=Tguess+TE(IX,BCPOS(IB))
          ENDDO
          Tguess=Tguess/(BCEND(IB)-BCSTART(IB)+1)
          if(diagno.ge.2) write(*,*) 'Initial Tguess:',Tguess/ev
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
          DO IX=BCSTART(IB), BCEND(IB)
cxpb code segment analog to b2tfhe
            if(ncall.eq.0) then
             vedia(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1) =          !xpb
     &        ti(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))*               !xpb
     &        bb(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),2)/             !xpb
     &        (qe*hx(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))))           !xpb
     &        *(2.0_R8/bb(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),3)**3) !xpb
     &        *(wbbr(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),3)-         !xpb
     &          wbbl(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),3))         !xpb
     &        *facdrift(ix,bcpos(ib))                                    !xpb
             veecrb(ix,bcpos(ib),1) = 1.0_R8/hx(ix,bcpos(ib))*           !xpb
     1          bb(ix,bcpos(ib),2)/bb(ix,bcpos(ib),3)**2*                !xpb
     2         (po(rightix(ix,bcpos(ib)),rightiy(ix,bcpos(ib)))+
     3          po(leftix(ix,bcpos(ib)),leftiy(ix,bcpos(ib))))*0.5_R8*
     &         fac_ExB(ix,bcpos(ib))                                     !xpb
            endif
!xpb            if(no_current.eq.1) then
!xpb               ehy = 0.0_R8
!xpb            else
            ehy = -(po(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))-
     &                po(ix,bcpos(ib)))+
     &         (ne(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))*
     &          te(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))-
     &          ne(ix,bcpos(ib))*Tguess)*
     &         (2.0_R8/qe)/(ne(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))
     &                     +ne(ix,bcpos(ib)))
!xpb            endif
            if (conduction_only.eq.1) then
              flo = 0.0_R8
            else
              prefactor=
     &         (vol(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))+
     &          vol(ix,bcpos(ib)))
     &         /(hy1(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))+
     &           hy1(ix,bcpos(ib)))
     &         *(ne(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))                    !xpb
     &          +ne(ix,bcpos(ib)))/2.0_R8                                       !xpb
              flo_vedia=prefactor*vedia(ix,bcpos(ib),1)
              flo_veecrb=prefactor*veecrb(ix,bcpos(ib),1)
              flo = (1.5_R8+BoRiS)*
     .         fne_32(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1)+ 
     1         2.5_R8*fne_52(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1)+
     &         chve(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1)*
     &        (ne(ix,bcpos(ib))+
     &         ne(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))))/2.0_R8+
     &         2.5_R8*flo_vedia+1.5_R8*flo_veecrb
            endif
            if (no_hybr_e.eq.0) then
              cone0 = hybr(flo,
     &         chce(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1))
              fy=flo*(te(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))+
     &         Tguess)/2.0_R8-cone0*
     &         (te(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))-Tguess)+
     &         calf(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1)*ehy*
     &         harm(te(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))),Tguess)
            else
              cone0 = chce(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1)
              t2 = dim ( cone0,
     .         abs(flo) *
     .         max(0.0_R8, 
     .          min(Tguess,
     .              te(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))))/
     .         max(cutlo, Tguess+
     .          te(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))) )
              fy = 
     .         upwind(flo, Tguess,
     .                te(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))))
     .         - t2 * (te(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))-
     .                 Tguess) +
     &         calf(ix,bcpos(ib),1)*
     &         harm(te(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))),
     &              Tguess)*ehy
            endif
            fyt=fyt + fy
          ENDDO
          if(diagno.ge.2)
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
          if(diagno.ge.1)
     >     write(*,'(a,f15.6,i4,1pe14.6)') ' tebound:',T1/ev,ibnd,fyt
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)= VBIG*NE00*T1*SY(IX,BCPOS(IB))
            she0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=6, NOT AVAILABLE FOR W,E')
        ENDIF
#else
        CALL XERRAB ('B2STBC -- BCENE=6, NOT IMPLEMENTED.')
#endif
        GOTO 30000
cghgh

* -- BCENE=7 -- PRESCRIBE THE ELECTRON TEMPERATURE AS A FUNCTION OF
*                 OTHER PLASMA PARAMETERS

20700   CONTINUE
        CALL XERRAB ('B2STBC -- BCENE=7, NOT YET AVAILABLE')
        GOTO 30000

* -- BCENE=8 -- PRESCRIBE THE TOTAL ELECTRON HEAT FLUX WITH CONSTANT
*                 FLUX DENSITY

20800   CONTINUE
        if(ncall.eq.0) then
#ifdef ASTRA
       ENEPAR(IB,1) = fheb_astra                                !srv 28.07.11
#endif
          write(*,*) 'total electron energy flux',ENEPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2     BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          US=0.0_R8
          DO IX=BCSTART(IB), BCEND(IB)
            US=US+SY(IX,BCPOS(IB))
          end do
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)=ENEPAR(IB,1)/US*SY(IX,BCPOS(IB))
            she0(IX,BCPOS(IB),1)=0.0_R8
	      if(mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(ix,bcpos(ib),0) = she0(ix,bcpos(ib),0) -
     &                              fhePSch(ix,topiy(ix,bcpos(ib)),1)
              she0(ix,bcpos(ib),1) = she0(ix,bcpos(ib),1) +
     &            1.5_R8*vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &           veecrb(ix,topiy(ix,bcpos(ib)),1)*ne(ix,bcpos(ib))
	      endif                                                !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          US=0.0_R8
          DO IX=BCSTART(IB), BCEND(IB)
            US=US+SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          end do
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)=ENEPAR(IB,1)/US*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            she0(IX,BCPOS(IB),1)=0.0_R8
	      if(mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(ix,bcpos(ib),0) = she0(ix,bcpos(ib),0) +
     &                            fhePSch(ix,topiy(ix,bcpos(ib)),1)
              she0(ix,bcpos(ib),1) = she0(ix,bcpos(ib),1) -
     &          1.5_R8*vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &          veecrb(ix,topiy(ix,bcpos(ib)),1)*ne(ix,bcpos(ib))
	      endif                                                !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          US=0.0_R8
          DO IY=BCSTART(IB), BCEND(IB)
            US=US+SX(BCPOS(IB),IY)
          end do
          DO IY=BCSTART(IB), BCEND(IB)
            she0(BCPOS(IB),IY,0)=ENEPAR(IB,1)/US*SX(BCPOS(IB),IY)
            she0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          US=0.0_R8
          DO IY=BCSTART(IB), BCEND(IB)
            US=US+SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          end do
          DO IY=BCSTART(IB), BCEND(IB)
            she0(BCPOS(IB),IY,0)=ENEPAR(IB,1)/US*
     1       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            she0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=8, NOT AVAILABLE')
        ENDIF
        GOTO 30000

* -- BCENE=9 -- PRESCRIBE THE DECAY LENGTH FOR THE ELECTRON TEMPERATURE

20900   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'electron temperature decay length',ENEPAR(IB,1),
     1     ' on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2     ' to ',BCEND(IB)
          call xertst(ENEPAR(IB,1).ne.0.0_R8,
     1     'BCENE = 9, ENEPAR(,1) = 0')
!          call xertst(all(BCCON(0:NS-1,IB).EQ.9),
!     1     'Should use BCCON(IS,IB).EQ.9 together with BCENE=9 ')       !srv 16.07.12
          write(*,*) 'Please be assured that this boundary condition',
     &               'is consistent with density boundary condition'     !srv 21.09.12
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)=0.0_R8
            she0(IX,BCPOS(IB),1)=-hce0(ix,BCPOS(IB))/ENEPAR(IB,1)*
     1       SY(IX,BCPOS(IB))
            DO IS=0,NS-1
              she0(IX,BCPOS(IB),0)=she0(IX,BCPOS(IB),0) +
     .         rza(ix,bcpos(ib),is)*2.5_R8*TE(IX,BCPOS(IB))*
     .         SNA0_NO_MDF(IX,BCPOS(IB),0,IS)                          !srv 23.09.08
              she0(IX,BCPOS(IB),1)=she0(IX,BCPOS(IB),1) +
     .         rza(ix,bcpos(ib),is)*2.5_R8*
     .         SNA0_NO_MDF(IX,BCPOS(IB),1,IS)*na(IX,BCPOS(IB),IS)      !srv 23.09.08
            ENDDO
	      if(mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(ix,bcpos(ib),0) = she0(ix,bcpos(ib),0) -
     &                              fhePSch(ix,topiy(ix,bcpos(ib)),1)
              she0(ix,bcpos(ib),1) = she0(ix,bcpos(ib),1) +
     &            1.5_R8*vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &           veecrb(ix,topiy(ix,bcpos(ib)),1)*ne(ix,bcpos(ib))
	      endif                                                !srv 12.09.08 }
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)=0.0_R8
            she0(IX,BCPOS(IB),1)=-hce0(ix,BCPOS(IB))/ENEPAR(IB,1)*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            DO IS=0,NS-1
              she0(IX,BCPOS(IB),0)=she0(IX,BCPOS(IB),0) +
     .         rza(ix,bcpos(ib),is)*2.5_R8*TE(IX,BCPOS(IB))*
     .         SNA0_NO_MDF(IX,BCPOS(IB),0,IS)                           !srv 23.09.08
              she0(IX,BCPOS(IB),1)=she0(IX,BCPOS(IB),1) +
     .         rza(ix,bcpos(ib),is)*2.5_R8*
     .         SNA0_NO_MDF(IX,BCPOS(IB),1,IS)*na(IX,BCPOS(IB),IS)       !srv 23.09.08
            ENDDO
	      if(mdf_fhe.ne.0) then                                !srv 12.09.08 {
              she0(ix,bcpos(ib),0) = she0(ix,bcpos(ib),0) +
     &                              fhePSch(ix,bcpos(ib),1)
              she0(ix,bcpos(ib),1) = she0(ix,bcpos(ib),1) -
     &           1.5_R8*vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &           veecrb(ix,bcpos(ib),1)*ne(ix,bcpos(ib))
	      endif                                                !srv 12.09.08 }
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=9, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 30000

* -- BCENE=10 -- Feedback option for core

21000   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'electron heat flux feedback bc on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),
     2     ' from ',BCSTART(IB),' to ',BCEND(IB)
          call xertst(lfeedback, 'LFEEDBACK must be turned on!')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            if (region(ix,bcpos(ib),0).eq.1) then
              she0(IX,BCPOS(IB),0)=cbshe(0,coreregno)*SY(IX,BCPOS(IB))
            elseif (region(ix,bcpos(ib),0).eq.5) then
              she0(IX,BCPOS(IB),0)=cbshe(0,coreregn2)*SY(IX,BCPOS(IB))
            endif
            she0(IX,BCPOS(IB),1)=0.0_R8
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=10, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 30000

* -- BCENE=11 -- FROM B2STBC_SPB                                       !srv 01.02.09 {

21100   CONTINUE
        if(ncall.eq.0) then   !srv added accumulation in order to account twice b.c. for corner cells
          write(*,*) ' FROM B2STBC_SPB with enepar1 =',ENEPAR(IB,1),
     1     ' and enepar2 =',ENEPAR(IB,2),
     2     ' on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     3     ' to ',BCEND(IB)
!          call xertst(ENEPAR(IB,1).ne.0.0_R8,
!     1     'BCENE = 11, ENEPAR(,1) = 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
!            vte = sqrt(bcint(te,nx,ny,ix,ib)/me)
            vte = sqrt(te(IX,BCPOS(IB))/me)
            s0 = vol(IX,BCPOS(IB))/hy1(IX,BCPOS(IB))
            if(mdf_fhe.ne.0) then                                    !srv 21.03.09 {
             she0(IX,BCPOS(IB),0) = she0(IX,BCPOS(IB),0) - 
     &          fhePSch(ix,topiy(ix,bcpos(ib)),1)
	      endif                                                    !srv 21.03.09 }
            she0(IX,BCPOS(IB),1)= she0(IX,BCPOS(IB),1)+ENEPAR(IB,1)*
     &       max(0.0e0_R8,-floe_noc(ix,topiy(ix,bcpos(ib)),1))       !srv 28.10.02
            she0(IX,BCPOS(IB),3)= she0(IX,BCPOS(IB),3)+ENEPAR(IB,2)*
     &        vte*s0
!     &       vte*SY(IX,BCPOS(IB))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
!            vte = sqrt(bcint(te,nx,ny,ix,ib)/me)
            vte = sqrt(te(IX,BCPOS(IB))/me)
            s0 = vol(IX,BCPOS(IB))/hy1(IX,BCPOS(IB))
            if(mdf_fhe.ne.0) then                                    !srv 21.03.09 {
             she0(IX,BCPOS(IB),0) = she0(IX,BCPOS(IB),0) + 
     &          fhePSch(ix,bcpos(ib),1)
	      endif                                                    !srv 21.03.09 }
            she0(IX,BCPOS(IB),1)=she0(IX,BCPOS(IB),1)+ENEPAR(IB,1)*
     &       max(0.0e0_R8,floe_noc(ix,bcpos(ib),1))                  !srv 28.10.02
            she0(IX,BCPOS(IB),3)= she0(IX,BCPOS(IB),3)+ENEPAR(IB,2)*
     &        vte*s0
!     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=11, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 30000                                                    !srv 01.02.09 }

* -- BCENE=12 -- SHEATH CONDITIONS, EL ENERGY TRANSMISSION COEFFICIENT

21200   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'sheath (12) el on ',BCCHAR(IB),
     1     ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(ix,bcpos(ib),3) = she0(ix,bcpos(ib),3)+
     1       ENEPAR(IB,1)*max(0.0_R8,
     1        -fne(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1))/
     2       ne(ix,bcpos(ib))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(ix,bcpos(ib),3) = she0(ix,bcpos(ib),3)+
     1             ENEPAR(IB,1)*max(0.0_R8,fne(ix,bcpos(ib),1))/
     2       ne(ix,bcpos(ib))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            she0(bcpos(ib),iy,3) = she0(bcpos(ib),iy,3)+
     1       ENEPAR(IB,1)*max(0.0_R8,
     1        -fne(rightix(bcpos(ib),iy),rightiy(bcpos(ib),iy),0))/
     2       ne(bcpos(ib),iy)
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            she0(bcpos(ib),iy,3) = she0(bcpos(ib),iy,3)+
     1             ENEPAR(IB,1)*max(0.0_R8,fne(bcpos(ib),iy,0))/
     2       ne(bcpos(ib),iy)
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=12, INVALID OPTION')
        END IF
        GOTO 30000

* -- BCENE=13 -- PRESCRIBE THE ELECTRON ENERGY FLUX PER UNIT AREA
*             -- PROPORTIONAL TO TEMPERATURE

21300   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'electron energy flux per area proportional to ',
     1     'temperature ',ENEPAR(IB,1),' on ',BCCHAR(IB),' at ',
     2     BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)=ENEPAR(IB,1)*TE(IX,BCPOS(IB))*
     1       SY(IX,BCPOS(IB))
            she0(IX,BCPOS(IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)=ENEPAR(IB,1)*TE(IX,BCPOS(IB))*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            she0(IX,BCPOS(IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            she0(BCPOS(IB),IY,0)=ENEPAR(IB,1)*TE(BCPOS(IB),IY)*
     1       SX(BCPOS(IB),IY)
            she0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            she0(BCPOS(IB),IY,0)=ENEPAR(IB,1)*TE(BCPOS(IB),IY)*
     1       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            she0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=13, NOT YET AVAILABLE')
        END IF
        GOTO 30000

* -- BCENE=14 -- leakage option for electron energy

21400   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'electron energy leakage',ENEPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2     BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            SHE0(IX,BCPOS(IB),0)=0.0_R8
            if(sound_model.eq.0) then
              CS=SQRT(bcint(PR,nx,ny,IX,IB)/bcint(RO,nx,ny,IX,IB))   !srv 01.02.09
            else
               zarg = 0.0_R8
               do is = 0, ns-1
                  zarg = zarg + bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .                 (bcint(rza(-1,-1,is),nx,ny,IY,IB)*
     .                  bcint(te,nx,ny,IY,IB)+bcint(ti,nx,ny,IY,IB))
     .                 /(am(is)*bcint(ni(-1,-1,0),nx,ny,IY,IB)*mp)
               enddo
               CS=SQRT(zarg)
            endif
            SHE0(IX,BCPOS(IB),1)=
     1       vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*cs*ENEPAR(IB,1)
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            SHE0(IX,BCPOS(IB),0)=0.0_R8
            if(sound_model.eq.0) then
              CS=SQRT(bcint(PR,nx,ny,IX,IB)/bcint(RO,nx,ny,IX,IB))   !srv 01.02.09
            else
               zarg = 0.0_R8
               do is = 0, ns-1
                  zarg = zarg + bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .                 (bcint(rza(-1,-1,is),nx,ny,IY,IB)*
     .                  bcint(te,nx,ny,IY,IB)+bcint(ti,nx,ny,IY,IB))
     .                 /(am(is)*bcint(ni(-1,-1,0),nx,ny,IY,IB)*mp)
               enddo
               CS=SQRT(zarg)
            endif
            SHE0(IX,BCPOS(IB),1)=
     1       vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*cs*ENEPAR(IB,1)
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=14, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 30000

* -- BCENE=15 -- FROM B2STBC_SPB                                     !srv 01.02.09 {

21500   CONTINUE
        if(bcene_15_style.eq.0) then
         if(ncall.eq.0) then   !srv added accumulation in order to account twice b.c. for corner cells
           write(*,*) 'bcene_15_style = 0 with GAMMAE =',GAMMAE,
     1      ' on ',BCCHAR(IB),
     2      ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
         endif
         IF(BCCHAR(IB).EQ.'W') THEN
           DO IY=BCSTART(IB), BCEND(IB)
!            CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))
            CS=SQRT(PR(bcpos(ib),iy)/RO(bcpos(ib),iy))
            t0 = max(0.0e0_R8,                                        !srv 05.05.09
     &       (abs(pbs(rightix(bcpos(ib),iy),iy,0))*ne(bcpos(ib),iy)*cs +   !srv 25.12.06 11.01.13
     &        fch_p(rightix(bcpos(ib),iy),iy,0)/qe))                  !srv 25.12.06 11.03.09 }
	      if(pot_eq.eq.1) then                                      !srv 27.02.13 {
             she0(bcpos(ib),iy,1) = she0(bcpos(ib),iy,1)-             !srv 11.03.09 {
     &	        (1.0e0_R8+GAMMAE)/(1.0e0_R8-GAMMAE)*t0            !srv 19.04.13
             she0(bcpos(ib),iy,0) = she0(bcpos(ib),iy,0) - 
     &                               qe*po(bcpos(ib),iy)*t0           !srv 07.11.12
	      else                                                      !srv 27.02.13
             she0(bcpos(ib),iy,1) = she0(bcpos(ib),iy,1)-             !srv 11.03.09 {
     &	                           (1.0e0_R8+3.1_R8)*t0           !srv 19.04.13
             she0(bcpos(ib),iy,0) = she0(bcpos(ib),iy,0) 
            endif                                                     !srv 27.02.13 }
            if(mdf_fhe.ne.0) then                                     !srv 01.06.06 {
              she0(bcpos(ib),iy,0) = she0(bcpos(ib),iy,0) - 
     &           fhePSch(rightix(bcpos(ib),iy),iy,0)
	      endif                                                     !srv 01.06.06 }
           END DO
         ELSEIF(BCCHAR(IB).EQ.'E') THEN
           DO IY=BCSTART(IB), BCEND(IB)
!            CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))
             CS=SQRT(PR(bcpos(ib),iy)/RO(bcpos(ib),iy))
             t0 = max(0.0e0_R8,                                        !srv 05.05.09
     &        (abs(pbs(bcpos(ib),iy,0))*ne(bcpos(ib),iy)*cs -          !srv 25.12.06 11.01.13
     &         fch_p(bcpos(ib),iy,0)/qe))                              !srv 25.12.06 11.03.09 }
	       if(pot_eq.eq.1) then                                      !srv 27.02.13 {
              she0(bcpos(ib),iy,1) = she0(bcpos(ib),iy,1)-               !srv 11.03.09 {
     &	        (1.0e0_R8+GAMMAE)/(1.0e0_R8-GAMMAE)*t0            !srv 19.04.13
              she0(bcpos(ib),iy,0) = she0(bcpos(ib),iy,0) - 
     &                              qe*po(bcpos(ib),iy)*t0             !srv 07.11.12
	       else                                                      !srv 27.02.13
              she0(bcpos(ib),iy,1) = she0(bcpos(ib),iy,1)-             !srv 11.03.09 {
     &	                          (1.0e0_R8+3.1_R8)*t0             !srv 19.04.13
              she0(bcpos(ib),iy,0) = she0(bcpos(ib),iy,0)
             endif                                                     !srv 27.02.13 }
             if(mdf_fhe.ne.0) then                                     !srv 01.06.06 {
              she0(bcpos(ib),iy,0) = she0(bcpos(ib),iy,0) + 
     &           fhePSch(bcpos(ib),iy,0)
	       endif                                                     !srv 01.06.06 }
           end do
         ELSE
           CALL XERRAB ('B2STBC -- BCENE=15, not aviable for S and N')
         END IF

        elseif(bcene_15_style.eq.1) then
         if(ncall.eq.0) then   !srv added accumulation in order to account twice b.c. for corner cells
           write(*,*) 'bcene_15_style = 1 with GAMMAE =',GAMMAE,
     1      ' on ',BCCHAR(IB),
     2      ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
         endif
         IF(BCCHAR(IB).EQ.'W') THEN
           DO IY=BCSTART(IB), BCEND(IB)
!             CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))
             s1=abs(pbs(rightix(bcpos(ib),iy),rightiy(bcpos(ib),iy),0)) !srv 18.04.13 {
             vte=sqrt(te(bcpos(ib),iy)/me)
             t0 =s1*ne(bcpos(ib),iy)*sqrt(1.0_R8/(2.0_R8*pi))*vte*
     &           expu2(max(-50.0_R8,
     &          -qe*(po(bcpos(ib),iy)-potpar(ib,2))/te(bcpos(ib),iy)))
	       if(pot_eq.eq.1) then                                      !srv 27.02.13 {
              she0(bcpos(ib),iy,1) = she0(bcpos(ib),iy,1)-
     &                   (1.0e0_R8+GAMMAE)*t0
              she0(bcpos(ib),iy,0) = she0(bcpos(ib),iy,0) - 
     &                   (1.0e0_R8-GAMMAE)*qe*po(bcpos(ib),iy)*t0
	       else                                                      !srv 27.02.13
              cs=sqrt(pr(bcpos(ib),iy)/ro(bcpos(ib),iy))
              she0(bcpos(ib),iy,1) = she0(bcpos(ib),iy,1)-
     &	                s1*(1.0e0_R8+3.1_R8)*ne(bcpos(ib),iy)*cs
              she0(bcpos(ib),iy,0) = she0(bcpos(ib),iy,0)              !srv 18.04.13 }
             endif                                                     !srv 27.02.13 }
             if(mdf_fhe.ne.0) then                                     !srv 01.06.06 {
              she0(bcpos(ib),iy,0) = she0(bcpos(ib),iy,0) - 
     &           fhePSch(rightix(bcpos(ib),iy),iy,0)
	       endif                                                     !srv 01.06.06 }
           END DO
         ELSEIF(BCCHAR(IB).EQ.'E') THEN
           DO IY=BCSTART(IB), BCEND(IB)
!             CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))
             s1=abs(pbs(bcpos(ib),iy,0))                               !srv 18.04.13 { 
             vte=sqrt(te(bcpos(ib),iy)/me)
             t0 =s1*ne(bcpos(ib),iy)*sqrt(1.0_R8/(2.0_R8*pi))*vte*
     &           expu2(max(-50.0_R8,
     &          -qe*(po(bcpos(ib),iy)-potpar(ib,2))/te(bcpos(ib),iy)))
	       if(pot_eq.eq.1) then                                      !srv 27.02.13 {
              she0(bcpos(ib),iy,1) = she0(bcpos(ib),iy,1)-
     &	            (1.0e0_R8+GAMMAE)*t0
              she0(bcpos(ib),iy,0) = she0(bcpos(ib),iy,0) - 
     &                   (1.0e0_R8-GAMMAE)*qe*po(bcpos(ib),iy)*t0 
	       else                                                      !srv 27.02.13
              cs=sqrt(pr(bcpos(ib),iy)/ro(bcpos(ib),iy))
              she0(bcpos(ib),iy,1) = she0(bcpos(ib),iy,1)-
     &	       s1*(1.0e0_R8+3.1_R8)*ne(bcpos(ib),iy)*cs
              she0(bcpos(ib),iy,0) = she0(bcpos(ib),iy,0)              !srv 18.04.13 } 
             endif                                                     !srv 27.02.13 }
             if(mdf_fhe.ne.0) then                                     !srv 01.06.06 {
              she0(bcpos(ib),iy,0) = she0(bcpos(ib),iy,0) + 
     &           fhePSch(bcpos(ib),iy,0)
	       endif                                                     !srv 01.06.06 }
           end do
         ELSE
           CALL XERRAB ('B2STBC -- BCENE=15, not aviable for S and N')
         END IF
        else
           CALL XERRAB ('B2STBC -- BCENE=15, incoorrect bcene_15_style')
	  endif
        GOTO 30000                                                   !srv 01.02.09 }

* -- BCENE=16 -- Feedback boundary condition with constant temperature

21600   CONTINUE

        if(bcene_16_style.eq.0) then
	   
         if(ncall.eq.0) then
           write(*,*) 'electron temperature feedback by flux',
     1      ENEPAR(IB,1),' with ',ENEPAR(IB,2),' on ', 
     1      BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2      ' to ',BCEND(IB)
           call xertst(enepar(ib,1).ne.0, 'faulty parameter enepar1')
           call xertst(enepar(ib,2).gt.0, 'faulty parameter enepar2')
         endif
         IF(BCCHAR(IB).EQ.'S') THEN
           totflux = 0.0_R8
           t0 = 0.0_R8
           weight = 0.0_R8
           DO IX=BCSTART(IB), BCEND(IB)
             bc_type16_ix=ix
             bc_type16_iy=BCPOS(IB)
             do bc_type16_i=1,bc_type16_ref      ! how many steps away from the boundary?
               bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
               bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
               bc_type16_ix=bc_type16_ixx
               call xertst (-2.lt.bc_type16_ix.and.bc_type16_ix.lt.nx+1,
     1          'stepped off x-grid in type 16 bc')
               call xertst (-2.lt.bc_type16_iy.and.bc_type16_iy.lt.ny+1,
     1          'stepped off y-grid in type 16 bc')
             enddo
             totflux=totflux+
     .        fhe(bc_type16_ix,bc_type16_iy,1)
             t0 = t0 +
     .        te(ix,bcpos(ib))*vol(ix,bcpos(ib))*ne(ix,bcpos(IB))
             weight = weight + vol(ix,bcpos(ib))*ne(ix,bcpos(IB))
           ENDDO
           t0 = t0 / weight
           if(diagno.ge.2) write(*,*) 'Totflux :',totflux
           T1 = (1.0_R8 - 
     .      enepar(ib,2)*(totflux-enepar(ib,1))/enepar(ib,1)) * t0 
           if(diagno.ge.2)
     1      write(*,*) 'Old/New electron temp :',t0/qe, t1/qe
           DO IX=BCSTART(IB), BCEND(IB)
             she0(IX,BCPOS(IB),0)= VBIG*NE00*T1*SY(IX,BCPOS(IB))
             she0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
           ENDDO
         ELSE
           CALL XERRAB ('B2STBC -- BCENE=16, NOT AVAILABLE FOR N,W,E')
         ENDIF
         GOTO 30000

        elseif(bcene_16_style.eq.1)then

          if(ncall.eq.0) then
		 write(*,*) 'Given integrated electron heat flux',
     &                  ' with constant Te ' 
!          write all boundaries with this type of boundary condition
           DO IBW=IB,NBC
		  if(BCCHAR(IB).EQ.'S'.and.BCENE(IBW).EQ.BCENE(IB))THEN	       
		   write(*,*) ' on ',BCCHAR(IBW),' at ',BCPOS(IBW),' from ',
     2          BCSTART(IBW),' to ', BCEND(IBW)
              call xertst(enepar(ibw,1).ne.0,'faulty parameter enepar1')
	      endif
	     ENDDO
#ifdef ASTRA
           SENEPAR = fheb_astra
#else
           SENEPAR = 0.0_R8
           DO IBW=IB,NBC
		  if(BCCHAR(IB).EQ.'S'.and.BCENE(IBW).EQ.BCENE(IB))THEN	       
               SENEPAR = SENEPAR + ENEPAR(IBW,1)
	      endif
	     ENDDO
#endif
		 write(*,*) ' Summary value is ', SENEPAR

	     TES = TE((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB))
          endif
	    if(.not.once_bcene_16) then
	     IF(BCCHAR(IB).EQ.'S') THEN
	      totflux = 0.0_R8
            DO IBW=IB,NBC
             IF(BCCHAR(IBW).EQ.'S'.AND.BCENE(IBW).EQ.BCENE(IB)) THEN
	        if(mdf_fhe.ne.0) then
               DO IX=BCSTART(IBW), BCEND(IBW)
	          totflux = totflux + 
     &          fhe_mdf(topix(ix,BCPOS(IBW)),topiy(ix,BCPOS(IBW)),1)
     &       +fac_ExB(ix,BCPOS(IBW))*                                    !srv 11.04.13 {
     &      ((po(ix,BCPOS(IBW))-po(leftix(ix,BCPOS(IBW)),BCPOS(IBW)))/
     &       (hx(ix,BCPOS(IBW))+hx(leftix(ix,BCPOS(IBW)),BCPOS(IBW)))+
     &       (po(rightix(ix,BCPOS(IBW)),BCPOS(IBW))-po(ix,BCPOS(IBW)))/
     &       (hx(rightix(ix,BCPOS(IBW)),BCPOS(IBW))+hx(ix,BCPOS(IBW))))*
     &        ne(ix,BCPOS(IBW))*te(ix,BCPOS(IBW))*SY(IX,BCPOS(IBW))*
     &        bb(ix,BCPOS(IBW),2)/bb(ix,BCPOS(IBW),3)**2                 !srv 11.04.13 }
               ENDDO
	        else
               DO IX=BCSTART(IBW), BCEND(IBW)
	          totflux = totflux + 
     &          fhe(topix(ix,BCPOS(IBW)),topiy(ix,BCPOS(IBW)),1) 
               ENDDO
	        endif
	       ENDIF
	      ENDDO
            TES = TE((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB)) * 
     &         (1.0_R8-bc_ref_te*(totflux-SENEPAR)/                      !srv 11.04.13
     &           max(abs(totflux),SENEPAR))
             write(*,"(' BCENE=16: TES=',es14.7,' totflux=',es14.7,
     *       ' SENEPAR=',es14.7)") TES/ev, totflux, SENEPAR
            DO IBW=IB,NBC
             IF(BCCHAR(IBW).EQ.'S'.AND.BCENE(IBW).EQ.BCENE(IB)) THEN
              DO IX=BCSTART(IBW), BCEND(IBW)
               she0(IX,BCPOS(IBW),0)= VBIG*NE00*TES*SY(IX,BCPOS(IBW))
               she0(IX,BCPOS(IBW),1)=-VBIG*NE00*SY(IX,BCPOS(IBW))
              ENDDO
	       ENDIF
	      ENDDO
		  once_bcene_16 = .true.
           ELSE
             CALL XERRAB
     1         ('B2STBC -- BCENE=16, NOT AVAILABLE FOR N,W,E')
           ENDIF
	    endif                                                          !srv 16.02.12 }
          GOTO 30000

	  else
	   write(*,*) ' Incorrect style of boundary condition'
	   stop
	  endif

* -- BCENE=17 -- Feedback boundary condition with constant shared temperature

21700   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'shared temperature feedback by flux',
     1     ENEPAR(IB,1),' with ',ENEPAR(IB,2),' on ', 
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2     ' to ',BCEND(IB)
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
          DO IX=BCSTART(IB), BCEND(IB)
            bc_type16_ix=ix
            bc_type16_iy=BCPOS(IB)
            do bc_type16_i=1,bc_type16_ref      ! how many steps away from the boundary?
              bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
              bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
              bc_type16_ix=bc_type16_ixx
              call xertst (-2.lt.bc_type16_ix.and.bc_type16_ix.lt.nx+1,
     1         'stepped off x-grid in type 16 bc')
              call xertst (-2.lt.bc_type16_iy.and.bc_type16_iy.lt.ny+1,
     1         'stepped off y-grid in type 16 bc')
            enddo
            totflux=totflux+fhe(bc_type16_ix,bc_type16_iy,1)+
     .                      fhi(bc_type16_ix,bc_type16_iy,1)
            if (bc_type16_kinetic_energy.eq.1) then
              do is = 0, ns-1
                totflux=totflux+(1.0_R8-BoRiS)*
     .           fna(bc_type16_ix,bc_type16_iy,1,is)*0.5_R8*
     .           (kinrgy(bc_type16_ix,bc_type16_iy,is)+
     .            kinrgy(bottomix(bc_type16_ix,bc_type16_iy),
     .                   bottomiy(bc_type16_ix,bc_type16_iy),is))
              enddo
            endif
            t0 = t0 +
     .       te(ix,bcpos(ib))*vol(ix,bcpos(ib))*ne(ix,bcpos(IB)) +
     .       ti(ix,bcpos(ib))*vol(ix,bcpos(ib))*ni(ix,bcpos(IB),0)
            weight = weight +
     1       vol(ix,bcpos(ib))*ne(ix,bcpos(IB)) +
     2       vol(ix,bcpos(ib))*ni(ix,bcpos(IB),0)
          ENDDO
          t0 = t0 / weight
          if(diagno.ge.2) write(*,*) 'Totflux :',totflux
          T1 = (1.0_R8 - 
     .     enepar(ib,2)*
     1       (totflux*enepar(ib,1)/(enepar(ib,1)+enipar(ib,1))-
     2         enepar(ib,1))/enepar(ib,1) -
     .     enipar(ib,2)*
     1       (totflux*enipar(ib,1)/(enepar(ib,1)+enipar(ib,1))-
     2         enipar(ib,1))/enipar(ib,1)
     .    ) * t0 
          if(diagno.ge.2)
     1     write(*,*) 'Old/New shared temp :',t0/qe, t1/qe
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)= VBIG*NE00*T1*SY(IX,BCPOS(IB))
            she0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
            shi0(IX,BCPOS(IB),0)= VBIG*NE00*T1*SY(IX,BCPOS(IB))
            shi0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=17, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 30000

* -- BCENE=18 -- FRACTIONAL DROP CONDITION (FOR EDGE2D COMPARISON)

21800   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'electron temperature fractional drop',
     1     ENEPAR(IB,1),' on ',BCCHAR(IB),' at ',BCPOS(IB),
     2     ' from ',BCSTART(IB),' to ',BCEND(IB)
          call xertst(ENEPAR(IB,1).lt.1.0_R8,
     1     'BCENE = 18, ENEPAR(,1) => 1')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            T1 = (1.0_R8 - ENEPAR(IB,1)) *
     .       TE(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)))
            she0(IX,BCPOS(IB),0)= VBIG*NE00*T1*SY(IX,BCPOS(IB))
            she0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            T1 = (1.0_R8 - ENEPAR(IB,1)) *
     .       TE(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            she0(IX,BCPOS(IB),0)= VBIG*NE00*T1*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            she0(IX,BCPOS(IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            T1 = (1.0_R8 - ENEPAR(IB,1)) *
     .       TE(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))
            she0(BCPOS(IB),IY,0)= VBIG*NE00*T1*SX(BCPOS(IB),IY)
            she0(BCPOS(IB),IY,1)=-VBIG*NE00*SX(BCPOS(IB),IY)
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            T1 = (1.0_R8 - ENEPAR(IB,1)) *
     .       TE(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            she0(BCPOS(IB),IY,0)= VBIG*NE00*T1*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            she0(BCPOS(IB),IY,1)=-VBIG*NE00*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=18, INVALID OPTION')
        ENDIF
        GOTO 30000

* -- BCENE=19 -- PRESCRIBE THE DECAY LENGTH FOR THE ELECTRON TEMPERATURE

21900   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'electron temperature decay length (new)',
     .     ENEPAR(IB,1), ' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2     BCEND(IB)
          call xertst(ENEPAR(IB,1).ne.0.0_R8,
     1     'BCENE = 19, ENEPAR(,1) = 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            she0(IX,BCPOS(IB),0)=0.0_R8
            she0(IX,BCPOS(IB),1)=-hce0(IX,BCPOS(IB))/ENEPAR(IB,1)*
     1       SY(IX,BCPOS(IB))
            t0=0.0_R8
            DO IS=0,NS-1
               t0=t0-rza(ix,bcpos(ib),is)*(
     .              1.5_R8*fna_32(topix(ix,bcpos(ib)),
     .                            topiy(ix,bcpos(ib)),1,is)+
     .              2.5_R8*fna_52(topix(ix,bcpos(ib)),
     .                            topiy(ix,bcpos(ib)),1,is))
            ENDDO
            if(t0.gt.0.0_R8) then
               she0(IX,BCPOS(IB),0)=she0(IX,BCPOS(IB),0) + 
     .              t0 * TE(IX,BCPOS(IB))
            else
               she0(IX,BCPOS(IB),1)=she0(IX,BCPOS(IB),1) + t0
            endif
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
             she0(IX,BCPOS(IB),0)=0.0_R8
             she0(IX,BCPOS(IB),1)=-hce0(IX,BCPOS(IB))/ENEPAR(IB,1)*
     1            SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
             t0=0.0_R8
             DO IS=0,NS-1
                t0=t0+rza(ix,bcpos(ib),is)*(
     .               1.5_R8*fna_32(ix,bcpos(ib),1,is)+
     .               2.5_R8*fna_52(ix,bcpos(ib),1,is))
             ENDDO
            if(t0.gt.0.0_R8) then
               she0(IX,BCPOS(IB),0)=she0(IX,BCPOS(IB),0) + 
     .              t0 * TE(IX,BCPOS(IB))
             else
                she0(IX,BCPOS(IB),1)=she0(IX,BCPOS(IB),1) + t0
             endif
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=19, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 30000

* -- BCENE=20 -- Feedback boundary condition with constant temperature

22000   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'electron temperature feedback by flux',
     1     ENEPAR(IB,1),' with ',ENEPAR(IB,2),' on ', 
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2     ' to ',BCEND(IB)
          call xertst(enepar(ib,1).ne.0, 'faulty parameter enepar1')
          call xertst(enepar(ib,2).gt.0, 'faulty parameter enepar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          totflux = 0.0_R8
          t0 = 0.0_R8
          weight = 0.0_R8
	  teave = 0.0_R8
	  neave = 0.0_R8
          DO IX=BCSTART(IB), BCEND(IB)
            bc_type16_ix=ix
            bc_type16_iy=BCPOS(IB)
            do bc_type16_i=1,bc_type16_ref      ! how many steps away from the boundary?
              bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
              bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
              bc_type16_ix=bc_type16_ixx
              call xertst (-2.lt.bc_type16_ix.and.bc_type16_ix.lt.nx+1,
     1         'stepped off x-grid in type 16 bc')
              call xertst (-2.lt.bc_type16_iy.and.bc_type16_iy.lt.ny+1,
     1         'stepped off y-grid in type 16 bc')
            enddo
            totflux=totflux+
     .       fhe(bc_type16_ix,bc_type16_iy,1)
            t0 = t0 +
     .       te(ix,bcpos(ib))*vol(ix,bcpos(ib))*ne(ix,bcpos(IB))
            weight = weight + vol(ix,bcpos(ib))*ne(ix,bcpos(IB))
            bc_type20_ix=ix
            bc_type20_iy=BCPOS(IB)
            do bc_type20_i=1,bc_type20_ref      ! how many steps away from the boundary?
              bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
              bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
              bc_type20_ix=bc_type20_ixx
              call xertst (-2.lt.bc_type20_ix.and.bc_type20_ix.lt.nx+1,
     1         'stepped off x-grid in type 20 bc')
              call xertst (-2.lt.bc_type20_iy.and.bc_type20_iy.lt.ny+1,
     1         'stepped off y-grid in type 20 bc')
            enddo
	    teave=teave+te(bc_type20_ix,bc_type20_iy)*
     1	     ne(bc_type20_ix,bc_type20_iy)*
     &          vol(bc_type20_ix,bc_type20_iy)
	    neave=neave+ne(bc_type20_ix,bc_type20_iy)*
     1	     vol(bc_type20_ix,bc_type20_iy)
          ENDDO
          t0 = t0 / weight
	  teave=teave/neave
          if(diagno.ge.2) write(*,*) 'Totflux :',totflux
          T1 = (1.0_R8 - 
     .     enepar(ib,2)*(totflux-enepar(ib,1))/enepar(ib,1)) * t0 
          if(diagno.ge.2)
     1     write(*,*) 'Old/New electron temp :',t0/qe, t1/qe
          DO IX=BCSTART(IB), BCEND(IB)
            bc_type20_ix=ix
            bc_type20_iy=BCPOS(IB)
            do bc_type20_i=1,bc_type20_ref      ! how many steps away from the boundary?
	      bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
	      bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
	      bc_type20_ix=bc_type20_ixx
	      call xertst (-2.lt.bc_type20_ix.and.bc_type20_ix.lt.nx+1,
     1	       'stepped off x-grid in type 20 bc')
	      call xertst (-2.lt.bc_type20_iy.and.bc_type20_iy.lt.ny+1,
     1	       'stepped off y-grid in type 20 bc')
	    enddo
	    she0(IX,BCPOS(IB),0)= VBIG*NE00*SY(IX,BCPOS(IB))*T1*
     1	     te(bc_type20_ix,bc_type20_iy)/teave
            she0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=20, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 30000

* -- BCENE=21 -- 

22100   CONTINUE
        CALL XERRAB ('B2STBC -- BCENE=21, NOT YET AVAILABLE')
        GOTO 30000

* -- BCENE=22 -- THE MODIFIED DECAY LENGTH FOR THE ELECTRON TEMPERATURE   !srv 23.07.10 {

22200   CONTINUE
        if(ncall.eq.0) then
          write(*,*) ' BCENE=22 with ENEPAR(IB,1)=',ENEPAR(IB,1),
     1     ' on',BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2     ' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            vte = sqrt(te(IX,BCPOS(IB))/me)
            she0(IX,BCPOS(IB),0)=VBIG*NE00*EV*
     &           CHCE(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)),1)*
     &             TE(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)))

            she0(IX,BCPOS(IB),1)=- VBIG*NE00*EV*
     &           (CHCE(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)),1)-
     &            ENEPAR(IB,1)*SY(IX,BCPOS(IB))*VTE*NE(IX,BCPOS(IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            vte = sqrt(te(IX,BCPOS(IB))/me)
            she0(IX,BCPOS(IB),0)=VBIG*NE00*EV*
     &           CHCE(IX,BCPOS(IB),1)*
     &             TE(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))

            she0(IX,BCPOS(IB),1)=- VBIG*NE00*EV*
     &           (CHCE(IX,BCPOS(IB),1)-ENEPAR(IB,1)*
     &            SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))*
     &            VTE*NE(IX,BCPOS(IB)))
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENE=22, NOT AVAILABLE FOR W,E')
        END IF                                                          !srv 23.07.10 }
        GOTO 30000

* -- BOUNDARY CONDITIONS FOR ION ENERGY EQUATIONS --

30000   continue
        if(BCENI(IB).eq.0) goto 40000
	GOTO(30100,30200,30300,30400,30500,30600,30700,30800,30900,31000,
     1     31100,31200,31300,31400,31500,31600,31700,31800,31900,32000,
     2     32100,32200,32300,32400), BCENI(IB)                         !srv 01.02.09 19.04.13
        write(*,*) BCENI(IB)
        call xerrab ('BCENI out of bounds')

* -- BCENI=1 -- PRESCRIBE THE VALUE OF THE ION TEMPERATURE

30100   CONTINUE
        if(ncall.eq.0) then
#ifdef ASTRA
          IF(BCCHAR(IB).EQ.'S') THEN
	      ENIPAR(IB,1) = tib_astra*1.0e3_R8                           !srv 14.06.12
          ENDIF
#endif
          write(*,*) 'specified ion temperature',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2     BCEND(IB)
          call xertst(ENIPAR(IB,1).gt.0.0_R8,
     1     'BCENI = 1, ENIPAR(,1) !> 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)= VBIG*NE00*ENIPAR(IB,1)*EV*
     .       SY(IX,BCPOS(IB))
            shi0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)= VBIG*NE00*ENIPAR(IB,1)*EV*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            shi0(IX,BCPOS(IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            shi0(BCPOS(IB),IY,0)= UBIG*NE00*ENIPAR(IB,1)*EV*
     .       SX(BCPOS(IB),IY)
            shi0(BCPOS(IB),IY,1)=-UBIG*NE00*SX(BCPOS(IB),IY)
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            shi0(BCPOS(IB),IY,0)= UBIG*NE00*ENIPAR(IB,1)*EV*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            shi0(BCPOS(IB),IY,1)=-UBIG*NE00*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          end do
        ELSEIF(BCCHAR(IB).EQ.'X') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)= UBIG*NE00*
     .       ENIPAR(IB,1)*EV*VOL(IX,BCPOS(IB))
            shi0(IX,BCPOS(IB),1)=-UBIG*NE00*
     .       VOL(IX,BCPOS(IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'Y') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            shi0(BCPOS(IB),IY,0)= UBIG*NE00*
     .       ENIPAR(IB,1)*EV*VOL(BCPOS(IB),IY)
            shi0(BCPOS(IB),IY,1)=-UBIG*NE00*
     .       VOL(BCPOS(IB),IY)
          end do
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCENI=1, INVALID OPTION.')
        END IF
        GOTO 40000

* -- BCENI=2 -- PRESCRIBE THE GRADIENT OF THE ION TEMPERATURE

30200   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'ion temperature gradient',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2           BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)= VBIG*NE00*EV*SY(IX,BCPOS(IB))*
     .       (-ENIPAR(IB,1)*0.5_R8*
     .         (HY1(IX,BCPOS(IB))+
     .          HY1(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)))) +
     .       TI(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB))))
            shi0(IX,BCPOS(IB),1)=-VBIG*NE00*EV*SY(IX,BCPOS(IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)= VBIG*NE00*EV*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))*
     .       (ENIPAR(IB,1)*0.5_R8*
     .         (HY1(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))+
     .          HY1(IX,BCPOS(IB))) +
     .       TI(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB))))
            shi0(IX,BCPOS(IB),1)=-VBIG*NE00*EV*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            shi0(BCPOS(IB),IY,0)= UBIG*NE00*EV*SX(BCPOS(IB),IY)*
     .       (-ENIPAR(IB,1)*0.5_R8*
     .         (HX(BCPOS(IB),IY)+
     .          HX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))) +
     .       TI(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)))
            shi0(BCPOS(IB),IY,1)=-UBIG*NE00*EV*SX(BCPOS(IB),IY)
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            shi0(BCPOS(IB),IY,0)= UBIG*NE00*EV*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))*
     .       (ENIPAR(IB,1)*0.5_R8*
     .         (HX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))+
     .          HX(BCPOS(IB),IY)) +
     .       TI(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)))
            shi0(BCPOS(IB),IY,1)=-UBIG*NE00*EV*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          end do
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCENI=2, INVALID OPTION.')
        END IF
        GOTO 40000

* -- BCENI=3 -- SHEATH CONDITIONS, ION ENERGY TRANSMISSION COEFFICIENT

30300   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'sheath ion with MOMPAR(ismain,IB,2) =',
     1     MOMPAR(ismain,IB,2),' and MOMPAR(ismain,IB,1) =',
     2     MOMPAR(ismain,IB,1)
          write(*,*) ' and enipar(ib,1) =',enipar(ib,1),               !srv 24.02.10 {
     1     ' and enipar(ib,2) =',enipar(ib,2),
     2     ' on ',BCCHAR(IB),' at ',BCPOS(IB),
     3     ' from ',BCSTART(IB),' to ',BCEND(IB)                       !srv 24.02.10 }
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
              if(sound_model.eq.0) then
                CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))  !srv 01.02.09
              else
                zarg = 0.0_R8
                do is = 0, ns-1
                   zarg = zarg + bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .                  (bcint(rza(-1,-1,is),nx,ny,IY,IB)*
     .                   bcint(te,nx,ny,IY,IB)+bcint(ti,nx,ny,IY,IB))
     .                  /(am(is)*
     .                  bcint(ni(-1,-1,0),nx,ny,IY,IB)*mp)
                enddo
                CS=SQRT(zarg)
              endif
              UUOUT=MOMPAR(ismain,IB,1)*CS*
     .         abs(wbbr(bcpos(ib),IY,0)/wbbr(bcpos(ib),iy,3))
              UPOUT=max(0.0_R8,UUOUT/
     .         abs(wbbr(bcpos(ib),IY,0)/wbbr(bcpos(ib),iy,3))+
     .         abs(wbbr(bcpos(ib),iy,2)/wbbr(bcpos(ib),iy,0))*
     .         ud(bcpos(ib),iy,1))
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  shi0(BCPOS(IB),IY,0)=shi0(BCPOS(IB),IY,0) -
     .             ENIPAR(IB,2)*0.5_R8*UUOUT*AM(IS)*mp*
     .             bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .             (UPOUT**2+uadia(bcpos(ib),iy,0,is)**2)*               !xpb
     .             SX(BCPOS(IB),IY)
                  if(enipar(ib,1).gt.0.0_R8) then                        !xpb
                    shi0(BCPOS(IB),IY,1)=shi0(BCPOS(IB),IY,1) - UUOUT *
     .               ENIPAR(IB,1)*bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .               SX(BCPOS(IB),IY)
                  else     ! xpb new pressure linearisation as per B2.5-style
                    shi0(BCPOS(IB),IY,3)=SHI0(BCPOS(IB),IY,3) - UUOUT *  !xpb
     .               (-ENIPAR(IB,1))*bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .               SX(BCPOS(IB),IY)/ni(bcpos(ib),iy,0)                 !xpb
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
                  UUTEST=
     .             UU(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY),IS)*
     .             (1.0_R8+
     .              HX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))/
     .      HX(RIGHTIX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),
     .         RIGHTIY(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))))-
     .      UU(RIGHTIX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),
     .         RIGHTIY(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),IS)*
     .              HX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))/
     .      HX(RIGHTIX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),
     .         RIGHTIY(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)))
                  UUOUT=MAX(MOMPAR(IS,IB,1)*CS*
     .             abs(wbbr(bcpos(ib),iy,0)/wbbr(bcpos(ib),iy,3)),
     .             -UUTEST)
                  UPOUT=max(0.0_R8,UUOUT/
     .             abs(wbbr(bcpos(ib),iy,0)/wbbr(bcpos(ib),iy,3))+
     .             abs(wbbr(bcpos(ib),iy,2)/wbbr(bcpos(ib),iy,0))*
     .             ud(bcpos(ib),iy,is))
                  shi0(BCPOS(IB),IY,0)=shi0(BCPOS(IB),IY,0) -
     .             ENIPAR(IB,2)*0.5_R8*UUOUT*AM(IS)*mp*
     .             bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .             (UPOUT**2+UD(BCPOS(IB),IY,IS)**2)*SX(BCPOS(IB),IY)
                  if(enipar(ib,1).gt.0.0_R8) then                        !xpb
                    shi0(BCPOS(IB),IY,1)=shi0(BCPOS(IB),IY,1) - UUOUT *
     .               ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .               SX(BCPOS(IB),IY)
                  else     ! xpb new pressure linearisation as per B2.5-style
                    shi0(BCPOS(IB),IY,3)=SHI0(BCPOS(IB),IY,3) - UUOUT *  !xpb
     .               (-ENIPAR(IB,1))*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .               SX(BCPOS(IB),IY)/ni(bcpos(ib),iy,0)                 !xpb
                  endif                                                  !xpb
                endif
              end do
            ENDIF
	      if(mdf_fhi.ne.0) then                                !srv 12.09.08 {
             shi0(BCPOS(IB),IY,0) = shi0(BCPOS(IB),IY,0) -
     &                  fhiPSch(rightix(BCPOS(IB),iy),iy,0)
	      endif                                                !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            IF(MOMPAR(ismain,IB,2).LT.0.5_R8) THEN
              if(sound_model.eq.0) then
                CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))  !srv 01.02.09
              else
                zarg = 0.0_R8
                do is = 0, ns-1
                  zarg = zarg + bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .             (bcint(rza(-1,-1,is),nx,ny,IY,IB)*
     .              bcint(te,nx,ny,IY,IB)+bcint(ti,nx,ny,IY,IB))
     .             /(am(is)*bcint(ni(-1,-1,0),nx,ny,IY,IB)*mp)
                enddo
                CS=SQRT(zarg)
              endif
              UUOUT=MOMPAR(ismain,IB,1)*CS*
     .         abs(wbbl(bcpos(ib),IY,0)/wbbl(bcpos(ib),iy,3))
              UPOUT=max(0.0_R8,UUOUT/
     .         abs(wbbl(bcpos(ib),IY,0)/wbbl(bcpos(ib),iy,3))-
     .         abs(wbbl(bcpos(ib),iy,2)/wbbl(bcpos(ib),iy,0))*
     .         ud(bcpos(ib),iy,1))
              DO IS=0,NS-1
                if (.not.is_neutral(is)) then
                  shi0(BCPOS(IB),IY,0)=shi0(BCPOS(IB),IY,0) -
     .             ENIPAR(IB,2)*0.5_R8*UUOUT*AM(IS)*mp*
     .             bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .             (UPOUT**2+UD(BCPOS(IB),IY,IS)**2)*
     .             SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))         !xpb
                  if(enipar(ib,1).gt.0.0_R8) then                        !xpb
                    shi0(BCPOS(IB),IY,1)=shi0(BCPOS(IB),IY,1) - UUOUT *
     .               ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .               SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
                  else     ! xpb new pressure linearisation as per B2.5-style
                    shi0(BCPOS(IB),IY,3)=SHI0(BCPOS(IB),IY,3) - UUOUT *  !xpb
     .               (-ENIPAR(IB,1))*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .                SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))/
     .                ni(bcpos(ib),iy,0)
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
                  UUTEST=
     .             UU(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY),IS)*
     1             (1.0_R8+
     2              HX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))/
     .       HX(LEFTIX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),
     .          LEFTIY(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))))-
     .       UU(LEFTIX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),
     .          LEFTIY(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),IS)*
     .              HX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))/
     .       HX(LEFTIX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),
     .          LEFTIY(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)))
                  UUOUT=MAX(MOMPAR(IS,IB,1)*CS*
     .             abs(wbbl(bcpos(ib),iy,0)/wbbl(bcpos(ib),iy,3)),
     .             UUTEST)
                  UPOUT=max(0.0_R8,UUOUT/
     .             abs(wbbl(bcpos(ib),iy,0)/wbbl(bcpos(ib),iy,3))-
     .             abs(wbbl(bcpos(ib),iy,2)/wbbl(bcpos(ib),iy,0))*
     .             ud(bcpos(ib),iy,is))
                  shi0(BCPOS(IB),IY,0)=shi0(BCPOS(IB),IY,0) -
     .             ENIPAR(IB,2)*0.5_R8*UUOUT*AM(IS)*mp*
     .             bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .             (UPOUT**2+UD(BCPOS(IB),IY,IS)**2)*
     .             SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))         !xpb
                  if(enipar(ib,1).gt.0.0_R8) then                        !xpb
                    shi0(BCPOS(IB),IY,1)=shi0(BCPOS(IB),IY,1) - UUOUT *
     .              ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .              SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
                  else     ! xpb new pressure linearisation as per B2.5-style
                    shi0(BCPOS(IB),IY,3)=SHI0(BCPOS(IB),IY,3) - UUOUT *  !xpb
     .               (-ENIPAR(IB,1))*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .                 SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))/
     .               ni(bcpos(ib),iy,0)
                  endif                                                  !xpb
                endif
              end do
            ENDIF
	      if(mdf_fhi.ne.0) then                                !srv 12.09.08 {
             shi0(BCPOS(IB),IY,0) = shi0(BCPOS(IB),IY,0) + 
     &                           fhiPSch(BCPOS(IB),iy,0)
	      endif                                                !srv 12.09.08 }
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=3, NO SHEATH IN RADIAL DIR.')
        END IF
        GOTO 40000

* -- BCENI=4 -- PRESCRIBE THE VALUE OF THE ION TEMPERATURE, WEAKLY
*                 A MIXED BOUNDARY CONDITION

30400   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'weak ion temperature',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2     BCEND(IB)
          call xertst(ENIPAR(IB,1).gt.0.0_R8,
     1     'BCENI = 4, ENIPAR(,1) !> 0')
          call xertst(ENIPAR(IB,2).gt.0.0_R8,                !srv 02.02.11
     1     'BCENI = 4, ENIPAR(,2) !> 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)= ENIPAR(IB,2)
     .       *ENIPAR(IB,1)*SY(IX,BCPOS(IB))
            shi0(IX,BCPOS(IB),1)=-ENIPAR(IB,2)/EV*SY(IX,BCPOS(IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)= ENIPAR(IB,2)*ENIPAR(IB,1)*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            shi0(IX,BCPOS(IB),1)=-ENIPAR(IB,2)/EV*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            shi0(BCPOS(IB),IY,0)= ENIPAR(IB,2)*
     1       ENIPAR(IB,1)*SX(BCPOS(IB),IY)
            shi0(BCPOS(IB),IY,1)=-ENIPAR(IB,2)/EV*SX(BCPOS(IB),IY)
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            shi0(BCPOS(IB),IY,0)= ENIPAR(IB,2)*ENIPAR(IB,1)*
     1       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            shi0(BCPOS(IB),IY,1)=-ENIPAR(IB,2)/EV*
     1       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=4, NOT YET AVAILABLE')
        END IF
        GOTO 40000

* -- BCENI=5 -- PRESCRIBE THE ION ENERGY FLUX PER UNIT AREA

30500   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'ion energy flux per area',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2     BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)=ENIPAR(IB,1)*SY(IX,BCPOS(IB))
            shi0(IX,BCPOS(IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)=ENIPAR(IB,1)*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            shi0(IX,BCPOS(IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            shi0(BCPOS(IB),IY,0)=ENIPAR(IB,1)*SX(BCPOS(IB),IY)
            shi0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            shi0(BCPOS(IB),IY,0)=ENIPAR(IB,1)*
     1       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            shi0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=5, NOT YET AVAILABLE')
        END IF
        GOTO 40000

* -- BCENI=6 -- PRESCRIBE THE TOTAL ION ENERGY FLUX FOR A
*                 CONSTANT ION TEMPERATURE

cghgh
30600   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'total ion energy flux',ENIPAR(IB,1),
     1     ' with const Ti on ',BCCHAR(IB),' at ',BCPOS(IB),
     2     ' from ',BCSTART(IB),' to ',BCEND(IB)
        endif
#ifdef WANT_THIS
        IF(BCCHAR(IB).EQ.'N') THEN
          Tguess=0.0_R8
          DO IX=BCSTART(IB), BCEND(IB)
            Tguess=Tguess+
     +       TI(bottomix(IX,BCPOS(IB)),bottomiy(ix,bcpos(ib)))
          ENDDO
          Tguess=Tguess/(BCEND(IB)-BCSTART(IB)+1)
          if(diagno.ge.2) write(*,*) 'Initial Tguess:',Tguess/ev
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
          DO IX=BCSTART(IB), BCEND(IB)
cxpb code segment analog to b2tfhi
            if(ncall.eq.0) then
             do is = 0, ns-1
              if(.not.is_neutral(is)) then
               vadia(ix,bcpos(ib),1,is) = 
     .          ti(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)))*
     .          bb(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)),2)
     &          /(qe*
     &            rza(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)),is)*
     *            hx(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib))))
     &          *(2.0_R8/
     &          bb(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)),3)**3)
     &          *(wbbr(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)),3)
     &          -wbbl(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)),3))
     &          *facdrift(ix,bcpos(ib))
               vaecrb(ix,bcpos(ib),1,is) = 1.0_R8/hx(ix,bcpos(ib))*
     1          bb(ix,bcpos(ib),2)/bb(ix,bcpos(ib),3)**2*0.5_R8*
     2         (po(rightix(ix,bcpos(ib)),rightiy(ix,bcpos(ib)))+
     2          po( leftix(ix,bcpos(ib)), leftiy(ix,bcpos(ib))))*
     &         fac_ExB(ix,bcpos(ib))
              else
               vadia(ix,bcpos(ib),1,is)=0.0_R8
               vaecrb(ix,bcpos(ib),1,is)=0.0_R8
              endif
             enddo
            endif
            if (conduction_only.eq.1) then
              flo = 0.0_R8
            else
             do is=0,ns-1
               prefactor=
     &          (vol(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)))+
     &           vol(ix,bcpos(ib)))
     &          /(hy1(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)))+
     &            hy1(ix,bcpos(ib)))
     &          *(na(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)),is)
     &           +na(ix,bcpos(ib),is))/2.0_R8
               cum_flo_vadia=cum_flo_vadia+
     1          prefactor*vadia(ix,bcpos(ib),1,is)
               cum_flo_vaecrb=cum_flo_vaecrb+
     1          prefactor*vaecrb(ix,bcpos(ib),1,is)
             enddo
              flo = (1.5_R8+BoRiS)*fni_32(ix,bcpos(ib),1)+
     1         2.5_R8*fni_52(ix,bcpos(ib),1)+
     &         chvi(ix,bcpos(ib),1)*
     &         (ni(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)),0)+
     &          ni(ix,bcpos(ib),0))/2.0_R8+
     1         2.5_R8*cum_flo_vadia+1.5_R8*cum_flo_vaecrb
            endif
            if (no_hybr_i.eq.0) then
              coni0 = hybr(flo,chci(ix,bcpos(ib),1))
              fy = flo*
     .         (ti(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)))
     &         +Tguess)/2.0_R8-coni0*
     1         (Tguess-ti(bottomix(ix,bcpos(ib)),
     2                    bottomiy(ix,bcpos(ib))))
            else
              coni0 = chci(ix,bcpos(ib),1)
              t2 = dim ( coni0,
     .         abs(flo) *
     .         max(0.0_R8, 
     .             min(Tguess,ti(bottomix(ix,bcpos(ib)),
     .                           bottomiy(ix,bcpos(ib)))))/
     .         max(cutlo, 
     .                 Tguess+ti(bottomix(ix,bcpos(ib)),
     .                           bottomiy(ix,bcpos(ib)))) )
              fy = 
     .         upwind(flo, 
     .                ti(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib))), 
     .                Tguess)
     .         - t2 * (Tguess-ti(bottomix(ix,bcpos(ib)),
     .                           bottomiy(ix,bcpos(ib))))
            endif
            fyt=fyt + fy
          ENDDO
          if(diagno.ge.2)
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
          if(diagno.ge.1)
     >     write(*,'(a,f15.6,i4,1pe14.6)') ' tibound:',T1/ev,ibnd,fyt
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)= VBIG*NE00*T1*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            shi0(IX,BCPOS(IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'S') THEN
          Tguess=0.0_R8
          DO IX=BCSTART(IB), BCEND(IB)
            Tguess=Tguess+TI(IX,BCPOS(IB))
          ENDDO
          Tguess=Tguess/(BCEND(IB)-BCSTART(IB)+1)
          if(diagno.ge.2) write(*,*) 'Initial Tguess:',Tguess/ev
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
          DO IX=BCSTART(IB), BCEND(IB)
cxpb code segment analog to b2tfhi
            if(ncall.eq.0) then
             do is = 0, ns-1
              if(.not.is_neutral(is)) then
               vadia(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1,is) =
     &          ti(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))*
     &          bb(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),2)/
     &          (qe*rza(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),is)*
     *           hx(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))))
     &          *(2.0_R8/
     &          bb(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),3)**3)
     &          *(wbbr(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),3)-
     &            wbbl(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),3))
     &          *facdrift(ix,bcpos(ib))
               vaecrb(ix,bcpos(ib),1,is) = 1.0_R8/hx(ix,bcpos(ib))*
     1          bb(ix,bcpos(ib),2)/bb(ix,bcpos(ib),3)**2*0.5_R8*
     2         (po(rightix(ix,bcpos(ib)),rightiy(ix,bcpos(ib)))+
     2          po( leftix(ix,bcpos(ib)), leftiy(ix,bcpos(ib))))*
     &         fac_ExB(ix,bcpos(ib))
              else
               vadia(ix,bcpos(ib),1,is)=0.0_R8
               vaecrb(ix,bcpos(ib),1,is)=0.0_R8
              endif
             enddo
            endif
            if (conduction_only.eq.1) then
              flo = 0.0_R8
            else
             do is=0,ns-1
               prefactor=
     &          (vol(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))+
     &           vol(ix,bcpos(ib)))
     &          /(hy1(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))+
     &            hy1(ix,bcpos(ib)))
     &          *(na(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),is)
     &           +na(ix,bcpos(ib),is))/2.0_R8
               cum_flo_vadia=cum_flo_vadia+
     1          prefactor*vadia(ix,bcpos(ib),1,is)
               cum_flo_vaecrb=cum_flo_vaecrb+
     1          prefactor*vaecrb(ix,bcpos(ib),1,is)
             enddo
             flo = (1.5_R8+BoRiS)*
     1         fni_32(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1)+
     1         2.5_R8*fni_52(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1)+
     &         chvi(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1)*
     &         (ni(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),0)+
     &          ni(ix,bcpos(ib),0))/2.0_R8+
     1         2.5_R8*cum_flo_vadia+1.5_R8*cum_flo_vaecrb
            endif
            if (no_hybr_i.eq.0) then
              coni0 = hybr(flo,
     &         chci(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1))
              fy = flo*(ti(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))
     &                 +Tguess)/2.0_R8-coni0*
     1         (ti(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))-Tguess)
            else
              coni0 = chci(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1)
              t2 = dim ( coni0,
     .         abs(flo) *
     .         max(0.0_R8, 
     .             min(Tguess,ti(topix(ix,bcpos(ib)),
     .                           topiy(ix,bcpos(ib)))))/
     .         max(cutlo, 
     .                 Tguess+ti(topix(ix,bcpos(ib)),
     .                           topiy(ix,bcpos(ib)))) )
              fy = 
     .         upwind(flo, Tguess,
     .                ti(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))))
     .         - t2 * (ti(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))-
     .                 Tguess)
            endif
            fyt=fyt + fy
          ENDDO
          if(diagno.ge.2)
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
          if(diagno.ge.1)
     >     write(*,'(a,f15.6,i4,1pe14.6)') ' tibound:',T1/ev,ibnd,fyt

          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)= VBIG*NE00*T1*SY(IX,BCPOS(IB))
            shi0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=6, NOT AVAILABLE FOR W,E')
        ENDIF
#else
        CALL XERRAB ('B2STBC -- BCENI=6, NOT IMPLEMENTED.')
#endif
        GOTO 40000
cghgh

* -- BCENI=7 -- PRESCRIBE THE ION TEMPERATURE AS A FUNCTION OF OTHER
*                 PLASMA PARAMETERS

30700   CONTINUE
        CALL XERRAB ('B2STBC -- BCENI=7, NOT YET AVAILABLE')
        GOTO 40000

* -- BCENI=8 -- PRESCRIBE THE TOTAL ION HEAT FLUX WITH CONSTANT
*                 FLUX DENSITY

30800   CONTINUE
        if(ncall.eq.0) then
#ifdef ASTRA
          ENIPAR(IB,1) = fhib_astra                              !srv 28.07.11
#endif
          write(*,*) 'total ion energy flux',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2     BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          US=0.0_R8
          DO IX=BCSTART(IB), BCEND(IB)
            US=US+SY(IX,BCPOS(IB))
          end do
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)=ENIPAR(IB,1)/US*SY(IX,BCPOS(IB))
            shi0(IX,BCPOS(IB),1)=0.0_R8
	      if(mdf_fhi.ne.0) then                                !srv 12.09.08 {
             shi0(ix,bcpos(ib),0) = shi0(ix,bcpos(ib),0) -
     &                           fhiPSch(ix,topiy(ix,bcpos(ib)),1)
             do is=0,ns-1
              shi0(ix,bcpos(ib),1) = shi0(ix,bcpos(ib),1) +
     &            1.5_R8*vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &           vaecrb(ix,topiy(ix,bcpos(ib)),1,is)*
     &           na(ix,bcpos(ib),is)
             end do
	      endif                                                !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          US=0.0_R8
          DO IX=BCSTART(IB), BCEND(IB)
            US=US+SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          end do
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)=ENIPAR(IB,1)/US*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            shi0(IX,BCPOS(IB),1)=0.0_R8
	      if(mdf_fhi.ne.0) then                                !srv 12.09.08 {
             shi0(ix,bcpos(ib),0) = shi0(ix,bcpos(ib),0) +
     &                           fhiPSch(ix,bcpos(ib),1)
             do is=0,ns-1
              shi0(ix,bcpos(ib),1) = shi0(ix,bcpos(ib),1) -
     &          1.5_R8*vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &          vaecrb(ix,bcpos(ib),1,is)*na(ix,bcpos(ib),is)
             end do
	      endif                                                !srv 12.09.08 }
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          US=0.0_R8
          DO IY=BCSTART(IB), BCEND(IB)
            US=US+SX(BCPOS(IB),IY)
          end do
          DO IY=BCSTART(IB), BCEND(IB)
            shi0(BCPOS(IB),IY,0)=ENIPAR(IB,1)/US*SX(BCPOS(IB),IY)
            shi0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          US=0.0_R8
          DO IY=BCSTART(IB), BCEND(IB)
            US=US+SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          end do
          DO IY=BCSTART(IB), BCEND(IB)
            shi0(BCPOS(IB),IY,0)=ENIPAR(IB,1)/US*
     1       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            shi0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=8, NOT AVAILABLE')
        ENDIF
        GOTO 40000

* -- BCENI=9 -- PRESCRIBE THE DECAY LENGTH FOR THE ION TEMPERATURE

30900   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'ion temperature decay length',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2     BCEND(IB)
          call xertst(ENIPAR(IB,1).ne.0.0_R8,
     1     'BCENI = 9, ENIPAR(,1) = 0')
!          call xertst(all(BCCON(0:NS-1,IB).EQ.9),
!     1     'Should use BCCON(IS,IB).EQ.9 together with BCENI=9 ')       !srv 16.07.12
          write(*,*) 'Please be assured that this boundary condition',
     &               'is consistent with density boundary condition'     !srv 21.09.12
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)=0.0_R8
            shi0(IX,BCPOS(IB),1)=-hci0(IX,BCPOS(IB))/ENIPAR(IB,1)*
     1       SY(IX,BCPOS(IB))
            DO IS=0,NS-1
cxpb          if (.not.is_neutral(is)) then   ! dpc --- need to check
                shi0(IX,BCPOS(IB),0)=shi0(IX,BCPOS(IB),0) +
     .           2.5_R8*TI(IX,BCPOS(IB))*SNA0_NO_MDF(IX,BCPOS(IB),0,IS) !srv 23.09.08
                shi0(IX,BCPOS(IB),1)=shi0(IX,BCPOS(IB),1) +
     .           2.5_R8*SNA0_NO_MDF(IX,BCPOS(IB),1,IS)*                 !srv 23.09.08
     &                           na(IX,BCPOS(IB),IS)
cxpb          endif
            ENDDO
 	      if(mdf_fhi.ne.0) then                                !srv 12.09.08 {
             shi0(ix,bcpos(ib),0) = shi0(ix,bcpos(ib),0) -
     &                              fhiPSch(ix,topiy(ix,bcpos(ib)),1)
             do is=0,ns-1
              shi0(ix,bcpos(ib),1) = shi0(ix,bcpos(ib),1) +
     &          1.5_R8*vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &          vaecrb(ix,topiy(ix,bcpos(ib)),1,is)*na(ix,bcpos(ib),is)
             end do
	      endif                                                !srv 12.09.08 }
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)=0.0_R8
            shi0(IX,BCPOS(IB),1)=-hci0(IX,BCPOS(IB))/ENIPAR(IB,1)*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            DO IS=0,NS-1
cxpb          if (.not.is_neutral(is)) then   ! dpc --- need to check
cxpb            shi0(IX,BCPOS(IB),0)=shi0(IX,BCPOS(IB),0) +
cxpb     .           rza(ix,bcpos(ib),is)*2.5_R8*TI(IX,BCPOS(IB))*
cxpb     .           SNA0_NO_MDF(IX,BCPOS(IB),0,IS)
cxpb            shi0(IX,BCPOS(IB),1)=shi0(IX,BCPOS(IB),1) +
cxpb     .           rza(ix,bcpos(ib),is)*2.5_R8*SNA0_NO_MDF(IX,BCPOS(IB),1,IS)*
cxpb     .           na(IX,BCPOS(IB),IS)
cxpb          endif
                 shi0(IX,BCPOS(IB),0)=shi0(IX,BCPOS(IB),0) +
     .           2.5_R8*TI(IX,BCPOS(IB))*SNA0_NO_MDF(IX,BCPOS(IB),0,IS) !srv 23.09.08
                shi0(IX,BCPOS(IB),1)=shi0(IX,BCPOS(IB),1) +
     .           2.5_R8*SNA0_NO_MDF(IX,BCPOS(IB),1,IS)*                 !srv 23.09.08
     &                           na(IX,BCPOS(IB),IS)
            ENDDO
	      if(mdf_fhi.ne.0) then                                !srv 12.09.08 {
              shi0(ix,bcpos(ib),0) = shi0(ix,bcpos(ib),0) +
     &                              fhiPSch(ix,bcpos(ib),1)
             do is=0,ns-1
               shi0(ix,bcpos(ib),1) = shi0(ix,bcpos(ib),1) -
     &          1.5_R8*vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*
     &          vaecrb(ix,bcpos(ib),1,is)*na(ix,bcpos(ib),is)
             end do
	      endif                                                !srv 12.09.08 }
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=9, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 40000

* -- BCENI=10 -- Feedback option for core

31000   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'ion heat flux feedback bc on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),
     2     ' from ',BCSTART(IB),' to ',BCEND(IB)
          call xertst(lfeedback, 'LFEEDBACK must be turned on!')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            if (region(ix,bcpos(ib),0).eq.1) then
              shi0(IX,BCPOS(IB),0)=
     1         cbshi(0,ismain,coreregno)*SY(IX,BCPOS(IB))
            elseif (region(ix,bcpos(ib),0).eq.5) then
              shi0(IX,BCPOS(IB),0)=
     1         cbshi(0,ismain,coreregn2)*SY(IX,BCPOS(IB))
            endif
            shi0(IX,BCPOS(IB),1)=0.0_R8
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=10, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 40000

* -- BCENI=11 -- SHEATH CONDITIONS, ION ENERGY TRANSMISSION COEFFICIENT

31100   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'sheath (11) ion on ',BCCHAR(IB),
     1     ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            DO IS=0,NS-1
              if (.not.is_neutral(is)) then
                ZBAR=bcint(na(-1,-1,IS),nx,ny,IY,IB)/
     .           bcint(NE,nx,ny,IY,IB)
                ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .           bcint(rz2(-1,-1,IS),nx,ny,IY,IB)*ZBAR*
     .           bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                CS=SQRT(ZARG)
                UUTEST=
     .           UU(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY),IS)*
     .           (1.0_R8+
     .            HX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))/
     .    HX(RIGHTIX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),
     .       RIGHTIY(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))))-
     .    UU(RIGHTIX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),
     .       RIGHTIY(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),IS)*
     .            HX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))/
     .    HX(RIGHTIX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)),
     .       RIGHTIY(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)))
                UUOUT=MAX(MOMPAR(IS,IB,1)*CS*
     .           abs(wbbr(bcpos(ib),iy,0)/wbbr(bcpos(ib),iy,3)),
     .           -UUTEST)
                UPOUT=max(0.0_R8,UUOUT/
     .           abs(wbbr(bcpos(ib),iy,0)/wbbr(bcpos(ib),iy,3))+
     .           abs(wbbr(bcpos(ib),iy,2)/wbbr(bcpos(ib),iy,0))*
     .           ud(bcpos(ib),iy,is))
                shi0(BCPOS(IB),IY,0)=shi0(BCPOS(IB),IY,0) -
     .           ENIPAR(IB,2)*0.5_R8*UUOUT*AM(IS)*mp*
     .           bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .           (UPOUT**2+UD(BCPOS(IB),IY,IS)**2)*SX(BCPOS(IB),IY)
                shi0(BCPOS(IB),IY,1)=shi0(BCPOS(IB),IY,1) - UUOUT *
     .           ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .           SX(BCPOS(IB),IY)
              endif
            end do
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            DO IS=0,NS-1
              if (.not.is_neutral(is)) then
                ZBAR=bcint(na(-1,-1,IS),nx,ny,IY,IB)/
     .           bcint(NE,nx,ny,IY,IB)
                ZARG=(GAMMAI*bcint(TI,nx,ny,IY,IB)+
     .           bcint(rz2(-1,-1,IS),nx,ny,IY,IB)*ZBAR*
     .           bcint(TE,nx,ny,IY,IB))/(AM(IS)*mp)
                CS=SQRT(ZARG)
                UUTEST=UU(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY),IS)*
     1        (1.0_R8+HX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))/
     .         HX(LEFTIX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),
     .            LEFTIY(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))))-
     .         UU(LEFTIX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),
     .            LEFTIY(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),IS)*
     .                HX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))/
     .         HX(LEFTIX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)),
     .            LEFTIY(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)))
                UUOUT=MAX(MOMPAR(IS,IB,1)*CS*
     .           abs(wbbl(bcpos(ib),iy,0)/wbbl(bcpos(ib),iy,3)),
     .           UUTEST)
                UPOUT=max(0.0_R8,UUOUT/
     .           abs(wbbl(bcpos(ib),iy,0)/wbbl(bcpos(ib),iy,3))-
     .           abs(wbbl(bcpos(ib),iy,2)/wbbl(bcpos(ib),iy,0))*
     .           ud(bcpos(ib),iy,is))
                shi0(BCPOS(IB),IY,0)=shi0(BCPOS(IB),IY,0) -
     .           ENIPAR(IB,2)*0.5_R8*UUOUT*AM(IS)*mp*
     .           bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .           (UPOUT**2+UD(BCPOS(IB),IY,IS)**2)*
     .           SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))             !xpb
                shi0(BCPOS(IB),IY,1)=shi0(BCPOS(IB),IY,1) - UUOUT *
     .           ENIPAR(IB,1)*bcint(na(-1,-1,IS),nx,ny,IY,IB)*
     .           SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
              endif
            end do
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=11, NO SHEATH IN RADIAL DIR.')
        END IF
        GOTO 40000

* -- BCENI=12 -- SHEATH CONDITIONS, ION ENERGY TRANSMISSION COEFFICIENT

31200   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'sheath (12) ion on ',BCCHAR(IB),
     1     ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            DO IS=0,NS-1
              shi0(ix,bcpos(ib),3) = shi0(ix,bcpos(ib),3)+
     1         ENIPAR(IB,1)*max(0.0_R8,
     1          -fna(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),1,is))/
     2         ni(ix,bcpos(ib),0)
            end do
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            DO IS=0,NS-1
              shi0(ix,bcpos(ib),3) = shi0(ix,bcpos(ib),3)+
     1         ENIPAR(IB,1)*max(0.0_R8,fna(ix,bcpos(ib),1,is))/
     2         ni(ix,bcpos(ib),0)
            end do
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            DO IS=0,NS-1
              shi0(bcpos(ib),iy,3) = shi0(bcpos(ib),iy,3)+
     1         ENIPAR(IB,1)*max(0.0_R8,
     1          -fna(rightix(bcpos(ib),iy),rightiy(bcpos(ib),iy),0,is))/
     2         ni(bcpos(ib),iy,0)
            end do
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            DO IS=0,NS-1
              shi0(bcpos(ib),iy,3) = shi0(bcpos(ib),iy,3)+
     1         ENIPAR(IB,1)*max(0.0_R8,fna(bcpos(ib),iy,0,is))/
     2         ni(bcpos(ib),iy,0)
            end do
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=12, INVALID OPTION')
        END IF
        GOTO 40000

* -- BCENI=13 -- PRESCRIBE THE ION ENERGY FLUX PER UNIT AREA
*             -- PROPORTIONAL TO TEMPERATURE

31300   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'ion energy flux per area proportional to ',
     1     'temperature ',ENIPAR(IB,1),' on ',BCCHAR(IB),' at ',
     2     BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)=ENIPAR(IB,1)*TI(IX,BCPOS(IB))*
     1       SY(IX,BCPOS(IB))
            shi0(IX,BCPOS(IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)=ENIPAR(IB,1)*TI(IX,BCPOS(IB))*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            shi0(IX,BCPOS(IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            shi0(BCPOS(IB),IY,0)=ENIPAR(IB,1)*TI(BCPOS(IB),IY)*
     1       SX(BCPOS(IB),IY)
            shi0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            shi0(BCPOS(IB),IY,0)=ENIPAR(IB,1)*TI(BCPOS(IB),IY)*
     1       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            shi0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=13, NOT YET AVAILABLE')
        END IF
        GOTO 40000

* -- BCENI=14 -- leakage option for ion energy

31400   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'ion energy leakage',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2     ' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            SHI0(IX,BCPOS(IB),0)=0.0_R8
            if(sound_model.eq.0) then
             CS=SQRT(bcint(PR,nx,ny,IX,IB)/bcint(RO,nx,ny,IX,IB))  !srv 01.02.09
            else
               zarg = 0.0_R8
               do is = 0, ns-1
                  zarg = zarg + bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .                 (bcint(rza(-1,-1,is),nx,ny,IY,IB)*
     .                  bcint(te,nx,ny,IY,IB)+bcint(ti,nx,ny,IY,IB))
     .                 /(am(is)*bcint(ni(-1,-1,0),nx,ny,IY,IB)*mp)
               enddo
               CS=SQRT(zarg)
            endif
            SHI0(IX,BCPOS(IB),1)=
     1       vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*cs*ENIPAR(IB,1)
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            SHI0(IX,BCPOS(IB),0)=0.0_R8
            if(sound_model.eq.0) then
              CS=SQRT(bcint(PR,nx,ny,IX,IB)/bcint(RO,nx,ny,IX,IB))  !srv 01.02.09
            else
               zarg = 0.0_R8
               do is = 0, ns-1
                  zarg = zarg + bcint(na(-1,-1,is),nx,ny,IY,IB)*
     .                 (bcint(rza(-1,-1,is),nx,ny,IY,IB)*
     .                  bcint(te,nx,ny,IY,IB)+bcint(ti,nx,ny,IY,IB))
     .                 /(am(is)*bcint(ni(-1,-1,0),nx,ny,IY,IB)*mp)
               enddo
               CS=SQRT(zarg)
            endif
            SHI0(IX,BCPOS(IB),1)=
     1       vol(ix,bcpos(ib))/hy1(ix,bcpos(ib))*cs*ENIPAR(IB,1)
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=14, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 40000

* -- BCENI=15 -- from b2stbc_spb                              !srv 01.02.09 {

31500   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'from b2stbc_spb for BCENI with enipar1 =',
     &     ENIPAR(IB,1),' and enipar2 =', ENIPAR(IB,2),' on ',
     &     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     &     ' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
!            SHI0(BCPOS(IB),IY,0)=0.0_R8
!            CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))
            CS=SQRT(PR(BCPOS(IB),IY)/RO(BCPOS(IB),IY))
	      t0 = 0.0_R8
            do is = 0, ns-1
             if(.not.is_neutral(is)) t0 = t0 + na(BCPOS(IB),iy,is)       !srv 11.03.09
            enddo
            SHI0(BCPOS(IB),IY,1) = SHI0(BCPOS(IB),IY,1) +                !srv 11.03.09
     &           ENIPAR(IB,1)*t0*cs*abs(pbs(rightix(BCPOS(IB),iy),iy,0)) !srv 26.11.02 28.01.03
            if(mdf_fhi.ne.0) then                                        !srv 13.10.06 {
             shi0(BCPOS(IB),iy,0) = shi0(BCPOS(IB),iy,0) - 
     &               fhiPSch(rightix(BCPOS(IB),iy),iy,0)
	      endif                                                         !srv 13.10.06 }
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
!            SHI0(BCPOS(IB),IY,0)=0.0_R8
!            CS=SQRT(bcint(PR,nx,ny,IY,IB)/bcint(RO,nx,ny,IY,IB))
            CS=SQRT(PR(BCPOS(IB),IY)/RO(BCPOS(IB),IY))
	      t0 = 0.0_R8
            do is = 0, ns-1
             if(.not.is_neutral(is)) t0 = t0 + na(BCPOS(IB),iy,is)   !srv 11.03.09
            enddo
            SHI0(BCPOS(IB),IY,1) = SHI0(BCPOS(IB),IY,1) +            !srv 11.03.09
     &           ENIPAR(IB,1)*t0*cs*abs(pbs(BCPOS(IB),iy,0))         !srv 26.11.02 28.01.03
            if(mdf_fhi.ne.0) then                                    !srv 13.10.06 {
             shi0(BCPOS(IB),iy,0) = shi0(BCPOS(IB),iy,0) + 
     &               fhiPSch(BCPOS(IB),iy,0)
	      endif                                                    !srv 13.10.06 }
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=15, NOT AVAILABLE FOR S,N')
        ENDIF
        GOTO 40000                                                   !srv 01.02.09 }

* -- BCENI=16 -- Feedback boundary condition with constant temperature

31600   CONTINUE
        if(bceni_16_style.eq.0) then 
         if(ncall.eq.0) then
           write(*,*) 'ion temperature feedback by flux',
     1      ENIPAR(IB,1),' with ',ENIPAR(IB,2),' on ',
     1      BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2      ' to ',BCEND(IB)
           call xertst(enipar(ib,1).ne.0, 'faulty parameter enipar1')
           call xertst(enipar(ib,2).gt.0, 'faulty parameter enipar2')
         endif
         IF(BCCHAR(IB).EQ.'S') THEN
           totflux = 0.0_R8
           t0 = 0.0_R8
           weight = 0.0_R8
           DO IX=BCSTART(IB), BCEND(IB)
             bc_type16_ix=ix
             bc_type16_iy=BCPOS(IB)
             do bc_type16_i=1,bc_type16_ref      ! how many steps away from the boundary?
               bc_type16_ixx=topix(bc_type16_ix,bc_type16_iy)
               bc_type16_iy=topiy(bc_type16_ix,bc_type16_iy)
               bc_type16_ix=bc_type16_ixx
               call xertst (-2.lt.bc_type16_ix.and.bc_type16_ix.lt.nx+1,
     1          'stepped off x-grid in type 16 bc')
               call xertst (-2.lt.bc_type16_iy.and.bc_type16_iy.lt.ny+1,
     1          'stepped off y-grid in type 16 bc')
             enddo
             totflux = totflux + 
     .        fhi(bc_type16_ix,bc_type16_iy,1)
             if (bc_type16_kinetic_energy.eq.1) then
               do is = 0, ns-1
                 totflux=totflux+(1.0_R8-BoRiS)*
     .            fna(bc_type16_ix,bc_type16_iy,1,is)*0.5_R8*
     .            (kinrgy(bc_type16_ix,bc_type16_iy,is)+
     .             kinrgy(bottomix(bc_type16_ix,bc_type16_iy),
     .                    bottomiy(bc_type16_ix,bc_type16_iy),is))
               enddo
             endif
             t0 = t0 + 
     .        ti(ix,bcpos(ib))*vol(ix,bcpos(ib))*ni(ix,bcpos(IB),0)
             weight = weight + vol(ix,bcpos(ib))*ni(ix,bcpos(IB),0)
           ENDDO
           t0 = t0 / weight
           if(diagno.ge.2) write(*,*) 'Totflux :',totflux
           T1 = (1.0_R8 - 
     .      enipar(ib,2)*(totflux-enipar(ib,1))/enipar(ib,1)) * t0 
           if(diagno.ge.2)
     1      write(*,*) 'Old/New ion temp :',t0/qe, t1/qe
           DO IX=BCSTART(IB), BCEND(IB)
             shi0(IX,BCPOS(IB),0)= VBIG*NE00*T1*SY(IX,BCPOS(IB))
             shi0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
           ENDDO
         ELSE
           CALL XERRAB ('B2STBC -- BCENI=16, NOT AVAILABLE FOR N,W,E')
         ENDIF
         GOTO 40000

        elseif(bceni_16_style.eq.1)then                                 !srv 16.02.12   {

          if(ncall.eq.0) then
	       write(*,*) 'Given integrated ion heat flux',
     &                  ' with constant Ti '
!          write all boundaries with this type of boundary condition
           DO IBW=IB,NBC
		if(BCCHAR(IB).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN	       
		   write(*,*) ' on ',BCCHAR(IBW),' at ',BCPOS(IBW),' from ',
     2          BCSTART(IBW),' to ', BCEND(IBW)
              call xertst(enipar(ibw,1).ne.0,'faulty parameter enipar1')
	      endif
	     ENDDO
#ifdef ASTRA
           SENIPAR = fhib_astra
#else
           SENIPAR = 0.0_R8
           DO IBW=IB,NBC
	      if(BCCHAR(IB).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN	       
               SENIPAR = SENIPAR + ENIPAR(IBW,1)
	      endif
	     ENDDO
#endif
		 write(*,*) ' Summary value is ', SENIPAR

	     TIS = TI((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB))
          endif
	    if(.not.once_bceni_16) then
           IF(BCCHAR(IB).EQ.'S') THEN
	      totflux = 0.0_R8
            DO IBW=IB,NBC
             IF(BCCHAR(IBW).EQ.'S'.AND.BCENI(IBW).EQ.BCENI(IB)) THEN
	        if(mdf_fhi.ne.0) then
               DO IX=BCSTART(IBW), BCEND(IBW)
	          totflux = totflux + 
     &          fhi_mdf(topix(ix,BCPOS(IBW)),topiy(ix,BCPOS(IBW)),1)
     &       +fac_ExB(ix,BCPOS(IBW))*                                    !srv 11.04.13 {
     &      ((po(ix,BCPOS(IBW))-po(leftix(ix,BCPOS(IBW)),BCPOS(IBW)))/
     &       (hx(ix,BCPOS(IBW))+hx(leftix(ix,BCPOS(IBW)),BCPOS(IBW)))+
     &       (po(rightix(ix,BCPOS(IBW)),BCPOS(IBW))-po(ix,BCPOS(IBW)))/
     &       (hx(rightix(ix,BCPOS(IBW)),BCPOS(IBW))+hx(ix,BCPOS(IBW))))*
     &        ni(ix,BCPOS(IBW),1)*ti(ix,BCPOS(IBW))*SY(IX,BCPOS(IBW))*
     &        bb(ix,BCPOS(IBW),2)/bb(ix,BCPOS(IBW),3)**2                 !srv 11.04.13 }
              ENDDO
	        else
               DO IX=BCSTART(IBW), BCEND(IBW)
	          totflux = totflux + 
     &          fhi(topix(ix,BCPOS(IBW)),topiy(ix,BCPOS(IBW)),1) 
               ENDDO
	        endif
	       ENDIF
	      ENDDO
            TIS = TI((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB)) * 
     &         (1.0_R8-bc_ref_ti*(totflux-SENIPAR)/                      !srv 11.04.13
     &           max(abs(totflux),SENIPAR))
             write(*,"(' BCENI=16: TIS=',es14.7,' totflux=',es14.7,
     *       ' SENIPAR=',es14.7)") TIS/ev, totflux, SENIPAR
            DO IBW=IB,NBC
             IF(BCCHAR(IBW).EQ.'S'.AND.BCENE(IBW).EQ.BCENE(IB)) THEN
              DO IX=BCSTART(IBW), BCEND(IBW)
               shi0(IX,BCPOS(IBW),0)= VBIG*NE00*TIS*SY(IX,BCPOS(IBW))
               shi0(IX,BCPOS(IBW),1)=-VBIG*NE00*SY(IX,BCPOS(IBW))
              ENDDO
	       ENDIF
	      ENDDO
		  once_bceni_16 = .true.
           ELSE
             CALL XERRAB
     1         ('B2STBC -- BCENI=16, NOT AVAILABLE FOR N,W,E')
           ENDIF
	    endif                                                          !srv 16.02.12 }
          GOTO 40000
	  else
	   write(*,*) ' Incorrect style of boundary condition'
	   stop
	  endif
        

* -- BCENI=17 -- Feedback boundary condition with constant temperature

31700   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'shared temperature feedback by flux',
     1     ENIPAR(IB,1),' with ',ENIPAR(IB,2),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2     ' to ',BCEND(IB)
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
          CALL XERRAB ('B2STBC -- BCENI=16, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 40000

* -- BCENI=18 -- FRACTIONAL DROP CONDITION (FOR EDGE2D COMPARISON)

31800   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'Ion temperature fractional drop',
     1     ENIPAR(IB,1),' on ',BCCHAR(IB),' at ',BCPOS(IB),
     2     ' from ',BCSTART(IB),' to ',BCEND(IB)
          call xertst(ENIPAR(IB,1).lt.1.0_R8,
     1     'BCENI = 18, ENIPAR(,1) => 1')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            T1 = (1.0_R8 - ENIPAR(IB,1)) *
     .       TI(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)))
            shi0(IX,BCPOS(IB),0)= VBIG*NE00*T1*SY(IX,BCPOS(IB))
            shi0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            T1 = (1.0_R8 - ENIPAR(IB,1)) *
     .       TI(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            shi0(IX,BCPOS(IB),0)= VBIG*NE00*T1*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            shi0(IX,BCPOS(IB),1)=-VBIG*NE00*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            T1 = (1.0_R8 - ENIPAR(IB,1)) *
     .       TI(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))
            shi0(BCPOS(IB),IY,0)= VBIG*NE00*T1*SX(BCPOS(IB),IY)
            shi0(BCPOS(IB),IY,1)=-VBIG*NE00*SX(BCPOS(IB),IY)
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            T1 = (1.0_R8 - ENIPAR(IB,1)) *
     .       TI(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            shi0(BCPOS(IB),IY,0)= VBIG*NE00*T1*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            shi0(BCPOS(IB),IY,1)=-VBIG*NE00*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=18, INVALID OPTION')
        ENDIF
        GOTO 40000

* -- BCENI=19 -- PRESCRIBE THE DECAY LENGTH FOR THE ION TEMPERATURE

31900   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'ion temperature decay length (new)',ENIPAR(IB,1),
     .     ' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2     BCEND(IB)
          call xertst(ENIPAR(IB,1).ne.0.0_R8,
     1     'BCENI = 19, ENIPAR(,1) = 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            shi0(IX,BCPOS(IB),0)=0.0_R8
            shi0(IX,BCPOS(IB),1)=-hci0(IX,BCPOS(IB))/ENIPAR(IB,1)*
     1       SY(IX,BCPOS(IB))
            t0=0.0_R8
            DO IS=0,NS-1
               t0=t0-
     .              1.5_R8*fna_32(topix(ix,bcpos(ib)),
     .                            topiy(ix,bcpos(ib)),1,is)-
     .              2.5_R8*fna_52(topix(ix,bcpos(ib)),
     .                            topiy(ix,bcpos(ib)),1,is)
            ENDDO
            if(t0.gt.0.0_R8) then
               shi0(IX,BCPOS(IB),0)=shi0(IX,BCPOS(IB),0) + 
     .              t0 * TI(IX,BCPOS(IB))
            else
               shi0(IX,BCPOS(IB),1)=shi0(IX,BCPOS(IB),1) + t0
            endif
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
             shi0(IX,BCPOS(IB),0)=0.0_R8
             shi0(IX,BCPOS(IB),1)=-hci0(IX,BCPOS(IB))/ENIPAR(IB,1)*
     1            SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
             t0=0.0_R8
             DO IS=0,NS-1
                t0=t0+
     .               1.5_R8*fna_32(ix,bcpos(ib),1,is)+
     .               2.5_R8*fna_52(ix,bcpos(ib),1,is)
             ENDDO
            if(t0.gt.0.0_R8) then
               shi0(IX,BCPOS(IB),0)=shi0(IX,BCPOS(IB),0) + 
     .              t0 * TI(IX,BCPOS(IB))
             else
                shi0(IX,BCPOS(IB),1)=shi0(IX,BCPOS(IB),1) + t0
             endif
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=19, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 40000

* -- BCENI=20 -- Feedback boundary condition with constant temperature scaled by temperature on the ring bc_type20_ref away

32000   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'ion temperature feedback by flux (neo)',
     1     ENIPAR(IB,1),' with ',ENIPAR(IB,2),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2     ' to ',BCEND(IB)
          call xertst(enipar(ib,1).ne.0, 'faulty parameter enipar1')
          call xertst(enipar(ib,2).gt.0, 'faulty parameter enipar2')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          totflux = 0.0_R8
          t0 = 0.0_R8
          weight = 0.0_R8
	  tiave = 0.0_R8
	  niave = 0.0_R8
          DO IX=BCSTART(IB), BCEND(IB)
            bc_type16_ix=ix
            bc_type16_iy=BCPOS(IB)
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
     .       fhi(bc_type16_ix,bc_type16_iy,1)
            t0 = t0 + 
     .       ti(ix,bcpos(ib))*vol(ix,bcpos(ib))*ni(ix,bcpos(IB),0)
            weight = weight + vol(ix,bcpos(ib))*ni(ix,bcpos(IB),0)
            bc_type20_ix=ix
            bc_type20_iy=BCPOS(IB)
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
     1	     ni(bc_type20_ix,bc_type20_iy,0)*
     &          vol(bc_type20_ix,bc_type20_iy)
	    niave=niave+ni(bc_type20_ix,bc_type20_iy,0)*
     1	     vol(bc_type20_ix,bc_type20_iy)
          ENDDO
          t0 = t0 / weight
	  tiave=tiave/niave
          if(diagno.ge.2) write(*,*) 'Totflux :',totflux
          T1 = (1.0_R8 - 
     .     enipar(ib,2)*(totflux-enipar(ib,1))/enipar(ib,1)) * t0 
          if(diagno.ge.2)
     1     write(*,*) 'Old/New ion temp :',t0/qe, t1/qe
          DO IX=BCSTART(IB), BCEND(IB)
            bc_type20_ix=ix
            bc_type20_iy=BCPOS(IB)
            do bc_type20_i=1,bc_type20_ref      ! how many steps away from the boundary?
	      bc_type20_ixx=topix(bc_type20_ix,bc_type20_iy)
	      bc_type20_iy=topiy(bc_type20_ix,bc_type20_iy)
	      bc_type20_ix=bc_type20_ixx
	      call xertst (-2.lt.bc_type20_ix.and.bc_type20_ix.lt.nx+1,
     1	       'stepped off x-grid in type 20 bc')
	      call xertst (-2.lt.bc_type20_iy.and.bc_type20_iy.lt.ny+1,
     1	       'stepped off y-grid in type 20 bc')
	    enddo
	    shi0(IX,BCPOS(IB),0)= VBIG*NE00*SY(IX,BCPOS(IB))*T1*
     1	     ti(bc_type20_ix,bc_type20_iy)/tiave
	    shi0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=20, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 40000

* -- BCENI=21 -- from b2stbc_spb                              !srv 01.02.09 {

32100   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'from b2stbc_spb for BCENI with enipar1 =',
     &     ENIPAR(IB,1),' and enipar2 =', ENIPAR(IB,2),' on ',
     &     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     &     ' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            s0 = vol(IX,BCPOS(IB))/hy1(IX,BCPOS(IB))
            if(mdf_fhi.ne.0) then                                    !srv 21.03.09 {
             shi0(IX,BCPOS(IB),0) = shi0(IX,BCPOS(IB),0) - 
     &          fhiPSch(ix,topiy(ix,bcpos(ib)),1)
	      endif                                                    !srv 21.03.09 }
            do is = 0, ns-1
	       if(.not.is_neutral(is)) then                            !srv 19.03.09 {
              t0 = na(ix,BCPOS(IB),is)
              csb = sqrt((rza(ix,BCPOS(IB),is)*te(ix,BCPOS(IB))+
     &                    ti(ix,BCPOS(IB)))/(am(is)*mp))
	        shi0(IX,BCPOS(IB),1)=shi0(IX,BCPOS(IB),1)+ENIPAR(IB,2)*
     &                            t0*csb*s0
!     &                            t0*csb*SY(IX,BCPOS(IB))
             endif
            enddo
            SHI0(IX,BCPOS(IB),1)= SHI0(IX,BCPOS(IB),1) - 
     &                 max(0.0e0_R8,-floi_noc(ix,topiy(ix,BCPOS(IB)),1)) !srv 28.10.02
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            s0 = vol(IX,BCPOS(IB))/hy1(IX,BCPOS(IB))
            if(mdf_fhi.ne.0) then                                    !srv 21.03.09 {
             shi0(IX,BCPOS(IB),0) = shi0(IX,BCPOS(IB),0) + 
     &          fhiPSch(ix,bcpos(ib),1)
	      endif                                                    !srv 21.03.09 }
            do is = 0, ns-1
	       if(.not.is_neutral(is)) then                            !srv 19.03.09 }
              t0 = na(ix,BCPOS(IB),is)
              csb = sqrt((rza(ix,BCPOS(IB),is)*te(ix,BCPOS(IB))+
     &                    ti(ix,BCPOS(IB)))/(am(is)*mp))
	        shi0(IX,BCPOS(IB),1)=shi0(IX,BCPOS(IB),1)+ENIPAR(IB,2)*
     &          t0*csb*s0
!     &          t0*csb*SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
             endif                                                       !srv 19.03.09 {
            enddo
            SHI0(IX,BCPOS(IB),1)=SHI0(IX,BCPOS(IB),1)-
     &           max(0.0e0_R8,floi_noc(ix,BCPOS(IB),1))              !srv 28.10.02
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=21, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 40000                                                   !srv 01.02.09 }

* -- BCENI=22 -- THE MODIFIED DECAY LENGTH FOR THE ION TEMPERATURE   !srv 23.07.10 {

32200   CONTINUE
        if(ncall.eq.0) then
          write(*,*) ' BCENI=22 with ENIPAR(IB,1) =',
     &     ENIPAR(IB,1),' on ', BCCHAR(IB),' at ',BCPOS(IB),
     &     ' from ',BCSTART(IB),' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
	      t0 = 0.0_R8
            do is = 0, ns-1
	       if(.not.is_neutral(is)) then                            !srv 19.03.09 {
              csb = sqrt((rza(ix,BCPOS(IB),is)*te(ix,BCPOS(IB))+
     &                    ti(ix,BCPOS(IB)))/(am(is)*mp))
	        t0 = t0 + csb*na(ix,BCPOS(IB),is)
             endif                                                   !srv 19.03.09 }
            enddo
            shi0(IX,BCPOS(IB),0)=VBIG*NE00*EV*
     &           CHCI(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)),1)*
     &             TI(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)))
            shi0(IX,BCPOS(IB),1)=- VBIG*NE00*EV*
     &           (CHCI(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)),1)-
     &            ENIPAR(IB,1)*SY(IX,BCPOS(IB))*T0)
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
	      t0 = 0.0_R8
            do is = 0, ns-1
	       if(.not.is_neutral(is)) then                            !srv 19.03.09 {
              csb = sqrt((rza(ix,BCPOS(IB),is)*te(ix,BCPOS(IB))+
     &                    ti(ix,BCPOS(IB)))/(am(is)*mp))
	        t0 = t0 + csb*na(ix,BCPOS(IB),is)
             endif                                                   !srv 19.03.09 }
            enddo
            shi0(IX,BCPOS(IB),0)=VBIG*NE00*EV*
     &           CHCI(IX,BCPOS(IB),1)*
     &             TI(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            shi0(IX,BCPOS(IB),1)=- VBIG*NE00*EV*
     &           (CHCI(IX,BCPOS(IB),1)-ENIPAR(IB,1)*                 !srv 02.02.11
     &            SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))*T0)
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCENI=22, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 40000                                                   !srv 23.07.10 }

* -- BCENI=23 -- PRESCRIBE THE POLOIDALLY AVERAGE VALUE OF THE ION TEMPERATURE
* AND INTRODUCE A POLOIDAL VARIATION AS CLOSE AS POSSIBLE TO NEOCLASSICAL SOLUTION.
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH CORRESPONDING
* CONDITION ON IONS DENSITY

32300   CONTINUE
        if(ncall.eq.0) then
#ifdef ASTRA
          IF(BCCHAR(IB).EQ.'S') THEN
	      ENIPAR(IB,1) = tib_astra*1.0e3_R8                           !srv 14.06.12
          ENDIF
#endif
          write(*,*) 'specified ion temperature',ENIPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2     BCEND(IB)
          call xertst(ENIPAR(IB,1).gt.0.0_R8,
     1     'BCENI = 23, ENIPAR(,1) !> 0')
        endif
        IF(BCCHAR(IB).EQ.'S' .and. 
     $   (mod(region(BCSTART(IB),bcpos(ib),0),4).eq.1)) THEN
!          DO IX=BCSTART(IB), BCEND(IB)
!            shi0(IX,BCPOS(IB),0)= VBIG*NE00*ENIPAR(IB,1)
!     .       *EV*SY(IX,BCPOS(IB))
!            shi0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
!          end do
	    t0=0.0e0_R8   ! Volume                                              !lk 12.12.12 { 
	    t1=0.0e0_R8   ! Poloidally volume averaged ions temperature in the cell next to boundary
	    t2=0.0e0_R8   ! Poloidally volume averaged electrons density in the cell next to boundary 
	    t3=0.0e0_R8   ! Poloidally volume averaged electrons density in the boundary cell
          do ix = -1,nx                                               
           if (mod(region(ix,bcpos(ib),0),4).eq.1) then                      
	      t0=t0 + vol(ix,0)                        
	      t1=t1 + ti(ix,0)*vol(ix,0)                        
	      t2=t2 + ne(ix,0)*vol(ix,0)                        
	      t3=t3 + ne(ix,-1)*vol(ix,0)                        
           endif                                                      
	    enddo
		t1=t1/t0                                                       
		t2=t2/t0                                                       
		t3=t3/t0                                                       
          do IX=BCSTART(IB), BCEND(IB)                                               
		 if(mod(region(ix,bcpos(ib),0),4).eq.1) then
            shi0(IX,BCPOS(IB),0)= shi0(IX,BCPOS(IB),0)+VBIG*NE00
     &      *(ENIPAR(IB,1)*EV+(ti(ix,0)-t1)*t1*t2/t3/ENIPAR(IB,1)/EV)
     &      *SY(IX,BCPOS(IB))
            shi0(IX,BCPOS(IB),1)=-VBIG*NE00*SY(IX,BCPOS(IB))
	     endif
	    enddo
	                                                     !lk 12.12.12 } 
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCENI=23 IS AVAILABLE FOR CORE BOUNDARY ONLY.')
        END IF
        GOTO 40000

* -- BCENI=24 -- Feedback boundary condition with prescribed total ions flux,
*                constant poloidally average ions temperature and a 
*        POLOIDAL VARIATION AS CLOSE AS POSSIBLE TO NEOCLASSICAL SOLUTION.
* IT IS RECOMMENDED TO USE THIS BOUNDARY CONDITION TOGETHER WITH CORRESPONDING
* CONDITION ON IONS DENSITY

32400   CONTINUE
!IYS 04.03.2013   {
          if(ncall.eq.0) then
	       write(*,*) 'Given integrated ion heat flux',
     &                  ' with constant Ti '
!          write all boundaries with this type of boundary condition
           DO IBW=IB,NBC
		if(BCCHAR(IB).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN	       
		   write(*,*) ' on ',BCCHAR(IBW),' at ',BCPOS(IBW),' from ',
     2          BCSTART(IBW),' to ', BCEND(IBW)
              call xertst(enipar(ibw,1).ne.0,'faulty parameter enipar1')
	      endif
	     ENDDO
#ifdef ASTRA
           SENIPAR = fhib_astra
#else
           SENIPAR = 0.0_R8
           DO IBW=IB,NBC
	      if(BCCHAR(IB).EQ.'S'.and.BCENI(IBW).EQ.BCENI(IB))THEN	       
               SENIPAR = SENIPAR + ENIPAR(IBW,1)
	      endif
	     ENDDO
#endif
		 write(*,*) ' Summary value is ', SENIPAR

	     TIS = TI((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB))
          endif
	    if(.not.once_bceni_16) then
           IF(BCCHAR(IB).EQ.'S' .and. 
     $   (mod(region(BCSTART(IB),bcpos(ib),0),4).eq.1)) THEN
	      totflux = 0.0_R8
	      t0=0.0e0_R8   ! Volume                                              !lk 12.12.12 { 
	      t1=0.0e0_R8   ! Poloidally volume averaged ions temperature in the cell next to boundary
              TIS = 0.0e0_R8 ! Poloidally volume averaged ions temperature in the boundary cell
	      t2=0.0e0_R8   ! Poloidally volume averaged electrons density in the cell next to boundary 
	      t3=0.0e0_R8   ! Poloidally volume averaged electrons density in the boundary cell
            DO IBW=IB,NBC
             IF(BCCHAR(IBW).EQ.'S'.AND.BCENI(IBW).EQ.BCENI(IB)) THEN
	        if(mdf_fhi.ne.0) then
               DO IX=BCSTART(IBW), BCEND(IBW)
	          totflux = totflux + 
     &          fhi_mdf(topix(ix,BCPOS(IBW)),topiy(ix,BCPOS(IBW)),1) 
     &       +fac_ExB(ix,BCPOS(IBW))*                                    !srv 11.04.13 {
     &      ((po(ix,BCPOS(IBW))-po(leftix(ix,BCPOS(IBW)),BCPOS(IBW)))/
     &       (hx(ix,BCPOS(IBW))+hx(leftix(ix,BCPOS(IBW)),BCPOS(IBW)))+
     &       (po(rightix(ix,BCPOS(IBW)),BCPOS(IBW))-po(ix,BCPOS(IBW)))/
     &       (hx(rightix(ix,BCPOS(IBW)),BCPOS(IBW))+hx(ix,BCPOS(IBW))))*
     &        ni(ix,BCPOS(IBW),1)*ti(ix,BCPOS(IBW))*SY(IX,BCPOS(IBW))*
     &        bb(ix,BCPOS(IBW),2)/bb(ix,BCPOS(IBW),3)**2                 !srv 11.04.13 }
	          t0=t0 + vol(ix,0)                        
	          t1=t1 + ti(ix,0)*vol(ix,0)                        
	          TIS=TIS + ti(ix,-1)*vol(ix,0)                        
	          t2=t2 + ne(ix,0)*vol(ix,0)                        
	          t3=t3 + ne(ix,-1)*vol(ix,0)                        
               ENDDO
	        else
               DO IX=BCSTART(IBW), BCEND(IBW)
	          totflux = totflux + 
     &          fhi(topix(ix,BCPOS(IBW)),topiy(ix,BCPOS(IBW)),1) 
	          t0=t0 + vol(ix,0)                        
	          t1=t1 + ti(ix,0)*vol(ix,0)                        
	          TIS=TIS + ti(ix,-1)*vol(ix,0)                        
	          t2=t2 + ne(ix,0)*vol(ix,0)                        
	          t3=t3 + ne(ix,-1)*vol(ix,0)                        
               ENDDO
	        endif
	       ENDIF
	      ENDDO
!            TIS = TI((BCEND(IB)+BCSTART(IB))/2,BCPOS(IB)) * 
!     &         (1.0_R8-bc_ref*(totflux-SENIPAR)/
!     &           max(abs(totflux),SENIPAR))
		t1=t1/t0                                                       
		t2=t2/t0                                                       
		t3=t3/t0                                                       
                TIS=TIS/t0* 
     &         (1.0_R8-bc_ref_ti*(totflux-SENIPAR)/
     &           max(abs(totflux),SENIPAR))
             write(*,"(' BCENI=24: TIS=',es14.7,' totflux=',es14.7,
     *       ' SENIPAR=',es14.7)") TIS/ev, totflux, SENIPAR
            DO IBW=IB,NBC
             IF(BCCHAR(IBW).EQ.'S'.AND.BCENE(IBW).EQ.BCENE(IB)) THEN
              DO IX=BCSTART(IBW), BCEND(IBW)
               shi0(IX,BCPOS(IBW),0)= VBIG*NE00*
     $           (TIS+(ti(ix,0)-t1)*t1*t2/t3/TIS)*SY(IX,BCPOS(IBW))
               shi0(IX,BCPOS(IBW),1)=-VBIG*NE00*SY(IX,BCPOS(IBW))
              ENDDO
	       ENDIF
	      ENDDO
		  once_bceni_16 = .true.
           ELSE
             CALL XERRAB
     1         ('B2STBC -- BCENI=24 IS AVAILABLE FOR CORE ONLY')
           ENDIF
	    endif                                                          !srv 16.02.12 }
          GOTO 40000

* -- BOUNDARY CONDITIONS FOR POTENTIAL EQUATION --

40000   CONTINUE
        if(BCPOT(IB).eq.0) goto 50000
	GOTO(40100,40200,40300,40400,40500,
     1	 40600,40700,40800,40900,41000,41100),BCPOT(IB)
        write(*,*) BCPOT(IB)
        call xerrab ('BCPOT out of bounds')

* -- BCPOT=1 -- PRESCRIBE THE VALUE OF THE POTENTIAL

40100   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'specified potential',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2     ' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            sch0(IX,BCPOS(IB),0)= VBIG*POTPAR(IB,1)*SY(IX,BCPOS(IB))
            sch0(IX,BCPOS(IB),1)=-VBIG*SY(IX,BCPOS(IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            sch0(IX,BCPOS(IB),0)= VBIG*POTPAR(IB,1)*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            sch0(IX,BCPOS(IB),1)=-VBIG*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            sch0(BCPOS(IB),IY,0)= UBIG*POTPAR(IB,1)*SX(BCPOS(IB),IY)
            sch0(BCPOS(IB),IY,1)=-UBIG*SX(BCPOS(IB),IY)
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            sch0(BCPOS(IB),IY,0)= UBIG*POTPAR(IB,1)*
     1       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            sch0(BCPOS(IB),IY,1)=-UBIG*
     1       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          end do
        ELSEIF(BCCHAR(IB).EQ.'X') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            sch0(IX,BCPOS(IB),0)= UBIG*
     .       POTPAR(IB,1)*VOL(IX,BCPOS(IB))
            sch0(IX,BCPOS(IB),1)=-UBIG*
     .       VOL(IX,BCPOS(IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'Y') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            sch0(BCPOS(IB),IY,0)= UBIG*
     .       POTPAR(IB,1)*VOL(BCPOS(IB),IY)
            sch0(BCPOS(IB),IY,1)=-UBIG*
     .       VOL(BCPOS(IB),IY)
          end do
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCPOT=1, INVALID OPTION.')
        ENDIF
        GOTO 50000

* -- BCPOT=2 -- PRESCRIBE THE GRADIENT OF THE POTENTIAL

40200   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'potential gradient',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2     ' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            sch0(IX,BCPOS(IB),0)= VBIG*SY(IX,BCPOS(IB))*
     .       (-POTPAR(IB,1)*0.5_R8*
     .       (HY1(IX,BCPOS(IB))+
     .        HY1(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB)))) +
     .       PO(TOPIX(IX,BCPOS(IB)),TOPIY(IX,BCPOS(IB))))
            sch0(IX,BCPOS(IB),1)=-VBIG*SY(IX,BCPOS(IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            sch0(IX,BCPOS(IB),0)= VBIG*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))*
     .       (POTPAR(IB,1)*0.5_R8*
     .       (HY1(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))+
     .        HY1(IX,BCPOS(IB))) +
     .       PO(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB))))
            sch0(IX,BCPOS(IB),1)=-VBIG*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            sch0(BCPOS(IB),IY,0)= UBIG*SX(BCPOS(IB),IY)*
     .       (-POTPAR(IB,1)*0.5_R8*
     .       (HX(BCPOS(IB),IY)+
     .        HX(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY))) +
     .        PO(RIGHTIX(BCPOS(IB),IY),RIGHTIY(BCPOS(IB),IY)))
            sch0(BCPOS(IB),IY,1)=-UBIG*SX(BCPOS(IB),IY)
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            sch0(BCPOS(IB),IY,0)= UBIG*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))*
     .       (POTPAR(IB,1)*0.5_R8*
     .       (HX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))+
     .        HX(BCPOS(IB),IY)) +
     .       PO(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY)))
            sch0(BCPOS(IB),IY,1)=-UBIG*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          end do
        ELSE
          CALL XERRAB
     1     ('B2STBC -- BCPOT=2, INVALID OPTION.')
        END IF
        GOTO 50000

* -- BCPOT=3 -- SHEATH CONDITIONS, ELECTRON ENERGY TRANSMISSION

40300   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'sheath potential bc on ',BCCHAR(IB),
     1     ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            s1 = abs(pbs(rightix(bcpos(ib),iy),rightiy(bcpos(ib),iy),0))
            vte = sqrt(bcint(te,nx,ny,iy,ib)/me)
            t0 = 0
            do is = 0, ns-1
              t0 = t0-bcint(rza(-1,-1,is),nx,ny,iy,ib)*
     1         fna(rightix(bcpos(ib),iy),rightiy(bcpos(ib),iy),0,is)
            enddo
            fchi = qe*t0
            if(t0.lt.0 .and. bottomiy(bcpos(ib),iy).ne.-2 .and.
     &       topiy(bcpos(ib),iy).ne.ny+1) then
              wrong_flow = .true.                                     !xpb
              if(diagno.ge.0) then
                write(*,*) 'b2stbc_west:bcpos,iy,fchi',
     1           bcpos(ib),iy,fchi
                write(*,*) 
     .         rza(rightix(bcpos(ib),iy),rightiy(bcpos(ib),iy),ismain),
     &         fna(rightix(bcpos(ib),iy),rightiy(bcpos(ib),iy),0,ismain)
                write(*,*) 'Bx,Bz,Uadia,-Ua =',
     &           wbbr(bcpos(ib),iy,0),wbbr(bcpos(ib),iy,2),
     &           uadia(bcpos(ib),iy,0,ismain),-ua(bcpos(ib),iy,ismain)
               endif
c              fchi=0.0_R8
            endif
            fche = sqrt(1/(2*pi))*qe*bcint(ne,nx,ny,iy,ib)*vte*s1*
     &       expu2(-qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
     .       bcint(te,nx,ny,iy,ib))
            if (secmodel.eq.1) then
              if (ne(bcpos(ib),iy)/(rza(bcpos(ib),iy,ismain)*
     .            na(bcpos(ib),iy,ismain)).gt.1.1_R8) then
                write(*,*) 'Large ne contribution from impurities!'
                write(*,*) 'Kinetic s.e.c. model may not be valid!'
              endif
              call eseec0 (bcint(te,nx,ny,iy,ib),
     &         bcint(ti,nx,ny,iy,ib),bcint(ne,nx,ny,iy,ib),
     &         bcint(po,nx,ny,iy,ib),wbbr(bcpos(ib),iy,0),
     &         wbbr(bcpos(ib),iy,3),0,bcpos(ib),iy,
     &         nx,ny,ns,ismain,seec,diagno)                           !cfs
            else
              seec = gammae
            endif
            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/bcint(te,nx,ny,iy,ib)
*    ..compute charge source
            sch0(bcpos(ib),iy,0) = sch0(bcpos(ib),iy,0)+
     &       (1.0_R8-seec)*fche-fchi+t0*po(bcpos(ib),iy)
            sch0(bcpos(ib),iy,1) = sch0(bcpos(ib),iy,1)-t0
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            s1 = abs(pbs(bcpos(ib),iy,0))
            vte = sqrt(bcint(te,nx,ny,iy,ib)/me)
            t0 = 0
            do is = 0, ns-1
              t0 = t0+
     &         bcint(rza(-1,-1,is),nx,ny,iy,ib)*fna(bcpos(ib),iy,0,is)
            enddo
            fchi = qe*t0
            if(t0.lt.0 .and. bottomiy(bcpos(ib),iy).ne.-2 .and.
     &       topiy(bcpos(ib),iy).ne.ny+1) then
              wrong_flow = .true.                                     !xpb
              if(diagno.ge.0) then
                write(*,*) 'b2stbc_east:bcpos,iy,fchi',
     1           bcpos(ib),iy,fchi
                write(*,*) 
     &           rza(bcpos(ib),iy,ismain),fna(bcpos(ib),iy,0,ismain)
                write(*,*) 'Bx,Bz,Uadia,Ua =',
     &           wbbl(bcpos(ib),iy,0),wbbl(bcpos(ib),iy,2),
     &           uadia(bcpos(ib),iy,0,ismain),ua(bcpos(ib),iy,ismain)
               endif
c              fchi=0.0_R8
            endif
            fche = sqrt(1/(2*pi))*qe*bcint(ne,nx,ny,iy,ib)*vte*s1*
     &       expu2(-qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
     .       bcint(te,nx,ny,iy,ib))
            if (secmodel.eq.1) then
              if (ne(bcpos(ib),iy)/(rza(bcpos(ib),iy,ismain)*
     .            na(bcpos(ib),iy,ismain)).gt.1.1_R8) then
                write(*,*) 'Large ne contribution from impurities!'
                write(*,*) 'Kinetic s.e.c. model may not be valid!'
              endif
              call eseec0 (bcint(te,nx,ny,iy,ib),
     &         bcint(ti,nx,ny,iy,ib),bcint(ne,nx,ny,iy,ib),
     &         bcint(po,nx,ny,iy,ib),wbbl(bcpos(ib),iy,0),
     &         wbbl(bcpos(ib),iy,3),0,bcpos(ib),iy,
     &         nx,ny,ns,ismain,seec,diagno)                           !cfs
            else
              seec = gammae
            endif
            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/bcint(te,nx,ny,iy,ib)
*    ..compute charge source
            sch0(bcpos(ib),iy,0) = sch0(bcpos(ib),iy,0)+
     &       (1.0_R8-seec)*fche-fchi+t0*po(bcpos(ib),iy)
            sch0(bcpos(ib),iy,1) = sch0(bcpos(ib),iy,1)-t0
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=3, NO SHEATH IN RADIAL DIR.')
        END IF
        GOTO 50000

* -- BCPOT=4 -- PRESCRIBE THE VALUE OF THE POTENTIAL,
*                 WEAKLY A MIXED BOUNDARY CONDITION

40400   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'weak potential',POTPAR(IB,1),' on ',BCCHAR(IB),
     1     ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
          call xertst(POTPAR(IB,2).gt.0.0_R8,
     1     'BCPOT = 4, POTPAR(,2) !> 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            sch0(IX,BCPOS(IB),0)= POTPAR(IB,2)*POTPAR(IB,1)*
     1       SY(IX,BCPOS(IB))
            sch0(IX,BCPOS(IB),1)=-POTPAR(IB,2)*SY(IX,BCPOS(IB))
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            sch0(IX,BCPOS(IB),0)= POTPAR(IB,2)*POTPAR(IB,1)*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            sch0(IX,BCPOS(IB),1)=-POTPAR(IB,2)*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            sch0(BCPOS(IB),IY,0)= POTPAR(IB,2)*POTPAR(IB,1)*
     1       SX(BCPOS(IB),IY)
            sch0(BCPOS(IB),IY,1)=-POTPAR(IB,2)*SX(BCPOS(IB),IY)
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            sch0(BCPOS(IB),IY,0)= POTPAR(IB,2)*
     1       POTPAR(IB,1)*SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            sch0(BCPOS(IB),IY,1)=-POTPAR(IB,2)*
     1       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=4, NOT YET AVAILABLE')
        END IF
        GOTO 50000

* -- BCPOT=5 -- PRESCRIBE THE CURRENT FLUX DENSITY PER UNIT AREA

40500   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'current flux density per area',POTPAR(IB,1),
     1     ' on ',BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),
     2     ' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            sch0(IX,BCPOS(IB),0)=POTPAR(IB,1)*SY(IX,BCPOS(IB))
            sch0(IX,BCPOS(IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            sch0(IX,BCPOS(IB),0)=POTPAR(IB,1)*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            sch0(IX,BCPOS(IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            sch0(BCPOS(IB),IY,0)=POTPAR(IB,1)*SX(BCPOS(IB),IY)
            sch0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            sch0(BCPOS(IB),IY,0)=POTPAR(IB,1)*
     1       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            sch0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=5, NOT YET AVAILABLE')
        END IF
        GOTO 50000

* -- BCPOT=6 -- PRESCRIBE THE TOTAL CURRENT FLUX DENSITY FOR A
*                 CONSTANT POTENTIAL

40600   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'total current flux density',POTPAR(IB,1),
     1     ' with const potential on ',BCCHAR(IB),' at ',BCPOS(IB),
     2     ' from ',BCSTART(IB),' to ',BCEND(IB)
        endif
        CALL XERRAB ('B2STBC -- BCPOT=6, NOT YET AVAILABLE')
        GOTO 50000

* -- BCPOT=7 -- PRESCRIBE THE POTENTIAL AS A FUNCTION OF
*                 OTHER PLASMA PARAMETERS

40700   CONTINUE
        CALL XERRAB ('B2STBC -- BCPOT=7, NOT YET AVAILABLE')
        GOTO 50000

* -- BCPOT=8 -- PRESCRIBE THE TOTAL ELECTRON CURRENT WITH CONSTANT
*                 FLUX DENSITY

40800   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'total current',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',
     2     BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          US=0.0_R8
          DO IX=BCSTART(IB), BCEND(IB)
            US=US+SY(IX,BCPOS(IB))
          end do
          DO IX=BCSTART(IB), BCEND(IB)
            sch0(IX,BCPOS(IB),0)=POTPAR(IB,1)/US*SY(IX,BCPOS(IB))
            sch0(IX,BCPOS(IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          US=0.0_R8
          DO IX=BCSTART(IB), BCEND(IB)
            US=US+SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
          end do
          DO IX=BCSTART(IB), BCEND(IB)
            sch0(IX,BCPOS(IB),0)=POTPAR(IB,1)/US*
     .       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))
            sch0(IX,BCPOS(IB),1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          US=0.0_R8
          DO IY=BCSTART(IB), BCEND(IB)
            US=US+SX(BCPOS(IB),IY)
          end do
          DO IY=BCSTART(IB), BCEND(IB)
            sch0(BCPOS(IB),IY,0)=POTPAR(IB,1)/US*SX(BCPOS(IB),IY)
            sch0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          US=0.0_R8
          DO IY=BCSTART(IB), BCEND(IB)
            US=US+SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
          end do
          DO IY=BCSTART(IB), BCEND(IB)
            sch0(BCPOS(IB),IY,0)=POTPAR(IB,1)/US*
     .       SX(LEFTIX(BCPOS(IB),IY),LEFTIY(BCPOS(IB),IY))
            sch0(BCPOS(IB),IY,1)=0.0_R8
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=8, NOT AVAILABLE')
        ENDIF
        GOTO 50000

* -- BCPOT=9 -- PRESCRIBE THE DECAY LENGTH FOR
*                 THE POTENTIAL

40900   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'potential decay length',POTPAR(IB,1),' on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),
     2     ' from ',BCSTART(IB),' to ',BCEND(IB)
          call xertst(POTPAR(IB,1).ne.0.0_R8,
     1     'BCPOT = 9, POTPAR(,1) = 0')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            sch0(IX,BCPOS(IB),0)=0.0_R8
            sch0(IX,BCPOS(IB),1)=-sig0(ix,BCPOS(IB))/POTPAR(IB,1)*
     1       SY(IX,BCPOS(IB))/qe/ne(ix,bcpos(ib))
          ENDDO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            sch0(IX,BCPOS(IB),0)=0.0_R8
            sch0(IX,BCPOS(IB),1)=-sig0(ix,BCPOS(IB))/POTPAR(IB,1)*
     1       SY(BOTTOMIX(IX,BCPOS(IB)),BOTTOMIY(IX,BCPOS(IB)))/qe/
     2       ne(ix,bcpos(ib))
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=9, NOT AVAILABLE FOR W,E')
        ENDIF
        GOTO 50000

* -- BCPOT=10 -- Feedback option for core
*                NOT YET TESTED!!!!!!!!!

41000   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'current feedback bc on ',
     1     BCCHAR(IB),' at ',BCPOS(IB),
     2     ' from ',BCSTART(IB),' to ',BCEND(IB)
          write(*,*) 'NOT YET TESTED !!!!!!!!!!'
          call xertst(lfeedback, 'LFEEDBACK must be turned on!')
        endif
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB), BCEND(IB)
            if (region(ix,bcpos(ib),0).eq.1) then
              sch0(IX,BCPOS(IB),0)=cbsch(0,coreregno)*SY(IX,BCPOS(IB))
            elseif (region(ix,bcpos(ib),0).eq.5) then
              sch0(IX,BCPOS(IB),0)=cbsch(0,coreregn2)*SY(IX,BCPOS(IB))
            endif
            sch0(IX,BCPOS(IB),1)=0.0_R8
          ENDDO
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=10, NOT AVAILABLE FOR N,W,E')
        ENDIF
        GOTO 50000

* -- BCPOT=11 -- SHEATH CONDITIONS, ELECTRON ENERGY TRANSMISSION FROM b2stbc_spb !srv 01.02.09 {

41100   CONTINUE
        if(ncall.eq.0) then
          write(*,*) 'sheath potential bc from b2stbc_spb',
     &     ' with POTPAR(IB,2) =',POTPAR(IB,2),' gammae =', gammae,
     &     'on ',BCCHAR(IB),
     1     ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
        endif
        IF(BCCHAR(IB).EQ.'W') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            s1 = abs(pbs(rightix(bcpos(ib),iy),rightiy(bcpos(ib),iy),0))
!            vte = sqrt(bcint(te,nx,ny,iy,ib)/me)
            vte = sqrt(te(bcpos(ib),iy)/me)
            wrong_flow=.false.
*    ..set boundary condition in explicit form (the same form as fna)
            CS=SQRT(PR(BCPOS(IB),IY)/RO(BCPOS(IB),IY))
            t0 = 0.0e0_R8
            do is = 0, ns-1
	       if(istyle_fchi.eq.1 .and. BCCON(IS,IB).eq.14) then              !srv 14.02.05 {
              t0 = t0-rza(bcpos(ib),iy,is)*
     *             CONPAR(IS,IB,1)*cs*s1*na(bcpos(ib),iy,is)
	       else 
              t0 = t0-rza(bcpos(ib),iy,is)*
     &             fna(rightix(bcpos(ib),iy),rightiy(bcpos(ib),iy),0,is)     !srv 26.11.02
	       endif                                                           !srv 14.02.05 }
            enddo
            fchi = qe*t0
!            t0 = 0
!            do is = 0, ns-1
!              t0 = t0-bcint(rza(-1,-1,is),nx,ny,iy,ib)*
!     1         fna(rightix(bcpos(ib),iy),rightiy(bcpos(ib),iy),0,is)
!            enddo
!            fchi = qe*t0
!            if(t0.lt.0 .and. bottomiy(bcpos(ib),iy).ne.-2 .and.
!     &       topiy(bcpos(ib),iy).ne.ny+1) then
!              wrong_flow = .true.                                     !xpb
!              if(diagno.ge.0) then
!                write(*,*) 'b2stbc_west:bcpos,iy,fchi',
!     1           bcpos(ib),iy,fchi
!                write(*,*) 
!     .         rza(rightix(bcpos(ib),iy),rightiy(bcpos(ib),iy),ismain),
!     &         fna(rightix(bcpos(ib),iy),rightiy(bcpos(ib),iy),0,ismain)
!                write(*,*) 'Bx,Bz,Uadia,-Ua =',
!     &           wbbr(bcpos(ib),iy,0),wbbr(bcpos(ib),iy,2),
!     &           uadia(bcpos(ib),iy,0,ismain),-ua(bcpos(ib),iy,ismain)
!               endif
c              fchi=0.0_R8
!            endif
!            fche = sqrt(1/(2*pi))*qe*bcint(ne,nx,ny,iy,ib)*vte*s1*
!     &       expu2(-qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
!     .       bcint(te,nx,ny,iy,ib))
            fche = sqrt(1/(2*pi))*qe*ne(bcpos(ib),iy)*vte*s1*
     &       expu2(-qe*(po(bcpos(ib),iy)-potpar(ib,2))/te(bcpos(ib),iy))
            if (secmodel.eq.1) then
              if (ne(bcpos(ib),iy)/(rza(bcpos(ib),iy,ismain)*
     .            na(bcpos(ib),iy,ismain)).gt.1.1_R8) then
                write(*,*) 'Large ne contribution from impurities!'
                write(*,*) 'Kinetic s.e.c. model may not be valid!'
              endif
              call eseec0 (bcint(te,nx,ny,iy,ib),
     &         bcint(ti,nx,ny,iy,ib),bcint(ne,nx,ny,iy,ib),
     &         bcint(po,nx,ny,iy,ib),wbbr(bcpos(ib),iy,0),
     &         wbbr(bcpos(ib),iy,3),0,bcpos(ib),iy,
     &         nx,ny,ns,ismain,seec,diagno)                           !cfs
            else
              seec = gammae
            endif
!            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/bcint(te,nx,ny,iy,ib)
            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/te(bcpos(ib),iy)
*    ..compute charge source
            sch0(bcpos(ib),iy,0) = sch0(bcpos(ib),iy,0)+
     &       (1.0_R8-seec)*fche-fchi+t0*po(bcpos(ib),iy)
            sch0(bcpos(ib),iy,1) = sch0(bcpos(ib),iy,1)-t0
          end do
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=BCSTART(IB), BCEND(IB)
            s1 = abs(pbs(bcpos(ib),iy,0))
!            vte = sqrt(bcint(te,nx,ny,iy,ib)/me)
!            t0 = 0
!            do is = 0, ns-1
!              t0 = t0+
!     &         bcint(rza(-1,-1,is),nx,ny,iy,ib)*fna(bcpos(ib),iy,0,is)
!            enddo
!            fchi = qe*t0
!            if(t0.lt.0 .and. bottomiy(bcpos(ib),iy).ne.-2 .and.
!     &       topiy(bcpos(ib),iy).ne.ny+1) then
!              wrong_flow = .true.                                     !xpb
!              if(diagno.ge.0) then
!                write(*,*) 'b2stbc_east:bcpos,iy,fchi',
!     1           bcpos(ib),iy,fchi
!                write(*,*) 
!     &           rza(bcpos(ib),iy,ismain),fna(bcpos(ib),iy,0,ismain)
!                write(*,*) 'Bx,Bz,Uadia,Ua =',
!     &           wbbl(bcpos(ib),iy,0),wbbl(bcpos(ib),iy,2),
!     &           uadia(bcpos(ib),iy,0,ismain),ua(bcpos(ib),iy,ismain)
!               endif
c              fchi=0.0_R8
!            endif
!            fche = sqrt(1/(2*pi))*qe*bcint(ne,nx,ny,iy,ib)*vte*s1*
!     &       expu2(-qe*(bcint(po,nx,ny,iy,ib)-potpar(ib,2))/
!     .       bcint(te,nx,ny,iy,ib))
            vte = sqrt(te(bcpos(ib),iy)/me)
            wrong_flow=.false.
            CS=SQRT(PR(BCPOS(IB),IY)/RO(BCPOS(IB),IY))
            t0 = 0.0e0_R8
            do is = 0, ns-1
	       if(istyle_fchi.eq.1 .and. BCCON(IS,IB).eq.14) then              !srv 14.02.05 {
*    ..set boundary condition in explicit form (the same form as fna)
              t0 = t0-rza(bcpos(ib),iy,is)*
     &                CONPAR(IS,IB,1)*cs*s1*na(bcpos(ib),iy,is)
	       else                                                            !srv 14.02.05 {
              t0 = t0+rza(bcpos(ib),iy,is)*fna(bcpos(ib),iy,0,is)            !srv 26.11.02
	       endif                                                           !srv 14.02.05 }
            enddo
            fchi = qe*t0
            fche = sqrt(1/(2*pi))*qe*ne(bcpos(ib),iy)*vte*s1*
     &       expu2(-qe*(po(bcpos(ib),iy)-potpar(ib,2))/te(bcpos(ib),iy))
            if (secmodel.eq.1) then
              if (ne(bcpos(ib),iy)/(rza(bcpos(ib),iy,ismain)*
     .            na(bcpos(ib),iy,ismain)).gt.1.1_R8) then
                write(*,*) 'Large ne contribution from impurities!'
                write(*,*) 'Kinetic s.e.c. model may not be valid!'
              endif
              call eseec0 (bcint(te,nx,ny,iy,ib),
     &         bcint(ti,nx,ny,iy,ib),bcint(ne,nx,ny,iy,ib),
     &         bcint(po,nx,ny,iy,ib),wbbl(bcpos(ib),iy,0),
     &         wbbl(bcpos(ib),iy,3),0,bcpos(ib),iy,
     &         nx,ny,ns,ismain,seec,diagno)                           !cfs
            else
              seec = gammae
            endif
!            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/bcint(te,nx,ny,iy,ib)
            t0 = max(fchi,(1.0_R8-seec)*fche)*qe/te(bcpos(ib),iy)
*    ..compute charge source
            sch0(bcpos(ib),iy,0) = sch0(bcpos(ib),iy,0)+
     &       (1.0_R8-seec)*fche-fchi+t0*po(bcpos(ib),iy)
            sch0(bcpos(ib),iy,1) = sch0(bcpos(ib),iy,1)-t0
          end do
        ELSE
          CALL XERRAB ('B2STBC -- BCPOT=3, NO SHEATH IN RADIAL DIR.')
        END IF
        GOTO 50000                                                    !srv 01.02.09 }

50000	continue

        if(istyle_cur_contr_on_S_and_N.eq.0) then                       !srv 01.02.09
!xpb
!xpb  Adding fchy_dia contributions on south and north boundaries
!xpb  following "sv 01.07.99" coding (see B2.5-style implementation)
!xpb  Special attention is paid to the case where there is a discontinuity in B
!xpb
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=BCSTART(IB),BCEND(IB)
            if(maxval(facdrift) .ne. 0.0_R8) then
              if(mod(region(ix,bcpos(ib),0),4).eq.1) then
                if(periodic_bc.ne.1 .or. rightix(ix,bcpos(ib)).eq.ix+1
     >           .or. nncut.eq.0) then
              fchy_dia = -praverage(1)*facdrift(ix,bcpos(ib))
     &         *b2stbc_fchy_dia                                   !srv 01.07.99
     &         *(OnedBsqc(rightix(topix(ix,bcpos(ib)),
     1                            topiy(ix,bcpos(ib))),
     &                    rightiy(topix(ix,bcpos(ib)),
     1                            topiy(ix,bcpos(ib))))
     &          -OnedBsqc(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))))      !srv
                else
      if(ncall.eq.0)
     > write(*,*) 'b2stbc: special treatment for cell ix,iy = ',
     . ix,bcpos(ib)
              fchy_dia = -praverage(1)*facdrift(ix,bcpos(ib))
     &         *b2stbc_fchy_dia  !srv 01.07.99
     &           *(1.0_R8/
     1             wbbv(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),3)**2
     &             -OnedBsqc(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))))   !xpb
                endif
              else
                if(periodic_bc.ne.1 .or. nncut.eq.0 .or.
     &     rightix(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))).eq.ix+1)
     >           then
                  fchy_dia = -bcint(wdia,nx,ny,ix,ib)
     .             *(vol(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))+
     .               vol(ix,bcpos(ib)))
     .             /(gs(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),2)+
     .               gs(ix,bcpos(ib),2))
     .      *(OnedBsqc(rightix(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))),
     .                 rightiy(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))))
     .             -OnedBsqc(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))))
     .             *facdrift(ix,bcpos(ib))*b2stbc_fchy_dia
                else
      if (ncall.eq.0)
     > write(*,*) 'b2stbc: special treatment for cell ix,iy = ',
     . ix,bcpos(ib)
                 fchy_dia = -bcint(wdia,nx,ny,ix,ib)
     .             *(vol(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)))+
     .               vol(ix,bcpos(ib)))
     .             /(gs(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),2)+
     .               gs(ix,bcpos(ib),2))
     .     *(1.0_R8/wbbv(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),3)**2
     .     -OnedBsqc(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))))
     .             *facdrift(ix,bcpos(ib))*b2stbc_fchy_dia
                endif
              endif
            else
              if(b2stbc_neoclassical .ne. 0.0_R8 .and.
     &           mod(region(ix,bcpos(ib),0),4).eq.1) then
                if(periodic_bc.ne.1 .or. nncut.eq.0 .or.
     &     rightix(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))).eq.ix+1)
     >           then
                  fchy_dia = -praverage(1)
     .    *(OnedBsqc(rightix(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))),
     .               rightiy(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))))
     .              -OnedBsqc(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))))
     .              *b2stbc_neoclassical
                else
                  fchy_dia = -praverage(1)
     .     *(1.0_R8/wbbv(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib)),3)**2
     .              -OnedBsqc(topix(ix,bcpos(ib)),topiy(ix,bcpos(ib))))
     .              *b2stbc_neoclassical
                endif
              else
                fchy_dia = 0.0_R8
              endif
            endif
            if(mod(region(ix,bcpos(ib),0),4).ne.1.and.
     1       b2stbc_fchy_dia_coreonly.eq.1) fchy_dia = 0.0_R8
            sch0(ix,bcpos(ib),0) = sch0(ix,bcpos(ib),0) + fchy_dia
          END DO
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=BCSTART(IB),BCEND(IB)
            if(maxval(facdrift) .ne. 0.0_R8) then
              if(region(ix,bcpos(ib),0).eq.5) then
            fchy_dia = -praverage(2)*facdrift(ix,bcpos(ib))*
     &       b2stbc_fchy_dia                                         !srv 01.07.99
     &       *(OnedBsqc(rightix(ix,bcpos(ib)),rightiy(ix,bcpos(ib)))-
     1         OnedBsqc(ix,bcpos(ib)))                               !srv
              else
                fchy_dia = -bcint(wdia,nx,ny,ix,ib)
     .           *(vol(ix,bcpos(ib))+
     .             vol(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib))))
     .           /(gs(ix,bcpos(ib),2)+
     .             gs(bottomix(ix,bcpos(ib)),bottomiy(ix,bcpos(ib)),2))
     .           *(OnedBsqc(rightix(ix,bcpos(ib)),rightiy(ix,bcpos(ib)))
     .            -OnedBsqc(ix,bcpos(ib)))*facdrift(ix,bcpos(ib))
     &            *b2stbc_fchy_dia
              endif
            else
              if(b2stbc_neoclassical .ne. 0.0_R8 .and.
     &           region(ix,bcpos(ib),0).eq.5) then
                fchy_dia = -praverage(2)
     .     *(OnedBsqc(rightix(ix,bcpos(ib)),rightiy(ix,bcpos(ib)))
     .      -OnedBsqc(ix,bcpos(ib)))
     .     *b2stbc_neoclassical
              else
                fchy_dia = 0.0_R8
              endif
            endif
            if(mod(region(ix,bcpos(ib),0),4).ne.1.and.
     1       b2stbc_fchy_dia_coreonly.eq.1) fchy_dia = 0.0_R8
            sch0(ix,bcpos(ib),0) = sch0(ix,bcpos(ib),0) - fchy_dia
          END DO
        ENDIF
!xpb
        else                                                         !srv 01.02.09 {
         IF(BCCHAR(IB).EQ.'S') THEN
           DO IX=BCSTART(IB),BCEND(IB)
             if(mod(region(ix,bcpos(ib),0),4).eq.1) then
*        for inner surfaces
              if(maxval(facdrift) .ne. 0.0_R8) then
                fchy_dia =                                             !srv 01.07.99
     &          -pzaverage*facdrift(ix,bcpos(ib))*phm0       ! *signmf        !srv 24.03.99 01.07.99 21.04.00 18.01.02 2
     &          *bzhzaverage                                           !lk 20.11.07 }
     &          *(OnedBsqc(rightix(ix,topiy(ix,bcpos(ib))),
     &                     rightiy(ix,topiy(ix,bcpos(ib))))-
     &            OnedBsqc(ix,topiy(ix,bcpos(ib))))	                       !srv 01.07.99
                sch0(ix,bcpos(ib),0) = sch0(ix,bcpos(ib),0) + fchy_dia

	          fchy_inertia = 0.0_R8                                         !srv 25.06.03
	          do is = 0, ns-1                                                 !srv 25.06.03
	           if(.not.is_neutral(is))then                                  !srv 25.06.03
                  fchy_inertia = fchy_inertia                                   !srv 07.04.00 25.06.03
     &             -0.5e0_R8*am(is)*mp*facdrift(ix,bcpos(ib))*phm1  ! *signmf                 !srv 07.04.00 21.04.00 21.10.02
     &             *naaverage(is)                                              !lk 20.11.07                       !srv 07.04.00
     &             *ua2average(is)                                            !srv 07.04.00 15.02.02
     &             *bzhzaverage                                               !lk 20.11.07
     &             *(OnedBsqc(rightix(ix,topiy(ix,bcpos(ib))),
     &                        rightiy(ix,topiy(ix,bcpos(ib))))-
     &               OnedBsqc(ix,topiy(ix,bcpos(ib))))                        !srv 07.04.00
	           endif
	          enddo
                sch0(ix,bcpos(ib),0) = sch0(ix,bcpos(ib),0)+fchy_inertia
               endif
             else
*        for PR region
               sch0(ix,bcpos(ib),0) = sch0(ix,bcpos(ib),0) + 
     &             fchdia(ix,topiy(ix,bcpos(ib)),1)+
     &             fchin(ix,topiy(ix,bcpos(ib)),1)+
     &             fchvispar(ix,topiy(ix,bcpos(ib)),1)+
     &             fchvisq(ix,topiy(ix,bcpos(ib)),1)+
     &             fchvisper(ix,topiy(ix,bcpos(ib)),1)+
     &             fchinert(ix,topiy(ix,bcpos(ib)),1)
	       endif
           END DO
         ELSEIF(BCCHAR(IB).EQ.'N') THEN
           DO IX=BCSTART(IB),BCEND(IB)
             sch0(ix,bcpos(ib),0) = sch0(ix,bcpos(ib),0) - 
     &             fchdia(ix,bcpos(ib),1)-
     &             fchin(ix,bcpos(ib),1)-
     &             fchvispar(ix,bcpos(ib),1)-
     &             fchvisq(ix,bcpos(ib),1)-
     &             fchvisper(ix,bcpos(ib),1)-
     &             fchinert(ix,bcpos(ib),1)
           END DO
         ENDIF
	  endif                                                        !srv 01.02.09 } 

      enddo

* ---------------------------------------------------------------------
*  INTERIOR BOUNDARY AT DEAD REGIONS
*  ---------------------------------------------------------------------
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            DO IS=0, NS-1
              SNA0(IX,IY,0,IS)= UBIG*NIISO(IS)
              SNA0(IX,IY,1,IS)=-UBIG
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

	if (iout.eq.1 .or. iout_b2npmo.eq.1) then                       !srv 11.09.09 {
       do is = 0, ns-1
	  write (chns,'(i2.2)') is
!	  do k=0,3
!	   write (chk,'(i1)') k
!	   call my_out(70,nx,ny,smo0(-1,-1,k,is),
!     &               'b2stbc_phys_smo'//chk//chns)
!	  enddo
	 call my_out(70,nx,ny,smo0(:,:,0,is)+
     &        smo0(:,:,1,is)*ua(:,:,is)+
     &        smo0(:,:,2,is)*na(:,:,is)*mp*am(is)+
     &        smo0(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is),
     &        'b2stbc_phys_smo'//chns)
	 enddo
	endif
	if (iout.eq.1 .or. (iout_b2npc9.eq.1 .or. iout_b2npco.eq.1)) then
       do is = 0, ns-1
	  write (chns,'(i2.2)') is
!	  do k=0,1
!	   write (chk,'(i1)') k
!	   call my_out(70,nx,ny,sna0(-1,-1,k,is),
!     &               'b2stbc_phys_sna'//chk//chns)
!	  enddo
	 call my_out(70,nx,ny,sna0(:,:,0,is)+
     &             sna0(:,:,1,is)*na(:,:,is),'b2stbc_phys_sna'//chns)
	 enddo
	endif
	if (iout.eq.1 .or. (iout_b2npht.eq.1 .or. iout_b2nph9.eq.1)) then
!	 do k=0,3
!	  write (chk,'(i1)') k
!	  call my_out(70,nx,ny,she0(-1,-1,k),'b2stbc_phys_she'//chk)
!	  call my_out(70,nx,ny,shi0(-1,-1,k),'b2stbc_phys_shi'//chk)
!	 enddo
	 call my_out(70,nx,ny,she0(:,:,0)+she0(:,:,1)*te+
     &   she0(:,:,2)*ne + she0(:,:,3)*te*ne,'b2stbc_phys_she')
	 call my_out(70,nx,ny,shi0(:,:,0)+shi0(:,:,1)*ti+
     &   shi0(:,:,2)*ni(:,:,0)+shi0(:,:,3)*ti*ni(:,:,0),
     &   'b2stbc_phys_shi')
	endif
	if (iout.eq.1 .or. (iout_b2npp7.eq.1 .or. iout_b2nppo.eq.1)) then
	 do k=0,3
	  write (chk,'(i1)') k
	  call my_out(70,nx,ny,sch0(-1,-1,k),'b2stbc_phys_sch'//chk)
	 enddo
!	 call my_out(70,nx,ny,sch0(:,:,0)+sch0(:,:,1)*po+
!     &   sch0(:,:,2)*ne+sch0(:,:,3)*po*ne(:,:),'b2stbc_phys_sch')
	endif                                                           !srv 11.09.09 }

* ..return
      if(ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbc_phys

      end subroutine b2stbc_phys

*-----------------------------------------------------------------------
*.specification
      subroutine b2stbc_fb (nx, ny, ns, ismain, BoRiS,
     1  vol, gs, na, ne, ni, te, ti, fna, fhe, fhi, fch, kinrgy,
     &  main_call, diagno)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_feedback
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, diagno
      real (kind=R8) :: BoRiS,
     *  vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  te(-1:nx, -1:ny), ti(-1:nx, -1:ny)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1)
      real (kind=R8) ::
     *  kinrgy(-1:nx,-1:ny,0:ns-1)
      logical main_call                                                  !xpb
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBC_FB provides linearised expressions for the sources of
*     density, momentum, electron heat, atom heat, and charge, due to
*     the boundary conditions.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, coreregno, coreregn2, 
     1  ncallfeedback, pfrregno1, pfrregno2, isfeedback, 
     2  ixl1, ixl2, ixr1, ixr2, solregno,
     3  fhiycore_kinetic_energy
      real (kind=R8) ::
     &  fnaycore, fnaycurr, fnaycore_alpha,
     1  fheycore, fheycurr, fhiycore, fhiycurr, fchycore, fchycurr,
     2  ndes, ndes_sol, ntot, cbsnainc, cbsnafac, den(0:1), rcbsna(0:1),
     3  amtol, neutral_sources_rescale, nesepm_sol, nesepm_pfr,
     4  nesepm, nesepm_curr, nesepm_alpha, fheycore_alpha,
     5  fhiycore_alpha, fchycore_alpha, private_flux_puff,
     6  nesepm_overshoot, save_gas_puff_while_overshot,
     7  save_pfr1_puff_while_overshot, save_pfr2_puff_while_overshot,
     8  save_core_flux_while_overshot, volrec_sol, volrec_alpha, 
     9  volrec_beta, volrec_overshoot, volrec_curr
      integer jsep, jxa, jxi
      logical zeroslope, feedback_namelist_used, file_ok
      real (kind=R8) ::
     1  private_flux_area, target_flux
      save ncall,
     &  fnaycore, fnaycore_alpha, fhiycore_kinetic_energy,
     1  fheycore, fhiycore, fchycore, nesepm_sol, solregno,
     1  coreregno, coreregn2, pfrregno1, pfrregno2, isfeedback,
     1  ncallfeedback, ndes, cbsnafac, zeroslope, den, rcbsna,
     1  nesepm, nesepm_alpha, jsep, jxa, jxi, ixl1, ixl2, ixr1, ixr2,
     1  fheycore_alpha, fhiycore_alpha, fchycore_alpha, ndes_sol,
     1  private_flux_puff, private_flux_area, nesepm_pfr,
     1  feedback_namelist_used, neutral_sources_rescale,
     1  nesepm_overshoot, save_gas_puff_while_overshot,
     1  save_pfr1_puff_while_overshot, save_pfr2_puff_while_overshot,
     1  save_core_flux_while_overshot, 
     1  volrec_sol, volrec_alpha, volrec_beta, volrec_overshoot
*   ..procedures
      intrinsic abs, max, mod, nint
      external subini, subend, xertst, ipgetr, ipgeti
      external b2xvsg
*   ..initialisation
      parameter (amtol = 1.0e-2_R8)
      data ncall/0/,
     1  fnaycore/-1.0e30_R8/, fnaycore_alpha/1.0e-1_R8/,
     1  fheycore/0.0_R8/, fheycore_alpha/1.0e-1_R8/,
     1  fhiycore/0.0_R8/, fhiycore_alpha/1.0e-1_R8/,
     1  fhiycore_kinetic_energy/0/, ndes_sol/0.0_R8/,
     1  fchycore/-1.0e30_R8/, fchycore_alpha/1.0e-1_R8/,
     1  coreregno/1/, coreregn2/4/, ncallfeedback /0/, ndes/0.0_R8/,
     1  pfrregno1/0/, pfrregno2/2/, isfeedback/0/, solregno/3/,
     2  cbsnafac/0.01_R8/, zeroslope/.false./,  nesepm_pfr/0.0_R8/,
     3  nesepm/0.0_R8/, nesepm_sol/0.0_R8/, nesepm_alpha/1.0e-3_R8/,
     4  private_flux_puff/0.0_R8/, neutral_sources_rescale/1.0_R8/, 
     2  feedback_namelist_used/.false./, nesepm_overshoot/0.0_R8/,
     3  save_gas_puff_while_overshot/0.0_R8/,
     4  save_core_flux_while_overshot/0.0_R8/,
     5  save_pfr1_puff_while_overshot/0.0_R8/,
     6  save_pfr2_puff_while_overshot/0.0_R8/,
     7  volrec_sol/0.0_R8/, volrec_alpha/1.0e-3_R8/,
     8  volrec_beta/1.0_R8/, volrec_overshoot/0.0_R8/
#ifdef COMPLICATED
      integer i
      real (kind=R8) ::
     .  teave1, teave2, neave1, neave2, tenext,
     .  tiave1, tiave2, niave1, niave2, tinext
#endif
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbc_fb')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2stbc_fnaycore', fnaycore)
       call ipgetr ('b2stbc_fnaycore_alpha', fnaycore_alpha)
       call xertst (0.0_R8.le.fnaycore_alpha,
     1   'faulty internal parameter fnaycore_alpha')
       call ipgetr ('b2stbc_fheycore', fheycore)
       call xertst (0.0_R8.le.fheycore, 
     1   'faulty internal parameter fheycore')
       call ipgetr ('b2stbc_fheycore_alpha', fheycore_alpha)
       call xertst (0.0_R8.le.fheycore_alpha,
     1   'faulty internal parameter fheycore_alpha')
       call ipgetr ('b2stbc_fhiycore', fhiycore)
       call xertst (0.0_R8.le.fhiycore, 
     1   'faulty internal parameter fhiycore')
       call ipgetr ('b2stbc_fhiycore_alpha', fhiycore_alpha)
       call xertst (0.0_R8.le.fhiycore_alpha,
     1   'faulty internal parameter fhiycore_alpha')
       call ipgeti ('b2tsbc_fhiycore_kinetic_energy',
     1  fhiycore_kinetic_energy)
       call ipgetr ('b2stbc_fchycore', fchycore)
       call ipgetr ('b2stbc_fchycore_alpha', fchycore_alpha)
       call xertst (0.0_R8.le.fchycore_alpha,
     1   'faulty internal parameter fchycore_alpha')
       call ipgetr ('b2stbc_nesepm', nesepm)
       call xertst (0.0_R8.le.nesepm, 
     1   'faulty internal parameter nesepm')
       call ipgetr ('b2stbc_nesepm_sol', nesepm_sol)
       call xertst (0.0_R8.le.nesepm_sol, 
     1   'faulty internal parameter nesepm_sol')
       call ipgetr ('b2stbc_nesepm_pfr', nesepm_pfr)
       call xertst (0.0_R8.le.nesepm_pfr, 
     1   'faulty internal parameter nesepm_pfr')
       call ipgetr ('b2stbc_volrec_sol', volrec_sol)
       call xertst (0.0_R8.le.volrec_sol, 
     1   'faulty internal parameter volrec_sol')
       call ipgetr ('b2stbc_nesepm_alpha', nesepm_alpha)
       call xertst (0.0_R8.le.nesepm_alpha,
     1   'faulty internal parameter nesepm_alpha')
       call ipgetr ('b2stbc_nesepm_overshoot', nesepm_overshoot)
       call xertst (1.0_R8.le.nesepm_overshoot.or.
     1  0.0_R8.eq.nesepm_overshoot,
     1  'faulty internal parameter nesepm_overshoot')
       call ipgetr ('b2stbc_volrec_alpha', volrec_alpha)
       call xertst (0.0_R8.le.volrec_alpha,
     1   'faulty internal parameter volrec_alpha')
       call ipgetr ('b2stbc_volrec_beta', volrec_beta)
       call ipgetr ('b2stbc_volrec_overshoot', volrec_overshoot)
       call xertst (1.0_R8.le.volrec_overshoot.or.
     1  0.0_R8.eq.volrec_overshoot,
     1  'faulty internal parameter volrec_overshoot')
       call ipgeti ('b2stbc_coreregno', coreregno)
       call ipgeti ('b2stbc_coreregn2', coreregn2)
       call ipgeti ('b2stbc_solregno', solregno)
       call ipgetr ('b2stbc_ndes_sol', ndes_sol)
       call xertst (0.0_R8.le.ndes_sol, 
     .  'faulty internal parameter ndes_sol')
       if(fnaycore.ne.-1.0e30_R8.or.fheycore.gt.0.0_R8.or.
     &    fhiycore.gt.0.0_R8.or.fchycore.ne.-1.0e30_R8.or.
     2    nesepm.gt.0.0_R8) then
         call xertst (cbirso.le.coreregno.and.
     1          coreregno.lt.cbirso+cbnrso-1,
     2    'faulty internal parameter coreregno')
         if (nncut.ge.2)
     >    call xertst (cbirso.le.coreregn2.and.
     1    coreregn2.lt.cbirso+cbnrso-1,
     2    'faulty internal parameter coreregn2')
       endif
       if(nesepm_sol.gt.0.0_R8.or.volrec_sol.gt.0.0_R8.or.
     1    ndes_sol.gt.0.0_R8) then
         call xertst (cbirno.le.solregno.and.
     1    solregno.le.cbirno+cbnrno-1,
     2    'faulty internal parameter solregno')
         call xertst(
     .    (ndes_sol.gt.0.0_R8 .and. volrec_sol.eq.0.0_R8 .and. 
     .     nesepm_sol.eq.0.0_R8) .or.
     .    (ndes_sol.eq.0.0_R8 .and. volrec_sol.gt.0.0_R8 .and. 
     .     nesepm_sol.eq.0.0_R8) .or.
     .    (ndes_sol.eq.0.0_R8 .and. volrec_sol.eq.0.0_R8 .and. 
     .     nesepm_sol.gt.0.0_R8),
     .    'Simultaneous SOL feedback schemes not allowed!')
       endif
       if(nesepm.gt.0.0_R8 .or. nesepm_sol.gt.0.0_R8 .or. 
     &    nesepm_pfr.gt.0.0_R8) then
          call get_jsep(nx,ny,jxi,jxa,jsep)
       endif
       call ipgeti ('b2stbc_ncallfeedback', ncallfeedback)
       call ipgetr ('b2stbc_ndes', ndes)
       call xertst (0.0_R8.le.ndes, 'faulty internal parameter ndes')
       call ipgetr ('b2stbc_cbsnafac', cbsnafac)
       call xertst (0.0_R8.lt.cbsnafac, 
     1   'faulty internal parameter cbsnafac')
       call ipgetr ('b2stbc_private_flux_puff', private_flux_puff)
       call xertst (0.0_R8.le.private_flux_puff,
     1  'faulty internal parameter private_flux_puff')
       call ipgeti ('b2stbc_pfrregno1', pfrregno1)
       call ipgeti ('b2stbc_pfrregno2', pfrregno2)
       call ipgeti ('b2stbc_isfeedback', isfeedback)
       call xertst (0.le.isfeedback.and.isfeedback.lt.ns,
     1  'faulty internal parameter isfeedback')
       if (ndes.gt.0.0_R8.or.private_flux_puff.gt.0.0_R8.or.
     &     nesepm_pfr.gt.0.0_R8) then
         call xertst (cbirso.le.pfrregno1.and.
     1    pfrregno1.le.cbirso+cbnrso-1,
     2    'faulty internal parameter pfrregno1')
         call xertst (cbirso.le.pfrregno2.and.
     1    pfrregno2.le.cbirso+cbnrso-1,
     2    'faulty internal parameter pfrregno2')
         if (pfrregno1.gt.cbirso) then
           ixl1 = nint(cbrbrk(pfrregno1-1)*nx)
         else
           ixl1 = 0
         endif
         ixr1 = nint(cbrbrk(pfrregno1)*nx)
         if (pfrregno2.gt.cbirso) then
           ixl2 = nint(cbrbrk(pfrregno2-1)*nx)
         else
           ixl2 = 0
         endif
         ixr2 = nint(cbrbrk(pfrregno2)*nx)
         call xertst (pfrregno1.eq.pfrregno2 .or.
     &    rightix(ixr1,-1).eq.ixl2 .or. rightix(ixr2,-1).eq.ixl1,
     .    'pfrregno1 and pfrregno2 are not connected!')
         call xertst(
     .    (ndes.gt.0.0_R8 .and. private_flux_puff.eq.0.0_R8 .and. 
     .     nesepm_pfr.eq.0.0_R8) .or.
     .    (ndes.eq.0.0_R8 .and. private_flux_puff.gt.0.0_R8 .and. 
     .     nesepm_pfr.eq.0.0_R8) .or.
     .    (ndes.eq.0.0_R8 .and. private_flux_puff.eq.0.0_R8 .and. 
     .     nesepm_pfr.gt.0.0_R8),
     .    'Simultaneous PFR feedback schemes not allowed!')
       endif
       call ipgetr ('b2mndr_rescale_neutrals_sources',
     *       neutral_sources_rescale)
       call xertst (0.0_R8.le.neutral_sources_rescale,
     1   'faulty input neutral_sources_rescale')
       inquire(file='../b2.feedback_save.parameters',exist=file_ok)
       if(file_ok) then
         open(99,file='../b2.feedback_save.parameters')
         call read_b2mod_feedback(99)
         close(99)
       endif
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
*    ..test sign of ni, ne, te, ti
       call b2xvsg (n2, ni(-1,-1,0), 1, 'ni0', '.gt.')
       call b2xvsg (n2, ni(-1,-1,1), 1, 'ni1', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
      endif

* DPC special treatment of core boundary condition
      if(fnaycore.ne.-1.0e30_R8) then
        fnaycurr=0.0_R8
        do ix=-1,nx
          if(mod(region(ix,-1,0),4).eq.1) then
            fnaycurr=fnaycurr+
     1       fna(topix(ix,-1),topiy(ix,-1),1,isfeedback)
          endif
        enddo
        if(ncall.eq.0) then
          if(saved_cbsna_core.ne.0.0_R8) then
            cbsna(0,isfeedback,coreregno)=saved_cbsna_core
            if(nncut.ge.2)
     1       cbsna(0,isfeedback,coreregn2)=saved_cbsna_core
          endif
          feedback_namelist_used=.true.
        endif
        cbsna(0,isfeedback,coreregno)=cbsna(0,isfeedback,coreregno)*
     *   (1.0_R8+fnaycore_alpha*
     1   (fnaycore-fnaycurr)/(abs(fnaycore)+abs(fnaycurr)+1.0e16_R8))
        if(nncut.ge.2)
     1   cbsna(0,isfeedback,coreregn2)=cbsna(0,isfeedback,coreregn2)*
     *    (1.0_R8+fnaycore_alpha*
     1    (fnaycore-fnaycurr)/(abs(fnaycore)+abs(fnaycurr)+1.0e16_R8))
        saved_cbsna_core=cbsna(0,isfeedback,coreregno)
        if(diagno.ge.1) then
          if(nncut.lt.2)
     >     write(*,*) 'fnaycurr= ',fnaycurr,
     1     ' ; target fnaycore= ',fnaycore,
     *     ' ; rescaling cbsna to ',cbsna(0,isfeedback,coreregno)
          if(nncut.ge.2)
     >     write(*,*) 'fnaycurr= ',fnaycurr,
     1     ' ; target fnaycore= ',fnaycore,
     *     ' ; rescaling cbsna to ',cbsna(0,isfeedback,coreregno),
     *     cbsna(0,isfeedback,coreregn2)
        endif
      endif
      if(fheycore.gt.0.0_R8) then
#ifdef COMPLICATED
        teave1=0.0_R8
        neave1=0.0_R8
        teave2=0.0_R8
        neave2=0.0_R8
#endif
        fheycurr=0.0_R8
        do ix=-1,nx
          if(mod(region(ix,-1,0),4).eq.1) then
#ifdef COMPLICATED
            teave1=teave1+te(ix,-1)*ne(ix,-1)*vol(ix,-1)
            neave1=neave1+ne(ix,-1)*vol(ix,-1)
            teave2=teave2+te(topix(ix,-1),topiy(ix,-1))*
     .                    ne(topix(ix,-1),topiy(ix,-1))*
     .                   vol(topix(ix,-1),topiy(ix,-1))
            neave2=neave2+ne(topix(ix,-1),topiy(ix,-1))*
     .                   vol(topix(ix,-1),topiy(ix,-1))
#endif
            fheycurr=fheycurr+fhe(topix(ix,-1),topiy(ix,-1),1)
          endif
        enddo
#ifdef COMPLICATED
        teave1=teave1/neave1
        teave2=teave2/neave2
        tenext=teave1+abs(teave1-teave2)*
     1   (fheycore-max(0.0_R8,fheycurr))/fheycore*0.5_R8
        if(ncall.lt.ncallfeedback.and.main_call.and.diagno.ge.1) then
          write(*,*) ncall,': CORETEFEEDBACK (-)',fheycore,fheycurr,
     1     teave1/ev,teave2/ev,tenext/ev
        else
          if(diagno.ge.1)
     >     write(*,*) ncall,': CORETEFEEDBACK (+)',fheycore,fheycurr,
     1     teave1/ev,teave2/ev,tenext/ev
          call xertst (0.lt.tenext,
     1     'faulty feedback for core electron temperature')
          do i=0,7
            cbshe(i,coreregno)=0.0_R8
            if(nncut.ge.2) cbshe(i,coreregn2)=0.0_R8
          enddo
          cbshe(0,coreregno)=tenext*1.0e30_R8
          if(nncut.ge.2) cbshe(0,coreregn2)=tenext*1.0e30_R8
          cbshe(1,coreregno)=-1.0e30_R8
          if(nncut.ge.2) cbshe(1,coreregn2)=-1.0e30_R8
        endif
#else
        if(ncall.eq.0) then
          if(saved_cbshe_core.ne.0.0_R8) then
            cbshe(0,coreregno)=saved_cbshe_core
            if(nncut.ge.2) cbshe(0,coreregn2)=saved_cbshe_core
          endif
          feedback_namelist_used=.true.
        endif
        cbshe(0,coreregno)=cbshe(0,coreregno)*
     *       (1.0_R8+fheycore_alpha*fheycore/fheycurr)/
     *       (1.0_R8+fheycore_alpha)
        if(nncut.ge.2) cbshe(0,coreregn2)=cbshe(0,coreregn2)*
     *       (1.0_R8+fheycore_alpha*fheycore/fheycurr)/
     *       (1.0_R8+fheycore_alpha)
        saved_cbshe_core=cbshe(0,coreregno)
        if(diagno.ge.1)
     >   write(*,*) 'Rescaling cbshe to ',cbshe(0,coreregno)
#endif
      endif
      if(fhiycore.gt.0.0_R8) then
#ifdef COMPLICATED
        tiave1=0.0_R8
        niave1=0.0_R8
        tiave2=0.0_R8
        niave2=0.0_R8
#endif
        fhiycurr=0.0_R8
        do ix=-1,nx
          if(mod(region(ix,-1,0),4).eq.1) then
#ifdef COMPLICATED
            tiave1=tiave1+ti(ix,-1)*ni(ix,-1,0)*vol(ix,-1)
            niave1=niave1+ni(ix,-1,0)*vol(ix,-1)
            tiave2=tiave2+ti(topix(ix,-1),topiy(ix,-1))*
     .                    ni(topix(ix,-1),topiy(ix,-1),0)*
     .                   vol(topix(ix,-1),topiy(ix,-1))
            niave2=niave2+ni(topix(ix,-1),topiy(ix,-1),0)*
     .                   vol(topix(ix,-1),topiy(ix,-1))
#endif
            fhiycurr=fhiycurr+fhi(topix(ix,-1),topiy(ix,-1),1)
            if (fhiycore_kinetic_energy.eq.1) then
              do is = 0, ns-1
                fhiycurr=fhiycurr+(1.0_R8-BoRiS)*
     .           fna(topix(ix,-1),topiy(ix,-1),1,is)*0.5_R8*
     .           (kinrgy(topix(ix,-1),topiy(ix,-1),is)+
     .            kinrgy(ix,-1,is))
              enddo
            endif
          endif
        enddo
#ifdef COMPLICATED
        tiave1=tiave1/niave1
        tiave2=tiave2/niave2
        tinext=tiave1+abs(tiave1-tiave2)*
     1         (fhiycore-max(0.0_R8,fhiycurr))/fhiycore*0.5_R8
        if(ncall.lt.ncallfeedback.and.main_call.and.diagno.ge.1) then
          write(*,*) ncall,': CORETIFEEDBACK (-)',fhiycore,fhiycurr,
     1     tiave1/ev,tiave2/ev,tinext/ev
        else
          if(diagno.ge.1)
     >     write(*,*) ncall,': CORETIFEEDBACK (+)',fhiycore,fhiycurr,
     1     tiave1/ev,tiave2/ev,tinext/ev
          call xertst (0.0_R8.lt.tinext,
     1           'faulty feedback for core ion temperature')
          do is=0,ns-1
            do i=0,7
              cbshi(i,is,coreregno)=0.0_R8
              if(nncut.ge.2) cbshi(i,is,coreregn2)=0.0_R8
            enddo
            cbshi(0,is,coreregno)=tinext*1.0e30_R8
            if(nncut.ge.2) cbshi(0,is,coreregn2)=tinext*1.0e30_R8
            cbshi(1,is,coreregno)=-1.0e30_R8
            if(nncut.ge.2) cbshi(1,is,coreregn2)=-1.0e30_R8
          enddo
        endif
#else
        if(ncall.eq.0) then
          if(saved_cbshi_core.ne.0.0_R8) then
            cbshi(0,0:ns-1,coreregno)=saved_cbshi_core
            if(nncut.ge.2) cbshi(0,0:ns-1,coreregn2)=saved_cbshi_core
          endif
          feedback_namelist_used=.true.
        endif
        cbshi(0,0:ns-1,coreregno)=
     *       cbshi(0,0:ns-1,coreregno)*
     *       (1.0_R8+fhiycore_alpha*fhiycore/fhiycurr)/
     *       (1.0_R8+fhiycore_alpha)
        if(nncut.ge.2) cbshi(0,0:ns-1,coreregn2)=
     *       cbshi(0,0:ns-1,coreregn2)*
     *       (1.0_R8+fhiycore_alpha*fhiycore/fhiycurr)/
     *       (1.0_R8+fhiycore_alpha)
        saved_cbshi_core=cbshi(0,0,coreregno)
        if(diagno.ge.1)
     >   write(*,*) 'Rescaling cbshi to ',cbshi(0,0:ns-1,coreregno)
#endif
      endif
cxpb current feedback not yet tested!
      if(fchycore.ne.-1.0e30_R8) then
        fchycurr=0.0_R8
        do ix=-1,nx
          if(mod(region(ix,-1,0),4).eq.1)
     >     fchycurr=fchycurr+fch(topix(ix,-1),topiy(ix,-1),1)
        enddo
        if(ncall.eq.0) then
          if(saved_cbsch_core.ne.0.0_R8) then
            cbsch(0,coreregno)=saved_cbsch_core
            if(nncut.ge.2) cbsch(0,coreregn2)=saved_cbsch_core
          endif
          feedback_namelist_used=.true.
        endif
        cbsch(0,coreregno)=
     *       cbsch(0,coreregno)*
     *       (1.0_R8+fchycore_alpha*fchycore/fchycurr)/
     *       (1.0_R8+fchycore_alpha)
        if(nncut.ge.2) cbsch(0,coreregn2)=
     *       cbsch(0,coreregn2)*
     *       (1.0_R8+fchycore_alpha*fchycore/fchycurr)/
     *       (1.0_R8+fchycore_alpha)
        saved_cbsch_core=cbsch(0,coreregno)
        if(diagno.ge.1)
     >   write(*,*) 'Rescaling cbsch to ',cbsch(0,coreregno)
      endif
      if(ndes.gt.0.0_R8) then
        if (ncall.eq.0) then
          if (saved_cbsna_pfr1.ne.0.0_R8.or.
     &        saved_cbsna_pfr2.ne.0.0_R8) then
            cbsna(0,isfeedback,pfrregno1)=saved_cbsna_pfr1
            if (pfrregno1.ne.pfrregno2)
     >       cbsna(0,isfeedback,pfrregno2)=saved_cbsna_pfr2
          endif
          feedback_namelist_used=.true.
        endif
        ntot = 0.0_R8
        do iy=-1,ny
          do ix=-1,nx
            do is = 0, ns-1
              if (zn(is).eq.zn(isfeedback).and.
     &            abs(am(is)-am(isfeedback)).lt.amtol)
     >       ntot = ntot + na(ix,iy,is)*vol(ix,iy)
            enddo
          enddo
        enddo
        cbsnainc=(ndes-ntot)*cbsnafac
        if(ncall.gt.2) then
          if((ntot-den(1))*(den(1)-den(0)).lt.0.0_R8) then
            zeroslope=.true.
            rcbsna(0)=cbsna(0,isfeedback,pfrregno1)
            if (pfrregno1.ne.pfrregno2)
     >       rcbsna(1)=cbsna(0,isfeedback,pfrregno2)
            if (pfrregno1.ne.pfrregno2) then
              if (diagno.ge.1) write(*,*) 
     1         'ZEROSLOPE:ncall,cbsna(0,',isfeedback,',ireg)',ncall,
     1          cbsna(0,isfeedback,pfrregno1),
     1          cbsna(0,isfeedback,pfrregno2)
            else
              if (diagno.ge.1) write(*,*) 
     1         'ZEROSLOPE:ncall,cbsna(0,',isfeedback,',ireg)',ncall,
     1          cbsna(0,isfeedback,pfrregno1)
            endif
          endif
          if(zeroslope.and.(ntot-ndes)*(den(1)-ndes).lt.0) then
            cbsna(0,isfeedback,pfrregno1)=rcbsna(0)
            if (pfrregno1.ne.pfrregno2)
     >       cbsna(0,isfeedback,pfrregno2)=rcbsna(1)
            cbsnainc=0
            zeroslope=.false.
            if (pfrregno1.ne.pfrregno2) then
              if(diagno.ge.1) write(*,*)
     >         'ZEROCROSS:ncall,cbsna(0,',isfeedback,',ireg)',ncall,
     1         cbsna(0,isfeedback,pfrregno1),
     1         cbsna(0,isfeedback,pfrregno2)
            else
              if (diagno.ge.1) write(*,*) 
     >         'ZEROCROSS:ncall,cbsna(0,',isfeedback,',ireg)',ncall,
     1         cbsna(0,isfeedback,pfrregno1)
            endif
          endif
        endif
        if(ncall.lt.ncallfeedback.and.diagno.ge.1.and.main_call) then
          write(*,*) ncall,': PFLUXPUFFFEEDBACK (-)',ntot,ndes,cbsnainc
        else
          if(diagno.ge.1)
     >     write(*,*)
     1     ncall,': PFLUXPUFFFEEDBACK (+)',ntot,ndes,cbsnainc
          cbsna(0,isfeedback,pfrregno1)= neutral_sources_rescale * 
     .     max(0.0_R8,cbsna(0,isfeedback,pfrregno1)+cbsnainc)
          if (pfrregno1.ne.pfrregno2)
     >     cbsna(0,isfeedback,pfrregno2)= neutral_sources_rescale * 
     .     max(0.0_R8,cbsna(0,isfeedback,pfrregno2)+cbsnainc)
        endif
        if (pfrregno1.ne.pfrregno2) then
          if (diagno.ge.1)
     >     write(*,*) 'cbsna(0,',isfeedback,',',pfrregno1,'),',
     .                'cbsna(0,',isfeedback,',',pfrregno2,') ',
     1     cbsna(0,isfeedback,pfrregno1),cbsna(0,isfeedback,pfrregno2)
        else
          if (diagno.ge.1)
     >     write(*,*) 'cbsna(0,',isfeedback,',',pfrregno1,'),',
     1     cbsna(0,isfeedback,pfrregno1)
        endif
        den(0)=den(1)
        den(1)=ntot
      endif
      if(nesepm.gt.0.0_R8) then
        if(ncall.eq.0) then
          if(saved_cbsna_core.ne.0.0_R8) then
            cbsna(0,isfeedback,coreregno)=saved_cbsna_core
            if(nncut.ge.2) 
     >      cbsna(0,isfeedback,coreregn2)=saved_cbsna_core
          endif
          feedback_namelist_used=.true.
        endif
        nesepm_curr = 0.5_R8 *
     *   (ne(jxa,jsep)+ne(topix(jxa,jsep),topiy(jxa,jsep)))
        if(save_core_flux_while_overshot.ne.0.0_R8) then
          cbsna(0,isfeedback,coreregno)=save_core_flux_while_overshot
          if(nncut.ge.2)
     >    cbsna(0,isfeedback,coreregn2)=save_core_flux_while_overshot
          save_core_flux_while_overshot=0.0_R8
        endif
        cbsna(0,isfeedback,coreregno)=cbsna(0,isfeedback,coreregno)*
     *   (1.0_R8+nesepm_alpha*nesepm/nesepm_curr)/
     *   (1.0_R8+nesepm_alpha)
        if(nncut.ge.2)
     *   cbsna(0,isfeedback,coreregn2)=cbsna(0,isfeedback,coreregn2)*
     *   (1.0_R8+nesepm_alpha*nesepm/nesepm_curr)/
     *   (1.0_R8+nesepm_alpha)
        saved_cbsna_core=cbsna(0,isfeedback,coreregno)
        if(nesepm_overshoot.ge.1.0_R8.and.
     1   nesepm_curr.gt.nesepm*nesepm_overshoot) then
          save_core_flux_while_overshot=cbsna(0,ismain,coreregno)
          cbsna(0,isfeedback,coreregno)=0.0_R8
          if (nncut.ge.2) cbsna(0,isfeedback,coreregn2)=0.0_R8
          if (diagno.ge.1) write(*,*) 'Overshot density! ',
     1     nesepm_curr,nesepm*nesepm_overshoot
        endif
        if(diagno.ge.1) then
          if(nncut.lt.2)
     >     write(*,'(a,1p,g20.10,a,g20.10)')
     1     'nesepm_curr= ',nesepm_curr,
     *   ' ; rescaling cbsna to ',cbsna(0,isfeedback,coreregno)
          if(nncut.ge.2)
     >     write(*,'(a,1p,g20.10,a,2g20.10)')
     1     'nesepm_curr= ',nesepm_curr,
     *   ' ; rescaling cbsna to ',cbsna(0,isfeedback,coreregno),
     *     cbsna(0,isfeedback,coreregn2)
        endif
      endif
      if (nesepm_sol.gt.0.0_R8.or.ndes_sol.gt.0.0_R8) then
c nesepm contains the desired separatrix density
        if (ncall.eq.0) then
c do this only once at the beginning
          if (saved_cbsna_sol.ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=saved_cbsna_sol
     .       *neutral_sources_rescale
c cbsna(0,isfeedback,solregno) contains the boundary condition we are
c changing under feedback control
          elseif (cbsna(0,isfeedback,solregno).ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=cbsna(0,isfeedback,solregno)
     .       *neutral_sources_rescale
          else
            cbsna(0,isfeedback,solregno)=
     1       1.0e18_R8*neutral_sources_rescale
            write(*,*) 'Warning!  ',
     1       'The gas puff was started with an initial rate of 1e18!'
          endif
          feedback_namelist_used=.true.
c make sure that we save the feedback value before the code stops
        endif
        ntot = 0.0_R8
        nesepm_curr = 0.0_R8
        if (nesepm_sol.gt.0.0_R8) then
          nesepm_curr = 0.5_R8 *
     *     (ne(jxa,jsep)+ne(topix(jxa,jsep),topiy(jxa,jsep)))
        else if (ndes_sol.gt.0.0_R8) then
          do iy = -1, ny
            do ix = -1, nx
              do is = 0, ns-1
                if (zn(is).eq.zn(isfeedback).and.
     &              abs(am(is)-am(isfeedback)).lt.amtol)
     >            ntot = ntot + na(ix,iy,is)*vol(ix,iy)
              enddo
            enddo
          enddo
        endif
c the current value of the separatrix density
        if(save_gas_puff_while_overshot.ne.0.0_R8) then
          cbsna(0,isfeedback,solregno)=save_gas_puff_while_overshot
          save_gas_puff_while_overshot=0.0_R8
        endif
        if (nesepm_sol.gt.0.0_R8) then
          cbsna(0,isfeedback,solregno)= 
     .     cbsna(0,isfeedback,solregno)*
     *      (1.0_R8+nesepm_alpha*nesepm_sol/nesepm_curr)/
     *      (1.0_R8+nesepm_alpha)
        else if (ndes_sol.gt.0.0_R8) then
          cbsna(0,isfeedback,solregno)= 
     .     cbsna(0,isfeedback,solregno)*
     *      (1.0_R8+nesepm_alpha*ndes_sol/ntot)/
     *      (1.0_R8+nesepm_alpha)
        endif
c the key part: update the gas puff based on the ration of desired
c value to the current value; nesepm_alpha is used to stabilise the
c system and defaults to 1e-3
        saved_cbsna_sol=cbsna(0,isfeedback,solregno)/
     .   neutral_sources_rescale
c we save the value
        if (nesepm_overshoot.ge.1.0_R8 .and.
     1    ((nesepm_sol.gt.0.0_R8 .and. 
     2      nesepm_curr.gt.nesepm_sol*nesepm_overshoot) .or.
     3     (ndes_sol.gt.0.0_R8 .and.
     4      ntot.gt.ndes_sol*nesepm_overshoot))) then
c if we have overshot the density by more than
c nesepm_sol*nesepm_overshoot then we force the gas puff to 0
          save_gas_puff_while_overshot=cbsna(0,isfeedback,solregno)
          cbsna(0,isfeedback,solregno)=0.0_R8
          if (diagno.ge.1 .and. nesepm_sol.gt.0.0_R8) 
     >     write(*,*) 'Overshot density! ',
     1     nesepm_curr,nesepm_sol*nesepm_overshoot
          if (diagno.ge.1 .and. ndes_sol.gt.0.0_R8) 
     >     write(*,*) 'Overshot density! ',
     1     ntot,ndes_sol*nesepm_overshoot
        endif
        if (diagno.ge.1 .and. nesepm_sol.gt.0.0_R8)
     >   write(*,'(a,1p,g20.10,a,g20.10,a,g20.10)')
     1    ' nesepm_curr= ',nesepm_curr,
     *    ' target = ',nesepm_sol,
     *    ' ; rescaling cbsna to ',cbsna(0,isfeedback,solregno)
        if (diagno.ge.1 .and. ndes_sol.gt.0.0_R8)
     >   write(*,'(a,1p,g20.10,a,g20.10,a,g20.10)')
     1    ' ntot= ',ntot,
     *    ' target = ',ndes_sol,
     *    ' ; rescaling cbsna to ',cbsna(0,isfeedback,solregno)
      endif
      if (nesepm_pfr.gt.0.0_R8) then
        if (ncall.eq.0) then
          private_flux_area=0
          if (saved_cbsna_pfr1.ne.0.0_R8.or.
     &        saved_cbsna_pfr2.ne.0.0_R8) then
            cbsna(0,isfeedback,pfrregno1)=saved_cbsna_pfr1
            if (pfrregno1.ne.pfrregno2)
     >       cbsna(0,isfeedback,pfrregno2)=saved_cbsna_pfr2
          endif
          feedback_namelist_used=.true.
        endif
        nesepm_curr = 0.5_R8 *
     *   (ne(jxa,jsep)+ne(topix(jxa,jsep),topiy(jxa,jsep)))
        if(save_pfr1_puff_while_overshot.ne.0.0_R8) then
          cbsna(0,isfeedback,pfrregno1)=save_pfr1_puff_while_overshot
          save_pfr1_puff_while_overshot=0.0_R8
        endif
        if(save_pfr2_puff_while_overshot.ne.0.0_R8) then
          cbsna(0,isfeedback,pfrregno2)=save_pfr2_puff_while_overshot
          save_pfr2_puff_while_overshot=0.0_R8
        endif
        cbsna(0,isfeedback,pfrregno1)= 
     .   cbsna(0,isfeedback,pfrregno1)*
     *   (1.0_R8+nesepm_alpha*nesepm_pfr/nesepm_curr)/
     *   (1.0_R8+nesepm_alpha)
        cbsna(0,isfeedback,pfrregno2)= 
     .   cbsna(0,isfeedback,pfrregno2)*
     *   (1.0_R8+nesepm_alpha*nesepm_pfr/nesepm_curr)/
     *   (1.0_R8+nesepm_alpha)
        saved_cbsna_pfr1=cbsna(0,isfeedback,pfrregno1)/
     .   neutral_sources_rescale
        saved_cbsna_pfr2=cbsna(0,isfeedback,pfrregno2)/
     .   neutral_sources_rescale
        if(nesepm_overshoot.ge.1.0_R8.and.
     1   nesepm_curr.gt.nesepm_pfr*nesepm_overshoot) then
          save_pfr1_puff_while_overshot=cbsna(0,isfeedback,pfrregno1)
          save_pfr2_puff_while_overshot=cbsna(0,isfeedback,pfrregno2)
          cbsna(0,isfeedback,pfrregno1)=0.0_R8
          cbsna(0,isfeedback,pfrregno2)=0.0_R8
          if (diagno.ge.1) write(*,*) 'Overshot density! ',
     1     nesepm_curr,nesepm_pfr*nesepm_overshoot
        endif
        if (diagno.ge.1)
     >   write(*,'(a,1p,g20.10,a,g20.10,a,g20.10)')
     1   'nesepm_curr= ',nesepm_curr,
     *    ' target = ',nesepm_pfr,
     *    ' ; rescaling cbsna to ',cbsna(0,isfeedback,pfrregno1),
     *    cbsna(0,isfeedback,pfrregno2)
      endif
      if(private_flux_puff.gt.0.0_R8.or.nesepm_pfr.gt.0.0_R8) then
        if(ncall.eq.0) then
          private_flux_area=0.0_R8
          if (saved_cbsna_pfr1.ne.0.0_R8.or.
     &        saved_cbsna_pfr2.ne.0.0_R8) then
            cbsna(0,isfeedback,pfrregno1)=saved_cbsna_pfr1
            if (pfrregno1.ne.pfrregno2)
     >       cbsna(0,isfeedback,pfrregno2)=saved_cbsna_pfr2
          endif
          feedback_namelist_used=.true.
          do ix=-1,nx
            if((mod(region(ix,-1,0),4).eq.3 .or.
     &         (mod(region(ix,-1,0),4).eq.0 .and. region(ix,-1,0).ne.0))
     &       .and.((ix.ge.ixl1 .and. ix.le.ixr1) .or.
     &             (ix.ge.ixl2 .and. ix.le.ixr2)))
     >       private_flux_area=private_flux_area+
     2        gs(topix(ix,-1),topiy(ix,-1),1)
          enddo
          if(diagno.ge.1) then
            write(*,*) 'private_flux_area = ',private_flux_area
            if (pfrregno1.ne.pfrregno2) then
              write(*,*) 'Original cbsna(0,isfeedback,pfrregno1&2) = ',
     1         cbsna(0,isfeedback,pfrregno1),
     2         cbsna(0,isfeedback,pfrregno2)
            else
              write(*,*) 'Original cbsna(0,isfeedback,pfrregno1) = ',
     1         cbsna(0,isfeedback,pfrregno1)
            endif
          endif
        endif
        if (nesepm_pfr.gt.0.0_R8) then
          nesepm_curr = 0.5_R8 *
     *     (ne(jxa,jsep)+ne(topix(jxa,jsep),topiy(jxa,jsep)))
          cbsna(0,isfeedback,pfrregno1)= 
     .     cbsna(0,isfeedback,pfrregno1)*
     *     (1.0_R8+nesepm_alpha*nesepm_pfr/nesepm_curr)/
     *     (1.0_R8+nesepm_alpha)
          saved_cbsna_pfr1=cbsna(0,isfeedback,pfrregno1)/
     .     neutral_sources_rescale
          if (pfrregno1.ne.pfrregno2) then
            cbsna(0,isfeedback,pfrregno2)= 
     .       cbsna(0,isfeedback,pfrregno2)*
     *       (1.0_R8+nesepm_alpha*nesepm_pfr/nesepm_curr)/
     *       (1.0_R8+nesepm_alpha)
            saved_cbsna_pfr2=cbsna(0,isfeedback,pfrregno2)/
     .       neutral_sources_rescale
          endif
        elseif (private_flux_puff.gt.0.0_R8) then
          target_flux=0.0_R8
          do ix=-1,nx
            if (leftix(ix,-1).eq.-2.and.
     &            ((ix.ge.ixl1 .and. ix.le.ixr1) .or.
     &             (ix.ge.ixl2 .and. ix.le.ixr2))) then
              do iy=-1,ny
                target_flux=target_flux
     1           -rza(rightix(ix,iy),rightiy(ix,iy),ismain)
     2           *fna(rightix(ix,iy),rightiy(ix,iy),0,ismain)
              enddo
            elseif (rightix(ix,-1).eq.nx+1.and.
     &            ((ix.ge.ixl1 .and. ix.le.ixr1) .or.
     &             (ix.ge.ixl2 .and. ix.le.ixr2))) then
              do iy=-1,ny
                target_flux=target_flux
     1           +rza(ix,iy,ismain)*fna(ix,iy,0,ismain)
              enddo
            endif
          enddo
          if(diagno.ge.1) write(*,*) 'target_flux = ',target_flux
          cbsna(0,isfeedback,pfrregno1)= neutral_sources_rescale*
     .     target_flux*private_flux_puff/private_flux_area
          if (pfrregno1.ne.pfrregno2)
     >     cbsna(0,isfeedback,pfrregno2)= neutral_sources_rescale*
     .     target_flux*private_flux_puff/private_flux_area
        endif
        if (pfrregno1.ne.pfrregno2) then
          if (diagno.ge.1) 
     >     write(*,*) 'cbsna(0,isfeedback,pfrregno1&2) = ',
     1     cbsna(0,isfeedback,pfrregno1),cbsna(0,isfeedback,pfrregno2)
        else
          if (diagno.ge.1) 
     >     write(*,*) 'cbsna(0,isfeedback,pfrregno1) = ',
     1     cbsna(0,isfeedback,pfrregno1)
        endif
      endif
      if (volrec_sol.gt.0.0_R8) then
        if (ncall.eq.0) then
          if (saved_cbsna_sol.ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=saved_cbsna_sol
     .       *neutral_sources_rescale
          elseif (cbsna(0,isfeedback,solregno).ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=cbsna(0,isfeedback,solregno)
     .       *neutral_sources_rescale
          else
            cbsna(0,isfeedback,solregno)=
     1       1.0e18_R8*neutral_sources_rescale
            write(*,*) 'Warning!  ',
     1       'The gas puff was started with an initial rate of 1e18!'
          endif
          feedback_namelist_used=.true.
        endif
        volrec_curr = cum_volrec
        if(volrec_curr.le.0.0_R8) then
          write(*,*) 'volrec_curr <= 0 ', volrec_curr
        else
          if(save_gas_puff_while_overshot.ne.0.0_R8) then
            cbsna(0,isfeedback,solregno)=save_gas_puff_while_overshot
            save_gas_puff_while_overshot=0.0_R8
          endif
          cbsna(0,isfeedback,solregno)= 
     .     cbsna(0,isfeedback,solregno)*
     *     (1.0_R8+volrec_alpha*(volrec_sol/volrec_curr)**volrec_beta)/
     *     (1.0_R8+volrec_alpha)
          saved_cbsna_sol=cbsna(0,isfeedback,solregno)/
     .     neutral_sources_rescale
          if(volrec_overshoot.ge.1.0_R8.and.
     1     volrec_curr.gt.volrec_sol*volrec_overshoot) then
            save_gas_puff_while_overshot=cbsna(0,isfeedback,solregno)
            cbsna(0,isfeedback,solregno)=0.0_R8
            if (diagno.ge.1) write(*,*) 'Overshot density! ',
     1       volrec_curr,volrec_sol*volrec_overshoot
          endif
        endif
        if (diagno.ge.1)
     >   write(*,*) 'volrec_curr= ',volrec_curr,
     *    ' target = ',volrec_sol,
     *    ' ; rescaling cbsna to ',cbsna(0,isfeedback,solregno)
      endif

* ..return
      if(ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      if(ncall.gt.1) return

* ..finalise feedback
      entry b2stbc_fb_final
      if(feedback_namelist_used) then
        call write_b2mod_feedback(99)
      endif
      return

*-----------------------------------------------------------------------
*.end b2stbc_fb

      end subroutine b2stbc_fb

*-----------------------------------------------------------------------
*.specification
      function bcint(fun,nx,ny,il,ib)
      use b2mod_types
      use b2mod_geo
      use b2mod_indirect
      use b2mod_boundary_namelist
      implicit none
*   ..output result (undefined upon entry)
      real (kind=R8) :: bcint
*   ..input arguments (unchanged on exit)
      integer nx, ny, il, ib
      real (kind=R8) :: fun(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     BCINT provides interpolated face-centred values for cell-centred
*     quantities (array fun), for use in boundary conditions.
*     The interpolation is done according to distance away from the
*     cell face. The boundary index is ib, and the cell index along
*     the boundary is il.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix, iy
*   ..procedures
      real (kind=R8) ::
     *  hy1
      external xertst, xerrab
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*  ..test input parameters
      call xertst (0.le.nx.and.0.le.ny, 'faulty parameter nx, ny')
      call xertst (0.lt.ib.and.nbc.ge.ib, 'faulty parameter ib')
      call xertst (bcstart(ib).le.il.and.bcend(ib).ge.il,
     .  'faulty parameter il')

* ..main computation
      if (bcchar(ib).eq.'S') then
        bcint =
     .   (fun(il,bcpos(ib))*
     .    hy1(topix(il,bcpos(ib)),topiy(il,bcpos(ib)))+
     .    fun(topix(il,bcpos(ib)),topiy(il,bcpos(ib)))*
     .    hy1(il,bcpos(ib)))
     .  /(hy1(il,bcpos(ib))+
     .    hy1(topix(il,bcpos(ib)),topiy(il,bcpos(ib))))
      elseif (bcchar(ib).eq.'N') then
        bcint =
     .   (fun(il,bcpos(ib))*
     .    hy1(bottomix(il,bcpos(ib)),bottomiy(il,bcpos(ib)))+
     .    fun(bottomix(il,bcpos(ib)),bottomiy(il,bcpos(ib)))*
     .    hy1(il,bcpos(ib)))
     .  /(hy1(il,bcpos(ib))+
     .    hy1(bottomix(il,bcpos(ib)),bottomiy(il,bcpos(ib))))
      elseif (bcchar(ib).eq.'E') then
        bcint =
     .   (fun(bcpos(ib),il)*
     .    hx(leftix(bcpos(ib),il),leftiy(bcpos(ib),il))+
     .    fun(leftix(bcpos(ib),il),leftiy(bcpos(ib),il))*
     .    hx(bcpos(ib),il))
     .  /(hx(bcpos(ib),il)+
     .    hx(leftix(bcpos(ib),il),leftiy(bcpos(ib),il)))
      elseif (bcchar(ib).eq.'W') then
        bcint =
     .   (fun(bcpos(ib),il)*
     .    hx(rightix(bcpos(ib),il),rightiy(bcpos(ib),il))+
     .    fun(rightix(bcpos(ib),il),rightiy(bcpos(ib),il))*
     .    hx(bcpos(ib),il))
     .  /(hx(bcpos(ib),il)+
     .    hx(rightix(bcpos(ib),il),rightiy(bcpos(ib),il)))
      else
        call xerrab ('bcint : Unknown orientation!')
      endif

* ..return
      return

*-----------------------------------------------------------------------
*.end bcint

      end function bcint
