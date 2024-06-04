*-----------------------------------------------------------------------
*.specification

      subroutine b2sral (nx, ny, ns,
     &  nscx, nscxmax, nsmin, nsmax, iscx, ismain, ismain0,
     &  dtim, BoRiS, facdrift, fac_ExB, fac_vis,                         !xpb
     &  vol, hx, hy, hy1, hz, qz, qc, qs, gs, pbs, bb, lnlam,            !srv 09.01.01
     &  na, ua, uadia, vedia, vadia, wadia, veecrb, vaecrb, ve, wedia,
     &  te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,                  !srv 28.10.02
     &  fna, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,        !srv 13.10.06 25.12.06 23.10.09
     &  fhe, fhi, fhm, fht, fnaPSch, fhePSch, fhiPSch, fch,
     &  fchanml_a, fchinert_a, fchvispar_a, fchvisper_a, fchvisq_a,      !srv 08.09.21
     &  fchin_a, fchdia, fchin, fch_p, fchvispar, fchvisper,             !srv 21.10.02, IYS 17.10.2023
     &  fchvisq, fchinert, fchanml, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  cdna, cdpa, cvsa_cl, cvsa_drho,                                  !som 18.08.21
     &  chci_al_ast, chci_a, cvsa_hAdp_albe, cvsa_hBdp_al,               !som 20.10.21
     &  cvsa_uAdp_albe, cvsa_uBdp_al,                                    !som 26.04.24
     &  cvla, chce, chve, chci, chvi, calf, csigin,                      !xpb
     &  rlsa, rlra, rlqa, rlcx, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &  rza, rz2, rpt, rpi,
     &  sna, smo, smq, she, shi, sch, sne,
     &  wrong_flow, main_call)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_tallies
      use b2mod_sources
     &     , only: rcxna, rcxhi, rcxmo
      use b2mod_input_profile
      use b2mod_numerics_namelist
      use b2mod_boundary_namelist
      use b2mod_time
      use b2mod_anomalous_transport
      use b2mod_work                                                     !srv 02.04.18 {
     &     , only: avisper_a, bvisper_a, cvisper_a
      use b2mod_plasma
     &     , only: ne2                                                   !srv 02.04.18 }
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_external
      use b2mod_b2_to_astra
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns,
     *  nscx, nscxmax, iscx(0:nscxmax-1), ismain, ismain0, nsmin, nsmax
      real (kind=R8) ::                                                  !xpb
     *  dtim, BoRiS, facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny),        !xpb
     &  fac_vis(-1:nx,-1:ny)
*     (geometry)
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny),
     &  hy(-1:nx,-1:ny), hy1(-1:nx,-1:ny), hz(-1:nx,-1:ny),
     &  pbs(-1:nx,-1:ny,0:1),
     &  qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny,0:1), qs(-1:nx,-1:ny,0:1),
     &  bb(-1:nx,-1:ny,0:3), gs(-1:nx,-1:ny,0:2)
*     (state)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1), ve(-1:nx,-1:ny),
     1  uadia(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     1  vadia(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     2  vedia(-1:nx,-1:ny,0:1,0:1), veecrb(-1:nx,-1:ny,0:1,0:1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  kinrgy(-1:nx,-1:ny,0:ns-1), lnlam(-1:nx,-1:ny),
     &  wedia(-1:nx,-1:ny,0:1,0:1),
     &  wadia(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  vaecrb(-1:nx,-1:ny,0:1,0:1,0:ns-1)
*     (fluxes)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  fna_32(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  fna_52(-1:nx,-1:ny,0:1,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1,0:1),
     &  fhm(-1:nx,-1:ny,0:1,0:1,0:ns-1), fht(-1:nx,-1:ny,0:1,0:1),
     &  fhi(-1:nx,-1:ny,0:1,0:1), fch(-1:nx,-1:ny,0:1,0:1),
     &  fch_p(-1:nx,-1:ny,0:1,0:1),
     &  fni_32(-1:nx,-1:ny,0:1,0:1), fni_52(-1:nx,-1:ny,0:1,0:1),
     *  fne_32(-1:nx,-1:ny,0:1,0:1), fne_52(-1:nx,-1:ny,0:1,0:1),
     &  floe_noc(-1:nx,-1:ny,0:1,0:1), floi_noc(-1:nx,-1:ny,0:1,0:1),    !srv 28.10.02
     &  fna_fcor(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  fnaPSch(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  fna_mdf(-1:nx,-1:ny,0:1,0:1,0:ns-1),                             !srv 13.10.06 25.12.06
     &  fna_nodrift(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  fna_he(-1:nx,-1:ny,0:1,0:1,0:ns-1),                              !srv 25.02.00 24.10.02
     &  fhePSch(-1:nx,-1:ny,0:1,0:1), fhiPSch(-1:nx,-1:ny,0:1,0:1),      !srv 13.10.06
     &  fhe_mdf(-1:nx,-1:ny,0:1,0:1), fhi_mdf(-1:nx,-1:ny,0:1,0:1),      !srv 23.10.09
     &  fna_eir(-1:nx,-1:ny,0:1,0:1, 0:ns-1),
     &  fne_eir(-1:nx,-1:ny,0:1,0:1),
     &  fhe_eir(-1:nx,-1:ny,0:1,0:1), fhi_eir(-1:nx,-1:ny,0:1,0:1)
*     (transport coefficients)                                           !xpb
      real (kind=R8) ::                                                  !xpb
     1  chce(-1:nx,-1:ny,0:1,0:1), chve(-1:nx,-1:ny,0:1,0:1),            !xpb
     2  chci(-1:nx,-1:ny,0:1,0:1), chvi(-1:nx,-1:ny,0:1,0:1),            !xpb
     3  calf(-1:nx,-1:ny,0:1,0:1), cdna(-1:nx,-1:ny,0:1,0:1,0:ns-1),     !xpb
     4  cvsa_cl(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     5  cvla(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  cvsa_drho(-1:nx,-1:ny,0:1,0:1,0:ns-1),                           !som 18.08.21
     &  chci_al_ast(-1:nx,-1:ny,0:1,0:1,0:nspecies-1),                   !som 20.10.21
     &  chci_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),                              !som 20.10.21
     &  cvsa_hAdp_albe(-1:nx,-1:ny,0:1,0:1,0:nspecies-1,0:nspecies-1),   !som 20.10.21
     &  cvsa_hBdp_al(-1:nx,-1:ny,0:1,0:1,0:nspecies-1),                  !som 20.10.21
     &  cvsa_uAdp_albe(-1:nx,-1:ny,0:1,0:1,0:nspecies-1,0:nspecies-1),   !som 26.04.24
     &  cvsa_uBdp_al(-1:nx,-1:ny,0:1,0:1,0:nspecies-1),                  !som 26.04.24
     6  cdpa(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  csigin(-1:nx,-1:ny,0:1,0:1,0:ns-1,0:ns-1)
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
     &  fchanml_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),                           !srv 08.09.21 {
     &  fchinert_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  fchvispar_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  fchvisper_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  fchvisq_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),                           !srv 08.09.21 }
     &  fchin_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),                             !IYS 17.10.2023
     &  fchvispar(-1:nx,-1:ny,0:1,0:1),                                  !srv 21.10.02 {
     &  fchvisper(-1:nx,-1:ny,0:1,0:1), fchin(-1:nx,-1:ny,0:1,0:1),
     &  fchvisq(-1:nx,-1:ny,0:1,0:1), fchdia(-1:nx,-1:ny,0:1,0:1),
     &  fchinert(-1:nx,-1:ny,0:1,0:1), fchanml(-1:nx,-1:ny,0:1,0:1)      !srv 21.10.02 }
      real (kind=R8) ::
     *  rza(-1:nx,-1:ny,0:ns-1), rz2(-1:nx,-1:ny,0:ns-1),
     &  rpt(-1:nx,-1:ny,0:ns-1), rpi(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *  sna(-1:nx,-1:ny,0:1,0:ns-1), smo(-1:nx,-1:ny,0:3,0:ns-1),
     &  smq(-1:nx,-1:ny,0:3,0:ns-1), she(-1:nx,-1:ny,0:3),
     &  shi(-1:nx,-1:ny,0:3), sch(-1:nx,-1:ny,0:3), sne(-1:nx,-1:ny,0:1)
      logical wrong_flow                                                 !xpb
      logical main_call                                                  !dpc
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SRAL provides linearised expressions for the sources of
*     density, parallel momentum, electron heat, atom heat, and charge.
*     Additionally, an unbalanced electron particle source can be
*     provided.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, k, is, ix, iy, sources_inputfile
      integer guard_flows, numerics_namelist
      integer use_astra
      integer, save :: elm_count
      integer, save :: style = 2                                         !srv 28.06.08
      integer, save :: iout = 0, iout_b2wdat = 0                         !srv 30.06.08, IYS 06.11.2017
      integer, save :: style_b2mndt = 1                                  !srv 07.07.21
      integer, save :: use_eirene =0, save_f30 = 0, save_f31 = 0         !srv 14.07.10
      integer, save :: no_fluid_neutral_sources = 0
      character chns*3, chk*1                                            !srv 30.06.08
      real (kind=R8) ::
     *  snatot(0:ns-1),                                                  !srv 18.05.09
     1  rsa(-1:nx,-1:ny,0:ns-1),
     2  rra(-1:nx,-1:ny,0:ns-1),
     3  rcx(-1:nx,-1:ny,0:ns-1),
     4  smq0(-1:nx,-1:ny,0:3,0:ns-1), smq0_ion(-1:nx,-1:ny,0:3,0:ns-1),  !srv 11.09.09
     4  smq0_rec(-1:nx,-1:ny,0:3,0:ns-1), wrk(-1:nx,-1:ny),
     5  rqa(-1:nx,-1:ny,0:ns-1), rrd(-1:nx,-1:ny,0:ns-1),
     6  rbr(-1:nx,-1:ny,0:ns-1), fni(-1:nx,-1:ny,0:1,0:1),
     7  fne(-1:nx,-1:ny,0:1,0:1),
     8  sna0(-1:nx,-1:ny,0:1,0:ns-1), sna0_ion(-1:nx,-1:ny,0:1,0:ns-1), !srv 11.09.09
     8  sna0_rec(-1:nx,-1:ny,0:1,0:ns-1), sne0(-1:nx,-1:ny,0:1),
     9  shi0(-1:nx,-1:ny,0:3), shi0_ion(-1:nx,-1:ny,0:3),
     &  shi0_rec(-1:nx,-1:ny,0:3), she0(-1:nx,-1:ny,0:3),
     &  smo0(-1:nx,-1:ny,0:3,0:ns-1), sch0(-1:nx,-1:ny,0:3)
      real (kind=R8), allocatable :: sna0_cx(:,:,:,:), smo0_cx(:,:,:,:), !IYS 06.11.2017
     &  shi0_cx(:,:,:)
      logical catch_up
      character*32 hlp_frm
      save ncall, sources_inputfile, guard_flows, numerics_namelist
      save use_astra
*   ..procedures
      real (kind=R8) ::
     *  smin, smax
      external xertst, sfill, b2saxpy, smin, smax
      external b2xvsg, b2xvff, b2xvps, b2xpfe, b2xpfi, b2sqel, b2sqcx,
     &  b2stbc, b2stbm, b2stbr, b2stel, b2stcx, ipgeti
*   ..initialisation
      data ncall/0/, sources_inputfile/0/, guard_flows/2/, use_astra/0/,
     *  numerics_namelist/0/

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
      call xertst (0.le.ismain0.and.ismain0.lt.ns,
     . 'invalid main neutral species index ismain0')
      call xertst (is_neutral(ismain0).or.ismain.eq.ismain0,
     . 'invalid main neutral species ismain0; must be neutral')
*   ..test facdrift (not on first call)
      if (ncall.gt.0) then
        call xertst (0.0_R8.le.minval(facdrift) .and.
     &       1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
        call xertst (0.0_R8.le.minval(fac_ExB) .and.
     &       1.0_R8.ge.maxval(fac_ExB), 'faulty argument fac_ExB')
       endif
*   ..inputfile for profile sources
      if (ncall.eq.0) then
        call ipgeti ('b2sral_inputfile', sources_inputfile)
        call alloc_input_profile(nx,ny,ns)
        call ipgeti ('b2sral_style', style)                              !srv 28.07.08
        call ipgeti ('b2sral_iout', iout)                                !srv 30.06.08
        call ipgeti ('b2mndr_astra', use_astra)
        call ipgeti ('b2mndr_eirene', use_eirene)                        !srv 14.07.10 {
        call ipgeti ('eirene_savef30', save_f30)
        call ipgeti ('eirene_savef31', save_f31)                         !srv 14.07.10 }
        call ipgeti ('b2sral_no_fluid_neutral_sources',
     .                                       no_fluid_neutral_sources)
        call ipgeti ('b2wdat_iout', iout_b2wdat)                         !iys 06.11.17
        call ipgeti ('b2mndt_style', style_b2mndt)                       !srv 07.07.21
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol, hx, hy, hy1
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
       call b2xvsg (n2, hy1, 1, 'hy1', '.gt.')
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
      call b2xpfe (nx, ny, ns, ns_ext, qe,
     &  rza, za_ext, fna, fa_ext, fch, fne)
      call b2xpfi (nx, ny, ns, fna, fni)

*   ..test facdrift and fac_ExB
      if (ncall.eq.0) then                                               !xpb
        call ipgeti ('b2news_guard_flows', guard_flows)                  !xpb
        call ipgeti ('b2srdt_numerics_namelist', numerics_namelist)
        if(numerics_time_mod.gt.0.0_R8) then
          elm_count=int(tim/numerics_time_mod)
        else
          elm_count=0
        endif
        if(numerics_time_mod.gt.0.0_R8) then
          catch_up=mod(tim,numerics_time_mod).ge.numerics_time_switch
     &      .and.numerics_time_switch.gt.0.0_R8
        else if (numerics_time_switch.gt.0.0_R8) then
          catch_up=tim.ge.numerics_time_switch
        else
          catch_up=.false.
        end if
        do while(catch_up)
          write(*,*) 'b2sral: catching up ...'
          call read_b2mod_numerics_namelist(nx,ny,ns,nsmin,nsmax,
     &                                     .false.)
          call write_b2mod_numerics_namelist
          if(numerics_time_mod.gt.0.0_R8) then
            catch_up=mod(tim,numerics_time_mod).ge.numerics_time_switch
     &        .and.numerics_time_switch.gt.0.0_R8
          else if (numerics_time_switch.gt.0.0_R8) then
            catch_up=tim.ge.numerics_time_switch
          else
            catch_up=.false.
          endif
        enddo
      else
        if(numerics_time_mod.gt.0.0_R8) then
          write(*,'(a,1p,3e12.4,2i6)')
     .     'ELM: ',tim,mod(tim,numerics_time_mod),
     1     numerics_time_switch,int(tim/numerics_time_mod),elm_count
          catch_up=
     1    (mod(tim,numerics_time_mod).ge.numerics_time_switch.and.
     1     numerics_time_switch.gt.0.0_R8).or.
     2     int(tim/numerics_time_mod).gt.elm_count
        else if (numerics_time_switch.gt.0.0_R8) then
          catch_up=tim.ge.numerics_time_switch
          if(catch_up) write(*,'(a,a)')
     .     'Switching to next numerics namelist file ',
     .      trim(numerics_filename)
        else
          catch_up=.false.
        end if
        if(catch_up) then
          if(numerics_time_mod.gt.0.0_R8)
     >     elm_count=int(tim/numerics_time_mod)
          call read_b2mod_numerics_namelist(nx,ny,ns,nsmin,nsmax,
     &                                     .false.)
          call write_b2mod_numerics_namelist
        endif
      endif                                                              !xpb
      call xertst (0.0_R8.le.minval(facdrift) .and.
     &     1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
      call xertst (0.0_R8.le.minval(fac_ExB) .and.
     &     1.0_R8.ge.maxval(fac_ExB), 'faulty argument fac_ExB')

      if(use_eirene.ne.0 .or. save_f31.ne.0) then                        !srv 14.07.10 {
*   ..recompute currents
       call b2tanml (nx, ny, ns, ismain, vol, hx, hy, qs,                !srv 20.09.06 }!som 02.11.21
     &               csig_an, po, ne, na, fchanml_a, fchanml)            !srv 08.09.21
       call b2tvspa (nx, ny, ns, ismain, vol, hx, hy, hy1, qs, bb,       !srv 19.02.21
     &               na, ua, ti, cvsa_cl, cvsa_drho, fac_vis,            !som 18.08.21
     &               cvsa_uAdp_albe, cvsa_uBdp_al,                       !som 26.04.24
     &               fchvispar_a, fchvispar)
       call b2tiner (nx, ny, ns, ismain, vol, gs, hx, hy, qs, bb,
     &               na, ua, facdrift, fchinert_a, fchinert)
       call b2tvspr (nx, ny, ns, ismain, fac_vis,                        !srv 10.10.17 19.02.21
     &               vol, hy1, po, ti, na,
     &               avisper_a, bvisper_a, cvisper_a,
     &               fchvisper_a, fchvisper)
       call b2tvsq  (nx, ny, ns, ismain, fac_vis,                        !srv 10.10.17 19.02.21
     &               ne2, lnlam, te, ti, ni(-1,-1,1), chci, na,          !som 21.10.21
     &               chci_al_ast, chci_a, cvsa_hAdp_albe, cvsa_hBdp_al,  !som 20.10.21
     &               fchvisq_a, fchvisq)
      call b2tinnt (nx, ny, ns, nscx, iscx,
     &              facdrift, fac_ExB,
     &              vol, gs, hx, hy, hy1, hz, qz, qc, qs, bb,
     &              te, po, ti, na, fna, csigin, fchin_a, fchin)         !IYS 17.10.2023
csw 30sep2013, added cdpa in call, uadia
       call b2tfrn (nx, ny, ns, guard_flows,                             !srv 19.02.21
     &  facdrift, fac_ExB, rza, rz2,                                     !srv 13.01.17
     &  vol, hx, hy, gs, qc, qs, bb, pbs,
     &  ua, na, ne, te, ti, po, fch_p,
     &  fchvispar_a, fchvisper_a, fchinert_a, fchanml_a, fchvisq_a,      !srv 08.09.21
     &  fchin_a, fne, cvla, cdna, cdpa, chce, chci,
     &  veecrb, vaecrb, wadia, ve, wedia,                                !srv 05.06.18
     &  fna_eir, fne_eir, fhe_eir, fhi_eir)
      endif                                                              !srv 14.07.10 }

* ..compute source coefficients
*   ..compute electron rate coefficients
      call b2sqel (nx, ny, ns, ismain, ev, te,
     &  rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &  rsa, rra, rqa, rrd, rbr, rza, rz2, rpt, rpi)
*   ..compute standard contributions due to boundaries
      wrong_flow = .false.                                               !xpb
      if(style.eq.0 .or. style.eq.2) then                                !srv 28.07.08 03.11.17
       if(ncall.eq.0) write(*,*) ' style=',style,' b2stbc called'
       call b2stbc (nx, ny, ns, ismain, ismain0,
     &  BoRiS, facdrift, fac_ExB,                                        !xpb
     &  na, ua, uadia, vedia, vadia, veecrb, vaecrb,
     .  te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,
     1  fna, fna_mdf, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fhe, fhi, fhe_mdf, fhi_mdf, fhm, fch, fne, fni, fna_fcor,        !srv 16.02.12
     &  fna_nodrift, fna_he, fna_eir, fnaPSch, fhePSch, fhiPSch,         !srv 25.02.00 13.10.06
     &  fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,   !srv 21.10.02
     &  fchanml,
     2  cdna, chce, chve, chci, chvi, calf,                              !xpb
     &  sna, smo, she, shi, sch, sne, wrong_flow, main_call)             !xpb
       if(wrong_flow) write(*,*)
     &  'b2sral [1]: wrong_flow returned from b2stbc'
      elseif(style.eq.1) then
!! FIXME : Need to add sne in argument list here ???
       if(ncall.eq.0) write(*,*) ' style=',style,' b2stbc_spb called'
       call b2stbc_spb (nx, ny, ns, facdrift,                            !srv 09.01.01
     &  na, ua, uadia, vaecrb, vadia, te, ti, po, ne,                    !srv 14.10.99 12.05.00 09.01.01 28.10.02
     1  ni, floe_noc, floi_noc, fna, fna_fcor, fnaPSch, fhe, fhi,        !srv 13.10.06
     &  fhePSch, fhiPSch, fne, fni, fch,                                 !srv 25.02.00 21.10.02
     &  fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,   !srv 21.10.02
     &  sna, smo, she, shi, sch, wrong_flow)                             !srv 09.07.99 01.07.05
       if(wrong_flow) write(*,*)
     &  'b2sral [2]: wrong_flow returned from b2stbc_spb'
      else
       call xertst (0.le.style.and.style.le.2,                           !srv 03.11.17
     .  'invalid switch b2sral_style')
      endif                                                              !srv 28.07.08 }
      if (use_astra.ne.0) then
        se  =0.0_R8                                                      !srv 28.07.11
        seh =0.0_R8                                                      !srv 11.06.13
        sih =0.0_R8                                                      !srv 11.06.13
        sheNB_astra=0.0_R8                                               !srv 21.06.12 {
        shiNB_astra=0.0_R8
        do iy = -1,ny
          do ix = -1,nx
            if(on_closed_surface(ix,iy)) then                            !srv 23.10.12
              sheNb_astra(ix,iy,0)=vol(ix,iy)*PeNB_astra(iy)*1.0e6_R8
              shiNb_astra(ix,iy,0)=vol(ix,iy)*PiNB_astra(iy)*1.0e6_R8
            endif
          enddo
        enddo
        call b2saxpy (n2*4, 1.0_R8, sheNb_astra, 1, she, 1)
        call b2saxpy (n2*4, 1.0_R8, shiNb_astra, 1, shi, 1)              !srv 21.06.12 }
      endif
*   ..initialise smq to 0
      call sfill (n2*4*ns, 0.0_R8, smq, 1)
*   ..compute special contributions (boundaries and recycling)
      call b2stbr (nx, ny, ns,
     &  nscx, nscxmax, iscx, ismain, dtim, BoRiS,
     &  vol, hx, hy, hy1, hz, gs, qz, pbs,
     &  na, ua, te, ti, po, ne, ni,
     &  fna, fna_mdf, fhe, fhi, fht, fch, fne, fni, rsa, bb,             !iyv 07.03.18
     &  sna0, smo0, she0, shi0, sch0, sne0, main_call)
*   ..add to previous contributions
      call b2saxpy (n2*2*ns, 1.0_R8, sna0, 1, sna, 1)
      call b2saxpy (n2*4*ns, 1.0_R8, smo0, 1, smo, 1)
      call b2saxpy (n2*2, 1.0_R8, sne0, 1, sne, 1)
      call b2saxpy (n2*4, 1.0_R8, she0, 1, she, 1)
      call b2saxpy (n2*4, 1.0_R8, shi0, 1, shi, 1)
      call b2saxpy (n2*4, 1.0_R8, sch0, 1, sch, 1)
      do is=0,ns-1                                                       !srv 18.05.09 {
        snatot(is)=0.0_R8
        do ix=-1,nx
          do iy=-1,ny
            snatot(is)=snatot(is)+sna(ix,iy,0,is)+
     &       sna(ix,iy,1,is)*na(ix,iy,is)
          enddo
        enddo
      enddo
      write(hlp_frm,'(a,i3,a)') '(a,1p,',ns,'(1x,1e14.7))'
      if(main_call.or.wrong_flow)                                   !xpb !srv 18.05.09 }
     > write(*,hlp_frm) 'snatot',(snatot(is),is=0,ns-1)
*   ..compute more special contributions
      call b2stbm (nx, ny, ns, ismain, dtim,
     &  vol, gs, hx, hy, qz, qc, pbs, bb,
     &  na, ua, te, ti, po, ne, ni, fna, fhe, fhi, fch, fne, fni,
     &  sna0, smo0, she0, shi0, sch0, sne0, main_call)
*   ..add to previous contributions
      call b2saxpy (n2*2*ns, 1.0_R8, sna0, 1, sna, 1)
      call b2saxpy (n2*4*ns, 1.0_R8, smo0, 1, smo, 1)
      call b2saxpy (n2*2, 1.0_R8, sne0, 1, sne, 1)
      call b2saxpy (n2*4, 1.0_R8, she0, 1, she, 1)
      call b2saxpy (n2*4, 1.0_R8, shi0, 1, shi, 1)
      call b2saxpy (n2*4, 1.0_R8, sch0, 1, sch, 1)
*   ..compute sources by inputfiles (heating,...)
      if (sources_inputfile.eq.1) then
        call source_input(nx,ny,ns,
     &   sna0,smo0,she0,shi0,sch0,sne0,vol,ncall)
*   ..add to previous contributions
        call b2saxpy (n2*2*ns, 1.0_R8, sna0, 1, sna, 1)
        call b2saxpy (n2*4*ns, 1.0_R8, smo0, 1, smo, 1)
        call b2saxpy (n2*2, 1.0_R8, sne0, 1, sne, 1)
        call b2saxpy (n2*4, 1.0_R8, she0, 1, she, 1)
        call b2saxpy (n2*4, 1.0_R8, shi0, 1, shi, 1)
        call b2saxpy (n2*4, 1.0_R8, sch0, 1, sch, 1)
      endif

*   ..compute sources due to electron-atom processes
*     (parallel momentum source must be pure friction)
      call b2stel (nx, ny, ns, ismain, BoRiS, bb, vol, hx, hy, hz, na,   !srv 09.01.01
     & ua, wadia, vadia, vaecrb, te, ti, po, ne, ni, facdrift, fac_ExB,  !srv 07.07.09
     & rsa, rra, rqa, rrd, rbr, rza, rz2, rpt, rpi, rlsa, rlra, rlqa,
     & sna0_ion, sna0_rec, smq0_ion, smq0_rec, she0, shi0_ion, shi0_rec) !srv 11.09.09 {
*   ..add to previous contributions
      call b2saxpy (n2*2*ns, 1.0_R8, sna0_ion, 1, sna, 1)
      call b2saxpy (n2*2*ns, 1.0_R8, sna0_rec, 1, sna, 1)
      call b2saxpy (n2*4*ns, 1.0_R8, smq0_ion, 1, smq, 1)
      call b2saxpy (n2*4*ns, 1.0_R8, smq0_rec, 1, smq, 1)
      call b2saxpy (n2*4, 1.0_R8, she0, 1, she, 1)
      call b2saxpy (n2*4, 1.0_R8, shi0_ion, 1, shi, 1)
      call b2saxpy (n2*4, 1.0_R8, shi0_rec, 1, shi, 1)                   !srv 11.09.09 }
*     (no contribution to sch or sne)
*   ..compute sources due to charge exchange
      rcxnareg=0.0_R8
      rcxhireg=0.0_R8
      rcxmoreg=0.0_R8
      rcxna=0.0_R8
      rcxhi=0.0_R8
      rcxmo=0.0_R8

      if (iout_b2wdat .eq. 4) then  ! IYS 06.11.2017
        if (.not. allocated(shi0_cx)) then
          write(*,*) 'Allocating arrays for CX'
          allocate (shi0_cx(-1:nx,-1:ny,0:3))
          allocate (sna0_cx(-1:nx,-1:ny,0:1,0:ns-1))
          allocate (smo0_cx(-1:nx,-1:ny,0:3,0:ns-1))
          call sfill (n2*2*ns, 0.0_R8, sna0_cx, 1)
          call sfill (n2*4*ns, 0.0_R8, smo0_cx, 1)
          call sfill (n2*4, 0.0_R8, shi0_cx, 1)
        endif
      endif  ! IYS 06.11.2017

      do k = 0, nscx-1
*    ..compute charge exchange rate coefficients
        call b2sqcx (nx, ny, ns, ev, am(iscx(k)), ti,
     &   rlcx(-1,-1,0,0,k), rcx)
*    ..compute sources due to charge exchange
*      (parallel momentum source must be pure friction)
        call b2stcx (nx, ny, ns, iscx(k), ismain, hz, BoRiS, vol, na,    !srv 09.01.01
     &    ua, wadia, vadia, vaecrb, ti, ni, rcx, sna0, smq0, shi0)       !srv 07.07.09
*    ..add to previous contributions
        call b2saxpy (n2*2*ns, 1.0_R8, sna0, 1, sna, 1)
        call b2saxpy (n2*4*ns, 1.0_R8, smq0, 1, smq, 1)
        call b2saxpy (n2*4, 1.0_R8, shi0, 1, shi, 1)
*      (no contribution to she, sne or sch)
        if (iout_b2wdat .eq. 4) then  ! IYS 06.11.2017
          write(*,*) 'Summing CX sources over hydrogenic neutrals'
          call b2saxpy (n2*2*ns, 1.0_R8, sna0, 1, sna0_cx, 1)
          call b2saxpy (n2*4*ns, 1.0_R8, smq0, 1, smo0_cx, 1)
          call b2saxpy (n2*4, 1.0_R8, shi0, 1, shi0_cx, 1)
        endif                         ! IYS 06.11.2017
      enddo

      if (iout_b2wdat .eq. 4) then  ! IYS 06.11.2017
        write(*,*) 'Writing CX sources'
        do is = 0,ns-1
          write (chns,'(i3.3)') is
          wrk(:,:) = sna0_cx(:,:,0,is)+sna0_cx(:,:,1,is)*na(:,:,is)
          call my_out(70,nx,ny,wrk,'b2stcx_sna_'//chns)
          wrk(:,:) = shi0_cx(:,:,0)+shi0_cx(:,:,1)*ti+
     &               shi0_cx(:,:,2)*ni(:,:,0)+
     &               shi0_cx(:,:,3)*ti*ni(:,:,0)
          call my_out(70,nx,ny,wrk,'b2stcx_shi')
          wrk(:,:) = smo0_cx(:,:,0,is)+
     &               smo0_cx(:,:,1,is)*ua(:,:,is)+
     &               smo0_cx(:,:,2,is)*na(:,:,is)*mp*am(is)+
     &               smo0_cx(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is)
          call my_out(70,nx,ny,wrk,'b2stcx_smq'//chns)
        enddo
        if (allocated(shi0_cx)) then
          write(*,*) 'Deallocating arrays for cx'
          deallocate (shi0_cx)
          deallocate (sna0_cx)
          deallocate (smo0_cx)
        endif
      endif

      if(iout.ne.0) then                                                 !srv 30.06.08 {
        do is = 0,ns-1
          write (chns,'(i3.3)') is
          do k=0,1
            write (chk,'(i1)') k
            call my_out(70,nx,ny,sna(-1,-1,k,is),
     &       'b2sral_sna'//chk//chns)
          enddo
          do k=0,3
            write (chk,'(i1)') k
            call my_out(70,nx,ny,smo(-1,-1,k,is),
     &       'b2sral_smo'//chk//chns)
            call my_out(70,nx,ny,smq(-1,-1,k,is),
     &       'b2sral_smq'//chk//chns)
          enddo
        enddo
        do k=0,3
          write (chk,'(i1)') k
          call my_out(70,nx,ny,she(-1,-1,k),'b2sral_she'//chk)
          call my_out(70,nx,ny,shi(-1,-1,k),'b2sral_shi'//chk)
          call my_out(70,nx,ny,sch(-1,-1,k),'b2sral_sch'//chk)
          if (use_astra.ne.0) then
            call my_out(70,nx,ny,sheNb_astra(-1,-1,k),
     &       'b2sral_sheNb_astra'//chk)
            call my_out(70,nx,ny,shiNb_astra(-1,-1,k),
     &       'b2sral_shiNb_astra'//chk)
          endif
        enddo
        do k=0,1
          write (chk,'(i1)') k
          call my_out(70,nx,ny,sne(-1,-1,k),'b2sral_sne'//chk)
        enddo
      endif                                                              !srv 30.06.08 }

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
