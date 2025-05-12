      module b2mod_driver
      use b2mod_types
      use b2mod_version
      use b2mod_tallies
      use b2mod_constants
      use b2mod_time
      use b2mod_diag
      use b2mod_elements
      use b2us_feedback
     & , only : write_b2us_feedback, init_feedback, dealloc_feedback,
     &          fb_target, fb_prev, fb_current, cum_volrec, fb_const,
     &          charge_frac, saved_fb_actuator, fb_rescale, dt_prev,
     &          feedback_namelist_used
      use b2mod_sputter
     & , only : sput_dst, dealloc_sputter_data, dealloc_b2mod_sputter
      use b2mod_solpstop
      use b2mod_ysmp_sdrv
     & , only : dealloc_b2mod_ysmp_sdrv
      use b2mod_wall
     & , only : dealloc_b2mod_wall
      use b2mod_neoclassical
     & , only : dealloc_b2mod_neoclassical
      use b2mod_facdrift_exb
     & , only : dealloc_b2mod_facdrift_exb
      use b2mod_first_flight
      use b2mod_input_profile
      use b2mod_user_namelist
      use b2mod_ma28_for_us
      use b2mod_transport_models
      use b2mod_neutr_src_scaling
      use b2mod_neutrals_namelist
      use b2mod_transport_namelist
      use b2mod_transport_disruption
     & , only : dealloc_b2mod_transport_disruption
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_b2cmrc
      use b2mod_boundary_namelist
      use b2mod_mwti
      use b2mod_file
      use b2mod_batch_average
     & , only : batch_av_all_init, batch_av_all_save, batch_av_all_fin
      use b2mod_trace
      use b2mod_ipmain
      use b2mod_openmp
      use b2mod_running_average
     & , only : run_av_init, run_av_get_plasma, run_av_save,
     & run_av_fin
#ifndef FIXED_POINT
      use b2mod_balance !djm Jan2017
     & , only : balance_netcdf, alloc_b2mod_balance, write_balance,
     &          dealloc_b2mod_balance, dealloc_b2mod_eirene_sources,
     &          balance_average, update_average_balance
#endif
      use b2mod_equation_sources                                        !srv 26.02.18
      use b2mod_boundary_sources                                        !srv 26.02.18
      use b2mod_switches
      use b2us_map
      use b2us_geo
      use b2us_plasma
      use b2us_io
      use b2mod_par_opt
      use b2mod_b2plot
     , , only : wklng, nlimi, alloc_b2mod_b2plot, dealloc_b2mod_b2plot

#ifdef B25_EIRENE
      use eirmod_parmmod
     * , only : ndx, ndy, nstra, nfl, natm, nmol, nion, npls,
     *          nsrfs, n1st, n2nd, n3rd, nreac, 
     *          eirene_set_parmmod
      use eirmod_cinit
     * , only : eirene_init_cinit, fort_lc
      use eirmod_ccoupl
     , , only : eirene_alloc_ccoupl, eirene_dealloc_ccoupl
      use eirmod_clogau
     , , only : eirene_alloc_clogau
      use eirmod_ctrcei
     , , only : eirene_alloc_ctrcei, eirene_dealloc_ctrcei
      use eirmod_comprt
     , , only : iunout
      use eirmod_braeir
     , , only : eirene_alloc_braeir, eirene_dealloc_braeir
      use eirmod_eirbra
     , , only : eirene_alloc_eirbra, eirene_dealloc_eirbra
      use eirmod_infcop
     , , only : eirene_dealloc_infcop
      use eirmod_extrab25
      use eirmod_wneutrals
     , , only : eirene_wneutrals_init
      use b2mod_eirene
     , , only : eirene_mc_final
#ifdef NAGFOR
      use f90_unix_proc ! IGNORE
#endif
#else
      use b2mod_braeir
      use b2mod_eirbra
#endif
      use b2mod_b2_to_astra
     , , only : write_b2_to_astra
#ifdef JET
      use b2mod_ppout
#endif
#ifdef USE_IMPGYRO
      use b2mod_impgyro !yamoto
#endif
#ifdef IMAS
      use helper
#if AL_MAJOR_VERSION > 4
      use ids_routines  ! IGNORE
     , , only : imas_open, al_build_uri_from_legacy_parameters,
     ,          OPEN_PULSE, STRMAXLEN, MDSPLUS_BACKEND
#else
      use ids_routines  ! IGNORE
     , , only : imas_open_env
#endif
      use ids_routines  ! IGNORE
     , , only : ids_deallocate, ids_get
      use ids_schemas   ! IGNORE
     , , only : ids_edge_profiles, ids_edge_sources, ids_edge_transport,
     ,          ids_radiation, ids_dataset_description, ids_equilibrium
#if ( IMAS_MINOR_VERSION > 8 || IMAS_MAJOR_VERSION > 3 )
      use ids_schemas   ! IGNORE
     , , only : ids_real, ids_real_invalid
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
      use ids_schemas   ! IGNORE
     , , only : ids_summary
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
      use ids_schemas   ! IGNORE
     , , only : ids_numerics
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
      use ids_schemas   ! IGNORE
     , , only : ids_divertors
#endif
      use b2mod_ual
     , , only : put_ids_edge, new_ids_edge, close_ual,
     ,          dealloc_ids_edge, dealloc_batch_edge
      use b2mod_ual_io
     , , only : b25_process_ids, b25_av_ids,
     ,          num_time_slices, num_batch_slices
#endif
      use b2mod_subsys
      use b2mod_dimensions
      implicit none
*   ..common blocks
*   ..local variables
      integer ncall, ntim, mvinc, mvnum, jsep,
     &  isfb, ix1, ix2, ix3, ix4, ix5, iy1, iy2, iy3, iy4, iy5,
     &  nxtl, nxtr
      save isfb,
     &     ix1, ix2, ix3, ix4, ix5, iy1, iy2, iy3, iy4, iy5
      integer, save :: inverse_ua = 0, iout = 0                      !srv 11.01.13 25.01.16
      integer, save :: old_style = 0                                 !pb 29.07.16
      integer, save :: ird_aver = 0
      integer, save :: shot = 0, run = 0, ids_save = 0, ids_av = 0
#ifdef _OPENMP
      real (kind=R8) :: cpuval, cpustart, cpuinit
#else
      real (kind=R4) :: cpuval, cpustart, cpuinit
#endif
      real (kind=R8) ::
     *  dtim, etim,
     &  b2mndr_cpu,
     &  b2mndr_elapsed, elapsedval, elapsedstart, elapsedinit,
     &  density_rescale, ne_wanted,
     &  ne_wanted_time, ne_wanted_next_time, ne_wanted_mod_time, factor,
     &  ixfb, iyfb, ax, ay, scl
      real (kind=R8), allocatable ::
     *  rza0(:,:), rz20(:,:), rpt0(:,:), rpi0(:,:)
      save ne_wanted, ne_wanted_time, ne_wanted_next_time,
     &  ne_wanted_mod_time, ixfb, iyfb
      character idout0*32, idout1*32, savefile*16, savefilestatus*3
      character, save :: savefile_p*14 = 'b2fplasmf.xxxx'
      character chns*3                                               !srv 25.01.16
      integer itim, ierr, edition, b2time, coronal,
     *  tally, ismain, ismain0, natmi, no_solve,
     *  external_species, idx
      integer, save :: edition_p = 0, plasinc = 1, plasnum = 0
      integer, save :: itim_plas = 0, nwrit_plas = 0, ntim_save = 0,
     *  ntim_batch = 500, ibatch_av_all = 0, edition_run = 0,
     *  edition_batch = 0, iav_run = 0, iav_cont = 1, ntim_run = 1000
      logical lwti, lwav, lrav, lrer, continued, ids_replace
      character(22) :: batch_name, aver_name, run_av_name
      real (kind=R8) ::
     *  cputarget, cpuincrement,
     1  na_eps, te_eps, ti_eps, tn_eps, po_eps, ua_eps, kt_eps, zt_eps,
     2  sput_frc, sput_phys,
     3  delta_min, dt_change_inc, delta_max, dt_change_dec,
     4  dt_min, dt_max, max_delta
#ifndef NO_CDF
#     include <netcdf.inc>
      integer ncid
#endif
#ifdef CATALYST
      logical use_catalyst
#endif
#ifdef JET
      real (kind=R8) ::
     *  delta_tran_time, save_tran_time
      integer tran_num
      character*6 tran_file_name
#endif
#ifdef IMAS
      type(ids_equilibrium) :: equilibrium
      type(ids_edge_profiles) :: edge_profiles
      type(ids_edge_profiles) :: batch_profiles
      type(ids_edge_profiles) :: old_edge_profiles
      type(ids_edge_profiles) :: old_batch_profiles
      type(ids_edge_sources) :: edge_sources
      type(ids_edge_sources) :: batch_sources
      type(ids_edge_transport) :: edge_transport
      type(ids_radiation) :: radiation
      type(ids_dataset_description) :: description
      type(ids_dataset_description) :: old_description
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
      type(ids_summary) :: summary
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
      type(ids_numerics) :: numerics
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
      type(ids_divertors) :: divertors
#endif
      integer :: batch_index
      real (IDS_real) :: old_start_time, old_end_time, ids_end_time
#if AL_MAJOR_VERSION < 5
      character treename*24
#endif
      character username*24, database*24, version*24
      character old_imas_version*132, imas_version*132
      character ids_path*256
      logical new_eq_ggd
#endif
      real (kind=R8) ::
     *  delta_cdfmovie_time, save_cdfmovie_time
      real (kind=R8) ::
     *  delta_plasma_time, save_plasma_time
      real (kind=R8) ::
     *  delta_ids_time, save_ids_time
      real (kind=R8) ::
     *  te_hot, ne_hot_frac
      logical quitexist, quit, quitexist_
      logical write_save, write_ids, write_av
      logical process_ids, process_av
      logical, save :: summary_av = .false.
      real (kind=R8) ::
     *  stack_dtim(20)
      integer stack_ntim(20)
      integer stack_ptr
      real (kind=R8) ::
     *  min_areshe, min_areshi, min_aresco, res_quit, res_max
      logical quit_residual, test_residual
      real (kind=R8) :: hzconst                                      !srv 30.01.10

*   ..initialisation
      save ncall, ntim, mvinc, mvnum, dtim, etim, idout0, idout1,
     1  b2mndr_cpu, edition, cputarget, cpuincrement,
     2  ismain, ismain0, savefilestatus, b2time, coronal,
     3  natmi, tally,
     4  density_rescale,
     .  na_eps, te_eps, ti_eps, tn_eps, po_eps, ua_eps, kt_eps, zt_eps,
     5  sput_frc, sput_phys, delta_min, dt_change_inc, delta_max,
     6  dt_change_dec, no_solve, dt_min, dt_max,
     7  te_hot, ne_hot_frac, external_species,
     8  b2mndr_elapsed
#ifdef JET
      save delta_tran_time, save_tran_time, tran_num
#endif
      save delta_cdfmovie_time, save_cdfmovie_time
      save delta_plasma_time, save_plasma_time
      save delta_ids_time, save_ids_time
      save min_areshe, min_areshi, min_aresco, res_quit, res_max
      data ncall/0/, ntim/1/, mvinc/1/, mvnum/0/, dtim/1.0_R8/,
     1     etim/0.0_R8/, ismain/1/, ismain0/0/,
     2     te_hot/0.0_R8/, ne_hot_frac/0.0_R8/
      data idout0/'pgnl;pgmm;pzmm'/, idout1/'pzmm'/
      data b2mndr_cpu/0.0_R8/,
     &     b2mndr_elapsed/0.0_R8/
      data savefile/'plasmastate.xxxx'/,
     1     cpuincrement/3600.0_R8/, savefilestatus/'new'/, b2time/1/,
     2     coronal/0/,
     *     tally/1/,
     1     density_rescale/1.0_R8/, sput_frc/0.0_R8/, sput_phys/0.0_R8/,
     *     no_solve/0/,
     2     na_eps/1.0e19_R8/,
     3     te_eps/1.0e+1_R8/,
     4     ti_eps/1.0e+1_R8/,
     4     tn_eps/1.0e+1_R8/,
     5     po_eps/1.0e+1_R8/,
     5     kt_eps/1.0e-2_R8/,
     5     zt_eps/1.0e-2_R8/,
     6     ua_eps/1.0e+4_R8/,
     7     delta_min/0.0_R8/,
     8     dt_change_inc/1.0_R8/,
     9     delta_max/0.0_R8/,
     1     dt_change_dec/1.0_R8/,
     2     dt_min/1.0e-30_R8/,
     3     dt_max/1.0e+01_R8/,
     5     external_species/0/
#ifdef JET
      data delta_tran_time/0.0_R8/, save_tran_time/0.0_R8/,
     1     tran_num/0/
#endif
      data delta_cdfmovie_time/0.0_R8/, save_cdfmovie_time/0.0_R8/
      data delta_plasma_time/0.0_R8/, save_plasma_time/0.0_R8/
      data delta_ids_time/0.0_R8/, save_ids_time/0.0_R8/
      data min_areshe/0.0_R8/, min_areshi/0.0_R8/, min_aresco/0.0_R8/,
     &     res_quit/1.0e-7_R8/
      data ne_wanted/0.0_R8/, ne_wanted_time/0.0_R8/,
     &     ne_wanted_mod_time/0.0_R8/

      contains

      subroutine b2mndr_0 (ninp, nout, ns, ns0, switch, geo,
     &  mpg, state, state_ext, state_avg)
      use b2mod_ad
     & , only : old_erosion, old_deposition
      implicit none
*   ..input arguments (unchanged on exit)
#ifdef JET
      integer ninp(0:6), nout(0:11)
#else
      integer ninp(0:6), nout(0:10)
#endif
      integer nCv, nFc, nVx, ns, ns0
      type(switches), intent (inout)   :: switch
      type(mapping), intent (inout)    :: mpg
      type(geometry), intent (inout)   :: geo
      type(B2State), intent (inout)    :: state
      type(B2StateExt), intent (inout) :: state_ext
      type(B2Average), intent (inout)  :: state_avg
*   ..output arguments (unspecified on entry)
*     (none)
*   ..local variables
      real (kind=R8) ::
     1  rescale_sa(0:nsdecl-1), rescale_ra(0:nsdecl-1),
     2  rescale_qa(0:nsdecl-1), rescale_cx(0:nsdecl-1),
     3  rescale_rd(0:nsdecl-1), rescale_br(0:nsdecl-1)
      integer atomic_physics_rescale_flag
      integer iCv, is, istra, k, idum(0:9)
#ifdef WG_TODO
      integer ix, iy
#endif
#ifdef B25_EIRENE
      integer iflag
#endif
      character*256 filename
      logical ex, found
#ifdef IMAS
      integer l, m
      character(len=256) :: home_dir, imasdir, zeile, path
      character(len=24) :: device_env
      logical absolute_path, not_default
#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
#endif
*   ..initialisation
      data atomic_physics_rescale_flag/0/
      namelist /atomic_physics_rescale/
     & rescale_sa, rescale_ra, rescale_qa, rescale_cx,
     & rescale_rd, rescale_br
*   ..procedures
      real (kind=R8) :: rratio, smax
      logical streql
      external rratio, smax, streql
#ifdef IMAS
      character*16 usrnam
      external usrnam
      intrinsic index
#endif
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MNDR is the main driver routine for the B2 code.
*
*
*  3. description (see also routine b2cdca)
*
*     For a detailed description of the equations and the method of
*     solution, see the references cited in the main program unit.
*
*     The routine b2mndr controls the timestepping procedure in a B2
*     calculation. An outline of the procedure follows.
*     ..initialise files and common
*     ..do while ()
*      ..perform one timestep
*     ..enddo
*     ..display result of calculation
*     ..save result of calculation
*
*     We have not implemented any adaptive procedure to determine a
*     timestep or to terminate the iteration; the number of iterations
*     and the timestep are predetermined. It is, however, easy to
*     continue a computation starting from a previous final state, in
*     case such a state is judged to be not well converged.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:6) integer array, input.
*     ninp specifies the input unit numbers, as follows.
*     ninp(0): formatted; provides physics and numerics parameters.
*     ninp(1): un*formatted; provides the geometry.
*     ninp(2): un*formatted; provides default physics parameters.
*     ninp(3): un*formatted; provides the initial state.
*     ninp(4): un*formatted; provides atomic rate data.
*     ninp(5): un*formatted; provides source term data.
*     ninp(6): un*formatted; provides the initial differentiated state.
*
*     nout - (0:11) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides summary print output.
*     nout(1): un*formatted; provides actual physics parameters.
*     nout(2): un*formatted; provides the final state.
*     nout(3): un*formatted; may provide data for movie.
*     nout(4): un*formatted; provides data on the numerical evolution.
*     nout(5): un*formatted; may provide detailed convergence data.
*     nout(9): un*formatted; provides the final differentiated state.
*     nout(10): un*formatted; provides data on the numerical evolution of differentiated quantities.
*     nout(11): formatted; postprocessor output. [JET]
*
*     nx, ny - integer, input.
*     nx and ny specify the number of interior cells along the first
*     and the second coordinate, respectively. The total number of
*     cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
*     It will hold that 0.le.nx and 0.le.ny.
*
*     ns - integer, input.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*     ns0 - integer, input.
*     ns0 specifies the number of atomic species in the data file
*     for the initial state.  The initial data for species (0:ns0-1)
*     are taken from that file, while the initial densities for species
*     (ns0:ns-1) are set to some very small value.
*     It will hold that 1.le.ns0.le.ns.
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations


*-----------------------------------------------------------------------
*.documentation-internal
*
*     The workspace array w is used as a stack for storage of dynamic
*     real*8 arrays. An array 'name' in w(0:*) will be stored starting at
*     location w(l'name').
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     /b2cmpa/ specifies basic physics parameters.
*      me, mp, ev, qe, eps0, mu0 - real*8 constant.
*     These parameters are included with /b2cmpa/.
*       me = electron mass (kg)
*       mp = proton mass (kg)
*       ev = 1 electron-Volt (Joule)
*       qe = elementary charge (Coulomb)
*       eps0 = vacuum dielectric permittivity (SI)
*       mu0 = vacuum magnetic permeability (SI)
*      nsdecl - integer constant.
*     nsdecl is an upper bound for the number of species (ns), and is
*     used in dimension specifications.
*      zamin, zamax, zn, am - (0:nsdecl-1) real*8 array.
*     For (is) in (0:ns-1), zn(is) specifies the nuclear charge of
*     species (is), am(is) specifies its atomic mass, and zamin(is) and
*     zamax(is) specify its atomic charge range. All are dimensionless.
*     It will hold that 0.le.zamin(is).le.zamax(is), 0.lt.zn(is) and
*     0.lt.am(is). It is required that neutral species not be bundled
*     with charged species.
*     For ns.le.is.lt.nsdecl, zamin(is), zamax(is), zn(is), am(is)
*     undefined.
*
*   ..transport coefficients, parameters
*
*     /b2cmpt/ - common block.
*     (Note: the integer constant nsdecl, which comes with /b2cmpa/,
*     is used instead of (ns) in dimension specifications.)
*      cfdf0 - (0:7,0:ns-1) real*8 array.
*     cfdf0 contains some auxiliary coefficients. For 0.le.is.lt.ns,
*     the local expression
*       df0(,is) = cfdf0(0,is)+cfdf0(1,is)/(ni(,0)/1.0e20)+
*         cfdf0(2,is)*sqrt(ti()/(am(is)*mp))/
*          (cfdf0(3,is)*ni(,0)+cfdf0(4,is)*ne())
*     is used in the definition of the anomalous transport coefficients.
*     df0 is meant to represent a diffusivity. If cfdf0(2,is).eq.0 then
*     the corresponding term in df0 is understood as 0 independent of
*     cfdf0(3:4,is).
*      cfdna - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfdna(,is) are the coefficients in the
*     local expression
*       difna(,is) = cfdna(0,is)+cfdna(1,is)/(ne()/1.0e20)+
*         cfdna(2,is)*df0(,is)+
*         cfdna(3,is)*te()/ev/abs(bb(,3))/16
*     that specifies the anomalous diffusivity with respect to
*     grad.na; the transport equation for species (is) contains a term
*     difna(,is)*'grad.na'.
*     If cfdna(7,is) is not zero, the dna0 array is flux-scaled with
*     a scaling factor of cfdna(7,is)*(hy1()/hy1(omp))**2
*     where hy1 is the radial width of the cell and omp refers
*     to the cell from the sam flux tube intersecting the outer midplane.
*     cfdna(4:6,) are not used in the default version of the code.
*      cfdpa - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfdpa(,is) are the coefficients in the
*     local expression
*       difpa(,is) = cfdpa(0,is)+cfdpa(1,is)/(ne()/1.0e20)+
*         cfdpa(2,is)*df0(,is)+
*         cfdpa(3,is)*te()/ev/abs(bb(,3))/16
*     that specifies the anomalous diffusivity with respect to
*     grad.pa; the transport equation for species (is) contains a term
*     (1/(rza(,is)*te()+ti()))*difpa(,is)*'grad.pa'.
*     If cfdpa(7,is) is not zero, the dpa0 array is flux-scaled with
*     a scaling factor of cfdpa(7,is)*(hy1()/hy1(omp))**2
*     where hy1 is the radial width of the cell and omp refers
*     to the cell from the sam flux tube intersecting the outer midplane.
*     cfdpa(4:6,) are not used in the default version of the code.
*      cfvla - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfvla(,is) are the coefficients in the
*     local expression
*       vlay(,is) = cfvla(0,is)+cfvla(1,is)/(ne()/1.0e20)+
*         cfvla(2,is)*df0(,is)+
*         cfvla(3,is)*te()/ev/abs(bb(,3))/16
*     that specifies the anomalous y-velocity of species (is); the
*     transport equation for species (is) contains a term
*     vlay(,is)*na(,is).
*     If cfvla(7,is) is not zero, the vla0 array is flux-scaled with
*     a scaling factor of cfvla(7,is)*(hy1()/hy1(omp))**2
*     where hy1 is the radial width of the cell and omp refers
*     to the cell from the sam flux tube intersecting the outer midplane.
*     cfvla(4:6,) are not used in the default version of the code.
*      cfvsa - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfvsa(,is) are the coefficients in the
*     local expression
*       dfvsa(,is) = cfvsa(0,is)+cfvsa(1,is)/(ne()/1.0e20)+
*         cfvsa(2,is)*df0(,is)+
*         cfvsa(3,is)*te()/ev/abs(bb(,3))/16
*     that specifies the diffusivity corresponding to the anomalous
*     transverse viscosity, so that this viscosity is
*       vsa(,is) = mp*am(is)*na(,is)*dfvsa
*     If cfvsa(7,is) is not zero, the vsa0 array is flux-scaled with
*     a scaling factor of cfvsa(7,is)*(hy1()/hy1(omp))**2
*     where hy1 is the radial width of the cell and omp refers
*     to the cell from the sam flux tube intersecting the outer midplane.
*     cfvsa(4:6,) are not used in the default version of the code.
*      cfhci - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfhci(,is) are the coefficients in the
*     local expression
*       dfhca(,is) = cfhci(0,is)+cfhci(1,is)/(ne()/1.0e20)+
*         cfhci(2,is)*df0(,is)+
*         cfhci(3,is)*te()/ev/abs(bb(,3))/16
*     that specifies the anomalous atom heat conductivity associated
*     with species (is), so that the atom heat conduction coefficient
*     contains a term dfhca(,is)*na(,is).
*     If cfhci(7,is) is not zero, the hci arrays are flux-scaled with
*     a scaling factor of cfhci(7,is)*(hy1()/hy1(omp))**2
*     where hy1 is the radial width of the cell and omp refers
*     to the cell from the sam flux tube intersecting the outer midplane.
*     cfhci(4:6,) are not used in the default version of the code.
*      cfhce - (0:7) real*8 array.
*     cfhce() are the coefficients in the local expression
*       dfhce() = cfhce(0)+cfhce(1)/(ne()/1.0e20)+
*         cfhce(3)*te()/ev/abs(bb(,3))/16
*     that specifies the anomalous electron heat conductivity, so that
*     the electron heat conduction coefficient is dfhce()*ne().
*     If cfhce(7) is not zero, the hce0 array is flux-scaled with
*     a scaling factor of cfhce(7)*(hy1()/hy1(omp))**2
*     where hy1 is the radial width of the cell and omp refers
*     to the cell from the sam flux tube intersecting the outer midplane.
*     cfhce(2) and (4:6) are not used in the default version of the code.
*      cfsig - (0:7) real*8 array.
*     cfsig() are the coefficients in the local expression
*       dfsig() = cfsig(0)+cfsig(1)/(ne()/1.0e20)+
*         cfsig(3)*te()/ev/abs(bb(,3))/16+
*         cfsig(4)*exp(-cfsig(5)*(y/L)**2)
*     that specifies the diffusivity corresponding to the anomalous
*     electrical conductivity, which is
*       sig() = qe*ne()*dfsig().
*     L is the local radial width of the plasma and
*     y is the radial distance from the inner (South) boundary.
*     If cfsig(7) is not zero, the sig0 array is flux-scaled with
*     a scaling factor of cfsig(7)*(hy1()/hy1(omp))**2
*     where hy1 is the radial width of the cell and omp refers
*     to the cell from the sam flux tube intersecting the outer midplane.
*     cfsig(2) and (6) are not used in the default version of the code.
*      cfalf - (0:7) real*8 array.
*     cfalf() are the coefficients in the local expression
*       dfalf() = cfalf(0)+cfalf(1)/(ne()/1.0e20)+
*         cfalf(3)*te()/ev/abs(bb(,3))/16
*     that specifies the diffusivity corresponding to the anomalous
*     thermo-electric coefficient, which is
*       alf() = ne()*sqrt(qe/te())*dfalf().
*     If cfalf(7) is not zero, the alf0 array is flux-scaled with
*     a scaling factor of cfalf(7)*(hy1()/hy1(omp))**2
*     where hy1 is the radial width of the cell and omp refers
*     to the cell from the sam flux tube intersecting the outer midplane.
*     cfalf(2) and (4:6) are not used in the default version of the code.
*      cflim - (0:7) real*8 array.
*     cflim contains coefficients for flux limiting. cflim(0)=cflme
*     specifies the flux limit factor for poloidal electron heat
*     transport. cflim(1)=cflmi specifies the flux limit factor for
*     poloidal all atom heat transport. cflim(2)=cflmv specifies the
*     flux limit factor for poloidal atom momentum transport.
*     cflim(3)=cflal specifies the flux limit factor for the
*     thermo-electric coefficient.
*     cflim(4)=cflab specifies the flux limit factor for the
*     friction force.
*     Zero settings turn off these flux limits.
*     cflim(5:7) are not used in the default version of the code.
*
*     /b2cmrc/ contains a table of atomic rate coefficients.
*      rtnt - integer constant.
*     rtnt specifies the number of intervals in the discretisation
*     of the temperature range. The table will contain data for
*     (rtnt+1) values of the temperature variable.
*      rtnn - integer constant.
*     rtnn specifies the number of intervals in the discretisation
*     of the density range. The table will contain data for (rtnn+1)
*     values of the density variable.
*      rtns - integer.
*     rtns specifies the total number of atomic species; i.e.
*       (sum inuc : 0.le.inuc.lt.nnuc : izhi(inuc)-izlo(inuc)+1) ,
*     where inuc is a nuclear species index, izlo(inuc) is the lowest
*     and izhi(inuc) is the highest atomic charge considered in B2 for
*     that nuclear species.
*     rtns equals the variable ns used elsewhere in B2.
*      rtdumm - integer.
*     rtdumm is a dummy integer variable, introduced in order to obtain
*     an even-length sequence of integer variables in the common block.
*      rtzmin, rtzmax, rtzn - (0:rtns-1) real*8 array.
*     For 0.le.is.lt.rtns, rtzmin(is) and rtzmax(is) specify the atomic
*     charge range, and rtzn(is) specifies the nuclear charge
*     of atomic species (is).
*     The izhi(inuc)-izlo(inuc)+1 charge states for each nuclear
*     species are taken in natural order, starting at izlo(inuc); the
*     highest state of nuclear species (inuc) is followed by the lowest
*     state of (inuc+1). rtzmin, rtzmax and rtzn are equal to the arrays
*     zamin, zamax and zn used elsewhere in B2.
*      rtt - (0:rtnt) real*8 array.
*     rtt(0:rtnt) specifies the sequence of values of the temperature
*     variable, expressed in ev, for which the rate coefficients are
*     evaluated. rtt(0:rtnt) will be positive and increasing.
*      rtn - (0:rtnn) real*8 array.
*     rtn(0:rtnn) specifies the sequence of values of the density
*     variable, expressed in m**-3, for which the rate coefficients
*     are evaluated. rtn(0:rtnn) will be positive and increasing.
*      rtlt - (0:rtnt) real*8 array.
*     rtlt(0:rtnt) holds log(rtt(0:rtnt)).
*      rtln - (0:rtnn) real*8 array.
*     rtln(0:rtnn) holds log(rtn(0:rtnn)).
*      rtlsa - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlsa(i,j,k) specifies
*     logu(rtrsa), where rtrsa is the ionisation rate coefficient at
*     electron temperature rtt(i), electron density rtn(j), for the
*     process k->k+1. rtrsa has dimension m**3/s. The function logu is
*     the same as (log) except that logu(0.0) is equal to some huge
*     negative number.
*      rtlra - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlra(i,j,k) specifies
*     logu(rtrra), where rtrra is the recombination rate coefficient at
*     electron temperature rtt(i) and electron density rtn(j), for the
*     process k->k-1. rtrra has dimension m**3/s.
*      rtlqa - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlqa(i,j,k) specifies the
*     logu(rtrqa), where rtrqa is the electron heat loss rate coefficient
*     at electron temperature rtt(i) and electron density rtn(j) for
*     processes starting from charge state k, excluding the contributions
*     coming from potential energy of recombination and other recombination
*     process contributions. In short, it includes the processes not
*     accounted for in Eirene.
*     rtqa has dimension eV*m**3/s.
*      rtlqr - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtqr(i,j,k) specifies the
*     logu(rtrqr), where rtrqr is the electron heat loss rate
*     coefficient at electron temperature rtt(i)
*     and electron density rtn(j) for processes starting from charge
*     state k. rtqr only includes the contribution of recombination
*     which replaces AMJUEL when Eirene is not used. rtqr has
*     dimension eV*m**3/s.
*      rtlcx - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlcx(i,j,k) specifies
*     logu(rtrcx), where rtrcx is the charge exchange rate coefficient
*     at neutral hydrogen temperature/amu rtt(i) and electron density
*     rtn(j), for the process k->k-1. rtrcx has dimension m**3/s.
*      rtlza - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlza(i,j,k) specifies
*     the effective charge state at
*     electron temperature rtt(i), electron density rtn(j), for the
*     superstage k. rtlza has dimension e.
*      rtlz2 - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlz2(i,j,k) specifies
*     the effective squared charge at
*     electron temperature rtt(i), electron density rtn(j), for the
*     superstage k. rtlz2 has dimension e**2.
*      rtlpt - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlpt(i,j,k) specifies
*     the cumulative ionisation potential at
*     electron temperature rtt(i), electron density rtn(j),
*     to reach superstage k. rtlpt has dimension eV.
*      rtlpi - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlpi(i,j,k) specifies
*     the effective ionisation potential at
*     electron temperature rtt(i), electron density rtn(j),
*     to reach superstage k from superstage k-1. rtlpi has dimension eV.
*
*     /b2cmsr/ specifies source terms from neutral gas model.
*     (Presently empty--to be used in conjunction with a Monte Carlo
*     neutral gas code.)
*
*     Description of some local variables (including arrays in w(0:*)):
*
*     cvX, cvY - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), cvX(iCv) specifies the
*     Cartesian X-coordinate of the center of the cell and cvX(iCv) specifies
*     the Cartesian Y-coordinate of the center of the cell.
*
*     vxX, vxY - (1:nVx) real*8 array.
*     vxX(iVx) specify the horizontal coordinate and vxY(iVx)
*     specify the vertical coordinates (in the physical plane) of
*     the vertex with number iVx.
*
*     cvBb - (1:nCv,0:3) real*8 array.
*     For (iCv) in (1:nCv), cvBb(iCv,0:3) specifies the
*     magnetic field at the center of the (iCv) cell.
*     cvBb(iCv,0:2) are the (x,y,z)-components and cvBb(iCv,3) is
*     the absolute magnetic field strength.
*     (cvBb(,3)=sqrt(cvBb(,0)**2+cvBb(,1)**2+cvBb(,,2)**2.)
*     It will hold that 0.le.cvBb(,0:3), because the direction of the
*     magnetic field follows from cvEb (see below).
*
*     cvEb - (1:nCv,0:2) real*8 array.
*     cvEb(iCv,0:2) specifies the direction of the magnetic field at
*     the center of the cell iCv in the (R,Z,toroidal) axis system.
*     (sqrt(cvEb(iCv,0)**2+cvEb(iCv,1)**2+cvEb(iCv,2)**2)=1)
*
*     cvVol - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), cvVol(iCv) specifies the volume of
*     the (iCv) cell. (In the case of toroidal symmetry this volume
*     includes a factor 2*pi*R, while in the case of cylindrical
*     symmetry it includes a factor 2*pi*Z.)
*     It will hold that 0.lt.cvVol().
*
*     cvHz - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), cvHz(iCv) specifies the length of
*     the (iCv) cell in the z-direction, measured at the cell center.
*
*     cvSz - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), cvSz(iCv) specifies the area of
*     the (iCv) cell perpendicular to the z-direction.
*
*     fcQgam - (1:nFc,0:1) real*8 array.
*     For iFc in (1:nFc), fcQgam(iFc,0) specifies cos(t) and
*     fcQgam(iFc,1) specifies sin(t), where t is the included angle
*     between the cell-connector vector and the normal vector
*     of face iFc. It will hold that 0.lt.fcQgam(,).le.1.
*
*     fcQalf - (1:nFc,0:1) real*8 array.
*     For iFc in (1:nFc), fcQalf(iFc,0) specifies cos(t) and
*     fcQalf(iFc,1) specifies sin(t), where t is the included angle
*     between the magnetic field vector and the face normal.
*     If Bz is negative, fcQalf(iFc,1) is modified to be -sin(t).
*
*     fcQbet - (1:nFc,0:1) real*8 array.
*     For iFc in (1:nFc), fcQbet(iFc,0) specifies cos(t) and
*     fcQbet(iFc,1) specifies sin(t), where t is the included angle
*     between the cell-connector and the magnetic field vector.
*     If Bz is negative, fcQbet(iFc,0) is modified to be -cos(t).
*
*     fcPbs - (1:nFc) real*8 array.
*     For iFc in (1:nFc), fcPbs(iFc) specifies the product
*     (fcBx/fcBb)*sx or (magnetic field pitch)*(surface area),
*     namely the parallel contact area.
*
*     fcPbshz - (1:nFc) real*8 array.
*     For iFc in (1:nFc), fcPbshz(iFc) specifies the product
*     (fcBx/fcBb)*sx*hz or (magnetic field pitch)*(surface area)*(hz).
*
*     fcS - (1:nFc,1:2) real*8 array.
*     For iFc in (1:nFc), fcHc(iFc,:) specifies the total area of face iFc.
*
*     na - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), na(iCv,is) specifies
*     the density of atomic species (is) on the (iCv) cell.
*     It will hold that 0.lt.na(:).
*
*     ua - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), ua(iCv,is) specifies
*     the parallel velocity of atomic species (is) on the (iCv) cell.
*     The sign of the parallel velocity indicates whether it is flowing
*     in the same direction as the local poloidal magnetic field.
*
*     te - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), te(iCv) specifies the electron
*     temperature on the (iCv) cell. It will hold that 0.lt.te(:).
*
*     ti - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), ti(iCv) specifies the all atom
*     temperature on the (iCv) cell. It will hold that 0.lt.ti(:).
*
*     tn - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), tn(iCv) specifies temperature of the
*     hydrogenic neutrals on the (iCv) cell. It will hold that 0.lt.tn(:).
*     When no separate fluid neutral energy equation is solved, tn = ti.
*
*     po - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), po(iCv) specifies the electric
*     potential on the (iCv) cell.
*
*     kt - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), kt(iCv) specifies the turbulent kinetic
*     energy on the (iCv) cell.
*
*     zt - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), kt(iCv) specifies the enstrophy
*     on the (iCv) cell.
*
*     For all the fluxes below, the sign convention is as follows.
*     Fluxes are positive when they flow in the local x-direction,
*     where x is the connector line from C1 to C2.
*
*     fna - (1:nFc,0:1,0:ns-1) real*8 array.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), fna(iFc,0,is) specifies
*     the poloidal component and fna(iFc,1,is) specifies
*     the radial component of the flux of atoms of species (is) through the
*     face iFc.
*
*     fhe - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fhe(iFc,0) specifies the poloidal component and
*     fhe(iFc,1) specifies the radial component of the electron heat flux
*     through the face iFc.
*
*     fhi - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fhi(iFc,0) specifies the poloidal component and
*     fhi(iFc,1) specifies the radial component of the all atom heat flux
*     through the face iFc.
*
*     fch - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fch(iFc,0) specifies the poloidal component and
*     fch(iFc,1) specifies the radial component of the electric
*     current through the face iFc.
*
*     uadia - (1:nFc,0:1,0:ns-1) real*8 array.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), uadia(iFc,0,is) specifies
*     the total effective drift velocity of species (is) on the (iFc) face
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but projected onto the
*     poloidal direction.
*     For (iFc,is) in (1:nFc,0:1,0:ns-1), uadia(iFc,1,is) specifies
*     the total effective drift velocity of species (is) on the (iFc) face
*     in the radial direction.
*
*     na0, ua0, te0, ti0, po0, fna0, fhe0, fhi0, fch0 - real*8 array.
*     state of the plasma at the start of the present timestep.
*
*     nap, uap, tep, tip, pop, fnap, fhep, fhip, fchp - real*8 array.
*     state of the plasma at the start of the previous timestep.
*
*     ntim - integer.
*     ntim specifies the number of timesteps to be performed.
*
*     itim - integer.
*     itim specifies the timestep counter. It is initialised to 0 and
*     incremented by 1 following each successful call to b2mndt.
*
*     dtim - real.
*     dtim specifies the timestep for the present call to b2mndt.
*     It will hold that 0.lt.dtim.
*
*     stim - real.
*     stim specifies the initial time --- default 0.
*     If set to a positive or zero value, this overwrites the time
*     value read from b2fstati.
*     If set to a negative value, the run continues from the time
*     read in b2fstati.
*
*     etim - real.
*     etim specifies the end time --- default 0.
*     Only active if etim > stim
*
*     tim - real.
*     tim is the simulation time. It will hold that 0.le.tim.
*
*     idout0, idout1 - character*32.
*     idout0 and idout1 specify the desired selection of output
*     segments; idout0 for output before the calculation and idout1
*     for output at the conclusion of the calculation. The generic
*     description of either of these 'idout' variables follows.
*     idout specifies the desired selection of output segments for
*     printed and graphical output. Associated with each segment is a
*     four-character word, which is identified on the printed output.
*     (The same word may be associated with several segments.)
*     A particular segment will be produced only if its associated word
*     has a match, as defined by the routine strmas, in idout.
*
*     mvinc - integer.
*     mvinc specifies the number of timesteps between movie output.
*
*     mvnum - integer.
*     mvnum specifies the maximum number of instances at which movie
*     data will be output.
*
*     ismain - integer.
*     ismain identifies the index of the main plasma species.
*     This cannot be a neutral species.
*
*     ismain0 - integer.
*     ismain0 identifies the index of the main neutral species.
*     This can be equal to the main plasma species index only for
*     cases where the plasma solution does not follow neutrals.
*     Otherwise, cannot be a charged species.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mndr_0')
*   ..test input data
*     (done in b2mnds)
      nCv = mpg%nCv
      nFc = mpg%nFc
      nVx = mpg%nVx
*   ..initialise files and common
      call b2mnds (ninp, nout, nCv, nFc, ns, ns0, switch)
*   ..set internal parameters on first (and only) call
*     (Note that the internal parameters are only defined following
*     the call to b2mnds.)
      if (ncall.eq.0) then
#ifndef FIXED_POINT
       balance_netcdf  = 0 ! Default: do not create balance arrays
       balance_average = 0 ! Default: do not perform averages of balance arrays
       call ipgeti ('balance_netcdf', balance_netcdf)                   !djm Jan2017
       call ipgeti ('balance_average', balance_average)
       if (balance_netcdf.ne.0) call alloc_b2mod_balance(nCv,nFc,ns)
#endif
       if(switch%boundary_sources.ne.0)
     &    call alloc_b2mod_boundary_sources(nCv,ns)
       if(switch%equation_sources.ne.0)
     &    call alloc_b2mod_equation_sources(nCv,ns)                   !srv 26.02.18 }
       call ipgeti ('b2mndr_ntim', ntim)
       call ipgetr ('b2mndr_dtim', dtim)
       call xertst (0.le.ntim.and.0.0_R8.lt.dtim,
     &   'faulty internal parameter ntim, dtim')
       call ipgeti ('b2mndr_inverse_ua', inverse_ua)                    !srv 11.01.13
       call ipgetr ('b2mndr_etim', etim)
       call ipgeti ('b2news_no_solve', no_solve)
       call xertst (-15.le.no_solve.and.no_solve.le.1,
     &   'faulty input no_solve')
       if (no_solve.eq.1)
     >  write(6,*) 'No_solve switch activated. This is a dummy run!'
       call ipgeti ('b2mndr_ismain', ismain)
       call xertst (0.le.ismain.and.ismain.lt.ns,
     .  'invalid main plasma species index ismain')
       call xertst (.not.is_neutral(ismain),
     .  'invalid main plasma species ismain; must not be neutral')
       found = .false.
       is = ismain
       do while (is.ge.0 .and. .not.found)
         if (is_neutral(is) .and. nint(zn(is)).eq.nint(zn(ismain)) .and.
     &       nint(am(is)).eq.nint(am(ismain))) then
           ismain0 = is
           found = .true.
         end if
         is = is - 1
       end do
       if (.not.found.and.ismain.ne.1) ismain0 = ismain
       call ipgeti ('b2mwti_ismain0', ismain0)
       call xertst (0.le.ismain0.and.ismain0.lt.ns,
     .   'invalid main neutral species index ismain0')
       call xertst (is_neutral(ismain0).or.ismain.eq.ismain0,
     .   'invalid main neutral species index ismain0')
       call ipgetr ('b2mndr_min_areshe', min_areshe)                     !jwk
       call xertst (min_areshe.ge.0.0_R8,                                !jwk
     &      'faulty internal parameter min_areshe')                      !jwk
       call ipgetr ('b2mndr_min_areshi', min_areshi)                     !jwk
       call xertst (min_areshi.ge.0.0_R8,                                !jwk
     &      'faulty internal parameter min_areshi')                      !jwk
       call ipgetr ('b2mndr_min_aresco', min_aresco)                     !jwk
       call xertst (min_aresco.ge.0.0_R8,                                !jwk
     &      'faulty internal parameter min_aresco')                      !jwk
       if (min_areshe.gt.0.0_R8.or.min_areshi.gt.0.0_R8.or.              !jwk
     &     min_aresco.gt.0.0_R8) then                                    !jwk
          test_residual=.true.                                           !jwk
       else                                                              !jwk
          test_residual=.false.                                          !jwk
       endif                                                             !jwk
       call ipgetr ('b2mndr_res_quit', res_quit)                         !sc
       call xertst (res_quit.ge.0.0_R8,                                  !sc
     &      'faulty internal parameter res_quit')                        !sc
       call ipgeti ('b2mndr_plasinc', plasinc)
       call ipgeti ('b2mndr_plasnum', plasnum)
       call ipgeti ('b2mndr_mvinc', mvinc)
       call ipgeti ('b2mndr_mvnum', mvnum)
       call xertst (1.le.mvinc.and.0.le.mvnum,
     &   'faulty internal parameter mvinc, mvnum')
       call xertst (1.le.plasinc.and.0.le.plasnum,
     &   'faulty internal parameter plasinc,plasnum')
       call ipgetc ('b2mndr_idout0', idout0)
       call ipgetc ('b2mndr_idout1', idout1)
       call ipgetr ('b2mndr_cpu', b2mndr_cpu)
       call xertst (b2mndr_cpu.ge.0.0_R8,
     1  'faulty internal parameter cpu')
       call ipgetr ('b2mndr_elapsed', b2mndr_elapsed)
       call xertst (b2mndr_elapsed.ge.0.0_R8,
     1  'faulty internal parameter elapsed')
       call ipgetr ('b2mndr_savecpu', cpuincrement)
       call xertst (cpuincrement.ge.0.0_R8,
     1  'faulty internal parameter cpuincrement')
       cputarget=cpuincrement

#ifdef CATALYST
       use_catalyst = .false.
       filename='coproc.py'
       call find_file(filename,use_catalyst)
       if (use_catalyst .eqv. .false.)
     >  write(*,*) 'File coproc.py missing ',
     > 'from run directory - Catalyst will not start.'
       if (use_catalyst .eqv. .true.) then
          call coprocessorinitializewithpython
     >     (filename,len_trim(filename))
       endif
#endif

#ifdef JET
       call ipgetr ('b2mndr_trantim', delta_tran_time)
       if (delta_tran_time.gt.0.0_R8)
     >  write(*,*) 'Saving tran files every ',delta_tran_time,' seconds'
       call xertst (delta_tran_time.eq.0.0_R8.or.
     &              delta_tran_time.ge.dtim,
     &  'b2mndr_trantim too short compared to timestep !')
#endif
       call ipgetr ('b2mndr_cdfmovietim', delta_cdfmovie_time)
       if (delta_cdfmovie_time.gt.0.0_R8)
     >  write(*,*) 'Saving CDF movie files every ',
     1  delta_cdfmovie_time,' simulation seconds'
       call xertst (delta_cdfmovie_time.eq.0.0_R8.or.
     &              delta_cdfmovie_time.ge.dtim,
     &  'b2mndr_cdfmovietim too short compared to timestep !')
       call ipgetr ('b2mndr_plasmatim', delta_plasma_time)
       if (delta_plasma_time.gt.0.0_R8)
     >  write(*,*) 'Saving plasma state files every ',
     1  delta_plasma_time,' simulation seconds'
       call xertst (delta_plasma_time.eq.0.0_R8.or.
     &              delta_plasma_time.ge.dtim,
     &  'b2mndr_plasmatim too short compared to timestep !')
       call ipgetr ('b2mndr_ids_time', delta_ids_time)
       if (delta_ids_time.gt.0.0_R8)
#ifdef IMAS
     >  write(*,*) 'Writing IDS every ',
     1  delta_ids_time,' simulation seconds'
       call xertst (delta_ids_time.eq.0.0_R8.or.
     &              delta_ids_time.ge.dtim,
     &  'b2mndr_ids_time too short compared to timestep !')
#else
     >  write(*,*)
     1   'Not compiled with IMAS option, IDS will not be written.'
#endif
       call ipgeti ('b2mndr_ntim_save', ntim_save)
       call ipgeti ('b2mndt_av', iav_run)
       call ipgeti ('b2mndt_av_continue', iav_cont)
       call ipgeti ('b2mndt_av_ntim_batch', ntim_batch)
       call ipgeti ('b2mndt_av_ntim_run', ntim_run)
       call ipgeti ('b2mndt_av_batch_all', ibatch_av_all)
       if (ntim_save.gt.0)
     >  write(*,*) 'Saving plasma state files every ',
     1  ntim_save,' iterations'
       call ipgeti ('b2mndr_b2time', b2time)
       call xertst (0.le.b2time,
     &   'faulty internal parameter b2time')
       call ipgetr ('b2mndr_density_rescale', density_rescale)
csw dec2013 allow for negative values; this is tweak to scale a specific density, e.g. Helium,
c   see below --> we need a better solution for this
c       call xertst (0.0_R8.lt.density_rescale,
c     &   'faulty internal parameter density_rescale')
       call ipgeti ('b2mndr_tally', tally)
       call xertst (0.le.tally,
     &   'faulty internal parameter tally')
       call ipgetr ('b2mndr_na_eps',na_eps)
       call xertst (0.0_R8.lt.na_eps,'faulty internal parameter na_eps')
       call ipgetr ('b2mndr_te_eps',te_eps)
       call xertst (0.0_R8.lt.te_eps,'faulty internal parameter te_eps')
       call ipgetr ('b2mndr_ti_eps',ti_eps)
       call xertst (0.0_R8.lt.ti_eps,'faulty internal parameter ti_eps')
       call ipgetr ('b2mndr_tn_eps',tn_eps)
       call xertst (0.0_R8.lt.tn_eps,'faulty internal parameter tn_eps')
       call ipgetr ('b2mndr_po_eps',po_eps)
       call xertst (0.0_R8.lt.po_eps,'faulty internal parameter po_eps')
       call ipgetr ('b2mndr_ua_eps',ua_eps)
       call xertst (0.0_R8.lt.ua_eps,'faulty internal parameter ua_eps')
       call ipgetr ('b2mndr_kt_eps',kt_eps)
       call xertst (0.0_R8.lt.kt_eps,'faulty internal parameter kt_eps')
       call ipgetr ('b2mndr_zt_eps',zt_eps)
       call xertst (0.0_R8.lt.zt_eps,'faulty internal parameter zt_eps')
       call ipgetr ('b2mndr_delta_min',delta_min)
       call xertst (0.0_R8.le.delta_min,
     1  'faulty internal parameter delta_min')
       call ipgetr ('b2mndr_dt_change_inc',dt_change_inc)
       call xertst (0.0_R8.le.dt_change_inc,
     1  'faulty internal parameter dt_change_inc')
       call ipgetr ('b2mndr_delta_max',delta_max)
       call xertst (0.0_R8.le.delta_max,
     1  'faulty internal parameter delta_max')
       call ipgetr ('b2mndr_dt_change_dec',dt_change_dec)
       call xertst (0.0_R8.le.dt_change_dec,
     1  'faulty internal parameter dt_change_dec')
       call ipgetr ('b2mndr_dt_min',dt_min)
       call xertst (0.0_R8.le.dt_min,'faulty internal parameter dt_min')
       call ipgetr ('b2mndr_dt_max',dt_max)
       call xertst (0.0_R8.le.dt_max,'faulty internal parameter dt_max')
* modify rate coefficients for coronal model
       allocate(
     *  rza0(nCv,0:ns-1), rz20(nCv,0:ns-1),
     *  rpt0(nCv,0:ns-1), rpi0(nCv,0:ns-1))
       call ipgeti ('b2mndr_coronal_model', coronal)
       if(coronal.ne.0) then
          write(*,*)
     1     'Modifying rate coefficients to suppress density dependence'
          do is=1,rtnn
             rtlsa(:,is,:)=rtlsa(:,0,:)
             rtlra(:,is,:)=rtlra(:,0,:)
             rtlqa(:,is,:)=rtlqa(:,0,:)
             rtlqr(:,is,:)=rtlqr(:,0,:)                        !mb 22.05.18
             rtlcx(:,is,:)=rtlcx(:,0,:)
             rtlrd(:,is,:)=rtlrd(:,0,:)
             rtlbr(:,is,:)=rtlbr(:,0,:)
             rtlza(:,is,:)=rtlza(:,0,:)
             rtlz2(:,is,:)=rtlz2(:,0,:)
             rtlpt(:,is,:)=rtlpt(:,0,:)
             rtlpi(:,is,:)=rtlpi(:,0,:)
          enddo
       endif
       call ipgeti ('b2mndr_atomic_physics_rescale',
     1      atomic_physics_rescale_flag)
* ..rescale atomic physics coefficients if desired
       if(atomic_physics_rescale_flag.ne.0) then
          write(*,*) 'b2mndr_atomic_physics_rescale',
     1      atomic_physics_rescale_flag
          rescale_sa=1.0_R8
          rescale_ra=1.0_R8
          rescale_qa=1.0_R8
          rescale_cx=1.0_R8
          rescale_rd=1.0_R8
          rescale_br=1.0_R8
          filename='b2.atomic_physics_rescale.parameters'
          call find_file(filename,ex)
          open(99,file=filename)
          write(*,*) 'Reading ',trim(filename)
#ifdef USE_MPI
          read(99,atomic_physics_rescale,err=999)
#else
          read(99,atomic_physics_rescale)
#endif
          write(*,atomic_physics_rescale)
          close(99)
          write(*,*) 'Read ',trim(filename)
          call xertst (min(
     1         minval(rescale_sa(0:rtns-1)),
     2         minval(rescale_ra(0:rtns-1)),
     3         minval(rescale_qa(0:rtns-1)),
     4         minval(rescale_cx(0:rtns-1)),
     5         minval(rescale_rd(0:rtns-1)),
     6         minval(rescale_br(0:rtns-1)))
     7         .gt.0.0_R8,
     8    'faulty value set in b2.atomic_physics_rescale.parameters')
          do is = 0, rtns-1
             rtlsa(:,:,is)=rtlsa(:,:,is)+log(rescale_sa(is))
             rtlra(:,:,is)=rtlra(:,:,is)+log(rescale_ra(is))
             rtlqa(:,:,is)=rtlqa(:,:,is)+log(rescale_qa(is))
             rtlqr(:,:,is)=rtlqr(:,:,is)+log(rescale_qa(is))
             rtlcx(:,:,is)=rtlcx(:,:,is)+log(rescale_cx(is))
             rtlrd(:,:,is)=rtlrd(:,:,is)+log(rescale_rd(is))
             rtlbr(:,:,is)=rtlbr(:,:,is)+log(rescale_br(is))
          enddo
       endif
       call ipgeti ('b2stbr_sput_dst', sput_dst)
       call ipgetr ('b2stbr_sput_frc', sput_frc)
       call ipgetr ('b2stbr_sput_phys', sput_phys)
       call xertst (0.0_R8.le.sput_phys,
     1  'faulty input b2stbr_sput_phys')
       if(sput_frc.ne.0.0_R8 .or. sput_phys.ne.0.0_R8) then
         call xertst(sput_dst.lt.ns,
     1    'faulty input parameter b2stbr_sput_dst')
       endif
       call ipgeti ('b2stbr_plate_model', plate_model)
       call ipgetr ('b2stbr_plate_temp', plate_temp)
       call xertst (0.0_R8.lt.plate_temp, 'faulty b2stbr_plate_temp')
       call ipgetr ('b2stbr_plate_thick', plate_thick)
       call xertst (0.0_R8.le.plate_thick, 'faulty b2stbr_plate_thick')
* change atomic physics for coupling to Eirene    ! DPC 1999.06.28
cdpc-2002.03.31       if(use_eirene.eq.0) neutral_rescale=1.0_R8
       if(switch%use_eirene.eq.0 .and.
     &  .not.(switch%neutral_sources_rescale.eq.1.0_R8)) write (*,*)
     &  ' WARNING! Using fluid neutrals with sources rescale /= 1!'
       if(switch%use_eirene.ne.0) then
c when using kinetic neutrals, volume recombination is fully handled on the Eirene side.
c The recombination rate depends on the ion density and not on the neutral density, so the neutral rescaling is not sufficient to turn it off on the B2.5 side
c and hence the reaction rates have to be rescaled here.
         if (switch%neutral_sources_rescale.ne.1.0_R8) then
           write(*,*) 'Adding ',log(switch%neutral_sources_rescale),
     *         ' to relevant rtlra and rtlqr'
           do is=0,ns-2
             scl = switch%neutral_sources_rescale
             do istra=1,nstrai
               if (streql(crcstra(istra),'V').and.
     *           (species_start(istra).eq.is)) then
                 scl = scl*strasclfl(istra)
               endif
             enddo
             if(is_neutral(is) .and. lnext(is,is+1)) then
                rtlra(:,:,is+1)=rtlra(:,:,is+1)+log(scl)
                write(*,*) is,zamax(is),rtlra(0,0,is+1)
                rtlqr(:,:,is+1)=rtlqr(:,:,is+1)+log(scl)      !mb 22.05.18
                write(*,'(i3,1p,2e16.8)')
     *             is,zamax(is),rtlqr(0,0,is+1)
             endif
           enddo
           if (switch%spatial_hybrid.eq.1) then
             write(*,*) 'Neutral sources rescaling has to be turned '//
     *                  'off for hybrid fluid-kinetic model!'
             write(0,*) 'Neutral sources rescaling has to be turned '//
     *                  'off for hybrid fluid-kinetic model!'
           endif
         elseif (switch%spatial_hybrid.eq.0) then
           write(*,*) 'Neutral sources rescaling not turned on '//
     *                'while coupled to Eirene!'
           write(0,*) 'Neutral sources rescaling not turned on '//
     *                'while coupled to Eirene!'
         endif
       endif
       call ipgetr ('b2mndr_dpc_mod_rates_te_hot',te_hot)
       call ipgetr ('b2mndr_dpc_mod_rates_ne_hot_frac',ne_hot_frac)
       if(te_hot.gt.0.0_R8.and.ne_hot_frac.gt.0.0_R8)
     &  call dpc_mod_rates(te_hot, ne_hot_frac)
       call ipgeti ('b2mndr_external_species', external_species)
       if (external_species.ne.0) then
        call xerrab ('WG -- External species not yet supported')
cwdk for now, dummy external state is created in b2mod_main
#ifdef WG_TODO
        call read_b2mod_external(nx,ny,ns,ismain)
       else
        call alloc_b2mod_external_dummy(nx,ny)
#endif
       endif
*   ..find nscx, iscx
*     (indices for neutral hydrogen species)
       state%rt%nscx = 0
       do is = 0, ns-1
        if (is_neutral(is).and.nint(zn(is)).eq.1) then
         call xertst (state%rt%nscx.lt.state%rt%nscxmax,
     1    'too many neutral hydrogenic species')
         state%rt%iscx(state%rt%nscx) = is
         state%rt%nscx = state%rt%nscx+1
        endif
       enddo
       do k = state%rt%nscx, state%rt%nscxmax-1
        state%rt%iscx(k) = -1
       enddo
       k = 0
       do is = 0, ns-1
         if (nint(zn(is)).eq.1.and.
     .      (nint(zamin(is)).eq.1.and.nint(zamax(is)).eq.1)) k = k + 1
       end do
       call xertst (k.eq.1.or.(switch%b2tfnb_xcur.eq.0.0_R8 .and.
     &                         switch%b2tfnb_ycur.eq.0.0_R8),
     .  'Currents assignment to ions only supported for single '//
     .  'hydrogenic species!')
*
       ixfb=rratio(nCv,8) !WG_TODO dummy values provided now; probably need icvfb instead...
       iyfb=rratio(nCv,2) !WG_TODO dummy values provided now
       isfb=ismain
       call ipgetr ('b2mndr_ixfb', ixfb)
       call ipgetr ('b2mndr_iyfb', iyfb)
       call ipgeti ('b2mndr_isfb', isfb)
       call ipgetr ('b2mndr_ne_wanted', ne_wanted)
       call ipgetr ('b2mndr_ne_wanted_time', ne_wanted_time)
       call ipgetr ('b2mndr_ne_wanted_mod_time', ne_wanted_mod_time)
       call xertst (0.le.isfb.and.isfb.lt.ns,
     &   'faulty internal parameter isfb')
csw 09oct2012
       call ipgeti ('b2mndr_old_style', old_style)             !pb 29.07.16
       call ipgeti ('b2mndr_av_read', ird_aver)                !pb 29.09.16
       call ipgeti ('b2mndr_ids_save', ids_save )
       call ipgeti ('b2mndr_ids_av', ids_av )
       call xertst ( 0.le.ids_av, 'Invalid ids_av value')
       summary_av =
     &  ids_av.gt.0 .and. ids_save.eq.0 .and. delta_ids_time.eq.0.0_R8
#ifdef IMAS
       if (ids_save.eq.1) then
         write(*,*) 'Writing IDS every iteration'
       else if (ids_save.gt.1) then
         write(*,*) 'Writing IDS every ',ids_save,' iterations'
       else if (ids_save.lt.0) then
         write(*,*) 'Writing IDS at the end of the run'
       end if
       if (ids_av.eq.1 .and. ntim_batch.gt.0) then
         write(*,*) 'Writing IDS for every averaged solution'
       else if (ids_av.gt.1 .and. ntim_batch.gt.0) then
         write(*,*) 'Writing IDS every ',ids_av,' averaged solutions'
       else if (ids_av.gt.0 .and. ntim_batch.eq.0) then
         ids_av = 0
         write(*,*) 'No batch averaging frequency provided'
         write(*,*) 'Cannot write averaged solution IDS'
       end if
       if (ids_save.ne.0 .or. ids_av.gt.0 .or. delta_ids_time.gt.0.0_R8)
     .  then
         path = ' '
         username = usrnam()
         call ipgetc('b2mndr_user', username )
         call xertst( .not.streql(username,' '),
     .    'User name not defined !')
         not_default = .not.streql(username, usrnam())
         home_dir = '/home/'//trim(username)
         database = 'solps-iter'
         write(version,'(i1)') IMAS_MAJOR_VERSION
#ifdef NO_GETENV
         device_env = ' '
         write(imas_version,'(i1,a1,i2,a1,i1)') IMAS_MAJOR_VERSION,'.',
     .                                          IMAS_MINOR_VERSION,'.',
     .                                          IMAS_MICRO_VERSION
#else
#ifdef USE_PXFGETENV
         CALL PXFGETENV ('HOME', 0, home_dir, lenval, ierror)
         CALL PXFGETENV ('DEVICE', 0, device_env, lenval, ierror)
         CALL PXFGETENV ('IMAS_VERSION', 0,
     .                              imas_version, lenval, ierror)
#else
         call get_environment_variable('HOME', status=ierror,
     .                                         length=lenval)
         if (ierror.eq.0) call get_environment_variable('HOME',
     .                                         value=home_dir)
         call get_environment_variable('DEVICE', status=ierror,
     .                                         length=lenval)
         if (ierror.eq.0) call get_environment_variable('DEVICE',
     .                                         value=device_env)
         call get_environment_variable('IMAS_VERSION', status=ierror,
     .                                         length=lenval)
         if (ierror.eq.0) call get_environment_variable('IMAS_VERSION',
     .                                         value=imas_version)
#endif
         if (.not.streql(device_env,' ')) database = device_env
         if (streql(database,'iter')) database = 'ITER'
#endif
         call ipgetc('b2mndr_device', database )
         call ipgetc('b2mndr_database', database )
         call xertst( .not.streql(database,' '),
     .    'Database not defined !')
         not_default = not_default.or.
     &                (.not.streql(database, device_env).and.
     &                 .not.streql(database,'solps-iter'))
         call ipgeti('b2mndr_pulse_number', shot )
         if (shot.eq.0) call ipgeti('b2mndr_shot_number', shot )
         call ipgeti('b2mndr_run_number', run )
#if AL_MAJOR_VERSION > 4
         call xertst( 0.le.shot, 'Invalid pulse number')
         call xertst( 0.le.run, 'Invalid run number')
#else
         call xertst( 0.lt.shot.and.shot.le.214748,
     .    'Invalid shot number')
         call xertst( 0.le.run.and.run.le.99999, 'Invalid run number')
#endif
         not_default = not_default.or.(shot.gt.0.and.run.ge.0)
         imasdir = trim(home_dir)//'/public/imasdb/'
     &           //trim(database)//'/'//trim(version)
#if AL_MAJOR_VERSION < 5
     &           //'/0'
#endif
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
         CALL PXFGETENV ('IMASDIR', 0, imasdir, lenval, ierror)
#else
         call get_environment_variable('IMASDIR',
     .    status=ierror, length=lenval)
         if (ierror.eq.0)
     >    call get_environment_variable('IMASDIR',value=imasdir)
#endif
#endif
         if (index(imasdir,'imasdb/'//trim(database)).eq.0) then
           l=index(imasdir,'imasdb/')
           m=index(imasdir(l+7:256),'/')
           write(imasdir,'(a)')
     >      imasdir(1:l+6)//trim(database)//trim(imasdir(m+l+6:256))
         end if
#if AL_MAJOR_VERSION < 5
         if (run.ge.10000) then
           m=len_trim(imasdir)
           write(imasdir,'(a)') imasdir(1:m-3)//trim(version)
     &       //'/'//int2str(run/10000)
         end if
#endif
#if AL_MAJOR_VERSION > 4
         call ipgetc('b2mndr_ids_path', path )
         not_default = not_default.or..not.streql(path,' ')
#else
         treename = 'ids'
#endif
         call strip_spaces(path)
         filename='shotnumber.history'
         call find_file(filename,ex)
         if (ex) then
           open(99,file=filename)
   11      continue
           read(99,'(a)',end=12) zeile
           goto 11
         else
           zeile = ' '
         end if
   12    continue
         if (ex) close(99)
         call strip_spaces(zeile)
         if (.not.streql(zeile,' ').or..not.streql(path,' ')) then
           if (.not.not_default) path = trim(zeile)
           if (streql(path,' ')) path = trim(imasdir)
           l=index(path,'$HOME')
           if (l.gt.0) then
             ids_path = trim(home_dir)//trim(path(l+5:256))
             path = ids_path
           end if
           l=index(path,'$SOLPSTOP')
           if (l.gt.0) then
             ids_path = trim(solpstop())//trim(path(l+9:256))
             path = ids_path
           end if
           l=index(path,'$IMASDIR')
           if (l.gt.0) then
             ids_path = trim(imasdir)//trim(path(l+8:256))
             path = ids_path
           end if
           l=index(path,trim(imasdir))
           if (l.eq.0) then
             m=index(path,'/')
           else
             m=index(path(l+len_trim(imasdir):256),'/')
           end if
           absolute_path = l.eq.1.or.(l.eq.0.and.m.eq.1)
           if (absolute_path) then
             ids_path = path
           else if (l.eq.0) then
             ids_path = trim(imasdir)//'/'//trim(path)
             path = ids_path
           end if
         end if
         if (streql(path,' ').or.streql(path,imasdir)) then
#if AL_MAJOR_VERSION > 4
           call xertst( 0.lt.shot, 'Invalid pulse number !')
           ids_path = trim(imasdir)//'/'//int2str(shot)//
     &                               '/'//int2str(run)
#else
           ids_path = imasdir
#endif
         end if
         if (.not.streql(zeile,' ').and..not.not_default) then
           write(*,'(a)')
     &      'IMAS data entry path is '//trim(ids_path)//
     &      ' (from shotnumber.history file)'
         else if (not_default) then
           write(*,'(a)')
     &      'IMAS data entry path is '//trim(ids_path)//
     &      ' (from b2mn.dat file settings)'
         else
           write(*,'(a)')
     &      'IMAS data entry path is '//trim(ids_path)//
     &      ' (from default path settings)'
         end if
         if (.not.streql(ids_path,' ')) then
           filename='shotnumber.history'
           call find_file(filename,ex)
           open(99,file=filename,position='append')
           write(99,'(a)') trim(ids_path)
           close(99)
         end if
       end if
#else
       if (ids_save.ne.0 .or. ids_av.gt.0 .or. delta_ids_time.gt.0.0_R8)
     >  write(*,*)
     1   'Not compiled with IMAS option, IDS will not be written.'
#endif
* ..return
      else
       call xerrab ('error, b2mndr has been called more than once')
      endif

* ..initialisation
#ifdef USE_IMPGYRO
*   ..call IMPGYRO
      call impgyro_init (nx, ny, ns, ismain) !yamoto
      write(*,*) 'b2ig: set ns_ext', state_ext%ns
#endif
      if(switch%use_eirene.ne.0) then
#ifdef B25_EIRENE
        IUNOUT=6
        CALL EIRENE_DEFAULTS_USR
        write(*,*) 'Calling FIND_PARAM'
        CALL EIRENE_FIND_PARAM(0)
        CALL EIRENE_SET_PARMMOD(1)
        write(*,*) 'Called FIND_PARAM'
        write(*,*) '   NDX = ',NDX
        write(*,*) '   NDY = ',NDY
        write(*,*) ' NSTRA = ',NSTRA
        write(*,*) '   NFL = ',NFL
        write(*,*) '  NATM = ',NATM
        write(*,*) '  NMOL = ',NMOL
        write(*,*) '  NION = ',NION
        write(*,*) 'NLIMPS = ',NLIMPS
      else
        ndx=DEF_NXD+DEF_NCUT*5
        ndy=DEF_NYD
        nfl=DEF_NFL
        natm=DEF_NATM
        nmol=DEF_NMOL
        nion=DEF_NION
        npls=DEF_NPLS
        nlim=DEF_NLIM
        nsts=DEF_NSTS
        nstra=DEF_NSTRA
        nsrfs=DEF_NSRFS
        N1ST=DEF_NYD+1
        N2ND=DEF_NXD+MAX(1,DEF_NCUT/2*3)+
     .       DEF_ISOEXTRA*MAX(0,DEF_NCUT/2-1)
        N3RD=DEF_NPHID
        NREAC=7
        IUNOUT=6
        CALL EIRENE_DEFAULTS_USR
        CALL EIRENE_SET_PARMMOD(1)
        CALL EIRENE_ALLOC_CLOGAU
        CALL EIRENE_INIT_CINIT
#else
        call xerrab
     &   ('Compile code with B25_EIRENE option for coupled runs')
#endif
      endif
      NFLAI = NFL
#ifdef B25_EIRENE
      CALL EIRENE_ALLOC_BRAEIR(NCV,NFC,NFL)
      CALL EIRENE_ALLOC_EIRBRA(NCV,NFL,NATM,NMOL,NION,NSTRA)
      CALL EIRENE_ALLOC_CCOUPL(1)
      CALL EIRENE_ALLOC_CTRCEI(1)
      CALL EIRENE_ALLOC_CTRCEI(2)
#else
      call alloc_b2mod_braeir(nCv,nFc,ns,nstra)
      call createB2SourceEir( nCv, ns, nstra, state % sr_eir )
#endif

*   ..read geometry and metric
      call read_b2fgmtry (ninp(1), mpg, geo)
      call init_geometry (mpg, geo, switch)

#ifdef WG_TODO
! can only be calculated after the indirect arrays have been read
      ix1=int(ixfb)
      call xertst(0.le.ix1.and.ix1.lt.nx, 'faulty ixfb (1)')
      iy1=int(iyfb)
      call xertst(0.le.iy1.and.iy1.lt.ny, 'faulty iyfb (1)')
      ix2=rightix(ix1,iy1)
      call xertst(0.le.ix2.and.ix2.lt.nx, 'faulty ixfb (2)')
      iy2=rightiy(ix1,iy1)
      call xertst(0.le.iy2.and.iy2.lt.ny, 'faulty iyfb (2)')
      ix3=topix(ix1,iy1)
      call xertst(0.le.ix3.and.ix3.lt.nx, 'faulty ixfb (3)')
      iy3=topiy(ix1,iy1)
      call xertst(0.le.iy3.and.iy3.lt.ny, 'faulty iyfb (3)')
      ix4=rightix(ix3,iy3)
      call xertst(0.le.ix4.and.ix4.lt.nx, 'faulty ixfb (4)')
      iy4=rightiy(ix3,iy3)
      call xertst(0.le.iy4.and.iy4.lt.ny, 'faulty iyfb (4)')
      ix5=topix(ix2,iy2)
      call xertst(ix4.eq.ix5, 'fault ixfb (5)')
      iy5=topiy(ix2,iy2)
      call xertst(iy4.eq.iy5, 'fault iyfb (5)')
      ax=ixfb-int(ixfb)
      ay=iyfb-int(iyfb)
      write(*,*) 'ix# ', ix1, ix2, ix3, ix4, ix5
      write(*,*) 'iy# ', iy1, iy2, iy3, iy4, iy5
      write(*,*) 'ax,ay ', ax, ay
#endif
!
      call alloc_b2mod_tallies(mpg%nnreg,ns)
#ifdef WG_TODO
cwdk eventually to be moved into init_geometry (or: defined through b2fgmtry)
*   ..find the poloidal number of midplane cell as a number of cell with min Bz
      if(nncut.eq.2) then                                            !srv 22.05.09 {
        nmdpl = rightcut(nncut)
        do ix = rightcut(nncut)+1, rightcut(nncut-1)-1
          if(abs(bb(ix,ny,2)) .lt. abs(bb(nmdpl,ny,2))) nmdpl = ix   !srv 11.01.13
        enddo
      elseif(nncut.eq.1) then
        nmdpl = leftcut(nncut)
        do ix = leftcut(nncut)+1, rightcut(nncut)-1
          if(abs(bb(ix,ny,2)) .lt. abs(bb(nmdpl,ny,2))) nmdpl = ix   !srv 11.01.13
        enddo
      else
        write(*,*) ' b2mndr: nncut = ', nncut
        nmdpl = (nx+2)/2
        do ix = -1, nx
          if(abs(bb(ix,ny,2)) .lt. abs(bb(nmdpl,ny,2))) nmdpl = ix   !srv 11.01.13
        enddo
      endif
      write(*,"(' Number of midplane cell =',i4)") nmdpl             !srv 22.05.09 }
*   ..compute length along x-direction                               !srv 13.02.17 {
      lengthx(:,:) = 0.0_R8
      lengthx(-1,:) = hx(-1,:)*0.5_R8
      do iy = -1, ny
       do ix = 0, nx
        if(on_closed_surface(ix,iy)) then
         if(rightix(ix,iy).ne.nx+1) then
          lengthx(ix,iy) = lengthx(leftix(ix,iy),iy) +
     *      (hx(ix,iy)+hx(rightix(ix,iy),rightiy(ix,iy)))*0.5_R8
         else
          lengthx(ix,iy) = hx(ix,iy)*0.5_R8
         endif
        else
         if(leftix(ix,iy).ne.-2) then
          lengthx(ix,iy) = lengthx(leftix(ix,iy),iy) +
     &     (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*0.5_R8
         else
          lengthx(ix,iy) = hx(ix,iy)*0.5_R8
         endif
        endif
       enddo
      enddo                                                          !srv 13.02.17 }
*   ..compute length along y-direction                               !srv 24.04.02 {
      lengthy(:,-1) = hy(:,-1)*0.5_R8
      do iy = 0, ny
       do ix = -1, nx
        lengthy(ix,iy) = lengthy(ix,iy-1) +
     *      (hy(ix,iy)+hy(ix,iy-1))*0.5_R8
       enddo
      enddo                                                          !srv 24.04.02 }
*
      len_surf=0.0d0                                                 !srv 02.05.09 {
      do ix=-1,nx
        do iy=1,ny
          len_surf(ix,iy)=len_surf(ix,iy-1)+hy(ix,iy-1)
        enddo
      enddo
*
      if (nncut.gt.0) then
        do ix=-1,nx
          crdy_surf(ix,-1:ny) = len_surf(ix,-1:ny) -
     .                          len_surf(ix,jsep+1)
        enddo
      else
        crdy_surf = len_surf
      end if
*
      do ix=-1,nx
        do iy = -1, ny-1
          crdy_c_cell(ix,iy) = (crdy_surf(ix,iy) +
     &                          crdy_surf(ix,iy+1))/2
        enddo
        iy = ny
        crdy_c_cell(ix,iy) = crdy_surf(ix,iy)
      enddo                                                          !srv 02.05.09 }
#endif
*
*   ..read initial state
      call read_b2fstate (ninp(3), mpg%nCv, mpg%nFc, ns, state)
      if (switch%solve_keps.gt.0 .and.
     &   .not.smax(nCv,state%pl%kt,1).gt.0.0_R8)
     .  state%pl%kt = state%pl%ti*0.1_R8
!      if (switch%solve_keps.gt.1 .and.
!     &   .not.smax(nCv,state%pl%zt,1).gt.0.0_R8)
!     .  state%pl%zt = state%pl%ti*1e5_R8
      call createB2Diagnostic( nCv, ns, mpg%nnreg, state % diag )
!WG_TODO      call init_state(mpg, geo, state)

*   initialize 2d batch averages
      if (ibatch_av_all > 0) call batch_av_all_init (nCv, ns)

*   initialize and/or read running averages
      if ((iav_run > 0) .or. (ird_aver > 0)) then
        lrer = (ird_aver == 0)
        call run_av_init (iav_cont, nCv, ns, .true., lrer, state_avg)
      end if

*   get averaged plasma profiles
      if (ird_aver > 0)
     &  call run_av_get_plasma (state%pl%na, state%pl%te,
     &             state%pl%ti, state%pl%ua, state%pl%po,
     &             state%pl%kt, state%pl%zt, nCv, ns, state_avg)

      write (*,*) 'Starting conversion to unstructured format...'

!WG_TODO: move this back to original location...
      if (switch%b2stbc_boundary_namelist.ge.1) then
        call initialize_b2mod_boundary(ns,mpg,0.0_R8)
      else
        call xerrab('b2mndr -- specify BCs using boundary_namelist')
      endif

      if(switch%b2mndr_iout.ne.0) then                                             !srv 25.01.16 {
       call my_out_us(70,nCv,0,state%dv%ne,'b2mndr_ne')
       do is = 0,ns-1
        write (chns,'(i3.3)') is
        call my_out_us(70,nCv,0,state%pl%na(1,is),'b2mndr_na'//chns)
        call my_out_us(70,nCv,0,state%pl%ua(1,is),'b2mndr_ua'//chns)
        call my_out_us(70,nFc,1,state%dv%uadia(1,0,is),
     .    'b2mndr_uadia_th'//chns)
        call my_out_us(70,nCv,0,state%dv%uadia(1,1,is),
     .    'b2mndr_uadia_r'//chns)
        call my_out_us(70,nFc,1,state%dv%fna(1,0,is),
     &    'b2mndr_fna_th'//chns)
        call my_out_us(70,nFc,1,state%dv%fna(1,1,is),
     &    'b2mndr_fna_r'//chns)
       enddo
       call my_out_us(70,nCv,0,state%pl%te,'b2mndr_te')
       call my_out_us(70,nCv,0,state%pl%te/ev,'b2mndr_te_eV')
       call my_out_us(70,nCv,0,state%pl%ti,'b2mndr_ti')
       call my_out_us(70,nCv,0,state%pl%ti/ev,'b2mndr_ti_eV')
       call my_out_us(70,nCv,0,state%pl%po,'b2mndr_po')
       call my_out_us(70,nFc,1,state%dv%fch(1,0),'b2mndr_fch_th')
       call my_out_us(70,nFc,1,state%dv%fch(1,1),'b2mndr_fch_r')
       call my_out_us(70,nFc,1,state%dv%fch_p(1,0),'b2mndr_fchx_p')
       call my_out_us(70,nFc,1,state%dv%fch_p(1,1),'b2mndr_fchy_p')
       call my_out_us(70,nFc,1,state%dv%fhe(1,0),'b2mndr_fhe_th')
       call my_out_us(70,nFc,1,state%dv%fhe(1,1),'b2mndr_fhe_r')
       call my_out_us(70,nFc,1,state%dv%fhi(1,0),'b2mndr_fhi_th')
       call my_out_us(70,nFc,1,state%dv%fhi(1,1),'b2mndr_fhi_r')
       call my_out_us(70,nFc,1,geo%fcQgam(1,0),'b2mndr_fcQgam0_cosg')
       call my_out_us(70,nFc,1,geo%fcQgam(1,1),'b2mndr_fcQgam1_sing')
       call my_out_us(70,nFc,1,geo%fcQalf(1,0),'b2mndr_fcQalf0_cosa')
       call my_out_us(70,nFc,1,geo%fcQalf(1,1),'b2mndr_fcQalf1_sina')
       call my_out_us(70,nFc,1,geo%fcQbet(1,0),'b2mndr_fcQbet0_cosb')
       call my_out_us(70,nFc,1,geo%fcQbet(1,1),'b2mndr_fcQbet1_sinb')
       call my_out_us(70,nFc,1,geo%fcHc(1,1),'b2mndr_fcHc1_0p5h')
       call my_out_us(70,nFc,1,geo%fcHc(1,2),'b2mndr_fcHc2_0p5h')
       call my_out_us(70,nFc,1,geo%fcHt,'b2mndr_fcHt_hti')
       call my_out_us(70,nCv,0,geo%cvVol,'b2mndr_cvVol')
       call my_out_us(70,nCv,0,geo%cvHz,'b2mndr_cvHz')
       call my_out_us(70,nCv,0,geo%cvEb(1,0),'b2mndr_cvEbX')
       call my_out_us(70,nCv,0,geo%cvEb(1,1),'b2mndr_cvEbY')
       call my_out_us(70,nCv,0,geo%cvEb(1,2),'b2mndr_cvEbZ')
       call my_out_us(70,nFc,1,geo%fcEb(1,0),'b2mndr_fcEbX')
       call my_out_us(70,nFc,1,geo%fcEb(1,1),'b2mndr_fcEbY')
       call my_out_us(70,nFc,1,geo%fcEb(1,2),'b2mndr_fcEbZ')
      endif                                                          !srv 25.01.16 }
      do is = 0,ns-1                                                 !srv 25.01.16 {
        do iCv=1,nCv
          state % pl % na(iCv,is) = min(max(state % pl % na(iCv,is),
     &      switch % b2mndr_na_min),switch % b2mndr_na_max)
        enddo
      enddo                                                          !srv 25.01.16 }
*   ..reset initial state if using exact solution
       call set_exact_solution(nCv, ns, switch, state % pl, state % dv)
*   ..test plasma state
      call b2xvps (nCv, nFc, ns, state % pl, state % dv)
*   ..compute log-log electron rate coefficients
      call b2spel (nCv, ns0, ev, state%pl%te, state%dv%ne, state%rt)
*   ..compute electron rate coefficients
      call xertst (0.le.ismain.and.ismain.lt.ns0,
     . 'invalid main plasma species index ismain')
      call b2sqel (nCv, ns0, ismain, switch, ev, state%pl%te, state%rt,
     .  state%rtw)
*   ..extend initial state to species (ns0:ns-1)
      if (ns0.lt.ns) then
       call sfill (nCv*(ns-ns0), switch % b2mndr_na_new, 
     .   state % pl % na(1,ns0), 1)
       call sfill (nCv*(ns-ns0), 0.0_R8, state % pl % ua(1,ns0), 1)
       call sfill (nCv*(ns-ns0), 0.0_R8, state % dv % kinrgy(1,ns0), 1)  !xpb
       call sfill (nFc*2*(ns-ns0), 0.0_R8, state%dv%uadia(1,0,ns0), 1)   !xpb
       call sfill (nFc*2*(ns-ns0), 0.0_R8, state%dv%vadia(1,0,ns0), 1)
       call sfill (nFc*2*(ns-ns0), 0.0_R8, state%dv%wadia(1,0,ns0), 1)
       call sfill (nFc*2*(ns-ns0), 0.0_R8, state%dv%vaecrb(1,0,ns0), 1)
       call sfill (nFc*2*(ns-ns0), 0.0_R8, state%dv%fna(1,0,ns0), 1)
       call sfill (nFc*2*(ns-ns0), 0.0_R8, state%dv%fna_mdf(1,0,ns0), 1)
       call sfill (nFc*2*(ns-ns0), 0.0_R8, state%dv%fna_fcor(1,0,ns0),1)
       call sfill (nFc*2*(ns-ns0), 0.0_R8, state%dv%fna_nodrift(1,0,ns0)
     .   , 1)
       call sfill (nFc*2*(ns-ns0), 0.0_R8, state%dv%fna_he(1,0,ns0), 1)
       call sfill (nFc*2*(ns-ns0), 0.0_R8, state%dv%fnaPSch(1,0,ns0), 1)
       call sfill (nFc*2*(ns-ns0), 0.0_R8, state%dv%fna_eir(1,0,ns0), 1)
       call sfill (nFc*2*(ns-ns0), 0.0_R8, state%dv%fna_32(1,0,ns0), 1)
       call sfill (nFc*2*(ns-ns0), 0.0_R8, state%dv%fna_52(1,0,ns0), 1)
      endif
*   ..rescale densities and associated fluxes
      if (density_rescale.ne.1.0_R8) then
csw 03oct2013 (this just quick and dirty for rescaling Helium
c   --> better way doing this is necessary)
        if(density_rescale < 0.d0) then
          do is=0,ns-1
           if (is_neutral(is).or.nint(zn(is)).ne.2) cycle
           state%pl%na(:,is) = state%pl%na(:,is)*(-density_rescale)
           state%dv%ne(:) = state%dv%ne(:) +
     .      state%pl%na(:,is)*state%rt%rza(:,is)*(-density_rescale-1.0)
           state%dv%fna(:,:,is) =
     .      state%dv%fna(:,:,is)*(-density_rescale)
           state%dv%fna_mdf(:,:,is) =
     .      state%dv%fna_mdf(:,:,is)*(-density_rescale)
           state%dv%fna_fcor(:,:,is) =
     .      state%dv%fna_fcor(:,:,is)*(-density_rescale)
           state%dv%fna_nodrift(:,:,is) =
     .      state%dv%fna_nodrift(:,:,is)*(-density_rescale)
           state%dv%fna_he(:,:,is) =
     .      state%dv%fna_he(:,:,is)*(-density_rescale)
           state%dv%fnaPSch(:,:,is) =
     .      state%dv%fnaPSch(:,:,is)*(-density_rescale)
           state%dv%fna_eir(:,:,is) =
     .      state%dv%fna_eir(:,:,is)*(-density_rescale)
           state%dv%fna_32(:,:,is) =
     .      state%dv%fna_32(:,:,is)*(-density_rescale)
           state%dv%fna_52(:,:,is) =
     .      state%dv%fna_52(:,:,is)*(-density_rescale)
           state%dv%kinrgy(:,is) =
     .      state%dv%kinrgy(:,is)*(-density_rescale)
           write(*,*) 'Density ',is,' rescaled by ',-density_rescale
          enddo
        else
         state%pl%na = state%pl%na * density_rescale
         state%dv%ne = (state%dv%ne-state_ext%ne) * density_rescale +
     .    state_ext%ne
         state%dv%fna = state%dv%fna * density_rescale
         state%dv%fna_mdf = state%dv%fna_mdf * density_rescale
         state%dv%fna_fcor = state%dv%fna_fcor * density_rescale
         state%dv%fna_nodrift = state%dv%fna_nodrift * density_rescale
         state%dv%fna_he = state%dv%fna_he * density_rescale
         state%dv%fnaPSch = state%dv%fnaPSch * density_rescale
         state%dv%fna_eir = state%dv%fna_eir * density_rescale
         state%dv%fna_32 = state%dv%fna_32 * density_rescale
         state%dv%fna_52 = state%dv%fna_52 * density_rescale
         state%dv%fne_eir = state%dv%fne_eir * density_rescale
         state%dv%fne_32 = state%dv%fne_32 * density_rescale
         state%dv%fne_52 = state%dv%fne_52 * density_rescale
         state%dv%fni_32 = state%dv%fni_32 * density_rescale
         state%dv%fni_52 = state%dv%fni_52 * density_rescale
         state%dv%fhe = state%dv%fhe * density_rescale
         state%dv%fhe_mdf = state%dv%fhe_mdf * density_rescale
         state%dv%fhePSch = state%dv%fhePSch * density_rescale
         state%dv%fhe_eir = state%dv%fhe_eir * density_rescale
         state%dv%fhi = state%dv%fhi * density_rescale
         state%dv%fhi_mdf = state%dv%fhi_mdf * density_rescale
         state%dv%fhiPSch = state%dv%fhiPSch * density_rescale
         state%dv%fhi_eir = state%dv%fhi_eir * density_rescale
         state%dv%fch = state%dv%fch * density_rescale
         state%dv%fch_32 = state%dv%fch_32 * density_rescale
         state%dv%fch_52 = state%dv%fch_52 * density_rescale
         state%dv%fch_p = state%dv%fch_p * density_rescale
         state%dv%fchdia = state%dv%fchdia * density_rescale
         state%dv%fchin = state%dv%fchin * density_rescale
         state%dv%fchvispar = state%dv%fchvispar * density_rescale
         state%dv%fchvisper = state%dv%fchvisper * density_rescale
         state%dv%fchvisq = state%dv%fchvisq * density_rescale
         state%dv%fchinert = state%dv%fchinert * density_rescale
         state%dv%floe_noc = state%dv%floe_noc * density_rescale
         state%dv%floi_noc = state%dv%floi_noc * density_rescale
         state%dv%kinrgy = state%dv%kinrgy * density_rescale
         write(*,*) 'All densities and fluxes rescaled by ',
     .    density_rescale
        endif
      endif
*   ..rescale neutrals and associated quantities
      if(ncall.eq.0.and.switch%neutral_rescale.ne.1.0_R8) then
        do is=0,ns-1
          if(is_neutral(is)) then
            state%pl%na(:,is)=state%pl%na(:,is)*switch%neutral_rescale
            state%dv%fna(:,:,is)=state%dv%fna(:,:,is)*
     .       switch%neutral_rescale
            state%dv%fna_mdf(:,:,is) =
     .       state%dv%fna_mdf(:,:,is)*switch%neutral_rescale
            state%dv%fna_fcor(:,:,is) =
     .       state%dv%fna_fcor(:,:,is)*switch%neutral_rescale
            state%dv%fna_nodrift(:,:,is) =
     .       state%dv%fna_nodrift(:,:,is)*switch%neutral_rescale
            state%dv%fna_he(:,:,is) =
     .       state%dv%fna_he(:,:,is)*switch%neutral_rescale
            state%dv%fnaPSch(:,:,is) =
     .       state%dv%fnaPSch(:,:,is)*switch%neutral_rescale
            state%dv%fna_eir(:,:,is) =
     .       state%dv%fna_eir(:,:,is)*switch%neutral_rescale
            state%dv%fna_32(:,:,is) =
     .       state%dv%fna_32(:,:,is)*switch%neutral_rescale
            state%dv%fna_52(:,:,is) =
     .       state%dv%fna_52(:,:,is)*switch%neutral_rescale
            state%dv%kinrgy(:,is)=state%dv%kinrgy(:,is)*
     .       switch%neutral_rescale
          endif
        enddo
        write(*,*) 'Rescaled neutrals by ',switch%neutral_rescale
      endif
*   ..if we rescaled or added species, recompute dependent quantities
*   ..until convergence
      if (ns0.lt.ns .or. density_rescale.ne.1.0_R8) then
        if (ns0.lt.ns) then
          do is = ns0, ns-1
            state%rt%rza(:,is) = (zamin(is)+zamax(is))/2.0_R8
            state%rt%rz2(:,is) = ((zamin(is)+zamax(is))/2.0_R8)**2
            state%rt%rpt(:,is) = (potmin(is)+potmax(is))/2.0_R8
            state%rt%rpi(:,is) = (potimin(is)+potimax(is))/2.0_R8
          enddo
          call b2xpne (nCv, ns, state % rt % rza, state % pl % na,
     .       state_ext % ne, state % dv % ne)
        endif
        max_delta = 1.0_R8
        is = 0
        do while (max_delta.gt.1.0e-6_R8 .and. is.lt.10)
          is = is + 1
          call b2scopy (nCv*ns, state%rt%rza, 1, rza0, 1)
          call b2scopy (nCv*ns, state%rt%rz2, 1, rz20, 1)
          call b2scopy (nCv*ns, state%rt%rpt, 1, rpt0, 1)
          call b2scopy (nCv*ns, state%rt%rpi, 1, rpi0, 1)
          call b2scopy (nCv, state%dv%ne, 1, state%psnc%ne, 1)
*    ..compute log-log electron rate coefficients
          call b2spel (nCv, ns, ev, state%pl%te, state%dv%ne, state%rt)
*   ..compute electron rate coefficients
          call b2sqel (nCv, ns, ismain, switch, ev, state%pl%te,
     .      state%rt, state%rtw)
          call b2xpne (nCv, ns, state%rt%rza, state%pl%na,
     .       state_ext%ne, state%dv%ne)
          max_delta = max(
     1     maxval(dabs((state%rt%rza-rza0)/(state%rt%rza+1.0_R8))),
     2     maxval(dabs((state%rt%rz2-rz20)/(state%rt%rz2+1.0_R8))),
     3     maxval(dabs((state%rt%rpt-rpt0)/(state%rt%rpt+1.0_R8))),
     4     maxval(dabs((state%rt%rpi-rpi0)/(state%rt%rpi+1.0_R8))),
     5     maxval(dabs((state%dv%ne - state%psnc%ne)/
     6       (state%dv%ne + na_eps))))
          write (*,*) 'Convergence of bundled rates = ',is,max_delta
        enddo
      endif
*   ..rewrite input file with inverted velocity
      if(inverse_ua.ne.0) then                                           !srv 04.02.15 {
        call cfopen (70,'b2fstati_with-ua','new','un*formatted')
        idum(0) = nCv
        idum(1) = nFc
        idum(2) = ns
        call cfverw (70, newversion)
        call cfwuin (70, 3, idum, 'nCv,nFc,ns')
        label = trim(label)//' Inverted velocities'
        call cfwuch (70, 120, label, 'label')
        call b2wuzd (70, newversion, ns, zamin, zamax, zn, am)
        state % pl % ua = -state % pl % ua
*     ..write plasma state
        call write_b2fstate (70, nCv, nFc, ns, state)
        close(70)
      endif                                                              !srv 04.02.15 }
*   ..compute initial ni
      call b2xpni (nCv, ns, state%pl%na, state%dv%ni)
      call b2xpnn (nCv, ns, state%pl%na, state%dv%nn)
*   ..initialise previous state and fluxes
      call getB2PlasmaSnapshot(state%pl, state%dv, state%psnc)
      call getB2PlasmaSnapshot(state%pl, state%dv, state%psnl)
!WG_TODO*   ..write output control
!WG_TODO      write (nout(0),'(2(/2x,a,2x,a))')
!WG_TODO     &  'idout0:', idout0, 'idout1:', idout1
!WG_TODO*   ..write geometric information
!WG_TODO      call b2wfgm (nout(0), idout0, nx, ny)
!WG_TODO*   ..write summary of initial state
!WG_TODO      call b2wfpi (nout(0), idout0, nx, ny, ns, ne, na, ua,
!WG_TODO     &  te, ti, po, fna, fhe, fhi, fch)

* ..main loop
      if (switch%b2mndr_stim.ge.0.0_R8) then
        continued = tim.eq.switch%b2mndr_stim
        tim = switch%b2mndr_stim
        run_start_time = switch%b2mndr_stim
      else
        continued = .true.
#ifdef IMAS
        run_start_time = IDS_REAL_INVALID
#else
        run_start_time = 0.0_R8
#endif
      end if
      call xertst (etim.eq.0.0_R8 .or. etim.gt.tim,
     & 'faulty input etim')
      if (etim.gt.tim) then
        run_end_time = etim
      else
#ifdef IMAS
        run_end_time = IDS_REAL_INVALID
#else
        run_end_time = 0.0_R8
#endif
      end if
      save_plasma_time = tim + delta_plasma_time
      save_ids_time = tim + delta_ids_time
      save_cdfmovie_time = tim + delta_cdfmovie_time
*   ..initialise tracing routines
      write(*,*) 'Calling alloc_b2mod_elements'
      call alloc_b2mod_elements(ns,zn,am)
      call alloc_b2mod_neutrals
     & (ns,natm,nmol,nion,nstraid,nsts,sput_dst)
      natmi = 0
      do is = 0, ns-1
        if (is_neutral(is)) natmi = natmi+1
      enddo
      call read_neutrals_namelist(ns, mpg, switch, .true.)
      call alloc_b2mod_user(geo,mpg,ns,nlim,nmol,switch)
      if(switch%b2tqna_transport_namelist.eq.1) then
        call alloc_transport_namelist(ns)
        call read_b2mod_transport_namelist
      endif
!csc  reading transpot input is needed here otherwise when optimization is running it will override the optimized parameters
      if (switch%b2tqna_transport_inputfile .eq. 1) then
        call read_transport_input(transport_ip_filename)
      endif
      call alloc_b2mod_numerics(nCv,ns,mpg%nnreg, switch)
      allocate(
     *   old_erosion(nwall,ntrack), old_deposition(nwall,ntrack))
      if (flag_optim .or. switch%b2optim_namelist .eq. 1) then
        call read_b2mod_par_opt(ns, mpg, switch)
        allocate(par_opt_phys(npar_opt))
        allocate(xsave(npar_opt))
        par_opt_phys(1:npar_opt) = x0(1:npar_opt)
        xsave(1:npar_opt) = x0(1:npar_opt)
      endif
      call xertst (natmi.le.DEF_NATM,
     .  'Increase DEF_NATM in b2mod_dimensions and recompile !')
      call alloc_b2mod_diag(nCv,nFc,ns,mpg%nnreg(0),nstra,natm,switch)
cxpb initialize the data necessary for boundary conditions
      call b2stbr_init (ns, switch, mpg)
cxpb read the Eirene output for 0-timestep runs and
cxpb for neutral fluxes used in b2stbc_phys
      if(switch%b2stbc_feedback.eq.1.or.lfeedback) then
        call init_feedback(nFc, ns, switch, mpg, geo)
      endif
!WG_TODO      call calc_mapping_rc (nx, ny, ns, mpg)
      nlimi = 0
#ifdef B25_EIRENE
      filename=fort_lc//'44'
      call find_file(filename,ex)
      CALL EIRENE_SETCON
      if (switch%use_eirene.ne.0) CALL EIRENE_WNEUTRALS_INIT
      if (ex.or.ids_save.ne.0.or.ids_av.gt.0.or.
     &    delta_ids_time.gt.0.0_R8) then
        call alloc_b2mod_b2plot(nCv,nFc,nVx,ns,state_ext%ns,wklng,
     .                          natm,nmol,nion,nlimps)
        if (switch%use_eirene.ne.0.and.ex) then
          open(unit=44,file=filename)
          call ngread(mpg,geo,ns,switch,ismain,iflag)
        endif
      endif
#else
      if (ids_save.ne.0.or.ids_av.gt.0.or.delta_ids_time.gt.0.0_R8)
     > call alloc_b2mod_b2plot(nCv,nFc,nVx,ns,state_ext%ns,wklng,
     .                         natm,nmol,nion,nlimps)
#endif
c new CDF movie option (1st frame at t=0)
      if(tim.ge.(save_cdfmovie_time-dtim/2.0_R8).and.
     1 delta_cdfmovie_time.gt.0.0_R8) then
        if(switch%use_eirene.eq.0) then
#ifndef NO_CDF
          call cdfmovie(ncid,nCv,mpg%nCi,ns,geo,state,switch)
#else
          write(*,*) 'Compile b2mndr without NO_CDF option for movies!'
#endif
        endif
        save_cdfmovie_time=save_cdfmovie_time+delta_cdfmovie_time
      endif
      idx = 0
      write (*,*) 'b2mndr_0, end: nCv = ',mpg%nCv

      !Initialize residual if using Crank-Nicolson scheme
      if (switch%b2mndt_ckn.eq.2.and.itim.eq.0) then
        !Dummy call to b2mndt
        !Temporarily set dummy switch to 1
        switch%b2mndt_dummy = 1
        !Temporarily set rxf to 0 to not update any physical quantities
        switch%b2mndt_rxf = 0.0_R8
        !Temporarily set inner iterations to 1
        switch%nstg(0) = 1
        switch%nstg(1) = 1
        switch%nstg(2) = 1
        call b2mndt (nout, nCv, nFc, nVx, ns,
     &       ismain, ismain0, state%rt%nscx, state%rt%nscxmax,
     &       state%rt%iscx, itim, dtim, ntim,
     &       switch, geo, mpg, state, state_ext, state_avg, ierr)
        !Read again rxf switch and inner iterations switches
        call ipgetr ('b2mndt_rxf', switch%b2mndt_rxf)
        call ipgeti ('b2mndt_nstg0', switch%nstg(0))
        call ipgeti ('b2mndt_nstg1', switch%nstg(1))
        call ipgeti ('b2mndt_nstg2', switch%nstg(2))
        !Set back dummy switch to default value of 0
        switch%b2mndt_dummy = 0
        !Update plasma snapshots with initialized residuals (without time-stepping sources)
        state%psnl%resco0 = state%dv%resco0
        state%psnl%resmo0 = state%dv%resmo0
        state%psnl%reshe0 = state%dv%reshe0
        state%psnl%reshi0 = state%dv%reshi0
        state%psnl%reshn0 = state%dv%reshn0
        state%psnl%reskt0 = state%dv%reskt0
        state%psnl%reszt0 = state%dv%reszt0
        endif

      call subend()
      return
#ifdef USE_MPI
  999 call xerrab
     . ('Error while reading ATOMIC_PHYSICS_RESCALE namelist!')
#endif
      end subroutine b2mndr_0

      subroutine b2mndr_1 (nout, ns, switch, geo, mpg, state,
     &  state_ext, state_avg, J)
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_batch_average
     & , only : e_she, e_shi, e_ua, ua_mean, e_kt, po_mean, te_mean,
     &          na_mean, ti_mean, sna_mean, kt_mean, she_mean, e_po,
     &          e_smo, shi_mean, naver, e_te, smo_mean, e_ti, e_sna,
     &          e_na
      use b2mod_recycle
     & , only : int0r, int0l, int3r, int3l, int6r, int6l, int2r, int2l,
     &          int5r, int5l, int1r, int1l, int4r, int4l, uold, mold,
     &          told, fna_mol
      use b2mod_math
     & , only : cutlo, cutll, b2mod_math_initialised, small_r4_constant
      use b2mod_facdrift_ExB
     & , only : ncall_drift, facdrift_scalar, fac_exb_scalar,
     &          fac_vis_scalar
#endif
      use b2mod_ad
#ifdef IMAS
      use b2mod_ual
     , , only : new_batch_edge, put_batch_edge
#endif
      implicit none
#ifdef JET
      integer nout(0:11), ns, idum(0:9)
#else
      integer nout(0:10), ns, idum(0:9)
#endif
#ifdef WG_TODO
      integer nx, ny
#endif
#ifdef WINDOWS
      logical(4) :: lresult
#endif
      type(mapping), intent (inout)    :: mpg
      type(geometry), intent (inout)   :: geo
      type(B2State), intent (inout)    :: state
      type(B2StateExt), intent (inout) :: state_ext
      type(B2Average), intent (inout)  :: state_avg
      type(switches), intent (inout)   :: switch
      real(kind=R8) :: J(nncf)
      integer nCv, nFc, nVx
      integer i, is, icf, iCv
#ifdef IMAS
      integer l, tmp_run, time_slice_index, num_batch_slices
#ifndef NAGFOR
      integer status
#endif
#if AL_MAJOR_VERSION > 4
      character(len=:), allocatable :: message
      character(len=STRMAXLEN) :: uri, uri_source, uri_dest
      character(len=256) :: olddir
#endif
#endif
      character*256 filename, fort31name
#if ( defined(B25_EIRENE) || defined(IMAS) )
      character*512 systemarg
#endif
#ifdef B25_EIRENE
#ifdef NAGFOR
      integer ierror, status
#endif
#endif
      character*2 ss
      logical ok, found
*   ..procedures
      real(kind=R8) :: epoch_seconds
      external epoch_seconds
#ifdef IMAS
      logical streql
      external streql
      intrinsic index
#endif
#if ( defined(CRAY) && defined(_OPENMP) )
      real(kind=R8) :: omp_get_wtime
      external omp_get_wtime
#endif

      call subini ('b2mndr_1')
*   ..no error so far
      ierr = 0
*   ..initialise counters
      itim = 0
      itim_plas = 0
      cpuval = 0.0
      elapsedval = 0.0_R8
      stack_ptr = 0
      first_time_step=.true.
      quit=.false.
      quitexist=.false.
      quitexist_=.false.                                            !srv 13.04.11
      quit_residual=.false.                                         !jwk
      res_max = 10.0_R8*res_quit
      continued = continued .and. ntim.gt.0
      ne_wanted_next_time=tim
#ifdef _OPENMP
      cpuinit = omp_get_wtime()
#else
      call cpu_time(cpuinit)
#endif
      elapsedinit=epoch_seconds()
cwdk  set some integers
      nCv = mpg%nCv
      nFc = mpg%nFc
      nVx = mpg%nVx

      write (*,*) 'b2mndr_1, start: nCv = ',nCv
#ifdef FIXED_POINT
! The FIXED POINT compilation option here and below transforms the time-stepping loop into
! a fixed-point type of iterator such that Tapenade can apply the two-phase checkpointing strategy
! and save only the final converged primal state instead of all the intermediate states. The AD pragma
! below tells Tapenade that a fixed-point loop is next and which variable is the 'state' variable.
! Note that when compiling the 'FIXED_POINT' version some parts of the standard time-stepping are 
! eliminated, otherwise Tapenade does not recognize it as a real fixed-point loop.
! $AD FP-LOOP "state"
      do while ((res_max.ge.res_quit) .and. (itim.lt.ntim) .and.
     &          (.not.quit))
#else
    1 continue
#endif
      write(*,'(1x,a,i9,1p,g14.7,i9,i3,1x,l1)')
     . 'b2mndr_00:itim,dtim,ntim,stack_ptr',
     1            itim,dtim,ntim,stack_ptr,quit

#ifdef _OPENMP
      cpustart = omp_get_wtime()
#else
      call cpu_time(cpustart)
#endif
      elapsedstart=epoch_seconds()
#ifndef FIXED_POINT
*   ..loop test
*!!   (Note: better loop control yet to be developed)
      if (first_time_step.or.
     &    (.not.quit.and.itim.lt.ntim.and.ierr.eq.0.and.
     &     (etim.gt.tim.or.etim.eq.0.0_R8))) then !{
#endif

#ifdef WG_TODO
!!!
!!! DPC experimental
!!!
        if( ne_wanted .gt. 0.0_R8 .and. tim .lt. ne_wanted_time .and.
     &        tim .ge. ne_wanted_next_time) then
           ne_wanted_next_time = ne_wanted_next_time +
     &          ne_wanted_mod_time
           factor=ne_wanted /
     &          (  ax*ay*ne(ix1,iy1) +
     &             (1.0_R8-ax)*ay*ne(ix2,iy2) +
     &             ax*(1.0_R8-ay)*ne(ix3,iy3) +
     &             (1.0_R8-ax)*(1.0_R8-ay)*ne(ix4,iy4) )
           write(*,'(a,1p,3g15.7)')
     &          'time, factor, created_electrons = ', tim, factor,
     &          sum(ne*vol)*(factor-1.0)/dtim
           na(:,:,isfb) = na(:,:,isfb) +
     &          ne(:,:) * (factor-1) / rza(:,:,isfb)
           ne(:,:) = ne(:,:) * factor
        endif
!!!
!!! end of experiment
!!!
#endif

*    ..set snapshot values. The behaviour is made to be backward compatible with ua0, uap etc in 3.0.9
c      nap, uap, etc from 3.0.9 -> psnl = Plasma SNapshot Last
c      na0, ua0, etc from 3.0.9 -> psnc = Plasma SNapshot Current
c      the code below copies ua -> psnl%ua (= uap) and ua -> psnc%ua (= ua0) etc, which is what currently happens in the b2mxnp and b2mxnu routines in 3.0.9
        call getB2PlasmaSnapshot(state%pl, state%dv, state%psnl)
        call getB2PlasmaSnapshot(state%pl, state%dv, state%psnc)

        call b2trcs
csw 09oct2012 set default potential already at this point (BCs)
        if (switch%pot_eq.eq.0) then                           !srv 13.07.05 {
          state%pl%po=0.0_R8
        else if (switch%pot_eq.eq.2) then
          state%pl%po=3.1_R8*state%pl%te/qe
        endif

*    ..perform one time step
        ok=.false.
#ifndef FIXED_POINT
        do while(.not.ok)
#endif
          write(*,'(1x,a,i9,1p,g14.7,i9,i3)')
     .     'b2mndr_ok:itim,dtim,ntim,stack_ptr',
     1                itim,dtim,ntim,stack_ptr
          call b2mndt (nout, nCv, nFc, nVx, ns,
     .      ismain, ismain0, state%rt%nscx, state%rt%nscxmax,
     .      state%rt%iscx, itim, dtim, ntim,
     .      switch, geo, mpg, state, state_ext, state_avg, ierr)
#ifndef FIXED_POINT
          if(ierr.eq.10) then
            stack_ptr=stack_ptr+1
            call xertst(stack_ptr.le.20, 'Too deep in b2mndr')
            stack_ntim(stack_ptr)=ntim-itim     ! dpc 2000.06.24 was ntim+10
            stack_dtim(stack_ptr)=dtim
            dtim=dtim/10.0_R8
            ntim=itim+10
            itim_plas=itim_plas-10
            call putB2PlasmaSnapshot( state%pl, state%dv, state%psnc )
          else
            ok=ierr.eq.0
          endif
        enddo
*    ..call cost function
        call b2usr_cost_function(nCv, nFc, nVx, ns, geo, mpg, state,
     &    state_ext, switch, J)
        do icf = 1, ncf
          write (ss,'(I1)') icf
          if (icf.gt.9) write (ss,'(I2)') icf
          write(*,*) 'Cost function value '//ss//': ',J(icf)
        end do
*    ..produce data for movie
        if (mod(itim,mvinc).eq.0.and.itim/mvinc.le.mvnum) then
          call b2mwmv (nout(3), nCv, state, itim, tim, dtim)
        endif
*    ..extra writes of b2fplasmf.xxxx
        if (mod(itim_plas,plasinc).eq.0.and.nwrit_plas.le.plasnum
     &    .and.plasnum.gt.0.and.plasinc.gt.0) then
          nwrit_plas = nwrit_plas + 1
          write(savefile_p(11:14),'(i4)') edition_p
          edition_p=edition_p+1
          do i=11,14
            if(savefile_p(i:i).eq.' ') savefile_p(i:i)='0'
          enddo
          call cfopen (nout(7),savefile_p,'new','formatted')
          write(*,*) 'Writing extra b2plasmf:',itim,itim_plas,
     &                plasinc,mod(itim_plas,plasinc),tim
          write(*,*) savefile_p
          call cfverw(nout(7), newversion)
          idum(0) = nCv
          idum(1) = nFc
          idum(2) = ns
          call cfwuin (nout(7), 3, idum, 'nCv,nFc,ns')
          call write_b2fplasma (nout(7), nCv, nFc, ns, state)
          close(nout(7))
        endif
!djm Jun2017
#ifndef FIXED_POINT
        if((balance_netcdf.ne.0).and.(ntim.gt.0)) then
          call update_average_balance(state,itim,ntim)
        endif
#endif
#endif
*    ..increment
        itim = itim+1
        itim_plas = itim_plas + 1
        if (no_solve.le.0) then
          tim = tim+dtim
        endif
        write(*,"(1x,a,i9,a,es14.6,a,i9,a,es10.2)") 'ITER ',itim,
     &    ' TIME ',tim,' NTIM ',ntim,' DTIM ',dtim
#ifndef FIXED_POINT
        write(*,*) 'delta na ',maxval(dabs((state%pl%na-state%psnc%na)/
     &    (state%pl%na+na_eps)))
        write(*,*) 'delta te ',maxval(dabs((state%pl%te-state%psnc%te)/
     &    (state%pl%te+te_eps*ev)))
        write(*,*) 'delta ti ',maxval(dabs((state%pl%ti-state%psnc%ti)/
     &    (state%pl%ti+ti_eps*ev)))
        write(*,*) 'delta tn ',maxval(dabs((state%pl%tn-state%psnc%tn)/
     &    (state%pl%tn+tn_eps*ev)))
        write(*,*) 'delta po ',maxval(dabs((state%pl%po-state%psnc%po)/
     &    (dabs(state%pl%po)+po_eps)))
        write(*,*) 'delta kt ',maxval(dabs((state%pl%kt-state%psnc%kt)/
     &    (dabs(state%pl%kt)+kt_eps*ev)))
        write(*,*) 'delta zt ',maxval(dabs((state%pl%zt-state%psnc%zt)/
     &    (dabs(state%pl%zt)+zt_eps*ev)))
        write(*,*) 'delta ua ',
     *   maxval(dabs((state%pl%na*state%pl%ua-
     &    state%psnc%na*state%psnc%ua)/
     &   (abs(state%pl%na*state%pl%ua)+na_eps*ua_eps)))
        if(delta_min.gt.0.0_R8.and.delta_max.gt.0.0_R8) then
          max_delta=max(
     &     maxval(dabs((state%pl%na-state%psnc%na)/
     &    (state%pl%na+na_eps))),
     &     maxval(dabs((state%pl%te-state%psnc%te)/
     &    (state%pl%te+te_eps*ev))),
     &     maxval(dabs((state%pl%ti-state%psnc%ti)/
     &    (state%pl%ti+ti_eps*ev))),
     &     maxval(dabs((state%pl%tn-state%psnc%tn)/
     &    (state%pl%tn+tn_eps*ev))),
     &     maxval(dabs((state%pl%po-state%psnc%po)/
     &    (abs(state%pl%po)+po_eps))),
     &     maxval(dabs((state%pl%kt-state%psnc%kt)/
     &    (state%pl%kt+kt_eps*ev))),
     &     maxval(dabs((state%pl%zt-state%psnc%zt)/
     &    (state%pl%zt+zt_eps*ev))),
     &     maxval(dabs((state%pl%na*state%pl%ua-
     &    state%psnc%na*state%psnc%ua)/
     &   (abs(state%pl%na*state%pl%ua)+na_eps*ua_eps))))
          if(max_delta.gt.delta_max) then
            write(*,*) 'b2mndr_delta: max_delta > delta_max'
            dtim=max(dt_min,dtim*dt_change_dec)
            write(*,*) 'b2mndr_dtim: dtim ',dtim
          elseif(max_delta.lt.delta_min) then
            write(*,*) 'b2mndr_delta: max_delta < delta_min'
            dtim=min(dt_max,dtim*dt_change_inc)
            write(*,*) 'b2mndr_dtim: dtim ',dtim
          else
            write(*,*) 'b2mndr_delta: max_delta OK'
            write(*,*) 'b2mndr_dtim: dtim ',dtim
          endif
        endif
        do while(itim.eq.ntim.and.stack_ptr.gt.0)
          write(*,'(a,i8,1p,g14.7,i8,i3)')
     .     'b2mndr_01:itim,dtim,ntim,stack_ptr',
     1                itim,dtim,ntim,stack_ptr
          ntim=stack_ntim(stack_ptr)+itim-1     ! dpc 2000.06.24 was stack_ntim(stack_ptr)
          dtim=stack_dtim(stack_ptr)
          stack_ptr=stack_ptr-1
          write(*,'(a,i8,1p,g14.7,i8,i3)')
     .     'b2mndr_02:itim,dtim,ntim,stack_ptr',
     1                itim,dtim,ntim,stack_ptr
        end do
*    ..produce time-dependent file
        if (switch%save_f31.lt.0) then
          write(fort31name,'(a,i8.8)') fort_lc//'31.', itim
          write(*,*) 'Saving ', fort31name, ' at time ', tim
          open(31,file=fort31name)
          call copy_background(nCv, nFc, ns, switch, mpg, geo,
     .                         state%pl, state%dv, state_ext)
          call write_f31(mpg)
          close(31)
        end if

#ifdef JET
        if(tim.ge.(save_tran_time-dtim/2.0_R8).and.
     &     delta_tran_time.gt.0.0_R8) then
          tran_num=tran_num+1
          if(tran_num.le.99) then
            write(tran_file_name,'(a4,i2)') 'tran',tran_num
            if(tran_file_name(5:5).eq.' ') tran_file_name(5:5)='0'
            call cfopen (nout(11),tran_file_name,'new','unformatted')
            call set_b2mod_ppout
            call write_b2mod_ppout(nout(11))
            close(nout(11))
            write(*,*) 'Saved ',tran_file_name
          else
            write(*,*) 'More than 99 tran files --- no more saved'
          endif
          save_tran_time=save_tran_time+delta_tran_time
        endif
#endif
c new CDF movie option
        if(tim.ge.(save_cdfmovie_time-dtim/2.0_R8).and.
     1    delta_cdfmovie_time.gt.0.0_R8) then
#ifndef NO_CDF
          call cdfmovie(ncid,nCv,mpg%nCi,ns,geo,state,switch)
#else
          write(*,*) 'Compile b2mndr without NO_CDF option for movies!'
#endif
          save_cdfmovie_time=save_cdfmovie_time+delta_cdfmovie_time
        endif
#endif
        lwti = .false.
        lwav = .false.
        if(b2time.gt.0) then
          if(mod(itim,b2time).eq.0.or.itim.eq.1) lwti = .true.
        endif
        if(ntim_batch.gt.0) then
          if(mod(itim,ntim_batch).eq.0) lwav = .true.
        endif
        call b2mwti(itim,tim,
#ifndef NO_CDF
     .              ntim,b2time,ntim_batch,
#endif
     .              nCv,nFc,ns,max(1,mpg%nXpt),geo,mpg,switch,
     .              state%pl,state%dv,state%co,state%rt,state%srw,
     .              state_ext,ismain,ismain0,lwti,lwav,.true.)
        if ((ibatch_av_all > 0) .and. lwav) then
          edition_batch = edition_batch + 1
          write (batch_name,'(a18,i4.4)')
     .           'batch_av/batch_av.',edition_batch
          call batch_av_all_save (batch_name,nCv,ns)
          write(*,*) 'Saved ',batch_name
        end if
        if ((iav_run > 0) .and. (ntim_run > 0)) then
          lrav = mod(itim,ntim_run) == 0
          write (*,*) ' itim,ntim_batch,lrav ', itim,ntim_run,lrav
          if (lrav) then
             edition_run = edition_run + 1
             write (run_av_name,'(a14,i4.4)')
     .              'run_av/run_av.',edition_run
             call run_av_save (run_av_name,nCv,ns,.true.,.false.,
     .                         state_avg)
             write(*,*) 'Saved ',run_av_name
          end if
        end if
#ifndef FIXED_POINT
*    ..tallies
        if(tally.gt.0) then
          if(mod(itim,tally).eq.0.or.itim.eq.1) then
            call tallies(ns,tim,mpg,geo,state,state_ext,switch)
          endif
        endif
        if (first_time_step) then
          if(.not.allocated(iz))
     >     call alloc_b2mod_diag(nCv,nFc,ns,mpg%nnreg(0),nstra,natm,
     .                           switch)
          call b2trci(geo,mpg,ns,nlim,nmol,switch)
          if (write_nml_user) call write_b2mod_user_namelist
        endif
        if (ank_tracing.ne.0.and.ntim.ne.0) then
          if (mod(itim,ank_tracing).eq.1.or.ank_tracing.eq.1)
     >     call b2trace(ns,nxtl,nxtr,ismain,dtim,
     .                  geo,mpg,state,switch,state_ext)
        endif
*    ..loop
#endif
        inquire(file='_quit',exist=quitexist_)                           !srv 18.05.09 13.04.11
        inquire(file='.quit',exist=quitexist)
#ifdef _OPENMP
        cpuval = omp_get_wtime()
#else
        call cpu_time(cpuval)
#endif
        elapsedval=epoch_seconds()
        if (b2mndr_cpu.gt.0.0_R8) then
          write(*,'(1x,3(a,es11.3))') 'b2-step-cpu =',cpuval-cpustart,
     &     ' elapsed ', elapsedval-elapsedstart,
     &     ' estimated remaining CPU time (s) =',
     &            min((ntim-itim)*real(cpuval-cpustart,R8),
     &      max(0.0_R8,b2mndr_cpu-real(cpuval-cpuinit,R8)))              !srv 18.05.09
        else if (b2mndr_elapsed.gt.0.0_R8) then
          write(*,'(1x,3(a,es11.3))') 'b2-step-cpu =',cpuval-cpustart,
     &     ' elapsed ', elapsedval-elapsedstart,
     &     ' estimated remaining elapsed time (s) =',
     &                min((ntim-itim)*(elapsedval-elapsedstart),
     &      max(0.0_R8,b2mndr_elapsed-(elapsedval-elapsedinit)))         !srv 18.05.09
        else if (etim.ne.0.0_R8 .and. etim.gt.switch%b2mndr_stim) then
          write(*,'(1x,3(a,es11.3))') 'b2-step-cpu =',cpuval-cpustart,
     &     ' elapsed ', elapsedval-elapsedstart,
     &     ' estimated remaining elapsed time (s) =',
     &                 min((ntim-itim)*(elapsedval-elapsedstart),
     &      max(0.0_R8,(etim-tim)/dtim-(elapsedval-elapsedstart)))       !srv 18.05.09
        else
          write(*,'(1x,3(a,es11.3))') 'b2-step-cpu =',cpuval-cpustart,
     &     ' elapsed ', elapsedval-elapsedstart,
     &     ' estimated remaining elapsed time (s) =',
     &                 (ntim-itim)*(elapsedval-elapsedstart)
        end if
        res_max = 0.0_R8
        do is = 0, ns-1
          if (any(solveco(is,:)))
     &       res_max = max(res_max,state%diag%aresco(is))
          if (any(solvemo(is,:)))
     &      res_max = max(res_max,state%diag%aresmo(is))
        end do
        if (any(solvepo))
     &    res_max = max(res_max,state%diag%arespo)
        if (any(solveee))
     &    res_max = max(res_max,state%diag%areshe)
        if (any(solveei))
     &    res_max = max(res_max,state%diag%areshi)
        if (any(solveen).and.switch%tn_style.eq.2)
     &    res_max = max(res_max,state%diag%areshn)
        if (any(solvekt).and.switch%solve_keps.ge.1)
     &    res_max = max(res_max,state%diag%areskt)
        if (any(solvezt).and.switch%solve_keps.eq.2)
     &    res_max = max(res_max,state%diag%areszt)
#ifndef FIXED_POINT
        if(test_residual.and.mod(itim,10).eq.0.and.no_solve.le.0) then   !jwk
          quit_residual=                                                 !jwk
     &        (state%diag%areshe*min_areshe).le.(min_areshe*min_areshe)  !jwk
     &        .and.
     &        (state%diag%areshi*min_areshi).le.(min_areshi*min_areshi)  !jwk
     &        .and.
     &        (state%diag%aresco(ismain)*
     &                           min_aresco).le.(min_aresco*min_aresco)  !jwk
        endif                                                            !jwk
#endif
        quit=quitexist_.or.quitexist.or.                                 !jwk !srv 13.04.11
     1      (cpuval-cpuinit.gt.b2mndr_cpu.and.b2mndr_cpu.gt.0.0_R8).or.
     1      (elapsedval-elapsedinit.gt.b2mndr_elapsed.and.
     &       b2mndr_elapsed.gt.0.0_R8)
#ifndef FIXED_POINT
     &      .or. quit_residual.or.
     2      (etim.gt.switch%b2mndr_stim.and.tim.ge.etim.and.
     &       etim.ne.0.0_R8).or.
     3      (no_solve.eq.1.and.itim.eq.ntim)
#ifdef _OPENMP
        write_save=(cpuval-cpuinit).gt.cputarget.and.
     &              cpuincrement.gt.0.0_R8
#else
        write_save=cpuval.gt.cputarget.and.cpuincrement.gt.0.0_R8
#endif
        write_save=write_save.or.
     .   (tim.ge.(save_plasma_time-dtim/2.0_R8).and.
     .    delta_plasma_time.gt.0.0_R8)
        if(ntim_save.gt.0.and.itim.gt.0) write_save=write_save.or.
     .   mod(itim,ntim_save).eq.0
        write_save=write_save.and..not.quit
        if(write_save) then
          write(savefile(13:16),'(i4)') edition
          edition=edition+1
          if(cpuval.gt.cputarget.and.cpuincrement.gt.0.0_R8)
     .     cputarget=cputarget+cpuincrement
          if(tim.ge.(save_plasma_time-dtim/2.0_R8).and.
     .       delta_plasma_time.gt.0.0_R8)
     .     save_plasma_time=save_plasma_time+delta_plasma_time
          if(edition.gt.9999) then
            edition=0
            savefilestatus='old'
          endif
          do i=13,16
            if(savefile(i:i).eq.' ') savefile(i:i)='0'
          enddo
          call cfopen (nout(6), savefile, savefilestatus, 'unformatted')
          call cfverw (nout(6), newversion)                          !xpb/dpc
*       ..test plasma state
          call b2xvps (nCv, nFc, ns, state % pl, state % dv)
          call b2xvsg (nCv*ns, state % dv % kinrgy, 1, 'kinrgy', '.ge.')
*       ..write plasma state
          call write_b2fstate (nout(6), nCv, nFc, ns, state)
          close(nout(6))
*       ..write additional save namelist data
          filename='b2.feedback_save.parameters'
          call find_file(filename,found)
          if (found) call write_b2us_feedback(99,savefile(13:16))
          filename='b2.neutrals_save.parameters'
          call find_file(filename,found)
          if (found)
     >     call write_b2mod_neutrals_save(99,savefile(13:16))
#ifdef WG_TODO
          filename='b2.sputter_save.parameters'
          call find_file(filename,found)
          if (found) call write_b2mod_sputter_namelist(99,nx,ny,ns,
     >     savefile(13:16))
          filename='b2.wall_save.parameters'
          call find_file(filename,found)
          if (found) call write_b2mod_wall_namelist(99,nx,ny,ns,
     >     savefile(13:16))
#endif
          filename='b2.transport_models_save.parameters'
          call find_file(filename,found)
          if (found) call write_b2mod_transport_models(99,
     >     savefile(13:16))
          if (ank_tracing.ne.0) call b2filef

#ifdef B25_EIRENE
          if (switch%use_eirene.ne.0) then
            call write_f44(savefile(13:16))
            filename=fort_lc//'13'
            inquire(file=filename,exist=found)
            if (found) then
#ifdef WINDOWS
              systemarg='copy '//trim(filename)//' '//
     &                           trim(filename)//'.'//savefile(13:16)//
     &        ' > nul'
#else
              systemarg='cp '//trim(filename)//' '//
     &                         trim(filename)//'.'//savefile(13:16)
#endif
#ifdef NAGFOR
              call system(systemarg, status, ierror)
#else
              call system(systemarg)
#endif
            endif
            filename=fort_lc//'14'
            inquire(file=filename,exist=found)
            if (found) then
#ifdef WINDOWS
              systemarg='copy '//trim(filename)//' '//
     &                           trim(filename)//'.'//savefile(13:16)//
     &        ' > nul'
#else
              systemarg='cp '//trim(filename)//' '//
     &                         trim(filename)//'.'//savefile(13:16)
#endif
#ifdef NAGFOR
              call system(systemarg, status, ierror)
#else
              call system(systemarg)
#endif
            endif
            filename=fort_lc//'15'
            inquire(file=filename,exist=found)
            if (found) then
#ifdef WINDOWS
              systemarg='copy '//trim(filename)//' '//
     &                           trim(filename)//'.'//savefile(13:16)//
     &        ' > nul'
#else
              systemarg='cp '//trim(filename)//' '//
     &                         trim(filename)//'.'//savefile(13:16)
#endif
#ifdef NAGFOR
              call system(systemarg, status, ierror)
#else
              call system(systemarg)
#endif
            endif
          endif
#endif
          write(*,*) 'Plasmastate savefile ',savefile,' saved'
          if (iav_run > 0) then
            write (aver_name,'(a17,a4)')
     .              'plasmastate_aver.',savefile(13:16)
            call run_av_save (aver_name,nCv,ns,.true.,.true.,state_avg)
            write (*,*) aver_name,' saved '
          end if
        endif
        ids_replace = .false.
        write_ids = ids_save.gt.0 .and. ntim.gt.0
        if (write_ids) write_ids = mod(itim,ids_save).eq.0
        write_ids = write_ids .or.
     .   (tim.ge.(save_ids_time-dtim/2.0_R8).and.
     .    delta_ids_time.gt.0.0_R8)
        process_ids = write_ids
        write_av = ids_av.gt.0 .and. ntim.gt.0 .and. ntim_batch.gt.0
        if (write_av) write_av = mod(itim,ids_av*ntim_batch).eq.0
        process_av = write_av
        if (write_ids .or. write_av) then
#ifdef IMAS
          status = 0
          if (idx.eq.0.and.continued) then
#if AL_MAJOR_VERSION > 4
            uri =
     .       'imas:mdsplus?path='//trim(ids_path)
            write(*,'(2a)')
     .        "Checking if IMAS data entry already exists : ",
     .         trim(ids_path)
            call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
            write(*,'(2a,2i8)')
     .       "Checking if IMAS data entry already exists : ",
     .        trim(database), shot, run
            call imas_open_env(treename, shot, run, idx,
     .       username, database, version, status)
#endif
            if (status.ne.0) then
              if (database.eq.'ITER'.or.
     .            index(ids_path,'imasdb/ITER').gt.0) then
                write(*,*) "Did not find ITER database IDS file."
                write(*,*) "Checking if old ''iter'' case exists."
#if AL_MAJOR_VERSION > 4
                l=index(ids_path,'imasdb/ITER')
                write(olddir,'(a)')
     .           ids_path(1:l-1)//'imasdb/iter'//ids_path(l+11:256)
                uri = 'imas:mdsplus?path='//trim(olddir)
                call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
                call imas_open_env(treename, shot, run, idx,
     .           username, 'iter', version, status)
#endif
                if (status.eq.0) then
                  database = 'iter'
                  write(*,*) "Old database case found."
                  write(*,*) "Will be rewritten in new location."
                end if
              else if (database.eq.'iter'.or.
     .                 index(ids_path,'imasdb/iter').gt.0) then
#if AL_MAJOR_VERSION > 4
                l=index(ids_path,'imasdb/iter')
                write(olddir,'(a)')
     .           ids_path(1:l-1)//'imasdb/ITER'//ids_path(l+11:256)
                uri = 'imas:mdsplus?path='//trim(olddir)
                call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
                call imas_open_env(treename, shot, run, idx,
     .           username, 'ITER', version, status)
#endif
                database = 'ITER'
              end if
            end if
            !! If this is a time continuation run, append the new data to the IDS
            if (status.eq.0 .and. idx.ne.0) then
              write(*,*) "Reading old IMAS data entry: ",
     .         trim(database), shot, run
              call ids_get( idx, "equilibrium", equilibrium, status)
              call ids_get( idx, "edge_profiles", old_edge_profiles,
     .          status)
              if (status.ne.0) then
                write(*,*) 'Error opening old edge_profiles IDS ! '//
     .                     'Will create a new one.'
                idx = 0
                continued = .false.
              else
                num_time_slices = size(old_edge_profiles%time)
                if (num_time_slices.gt.0) then
                  ids_end_time = old_edge_profiles%time(num_time_slices)
                else
                  ids_end_time = IDS_REAL_INVALID
                end if
                call ids_deallocate( old_edge_profiles )
                call ids_get( idx, "edge_profiles/1",
     .            old_batch_profiles, status)
                if (status.ne.0) then
                  num_batch_slices = 0
                else
                  num_batch_slices = size(old_batch_profiles%time)
                  call ids_deallocate( old_batch_profiles )
                end if
                old_start_time = 0.0_IDS_real
                old_end_time = IDS_REAL_INVALID
                old_imas_version = 'x.xx.x'
                call ids_get( idx, "dataset_description",
     .                                  old_description, status)
                if ( status.ne.0 ) then
                  write(*,*)
     .              'Error opening old dataset_description IDS !'
                else if (associated(old_description%ids_properties%
     &                              version_put%data_dictionary)) then
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
                  old_start_time = description%simulation%time_begin
                  old_end_time = description%simulation%time_end
#endif
                  old_imas_version = old_description%ids_properties%
     &                               version_put%data_dictionary(1)
                  call ids_deallocate( old_description )
                end if
                continued = continued .and.
     &           (ids_end_time.lt.tim .and.
     &            ids_end_time.ne.IDS_REAL_INVALID)
                continued = continued .or.
     &            run_start_time.ge.ids_end_time
                if (continued.or.database.eq.'iter'.or.
     &              index(ids_path,'imasdb/iter').gt.0) then
                  if (.not.streql(old_imas_version,imas_version).or.
     &                database.eq.'iter'.or.
     &                index(ids_path,'imasdb/iter').gt.0) then
                    if (.not.streql(old_imas_version,imas_version)) then
                      write(*,*)
     &                 'Old IMAS data entry was written '//
     &                 'using Data Dictionary version '//
     &                  trim(old_imas_version)//'.'
                      write(*,*)
     &                 'Recreating using Data Dictionary version '//
     &                  trim(imas_version)//'.'
                    end if
                    if (database.eq.'iter'.or.
     &                  index(ids_path,'imasdb/iter').gt.0) write(*,*)
     &               'IDS file will be moved to ITER database.'
                    call close_ual(idx)
cxpb Copy the IDS to a temporary location with the new DD and then bring it back
                    tmp_run = run
                    if (database.ne.'iter'.and.
     &                  index(ids_path,'imasdb/iter').eq.0)
     &               tmp_run = run + 1000
#if AL_MAJOR_VERSION > 4
                    call al_build_uri_from_legacy_parameters
     &               ( MDSPLUS_BACKEND, shot, run, trim(username),
     &                 trim(database), int2str(IMAS_MAJOR_VERSION), '',
     &                 uri_source, status )
                    if (database.eq.'iter'.or.
     &                  index(ids_path,'imasdb/iter').gt.0) then
                      call al_build_uri_from_legacy_parameters
     &                 ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),
     &                   'ITER', int2str(IMAS_MAJOR_VERSION), '',
     &                   uri_dest, status )
                    else
                      call al_build_uri_from_legacy_parameters
     &                 ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),
     &                   trim(database), int2str(IMAS_MAJOR_VERSION), '',
     &                   uri_dest, status )
                    end if
#if IMAS_MAJOR_VERSION < 4
                    write(systemarg,'(a,a,a,a,a)')
     &                'idscp --set-dataset-version'//
     &                     ' -s ',trim(uri_source),
     &                     ' -d ',trim(uri_dest),' --dd-update'
#else
                    write(systemarg,'(a,a,a,a,a)')
     &                'idscp -s ',trim(uri_source),
     &                     ' -d ',trim(uri_dest),' --dd-update'
#endif
#else
#if ( IMAS_MINOR_VERSION > 31 || IMAS_MAJOR_VERSION > 3 )
                    write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)')
     &                'idscp --setDatasetVersion'//
     &                     ' -si ',shot,' -ri ',run,
     &                     ' -so ',shot,' -ro ',tmp_run,
     &                     ' -d ',trim(database),' -u ',trim(username)
#else
                    write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)')
     &                'idscp -si ',shot,' -ri ',run,
     &                     ' -so ',shot,' -ro ',tmp_run,
     &                     ' -d ',trim(database),' -u ',trim(username)
#endif
                    if (database.eq.'iter'.or.
     &                  index(ids_path,'imasdb/iter').gt.0)
     &               systemarg = trim(systemarg)//' -do ITER'
#endif
#ifdef NAGFOR
                    call system(systemarg, status, ierror)
#else
                    call system(systemarg)
#endif
                    if (database.ne.'iter'.and.
     &                  index(ids_path,'imasdb/iter').eq.0) then
#if AL_MAJOR_VERSION > 4
                      call al_build_uri_from_legacy_parameters
     &                 ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),
     &                   trim(database), int2str(IMAS_MAJOR_VERSION), '',
     &                   uri_source, status )
                      call al_build_uri_from_legacy_parameters
     &                 ( MDSPLUS_BACKEND, shot, run, trim(username),
     &                   trim(database), int2str(IMAS_MAJOR_VERSION), '',
     &                   uri_dest, status )
                      write(systemarg,'(a,a,a,a)')
     &                  'idscp -s ',trim(uri_source),
     &                       ' -d ',trim(uri_dest)
#else
#if ( IMAS_MINOR_VERSION > 31 || IMAS_MAJOR_VERSION > 3 )
                      write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)')
     &                 'idscp --setDatasetVersion'//
     &                    ' -si ',shot,' -ri ',tmp_run,
     &                    ' -so ',shot,' -ro ',run,
     &                    ' -d ',trim(database),' -u ',trim(username)
#else
                      write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)')
     &                 'idscp -si ',shot,' -ri ',tmp_run,
     &                    ' -so ',shot,' -ro ',run,
     &                    ' -d ',trim(database),' -u ',trim(username)
#endif
#endif
#ifdef NAGFOR
                      call system(systemarg, status, ierror)
#else
                      call system(systemarg)
#endif
                    end if
                    if (database.eq.'iter') database = 'ITER'
                    if (index(ids_path,'imasdb/iter').gt.0) then
                      l=index(ids_path,'imasdb/iter')
                      write(ids_path(l+7:l+10),'(a4)') 'ITER'
                    end if
#if AL_MAJOR_VERSION > 4
                    uri = 'imas:mdsplus?path='//trim(ids_path)
                    call imas_open( uri, OPEN_PULSE,
     .               idx, status, message )
#else
                    call imas_open_env(treename, shot, run, idx,
     .               username, database, version, status)
#endif
                  end if
                  if (continued) then
                    if (process_ids) then
                      write(*,*)
     .                 "Appending a new time slice at t = ", tim, " s."
                      num_time_slices = num_time_slices + 1
                    end if
                    if (process_av) then
                      write(*,*)
     .                 "Appending a new averaged time slice at t = ",
     .                  tim, " s."
                      num_batch_slices = num_batch_slices + 1
                    end if
                  end if
                  if (process_ids) then
                    time_slice_index = num_time_slices
                    call B25_process_ids(
     &               geo, mpg, state, state_ext, state_avg, switch,
     &               edge_profiles, edge_sources, edge_transport,
     &               radiation, description, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &               summary,
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
     &               numerics,
#endif
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
     &               old_start_time, run_end_time,
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
     &               divertors,
#endif
     &               tim, dtim, shot, run, database, version,
     &               new_eq_ggd,
     &               time_slice_index, num_time_slices )
                    process_ids = .false.
                  end if
                  if (process_av) then
                    batch_index = num_batch_slices
                    call B25_av_ids( summary_av,
     &               geo, mpg, state, state_ext, state_avg, switch,
     &               batch_profiles, batch_sources,
     &               description, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &               summary,
#endif
     &               tim, shot, run, database, version, new_eq_ggd,
     &               batch_index, num_batch_slices )
                    process_av = .false.
                  end if
                else
                  write(*,*)
     .             'Not a time continuation, edge IDS components '//
     .             'will be replaced !'
                  ids_replace = .true.
                end if
              end if
            else
              write(*,*)
     .         "No previous IMAS data entry found, "//
     .         "a new one will be created"
              idx = 0
              if (database.eq.'iter') database = 'ITER'
              if (index(ids_path,'imasdb/iter').gt.0) then
                l=index(ids_path,'imasdb/iter')
                write(ids_path(l+7:l+10),'(a4)') 'ITER'
              end if
            end if
          end if
          if ((process_ids .or. process_av) .and.
     &        (status.ne.0 .or. idx.eq.0)) then
            if (status.ne.0) idx = 0
            if (process_ids) then
              call B25_process_ids(
     &          geo, mpg, state, state_ext, state_avg, switch,
     &          edge_profiles, edge_sources, edge_transport,
     &          radiation, description, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &          summary,
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
     &          numerics,
#endif
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
     &          run_start_time, run_end_time,
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
     &          divertors,
#endif
     &          tim, dtim, shot, run, database, version, new_eq_ggd )
              old_start_time = run_start_time
              process_ids = .false.
            end if
            if (process_av) then
              call B25_av_ids( summary_av,
     &         geo, mpg, state, state_ext, state_avg, switch,
     &         batch_profiles, batch_sources,
     &         description, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &         summary,
#endif
     &         tim, shot, run, database, version, new_eq_ggd )
              process_av = .false.
            end if
          else if (process_ids .or. process_av) then
            if (process_ids) then
              write(*,*)
     .         "Appending a new time slice at t = ", tim, " s."
              num_time_slices = num_time_slices + 1
              time_slice_index = num_time_slices
              call B25_process_ids(
     &          geo, mpg, state, state_ext, state_avg, switch,
     &          edge_profiles, edge_sources, edge_transport,
     &          radiation, description, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &          summary,
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
     &          numerics,
#endif
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
     &          old_start_time, run_end_time,
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
     &          divertors,
#endif
     &          tim, dtim, shot, run, database, version, new_eq_ggd,
     &          time_slice_index, num_time_slices )
              process_ids = .false.
            end if
            if (process_av) then
              write(*,*)
     .         "Appending a new averaged time slice at t = ",
     .          tim, " s."
              num_batch_slices = num_batch_slices + 1
              batch_index = num_batch_slices
              call B25_av_ids( summary_av,
     &          geo, mpg, state, state_ext, state_avg, switch,
     &          batch_profiles, batch_sources,
     &          description, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &          summary,
#endif
     &          tim, shot, run, database, version, new_eq_ggd,
     &          batch_index, num_batch_slices )
              process_av = .false.
            end if
          end if
          if (ids_replace) then
            if (write_ids)
     >       call new_ids_edge( edge_profiles, edge_sources,
     &        edge_transport, radiation, description, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &        summary,
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
     &        numerics,
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
     &        divertors,
#endif
#if IMAS_MAJOR_VERSION > 3
     &        uri,
#endif
     &        idx, new_eq_ggd )
            new_eq_ggd = new_eq_ggd .and. .not.write_ids
            if (write_av)
     &       call new_batch_edge( idx,
#if IMAS_MAJOR_VERSION > 3
     &        uri,
#endif
     &        batch_profiles, batch_sources, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &        summary,
#endif
     &        description, summary_av, new_eq_ggd )
          else
            if (database.eq.'iter') database = 'ITER'
            if (write_ids)
     >       call put_ids_edge( edge_profiles, edge_sources,
     &        edge_transport, radiation, description, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &        summary,
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
     &        numerics,
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
     &        divertors,
#endif
#if AL_MAJOR_VERSION > 4
     &        ids_path,
#else
     &        treename, shot, run, username, database, version,
#endif
     &        idx, new_eq_ggd )
            new_eq_ggd = new_eq_ggd .and. .not.write_ids
            if (write_av)
     &       call put_batch_edge(
#if AL_MAJOR_VERSION > 4
     &        ids_path,
#else
     &        treename, shot, run, username, database, version,
#endif
     &        idx, batch_profiles, batch_sources, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &        summary,
#endif
     &        description, summary_av, new_eq_ggd )
          end if
#endif
          if (.not.continued) continued = .true.
          if (tim.ge.(save_ids_time-dtim/2.0_R8).and.
     .        delta_ids_time.gt.0.0_R8)
     .     save_ids_time=save_ids_time+delta_ids_time
        end if
        first_time_step=.false.

#ifdef CATALYST
        if (use_catalyst .eqv. .true.) then
           call coprocessor(nVx, nFc, nCv, itim, tim, geo, state)
        endif
#endif
#endif
#ifdef FIXED_POINT
        write(*,*) 'MAX RESIDUAL ',res_max
        primal_iterations = itim
*   ..end loop
      end do
*   ..call cost function
      call b2usr_cost_function(nCv, nFc, nVx, ns, geo, mpg, state,
     &    state_ext, switch, J)
      do icf = 1, ncf
        write (ss,'(I1)') icf
        if (icf.gt.9) write (ss,'(I2)') icf
        write(*,*) 'Cost function value '//ss//': ',J(icf)
      end do
#else
        goto 1
*   ..end loop
      endif  !}
      primal_iterations = itim
#endif
      call subend ()
      return
      end subroutine b2mndr_1

      subroutine b2mndr_2(nout, ns, switch, geo, mpg,
     &                    state, state_ext, state_avg)
      use b2mod_switches
      use b2mod_ad
     & , only : old_erosion, old_deposition
      implicit none
#ifdef JET
      integer nout(0:11), ns, idum(0:9)
#else
      integer nout(0:10), ns, idum(0:9)
#endif
#ifdef WG_TODO
      integer nx, ny
#endif
      type(switches), intent (inout)   :: switch
      type(mapping), intent (inout)    :: mpg
      type(geometry), intent (inout)   :: geo
      type(B2State), intent (inout)    :: state
      type(B2StateExt), intent (inout) :: state_ext
      type(B2Average), intent (inout)  :: state_avg
      integer nCv, nFc, is                                          !srv 29.11.18
#ifndef NO_CDF
      integer nVx
#endif
#ifdef WG_TODO
      character*256 filename
      logical found
#endif
      logical ids_done
#ifdef IMAS
      integer l, tmp_run, time_slice_index, status
      character*512 systemarg
#if AL_MAJOR_VERSION > 4
      character(len=:), allocatable :: message
      character(len=STRMAXLEN) :: uri, uri_source, uri_dest
      character(len=256) :: olddir
#endif
      logical streql
      external streql
      intrinsic index
#endif

      call subini ('b2mndr_2')
      if(quitexist)
     1 write(*,*) 'stopping because .quit seen'
      if(quitexist_) write(*,*) 'stopping because _quit seen'            !srv 13.04.11
      if(cpuval.gt.b2mndr_cpu.and.b2mndr_cpu.gt.0.0_R8)
     1 write(*,'(1x,a,1p,2(e12.5,1x))')
     . 'stopping because cpu > limit: ',cpuval,b2mndr_cpu
      if(itim.ge.ntim) write(*,'(1x,a,2(i8,1x))')
     1 'stopping because itim >= ntim: ',itim,ntim
      if(etim.gt.switch%b2mndr_stim.and.tim.ge.etim.and.etim.ne.0.0_R8)
     1 write(*,'(1x,a,1p,2(e12.5,1x))')
     . 'stopping because tim >= etim: ',tim,etim
      if(etim.eq.switch%b2mndr_stim.and.switch%b2mndr_stim.gt.0.0_R8)
     & write(*,'(1x,a,1p,2(e12.5,1x))')
     1 'stopping because etim = stim: ',etim,switch%b2mndr_stim
      if(elapsedval-elapsedinit.gt.b2mndr_elapsed.and.
     &   b2mndr_elapsed.gt.0.0_R8) write(*,'(1x,a,1p,2(e12.5,1x))')
     1 'stopping because elapsed CPU time > limit: ',
     .  elapsedval-elapsedinit,b2mndr_elapsed
      if(ierr.ne.0)
     1 write(*,*) 'stopping because ierr <> 0: ',ierr
      if(quit_residual)                                                  !jwk
     1 write(*,*) 'stopping because min_residual criteria satisfied'     !jwk
#ifdef FIXED_POINT
      if(res_max.le.res_quit)
     1 write(*,*) 'stopping because res_quit criteria satisfied'
#endif

      nCv = mpg%nCv
      nFc = mpg%nFc
#ifndef NO_CDF
      nVx = mpg%nVx
#endif

#ifdef USE_IMPGYRO
* finalise impgyro !yamoto
      call impgyro_last !yamoto
#endif
* finalise feedback
      call write_b2us_feedback(99,'    ')
* finalise transport_models
      call write_b2mod_transport_models(99,'    ')
#ifdef WG_TODO
* finalise sputtering
      call write_b2mod_sputter_namelist(99, nx, ny, ns, '    ')
* finalise time-dependent wall data
      filename='b2.wall_save.parameters'
      call find_file(filename,found)
      call wall_final(nx,ny,ns,found)

* ..display and save the result of the calculation
*   ..write summary of final state
      call b2wfpi (nout(0), idout1, nx, ny, ns, ne, na, ua,
     &  te, ti, po, fna, fhe, fhi, fch)
      if (old_style.gt.0) call b2waos (nout(0), nx, ny, ns,              !pb 29.07.16
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
#endif
*   ..test plasma state
      call b2xvps (nCv, nFc, ns, state%pl, state%dv)
      call b2xvsg (nCv*ns, state%dv%kinrgy, 1, 'kinrgy', '.ge.')
*   ..save final state
      call write_b2fstate (nout(2), nCv, nFc, ns, state)
*   ..write plasma state
      call cfopen(nout(7),'b2fplasma','new','unformatted')
      call cfverw(nout(7), newversion)
      idum(0) = nCv
      idum(1) = nFc
      idum(2) = ns
      call cfwuin (nout(7), 3, idum, 'nCv,nFc,ns')
      call write_b2fplasma (nout(7), nCv, nFc, ns, state)
      close(nout(7))

*   ..mms steady-state error analysis
      if(switch%use_mms.eq.1.and.switch%iout_b2wdat.eq.10) then                 !srv 26.02.18 {
        call calc_err(nCv, ns, state%pl, switch, itim, mpg)
      endif                                                              !srv 26.02.18 }
      if(switch%b2mndr_iout.ne.0) then
        if(switch%use_mms.eq.0) then                                            !srv 29.11.18 {
         do is = 0,ns-1
          write (chns,'(i3.3)') is
          call my_out_us(70,nCv,0,state%pl%na(1,is),
     &      'b2mndr_na'//chns//'out')
          call my_out_us(70,nCv,0,state%pl%ua(1,is),
     &      'b2mndr_ua'//chns//'out')
         enddo
         call my_out_us(70,nCv,0,state%pl%te,'b2mndr_te_out')
         call my_out_us(70,nCv,0,state%pl%te/ev,'b2mndr_te_eV_out')
         call my_out_us(70,nCv,0,state%pl%ti,'b2mndr_ti_out')
         call my_out_us(70,nCv,0,state%pl%ti/ev,'b2mndr_ti_eV_out')
         call my_out_us(70,nCv,0,state%pl%po,'b2mndr_po_out')
       endif
      endif                                                              !srv 29.11.18 }
      call b2wdat (ns, mpg, geo, state)                                  !srv 05.10.01
#ifdef WG_TODO
      if (switch%use_astra.ne.0) then
       call cfopen (nout(8),'b2_to_astra.dat','replace','un*formatted')
       call write_b2_to_astra(nx, ny, ns, ismain, ismain0, nout(8))      !srv 28.07.11 19.09.11
       close(nout(8))
      end if
#endif
*   ..write out fort.30 and fort.31 files if not yet done
      if (switch%use_eirene.eq.0) then
!WG_TODO       if (switch%save_f30.ne.0) call write_f30(nx, ny)
       if (switch%save_f31.ne.0) then
        write(*,*) 'Saving '//trim(fort_lc)//'31 at end of run'
        open(31,file=trim(fort_lc)//'31')
        call copy_background(nCv, nFc, ns, switch, mpg, geo,
     &                       state%pl, state%dv, state_ext)
        call write_f31(mpg)
        close(31)
       end if
      end if
*   ..produce data for last frame in movie if not yet done
      if (mod(itim,mvinc).ne.0.and.(itim/mvinc+1).le.mvnum) then
        call b2mwmv (nout(3), nCv, state, itim, tim, dtim)
      endif
#ifdef WG_TODO
*   ..produce last time-dependent file if not yet done
#ifdef JET
      if(tim.gt.save_tran_time-delta_tran_time+dtim.and.
     &   delta_tran_time.gt.0.0_R8) then
        tran_num=tran_num+1
        if(tran_num.le.99) then
          write(tran_file_name,'(a4,i2)') 'tran',tran_num
          if(tran_file_name(5:5).eq.' ') tran_file_name(5:5)='0'
          call cfopen (nout(11),tran_file_name,'new','unformatted')
          call set_b2mod_ppout
          call write_b2mod_ppout(nout(11))
          close(nout(11))
          write(*,*) 'Saved ',tran_file_name
        else
          write(*,*) 'More than 99 tran files --- no more saved'
        endif
      endif
#endif
#endif
*   ..produce last CDF movie frame if not yet done
      if(tim.gt.save_cdfmovie_time-delta_cdfmovie_time+dtim.and.
     1  delta_cdfmovie_time.gt.0.0_R8) then
#ifndef NO_CDF
        call cdfmovie(ncid,nCv,mpg%nCi,ns,geo,state,switch)
#else
        write(*,*) 'Compile b2mndr without NO_CDF option for movies!'
#endif
      endif
*   ..produce 2d profiles if not yet done
      if(b2time.gt.0) then
        if((mod(itim,b2time).ne.0.and.itim.gt.1).or.itim.eq.0) then
          call b2mwti(itim,tim,
#ifndef NO_CDF
     *                ntim,b2time,ntim_batch,
#endif
     *                nCv,nFc,ns,max(1,mpg%nXpt),geo,mpg,switch,
     *                state%pl,state%dv,state%co,state%rt,state%srw,
     *                state_ext,ismain,ismain0,.true.,.false.,.false.)
        endif
        call dealloc_b2mod_mwti
      endif
*   ..write out last tallies if not yet done
      if(tally.gt.0) then
        if(mod(itim,tally).ne.0.and.itim.gt.1) then
          call tallies(ns,tim,mpg,geo,state,state_ext,switch)
        endif
      endif
#ifndef FIXED_POINT
*   ..create NetCDF file containing everything for balance postprocessing routines !djm Jan2017
      if(balance_netcdf.ne.0) then
#ifndef NO_CDF
        if (balance_average .ne. 0) then
          call update_average_balance(state,itim,-1) ! scale quantities by 1/itim if averaged
        endif
        call write_balance(nCv,nFc,nVx,ns,geo,mpg)
#else
        write(*,*)
     >   'NO_CDF is set at compilation, cannot write balance.nc'
#endif
      endif
#endif
*   ..write out last IDS if not yet done
      ids_done = .false.
      write_ids = .false.
      ids_replace = .false.
      if (ids_save.gt.0) then
        ids_done = mod(itim,ids_save).eq.0
        write_ids = mod(itim,ids_save).ne.0
      end if
      if (delta_ids_time.gt.0.0_R8) then
        ids_done = ids_done .or.
     .    tim.le.save_ids_time-delta_ids_time+dtim
        write_ids = write_ids .or.
     .    tim.gt.save_ids_time-delta_ids_time+dtim
      end if
      write_ids = write_ids .or. ids_save.lt.0
      write_ids = write_ids .and. ntim.gt.0 .and. .not.ids_done
      if (write_ids) then
#ifdef IMAS
        status = 0
        if (idx.ne.0) then
          write(*,*) "Appending a new time slice at t = ", tim, " s."
          num_time_slices = num_time_slices + 1
          time_slice_index = num_time_slices
          call B25_process_ids(
     &      geo, mpg, state, state_ext, state_avg, switch,
     &      edge_profiles, edge_sources, edge_transport,
     &      radiation, description, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &      summary,
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
     &      numerics,
#endif
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
     &      old_start_time, run_end_time,
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
     &      divertors,
#endif
     &      tim, dtim, shot, run, database, version, new_eq_ggd,
     &      time_slice_index, num_time_slices )
        else if (continued) then
#if AL_MAJOR_VERSION > 4
          uri = 'imas:mdsplus?path='//trim(ids_path)
          write(*,'(2a)')
     .     "Checking if IMAS data entry already exists : ",
     .      trim(ids_path)
          call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
          write(*,'(2a,2i8)')
     .     "Checking if IMAS data entry already exists : ",
     .      trim(database), shot, run
          call imas_open_env(treename, shot, run, idx,
     .      username, database, version, status)
#endif
          if (status.ne.0) then
            if (database.eq.'ITER'.or.
     .          index(ids_path,'imasdb/ITER').gt.0) then
              write(*,*) "Did not find ITER database IDS file."
              write(*,*) "Checking if old ''iter'' case exists."
#if AL_MAJOR_VERSION > 4
              l=index(ids_path,'imasdb/ITER')
              write(olddir,'(a)')
     .         ids_path(1:l-1)//'imasdb/iter'//ids_path(l+11:256)
              ids_path = trim(olddir)
              uri = 'imas:mdsplus?path='//trim(ids_path)
              call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
              call imas_open_env(treename, shot, run, idx,
     .         username, 'iter', version, status)
#endif
              if (status.eq.0) then
                database = 'iter'
                write(*,*) "Old database case found."
                write(*,*) "Will be rewritten in new location."
              end if
            else if (database.eq.'iter'.or.
     .               index(ids_path,'imasdb/iter').gt.0) then
#if AL_MAJOR_VERSION > 4
              l=index(ids_path,'imasdb/ITER')
              write(olddir,'(a)')
     .         ids_path(1:l-1)//'imasdb/iter'//ids_path(l+11:256)
              ids_path = trim(olddir)
              uri = 'imas:mdsplus?path='//trim(ids_path)
              call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
              call imas_open_env(treename, shot, run, idx,
     .         username, 'ITER', version, status)
#endif
              database = 'ITER'
            end if
          end if
          !! If this is a time continuation run, append the new data to the IDS
          if (status.eq.0 .and. idx.ne.0) then
            write(*,*) "Reading old IMAS data entry: ",
     .       trim(database), shot, run
            call ids_get( idx, "equilibrium", equilibrium, status)
            call ids_get( idx, "edge_profiles", old_edge_profiles,
     .        status)
            if (status.ne.0) then
              write(*,*) 'Error opening old edge_profiles IDS ! '//
     .                   'Will create a new one.'
              idx = 0
              continued = .false.
            else
              num_time_slices = size(old_edge_profiles%time)
              if (num_time_slices.gt.0) then
                ids_end_time = old_edge_profiles%time(num_time_slices)
              else
                ids_end_time = IDS_REAL_INVALID
              end if
              call ids_deallocate( old_edge_profiles )
              old_start_time = 0.0_IDS_real
              old_end_time = IDS_REAL_INVALID
              old_imas_version = 'x.xx.x'
              call ids_get( idx, "dataset_description",
     .                                old_description, status)
              if ( status.ne.0 ) then
                write(*,*) 'Error opening old dataset_description IDS !'
              else if (associated(old_description%ids_properties%
     &                            version_put%data_dictionary)) then
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
                old_start_time = description%simulation%time_begin
                old_end_time = description%simulation%time_end
#endif
                old_imas_version = old_description%ids_properties%
     &                             version_put%data_dictionary(1)
                call ids_deallocate( old_description )
              end if
              continued = continued .and.
     &         (ids_end_time.lt.tim .and.
     &          ids_end_time.ne.IDS_REAL_INVALID)
              continued = continued .or.
     &         (run_start_time.ge.ids_end_time .and.
     &          ids_end_time.ne.IDS_REAL_INVALID)
              if (continued.or.database.eq.'iter'.or.
     &            index(ids_path,'imasdb/iter').gt.0) then
                if (.not.streql(old_imas_version,imas_version).or.
     &              database.eq.'iter'.or.
     &              index(ids_path,'imasdb/iter').gt.0) then
                  if (.not.streql(old_imas_version,imas_version)) then
                    write(*,*)
     &               'Old IMAS data entry was written '//
     &               'using Data Dictionary version '//
     &                trim(old_imas_version)//'.'
                    write(*,*)
     &               'Recreating using Data Dictionary version '//
     &                trim(imas_version)//'.'
                  end if
                  if (database.eq.'iter'.or.
     &                index(ids_path,'imasdb/iter').gt.0) write(*,*)
     &             'IDS file will be moved to ITER database.'
                  call close_ual(idx)
cxpb Copy the IDS to a temporary location with the new DD and then bring it back
                  tmp_run = run
                  if (database.ne.'iter'.and.
     &                index(ids_path,'imasdb/iter').eq.0)
     &             tmp_run = run + 1000
#if AL_MAJOR_VERSION > 4
                  call al_build_uri_from_legacy_parameters
     &             ( MDSPLUS_BACKEND, shot, run, trim(username),
     &               trim(database), int2str(IMAS_MAJOR_VERSION), '',
     &               uri_source, status )
                  if (database.eq.'iter'.or.
     &                index(ids_path,'imasdb/iter').gt.0) then
                    call al_build_uri_from_legacy_parameters
     &               ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),
     &                 'ITER', int2str(IMAS_MAJOR_VERSION), '',
     &                 uri_dest, status )
                  else
                    call al_build_uri_from_legacy_parameters
     &               ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),
     &                 trim(database), int2str(IMAS_MAJOR_VERSION), '',
     &                 uri_dest, status )
                  end if
#if IMAS_MAJOR_VERSION < 4
                  write(systemarg,'(a,a,a,a,a)')
     &              'idscp --set-dataset-version'//
     &                   ' -s ',trim(uri_source),
     &                   ' -d ',trim(uri_dest),' --dd-update'
#else
                  write(systemarg,'(a,a,a,a,a)')
     &              'idscp -s ',trim(uri_source),
     &                   ' -d ',trim(uri_dest),' --dd-update'
#endif
#else
#if ( IMAS_MINOR_VERSION > 31 || IMAS_MAJOR_VERSION > 3 )
                  write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)')
     &              'idscp --setDatasetVersion'//
     &                   ' -si ',shot,' -ri ',run,
     &                   ' -so ',shot,' -ro ',tmp_run,
     &                   ' -d ',trim(database),' -u ',trim(username)
#else
                  write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)')
     &              'idscp -si ',shot,' -ri ',run,
     &                   ' -so ',shot,' -ro ',tmp_run,
     &                   ' -d ',trim(database),' -u ',trim(username)
#endif
                  if (database.eq.'iter'.or.
     &                index(ids_path,'imasdb/iter').gt.0)
     &             systemarg = trim(systemarg)//' -do ITER'
#endif
#ifdef NAGFOR
                  call system(systemarg, status, ierror)
#else
                  call system(systemarg)
#endif
                  if (database.ne.'iter'.and.
     &                index(ids_path,'imasdb/iter').eq.0) then
#if AL_MAJOR_VERSION > 4
                    call al_build_uri_from_legacy_parameters
     &               ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),
     &                 trim(database), int2str(IMAS_MAJOR_VERSION), '',
     &                 uri_source, status )
                    call al_build_uri_from_legacy_parameters
     &               ( MDSPLUS_BACKEND, shot, run, trim(username),
     &                 trim(database), int2str(IMAS_MAJOR_VERSION), '',
     &                 uri_dest, status )
                    write(systemarg,'(a,a,a,a)')
     &               'idscp -s ',trim(uri_source),
     &                    ' -d ',trim(uri_dest)
#else
#if ( IMAS_MINOR_VERSION > 31 || IMAS_MAJOR_VERSION > 3 )
                    write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)')
     &               'idscp --setDatasetVersion'//
     &                  ' -si ',shot,' -ri ',tmp_run,
     &                  ' -so ',shot,' -ro ',run,
     &                  ' -d ',trim(database),' -u ',trim(username)
#else
                    write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)')
     &               'idscp -si ',shot,' -ri ',tmp_run,
     &                  ' -so ',shot,' -ro ',run,
     &                  ' -d ',trim(database),' -u ',trim(username)
#endif
#endif
#ifdef NAGFOR
                    call system(systemarg, status, ierror)
#else
                    call system(systemarg)
#endif
                  end if
                  if (database.eq.'iter') database = 'ITER'
                  if (index(ids_path,'imasdb/iter').gt.0) then
                    l=index(ids_path,'imasdb/iter')
                    write(ids_path(l+7:l+10),'(a4)') 'ITER'
                  end if
#if AL_MAJOR_VERSION > 4
                  uri = 'imas:mdsplus?path='//trim(ids_path)
                  call imas_open( uri, OPEN_PULSE,
     .             idx, status, message )
#else
                  call imas_open_env(treename, shot, run, idx,
     .             username, database, version, status)
#endif
                end if
                if (continued) then
                  write(*,*)
     .             "Appending a new time slice at t = ", tim, " s."
                  num_time_slices = num_time_slices + 1
                end if
                time_slice_index = num_time_slices
                call B25_process_ids(
     &            geo, mpg, state, state_ext, state_avg, switch,
     &            edge_profiles, edge_sources, edge_transport,
     &            radiation, description, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &            summary,
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
     &            numerics,
#endif
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
     &            old_start_time, run_end_time,
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
     &            divertors,
#endif
     &            tim, dtim, shot, run, database, version, new_eq_ggd,
     &            time_slice_index, num_time_slices )
              else
                write(*,'(a)')
     .           'Not a time continuation, edge IDS components '//
     .           'will be replaced !'
                ids_replace = .true.
              end if
            end if
          else
            write(*,*)
     .       "No previous IMAS data entry found, "//
     .       "a new one will be created"
            idx = 0
            if (database.eq.'iter') database = 'ITER'
            if (index(ids_path,'imasdb/iter').gt.0) then
              l=index(ids_path,'imasdb/iter')
              write(ids_path(l+7:l+10),'(a4)') 'ITER'
            end if
          end if
        else
          write(*,*) "Creating a new IMAS data entry"
          idx = 0
          if (database.eq.'iter') database = 'ITER'
          if (index(ids_path,'imasdb/iter').gt.0) then
            l=index(ids_path,'imasdb/iter')
            write(ids_path(l+7:l+10),'(a4)') 'ITER'
          end if
        end if
        if (status.ne.0 .or. idx.eq.0 .or. ids_replace) then
          if (status.ne.0) idx = 0
          call B25_process_ids(
     &      geo, mpg, state, state_ext, state_avg, switch,
     &      edge_profiles, edge_sources, edge_transport,
     &      radiation, description, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &      summary,
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
     &      numerics,
#endif
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
     &      run_start_time, run_end_time,
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
     &      divertors,
#endif
     &      tim, dtim, shot, run, database, version, new_eq_ggd )
        end if
        if (ids_replace) then
          call new_ids_edge( edge_profiles, edge_sources,
     &      edge_transport, radiation, description, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &      summary,
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
     &      numerics,
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
     &      divertors,
#endif
#if IMAS_MAJOR_VERSION > 3
     &      uri,
#endif
     &      idx, new_eq_ggd )
        else
          call put_ids_edge( edge_profiles, edge_sources,
     &      edge_transport, radiation, description, equilibrium,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &      summary,
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
     &      numerics,
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
     &      divertors,
#endif
#if AL_MAJOR_VERSION > 4
     &      ids_path,
#else
     &      treename, shot, run, username, database, version,
#endif
     &      idx, new_eq_ggd )
        end if
        call dealloc_ids_edge(
     &      edge_profiles, edge_sources, edge_transport,
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
     &      numerics,
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
     &      divertors,
#endif
     &      radiation )
        call dealloc_batch_edge( batch_profiles, batch_sources,
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
     &      summary,
#endif
     &      description )
        call close_ual(idx)
#endif
      end if

#ifdef CATALYST
      if (use_catalyst .eqv. .true.) call coprocessorfinalize()
#endif

* finalise eirene_mc
      if(switch%use_eirene.ne.0.and.itim.gt.0) then
#ifdef B25_EIRENE
         call eirene_mc_final(nCv,nFc,ns,switch,mpg,geo,state,state_ext)
#else
         call xerrab ('b2mndr not compiled with the B25_EIRENE option')
#endif
      endif

      if (iav_run > 0) then
        call run_av_save ('b2favere',nCv,ns,.true.,.true.,state_avg)
        call run_av_fin (state_avg)
      end if

      if (ibatch_av_all > 0) call batch_av_all_fin

*  .. finalize tracing routines
      call b2trcf

* ..return
      ncall = ncall+1
      if ((.not.flag_optim).and.(switch%b2optim_namelist.eq.1)) then
        deallocate(par_opt_phys)
        deallocate(xsave)
      endif
      deallocate(old_erosion, old_deposition)
      deallocate(rza0, rz20, rpt0, rpi0)
      call dealloc_b2mod_diag
!WG_RM      call dealloc_b2mod_ma28
      call dealloc_b2mod_wall
      call dealloc_b2mod_layer
      call dealloc_b2mod_b2cmrc
      call dealloc_sputter_data
      call dealloc_b2mod_sputter
      call dealloc_b2mod_tallies
#ifndef FIXED_POINT
      call dealloc_b2mod_balance
#endif
      call dealloc_input_profile
      call dealloc_b2mod_elements
      call dealloc_b2mod_boundary
      call dealloc_b2mod_neutrals
      call dealloc_b2mod_numerics
      call dealloc_b2mod_user
      call dealloc_b2mod_par_opt
      call dealloc_feedback
!WG_RM      call dealloc_b2mod_external
      call dealloc_b2mod_wall_init
      call dealloc_b2mod_ysmp_sdrv
      call dealloc_b2mod_facdrift_ExB
      call dealloc_transport_namelist
      call dealloc_b2mod_first_flight
      call dealloc_b2mod_neoclassical
!WG_RM      call dealloc_b2mod_ma28_for_7diag
!WG_RM      call dealloc_b2mod_ma28_for_9diag
      call dealloc_b2mod_ma28_for_us
#ifndef FIXED_POINT
      call dealloc_b2mod_eirene_sources
#endif
!WG_TODO      call dealloc_b2mod_external_dummy
      call dealloc_b2mod_transport_models
      call dealloc_b2mod_transport_disruption
      if(switch%boundary_sources.ne.0)
     &  call dealloc_b2mod_boundary_sources     !srv 26.02.18
      if(switch%equation_sources.ne.0)
     &  call dealloc_b2mod_equation_sources     !srv 26.02.18
      call dealloc_b2mod_b2plot
#ifdef B25_EIRENE
      CALL EIRENE_DEALLOC_BRAEIR
      CALL EIRENE_DEALLOC_EIRBRA
      CALL EIRENE_DEALLOC_CCOUPL
      CALL EIRENE_DEALLOC_CTRCEI
      CALL EIRENE_DEALLOC_INFCOP
#else
      call dealloc_b2mod_braeir
      call destroyB2SourceEir( state % sr_eir )
#endif
      call subend ()
      return
      end subroutine b2mndr_2

*-----------------------------------------------------------------------
*.end b2mndr
c
      subroutine dpc_mod_rates(te_hot, ne_hot_frac)
      use b2mod_types
      use b2mod_b2cmrc

      implicit none

      real (R8) :: te_hot
      real (R8) :: ne_hot_frac
      real (R8) :: rlte, fxte
      real (R8) :: rate_sa, rate_ra, rate_qa, rate_rd, rate_br,
     *             rate_za, rate_z2, rate_pt, rate_pi, rate_qr
!!    real (R8) :: rate_cx
      integer ite, ine, ise, ixte

      write(*,*) 'dpc_mod_rates: ', te_hot, ne_hot_frac

      ite = 0
      rlte = max(rtlt(0),min(rtlt(rtnt),log(te_hot)))

    1 if (ite.lt.rtnt-1) then
        if (rtlt(ite+1).le.rlte) then
          ite = ite+1
          goto 1
        endif
      endif
    2 if (0.lt.ite) then
        if (rlte.lt.rtlt(ite)) then
          ite = ite-1
          goto 2
        endif
      endif
*     ..compute coefficient for linear interpolation
      ixte = ite
      fxte = (rlte-rtlt(ite))/(rtlt(ite+1)-rtlt(ite))
      write(*,*) ixte
      write(*,*) rtn(rtnn/2)
      write(*,*) rtt(:)
      write(*,*) exp(rtlsa(:,rtnn/2,0))
      do ise= 0, rtns-1
        do ine= 0, rtnn
          rate_sa = exp((1.0_R8-fxte)*rtlsa(ixte,ine,ise) +
     1     fxte*rtlsa(ixte+1,ine,ise)) * ne_hot_frac
          rate_ra = exp((1.0_R8-fxte)*rtlra(ixte,ine,ise) +
     1     fxte*rtlra(ixte+1,ine,ise)) * ne_hot_frac
          rate_qa = exp((1.0_R8-fxte)*rtlqa(ixte,ine,ise) +
     1     fxte*rtlqa(ixte+1,ine,ise)) * ne_hot_frac
          rate_qr = exp((1.0_R8-fxte)*rtlqr(ixte,ine,ise) +
     1     fxte*rtlqr(ixte+1,ine,ise)) * ne_hot_frac
          rate_rd = exp((1.0_R8-fxte)*rtlrd(ixte,ine,ise) +
     1     fxte*rtlrd(ixte+1,ine,ise)) * ne_hot_frac
          rate_br = exp((1.0_R8-fxte)*rtlbr(ixte,ine,ise) +
     1     fxte*rtlbr(ixte+1,ine,ise)) * ne_hot_frac
!!        rate_cx = exp((1.0_R8-fxte)*rtlcx(ixte,ine,ise) +
!!   1     fxte*rtlcx(ixte+1,ine,ise)) * ne_hot_frac
          rate_za = ((1.0_R8-fxte)*rtlza(ixte,ine,ise) +
     1     fxte*rtlza(ixte+1,ine,ise)) * ne_hot_frac
          rate_z2 = ((1.0_R8-fxte)*rtlz2(ixte,ine,ise) +
     1     fxte*rtlz2(ixte+1,ine,ise)) * ne_hot_frac
          rate_pt = ((1.0_R8-fxte)*rtlpt(ixte,ine,ise) +
     1     fxte*rtlpt(ixte+1,ine,ise)) * ne_hot_frac
          rate_pi = ((1.0_R8-fxte)*rtlpi(ixte,ine,ise) +
     1     fxte*rtlpi(ixte+1,ine,ise)) * ne_hot_frac
          do ite= 0, rtnt
            rtlsa(ite,ine,ise) =
     1       log(exp(rtlsa(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_sa)
            rtlra(ite,ine,ise) =
     1       log(exp(rtlra(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_ra)
            rtlqa(ite,ine,ise) =
     1       log(exp(rtlqa(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_qa)
            rtlqr(ite,ine,ise) =
     1       log(exp(rtlqr(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_qr)
            rtlrd(ite,ine,ise) =
     1       log(exp(rtlrd(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_rd)
            rtlbr(ite,ine,ise) =
     1       log(exp(rtlbr(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_br)
!!          rtlcx(ite,ine,ise) =
!!   1       log(exp(rtlcx(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_cx)
            rtlza(ite,ine,ise) =
     1       rtlza(ite,ine,ise)*(1.0_R8-ne_hot_frac)+rate_za
            rtlz2(ite,ine,ise) =
     1       rtlz2(ite,ine,ise)*(1.0_R8-ne_hot_frac)+rate_z2
            rtlpt(ite,ine,ise) =
     1       rtlpt(ite,ine,ise)*(1.0_R8-ne_hot_frac)+rate_pt
            rtlpi(ite,ine,ise) =
     1       rtlpi(ite,ine,ise)*(1.0_R8-ne_hot_frac)+rate_pi
          enddo
        enddo
      enddo
      write(*,*) exp(rtlsa(:,rtnn/2,0))
      return
      end subroutine dpc_mod_rates

      end module b2mod_driver
