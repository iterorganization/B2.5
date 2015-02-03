*-----------------------------------------------------------------------
*.specification

      subroutine b2mndr (ninp, nout, nx, ny, ns, ns0)
      use b2mod_types
      use b2mod_version
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_work
      use b2mod_tallies
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_time
      use b2mod_diag
     * , only : ank_tracing
      use b2mod_elements
      use b2mod_transport_models
      use b2mod_neutrals_namelist
      use b2mod_neutr_src_scaling
      use b2mod_b2cmfs
      use b2mod_b2cmgs
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      use b2mod_b2cmrc
      use b2mod_boundary_namelist
#ifdef ASTRA
      use b2mod_b2_to_astra
#endif	                       
#ifdef JET
      use b2mod_ppout
#endif
      implicit none
*   ..input arguments (unchanged on exit)
#ifdef JET
      integer ninp(0:5), nout(0:8), nx, ny, ns, ns0
#else
      integer ninp(0:5), nout(0:8), nx, ny, ns, ns0
#endif
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
#     include "DIMENSIONS.F"
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MNDR is the main driver routine for the b2 code.
*
*
*  3. description (see also routine b2cdca)
*
*     For a detailed description of the equations and the method of
*     solution, see the references cited in the main program unit.
*
*     The routine b2mndr controls the timestepping procedure in a b2
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
*     ninp - (0:5) integer array, input.
*     ninp specifies the input unit numbers, as follows.
*     ninp(0): formatted; provides physics and numerics parameters.
*     ninp(1): un*formatted; provides the geometry.
*     ninp(2): un*formatted; provides default physics parameters.
*     ninp(3): un*formatted; provides the initial state.
*     ninp(4): un*formatted; provides atomic rate data.
*     ninp(5): un*formatted; provides source term data.
*
*     nout - (0:5) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides summary print output.
*     nout(1): un*formatted; provides actual physics parameters.
*     nout(2): un*formatted; provides the final state.
*     nout(3): un*formatted; may provide data for movie.
*     nout(4): un*formatted; provides data on the numerical evolution.
*     nout(5): un*formatted; may provide detailed convergence data.
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

*   ..common blocks
*   ..local variables
      integer ncall, ntim, mvinc, mvnum, ix, iy, is, ib,js,jsep,jxi,jxa
      integer, save :: save_f30 = 0, save_f31 = 0, background_new=0  !srv 01.03.10 14.07.10
      integer, save :: inverse_ua = 0                                !srv 11.01.13
      integer, save :: pot_eq=1
      real (kind=R8), save :: fac_hz = 0.0e0_R8                      !srv 09.01.01 09.06.08
      real (kind=R8) ::
     *  na0eps, dtim, stim, etim, b2mndr_cpu, cpuval, cpustart,
     &  density_rescale, BoRiS
      character idout0*32, idout1*32, savefile*16, savefilestatus*3
      integer n2, itim, ierr, edition, i, b2time, coronal, use_eirene,
     *  tally, ismain, ismain0, plate_model, natmi, sput_dst, no_solve
      real (kind=R8) ::
     *  cputarget, cpuincrement, neutral_rescale, 
     *  plate_temp, plate_thick, neutral_sources_rescale,
     1  na_eps, te_eps, ti_eps, po_eps, ua_eps, sput_frc, sput_phys,
     2  delta_min, dt_change_inc, delta_max, dt_change_dec,
     3  dt_min, dt_max, max_delta
#ifdef JET
      real (kind=R8) ::
     *  delta_tran_time, save_tran_time
      integer tran_num
      character*6 tran_file_name
#endif
      real (kind=R8) ::
     *  delta_cdfmovie_time, save_cdfmovie_time
      real (kind=R8) ::
     *  delta_plasma_time, save_plasma_time
      logical quitexist, quit, quitexist_
	logical, save :: outputexist = .false.
      real (kind=R8) ::
     *  stack_dtim(20)
      integer stack_ntim(20)
      integer stack_ptr
      logical ok, found
      real (kind=R8) ::                                                  !jwk
     *  min_areshe, min_areshi, min_aresco                               !jwk
      logical quit_residual, test_residual                               !jwk
      real (kind=R8) :: amtol
      parameter (amtol=1.0e-2_R8)
      real (kind=R8) :: 
     *  rsa(-1:nx,-1:ny,0:ns-1), rra(-1:nx,-1:ny,0:ns-1),
     &  rqa(-1:nx,-1:ny,0:ns-1), rrd(-1:nx,-1:ny,0:ns-1),
     &  rbr(-1:nx,-1:ny,0:ns-1),
     *  rza0(-1:nx,-1:ny,0:ns-1), rz20(-1:nx,-1:ny,0:ns-1),
     *  rpt0(-1:nx,-1:ny,0:ns-1), rpi0(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) :: hzconst                                      !srv 30.01.10 
csw 18feb2014 added for correct implementation of boundary conditions
      integer, save :: l_fbound=0, l_enebound=0, l_timbound=0
*   ..procedures
      logical lnext
      real (kind=R8) ::
     *  x05baf
      external x05baf
      intrinsic mod, abs, sqrt
      external subini, subend, xertst, xerrab, ipgeti, ipgetr, 
     &  cfrure, cfwure, scopy
      external b2rugm, b2rups, b2wfgm, b2wfpi, b2wups, b2mnds, b2mndt,
     &  b2mwmv, b2mxnp, b2mxnu, b2xgbs, intfaceh1, intfacer1, intfacev1, !xpb
     &  my_intcor, b2stbc_fb_final, b2spel, b2sqel, b2xvps, b2xvsg
      lnext(is,js) = zamax(is).lt.zamin(js).and.zn(is).eq.zn(js).and.
     &  abs(am(is)-am(js)).lt.amtol
*   ..initialisation
      save ncall, ntim, mvinc, mvnum, dtim, stim, etim, idout0, idout1,
     1  na0eps, b2mndr_cpu, edition, cputarget, cpuincrement, ismain,
     2  ismain0, savefilestatus, b2time, coronal, use_eirene, natmi, 
     3  sput_dst, neutral_rescale, neutral_sources_rescale, tally, 
     4  BoRiS, density_rescale, na_eps, te_eps, ti_eps, po_eps, ua_eps, 
     5  sput_frc, sput_phys, delta_min, dt_change_inc, delta_max, 
     6  dt_change_dec, no_solve, dt_min, dt_max, plate_temp, 
     7  plate_thick, plate_model
#ifdef JET
      save delta_tran_time, save_tran_time, tran_num
#endif
      save delta_cdfmovie_time, save_cdfmovie_time
      save delta_plasma_time, save_plasma_time
      save min_areshe, min_areshi, min_aresco                            !jwk
      data ncall/0/, ntim/1/, mvinc/1/, mvnum/0/, dtim/1.0_R8/,
     1     stim/0.0_R8/, etim/0.0_R8/, ismain/1/, ismain0/0/, 
     2     plate_model/0/
      data idout0/'pgnl;pgmm;pzmm'/, idout1/'pzmm'/
      data na0eps/1.0e-30_R8/, b2mndr_cpu/0.0_R8/
      data savefile/'plasmastate.xxxx'/,
     1     cpuincrement/3600.0_R8/,savefilestatus/'new'/,b2time/1/,
     2     coronal/0/, use_eirene/0/, neutral_rescale/1.0_R8/,
     *     neutral_sources_rescale/1.0_R8/,tally/1/, sput_dst/-1/,
     1     density_rescale/1.0_R8/, sput_frc/0.0_R8/, sput_phys/0.0_R8/,
     *     plate_temp/300.0_R8/, plate_thick/0.00_R8/, no_solve/0/,
     2     na_eps/1.0e19_R8/,
     3     te_eps/1.0e+1_R8/,
     4     ti_eps/1.0e+1_R8/,
     5     po_eps/1.0e+1_R8/,
     6     ua_eps/1.0e+4_R8/,
     7     delta_min/0.0_R8/, 
     8     dt_change_inc/1.0_R8/, 
     9     delta_max/0.0_R8/, 
     1     dt_change_dec/1.0_R8/,
     2     dt_min/1.0e-30_R8/,
     3     dt_max/1.0e+01_R8/,
     4     BoRiS/0.0_R8/
#ifdef JET
      data delta_tran_time/0.0_R8/, save_tran_time/0.0_R8/,
     1     tran_num/0/
#endif
      data delta_cdfmovie_time/0.0_R8/, save_cdfmovie_time/0.0_R8/
      data delta_plasma_time/0.0_R8/, save_plasma_time/0.0_R8/
      data min_areshe/0.0_R8/, min_areshi/0.0_R8/, min_aresco/0.0_R8/    !jwk
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
*     /b2cmgs/ specifies the singular regions in the geometry.
*      nlrmx - integer constant.
*     nlrmx is an upper bound for nlreg, and is used in dimension
*     specifications.
*      nlreg - integer.
*     nlreg specifies the number of singular regions. It will hold that
*     0.le.nlreg.le.nlrmx.
*      nlxlo, nlxhi, nlylo, nlyhi - (0:nlrmx-1) integer array.
*     For 0.le.i.lt.nlreg, nlxlo(i), nlxhi(i), nlylo(i) and nlyhi(i) are
*     the bounds of the i-th singular region. Mesh locations (ix,iy)
*     in (nlxlo(i):nlxhi(i),nlylo(i):nlyhi(i)) make up the region, which
*     will be included in (-1:nx,-1:ny).
*     For nlreg.le.i.lt.nlrmx, nlxlo(i), nlxhi(i), nlylo(i) and nlyhi(i)
*     are unspecified.
*      nlloc - (0:nlrmx-1) integer array.
*     For 0.le.i.lt.4, nlloc(i) specifies the class of quantities for
*     which the i-th region is special, as follows:
*       nlloc(i) = 0: cell centers, indices (-1:nx,-1:ny).
*       nlloc(i) = 1: cell x-faces, indices (0:nx,-1:ny).
*       nlloc(i) = 2: cell y-faces, indices (-1:nx,0:ny).
*       nlloc(i) = 3: cell corners, indices (0:nx,0:ny).
*     For nlreg.le.i.lt.nlrmx, nlloc(i) is unspecified.
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
*     /b2cmpb/ specifies physics parameters for boundary conditions.
*     For this purpose the boundary is divided up into regions, and for
*     each region a set of coefficients defines the local conditions.
*     For detailed description of /b2cmpb/ see the code, especially the
*     routine b2stbc.
*      cbnrmx - integer constant.
*     cbnrmx is an upper bound for the number of boundary regions nreg,
*     (nreg.eq.cbnrso+cbnrno+cbnrwe+cbnrea) and is used in dimension
*     specifications.
*      cbirso, cbirno, cbirwe, cbirea - integer.
*     cbirso(-no,-we,-ea) is the index of the first region on the
*     south (north, west, east) boundary.
*      cbnrso, cbnrno, cbnrwe, cbnrea - integer.
*     cbnrso(-no,-we,-ea) is the number of regions on the south
*     (north, west, east) boundary.
*     It will hold that 1.le.cbnrso, 1.le.cbnrno, 1.le.cbnrwe, and
*     1.le.cbnrea; and 0.eq.cbirso, cbirso+cbnrso.eq.cbirno,
*     cbirno+cbnrno.eq.cbirwe, cbirwe+cbnrwe.eq.cbirea.
*      cbrbrk - (0:cbnrmx-1) real*8 array.
*     cbrbrk defines the boundaries between the boundary regions.
*      cbsna - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.is.lt.ns and 0.le.ireg.lt.nreg, the eight coefficients
*     cbsna(0:7,is,ireg) specify standard form boundary conditions for
*     the continuity equation for species (is) on the (ireg) boundary
*     region.
*      cbsmo - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.is.lt.ns and 0.le.ireg.lt.nreg, the eight coefficients
*     cbsmo(0:7,is,ireg) specify standard form boundary conditions for
*     the momentum balance equation for species (is) on the (ireg)
*     boundary region.
*      cbshi - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.is.lt.ns and 0.le.ireg.lt.nreg, the eight coefficients
*     cbshi(0:7,is,ireg) specify standard form boundary conditions for
*     the contribution to that atom heat equation by species (is) on
*     the (ireg) boundary region.
*      cbshe - (0:7,0:cbnrmx-1) real*8 array.
*     For 0.le.ireg.lt.nreg, the eight coefficients cbshe(0:7,ireg)
*     specify standard form boundary conditions for the electron heat
*     equation on the (ireg) boundary region.
*      cbsch - (0:7,0:cbnrmx-1) real*8 array.
*     For 0.le.ireg.lt.nreg, the eight coefficients cbsch(0:7,ireg)
*     specify standard form boundary conditions for the electric
*     potential equation on the (ireg) boundary region.
*      cbrec - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.is.lt.ns-1 and 0.le.ireg.lt.nreg, the eight coefficients
*     cbrec(0:7,is,ireg) specify the recycling process for species (is)
*     on the boundary region ireg. The exact interpretation of these
*     coefficients is defined by the routine that computes the source
*     terms due to the recycling process (b2stbr). Routine b2stbr is
*     problem dependent.
*      cbmsa - (0:7,0:cbnrmx-1) real*8 array.
*      cbmsc - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.ireg.lt.nreg, the eight coefficients cbmsa(0:7,ireg)
*     specify miscellaneous processes on the boundary region ireg.
*     For 0.le.is.lt.ns-1 and 0.le.ireg.lt.nreg, the eight coefficients
*     cbmsc(0:7,is,ireg) specify miscellaneous processes for species
*     (is) on the boundary region ireg. The interpretation of these
*     coefficients is defined by the routine that computes the source
*     terms due to miscellaneous processes (b2stbm). Routine b2stbm is
*     application dependent.
*
*   ..transport coefficients, parameters
*
*     /b2cmpt/ - common block.
*     (Note: the integer constant nsdecl, which comes with /b2cmpa/,
*     is used instead of (ns) in dimension specifications.)
*      cfdf0 - (0:7,0:ns-1) real*8 array.
*     cfdf0 contains some auxiliary coefficients.  For 0.le.is.lt.ns,
*     the local expression
*       df0(,,is) = cfdf0(0,is)+cfdf0(1,is)/(ni(,,0)/1.0e20)+
*         cfdf0(2,is)*sqrt(ti(,)/(am(is)*mp))/
*          (cfdf0(3,is)*ni(,,0)+cfdf0(4,is)*ne(,))
*     is used in the definition of the anomalous transport coefficients.
*     df0 is meant to represent a diffusivity.  If cfdf0(2,is).eq.0 then
*     the corresponding term in df0 is understood as 0 independent of
*     cfdf0(3:4,is).
*      cfdna - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfdna(,is) are the coefficients in the
*     local expression
*       difna(,,is) = cfdna(0,is)+cfdna(1,is)/(ne(,)/1.0e20)+
*         cfdna(2,is)*df0(,)+
*         cfdna(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous diffusivity with respect to
*     grad.na; the transport equation for species (is) contains a term
*     difna(,,is)*'grad.na'.
*     If cfdna(7,is) is not zero, the dna0 array is flux-scaled with
*     a scaling factor of cfdna(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfdna(4:6,) are not used in the default version of the code.
*      cfdpa - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfdpa(,is) are the coefficients in the
*     local expression
*       difpa(,,is) = cfdpa(0,is)+cfdpa(1,is)/(ne(,)/1.0e20)+
*         cfdpa(2,is)*df0(,,is)+
*         cfdpa(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous diffusivity with respect to
*     grad.pa; the transport equation for species (is) contains a term
*     (1/(rza(,,is)*te(,)+ti(,)))*difpa(,,is)*'grad.pa'.
*     If cfdpa(7,is) is not zero, the dpa0 array is flux-scaled with
*     a scaling factor of cfdpa(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfdpa(4:6,) are not used in the default version of the code.
*      cfvla - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfvla(,is) are the coefficients in the
*     local expression
*       vlay(,,is) = cfvla(0,is)+cfvla(1,is)/(ne(,)/1.0e20)+
*         cfvla(2,is)*df0(,,is)+
*         cfvla(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous y-velocity of species (is); the
*     transport equation for species (is) contains a term
*     vlay(,,is)*na(,,is).
*     If cfvla(7,is) is not zero, the vla0 array is flux-scaled with
*     a scaling factor of cfvla(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfvla(4:6,) are not used in the default version of the code.
*      cfvsa - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfvsa(,is) are the coefficients in the
*     local expression
*       dfvsa(,,is) = cfvsa(0,is)+cfvsa(1,is)/(ne(,)/1.0e20)+
*         cfvsa(2,is)*df0(,,is)+
*         cfvsa(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the diffusivity corresponding to the anomalous
*     transverse viscosity, so that this viscosity is
*       vsa(,,is) = mp*am(is)*na(,,is)*dfvsa
*     If cfvsa(7,is) is not zero, the vsa0 array is flux-scaled with
*     a scaling factor of cfvsa(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfvsa(4:6,) are not used in the default version of the code.
*      cfhci - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfhci(,is) are the coefficients in the
*     local expression
*       dfhca(,,is) = cfhci(0,is)+cfhci(1,is)/(ne(,)/1.0e20)+
*         cfhci(2,is)*df0(,,is)+
*         cfhci(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous atom heat conductivity associated
*     with species (is), so that the atom heat conduction coefficient
*     contains a term dfhca(,,is)*na(,,is).
*     If cfhci(7,is) is not zero, the hci arrays are flux-scaled with
*     a scaling factor of cfhci(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfhci(4:6,) are not used in the default version of the code.
*      cfhce - (0:7) real*8 array.
*     cfhce() are the coefficients in the local expression
*       dfhce(,) = cfhce(0)+cfhce(1)/(ne(,)/1.0e20)+
*         cfhce(3)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous electron heat conductivity, so that
*     the electron heat conduction coefficient is dfhce(,)*ne(,).
*     If cfhce(7) is not zero, the hce0 array is flux-scaled with
*     a scaling factor of cfhce(7)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfhce(2) and (4:6) are not used in the default version of the code.
*      cfsig - (0:7) real*8 array.
*     cfsig() are the coefficients in the local expression
*       dfsig(,) = cfsig(0)+cfsig(1)/(ne(,)/1.0e20)+
*         cfsig(3)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the diffusivity corresponding to the anomalous
*     electrical conductivity, which is
*       sig(,) = qe*ne(,)*dfsig(,).
*     If cfsig(7) is not zero, the sig0 array is flux-scaled with
*     a scaling factor of cfsig(7)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfsig(2) and (4:6) are not used in the default version of the code.
*      cfalf - (0:7) real*8 array.
*     cfalf() are the coefficients in the local expression
*       dfalf(,) = cfalf(0)+cfalf(1)/(ne(,)/1.0e20)+
*         cfalf(3)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the diffusivity corresponding to the anomalous
*     thermo-electric coefficient, which is
*       alf(,) = ne(,)*sqrt(qe/te(,))*dfalf(,).
*     If cfalf(7) is not zero, the alf0 array is flux-scaled with
*     a scaling factor of cfalf(7)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfalf(2) and (4:6) are not used in the default version of the code.
*      cflim - (0:7) real*8 array.
*     cflim contains coefficients for flux limiting. cflim(0)=cflme
*     specifies the flux limit factor for poloidal electron heat
*     transport. cflim(1)=cflmi specifies the flux limit factor for
*     poloidal all atom heat transport. cflim(2)=cflmv specifies the
*     flux limit factor for poloidal atom momentum transport.
*     Zero settings turn off these flux limits.
*     cflim(3:7) are not used in the default version of the code.
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
*     and izhi(inuc) is the highest atomic charge considered in b2 for
*     that nuclear species.
*     rtns equals the variable ns used elsewhere in b2.
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
*     zamin, zamax and zn used elsewhere in b2.
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
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlqa(i,j,k) specifies
*     logu(rtrqa), where rtrqa is the electron heat loss rate
*     coefficient at electron temperature rtt(i) and electron density
*     rtn(j), for processes starting from atomic species k. rtrqa has
*     dimension ev*m**3/s.
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
*     crx, cry - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), crx(ix,iy,0:3) specify the
*     horizontal coordinates and cry(ix,iy,0:3) specify the vertical
*     coordinates (in the physical plane) of the four corners of the
*     (ix,iy) cell, mapped as follows:
*     (,,0): lower left corner,  (,,1): lower right corner,
*     (,,2): upper left corner,  (,,3): upper right corner.
*     One expects the following identities to hold away from the
*     boundaries:
*     crx(ix,iy,0)=crx(ix-1,iy,1)=crx(ix,iy-1,2)=crx(ix-1,iy-1,3),
*     cry(ix,iy,0)=cry(ix-1,iy,1)=cry(ix,iy-1,2)=cry(ix-1,iy-1,3),
*     However, there may be cuts in the geometry, in which case these
*     identities need not hold.
*
*     bb - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), bb(ix,iy,0:3) specifies the
*     magnetic field at the center of the (ix,iy) cell.
*     bb(ix,iy,0:2) are the (x,y,z)-components and bb(ix,iy,3) is
*     the absolute magnetic field strength.
*     (bb(,,3)=sqrt(bb(,,0)**2+bb(,,1)**2+bb(,,2)**2.)
*     It will hold that 0.lt.bb(,,3).
*
*     vol - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), vol(ix,iy) specifies the volume of
*     the (ix,iy) cell. (In the case of toroidal symmetry this volume
*     includes a factor 2*pi*crx, while in the case of cylindrical
*     symmetry it includes a factor 2*pi*cry.)
*     It will hold that 0.lt.vol(,).
*
*     hx, hy - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), hx(ix,iy) and hy(ix,iy) specify the
*     x-diameter and the y-diameter, respectively, of the (ix,iy) cell.
*     It will hold that 0.lt.hx(,) and 0.lt.hy(,).
*
*     qz - (-1:nx,-1:ny,0:1) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), qz(ix,iy,0) specifies sin(t) and
*     qz(ix,iy,1) specifies cos(t), where t is the complementary angle
*     between the x-direction and the y-direction on the (ix,iy) cell.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. (qz(ix,iy,0) is the inner product of the x and y unit
*     vectors.)
*     It will hold that -1.lt.qz(,,0).lt.1 and 0.lt.qz(,,1).le.1.
*
*     qc - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qc(ix,iy) specifies cos(t), where 
*     t is the complementary angle between the x-direction of the 
*     (ix,iy) cell and its left face.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. 
*     It will hold that 0.lt.qc(,).le.1.
*
*     gs - (-1:nx,-1:ny,0:2) real*8 array.
*     (Let gsx(,)=gs(,,0), gsy(,)=gs(,,1) and gsz(,)=gs(,,2).)
*     For (ix,iy) in (0:nx,-1:ny), gsx(ix,iy) specifies the area of the
*     face between the (ix-1,iy) and the (ix,iy) cell. For ix.eq.-1,
*     gsx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), gsy(ix,iy) specifies the area of the
*     face between the (ix,iy-1) and the (ix,iy) cell. For iy.eq.-1,
*     gsy(-1:nx,iy) will hold 0.
*     For (ix,iy) in (-1:nx,-1:ny), gsz(ix,iy) specifies the area of the
*     (ix,iy) cell normal to the third (ignorable) coordinate.
*
*     na - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), na(ix,iy,is) specifies
*     the density of atomic species (is) on the (ix,iy) cell.
*     It will hold that 0.lt.na(,,).
*
*     ua - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), ua(ix,iy,is) specifies
*     the parallel velocity of atomic species (is) on the (ix,iy) cell.
*
*     uadia - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), uadia(ix,iy,0,is) specifies
*     the diamagnetic drift velocity of species (is) on the (ix,iy) cell, and
*     uadia(ix,iy,1,is) the radial drift velocity.
* 
*     te - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), te(ix,iy) specifies the electron
*     temperature on the (ix,iy) cell. It will hold that 0.lt.te(,).
*
*     ti - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), ti(ix,iy) specifies the all atom
*     temperature on the (ix,iy) cell. It will hold that 0.lt.ti(,).
*
*     po - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), po(ix,iy) specifies the electric
*     potential on the (ix,iy) cell.
*
*     fna - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let fnax(,,)=fna(,,0,) and fnay(,,)=fna(,,1,).)
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), fnax(ix,iy,is) specifies
*     the flux of atoms of species (is) through the face between the
*     (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fnax(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), fnay(ix,iy,is) specifies
*     the flux of atoms of species (is) through the face between the
*     (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fnay(-1:nx,iy,0:ns-1) will hold 0.
*
*     fhe - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhex(,)=fhe(,,0) and fhey(,)=fhe(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhex(ix,iy) specifies the electron
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For ix.eq.-1, fhex(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhey(ix,iy) specifies the electron
*     heat flux through the face between the (ix,iy-1) and (ix,iy)
*     cells. For iy.eq.-1, fhey(-1:nx,iy) will hold 0.
*
*     fhi - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhix(,)=fhi(,,0) and fhiy(,)=fhi(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhix(ix,iy) specifies the all atom
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For ix.eq.-1, fhix(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhiy(ix,iy) specifies the all atom
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For iy.eq.-1, fhiy(-1:nx,iy) will hold 0.
*
*     fch - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchx(,)=fch(,,0) and fchy(,)=fch(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchx(ix,iy) specifies the electric
*     current through the face between the (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fchx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchy(ix,iy) specifies the electric
*     current through the face between the (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fchy(-1:nx,iy) will hold 0.
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
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mndr')
*   ..test input data
*     (done in b2mnds)
*   ..initialise files and common
      call b2mnds (ninp, nout, nx, ny, ns, ns0)
*   ..set internal parameters on first (and only) call
*     (Note that the internal parameters are only defined following
*     the call to b2mnds.)
      if (ncall.eq.0) then
       call ipgeti ('b2mndr_inverse_ua', inverse_ua)                 !srv 11.01.13
       call ipgeti ('b2mndr_ntim', ntim)
       call ipgetr ('b2mndr_dtim', dtim)
       call xertst (0.le.ntim.and.0.0_R8.lt.dtim,
     &   'faulty internal parameter ntim, dtim')
       call ipgetr ('b2mndr_stim', stim)
       call ipgetr ('b2mndr_etim', etim)
       call xertst (stim.le.0.0_R8 .or. etim.gt.stim, 
     &   'faulty inputs stim,etim')
       call ipgetr ('b2news_BoRiS', BoRiS)
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
         if (is_neutral(is) .and. zn(is).eq.zn(ismain) .and.
     &       am(is).eq.am(ismain)) then
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
       call ipgeti ('b2mndr_mvinc', mvinc)
       call ipgeti ('b2mndr_mvnum', mvnum)
       call xertst (1.le.mvinc.and.0.le.mvnum,
     &   'faulty internal parameter mvinc, mvnum')
       call ipgetc ('b2mndr_idout0', idout0)
       call ipgetc ('b2mndr_idout1', idout1)
       call ipgetr ('b2mndr_na0eps', na0eps)
       call xertst (na0eps.gt.0.0_R8,
     &  'faulty internal parameter na0eps')
       call ipgetr ('b2mndr_cpu', b2mndr_cpu)
       call xertst (b2mndr_cpu.ge.0.0_R8,
     1  'faulty internal parameter cpu')
       call ipgetr ('b2mndr_savecpu', cpuincrement)
       call xertst (cpuincrement.ge.0.0_R8,
     1  'faulty internal parameter cpuincrement')
       cputarget=cpuincrement
#ifdef JET
       call ipgetr ('b2mndr_trantim', delta_tran_time)
       if (delta_tran_time.gt.0.0_R8) 
     >  write(*,*) 'Saving tran files every ',delta_tran_time,' seconds'
#endif
       call ipgetr ('b2mndr_cdfmovietim', delta_cdfmovie_time)
       if (delta_cdfmovie_time.gt.0.0_R8)
     >  write(*,*) 'Saving CDF movie files every ',
     1  delta_cdfmovie_time,' simulation seconds'
       call ipgetr ('b2mndr_plasmatim', delta_plasma_time)
       if (delta_plasma_time.gt.0.0_R8)
     >  write(*,*) 'Saving plasma state files every ',
     1  delta_plasma_time,' simulation seconds'
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
       call ipgetr ('b2mndr_po_eps',po_eps)
       call xertst (0.0_R8.lt.po_eps,'faulty internal parameter po_eps')
       call ipgetr ('b2mndr_ua_eps',ua_eps)
       call xertst (0.0_R8.lt.ua_eps,'faulty internal parameter ua_eps')
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
       call ipgeti ('b2mndr_coronal_model', coronal)
       if(coronal.ne.0) then
          write(*,*)
     1     'Modifying rate coefficients to suppress density dependence'
          do is=1,rtnn
             rtlsa(:,is,:)=rtlsa(:,0,:)
             rtlra(:,is,:)=rtlra(:,0,:)
             rtlqa(:,is,:)=rtlqa(:,0,:)
             rtlcx(:,is,:)=rtlcx(:,0,:)
             rtlrd(:,is,:)=rtlrd(:,0,:)
             rtlbr(:,is,:)=rtlbr(:,0,:)
             rtlza(:,is,:)=rtlza(:,0,:)
             rtlz2(:,is,:)=rtlz2(:,0,:)
             rtlpt(:,is,:)=rtlpt(:,0,:)
             rtlpi(:,is,:)=rtlpi(:,0,:)
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
       call ipgeti ('b2mndr_eirene',use_eirene)
       call ipgeti ('eirene_savef30', save_f30)                    !srv 01.03.10
       call ipgeti ('eirene_savef31', save_f31)                    !srv 01.03.10
       call ipgeti ('eirene_background_new', background_new)       !srv 14.07.10
       call ipgetr ('b2mndr_rescale_neutrals',neutral_rescale)
       call ipgetr ('b2mndr_rescale_neutrals_sources',
     *       neutral_sources_rescale)
c-dpc-2002.03.31       if(use_eirene.eq.0) neutral_rescale=1.0_R8
       if(use_eirene.eq.0) neutral_sources_rescale=1.0_R8
       if(use_eirene.ne.0) then
         write(*,*) 'Adding ',log(neutral_sources_rescale),
     *         ' to relevant rtlra'
         do is=0,ns-2
           if(is_neutral(is) .and. lnext(is,is+1)) then
              rtlra(:,:,is+1)=rtlra(:,:,is+1)+
     *             log(neutral_sources_rescale)
              write(*,*) is,zamax(is),rtlra(0,0,is+1)
           endif
         enddo
       endif
       call ipgetr ('b2mndr_hz', fac_hz)				               !srv 09.01.01 09.06.08 {
       call xertst (0.0.le.fac_hz .and. 1.0.ge.fac_hz,
     &  'faulty parameter b2mndr_hz')					               !srv 09.01.01 09.06.08 }
csw 09oct2012
       call ipgeti ('b2news_poteq', pot_eq) 
csw 21feb2014
       call ipgeti ('b2news_l_fbound',l_fbound)
       call ipgeti ('b2news_l_enebound',l_enebound)
       call ipgeti ('b2news_l_timbound',l_timbound)
* ..return
      else
       call xerrab ('error, b2mndr has been called more than once')
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)

* ..initialisation
*   ..read geometry and metric
      call b2rugm (ninp(1), nx, ny, crx, cry, fpsi,
     &  ffbz, bb, vol, hx, hy, qz, qc, gs, pbs, wbbl, wbbr, wbbv)
      call alloc_b2mod_tallies(nnreg,ns)
      call get_jsep(nx,ny,jxi,jxa,jsep)
      if (bb(jxa,jsep,2).ne.0.0_R8) then
        signmf = bb(jxa,jsep,2)/abs(bb(jxa,jsep,2))
      else
        signmf = 0.0_R8
        call xerrab ('error in b2mndr , signmf = 0 in b2mndr')         !srv 09.10.09
      endif
      write(*,*) 'conn'
      do iy=-1,ny
         conn(iy)=0.0_R8
         do ix=-1,nx
            if(mod(region(ix,iy,0),4).eq.1) conn(iy)=conn(iy)+
     1           hx(ix,iy)*bb(ix,iy,3)/abs(bb(ix,iy,0))
         enddo
         write(*,*) iy,conn(iy)
      enddo
*  .. calculate 1/B/B                                                    !srv {
      do iy=-1,ny
        do ix=-1,nx
          OnedBsq(ix,iy)=1.0_R8/bb(ix,iy,3)**2
        enddo 
      enddo
*  .. interpolate OnedBsq for corners of cell
      call my_intcor (nx,ny,vol,OnedBsq,OnedBsqc)                        !srv }
*  .. interpolate magnetic field to cell faces                           !xpb
c!!!This interpolation scheme does not guarantee that bb(3)=sqrt(bb(0:2)**2)!!!
*   ..extrapolate magnetic field to edges
      call my_intcor1(nx,ny,2,vol,bb,wbbc)                               !xpb
*   ..compute extrapolated wbbc(:,:,3) component
      wbbc(:,:,3)=sqrt(wbbc(:,:,0)**2+wbbc(:,:,1)**2+wbbc(:,:,2)**2)
*   ..compute hz                                                     !srv 09.01.01  09.06.08 {
c      hz = 1.0e0_R8*(1.0e0_R8-fac_hz)+fac_hz*vol/hx/hy
!      hz = 1.0e0_R8*(1.0e0_R8-fac_hz)+fac_hz*vol/gs(:,:,2)           !srv 09.01.01
	hzconst = 2.0_R8*pi*6.0_R8   ! for ITER for testing
	hzconst = 2.0_R8*pi*1.6_R8   ! for ASDEX-UPGRADE for testing
	hzconst = 1.0_R8                                               !srv 30.01.10
      hz = hzconst*(1.0e0_R8-fac_hz)+fac_hz*vol/gs(:,:,2)            !srv 16.06.09
*   ..compute pbshz (magnetic field pitch * surface area  at cell faces * hz)
c      wrk0 = hz*vol
c      call b2xgbs (nx, ny, bb, wrk0, hx, hy, qz, pbshz)
      do iy=-1,ny                                                    !srv 11.03.08 {
        do ix=-1,nx
	   if(leftix(ix,iy).ne.-2) then
          pbshz(ix,iy,0)=pbs(ix,iy,0)*
     &      (hz(ix,iy)*vol(leftix(ix,iy),leftiy(ix,iy))
     &      +hz(leftix(ix,iy),leftiy(ix,iy))*vol(ix,iy))
     &      /(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy)) 
         else
	    pbshz(ix,iy,0)=0.0e0_R8 
	   endif             
        enddo 
      enddo    
      do iy=-1,ny
        do ix=-1,nx
	   if(bottomiy(ix,iy).ne.-2) then
          pbshz(ix,iy,1)=pbs(ix,iy,1)*
     &      (hz(ix,iy)*vol(bottomix(ix,iy),bottomiy(ix,iy))
     &      +hz(bottomix(ix,iy),bottomiy(ix,iy))*vol(ix,iy))
     &      /(vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy)) 
         else
	    pbshz(ix,iy,1)=0.0e0_R8 
	   endif             
        enddo 
      enddo                                                          !srv 11.03.08 }
*   ..find the poloidal number of midplane cell as a number of cell with min Bz                         
         if(nncut.eq.2) then                                         !srv 22.05.09 {
	    nmdpl = rightcut(nncut)
	    do ix = rightcut(nncut)+1, rightcut(nncut-1)-1
	     if(abs(bb(ix,ny,2)) .lt. abs(bb(nmdpl,ny,2))) nmdpl = ix  !srv 11.01.13
	    enddo
	   elseif(nncut.eq.1) then
	    nmdpl = leftcut(nncut)
	    do ix = leftcut(nncut)+1, rightcut(nncut)-1
	     if(abs(bb(ix,ny,2)) .lt. abs(bb(nmdpl,ny,2))) nmdpl = ix  !srv 11.01.13
	    enddo
	   else
	    write(*,*) ' b2mndr: nncut = ', nncut
	    nmdpl = (nx+2)/2
	    do ix = -1, nx
	     if(abs(bb(ix,ny,2)) .lt. abs(bb(nmdpl,ny,2))) nmdpl = ix  !srv 11.01.13
	    enddo
	   endif 
	   write(*,"(' Number of midplane cell =',i4)") nmdpl          !srv 22.05.09 }
*   ..compute length along y-direction                               !srv 24.04.02 {
      lengthy(:,-1) = hy(:,-1)*0.5e0_R8
      do iy = 0, ny
       do ix = -1, nx
        lengthy(ix,iy) = lengthy(ix,bottomiy(ix,iy)) + 
     *      (hy(ix,iy)+hy(bottomix(ix,iy),bottomiy(ix,iy)))*0.5e0_R8
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
	do ix=-1,nx
       crdy_surf(ix,-1:ny) = len_surf(ix,-1:ny) - len_surf(ix,topcut(1))
	enddo
*
	do ix=-1,nx
	 do iy = -1, ny-1
        crdy_c_cell(ix,iy) = (crdy_surf(ix,iy) + 
     &                        crdy_surf(topix(ix,iy),topiy(ix,iy)))/2
	 enddo
	 iy = ny
       crdy_c_cell(ix,iy) = crdy_surf(ix,iy)
	enddo                                                          !srv 02.05.09 }
*
*   ..compute Bz/B**2 - Bz/<B**2>
      call b2xbzb (nx, ny, bb, vol, bzb)                          !srv 13.10.06  09.06.08 23.07.10}
*   ..read initial state
      call b2rups (ninp(3), nx, ny, ns0, ne, na, ua, uadia, te, ti, po,  !xpb
     & fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy, 
     & b2fstati_version)
*   ..test plasma state and fluxes
      call b2xvps (nx, ny, ns0, 
     & ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
      call b2xvsg (n2*ns, kinrgy, 1, 'kinrgy', '.ge.')
*   ..compute log-log electron rate coefficients
      call b2spel (nx, ny, ns0, ev, te, ne, 
     &     rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi)
*   ..compute electron rate coefficients
      call b2sqel (nx, ny, ns0, ismain, ev, te, 
     &    rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &    rsa, rra, rqa, rrd, rbr, rza, rz2, rpt, rpi, wrk0)
*   ..extend initial state to species (ns0:ns-1)
      if (ns0.lt.ns) then
       call sfill (n2*(ns-ns0), na0eps, na(-1,-1,ns0), 1)
       call sfill (n2*(ns-ns0), 0.0_R8, ua(-1,-1,ns0), 1)
       call sfill (n2*(ns-ns0), 0.0_R8, kinrgy(-1,-1,ns0), 1)            !xpb
       call sfill (n2*2*(ns-ns0), 0.0_R8, uadia(-1,-1,0,ns0), 1)         !xpb
       call sfill (n2*2*(ns-ns0), 0.0_R8, fna(-1,-1,0,ns0), 1)
      endif
*   ..rescale densities
      if (density_rescale.ne.1.0_R8) then
csw 03oct2013 (this just quick and dirty for rescaling Helium, species 3 and 4
c   --> better way doing this is necessary)
        if(density_rescale < 0.d0) then
          do is=3,4
           na(:,:,is)=na(:,:,is)*(-density_rescale)
           fna(:,:,0:1,is) = fna(:,:,0:1,is)*(-density_rescale)
           kinrgy(:,:,is)=kinrgy(:,:,is)*(-density_rescale)
           write(*,*) 'Density ',is,' rescaled by ',-density_rescale
          enddo
        else
         na = na * density_rescale
         ne = ne * density_rescale
         fna = fna * density_rescale
         kinrgy = kinrgy * density_rescale
         write(*,*) 'All densities rescaled by ',density_rescale
        endif
      endif
*   ..rescale neutrals
      if(ncall.eq.0.and.neutral_rescale.ne.1.0_R8) then
        do is=0,ns-1
          if(is_neutral(is)) then
            na(:,:,is)=na(:,:,is)*neutral_rescale
            fna(:,:,0:1,is)=fna(:,:,0:1,is)*neutral_rescale
            kinrgy(:,:,is)=kinrgy(:,:,is)*neutral_rescale
          endif
        enddo
        write(*,*) 'Rescaled neutrals by ',neutral_rescale
      endif
*   ..if we rescaled or added species, recompute dependent quantities
*   ..until convergence
      if (ns0.lt.ns .or. density_rescale.ne.1.0_R8) then
        if (ns0.lt.ns) then
          do is = ns0, ns-1
            rza(-1:nx,-1:ny,is) = (zamin(is)+zamax(is))/2.0_R8
            rz2(-1:nx,-1:ny,is) = ((zamin(is)+zamax(is))/2.0_R8)**2
            rpt(-1:nx,-1:ny,is) = (potmin(is)+potmax(is))/2.0_R8
            rpi(-1:nx,-1:ny,is) = (potimin(is)+potimax(is))/2.0_R8
          enddo
          call b2xpne (nx, ny, ns, rza, na, ne)
        endif
        max_delta = 1.0_R8
        is = 0
        do while (max_delta.gt.1.0e-6_R8 .and. is.lt.10)
          is = is + 1
          call scopy (n2*ns, rza, 1, rza0, 1)
          call scopy (n2*ns, rz2, 1, rz20, 1)
          call scopy (n2*ns, rpt, 1, rpt0, 1)
          call scopy (n2*ns, rpi, 1, rpi0, 1)
          max_delta = max(
     1     maxval(abs((rza-rza0)/(rza+1.0_R8))),
     2     maxval(abs((rz2-rz20)/(rz2+1.0_R8))),
     3     maxval(abs((rpt-rpt0)/(rpt+1.0_R8))),
     4     maxval(abs((rpi-rpi0)/(rpi+1.0_R8))))
*    ..compute log-log electron rate coefficients
          call b2spel (nx, ny, ns, ev, te, ne, 
     &     rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi)
*   ..compute electron rate coefficients
          call b2sqel (nx, ny, ns, ismain, ev, te, 
     &     rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &     rsa, rra, rqa, rrd, rbr, rza, rz2, rpt, rpi, wrk0)
          call b2xpne (nx, ny, ns, rza, na, ne)
          write (*,*) 'Convergence of bundled rates = ',is,max_delta
        enddo
      endif
*   ..initialise previous state and fluxes
      call scopy (n2*ns, na, 1, na0, 1)
      call scopy (n2, ne, 1, ne0, 1)
      call scopy (n2*ns, ua, 1, ua0, 1)
      call scopy (n2, te, 1, te0, 1)
      call scopy (n2, ti, 1, ti0, 1)
      call scopy (n2, po, 1, po0, 1)
      call scopy (n2*2*ns, fna, 1, fna0, 1)
      call scopy (n2*2, fhe, 1, fhe0, 1)
      call scopy (n2*2, fhi, 1, fhi0, 1)
      call scopy (n2*2, fch, 1, fch0, 1)
      call scopy (n2*ns, kinrgy, 1, kinrgy0, 1)
      call scopy (n2*ns, na, 1, nap, 1)
      call scopy (n2, ne, 1, nep, 1)
      call scopy (n2*ns, ua, 1, uap, 1)
      call scopy (n2, te, 1, tep, 1)
      call scopy (n2, ti, 1, tip, 1)
      call scopy (n2, po, 1, pop, 1)
      call scopy (n2*2*ns, fna, 1, fnap, 1)
      call scopy (n2*2, fhe, 1, fhep, 1)
      call scopy (n2*2, fhi, 1, fhip, 1)
      call scopy (n2*2, fch, 1, fchp, 1)
      call scopy (n2*ns, kinrgy, 1, kinrgyp, 1)
*   ..compute initial ni
      call b2xpni (nx, ny, ns, na, ni)
*   ..compute ni0
      call b2xpni (nx, ny, ns, na0, ni0)
*   ..write output control
      write (nout(0),'(2(/2x,a,2x,a))')
     &  'idout0:', idout0, 'idout1:', idout1
*   ..write geometric information
      call b2wfgm (nout(0), idout0, nx, ny)
*   ..write summary of initial state
      call b2wfpi (nout(0), idout0, nx, ny, ns, ne, na, ua,
     &  te, ti, po, fna, fhe, fhi, fch)

csw 21feb2014, setup some additional work variables for boundary conds
      iffboundx=0
      iffboundy=0
      if(l_fbound .ne. 0) then
       do ix=-1,nx
        do iy=-1,ny
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            iffboundy(ix,iy)=1
          endif
        enddo
       enddo
       do iy=-1,ny,ny+1
         iffboundx(:,iy)=1
       enddo 
       ! core boundary special treatment
       do ib=1,nbc
         if(bcchar(ib).eq.'S') then
           do ix=bcstart(ib), bcend(ib)
             if(mod(region(ix,bcpos(ib),0),4).eq.1) then
               iy=bcpos(ib)
               iffboundx(ix,iy)=1
             endif
           enddo
         endif
       enddo
      endif

      ifenebound=0
      if(l_enebound .ne. 0) then
       do ix=-1,nx
        do iy=-1,ny
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            ifenebound(ix,iy)=1
          endif
        enddo
       enddo
       do iy=-1,ny,ny+1
         ifenebound(:,iy)=1
       enddo 
       ! core boundary special treatment
       do ib=1,nbc
         if(bcchar(ib).eq.'S') then
           if(mod(region(bcstart(ib),bcpos(ib),0),4).eq.1) then
             iy=bcpos(ib)
             do ix = bcstart(ib), bcend(ib)
               ifenebound(ix,iy)=1
             enddo
           endif
         endif
       enddo
      endif

      iftimbound=0
      if(l_timbound .ne. 0) then
       do ix=-1,nx
        do iy=-1,ny
          iftimbound(ix,iy)=ifenebound(ix,iy)
        enddo
       enddo
      endif




* ..main loop
*   ..initialise counters
      itim = 0
      if (stim.ge.0.0_R8) tim = stim
      call xertst (etim.eq.0.0_R8 .or. etim.gt.tim,
     & 'faulty input etim')
      save_plasma_time = tim + delta_plasma_time
*   ..initialise tracing routines
      write(*,*) 'Calling alloc_b2mod_elements'
      call alloc_b2mod_elements(ns,zn,am)
      call alloc_b2mod_neutrals(ns,nstraid,sput_dst)
      natmi = 0
      do is = 0, ns-1
        if (is_neutral(is)) natmi = natmi+1
      enddo
      call xertst (natmi.le.DEF_NATM,
     .  'Increase DEF_NATM in DIMENSIONS.F and recompile !')
c new CDF movie option (1st frame at t=0)
      if(tim.ge.save_cdfmovie_time.and.
     1 delta_cdfmovie_time.gt.0.0_R8) then
        if(use_eirene.eq.0) then
#ifndef NO_CDF
          call cdfmovie(nx,ny,ns)
#else
          write(*,*) 'Compile b2mndr without NO_CDF option for movies!'
#endif
        endif
        save_cdfmovie_time=save_cdfmovie_time+delta_cdfmovie_time
      endif
*   ..no error so far
      ierr = 0
      quit=.false.
      quitexist=.false.
      quitexist_=.false.                                            !srv 13.04.11
      quit_residual=.false.                                              !jwk
      cpuval = 0.0_R8
      stack_ptr=0
    1 continue
      cpustart=x05baf()
*   ..loop test
*!!   (Note: better loop control yet to be developed)
      if (.not.quit.and.itim.lt.ntim.and.ierr.eq.0.and.
     &    (etim.gt.tim.or.etim.eq.0.0_R8)) then
*    ..guess the next state
       call b2mxnp (n2*ns, nap, na0, na)
       call b2mxnu (n2*ns, uap, ua0, ua)
       call b2mxnp (n2, nep, ne0, ne)
       call b2mxnp (n2, tep, te0, te)
       call b2mxnp (n2, tip, ti0, ti)
       call b2mxnp (n2, kinrgyp, kinrgy0, kinrgy)
       call b2mxnu (n2, pop, po0, po)
       call b2mxnu (n2*2*ns, fnap, fna0, fna)
       call b2mxnu (n2*2, fhep, fhe0, fhe)
       call b2mxnu (n2*2, fhip, fhi0, fhi)
       call b2mxnu (n2*2, fchp, fch0, fch)
       call b2trcs (nx, ny, ns, ismain, dtim, BoRiS)
csw 09oct2012 set default potential already at this point (BCs)
       if (pot_eq.eq.0) then                                           !srv 13.07.05 {
         po   =0.0e0_R8
       else if (pot_eq.eq.2) then
         po(:,:)=3.1_R8*te(:,:)/qe
       endif

*    ..perform one time step
       ok=.false.
       do while(.not.ok)
         write(*,"(1x,a,i9,e10.2,i9,i3)")
     &  'b2mndr_ok:itim,dtim,ntim,stack_ptr', itim,dtim,ntim,stack_ptr
          call b2mndt (nout, nx, ny, ns, ismain, 
     &     itim, dtim, ierr)

          if(ierr.eq.10) then
             stack_ptr=stack_ptr+1
             if(stack_ptr.gt.20) then
                write(*,*) 'Too deep in b2mndr'
                stop
             endif
             stack_ntim(stack_ptr)=ntim-itim     ! dpc 2000.06.24 was ntim+10
             stack_dtim(stack_ptr)=dtim
             dtim=dtim/10.0_R8
             ntim=itim+10
             na=na0
             ne=ne0
             ua=ua0
             te=te0
             ti=ti0
             po=po0
             fna=fna0
             fhe=fhe0
             fhi=fhi0
             fch=fch0
             kinrgy=kinrgy0
          else
             ok=ierr.eq.0
          endif
       enddo
*    ..produce data for movie
       if (mod(itim,mvinc).eq.0.and.itim/mvinc.lt.mvnum) then
        call b2mwmv (nout(3), nx, ny, ns, itim, tim, dtim)
       endif
*    ..increment
       itim = itim+1
       if (no_solve.le.0) then
         tim = tim+dtim
       endif
       inquire(file='_out',exist=outputexist)                     !srv 18.05.09
	 if(outputexist) then
        write(*,"(1x,a,i9,a,es14.6,a,i9,a,es10.2)") 'ITER ',itim,
     &     ' TIME ',tim, ' NTIM ', ntim,' DTIM ',dtim
        write(*,*) 'delta na ',maxval(abs((na-na0)/(na+na_eps)))
        write(*,*) 'delta te ',maxval(abs((te-te0)/(te+te_eps*ev)))
        write(*,*) 'delta ti ',maxval(abs((ti-ti0)/(ti+ti_eps*ev)))
        write(*,*) 'delta po ',maxval(abs((po-po0)/(abs(po)+po_eps)))
        write(*,*) 'delta ua ',
     *   maxval(abs((na*ua-na0*ua0)/(abs(na*ua)+na_eps*ua_eps)))
	 endif
       if(delta_min.gt.0.0_R8.and.delta_max.gt.0.0_R8) then
         max_delta=max(
     1    maxval(abs((na-na0)/(na+na_eps))),
     2    maxval(abs((te-te0)/(te+te_eps*ev))),
     3    maxval(abs((ti-ti0)/(ti+ti_eps*ev))),
     4    maxval(abs((po-po0)/(abs(po)+po_eps))),
     5    maxval(abs((na*ua-na0*ua0)/(abs(na*ua)+na_eps*ua_eps))))
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
       if(itim.eq.ntim.and.stack_ptr.gt.0) then
          ntim=stack_ntim(stack_ptr)+itim     ! dpc 2000.06.24 was stack_ntim(stack_ptr)
          dtim=stack_dtim(stack_ptr)
          stack_ptr=stack_ptr-1
       endif
*    ..produce time dependent file
#ifdef JET
       if(tim.ge.save_tran_time.and.delta_tran_time.gt.0.0_R8) then
         tran_num=tran_num+1
         if(tran_num.le.99) then
           write(tran_file_name,'(a4,i2)') 'tran',tran_num
           if(tran_file_name(5:5).eq.' ') tran_file_name(5:5)='0'
           call cfopen (nout(8),tran_file_name,'new','unformatted')
           call set_b2mod_ppout
           call write_b2mod_ppout(nout(8))
           close(nout(8))
           write(*,*) 'Saved ',tran_file_name
         else
           write(*,*) 'More than 99 tran files --- no more saved'
         endif
         save_tran_time=save_tran_time+delta_tran_time
       endif
#endif
c new CDF movie option
       if(tim.ge.save_cdfmovie_time.and.
     1  delta_cdfmovie_time.gt.0.0_R8) then
#ifndef NO_CDF
         call cdfmovie(nx,ny,ns)
#else
         write(*,*) 'Compile b2mndr without NO_CDF option for movies!'
#endif
         save_cdfmovie_time=save_cdfmovie_time+delta_cdfmovie_time
       endif
       if(b2time.gt.0) then
          if(mod(itim,b2time).eq.0) then
             call b2mwti(itim,tim,nx,ny,ns,ismain,ismain0,BoRiS)
          endif
       endif
*    ..tallies
       if(tally.gt.0) then
         if(mod(itim,tally).eq.0) then
            call tallies(nx,ny,ns,BoRiS)
         endif
       endif
       if (itim.eq.1) call b2trci(nx,ny,ns,ismain,dtim,BoRiS)
       call b2trace(nx,ny,ns,ismain,dtim,BoRiS)
c-rfs       write(*,*) 'fhe(nx/2,2,0)/gs(nx/2,2,0)/qc(nx/2,2)',
c-rfs     1  fhe(nx/2,2,0)/gs(nx/2,2,0)/qc(nx/2,2)
*    ..loop
       inquire(file='_quit',exist=quitexist_)                            !srv 18.05.09 13.04.11
       inquire(file='.quit',exist=quitexist)
       cpuval=x05baf()
!	 if(outputexist)
       write(*,"(2(a,es11.3))") 'b2-step-cpu =',cpuval-cpustart,
     & ' remaining time (s) =',
     &            min((ntim-itim)*(cpuval-cpustart),
     &                b2mndr_cpu-itim*(cpuval-cpustart))                 !srv 18.05.09
       if(test_residual.and.mod(itim,10).eq.0.and.no_solve.le.0) then    !jwk
          quit_residual=                                                 !jwk
     &        (areshe*min_areshe).le.(min_areshe*min_areshe).and.        !jwk
     &        (areshi*min_areshi).le.(min_areshi*min_areshi).and.        !jwk
     &        (aresco(ismain)*min_aresco).le.(min_aresco*min_aresco)     !jwk
       endif                                                             !jwk
       quit=quitexist_.or.quitexist.or.quit_residual.or.                 !jwk !srv 13.04.11
     1      (cpuval.gt.b2mndr_cpu.and.b2mndr_cpu.gt.0).or.
     2      (etim.gt.stim.and.tim.ge.etim.and.etim.ne.0.0_R8).or.
     3      (no_solve.eq.1.and.itim.eq.ntim)
       if(.not.quit.and.
     .   ((cpuval.gt.cputarget.and.cpuincrement.gt.0.0_R8).or.
     .    (tim.ge.save_plasma_time.and.delta_plasma_time.gt.0.0_R8))) 
     >   then
          write(savefile(13:16),'(i4)') edition
          edition=edition+1
          if(cpuval.gt.cputarget.and.cpuincrement.gt.0.0_R8) 
     .      cputarget=cputarget+cpuincrement
          if(tim.ge.save_plasma_time.and.delta_plasma_time.gt.0.0_R8) 
     .      save_plasma_time=save_plasma_time+delta_plasma_time
          if(edition.gt.9999) then
             edition=0
             savefilestatus='old'
          endif
          do i=13,16
             if(savefile(i:i).eq.' ') savefile(i:i)='0'
          enddo
          call cfopen (nout(6), savefile, savefilestatus, 'unformatted') 
          call cfverw (nout(6), newversion)			     !xpb/dpc
*       ..test plasma state
          call b2xvps (nx, ny, ns, 
     .     ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
          call b2xvsg (n2*ns, kinrgy, 1, 'kinrgy', '.ge.')
*       ..write plasma state
          call b2wups (nout(6), nx, ny, ns, ne, na, ua, uadia, te,   !xpb
     &     ti, po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy)
          close(nout(6))
          write(*,*) 'Plasmastate savefile ',savefile,' saved'
!          open (70, file='b2fstate_add.'//savefile(13:16),           !srv 14.07.08 {
!     &             form='unformatted')
!	     write(70) vaecrb, vadia, floe_noc, floi_noc, fna, fna_fcor,
!     &     fna_nodrift, fna_he, fhe, fhi, fch,
!     &     fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq,fchinert,
!     &     fnaPSch, fna_mdf, fhePSch, fhiPSch, fhe_mdf, fhi_mdf      !srv 29.10.09
!	    close(70)                                                  !srv 14.07.08 }
          if (ank_tracing.ne.0) call b2filef(ns)
       endif
       goto 1
*   ..end loop
      endif

      if(quitexist)
     1 write(*,*) 'stopping because .quit seen'
      if(quitexist_) write(*,*) 'stopping because _quit seen'         !srv 13.04.11
      if(cpuval.gt.b2mndr_cpu.and.b2mndr_cpu.gt.0.0_R8)
     1 write(*,*) 'stopping because cpu > limit: ',cpuval,b2mndr_cpu
      if(itim.ge.ntim)
     1 write(*,*) 'stopping because itim >= ntim: ',itim,ntim
      if(etim.gt.stim.and.tim.ge.etim.and.etim.ne.0.0_R8)
     1 write(*,*) 'stopping because tim >= etim: ',tim,etim
      if(etim.eq.stim.and.stim.gt.0.0_R8)
     1 write(*,*) 'stopping because etim = stim: ',etim,stim
      if(ierr.ne.0)
     1 write(*,*) 'stopping because ierr <> 0: ',ierr
      if(quit_residual)                                                   !jwk
     1 write(*,*) 'stopping because min_residual criteria satisfied'      !jwk

      call b2wdat (nx, ny, ns, ismain, BoRiS)                             !srv 05.10.01 13.01.13 {
#ifdef ASTRA
      call write_b2_to_astra(nx, ny, ns, ismain, nout(8))                 !srv 28.07.11 19.09.11
#endif
*   ..write out fort.30 and fort.31 files                                 !srv 01.03.10 {
      if (use_eirene.eq.0) then
       if (save_f30.ne.0) call write_f30(nx, ny)
       if (save_f31.ne.0) then
        if(background_new.eq.0) call copy_background(nx, ny, ns, BoRiS)
        if(background_new.ne.0) call copy_background_(nx, ny, ns, BoRiS)
        call write_f31(nx, ny)
       end if
      end if                                                              !srv 01.03.10 } 13.01.13 }

* finalise feedback in b2stbc
      call b2stbc_fb_final
* finalise transport_models
      call write_b2mod_transport_models
      if(use_eirene.ne.0.and.itim.gt.0) then
#ifdef USE_EIRENE
         call eirene_mc_final(nx,ny,ns,BoRiS)
#else
         call xerrab ('b2mndr not compiled with the USE_EIRENE option')
#endif
      endif
* finalise sputtering in b2stbr
      call b2stbr_sputter_final(nx, ny, ns)                        !srv 25.10.08
* finalise time-dependent wall data
      call wall_final(nx,ny,ns,plate_model,plate_temp,plate_thick) !srv 25.10.08

* ..display and save the result of the calculation
*   ..write summary of final state
      call b2wfpi (nout(0), idout1, nx, ny, ns, ne, na, ua,
     &  te, ti, po, fna, fhe, fhi, fch)
*   ..test plasma state
      call b2xvps (nx, ny, ns, 
     .  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
      call b2xvsg (n2*ns, kinrgy, 1, 'kinrgy', '.ge.')
*   ..save final state
      call b2wups (nout(2), nx, ny, ns, ne, na, ua, uadia, te,            !xpb
     &  ti, po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy)

!      open (70, file='b2fstate_add',form='unformatted')                   !srv 11.10.02 {
!	 write(70) vaecrb, vadia, floe_noc, floi_noc, fna, fna_fcor,
!     &    fna_nodrift, fna_he, fhe, fhi, fch,
!     &    fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq,fchinert,
!     &    fnaPSch, fna_mdf, fhePSch, fhiPSch, fhe_mdf, fhi_mdf       !srv 13.10.06 !srv 29.10.09  
!	close(70)                                                      !srv 11.10.02 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mndr

      end subroutine b2mndr
c
