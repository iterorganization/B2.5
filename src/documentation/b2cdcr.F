*-----------------------------------------------------------------------
*.specification

      subroutine b2cdcr ()
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2CDCR contains a brief statement of the purpose of many of the
*     subroutines of the B2.5 code. No computation is performed in this
*     routine.
*
*-----------------------------------------------------------------------
*.text
*
*  0. Overview and naming conventions.
*
*     B2A*, B2B* - Auxiliary programs and their subroutines.
*       B2AB* - Preparation of databases for advanced fluid neutral BCs
*       B2AG* - Preparation of a geometry.
*       B2AH* - Preparation of default physics parameters.
*       B2AI* - Preparation of an initial state.
*       B2AR* - Preparation of atomic rate coefficients.
*       B2AW* - Reserved for preparation of wall data.
*     B2C* - Documentation and 'include' segments.
*       B2CD* - Documentation.
*       B2CM* - 'include' segments: parameters and common blocks.
*     B2D* - Reserved for DEGAS interface.
*     B2E* - Reserved for EIRENE interface.
*     B2F*-B2L* - Not used.
*     B2M* - Core routines for main program: driver and timestepping.
*       B2MN_INIT - Main program initialisation.
*       B2MN_STEP - Main program call.
*       B2MN_FIN - Main program finalization.
*       B2MNDR_0 - Main driver routine initialisation.
*       B2MNDR_1 - Main driver routine.
*       B2MNDR_2 - Main driver routine finalisation.
*       B2MNDS - Initialisation.
*       B2MNDT - Time stepping.
*       B2MW* - Monitoring of the calculation.
*       B2MX* - Auxiliaries for b2m*.
*     B2N* - Driver routines for inner iteration.
*       B2NEWS_ - Controls inner iteration.
*       B2NEWS_M - Controls inner iteration for time-dependent mode.
*       B2NP* - Iteration for individual equations.
*       B2NX* - Auxiliary routines for b2n*.
*     B2O*-B2Q* - Reserved for alternative numerical core.
*     B2R* - Read from an external file.
*       B2RF* - Read formatted.
*       B2RU* - Read un*formatted.
*     B2S* - Evaluate source term coefficients.
*       B2SI* - Evaluate internal source terms.
*       B2SP* - Evaluate log-log linearised rate coefficients.
*       B2SQ* - Evaluate rate coefficients.
*       B2SR* - Evaluate or modify total source term coefficients.
*       B2ST* - Evaluate partial source term coefficients.
*     B2T* - Evaluate transport coefficients and fluxes.
*       B2TF* - Evaluate fluxes.
*       B2TL* - Apply flux limits to transport coefficients.
*       B2TQ* - Evaluate differential transport coefficients.
*       B2TR* - Evaluate cell-integrated transport coefficients.
*       B2TX* - Auxiliary routines for b2t*.
*     B2U* - Inner routines for the inner iteration.
*       B2UP* - Update some quantity.
*       B2UR* - Evaluate a residual.
*       B2US* - Solve a correction equation.
*       B2UX* - Auxiliaries related to b2u*.
*     B2US_* - Routines and modules related to unstructured format.
*     B2V* - Not used.
*     B2W* - Write to an external file.
*       B2WF* - Write formatted.
*       B2WG* - Write graphical.
*       B2WS* - Write special purpose.
*       B2WU* - Write un*formatted.
*     B2X* - Basic auxiliary routines.
*       B2XG* - Auxiliaries for geometry.
*       B2XP* - Auxiliaries for plasma state.
*       B2XV* - Test validity of some quantities.
*       B2XX* - Miscellaneous auxiliaries.
*     B2Y* - Postprocessing.
*       B2YG* - Display the geometry and magnetic field.
*       B2YH* - Display the table of physics parameters.
*       B2YI* - Display the plasma state.
*       B2YM* - Display the evolution of the plasma state.
*       B2YN* - Display the progress of the inner iterations.
*       B2YP* - Comprehensive postprocessing of a calculation.
*       B2YQ* - Quick display of the evolution of a calculation.
*       B2YR* - Display the table of atomic rate coefficients.
*       B2YS* - Display the externally imposed source terms.
*     B2Z* - Reserved for application-dependent postprocessing.
*
*     Routines that end in _st are remnants of the structured code versions.
*
*
*  1. Purpose of individual B2.5 routines.
*
*     B2AB (main program) prepares databases for the advanced fluid neutral
*     (AFN) boundary conditions.  This program calls the driver routine
*     b2abdr for the actual computation.
*
*     B2ABDR computes and writes to an external file the databases for the
*     advanced fluid neutral (AFN) boundary conditions.
*
*     B2AG (main program) prepares a geometry and magnetic field for
*     B2.5. This program calls the driver routine b2agdr for the actual
*     computation.
*
*     B2AGDR computes and writes to an external file the specification
*     of a B2.5 geometry and magnetic field.
*
*     B2AGBB computes the magnetic field from the flux functions.
*
*     B2AGBF computes computes the three components of the magnetic field
*     and the total magnetic field strength, all at cell faces.
*     Not converted for WG.
*
*     B2AGBS computes a correction to the pbsx array to ensure the
*     constancy of magnetic flux along flux tubes.
*     Not converted for WG.
*
*     B2AGCF computes the locations of the cell vertices, and computes
*     the flux functions at those locations.
*     Not converted for WG.
*
*     B2AGFS obtains a description of the flux surface structure.
*
*     B2AGGS marks a singular region in the geometry.
*
*     B2AGMT computes the metric coefficients given the locations of
*     cell vertices.
*     Not converted for WG.
*
*     B2AGMX projects the metric coefficients from a fine mesh to a
*     coarsening of that mesh.
*
*     B2AGSP evaluates the coordinates and the flux functions at one
*     single point on the mesh. This routine is application dependent,
*     and defines the specific geometry and magnetic field for a B2.5
*     calculation.
*     Not converted for WG.
*
*     B2AGX0 reads the mesh dimensions from a formatted file.
*
*     B2AH (main program) prepares default physics parameters for B2.5.
*     This program calls the driver routine b2ahdr for the actual
*     computation.
*
*     B2AHDR prepares and writes to an external file a table of default
*     physics parameters for B2.5.
*
*     B2AHX0 reads the number of species from a formatted file.
*
*     B2AI (main program) prepares an initial plasma state for B2.5.
*     This program calls the driver routine b2aidr for the actual
*     computation.
*
*     B2AIDR prepares and writes to an external file an initial plasma
*     state for B2.5.
*
*     B2AIX0 reads the mesh dimensions and the number of species from a
*     formatted file.
*
*     B2AR (main program) prepares a table of atomic rate coefficients
*     for B2.5. This program calls the driver routine b2ardr for the
*     actual computation.
*
*     B2ARDR computes and writes to an external file a table of atomic
*     rate coefficients for B2.5. The coefficients are obtained using
*     either the ADPAK or the STRAHL set of routines.
*
*     B2ARIN reads from a formatted file the list of nuclear species
*     and atomic states for which the rate coefficients are to be
*     computed.
*
*     B2CDCA provides some general documentation of the B2.5 code.
*     No computation is performed in this routine.
*
*     B2CDCR contains a brief statement of the purpose of many of the
*     subroutines of the B2.5 code. No computation is performed in this
*     routine.
*
*     B2CDCV contains a description of some of the variables used in
*     the B2.5 code. No computation is performed in this routine.
*
*     B2MN is the main program for the B2.5 code. The code solves a
*     system of fluid equations that models a multi-species plasma
*     in a two-dimensional geometry. This program calls the driver
*     routines b2mndr_[0-2] for the actual computation.
*
*     B2MNDR is the main driver routine for the B2.5 code.
*
*     B2MNDS performs initialisation functions for a B2.5 calculation.
*
*     B2MNDT performs one implicit timestep for the B2.5 system of
*     equations.
*
*     B2MWMV writes to an un*formatted file a selection of data
*     describing the state of the plasma.
*     Not converted for WG.
*
*     B2MWIT writes to an un*formatted file detailed monitoring output
*     for the inner iterations of the B2.5 solution process.
*     Not converted for WG.
*
*     B2MWQ0 writes a header for the output from b2mwqt.
*
*     B2MWQT writes brief data concerning the numerical evolution.
*
*     B2MXAC computes the (maximum-)norms of the corrections.
*
*     B2MXAR computes scaled (L1-)norms of the residuals.
*
*     B2MXNP extrapolates a positive quantity to the next time.
*
*     B2MXNU extrapolates an indefinite quantity to the next time.
*
*     B2MXZR sets the residuals to zero along the boundaries of the
*     mesh and on the special regions.
*
*     B2NEWS_ performs one inner iteration for the B2.5 equations.
*     This routine defines the core structure of the solution
*     algorithm.
*
*     B2NEWS_m performs one inner iteration for the B2.5 equations,
*     aiming at a correct time-dependent behavior.
*     This routine defines the core structure of the solution
*     algorithm.
*
*     B2NPCO performs one iteration of a relaxation process for the
*     solution of the continuity equation for one atomic species.
*
*     B2NPHT performs one iteration of a relaxation process for the
*     solution of the equations for internal electron and atom energy.
*
*     B2NPMO performs one iteration of a relaxation process for the
*     solution of the B2.5 equations for parallel velocity.
*
*     B2NPPO performs one iteration of a relaxation process for the
*     solution of the charge conservation equation, which governs the
*     electric potential.
*
*     B2NXCM computes the adjustment to the coefficients pcca that
*     comes from the combined parallel momentum equation solver.
*
*     B2NXDP computes the correction to the electric potential following
*     an update to the electron density.
*
*     B2NXDU computes the artificial velocity shift that is
*     required to make the total numerical friction force cancel.
*
*     B2NXDV computes the correction to the parallel momentum source
*     that is required to make the total friction force cancel.
*
*     B2NXFC computes the transport terms in the parallel momentum
*     correction equation.
*
*     B2NXFX computes an expression that occurs in the electron-atom
*     thermal force.
*
*     B2NXST accumulates the parallel momentum source linearisation
*     over the distinct atomic species.
*
*     B2PRNTI prints out the residuals in integer format.
*     Not converted for WG.
*
*     B2RFCP reads from a formatted file some parameters in the
*     physics common blocks.
*
*     B2RFLB reads from a formatted file a label for one run of one
*     of the programs that make up B2.5.
*
*     B2RUCP reads from an un*formatted file the parameters in the
*     physics common blocks.
*
*     B2RUGM reads from an un*formatted file the geometry and magnetic
*     field.
*
*     B2RUPS reads the state of the plasma from an un*formatted file.
*
*     B2RURC reads a table of atomic rate data from an un*formatted
*     file.
*
*     B2RUSR reads a table of source term data from an un*formatted
*     file.
*
*     B2SCDT modifies the source coefficients to account for the
*     implicit time discretisation for the isb-th density equation.
*
*     B2SHDT modifies the source coefficients to account for the
*     implicit time discretisation for Te and Ti equations.
*
*     B2SIAN computes contribution to the parallel momentum source due
*     to the anomalous current.
*
*     B2SIAV computes a linearised representation of the additional
*     viscosity.
*
*     B2SICF computes centrifugal force.
*
*     B2SIFR computes a linearised representation of the friction
*     force between one species and all others.
*
*     B2SIFR_ computes a linearised representation of the friction
*     force between one species and all others.
*
*     B2SIFRTF computes a linearised representation of the thermal
*     and friction forces.
*
*     B2SIGP computes the pressure gradient term on the right hand
*     of the parallel momentum balance equation for one species.
*
*     B2SIHS_ computes the source coefficients for a set of terms in the
*     heat equations. These terms specify the contributions of the
*     divergence of the flow field, of viscous heating, and of the
*     electron-atom and atom-atom friction.
*
*     B2SIKT computes the sources in the energy equations due to the
*     k-eps model.
*
*     B2SMDT modifies the source coefficients to account for the
*     implicit time discretisation for parallel momentum equations.
*
*     B2SPCX computes on all cells of the mesh a log-log representation
*     of the rate coefficients for charge exchange.
*
*     B2SPEL computes on all cells of the mesh a log-log representation
*     of the rate coefficients for ionisation, recombination, and
*     electron heat loss.
*
*     B2SQCX computes on all cells of the mesh the rate coefficients
*     for charge exchange with one neutral hydrogen species.
*
*     B2SQEL computes on all cells of the mesh the rate coefficients
*     for ionisation, recombination, and electron heat loss.
*
*     B2SRAL provides linearised expressions for the sources of
*     density, parallel momentum, electron heat, atom heat, and charge.
*
*     B2SRDT modifies the source coefficients to account for the
*     implicit time discretisation.
*
*     B2SRSM smoothes expressions the sources strengths of
*     density, parallel momentum, electron heat, atom heat, and charge.
*
*     B2SRST modifies the source coefficients to provide numerical
*     stabilisation.
*
*     B2STBC provides linearised expressions for the sources of
*     density, parallel momentum, electron heat, atom heat, and charge,
*     due to the boundary conditions.
*
*     B2STBC_FB provides linearised expressions for the sources of
*     density, parallel momentum, electron heat, atom heat, and charge,
*     due to feedback boundary conditions after the normal boundary
*     conditions have been applied.
*
*     B2STBC_PHYS provides linearised expressions for the sources of
*     density, parallel momentum, electron heat, atom heat, and charge,
*     due to the boundary conditions.
*     Here the boundary conditions are expressed by use of the old B2
*     namelist format.
*
*     B2STBM provides linearised expressions for the sources of
*     density, parallel momentum, and heat, due to miscellaneous
*     application specific processes.
*
*     B2STBR provides linearised expressions for the sources of
*     density, parallel momentum, and heat, due to special boundary
*     conditions and the neutral gas recycling model.
*
*     B2STBR_PHYS provides linearized expressions for the sources of
*     density, parallel momentum and heat, due to special boundary
*     conditions and the neutral gas recycling model, using the neutrals
*     namelist. Also the advanced fluid neutral (AFN) boundary conditions are computed here.
*
*     B2STCX provides linearised expressions for the sources of atom
*     density and atom parallel momentum, due to charge exchange.
*
*     B2STEL provides linearised expressions for the sources of
*     density, parallel momentum, electron heat and atom heat, due
*     to the electron-atom collision processes of ionisation,
*     recombination, and radiation.
*
*     B2TANML computes the X and Y components of the anomalous current.
*
*     B2TBFL computes the total flux through the boundaries.
*
*     B2TCPA computes the parallel current between neighboring cells.
*
*     B2TDIA computes the X and Y components of the divergent part of
*     the diamagnetic current.
*
*     B2TEPSCH computes the X and Y components of the Pfirsch-Schlueter
*     heat flux.
*
*     B2TFCC adds contribution from cdpa to floe, cone, floi, coni
*
*     B2TFCH computes the electric current between neighboring cells,
*     and the coefficients that specify its linearisation.
*
*     B2TFED computes the electron drift velocities.
*
*     B2TFHE_ computes the electron heat flux between neighboring cells,
*     and the coefficients that specify its linearisation.
*
*     B2TFHI_ computes the all atom heat flux between neighboring cells,
*     and the coefficients that specify its linearisation.
*
*     B2TFNB computes the particle flux of one atomic species between
*     neighboring cells, and the coefficients that specify its
*     linearisation.
*
*     B2TFRN computes the particle and energy fluxes between
*     neighbouring cells which are to be provided to Eirene.

*     B2TINER computes the inertial current.
*
*     B2TINNT computes the ion-neutral current.
*
*     B2TLC0 computes the neutral diffusion flux limit.       
*
*     B2TLHE applies a flux limit to the electron heat and charge
*     transport coefficients.
*
*     B2TLHI applies a flux limit to the all atom heat transport
*     coefficients.
*
*     B2TLH0 computes the conductive neutral heat transport flux limit.
*
*     B2TLMV applies a flux limit to the viscosity.
*
*     B2TLNL computes the Coulomb logarithm according to Braginskii or 
*     Wesson formulas.
*
*     B2TLV0 computes the neutral viscosity flux limit.
*
*     B2TQCA computes the classical contributions to the atom
*     transport coefficients.
*
*     B2TQCE computes the classical contributions to the electron
*     transport coefficients.
*
*     B2TQIN computes the classical ion-neutral perpendicular conductivity.
*
*     B2TQNA computes the anomalous contributions to the transport
*     coefficients.
*
*     B2TRAL computes the transport coefficients, including geometric
*     factors.
*
*     B2TRCL computes the classical transport coefficients, including
*     geometric factors.
*
*     B2TREQ computes the cell-integrated energy equipartition
*     multiplier.
*
*     B2TRNO computes the anomalous transport coefficients, including
*     geometric factors.
*
*     B2TRQL computes the flux limit velocities.
*
*     B2TTIA computes the classical Braginskii time of the ion-atom collision.
*
*     B2TVSKT computes the current due to RS of kt-model.
*
*     B2TVSPA computes the X and Y components of the parallel viscosity
*     current.
*
*     B2TVSPR computes the X and Y components of the perpendicular viscosity
*     current.
*
*     B2TVSQ computes the X and Y components of the heat viscosity
*     current.
*
*     B2TXCC computes normal conductances between neighbouring cells, 
*     given the conductivities.
*
*     B2TXCV computes viscous conductances between neighbouring
*     staggered cells, given the viscosities (not used anymore).
*
*     B2TXCX computes conductances between neighbouring cells in the
*     poloidal direction, given the conductivities.
*
*     B2TXCY computes conductances between neighbouring cells in the
*     radial direction, given the conductivities.
*
*     B2TXFX computes conductances between neighbouring cells in the
*     poloidal direction, given the flux-limited conductivities.
*     Not converted for WG.
*
*     B2TXFY computes conductances between neighbouring cells in the
*     radial direction, given the flux-limited conductivities.
*     Not converted for WG.
*
*     B2TXSX computes the integrals over cell x-faces of a function
*     that is provided on the cells.
*
*     B2TXSY computes the integrals over cell y-faces of a function
*     that is provided on the cells.
*
*     B2UPCO corrects the state of the plasma towards the solution of
*     continuity equation for one atomic species.
*
*     B2UPHT adjusts the state of the plasma towards the solution
*     of the heat equations.
*
*     B2UPPO corrects the electric potential po.
*
*     B2URMO computes the residual of the parallel momentum balance
*     equation for one atomic species.
*
*     B2URSC computes the residual of a scalar conservation equation.
*
*     B2URSD computes the residual of a scalar conservation equation.
*
*     B2USCO computes a correction towards the solution of the
*     continuity equation for one atomic species.
*
*     B2USHT computes a correction to the fields te and ti.
*
*     B2USMO Computes a correction towards the solution of the
*     parallel momentum balance equation for one atomic species.
*
*     B2USPO computes a correction towards the solution of the
*     potential equation.
*
*     B2UXUS solves a system of equations on an unstructured mesh. This routine
*     is an interface to an iterative or a direct solver; the choice
*     depends upon an internal parameter.
*
*     B2WAOS prints out the results of a B2.5 calculation in old-fashioned
*     style.
*     Not converted for WG.
*
*     B2WDAT is an output routine that writes out B2.5 arrays into ASCII
*     data files.
*
*     B2WFCP writes to a formatted file the parameters in the physics
*     common blocks.
*
*     B2WFGM writes to a formatted file a summary of the B2.5 geometry
*     and magnetic field.
*
*     B2WFPI writes a summary of the plasma state to a formatted file.
*     This routine is an interface to b2wfpj.
*
*     B2WFPJ writes a summary of the plasma state to a formatted file.
*     This routine does not make use of any geometric information.
*
*     B2WFPK writes supplementary information about the plasma state to
*     a formatted file.
*
*     B2WFPS writes a summary of the plasma state to a formatted file.
*
*     B2WFR0 writes to a formatted file a text containing credits for
*     the ADPAK atomic physics rate coefficients.
*
*     B2WFR1 writes to a formatted file a text containing credits for
*     the STRAHL atomic physics rate coefficients.
*
*     B2WFRC writes to a formatted file a table of atomic rate
*     coefficients that has been prepared for B2.5.
*
*     B2WGCR provides graphical output of the residuals and the
*     corrections.
*
*     B2WGDT provides graphical output of the change in the state.
*
*     B2WGGM produces graphical displays of the B2.5 geometry and
*     magnetic field.
*
*     B2WGPI provides graphical output of the plasma state on a uniform
*     rectangular mesh. This routine is an interface to b2wgpj.
*
*     B2WGPJ provides graphical output of the plasma state on a uniform
*     cartesian mesh.
*
*     B2WGPS provides graphical output of the plasma state.
*     This routine is an interface to b2wgst.
*
*     B2WGST provides graphical output of the plasma state.
*
*     B2WGRC displays a table of atomic rate coefficients that has been
*     prepared for B2.5.
*
*     B2WGX0 displays the B2.5 mesh in a cartesian representation.
*
*     B2WGX1 displays the B2.5 mesh on the physical plane.
*
*     B2WGX2 produces contour plots on a sequence of cartesian domains.
*
*     B2WGX3 produces contour plots on the physical domain.
*
*     B2WGX4 produces log-contour plots on a sequence of cartesian
*     domains.
*
*     B2WGX5 produces log-contour plots on the physical domain.
*
*     B2WGXB plots an outline of the B2.5 mesh on the physical domain.
*
*     B2WGXP produces contour plots on a cartesian domain.
*
*     B2WGXQ produces log-contour plots on a cartesian domain.
*
*     B2WGYP performs the core NCAR operations to produce a contour
*     plot.
*
*     B2WGYQ performs the core NCAR operations to produce a log-contour
*     plot.
*
*     B2WSAV writes the state of the plasma to an unformatted file,
*     in the style required for input into the 'old B2' code.
*     The routine operates as an interface to b2wsax.
*
*     B2WSAX writes the state of the plasma to an unformatted file,
*     in the style required for input into the 'old B2' code.
*
*     B2WSCP provides special purpose printed output of data in the
*     physics common blocks.
*
*     B2WSPS provides special purpose printed output concerning the
*     plasma state. This routine is an interface to b2wsst.
*
*     B2WSST provides special purpose printed output concerning the
*     plasma state.
*
*     B2WUCP writes to an un*formatted file the parameters in the
*     physics common blocks.
*
*     B2WUGM writes to an un*formatted file the geometry and magnetic
*     field.
*
*     B2WUPS writes the state of the plasma to an un*formatted file.
*
*     B2WURC writes a table of atomic rate coefficients to an
*     un*formatted file.
*
*     B2XBZB calculates Bz/B**2 - Bz/<B**2>.
*
*     B2XGBS computes the product (magnetic field pitch)*(surface area)
*     at cell faces.
*     Obsolete for WG.
*
*     B2XGXY computes cartesian coordinates suitable for graphical and
*     printed output.
*
*     B2XPFE computes the electron particle flux, fne:
*       fne(iFc,0:1) = (sum is :: rza(iFc,is)*fna(iFc,0:1,is))-(1/qe)*fch(iFc,0:1).
*
*     B2XPFI computes the total atom particle flux, fni:
*       fni(iFc,0:1) = (sum is :: fna(iFc,0:1,is)) .
*
*     B2XPFZ computes the total ion particle flux, fni:
*       fni(iFc,0:1) = (sum is, .not.is_neutral(is) :: fna(iFc,0:1,is)) .
*
*     B2XPIC computes the poloidal projection of the ion parallel
*     current on cell centers.
*
*     B2XPIF computes the ion parallel current on cell faces.
*
*     B2XPMC computes the mach number of the parallel flow.
*
*     B2XPNE computes the electron density, ne:
*       ne(iCv) = (sum is :: rza(iCv,is)*na(iCv,is)) .
*
*     B2XPNI computes the total atom density, ni0:
*       ni(iCv,0) = (sum is :: na(iCv,is)) .
*     and the total ion density, ni1:
*       ni(iCv,1) = (sum is :: na(iCv,is) if is_neutral(is).eqv..false.) .
*
*     B2XPNM computes the atom density weighted by rz2/am, ne2m:
*       ne2m(,) = (sum is :: (rz2(,,is)/am(is))*na(,,is)) .
*
*     B2XPNR computes the atom density weighted by 1/sqrt(am), nirm:
*       nirm(,) = (sum is :: na(,,is)/sqrt(am(is))) .
*
*     B2XPPB computes the partial pressure, pb = nb*(rzb*te+ti).
*
*     B2XPPE computes the electron pressure, pe = ne*te.
*
*     B2XPPR computes the total pressure, pr:
*       pr(iCv) = ne(iCv)*te(iCv) + (sum is :: na(iCv,is)*ti(iCv)) .
*
*     B2XPPZ computes the plasma pressure, pz:
*       pz(iCv) = ne(iCv)*te(iCv) + (sum is :: na(iCv,is)*ti(iCv), not including neutrals).
*
*     B2XPRO computes the mass density, ro:
*       ro(iCv) = (sum is :: am(is)*mp*na(iCv,is)) .
*
*     B2XPRZ computes the ion mass density, rz:
*       rz(iCv) = (sum is :: am(is)*mp*na(iCv,is), including only charged species).
*
*     B2XPVE computes the parallel electron velocity on cell centers.
*
*     B2XPXC interpolates a function from x-faces to cell centers.
*     Obsolete for WG.
*
*     B2XVCP tests the validity of parameters in the physics common
*     blocks.
*
*     B2XVFF verifies the zero entries in the edge locations of some
*     quantity defined on faces of the B2.5 mesh. This routine is an
*     interface to b2xvfx and b2xvfy.
*     Obsolete for WG.
*
*     B2XVFV verifies the zero entries in the edge locations of some
*     quantity that is defined in the manner of a viscosity on the
*     cells and corners of the B2.5 mesh (not used anymore).
*     Obsolete for WG.
*
*     B2XVFX verifies the zero entries in the edge locations of some
*     quantity defined on x-faces of the B2.5 mesh.
*     Obsolete for WG.
*
*     B2XVFY verifies the zero entries in the edge locations of some
*     quantity defined on y-faces of the B2.5 mesh.
*     Obsolete for WG.
*
*     B2XVGM tests the validity of the geometric parameters.
*     Not converted for WG.
*
*     B2XVPS tests the validity of the state of the plasma.
*
*     B2XVSG verifies the sign of some vector quantity.
*
*     B2XXGS sets a mesh function to a constant value on the singular
*     regions of the B2.5 geometry.
*     Not converted for WG.
*
*     B2XXID returns an identification of the current version of B2.5.
*
*     B2XXMM computes the locations of the minimum and maximum entries
*     of a quantity defined on the B2.5 mesh.
*     Not converted for WG.
*
*     B2XPEN computes contributions to the total energy fluxes,
*     as well as the total energy fluxes
*     The contributions from external species are not included and must
*     be added explicitly in the calling routine when needed.
*
*     B2XZEF computes zeff:
*       zeff(iCv) = (sum is :: rz2(iCv,is)*na(iCv,is))/ne(iCv) .
*
*     B2XZMF computes zmff:
*       zmff(iCv) = (sum over is /= ismain :: rz2(iCv,is)*na(iCv,is))/
*                  na(iCv,ismain)*rz2(iCv,ismain)
*
*     B2XZDD zeroes out the array elements corresponding to dead cells
*      and core boundary cells.
*
*     B2XEHX computes the modified poloidal electrical field on all cell
*     faces.
*
*     B2XEHY computes the modified radial electrical field on all cell
*     faces.
*
*     B2XPFN computes the total particle flux of hydrogenic atoms, fnn:
*       fnn(iFc,0:1) = (sum is, is_neutral(is).and.zn(is).eq.1 :: fna(iFc,0:1,is)).
*
*     B2XPNN computes the total density of hydrogenic neutrals, nn:
*       nn = (sum is, is_neutral(is).and.zn(is).eq.1 :: na(,is)).
*
*     B2YE (main program) is a postprocessor for B2.5. It calculates
*     statistical errors based on averaged quantities collected during
*     the averaging phase of a B2.5 calculation.
*
*     B2YG (main program) is a postprocessor for B2.5. It reads and
*     displays the geometry and magnetic field. This program calls the
*     driver routine b2ygdr for the actual computation.
*
*     B2YGDR is a postprocessor for B2.5. It reads and displays the
*     geometry and magnetic field.
*
*     B2YH (main program) is a postprocessor for B2.5. It reads and
*     displays the table of physics parameters. This program calls the
*     driver routine b2yhdr for the actual computation.
*
*     B2YHDR is a postprocessor for B2.5. It reads and displays the table
*     of physics parameters.
*
*     B2YI (main program) is a postprocessor for B2.5. It reads and
*     displays the plasma state. This program calls the driver routine
*     b2yidr for the actual computation.
*
*     B2YIDR is a postprocessor for B2.5. It reads and displays the
*     plasma state. This routine does not make use of any metric
*     information.
*
*     B2YM (main program) is a postprocessor for B2.5. It is meant to
*     produce movie output. This program calls the driver routine b2ymdr
*     for the actual computation.
*
*     B2YMDR is a postprocessor for B2.5. It is meant to produce movie
*     output.
*
*     B2YN (main program) is a postprocessor for B2.5. It displays the
*     progress of the inner iterations. This program calls the driver
*     routine b2yndr for the actual computation.
*
*     B2YNDR is a postprocessor for B2.5. It displays the progress of
*     the inner iterations.
*
*     B2YP (main program) is a postprocessor for B2.5. It displays the
*     plasma state. This program calls the driver routine b2ypdr for
*     the actual computation.
*
*     B2YPDR is a postprocessor for B2.5. It displays the plasma state,
*     making use also of the geometric information.
*
*     B2YQ (main program) is a postprocessor for B2.5. It provides a
*     quick display of the evolution data. This program calls the driver
*     routine b2yqdr for the actual computation.
*
*     B2YQDR is a postprocessor for B2.5. It provides a quick display of
*     the evolution data.
*
*     B2YR (main program) is a postprocessor for B2.5. It reads and
*     displays a table of atomic rate coefficients. This program calls
*     the driver routine b2yrdr for the actual computation.
*
*     B2YRDR is a postprocessor for B2.5. It reads and displays a table
*     of atomic rate coefficients.
*
*
*  2. Purpose of some additional routines employed with B2.5.
*
*   ..2a program management and error handling
*
*     PRGINI performs any tasks required at the start of program
*     execution.
*
*     PRGEND performs any tasks required at the end of program
*     execution.
*
*     SUBSYS is a group of routines (entry points in SUBSYS) that
*     perform certain functions related to workspace management and
*     error handling. The entry points are subini, subend, and subxnd.
*     SUBINI performs any tasks required at the start of subprogram
*     execution.
*     SUBEND performs any tasks required at the end of subprogram
*     execution.
*     SUBXND is employed to build down the stack when the program
*     is terminated because of an error condition.
*
*     IPMAIN is a group of routines that
*     operate on a collection of names and values representing
*     "internal parameters" in a code. The entry points are ipsetc,
*     ipgetc, and ipfini.
*     IPSETC saves the value of an internal parameter.
*     IPGETC returns the value of an internal parameter.
*     IPFINI writes an error message for any parameter that has been
*     set but not read.
*
*     IPGETI is an interface to ipgetc; it retrieves an integer value
*     ival from the table of internal parameters, or returns ival
*     unchanged if the name is not found or if the corresponding entry
*     is a blank string. The internal parameter may have been set by an
*     earlier call to ipseti or to ipsetc.
*
*     IPGETR is an interface to ipgetc; it retrieves a real*8 value
*     rval from the table of internal parameters, or returns rval
*     unchanged if the name is not found or if the corresponding entry
*     is a blank string. The internal parameter may have been set by an
*     earlier call to ipsetr or to ipsetc.
*
*     IPSETI is an interface to ipsetc; it saves an integer value in
*     the table of internal parameters.
*
*     IPSETR is an interface to ipsetc; it saves a real*8 value in the
*     table of internal parameters.
*
*     LWIMAI is a group of routines (entry points in lwimai) that
*     maintain a base address for free space in the workspace area
*     b2mod_stack%wi. The entry points are lwiset, lwiget, and lwimax.
*     LWISET frees workspace by resetting the base address.
*     LWIGET reserves workspace by incrementing the base address.
*     LWIMAX returns the maximum workspace use.
*
*     LWMAIN is a group of routines (entry points in lwmain) that
*     maintain a base address for free space in the workspace area
*     b2mod_stack%w. The entry points are lwset, lwget, and lwmax.
*     LWSET frees workspace by resetting the base address.
*     LWGET reserves workspace by incrementing the base address.
*     LWMAX returns the maximum workspace use.
*
*     XERGET gets the unit number for error messages.
*
*     XERSET sets the unit number for error messages.
*
*     XERRAB writes an error message and the call chain, and terminates
*     execution of the program.
*
*     XERTST tests a condition; if violated, calls xerrab to terminate
*     execution.
*
*   ..2b File management and un*formatted input/output
*
*     CFENDF tests whether an un*formatted file is positioned at
*     end-of-file.
*
*     CFOPEN opens an external file. The file may be formatted or
*     unformatted.
*!!   This routine may contain system dependencies.
*
*     CFRUCH reads a character array from an un*formatted file.
*
*     CFRUIN reads an integer array from an un*formatted file.
*
*     CFRURE reads a real*8 array from an un*formatted file.
*
*     CFWUCH writes a character array to an un*formatted file.
*
*     CFWUIN writes an integer array to an un*formatted file.
*
*     CFWURE writes a real*8 array to an un*formatted file.
*
*   ..2c Low-level auxiliary routines
*
*     BLAS and BLAS-like routines (basic linear algebra subprograms):
*       IFILL (BLAS-like) fills an integer vector with a constant value.
*       SMIN (BLAS-like) computes the minimum value in a vector.
*       SMAX (BLAS-like) computes the maximum value in a vector.
*       SAMAX (BLAS-like) computes the maximum absolute value in
*       a vector.
*       B2SSUM (BLAS-like) sums the elements of a vector.
*       B2SASUM (BLAS) sums the absolute values of elements of a vector.
*       SFILL (BLAS-like) fills a vector with a constant value.
*       SBIAS (BLAS-like) adds a constant value to a vector.
*       B2SSCAL (BLAS) scales a vector by a constant.
*       B2SCOPY (BLAS) copies a vector to another.
*       B2SAXPY (BLAS) adds a constant-times-vector to a vector.
*
*     ESSEQ determines whether two real*8 numbers differ by no more than
*     a small multiple of their machine precision.
*
*     EVAL returns its real*8 argument through the function name.  It may
*     be used to force a compiler to store an intermediate operand.
*
*     RATIO evaluates very carefully the real*8 ratio of two integers.
*
*     SCAL10 is meant to compute a suitable scale for graphical output.
*
*     CHCASE converts a character variable to either all uppercase or
*     all lowercase.
*
*     STRCOPY copies the first l0 characters of string str1 into str0.
*
*     STREQL checks for equality between two strings, ignoring case
*     shift and trailing blanks.
*
*     STRMAS checks for a certain matching between two strings.
*
*   ..2d Auxiliary routines for discretised PDEs
*
*     ILU5F computes the ILU factorisation of a 2d 5-point matrix.
*
*     ILU5G computes the ILU factorisation of a 2d 5-point matrix and
*     applies the factors for one iteration of an update equation.
*
*     ILU5S solves an equation of the form L*U*x=r, where L and U are
*     the factors in an ILU factorisation of a 2d 5-point matrix.
*
*     ILUTER performs a number of iterations of the Stone strongly
*     implicit procedure for a system of equations of five-point form.
*
*     SIP5F computes the SIP-ILU factorisation of a 2d 5-point matrix.
*
*     SIP5G computes the SIP-ILU factorisation of a 2d 5-point matrix,
*     and applies the factors for one iteration of an update equation.
*
*     SLV5PT solves a five-point system of equations by a direct
*     method.
*
*     STABEQ adjusts the diagonal of a five-point matrix to satisfy
*     a property related to diagonal dominance.
*
*   ..2e System-dependent routines
*
*     DAYTIM returns the date and time in a character*17 format.
*!!   This routine may contain system dependencies.
*
*     JOBNAM returns an identification of the present job.
*!!   This routine may contain system dependencies.
*
*     MACHEPS returns the machine precision.
*!!   This routine may contain system dependencies.
*
*     MACHSFR returns the safe numerical range (a tiny positive machine
*     value e0 for which -e0, 1/e0, -1/e0 are all computable).
*!!   This routine may contain system dependencies.
*
*     PRVRT sets all elements of a real*8 array to a value that is
*     expected to cause an error when used in an arithmetic operation.
*     This is used to protect against certain programming errors.
*!!   This routine may contain system dependencies.
*
*     PRVRTI sets all elements of an integer array to a value that is
*     expected to cause an error when used in an arithmetic operation.
*     This is used to protect against certain programming errors.
*!!   This routine may contain system dependencies.
*
*     SYSINI performs any system dependent calls that may be required
*     at the start of program execution.
*!!   This routine may contain system dependencies.
*
*     SYSEND performs any system dependent calls that may be required
*     at the finish of program execution.
*!!   This routine may contain system dependencies.
*
*     USRNAM returns the name or other appropriate identification of
*     the submitter of the present job.
*!!   This routine may contain system dependencies.
*
*   ..2f Basic graphics routines
*
*     GNOPEN opens a graphics file for use with GKS.
*
*     GWSREP defines the workstation representation for GKS. It should
*     normally be called just once for each workstation, immediately
*     after opening the workstation.
*
*     PAGDEF sets GKS and NCAR properties to default values. It should
*     normally be called at the start of each new plot.
*
*     PAGINI starts a new plot.
*
*     PAGEND terminates a plot.
*
*     CPMPXY is a substitute for a routine of the same name in the
*     CONPACK sub-package of NCARG graphics. It defines a mapping from
*     rectangular array index coordinates in (1:*,1:*) to the plotting
*     window.
*
*     GNMAP2 (together with the alternate entry GNDFM2) is a service
*     routine for use with the CONPACK routines in the NCARG graphics
*     package. It defines a rectangular map for CONPACK. The map is
*     specified by a call to the entry gndfm2.
*
*     GNMAP3 (together with the alternate entry GNDFM3) is a service
*     routine for use with the CONPACK routines in the NCARG graphics
*     package. It defines a curvilinear, logically rectangular map for
*     CONPACK. The map is specified by a call to the entry gndfm3.
*
*     GNCOLR is for use by the routine ARSCAM in the AREAS utility of
*     NCAR graphics. It defines the filling of a polygon.
*
*     LBFILL is a substitute for a routine of the same name in the
*     LABELBAR utility of NCAR graphics. It fills a specified area.
*     LBFILX and LBFILQ are alternate entries in LBFILL. LBFILX is
*     called to define the filling style and LBFILQ is called to query
*     the filling style.
*
*   ..2g Atomic physics and surface physics
*
*     RATADP uses the ADPAK code to prepare tables of rate coefficients
*     for ionisation, recombination, electron energy loss, and charge
*     exchange; as function of temperature and density.
*
*     RATSTR uses the STRAHL code to prepare tables of rate coefficients
*     for ionisation, recombination, electron energy loss, and charge
*     exchange; as function of temperature and density.
*
*     RATPRT prints an excerpt from a table of atomic rate
*     coefficients.
*
*     RATTAI computes a smoothed set of rate coefficients to represent
*     the effect of a power-law tail on the distribution function.
*
*     RATWEI uses WEISHEIT data to improve tables of rate coefficients
*     for ionisation, recombination, electron energy loss, and charge
*     exchange; as function of temperature and density.
*     The WEISHEIT dataset is incomplete.  The call to ratwei must
*     be preceded by a call to one of the comprehensive routines.
*
*-----------------------------------------------------------------------
*.end b2cdcr

      end subroutine b2cdcr
