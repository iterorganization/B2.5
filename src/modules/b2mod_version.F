      module b2mod_version
      character*10, save :: newversion ='03.000.007'

      character*10, save ::
     1 b2fstati_version = '00.000.000',
     2 b2fstate_version = '00.000.000',
     2 b2frates_version = '00.000.000',
     3 b2fmovie_version = '00.000.000',
     3 b2ftrack_version = '00.000.000',
     4 plasmastate_version = '00.000.000',
     4 b2fplasma_version = '00.000.000',
     5 b2fgmtry_version = '00.000.000',
     6 b2fparam_version = '00.000.000',
     7 boundary_version = '00.000.000',
     8 neutrals_version = '00.000.000',
     9 b2faveri_version = '00.000.000'

      contains

      subroutine cfverr(nget, readversion)
      implicit none
      integer, intent(in) :: nget
      character*10 version_in, readversion
      character*12 chf
      character*7 label
      logical streql
      external streql

      inquire (nget,form=chf,err=10)
      if (streql(chf,'formatted')) then
        read(nget,'(a,a)',err=10) label,version_in
      else
        read(nget,err=10) label,version_in
      endif
      if(label.eq.'VERSION') then
        readversion=version_in
        write(*,'(a,a)') 'VERSION set to: ',readversion
        return
      endif
   10 continue
      rewind(nget)
      write(*,'(a,a)') 'VERSION defaulted to : ',newversion
      readversion=newversion
      return
      end subroutine cfverr

      subroutine cfverw(nput, writeversion)
      implicit none
      integer, intent(in) :: nput
      character*10, intent(in) :: writeversion
      character*12 chf
      character*32 get_B25_hash
      logical streql
      external streql, get_B25_hash

      inquire (nput,form=chf,err=20)
      if (streql(chf,'formatted')) then
        write(nput,'(a,a,1x,a)') 'VERSION',
     .                          writeversion,get_B25_hash()
      else
        write(nput) 'VERSION',writeversion,' ',get_B25_hash()
      endif
      write(*,'(4a)') 'VERSION: ',writeversion,
     .                ' GIT hash: ',get_B25_hash()
   20 continue
      return
      end subroutine cfverw

      end module b2mod_version
cxpb Version 1.1 corresponds to adding a version number and the UADIA array
cxpb to the 'b2fstat*' files. Implemented 20.01.2000.
cxpb Version 1.1 also adds the 'periodic_bc' variable to the geometry file.
cxpb Finalised 14.07.2000.
cdpc Version 01.001.001 corresponds to adding b2stbr_sna,smo,she,shi,sch to
cdpc b2mod_sources
cdpc Version 01.001.002 corresponds to adding the anomalous transport
cdpc (D,chi...) quantities to b2fplasma [2001.07.18]
cdpc Version 01.001.003 corresponds to adding the heat sources
cdpc from b2sihs to b2fplasma [2002.08.07]
cxpb Version 01.001.004 corresponds to adding the parallel momentum sources
cxpb from atomic processes and the total energy flows to b2fplasma [2002.08.13]
cdpc Version 01.001.005: added saving of fch_32 and fch_52 to b2fstat[ie]
cdpc Version 01.001.006: added saving of b2npmo_sma[fgv] to b2fplasma
cdpc Version 01.001.007: added saving of csigin & cvsa_cl to b2fplasma
cdpc Version 01.001.008: added saving of fmo to b2fplasma
cxpb Version 01.001.009: added saving of kinrgy to b2fstat[ie]
cxpb Version 01.001.010: added saving of wadia and vaecrb to b2fplasma
cdpc Version 01.001.011: added saving of b2stbc_s?? to b2fplasma
cdpc Version 01.001.012: added saving of facdrift and fac_ExB to b2fplasma
cdpc Version 01.001.013: added saving of ext sources
cdpc Version 01.001.014: added saving of additional current arrays
cdpc Version 01.001.015: added saving of fac_vis
cxpb Version 01.001.016: added rtlrd line radiation rate array to b2frates
cdpc Version 01.001.017: added neutral flux limiter diagnostic arrays
cxpb Version 01.001.018: added simulation time to b2fstate (per ANK request)
cxpb Version 01.001.019: added rtlbr bremsstrahlung rate array to b2frates
cxpb Version 01.001.020: added viscosity flux limiter diagnostic array
cxpb Version 01.001.021: added electron and ion heat flux limit arrays
cxpb Version 01.001.022: added qc array to b2fgmtry and b2fplasma
cxp                      added fch_p parallel current array to b2fstat[ie]
cxpb Version 01.001.023: added sclrtio, sclstra arrays to b2fplasma
cxpb Version 01.001.024: added rtlza, rtlz2, rtlpt and rtlpi arrays to b2frates
cxpb                     added saving of ne to b2fstat[ie]
cxpb                     new species data format: zamin, zamax, zn, am
cxpb                     added version number to b2fparam and b2fpardf
cxpb                     added nstra to b2fparam
cxpb                     made ppout a permanent fixture of b2frates
cxpb Version 01.001.025: added pbs computed from basis mesh contributions to b2fgmtry
cxpb Version 01.001.026: added sne, b2stbr_sne, b2stbc_sne and ext_sne to b2fplasma
cxpb Version 01.001.027: added saving of b2stbm_s?? to b2fplasma
cxpb Version 01.001.028: new style for boundary and neutrals namelist to allow for
cxpb                     full vacuum vessel grid coverage
cxpb Version 01.001.029: split the thermo-electric coefficient into a classical and
cxpb                     and an anomalous piece
cdpc Version 01.001.030: added saving vsal to b2fplasma
cdpc Version 01.001.031: added saving smpr, smth, smfr to b2fplasma
cxpb

csrv                     B2solps5.2 starts 2008
csrv Version 01.001.025  fch_p(:,:) replaced by fch_p(:,:,0:1) and many otheres
c                        changes have been made. Sorry for absence of details.
csrv Version 01.001.026  from 28.02.2009.
csrv Version 01.001.027  09.03.09 inserted area_fix=3 in b2news_
c                                 corrected b2stbc_bas and b2stbc_phys
c                        11.03.09 corrected boundary conditions bceni=15 in b2stbc_phys
c    Version 01.001.028  05.05.09 corrected boundary conditions bcene=15 in b2stbc_phys
c    Version 01.001.029  07.05.09 modified boundary conditions bcmom=13 in b2stbc_phys
c                                 and signs of coefficients in b2.transport.parameters
c                                 for bcmom=13
c    Version 01.001.030  08.05.09 corrected boundary conditions bcene=3 in b2stbc_phys
c    Version 01.001.031  11.05.09 introduced subroutine b2tinnt for calculating
c                                 ion-neutrals current
c    Version 01.001.032  15.05.09 exchanged subroutine set_transport_srv for case(5)
c    Version 01.001.033  22.05.09 introduced case(9) instead of case(5) with new barrier
c                                 function for h-regime
c    Version 01.001.034  10.06.09 exchanged definition of sig0(ix,iy) in b2tqna
c    Version 01.001.035  16.06.09 corrected definition of alfe0(ix,iy,0) in b2tfhe_
c                                 deleted SPb style for csigin in b2trcl.F
c                                 added factor hz for b2sigp_style=0
c                                 exchanged the interface of b2tinnt.F
c    Version 01.001.036  01.07.09 introduced b2tlnl.F for the Coulomb logarithm
c    Version 01.001.037  07.07.09 added the heat terms into b2stel.F, b2stcx.F, b2sihs_.F
c                                 changed interface of mapx.F and mapxy.F
c    Version 01.001.038  08.07.09 re-wrote b2tinnt.F and changed b2tqin.F
c    Version 01.001.039  15.07.09 b2agmx.F replaced by b2agmx.F from B2solps5.0
c                                 modified subroutine open_file.F
c    Version 01.001.040  11.09.09 done new diagnostics output, introduced new arrays for
c                                 storing physical sources in r.h.s.
c                        18.09.09 cleared guard cells for shei in b2nph9.F and b2npht.F
c    Version 01.001.041  22.09.09 changed the calculation of uesx in b2sihs_.F
c    Version 01.001.042  25.09.09 added bccon=16 in b2stbc_phys.F
c    Version 01.001.043  07.10.09 deleted bccon=15 in b2stbc_phys.F. Use bccon=5.
c    Version 01.001.044  14.10.09 changed status of the following files b2mn.prt, b2fparam,
c                                 b2fstate, b2fmovie, b2ftrace, b2ftrack.
c    Version 01.001.045  20.10.09 changed logical condition in b2tral
c                                 added term in flo_nodrift and flo2dia in b2tfnb
c                        21.10.09 introduced internal key 'b2tfhe_mode_ehx'
c                        23.10.09 modified energy recycling terms in b2stbr_phys
c    Version 01.001.046  29.10.09 changed contents of b2fstati_add and b2fstate_add
c    Version 01.001.047  05.11.09 added termhz to flo2dia in b2tfnb
c                        12.11.09 corrected b2mod_input_profile
c                        17.11.09 corrected b2gxdr_dummy in b2agfs.F
c    Version 01.001.048  23.11.09 updated module b2mod_constants and subroutines rat*.F
c                        03.12.09 inserted internal key 'b2trcl_lluciani'
c    Version 01.001.049  22.12.09 in order to copy the data of first n=min(ns0,ns1) species
c                                 from the plasma state file b2fstate for ns0 species
c                                 into the initial plasma state file b2fstati for ns1 species
c                                 the special internal key 'b2aidr_read_b2fstate'
c                                 has been added in b2aidr.F
c                                 Grid size, nuclear and atom charge must be the same
c                                 in both files.
c    Version 01.001.050  30.12.09 modified term in the electron heat equation b2tfhe_.F
c                                 deleted statement use b2mod_ma28 in b2co.F and b2fu.F
c    Version 01.001.051  30.01.10 changed value of internal key b2tfnb_xfrhiehz,
c                                 b2nph9_style_shei by default. The name of key b2tfhe_no_alfTeEh
c                                 replaced by the name b2tfhe_alfTeEh
c                                 corrected the long line in b2tfnb which produced small errors
c                                 for working in B2solps5.0 style
c    Version 01.001.052  15.02.10 corrected the name of key in ipgeti in b2npht.F
c                                 uncommented the lines in b2npc9.F for guard_flows.eq.2
c                                 introduced key 'b2news_re_compute_pb'
c    Version 01.001.053  24.02.10 added output fmo in b2wdat.F
c                                 added the diagnostic output in b2stbc_phys.F
c                        01.03.10 added in b2mndr.F the output of fort.30 and fort.31 for Eirene
c                                 corrected the diagnostic output of fhe and fhi
c    Version 01.001.054  12.04.10 corrected the calculation of thermal force in b2sifr_.F
c    Version 01.001.055  19.04.10 added the calculation of anomalous velocity for H-regime in b2tqna.F
c                                 added additional diagnostic into b2stbr.F for B2-Eirene case
c                                 corrected output of dab2,dmb2 and tab2,tmb2 in b2wdat.F
c    Version 01.001.056  06.07.10 modified the condition of solving of equations for Ua in b2npmo.F
c                        12.07.10 corrected the calculation of transport parameters for H-regime
c                                 for UDND and LDND cases in b2tqna.F
c                        14.07.10 added b2tfrn.F for evaluating of special fluxes for B2-Eirene coupling
c                                 changed copy_background.F for passing new special fluxes to Eirene
c                        16.07.10 corrected diagnostic output for b2srdt.F
c                        23.07.10 Corrected diagnostic calculation of total energy, total current,
c                                 total momentum and power.
c                                 Changed numerical realization of bccon=16.
c                                 Implemented the modified decay length boundary condition
c                                 for the electron temperature bcene=22.
c                                 Implemented the modified decay length boundary condition
c                                 for the ion temperature bceni=22.
c    Version 01.001.057  28.07.10 Corrected my_intcor.F
c                                 Corrected b2tqna.F by adding of attribute save for variable mintopcut
c                                 Corrected default value of eirene_fixmap in the eirene_f30f31.F
c    Version 01.001.058  21.10.10 Released pressure correction procedure for cases without and with drifts
c                                 with modificated particle anf heating fluxes.
c                        27.10.10 commented parts which use files b2fstati_add and b2fstate_add
c                                 commented key 'b2news_re_compute_pb'
c    Version 01.001.059  18.11.10 replaced b2ardr.F by new from xpb
c                        02.12.10 introduced Zhdanov expression for fke(z) in b2tqce.F
c    Version 01.001.060  20.01.11 inserted factor 5/3 before cvla into fna_he
c    Version 01.001.061  27.01.11 changed the applying of correction for boundary cells. See b2npmo.F
c    Version 01.001.062  02.02.11 corrected boundary condition BCENI==22 at 'N'
c    Version 01.001.063  15.02.11 updated some common subroutines which are used B2solps5.0 (release 3066
c                                 and B2solps5.2
c                                 removed the heat terms into b2stel.F, b2stcx.F, b2sihs_.F which were added in 07.07.09
c                                 corrected the coefficient before csigin(,,0,,) in b2tfhe_.F from x-component of ion-neutral current
c                                 added calculation ion-neutral current in b2tinnt.F for B2-Eirene case.
c    Version 01.001.064  05.03.11 changed calculation of ue in b2sihs_.F
c                                 changed condition of stopping B2mn using ether file .quit or _quit 13.04.11
c                                 is = 1 has been replaced by is = ismain   12.05.11
c    Version 01.001.065  24.05.11 changed case of lluciani=3 in b2trcl.F
c    Version 01.001.066  30.05.11 added currents in radial part of ion heat flux in b2tfhi_.F
c                                 corrected output in b2urmo.F and its interface too.
c                                 replaced default value area_fix=0 by area_fix=3
c    Version 01.001.067  02.06.11 introduced restriction of diamagnetic velocities of impurity ions in b2tfnb.F
c                        09.06.11 introduced new numerical manner of calculation fnbPSch(ix,iy,0) and fnbPSch(ix,iy,1)
c                        15.07.11 replaced 1.e-3 by 3.e-3 in Eirene sputer.F There was very strange "access violation" in this place
c                        22.07.11 corrected grlib.F (dummy) and masr3.F
c    Version 01.001.068  28.07.11 added coupling B2 and ASTRA. For coupling it is necessary to use pre-processor variable ASTRA
c    Version 01.001.069  09.09.11 changed output in b2fstate, b2fplasma and input from b2fstati fch_p. In both cases it is used array fch_p(:,:)
c                        20.09.11 modyfied b2tlnl.F
c                        21.09.11 introduced nwtf in b2plot.F instead of 97
c    Version 01.001.070  17.10.11 modyfied Braginsky formula in b2tlnl.F Introduced restriction lnlam >= 5
c                        03.11.11 changed calculation of chci_koeff in b2trcl.F
c                                 added the usage of -1 cell in b2mod_astra_to_b2.F
c                        15.11.11 changed b2mwti.F and b2mndr.F in according to request of Bonnin
c                        29.11.11 introduced BCCON=18 for B2-ASTRA coupling
c    Version 01.001.071  11.01.12 inserted call b2stbr_neutr_scl(nx, ny, ns) into b2news_.F
c                                 corrected the contribution from friction between charged species in b2sihs_.F
c                        17.01.12 added internal key 'b2trcl_conductive_limit' in b2tfhe_.F and b2tfhi_.F
c                                 corrected default value of lluciani in b2siav.F
c                        27.01.12 added diagnostic keys into b2tfrn.F
c                        07.02.12 corrected diagnostic output in b2sifr_.F
c    Version 01.001.072  15.02.12 reinstalled keys b2tfnb_no_Ptncr_x and b2tfnb_no_Ptncr_y
c                        15.02.12 corrected calculation of kinetic energy in b2tfnb.F and fluxes in b2tfrn.F
c                        16.02.12 replaced the number BCCON=16 by BCCON=15
c                        16.02.12 added another style of BCCON=16 feedback boundary condition for given integrated particle
c                                 flux with constant density
c                        16.02.12 added another style of BCENE=16 feedback boundary condition for given integrated electron
c                                 heat flux with constant electron temperature
c                        16.02.12 added another style of BCENI=16 feedback boundary condition for given integrated ion heat
c                                 flux with constant ion temperature
c                        20.02.12 added the contribution to the integrated source se from Eirene
c                        12.03.12 added check that we have some particle diffusion in b2tral.F
c                        14.03.12 added diagnostic output in subroutine write_b2_to_astra
c                        23.03.12 added diagnostic output in subroutine b2stbr.F
c                        26.03.12 introduced restriction for argument of exp in subroutine recycle from file b2stbr_phys.F
c                        28.03.12 added the procedure of smoothing of Eirene sources suggested by Senichenkov
c                        12.04.12 added additional diagnostics for particles of fluxes
c    Version 01.001.073  12.04.12 deleted the modified fluxes from recycling boundary conditions
c                        21.05.12 modified the algorithm of searching radial midpalne coordinates in module b2mod_b2_to_astra
c                        14.06.12 added setting constatnt electron and ion temperature from ASTRA
c                        15.06.12 replaced natm by nnatmi in b2stbr.F
c                        20.06.12 added restriction of usage BCCON=8 in case of coupling B2-ASTRA
c                        21.06.12 added NB sources in case of coupling B2-ASTRA
c                        16.07.12 added test of correct using boundary conditions BCENE=9, BCENI=9, BCMOM=10 along with BCCON=9
c                                 ehxp is set up zero subroutine b2nxfx when electric potential is zero (pot_eq=0)
c                        26.07.12 corrected contribution of ASTRA sources in b2sral
c                        21.09.12 test of correct using boundary conditions BCENE=9, BCENI=9, BCMOM=10 along with BCCON=9
c                                 has been replaced by warning message.
c                        26.09.12 introduced the freezing of na, te, ti when 'b2tfhe_mdf','b2tfhi_mdf' and 'b2tfnb_mdf' are '1'
c                        10.10.12 for compatibility with gfortran 'binary' has been replaced by 'unformatted' in open statement
c                        23.10.12 corrected contribution of PeNB_astra and PiNB_astra to she and shi
c                        01.11.12 corrected boundary condition bccon=19
c    Version 01.001.074  07.11.12 corrected boundary condition bcene=15
c    Version 01.001.075  22.11.12 corrected multiplication of fchvispar by factdrift in b2tfhe_.F
c                                 replaced lines of calculation fchinert in b2tfhe_.F by calling b2tiner.F
c                        20.12.12 redone calculating intergated fluxes for B2-ASTRA case.
c    Version 01.001.076  11.01.13 deleted factor signmf in drifts terms in b2tfnb.F
c                                 replaced bb(:,:,3) by bb(:,:,2) in b2sihs_.F and b2tinnt.F
c    Version 01.001.077  11.01.13 changes in subroutines have been done in order to change sign of Bx and/or Bz automatically
c                        27.02.13 eliminated influence of Ehx, Ehy and currents when pot_eq = 2
c                        11.03.13 changes in subroutines have been done in order to have compatibility with gfortran
c                        13.04.13 changes in subroutine b2wdat.F have been done in order to have compatibility with gfortran
c    Version 01.001.078  11.04.13 modified bcene=16 and bceni=16 for B2-ASTRA case
c                        18.04.13 corrected bcene=15 and made two variants of implementation.
c                        19.04.13 changed the applying of alfa3 to anomalous transport coefficients in H-regime and introduced alfa4
c                                 corrected calculation of qe at the plates and modified output for B2-ASTRA case
c                                 implemented option in b2tqce.F in order to set sigx=0 and alfx=0 in guard cells at north and south except of core region
c                                 implemented option in b2tiner.F in order to calculate inerthia current using different approximation
c                                 inserted boundary conditions BCCON = 21, 22, 23 suggested by L. Kaveeva and I. Senichenkov
c                                 inserted boundary conditions BCMOM = 16 suggested by L. Kaveeva and I. Senichenkov
c                                 inserted boundary conditions BCENI = 23 and 24 (with modified flux) suggested by L. Kaveeva and I. Senichenkov
c                                 changed subroutine intcor.F Complicated operator IF was splitted by two in order that debugger not complain about index
c                        11.06.13 added calculation of integrated electron and ion heat sources due to atomic processes for B2-ASTRA coupling
c
c    Versions 02.00x.0xx  Reserved for SOLPS6.0
c
c    Version 03.000.000  09.02.15 Inaugural SOLPS-ITER version
c                                 Added smcf centrifugal force array to b2fplasma
c    Version 03.000.001  01.06.15 Adding the POB array to fort.44 per Juelich code
c                        11.06.15 Making the default switches follow the 5.2 code
c                                 base physics as opposed to 5.0
c    Version 03.000.002  29.06.15 Changed default value to L_NEUTRAD to 0
c                                 Brings in neutral radiation from Eirene computation
c                                 Changed default value of sheath_drift_fix to 1
c    Version 03.000.003  13.07.15 Added modified and Pfirsch-Schlueter fluxes to b2fplasma
c                                 Added fna_fcor, fna_he, fchvisq, fchinert to b2fplasma
c    Version 03.000.004  16.10.15 Changed to default BCPOT=12 amd 13 for free boundaries with currents
c                                 Changed format and fitting formulae for bulk material properties files
c    Version 03.000.005  25.04.16 Added various quantities to b2fstati/e to avoid restart
c                                 effects when using various advanced features:
c                                 fna_mdf, fhe_mdf, fhi_mdf,
c                                 fna_fcor, fna_nodrift, fna_he,
c                                 fnaPSch, fhePSch, fhiPSch,
c                                 fna_eir, fne_eir, fhe_eir, fhi_eir,
c                                 fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
c                                 fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
c                                 vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc
c    Version 03.000.006  29.06.18 New default model b2sigp_style.eq.2
c                                 Added redef_pbs switch to b2fgmtry file
c                                 Added cflim(3:4) flux limit factors for
c                                 thermo-electric coefficient and friction terms
c                                 Added flux limit coefficients to b2mod_transport
c    Version 03.000.007  01.07.19 Added fix_recomb and zmax_recomb to b2frates
c                                 Replaced redef_pbs by redef_gmtry in b2fgmtry
c                                 Added wbbl to b2fgmtry
c                                 Added (corrected) fht, fhj, fhm, fhp to b2fplasma
c
