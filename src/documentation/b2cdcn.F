*-----------------------------------------------------------------------
*.specification

      subroutine b2cdcn ()
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2CDCN provides information about namelists of the b2 code.
*     No computation is performed in this routine.
*
*     The namelists are listed in no particular order. Each description
*     will begin with the file from which the code expects to read the 
*     namelist and its general purpose. Then, each namelist variable
*     will be listed with its use and size.
*
*-----------------------------------------------------------------------
*.text
*
* NAMELIST /NEUTRALS/
*  Found in b2.neutrals.parameters.
*    This namelist is used to either set parameters for the B2 neutral fluid 
*    model or for the Eirene Monte-Carlo neutrals run.
*    When coupling to Eirene, it is essential that the neutral boundaries
*    conditions declared here be in the same exact order than they are in
*    the Eirene input file.
*    For coupled cases, some of the parameters will become active only on the 
*    second call to Eirene, while values inherited from the Eirene input file 
*    will be used for the first call. In such cases, use of the switch 
*    'eirene_repeat_first_call' is recommended.
*    This namelist is read if either 'b2stbr_neutrals_namelist' or 
*    'b2mndr_eirene' are set to 1.
*    Some of the arrays are dimensioned by NTRACK, which is the maximum
*    number of tracked species in the mixed materials model. This parameter
*    is set in b2mod_layer and currently holds the value 4.
*  NSTRAI - integer. Number of neutral sources, or 'strata'. Must not be larger
*     than DEF_NSTRA from $(SOLPSTOP)/include(.local)/DIMENSIONS.F file.
*     Defaults to 0.
*     Need not include the time-dependent stratum. If a time-dependent stratum
*     is used, is incremented internally as needed. The incremented value is 
*     referred to as NSTRAT below.
*  RCPOS - integer array of length (NSTRAT). Position in the B2 grid of the 
*     of strata. Similar use as BCPOS from /BOUNDARY/ namelist. Defaults to -2.
*  RCSTART - integer array of length (NSTRAT). Start coordinate on the B2 grid 
*     the strata. Similar use as BCSTART from /BOUNDARY/ namelist. Defaults to -2.
*  RCEND - integer array of length (NSTRAT). End coordinate on the B2 grid of
*     the strata. Similar use as BCEND from /BOUNDARY/ namelist. Defaults to -2.
*  RC_LIST_SIZE - integer array of length (NSTRAT). Contains the size of the 
*     recycling boundary lists. Similar use as BC_LIST_SIZE from /BOUNDARY/ 
*     namelist. Defaults to 0.
*  RC_LIST_X - integer array of length (2*(NXD+NYD),NSTRAT). Contains the 
*     X-coordinate of the cells where recycling boundaries are applied. 
*     Similar use as BC_LIST_X from /BOUNDARY/ namelist. Defaults to -2.
*  RC_LIST_Y - integer array of length (2*(NXD+NYD),NSTRAT). Contains the 
*     Y-coordinate of the cells where recycling boundaries are applied. 
*     Similar use as BC_LIST_Y from /BOUNDARY/ namelist. Defaults to -2.
*  TARGSP - integer array of size (NSTRAT,NTRACK). Identifies the base 
*     material(s) of this stratum wall. The number corresponds to the B2 
*     species index produced by sputtering. Defaults to 'b2stbr_sput_dst'. 
*     Check b2cdci for details. If the bulk species is not sputtered,
*     should contain -1.
*  CHEMSP - logical array of length NSTRAT. Indicates whether chemical 
*     sputtering is allowed from this wall stratum. Defaults to .false..
*  RECYC - real*8 array of size (0:NS-1,NSTRAT). Stores the particle recycling 
*     coefficients of species (is) on stratum (istra). Defaults to 1.0 for 
*     non-carbon species and 0.0 for carbon species. Particles recycle into
*     the neutral species associated to their homonuclear sequence.
*     Applies to B2 neutral fluid species.
*     Also multiplies Eirene recycling fluxes if 'eirene_use_recyceir' is set 
*     to 0 (see b2cdci for details).
*  MRECYC - real*8 array of size (0:NS-1,NSTRAT). Stores the parallel momentum
*     recycling coefficients of species (is) on stratum (istra). Defaults to 0.0. 
*     The parallel momentum is carried back by the recycled neutrals.
*     Applies only to B2 neutral fluid species.
*  ERECYC - real*8 array of size (0:NS-1,NSTRAT). Stores the energy recycling 
*     coefficients of species (is) on stratum (istra). Defaults to 0.3 for 
*     non-carbon species and 0.0 for carbon species. The energy is carried back
*     by the recycled neutrals.
*     Applies only to B2 neutral fluid species.
*  RCION - real*8 array of size (0:NS-1,NSTRAT). Stores the particle-into-ion 
*     recycling coefficients of species (is) on stratum (istra). Defaults to 
*     0.0. Particles recycle into the next ionised species associated to their 
*     homonuclear sequence. Applies only to B2 neutral fluid species.
*  RECYCEIR - real*8 array of size (NSTRAT). Multiplier to the Eirene recycling
*     fluxes from stratum (istra). Only used if 'eirene_use_recyceir' is set 
*     to 1 (default). Defaults to 1.0.
*  USERFLUXPARM - real*8 array of size (NSTRAT,2). The element (istra,1) 
*     contains the strength of constant gas puff strata (type 'C') for Eirene 
*     in particles/second. Defaults to 0.
*  CRCSTRA - character*1 array of length (NSTRAT). Contains the type of stratum
*     for Eirene. Defaults to ' '. Possible options include:
*     'N','S','W','E' - topological mesh boundaries (same use as BCCHAR in 
*           /BOUNDARY/ namelist)
*     'V' - volume recombination source (related RCPOS, RCSTART, RCEND and
*           RC_LIST variables are ignored).
*     'C' - constant or feedback gas puff source (related RCPOS, RCSTART,
*           RCEND and RC_LIST variables are ignored). Gas puffs for B2 fluid 
*           neutrals are handled via b2.boundary.parameters.
*           Unless specified otherwise by use of 'eirene_nesepm_istra', it is
*           the first 'C' stratum that is used for feedback puff schemes.
*           See b2cdci for details.
*     'T' - time-dependent source for Eirene (related RCPOS, RCSTART, RCEND 
*           and RC_LIST variables are ignored). Long-lived Eirene neutrals  
*           are stored in this stratum. See EIRENE_STEP_DT below.
*  RF_NEUT - real*8 array of size (4). Relaxation parameters multiplying the 
*           Eirene particle, parallel momentum, electron energy and ion energy
*           sources respectively before use in B2. Defaults to 1.0.
*  PHYS_SPUT - real*8 array of size (0:NS-1,NSTRAT). Stores the physical 
*     sputtering multiplier for B2 species (is) sputtering on the wall of 
*     stratum (istra). Must be non-zero for physical sputtering to occur.
*     Only applies if using B2 fluid neutral model. Defaults to 0.0.
*  CHEM_SPUT - real*8 array of size (0:NS-1,NSTRAT). Stores the chemical 
*     sputtering multiplier for B2 species (is) sputtering on the wall of 
*     stratum (istra). Must be non-zero for chemical sputtering to occur.
*     Only applies if using B2 fluid neutral model. Defaults to 0.0.
*  EIRENE_STEP_CPU - real*8. Length of CPU time devoted to Eirene calls after 
*     the first one (in s). Defaults to the value given in input.dat.
*  EIRENE_STEP_DT - real*8. Physical time (in s) the Eirene particles are to 
*     be followed before being passed to the time-dependent stratum.
*     Defaults to 1.0e-3 seconds.
*  EIRENE_MOD - integer. Frequency of Eirene calls. Eirene is called every 
*     EIRENE_MOD full B2 iterations. Defaults to 1.
*  VOLRECSTART - real*8 array of size (NSTRAT). Initial volume recombination 
*     rate from stratum (istra) to be used for the volume recombination 
*     scaling regulation algorithm. Rendered obsolete by 'eirene_dpc_fix'. 
*     Defaults to 1.e21 particles/second. 
*  VOLRECINC - real*8. Maximum rate of increase of the volume recombination 
*     source strength. The volume recombination regulation scaling scheme is 
*     activated if VOLRECINC is neither 0 nor 1. 
*     Rendered obsolete by 'eirene_dpc_fix'.
*  VOLRECWT - real*8. Weight of new volume recombination strength relative to 
*     that of previous iteration in the volume recombination regulation 
*     scaling scheme. Defaults to 0.1.
*     Rendered obsolete by 'eirene_dpc_fix'.
*  SPECIES_START - integer array of size (NSTRAT). Specifies the index of the 
*     first B2 species involved in stratum (istra). Defaults to 0.
*  SPECIES_END - integer array of size (NSTRAT). Specifies the index of the 
*     last B2 species involved in stratum (istra). Defaults to ns-1.
*  NEUTRALS_FILENAME - character*256. Name of the next file to use for reading 
*     a new /NEUTRALS/ namelist. Defaults to 'b2.neutrals.parameters'.
*  NEUTRALS_TIME_MOD - real*8. When the B2 run simulation time, in seconds, 
*     modulo(NEUTRALS_TIME_MOD), exceeds NEUTRALS_TIME_SWITCH, reads the new 
*     namelist from NEUTRALS_FILENAME. Also switches to the new namelist if 
*     the ELM count (here time/NEUTRALS_TIME_MOD) changes. Defaults to 0.0.
*     Only active if NEUTRALS_TIME_MOD is greater than 0.
*  NEUTRALS_TIME_SWITCH - real*8. Time (in seconds) within an ELM cycle after 
*     which a new NEUTRALS namelist is read. Defaults to 0.0.
*  L_NEUTRAD - integer. If l_neutrad >= 0, then the radiated power due to the
*     neutrals atoms is taken directly from the Eirene calculation, instead
*     of being recomputed by B2. Defaults to 0.
*  L_NEUTFLUX - integer. If l_neutflux >=0, then correct treatment of the
*     incident fluxes in B2 and b2plot; if <0, then old (approximate) treatment
*  LSTRASCL - integer array of size (NSTRAT,0:natm). Indicates with which 
*    Eirene atomic species to scale the Eirene stratum (istra) 
*    (0 means no scaling, default). Atomic species 0 stands for electrons.
*  B2ESPCR - integer array of size (0:NS-1). Contains the index of the Eirene 
*    atomic species corresponding to the B2 species (is). Defaults to the
*    ordering of the B2 isonuclear sequences.
*  EB2SPCR - integer array of size (NATM). Contains the index of the B2 
*    neutral fluid species corresponding to the Eirene atomic species (iatm).
*    Defaults to the first B2 species of each isonuclear sequence.
*  LMOLSCL - integer array of size (NMOL). Contains the index of the Eirene 
*    atomic species with which the Eirene molecular species (IMOL) should be 
*    scaled. Defaults to 0.
*  MLCMP - integer array of size (NATM,NMOL). Contains the number of atoms 
*    from Eirene atomic species IATM within each molecule of Eirene molecular 
*    species IMOL. Defaults to 0.
*  LIONSCL - integer array of size (NION). Contains the index of the Eirene 
*    atomic species with which the Eirene test ion species (IION) should be 
*    scaled. Defaults to 0.
*  LCNS - integer array of size (NSTS). Contains the indices of Eirene 
*    non-standard surfaces (block 3A in input.dat) corresponding to core 
*    boundaries. Defaults to 0.
*  LTNS - integer array of size (NSTS). Contains the indices of Eirene 
*    non-standard surfaces (block 3A in input.dat) corresponding to target 
*    boundaries. Defaults to 0.
*  LSNS - integer array of size (NSTRAT). Contains the indices of Eirene 
*    surfaces related to the recycling sources. Defaults to 0.
*  KSNS - integer array of size (NSTRAT). Contains the number of parts
*    for each Eirene recycling stratum. Defaults to 0. 
*  GPFC - real*8 array of size (NATM,NSTRAT). Specifies the fraction of Eirene 
*    atomic species (iatm) in one particle puffed from stratum (istra).
*    Defaults to 0.0.
*  DBG_EIR_MC - integer. Debug output control for eirene_mc routine. See code 
*    for usage. Defaults to 0.
*  DEBUG_FLAGS - integer array of size (100). Non-zero elements yield 
*    additional print-out data for debugging B2-Eirene coupling. See code for 
*    specific uses. Many slots still available for user-specific needs.
*    Defaults to 0.
*  NEUT_SCL_LIM - real*8. Maximum number by which Eirene neutral sources may 
*    be scaled in either direction within the ank-mods scheme. See explaining
*    text in $SOLPSTOP/doc/Source_Scaling_in_B2.pdf for details. 
*    Defaults to 2.0.
*  TRACK_INDEX - integer array of size (0:NS-1). Specifies the mixed material 
*    species index related to B2 species (is). Default is 1 for species 
*    sput_dst, and 0 for others.
*  TRACK_CHEM_SPUT - logical array of size (NTRACK). Indicates whether mixed 
*    material species (itrack) participates in chemical sputtering.
*    Defaults to .true. for first tracked species if sput_dst is carbon,
*    and to .false. for all other cases. 
*  CHEMICAL_EROSION_REDEP_FAC - real*8. Multiplier to the chemical sputtering 
*    coefficient of carbon for computing the rate used for redeposited carbon. 
*    Defaults to 1.0.
*  CHEMICAL_EROSION_BE_FAC - logical. Indicates whether the presence of 
*    beryllium is to impact on the chemical sputtering yield of carbon. 
*    Defaults to .false.
*  CHEMICAL_EROSION_BE_FAC_A - real*8. Defaults to 0.2. The chemical 
*    sputtering yield of carbon is multiplied by 
*    (1.0-C/2*(tanh((frac-A)/B)-tanh((-A)/B))) 
*    where frac is the fractional content of Be in the surface layer material.
*  CHEMICAL_EROSION_BE_FAC_B - real*8. Defaults to 0.05. see above.
*  CHEMICAL_EROSION_BE_FAC_C - real*8. Defaults to 0.9. See above.
*  N_SPCSRF - integer. Number of special surfaces groups for diagnostics. 
*    Defaults to 0.
*  L_SPCSRF - integer array of length (NLIM+NSTS). List of surface segments 
*    (non-default standard surfaces [NDSS] or additional surfaces in Eirene 
*    notation) included in the groups. Negative numbers correspond to NDSS.
*    Defaults to 0.
*  SPS_ID - character*8 array of length (N_SPCSRF). Labels of groups of 
*    special surfaces.
*  I_SPCSRF - integer array of length (N_SPCSRF). Index of the first Eirene 
*    surface belonging to a special surface group in the L_SPCSRF list.
*    Defaults to 0.
*  J_SPCSRF - integer array of length (N_SPCSRF). Index of the last Eirene 
*    surface belonging to a special surface group in the L_SPCSRF list.
*    Defaults to 0.
*  SPS_ABSR - real*8 array of length (N_SPCSRF). Particle absorption on the 
*    surface. If positive, will supercede the setting from the Eirene input
*    file: RECYCT = 1-SPS_ABSR. Ignored if negative. Defaults to -1.0.
*  SPS_TRNO - real*8 array of length (N_SPCSRF). Surface transparency in 
*    positive direction. If positive, will supercede the setting from the
*    Eirene input file: TRANSP(1,) = SPS_TRNO. Ignored if negative.
*    Defaults to -1.0.
*  SPS_TRNI - real*8 array of length (N_SPCSRF). Surface transparency in 
*    negative direction. If positive, will supercede the setting from the
*    Eirene input file: TRANSP(2,) = SPS_TRNI. If negative, the setting
*    from SPS_TRNO is used. Defaults to -1.0.
*  SPS_MTRI - real*8 array of length (N_SPCSRF). Surface material in Eirene 
*    notation (e.g., 1206 for C). If positive, will supercede the setting
*    from the Eirene input file: ZNML = SPS_MTRI. Ignored if negative.
*    Defaults to 0.
*  SPS_MTRL - character*8 array of length (N_SPCSRF). Surface material in 
*    human notation (e.g., 'C').
*  SPS_TMPR - real*8 array of length (N_SPCSRF). Surface temperature (in eV).
*    If set larger to 1.e10, will be ignored. Otherwise, will supercede the
*    setting from the Eirene input file: EWALL = SPS_TMPR. Defaults to 1.e15.
*  SPS_SPPH - real*8 array of length (N_SPCSRF). Fudge factor for physical 
*    sputtering. If positive, will supercede the setting from the Eirene
*    input file: RECYCS = SPS_SPPH. Ignored if negative. Defaults to -1.0.
*  SPS_SPCH - real*8 array of length (N_SPCSRF). Fudge factor for chemical 
*    sputtering. If positive, will supercede the setting from the Eirene
*    input file: RECYCC = SPS_SPCH. Ignored if negative. Defaults to -1.0.
*  SPS_SGRP - integer array of length (N_SPCSRF). Sputtered particle species
*    flag for chemical sputtering. If positive, will supercede the setting
*    from the Eirene input file: ISRC = SPS_SGRP. Ignored if negative.
*    Defaults to -1.
*  WRITE_NML_NEUT - logical. If .true. (default), writes the content of the
*    namelist to stdout after it has been read.
*  TIME_DEP_PUFF - logical array of length (NSTRAT). Indicates whether a 
*    stratum is a time-dependent gas puff. Should only be .true. for strata
*    defined as gas puffs. Defaults to .false. for all strata.
*  NGPDATA - integer data of size (NSTRAT). Number of points over which the
*    time profile of the strength of gas puff (istra) is given. Defaults to 0.
*    Current maximum value is 100.
*  GPDATA - real*8 data of size (NGPDATA,2,NSTRAT). 
*     For (i,ik,istra) in (1:NGPDATA(istra),1:2,1:NSTRAT), 
*     GPDATA(i,1,istra) contains the time point (i) for stratum (istra) (in s).
*     GPDATA(i,2,istra) contains the gas puff strength of stratum (istra) at
*       time point (i) (in particles/s).
*     The gas puff strength before the first time point is given by 
*       USERFLUXPARM(1,istra).
*     The gas puff strength after the last time point is given by the
*       GPDATA value of the last time point.
*     Otherwise, the gas puff strength is linearly interpolated from the
*       given data.
*     Defaults to 0.0.
*  
* The following are active only if NO_B2_CHEM_SPUT is not defined among the 
* compiler options and are used to override the Eirene chemical sputtering
* settings:
*
*  CHEMICAL_SPUTTER_YIELD - real*8 array of size (0:NLIM+NSTS). 
*    Passed to Eirene. Chemical sputter yield of wall surface (ilim).
*    Defaults to 0.0.
*  FCHAR_CHEMICAL - real*8. Passed to Eirene. Nuclear charge of atomic species 
*    causing the sputtering. Defaults to 0 [no chemical sputtering].
*  IGASS_CHEMICAL - integer. Passed to Eirene. Eirene atomic species index of 
*    the sputtered particle. Defaults to 0. If igass_chemical > natmi, the
*    data from chemical_sputter_yield is not used and the yield from the Eirene
*    surface blocks is used instead.
*  ITSPUT_CHEMICAL - integer. Passed to Eirene. Eirene type index of the 
*    sputtered particle. Atoms = 1, Ions = 4. Defaults to 0, meaning
*    1 eV atom chemical sputtering.
*  ISSPUT_CHEMICAL - integer. Passed to Eirene. Mass*100 + nuclear charge of
*    the sputtered particle. Defaults to 0, internally changed to 1206 = carbon.
*
*  NAMELIST /WALL_SAVE/ 
*   Found in b2.wall_save.parameters.
*   Read if present in the run directory.
*   Automatically generated by B2 but can be modified by user.
*   Contains data describing the state of the walls lining the vacuum vessel.
*   These data are used for smooth restarts and are not really intended for 
*   human use.
*  NDEPTH_NML - integer. Dimension NDEPTH used for the arrays within this 
*    namelist. Represents the number of depth layer discretising the wall 
*    elements for the wall model. If using the 0-D model or the 
*    time-independent 1-D model, will contain 1 (default). Should not exceed
*    the value of the parameter NDEPTH declared in b2mod_wall.F.
*  IMAPX - integer array of size (NWALL). Indicates the (ix) position in the 
*    grid of wall element (iwall). Default ordering of the wall elements 
*    proceeds clockwise in physical space from the lower left corner.
*  IMAPY - integer array of size (NWALL). Indicates the (iy) position in the 
*    grid of wall element (iwall). Default ordering of the wall elements 
*    proceeds clockwise in physical space from the lower left corner.
*  XYMAP - integer array of size (-1:NX,-1:NY). XYMAP(ix,iy) contains the wall 
*    index (iwall) of the wall element located a grid locaion (ix,iy). 
*    If there is no wall element at this position, contains 0. 
*    Default ordering of the wall elements proceeds clockwise in physical 
*    space from the lower left corner.
*  SURFACE_MATERIAL_NAME - character*6 of size (NWALL). Contains the filename 
*    from which to extract the surface material properties for wall element 
*    (iwall). The files are to be found in 
*    $(SOLPSTOP)/data(.local)/Surface_properties/. Defaults to 'C'.
*  COATING_MATERIAL_NAME - character*6 of size (NWALL). Contains the filename 
*     from which to extract the bulk material properties of the eventual 
*     coating for wall element (iwall). Defaults to ' '.
*  BULK_MATERIAL_NAME - character*6 of size (NWALL). Contains the filename from
*     which to extract the bulk material properties for wall element (iwall). 
*     The files are to be found in 
*     $(SOLPSTOP)/data(.local)/Bulk_properties/. Defaults to 'C'.
*  LAYER_ALLOYS - character*6 of size (NALLOYS). Contains the filename from
*     which to extract the material properties for alloy (nalloy) which may
*     bw present in mixed materials deposited layers. Not yet operational.
*     The files are to be found in $(SOLPSTOP)/data(.local)/Bulk_properties/ 
*     and $(SOLPSTOP)/data(.local)/Surface_properties/.
*  TARGET_TEMP - real*8 array of size (NWALL,NDEPTH). Contains the temperature 
*     (in Kelvin) of wall element (iwall) at depth layer (idepth).
*     If plate_option.eq.1, will be set to backplate_temp(iwall).
*     If plate_option.eq.2 and empty, will be set to equilibrium 1-D profile
*     deduced from plasma incoming fluxes. 
*     If plate-option.eq.3, must be set.
*  INERTIAL_COOLING - logical array of size (NWALL). Indicates whether wall
*     element (iwall) is inertially cooled instead of actively cooled.
*     Defaults to .false..
*  BACKPLATE_TEMP - real*8 array of size (NWALL). Contains the temperature 
*     (in Kelvin) maintained by cooling at the back end of wall element 
*     (iwall). Defaults to 'b2stbr_plate_temp'.
*  PLATE_THICKNESS - real*8 array of size (NWALL). Contains the thickness (in 
*     meters) of wall element (iwall). Defaults to 'b2stbr_plate_thick'.
*  COATING_THICKNESS - real*8 array of size (NWALL). Contains the thickness
*      (in meters) of the eventual coating on wall element (iwall). 
*      Defaults to 0.
*  PLATE_TIME_FACTOR - real*8 array of size (NWALL). Multiplier to the time 
*     for the equations for temperature and composition evolution of wall 
*     element (iwall). Defaults to 1.0.
*  DEPOSITION - real*8 array of size (NWALL,NTRACK). Contains the amount of 
*     deposited material (in atoms) from species (itrack) onto wall element 
*     (iwall). Defaults to 0.0.
*  EROSION - real*8 array of size (NWALL,NTRACK). Contains the amount of eroded
*     material (in atoms) of species (itrack) from wall element (iwall).
*     Defaults to 0.0.
*  CHEMICAL_SPUTTERING - real*8 array of size (NWALL,0:NS-1). Contains the 
*     chemical sputter yield of 'sput_dst' or 'sput_dst_bulk' species on wall 
*     element (iwall) caused by B2 species (is). Defaults to 0.0.
*  PHYSICAL_SPUTTERING - real*8 array of size (NWALL,0:NS-1). Contains the 
*     physical sputter yield of 'sput_dst' or 'sput_dst_bulk' species on wall 
*     element (iwall) caused by B2 species (is). Defaults to 0.0.
*  PHYSICAL_SPUTTERING_ENERGY - real*8 array of size (NWALL,0:NS-1). Contains 
*     the fraction of returned energy carried by sputtered particles of species
*     'sput_dst' or 'sput_dst_bulk' species from wall element (iwall) caused by
*     B2 species (is). Defaults to 0.0.
*  RES_SPUTTERING - real*8 array of size (NWALL,0:NS-1). Contains the RES 
*     sputter yield of 'sput_dst' or 'sput_dst_bulk' species on wall element 
*     (iwall) caused by B2 species (is). Defaults to 0.0.
*  THERMAL_EVAPORATION - real*8 array of size (NWALL,0:NS-1). Contains the  
*     rate of thermal evaporation of species (is) (in particles/second) 
*     from wall element (iwall). Defaults to 0.0.
*  BACKSCATTERING - real*8 array of size (NWALL,0:NS-1). Contains the 
*     backscattering fraction for incoming B2 species (is) onto wall element 
*     (iwall).  Defaults to 0.0.
*  BACKSCATTERING_ENERGY - real*8 array of size (NWALL,0:NS-1). Contains the 
*     backscattered energy fraction for incoming B2 species (is) onto wall 
*     element (iwall).  Defaults to 0.0.
*  PLATE_TIME - real*8 array of size (NWALL). Indicates how much simulation 
*     time has elapsed for wall element (iwall) (in seconds). Defaults to 0.0.
*  PLATE_AREA - real*8 array of size (NWALL). Indicates the wall area (in 
*     square meters) for wall element (iwall). Defaults to the area computed
*     from the B2 grid information.
*  MONOLAYER_DEPOSITION - real*8 array of size (NWALL,NTRACK). Contains the 
*     amount of deposited material (in monolayers) from species (itrack) onto 
*     wall element (iwall). Defaults to 0.0.
*  MONOLAYER_EROSION - real*8 array of size (NWALL,NTRACK). Contains the amount
*     of eroded material (in monolayers) of species (itrack) from wall element 
*     (iwall). Defaults to 0.0.
*  LAYER_NCONSTITUENTS - integer array of size (NWALL). Indicates the number of
*     elemental constituents within the surface layer of wall element (NWALL).
*     Defaults to 1.
*  LAYER_NZCONSTITUENTS - integer array of size (NWALL,6+NTRACK). Contains the 
*     atomic numbers Z of the various elements present within the surface layer
*     of wall element (iwall). Defaults to 6 for the first value, 0 otherwise.
*  LAYER_NRELCONSTITUENTS - real*8 array of size (NWALL,6+NTRACK). Contains the
*     relative atomic abundances of the various elements present within the 
*     surface layer of wall element (iwall). Defaults to 1.0 for the first 
*     value, 0.0 otherwise.
*
* NAMELIST /B2MD_NAMELIST/ 
*   Found in b2md.dat and read by b2md.
*   Contains data useful for cataloguing the run (done from 'save_mds' and 
*   'resave_mds' scripts).
*  EXP - character*128. Name of the experiment being modelled. 
*     Defaults to 'NOT_SET'.
*  SHOT - integer. Shot number identifying the run. Defaults to the last 
*     number found in shotnumber.history, or 0 if the file is not found.
*  TIME - real*8. Time point of the experimental shot being simulated.
*     Defaults to 0.0.
*  COMMENT - character*128. Label for the run. Defaults to 'NOT_SET'.
*  TIMEDEP - logical. If .true. (default), saves data from b2time.nc.
*  SNAPSHOT - logical. If .true. (default), saves data from b2fplasma.
*  TALLIES - logical. If .true. (default), saves data from b2tallies.nc.
*  MOVIES - logical. If .true. (default), saves data from b2movies.nc.
*  OVERWRITE_SHOTNUMBER - integer. Indicate the shot number to overwrite (to be
*     used only when updating an already saved run with 'resave_mds' script). 
*     Defaults to the last number found in shotnumber.history, or 0 if the file 
*     is not found.
*
* NAMELIST /BOUNDARY/ 
*   Found in b2.boundary.parameters.
*   Contains data used to specify the boundary conditions.
*   This namelist is read if 'b2stbc_boundary_namelist' is set to 1.
*   Necessary for double-null geometries or impinging limiting structures.
*   It is recommended but not strictly necessary that the boundary conditions
*   declared here match in number and order any neutral boundary conditions
*   declared in the b2.neutrals.parameters file.
*   This description is copied from the SOLPS manual.
*  NBC - integer. Number of boundary segments. Defaults to 0.
*  BCCHAR - character*1 array of length (NBC). Specifies the nature of
*   the boundary segment (default ' '): 
*   'N' = ``North'' boundary
*   'S' = ``South'' boundary
*   'W' = ``West'' boundary
*   'E' = ``East'' boundary
*   'X' = ``X'' boundary used for specifying a fixed value on a row of cells
*   'Y' = ``Y'' boundary used for specifying a fixed value on a column of cells
*  CONPAR - real*8 array of size (0:NS-1,NBC,3). Contains parameters helping
*     to define the boundary conditions for the continuity equation of species
*     (is). See description of BCCON below for details. Defaults to 0.0.
*  MOMPAR - real*8 array of size (0:NS-1,NBC,2). Contains parameters helping
*     to define the boundary conditions for the parallel momentum equation of
*     species (is). See description of BCMOM below for details. Defaults to 0.0.
*  ENEPAR - real*8 array of size (NBC,2). Contains parameters helping
*     to define the boundary conditions for the electron energy equation.
*     See description of BCENE below for details. Defaults to 0.0.
*  ENIPAR - real*8 array of size (NBC,2). Contains parameters helping
*     to define the boundary conditions for the ion energy equation. 
*     See description of BCENI below for details. Defaults to 0.0.
*  POTPAR - real*8 array of size (NBC,2). Contains parameters helping
*     to define the boundary conditions for the potential equation.
*     See description of BCPOT below for details. Defaults to 0.0.
*  BCPOS - integer array, length (NBC), for North, South
*     or X boundary conditions, it specifies the row index; for West, East
*     and Y boundary conditions, it specifies the column index. 
*     Defaults to -2.
*  BCSTART - integer array, length (NBC), for North,
*     South or X boundary conditions, it specifies the start column index;
*     for West, East and Y boundary conditions, it specifies the start row
*     index. Defaults to -2.
*  BCEND - integer array, length (NBC), for North, South
*     or X boundary conditions, it specifies the end column index; for
*     West, East and Y boundary conditions, it specifies the end row index.
*     Defaults to -2.
*  BC_LIST_SIZE - integer array of length (NBC). Contains the size of the 
*     list of cells where a boundary condition is applied. Defaults to 0.
*  BC_LIST_X - integer array of length (2*(NXD+NYD),NBC). Contains the 
*     X-coordinate of the cells where boundaries conditions are applied. 
*     Defaults to -2.
*  BC_LIST_Y - integer array of length (2*(NXD+NYD),NBC). Contains the 
*     Y-coordinate of the cells where boundaries conditions are applied. 
*     Defaults to -2.
*  BCCON - integer array, length (0:NS-1,NBC), specifying
*     the type of density boundary condition for each segment and species
*     (fastest varying index is species); makes use of CONPAR to
*     specify additional information, as indicated:
*     0 : default, no boundary condition is applied
*     1 : prescribe the value of the density, CONPAR(,,1)
*     specifies the required density in m^{-3}
*     2 : prescribe the gradient of the density, CONPAR(,,1)
*     specifies the required density gradient in $m^{-4}$
*     3 : sheath conditions, CONPAR(,,1) not used (zero
*     gradient is used)
*     4 : prescribe the value of the density, weakly a mixed boundary
*     condition, CONPAR(,,1) specifies the required density in
*     $m^{-3}$ and CONPAR(,,2) specifies the ``strength'' of the
*     boundary condition
*     5 : prescribe the particle flux per unit area,
*     CONPAR(,,1) specifies the required particle flux density in
*     $m^{-2} s^{-1}$
*     6 : prescribe the total particle flux for a constant density,
*     CONPAR(,,1) specifies the particle flux in
*     $s^{-1}$
*     7 : prescribe the density as a function of other plasma
*     parameters [not yet available]
*     8 : prescribe the total particle flux with constant flux
*     density, CONPAR(,,1) specifies the particle flux in
*     $s^{-1}$
*     9 : prescribe the decay length for the density,
*     CONPAR(,,1) specifies the gradient length in $m$
*     (should use type 15 instead when drifts are turned on)
*     10 : leakage option for density, CONPAR(,,1) specifies
*     the leakage factor, $\alpha$ in $\Gamma_{\mbox{loss}} = \alpha C_s n_a$
*     11 : particle flux feedback boundary condition,
*     CONPAR(,,1) not used, derived from CBSNA(0,IS,IREG).
*     The species used must be declared using the 'b2stbc_isfeedback' switch.
*     12 : particle density feedback boundary condition,
*     CONPAR(,,1) not used, derived from CBSNA(0,IS,IREG).
*     The species used must be declared using the 'b2stbc_isfeedback' switch.
*     13 : particle density to achieve specified total flux,
*     CONPAR(,,1) is the specified flux crossing the flux surface type13_ref 
*     steps from the boundary,
*     CONPAR(,,2) is the strength of the feedback,
*     CONPAR(,,3) is set internally to match the re-entering flux of ionised
*     neutrals that crossed the core boundary (when running with Eirene
*     and the 'ionising_core' option).
*     The feedback scheme can be further tweaked with the switches
*     type13_norm and type13_fac. See code for details.
*     14 : sound speed velocity flux,
*     CONPAR(,,1) is a multiplier to the outgoing sound speed $C_s$.
*     15 : prescribe the decay length for the density with drifts switched on,
*     CONPAR(,,1) specifies the gradient length in $m$
*     16 : particle density to achieve specified total flux, used with ASTRA 
*     coupling. The total desired flux is summed over all BCCON=16 core boundaries.
*     CONPAR(,,1) specifies the desired particle flux in $s^{-1}$.
*     18 : prescribe total main ion particle flux, used with ASTRA coupling.
*     19 : particle flux feedback boundary condition, flux is summed over neutrals
*     and ions, for coupling with ASTRA. The total desired flux is summed over all 
*     BCCON=19 core boundaries. This boundary condition type is applied to ions
*     in their highest ionisation state.
*     CONPAR(,,1) specifies the desired particle flux in $s^{-1}$.
*     20 : constant density feedback condition.
*     CONPAR(,,1) specifies the desired density in $m^{-3}$.
*     21 : prescribe the value of the density and add a density perturbation 
*     to get a solution which is as close as possible to neoclassical theory.
*     It is recommended to use this boundary condition together with
*     corresponding condition on ion temperature (BCENI=24). 
*     The total desired density is summed over all BCCON=21 core boundaries.
*     CONPAR(,,1) specifies the desired density in $m^{-3}$.
*     22 : Feedback boundary condition: given total particle flux with constant
*     average density. The poloidal density variation is added to make the solution
*     as close as possible to the neoclassical value derived for pure H/D/T plasma.
*     It is recommended to use this boundary condition with corresponding boundary
*     condition on ion temperature (BCENI=23,24).
*     The total desired flux is summed over all BCCON=22 core boundaries.
*     CONPAR(,,1) specifies the desired particle flux in $s^{-1}$.
*     23 : Feedback boundary condition: given sum of integrated neutrals and
*     main ion particle fluxes with constant average density.
*     The poloidal density variation is added to make the solution
*     as close as possible to the neoclassical value derived for pure H/D/T plasma.
*     It is recommended to use this boundary condition with corresponding
*     boundary condition on ion temperature (BCENI=23,24).
*     The total desired flux is summed over all BCCON=23 core boundaries.
*     CONPAR(,,1) specifies the desired particle flux in $s^{-1}$.
*     24 : constant density feedback scaled by density on the ring bc_type21_ref away
*     CONPAR(,,1) specifies the desired density in $m^{-3}$.
*     CONPAR(,,2) is the strength of the feedback
*     25 : Feedback boundary condition: prescribe the average value of the density
*     and add a density perturbation from neighbouring radial cell.
*     This boundary condition is suitable for any species of a
*     multi-species plasma (i.e. in the case when the condition used
*     in BCCON=21 fails).
*     It is recommended to use this boundary condition together with the
*     corresponding condition on ion temperature (BCENI=26,27).
*     CONPAR(,,1) specifies the desired average density in $m^{-3}$.
*     26 : Feedback boundary condition: prescribe the total ion flux and find
*     the average density.
*     A density perturbation is taken from the neighbouring radial cell.
*     This boundary condition is suitable for any species of a
*     multi-species plasma (i.e. in the case when the condition used
*     in BCCON=22 fails).
*     It is recommended to use this boundary condition together with
*     a corresponding condition on the ion temperature (BCENI=26,27).
*     The total desired flux is summed over all BCCON=26 core boundaries.
*     CONPAR(,,1) specifies the desired particle flux in $s^{-1}$.
*     27 : Feedback boundary condition: prescribe the particle flux sum
*     for all neutrals and ions belonging to a given isonuclear sequence
*     and find the average density of the highest ionization stage.
*     A density perturbation is taken from the neighbouring radial cell.
*     It is recommended to use this boundary condition together with
*     corresponding condition on ion temperature (BCENI=25,26).
*     The total desired flux is summed over all BCCON=27 core boundaries.
*     CONPAR(,,1) specifies the desired particle flux in $s^{-1}$.
*  BCMOM - integer array, length NS * NBC, specifying
*     the type of parallel momentum or velocity boundary condition for each
*     segment and species (fastest varying index is species); makes use of
*     MOMPAR to specify additional information, as indicated
*     0 : default, no boundary condition is applied
*     1 : prescribe the value of the parallel velocity,
*     MOMPAR(,,1) specifies the parallel velocity in $m.s^{-1}$
*     2 : prescribe the gradient of the parallel velocity,
*     MOMPAR(,,1) specifies the parallel velocity gradient in $s^{-1}$
*     3 : sheath conditions, mach number as input, if
*     MOMPAR(,,2) < 0.5, then the velocity is set to exactly
*     MOMPAR(,,1) $\times C_{s, \mbox{collective}}$, otherwise the
*     velocity is set to be at least MOMPAR(,,1) $\times
*     C_{s,\mbox{species}}$
*     4 : prescribe the value of the velocity, weakly a mixed
*     boundary condition, MOMPAR(,,1) specifies the parallel
*     velocity in $m.s^{-1}$ and MOMPAR(,,2) specifies the
*     ``strength'' of the boundary condition
*     5 : prescribe the parallel momentum flux per unit area,
*     MOMPAR(,,1) specifies the parallel momentum flux density in $N.m^{-2}$
*     6 : prescribe the total parallel momentum flux for a constant
*     parallel velocity [not yet available]
*     7 : prescribe the parallel momentum as a function of other plasma
*     parameters [not yet available]
*     8 : special : limited shear, imposes zero gradient for the Mach number.
*     [[[Eventually intended to have MOMPAR(,,1) specify the
*     gradient of the Mach number in $m^{-1}$]]]
*     9 : prescribe the total parallel momentum flux with constant flux
*     density, MOMPAR(,,1) specifies the parallel momentum flux in $N$
*     10 : prescribe the decay length for the parallel momentum,
*     MOMPAR(,,1) specifies the decay length in $m$
*     11 : Rozhansky viscosity condition for the parallel momentum,
*     MOMPAR(,,1) is not used
*     12 : Condition from b2stbc_spb for the parallel momentum
*     13 : sheath boundary condition from b2stbc_spb for the parallel momentum
*     14 : Condition from b2stbc_spb for the parallel momentum
*     15 : Prescribe the value of the parallel velocity as a function of
*     poloidal coordinate
*     16 : Prescribe the average value of the parallel velocity
*     MOMPAR(,,1) specifies the parallelv elocity in $m.s^{-1}$
*     BCENE - integer array, length NBC, specifying
*     the type of electron energy or temperature boundary condition for
*     each segment; makes use of ENEPAR to specify additional
*     information, as indicated
*     0 : default, no boundary condition is applied
*     1 : prescribe the value of the electron temperature,
*     ENEPAR(,1) specifies the temperature in $eV$
*     2 : prescribe the gradient of the electron temperature,
*     ENEPAR(,1) specifies the temperature gradient in $eV.m^{-1}$
*     3 : sheath conditions, electron energy transmission,
*     ENEPAR(,1) specifies an additional contribution to the energy
*     transmission coefficient in addition to that of the potential
*     difference [the sound speed used depends on settings of
*     MOMPAR(,ISMAIN,2)]
*     4 : prescribe the value of the electron temperature, weakly a
*     mixed boundary condition, ENEPAR(,1) specifies the
*     temperature in $eV$ and ENEPAR(,2) specifies the
*     ``strength'' of the boundary condition
*     5 : prescribe the electron energy flux per unit area,
*     ENEPAR(,1) specifies the energy flux density in $W.m^{-2}$
*     6 : prescribe the total electron energy flux for a constant
*     electron temperature, ENEPAR(,1) specifies the energy flux
*     in $W$ (no longer supported, use types [16] or [17] instead)
*     7 : prescribe the electron temperature as a function of
*     other plasma parameters [not yet available]
*     8 : prescribe the total electron heat flux with constant flux
*     density, ENEPAR(,1) specifies the energy flux in $W$
*     9 : prescribe the decay length for the electron temperature,
*     ENEPAR(,1) specifies the decay length in $m$ (can also use type [19] 
*     instead)
*     10 : feedback option for core, ENEPAR(,1) not used,
*     derived from CBSHE(0,coreregno)
*     11 : not used
*     12 : sheath conditions, electron energy transmission coefficient,
*     ENEPAR(,1) specifies an energy transmission factor,
*     $\delta_e$ in $Q_e = \delta_e \Gamma_e T_e$
*     13 : prescribe the electron energy flux per unit area
*     proportional to temperature, ENEPAR(,1) specifies the energy
*     flux density per temperature in $W.m^{-2}.J^{-1}$ (the temperature
*     here in $J$)
*     14 : leakage option for electron energy, ENEPAR(,1) specifies
*     the leakage factor, $\alpha$ in $\Gamma_{\mbox{loss}} = \alpha C_s n_e$
*     15 : not used
*     16 : Feedback boundary condition with constant temperature,
*     ENEPAR(,1) specifies the power flux in $W$ across the flux surface
*     with index 'b2stbc_type16_ref' (default=-1), ENEPAR(,2)
*     should be something like 0.1 and specifies the strength of the
*     feedback. Also see type [17] below.
*     Available if bcene_16_style.eq.0 (default).
*     If bcene_16_style.eq.1, integrated electron heat flux with constant Te,
*     summed over all core boundaries with BCENE=16.
*     17 : Feedback boundary condition with constant shared temperature
*     for both electrons and ions, with ENEPAR(,1) + ENIPAR(,1) giving 
*     the total power flux in $W$ across the flux surface
*     with index 'b2stbc_type16_ref' (default=-1), ENEPAR(,2)
*     should be something like 0.1 and specifies the strength of the
*     feedback. Replaces [16] for high densities and large values of 
*     'b2stbc_type16_ref'.
*     18 : Fractional drop condition. Not yet working.
*     19 : same as [9] but to be used when simultaneously setting BCCON=1 on 
*     the same boundary.
*     20 : Feedback boundary condition with constant temperature,
*     as per type [16] but with type20_ switches and a different
*     feedback scheme.      
*     21 : constant temperature feedback scaled by temperature on the ring 
*     bc_type21_ref away
*     ENEPAR(,,1) specifies the desired electron temperature in $eV$.
*     ENEPAR(,,2) is the strength of the feedback
*     22 : modified decay length for the electron temperature
*  BCENI - integer array, length NBC, specifying the
*     type of ion energy or temperature boundary condition for each
*     segment; makes use of ENIPAR to specify additional
*     information, as indicated
*     0 : default, no boundary condition is applied
*     1 : prescribe the value of the ion temperature,
*     ENIPAR(,1) specifies the temperature in $eV$
*     2 : prescribe the gradient of the ion temperature,
*     ENIPAR(,1) specifies the temperature gradient in $eV.m^{-1}$
*     3 : sheath conditions, ion energy transmission,
*     ENIPAR(,1) specifies the contribution to the energy
*     transmission coefficient [the sound speed used depends on settings
*     of MOMPAR(,ISMAIN,2)]
*     4 : prescribe the value of the ion temperature, weakly a mixed
*     boundary condition, ENIPAR(,1) specifies the temperature in
*     $eV$ and ENIPAR(,2) specifies the ``strength'' of the
*     boundary condition
*     5 : prescribe the ion energy flux per unit area,
*     ENIPAR(,1) specifies the energy flux density in $W.m^{-2}$
*     6 : prescribe the total ion energy flux for a constant ion
*     temperature, ENIPAR(,1) specifies the energy flux in $W$
*     (no longer supported, use types [16] or [17] instead)
*     7 : prescribe the ion temperature as a function of other plasma
*     parameters [not yet available]
*     8 : prescribe the total ion heat flux with constant flux
*     density, ENIPAR(,1) specifies the energy flux in $W$
*     9 : prescribe the decay length for the ion temperature,
*     ENIPAR(,1) specifies the decay length in $m$ (can also use type [19] 
*     instead)
*     10 : feedback option for core, ENIPAR(,1) not used,
*     derived from CBSHI(0,ISMAIN,coreregno)
*     11 : sheath conditions, ion energy transmission coefficient,
*     ENIPAR(,1) specifies an energy transmission factor,
*     $\delta_i$ in $Q_i = \delta_i T_i \sum_a n_a C_{s,a}$
*     12 : sheath conditions, ion energy transmission coefficient,
*     ENIPAR(,1) specifies an energy transmission factor,
*     $\delta_i$ in $Q_i = \delta_i T_i \sum_a \Gamma_a$
*     13 : prescribe the ion energy flux per unit area proportional
*     to temperature, ENIPAR(,1) specifies the energy flux density
*     per temperature in $W.m^{-2}.J^{-1}$ (the temperature here in $J$)
*     14 : leakage option for ion energy, ENIPAR(,1) specifies
*     the leakage factor, $\alpha$ in $\Gamma_{\mbox{loss}} = \alpha C_s n_a$
*     15 : from b2stbc_spb
*     16 : Feedback boundary condition with constant temperature,
*     ENIPAR(,1) specifies the power flux in $W$ across the flux surface
*     with index 'b2stbc_type16_ref' (default=-1), ENIPAR(,2)
*     should be something like 0.1 and specifies the strength of the
*     feedback. Also see type [17] below.
*     Available if bceni_16_style.eq.0 (default).
*     If bceni_16_style.eq.1, integrated ion heat flux with constant Ti,
*     summed over all core boundaries with BCENI=16.
*     ENIPAR(,1) specifies the power flux in $W$
*     17 : Feedback boundary condition with constant shared temperature
*     for both electrons and ions, see BCENE=17 above for description.
*     18 : Fractional drop condition. Not yet working.
*     19 : same as [9] but to be used when simultaneously setting BCCON=1 on 
*     the same boundary.
*     20 : Feedback boundary condition with constant temperature,
*     as per type [16] but with type20_ switches and a different
*     feedback scheme.      
*     21 : from b2stbc_spb
*     22 : modified decay length for the ion temperature
*     23 : Prescribe the poloidally averaged value of the ion temperature
*     and introduce a poloidal variation as close as possible to neoclassical 
*     solution. It is recommended to use this boundary condition together with
*     corresponding condition on ion density (BCCON=21,22,23).
*     The average is taken over all core boundaries with BCENI=23.
*     ENIPAR(,1) specifies the temperature in $eV$
*     24 : Feedback boundary condition with prescribed total ion flux,
*     constant poloidally averaged ion temperature and a
*     poloidal variation as close as possible to neoclassical solution.
*     It is recommended to use this boundary condition together with
*     corresponding condition on ion density (BCCON=21,22,23).
*     The flux is summed over all core boundaries with BCENI=24.
*     ENIPAR(,1) specifies the energy flux in $W$
*     25 : constant temperature feedback scaled by temperature on
*     the ring bc_type21_ref away
*     ENIPAR(,1) specifies the desired ion temperature in $eV$.
*     ENIPAR(,2) is the strength of the feedback
*     26 : prescribe the poloidally averaged value of the ion temperature
*     and introduce a poloidal variation in a simplified manner.
*     This boundary condition is suitable for any plasma composition
*     (i.e. when BCENI=23 fails).
*     It is recommended to use this boundary condition together with corresponding
*     condition on ion density (BCCON = 25,26,27).
*     The average is taken over all core boundaries with BCENI=26.
*     ENIPAR(,1) specifies the temperature in $eV$
*     27 : Feedback boundary condition with prescribed total ion heat flux,
*     constant poloidally averaged ion temperature and a
*     poloidal variation in a simplified manner.
*     This boundary condition is suitable for any plasma composition
*     (i.e. when BCENI=24 fails).
*     It is recommended to use this boundary condition together with corresponding
*     condition on ion density (BCCON = 25,26,27).
*     The flux is summed over all core boundaries with BCENI=27.
*     ENIPAR(,1) specifies the energy flux in $W$
*  BCPOT - integer array, length NBC, specifying the
*     type of electric potential or current boundary condition for each
*     segment; makes use of POTPAR to specify additional
*     information, as indicated
*     0 : default, no boundary condition is applied
*     1 : prescribe the value of the potential, POTPAR(,1)
*     specifies the potential in $V$
*     2 : prescribe the gradient of the potential, POTPAR(,1)
*     specifies the potential gradient in $V.m^{-1}$
*     3 : sheath conditions, 
*     POTPAR(,2) used for biasing [see code for details]
*     4 : prescribe the value of the potential weakly a mixed
*     boundary condition, POTPAR(,1) specifies the potential in
*     $V$ and POTPAR(,2) specifies the ``strength'' of the
*     boundary condition
*     5 : prescribe the current flux density per unit area,
*     POTPAR(,1) specifies the electric current flux density in
*     $A.m^{-2}$
*     6 : prescribe the total current flux density for a
*     constant potential [not yet available]
*     7 : prescribe the potential as a function of other
*     plasma parameters [not yet available]
*     8 : prescribe the total electric current with constant flux
*     density, POTPAR(,1) specifies the electric current in $A$
*     9 : prescribe the decay length for the potential,
*     POTPAR(,1) specifies the decay length in $m$
*     10 : feedback option for core [not yet tested!!!!!!!!!] (based
*     on using CBSCH(0,coreregno))
*     11 : sheath conditions, electron energy transmission from b2stbc_spb
*     POTPAR(,2) specifies the bias potential in $V$
*     12 : Imposes the currents due to drifts for the South core boundary.
*     Must be used in conjunction with istyle_cur_contr_on_S_and_N.eq.2
*     13 : Imposes the currents due to drifts for the South private flux
*     and North boundaries.
*     Must be used in conjunction with istyle_cur_contr_on_S_and_N.eq.2
*     16 : Constant electric potential feedback on imposed total current
*     The current prescribed is given by the sum of the POTPAR(IB,1) (in A)
*     over all the BCPOT=16 boundaries. 
*     Still experimental, will not work for drift cases.
*     21 : constant potential feedback scaled by potential on
*     the ring bc_type21_ref away.
*     POTPAR(,,1) specifies the desired potential in $V$.
*     POTPAR(,,2) is the strength of the feedback
*  GAMMAI - real*8. Ratio of specific heats (adiabatic coefficient).
*     Defaults to 5/3.
*  GAMMAE - real*8. Secondary electron emission coefficient. Defaults to 0.5.
*  LBNDUSR - logical. If .true. will also call Bas' boundary condition 
*     routine after the end of the physics boundary condition routine
*     (governed by the data from b2ah.dat and b2mn.dat). 
*     Defaults to .false..
*  LFEEDBACK - logical. Indicates whether a feedback scheme is used.
*     Obsolete. Superceded by 'b2stbc_feedback'. Defaults to .false..
*  NNISO - integer. Number of dead (or isolated) regions. Defaults to 0.
*  NIISO - real*8 array of size (0:NS-1). Density of species (is) in (m-3)
*     to impose in isolated regions.
*  TEISO - real*8. Electron temperature (in eV) to impose in isolated regions.
*     Defaults to 1.0.
*  TIISO - real*8. Ion temperature (in eV) to impose in isolated regions.
*     Defaults to 1.0.
*  PHIISO - real*8. Electric potential (in V) to impose in isolated regions.
*     Defaults to 0.0.
*  NXISO1 - integer array of size NNISO. Column number of bottom left corner 
*     of the dead region (iiso). Defaults to -2.
*  NXISO2 - integer array of size NNISO. Column number of top right corner 
*     of the dead region (iiso). Defaults to -2.
*  NYISO1 - integer array of size NNISO. Row number of bottom left corner 
*     of the dead region (iiso). Defaults to -2.
*  NYISO2 - integer array of size NNISO. Row number of top right corner 
*     of the dead region (iiso). Defaults to -2.
*  BOUNDARY_FILENAME - character*256. Name of the next file to use for reading 
*     a new /BOUNDARY/ namelist. Defaults to 'b2.boundary.parameters'.
*  BOUNDARY_TIME_MOD - real*8. When the B2 run simulation time, in seconds,
*     modulo(BOUNDARY_TIME_MOD), exceeds BOUNDARY_TIME_SWITCH, reads the new 
*     namelist from BOUNDARY_FILENAME. Also switches to the new namelist as
*     the ELM count (here time/BOUNDARY_TIME_MOD) changes.
*     Only active if BOUNDARY_TIME_MOD is greater than 0. Defaults to 0.0.
*  BOUNDARY_TIME_SWITCH - real*8. Time (in seconds) within an ELM cycle after 
*     which a new /BOUNDARY/ namelist is read. Defaults to 0.0.
*  LCBS - integer array of size (NBC). Indices of the core boundary segments
*     in B2 for passing to EIRENE. Defaults to the list of 'S' boundaries in
*     regions 1 and 5.
*  WRITE_NML_BND - logical. If .true. (default), writes the content of the
*    namelist to stdout after it has been read.
*
*  NAMELIST /FEEDBACK/ 
*   Found in b2.feedback_save.parameters.
*   Read if present in the run directory.
*   Used internally by B2 to regulate feedback loops and restarts.
*   Not really intended for human use.
*   See b2cdci and SOLPS manual for description of the feedback schemes.
*  SAVED_CBSHE_CORE - real*8. Last value used for the core electron energy 
*     radial flux feedback. Defaults to 0.0.
*  SAVED_CBSHI_CORE - real*8. Last value used for the core ion energy radial 
*     flux feedback. Defaults to 0.0.
*  SAVED_CBSNA_CORE - real*8. Last value used for the core particle radial flux
*     feedback. Corresponds to 'isfeedback' B2 species. Defaults to 0.0.
*  SAVED_CBSCH_CORE - real*8. Last value used for the core radial electric 
*     current feedback. Defaults to 0.0.
*  SAVED_CBSNA_PFR1 - real*8. Last value used for the radial inner PFR boundary
*     flux feedback. Corresponds to 'isfeedback' B2 species. Defaults to 0.0.
*  SAVED_CBSNA_PFR2 - real*8. Last value used for the radial outer PFR boundary
*     flux feedback. Corresponds to 'isfeedback' B2 species. Defaults to 0.0.
*  SAVED_CBSNA_SOL - real*8. Last value used for the radial particle flux 
*     feedback in the SOL. Corresponds to 'isfeedback' B2 species. Defaults 
*     to 0.0.
*
*  NAMELIST /FEEDBACK_CONTROL/
*   Found in b2.feedback_control.parameters.
*   Allows for more complex feedback schemes.
*  VACUUM_COMMUNICATION - Integer. Defaults to 0. If >0, allows for a
*   communication of particle fluxes across vacuum regions. This option only
*   applies to neutrals. The density boundary condition is based on the
*   difference between the average pressure and the local pressure.
*  VACUUM_COMMUNICATION_NREG - Integer array of size (NVAC). Defaults to 0.
*   Number of communicating vacuum regions.
*  VACUUM_COMMUNICATION_METHOD - Integer array of size (NVAC). Defaults to 0.
*   Option for resorbing the pressure difference.
*   1: Try to set a flux. Corr = (beta*ave_pressure - pressure)/temp * alpha
*   2: Try to set a density based on pressure equality.
*      Corr = alpha * beta * pressure / Ti
*  VACUUM_COMMUNICATION_IY - Integer array of size (NVACREG,NVAC).
*   Radial index of the ring on which the neutral pressure is computed for
*   region IREG. Defaults to -2.
*  VACUUM_COMMUNICATION_IX1 - Integer array of size (NVACREG,NVAC).
*   Poloidal lower bound of the range over which the neutral pressure is
*   computed for region IREG. Defaults to -2.
*  VACUUM_COMMUNICATION_IX2 - Integer array of size (NVACREG,NVAC).
*   Poloidal upper bound of the range over which the neutral pressure is
*   computed for region IREG. Defaults to -2.
*  VACUUM_COMMUNICATION_ALPHA - Real*8 array of size (0:NSPECIES-1,NVAC).
*    Parameter for setting the pressure correction. See above. Defaults to 0.0.
*  VACUUM_COMMINICATION_BETA - Real*8 array of size (0:NSPECIES-1,NVAC).
*    Parameter for setting the pressure correction. See above. Defaults to 1.0.
*  NA_FEEDBACK_TARGET - Real*8 array of size (0:NSPECIES-1). Sets the target
*    density of species (ISPECIES) for the feedback scheme. Defaults to 0.0.
*  NA_FEEDBACK_TIME - Real*8 array of size (0:NSPECIES-1). Sets the time
*    of reference (in s) for the feedback of species (ISPECIES).
*    Defaults to 0.0.
*  NA_FEEDBACK_CHOICE - Integer array of size (0:NSPECIES-1). Defaults to 0.
*    Choice of quantity on which the feedback is computed:
*    0: no action
*    1: local species density (averaged over the rectangle of cells
*       [IX1:IX2,IY1:IY2]). Summed over all charge states of that species.
*    2: local electron density (averaged over the rectangle of cells
*       [IX1:IX2,IY1:IY2]).
*    3: outer midplane separatrix electron density.
*    4: total particle content for that species.
*    5: total ion content for that species (not including neutrals).
*    6: neutral particle flux through the core boundary.
*    7: relative average concentration of this species at the separatrix.
*  NA_FEEDBACK_OPTION - Integer array of size (0:NSPECIES-1). Defaults to 0.
*    Option for computing the new fedback quantity.
*    0: no action
*    1: rescale slowed by na_feedback_alpha
*    2: pure rescale
*    3: rescaling slowed by tanh_log
*    4: rescale done according to SOLPS4 formula
*       The target waveform for the particle content is
*        N = C + V*(time-T)
*       and the current puffing rate S is adjusted
*        S --> max(0, min(X,S + D)), where
*        D = F*((N - <N>)/dt + (<N>_prev - <N>)/dt_prev)
*    5: rescale done according to SOLPS4 formula: N = C*exp((time-T)*V)
*    6: rescale slowed by na_feedback_alpha (SOLPS4 style)
*  NA_FEEDBACK_ACTUATOR - Integer array of size (0:NSPECIES-1). Defaults to 0.
*    Choice for the actuator used for the feedback.
*    0: no action
*    1: gas puff via boundary condition
*    2: rescale na
*    3: core boundary flux condition
*  NA_FEEDBACK_ALPHA - Real*8 array of size (0:NSPECIES-1). Defaults to 0.001.
*    Factor by which the rescaling is slowed. Rescaling factor is :
*    Option 1: (1 + alpha*target/current) / (1 + alpha)
*    Option 3: 2**(tanh(log(x)/beta)*log(alpha)/log(2))
*    Options 4 & 5: Corresponds to parameter F
*  NA_FEEDBACK_BETA - Real*8 array of size (0:NSPECIES-1). Defaults to 1.0.
*    Factor by which the rescaling is slowed. See above.
*    Options 4 & 5: Corresponds to parameter V (ffb_rtvn)
*  NA_FEEDBACK_CONST - Real*8 array of size (0:NSPECIES-1). Defaults to 0.0.
*    Options 4 & 5: Corresponds to parameter C. If negative, C is
*                   computed as the initial total particle content of the
*                   sequence.
*  NA_FEEDBACK_IX1 - Integer array of size (0:NSPECIES-1). Lower poloidal
*    bound for the region of which the density is being averaged.
*    Defaults to -2.
*  NA_FEEDBACK_IX2 - Integer array of size (0:NSPECIES-1). Upper poloidal
*    bound for the region of which the density is being averaged.
*    Defaults to -2.
*  NA_FEEDBACK_IY1 - Integer array of size (0:NSPECIES-1). Lower radial
*    bound for the region of which the density is being averaged.
*    Defaults to -2.
*  NA_FEEDBACK_IY2 - Integer array of size (0:NSPECIES-1). Upper poloidal
*    bound for the region of which the density is being averaged.
*    Defaults to -2.
*  NA_FEEDBACK_IB - Integer array of size (0:NSPECIES-1). Defaults to -1.
*    Index of boundary condition through which the feedback is being applied.
*  NA_FEEDBACK_PUFF_MIN - Real*8 array of size (0:NSPECIES-1). Defaults to 0.0.
*    Minimum gas puff being applied.
*  NA_FEEDBACK_PUFF_MAX - Real*8 array of size (0:NSPECIES-1). Defaults to 0.0.
*    Maximum gas puff being applied.
*  NA_FEEDBACK_OVERSHOOT - Real*8 array of size (0:NSPECIES-1).
*    If the density is larger than target*overshoot, the gas puff is
*    turned off. Defaults to 0.0.
*
*  NAMELIST /PROFILE/
*   Found in b2.sources.profile.
*   Allows the user to set external sources to use within the code.
*   Unless specifying a divertor heat source, the source will be
*   located in the core and main SOL regions only.
*   Used if 'b2sral_inputfile' is set to 1.
*  NSDATA - integer array of size (NKIND_DATA,NKIND_SOURCE,0:NS). Number of 
*     points over which the source profile of (kind_data,kind_source,is) is 
*     defined. Should not exceed NY+2. Defaults to 0.
*     If KIND_DATA=1, the data is expressed as a profile in physical distance 
*     from the separatrix (in metres) along the outer midplane.
*     The user can change this default reference location by use of the 
*     'set_transport_i[xy]ref' switches. 
*     If KIND_DATA=2, the data is expressed as a profile in flux distance from 
*     the separatrix (not yet available).
*  SDATA - real*8 data of size (2,NY+2,NKIND_SOURCE,0:NS). 
*     For (i,ir,ik,is) in (1:2,1:NY+2,1:NKIND_SOURCE,0:NS), 
*     SDATA(1,ir,:,:) contains the radial location of the profile point (ir).
*     SDATA(2,ir,:,:) contains the source profile value at point (ir).
*     KIND_SOURCE=1 means particle source of species (is) (in particles/m3)
*     KIND_SOURCE=2 means parallel momentum source for species (is) (in kg.m/s/m3)
*     KIND_SOURCE=3 means electron heat source (in Watts/m3)
*     KIND_SOURCE=4 means ion heat source (in Watts/m3)
*     KIND_SOURCE=5 means electric charge source (in Coulombs/m3)
*     KIND_SOURCE=6 means non-ambipolar electron particle source (in e/m3)
*     Defaults to 0.0.
*  NXDATA - integer array of size (NKIND_DATA,NKIND_SOURCE,0:NS). Number of 
*     points over which the axial profile of (kind_data,kind_source,is) is 
*     defined. Should not exceed NX+2. Defaults to 0.
*     If KIND_DATA=1, the data is expressed as a profile in physical distance,
*     here connection length, rescaled from 0.0 to 1.0.
*     For closed field lines, the reference location for the zero of distance 
*     is set by use of the 'set_transport_i[xy]ref' switches. 
*     If KIND_DATA=2, the data is expressed as a profile in (ix) cell index,
*     again normalized from 0.0 to 1.0 to match the [0:nx-1] interval.
*  XDATA - real*8 data of size (2,NY+2,NKIND_SOURCE,0:NS). 
*     Multiplier to the poloidal source profile in the axial direction.
*     Same convention for KIND_SOURCE as above.
*     Defaults to 1.0.
*  DIVHEAT - real*8. Additional divertor ion heat source (in Watts/m3)
*     Defaults to 0.0.
*  SOURCES_FILENAME - character*256. Name of the next file to use for reading
*     a new /PROFILE/ namelist. Defaults to 'b2.sources.profile'.
*     Quantities not present in the new file will be inherited from the old one.
*  SOURCES_TIME_MOD - real*8. When the B2 run simulation time, in seconds,
*     modulo(SOURCES_TIME_MOD),exceeds SOURCES_TIME_SWITCH, reads the new 
*     namelist from SOURCES_FILENAME. Also switches to the new namelist if 
*     the ELM count (here time/SOURCES_TIME_MOD) changes. Defaults to 0.0.
*     Only active if SOURCES_TIME_MOD is greater than 0.
*  SOURCES_TIME_SWITCH - real*8. Time (in seconds) within an ELM cycle after 
*     which a new /PROFILE/ namelist is read. Defaults to 0.0.
*
*  NAMELIST /TRANSPORT/
*   Found in b2.transport.inputfile.
*   Allows the user to set transport coefficients profiles.
*   Only used is 'b2tqna_inputfile' is set to 1.
*   Overwrites b2ah.dat, b2mn.dat and b2.transport.parameters.
*   See SOLPS manual for details.
*  NDATA - integer array of size (NKIND_DATA,NCOEF,0:NS). Number of points over
*     which the source profile of (kind_data,kind_coef,is) is defined. Should 
*     not exceed NY+2. Defaults to 0.
*     If KIND_DATA=1, the data is expressed as a profiles in physical distance 
*     from the separatrix (in metres).
*     If KIND_DATA=2, the data is expressed as a profile in flux distance from
*     the separatrix (not yet available).
*  TDATA - real*8 data of size (3,NY+2,NKIND_COEFF,0:NS).  Defaults to 0.0.
*     For (i,ir,ik,is) in (1:3,1:NY+2,1:NKIND_COEFF,0:NS), 
*     TDATA(1,ir,:,:) contains the radial location of the profile point (ir).
*     By default, these are measured at the outer midplane, as distance to the separatrix
*     in metres. The user can change this by setting the 'set_transport_i[xy]ref'
*     switches to choose a different location for the distance reference. 
*     TDATA(2,ir,:,:) contains the transport profile value at point (ir).
*     TDATA(3,ir,:,:) contains the ELM transport profile value at point (ir). 
*     KIND_COEFF=1 means density-driven particle diffusivity for species (is)
*     KIND_COEFF=2 means pressure-driven particle diffusivity for species (is)
*     KIND_COEFF=3 means ion heat diffusivity for species (is)
*     KIND_COEFF=4 means electron heat diffusivity
*     KIND_COEFF=5 means poloidal pinch velocity for species (is)
*     KIND_COEFF=6 means radial pinch velocity for species (is)
*     KIND_COEFF=7 means viscosity for species (is)
*     KIND_COEFF=8 means field-driven radial current conductivity
*     KIND_COEFF=9 means temperature-driven radial current conductivity
*  ADDSPEC - integer array of size (NS,NKIND_COEFF,0:NS). If 
*     ADDSPEC(is,ikind,spec).ge.0, then the transport coefficient profile of 
*     type (ikind) from species (spec) is also used for species with index 
*     ADDSPEC(is,ikind,spec). Defaults to -5.
*  TRANSPORT_FILENAME - character*256. Name of the next file to use for reading
*     a new /TRANSPORT/ namelist. Defaults to 'b2.transport.inputfile'.
*  TRANSPORT_TIME_MOD - real*8. When the B2 run simulation time, in seconds,
*     modulo(TRANSPORT_TIME_MOD),exceeds TRANSPORT_TIME_SWITCH, reads the new 
*     namelist from TRANSPORT_FILENAME. Also switches to the new namelist if 
*     the ELM count (here time/TRANSPORT_TIME_MOD) changes. Defaults to 0.0.
*     Only active if TRANSPORT_TIME_MOD is greater than 0.
*  TRANSPORT_TIME_SWITCH - real*8. Time (in seconds) within an ELM cycle after 
*     which a new /TRANSPORT/ namelist is read. Defaults to 0.0.
*  REGION_FLAGS - logical array of size (NREG,NKIND_COEFF). If 
*     region_flags(ireg,ikind) is .true. (default), then the transport 
*     parameters profiles for region (ireg) and kind (ikind) are used.
*  NO_PFLUX - logical. If .true., transport coefficients profiles are not 
*     implemented in the private flux regions. Defaults to .false..
*  POLOIDAL_SCALING - logical array of size (10).  If .true., then the 
*     transport coefficients profiles from the current
*     b2.transport.inputfile are increased by a factor of 1.0+Gaussian
*     where Gaussian is a Gaussian profile in the poloidal direction of
*     amplitude SCALING_STRENGTH extending from SCALING_IX_BEGIN to 
*     SCALING_IX_END inclusively. Defaults to .false..  The profile has
*     a decay length of SCALING_WIDTH (in units of the number of poloidal
*     cells).
*  SCALING_STRENGTH - real*8 array of size 10. See above. Defaults to zero.
*  SCALING_WIDTH - real*8 array of size 10. See above. Defaults to about
*     1/3 of the interval over which the scaling is to be done.
*  SCALING_IX_BEGIN - integer array of size 10. See above. Defaults to -2.
*  SCALING_IX_END - integer array of size 10. See above. Defaults to -2.
*  ELM_TIME_BEGIN - real*8. Time (in seconds) modulo ELM_TIME_PERIOD at
*     which the ELM phase begins and the ELM data must be used.
*     Defaults to 0.0.
*  ELM_TIME_END - real*8. Time (in seconds) modulo ELM_TIME_PERIOD at
*     which the ELM phase ends and the ELM data is no longer used.
*     Defaults to 0.0.
*  ELM_TIME_PERIOD - real*8. Indicates the real frequency of simulated ELMs.
*     See above for usage. Defaults to 0.0.
*     If zero, no ELM profiles are used.
*  ELM_IX_BEGIN - integer. Poloidal position at which the ELM profile
*     starts to be applied. Defaults to -2.
*  ELM_IX_END - integer. Poloidal position at which the ELM profile
*     ends being applied. Defaults to -2.
*
*  NAMELIST /NEUTRALS_SAVE/
*    Found in b2.neutrals_save.parameters.
*    Read if file is present.
*    Used internally for regulation and rescaling of volume 
*    recombination sources and restarts. Superceded by use of
*    'eirene_dpc_fix' switch. See code in eirene_mc and
*    description in b2cdci for details. 
*  SAVED_VOLREC - real*8 array of size (NSTRAT). Contains the last value of the
*     strength of volume recombination sources from stratum (istra). 
*     Defaults to 0.0.
*
*  NAMELIST /NUMERICS/ 
*    Found in b2.numerics.parameters.
*    Used to specifiy time multipliers per equation and per region
*    for the B2 solver.
*    Active is 'b2srdt_numerics_namelist' is set to 1 or greater.
*    It is recommended to use the same multiplier for all equations within 
*    one region.
*  DTCO - real*8 array of size (0:NS-1,0:NREG). Multiplier to the time used in 
*     solving the continuity equation of species (is) in region (ireg).
*     Defaults to 1.0.
*  DTMO - real*8 array of size (0:NS-1,0:NREG). Multiplier to the time used in 
*     solving the parallel momentum equation of species (is) in region (ireg).
*     Defaults to 1.0.
*  DTEE - real*8 array of size (0:NREG). Multiplier to the time used in solving
*     the electron heat equation in region (ireg).
*     Defaults to 1.0.
*  DTEI - real*8 array of size (0:NREG). Multiplier to the time used in solving
*     the ion heat equation in region (ireg).
*     Defaults to 1.0.
*  SOLVECO - logical array of size (0:NS-1,0:NREG). Indicates whether the 
*     continuity equation for species (is) is to be solved in region (ireg).
*     Subservient to 'no_solve' switch. Defaults to .true..
*  SOLVEMO - logical array of size (0:NS-1,0:NREG). Indicates whether the 
*     parallel momentum equation for species (is) is to be solved in region (ireg).
*     Subservient to 'no_solve' switch. Defaults to .true..
*  SOLVEPO - logical array of size (0:NREG). Indicates whether the 
*     potential equation is to be solved in region (ireg).
*     Subservient to 'no_solve' switch. Defaults to .true..
*  SOLVEEE - logical array of size (0:NREG). Indicates whether the 
*     electron energy equation is to be solved in region (ireg).
*     Subservient to 'no_solve' switch. Defaults to .true..
*  SOLVEEI - logical array of size (0:NREG). Indicates whether the 
*     ion energy equation is to be solved in region (ireg).
*     Subservient to 'no_solve' switch. Defaults to .true..
*  TIME_FACTOR_REQUIRED - real*8. Minimum time scale of evolution allowed for 
*     all equations. Only active is 'b2srsm_enable' is set to 1. 
*     Defaults to 0.1.
*  CORE_DT_SUPPRESSION - real*8. De-multiplier to the timestep in the core. 
*     Only active if less than 1. Should be larger than 0.
*     Applies fully to the innermost core ring of cells (IY .eq. -1).
*     See CORE_DT_FACTOR for further use. Defaults to 1.0.
*  CORE_DT_FACTOR - real*8. Multiplier to the timestep in the core. 
*     Only active is less than 1. Should be larger than 0. Defaults to 1.0.
*     Multiplies each successive core ring of cells (increasing IY)
*     by CORE_DT_FACTOR, until the local time step multiplier is equal to 1.
*  NUMERICS_FILENAME - character*256. Name of the next file to use for reading 
*     a new /NUMERICS/ namelist. Defaults to 'b2.numerics.namelist'.
*  NUMERICS_TIME_MOD - real*8. When the B2 run simulation time, in seconds,
*     modulo(NUMERICS_TIME_MOD), exceeds NUMERICS_TIME_SWITCH, reads the new 
*     namelist from NUMERICS_FILENAME. Also switches to the new namelist as the
*     ELM count (here time/NUMERICS_TIME_MOD) changes. Defaults to 0.0.
*     Only active if NUMERICS_TIME_MOD is greater than 0.
*  NUMERICS_TIME_SWITCH - real*8. Time (in seconds) within an ELM cycle after 
*     which a new /NUMERICS/ namelist is read. Defaults to 0.0.
*  WRITE_NML_NUM - logical. If .true. (default), writes the content of the
*    namelist to stdout after it has been read.
*
*  NAMELIST /TRANSPORT_MODELS/
*    Found in b2.transport_models_save.parameters.
*    Read if file is present.
*  ETA_HCE_MULT - real*8 array of size (-1:NY). Used by the user specified
*     set_transport_eta transport model. See code for details. Defaults to 1.0.
*
*  NAMELIST /TRANSPORT/ 
*    Found in b2.transport.parameters.
*    Allows user to choose between different anomalous transport models.
*    Read if 'b2mod_transport_namelist' is set to 1.
*    Overwrite data from b2ah.dat and b2mn.dat.
*    Does not apply to fluid neutral species.
*  FLAG_DNA - integer. Transport model flag for density-driven diffusion.
*  FLAG_DPA - integer. Transport model flag for pressure-driven dffusion.
*  FLAG_VLA - integer. Transport model flag for pinch velocity.
*  FLAG_VSA - integer. Transport model flag for viscosity.
*  FLAG_HCI - integer. Transport model flag for ion heat diffusivity.
*  FLAG_HCE - integer. Transport model flag for electron heat diffusivity.
*  FLAG_SIG - integer. Transport model flag for field-driven current radial 
*     conductivity.
*  FLAG_ALF - integer. Transport model flag for temperature-driven current 
*     radial conductivity.
*     All flags follow:
*     FLAG=0: Use model from b2ah.dat and b2mn.dat. Default.
*     FLAG=1: Constant, value set to PARM.
*     FLAG=2: 1/N model, multiplied by PARM.
*     FLAG=3: Bohm model, multiplied by PARM.
*     FLAG=4: flux-scaled model, multiplied by PARM.
*  PARM_DNA - real*8 array of size (0:NS-1). 
*     Parameter for the density-driven particle diffusion coefficient for 
*     species (is).
*  PARM_DPA - real*8 array of size (0:NS-1).
*     Parameter for the pressure-driven particle diffusion coefficient for 
*     species (is).
*  PARM_VLA - real*8 array of size (0:NS-1).
*     Parameter for the anomalous radial pinch velocity for species (is).
*  PARM_VSA - real*8 array of size (0:NS-1).
*     Parameter for the viscosity for species (is).
*  PARM_HCI - real*8 array of size (0:NS-1).
*     Parameter for the heat diffusivity coefficient for species (is).
*  PARM_HCE - real*8. Parameter for the electron heat diffusivity coefficient.
*  PARM_SIG - real*8. Parameter for the anomalous radial field-driven current 
*     conductivity.
*  PARM_ALF - real*8. Parameter for the anomalous radial temperature-driven 
*     current conductivity.
*  TRANSPORT_FILENAME - character*256. Name of the next file to use for 
*     reading a new /TRANSPORT/ namelist. Defaults to 'b2.transport.parameters'.
*  TRANSPORT_TIME_MOD - real*8. When the B2 run simulation time, in seconds,
*     modulo(TRANSPORT_TIME_MOD), exceeds TRANSPORT_TIME_SWITCH, reads the 
*     new namelist from TRANSPORT_FILENAME. Also switches to the new namelist
*     if the ELM count (here time/TRANSPORT_TIME_MOD) changes. Defaults to 0.0.
*     Only active if TRANSPORT_TIME_MOD is greater than 0.
*  TRANSPORT_TIME_SWITCH - real*8. Time (in seconds) within an ELM cycle after 
*     which a new /TRANSPORT/ namelist is read. defaults to 0.0.
*  CFLME - real*8. Multiplier to the electron heat flux limit. Default 
*     inherited from b2mn.dat.
*  CFLMI - real*8. Multiplier to the ion heat flux limit. Default inherited 
*     from b2mn.dat.
*  CFLMV - real*8. Multiplier to the viscous heat flux limit. Default 
*     inherited from b2mn.dat.
*    For all flux limit multipliers, setting them to zero turns off the 
*    corresponding flux limiter.
*  WRITE_NML_TRANSP - logical. If .true. (default), writes the content of the
*    namelist to stdout after it has been read.
*  The following four parameters are only activated if 'b2tqna_user_transport'
*     is set to '6'.
*  VOUT_CNV  - real*8 array of size (0:NS-1).
*    Value of the "blob" convection velocity for species (is) at the outer 
*    grid edge (in m/s). Defaults to 0.0.
*  PW0_CNV  - real*8 array of size (0:NS-1). Defaults to 1.0.
*    Exponent in radial profile of the "blob" velocity for species (is).
*  PW1_CNV - real*8 array of size (0:NS-1). Defaults to 5.0.
*    First exponent in poloidal profile of the "blob" velocity for  
*    species (is).
*  PW2_CNV - real*8 array of size (0:NS-1). Defaults to 2.0.
*    Second exponent in a poloidal profile of the "blob" velocity for  
*    species (is).
*
*  NAMELIST /USER/
*    Found in b2.user.parameters.
*    Contains additional data for user-specific diagnostics. 
*    This version provided by Andrei Kukushkin for ITER cases.
*  LHETRGTS - integer array of size (NLIM). List of surface indices (EIRENE 
*              notation) which are used for calculation of helium enrichment.
*              Defaults to -2 for the firt element, -3 for the second element
*              and 0 otherwise.
*  LPFRB_I - integer. B2 x-cell indices corresponding to the bottom of the
*              bypass between the inner and outer divertor (counted
*              from the bottom, if non-positive) (inner side). Defaults to -2.
*  LPFRB_O - integer. B2 x-cell indices corresponding to the bottom of the
*              bypass between the inner and outer divertor (counted
*              from the bottom, if non-positive) (outer side). Defaults to -2.
*  LPFRT_I - integer. B2 x-cell indices corresponding to the top of the
*              bypass between the inner and outer divertor (counted
*              from the x-point, if non-positive) (inner side). Defaults to -2.
*  LPFRT_O - integer. B2 x-cell indices corresponding to the top of the
*              bypass between the inner and outer divertor (counted
*              from the x-point, if non-positive) (outer side). Defaults to -2.
*  J_HE_AT - integer. Species index of the helium atoms in Eirene. The code
*              attempts to find a match by default.
*  J_NE_AT - integer. Species index of the neon atoms in Eirene. The code
*              attempts to find a match by default.
*  J_H_AT - integer. Species index of the hydrogen atoms in Eirene. The code
*              attempts to find a match by default.
*  L_H_MOL - integer. Species index of the hydrogen molecules in Eirene. The code
*              attempts to find a match by default.
*  FUSION_POWER - real*8. Fusion power occuring in core (including neutrons, 
*              in Megawatts). Defaults to 0.0.
*  SPMP_HE_TO_D - real*8. Ratio of He to DT pumping speeds (typically, 0.8).
*              Defaults to 1.0.
*  LPFRS_PMP - integer. Location of the pump.
*              0 no pump at all (default), 1 - lower PFR, 2 - lower outer
*  NPFRGRP - integer. Actual number of surface groups for PFR flows. Defaults to 0.
*  LPFRGRP - integer array of size (NLIM). List of surface segments for PFR 
*              flows.
*  IPFRGRP - integer array of size (NPFRGRP). First positions in this list 
*              for each group. Defaults to 0.
*  JPFRGRP - integer array of size (NPFRGRP). Last positions in this list 
*              for each group. Defaults to 0.
*  GPFRGRP - character*8 array of size (NPFRGRP). Labels for surface groups 
*              for PFR flows. Defaults to ' '.
*  NNTRGRP - integer. Actual number of surface groups for neutral data.
*  LNTRGRP - integer array of size (NLIM). List of surface segments for 
*              neutral data.
*  INTRGRP - integer array of size (NNTRGRP). First positions in this list 
*              for each group.
*  JNTRGRP - integer array of size (NNTRGRP). Last positions in this list 
*              for each group.
*  GNTRGRP - character*8 array of size (NNTRGRP). Labels for the neutral data 
*              surface groups.
*  SPMP_NOM - real*8. Nominal pumping speed. Defaults to 0.0.
*  USER_FILENAME - character*80. Filename where /USER/ namelist is stored. 
*              Defaults to 'b2.user.parameters'.
*
*  NAMELIST /SPUTTER_SAVE/ 
*    Found in b2.sputter_save.parameters.
*    Read if file is present.
*    Used to regulate the temporal evolution of sputtering and restarts.
*    Not intended for human use.
*  NX_SP - integer. Array dimension NX used in this namelist. Defaults to NX.
*  NY_SP - integer. Array dimension NY used in this namelist. Defaults to NY.
*  NS_SP - integer. Array dimension NS used in this namelist. Defaults to NS.
*  SPUTTER_YIELD - real*8 array of size (-1:NX,-1:NY,0:NS-1,1:2). Contains
*     the physical sputtering yield of species (is) at position (ix,iy) in 
*     element (:,:,:,1) and the chemical sputtering yield in element (:,:,:,2).
*     Defaults to 0.0.
*  SPUTTER_YIELD2 - real*8 array of size (-1:NX,-1:NY,0:NS-1,1:2). Contains
*     the energy physical sputtering yield of species (is) at position (ix,iy) 
*     in element (:,:,:,1) and the energy chemical sputtering yield in element 
*     (:,:,:,2). Defaults to 0.0.
*
*  NAMELIST /ATOMIC_PHYSICS_RESCALE/
*    Found in b2.atomic_physics_rescale.parameters.
*    Read if 'b2mndr_atomic_physics_rescale' is set to 1.
*    Contains rescaling multipliers for atomic physics rates.
*  RESCALE_SA - real*8 array of size(0:NS-1). Scaling factors for rtlsa:
*    ionisation rates of species (is). Defaults to 1.0.
*  RESCALE_RA - real*8 array of size(0:NS-1). Scaling factors for rtlra:
*    recombination rates of species (is). Defaults to 1.0.
*  RESCALE_QA - real*8 array of size(0:NS-1). Scaling factors for rtlqa:
*    electron cooling rates of species (is). Defaults to 1.0.
*  RESCALE_CX - real*8 array of size(0:NS-1). Scaling factors for rtlcx:
*    charge exchange rates of species (is). Defaults to 1.0.
*  RESCALE_RD - real*8 array of size(0:NS-1). Scaling factors for rtlrd:
*    line radiation rates of species (is). Defaults to 1.0.
*  RESCALE_BR - real*8 array of size(0:NS-1). Scaling factors for rtlbr:
*    bremsstrahlung radiation rates of species (is). Defaults to 1.0.
*
*-----------------------------------------------------------------------
*.end b2cdcn

      end subroutine b2cdcn
