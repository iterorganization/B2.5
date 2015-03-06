c VERSION: 01.07.2000 23:24
c
c*** Control parameters for transfer from B2 to Eirene
c*** This allows to input some of more often varied Eirene parameters
c*** to be specified in the B2 input file in a convenient way.
c*** 
c*** DELTAT and EIRENE_STEP_CPU are left in the dummy variable list in
c*** the eirsrt routine just to minimize the number of possible errors.
c***>>  Control parameters:
c*** DEF_NSPCSRFG:  max. number of groups of surfaces for which the data
c***                can be transferred (e. g., pump duct, e.t.c)
c*** bn_spcsrf:  actual number of such groups as input from B2
c*** bl_spcsrf:  list of surface segments (non-default standard surfaces
c***            or additional surfaces in Eirene notation) included
c***            in the groups. Negative numbers correspond to NDSS.
c*** bi_spcsrf:  starting position of each group in the list
c*** bj_spcsrf:  ending position of each group in the list
c***>>  Eirene parameters for the groups:
c*** bsps_absr:  particle absorption on the surface
c*** bsps_trno:  surface transparency in positive direction
c*** bsps_trni:  surface transparency in negative direction
c*** bsps_mtri:  surface material in Eirene notation (e.g., 1206 for C)
c*** bsps_tmpr:  surface temperature
c*** bsps_spph:  fudge-factor for physical sputtering
c*** bsps_spch:  fudge-factor for chemical sputtering
c*** bsps_sgrp:  surface group for chemical sputtering
c*** bsps_mtrl:  surface material in human notation (e.g., C)
c*** bsps_id:    identifier of the group of surfaces (assigned in DG)
     -  bsps_absr(DEF_NSPCSRFG), bsps_trno(DEF_NSPCSRFG), 
     ,  bsps_trni(DEF_NSPCSRFG), bsps_mtri(DEF_NSPCSRFG),
     ,  bsps_tmpr(DEF_NSPCSRFG), bsps_spph(DEF_NSPCSRFG),
     ,  bsps_spch(DEF_NSPCSRFG)
      integer nn_spcsrf
      parameter(nn_spcsrf=DEF_NLIM+DEF_NSTS)
      integer bn_spcsrf, bl_spcsrf(nn_spcsrf), bi_spcsrf(DEF_NSPCSRFG),
     ,  bj_spcsrf(DEF_NSPCSRFG), bsps_sgrp(DEF_NSPCSRFG)
      character*8 bsps_mtrl(DEF_NSPCSRFG), bsps_id(DEF_NSPCSRFG)
      common/spcsrfcom/
     r  bsps_absr,bsps_trno,bsps_trni,bsps_mtri,bsps_tmpr,
     r  bsps_spph,bsps_spch,
     i  bn_spcsrf,bl_spcsrf,bi_spcsrf,bj_spcsrf,bsps_sgrp
      common/spcsrfcch/bsps_mtrl,bsps_id
