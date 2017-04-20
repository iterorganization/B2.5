c  version : 19.02.97 20:49
c
C     EIRENE TO BRAAMS COMMON
c
      COMMON /EIRBRA/
     /    SNI(0:NDXP,0:NDYP,NFL,NSTRA), 
     ,    SMO(0:NDXP,0:NDYP,NFL,NSTRA),
     ,    SEE(0:NDXP,0:NDYP,NSTRA),
     ,    SEI(0:NDXP,0:NDYP,NSTRA),
cdjm Jan2017
     ,    SNI_PAPL(0:NDXP,0:NDYP,NFL,NSTRA),
     ,    SNI_PMPL(0:NDXP,0:NDYP,NFL,NSTRA),
     ,    SNI_PIPL(0:NDXP,0:NDYP,NFL,NSTRA),
     ,    SNI_PPPL(0:NDXP,0:NDYP,NFL,NSTRA),
     ,    SMO_MAPL(0:NDXP,0:NDYP,NFL,NSTRA),
     ,    SMO_MMPL(0:NDXP,0:NDYP,NFL,NSTRA),
     ,    SMO_MIPL(0:NDXP,0:NDYP,NFL,NSTRA),
     ,    SMO_CPPV(0:NDXP,0:NDYP,NFL,NSTRA),
     ,    SEE_EAEL(0:NDXP,0:NDYP,NSTRA),
     ,    SEE_EMEL(0:NDXP,0:NDYP,NSTRA),
     ,    SEE_EIEL(0:NDXP,0:NDYP,NSTRA),
     ,    SEE_EPEL(0:NDXP,0:NDYP,NSTRA),
     ,    SEI_EAPL(0:NDXP,0:NDYP,NSTRA),
     ,    SEI_EMPL(0:NDXP,0:NDYP,NSTRA),
     ,    SEI_EIPL(0:NDXP,0:NDYP,NSTRA),
     ,    SEI_EPPL(0:NDXP,0:NDYP,NSTRA),
c-dpc     ,    DSNIDNI(0:NDXP,0:NDYP,NFL,NSTRA),
c-dpc     ,    DSMODUP(0:NDXP,0:NDYP,NFL,NSTRA),
c-dpc     ,    DSEEDTE(0:NDXP,0:NDYP,NSTRA), DSEIDTE(0:NDXP,0:NDYP,NSTRA),
c-dpc     ,    DSEEDTI(0:NDXP,0:NDYP,NSTRA), DSEIDTI(0:NDXP,0:NDYP,NSTRA),
     ,    FLUXH(10,NSTRA),  srcstrn(nstra),
     ,    VOLSUMN(nstra), VOLSUMM(nstra),
     ,    VOLSUMEI(nstra), VOLSUMEE(nstra),
     ,    flxspci(npls,nstra),srccrfc(natm,nstra)
c
c*** fluxh  : not used
c*** volsumn: plasma particle source due to volumetric recombination
c*** volsum[m,ee,ei]: the same for momentum and energy
c*** srcstrn: intensities of different neutral sources (fluxt)
c*** flxspci: fluxes of plasma ions to the recycling surfaces
c*** srccrfc: source correction factors (from infcop)
C
#ifndef NO_SAVE
      save /EIRBRA/                      ! dpc
#endif
