c  version : 27.07.2000 22:52
     .	dab2(0:ndxp,0:ndyp,natm,1),    dmb2(0:ndxp,0:ndyp,nmol,1),
     .	dib2(0:ndxp,0:ndyp,nion,1),    tab2(0:ndxp,0:ndyp,natm,1),
     .	tmb2(0:ndxp,0:ndyp,nmol,1),    tib2(0:ndxp,0:ndyp,nion,1),
     .	rfluxa(0:ndxp,0:ndyp,natm,1),  rfluxm(0:ndxp,0:ndyp,nmol,1),
     .	refluxa(0:ndxp,0:ndyp,natm,1), refluxm(0:ndxp,0:ndyp,nmol,1),
     .	pfluxa(0:ndxp,0:ndyp,natm,1),  pfluxm(0:ndxp,0:ndyp,nmol,1),
     .	pefluxa(0:ndxp,0:ndyp,natm,1), pefluxm(0:ndxp,0:ndyp,nmol,1),
     .	emiss(0:ndxp,0:ndyp,1,1),      emissmol(0:ndxp,0:ndyp,1,1),
     .	srcml(0:ndxp,0:ndyp,nmol,1),   
     .  edissml(0:ndxp,0:ndyp,nmol,0:nstra),
     .	wldnek(nlimps,0:nstra),        wldnep(nlimps,0:nstra),
     .	wldna(nlimps,natm,0:nstra),    ewlda(nlimps,natm,0:nstra),
     .	wldnm(nlimps,nmol,0:nstra),    ewldm(nlimps,nmol,0:nstra),
     .	wldra(nlimps,natm,0:nstra),    wldrm(nlimps,nmol,0:nstra),
     .	wldpp(nlimps,npls,0:nstra),    wldpa(nlimps,natm,0:nstra),
     .	wldpm(nlimps,nmol,0:nstra),    wldpeb(nlimps,0:nstra),
     .	wldspt(nlimps,0:nstra),
     .  wldspta(nlimps,natm,0:nstra),  wldsptm(nlimps,nmol,0:nstra),
     .  eneutrad(0:ndxp,0:ndyp,natm,0:nstra),
     .  wlarea(nlimps),                wlabsrp(nspz,nlimps),
     .  wldna_res(natm,nsgmx),     !incident atom flux
     .  wldnm_res(nmol,nsgmx),     !incident molecule flux
     .  ewlda_res(natm,nsgmx),     !incident energy flux of atoms
     .  ewldm_res(nmol,nsgmx),     !incident energy flux of molecules
     .  ewldea_res(natm,nsgmx),    !emitted energy flux of atoms
     .  ewldem_res(nmol,nsgmx),    !net incident energy flux due of molecules
     .  ewldrp_res(nsgmx),         !kinetic energy of reflected neutrals originated from ions
     .  ewldmr_res(nmol,nsgmx),    !energy due to recombination of atoms into molecules
     .  wldspt_res(nsgmx),         !flux of sputtered wall material
     .  wldspta_res(natm,nsgmx),   !sputtered flux for each type of emitted atom
     .  wldsptm_res(nmol,nsgmx),   !sputtered flux for each type of emitted molecule
     .  wlpump(nspz,nlimps),       ! pumped flux for each type of pumped particle
     .  wlpump_res(nspz,nsgmx),    ! pumped flux for each type of pumped particle
     .  sarea_res(nsgmx)           !area of the surface elements

CVK NEUTRAL FLUXES, SPATIALLY RESOLVED ON NON-DEFAULT STANDARD SURFACES (NDS)
C First dimension: index of atom or molecule. Second dimension is the index of
C surface elements, controlled by array 'eirdiag_nds_ind'
C

! See WNEUTRAL_FLUXES for the description
      INTEGER ::
     .                                  eirdiag_nds_ind(nsts+1),
     .                                  eirdiag_nds_typ(nsts),
     .                                  eirdiag_nds_srf(nsts),
     .                                  eirdiag_nds_start(nsts),
     .                                  eirdiag_nds_end(nsts)

C Meaning of the arrays which are written to the fort.44 file
C 
C  eirdiag_nds_typ    : type of the surface
C                       1: poloidal surface (Y=CONST)
C                       2: radial surface (X=CONST) 
C  eirdiag_nds_srf    : index of the corresponding radial or poloidal SURFACE
C                       on which the distribution is given 
C  eirdiag_nds_start  : index of the first CELL on the surface '..._srf' 
C  eirdiag_nds_end    : ... of the last CELL ...
C 
C  _sfr, _start and _end are given in B2 indexing. 
C  '_srf' corresponds to the indexing of B2 arrays with fluxes (CELL-1 for outer wall and target).
C
C  eirdiag_nds_ind(IS) : is the index AFTER WHICH the data for surface IS is started in arrays '*_res'.
C                        IF _ind(IS)<0 then this surface is skipped
C      
C
      CHARACTER(8) :: eirtxt(nspz) !ak+vk

      integer nnlimi,nnstsi,nnatmi,nnmoli,nnioni,nnplsi,isrftype(nlimps)
      integer lkindp(npls), lkindm(nmol), lkindi(nion)
c
      common /eirdiag/ dab2,dmb2,dib2,tab2,tmb2,tib2,
     1 rfluxa,rfluxm,refluxa,refluxm,
     2 pfluxa,pfluxm,pefluxa,pefluxm,
     3 emiss,emissmol,eneutrad,
     4 srcml,edissml,wldnek,wldnep,wldna,ewlda,wldnm,ewldm,wldra,wldrm,
     5 wldpp,wldpa,wldpm,wldpeb,wldspt,wldspta,wldsptm,
     6 wlarea, wlabsrp, wlpump,
     7 wldna_res, wldnm_res, ewlda_res, ewldm_res, wldspt_res,
     8 ewldea_res, ewldem_res, ewldrp_res, ewldmr_res, sarea_res,
     9 wldspta_res, wldsptm_res, wlpump_res,
     i nnlimi,nnstsi,nnatmi,nnmoli,nnioni,nnplsi,isrftype,
     2 lkindp,lkindm,lkindi,
     3 eirdiag_nds_ind,eirdiag_nds_typ,eirdiag_nds_srf,
     4 eirdiag_nds_start,eirdiag_nds_end,
     c eirtxt
c
#ifndef NO_SAVE
      save /eirdiag/
#endif
c
c*** Volume data:
c***	srcml	:   power loss due to molecules, including
c***	edissml :   power loss due to molecule dissociation
c*** Surface data:
c***	wldnek	:   heat transferred with neutrals
c***	wldnep	:   potential energy released by neutrals
c***	wldna	:   flux of atoms impinging onto the surface
c***	ewlda	:   their average energy
c***	wldnm	:   flux of molecules impinging onto the surface
c***	ewldm	:   their average energy
c***	wldra	:   flux of reflected atoms
c***	wldrm	:   flux of reflected molecules
c***	wldpp	:   flux of plasma ions impinging onto the surface
c***	wldpa	:   flux of resulting atoms
c***	wldpm	:   flux of resulting molecules
c***	wldpeb	:   power carried away by these atoms and molecules
c***	wldspt	:   flux of sputtered wall material
c***	isrftype:   surface type (iliin in Eirene)
c***    eneutrad:   power radiated due to neutrals (only atoms at the moment)
