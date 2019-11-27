      module b2mod_neutrals_namelist

c  version : 27.07.2002 00:00

cank-20020724 {
c* Parameters controlling the source scaling procedure are added:
c*  b2eatcr(i)    : Eirene atomic species corresponding to B2 species i
c*  b2espcr(i)    : Nuclear sequence index corresponding to B2 species i
c*  eb2atcr(i)    : B2 atomic species corresponding to Eirene atomic species i
c*  eb2spcr(i)    : First B2 species corresponding to sequence i
c*  lstrascl(i,j) : scale atomic species j in stratum i with this atomic
c*                  species (0 means no scaling). Species 0 is reserved
c*                  for electron energy (TBD).
c*  gpfc(i,j)     : fraction of atoms of kind i (Eirene numbering) in one
c*                  particle puffed in stratum number j
c*  latmscl(i)    : B2 sequence number with which atoms of kind i scale
c*  lmolscl(i)    : atomic species with which molecules of kind i scale
c*  lionscl(i)    : atomic species with which test ions of kind i scale
c*  mlcmp(i,j)    : number of atoms i in a molecule j
c*  lcns          : Eirene 'non-default surfaces' for the core boundary
c*  ncns          : number of these surfaces
c*  ltns          : Eirene 'non-default surfaces' for the targets
c*  ntns          : number of these surfaces
c*  l_neutrad     : if >=0, import radiation losses due to neutrals from Eirene
c*  l_neutflux    : if >=0, then correct treatment of the incident fluxes in
c*                  B2 and b2plot; if <0, then old (approximate) treatment
c*  dbg_eir_mc    : control of debugging output for eirene_mc
c*  neut_scl_lim  : maximum variation of neutral sources allowed during
c*                  source correction in eirene_mc
cank }
      use b2mod_types
      use b2mod_wall
      use b2mod_b2cmfs
      use b2mod_version
      use b2us_map
      use b2us_geo
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_extrab25
      use eirmod_parmmod
     &    , only: nlim, nsts, nlimps
#else
      use parmmod   ! IGNORE
      use braeir    ! IGNORE
      use eirdiag   ! IGNORE
      use ceirsrt   ! IGNORE
#endif
#else
      use b2mod_eirdiag
#endif
      implicit none
#     include "DIMENSIONS.F"
#ifndef B25_EIRENE
#     include "COUPLE/KOPPLDIM.F"
#endif
#ifndef SOLPS4_3
#     include "COUPLE/EIRSRT.F"
#endif
      integer nstraid, nsdmax, natmid, nmolid, nionid, nstsid, nlimpsd
      parameter (NSTRAID=DEF_NSTRA, NSDMAX=DEF_NSD)
      parameter (natmid=DEF_NATM, nmolid=DEF_NMOL, nionid=DEF_NION)
      parameter (nstsid=DEF_NSTS)
      parameter (nlimpsd=DEF_NLIM+DEF_NSTS)
#ifdef ALLOCATE_AND_NAMELIST
      real (kind=R8), allocatable, save ::
     *   RECYC(:,:), ERECYC(:,:), MRECYC(:,:), RCION(:,:), recyceir(:),
     *   USERFLUXPARM(:,:),
     *   volrecstart(:), saved_volrec(:),
     *   phys_sput(:,:), chem_sput(:,:),
     1   gpfc(:,:), gpdata(:,:,:),
     2   RC_FACE_ORI(:,:)
#ifndef NO_B2_CHEM_SPUT
#ifndef SOLPS4_3
      real (kind=R8), allocatable, save ::
     *   chemical_sputter_yield(:)
#endif
#endif
      integer, allocatable, save ::
     *   RCPOS(:),
     *   RCSTART(:),
     *   RCEND(:),
     *   ARCSTART(:),
     *   ARCEND(:),
     *   RC_LIST_SIZE(:),
     *   RC_LIST_X(:,:),
     *   RC_LIST_Y(:,:),
     *   NRCFACES(:),
     *   RC_FACES(:,:),
     *   TARGSP(:,:),
     *   species_start(:), species_end(:),
     *   b2eatcr(:), b2espcr(:), eb2atcr(:), eb2spcr(:),
     *   lstrascl(:,:), latmscl(:), lmolscl(:), lionscl(:),
     1   mlcmp(:,:), lcns(:), ltns(:), msns(:,:), ksns(:), lsns(:),
     2     ngpdata(:)
      real (kind=R8), allocatable, save ::
     *   RCFE(:), RCFI(:)      
      logical, allocatable, save ::
     *   CHEMSP(:), TIME_DEP_PUFF(:)
      character*1, allocatable, save ::
     *   CRCSTRA(:),
     *   ARCSTRA(:),
     *   RC_LIST_CHAR(:,:)
#else
      real (kind=R8), save :: RECYC(0:nsdmax-1,NSTRAID)
      real (kind=R8), save :: ERECYC(0:nsdmax-1,NSTRAID)
      real (kind=R8), save :: MRECYC(0:nsdmax-1,NSTRAID)
      real (kind=R8), save :: RCION(0:nsdmax-1,NSTRAID)
      real (kind=R8), save :: recyceir(NSTRAID)
      real (kind=R8), save :: USERFLUXPARM(NSTRAID,2)
      real (kind=R8), save :: volrecstart(NSTRAID)
      real (kind=R8), save :: saved_volrec(NSTRAID)
      real (kind=R8), save :: chem_sput(0:nsdmax-1,NSTRAID)
      real (kind=R8), save :: phys_sput(0:nsdmax-1,NSTRAID)
      real (kind=R8), save :: gpdata(100,2,NSTRAID)
      real (kind=R8), save :: RC_FACE_ORI(2*(DEF_NXD+DEF_NYD),NSTRAID)
#ifndef NO_B2_CHEM_SPUT
#ifndef SOLPS4_3
      real (kind=R8), save ::
     ,  chemical_sputter_yield(0:DEF_NLIM+DEF_NSTS)
#endif
#endif
      integer, save :: RCPOS(NSTRAID)
      integer, save :: RCSTART(NSTRAID)
      integer, save :: RCEND(NSTRAID)
      integer, save :: ARCSTART(NSTRAID)
      integer, save :: ARCEND(NSTRAID)
      integer, save :: RC_LIST_SIZE(NSTRAID)
      integer, save :: RC_LIST_X(2*(DEF_NXD+DEF_NYD),NSTRAID)
      integer, save :: RC_LIST_Y(2*(DEF_NXD+DEF_NYD),NSTRAID)
      character, save :: RC_LIST_CHAR(2*(DEF_NXD+DEF_NYD),NSTRAID)
      integer, save :: NRCFACES(NSTRAID)
      integer, save :: RC_FACES(2*(DEF_NXD+DEF_NYD),NSTRAID)
      integer, save :: TARGSP(NSTRAID,ntrack)
      logical, save :: CHEMSP(NSTRAID)
      logical, save :: TIME_DEP_PUFF(NSTRAID)
      integer, save :: ngpdata(NSTRAID)
      integer, save :: species_start(NSTRAID), species_end(NSTRAID)
      integer, save :: b2eatcr(0:nsdmax-1)
      integer, save :: b2espcr(0:nsdmax-1)
      integer, save :: eb2atcr(1:natmid)
      integer, save :: eb2spcr(1:natmid)
      integer, save :: lstrascl(NSTRAID,0:natmid)
      character*1, save :: CRCSTRA(NSTRAID)
      character*1, save :: ARCSTRA(NSTRAID)
      integer, save :: mlcmp(natmid,nmolid),
     ,  latmscl(natmid), lmolscl(nmolid), lionscl(nionid),
     ,  lcns(nstsid), ltns(nstsid),
     ,  lsns(nstraid*DEF_NSRFS), ksns(NSTRAID), msns(2,NSTRAID)
      real (kind=R8), save :: gpfc(natmid,nstraid)
      real (kind=R8), save :: RCFE(NSTRAID), RCFI(NSTRAID)      
#endif

      integer, save :: NSTRAI, nstrat, eirene_mod, ncns, ntns,
     .     dbg_eir_mc, debug_flags(100), neutrals_count
      real (kind=R8), save :: RF_NEUT(4)
      real (kind=R8), save :: eirene_step_cpu, eirene_step_dt,
     .     volrecinc, volrecwt
      character*256, save :: neutrals_filename
      real (kind=R8), save :: neutrals_time_mod
      real (kind=R8), save :: neutrals_time_switch
      real (kind=R8), save :: neut_scl_lim

#ifndef B25_EIRENE
#ifndef NO_B2_CHEM_SPUT
      real (kind=R8), save :: fchar_chemical
      integer, save :: igass_chemical, itsput_chemical, issput_chemical
#endif
#endif

      real(kind=R8), save ::
     1 chemical_erosion_redep_fac,
     2 chemical_erosion_be_fac_a,
     3 chemical_erosion_be_fac_b,
     4 chemical_erosion_be_fac_c
      logical, save ::
     1 chemical_erosion_be_fac

      namelist /neutrals_save/ saved_volrec

      private :: neutrals
      namelist /neutrals/ nstrai, rcpos, rcstart, rcend, targsp, chemsp,
     .  rcfe, rcfi,
     1  rc_list_size, rc_list_x, rc_list_y, rc_list_char,
     4  recyc, erecyc, mrecyc, rcion, recyceir, userfluxparm, crcstra,
     .  rf_neut, phys_sput, chem_sput,
#ifndef NO_B2_CHEM_SPUT
     .  chemical_sputter_yield, fchar_chemical, igass_chemical,
     .  itsput_chemical, issput_chemical,
#endif
     .  eirene_step_cpu, eirene_step_dt, eirene_mod,
     .  volrecstart, volrecinc, volrecwt,
     .  species_start, species_end,
     5  neutrals_filename, neutrals_time_mod, neutrals_time_switch,
     ,  lstrascl, l_neutrad, l_neutflux,
     ,  b2eatcr, b2espcr, eb2spcr, eb2atcr,
     ,  latmscl, lmolscl, mlcmp, lionscl,
     ,  lcns, ltns, lsns, ksns, gpfc, dbg_eir_mc, debug_flags,
     ,  neut_scl_lim, track_index, track_chem_sput,
     *  chemical_erosion_redep_fac,
     *  chemical_erosion_be_fac,
     *  chemical_erosion_be_fac_a,
     *  chemical_erosion_be_fac_b,
     *  chemical_erosion_be_fac_c,
     6  n_spcsrf, l_spcsrf, sps_id, i_spcsrf, j_spcsrf,
     7  sps_absr, sps_trno, sps_trni, sps_mtri, sps_tmpr,
     8  sps_spph, sps_spch, sps_sgrp, sps_mtrl,
     .  msns, write_nml_neut, time_dep_puff, ngpdata, gpdata
      logical, save :: write_nml_neut = .true.

csw 21feb2012
      integer, save ::
     . l_neutrad=0, l_neutflux=-1

      interface read_neutrals_namelist
        module procedure read_neutrals_namelist_st
        module procedure read_neutrals_namelist_us
      end interface

      contains

      subroutine alloc_b2mod_neutrals(nxd,nyd,nsd,
     &  natm_in,nmol_in,nion_in,nstra_in,nsts_in,sput_dst)
      use b2mod_elements
     & ,only : is_codes
      use b2mod_b2cmpa
      implicit none
      integer, intent(in) :: nxd, nyd, nsd,
     &  natm_in, nmol_in, nion_in, nstra_in, nsts_in, sput_dst
      integer is, isstart, isend
      integer is0, ispecies
      character*12 hlp_frm
      intrinsic nint

#ifdef ALLOCATE_AND_NAMELIST
      NSTRAID=NSTRAID_in
      NSTSID=nsts_in
      NATMID=natm_in
      allocate(RECYC(0:nsd-1,NSTRAID))
      allocate(ERECYC(0:nsd-1,NSTRAID))
      allocate(MRECYC(0:nsd-1,NSTRAID))
      allocate(RCION(0:nsd-1,NSTRAID))
      allocate(recyceir(NSTRAID))
      allocate(USERFLUXPARM(NSTRAID,2))
      allocate(volrecstart(NSTRAID))
      allocate(chem_sput(0:nsd-1,NSTRAID))
      allocate(phys_sput(0:nsd-1,NSTRAID))
      allocate(RCPOS(NSTRAID))
      allocate(RCSTART(NSTRAID))
      allocate(RCEND(NSTRAID))
      allocate(ARCSTART(NSTRAID))
      allocate(ARCEND(NSTRAID))
      allocate(RCFE(NSTRAID))
      allocate(RCFI(NSTRAID))
      allocate(RC_LIST_SIZE(NSTRAID))
      allocate(RC_LIST_X(2*(NXD+NYD),NSTRAID))
      allocate(RC_LIST_Y(2*(NXD+NYD),NSTRAID))
      allocate(RC_LIST_CHAR(2*(NXD+NYD),NSTRAID))
      allocate(NRCFACES(NSTRAID))
      allocate(RC_FACES(2*(NXD+NYD),NSTRAID))
      allocate(RC_FACE_ORI(2*(NXD+NYD),NSTRAID))
      allocate(TARGSP(NSTRAID,ntrack))
      allocate(CHEMSP(NSTRAID))
      allocate(TIME_DEP_PUFF(NSTRAID))
      allocate(NGPDATA(NSTRAID))
      allocate(GPDATA(100,2,NSTRAID))
      allocate(CRCSTRA(NSTRAID))
      allocate(ARCSTRA(NSTRAID))
      allocate(species_start(NSTRAID))
      allocate(species_end(NSTRAID))
      allocate(b2eatcr(0:nsd-1))
      allocate(b2espcr(0:nsd-1))
      allocate(eb2atcr(natm_in))
      allocate(eb2spcr(nspecies))
      allocate(lstrascl(nstra_in,0:natm_in))
      allocate(latmscl(natm_in))
      allocate(lmolscl(nmol_in))
      allocate(lionscl(nion_in))
      allocate(gpfc(natm_in,nstra_in))
      allocate(mlcmp(natm_in,nmol_in))
      allocate(lcns(nsts_in))
      allocate(ltns(nsts_in))
      allocate(lsns(nstra_in*DEF_NSRFS))
      allocate(ksns(nstra_in))
      allocate(msns(2,nstra_in))
#ifndef NO_B2_CHEM_SPUT
      allocate(chemical_sputter_yield(0:nlimps))
#endif
#else
#ifdef WG_TODO
      call xertst (nxd.le.DEF_NXD, 'faulty parameter nxd')
      call xertst (nyd.le.DEF_NYD, 'faulty parameter nyd')
      call xertst (nsd.le.DEF_NSD, 'faulty parameter nsd')
#endif
      write (*,*) ' natm_id,natmid ',natm_in, natmid
      call xertst (natm_in.le.NATMID, 'faulty parameter natm_in')
      call xertst (nmol_in.le.NMOLID, 'faulty parameter nmol_in')
      call xertst (nion_in.le.NIONID, 'faulty parameter nion_in')
      call xertst (nsts_in.le.NSTSID, 'faulty parameter nsts_in')
      call xertst (nstra_in.le.NSTRAID, 'faulty parameter nstraid_in')
#endif
      call xertst (natm_in.ge.nspecies, 'faulty parameter natm_in')
      CRCSTRA=' '
      ARCSTRA=' '
      SPS_ID=' '
      SPS_MTRL=' '
      NSTRAI=0
      RCPOS=-2
      RCSTART=-2
      RCEND=-2
      ARCSTART=-2
      ARCEND=-2
      RCFE = 0.5_R8
      RCFI = 3.0_R8
      RC_LIST_SIZE=0
      RC_LIST_X=-2
      RC_LIST_Y=-2
      RC_LIST_CHAR=' '
      TARGSP=-1
      TARGSP(:,1)=sput_dst
      CHEMSP=.false.
      TIME_DEP_PUFF=.false.
      N_SPCSRF=0
      L_SPCSRF=0
      I_SPCSRF=0
      J_SPCSRF=0
      SPS_ABSR=-1.0_R8
      SPS_TRNO=-1.0_R8
      SPS_TRNI=-1.0_R8
      SPS_MTRI= 0.0_R8
      SPS_TMPR= 1.0e15_R8
      SPS_SPPH=-1.0_R8
      SPS_SPCH=-1.0_R8
      SPS_SGRP=-1
      NGPDATA=0
      GPDATA=0.0_R8
cxpb The default recycling model is 100% particle recycling and
cxpb 30% energy recycling for all ionised species except carbon.
cwdk 0% momentum recycling, for backwards compatibility
      RECYC=0.0_R8
      ERECYC=0.0_R8
      MRECYC=0.0_R8
      do is = 1, nsd-1
        if (.not.is_neutral(is) .and. nint(zn(is)).ne.6) then
          RECYC(is,:)=1.0_R8
          ERECYC(is,:)=0.3_R8
        endif
      enddo
      RCION=0.0_R8
      recyceir=1.0_R8
      USERFLUXPARM=0.0_R8
      RF_NEUT=1.0_R8
      phys_sput=0.0
      chem_sput=0.0
#ifndef NO_B2_CHEM_SPUT
      chemical_sputter_yield=0.0_R8
      fchar_chemical=0
      igass_chemical=0
      itsput_chemical=0
      issput_chemical=0
#endif
      eirene_step_cpu=0.0_R8
      eirene_step_dt=1.0e-3_R8
      eirene_mod=1
      volrecstart(:)=1.0e21_R8
      volrecinc=0.0_R8
      volrecwt=0.1_R8
      species_start(:)=0
      species_end(:)=nsd-1
      neutrals_filename='b2.neutrals.parameters'
      neutrals_time_mod=0.0_R8
      neutrals_time_switch=0.0_R8

      lstrascl=0
      b2eatcr=0
      b2espcr=0
      eb2atcr=-1
      eb2spcr=-1
      latmscl=0
      lmolscl=0
      lionscl=0
      mlcmp=0
      lcns=0
      ltns=0
      ncns=0
      ntns=0
      ksns=0
      lsns=0
      msns=0
      gpfc=0.0_R8
      dbg_eir_mc=100                     !srv 08.06.13
      debug_flags=0
      neut_scl_lim=2

      chemical_erosion_redep_fac=1.0_R8
      chemical_erosion_be_fac=.false.
c the following gives a transition from 1 to 0.1 in the range 0.1 to 0.3
      chemical_erosion_be_fac_a=0.2_R8
      chemical_erosion_be_fac_b=0.05_R8
      chemical_erosion_be_fac_c=0.9_R8

      write(*,*) 'Sput_dst : ',sput_dst
      track_index=0
      track_species='  '
      track_chem_sput=.false.
      if (sput_dst.ge.0 .and. sput_dst.lt.nsd) then
        track_index(sput_dst)=1
        track_species(1)=is_codes(sput_dst)
        if (nint(zn(sput_dst)).eq.6) track_chem_sput(1)=.true.
      endif

      is0 = 0
      ispecies = 0
      do while (is0.lt.nsd)
        is = is0
        ispecies = ispecies + 1
        do while (zn(is).eq.zn(is0).and.
     &     abs(am(is)-am(is0)).lt.1.0e-3_R8.and.is.lt.nsd)
          is = is + 1
          if (is.eq.nsd) exit
        enddo
        isstart = is0
        isend = is - 1
        b2eatcr(isstart:isend)=ispecies
        b2espcr(isstart:isend)=ispecies
        eb2atcr(ispecies)=isstart
        eb2spcr(ispecies)=isstart
        latmscl(ispecies)=ispecies
        is0 = is
      enddo
      write(hlp_frm,'(a6,i3,a3)') '(1x,a,',nsd,'i4)'
      write(*,hlp_frm) 'b2espcr : ',b2espcr(0:nsd-1)
      write(*,hlp_frm) 'b2eatcr : ',b2eatcr(0:nsd-1)

      return
      end subroutine alloc_b2mod_neutrals

      subroutine alloc_b2mod_neutrals_save(NSTRAID_in)
      integer nstraid_in

#ifdef ALLOCATE_AND_NAMELIST
      NSTRAID=NSTRAID_in
      allocate(saved_volrec(NSTRAID))
#else
      call xertst (nstraid_in.le.nstraid, 'faulty parameter nstraid_in')
#endif
      return
      end subroutine alloc_b2mod_neutrals_save

      subroutine dealloc_b2mod_neutrals
      implicit none

#ifdef ALLOCATE_AND_NAMELIST
      if (.not.allocated(recyc)) return
      deallocate(RECYC)
      deallocate(ERECYC)
      deallocate(MRECYC)
      deallocate(RCION)
      deallocate(recyceir)
      deallocate(USERFLUXPARM)
      deallocate(volrecstart)
      deallocate(chem_sput)
      deallocate(phys_sput)
      deallocate(RCPOS)
      deallocate(RCSTART)
      deallocate(RCEND)
      deallocate(ARCSTART)
      deallocate(ARCEND)
      deallocate(RCFE)
      deallocate(RCFI)
      deallocate(RC_LIST_SIZE)
      deallocate(RC_LIST_X)
      deallocate(RC_LIST_Y)
      deallocate(RC_LIST_CHAR)
      deallocate(NRCFACES)
      deallocate(RC_FACES)
      deallocate(RC_FACE_ORI)
      deallocate(TARGSP)
      deallocate(CHEMSP)
      deallocate(TIME_DEP_PUFF)
      deallocate(NGPDATA)
      deallocate(GPDATA)
      deallocate(CRCSTRA)
      deallocate(ARCSTRA)
      deallocate(species_start)
      deallocate(species_end)
      deallocate(b2eatcr)
      deallocate(b2espcr)
      deallocate(eb2atcr)
      deallocate(eb2spcr)
      deallocate(lstrascl)
      deallocate(latmscl)
      deallocate(lmolscl)
      deallocate(lionscl)
      deallocate(gpfc)
      deallocate(mlcmp)
      deallocate(lcns)
      deallocate(ltns)
      deallocate(lsns)
      deallocate(ksns)
      deallocate(msns)
#ifndef NO_B2_CHEM_SPUT
      deallocate(chemical_sputter_yield)
#endif
#endif
      return
      end subroutine dealloc_b2mod_neutrals

      subroutine dealloc_b2mod_neutrals_save
#ifdef ALLOCATE_AND_NAMELIST
      deallocate(saved_volrec)
#endif
      return
      end subroutine dealloc_b2mod_neutrals_save

      subroutine read_neutrals_namelist_st(nx, ny, ns, dotest)
      use b2mod_wall
     & ,only : track_index, ntrack, track_species, track_chem_sput
      use b2mod_indirect
      use b2mod_elements
     & ,only : is_codes
      use b2mod_boundary_namelist
      use b2mod_b2cmpa
      implicit none
      integer nx, ny, ns
      logical dotest
      integer is, ic, ix, istra, i, j, k, n, ii(1), itrack
      logical found, side
      intrinsic mod
      real (kind=R8) hlp
      logical :: file_ok
      logical, allocatable :: old_style(:)
      character*16 hlp_frm
      character*260 filename

      filename=neutrals_filename
      call find_file(filename,file_ok)
      if(file_ok) then
        open(99,file=filename)
        call cfverr(99,neutrals_version)
        write(*,*) 'Reading namelist neutrals from ',trim(filename)
#ifdef USE_MPI
        read(99,neutrals,err=999)
#else
        read(99,neutrals)
#endif
        write(*,*) 'Read namelist neutrals from ',trim(filename)
        close(99)
      else
        call xerrab ('No '//trim(filename))
      endif
      if (nstrai.gt.0) then
cxpb  In Eirene facelift the time stratum is handled separately
cxpb  On the B2 side, use nstrat as nstrai + the time stratum
        if (crcstra(nstrai).eq.'T') then
          nstrat=nstrai
          nstrai=nstrai-1
        elseif (crcstra(nstrai+1).eq.'T') then
          nstrat=nstrai+1
        else
#ifdef B25_EIRENE
          nstrat=nstrai+1
          crcstra(nstrat)='T'
          rcpos(nstrat)=0
          rcstart(nstrat)=0
          rcend(nstrat)=0
#else
          nstrat=nstrai
#endif
        endif
        allocate(old_style(nstrat))
        do is = 1, nstrat
          if (rc_list_size(is).eq.0) then
            if (crcstra(is)=="A") then
                ! Automatic creation of boundary condition list from structures
              old_style(is) = .false.      
              call build_boundary_list( rcstart(is), rcend(is),
     $               rc_list_size(is),
     $               rc_list_x(:,is), rc_list_y(:,is),
     $               rc_list_char(:,is))
            else
              old_style(is) = .true.
              rc_list_size(is) = rcend(is) - rcstart(is) + 1
              rc_list_char(1:rc_list_size(is),is) = crcstra(is)
              if (crcstra(is).eq.'N' .or. crcstra(is).eq.'S') then
                do ix = rcstart(is), rcend(is)
                  rc_list_x(ix-rcstart(is)+1,is) =  ix
                end do
                rc_list_y(1:rc_list_size(is),is) = rcpos(is)
              else if (crcstra(is).eq.'E' .or. crcstra(is).eq.'W') then
                rc_list_x(1:rc_list_size(is),is) = rcpos(is)
                do ix = rcstart(is), rcend(is)
                  rc_list_y(ix-rcstart(is)+1,is) =  ix
                end do
              endif
            endif
          else
            old_style(is) = .false.
          end if
        end do
      else
        nstrat = 0
      endif
c
      if(chemical_erosion_redep_fac.ne.1.0_R8)
     * write(*,*)
     1 'Note: chemical_erosion_redep_fac = ',
     2 chemical_erosion_redep_fac
      if(chemical_erosion_be_fac)
     1 write(*,*) 'Note: suppression of chemical sputtering by Be ',
     2 'is active with parameters ',
     3 chemical_erosion_be_fac_a,
     4 chemical_erosion_be_fac_b,
     5 chemical_erosion_be_fac_c
c
      do istra = 1, nstrat
        do is = 1, nnatmi
          if (lstrascl(istra,is).lt.0) lstrascl(istra,is)=is
          if (lstrascl(istra,is).gt.nnatmi) lstrascl(istra,is)=is
        end do
      end do

*   ..set additional variables
      do is = 1, nstrat
        if (crcstra(is).ne.'N' .and. crcstra(is).ne.'S' .and.
     &      crcstra(is).ne.'W' .and. crcstra(is).ne.'E') targsp(is,:)=-1
      end do

      arcstra = crcstra
      arcstart = rcstart
      arcend = arcend

*   ..test input
      call xertst (0.le.minval(track_index).and.
     1 maxval(track_index).le.ntrack, 'faulty input track_index')
      call xertst (0.le.nstrai, 'faulty input nstrai')
      call xertst (0.le.nstrat, 'faulty input nstrat')
      call xertst (0.0_R8.lt.rf_neut(1).and.0.0_R8.lt.rf_neut(2).and.
     .             0.0_R8.lt.rf_neut(3).and.0.0_R8.lt.rf_neut(4),
     .  'faulty input rf_neut')
#ifndef NO_B2_CHEM_SPUT
      call xertst (0.le.issput_chemical,
     .  'faulty input issput_chemical')
#endif
      call xertst (0.0_R8.le.eirene_step_cpu,
     . 'faulty input eirene_step_cpu')
      call xertst (0.0_R8.lt.eirene_step_dt,
     . 'faulty input eirene_step_dt')
      call xertst (1.0_R8.lt.neut_scl_lim, 'faulty input neut_scl_lim')
      call xertst (0.0_R8.le.volrecinc, 'faulty input volrecinc')
      call xertst (0.0_R8.lt.volrecwt, 'faulty input volrecwt')
      call xertst (0.0_R8.le.neutrals_time_mod,
     . 'faulty input neutrals_time_mod')
      call xertst (0.0_R8.le.neutrals_time_switch,
     . 'faulty input neutrals_time_switch')
      do istra = 1, nstrat
        if (crcstra(istra).eq.'N') then
          if (old_style(istra)) then
            call xertst (1.le.rcpos(istra).and.rcpos(istra).le.ny+1,
     .       'faulty input rcpos (North)')
            call xertst (-1.le.rcstart(istra).and.
     .                   rcstart(istra).le.nx+1,
     .       'faulty input rcstart (North)')
            call xertst (-1.le.rcend(istra).and.rcend(istra).le.nx+1,
     .       'faulty input rcend (North)')
          else
            do ix = 1, rc_list_size(istra)
              call xertst (1.le.rc_list_y(ix,istra).and.
     .                     rc_list_y(ix,istra).le.ny+1,
     .                     'faulty input rc_list_y (North)')
              call xertst (-1.le.rc_list_x(ix,istra).and.
     .                     rc_list_x(ix,istra).le.nx+1,
     .                     'faulty input rc_list_x (North)')
            end do
          end if
        else if (crcstra(istra).eq.'S') then
          if (old_style(istra)) then
            call xertst (-1.le.rcpos(istra).and.rcpos(istra).le.ny-1,
     .       'faulty input rcpos (South)')
            call xertst (-1.le.rcstart(istra).and.
     .                   rcstart(istra).le.nx+1,
     .       'faulty input rcstart (South)')
            call xertst (-1.le.rcend(istra).and.rcend(istra).le.nx+1,
     .       'faulty input rcend (South)')
          else
            do ix = 1, rc_list_size(istra)
              call xertst (-1.le.rc_list_y(ix,istra).and.
     .                     rc_list_y(ix,istra).le.ny-1,
     .                     'faulty input rc_list_y (South)')
              call xertst (-1.le.rc_list_x(ix,istra).and.
     .                     rc_list_x(ix,istra).le.nx+1,
     .                     'faulty input rc_list_x (South)')
            end do
          end if
        else if (crcstra(istra).eq.'E') then
          if (old_style(istra)) then
            call xertst (1.le.rcpos(istra).and.rcpos(istra).le.nx+1,
     .       'faulty input rcpos (East)')
            call xertst (-1.le.rcstart(istra).and.
     .       rcstart(istra).le.ny+1,
     .       'faulty input rcstart (East)')
            call xertst (-1.le.rcend(istra).and.rcend(istra).le.ny+1,
     .       'faulty input rcend (East)')
          else
            do ix = 1, rc_list_size(istra)
              call xertst (1.le.rc_list_x(ix,istra).and.
     .                     rc_list_x(ix,istra).le.nx+1,
     .                     'faulty input rc_list_x (East)')
              call xertst (-1.le.rc_list_y(ix,istra).and.
     .                     rc_list_y(ix,istra).le.ny+1,
     .                     'faulty input rc_list_y (East)')
            end do
          end if
        else if (crcstra(istra).eq.'W') then
          if (old_style(istra)) then
            call xertst (-1.le.rcpos(istra).and.rcpos(istra).le.nx-1,
     .       'faulty input rcpos (West)')
            call xertst (-1.le.rcstart(istra).and.
     .       rcstart(istra).le.ny+1,
     .       'faulty input rcstart (West)')
            call xertst (-1.le.rcend(istra).and.rcend(istra).le.ny+1,
     .       'faulty input rcend (West)')
          else
            do ix = 1, rc_list_size(istra)
              call xertst (-1.le.rc_list_x(ix,istra).and.
     .                     rc_list_x(ix,istra).le.nx-1,
     .                     'faulty input rc_list_x (West)')
              call xertst (-1.le.rc_list_y(ix,istra).and.
     .                     rc_list_y(ix,istra).le.ny+1,
     .                     'faulty input rc_list_y (West)')
            end do
          end if
        endif
        if (old_style(istra))
     .   call xertst (rcstart(istra).le.rcend(istra),
     .   'faulty inputs rcstart > rcend')
        call xertst (maxval(targsp(istra,:)).lt.ns,
     .   'faulty inputs targsp > ns')
        call xertst (targsp(istra,1).ge.-1,
     .   'faulty inputs targsp(,1)')
        call xertst (0.0_R8.le.chemical_erosion_redep_fac,
     .   'faulty input chemical_erosion_redep_fac')
        call xertst (0.0_R8.lt.chemical_erosion_be_fac_a,
     .   'faulty input chemical_erosion_be_fac_a')
        call xertst (0.0_R8.lt.chemical_erosion_be_fac_b,
     .   'faulty input chemical_erosion_be_fac_b')
        call xertst (0.0_R8.le.chemical_erosion_be_fac_c
     &          .and.1.0_R8.ge.chemical_erosion_be_fac_c,
     .   'faulty input chemical_erosion_be_fac_c')
        call xertst (0.0_R8.le.volrecstart(istra),
     .   'faulty input volrecstart')
        call xertst (0.le.species_start(istra).and.
     &   species_start(istra).lt.ns, 'faulty input species_start')
        call xertst (0.le.species_end(istra).and.
     &   species_end(istra).lt.ns, 'faulty input species_end')
        call xertst (species_start(istra).le.species_end(istra),
     .   'faulty inputs species_start > species_end')
        do is = 0, ns-1
          call xertst (0.0_R8.le.recyc(is,istra), 'faulty input recyc')
          call xertst (0.0_R8.le.erecyc(is,istra),
     .     'faulty input erecyc')
          call xertst (0.0_R8.le.mrecyc(is,istra),
     .     'faulty input mrecyc')
          call xertst (0.0_R8.le.rcion(is,istra), 'faulty input rcion')
          call xertst (0.0_R8.le.phys_sput(is,istra),
     .     'faulty input phys_sput')
          call xertst (0.0_R8.le.chem_sput(is,istra),
     .     'faulty input chem_sput')
        enddo
        call xertst (0.0_R8.le.recyceir(istra), 'faulty input recyceir')
        do is = 1, nnatmi
          call xertst(0.0_R8.le.gpfc(is,istra), 'faulty input gpfc')
          call xertst(0.le.lstrascl(istra,is) .and.
     .                lstrascl(istra,is).le.nnatmi,
     .                'faulty input lstrascl')
        enddo
        if (time_dep_puff(istra)) then
          call xertst(crcstra(istra).eq.'C',
     .     'Time_dep_puff set for a non-puff stratum')
          call xertst(ngpdata(istra).gt.0,
     .     'Time_dep_puff data list size not declared')
          call xertst(ngpdata(istra).le.100,
     .     'Time_dep_puff data list size too large')
          do ix = 1, ngpdata(istra)
            call xertst(gpdata(ix,2,istra).ge.0.0_R8,
     .       'Negative gas puff strength declared!')
          end do
        endif
      enddo
      do is = 0, ns-1
        call xertst(1.le.b2espcr(is).and.b2espcr(is).le.nspecies,
     .   'faulty input b2espcr')
        if(nnatmi.eq.0) then
          call xertst(1.le.b2eatcr(is).and.b2eatcr(is).le.nspecies,
     .     'faulty input b2eatcr')
        else
          call xertst(1.le.b2eatcr(is).and.b2eatcr(is).le.nnatmi,
     .     'faulty input b2eatcr')
        endif
      enddo
      do is = 1, nnatmi
        call xertst(0.le.eb2atcr(is).and.eb2atcr(is).lt.ns,
     .   'faulty input eb2atcr')
        call xertst(0.le.latmscl(is), 'faulty input latmscl')
#ifdef B25_EIRENE
        if (latmscl(is).eq.0) then
          write(6,'(a,i3,a)') 'Atom ',is,
     .     ' in Eirene does not scale with any B2.5 sequence'
          write(6,'(a)')
     .     'Please provide a scaling using the latmscl array '//
     .     'in b2.neutral.parameters'
          write(0,'(a,i3,a)') 'Atom ',is,
     .     ' in Eirene does not scale with any B2.5 sequence'
          write(0,'(a)')
     .     'Please provide a scaling using the latmscl array '//
     .     'in b2.neutral.parameters'
        end if
#endif
      enddo
      do is = 1, nspecies
        call xertst(0.le.eb2spcr(is).and.eb2spcr(is).lt.ns,
     .   'faulty input eb2spcr')
      enddo
      do is = 1, nnmoli
        do i = 1, nnatmi
          call xertst(0.le.mlcmp(i,is), 'faulty input mlcmp')
        enddo
        call xertst(0.le.lmolscl(is), 'faulty input lmolscl')
#ifdef B25_EIRENE
        if (lmolscl(is).eq.0) then
          write(6,'(a,i3,a)') 'Molecule ',is,
     .     ' does not scale with any atom species'
          write(6,'(a)')
     .     'Please provide a scaling using the lmolscl array '//
     .     'in b2.neutral.parameters'
          write(0,'(a,i3,a)') 'Molecule ',is,
     .     ' does not scale with any atom species'
          write(0,'(a)')
     .     'Please provide a scaling using the lmolscl array '//
     .     'in b2.neutral.parameters'
        end if
        if (maxval(mlcmp(1:nnatmi,is)).eq.0) then
          write(6,'(a,i3,a)')
     .     'Composition of molecule ',is,' is not given'
          write(6,'(a,i3,a)')
     .     'If it contains species followed in the plasma, '//
     .     'please specify them using the mlcmp(:,',is,
     .     ') array from b2.neutrals.parameters'
          write(0,'(a,i3,a)')
     .     'Composition of molecule ',is,' is not given'
          write(0,'(a,i3,a)')
     .     'If it contains species followed in the plasma, '//
     .     'please specify them using the mlcmp(:,',is,
     .     ') array from b2.neutrals.parameters'
        end if
#endif
      enddo
      do is = 1, nnioni
        call xertst(0.le.lionscl(is), 'faulty input lionscl')
      end do
      do is = 1, nnstsi
        call xertst (0.le.lcns(is).and.lcns(is).le.nnstsi,
     &   'faulty input lcns')
        call xertst (0.le.ltns(is).and.ltns(is).le.nnstsi,
     &   'faulty input ltns')
      enddo
      do is = 1, nstrai
        call xertst (0.le.ksns(is), 'faulty input ksns')
        call xertst (ksns(is).le.DEF_NSRFS, 'faulty input ksns')
      enddo
      do is = 1, sum(ksns(1:nstrai))
        call xertst (lsns(is).le.nlim+nsts, 'faulty input lsns')
        call xertst (lsns(is).ge.-nsts, 'faulty input lsns')
        call xertst (lsns(is).ne.0, 'faulty input lsns')
      end do

      ii=minloc(lcns)
      if(lcns(ii(1)).gt.0) then !{
        ncns=nstsid
      else !}{
        ncns=ii(1)-1
      end if !}
      ii=minloc(ltns)
      if(ltns(ii(1)).gt.0) then !{
        ntns=nstsid
      else !}{
        ntns=ii(1)-1
      end if !}
      call xertst (0.le.ncns, 'faulty input ncns')
      call xertst (0.le.ntns, 'faulty input ntns')

      do j=1,nstrai !{
        if(lsns(j).gt.0 .and. lsns(j).le.nsts) lsns(j)=-lsns(j)
      end do !}
      n=0
      do j=1,nstrai !{
        if(ksns(j).gt.0) n=j
      end do !}
c*** if ksns was read, then recalculate the msns list
      k=1
      do j=1,n !{
        msns(1,j)=k
        msns(2,j)=ksns(j)
        k=k+ksns(j)
      end do !}
c<<<
      if(n.gt.0) then !{
        print *,'msns recalculated:'
        print '(3x,50i5)',(msns(1,j),j=1,n)
        print '(3x,50i5)',(msns(2,j),j=1,n)
      end if !}
c>>>

      do j=1,nstrai !{
        hlp=0.0_R8
        do i=1,nnatmi !{
          hlp=hlp+gpfc(i,j)
        end do !}
        if(hlp.gt.0.0_R8) then !{
          do i=1,nnatmi !{
            gpfc(i,j)=gpfc(i,j)/hlp
          end do !}
        end if !}
      end do !}

      if(neut_scl_lim.lt.1.0e-6_R8) neut_scl_lim=2.0_R8  ! to avoid log(=<0)

      track_species='  '
      do is=0,ns-1
        itrack=track_index(is)
        if(itrack.gt.0.and.itrack.le.ntrack) then
          if(track_species(itrack).eq.'  ') then
            track_species(itrack)=is_codes(is)
          endif
        endif
      enddo
      write(hlp_frm,'(a,i2,a)') '(a,',ntrack,'(a2,1x))'
      write(*,hlp_frm) 'track_species : Bulk ', track_species(1:ntrack)
c.. Test consistency with region boundaries
      if (.not.dotest) return
      do ic=1, nncut
        found = on_closed_surface(leftcut(ic)-1,-1)
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          found = rcend(istra).eq.leftcut(ic)-1 .and.
     &            crcstra(istra).eq.'S' .and. old_style(istra)
        end do
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          if (.not.old_style(istra) .and. crcstra(istra).eq.'S') then
            side = .false.
            ix = 1
            do while (.not.side .and. ix.le.rc_list_size(istra))
              side = rc_list_x(ix,istra).eq.leftcut(ic)-1 .and.
     .               rc_list_y(ix,istra).le.topcut(ic)
              ix = ix + 1
            end do
            ix = 1
            do while (side .and. ix.le.rc_list_size(istra))
              found = rc_list_x(ix,istra).eq.leftcut(ic) .and.
     .                rc_list_y(ix,istra).le.topcut(ic)
              ix = ix + 1
            end do
          end if
        end do
        call xertst (found,
     .   'Mismatch with neutrals boundaries (South) !')
        found = on_closed_surface(leftcut(ic),-1)
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          found = rcstart(istra).eq.leftcut(ic) .and.
     &            crcstra(istra).eq.'S' .and. old_style(istra)
        end do
        call xertst (found,
     .   'Mismatch with neutrals boundaries (South) !')
        found = on_closed_surface(rightcut(ic)-1,-1)
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          found = rcend(istra).eq.rightcut(ic)-1 .and.
     &            crcstra(istra).eq.'S' .and. old_style(istra)
        end do
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          if (.not.old_style(istra) .and. crcstra(istra).eq.'S') then
            side = .false.
            ix = 1
            do while (.not.side .and. ix.le.rc_list_size(istra))
              side = rc_list_x(ix,istra).eq.rightcut(ic)-1 .and.
     .               rc_list_y(ix,istra).le.topcut(ic)
              ix = ix + 1
            end do
            ix = 1
            do while (side .and. ix.le.rc_list_size(istra))
              found = rc_list_x(ix,istra).eq.rightcut(ic) .and.
     .                rc_list_y(ix,istra).le.topcut(ic)
              ix = ix + 1
            end do
          end if
        end do
        call xertst (found,
     .   'Mismatch with neutrals boundaries (South) !')
        found = on_closed_surface(rightcut(ic),-1)
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          found = rcstart(istra).eq.rightcut(ic) .and.
     &            crcstra(istra).eq.'S' .and. old_style(istra)
        end do
        call xertst (found,
     .   'Mismatch with neutrals boundaries (South) !')
      end do
c.. Test consistency with physical edges within the computational domain
      do ix=0,nx-1
        if (rightix(ix,-1).eq.nx+1 .and. region(ix,-1,0).ne.0) then
          found = isymm.eq.2
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            found = crcstra(istra).eq.'S' .and. rcend(istra).eq.ix-1
     .       .and. old_style(istra)
          end do
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            if (.not.old_style(istra) .and. crcstra(istra).eq.'S') then
              ic = 1
              do while (.not.found .and. ic.le.rc_list_size(istra))
                found = rc_list_x(ic,istra).eq.ix-1 .and.
     .                  rc_list_y(ic,istra).eq.-1
                ic = ic + 1
              end do
            end if
          end do
          call xertst (found,
     .     'Mismatch with right internal edge (South) !')
        endif
        if (rightix(ix,ny).eq.nx+1 .and. region(ix,ny,0).ne.0) then
          found = .false.
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            found = crcstra(istra).eq.'N' .and. rcend(istra).eq.ix-1
     .       .and. old_style(istra)
          end do
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            if (.not.old_style(istra) .and. crcstra(istra).eq.'N') then
              ic = 1
              do while (.not.found .and. ic.le.rc_list_size(istra))
                found = rc_list_x(ic,istra).eq.ix-1 .and.
     .                  rc_list_y(ic,istra).eq.ny
                ic = ic + 1
              end do
            end if
          end do
          call xertst (found,
     .     'Mismatch with right internal edge (North) !')
        endif
        if (leftix(ix,-1).eq.-2 .and. region(ix,-1,0).ne.0) then
          found = isymm.eq.2
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            found = crcstra(istra).eq.'S' .and. rcstart(istra).eq.ix+1
     .       .and. old_style(istra)
          end do
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            if (.not.old_style(istra) .and. crcstra(istra).eq.'S') then
              ic = 1
              do while (.not.found .and. ic.le.rc_list_size(istra))
                found = rc_list_x(ic,istra).eq.ix+1 .and.
     .                  rc_list_y(ic,istra).eq.-1
                ic = ic + 1
              end do
            end if
          end do
          call xertst (found,
     .     'Mismatch with left internal edge (South) !')
        endif
        if (leftix(ix,ny).eq.-2 .and. region(ix,ny,0).ne.0) then
          found = .false.
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            found = crcstra(istra).eq.'N' .and. rcstart(istra).eq.ix+1
     .       .and. old_style(istra)
          end do
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            if (.not.old_style(istra) .and. crcstra(istra).eq.'N') then
              ic = 1
              do while (.not.found .and. ic.le.rc_list_size(istra))
                found = rc_list_x(ic,istra).eq.ix+1 .and.
     .                  rc_list_y(ic,istra).eq.ny
                ic = ic + 1
              end do
            end if
          end do
          call xertst (found,
     .     'Mismatch with left internal edge (North) !')
        endif
      enddo
      return

 999  call xerrab ('Error while reading NEUTRALS namelist!')
      end subroutine read_neutrals_namelist_st

      subroutine read_neutrals_namelist_us(ns, m, dotest)
      use b2mod_wall
     & ,only : track_index, ntrack, track_species, track_chem_sput
      use b2mod_indirect
      use b2mod_elements
     & ,only : is_codes
      use b2mod_boundary_namelist
      use b2mod_b2cmpa
      use b2us_map
      implicit none
      type(mapping), intent(inout) :: m
      integer ns
      logical dotest
      integer is, ic, ix, istra, i, j, k, n, ii(1), itrack,
     >        noss, instra, iss, indss, idb, nrcl
      logical found, side
      intrinsic mod
      real (kind=R8) hlp
      logical :: file_ok
      logical, allocatable :: old_style(:)
      character*16 hlp_frm
      character*260 filename

      filename=neutrals_filename
      call find_file(filename,file_ok)
      if(file_ok) then
        open(99,file=filename)
        call cfverr(99,neutrals_version)
        write(*,*) 'Reading namelist neutrals from ',trim(filename)
#ifdef USE_MPI
        read(99,neutrals,err=999)
#else
        read(99,neutrals)
#endif
        write(*,*) 'Read namelist neutrals from ',trim(filename)
        close(99)
      else
        call xerrab ('No '//trim(filename))
      endif
      if (nstrai > 0) then
        
        m%nRc = nstrai
        m%mxnRc = count(m%fcLbl /= 0)
        call alloc_mapping_rc (m)
        m%rcFcP = 0
        
        open (newunit=idb,file='debug_neutrals.out')
        instra = 0
        do istra = 1, nstrai
          if (crcstra(istra) == 'A') then
            nrcl = 0
            m%rcFcP(istra,1) = instra +1 
            noss = rcend(istra) - rcstart(istra) + 1
            do iss = 1, noss
              indss = rcstart(istra)+iss-1
              write (*,*) ' surface index ',indss
! find all faces belonging to surface structure INDSS
              call find_faces (indss, instra, nrcl, 
     >                       m%mxnrc,m%rcFc, m%rcFcOr, m, idb)
            end do              ! iss	 
            m%rcFcP(istra,2) = nrcl
         elseif ((crcstra(istra) /= 'T') .and. (crcstra(istra) /= 'V')
     >      .and. (crcstra(istra) /= 'C')) then      
            call xerrab ('CRCSTRA /= ACTV ')
          endif
        end do
        close(idb)
      end if

      if (any(rc_list_size > 0))
     .  call xerrab 
     .  ('read_neutrals_namelist_us -- rc_list_size > 0')

      call rc_to_struct(m)

      if (nstrai.gt.0) then
cxpb  In Eirene facelift the time stratum is handled separately
cxpb  On the B2 side, use nstrat as nstrai + the time stratum
        if (crcstra(nstrai).eq.'T') then
          nstrat=nstrai
          nstrai=nstrai-1
        elseif (crcstra(nstrai+1).eq.'T') then
          nstrat=nstrai+1
        else
#ifdef B25_EIRENE
          nstrat=nstrai+1
          crcstra(nstrat)='T'
          rcpos(nstrat)=0
          rcstart(nstrat)=0
          rcend(nstrat)=0
#else
          nstrat=nstrai
#endif
        endif
#ifdef WG_TODO
        allocate(old_style(nstrat))
        do is = 1, nstrat
          if (rc_list_size(is).eq.0) then
            old_style(is) = .true.
            rc_list_size(is) = rcend(is) - rcstart(is) + 1
            if (crcstra(is).eq.'N' .or. crcstra(is).eq.'S') then
              do ix = rcstart(is), rcend(is)
                rc_list_x(ix-rcstart(is)+1,is) =  ix
              end do
              rc_list_y(1:rc_list_size(is),is) = rcpos(is)
            else if (crcstra(is).eq.'E' .or. crcstra(is).eq.'W') then
              rc_list_x(1:rc_list_size(is),is) = rcpos(is)
              do ix = rcstart(is), rcend(is)
                rc_list_y(ix-rcstart(is)+1,is) =  ix
              end do
            endif
          else
            old_style(is) = .false.
          end if
        end do
#endif
      else
        nstrat = 0
      endif

      arcstra = crcstra
      arcstart = rcstart
      arcend = arcend    
c
      if(chemical_erosion_redep_fac.ne.1.0_R8)
     * write(*,*)
     1 'Note: chemical_erosion_redep_fac = ',
     2 chemical_erosion_redep_fac
      if(chemical_erosion_be_fac)
     1 write(*,*) 'Note: suppression of chemical sputtering by Be ',
     2 'is active with parameters ',
     3 chemical_erosion_be_fac_a,
     4 chemical_erosion_be_fac_b,
     5 chemical_erosion_be_fac_c
c
      do istra = 1, nstrat
        do is = 1, nnatmi
          if (lstrascl(istra,is).lt.0) lstrascl(istra,is)=is
          if (lstrascl(istra,is).gt.nnatmi) lstrascl(istra,is)=is
        end do
      end do

*   ..set additional variables
      do is = 1, nstrat
        if (crcstra(is).ne.'N' .and. crcstra(is).ne.'S' .and.
     &      crcstra(is).ne.'W' .and. crcstra(is).ne.'E') targsp(is,:)=-1
      end do

*   ..test input
      call xertst (0.le.minval(track_index).and.
     1 maxval(track_index).le.ntrack, 'faulty input track_index')
      call xertst (0.le.nstrai, 'faulty input nstrai')
      call xertst (0.le.nstrat, 'faulty input nstrat')
      call xertst (0.0_R8.lt.rf_neut(1).and.0.0_R8.lt.rf_neut(2).and.
     .             0.0_R8.lt.rf_neut(3).and.0.0_R8.lt.rf_neut(4),
     .  'faulty input rf_neut')
#ifndef NO_B2_CHEM_SPUT
      call xertst (0.le.issput_chemical,
     .  'faulty input issput_chemical')
#endif
      call xertst (0.0_R8.le.eirene_step_cpu,
     . 'faulty input eirene_step_cpu')
      call xertst (0.0_R8.lt.eirene_step_dt,
     . 'faulty input eirene_step_dt')
      call xertst (1.0_R8.lt.neut_scl_lim, 'faulty input neut_scl_lim')
      call xertst (0.0_R8.le.volrecinc, 'faulty input volrecinc')
      call xertst (0.0_R8.lt.volrecwt, 'faulty input volrecwt')
      call xertst (0.0_R8.le.neutrals_time_mod,
     . 'faulty input neutrals_time_mod')
      call xertst (0.0_R8.le.neutrals_time_switch,
     . 'faulty input neutrals_time_switch')

      do istra = 1, nstrat
#ifdef WG_TODO
        if (crcstra(istra).eq.'N') then
          if (old_style(istra)) then
            call xertst (1.le.rcpos(istra).and.rcpos(istra).le.ny+1,
     .       'faulty input rcpos (North)')
            call xertst (-1.le.rcstart(istra).and.
     .                   rcstart(istra).le.nx+1,
     .       'faulty input rcstart (North)')
            call xertst (-1.le.rcend(istra).and.rcend(istra).le.nx+1,
     .       'faulty input rcend (North)')
          else
            do ix = 1, rc_list_size(istra)
              call xertst (1.le.rc_list_y(ix,istra).and.
     .                     rc_list_y(ix,istra).le.ny+1,
     .                     'faulty input rc_list_y (North)')
              call xertst (-1.le.rc_list_x(ix,istra).and.
     .                     rc_list_x(ix,istra).le.nx+1,
     .                     'faulty input rc_list_x (North)')
            end do
          end if
        else if (crcstra(istra).eq.'S') then
          if (old_style(istra)) then
            call xertst (-1.le.rcpos(istra).and.rcpos(istra).le.ny-1,
     .       'faulty input rcpos (South)')
            call xertst (-1.le.rcstart(istra).and.
     .                   rcstart(istra).le.nx+1,
     .       'faulty input rcstart (South)')
            call xertst (-1.le.rcend(istra).and.rcend(istra).le.nx+1,
     .       'faulty input rcend (South)')
          else
            do ix = 1, rc_list_size(istra)
              call xertst (-1.le.rc_list_y(ix,istra).and.
     .                     rc_list_y(ix,istra).le.ny-1,
     .                     'faulty input rc_list_y (South)')
              call xertst (-1.le.rc_list_x(ix,istra).and.
     .                     rc_list_x(ix,istra).le.nx+1,
     .                     'faulty input rc_list_x (South)')
            end do
          end if
        else if (crcstra(istra).eq.'E') then
          if (old_style(istra)) then
            call xertst (1.le.rcpos(istra).and.rcpos(istra).le.nx+1,
     .       'faulty input rcpos (East)')
            call xertst (-1.le.rcstart(istra).and.
     .       rcstart(istra).le.ny+1,
     .       'faulty input rcstart (East)')
            call xertst (-1.le.rcend(istra).and.rcend(istra).le.ny+1,
     .       'faulty input rcend (East)')
          else
            do ix = 1, rc_list_size(istra)
              call xertst (1.le.rc_list_x(ix,istra).and.
     .                     rc_list_x(ix,istra).le.nx+1,
     .                     'faulty input rc_list_x (East)')
              call xertst (-1.le.rc_list_y(ix,istra).and.
     .                     rc_list_y(ix,istra).le.ny+1,
     .                     'faulty input rc_list_y (East)')
            end do
          end if
        else if (crcstra(istra).eq.'W') then
          if (old_style(istra)) then
            call xertst (-1.le.rcpos(istra).and.rcpos(istra).le.nx-1,
     .       'faulty input rcpos (West)')
            call xertst (-1.le.rcstart(istra).and.
     .       rcstart(istra).le.ny+1,
     .       'faulty input rcstart (West)')
            call xertst (-1.le.rcend(istra).and.rcend(istra).le.ny+1,
     .       'faulty input rcend (West)')
          else
            do ix = 1, rc_list_size(istra)
              call xertst (-1.le.rc_list_x(ix,istra).and.
     .                     rc_list_x(ix,istra).le.nx-1,
     .                     'faulty input rc_list_x (West)')
              call xertst (-1.le.rc_list_y(ix,istra).and.
     .                     rc_list_y(ix,istra).le.ny+1,
     .                     'faulty input rc_list_y (West)')
            end do
          end if
        endif
        if (old_style(istra))
     .   call xertst (rcstart(istra).le.rcend(istra),
     .   'faulty inputs rcstart > rcend')
#endif
        call xertst (maxval(targsp(istra,:)).lt.ns,
     .   'faulty inputs targsp > ns')
        call xertst (targsp(istra,1).ge.-1,
     .   'faulty inputs targsp(,1)')
        call xertst (0.0_R8.le.chemical_erosion_redep_fac,
     .   'faulty input chemical_erosion_redep_fac')
        call xertst (0.0_R8.lt.chemical_erosion_be_fac_a,
     .   'faulty input chemical_erosion_be_fac_a')
        call xertst (0.0_R8.lt.chemical_erosion_be_fac_b,
     .   'faulty input chemical_erosion_be_fac_b')
        call xertst (0.0_R8.le.chemical_erosion_be_fac_c
     &          .and.1.0_R8.ge.chemical_erosion_be_fac_c,
     .   'faulty input chemical_erosion_be_fac_c')
        call xertst (0.0_R8.le.volrecstart(istra),
     .   'faulty input volrecstart')
        call xertst (0.le.species_start(istra).and.
     &   species_start(istra).lt.ns, 'faulty input species_start')
        call xertst (0.le.species_end(istra).and.
     &   species_end(istra).lt.ns, 'faulty input species_end')
        call xertst (species_start(istra).le.species_end(istra),
     .   'faulty inputs species_start > species_end')
        do is = 0, ns-1
          call xertst (0.0_R8.le.recyc(is,istra), 'faulty input recyc')
          call xertst (0.0_R8.le.erecyc(is,istra),
     .     'faulty input erecyc')
          call xertst (0.0_R8.le.mrecyc(is,istra),
     .     'faulty input mrecyc')
          call xertst (0.0_R8.le.rcion(is,istra), 'faulty input rcion')
          call xertst (0.0_R8.le.phys_sput(is,istra),
     .     'faulty input phys_sput')
          call xertst (0.0_R8.le.chem_sput(is,istra),
     .     'faulty input chem_sput')
        enddo
        call xertst (0.0_R8.le.recyceir(istra), 'faulty input recyceir')
        do is = 1, nnatmi
          call xertst(0.0_R8.le.gpfc(is,istra), 'faulty input gpfc')
          call xertst(0.le.lstrascl(istra,is) .and.
     .                lstrascl(istra,is).le.nnatmi,
     .                'faulty input lstrascl')
        enddo
        if (time_dep_puff(istra)) then
          call xertst(crcstra(istra).eq.'C',
     .     'Time_dep_puff set for a non-puff stratum')
          call xertst(ngpdata(istra).gt.0,
     .     'Time_dep_puff data list size not declared')
          call xertst(ngpdata(istra).le.100,
     .     'Time_dep_puff data list size too large')
          do ix = 1, ngpdata(istra)
            call xertst(gpdata(ix,2,istra).ge.0.0_R8,
     .       'Negative gas puff strength declared!')
          end do
        endif
      enddo
      do is = 0, ns-1
        call xertst(1.le.b2espcr(is).and.b2espcr(is).le.nspecies,
     .   'faulty input b2espcr')
        if(nnatmi.eq.0) then
          call xertst(1.le.b2eatcr(is).and.b2eatcr(is).le.nspecies,
     .     'faulty input b2eatcr')
        else
          call xertst(1.le.b2eatcr(is).and.b2eatcr(is).le.nnatmi,
     .     'faulty input b2eatcr')
        endif
      enddo
      do is = 1, nnatmi
        call xertst(0.le.eb2atcr(is).and.eb2atcr(is).lt.ns,
     .   'faulty input eb2atcr')
        call xertst(0.le.latmscl(is), 'faulty input latmscl')
#ifdef B25_EIRENE
        if (latmscl(is).eq.0) then
          write(6,'(a,i3,a)') 'Atom ',is,
     .     ' in Eirene does not scale with any B2.5 sequence'
          write(6,'(a)')
     .     'Please provide a scaling using the latmscl array '//
     .     'in b2.neutral.parameters'
          write(0,'(a,i3,a)') 'Atom ',is,
     .     ' in Eirene does not scale with any B2.5 sequence'
          write(0,'(a)')
     .     'Please provide a scaling using the latmscl array '//
     .     'in b2.neutral.parameters'
        end if
#endif
      enddo
      do is = 1, nspecies
        call xertst(0.le.eb2spcr(is).and.eb2spcr(is).lt.ns,
     .   'faulty input eb2spcr')
      enddo
      do is = 1, nnmoli
        do i = 1, nnatmi
          call xertst(0.le.mlcmp(i,is), 'faulty input mlcmp')
        enddo
        call xertst(0.le.lmolscl(is), 'faulty input lmolscl')
#ifdef B25_EIRENE
        if (lmolscl(is).eq.0) then
          write(6,'(a,i3,a)') 'Molecule ',is,
     .     ' does not scale with any atom species'
          write(6,'(a)')
     .     'Please provide a scaling using the lmolscl array '//
     .     'in b2.neutral.parameters'
          write(0,'(a,i3,a)') 'Molecule ',is,
     .     ' does not scale with any atom species'
          write(0,'(a)')
     .     'Please provide a scaling using the lmolscl array '//
     .     'in b2.neutral.parameters'
        end if
        if (maxval(mlcmp(1:nnatmi,is)).eq.0) then
          write(6,'(a,i3,a)')
     .     'Composition of molecule ',is,' is not given'
          write(6,'(a,i3,a)')
     .     'If it contains species followed in the plasma, '//
     .     'please specify them using the mlcmp(:,',is,
     .     ') array from b2.neutrals.parameters'
          write(0,'(a,i3,a)')
     .     'Composition of molecule ',is,' is not given'
          write(0,'(a,i3,a)')
     .     'If it contains species followed in the plasma, '//
     .     'please specify them using the mlcmp(:,',is,
     .     ') array from b2.neutrals.parameters'
        end if
#endif
      enddo
      do is = 1, nnioni
        call xertst(0.le.lionscl(is), 'faulty input lionscl')
      end do
      do is = 1, nnstsi
        call xertst (0.le.lcns(is).and.lcns(is).le.nnstsi,
     &   'faulty input lcns')
        call xertst (0.le.ltns(is).and.ltns(is).le.nnstsi,
     &   'faulty input ltns')
      enddo
      do is = 1, nstrai
        call xertst (0.le.ksns(is), 'faulty input ksns')
        call xertst (ksns(is).le.DEF_NSRFS, 'faulty input ksns')
      enddo
      do is = 1, sum(ksns(1:nstrai))
        call xertst (lsns(is).le.nlim+nsts, 'faulty input lsns')
        call xertst (lsns(is).ge.-nsts, 'faulty input lsns')
        call xertst (lsns(is).ne.0, 'faulty input lsns')
      end do

      ii=minloc(lcns)
      if(lcns(ii(1)).gt.0) then !{
        ncns=nstsid
      else !}{
        ncns=ii(1)-1
      end if !}
      ii=minloc(ltns)
      if(ltns(ii(1)).gt.0) then !{
        ntns=nstsid
      else !}{
        ntns=ii(1)-1
      end if !}
      call xertst (0.le.ncns, 'faulty input ncns')
      call xertst (0.le.ntns, 'faulty input ntns')

      do j=1,nstrai !{
        if(lsns(j).gt.0 .and. lsns(j).le.nsts) lsns(j)=-lsns(j)
      end do !}
      n=0
      do j=1,nstrai !{
        if(ksns(j).gt.0) n=j
      end do !}
c*** if ksns was read, then recalculate the msns list
      k=1
      do j=1,n !{
        msns(1,j)=k
        msns(2,j)=ksns(j)
        k=k+ksns(j)
      end do !}
c<<<
      if(n.gt.0) then !{
        print *,'msns recalculated:'
        print '(3x,50i5)',(msns(1,j),j=1,n)
        print '(3x,50i5)',(msns(2,j),j=1,n)
      end if !}
c>>>

      do j=1,nstrai !{
        hlp=0.0_R8
        do i=1,nnatmi !{
          hlp=hlp+gpfc(i,j)
        end do !}
        if(hlp.gt.0.0_R8) then !{
          do i=1,nnatmi !{
            gpfc(i,j)=gpfc(i,j)/hlp
          end do !}
        end if !}
      end do !}

      if(neut_scl_lim.lt.1.0e-6_R8) neut_scl_lim=2.0_R8  ! to avoid log(=<0)

      track_species='  '
      do is=0,ns-1
        itrack=track_index(is)
        if(itrack.gt.0.and.itrack.le.ntrack) then
          if(track_species(itrack).eq.'  ') then
            track_species(itrack)=is_codes(is)
          endif
        endif
      enddo
      write(hlp_frm,'(a,i2,a)') '(a,',ntrack,'(a2,1x))'
      write(*,hlp_frm) 'track_species : Bulk ', track_species(1:ntrack)
c.. Test consistency with region boundaries
      if (.not.dotest) return
#ifdef WG_TODO
      do ic=1, nncut
        found = on_closed_surface(leftcut(ic)-1,-1)
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          found = rcend(istra).eq.leftcut(ic)-1 .and.
     &            crcstra(istra).eq.'S' .and. old_style(istra)
        end do
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          if (.not.old_style(istra) .and. crcstra(istra).eq.'S') then
            side = .false.
            ix = 1
            do while (.not.side .and. ix.le.rc_list_size(istra))
              side = rc_list_x(ix,istra).eq.leftcut(ic)-1 .and.
     .               rc_list_y(ix,istra).le.topcut(ic)
              ix = ix + 1
            end do
            ix = 1
            do while (side .and. ix.le.rc_list_size(istra))
              found = rc_list_x(ix,istra).eq.leftcut(ic) .and.
     .                rc_list_y(ix,istra).le.topcut(ic)
              ix = ix + 1
            end do
          end if
        end do
        call xertst (found,
     .   'Mismatch with neutrals boundaries (South) !')
        found = on_closed_surface(leftcut(ic),-1)
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          found = rcstart(istra).eq.leftcut(ic) .and.
     &            crcstra(istra).eq.'S' .and. old_style(istra)
        end do
        call xertst (found,
     .   'Mismatch with neutrals boundaries (South) !')
        found = on_closed_surface(rightcut(ic)-1,-1)
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          found = rcend(istra).eq.rightcut(ic)-1 .and.
     &            crcstra(istra).eq.'S' .and. old_style(istra)
        end do
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          if (.not.old_style(istra) .and. crcstra(istra).eq.'S') then
            side = .false.
            ix = 1
            do while (.not.side .and. ix.le.rc_list_size(istra))
              side = rc_list_x(ix,istra).eq.rightcut(ic)-1 .and.
     .               rc_list_y(ix,istra).le.topcut(ic)
              ix = ix + 1
            end do
            ix = 1
            do while (side .and. ix.le.rc_list_size(istra))
              found = rc_list_x(ix,istra).eq.rightcut(ic) .and.
     .                rc_list_y(ix,istra).le.topcut(ic)
              ix = ix + 1
            end do
          end if
        end do
        call xertst (found,
     .   'Mismatch with neutrals boundaries (South) !')
        found = on_closed_surface(rightcut(ic),-1)
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          found = rcstart(istra).eq.rightcut(ic) .and.
     &            crcstra(istra).eq.'S' .and. old_style(istra)
        end do
        call xertst (found,
     .   'Mismatch with neutrals boundaries (South) !')
      end do
#endif
c.. Test consistency with physical edges within the computational domain
#ifdef WG_TODO
      do ix=0,nx-1
        if (rightix(ix,-1).eq.nx+1 .and. region(ix,-1,0).ne.0) then
          found = isymm.eq.2
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            found = crcstra(istra).eq.'S' .and. rcend(istra).eq.ix-1
     .       .and. old_style(istra)
          end do
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            if (.not.old_style(istra) .and. crcstra(istra).eq.'S') then
              ic = 1
              do while (.not.found .and. ic.le.rc_list_size(istra))
                found = rc_list_x(ic,istra).eq.ix-1 .and.
     .                  rc_list_y(ic,istra).eq.-1
                ic = ic + 1
              end do
            end if
          end do
          call xertst (found,
     .     'Mismatch with right internal edge (South) !')
        endif
        if (rightix(ix,ny).eq.nx+1 .and. region(ix,ny,0).ne.0) then
          found = .false.
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            found = crcstra(istra).eq.'N' .and. rcend(istra).eq.ix-1
     .       .and. old_style(istra)
          end do
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            if (.not.old_style(istra) .and. crcstra(istra).eq.'N') then
              ic = 1
              do while (.not.found .and. ic.le.rc_list_size(istra))
                found = rc_list_x(ic,istra).eq.ix-1 .and.
     .                  rc_list_y(ic,istra).eq.ny
                ic = ic + 1
              end do
            end if
          end do
          call xertst (found,
     .     'Mismatch with right internal edge (North) !')
        endif
        if (leftix(ix,-1).eq.-2 .and. region(ix,-1,0).ne.0) then
          found = isymm.eq.2
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            found = crcstra(istra).eq.'S' .and. rcstart(istra).eq.ix+1
     .       .and. old_style(istra)
          end do
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            if (.not.old_style(istra) .and. crcstra(istra).eq.'S') then
              ic = 1
              do while (.not.found .and. ic.le.rc_list_size(istra))
                found = rc_list_x(ic,istra).eq.ix+1 .and.
     .                  rc_list_y(ic,istra).eq.-1
                ic = ic + 1
              end do
            end if
          end do
          call xertst (found,
     .     'Mismatch with left internal edge (South) !')
        endif
        if (leftix(ix,ny).eq.-2 .and. region(ix,ny,0).ne.0) then
          found = .false.
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            found = crcstra(istra).eq.'N' .and. rcstart(istra).eq.ix+1
     .       .and. old_style(istra)
          end do
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            if (.not.old_style(istra) .and. crcstra(istra).eq.'N') then
              ic = 1
              do while (.not.found .and. ic.le.rc_list_size(istra))
                found = rc_list_x(ic,istra).eq.ix+1 .and.
     .                  rc_list_y(ic,istra).eq.ny
                ic = ic + 1
              end do
            end if
          end do
          call xertst (found,
     .     'Mismatch with left internal edge (North) !')
        endif
      enddo
#endif
      return

 999  call xerrab ('Error while reading NEUTRALS namelist!')
      end subroutine read_neutrals_namelist_us

      subroutine write_b2mod_neutrals_namelist
      implicit none
      if (write_nml_neut) write(*,neutrals)
      end subroutine write_b2mod_neutrals_namelist

      subroutine read_b2mod_neutrals_save(nread,edition,found)
      implicit none
      integer nread
      character*4 edition
      logical found
      character*256 filename
      integer istra
      logical file_ok

      saved_volrec(:)=0.0_R8

      if (edition.eq.'    ') then
        filename = 'b2.neutrals_save.parameters'
      else
        filename = 'b2.neutrals_save.parameters.'//edition
      endif
      call find_file(filename,file_ok)
      if(file_ok) then
        open(nread,file=trim(filename))
        write(6,*) 'Reading ',trim(filename)
#ifdef USE_MPI
        read(nread,neutrals_save,err=999)
#else
        read(nread,neutrals_save)
#endif
        write(6,*) 'Read ',trim(filename)
        close(nread)
      else
        write(*,*) 'No ', trim(filename)
      endif
      found = maxval(saved_volrec).gt.0.0_R8

*   ..test input
      do istra = 1, nstrat
        call xertst (0.0_R8.le.saved_volrec(istra),
     .   'faulty input saved_volrec')
      enddo

      return
 999  call xerrab ('Error while reading NEUTRALS_SAVE namelist!')
      end subroutine read_b2mod_neutrals_save

      subroutine write_b2mod_neutrals_save(nwrite,edition)
      implicit none
      integer nwrite
      character*4 edition
      character*256 filename

      if (edition.eq.'    ') then
        filename = 'b2.neutrals_save.parameters'
      else
        filename = 'b2.neutrals_save.parameters.'//edition
      endif
      if (maxval(saved_volrec).gt.0.0_R8) then
        open(nwrite,file=trim(filename))
        write(nwrite,neutrals_save)
        close(nwrite)
        write(*,*) 'Saved ',trim(filename)
      endif
      return

      end subroutine write_b2mod_neutrals_save



***********************************************************************


      subroutine calc_mapping_rc (nx, ny, ns, m)
      use b2mod_indirect
      implicit none
      integer :: nx, ny, ns, istra, istrac, instra, iout, ian, ien, i
      integer :: rcs(nstrat), rce(nstrat), rcl(nstrat)
      type (mapping), intent (inout) :: m

      m%nRc   = nstrai  ! number of recycling boundaries
      m%mxnRc = 0       ! number of cells with imposed recycling
      do istra = 1, nstrai
        rcs(istra) = rcstart(istra)
        rce(istra) = rcend(istra)
        ! exclude corner cells if present
        if (crcstra(istra).eq.'N' .or. crcstra(istra).eq.'S') then
          if (leftix (rcstart(istra),rcpos(istra)) .eq. -2)  
     >      rcs(istra) = rcs(istra) + 1
          if (rightix(rcend(istra),  rcpos(istra)) .eq. nx+1) 
     >      rce(istra) = rce(istra) - 1
          rcl(istra) = rce(istra) - rcs(istra) + 1
          m%mxnRc = m%mxnRc + rcl(istra)
        elseif (crcstra(istra).eq.'E' .or. crcstra(istra).eq.'W') then
          if (bottomiy(rcpos(istra),rcstart(istra)) .eq. -2)
     >      rcs(istra) = rcs(istra) + 1
          if (topiy   (rcpos(istra),rcend(istra)  ) .eq. ny+1)
     >      rce(istra) = rce(istra) - 1
          rcl(istra) = rce(istra) - rcs(istra) + 1
          m%mxnRc = m%mxnRc + rcl(istra)
        elseif (crcstra(istra).eq.'V' .or. crcstra(istra).eq.'T' .or.
     >          crcstra(istra).eq.'C') then
          rcl(istra) = 0
        else
          call xerrab ('Wrong CRCSTRA or not yet available for WG')
        endif
      enddo

      call alloc_mapping_rc (m)

      instra = 0 ! index in rcCv
      do istra = 1, nstrai
        m%rcCvP(istra,1) = instra + 1
        m%rcCvP(istra,2) = rcl(istra)
        if (crcstra(istra).eq.'S') then
          do istrac = 1, rcl(istra)
            instra = instra + 1
            m%rcCv(instra,1) = 
     >        m%imapCv(rcs(istra)+istrac-1,rcpos(istra))      ! number of the guard cell
            m%rcCv(instra,2) = min(m%fcCv(m%cvFc(             ! number of the corresponding domain cell
     >        m%cvFcP(m%imapCv(rcs(istra)+istrac-1,rcpos(istra)),1)),1),
     >        m%fcCv(m%cvFc(
     >        m%cvFcP(m%imapCv(rcs(istra)+istrac-1,rcpos(istra)),1)),2))
            m%rcFcOr(instra) = -1.0_R8
          enddo
        elseif (crcstra(istra).eq.'N') then
          do istrac = 1, rcl(istra)
            instra = instra + 1
            m%rcCv(instra,1) = 
     >        m%imapCv(rcs(istra)+istrac-1,rcpos(istra))      ! number of the guard cell
            m%rcCv(instra,2) = min(m%fcCv(m%cvFc(             ! number of the corresponding domain cell
     >        m%cvFcP(m%imapCv(rcs(istra)+istrac-1,rcpos(istra)),1)),1),
     >        m%fcCv(m%cvFc(
     >        m%cvFcP(m%imapCv(rcs(istra)+istrac-1,rcpos(istra)),1)),2))
            m%rcFcOr(instra) = 1.0_R8
          enddo
        elseif (crcstra(istra).eq.'W') then
          do istrac = 1, rcl(istra)
            instra = instra + 1
           m%rcCv(instra,1) = 
     >        m%imapCv(rcpos(istra),rcs(istra)+istrac-1)      ! number of the guard cell
            m%rcCv(instra,2) = min(m%fcCv(m%cvFc(             ! number of the corresponding domain cell
     >        m%cvFcP(m%imapCv(rcpos(istra),rcs(istra)+istrac-1),1)),1),
     >        m%fcCv(m%cvFc(
     >        m%cvFcP(m%imapCv(rcpos(istra),rcs(istra)+istrac-1),1)),2))
            m%rcFcOr(instra) = -1.0_R8
          enddo
        elseif (crcstra(istra).eq.'E') then 
          do istrac = 1, rcl(istra)
            instra = instra + 1
            m%rcCv(instra,1) = 
     >        m%imapCv(rcpos(istra),rcs(istra)+istrac-1)      ! number of the guard cell
            m%rcCv(instra,2) = min(m%fcCv(m%cvFc(             ! number of the corresponding domain cell
     >        m%cvFcP(m%imapCv(rcpos(istra),rcs(istra)+istrac-1),1)),1),
     >        m%fcCv(m%cvFc(
     >        m%cvFcP(m%imapCv(rcpos(istra),rcs(istra)+istrac-1),1)),2))
            m%rcFcOr(instra) = 1.0_R8
          enddo
        endif
      enddo 

      open (newunit=iout,file='mapping_rc.testout')
      write (iout,*) 'nstrai =',m%nRc
      write (iout,'(a,/(10i6))') 'rcCvP(1) ',
     .      (m%rcCvP(istra,1),istra=1,m%nRc)
      write (iout,'(a,/(10i6))') 'rcCvP(2) ',
     .      (m%rcCvP(istra,2),istra=1,m%nRc)
      do istra = 1, m%nRc
        ian = m%rcCvP(istra,1)
        ien = m%rcCvP(istra,1) + m%rcCvP(istra,2) - 1
        write (iout,*) 'istra = ',istra
        write (iout,'(a,/(10i6))') 'rcCv(1) ',(m%rcCv(i,1),i=ian,ien)
        write (iout,'(a,/(10i6))') 'rcCv(2) ',(m%rcCv(i,2),i=ian,ien)
        write (iout,'(a,/(10f6.2))') 'rcFcOr  ',(m%rcFcOr(i),i=ian,ien)
      end do
      close (iout)
    
      return
      end subroutine calc_mapping_rc



***********************************************************************


      subroutine conv_rcparms (nx, ny, ns, m, g)
      use b2mod_indirect
      implicit none
      integer :: nx, ny, ns, istra, istrac, instra,
     .           noss, iss, indss, i, idb
      integer :: rcs(nstrat), rce(nstrat), rcl(nstrat)
      integer :: minlbl(nstrat), maxlbl(nstrat)
      character(1) :: lchar(nstrat)
      type (mapping), intent (inout) :: m
      type (geometry), intent (inout) :: g

      m%nRc   = nstrai  ! number of recycling boundaries
      m%mxnRc = 0       ! number of cells with imposed recycling

      do istra = 1, nstrai
        rcs(istra) = rcstart(istra)
        rce(istra) = rcend(istra)
        if (crcstra(istra).eq.'A') then
! find the number of faces for this boundary condition
cwdk need more complex test here; present scheme double counts;
cwdk also need to correct for degenerate faces?
          rcl(istra) = count((cflags(:,:,2:5) >= rcstart(istra)) .and.
     >                    (cflags(:,:,2:5) <= rcend(istra)))
          m%mxnRc = m%nCg
        else
        ! exclude corner cells if present
        if (crcstra(istra).eq.'N' .or. crcstra(istra).eq.'S') then
          if (leftix (rcstart(istra),rcpos(istra)) .eq. -2)  
     >      rcs(istra) = rcs(istra) + 1
          if (rightix(rcend(istra),  rcpos(istra)) .eq. nx+1) 
     >      rce(istra) = rce(istra) - 1
          rcl(istra) = rce(istra) - rcs(istra) + 1
          m%mxnRc = m%mxnRc + rcl(istra)
        elseif (crcstra(istra).eq.'E' .or. crcstra(istra).eq.'W') then
          if (bottomiy(rcpos(istra),rcstart(istra)) .eq. -2)
     >      rcs(istra) = rcs(istra) + 1
          if (topiy   (rcpos(istra),rcend(istra)  ) .eq. ny+1)
     >      rce(istra) = rce(istra) - 1
          rcl(istra) = rce(istra) - rcs(istra) + 1
          m%mxnRc = m%mxnRc + rcl(istra)
        elseif (crcstra(istra).eq.'V' .or. crcstra(istra).eq.'T' .or.
     >          crcstra(istra).eq.'C') then
          rcl(istra) = 0
        else
          call xerrab ('Wrong CRCSTRA or not yet available for WG')
        endif
        endif
      enddo

      call alloc_mapping_rc (m)
      open (newunit=idb,file='recyc.out')

      minlbl = huge(1)
      maxlbl = -huge(1)
      lchar = ' '
      lchar(1:nstrai) = crcstra(1:nstrai)

      instra = 0 ! index in rcCv
      do istra = 1, nstrai
        m%rcFcP(istra,1) = instra + 1
        m%rcFcP(istra,2) = rcl(istra)
        if (crcstra(istra).eq.'S') then
          do istrac = 1, rcl(istra)
            instra = instra + 1
            m%rcFc(instra) = m%cvFc(
     >         m%cvFcP(m%imapCv(rcs(istra)+istrac-1,rcpos(istra)),1))    ! number of the corresponding domain face
            m%rcFcOr(instra) = -1.0_R8
            minlbl(istra) = min(minlbl(istra),m%fcLbl(m%rcFc(instra)))
            maxlbl(istra) = max(maxlbl(istra),m%fcLbl(m%rcFc(instra)))
            lchar(istra) = 'A'
          enddo
          rcfe(istra) = 0.0_R8
          rcfi(istra) = 2.0_R8
        elseif (crcstra(istra).eq.'N') then
          do istrac = 1, rcl(istra)
            instra = instra + 1
            m%rcFc(instra) = m%cvFc(
     >         m%cvFcP(m%imapCv(rcs(istra)+istrac-1,rcpos(istra)),1))    ! number of the corresponding domain face
            m%rcFcOr(instra) = 1.0_R8
            minlbl(istra) = min(minlbl(istra),m%fcLbl(m%rcFc(instra)))
            maxlbl(istra) = max(maxlbl(istra),m%fcLbl(m%rcFc(instra)))
            lchar(istra) = 'A'
          enddo
          rcfe(istra) = 0.0_R8
          rcfi(istra) = 2.0_R8
        elseif (crcstra(istra).eq.'W') then
          do istrac = 1, rcl(istra)
            instra = instra + 1
            m%rcFc(instra) = m%cvFc(
     >         m%cvFcP(m%imapCv(rcpos(istra),rcs(istra)+istrac-1),1))    ! number of the corresponding domain face
            m%rcFcOr(instra) = -1.0_R8
            minlbl(istra) = min(minlbl(istra),m%fcLbl(m%rcFc(instra)))
            maxlbl(istra) = max(maxlbl(istra),m%fcLbl(m%rcFc(instra)))
            lchar(istra) = 'A'
          enddo
          rcfe(istra) = 0.5_R8
          rcfi(istra) = 3.0_R8
        elseif (crcstra(istra).eq.'E') then 
          do istrac = 1, rcl(istra)
            instra = instra + 1
            m%rcFc(instra) = m%cvFc(
     >         m%cvFcP(m%imapCv(rcpos(istra),rcs(istra)+istrac-1),1))    ! number of the corresponding domain face
            m%rcFcOr(instra) = 1.0_R8
            minlbl(istra) = min(minlbl(istra),m%fcLbl(m%rcFc(instra)))
            maxlbl(istra) = max(maxlbl(istra),m%fcLbl(m%rcFc(instra)))
            lchar(istra) = 'A'
          enddo
          rcfe(istra) = 0.5_R8
          rcfi(istra) = 3.0_R8
        elseif (crcstra(istra).eq.'A') then 
          minlbl(istra) = rcstart(istra)
          maxlbl(istra) = rcend(istra)
          lchar(istra) = 'A'
          rcl(istra) = 0
          noss = rcend(istra) - rcstart(istra) + 1
	  do iss = 1, noss
	    indss = rcstart(istra)+iss-1
            write (*,*) ' surface index ',indss
!  find all faces belonging to surface structure INDSS
            call find_faces (indss, instra, rcl(istra), 
     >                       m%mxnrc,m%rcFc, m%rcFcOr, m, idb)
	  end do  ! iss	 
          m%rcFcP(istra,2) = rcl(istra) 
          rcfe(istra) = 0.5_R8
          rcfi(istra) = 3.0_R8
        else
!          call xerrab ('conv_rcparms -- wrong CRCSTRA')
        endif
      enddo
      close (idb)

      arcstra(1:nstrai) = lchar(1:nstrai)
      where (minlbl == huge(1))
        minlbl = 0
      end where
      where (maxlbl == -huge(1))
        maxlbl = 0
      end where
      arcstart(1:nstrai) = minlbl(1:nstrai)
      arcend(1:nstrai) = maxlbl(1:nstrai)

      if (crcstra(nstrai) == 'T') then
        arcstart(nstrai) = 0
        arcend(nstrai) = 0
      end if
      return
      end subroutine conv_rcparms

***********************************************************************

      subroutine rc_to_struct (m)
      type(mapping), intent(inout) :: m
      integer :: istra, instra, istrac, ifc, ic1, ic2
      integer :: iout, ian, ien, i

      m%nRc = nstrai
!pb   m%mxnRc = sum(nrcfaces(1:nstrai))       ! number of cells with imposed stratum
      m%mxnRc = sum(m%rcFcP(1:nstrai,2))       ! number of cells with imposed stratum

      call alloc_mapping_rc(m)

      instra = 0 ! index in rcCv
      do istra = 1,m%nRc
        m%rcCvP(istra,1) = instra + 1
!pb     m%rcCvP(istra,2) = nrcfaces(istra)
        m%rcCvP(istra,2) = m%rcFcp(istra,2)
        m%rcFcP(istra,1) = instra + 1
!pb     m%rcFcP(istra,2) = nrcfaces(istra)
!pb     do istrac = 1, nrcfaces(istra)
        do istrac = 1, m%rcFcp(istra,2)
          instra = instra + 1
!pb       ifc = rc_faces(istrac,istra)
          ifc = m%rcFc(m%rcFcP(istra,1)+istrac-1)
          ic1 = min(m%fcCv(ifc,1), m%fcCv(ifc,2))
          ic2 = max(m%fcCv(ifc,1), m%fcCv(ifc,2))
          if ((ic1 <= m%nCi) .and. (ic2 >m%nCi)) then
            m%rcCv(instra,1) = ic2
            m%rcCv(instra,2) = ic1
!pb         m%rcFcOr(instra) = rc_face_ori(istrac,istra)
!pb         m%rcFc(instra) = ifc
          else
            call xerrab('rc_to_struct -- problem with boundary face') 
          end if
        end do
      end do

      open (newunit=iout,file='rc.testout')
      write (iout,*) 'nstrat =',m%nRc
      write (iout,'(a,/(10i6))') 'rcCvP(1) ',
     .      (m%rcCvP(istra,1),istra=1,m%nRc)
      write (iout,'(a,/(10i6))') 'rcCvP(2) ',
     .      (m%rcCvP(istra,2),istra=1,m%nRc)
      do istra = 1, m%nRc
        ian = m%rcCvP(istra,1)
        ien = m%rcCvP(istra,1) + m%rcCvP(istra,2) - 1
        write (iout,*) 'istra = ',istra
        write (iout,'(a,/(10i6))') 'rcCv(1) ',(m%rcCv(i,1),i=ian,ien)
        write (iout,'(a,/(10i6))') 'rcCv(2) ',(m%rcCv(i,2),i=ian,ien)
        write (iout,'(a,/(10f6.2))') 'rcFcOr  ',(m%rcFcOr(i),i=ian,ien)
      end do
      close (iout)

      end subroutine rc_to_struct



***********************************************************************

      subroutine write_rcparms (ns, m)
      use b2mod_b2cmpa
      implicit none
      integer, intent(in) :: ns
      type (mapping), intent (inout) :: m
      integer :: istra, iout, ian, ien, is, i
      character(1) :: ap =''''
      character(10) :: new_neutrals_version = '03.002.000'

      open(newunit=iout,file='b2us.neutrals.parameters')

      call cfverw(iout,new_neutrals_version)

      write (iout,'(1x,a)') '&NEUTRALS'
      write (iout,'(1x,a,es9.2)') 'EIRENE_STEP_CPU = ',eirene_step_cpu
      write (iout,'(1x,a,es9.2)') 'EIRENE_STEP_DT  = ',eirene_step_dt
      write (iout,'(1x,a,i0)') 'EIRENE_MOD = ',eirene_mod

      write (iout,'(1x,a,i0)') 'NSTRAI     = ',m%nRc
      write (iout,'(1x,a,(1x,T14,10(2x,4a1)))') 
     .       'CRCSTRA = ',(ap,arcstra(istra),ap,',',istra=1,m%nRc)

      write(iout,'(1x,a,(1x,T14,10(i5,a1)))')
     .       'RCSTART =  ',(arcstart(istra),',',istra=1,m%nRc)
      write(iout,'(1x,a,(1x,T14,10(i5,a1)))')
     .       'RCEND  =  ',(arcend(istra),',',istra=1,m%nRc)
      write(iout,'(1x,a,(1x,T14,10(f6.2,a1,1x)))')
     .       'RCFE  =  ',(rcfe(istra),',',istra=1,m%nRc)
      write(iout,'(1x,a,(1x,T14,10(f6.2,a1,1x)))')
     .       'RCFI  =  ',(rcfi(istra),',',istra=1,m%nRc)
 
      do istra = 1, m%nRc
        write (iout,'(1x,a,i2,a,(10(f6.2,a1)))') 
     .    'recyc(0,',istra,')= ',(recyc(is,istra),',',is=0,ns-1)
      end do

      do istra = 1, m%nRc
        write (iout,'(1x,a,i2,a,(10(f6.2,a1)))') 
     .    'erecyc(0,',istra,')= ',(erecyc(is,istra),',',is=0,ns-1)
      end do

      do istra = 1, m%nRc
        if (crcstra(istra) /= 'V') cycle
        write (iout,'(1x,a,i2,a,es10.2,a1)') 
     .    'volrecstart(',istra,')= ',volrecstart(istra),','
      end do
      write (iout,'(1x,a,es10.2,a1)') 'volrecinc= ',volrecinc,','
      write (iout,'(1x,a,es10.2,a1)') 'volrecwt= ',volrecwt,','

      write (iout,'(1x,a,(1x,T7,10(i5,a1,1x)))') 
     .   'lsns= ',(abs(lsns(is)),',',is =1,max(1,sum(ksns(1:nstrai))))
      write (iout,'(1x,a,(1x,T7,10(i5,a1,1x)))') 
     .   'ksns= ',(ksns(is),',',is =1,nstrai)
      if (nnstsi > 0) then
        write (iout,'(1x,a,(1x,T7,10(i5,a1,1x)))') 
     .   'lcns= ',(lcns(is),',',is =1,nnstsi)
        write (iout,'(1x,a,(1x,T7,10(i5,a1,1x)))') 
     .   'ltns= ',(ltns(is),',',is =1,nnstsi)
      end if

      write (iout,'(1x,a,(10(i5,a1,1x)))') 
     .    'b2espcr= ',(b2espcr(is),',',is=0,ns-1)
      write (iout,'(1x,a,(10(i5,a1,1x)))') 
     .    'eb2spcr= ',(eb2spcr(is),',',is=1,nspecies)
      if (nnmoli > 0)
     .  write (iout,'(1x,a,(10(i3,a1,1x)))') 
     .    'lmolscl(1)= ',(lmolscl(is),',',is=1,nnmoli)
      do is=1,nnmoli
        write (iout,'(1x,a,i2,a,(10(i5,a1,1x)))') 
     .    'mlcmp(1,',is,')= ',(mlcmp(i,is),',',i=1,nnatmi)
      end do
      do is=1,nnatmi
        write (iout,'(1x,a,i2,a,(10(i3,a1,1x)))') 
     .    'lstrascl(1,',is,')= ',(lstrascl(i,is),',',i=1,m%nRc)
      end do
      write (iout,'(1x,a,i5,a1)') 'l_neutrad= ',l_neutrad,','
      write (iout,'(1x,a,es10.2,a1)') 'neut_scl_lim= ',neut_scl_lim,','

      write (iout,'(1x,a,i2,a,(10(es10.2,a1,1x)))') 
     .    'userfluxparm(',1,',1)= ',(userfluxparm(i,1),',',i=1,m%nRc)
      write (iout,'(1x,a,i2,a,(10(es10.2,a1,1x)))') 
     .    'userfluxparm(',1,',2)= ',(userfluxparm(i,2),',',i=1,m%nRc)

#ifndef NO_B2_CHEM_SPUT
      write (iout,'(1x,a,(1x,T25,10(f6.2,a1,1x)))') 
     .    'chemical_sputter_yield= ',
     .    (chemical_sputter_yield(i),',',i=1,nlimps)
      write (iout,'(1x,a,es10.2,a1)') 
     .    'fchar_chemical= ',fchar_chemical,','
      write (iout,'(1x,a,i5,a1)') 
     .    'igass_chemical= ',igass_chemical,','
      write (iout,'(1x,a,i5,a1)') 
     .    'itsput_chemical= ',itsput_chemical,','
      write (iout,'(1x,a,i5,a1)') 
     .    'issput_chemical= ',issput_chemical,','
#endif

      do i = 1, ntrack
        write (iout,'(1x,a,i2,a,(10(i5,a1)))') 
     .    'targsp(1,',i,')= ',(targsp(is,i),',',is=1,m%nRc)
      end do
      write (iout,'(1x,a,(1x,T10,10(l3,a1,1x)))') 
     .   'chemsp= ',(chemsp(is),',',is =1,nstrai)
      do istra = 1, m%nRc
        write (iout,'(1x,a,i2,a,(10(f6.2,a1)))') 
     .    'mrecyc(0,',istra,')= ',(mrecyc(is,istra),',',is=0,ns-1)
      end do
      do istra = 1, m%nRc
        write (iout,'(1x,a,i2,a,(10(f6.2,a1)))') 
     .    'rcion(0,',istra,')= ',(rcion(is,istra),',',is=0,ns-1)
      end do
      write (iout,'(1x,a,(1x,T11,10(es10.2,a1,1x)))') 
     .   'recyceir= ',(recyceir(is),',',is =1,nstrai)
      write (iout,'(1x,a,4(es10.2,a1,1x))') 
     .   'rf_neut= ',(rf_neut(i),',',i =1,4)
      do istra = 1, m%nRc
        write (iout,'(1x,a,i2,a,(10(es10.2,a1)))') 
     .   'phys_sput(0,',istra,')= ',(phys_sput(is,istra),',',is=0,ns-1)
      end do
      do istra = 1, m%nRc
        write (iout,'(1x,a,i2,a,(10(es10.2,a1)))') 
     .   'chem_sput(0,',istra,')= ',(chem_sput(is,istra),',',is=0,ns-1)
      end do
      write (iout,'(1x,a,(10(i5,a1,1x)))') 
     .    'species_start= ',(species_start(is),',',is=1,m%nRc)
      write (iout,'(1x,a,(10(i5,a1,1x)))') 
     .    'species_end= ',(species_end(is),',',is=1,m%nRc)

      write (iout,'(1x,5a)') 'NEUTRALS_FILENAME=',ap,
     .                       trim(neutrals_filename),ap,','
      write (iout,'(1x,a,es10.2,a1)') 
     .    'neutrals_time_mod= ',neutrals_time_mod,','
      write (iout,'(1x,a,es10.2,a1)') 
     .    'neutrals_time_switch= ',neutrals_time_switch,','
      write (iout,'(1x,a,i5,a1)') 
     .     'l_neutflux= ',l_neutflux,','
      write (iout,'(1x,a,(10(i5,a1,1x),/))') 
     .    'b2eatcr= ',(b2eatcr(is),',',is=0,ns-1)
      if (nnatmi > 0) then
        write (iout,'(1x,a,(10(i5,a1,1x)))') 
     .    'eb2atcr= ',(eb2atcr(is),',',is=1,nnatmi)
        write (iout,'(1x,a,(10(i5,a1,1x)))') 
     .    'latmscl= ',(latmscl(is),',',is=1,nnatmi)
      end if
      if (nnioni > 0) 
     .   write (iout,'(1x,a,(10(i5,a1,1x)))') 
     .    'lionscl= ',(lionscl(is),',',is=1,nnioni)

      if (nnatmi > 0) then
        do istra = 1, m%nRc
          write (iout,'(1x,a,i2,a,(10(es10.2,a1)))') 
     .   'gpfc(1,',istra,')= ',(gpfc(i,istra),',',i=1,nnatmi)
        end do
      end if
      write (iout,'(1x,a,i5,a1)') 
     .    'dbg_eir_mc= ',dbg_eir_mc,','
      write (iout,'(1x,a,(1x,T14,10(i5,a1,1x)))') 
     .    'debug_flags= ',(debug_flags(i),',',i=1,100)
      write (iout,'(1x,a,(1x,T14,10(i5,a1,1x)))') 
     .    'track_index= ',(track_index(is),',',is=0,ns-1)
      write (iout,'(1x,a,(1x,T18,10(l3,a1,1x)))') 
     .    'track_chem_sput= ',(track_chem_sput(is),',',is=1,ntrack)

      write (iout,'(1x,a,es10.2,a1)') 
     .    'chemical_erosion_redep_fac= ',chemical_erosion_redep_fac,','
      write (iout,'(1x,a,l3,a1)') 
     .    'chemical_erosion_be_fac=   ',chemical_erosion_be_fac,','
      write (iout,'(1x,a,es10.2,a1)') 
     .    'chemical_erosion_be_fac_a= ',chemical_erosion_be_fac_a,','
      write (iout,'(1x,a,es10.2,a1)') 
     .    'chemical_erosion_be_fac_b= ',chemical_erosion_be_fac_b,','
      write (iout,'(1x,a,es10.2,a1)') 
     .    'chemical_erosion_be_fac_c= ',chemical_erosion_be_fac_c,','

      write (iout,'(1x,a,i5,a1)') 
     .    'n_spcsrf= ',n_spcsrf,','
      if (n_spcsrf > 0) then
        write (iout,'(1x,a,(10(i5,a1,1x)))') 
     .    'l_spcsrf= ',(l_spcsrf(i),',',i=1,n_spcsrf)
        write (iout,'(1x,a,(10(a1,a,a1,a1,1x)))') 
     .    'sps_id= ',(ap,sps_id(i),ap,',',i=1,n_spcsrf)
        write (iout,'(1x,a,(10(i5,a1,1x)))') 
     .    'i_spcsrf= ',(i_spcsrf(i),',',i=1,n_spcsrf)
        write (iout,'(1x,a,(10(i5,a1,1x)))') 
     .    'j_spcsrf= ',(j_spcsrf(i),',',i=1,n_spcsrf)
        write (iout,'(1x,a,(10(es10.2,a1,1x)))') 
     .    'sps_absr= ',(sps_absr(i),',',i=1,n_spcsrf)
        write (iout,'(1x,a,(10(es10.2,a1,1x)))') 
     .    'sps_trno= ',(sps_trno(i),',',i=1,n_spcsrf)
        write (iout,'(1x,a,(10(es10.2,a1,1x)))') 
     .    'sps_trni= ',(sps_trni(i),',',i=1,n_spcsrf)
        write (iout,'(1x,a,(10(es10.2,a1,1x)))') 
     .    'sps_mtri= ',(sps_mtri(i),',',i=1,n_spcsrf)
        write (iout,'(1x,a,(10(es10.2,a1,1x)))') 
     .    'sps_tmpr= ',(sps_tmpr(i),',',i=1,n_spcsrf)
        write (iout,'(1x,a,(10(es10.2,a1,1x)))') 
     .    'sps_spph= ',(sps_spph(i),',',i=1,n_spcsrf)
        write (iout,'(1x,a,(10(es10.2,a1,1x)))') 
     .    'sps_spch= ',(sps_spch(i),',',i=1,n_spcsrf)
        write (iout,'(1x,a,(10(i5,a1,1x)))') 
     .    'sps_sgrp= ',(sps_sgrp(i),',',i=1,n_spcsrf)
        write (iout,'(1x,a,(10(a1,a,a1,a1,1x)))') 
     .    'sps_mtrl= ',(ap,sps_mtrl(i),ap,',',i=1,n_spcsrf)
      end if

      do istra=1,m%nRc
        write (iout,'(1x,a,i2,a,(1x,T14,10(i5,a1,1x)))') 
     .   'msns(1,',istra,')= ',(msns(i,istra),',',i =1,2)
      end do
      write (iout,'(1x,a,l3,a1)') 
     .    'write_nml_neut= ',write_nml_neut,','
      write (iout,'(1x,a,(10(l3,a1,1x)))') 
     .    'time_dep_puff= ',(time_dep_puff(i),',',i=1,m%nRc)
      write (iout,'(1x,a,(10(i5,a1,1x)))') 
     .    'ngpdata= ',(ngpdata(i),',',i=1,m%nRc)
      do istra=1,m%nRc
        if (ngpdata(istra) > 0) then
          write (iout,'(1x,a,i2,a,(1x,T14,10(i5,a1,1x)))') 
     .     'gpdata(1,1,',istra,')= ',
     .     (gpdata(i,1,istra),',',i=1,ngpdata(istra))
          write (iout,'(1x,a,i2,a,(1x,T14,10(i5,a1,1x)))') 
     .     'gpdata(1,2,',istra,')= ',
     .    (gpdata(i,2,istra),',',i=1,ngpdata(istra))
        end if
      end do

      write (iout,'(1x,a)') '/'
      close (iout)

      end subroutine write_rcparms

      end module b2mod_neutrals_namelist
