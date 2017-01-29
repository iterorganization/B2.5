      module b2mod_neutrals_namelist

c  version : 27.07.2002 00:00

cank-20020724 {
c* Parameters controlling the source scaling procedure are added:
c*  b2espcr(i)    : Eirene atomic species corresponding to B2 species i
c*  eb2spcr(i)    : B2 atomic species corresponding to Eirene atomic species i
c*  lstrascl(i,j) : scale atomic species j in stratum i with this atomic
c*                  species (0 means no scaling). Species 0 is reserved
c*                  for electron energy (TBD). 
c*  gpfc(i,j)     : fraction of nuclei of kind i (Eirene atoms) in one
c*                  particle puffed in stratum number j
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
      use b2mod_user_namelist
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_extrab25
      use eirmod_parmmod
#else
      use parmmod   ! IGNORE
      use braeir    ! IGNORE
      use eirdiag   ! IGNORE
      use ceirsrt   ! IGNORE
#endif
#endif
      implicit none
#     include "DIMENSIONS.F"
#ifndef B25_EIRENE
#     include "COUPLE/KOPPLDIM.F"
#     include "COUPLE/EIRDIAG.F"
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
     1   gpfc(:,:), gpdata(:,:,:)
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
     *   RC_LIST_SIZE(:),
     *   RC_LIST_X(:,:),
     *   RC_LIST_Y(:,:),
     *   TARGSP(:,:),
     *   species_start(:), species_end(:),
     *   b2espcr(:), eb2spcr(:), lstrascl(:,:), lmolscl(:), lionscl(:),
     1   mlcmp(:,:), lcns(:), ltns(:), msns(:,:), ksns(:), lsns(:),
     2   ngpdata(:)
      logical, allocatable, save ::
     *   CHEMSP(:), TIME_DEP_PUFF(:)
      character*1, allocatable, save ::
     *   CRCSTRA(:)
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
#ifndef NO_B2_CHEM_SPUT
#ifndef SOLPS4_3
      real (kind=R8), save :: 
     ,  chemical_sputter_yield(0:DEF_NLIM+DEF_NSTS)
#endif
#endif
      integer, save :: RCPOS(NSTRAID)
      integer, save :: RCSTART(NSTRAID)
      integer, save :: RCEND(NSTRAID)
      integer, save :: RC_LIST_SIZE(NSTRAID)
      integer, save :: RC_LIST_X(2*(DEF_NXD+DEF_NYD),NSTRAID)
      integer, save :: RC_LIST_Y(2*(DEF_NXD+DEF_NYD),NSTRAID)
      integer, save :: TARGSP(NSTRAID,ntrack)
      logical, save :: CHEMSP(NSTRAID)
      logical, save :: TIME_DEP_PUFF(NSTRAID)
      integer, save :: ngpdata(NSTRAID)
      integer, save :: species_start(NSTRAID), species_end(NSTRAID)
      integer, save :: b2espcr(0:nsdmax-1)
      integer, save :: eb2spcr(1:natmid)
      integer, save :: lstrascl(NSTRAID,0:natmid)
      character*1, save :: CRCSTRA(NSTRAID)
      integer, save :: lmolscl(nmolid), mlcmp(natmid,nmolid),
     ,  lionscl(nionid),
     ,  lcns(nstsid), ltns(nstsid), lsns(nstraid), ksns(NSTRAID), 
     ,  msns(2,NSTRAID)
      real (kind=R8), save :: gpfc(natmid,nstraid)
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
     1  rc_list_size, rc_list_x, rc_list_y,
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
     ,  b2espcr, eb2spcr, lmolscl, mlcmp, lionscl,
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
      allocate(RC_LIST_SIZE(NSTRAID))
      allocate(RC_LIST_X(2*(NXD+NYD),NSTRAID))
      allocate(RC_LIST_Y(2*(NXD+NYD),NSTRAID))
      allocate(TARGSP(NSTRAID,ntrack))
      allocate(CHEMSP(NSTRAID))
      allocate(TIME_DEP_PUFF(NSTRAID))
      allocate(NGPDATA(NSTRAID))
      allocate(GPDATA(100,2,NSTRAID))
      allocate(CRCSTRA(NSTRAID))
      allocate(species_start(NSTRAID))
      allocate(species_end(NSTRAID))
      allocate(b2espcr(0:nsd-1))
      allocate(eb2spcr(1:natm_in))
      allocate(lstrascl(nstra_in,0:natm_in))
      allocate(lmolscl(nmol_in))
      allocate(lionscl(nion_in))
      allocate(gpfc(natm_in,nstra_in))
      allocate(mlcmp(natm_in,nmol_in))
      allocate(lcns(nsts_in))
      allocate(ltns(nsts_in))
      allocate(lsns(nstra_in))
      allocate(ksns(nstra_in))
      allocate(msns(2,nstra_in))
#ifndef NO_B2_CHEM_SPUT
      allocate(chemical_sputter_yield(0:nlimps))
#endif
#else
      call xertst (nxd.le.DEF_NXD, 'faulty parameter nxd')
      call xertst (nyd.le.DEF_NYD, 'faulty parameter nyd')
      call xertst (nsd.le.DEF_NSD, 'faulty parameter nsd')
      call xertst (natm_in.le.NATMID, 'faulty parameter natm_in')
      call xertst (nmol_in.le.NMOLID, 'faulty parameter nmol_in')
      call xertst (nion_in.le.NIONID, 'faulty parameter nion_in')
      call xertst (nsts_in.le.NSTSID, 'faulty parameter nsts_in')
      call xertst (nstra_in.le.NSTRAID, 'faulty parameter nstraid_in')
#endif
      call xertst (natm_in.ge.nspecies, 'faulty parameter natm_in')
      CRCSTRA=' '
      SPS_ID=' '
      SPS_MTRL=' '
      NSTRAI=0
      RCPOS=-2
      RCSTART=-2
      RCEND=-2
      RC_LIST_SIZE=0
      RC_LIST_X=-2
      RC_LIST_Y=-2
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
      b2espcr=0
      eb2spcr=-1
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
        b2espcr(isstart:isend)=ispecies
        eb2spcr(ispecies)=isstart
        is0 = is 
      enddo        
      write(hlp_frm,'(a6,i3,a3)') '(1x,a,',nsd,'i4)'
      write(*,hlp_frm) 'b2espcr : ',b2espcr(0:nsd-1)

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
      deallocate(RC_LIST_SIZE)
      deallocate(RC_LIST_X)
      deallocate(RC_LIST_Y)
      deallocate(TARGSP)
      deallocate(CHEMSP)
      deallocate(TIME_DEP_PUFF)
      deallocate(NGPDATA)
      deallocate(GPDATA)
      deallocate(CRCSTRA)
      deallocate(species_start)
      deallocate(species_end)
      deallocate(b2espcr)
      deallocate(eb2spcr)
      deallocate(lstrascl)
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

      subroutine read_neutrals_namelist(nx, ny, ns, dotest)
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
      logical hydrogen, helium, neon, found, side
      intrinsic mod
      real (kind=R8) hlp
      logical :: file_ok      
      logical, allocatable :: old_style(:)
      character*260 filename

      filename='../' // neutrals_filename
      inquire(file=filename,exist=file_ok)
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
        do is = 1, nspecies
          if (lstrascl(istra,is).lt.0) lstrascl(istra,is)=is
          if (lstrascl(istra,is).gt.nspecies) lstrascl(istra,is)=is
        end do
      end do

*   ..set additional variables
      hydrogen = .false.
      is = 0
      do while (.not.hydrogen .and. is.lt.ns)
        hydrogen = is_neutral(is) .and. nint(zn(is)).eq.1
        if (hydrogen) j_h_at = b2espcr(is)
        is = is + 1
      end do
      helium = .false.
      is = 0
      do while (.not.helium .and. is.lt.ns)
        helium = is_neutral(is) .and. nint(zn(is)).eq.2
        if (helium) j_he_at = b2espcr(is)
        is = is + 1
      end do
      neon = .false.
      is = 0
      do while (.not.neon .and. is.lt.ns)
        neon = is_neutral(is) .and. nint(zn(is)).eq.10
        if (neon) j_ne_at = b2espcr(is)
        is = is + 1
      end do
      do is = 1, nstrat
        if (crcstra(is).ne.'N' .and. crcstra(is).ne.'S' .and.
     &      crcstra(is).ne.'W' .and. crcstra(is).ne.'E') targsp(is,:)=-1
      end do

*   ..test input
      call xertst (0.le.minval(track_index).and.
     1 maxval(track_index).le.ntrack, 'faulty input track_index')
      call xertst (0.le.nstrai, 'faulty input nstrai')
      call xertst (0.le.nstrat, 'faulty input nstrat')
      call xertst (0.le.ncns, 'faulty input ncns')
      call xertst (0.le.ntns, 'faulty input ntns')
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
        do is = 1, nspecies
          call xertst(0.0_R8.le.gpfc(is,istra), 'faulty input gpfc')
          call xertst(0.le.lstrascl(istra,is) .and.
     .                lstrascl(istra,is).le.nspecies, 
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
      enddo
      if (hydrogen) then
        call xertst(nint(zn(eb2spcr(j_h_at))).eq.1,
     .   'faulty input j_h_at')
      end if
      if (helium) call xertst(nint(zn(eb2spcr(j_he_at))).eq.2,
     . 'faulty input j_he_at')
      if (neon) call xertst(nint(zn(eb2spcr(j_ne_at))).eq.10,
     . 'faulty input j_ne_at')
      do is = 1, nspecies
        call xertst(0.le.eb2spcr(is).and.eb2spcr(is).lt.ns,
     .   'faulty input eb2spcr')
      enddo
      do is = 1, nnmoli
        do i = 1, nspecies
          call xertst(0.le.mlcmp(i,is), 'faulty input mlcmp')
        enddo
#ifdef B25_EIRENE
        call xertst(0.le.lmolscl(is), 'faulty input lmolscl')
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
        if (maxval(mlcmp(1:nspecies,is)).eq.0) then
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
        if (hydrogen .and. user_namelist_used) call xertst
     .   (l_h_mol(is).eq.mlcmp(j_h_at,is), 
     .    'faulty inputs l_h_mol, mlcmp')
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
        call xertst (lsns(is).le.nlim+nsts, 'faulty input lsns')
        call xertst (lsns(is).ge.-nsts, 'faulty input lsns')
        call xertst (0.le.ksns(is), 'faulty input ksns')
      enddo
      
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
        do i=1,nspecies !{
          hlp=hlp+gpfc(i,j)
        end do !}
        if(hlp.gt.0.0_R8) then !{
          do i=1,nspecies !{
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
      write(*,*) 'track_species : Bulk ', track_species(1:ntrack)
c.. Test consistency with region boundaries
      if (.not.dotest) return
      do ic=1, nncut
        found = mod(region(leftcut(ic)-1,-1,0),4).eq.1
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
        found = mod(region(leftcut(ic),-1,0),4).eq.1
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          found = rcstart(istra).eq.leftcut(ic) .and.
     &            crcstra(istra).eq.'S' .and. old_style(istra)
        end do
        call xertst (found, 
     .   'Mismatch with neutrals boundaries (South) !')
        found = mod(region(rightcut(ic)-1,-1,0),4).eq.1
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
        found = mod(region(rightcut(ic),-1,0),4).eq.1
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
      end subroutine read_neutrals_namelist

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
        filename = '../b2.neutrals_save.parameters'
      else
        filename = 'b2.neutrals_save.parameters.'//edition
      endif
      inquire(file=trim(filename),exist=file_ok)
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
        filename = '../b2.neutrals_save.parameters'
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

      end module b2mod_neutrals_namelist
