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
c*  mlcmp(i,j)    : number of atoms i in a molecule j
c*  lcns          : Eirene 'non-default surfaces' for the core boundary
c*  ncns          : number of these surfaces
c*  ltns          : Eirene 'non-default surfaces' for the targets
c*  ntns          : number of these surfaces
c*  dbg_eir_mc    : control of debugging output for eirene_mc
c*  neut_scl_lim  : maximum variation of neutral sources allowed during
c*                  source correction in eirene_mc
cank }
      use b2mod_types
      use b2mod_wall
      use b2mod_user_namelist
#     include "DIMENSIONS.F"
#     include "COUPLE/KOPPLDIM.F"
#     include "COUPLE/EIRSRT.F"

      integer nstraid, nsdmax, natmid, nmolid, nstsid
      parameter (NSTRAID=DEF_NSTRA, NSDMAX=DEF_NSD)
      parameter (natmid=DEF_NATM, nmolid=DEF_NMOL, nstsid=DEF_NSTS)
#ifdef ALLOCATE_AND_NAMELIST
      real (kind=R8), allocatable, save ::
     *   RECYC(:,:), ERECYC(:,:), RCION(:,:), recyceir(:),
     *   USERFLUXPARM(:,:),
     *   volrecstart(:), saved_volrec(:), 
     *   phys_sput(:,:), chem_sput(:,:),
     1   gpfc(:,:)
      integer, allocatable, save ::
     *   RCPOS(:), 
     *   RCSTART(:), 
     *   RCEND(:),
     *   TARGSP(:,:),
     *   species_start(:), species_end(:),
     *   b2espcr(:), eb2spcr(:), lstrascl(:,:), lmolscl(:),
     1   mlcmp(:,:), lcns(:), ltns(:), msns(:,:), ksns(:), lsns(:)
      logical, allocatable, save ::
     *   CHEMSP(:)
      character*1, allocatable, save ::
     *   CRCSTRA(:)
#else
      real (kind=R8), save :: RECYC(0:nsdmax-1,NSTRAID)
      real (kind=R8), save :: ERECYC(0:nsdmax-1,NSTRAID)
      real (kind=R8), save :: RCION(0:nsdmax-1,NSTRAID)
      real (kind=R8), save :: recyceir(NSTRAID)
      real (kind=R8), save :: USERFLUXPARM(NSTRAID,2)
      real (kind=R8), save :: volrecstart(NSTRAID)
      real (kind=R8), save :: saved_volrec(NSTRAID)
      real (kind=R8), save :: chem_sput(0:nsdmax-1,NSTRAID)
      real (kind=R8), save :: phys_sput(0:nsdmax-1,NSTRAID)
      integer, save :: RCPOS(NSTRAID)
      integer, save :: RCSTART(NSTRAID)
      integer, save :: RCEND(NSTRAID)
      integer, save :: TARGSP(NSTRAID,ntrack)
      logical, save :: CHEMSP(NSTRAID)
      integer, save :: species_start(NSTRAID), species_end(NSTRAID)
      integer, save :: b2espcr(0:nsdmax-1)
      integer, save :: eb2spcr(1:natmid)
      integer, save :: lstrascl(NSTRAID,0:natmid)
      character*1, save :: CRCSTRA(NSTRAID)
      integer, save :: lmolscl(nmolid), mlcmp(natmid,nmolid),
     ,  lcns(nstsid), ltns(nstsid), lsns(nlmpgs), ksns(NSTRAID), 
     ,  msns(2,NSTRAID)
      real (kind=R8), save :: gpfc(natmid,nstraid)
#endif

      integer, save :: NSTRAI, eirene_mod, ncns, ntns, 
     .     dbg_eir_mc, debug_flags(100)
      real (kind=R8), save :: RF_NEUT(4)
      real (kind=R8), save :: eirene_step_cpu, eirene_step_dt,
     .     volrecinc, volrecwt
      character*256, save :: neutrals_filename
      real (kind=R8), save :: neutrals_time_mod
      real (kind=R8), save :: neutrals_time_switch
      real (kind=R8), save :: neut_scl_lim

      real(kind=R8), save ::
     1 chemical_erosion_redep_fac,
     2 chemical_erosion_be_fac_a,
     3 chemical_erosion_be_fac_b,
     4 chemical_erosion_be_fac_c
      logical, save ::
     1 chemical_erosion_be_fac

csw 21feb2012
      integer, save ::
     . l_neutrad=0

      contains

      subroutine alloc_b2mod_neutrals(nsd,NSTRAID_in,sput_dst)
      use b2mod_elements
     & ,only : is_codes
      use b2mod_b2cmpa
      implicit none
      integer nsd, nstraid_in, sput_dst
#ifndef ALLOCATE_AND_NAMELIST
      external xertst
#endif
#     include "DIMENSIONS.F"
#     include "COUPLE/KOPPLDIM.F"
#     include "COUPLE/BRAEIR.F"
      integer is, isstart, isend
      integer is0, ispecies
      intrinsic nint

#ifdef ALLOCATE_AND_NAMELIST
      NSTRAID=NSTRAID_in
      allocate(RECYC(0:nsd-1,NSTRAID))
      allocate(ERECYC(0:nsd-1,NSTRAID))
      allocate(RCION(0:nsd-1,NSTRAID))
      allocate(recyceir(NSTRAID))
      allocate(USERFLUXPARM(NSTRAID,2))
      allocate(volrecstart(NSTRAID))
      allocate(chem_sput(0:nsd-1,NSTRAID))
      allocate(phys_sput(0:nsd-1,NSTRAID))
      allocate(RCPOS(NSTRAID))
      allocate(RCSTART(NSTRAID))
      allocate(RCEND(NSTRAID))
      allocate(TARGSP(NSTRAID,ntrack))
      allocate(CHEMSP(NSTRAID))
      allocate(CRCSTRA(NSTRAID))
      allocate(species_start(NSTRAID))
      allocate(species_end(NSTRAID))
      allocate(b2espcr(0:nsd-1))
      allocate(eb2spcr(1:natmid))
      allocate(lstrascl(NSTRAID,0:natmid))
      allocate(lmolscl(nmolid))
      allocate(gpfc(natmid,NSTRAID))
      allocate(mlcmp(natmid,nmolid))
      allocate(lcns(nstsid))
      allocate(ltns(nstsid))
      allocate(lsns(nlmpgs))
      allocate(ksns(NSTRAID))
      allocate(msns(2,NSTRAID))
#else
      call xertst (nsd.le.nsdmax, 'faulty parameter nsd')
      call xertst (nstraid_in.le.nstraid, 'faulty parameter nstraid_in')
#endif
      CRCSTRA=' '
      SPS_ID=' '
      NSTRAI=0
      RCPOS=-2
      RCSTART=-2
      RCEND=-2
      TARGSP=-1
      TARGSP(:,1)=sput_dst
      CHEMSP=.false.
cxpb The default recycling model is 100% particle recycling and 
cxpb 30% energy recycling for all ionised species except carbon.
      RECYC=0.0_R8
      ERECYC=0.0_R8
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
      chemical_sputter_yield=0.0
      fchar_chemical=0
      igass_chemical=0
      itsput_chemical=0
      issput_chemical=0
#endif
      eirene_step_cpu=0.0_R8
      eirene_step_dt=1e-3_R8
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
      eb2spcr=0
      lmolscl=0
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
        track_chem_sput(1)=.true.
      endif

      is0 = 0
      ispecies = 0
      do while (is0.lt.nsd)
        is = is0
        ispecies = ispecies + 1
        do while (zn(is).eq.zn(is0).and.
     &     abs(am(is)-am(is0)).lt.1.0e-3_R8.and.is.lt.nsd)
          is = is + 1
        enddo
        isstart = is0
        isend = is - 1
        b2espcr(isstart:isend)=ispecies
        eb2spcr(ispecies)=isstart
        is0 = is 
      enddo        
      write(*,*) 'b2espcr : ',b2espcr(0:nsd-1)

      return
      end subroutine alloc_b2mod_neutrals

      subroutine alloc_b2mod_neutrals_save(NSTRAID_in)
      integer nstraid_in
#ifndef ALLOCATE_AND_NAMELIST
      external xertst
#endif

#ifdef ALLOCATE_AND_NAMELIST
      NSTRAID=NSTRAID_in
      allocate(saved_volrec(NSTRAID))
#else
      call xertst (nstraid_in.le.nstraid, 'faulty parameter nstraid_in')
#endif
      return
      end subroutine alloc_b2mod_neutrals_save

      subroutine read_b2mod_neutrals_namelist(nx, ny, ns, dotest)
      use b2mod_wall
     & ,only : track_index, ntrack, track_species, track_chem_sput
      use b2mod_indirect
      use b2mod_elements
     & ,only : is_codes
      use b2mod_boundary_namelist
      use b2mod_b2cmpa
#ifdef USE_EIRENE
      use eirmod_extrab25
      use eirmod_braeir
#endif
      implicit none
      integer nx, ny, ns
      logical dotest
#     include <DIMENSIONS.F>
#     include "COUPLE/KOPPLDIM.F"
      external xertst
      integer is, ic, istra, i, j, k, n, ii(1), itrack
      logical hydrogen, helium, found
      intrinsic mod
      real (kind=R8) hlp

      namelist /neutrals/ nstrai,rcpos,rcstart,rcend,targsp,chemsp,
     4  recyc,erecyc,rcion,recyceir,userfluxparm,crcstra,rf_neut,
     .  phys_sput, chem_sput,
#ifndef NO_B2_CHEM_SPUT
#ifdef USE_EIRENE
     .  chemical_sputter_yield,fchar_chemical,igass_chemical,
     .  itsput_chemical,issput_chemical,
#endif
#endif
     .  eirene_step_cpu,eirene_step_dt,eirene_mod,
     .  volrecstart, volrecinc, volrecwt,
     .  species_start, species_end,
     5  neutrals_filename,neutrals_time_mod,neutrals_time_switch,
     ,  lstrascl,b2espcr,eb2spcr,lmolscl,mlcmp,
     ,  lcns,ltns,lsns,ksns,gpfc,dbg_eir_mc,debug_flags,
     ,  neut_scl_lim, track_index, track_chem_sput,
     *  chemical_erosion_redep_fac,
     *  chemical_erosion_be_fac,
     *  chemical_erosion_be_fac_a,
     *  chemical_erosion_be_fac_b,
     *  chemical_erosion_be_fac_c,
     6  n_spcsrf, l_spcsrf, sps_id, i_spcsrf, j_spcsrf, 
     7  sps_absr, sps_trno, sps_trni, sps_mtri, sps_tmpr, 
     8  sps_spph, sps_spch, sps_sgrp,
csw 21feb2012
     .  l_neutrad,
csw 04sep2013
     .  msns
      logical :: file_ok
      character*260 filename

      filename='../' // neutrals_filename
      inquire(file=filename,exist=file_ok)
      if(file_ok) then
        open(99,file=filename)
        write(*,*) 'Reading namelist neutrals from ',trim(filename)
        read(99,neutrals)
        write(*,*) 'Read namelist neutrals from ',trim(filename)
      else
        write(*,*) 'No ',trim(filename)
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

*   ..test input
      call xertst (0.le.minval(track_index).and.
     1 maxval(track_index).le.ntrack, 'faulty input track_index')
      call xertst (0.le.nstrai, 'faulty input nstrai')
      call xertst (0.le.ncns, 'faulty input ncns')
      call xertst (0.le.ntns, 'faulty input ntns')
      call xertst (0.0_R8.lt.rf_neut(1).and.0.0_R8.lt.rf_neut(2).and.
     .             0.0_R8.lt.rf_neut(3).and.0.0_R8.lt.rf_neut(4), 
     .  'faulty input rf_neut')
#ifndef NO_B2_CHEM_SPUT
#ifdef USE_EIRENE
      call xertst (0.le.issput_chemical,
     .  'faulty input issput_chemical')
#endif
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
      do istra = 1, nstrai
        if (crcstra(istra).eq.'N') then
          call xertst (1.le.rcpos(istra).and.rcpos(istra).le.ny+1,
     .     'faulty input rcpos (North)')
          call xertst (-1.le.rcstart(istra).and.rcstart(istra).le.nx+1,
     .     'faulty input rcstart (North)')
          call xertst (-1.le.rcend(istra).and.rcend(istra).le.nx+1,
     .     'faulty input rcend (North)')
        else if (crcstra(istra).eq.'S') then
          call xertst (-1.le.rcpos(istra).and.rcpos(istra).le.ny-1,
     .     'faulty input rcpos (South)')
          call xertst (-1.le.rcstart(istra).and.rcstart(istra).le.nx+1,
     .     'faulty input rcstart (South)')
          call xertst (-1.le.rcend(istra).and.rcend(istra).le.nx+1,
     .     'faulty input rcend (South)')
        else if (crcstra(istra).eq.'E') then
          call xertst (1.le.rcpos(istra).and.rcpos(istra).le.nx+1,
     .     'faulty input rcpos (East)')
          call xertst (-1.le.rcstart(istra).and.rcstart(istra).le.ny+1,
     .     'faulty input rcstart (East)')
          call xertst (-1.le.rcend(istra).and.rcend(istra).le.ny+1,
     .     'faulty input rcend (East)')
        else if (crcstra(istra).eq.'W') then
          call xertst (-1.le.rcpos(istra).and.rcpos(istra).le.nx-1,
     .     'faulty input rcpos (West)')
          call xertst (-1.le.rcstart(istra).and.rcstart(istra).le.ny+1,
     .     'faulty input rcstart (West)')
          call xertst (-1.le.rcend(istra).and.rcend(istra).le.ny+1,
     .     'faulty input rcend (West)')
        endif
        call xertst (rcstart(istra).le.rcend(istra),
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
          call xertst (0.0_R8.le.rcion(is,istra), 'faulty input rcion')
          call xertst (0.0_R8.le.phys_sput(is,istra),
     .     'faulty input phys_sput')
          call xertst (0.0_R8.le.chem_sput(is,istra),
     .     'faulty input chem_sput')
        enddo
        call xertst (0.0_R8.le.recyceir(istra), 'faulty input recyceir')
        do is = 1, natmid
          call xertst(0.0_R8.le.gpfc(is,istra), 'faulty input gpfc')
          call xertst(0.le.lstrascl(istra,is), 'faulty input lstrascl')
        enddo
      enddo
      do is = 0, ns-1
        call xertst(1.le.b2espcr(is).and.b2espcr(is).le.natmid,
     .   'faulty input b2espcr')
      enddo
      if (hydrogen) then
        call xertst(nint(zn(eb2spcr(j_h_at))).eq.1,
     .   'faulty input j_h_at')
      end if
      if (helium) call xertst(nint(zn(eb2spcr(j_he_at))).eq.2,
     . 'faulty input j_he_at')
      do is = 1, natmid
        call xertst(0.le.eb2spcr(is).and.eb2spcr(is).lt.ns,
     .   'faulty input eb2spcr')
      enddo
      do is = 1, nmolid
        call xertst(0.le.lmolscl(is), 'faulty input lmolscl')
        do i = 1, natmid
          call xertst(0.le.mlcmp(i,is), 'faulty input mlcmp')
        enddo
        if (hydrogen .and. user_namelist_used) call xertst
     .   (l_h_mol(is).eq.mlcmp(j_h_at,is), 
     .    'faulty inputs l_h_mol, mlcmp')
      enddo
      do is = 1, nstsid
        call xertst (0.le.lcns(is).and.lcns(is).le.nsts, 
     &   'faulty input lcns')
        call xertst (0.le.ltns(is).and.ltns(is).le.nsts, 
     &   'faulty input ltns')
        call xertst (lsns(is).le.nlim+nsts, 'faulty input lsns')
      enddo
      do is = 1, nstraid
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

      do j=1,nstsid !{
        if(lsns(j).gt.0 .and. lsns(j).le.nlim) lsns(j)=-lsns(j)
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
        do i=1,natmid !{
          hlp=hlp+gpfc(i,j)
        end do !}
        if(hlp.gt.0.0_R8) then !{
          do i=1,natmid !{
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
c.. Test consistency with region boundaries
      if (.not.dotest) return
      do ic=1, nncut
        found = mod(region(leftcut(ic)-1,-1,0),4).eq.1
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          found = rcend(istra).eq.leftcut(ic)-1 .and.
     &            crcstra(istra).eq.'S'
        end do
        call xertst (found, 
     .   'Mismatch with region boundaries (South) !')
        found = mod(region(leftcut(ic),-1,0),4).eq.1
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          found = rcstart(istra).eq.leftcut(ic) .and.
     &            crcstra(istra).eq.'S'
        end do
        call xertst (found, 
     .   'Mismatch with region boundaries (South) !')
        found = mod(region(rightcut(ic)-1,-1,0),4).eq.1
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          found = rcend(istra).eq.rightcut(ic)-1 .and.
     &            crcstra(istra).eq.'S'
        end do
        call xertst (found, 
     .   'Mismatch with region boundaries (South) !')
        found = mod(region(rightcut(ic),-1,0),4).eq.1
        istra = 0
        do while (istra.lt.nstrai .and. .not.found)
          istra = istra + 1
          found = rcstart(istra).eq.rightcut(ic) .and.
     &            crcstra(istra).eq.'S'
        end do
        call xertst (found, 
     .   'Mismatch with region boundaries (South) !')
      end do
      if (nncut.ge.2) then
        do ic = 1, nniso
          found = .false.
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            found = crcstra(istra).eq.'S' .and. 
     &       (rcend(istra).eq.nxiso1(ic)-1.or.
     &        rcend(istra).eq.nxiso1(ic))
          end do
          call xertst (found, 
     .     'Mismatch with isolating boundaries (South) !')
          found = .false.
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            found = crcstra(istra).eq.'N' .and. 
     &       (rcend(istra).eq.nxiso1(ic)-1.or.
     &        rcend(istra).eq.nxiso1(ic))
          end do
          call xertst (found, 
     .     'Mismatch with isolating boundaries (North)  !')
          found = .false.
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            found = crcstra(istra).eq.'S' .and. 
     &       (rcstart(istra).eq.nxiso2(ic).or.
     &        rcstart(istra).eq.nxiso2(ic)+1)
          end do
          call xertst (found, 
     .     'Mismatch with isolating boundaries (South) !')
          found = .false.
          istra = 0
          do while (istra.lt.nstrai .and. .not.found)
            istra = istra + 1
            found = crcstra(istra).eq.'N' .and. 
     &       (rcstart(istra).eq.nxiso2(ic).or.
     &        rcstart(istra).eq.nxiso2(ic)+1)
          end do
          call xertst (found, 
     .     'Mismatch with isolating boundaries (North)  !')
        end do
      end if
    
      return
      
      entry write_b2mod_neutrals_namelist

      write(*,neutrals)
      return

      end subroutine read_b2mod_neutrals_namelist

      subroutine read_b2mod_neutrals_save
      integer istra
      external xertst

      namelist /neutrals_save/ saved_volrec
      logical :: file_ok

      saved_volrec(:)=0.0_R8

      inquire(file='../b2.neutrals_save.parameters',exist=file_ok)
      if(file_ok) then
        open(99,file='../b2.neutrals_save.parameters')
        read(99,neutrals_save)
        close(99)
      else
        write(*,*) 'No ../b2.neutrals_save.parameters'
      endif
      return

*   ..test input
      do istra = 1, nstrai
        call xertst (0.0_R8.le.saved_volrec(istra), 
     .   'faulty input saved_volrec')
      enddo
      
      end subroutine read_b2mod_neutrals_save

      subroutine write_b2mod_neutrals_save

      namelist /neutrals_save/ saved_volrec

      write(*,neutrals_save)
      return

      end subroutine write_b2mod_neutrals_save

      end module b2mod_neutrals_namelist
