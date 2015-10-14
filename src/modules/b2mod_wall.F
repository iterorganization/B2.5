      module b2mod_wall
      use b2mod_types
      use b2mod_layer
      implicit none

      integer ndepth, nwall_segments, nsegments, nline
      integer, save :: ntargsp
      parameter(
     *  ndepth=201,
     *  nline=15)
cxpb  Namelist data
      integer, save :: nbnd, ndepth_nml
#ifdef ALLOCATE_AND_NAMELIST
      integer, allocatable, save :: imapx(:), imapy(:)
      logical, allocatable :: inertial_cooling(:)
      real (kind=R8), allocatable, save ::
     *  target_temp(:,:),
     *  target_depth(:,:),
     *  target_width(:,:),
     *  target_vol(:,:),
     *  targetcr(:,:,:),
     *  targetcz(:,:,:),
     *  backplate_temp(:),
     *  plate_thickness(:),
     *  plate_time_factor(:),
     *  coating_thickness(:),
     *  deposition(:,:),
     *  monolayer_deposition(:,:),
     *  erosion(:,:),
     *  monolayer_erosion(:,:),
     *  plate_time(:),
     *  plate_area(:),
     *  chemical_sputtering(:,:),
     *  physical_sputtering(:,:),
     *  physical_sputtering_energy(:,:),
     *  RES_sputtering(:,:),
     *  thermal_evaporation(:,:),
     *  backscattering(:,:),
     *  backscattering_energy(:,:)

cxpb Data file material properties
      character*6, allocatable, save ::
     *  surface_material_name(:),
     *  coating_material_name(:),
     *  bulk_material_name(:)

      integer, allocatable, save ::
     *  xymap(:,:)

      integer, allocatable, save :: track_index(:)
#else
#     include "DIMENSIONS.F"
      integer, save ::
     *  imapx(nwall), imapy(nwall)
      integer, save :: plate_model, plate_option
      logical inertial_cooling(nwall)
      real (kind=R8), save ::
     *  plate_temp, plate_thick,
     *  target_temp(nwall,ndepth),
     *  target_depth(nwall,ndepth),
     *  target_width(nwall,ndepth),
     *  target_vol(nwall,ndepth),
     *  targetcr(nwall,ndepth,2),
     *  targetcz(nwall,ndepth,2),
     *  backplate_temp(nwall),
     *  plate_thickness(nwall),
     *  plate_time_factor(nwall),
     *  coating_thickness(nwall),
     *  deposition(nwall,ntrack),
     *  monolayer_deposition(nwall,ntrack),
     *  erosion(nwall,ntrack),
     *  monolayer_erosion(nwall,ntrack),
     *  plate_time(nwall),
     *  plate_area(nwall),
     *  chemical_sputtering
     .   (nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  physical_sputtering
     .   (nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  physical_sputtering_energy
     .   (nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  RES_sputtering(nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  thermal_evaporation
     .   (nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  backscattering(nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  backscattering_energy
     .   (nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION)

      character*6, save :: surface_material_name(nwall)
      character*6, save :: coating_material_name(nwall)
      character*6, save :: bulk_material_name(nwall)

      integer xymap(-1:DEF_NXD,-1:DEF_NYD)

      integer, save :: track_index(0:DEF_NSD-1)
#endif

      integer, allocatable, save :: imapw(:),
     *  imapreg(:), wall_begin(:), wall_end(:)

cxpb Data file material properties
      character*32, allocatable, save :: surface_long_name(:)
      integer, allocatable, save ::
     *  matsurf_nconstituents(:)
      integer, allocatable, save ::
     *  matsurf_nzconstituents(:,:),
     *  matsurf_nrelconstituents(:,:)
      real (kind=R8), allocatable, save ::
     *  matsurf_chemsput_erel(:),
     *  matsurf_blackbody_emiss(:),
     *  matsurf_vapprs_A(:),
     *  matsurf_vapprs_B(:),
     *  matsurf_stick_coef(:),
     *  matsurf_vapor_temp(:),
     *  matsurf_heat_sublim(:)

      character*32, allocatable, save :: coating_long_name(:)
      real (kind=R8), allocatable, save ::
     *  coating_target_density(:),
     *  coating_Kappa_KA(:),
     *  coating_Kappa_KB(:),
     *  coating_Kappa_A1(:),
     *  coating_Kappa_A2(:),
     *  coating_Kappa_A3(:),
     *  coating_Kappa_A4(:),
     *  coating_Kappa_perp_A1(:),
     *  coating_Kappa_perp_A2(:),
     *  coating_Kappa_perp_A3(:),
     *  coating_Kappa_perp_A4(:),
     *  coating_Cp_B1(:),
     *  coating_Cp_B2(:),
     *  coating_Cp_B3(:)

      character*32, allocatable, save :: bulk_long_name(:)
      real (kind=R8), allocatable, save ::
     *  matbulk_target_density(:),
     *  matbulk_Kappa_KA(:),
     *  matbulk_Kappa_KB(:),
     *  matbulk_Kappa_A1(:),
     *  matbulk_Kappa_A2(:),
     *  matbulk_Kappa_A3(:),
     *  matbulk_Kappa_A4(:),
     *  matbulk_Kappa_perp_A1(:),
     *  matbulk_Kappa_perp_A2(:),
     *  matbulk_Kappa_perp_A3(:),
     *  matbulk_Kappa_perp_A4(:),
     *  matbulk_Cp_B1(:),
     *  matbulk_Cp_B2(:),
     *  matbulk_Cp_B3(:)

      real(kind=R8), save ::         ! data from incatm file
     *  atomic_mass(92),             ! in amu
     1  surface_binding_energy(92),  ! in eV
     2  target_density(92)           ! in g/cm3

      integer, allocatable, save :: itrack_0(:)
      character*2, save :: track_species(ntrack)
      character*2, save :: track_species_save(ntrack)
      logical, save :: track_chem_sput(ntrack)

      real (kind=R8), allocatable :: r1(:), r2(:), z1(:), z2(:)

      namelist /wall_save/ ndepth_nml, track_species,
     *  imapx, imapy, xymap, surface_material_name, bulk_material_name,
     *  target_temp, backplate_temp, plate_thickness, plate_time_factor,
     *  deposition, erosion, chemical_sputtering, physical_sputtering,
     *  physical_sputtering_energy, RES_sputtering, thermal_evaporation,
     *  backscattering, backscattering_energy, plate_time, plate_area,
     *  monolayer_deposition, monolayer_erosion,
     *  layer_nconstituents, layer_nzconstituents,
     *  layer_nrelconstituents, layer_thermal_contact,
     *  coating_thickness, coating_material_name,
     *  layer_alloys, inertial_cooling

cxpb Mapping arrays:
cxpb A wall element numbered iwall can be found in the B2 grid at the position
cxpb   ( imapx(iwall) , imapy(iwall) )
cxpb Reversely, the B2 boundary cell ( ix, iy ) corresponds to the boundary
cxpb element: xymap(ix,iy)

      data plate_model/0/, plate_option/3/,
     .     plate_temp/300.0_R8/, plate_thick/0.00_R8/

      contains

      subroutine alloc_b2mod_wall_init(nx,ny,nncut)
      implicit none
      integer nx, ny, nncut
#ifndef ALLOCATE_AND_NAMELIST
      logical errdim
      external xertst
#endif
      external subini, subend

      call subini('b2mod_wall_init')
cxpb  Geometric information
#ifdef ALLOCATE_AND_NAMELIST
      allocate(imapx(nwall), imapy(nwall))
      allocate(plate_area(nwall))
      allocate(xymap(-1:nx,-1:ny))
#else
      errdim=.false.
      if(nx.gt.DEF_NXD) then
        write(*,*) 'alloc_b2mod_wall_init: nxd > def_nxd',nx,DEF_NXD
        errdim=.true.
      endif
      if(ny.gt.DEF_NYD) then
        write(*,*) 'alloc_b2mod_wall_init: nyd > def_nyd',ny,DEF_NYD
        errdim=.true.
      endif
      call xertst (.not.errdim, 'faulty def_nxd, def_nyd')
#endif

      nwall_segments=4*(2*max(nncut,1)+1)
      allocate(imapw(nwall), imapreg(nwall),
     *  wall_begin(nwall_segments), wall_end(nwall_segments))
      allocate(r1(nwall), r2(nwall), z1(nwall), z2(nwall))
      allocate(itrack_0(nwall))

      call subend()
      return
      end subroutine alloc_b2mod_wall_init

      subroutine alloc_b2mod_wall(ns,natm,nmol,nion)
      implicit none
      integer ns, natm, nmol, nion
      external subini, subend, xertst

      call subini('b2mod_wall')
#ifdef ALLOCATE_AND_NAMELIST
cxpb  Namelist data
      allocate(inertial_cooling(nwall))
      allocate(
     *  target_temp(nwall,ndepth),
     *  target_depth(nwall,ndepth),
     *  target_width(nwall,ndepth),
     *  target_vol(nwall,ndepth),
     *  targetcr(nwall,ndepth,2),
     *  targetcz(nwall,ndepth,2),
     *  backplate_temp(nwall),
     *  plate_thickness(nwall),
     *  plate_time_factor(nwall),
     *  coating_thickness(nwall),
     *  deposition(nwall,ntrack),
     *  monolayer_deposition(nwall,ntrack),
     *  erosion(nwall,ntrack),
     *  monolayer_erosion(nwall,ntrack),
     *  plate_time(nwall),
     *  chemical_sputtering(nwall,0:ns-1+natm+nmol+nion),
     *  physical_sputtering(nwall,0:ns-1+natm+nmol+nion),
     *  physical_sputtering_energy(nwall,0:ns-1+natm+nmol+nion),
     *  RES_sputtering(nwall,0:ns-1+natm+nmol+nion),
     *  thermal_evaporation(nwall,0:ns-1+natm+nmol+nion),
     *  backscattering(nwall,0:ns-1+natm+nmol+nion),
     *  backscattering_energy(nwall,0:ns-1+natm+nmol+nion))

      allocate(
     *  surface_material_name(nwall),
     *  coating_material_name(nwall),
     *  bulk_material_name(nwall))

      allocate(track_index(0:ns-1))
#endif

cxpb Data file material properties
      allocate(surface_long_name(nwall))
      allocate(matsurf_nconstituents(nwall))
      allocate(matsurf_chemsput_erel(nwall),
     *  matsurf_blackbody_emiss(nwall),
     *  matsurf_vapprs_A(nwall),
     *  matsurf_vapprs_B(nwall),
     *  matsurf_stick_coef(nwall),
     *  matsurf_vapor_temp(nwall),
     *  matsurf_heat_sublim(nwall))
      
      allocate (coating_long_name(nwall))
      allocate(
     *  coating_target_density(nwall),
     *  coating_Kappa_KA(nwall),
     *  coating_Kappa_KB(nwall),
     *  coating_Kappa_A1(nwall),
     *  coating_Kappa_A2(nwall),
     *  coating_Kappa_A3(nwall),
     *  coating_Kappa_A4(nwall),
     *  coating_Kappa_perp_A1(nwall),
     *  coating_Kappa_perp_A2(nwall),
     *  coating_Kappa_perp_A3(nwall),
     *  coating_Kappa_perp_A4(nwall),
     *  coating_Cp_B1(nwall),
     *  coating_Cp_B2(nwall),
     *  coating_Cp_B3(nwall))

      allocate(bulk_long_name(nwall))
      allocate(
     *  matbulk_target_density(nwall),
     *  matbulk_Kappa_KA(nwall),
     *  matbulk_Kappa_KB(nwall),
     *  matbulk_Kappa_A1(nwall),
     *  matbulk_Kappa_A2(nwall),
     *  matbulk_Kappa_A3(nwall),
     *  matbulk_Kappa_A4(nwall),
     *  matbulk_Kappa_perp_A1(nwall),
     *  matbulk_Kappa_perp_A2(nwall),
     *  matbulk_Kappa_perp_A3(nwall),
     *  matbulk_Kappa_perp_A4(nwall),
     *  matbulk_Cp_B1(nwall),
     *  matbulk_Cp_B2(nwall),
     *  matbulk_Cp_B3(nwall))

      call xertst(1.le.ndepth .and. ndepth.le.1001,
     . 'invalid dimension ndepth')

      ndepth_nml = 1
      inertial_cooling = .false.
      backplate_temp = 0.0_R8
      plate_thickness = 0.0_R8
      coating_thickness = 0.0_R8
      plate_time_factor = 0.0_R8
      surface_material_name = ' '
      coating_material_name = ' '
      bulk_material_name = ' '
      chemical_sputtering = 0.0_R8
      physical_sputtering = 0.0_R8
      physical_sputtering_energy = 0.0_R8
      RES_sputtering = 0.0_R8
      thermal_evaporation = 0.0_R8
      backscattering = 0.0_R8
      backscattering_energy = 0.0_R8
      deposition = 0.0_R8
      monolayer_deposition = 0.0_R8
      erosion = 0.0_R8
      monolayer_erosion = 0.0_R8
      target_temp = 0.0_R8
      plate_time = 0.0_R8

      call subend()
      return

      end subroutine alloc_b2mod_wall

      subroutine read_b2mod_wall_namelist(nread,ns,edition,found)
      implicit none
      integer nread, ns
      character*4 edition
      logical found
      character*256 filename
      logical file_ok
      external xertst

      if (edition.eq.'    ') then
        filename = '../b2.wall_save.parameters'
      else
        filename = 'b2.wall_save.parameters.'//edition
      endif
      track_species_save=track_species
      inquire(file=trim(filename),exist=file_ok)
      if(file_ok) then
        open(nread,file=trim(filename))
        write(*,*) 'Reading ', trim(filename)
        read(nread,wall_save)
        write(*,*) 'Read ', trim(filename)
        close(nread)
      else
        write(6,*) 'No ',trim(filename)
      endif
      found =
     .  layer_thermal_contact.ne.1.0_R8 .or.
     .  maxval(target_temp).gt.plate_temp .or.
     .  maxval(thermal_evaporation).gt.0.0_R8 .or.
     .  maxval(physical_sputtering).gt.0.0_R8 .or.
     .  maxval(chemical_sputtering).gt.0.0_R8 .or.
     .  maxval(RES_sputtering).gt.0.0_R8 .or.
     .  maxval(backscattering).gt.0.0_R8
      call xertst(ndepth_nml.le.ndepth .and. ndepth_nml.ge.1,
     . 'Invalid array size NDEPTH_NML!')
      return
      end subroutine read_b2mod_wall_namelist

      subroutine write_b2mod_wall_namelist(nwrite,nx,ny,ns,edition)
      implicit none
      integer nwrite, nx, ny, ns
      character*4 edition
      integer i, j, ix, iy, is, ibnd, itrack
      character*256 filename
      character*400 line
      character*17 numstr
      logical found, todo

      if (edition.eq.'    ') then
        filename = '../b2.wall_save.parameters'
      else
        filename = 'b2.wall_save.parameters.'//edition
      endif
      open(nwrite,file=trim(filename))
      write(nwrite,'(a10)') '&WALL_SAVE'
      do itrack=1,nalloys
        if(layer_alloys(itrack).ne.' ')
     .    write(nwrite,'(a13,i1,a4,a6,a2)')
     .    'LAYER_ALLOYS(',itrack,')= ''',
     .     layer_alloys(itrack),''','
      enddo
      write(nwrite,'(a11,i4,a1)') 'NDEPTH_NML=',ndepth_nml,','
      do itrack=1,ntargsp
        write(nwrite,'(a14,i2,a3,a2,a2)')
     .   'TRACK_SPECIES(',itrack,')=''',track_species(itrack),''','
      enddo
      if (layer_thermal_contact.ne.1.0_R8)
     > write(nwrite,'(a22,1es15.8,a1)')
     > 'LAYER_THERMAL_CONTACT=',layer_thermal_contact,','
      if (.not.(plate_model.ge.1.and.ndepth_nml.gt.1)) then
        do ibnd = 1, nbnd
          if (target_temp(ibnd,1).gt.plate_temp)
     >     write(nwrite,'(a12,i4,a4,1x,1es15.8,a1)')
     .     'TARGET_TEMP(',ibnd,',1)=',target_temp(ibnd,1),','
        enddo
      endif
      do is = 0, ns-1
        if (maxval(thermal_evaporation(1:nbnd,is)).gt.0.0_R8) then
          do i = 1, (nbnd-1)/nline+1
            if (i.eq.1) then
              write(line,'(a22,i3,a2)') 'THERMAL_EVAPORATION(1,',is,')='
            else
              line = ''
            endif
            do ibnd = (i-1)*nline+1, min(i*nline,nbnd)
              write(numstr,'(1x,1es15.8,a1)')
     .         thermal_evaporation(ibnd,is),','
              line = trim(line)//numstr
            enddo
            write(nwrite,'(a)') trim(line)
          enddo
        endif
        if (maxval(physical_sputtering(1:nbnd,is)).gt.0.0_R8) then
          do i = 1, (nbnd-1)/nline+1
            if (i.eq.1) then
              write(line,'(a22,i3,a2)') 'PHYSICAL_SPUTTERING(1,',is,')='
            else
              line = ''
            endif
            do ibnd = (i-1)*nline+1, min(i*nline,nbnd)
              write(numstr,'(1x,1e15.8,a1)')
     .         physical_sputtering(ibnd,is),','
              line = trim(line)//numstr
             enddo
           write(nwrite,'(a)') trim(line)
          enddo
          do i = 1, (nbnd-1)/nline+1
            if (i.eq.1) then
              write(line,'(a29,i3,a2)')
     .         'PHYSICAL_SPUTTERING_ENERGY(1,',is,')='
            else
              line = ''
            endif
            do ibnd = (i-1)*nline+1, min(i*nline,nbnd)
              write(numstr,'(1x,1es15.8,a1)')
     .         physical_sputtering_energy(ibnd,is),','
              line = trim(line)//numstr
            enddo
            write(nwrite,'(a)') trim(line)
          enddo
        endif
        if (maxval(chemical_sputtering(1:nbnd,is)).gt.0.0_R8) then
          do i = 1, (nbnd-1)/nline+1
            if (i.eq.1) then
              write(line,'(a22,i3,a2)')
     .         'CHEMICAL_SPUTTERING(1,',is,')='
            else
              line = ''
            endif
            do ibnd = (i-1)*nline+1, min(i*nline,nbnd)
              write(numstr,'(1x,1es15.8,a1)')
     .         chemical_sputtering(ibnd,is),','
              line = trim(line)//numstr
             enddo
           write(nwrite,'(a)') trim(line)
          enddo
        endif
        if (maxval(RES_sputtering(1:nbnd,is)).gt.0.0_R8) then
          do i = 1, (nbnd-1)/nline+1
            if (i.eq.1) then
              write(line,'(a17,i3,a2)')
     .         'RES_SPUTTERING(1,',is,')='
            else
              line = ''
            endif
            do ibnd = (i-1)*nline+1, min(i*nline,nbnd)
              write(numstr,'(1x,1es15.8,a1)')
     .         RES_sputtering(ibnd,is),','
              line = trim(line)//numstr
            enddo
            write(nwrite,'(a)') trim(line)
          enddo
        endif
        if (maxval(backscattering(1:nbnd,is)).gt.0.0_R8) then
          do i = 1, (nbnd-1)/nline+1
            if (i.eq.1) then
              write(line,'(a17,i3,a2)')
     .         'BACKSCATTERING(1,',is,')='
            else
              line = ''
            endif
            do ibnd = (i-1)*nline+1, min(i*nline,nbnd)
              write(numstr,'(1x,1es15.8,a1)')
     .         backscattering(ibnd,is),','
              line = trim(line)//numstr
            enddo
            write(nwrite,'(a)') trim(line)
          enddo
          do i = 1, (nbnd-1)/nline+1
            if (i.eq.1) then
              write(line,'(a24,i3,a2)')
     .         'BACKSCATTERING_ENERGY(1,',is,')='
            else
              line = ''
            endif
            do ibnd = (i-1)*nline+1, min(i*nline,nbnd)
              write(numstr,'(1x,1es15.8,a1)')
     .         backscattering_energy(ibnd,is),','
              line = trim(line)//numstr
            enddo
            write(nwrite,'(a)') trim(line)
          enddo
        endif
      enddo
      if (plate_model.ge.1) then
cxpb  Use explicit writing of non-default elements only
        do ix = -1, nx
          do iy = -1, ny
            if (xymap(ix,iy).ne.0)
     .       write(nwrite,'(a6,i4,a1,i4,a2,i4,a1)')
     .       'XYMAP(',ix,',',iy,')=',xymap(ix,iy),','
          end do
        end do
        do ibnd = 1, nbnd
          write(nwrite,'(a6,i4,a2,i4,a1,1x,a6,i4,a2,i4,a1)')
     .     'IMAPX(',ibnd,')=',imapx(ibnd),',',
     .     'IMAPY(',ibnd,')=',imapy(ibnd),','
          if (inertial_cooling(ibnd))
     .     write(nwrite,'(a17,i4,a2,l1,a1)')
     .     'INERTIAL_COOLING(',ibnd,')=',inertial_cooling(ibnd),','
          if (backplate_temp(ibnd).ne.plate_temp .and.
     .        backplate_temp(ibnd).ne.0.0_R8)
     .     write(nwrite,'(a15,i4,a2,1es15.8,a1)')
     .     'BACKPLATE_TEMP(',ibnd,')=',backplate_temp(ibnd),','
          if (plate_thickness(ibnd).ne.plate_thick .and.
     .        plate_thickness(ibnd).ne.0.0_R8)
     .     write(nwrite,'(a16,i4,a2,1es15.8,a1)')
     .     'PLATE_THICKNESS(',ibnd,')=',plate_thickness(ibnd),','
          if (coating_thickness(ibnd).gt.0.0_R8)
     .     write(nwrite,'(a18,i4,a2,1es15.8,a1)')
     .     'COATING_THICKNESS(',ibnd,')=',coating_thickness(ibnd),','
          if (plate_time_factor(ibnd).ne.1.0_R8 .and.
     .        plate_time_factor(ibnd).ne.0.0_R8)
     .     write(nwrite,'(a18,i4,a2,1es15.8,a1)')
     .     'PLATE_TIME_FACTOR(',ibnd,')=',plate_time_factor(ibnd),','
          if (trim(surface_material_name(ibnd)).ne.'C' .and.
     .             surface_material_name(ibnd).ne.' ')
     .     write(nwrite,'(a22,i4,a4,a6,a2)')
     .     'SURFACE_MATERIAL_NAME(',ibnd,')= ''',
     .      surface_material_name(ibnd),''','
          if (coating_thickness(ibnd).gt.0.0_R8)
     .     write(nwrite,'(a22,i4,a4,a6,a2)')
     .     'COATING_MATERIAL_NAME(',ibnd,')= ''',
     .      coating_material_name(ibnd),''','
          if (trim(bulk_material_name(ibnd)).ne.'C' .and.
     .             bulk_material_name(ibnd).ne.' ')
     .     write(nwrite,'(a19,i4,a4,a6,a2)')
     .     'BULK_MATERIAL_NAME(',ibnd,')= ''',
     .      bulk_material_name(ibnd),''','
          todo =
     .     layer_nconstituents(ibnd).ne.matsurf_nconstituents(ibnd)
          do i = 1, layer_nconstituents(ibnd)
            found = .false.
            j = 1
            do while (.not.found .and. j.le.matsurf_nconstituents(ibnd))
              found = layer_nzconstituents(ibnd,i).eq.
     .              matsurf_nzconstituents(ibnd,j)
              if (.not.found) j = j + 1
            enddo
            if (found) todo = todo.or.
     .          layer_nrelconstituents(ibnd,i).ne.
     .        matsurf_nrelconstituents(ibnd,j)
          end do
          if (todo) then
            write(nwrite,'(a20,i4,a2,i2,a1)')
     .       'LAYER_NCONSTITUENTS(',ibnd,')=',
     .        layer_nconstituents(ibnd),','
            do i = 1, layer_nconstituents(ibnd)
              write(nwrite,'(a21,i4,a1,i2,a4,i2,a1)')
     .         'LAYER_NZCONSTITUENTS(',ibnd,',',i,')=  ',
     .          layer_nzconstituents(ibnd,i),','
              write(nwrite,'(a23,i4,a1,i2,a2,1es15.8,a1)')
     .         'LAYER_NRELCONSTITUENTS(',ibnd,',',i,')=',
     .          layer_nrelconstituents(ibnd,i),','
            end do
          end if
          do itrack=1,ntargsp
            if (deposition(ibnd,itrack).ne.0.0_R8)
     .       write(nwrite,'(a11,i4,a1,i1,a2,1es15.8,a1)')
     .        'DEPOSITION(',ibnd,',',itrack,')=',
     1         deposition(ibnd,itrack),','
            if (monolayer_deposition(ibnd,itrack).ne.0.0_R8)
     .       write(nwrite,'(a,i4,a,i1,a,1es15.8,a)')
     .        'MONOLAYER_DEPOSITION(',ibnd,',',itrack,')=',
     1         monolayer_deposition(ibnd,itrack),','
            if (erosion(ibnd,itrack).ne.0.0_R8)
     .       write(nwrite,'(a8,i4,a1,i1,a2,1es15.8,a1)')
     .        'EROSION(',ibnd,',',itrack,')=',
     1         erosion(ibnd,itrack),','
            if (monolayer_erosion(ibnd,itrack).ne.0.0_R8)
     .       write(nwrite,'(a,i4,a,i1,a,1es15.8,a)')
     .        'MONOLAYER_EROSION(',ibnd,',',itrack,')=',
     1         monolayer_erosion(ibnd,itrack),','
          enddo
          if (plate_time(ibnd).ne.0.0_R8)
     .     write(nwrite,'(a11,i4,a2,1es15.8,a1)')
     .      'PLATE_TIME(',ibnd,')=',plate_time(ibnd),','
          if (plate_area(ibnd).ne.0.0_R8)
     .     write(nwrite,'(a11,i4,a2,1es15.8,a1)')
     .      'PLATE_AREA(',ibnd,')=',plate_area(ibnd),','
        enddo
        if (plate_model.ge.1.and.ndepth_nml.gt.1) then
          do ibnd = 1, nbnd
            do i = 1, (ndepth_nml-1)/nline+1
              if (i.eq.1) then
                write(line,'(a12,i4,a3,i4,a2)')
     .           'TARGET_TEMP(',ibnd,',1:',ndepth_nml,')='
              else
                line = ''
              endif
              do is = (i-1)*nline+1, min(i*nline,ndepth_nml)
                write(numstr,'(1x,1es15.8,a1)')
     .           target_temp(ibnd,is),','
                line = trim(line)//numstr
              enddo
              write(nwrite,'(a)') trim(line)
            enddo
          enddo
        endif
      endif
      write(nwrite,'(a1)') '/'
      close(nwrite)
      write(*,*) 'Saved ',trim(filename)
      return
      end subroutine write_b2mod_wall_namelist

      subroutine read_b2mod_wall(ntargsp)
      implicit none
      integer ntargsp
      integer nmax, i, ic, ibnd, ifile, itrack
      integer old_file, idone(100), idum(0:9)
      real (kind=R8) :: sum
      logical open_file, silent_open, done
      external open_file, silent_open, xertst
      character*2  adum
      character*24 filename, name(100)
      character*32 str_dummy
      character*256 filenamecat

      if(open_file(98,'incatm','data.local|modules/B2.5/Database')) then
        do ic=1,94
          read(98,*) adum, idum(0), atomic_mass(idum(0)),
     1     surface_binding_energy(idum(0)), target_density(idum(0))
        enddo
        close(98)
      endif

      nmax=1
      ifile=0
      write(*,*) 'Obtaining nominal surface properties'
      do ibnd = 1, nbnd
        filename='Surface_properties/'
     &   //trim(surface_material_name(ibnd))
        name(ifile+1) = filename
        idone(ifile+1) = ibnd
        done = .false.
        do i = 1, ifile
          if (trim(filename).eq.trim(name(i))) then
            done = .true.
            old_file = idone(i)
          end if
        end do
        if (.not.done) then
          ifile = ifile + 1
          if(open_file(99,filename,
     .       'data.local|modules/B2.5/Database')) then
            read(99,*) surface_long_name(ibnd)
            read(99,*) matsurf_nconstituents(ibnd)
            nmax = max(nmax,matsurf_nconstituents(ibnd))
            call xertst (nmax.le.6, 'Increase NMAX=6 in b2mod_wall!')
            read(99,*) (idum(i),i=1,matsurf_nconstituents(ibnd))   ! matsurf_nzconstituents(:,:) not yet allocated
            read(99,*) (idum(i),i=1,matsurf_nconstituents(ibnd))   ! matsurf_nrelconstituents(:,:) not yet allocated
            read(99,*) matsurf_chemsput_erel(ibnd)
            read(99,*) matsurf_blackbody_emiss(ibnd)
            read(99,*) matsurf_vapprs_A(ibnd)
            read(99,*) matsurf_vapprs_B(ibnd)
            read(99,*) matsurf_stick_coef(ibnd)
            read(99,*) matsurf_vapor_temp(ibnd)
            read(99,*) matsurf_heat_sublim(ibnd)
            close(99)
          endif
        else
          surface_long_name(ibnd) = surface_long_name(old_file)
          matsurf_nconstituents(ibnd) = matsurf_nconstituents(old_file)
          matsurf_chemsput_erel(ibnd) = matsurf_chemsput_erel(old_file)
          matsurf_blackbody_emiss(ibnd) =
     .     matsurf_blackbody_emiss(old_file)
          matsurf_vapprs_A(ibnd) = matsurf_vapprs_A(old_file)
          matsurf_vapprs_B(ibnd) = matsurf_vapprs_B(old_file)
          matsurf_stick_coef(ibnd) = matsurf_stick_coef(old_file)
          matsurf_vapor_temp(ibnd) = matsurf_vapor_temp(old_file)
          matsurf_heat_sublim(ibnd) = matsurf_heat_sublim(old_file)
        endif
      enddo
      allocate(matsurf_nzconstituents(nwall,nmax))
      allocate(matsurf_nrelconstituents(nwall,nmax))

      ifile=0
      do ibnd = 1, nbnd
        filename='Surface_properties/'
     &   //trim(surface_material_name(ibnd))
        done = .false.
        do i = 1, ifile
          if (trim(filename).eq.trim(name(i))) then
            done = .true.
            old_file = idone(i)
          end if
        end do
        if (.not.done) then
          ifile = ifile + 1
          if(silent_open(99,filename,
     .       'data.local|modules/B2.5/Database',filenamecat)) then
            read(99,*) str_dummy
            read(99,*) idum(0)
            read(99,*) (matsurf_nzconstituents(ibnd,i),i=1,idum(0))
            read(99,*) (matsurf_nrelconstituents(ibnd,i),i=1,idum(0))
            close(99)
          endif
        else
          matsurf_nzconstituents(ibnd,
     .      1:matsurf_nconstituents(old_file)) =
     .     matsurf_nzconstituents(old_file,
     .      1:matsurf_nconstituents(old_file))
          matsurf_nrelconstituents(ibnd,
     .      1:matsurf_nconstituents(old_file)) =
     .     matsurf_nrelconstituents(old_file,
     .      1:matsurf_nconstituents(old_file))
        endif
      enddo
c-xpb Normalise surface constituents
      do ibnd = 1, nbnd
        sum = 0.0_R8
        do i = 1, matsurf_nconstituents(ibnd)
          sum = sum + matsurf_nrelconstituents(ibnd,i)
        end do
        do i = 1, matsurf_nconstituents(ibnd)
          matsurf_nrelconstituents(ibnd,i) =
     .     matsurf_nrelconstituents(ibnd,i)/sum
        end do
      end do

      ifile=0
      write(*,*) 'Obtaining nominal bulk properties'
      do ibnd = 1, nbnd
        matbulk_Kappa_perp_A1(ibnd) = 0.0_R8
        matbulk_Kappa_perp_A2(ibnd) = 0.0_R8
        matbulk_Kappa_perp_A3(ibnd) = 0.0_R8
        matbulk_Kappa_perp_A4(ibnd) = 0.0_R8
        filename='Bulk_properties/'//trim(bulk_material_name(ibnd))
        name(ifile+1) = filename
        idone(ifile+1) = ibnd
        done = .false.
        do i = 1, ifile
          if (trim(filename).eq.trim(name(i))) then
            done = .true.
            old_file = idone(i)
          end if
        end do
        if (.not.done) then
          ifile = ifile + 1
          if(open_file(99,filename,
     .       'data.local|modules/B2.5/Database')) then
            read(99,*) bulk_long_name(ibnd)
            read(99,*) matbulk_target_density(ibnd)
            read(99,*) matbulk_Kappa_KA(ibnd)
            read(99,*) matbulk_Kappa_KB(ibnd)
            read(99,*) matbulk_Kappa_A1(ibnd)
            read(99,*) matbulk_Kappa_A2(ibnd)
            read(99,*) matbulk_Kappa_A3(ibnd)
            read(99,*) matbulk_Kappa_A4(ibnd)
            read(99,*) matbulk_Cp_B1(ibnd)
            read(99,*) matbulk_Cp_B2(ibnd)
            read(99,*) matbulk_Cp_B3(ibnd)
            read(99,*,err=98,end=98) matbulk_Kappa_perp_A1(ibnd)
            read(99,*) matbulk_Kappa_perp_A2(ibnd)
            read(99,*) matbulk_Kappa_perp_A3(ibnd)
            read(99,*) matbulk_Kappa_perp_A4(ibnd)
98          close(99)
            if (matbulk_Kappa_perp_A1(ibnd).eq.0.0_R8) then
              matbulk_Kappa_perp_A1(ibnd) = matbulk_Kappa_A1(ibnd)
              matbulk_Kappa_perp_A2(ibnd) = matbulk_Kappa_A2(ibnd)
              matbulk_Kappa_perp_A3(ibnd) = matbulk_Kappa_A3(ibnd)
              matbulk_Kappa_perp_A4(ibnd) = matbulk_Kappa_A4(ibnd)
            endif
          endif
        else
          bulk_long_name(ibnd) = bulk_long_name(old_file)
          matbulk_target_density(ibnd) = 
     .     matbulk_target_density(old_file) 
          matbulk_Kappa_KA(ibnd) = matbulk_Kappa_KA(old_file)
          matbulk_Kappa_KB(ibnd) = matbulk_Kappa_KB(old_file)
          matbulk_Kappa_A1(ibnd) = matbulk_Kappa_A1(old_file)
          matbulk_Kappa_A2(ibnd) = matbulk_Kappa_A2(old_file)
          matbulk_Kappa_A3(ibnd) = matbulk_Kappa_A3(old_file)
          matbulk_Kappa_A4(ibnd) = matbulk_Kappa_A4(old_file)
          matbulk_Kappa_perp_A1(ibnd) = matbulk_Kappa_perp_A1(old_file)
          matbulk_Kappa_perp_A2(ibnd) = matbulk_Kappa_perp_A2(old_file)
          matbulk_Kappa_perp_A3(ibnd) = matbulk_Kappa_perp_A3(old_file)
          matbulk_Kappa_perp_A4(ibnd) = matbulk_Kappa_perp_A4(old_file)
          matbulk_Cp_B1(ibnd) = matbulk_Cp_B1(old_file)
          matbulk_Cp_B2(ibnd) = matbulk_Cp_B2(old_file)
          matbulk_Cp_B3(ibnd) = matbulk_Cp_B3(old_file)
        endif
      enddo

      ifile=0
      do ibnd = 1, nbnd
        if (coating_thickness(ibnd).gt.0.0_R8) then
          if (ifile.eq.0) write(*,*) 'Obtaining coating properties'
          coating_Kappa_perp_A1(ibnd) = 0.0_R8
          coating_Kappa_perp_A2(ibnd) = 0.0_R8
          coating_Kappa_perp_A3(ibnd) = 0.0_R8
          coating_Kappa_perp_A4(ibnd) = 0.0_R8
          filename='Bulk_properties/'//trim(coating_material_name(ibnd))
          name(ifile+1) = filename
          idone(ifile+1) = ibnd
          done = .false.
          do i = 1, ifile
            if (trim(filename).eq.trim(name(i))) then
              done = .true.
              old_file = idone(i)
            end if
          end do
          if (.not.done) then
            ifile = ifile + 1
            if(open_file(99,filename,
     .         'data.local|modules/B2.5/Database')) then
              read(99,*) coating_long_name(ibnd)
              read(99,*) coating_target_density(ibnd)
              read(99,*) coating_Kappa_KA(ibnd)
              read(99,*) coating_Kappa_KB(ibnd)
              read(99,*) coating_Kappa_A1(ibnd)
              read(99,*) coating_Kappa_A2(ibnd)
              read(99,*) coating_Kappa_A3(ibnd)
              read(99,*) coating_Kappa_A4(ibnd)
              read(99,*) coating_Cp_B1(ibnd)
              read(99,*) coating_Cp_B2(ibnd)
              read(99,*) coating_Cp_B3(ibnd)
              read(99,*,err=97,end=97) coating_Kappa_perp_A1(ibnd)
              read(99,*) coating_Kappa_perp_A2(ibnd)
              read(99,*) coating_Kappa_perp_A3(ibnd)
              read(99,*) coating_Kappa_perp_A4(ibnd)
97            close(99)
              if (coating_Kappa_perp_A1(ibnd).eq.0.0_R8) then
                coating_Kappa_perp_A1(ibnd) = coating_Kappa_A1(ibnd)
                coating_Kappa_perp_A2(ibnd) = coating_Kappa_A2(ibnd)
                coating_Kappa_perp_A3(ibnd) = coating_Kappa_A3(ibnd)
                coating_Kappa_perp_A4(ibnd) = coating_Kappa_A4(ibnd)
              endif
            endif
          else
            coating_long_name(ibnd) = coating_long_name(old_file)
            coating_target_density(ibnd) = 
     .       coating_target_density(old_file) 
            coating_Kappa_KA(ibnd) = coating_Kappa_KA(old_file)
            coating_Kappa_KB(ibnd) = coating_Kappa_KB(old_file)
            coating_Kappa_A1(ibnd) = coating_Kappa_A1(old_file)
            coating_Kappa_A2(ibnd) = coating_Kappa_A2(old_file)
            coating_Kappa_A3(ibnd) = coating_Kappa_A3(old_file)
            coating_Kappa_A4(ibnd) = coating_Kappa_A4(old_file)
            coating_Kappa_perp_A1(ibnd) = 
     .       coating_Kappa_perp_A1(old_file)
            coating_Kappa_perp_A2(ibnd) = 
     .       coating_Kappa_perp_A2(old_file)
            coating_Kappa_perp_A3(ibnd) = 
     .       coating_Kappa_perp_A3(old_file)
            coating_Kappa_perp_A4(ibnd) = 
     .       coating_Kappa_perp_A4(old_file)
            coating_Cp_B1(ibnd) = coating_Cp_B1(old_file)
            coating_Cp_B2(ibnd) = coating_Cp_B2(old_file)
            coating_Cp_B3(ibnd) = coating_Cp_B3(old_file)
          endif
        endif
      enddo

      do ifile = 1, nalloys
        if(layer_alloys(ifile).ne.' ') then
          if (ifile.eq.0) write(*,*) 'Obtaining alloy properties'
          filename='Bulk_properties/'//trim(layer_alloys(ifile))
          if(open_file(99,filename,
     .       'data.local|modules/B2.5/Database')) then
            read(99,*) alloy_long_name(ifile)
            read(99,*) alloy_Kappa_KA(ifile)
            read(99,*) alloy_Kappa_KB(ifile)
            read(99,*) alloy_Kappa_A1(ifile)
            read(99,*) alloy_Kappa_A2(ifile)
            read(99,*) alloy_Kappa_A3(ifile)
            read(99,*) alloy_Kappa_A4(ifile)
            read(99,*) alloy_Cp_B1(ifile)
            read(99,*) alloy_Cp_B2(ifile)
            read(99,*) alloy_Cp_B3(ifile)
96          close(99)
          endif
          filename='Surface_properties/'//trim(layer_alloys(ifile))
          if(open_file(99,filename,
     .       'data.local|imodules/B2.5/Database')) then
            read(99,*) adum
            read(99,*) alloy_nconstituents(ifile)
            read(99,*) alloy_nzconstituents
     .       (ifile,1:alloy_nconstituents(ifile))
            read(99,*) alloy_nrelconstituents
     .       (ifile,1:alloy_nconstituents(ifile))
            read(99,*) alloy_chemsput_erel(ifile)
            read(99,*) alloy_blackbody_emiss(ifile)
            read(99,*) alloy_vapprs_A(ifile)
            read(99,*) alloy_vapprs_B(ifile)
            read(99,*) alloy_stick_coef(ifile)
            read(99,*) alloy_vapor_temp(ifile)
            read(99,*) alloy_heat_sublim(ifile)
            close(99)
          endif
        endif
      enddo

      do itrack = 1, ntargsp
        if (track_species(itrack).ne.'  ') then
          if (itrack.eq.1) 
     >     write(*,*) 'Obtaining tracked species properties'
          filename='Bulk_properties/'//trim(track_species(itrack))
          if(open_file(99,filename,
     .       'data.local|modules/B2.5/Database')) then
            read(99,*) species_long_name(itrack)
            read(99,*) species_target_density(itrack)
            read(99,*) species_Kappa_KA(itrack)
            read(99,*) species_Kappa_KB(itrack)
            read(99,*) species_Kappa_A1(itrack)
            read(99,*) species_Kappa_A2(itrack)
            read(99,*) species_Kappa_A3(itrack)
            read(99,*) species_Kappa_A4(itrack)
            read(99,*) species_Cp_B1(itrack)
            read(99,*) species_Cp_B2(itrack)
            read(99,*) species_Cp_B3(itrack)
            close(99)
          endif
          filename='Surface_properties/'//trim(track_species(itrack))
          if(open_file(99,filename,
     .       'data.local|modules/B2.5/Database')) then
            read(99,*) adum
            read(99,*) idum(0)
            read(99,*) idum(1:idum(0))
            read(99,*) idum(1:idum(0))
            read(99,*) adum
            read(99,*) species_blackbody_emiss(itrack)
            read(99,*) species_vapprs_A(itrack)
            read(99,*) species_vapprs_B(itrack)
            read(99,*) species_stick_coef(itrack)
            read(99,*) species_vapor_temp(itrack)
            read(99,*) species_heat_sublim(itrack)
            close(99)
          endif
        endif
      enddo

      return
      end subroutine read_b2mod_wall

#ifndef NO_CDF
      subroutine define_cdf_b2mod_wall(nnwall,nntrack)
      use b2mod_tallies
      implicit none
      integer nnwall,nntrack
#include <netcdf.inc>
      integer :: status, ncid, nwallid, ntrackid, ntwoid, tid
      integer :: vdims(3)

      write(*,*) 'Creating the definitions for b2wall.nc'

      status=nf_create('b2wall.nc',NF_CLOBBER,ncid)
      call check_cdf_status(status)

      write(*,*) 'element'
      status=nf_def_dim(ncid,'element',2,ntwoid)
      call check_cdf_status(status)
      write(*,*) 'nwall'
      status=nf_def_dim(ncid,'nwall',nnwall,nwallid)
      call check_cdf_status(status)
      write(*,*) 'ntrack'
      status=nf_def_dim(ncid,'ntrack',nntrack,ntrackid)
      call check_cdf_status(status)
      write(*,*) 'time'
      status=nf_def_dim(ncid,'time',NF_UNLIMITED,tid)
      call check_cdf_status(status)

      vdims(1)=ntwoid
      vdims(2)=ntrackid
      write(*,*) 'track_species'
      call create_cdf_var_with_attrs(
     *  ncid,'track_species',NF_CHAR,
     *  2,vdims,'','Tracked species')

      vdims(1)=nwallid
      write(*,*) 'plate_area'
      call create_cdf_var_with_attrs(
     *  ncid,'plate_area',NF_DOUBLE,
     *  1,vdims,'m^2','Plate area')

      vdims(1)=nwallid
      vdims(2)=tid
      write(*,*) 'plate_time'
      call create_cdf_var_with_attrs(
     *  ncid,'plate_time',NF_DOUBLE,
     *  2,vdims,'s','Time elapsed')

      vdims(1)=nwallid
      vdims(2)=ntrackid
      vdims(3)=tid
      write(*,*) 'erosion'
      call create_cdf_var_with_attrs(
     *  ncid,'erosion',NF_DOUBLE,
     *  3,vdims,'particles','Eroded particles (gross)')
      write(*,*) 'deposition'
      call create_cdf_var_with_attrs(
     *  ncid,'deposition',NF_DOUBLE,
     *  3,vdims,'particles','Deposited particles (gross)')
      write(*,*) 'monolayer_erosion'
      call create_cdf_var_with_attrs(
     *  ncid,'monolayer_erosion',NF_DOUBLE,
     *  3,vdims,'monolayers','Eroded particles (gross)')
      write(*,*) 'monolayer_deposition'
      call create_cdf_var_with_attrs(
     *  ncid,'monolayer_deposition',NF_DOUBLE,
     *  3,vdims,'monolayers','Deposited particles (gross)')

      status=nf_close(ncid)
      call check_cdf_status(status)

      write(*,*) 'Created the definitions for b2wall.nc'

      return
      end subroutine define_cdf_b2mod_wall
#endif

      end module b2mod_wall
