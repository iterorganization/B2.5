      module b2mod_wall
      use b2mod_types
      use b2mod_layer
#     include "DIMENSIONS.F"
      
      integer ndepth, nwall_segments, nsegments, nline
      parameter(
     *  nwall_segments=5*DEF_NCUT+1+1,
     *  ndepth=201,
     *  nline=15)
cxpb  Namelist data
      integer, save :: nbnd, ndepth_nml,
     *  imapx(nwall), imapy(nwall), imapw(nwall), imapreg(nwall),
     *  wall_begin(nwall_segments), wall_end(nwall_segments)
      logical inertial_cooling(nwall)
      real (kind=R8), save ::
     *  target_temp(nwall,ndepth),
     *  target_depth(nwall,ndepth),
     *  target_width(nwall,ndepth),
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

cxpb Data file material properties
      character*6, save ::
     *  surface_material_name(nwall)
      character*32, save :: surface_long_name(nwall)
      integer, save ::
     *  matsurf_nconstituents(nwall)
      integer, allocatable, save ::
     *  matsurf_nzconstituents(:,:),
     *  matsurf_nrelconstituents(:,:)
      real (kind=R8), save ::
     *  matsurf_chemsput_erel(nwall),
     *  matsurf_blackbody_emiss(nwall),
     *  matsurf_vapprs_A(nwall),
     *  matsurf_vapprs_B(nwall),
     *  matsurf_stick_coef(nwall),
     *  matsurf_vapor_temp(nwall),
     *  matsurf_heat_sublim(nwall)
      
      character*6, save :: coating_material_name(nwall)
      character*32, save :: coating_long_name(nwall)
      real (kind=R8), save ::
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
     *  coating_Cp_B3(nwall)

      character*6, save ::
     *  bulk_material_name(nwall)
      character*32, save :: bulk_long_name(nwall)
      real (kind=R8), save ::
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
     *  matbulk_Cp_B3(nwall)

      real(kind=R8), save ::         ! data from incatm file
     *  atomic_mass(92),             ! in amu
     1  surface_binding_energy(92),  ! in eV
     2  target_density(92)           ! in g/cm3

      integer, save :: itrack_0(nwall)
      integer, save :: track_index(0:DEF_NSD)
      character*2, save :: track_species(ntrack)
      logical, save :: track_chem_sput(ntrack)

      real (kind=R8) :: r1(nwall), r2(nwall), z1(nwall), z2(nwall)
#ifdef ALLOCATE_AND_NAMELIST
      integer, allocatable, save ::
     *  xymap(:,:)
#else
      integer xymap(-1:DEF_NXD,-1:DEF_NYD)
#endif

cxpb Mapping arrays:
cxpb A wall element numbered iwall can be found in the B2 grid at the position
cxpb   ( imapx(iwall) , imapy(iwall) )
cxpb Reversely, the B2 boundary cell ( ix, iy ) corresponds to the boundary
cxpb element: xymap(ix,iy)

      contains

      subroutine alloc_b2mod_wall(nx,ny)
      implicit none
      integer nx, ny
      logical errdim
      external subini, subend, xertst

      call subini('b2mod_wall')
#ifdef ALLOCATE_AND_NAMELIST
      allocate(xymap(-1:nx,-1:ny))
#else
      errdim=.false.
      if(nx.gt.DEF_NXD) then
        write(*,*) 'alloc_b2mod_wall: nxd > def_nxd',nx,DEF_NXD
        errdim=.true.
      endif
      if(ny.gt.DEF_NYD) then
        write(*,*) 'alloc_b2mod_wall: nyd > def_nyd',ny,DEF_NYD
        errdim=.true.
      endif
      call xertst (.not.errdim,
     1 'faulty def_nxd, def_nyd')
#endif
      call xertst(1.le.ndepth .and. ndepth.le.1001,
     . 'invalid dimension ndepth')
      call subend()
      return

      end subroutine alloc_b2mod_wall

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

      if(open_file(98,'incatm','data.local|data')) then
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
          if(open_file(99,filename,'data.local|data')) then
            read(99,*) surface_long_name(ibnd)
            read(99,*) matsurf_nconstituents(ibnd)
            nmax = max(nmax,matsurf_nconstituents(ibnd))
            call xertst (nmax.le.6, 'Increase NMAX=6 in b2mod_wall!')
            read(99,*) idum(1:matsurf_nconstituents(ibnd))   ! matsurf_nzconstituents(:,:) not yet allocated
            read(99,*) idum(1:matsurf_nconstituents(ibnd))   ! matsurf_nrelconstituents(:,:) not yet allocated
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
          if(silent_open(99,filename,'data.local|data',filenamecat)) 
     >     then
            read(99,*) str_dummy
            read(99,*) idum(0)
            read(99,*) matsurf_nzconstituents(ibnd,1:idum(0)) 
            read(99,*) matsurf_nrelconstituents(ibnd,1:idum(0)) 
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
          if(open_file(99,filename,'data.local|data')) then
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
            if(open_file(99,filename,'data.local|data')) then
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
          if(open_file(99,filename,'data.local|data')) then
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
          if(open_file(99,filename,'data.local|data')) then
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
          if(open_file(99,filename,'data.local|data')) then
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
          if(open_file(99,filename,'data.local|data')) then
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
