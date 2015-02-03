      subroutine init_wall(nx,ny,ns,BoRiS,
     .   plate_model, plate_option, plate_temp, plate_thick)
      use b2mod_types 
      use b2mod_geo
      use b2mod_rates
      use b2mod_plasma 
      use b2mod_sputter
      use b2mod_indirect
      use b2mod_wall
      use b2mod_layer
      use b2mod_b2cmpa
      use b2mod_b2cmpb
#     include "DIMENSIONS.F"
      implicit none

      integer nx, ny, ns
      integer plate_model, plate_option
      real (kind=R8) :: BoRiS, plate_temp, plate_thick

      integer ix, iy, is, ireg
      integer ibnd, ibndp, ibndm, i, j, itrack
      integer imapx_tmp(nwall), imapy_tmp(nwall)
      integer xymap_tmp(-1:nx,-1:ny)
      real (kind=R8) :: qin, qel, normal_r, normal_z, 
     *  n1_r, n1_z, n2_r, n2_z, dist, sina, cosa, 
     *  coshlf, sinhlf, tanhlf, nhalf_r, nhalf_z, sum,
     *  sput_frac(0:ntrack)
      logical wall_namelist_used, default_cooling
      save wall_namelist_used
      data wall_namelist_used /.false./
      character*400 line
      character*17 numstr
      character*2 track_species_save(ntrack)
      logical todo, found

      integer, save :: hd_linlog, output, sput_dst, sput_frac_flag
      real (kind=R8), save :: hd_ratio
      data hd_linlog/1/, hd_ratio/100.0_R8/, output/0/
      data sput_dst/-1/, sput_frac_flag/0/

      real (kind=R8) :: mass_density
      intrinsic mod, max, min, sqrt
      external mass_density
      external xertst, subini, subend, ipgeti, ipgetr

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

      call subini('init_wall')
      call alloc_b2mod_wall(nx,ny)
cxpb Establish numbering of boundary wall cells
cxpb We must avoid corner cells and dead cells
      nbnd = 0
      ndepth_nml = 1
      imapx_tmp= 0
      imapy_tmp= 0
      xymap_tmp= 0
      wall_begin = 0
      wall_end = 0
      nsegments = 0
      inertial_cooling = .false.
      backplate_temp = 0.0_R8
      plate_thickness = 0.0_R8
      coating_thickness = 0.0_R8
      plate_time_factor = 0.0_R8
      surface_material_name = ' '
      coating_material_name = ' '
      bulk_material_name = ' '
      itrack_0 = 0
      layer_alloys = ' '
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
      layer_nconstituents = 0
      layer_nzconstituents = 0
      layer_nrelconstituents = 0.0_R8
      layer_composition = 0.0_R8
      layer_thermal_contact = 1.0_R8
      
cxpb  We follow the same order as in b2plot
cxpb  Northern boundaries (excluding diverted island)
      ireg = cbirno
      do iy = ny, -1, -1
        do ix = -1, nx
          if (topiy(ix,iy).eq.ny+1 .and. region(ix,iy,0).ne.0 .and.
     &        region(ix,iy,0).ne.5 .and.
     &        leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            nbnd = nbnd + 1
            imapx_tmp(nbnd) = ix
            imapy_tmp(nbnd) = iy
            r1(nbnd)=crx(ix,iy,2)
            z1(nbnd)=cry(ix,iy,2)
            r2(nbnd)=crx(ix,iy,3)
            z2(nbnd)=cry(ix,iy,3)
            xymap_tmp(ix,iy)= nbnd
            plate_area(nbnd) = gs(ix,iy,1)
            if (xymap_tmp(ix,iy).ne.0 .and. 
     &          (xymap_tmp(leftix(ix,iy),leftiy(ix,iy)).eq.0 .or.
     &          (leftix(ix,iy).eq.ix .and. leftiy(ix,iy).eq.iy))) then
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd
            endif 
          endif
          if(iy.eq.ny) then
            if (cbrbrk(ireg)*nx.lt.ix+0.5_R8) then
              ireg = ireg + 1
            endif
            if (xymap_tmp(ix,iy).ne.0) imapreg(nbnd) = ireg
          endif
        enddo
      enddo
cxpb  Eastern boundaries
      ireg = cbirea+cbnrea-1
      do ix = nx, -1, -1
        do iy = ny, -1, -1
          if (rightix(ix,iy).eq.nx+1 .and. region(ix,iy,0).ne.0 .and.
     &        bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
            nbnd = nbnd + 1
            imapx_tmp(nbnd) = ix
            imapy_tmp(nbnd) = iy
            r1(nbnd)=crx(ix,iy,3)
            z1(nbnd)=cry(ix,iy,3)
            r2(nbnd)=crx(ix,iy,1)
            z2(nbnd)=cry(ix,iy,1)
            xymap_tmp(ix,iy)= nbnd
            plate_area(nbnd) = gs(ix,iy,0)
            if (xymap_tmp(ix,iy).ne.0 .and. 
     &          (xymap_tmp(topix(ix,iy),topiy(ix,iy)).eq.0 .or.
     &          (topix(ix,iy).eq.ix .and. topiy(ix,iy).eq.iy))) then
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd
            endif 
          endif
          if(ix.eq.nx) then
            if (ireg.gt.cbirea) then
              if (cbrbrk(ireg-1)*ny.gt.iy+0.5_R8) then
                ireg = ireg - 1
              endif
            endif
            if (xymap_tmp(ix,iy).ne.0) imapreg(nbnd) = ireg
          endif
        enddo
      enddo
cxpb  Southern boundaries (excluding core)
      ireg = cbirso+cbnrso-1
      do iy = -1, ny
        do ix = nx, -1, -1      
          if (bottomiy(ix,iy).eq.-2 .and. region(ix,iy,0).ne.0 .and.
     &        mod(region(ix,iy,0),4).ne.1 .and.
     &        leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            nbnd = nbnd + 1
            imapx_tmp(nbnd) = ix
            imapy_tmp(nbnd) = iy
            r1(nbnd)=crx(ix,iy,1)
            z1(nbnd)=cry(ix,iy,1)
            r2(nbnd)=crx(ix,iy,0)
            z2(nbnd)=cry(ix,iy,0)
            xymap_tmp(ix,iy)= nbnd
            plate_area(nbnd) = gs(topix(ix,iy),topiy(ix,iy),1)
            if (xymap_tmp(ix,iy).ne.0 .and. 
     &          (xymap_tmp(rightix(ix,iy),rightiy(ix,iy)).eq.0 .or.
     &          (rightix(ix,iy).eq.ix .and. rightiy(ix,iy).eq.iy))) then
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd
            endif 
          endif
          if(iy.eq.-1) then
            if (ireg.gt.cbirso) then
              if (cbrbrk(ireg-1)*nx.gt.ix+0.5_R8) then
                ireg = ireg - 1
              endif
            endif
            if (xymap_tmp(ix,iy).ne.0) imapreg(nbnd) = ireg
          endif
        enddo
      enddo
cxpb  Western boundaries
      ireg = cbirwe
      do ix = -1, nx
        do iy = -1, ny
          if (leftix(ix,iy).eq.-2 .and. region(ix,iy,0).ne.0 .and.
     &        bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
            nbnd = nbnd + 1
            imapx_tmp(nbnd) = ix
            imapy_tmp(nbnd) = iy
            r1(nbnd)=crx(ix,iy,0)
            z1(nbnd)=cry(ix,iy,0)
            r2(nbnd)=crx(ix,iy,2)
            z2(nbnd)=cry(ix,iy,2)
            xymap_tmp(ix,iy)= nbnd
            plate_area(nbnd) = gs(rightix(ix,iy),rightiy(ix,iy),0)
            if (xymap_tmp(ix,iy).ne.0 .and. 
     &          (xymap_tmp(bottomix(ix,iy),bottomiy(ix,iy)).eq.0 .or.
     &          (rightix(ix,iy).eq.ix .and. rightiy(ix,iy).eq.iy))) then
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd
            endif 
          endif
          if(ix.eq.-1) then
            if (cbrbrk(ireg)*ny.lt.iy+0.5_R8) then
              ireg = ireg + 1
            endif
            if (xymap_tmp(ix,iy).ne.0) imapreg(nbnd) = ireg
          endif
        enddo
      enddo
cxpb  Compute end segment positions
      do is = 2, nsegments
        wall_end(is-1) = wall_begin(is)-1
      enddo
      wall_end(nsegments) = nbnd
      do is = nsegments+1, nwall_segments
        wall_begin(is) = nbnd+1
      enddo
      call xertst (nbnd.le.2*DEF_NCUT*(DEF_NXD+2)+2*(DEF_NYD+2),
     & 'Increase size of B2 contribution to nwall in b2mod_wall')

      inquire(file='../b2.wall_save.parameters',
     1 exist=wall_namelist_used)
      if(wall_namelist_used) then
        track_species_save=track_species
        open(99,file='../b2.wall_save.parameters')
        write(*,*) 'Reading ../b2.wall_save.parameters'
        read(99,wall_save)
        close(99)
        write(*,*) 'Read ../b2.wall_save.parameters'
        call xertst(ndepth_nml.le.ndepth .and. ndepth_nml.ge.1,
     .   'Invalid array size NDEPTH_NML!')
        write(*,*) 'track_species : Bulk ', track_species(1:ntargsp)
        do is = 1, ntrack
          call xertst(track_species_save(is).eq.track_species(is),
     .     'inconsistent declaration of tracked species')
        enddo
        do ix = -1, nx
          do iy = -1, ny
            if (xymap(ix,iy).ne.0) then 
              call xertst(xymap(ix,iy).eq.xymap_tmp(ix,iy),
     .         'Inconsistent numbering between namelist and code!')
              call xertst(
     .           imapx(xymap(ix,iy)).eq.ix .and.
     .           imapy(xymap(ix,iy)).eq.iy,
     .       'Inconsistent boundary mapping in namelist WALL!')
            endif
          end do
        end do
      endif
      imapx(1:nbnd) = imapx_tmp(1:nbnd)
      imapy(1:nbnd) = imapy_tmp(1:nbnd)
      xymap(-1:nx,-1:ny) = xymap_tmp(-1:nx,-1:ny)
      call subend()
      return

      entry wall_data(nx,ny,ns,BoRiS,
     .   plate_model, plate_option, plate_temp, plate_thick)
      call subini('wall_data')
cxpb Fill up array data not given in namelist
      do ibnd = 1, nbnd
        if (backplate_temp(ibnd).eq.0.0_R8)
     .      backplate_temp(ibnd) = plate_temp
        if (plate_thickness(ibnd).eq.0.0_R8)
     .      plate_thickness(ibnd) = plate_thick
        if (plate_time_factor(ibnd).eq.0.0_R8)
     .      plate_time_factor(ibnd) = 1.0_R8
        if (len_trim(surface_material_name(ibnd)).eq.0)
     .      surface_material_name(ibnd)(1:1) = 'C'   ! Graphite
        if (len_trim(bulk_material_name(ibnd)).eq.0)
     .      bulk_material_name(ibnd)(1:1) = 'C'   ! Graphite
      enddo

cxpb Read material data
      call ipgeti('b2stbr_sput_dst', sput_dst)
      call xertst(sput_dst.lt.ns,
     .   'faulty input parameter sput_dst')
      if (sput_dst.lt.0) then
        write(*,*) 'No sputtering enabled : Target models turned off!'
        call subend()
        return
      endif
      call read_b2mod_wall(ntargsp)

cxpb For mixed materials, if unspecified, use normal material constituents
      do ibnd = 1, nbnd
        if (layer_nconstituents(ibnd).eq.0) then
          layer_nconstituents(ibnd)=
     1     matsurf_nconstituents(ibnd)
          layer_nzconstituents(ibnd,1:layer_nconstituents(ibnd))=
     1     matsurf_nzconstituents(ibnd,1:layer_nconstituents(ibnd))
          layer_nrelconstituents(ibnd,1:layer_nconstituents(ibnd))=
     1     matsurf_nrelconstituents(ibnd,1:layer_nconstituents(ibnd))
        end if
cxpb  normalise the layer constituents relative concentrations
        sum = 0.0_R8
        do i = 1, layer_nconstituents(ibnd)
          sum = sum + layer_nrelconstituents(ibnd,i)
        enddo
        layer_nrelconstituents(ibnd,1:layer_nconstituents(ibnd))=
     1   layer_nrelconstituents(ibnd,1:layer_nconstituents(ibnd))/sum
      end do

cxpb Compute target geometry data
      call ipgeti('b2stbr_output', output)
      call ipgeti('heatdiff1D_linlog', hd_linlog)
      call xertst(hd_linlog.eq.1 .or. hd_linlog.eq.2,
     .   'faulty parameter heatdiff1D_linlog')
      if (hd_linlog.eq.2) then
        call ipgetr('heatdiff1D_ratio', hd_ratio)
        call xertst(hd_ratio.gt.1.0_R8, 
     .   'faulty parameter heatdiff1D_ratio')
      endif
      call ipgeti ('b2stbr_sput_frac_flag', sput_frac_flag)
      do ibnd = 1, nbnd
c-xpb  Compute initial layer thickness (deposited - eroded)
c-xpb  and layer depth (sum of deposition), 
c-xpb  i.e. minimum depth over which material mixing occurs
        layer_thickness(ibnd)=0.0_R8
        layer_depth(ibnd)=0.0_R8
        do i=1,ndepth_nml-1
          if (hd_linlog.eq.1) then
            target_depth(ibnd,i) = plate_thickness(ibnd) * 
     .       float(i-1) / float(ndepth_nml - 1)
          else
            target_depth(ibnd,i) = plate_thickness(ibnd) * 
     .       (hd_ratio ** (float(i-1) / float(ndepth_nml - 1)) 
     .                                                - 1.0_R8)/
     .       (hd_ratio - 1.0_R8)
          endif
        enddo
        if (sput_frac_flag.ne.0) then  
c mixed materials with layers
          do itrack=1,ntargsp
            layer_thickness(ibnd)=layer_thickness(ibnd)+
     &       (monolayer_deposition(ibnd,itrack)-
     &        monolayer_erosion(ibnd,itrack))*monolayer_l(itrack)
            layer_depth(ibnd)=layer_depth(ibnd)+
     &        monolayer_deposition(ibnd,itrack)*monolayer_l(itrack)
          enddo
          call add_layer_thickness(nx,ny,ns,ibnd)
        else if (sput_frac_flag.eq.0) then
c no mixed materials
          layer_nconstituents(ibnd)=matsurf_nconstituents(ibnd)
          layer_nzconstituents(ibnd,1:matsurf_nconstituents(ibnd))=
     1      matsurf_nzconstituents
     .       (ibnd,1:matsurf_nconstituents(ibnd))
          layer_nrelconstituents(ibnd,1:matsurf_nconstituents(ibnd))=
     1      matsurf_nrelconstituents
     .       (ibnd,1:matsurf_nconstituents(ibnd))
          layer_mass_density(ibnd)=mass_density(ibnd,2)
          sput_frac(0)=1.0_R8
          sput_frac(1:ntrack)=0.0_R8
          call compose_layer(ibnd,sput_frac)
        endif
        if (ndepth_nml.gt.1) then
          target_depth(ibnd,ndepth_nml) = plate_thickness(ibnd)
        else
          target_depth(ibnd,ndepth_nml) = 0.0_R8
        endif
      enddo

cxpb Correct the target_width for non-parallel but neighbouring wall elements
      do ix = -1, nx
        do iy = -1, ny
          ibnd=xymap(ix,iy)
          if (ibnd.ne.0) then
            if (topiy(ix,iy).eq.ny+1) then        ! North boundary 
              ibndm=xymap( leftix(ix,iy), leftiy(ix,iy))
              ibndp=xymap(rightix(ix,iy),rightiy(ix,iy))
            else if (bottomiy(ix,iy).eq.-2) then  ! South boundary
              ibndp=xymap( leftix(ix,iy), leftiy(ix,iy))
              ibndm=xymap(rightix(ix,iy),rightiy(ix,iy))
            endif
            if (leftix(ix,iy).eq.-2) then         ! West boundary
              ibndm=xymap(bottomix(ix,iy),bottomiy(ix,iy))
              ibndp=xymap(   topix(ix,iy),   topiy(ix,iy))
            else if (rightix(ix,iy).eq.nx+1) then ! East boundary
              ibndp=xymap(bottomix(ix,iy),bottomiy(ix,iy))
              ibndm=xymap(   topix(ix,iy),   topiy(ix,iy))
            endif
            normal_r = -(z1(ibnd)-z2(ibnd))
            normal_z = -(r2(ibnd)-r1(ibnd))
            dist = sqrt((r2(ibnd)-r1(ibnd))**2+(z2(ibnd)-z1(ibnd))**2)
            target_width(ibnd,1:ndepth_nml) = dist
            normal_r = normal_r / dist
            normal_z = normal_z / dist
            if (ibndm.ne.0) then
              n1_r = -(z1(ibndm)-z2(ibndm))
              n1_z = -(r2(ibndm)-r1(ibndm))
              dist = sqrt((r2(ibndm)-r1(ibndm))**2+
     +                    (z2(ibndm)-z1(ibndm))**2)
              n1_r = n1_r / dist
              n1_z = n1_z / dist
              sina = normal_r * n1_z - normal_z * n1_r
              cosa = normal_r * n1_r + normal_z * n1_z
              cosa = max(-1.0_R8, min(1.0_R8, cosa))
              call xertst (cosa.ne.0.0_R8,
     .         'perpendicular neighbouring wall elements !')
              if(cosa.eq.1.0_R8) then
                coshlf = 1.0_R8
                sinhlf = 0.0_R8
              else if(sina.ge.0.0_R8) then
                coshlf = sqrt((1.0_R8 + cosa)/2.0_R8)
                sinhlf = sqrt((1.0_R8 - cosa)/2.0_R8)
              else
                coshlf = sqrt((1.0_R8 + cosa)/2.0_R8)
                sinhlf =-sqrt((1.0_R8 - cosa)/2.0_R8)
              endif
              tanhlf = sinhlf / coshlf
cxpb  rotate with an angle -hlf
              nhalf_r = normal_r * coshlf - normal_z * sinhlf
              nhalf_z = normal_r * sinhlf + normal_z * coshlf
              do i = 1, ndepth_nml
                target_width(ibnd,i) = max(0.0_R8, target_width(ibnd,i) 
     .           - tanhlf * target_depth(ibnd,i))
                if (target_width(ibnd,i).gt.0.0_R8) then
                  targetcr(ibnd,i,1) = r1(ibnd) 
     .             - target_depth(ibnd,i) * nhalf_r 
                  targetcz(ibnd,i,1) = z1(ibnd) 
     .             - target_depth(ibnd,i) * nhalf_z
                else
                  targetcr(ibnd,i,1) = (r1(ibnd) + r2(ibnd))/2.0_R8
     .             - target_depth(ibnd,i) * normal_r
                  targetcz(ibnd,i,1) = (z1(ibnd) + z2(ibnd))/2.0_R8
     .             - target_depth(ibnd,i) * normal_z
                endif 
              enddo
            else  ! case of edge cell with no neighbour on this side
              do i = 1, ndepth_nml
                targetcr(ibnd,i,1) = r1(ibnd)
     .           - target_depth(ibnd,i) * normal_r
                targetcz(ibnd,i,1) = z1(ibnd)
     .           - target_depth(ibnd,i) * normal_z
              enddo
            endif 
            if (ibndp.ne.0) then
              n2_r = -(z1(ibndp)-z2(ibndp))
              n2_z = -(r2(ibndp)-r1(ibndp))
              dist = sqrt((r2(ibndp)-r1(ibndp))**2+
     +                    (z2(ibndp)-z1(ibndp))**2)
              n2_r = n2_r / dist
              n2_z = n2_z / dist
              sina = normal_r * n2_z - normal_z * n2_r
              cosa = normal_r * n2_r + normal_z * n2_z
              cosa = max(-1.0_R8, min(1.0_R8, cosa))
              call xertst (cosa.ne.0.0_R8,
     .         'perpendicular neighbouring wall elements !')
              if(cosa.eq.1.0_R8) then
                coshlf = 1.0_R8
                sinhlf = 0.0_R8
              else if(sina.ge.0.0_R8) then
                coshlf = sqrt((1.0_R8 + cosa)/2.0_R8)
                sinhlf = sqrt((1.0_R8 - cosa)/2.0_R8)
              else
                coshlf = sqrt((1.0_R8 + cosa)/2.0_R8)
                sinhlf =-sqrt((1.0_R8 - cosa)/2.0_R8)
              endif
              tanhlf = sinhlf / coshlf
cxpb  rotate with an angle -hlf
              nhalf_r = normal_r * coshlf - normal_z * sinhlf
              nhalf_z = normal_r * sinhlf + normal_z * coshlf
              do i = 1, ndepth_nml
                target_width(ibnd,i) = max(0.0_R8, target_width(ibnd,i) 
     .           - tanhlf * target_depth(ibnd,i))
                if (target_width(ibnd,i).gt.0.0_R8) then
                  targetcr(ibnd,i,2) = r2(ibnd) 
     .             - target_depth(ibnd,i) * nhalf_r 
                  targetcz(ibnd,i,2) = z2(ibnd) 
     .             - target_depth(ibnd,i) * nhalf_z 
                else
                  targetcr(ibnd,i,2) = (r1(ibnd) + r2(ibnd))/2.0_R8
     .             - target_depth(ibnd,i) * normal_r
                  targetcz(ibnd,i,2) = (z1(ibnd) + z2(ibnd))/2.0_R8
     .             - target_depth(ibnd,i) * normal_z
                endif 
              enddo
            else  ! case of edge cell with no neighbour on this side
              do i = 1, ndepth_nml
                targetcr(ibnd,i,2) = r2(ibnd)
     .           - target_depth(ibnd,i) * normal_r
                targetcz(ibnd,i,2) = z2(ibnd)
     .           - target_depth(ibnd,i) * normal_z
              enddo
            endif 
          endif
        enddo
      enddo

cxpb Compute any missing temperature profiles if necessary
cxpb or change them according to plate_option setting
      if (ndepth_nml.gt.1.and.plate_model.ge.1) then
        do ibnd = 1, nbnd
          call xertst(plate_thickness(ibnd).gt.0.0_R8,
     &     'Incomplete PLATE_THICKNESS array!')
          if (plate_option.eq.1) then
cxpb Set constant temperature
            target_temp(ibnd,:) = backplate_temp(ibnd)
          else if (plate_option.eq.2 .or.
     &      target_temp(ibnd,ndepth_nml).eq.0.0_R8) then
cxpb Compute incoming heat flux
            if (bottomiy(imapx(ibnd),imapy(ibnd)).eq.-2) then
cxpb          South boundary
              qin = max(0.0_R8,
     &           -fhi(topix(imapx(ibnd),imapy(ibnd)),
     &                topiy(imapx(ibnd),imapy(ibnd)),1))/
     &            gs(topix(imapx(ibnd),imapy(ibnd)),
     &               topiy(imapx(ibnd),imapy(ibnd)),1)
              do is=0,ns-1
                qin=qin+max(0.0_R8,
     &             -fna(topix(imapx(ibnd),imapy(ibnd)),
     &                  topiy(imapx(ibnd),imapy(ibnd)),1,is))
     &             /gs(topix(imapx(ibnd),imapy(ibnd)),
     &                 topiy(imapx(ibnd),imapy(ibnd)),1)*
     &             (kinrgy(imapx(ibnd),imapy(ibnd),is)*(1.0_R8-BoRiS) + 
     &              rpt(imapx(ibnd),imapy(ibnd),is)*ev)
              enddo
              qel = max(0.0_R8,
     &           -fhe(topix(imapx(ibnd),imapy(ibnd)),
     &                topiy(imapx(ibnd),imapy(ibnd)),1))/
     &            gs(topix(imapx(ibnd),imapy(ibnd)),
     &               topiy(imapx(ibnd),imapy(ibnd)),1)
            else if (topiy(imapx(ibnd),imapy(ibnd)).eq.ny+1) then
cxpb          North boundary
              qin = max(0.0_R8,fhi(imapx(ibnd),imapy(ibnd),1))
     &           /gs(imapx(ibnd),imapy(ibnd),1)
              do is=0,ns-1
                qin=qin+max(0.0_R8,fna(imapx(ibnd),imapy(ibnd),1,is))
     &             /gs(imapx(ibnd),imapy(ibnd),1)*
     &             (kinrgy(imapx(ibnd),imapy(ibnd),is)*(1.0_R8-BoRiS) +  
     &              rpt(imapx(ibnd),imapy(ibnd),is)*ev)
              enddo
              qel = max(0.0_R8,fhe(imapx(ibnd),imapy(ibnd),1))
     &           /gs(imapx(ibnd),imapy(ibnd),1)
            else if (leftix(imapx(ibnd),imapy(ibnd)).eq.-2) then
cxpb          West boundary
              qin = max(0.0_R8,-fhi(rightix(imapx(ibnd),imapy(ibnd)),
     &                          rightiy(imapx(ibnd),imapy(ibnd)),0))
     &           /gs(rightix(imapx(ibnd),imapy(ibnd)),
     &               rightiy(imapx(ibnd),imapy(ibnd)),0)
              do is=0,ns-1
                qin=qin+max(0.0_R8,
     &             -fna(rightix(imapx(ibnd),imapy(ibnd)),
     &                  rightiy(imapx(ibnd),imapy(ibnd)),0,is))
     &              /gs(rightix(imapx(ibnd),imapy(ibnd)),
     &                  rightiy(imapx(ibnd),imapy(ibnd)),0)*
     &             (kinrgy(imapx(ibnd),imapy(ibnd),is)*(1.0_R8-BoRiS) +  
     &              rpt(imapx(ibnd),imapy(ibnd),is)*ev)
              enddo
              qel = max(0.0_R8,
     &           -fhe(rightix(imapx(ibnd),imapy(ibnd)),
     &                rightiy(imapx(ibnd),imapy(ibnd)),0))
     &           /gs(rightix(imapx(ibnd),imapy(ibnd)),
     &               rightiy(imapx(ibnd),imapy(ibnd)),0)
            else if (rightix(imapx(ibnd),imapy(ibnd)).eq.nx+1) then
cxpb          East boundary
              qin = max(0.0_R8,fhi(imapx(ibnd),imapy(ibnd),0))
     &           /gs(imapx(ibnd),imapy(ibnd),0)
              do is=0,ns-1
                qin = qin + 
     &             max(0.0_R8,fna(imapx(ibnd),imapy(ibnd),0,is))
     &             /gs(imapx(ibnd),imapy(ibnd),0)*
     &             (kinrgy(imapx(ibnd),imapy(ibnd),is)*(1.0_R8-BoRiS) +  
     &              rpt(imapx(ibnd),imapy(ibnd),is)*ev)
              enddo
              qel = max(0.0_R8,fhe(imapx(ibnd),imapy(ibnd),0))
     &           /gs(imapx(ibnd),imapy(ibnd),0)
            endif
cxpb Use 1-D equilibrium profile
            call heatdiff1D(ibnd,ns,qin+qel,0.0_R8,
     .       target_temp(ibnd,1:ndepth_nml),output)
          endif   ! plate_option
        enddo   !  ibnd
      endif   !  ndepth_nml
      call subend()
      return

* ..finalise time-dependent wall data
      entry wall_final(nx,ny,ns,plate_model,plate_temp,plate_thick)
      default_cooling = .true.
      do ibnd = 1, nbnd
        default_cooling = default_cooling .and. 
     &   .not.inertial_cooling(ibnd)
      enddo
      wall_namelist_used = wall_namelist_used .or. ntargsp.gt.0 .or.
     .  layer_thermal_contact.ne.1.0_R8 .or. plate_model.ge.1 .or.
     .  maxval(target_temp(1:nbnd,1)).gt.plate_temp .or.
     .  maxval(thermal_evaporation(1:nbnd,0:ns-1)).gt.0.0_R8 .or.
     .  maxval(physical_sputtering(1:nbnd,0:ns-1)).gt.0.0_R8 .or.
     .  maxval(chemical_sputtering(1:nbnd,0:ns-1)).gt.0.0_R8 .or.
     .  maxval(RES_sputtering(1:nbnd,0:ns-1)).gt.0.0_R8 .or.
     .  maxval(backscattering(1:nbnd,0:ns-1)).gt.0.0_R8 .or.
     .  .not.default_cooling
      if (.not.wall_namelist_used) return
      open(99,file='../b2.wall_save.parameters')
      write(99,'(a10)') '&WALL_SAVE'
      do itrack=1,nalloys
        if(layer_alloys(itrack).ne.' ') write(99,'(a13,i1,a4,a6,a2)')
     .    'LAYER_ALLOYS(',itrack,')= ''',
     .     layer_alloys(itrack),''','
      enddo
      write(99,'(a11,i4,a1)') 'NDEPTH_NML=',ndepth_nml,','
      do itrack=1,ntargsp
        write(99,'(a14,i2,a3,a2,a2)') 
     .   'TRACK_SPECIES(',itrack,')=''',track_species(itrack),''','
      enddo
      if (layer_thermal_contact.ne.1.0_R8) write(99,'(a22,1es15.8,a1)') 
     > 'LAYER_THERMAL_CONTACT=',layer_thermal_contact,','
      if (.not.(plate_model.ge.1.and.ndepth_nml.gt.1)) then
	do ibnd = 1, nbnd
	  if (target_temp(ibnd,1).gt.plate_temp) 
     >	   write(99,'(a12,i4,a4,1x,1es15.8,a1)') 
     .	   'TARGET_TEMP(',ibnd,',1)=',target_temp(ibnd,1),','
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
            write(99,'(a)') trim(line)
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
           write(99,'(a)') trim(line)
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
            write(99,'(a)') trim(line)
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
           write(99,'(a)') trim(line)
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
            write(99,'(a)') trim(line)
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
            write(99,'(a)') trim(line)
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
            write(99,'(a)') trim(line)
          enddo
        endif
      enddo
      if (plate_model.ge.1 .or. wall_namelist_used) then
cxpb  Use explicit writing of non-default elements only
        do ix = -1, nx
          do iy = -1, ny
            if (xymap(ix,iy).ne.0) write(99,'(a6,i4,a1,i4,a2,i4,a1)')
     .       'XYMAP(',ix,',',iy,')=',xymap(ix,iy),','     
          end do
        end do
        do ibnd = 1, nbnd
          write(99,'(a6,i4,a2,i4,a1,1x,a6,i4,a2,i4,a1)') 
     .     'IMAPX(',ibnd,')=',imapx(ibnd),',',
     .     'IMAPY(',ibnd,')=',imapy(ibnd),','
          if (inertial_cooling(ibnd)) 
     .     write(99,'(a17,i4,a2,l1,a1)')
     .     'INERTIAL_COOLING(',ibnd,')=',inertial_cooling(ibnd),','
          if (backplate_temp(ibnd).ne.plate_temp) 
     .     write(99,'(a15,i4,a2,1es15.8,a1)') 
     .     'BACKPLATE_TEMP(',ibnd,')=',backplate_temp(ibnd),','
          if (plate_thickness(ibnd).ne.plate_thick)
     .     write(99,'(a16,i4,a2,1es15.8,a1)') 
     .     'PLATE_THICKNESS(',ibnd,')=',plate_thickness(ibnd),','
          if (coating_thickness(ibnd).gt.0.0_R8)
     .     write(99,'(a18,i4,a2,1es15.8,a1)') 
     .     'COATING_THICKNESS(',ibnd,')=',coating_thickness(ibnd),','
          if (plate_time_factor(ibnd).ne.1.0_R8)
     .     write(99,'(a18,i4,a2,1es15.8,a1)') 
     .     'PLATE_TIME_FACTOR(',ibnd,')=',plate_time_factor(ibnd),','
          if (trim(surface_material_name(ibnd)).ne.'C') 
     .     write(99,'(a22,i4,a4,a6,a2)') 
     .     'SURFACE_MATERIAL_NAME(',ibnd,')= ''',
     .      surface_material_name(ibnd),''','
          if (coating_thickness(ibnd).gt.0.0_R8) 
     .     write(99,'(a22,i4,a4,a6,a2)') 
     .     'COATING_MATERIAL_NAME(',ibnd,')= ''',
     .      coating_material_name(ibnd),''','
          if (trim(bulk_material_name(ibnd)).ne.'C') 
     .     write(99,'(a19,i4,a4,a6,a2)') 
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
            write(99,'(a20,i4,a2,i2,a1)')
     .       'LAYER_NCONSTITUENTS(',ibnd,')=',
     .        layer_nconstituents(ibnd),','
            do i = 1, layer_nconstituents(ibnd)
              write(99,'(a21,i4,a1,i2,a4,i2,a1)')
     .         'LAYER_NZCONSTITUENTS(',ibnd,',',i,')=  ',
     .          layer_nzconstituents(ibnd,i),','
              write(99,'(a23,i4,a1,i2,a2,1es15.8,a1)')
     .         'LAYER_NRELCONSTITUENTS(',ibnd,',',i,')=',
     .          layer_nrelconstituents(ibnd,i),','
            end do
          end if
          do itrack=1,ntargsp
            if (deposition(ibnd,itrack).ne.0.0_R8)
     .       write(99,'(a11,i4,a1,i1,a2,1es15.8,a1)')
     .        'DEPOSITION(',ibnd,',',itrack,')=',
     1         deposition(ibnd,itrack),','
            if (monolayer_deposition(ibnd,itrack).ne.0.0_R8)
     .       write(99,'(a,i4,a,i1,a,1es15.8,a)') 
     .        'MONOLAYER_DEPOSITION(',ibnd,',',itrack,')=',
     1         monolayer_deposition(ibnd,itrack),','
            if (erosion(ibnd,itrack).ne.0.0_R8)
     .       write(99,'(a8,i4,a1,i1,a2,1es15.8,a1)')
     .        'EROSION(',ibnd,',',itrack,')=',
     1         erosion(ibnd,itrack),','
            if (monolayer_erosion(ibnd,itrack).ne.0.0_R8)
     .       write(99,'(a,i4,a,i1,a,1es15.8,a)') 
     .        'MONOLAYER_EROSION(',ibnd,',',itrack,')=',
     1         monolayer_erosion(ibnd,itrack),','
          enddo
          if (plate_time(ibnd).ne.0.0_R8)
     .     write(99,'(a11,i4,a2,1es15.8,a1)')
     .      'PLATE_TIME(',ibnd,')=',plate_time(ibnd),','
          if (plate_area(ibnd).ne.0.0_R8)
     .     write(99,'(a11,i4,a2,1es15.8,a1)')
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
              write(99,'(a)') trim(line)
            enddo
          enddo
        endif
      endif
      write(99,'(a1)') '/'
      close(99)
      return
        
      end subroutine init_wall

      subroutine add_layer_thickness(nx,ny,ns,ibnd)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      use b2mod_b2cmpa
      use b2mod_sputter
      use b2mod_indirect
      implicit none
      integer nx, ny, ns, nz, ibnd
      integer i, is, ix, iy, ibndm, ibndp, iend, itrack
      real (kind=R8) :: normal_r, normal_z, dist, 
     .  cosa, sina, coshlf, sinhlf, tanhlf,
     .  n1_r, n1_z, n2_r, n2_z, nhalf_r, nhalf_z, sum,
     .  sput_frac(0:ntrack)
      logical surface_only, found
      character*24 filename
      integer, save :: ncall, hd_linlog, sput_frac_flag
      real (kind=R8), save :: hd_ratio,
     .  sput_mixed_alpha, sput_mixed_beta
      data ncall/0/, hd_linlog/1/, hd_ratio/100.0_R8/,
     .     sput_mixed_alpha/1.0_R8/, sput_mixed_beta/1.0_R8/,
     2     sput_frac_flag/0/

      real (kind=R8) :: mass_density
      logical open_file
      intrinsic max, sqrt
      external ipgeti, ipgetr, xertst, xerrab
      external mass_density, open_file

      if (ncall.eq.0) then
        call ipgeti('heatdiff1D_linlog', hd_linlog)
        call xertst(hd_linlog.eq.1 .or. hd_linlog.eq.2,
     .   'faulty parameter heatdiff1D_linlog')
        if (hd_linlog.eq.2) then
          call ipgetr('heatdiff1D_ratio', hd_ratio)
          call xertst(hd_ratio.gt.1.0_R8, 
     .     'faulty parameter heatdiff1D_ratio')
        endif
        call ipgeti ('b2stbr_sput_frac_flag', sput_frac_flag)
        if(sput_frac_flag.ne.0) then
          call ipgetr ('b2stbr_sput_mixed_alpha', sput_mixed_alpha)
          call xertst (sput_mixed_alpha.ge.0.0_R8,
     .     'faulty parameter sput_mixed_alpha')
          call ipgetr ('b2stbr_sput_mixed_beta', sput_mixed_beta)
          call xertst (sput_mixed_beta.gt.0.0_R8,
     .     'faulty parameter sput_mixed_beta')
        endif
        ncall = ncall + 1
      endif
        
      if (layer_thickness(ibnd).eq.0.0_R8) then
        continue  ! nothing to do
      else if (layer_thickness(ibnd).gt.0.0_R8) then
cxpb  The deposited layer is included in the first depth element, raising it.
cxpb  Recall the depth is measured starting at 0 at the nominal surface and
cxpb  increases going away from the plasma.
        target_depth(ibnd,1) = -layer_thickness(ibnd)
        surface_only = .true.
      else if (layer_thickness(ibnd).lt.0.0_R8) then
cxpb  Shallow erosion case, does not cross eventual coating boundary
cxpb  The net erosion is layer_thickness but the mixed layer is
cxpb  layer_depth thick, which all must be kept within the surface 
cxpb  target element
        if ((layer_depth(ibnd)-layer_thickness(ibnd)).lt.
     .      target_depth(ibnd,2).and.
     .      ((layer_depth(ibnd)-layer_thickness(ibnd)).lt.
     .       coating_thickness(ibnd) .or.
     .       coating_thickness(ibnd).eq.0.0_R8)) then
          target_depth(ibnd,1) = -layer_thickness(ibnd)
          surface_only = .true.
        else if ((layer_depth(ibnd)-layer_thickness(ibnd)).lt.
     .       coating_thickness(ibnd) .or.
     .       coating_thickness(ibnd).eq.0.0_R8) then
cxpb  Deep erosion case, still within eventual coating boundary
cxpb  The first target element must be at least layer_depth thick
          target_depth(ibnd,1) = -layer_thickness(ibnd)
          target_depth(ibnd,2) = layer_depth(ibnd)-layer_thickness(ibnd)
          if (hd_linlog.eq.1) then
            target_depth(ibnd,2) = max(target_depth(ibnd,2),
     .       (plate_thickness(ibnd)-layer_thickness(ibnd))/
     .        float(ndepth_nml-1))
            do i = 3, ndepth_nml-1
              target_depth(ibnd,i) = target_depth(ibnd,2) +
     .         (plate_thickness(ibnd)-target_depth(ibnd,2)) * 
     .         float(i-2) / float(ndepth_nml - 2)
            end do
          else
            target_depth(ibnd,2) = max(target_depth(ibnd,2),
     .       (plate_thickness(ibnd)-layer_thickness(ibnd)) * 
     .         (hd_ratio ** (1.0_R8 / float(ndepth_nml - 1)) - 1.0_R8)/
     .         (hd_ratio - 1.0_R8))
            do i = 3, ndepth-1
              target_depth(ibnd,i) = target_depth(ibnd,2) +
     .         (plate_thickness(ibnd)-target_depth(ibnd,2)) * 
     .         (hd_ratio ** (float(i-2) / float(ndepth_nml - 2)) 
     .                                                      - 1.0_R8)/
     .         (hd_ratio - 1.0_R8)
            end do
          endif
          surface_only = .false.
        else
cxpb  Very deep erosion case, crosses coating boundary
cxpb  The first target element must be at least layer_depth thick
          target_depth(ibnd,1) = -layer_thickness(ibnd)
          target_depth(ibnd,2) = layer_depth(ibnd)-layer_thickness(ibnd)
          if (hd_linlog.eq.1) then
            target_depth(ibnd,2) = max(target_depth(ibnd,2),
     .       (plate_thickness(ibnd)-layer_thickness(ibnd))/
     .        float(ndepth_nml-1))
            do i = 3, ndepth-1
              target_depth(ibnd,i) = target_depth(ibnd,2) +
     .         (plate_thickness(ibnd)-target_depth(ibnd,2)) * 
     .         float(i-2) / float(ndepth_nml - 2)
            end do
          else
            target_depth(ibnd,2) = max(target_depth(ibnd,2),
     .       (plate_thickness(ibnd)-layer_thickness(ibnd)) * 
     .         (hd_ratio ** (1.0_R8 / float(ndepth_nml - 1)) - 1.0_R8)/
     .         (hd_ratio - 1.0_R8))
            do i = 3, ndepth-1
              target_depth(ibnd,i) = target_depth(ibnd,2) +
     .         (plate_thickness(ibnd)-target_depth(ibnd,2)) * 
     .         (hd_ratio ** (float(i-2) / float(ndepth_nml - 2)) 
     .                                                       - 1.0_R8)/
     .         (hd_ratio - 1.0_R8)
            end do
          endif
          surface_only = .false.
          write(*,*) 
     .     'Fully eroded coating. Read new bulk surface properties'
          filename='Surface_properties/'//trim(bulk_material_name(ibnd))
          if(open_file(99,filename,'data.local|data')) then
            read(99,*) surface_long_name(ibnd)
            read(99,*) matsurf_nconstituents(ibnd)
            read(99,*) matsurf_nzconstituents
     .       (ibnd,1:matsurf_nconstituents(ibnd))
            read(99,*) matsurf_nrelconstituents
     .       (ibnd,1:matsurf_nconstituents(ibnd))
            read(99,*) matsurf_chemsput_erel(ibnd)
            read(99,*) matsurf_blackbody_emiss(ibnd)
            read(99,*) matsurf_vapprs_A(ibnd)
            read(99,*) matsurf_vapprs_B(ibnd)
            read(99,*) matsurf_stick_coef(ibnd)
            read(99,*) matsurf_vapor_temp(ibnd)
            read(99,*) matsurf_heat_sublim(ibnd)
            close(99)
          else
            call xerrab ('Deep bulk material surface file not found!')
          endif
c-xpb Normalise surface constituents
          sum = 0.0_R8
          do i = 1, matsurf_nconstituents(ibnd)
            sum = sum + matsurf_nrelconstituents(ibnd,i)
          end do
          do i = 1, matsurf_nconstituents(ibnd)
            matsurf_nrelconstituents(ibnd,i) = 
     .       matsurf_nrelconstituents(ibnd,i)/sum
          end do
        endif  
      endif
cxpb  Modify the target_width, targetcr and targetcz arrays as needed
cxpb  Code similar to int_wall above
      if (surface_only .or. layer_thickness(ibnd).eq.0.0_R8) then
        iend=1
      else
        iend=ndepth_nml-1
      endif
      ix = imapx(ibnd)
      iy = imapy(ibnd)
      normal_r = -(z1(ibnd)-z2(ibnd))
      normal_z = -(r2(ibnd)-r1(ibnd))
      dist = sqrt((r2(ibnd)-r1(ibnd))**2+(z2(ibnd)-z1(ibnd))**2)
      normal_r = normal_r / dist
      normal_z = normal_z / dist
      if (topiy(ix,iy).eq.ny+1) then        ! North boundary 
        ibndm=xymap( leftix(ix,iy), leftiy(ix,iy))
        ibndp=xymap(rightix(ix,iy),rightiy(ix,iy))
      else if (bottomiy(ix,iy).eq.-2) then  ! South boundary
        ibndp=xymap( leftix(ix,iy), leftiy(ix,iy))
        ibndm=xymap(rightix(ix,iy),rightiy(ix,iy))
      endif
      if (leftix(ix,iy).eq.-2) then         ! West boundary
        ibndm=xymap(bottomix(ix,iy),bottomiy(ix,iy))
        ibndp=xymap(   topix(ix,iy),   topiy(ix,iy))
      else if (rightix(ix,iy).eq.nx+1) then ! East boundary
        ibndp=xymap(bottomix(ix,iy),bottomiy(ix,iy))
        ibndm=xymap(   topix(ix,iy),   topiy(ix,iy))
      endif
      if (ibndm.ne.0 .and. layer_thickness(ibnd).ne.0.0_R8) then
        n1_r = -(z1(ibndm)-z2(ibndm))
        n1_z = -(r2(ibndm)-r1(ibndm))
        dist = sqrt((r2(ibndm)-r1(ibndm))**2+
     +              (z2(ibndm)-z1(ibndm))**2)
        n1_r = n1_r / dist
        n1_z = n1_z / dist
        sina = normal_r * n1_z - normal_z * n1_r
        cosa = normal_r * n1_r + normal_z * n1_z
        cosa = max(-1.0_R8, min(1.0_R8, cosa))
        call xertst (cosa.ne.0.0_R8,
     .   'perpendicular neighbouring wall elements !')
        if(cosa.eq.1.0_R8) then
          coshlf = 1.0_R8
          sinhlf = 0.0_R8
        else if(sina.ge.0.0_R8) then
          coshlf = sqrt((1.0_R8 + cosa)/2.0_R8)
          sinhlf = sqrt((1.0_R8 - cosa)/2.0_R8)
        else
          coshlf = sqrt((1.0_R8 + cosa)/2.0_R8)
          sinhlf =-sqrt((1.0_R8 - cosa)/2.0_R8)
        endif
        tanhlf = sinhlf / coshlf
cxpb  rotate with an angle -hlf
        nhalf_r = normal_r * coshlf - normal_z * sinhlf
        nhalf_z = normal_r * sinhlf + normal_z * coshlf
        do i = 1, iend
          target_width(ibnd,i) = max(0.0_R8, 
     .     target_width(ibnd,i) 
     .           - tanhlf * target_depth(ibnd,i))
          if (target_width(ibnd,i).gt.0.0_R8) then
            targetcr(ibnd,i,1) = r1(ibnd) 
     .             - target_depth(ibnd,i) * nhalf_r 
            targetcz(ibnd,i,1) = z1(ibnd) 
     .             - target_depth(ibnd,i) * nhalf_z
          else
            targetcr(ibnd,i,1) = (r1(ibnd) + r2(ibnd))/2.0_R8
     .             - target_depth(ibnd,i) * normal_r
            targetcz(ibnd,i,1) = (z1(ibnd) + z2(ibnd))/2.0_R8
     .             - target_depth(ibnd,i) * normal_z
          endif 
        enddo
      else if (layer_thickness(ibnd).ne.0.0_R8) then ! case of edge cell with no neighbour on this side
        do i = 1, iend
          targetcr(ibnd,i,1) = r1(ibnd)
     .           - target_depth(ibnd,i) * normal_r
          targetcz(ibnd,i,1) = z1(ibnd)
     .           - target_depth(ibnd,i) * normal_z
        enddo
      endif
      if (ibndp.ne.0 .and. layer_thickness(ibnd).ne.0.0_R8) then
        n2_r = -(z1(ibndp)-z2(ibndp))
        n2_z = -(r2(ibndp)-r1(ibndp))
        dist = sqrt((r2(ibndp)-r1(ibndp))**2+
     +              (z2(ibndp)-z1(ibndp))**2)
        n2_r = n2_r / dist
        n2_z = n2_z / dist
        sina = normal_r * n2_z - normal_z * n2_r
        cosa = normal_r * n2_r + normal_z * n2_z
        cosa = max(-1.0_R8, min(1.0_R8, cosa))
        call xertst (cosa.ne.0.0_R8,
     .   'perpendicular neighbouring wall elements !')
        if(cosa.eq.1.0_R8) then
          coshlf = 1.0_R8
          sinhlf = 0.0_R8
        else if(sina.ge.0.0_R8) then
          coshlf = sqrt((1.0_R8 + cosa)/2.0_R8)
          sinhlf = sqrt((1.0_R8 - cosa)/2.0_R8)
        else
          coshlf = sqrt((1.0_R8 + cosa)/2.0_R8)
          sinhlf =-sqrt((1.0_R8 - cosa)/2.0_R8)
        endif
        tanhlf = sinhlf / coshlf
cxpb  rotate with an angle -hlf
        nhalf_r = normal_r * coshlf - normal_z * sinhlf
        nhalf_z = normal_r * sinhlf + normal_z * coshlf
        do i = 1, iend
          target_width(ibnd,i) = max(0.0_R8, target_width(ibnd,i) 
     .     - tanhlf * target_depth(ibnd,i))
          if (target_width(ibnd,i).gt.0.0_R8) then
            targetcr(ibnd,i,2) = r2(ibnd) 
     .             - target_depth(ibnd,i) * nhalf_r 
            targetcz(ibnd,i,2) = z2(ibnd) 
     .             - target_depth(ibnd,i) * nhalf_z 
          else
            targetcr(ibnd,i,2) = (r1(ibnd) + r2(ibnd))/2.0_R8
     .             - target_depth(ibnd,i) * normal_r
            targetcz(ibnd,i,2) = (z1(ibnd) + z2(ibnd))/2.0_R8
     .             - target_depth(ibnd,i) * normal_z
          endif 
        enddo
      else if (layer_thickness(ibnd).ne.0.0_R8) then  ! case of edge cell with no neighbour on this side
        do i = 1, ndepth_nml
          targetcr(ibnd,i,2) = r2(ibnd)
     .           - target_depth(ibnd,i) * normal_r
          targetcz(ibnd,i,2) = z2(ibnd)
     .           - target_depth(ibnd,i) * normal_z
        enddo
      endif    

c-xpb We now identify the track index of the bulk or coating species
c-xpb below the re-deposited layer.
      nz = matsurf_nconstituents(ibnd)
      found = .false.
      i = 1
      do while (i.le.nz .and. .not.found)
        is = 0
        do while (is.lt.ns .and. .not.found)
          if (is_neutral(is)) then
            found = nint(zn(is)).eq.matsurf_nzconstituents(ibnd,i)
          end if
          if (.not.found) is = is + 1
        end do
        if (.not.found) i = i + 1
      end do
      if (found) then
        itrack_0(ibnd)=track_index(is)
      else
        itrack_0(ibnd)=0
      endif

c-xpb  Code similar to b2stbr_phys.recycle in order to initialise the
c-xpb  layer quantities
c the base material ("0") is assumed to have an unnormalised weight 
c sput_mixed_alpha (default 1)
      sput_frac(0)=sput_mixed_alpha
c the remaining constituents each have an unnormalised weight equal
c to the number of monolayers deposited raised to the power 
c sput_mixed_beta (default 1)
      do itrack=1,ntargsp
        sput_frac(0)=sput_frac(0)+
     1   monolayer_deposition(ibnd,itrack)**sput_mixed_beta
      enddo
c the remaining constituents are then normalised by the sum of weights
      do itrack=1,ntargsp
        sput_frac(itrack)=
     1    monolayer_deposition(ibnd,itrack)**sput_mixed_beta/
     &    sput_frac(0)
      enddo
c the base material weight is calculated to ensure that the sum of 
c normalised weights is 1
      sput_frac(0)=1.0_R8
      do itrack=1,ntargsp
        sput_frac(0)=sput_frac(0)-sput_frac(itrack)
      enddo
c-xpb  Modify surface layer properties for backscattering
      layer_nconstituents(ibnd)=matsurf_nconstituents(ibnd)+ntargsp
c-xpb  Contributions from the bulk
      layer_nrelconstituents(ibnd,1:matsurf_nconstituents(ibnd))=
     1  matsurf_nrelconstituents
     .   (ibnd,1:matsurf_nconstituents(ibnd))*
     2  sput_frac(0)
c-xpb  Contributions from the deposited species
      do itrack=1,ntargsp
        layer_nzconstituents
     .     (ibnd,matsurf_nconstituents(ibnd)+itrack)=
     1     nint(zn(target_intcode(itrack)))
        layer_nrelconstituents
     .     (ibnd,matsurf_nconstituents(ibnd)+itrack)=
     .      sput_frac(itrack)
      enddo
      call compose_layer(ibnd,sput_frac)
c-xpb  Update layer mass density 
c-xpb  with bulk and deposited contributions
      layer_mass_density(ibnd)=mass_density(ibnd,2)*sput_frac(0)
      do itrack=1,ntargsp
        layer_mass_density(ibnd)=layer_mass_density(ibnd)+
     &    sput_frac(itrack)*(tgdns(itrack)*1.0e3_R8)
      enddo

      return
      end subroutine add_layer_thickness

      subroutine compose_layer(ibnd,sput_frac)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      use b2mod_b2cmpa
      use b2mod_sputter
      implicit none
      integer ibnd
      real (kind=R8) :: sput_frac(0:ntargsp)

      integer j, itrack, nz
      logical found
      real (kind=R8) :: mass_density
      external mass_density

      particle_inventory(ibnd,:)=0.0_R8
c-xpb  Compute layer composition
c-xpb  This will need refinement including alloy formation
      nz = matsurf_nconstituents(ibnd)
c-xpb  Contribution from the bulk
      layer_composition(ibnd,1:nz)=
     2  matsurf_nrelconstituents(ibnd,1:nz)*sput_frac(0)
      if (layer_nconstituents(ibnd).eq.nz) return
      if (ndepth_nml.gt.1) then
c-xpb  Inventory of bulk or coating particles
        if (coating_thickness(ibnd).gt.0.0_R8 .and.
     &      target_depth(ibnd,1).lt.coating_thickness(ibnd)) then
c-xpb  Coating case
          particle_inventory(ibnd,itrack_0(ibnd))=
     1     (coating_thickness(ibnd)-layer_depth(ibnd))*
     2      plate_area(ibnd)*mass_density(ibnd,2)/
     3     (atomic_mass(matsurf_nzconstituents(ibnd,1))*mp)
        else
c-xpb  Deep erosion or no coating case
          particle_inventory(ibnd,itrack_0(ibnd))=
     1     (plate_thickness(ibnd)-layer_depth(ibnd))*
     2      plate_area(ibnd)*mass_density(ibnd,2)/
     3     (atomic_mass(matsurf_nzconstituents(ibnd,1))*mp)
        endif
      endif
c-xpb  Additional contribution from deposited material
c-xpb  Must check if Z of deposited material already exists in bulk
      do itrack=1,ntargsp
        found = .false.
        j = 1
        do while (.not.found .and. j.le.nz)
          found = zn(target_intcode(itrack)).eq.
     .            matsurf_nzconstituents(ibnd,j)
          if (.not.found) j = j + 1
        enddo
        if (found) then
          layer_composition(ibnd,j)=
     .     layer_composition(ibnd,j)+sput_frac(itrack)
        else
          layer_composition(ibnd,nz+itrack)=
     .     layer_composition(ibnd,nz+itrack)+sput_frac(itrack) 
        endif
        if (sput_frac(itrack).gt.0.0_R8)
     >   particle_inventory(ibnd,itrack)= 
     .    particle_inventory(ibnd,itrack)+
     .    deposition(ibnd,itrack)-erosion(ibnd,itrack)
      enddo

      return
      end subroutine compose_layer

      function layer_Kappa(ibnd,T)
      use b2mod_types
      use b2mod_layer
      use b2mod_wall
      use b2mod_sputter
      implicit none
      integer ibnd
      integer nz, itrack
      real (kind=R8) :: T, layer_Kappa, a(4)
      real (kind=R8) :: bulk_value, coating_value, layer_value, frac

c-xpb First count the bulk contribution
      frac=layer_nrelconstituents(ibnd,1)/
     .   matsurf_nrelconstituents(ibnd,1)
      nz = matsurf_nconstituents(ibnd)
c-xpb The bulk or coating surface parameters
      if (coating_thickness(ibnd).gt.0.0_R8 .and.
     &    target_depth(ibnd,1).lt.coating_thickness(ibnd)) then
        a(1) = coating_Kappa_A1(ibnd)
        a(2) = coating_Kappa_A2(ibnd)
        a(3) = coating_Kappa_A3(ibnd)
        a(4) = coating_Kappa_A4(ibnd)
        coating_value = a(1) + a(2) * T + a(3) * T**2 + a(4) * T**3
        layer_Kappa = coating_value*frac
      else
        a(1) = matbulk_Kappa_A1(ibnd)
        a(2) = matbulk_Kappa_A2(ibnd)
        a(3) = matbulk_Kappa_A3(ibnd)
        a(4) = matbulk_Kappa_A4(ibnd)
        bulk_value = a(1) + a(2) * T + a(3) * T**2 + a(4) * T**3
        layer_Kappa = bulk_value*frac
      end if
      if (layer_nconstituents(ibnd).eq.nz) return
c-xpb Loop over deposited species contributions
      do itrack = 1, ntargsp
        frac = layer_nrelconstituents(ibnd,nz+itrack)
        a(1) = species_Kappa_A1(itrack)
        a(2) = species_Kappa_A2(itrack)
        a(3) = species_Kappa_A3(itrack)
        a(4) = species_Kappa_A4(itrack)
        layer_value = a(1) + a(2) * T + a(3) * T**2 + a(4) * T**3
        layer_Kappa = layer_Kappa + layer_value*frac
      end do

      return
      end function layer_Kappa
      
      function dlayer_KappabydT(ibnd,T)
      use b2mod_types
      use b2mod_layer
      use b2mod_wall
      use b2mod_sputter
      implicit none
      integer ibnd
      integer nz, itrack
      real (kind=R8) :: T, dlayer_KappabydT, a(2:4)
      real (kind=R8) :: bulk_value, coating_value, layer_value, frac

c-xpb First count the bulk contribution
      frac=layer_nrelconstituents(ibnd,1)/
     .   matsurf_nrelconstituents(ibnd,1)
      nz = matsurf_nconstituents(ibnd)
c-xpb The bulk or coating surface parameters
      if (coating_thickness(ibnd).gt.0.0_R8 .and.
     &    target_depth(ibnd,1).lt.coating_thickness(ibnd)) then
        a(2) = coating_Kappa_A2(ibnd)
        a(3) = coating_Kappa_A3(ibnd)
        a(4) = coating_Kappa_A4(ibnd)
        coating_value = a(2) + 2.0_R8 * a(3) * T + 3.0_R8 * a(4) * T**2
        dlayer_KappabydT = coating_value*frac
      else  
        a(2) = matbulk_Kappa_A2(ibnd)
        a(3) = matbulk_Kappa_A3(ibnd)
        a(4) = matbulk_Kappa_A4(ibnd)
        bulk_value = a(2) + 2.0_R8 * a(3) * T + 3.0_R8 * a(4) * T**2
        dlayer_KappabydT = bulk_value*frac
      end if
      if (layer_nconstituents(ibnd).eq.nz) return
c-xpb Loop over deposited species contributions
      do itrack = 1, ntargsp
        frac = layer_nrelconstituents(ibnd,nz+itrack)
        a(2) = species_Kappa_A2(itrack)
        a(3) = species_Kappa_A3(itrack)
        a(4) = species_Kappa_A4(itrack)
        layer_value = a(2) + 2.0_R8 * a(3) * T + 3.0_R8 * a(4) * T**2
        dlayer_KappabydT = dlayer_KappabydT + layer_value*frac
      end do

      return
      end function dlayer_KappabydT
      
      function layer_Cp(ibnd,T)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      use b2mod_sputter
      implicit none
      integer ibnd
      integer nz, itrack
      real (kind=R8) :: T, layer_Cp
      real (kind=R8) :: b(3)
      real (kind=R8) :: bulk_value, coating_value, layer_value, frac

c-xpb First count the bulk contribution
      frac=layer_nrelconstituents(ibnd,1)/
     .   matsurf_nrelconstituents(ibnd,1)
      nz = matsurf_nconstituents(ibnd)
c-xpb The bulk or coating surface parameters
      if (coating_thickness(ibnd).gt.0.0_R8 .and.
     &    target_depth(ibnd,1).lt.coating_thickness(ibnd)) then 
        b(1) = coating_Cp_B1(ibnd)
        b(2) = coating_Cp_B2(ibnd)
        b(3) = coating_Cp_B3(ibnd)
        coating_value = b(1) + b(2) * T - b(3) / T**2
        layer_Cp = coating_value * frac
      else
        b(1) = matbulk_Cp_B1(ibnd)
        b(2) = matbulk_Cp_B2(ibnd)
        b(3) = matbulk_Cp_B3(ibnd)
        bulk_value = b(1) + b(2) * T - b(3) / T**2
        layer_Cp = bulk_value * frac
      end if        
c-xpb Loop over deposited species contributions
      if (layer_nconstituents(ibnd).gt.nz) then
        do itrack = 1, ntargsp
          frac = layer_nrelconstituents(ibnd,nz+itrack)
          b(1) = species_Cp_B1(itrack)
          b(2) = species_Cp_B2(itrack)
          b(3) = species_Cp_B3(itrack)
          layer_value = b(1) + b(2) * T - b(3) / T**2
          layer_Cp = layer_Cp + layer_value * frac
        end do
      end if

c Converting from "Joules/gram/K" to "Joules/kg/K"
      layer_Cp = layer_Cp * 1.0e3_R8

      return
      end function layer_Cp
      
      function layer_emissivity(ibnd)
      use b2mod_types
      use b2mod_layer
      use b2mod_wall
      use b2mod_sputter
      implicit none
      integer ibnd
      integer nz, itrack
      real (kind=R8) :: layer_emissivity
      real (kind=R8) :: layer_value, frac

c-xpb First count the bulk contribution
      frac=layer_nrelconstituents(ibnd,1)/
     .   matsurf_nrelconstituents(ibnd,1)
      nz = matsurf_nconstituents(ibnd)
c-xpb The bulk or coating surface parameters
      layer_emissivity = matsurf_blackbody_emiss(ibnd)*frac
      if (layer_nconstituents(ibnd).eq.nz) return
c-xpb Loop over deposited species contributions
      do itrack = 1, ntargsp
        frac = layer_nrelconstituents(ibnd,nz+itrack)
        layer_value = species_blackbody_emiss(itrack)
        layer_emissivity = layer_emissivity + layer_value*frac
      end do

      return
      end function layer_emissivity
