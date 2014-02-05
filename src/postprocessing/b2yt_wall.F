      subroutine b2yt_wall(nx1,ny1,nx,ny,ns,ns1,is_neutral1,bundle_frac,
     *   plate_model,plate_temp,plate_thick)
      use b2mod_types
      use b2mod_geo2
      use b2mod_wall
      use b2mod_wall2
      use b2mod_layer
      use b2mod_indirect
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_sputter
#     include "DIMENSIONS.F"
      implicit none

      integer nx1, ny1, nx, ny, ns, ns1, nxtl, nxtr
      integer plate_model
      real (kind=R8) :: plate_temp, plate_thick
      real (kind=R8) :: bundle_frac(0:ns1-1,0:ns-1)
      logical is_neutral1(0:ns1-1)
      real (kind=R8), allocatable :: wall_data(:), wall_tmp(:) 

      integer ix, iy, iz, is, ibnd, ireg, i, i1, itrack
      integer ixa, ixb, iya, iyb, iza, ndepth1, nsm
      real (kind=R8) :: sum
      logical wall_namelist_used, layer_model
      save wall_namelist_used, layer_model
      data wall_namelist_used /.false./, layer_model/.false./
      character*400 line
      character*17 numstr
      logical found

      intrinsic mod, max, min, sqrt
      external xertst, subini, subend
      external convert_bundle_r, get_nxt

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

      call subini('b2yt_wall')
      call alloc_b2mod_wall(nx,ny)
cxpb Establish numbering of boundary wall cells
cxpb We must avoid corner cells and dead cells
      nbnd = 0
      ndepth_nml = 1
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
      layer_thermal_contact = 1.0_R8
      
      inquire(file='../b2.wall_save.parameters',
     1 exist=wall_namelist_used)
      if(wall_namelist_used) then
        open(99,file='../b2.wall_save.parameters')
        write(*,*) 'Reading ../b2.wall_save.parameters'
        read(99,wall_save)
        close(99)
        write(*,*) 'Read ../b2.wall_save.parameters'
        call xertst(ndepth_nml.le.ndepth .and. ndepth_nml.ge.1,
     .   'Invalid array size NDEPTH_NML!')
        ndepth1 = ndepth_nml
        call ipgeti ('b2ytdr_ndepth1', ndepth1)
        call xertst(ndepth1.le.ndepth .and. ndepth1.ge.1,
     .   'Invalid array size NDEPTH1!')
        call alloc_b2mod_wall2(nx1,ny1,ndepth1)
        if (plate_model.eq.0) then
          write(*,*) 
     .     'Using O-D time-independent plate temperature model'
        else if (plate_model.eq.1) then
          write(*,*) 'Using 1-D time-dependent plate temperature model'
          call xertst (ndepth1.gt.1, 
     .     'Inconsistent values of ndepth1 and plate_model')
        else if (plate_model.eq.2) then
          write(*,*) 'Using 2-D time-dependent plate temperature model'
          call xertst (ndepth1.gt.1, 
     .     'Inconsistent values of ndepth1 and plate_model')
        else
          call xerrab ('Invalid value for b2stbr_plate_model')
        endif
      endif

cxpb  We follow the same order as in b2plot
cxpb  First fill up the b2cmpb needed data
      cbirso = 0
      cbnrso = max(3*nncut,1)
      cbirno = cbirso+cbnrso
      cbnrno = max(nncut,1)
      cbirwe = cbirno+cbnrno
      cbnrwe = 1
      cbirea = cbirwe+cbnrwe
      cbnrea = 1
      call xertst (cbnrso+cbnrno+cbnrwe+cbnrea.le.cbnrmx,
     &  'faulty parameter cbnrmx')
      if (nncut.eq.1) then
        cbrbrk(cbirso) = real(leftcut1(1))/real(nx1)
        cbrbrk(cbirso+1) = real(rightcut1(1))/real(nx1) 
      else if (nncut.eq.2) then
        call get_nxt(nx,nxtl,nxtr)
        cbrbrk(cbirso) = real(leftcut1(1))/real(nx1)
        cbrbrk(cbirso+1) = real(rightcut1(2))/real(nx1) 
        cbrbrk(cbirso+2) = real(nxtl)/real(nx1) 
        cbrbrk(cbirso+3) = real(leftcut1(2))/real(nx1) 
        cbrbrk(cbirso+4) = real(rightcut1(2))/real(nx1) 
      endif
      cbrbrk(cbirso+cbnrso-1) = 2.00_R8
      if (nncut.eq.2) then
        cbrbrk(cbirno) = real(nxtl)/real(nx1) 
      endif
      cbrbrk(cbirno+cbnrno-1) = 2.00_R8
      cbrbrk(cbirwe:cbirwe+cbnrwe-1) = 2.00_R8
      cbrbrk(cbirea:cbirea+cbnrea-1) = 2.00_R8
cxpb  Northern boundaries (excluding diverted island)
      ireg = cbirno
      do iy = ny1, -1, -1
        do ix = -1, nx1
          if (topiy1(ix,iy).eq.ny1+1 .and. region1(ix,iy,0).ne.0 .and.
     &        region1(ix,iy,0).ne.5 .and.
     &        leftix1(ix,iy).ne.-2 .and. rightix1(ix,iy).ne.nx1+1) then
            nbnd = nbnd + 1
            imapx1(nbnd) = ix
            imapy1(nbnd) = iy
            xymap1(ix,iy)= nbnd
            plate_area1(nbnd) = gs1(ix,iy,1)
            if (xymap1(ix,iy).ne.0 .and. 
     &          xymap1(leftix1(ix,iy),leftiy1(ix,iy)).eq.0) then
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd
            endif 
          endif
          if(iy.eq.ny1) then
            if (cbrbrk(ireg)*nx1.lt.ix+0.5_R8) then
              ireg = ireg + 1
            endif
            if (xymap1(ix,iy).ne.0) imapreg(nbnd) = ireg
          endif
        enddo
      enddo
cxpb  Eastern boundaries
      ireg = cbirea+cbnrea-1
      do ix = nx1, -1, -1
        do iy = ny1, -1, -1
          if (rightix1(ix,iy).eq.nx1+1 .and. region1(ix,iy,0).ne.0 .and.
     &        bottomiy1(ix,iy).ne.-2 .and. topiy1(ix,iy).ne.ny1+1) then
            nbnd = nbnd + 1
            imapx1(nbnd) = ix
            imapy1(nbnd) = iy
            xymap1(ix,iy)= nbnd
            plate_area1(nbnd) = gs1(ix,iy,0)
            if (xymap1(ix,iy).ne.0 .and. 
     &          xymap1(topix1(ix,iy),topiy1(ix,iy)).eq.0) then
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd
            endif 
          endif
          if(ix.eq.nx1) then
            if (ireg.gt.cbirea) then
              if (cbrbrk(ireg-1)*ny1.gt.iy+0.5_R8) then
                ireg = ireg - 1
              endif
            endif
            if (xymap1(ix,iy).ne.0) imapreg(nbnd) = ireg
          endif
        enddo
      enddo
cxpb  Southern boundaries (excluding core)
      ireg = cbirso+cbnrso-1
      do iy = -1, ny1
        do ix = nx1, -1, -1      
          if (bottomiy1(ix,iy).eq.-2 .and. region1(ix,iy,0).ne.0 .and.
     &        mod(region1(ix,iy,0),4).ne.1 .and.
     &        leftix1(ix,iy).ne.-2 .and. rightix1(ix,iy).ne.nx1+1) then
            nbnd = nbnd + 1
            imapx1(nbnd) = ix
            imapy1(nbnd) = iy
            xymap1(ix,iy)= nbnd
            plate_area1(nbnd) = gs1(topix1(ix,iy),topiy1(ix,iy),1)
            if (xymap1(ix,iy).ne.0 .and. 
     &          xymap1(rightix1(ix,iy),rightiy1(ix,iy)).eq.0) then
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd
            endif 
          endif
          if(iy.eq.-1) then
            if (ireg.gt.cbirso) then
              if (cbrbrk(ireg-1)*nx1.gt.ix+0.5_R8) then
                ireg = ireg - 1
              endif
            endif
            if (xymap1(ix,iy).ne.0) imapreg(nbnd) = ireg
          endif
        enddo
      enddo
cxpb  Western boundaries
      ireg = cbirwe
      do ix = -1, nx1
        do iy = -1, ny1
          if (leftix1(ix,iy).eq.-2 .and. region1(ix,iy,0).ne.0 .and.
     &        bottomiy1(ix,iy).ne.-2 .and. topiy1(ix,iy).ne.ny1+1) then
            nbnd = nbnd + 1
            imapx1(nbnd) = ix
            imapy1(nbnd) = iy
            xymap1(ix,iy)= nbnd
            plate_area1(nbnd) = gs1(rightix1(ix,iy),rightiy1(ix,iy),0)
            if (xymap1(ix,iy).ne.0 .and. 
     &          xymap1(bottomix1(ix,iy),bottomiy1(ix,iy)).eq.0) then
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd
            endif 
          endif
          if(ix.eq.-1) then
            if (cbrbrk(ireg)*ny1.lt.iy+0.5_R8) then
              ireg = ireg + 1
            endif
            if (xymap1(ix,iy).ne.0) imapreg(nbnd) = ireg
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
 
cxpb Fill up array data not given in namelist
      do ibnd = 1, nwall
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
      call read_b2mod_wall(ntargsp)

cxpb For mixed materials, if unspecified, use normal material constituents
      do ibnd = 1, nwall
        if (layer_nconstituents(ibnd).eq.0) then
          layer_nconstituents(ibnd)=
     1     matsurf_nconstituents(ibnd)
          layer_nzconstituents(ibnd,1:layer_nconstituents(ibnd))=
     1     matsurf_nzconstituents(ibnd,1:layer_nconstituents(ibnd))
          layer_nrelconstituents(ibnd,1:layer_nconstituents(ibnd))=
     1     matsurf_nrelconstituents(ibnd,1:layer_nconstituents(ibnd))
        else
          layer_model = .true.
        end if
cxpb  normalise the layer constituents relative concentrations
        sum = 0.0_R8
        do i = 1, layer_nconstituents(ibnd)
          sum = sum + layer_nrelconstituents(ibnd,i)
        enddo
        if (sum.gt.0.0_R8)
     1   layer_nrelconstituents(ibnd,1:layer_nconstituents(ibnd))=
     2   layer_nrelconstituents(ibnd,1:layer_nconstituents(ibnd))/sum
      end do

cxpb Fill up the new wall arrays
      nsm = max(ns,ns1)
      allocate(wall_data(0:nsm-1+DEF_NATM+DEF_NMOL+DEF_NION))
      allocate(wall_tmp(0:nsm-1))
      do ibnd = 1, nbnd
        if (imapx1(ibnd).eq.-1) then
          ixa = -1
        else
          ixa = 2*nx-1-(4*nx1-1-2*imapx1(ibnd))*nx/(2*nx1)
        end if
        if (imapx1(ibnd).eq.nx1) then
          ixb = nx
        else
          ixb = (2*imapx1(ibnd)+2*nx1+1)*nx/(2*nx1)-nx
        end if
        if (imapy1(ibnd).eq.-1) then
          iya = -1
        else
          iya = 2*ny-1-(4*ny1-1-2*imapy1(ibnd))*ny/(2*ny1)
        end if
        if (imapy1(ibnd).eq.ny1) then
          iyb = ny
        else
          iyb = (2*imapy1(ibnd)+2*ny1+1)*ny/(2*ny1)-ny
        end if
cxpb Copy data arrays into new geometry
        do ix = ixa, ixb
          do iy = iya, iyb
            if (xymap(ix,iy).gt.0) then
              write(*,*) 'XYMAP(ix,iy) = ',xymap(ix,iy)
              surface_material_name1(ibnd)=
     .         surface_material_name(xymap(ix,iy))
              coating_material_name1(ibnd)=
     .         coating_material_name(xymap(ix,iy))
              bulk_material_name1(ibnd)=
     .         bulk_material_name(xymap(ix,iy))
              if (maxval(plate_thickness(1:nwall)).eq.plate_thick .and.
     &            minval(plate_thickness(1:nwall)).eq.plate_thick) then
                plate_thickness1(ibnd)=plate_thick
              else
                plate_thickness1(ibnd)=plate_thickness1(ibnd)+
     .           plate_thickness(xymap(ix,iy))*
     .           plate_area(xymap(ix,iy))/plate_area1(ibnd)
              endif
              coating_thickness1(ibnd)=coating_thickness1(ibnd)+
     .         coating_thickness(xymap(ix,iy))*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              if (maxval(plate_time_factor(1:nwall)).eq.1.0_R8 .and.
     &            minval(plate_time_factor(1:nwall)).eq.1.0_R8) then
                plate_time_factor1(ibnd)=1.0_R8
              else
                plate_time_factor1(ibnd)=plate_time_factor1(ibnd)+
     .           plate_time_factor(xymap(ix,iy))*
     .           plate_area(xymap(ix,iy))/plate_area1(ibnd)
              endif
              plate_time1(ibnd)=plate_time1(ibnd)+
     .         plate_time(xymap(ix,iy))*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              deposition1(ibnd,:)=deposition1(ibnd,:)+
     .         deposition(xymap(ix,iy),:)
              monolayer_deposition1(ibnd,:)=
     .         monolayer_deposition1(ibnd,:)+
     .         monolayer_deposition(xymap(ix,iy),:)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              erosion1(ibnd,:)=erosion1(ibnd,:)+
     .         erosion(xymap(ix,iy),:)
              monolayer_erosion1(ibnd,:)=monolayer_erosion1(ibnd,:)+
     .         monolayer_erosion(xymap(ix,iy),:)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION) =
     .         chemical_sputtering
     .         (xymap(ix,iy),ns:ns-1+DEF_NATM+DEF_NMOL+DEF_NION)
              wall_tmp(0:ns-1) = 
     .         chemical_sputtering(xymap(ix,iy),0:ns-1)
              call convert_bundle_r
     .         (ns1,ns,nsm,is_neutral1,bundle_frac,wall_tmp)
              chemical_sputtering1(ibnd,0:ns1-1) = 
     .         chemical_sputtering1(ibnd,0:ns1-1)+wall_tmp(0:ns1-1)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              chemical_sputtering1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)
     .         = chemical_sputtering1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)+
     .         wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION) =
     .         physical_sputtering
     .         (xymap(ix,iy),ns:ns-1+DEF_NATM+DEF_NMOL+DEF_NION)
              wall_tmp(0:ns-1) = 
     .         physical_sputtering(xymap(ix,iy),0:ns-1)
              call convert_bundle_r
     .         (ns1,ns,nsm,is_neutral1,bundle_frac,wall_tmp)
              physical_sputtering1(ibnd,0:ns1-1) = 
     .         physical_sputtering1(ibnd,0:ns1-1)+wall_tmp(0:ns1-1)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              physical_sputtering1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)
     .         = physical_sputtering1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)+
     .         wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION) =
     .         physical_sputtering_energy
     .         (xymap(ix,iy),ns:ns-1+DEF_NATM+DEF_NMOL+DEF_NION)
              wall_tmp(0:ns-1) = 
     .         physical_sputtering_energy(xymap(ix,iy),0:ns-1)
              call convert_bundle_r
     .         (ns1,ns,nsm,is_neutral1,bundle_frac,wall_tmp)
              physical_sputtering_energy1(ibnd,0:ns1-1) = 
     .         physical_sputtering_energy1(ibnd,0:ns1-1)+ 
     .         wall_tmp(0:ns1-1)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              physical_sputtering_energy1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)
     .         = physical_sputtering_energy1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)+
     .         wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION) =
     .         RES_sputtering
     .         (xymap(ix,iy),ns:ns-1+DEF_NATM+DEF_NMOL+DEF_NION)
              wall_tmp(0:ns-1) = 
     .         RES_sputtering(xymap(ix,iy),0:ns-1)
              call convert_bundle_r
     .         (ns1,ns,nsm,is_neutral1,bundle_frac,wall_tmp)
              RES_sputtering1(ibnd,0:ns1-1) = 
     .         RES_sputtering1(ibnd,0:ns1-1)+wall_tmp(0:ns1-1)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              RES_sputtering1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)
     .         = RES_sputtering1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)+
     .         wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION) =
     .         thermal_evaporation
     .         (xymap(ix,iy),ns:ns-1+DEF_NATM+DEF_NMOL+DEF_NION)
              wall_tmp(0:ns-1) = 
     .         thermal_evaporation(xymap(ix,iy),0:ns-1)
              call convert_bundle_r
     .         (ns1,ns,nsm,is_neutral1,bundle_frac,wall_tmp)
              thermal_evaporation1(ibnd,0:ns1-1) = 
     .         thermal_evaporation1(ibnd,0:ns1-1)+wall_tmp(0:ns1-1)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              thermal_evaporation1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)
     .         = thermal_evaporation1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)+
     .         wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION) =
     .         backscattering
     .         (xymap(ix,iy),ns:ns-1+DEF_NATM+DEF_NMOL+DEF_NION)
              wall_tmp(0:ns-1) = 
     .         backscattering(xymap(ix,iy),0:ns-1)
              call convert_bundle_r
     .         (ns1,ns,nsm,is_neutral1,bundle_frac,wall_tmp)
              backscattering1(ibnd,0:ns1-1) = 
     .         backscattering1(ibnd,0:ns1-1)+wall_tmp(0:ns1-1)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              backscattering1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)
     .         = backscattering1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)+
     .         wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION) =
     .         backscattering_energy
     .         (xymap(ix,iy),ns:ns-1+DEF_NATM+DEF_NMOL+DEF_NION)
              wall_tmp(0:ns-1) = 
     .         backscattering_energy(xymap(ix,iy),0:ns-1)
              call convert_bundle_r
     .         (ns1,ns,nsm,is_neutral1,bundle_frac,wall_tmp)
              backscattering_energy1(ibnd,0:ns1-1) = 
     .         backscattering_energy1(ibnd,0:ns1-1)+wall_tmp(0:ns1-1)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              backscattering_energy1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)
     .         = backscattering_energy1
     .         (ibnd,ns1:ns1-1+DEF_NATM+DEF_NMOL+DEF_NION)+
     .         wall_data(1:DEF_NATM+DEF_NMOL+DEF_NION)*
     .         plate_area(xymap(ix,iy))/plate_area1(ibnd)
              inertial_cooling1(ibnd)=.not.
     .         (.not.inertial_cooling1(ibnd).or.
     .          .not.inertial_cooling(xymap(ix,iy)))
              if (maxval(backplate_temp(1:nwall)).eq.plate_temp .and.
     &            minval(backplate_temp(1:nwall)).eq.plate_temp) then
                backplate_temp1(ibnd)=plate_temp
              else
                backplate_temp1(ibnd)=backplate_temp1(ibnd)+
     .           backplate_temp(xymap(ix,iy))*
     .           plate_area(xymap(ix,iy))/plate_area1(ibnd)
              endif
              do i = 1, layer_nconstituents(xymap(ix,iy))
                found = .false.
                do i1 = 1, layer_nconstituents1(ibnd)
                  if (layer_nzconstituents1(ibnd,i1).eq.
     .                layer_nzconstituents(xymap(ix,iy),i)) then
                    found = .true.
                    layer_nrelconstituents1(ibnd,i1)=
     .               layer_nrelconstituents1(ibnd,i1)+
     .               layer_nrelconstituents(xymap(ix,iy),i)*
     .               plate_area(xymap(ix,iy))/plate_area1(ibnd)
                  end if
                end do
                if (.not.found) then
                 layer_nconstituents1(ibnd)=layer_nconstituents1(ibnd)+1
                 layer_nzconstituents1(ibnd,layer_nconstituents1(ibnd))=
     .            layer_nzconstituents(xymap(ix,iy),i)
                 layer_nrelconstituents1(ibnd,
     .                                   layer_nconstituents1(ibnd))=
     .            layer_nrelconstituents(xymap(ix,iy),i)*
     .               plate_area(xymap(ix,iy))/plate_area1(ibnd)
                end if
              end do
              do iz = 1, ndepth1
                if (iz.eq.ndepth1) then
                  iza = ndepth_nml
                else
                  iza = ((iz-1)*ndepth_nml)/ndepth1+1
                endif
                if (maxval(target_temp(1:nwall,1:ndepth)).le.plate_temp)
     >           then
                  target_temp1(ibnd,iz)=plate_temp
                else
                  target_temp1(ibnd,iz)=target_temp1(ibnd,iz)+
     .             target_temp(xymap(ix,iy),iza)*
     .             plate_area(xymap(ix,iy))/plate_area(ibnd)
                endif
              end do
            end if
          end do
        end do          
      end do

      call subend()
      return

* ..finalise time-dependent wall data
      entry b2yt_wall_final(nx1,ny1,ns1,
     . plate_model,plate_temp,plate_thick)
      wall_namelist_used = wall_namelist_used .or. ntargsp.gt.0 .or.
     .  layer_thermal_contact.ne.1.0_R8 .or. plate_model.ge.1 .or.
     .  maxval(target_temp1(1:nbnd,1)).gt.plate_temp .or.
     .  maxval(thermal_evaporation1(1:nbnd,0:ns1-1)).gt.0.0_R8 .or.
     .  maxval(physical_sputtering1(1:nbnd,0:ns1-1)).gt.0.0_R8 .or.
     .  maxval(chemical_sputtering1(1:nbnd,0:ns1-1)).gt.0.0_R8 .or.
     .  maxval(RES_sputtering1(1:nbnd,0:ns1-1)).gt.0.0_R8 .or.
     .  maxval(backscattering1(1:nbnd,0:ns1-1)).gt.0.0_R8
      if (.not.wall_namelist_used) return
      open(99,file='../b2.wall_save.parameters.2')
      write(99,'(a10)') '&WALL_SAVE'
      do itrack=1,nalloys
        if(layer_alloys(itrack).ne.' ') write(99,'(a13,i1,a4,a6,a2)')
     .    'LAYER_ALLOYS(',itrack,')= ''',
     .     layer_alloys(itrack),''','
      enddo
      write(99,'(a11,i4,a1)') 'NDEPTH_NML=',ndepth_nml,','
      if (layer_thermal_contact.ne.1.0_R8) write(99,'(a22,1es15.8,a1)') 
     > 'LAYER_THERMAL_CONTACT=',layer_thermal_contact,','
      if (.not.(plate_model.ge.1.and.ndepth_nml.gt.1)) then
	do ibnd = 1, nbnd
        if (target_temp1(ibnd,1).gt.plate_temp) 
     >     write(99,'(a12,i4,a4,1x,1es15.8,a1)') 
     .     'TARGET_TEMP(',ibnd,',1)=',target_temp1(ibnd,1),','
        enddo
      endif
      do is = 0, ns1-1
        if (maxval(thermal_evaporation1(1:nbnd,is)).gt.0.0_R8) then
          do i = 1, (nbnd-1)/nline+1
            if (i.eq.1) then
              write(line,'(a22,i3,a2)') 'THERMAL_EVAPORATION(1,',is,')='
            else
              line = ''
            endif
            do ibnd = (i-1)*nline+1, min(i*nline,nbnd)
              write(numstr,'(1x,1es15.8,a1)') 
     .         thermal_evaporation1(ibnd,is),','
              line = trim(line)//numstr
            enddo
            write(99,'(a)') trim(line)
          enddo
        endif
        if (maxval(physical_sputtering1(1:nbnd,is)).gt.0.0_R8) then
          do i = 1, (nbnd-1)/nline+1
            if (i.eq.1) then
              write(line,'(a22,i3,a2)') 'PHYSICAL_SPUTTERING(1,',is,')='
            else
              line = ''
            endif
            do ibnd = (i-1)*nline+1, min(i*nline,nbnd)
              write(numstr,'(1x,1e15.8,a1)') 
     .         physical_sputtering1(ibnd,is),','
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
     .         physical_sputtering_energy1(ibnd,is),','
              line = trim(line)//numstr
            enddo
            write(99,'(a)') trim(line)
          enddo
        endif
        if (maxval(chemical_sputtering1(1:nbnd,is)).gt.0.0_R8) then
          do i = 1, (nbnd-1)/nline+1
            if (i.eq.1) then
              write(line,'(a22,i3,a2)') 
     .         'CHEMICAL_SPUTTERING(1,',is,')='
            else
              line = ''
            endif
            do ibnd = (i-1)*nline+1, min(i*nline,nbnd)
              write(numstr,'(1x,1es15.8,a1)') 
     .         chemical_sputtering1(ibnd,is),','
              line = trim(line)//numstr
            enddo
            write(99,'(a)') trim(line)
          enddo
        endif
        if (maxval(RES_sputtering1(1:nbnd,is)).gt.0.0_R8) then
          do i = 1, (nbnd-1)/nline+1
            if (i.eq.1) then
              write(line,'(a17,i3,a2)') 
     .         'RES_SPUTTERING(1,',is,')='
            else
              line = ''
            endif
            do ibnd = (i-1)*nline+1, min(i*nline,nbnd)
              write(numstr,'(1x,1es15.8,a1)') 
     .         RES_sputtering1(ibnd,is),','
              line = trim(line)//numstr
            enddo
            write(99,'(a)') trim(line)
          enddo
        endif
        if (maxval(backscattering1(1:nbnd,is)).gt.0.0_R8) then
          do i = 1, (nbnd-1)/nline+1
            if (i.eq.1) then
              write(line,'(a17,i3,a2)') 
     .         'BACKSCATTERING(1,',is,')='
            else
              line = ''
            endif
            do ibnd = (i-1)*nline+1, min(i*nline,nbnd)
              write(numstr,'(1x,1es15.8,a1)') 
     .         backscattering1(ibnd,is),','
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
     .         backscattering_energy1(ibnd,is),','
              line = trim(line)//numstr
            enddo
            write(99,'(a)') trim(line)
          enddo
        endif
      enddo
      if (plate_model.ge.1 .or. wall_namelist_used) then
cxpb  Use explicit writing of non-default elements only
        do ix = -1, nx1
          do iy = -1, ny1
            if (xymap1(ix,iy).ne.0) write(99,'(a6,i4,a1,i4,a2,i4,a1)')
     .       'XYMAP(',ix,',',iy,')=',xymap1(ix,iy),','     
          end do
        end do
        do ibnd = 1, nbnd
          write(99,'(a6,i4,a2,i4,a1,1x,a6,i4,a2,i4,a1)') 
     .     'IMAPX(',ibnd,')=',imapx1(ibnd),',',
     .     'IMAPY(',ibnd,')=',imapy1(ibnd),','
          if (inertial_cooling(ibnd))
     .     write(99,'(a17,i4,a2,l1,a1)')
     .     'INERTIAL_COOLING(',ibnd,')=',inertial_cooling1(ibnd),','
          if (backplate_temp1(ibnd).ne.plate_temp) 
     .     write(99,'(a15,i4,a2,1es15.8,a1)') 
     .     'BACKPLATE_TEMP(',ibnd,')=',backplate_temp1(ibnd),','
          if (plate_thickness1(ibnd).ne.plate_thick)
     .     write(99,'(a16,i4,a2,1es15.8,a1)') 
     .     'PLATE_THICKNESS(',ibnd,')=',plate_thickness1(ibnd),','
          if (coating_thickness1(ibnd).gt.0.0_R8)
     .     write(99,'(a18,i4,a2,1es15.8,a1)') 
     .     'COATING_THICKNESS(',ibnd,')=',coating_thickness1(ibnd),','
          if (plate_time_factor1(ibnd).ne.1.0_R8)
     .     write(99,'(a18,i4,a2,1es15.8,a1)') 
     .     'PLATE_TIME_FACTOR(',ibnd,')=',plate_time_factor1(ibnd),','
          if (trim(surface_material_name1(ibnd)).ne.'C') 
     .     write(99,'(a22,i4,a4,a6,a2)') 
     .     'SURFACE_MATERIAL_NAME(',ibnd,')= ''',
     .      surface_material_name1(ibnd),''','
          if (coating_thickness1(ibnd).gt.0.0_R8) 
     .     write(99,'(a22,i4,a4,a6,a2)') 
     .     'COATING_MATERIAL_NAME(',ibnd,')= ''',
     .      coating_material_name1(ibnd),''','
          if (trim(bulk_material_name1(ibnd)).ne.'C') 
     .     write(99,'(a19,i4,a4,a6,a2)') 
     .     'BULK_MATERIAL_NAME(',ibnd,')= ''',
     .      bulk_material_name1(ibnd),''','
          if (layer_model) then
            write(99,'(a20,i4,a2,i2,a1)')
     .       'LAYER_NCONSTITUENTS(',ibnd,')=',
     .        layer_nconstituents1(ibnd),','
            do i = 1, layer_nconstituents1(ibnd)
              write(99,'(a21,i4,a1,i2,a4,i2,a1)')
     .         'LAYER_NZCONSTITUENTS(',ibnd,',',i,')=  ',
     .          layer_nzconstituents1(ibnd,i),','
              write(99,'(a23,i4,a1,i2,a2,1es15.8,a1)')
     .         'LAYER_NRELCONSTITUENTS(',ibnd,',',i,')=',
     .          layer_nrelconstituents1(ibnd,i),','
            end do
          end if
          do itrack=1,ntrack
            if (deposition1(ibnd,itrack).ne.0.0_R8)
     .       write(99,'(a11,i4,a1,i1,a2,1es15.8,a1)')
     .        'DEPOSITION(',ibnd,',',itrack,')=',
     1        deposition1(ibnd,itrack),','
            if (monolayer_deposition1(ibnd,itrack).ne.0.0_R8)
     .       write(99,'(a,i4,a,i1,a,1es15.8,a)') 
     .        'MONOLAYER_DEPOSITION(',ibnd,',',itrack,')=',
     1         monolayer_deposition1(ibnd,itrack),','
            if (erosion1(ibnd,itrack).ne.0.0_R8)
     .       write(99,'(a8,i4,a1,i1,a2,1es15.8,a1)')
     .        'EROSION(',ibnd,',',itrack,')=',
     1         erosion1(ibnd,itrack),','
            if (monolayer_erosion1(ibnd,itrack).ne.0.0_R8)
     .       write(99,'(a,i4,a,i1,a,1es15.8,a)') 
     .        'MONOLAYER_EROSION(',ibnd,',',itrack,')=',
     1         monolayer_erosion1(ibnd,itrack),','
          enddo
          if (plate_time1(ibnd).ne.0.0_R8)
     .     write(99,'(a11,i4,a2,1es15.8,a1)')
     .      'PLATE_TIME(',ibnd,')=',plate_time1(ibnd),','
          if (plate_area1(ibnd).ne.0.0_R8)
     .     write(99,'(a11,i4,a2,1es15.8,a1)')
     .      'PLATE_AREA(',ibnd,')=',plate_area1(ibnd),','
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
     .           target_temp1(ibnd,is),','
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
        
      end subroutine b2yt_wall
