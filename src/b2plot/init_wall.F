      subroutine init_wall(nx,ny)
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
      use b2mod_b2cmfs
      use b2mod_external
      use b2mod_subsys
#ifndef ALLOCATE_AND_NAMELIST
#     include "DIMENSIONS.F"
#endif
      implicit none

      integer nx, ny, ns
      real (kind=R8) :: BoRiS

      integer ix, iy, is, ixn, iyn, idir, ireg
      integer ibnd, ibndm, i, j, ik, itrack, icoat, ikink, ikink0
      integer, allocatable, save :: imapx_tmp(:), imapy_tmp(:)
      integer, save :: debug_init_wall
      real (kind=R8) :: qin, qel, normal_r, normal_z, torl,
     *  xa, ya, xa2, ya2, xb, yb, xo, yo, nhalf_r, nhalf_z, nhalf,
     *  xap, yap, xaq, yaq, xbp, ybp, xbq, ybq, xbm, ybm, xint, yint,
     *  distance, distp, dint, sum,
     *  sput_frac(0:ntrack), wrk(ndepth_nml)
      character*16 hlp_frm
      logical wall_namelist_used, default_cooling, correct_quadrant
      logical error

      integer, save :: hd_linlog, output, sput_frac_flag,
     *  sput_phys_col, sput_phys_model, refl_model
      real (kind=R8), save :: hd_ratio, width
      data hd_linlog/1/, hd_ratio/100.0_R8/, output/0/
      data sput_frac_flag/0/, sput_phys_col/3/,
     *  sput_phys_model/1/, refl_model/1/, width/1.0_R8/
#ifdef DBG
      data debug_init_wall/1/
#else
      data debug_init_wall/0/
#endif

      real (kind=R8) :: mass_density
      logical overlap
      intrinsic mod, max, min, sqrt, abs, int
      external mass_density, overlap
      external xertst, ipgeti, ipgetr

      real (kind=R8) :: onemeps
      parameter (onemeps =  1.0_R8 - 1.0e-5_R8)

      call subini('init_wall')
      allocate(imapx_tmp(nwall), imapy_tmp(nwall))
      call ipgeti('init_wall_debug', debug_init_wall)
cxpb Establish numbering of boundary wall cells
cxpb We must avoid corner cells and dead cells
      nbnd = 0
      imapx_tmp= 0
      imapy_tmp= 0
      xymap_tmp= 0
      wall_begin = 0
      wall_end = 0
      nsegments = 0
      itrack_0 = 0

cxpb  We start from the lower left corner and go clockwise around the plasma
cxpb  keeping the plasma to our right (but we are moving in the guard cells!)
cxpb  Northern boundaries (excluding diverted island)
      iy = ny
      ix = -1
      ixn = rightix(ix,iy)
      iyn = rightiy(ix,iy)
      if (rightix(ixn,iyn).ne.nx+1) then
         idir = 2               ! 0 : left, 1: bottom, 2: right, 3: top
      else
! we have the special case that there are no cells on the top
         idir = 1
         ixn = bottomix(ixn,iyn)
         iyn = bottomiy(ixn,iyn)
         write(*,*) 'Top skipped'
      endif
      if(debug_init_wall.gt.0) then
        write(*,'(a,i3,1x,i3,1x,i1)')
     &   'init_wall(debug) start: ixn, iyn, idir ',
     &   ixn, iyn, idir
      endif
      nsegments = 1
      wall_begin(1) = 1
      do while (xymap_tmp(ixn,iyn).eq.0)
        nbnd = nbnd + 1
        imapx_tmp(nbnd) = ixn
        imapy_tmp(nbnd) = iyn
        xymap_tmp(ixn,iyn)= nbnd
        ix = ixn
        iy = iyn
        DIRVALUE: select case (idir)
          case(0)  ! left
            r1(nbnd)=crx(ix,iy,3)
            z1(nbnd)=cry(ix,iy,3)
            r2(nbnd)=crx(ix,iy,2)
            z2(nbnd)=cry(ix,iy,2)
            plate_area(nbnd) = gs(topix(ix,iy),topiy(ix,iy),1)
            ixn = leftix(ix,iy)
            iyn = leftiy(ix,iy)
            call check_ok(nx, ny, ixn, iyn, idir)
            if (leftix(ixn,iyn).eq.-2) then
              idir = 3
              ixn = leftix(topix(ix,iy),topiy(ix,iy))
              iyn = leftiy(topix(ix,iy),topiy(ix,iy))
              if(iyn.eq.ny) then
! special case for 1d linear device
                 write(*,*) 'Special treatment: 1d linear device'
                 iyn = iy+1
              endif
              wall_end(nsegments) = nbnd
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd+1
              if(debug_init_wall.gt.1) then
                write(*,*)
     &           'init_wall(debug): ',
     &           'switching from idir=0 (left) to idir=3 (top)'
              endif
            elseif (bottomiy(ixn,iyn).ne.-2) then
              idir = 1
              ixn = bottomix(leftix(ix,iy),leftiy(ix,iy))
              iyn = bottomiy(leftix(ix,iy),leftiy(ix,iy))
              wall_end(nsegments) = nbnd
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd+1
              if(debug_init_wall.gt.1) then
                write(*,*)
     &           'init_wall(debug): ',
     &           'switching from idir=0 (left) to idir=1 (bottom)'
              endif
            endif
          case(1)  ! bottom
            plate_area(nbnd) = gs(ix,iy,0)
            r1(nbnd)=crx(ix,iy,2)
            z1(nbnd)=cry(ix,iy,2)
            r2(nbnd)=crx(ix,iy,0)
            z2(nbnd)=cry(ix,iy,0)
            ixn = bottomix(ix,iy)
            iyn = bottomiy(ix,iy)
            call check_ok(nx, ny, ixn, iyn, idir)
            if (bottomiy(ixn,iyn).eq.-2) then
              idir = 0
              ixn = bottomix(leftix(ix,iy),leftiy(ix,iy))
              iyn = bottomiy(leftix(ix,iy),leftiy(ix,iy))
              wall_end(nsegments) = nbnd
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd+1
              if(debug_init_wall.gt.1) then
                write(*,*)
     &           'init_wall(debug): ',
     &           'switching from idir=1 (bottom) to idir=0 (left)'
              endif
            else if (leftix(ixn,iyn).eq.-2) then
              idir = 3
              wall_end(nsegments) = nbnd
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd+1
              if(debug_init_wall.gt.1) then
                write(*,*)
     &           'init_wall(debug): ',
     &           'switching from idir=1 (bottom) to idir=3 (top)'
              endif
            else if (rightix(ixn,iyn).ne.nx+1) then
              idir = 2
              ixn = rightix(bottomix(ix,iy),bottomiy(ix,iy))
              iyn = rightiy(bottomix(ix,iy),bottomiy(ix,iy))
              wall_end(nsegments) = nbnd
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd+1
              if(debug_init_wall.gt.1) then
                write(*,*)
     &           'init_wall(debug): ',
     &           'switching from idir=1 (bottom) to idir=2 (right)'
              endif
            elseif (iy.eq.iyn) then
! special case for 1d linear device
              write(*,*) 'Special treatment: 1d linear device'
              idir = 0
              ixn = bottomix(leftix(ix,iy),leftiy(ix,iy))
              iyn = -1
              wall_end(nsegments) = nbnd
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd+1
              if(debug_init_wall.gt.1) then
                write(*,*)
     &           'init_wall(debug): ',
     &           'switching from idir=1 (bottom) to idir=0 (left)'
              endif
            endif
          case(2)  ! right
            plate_area(nbnd) = gs(ix,iy,1)
            r1(nbnd)=crx(ix,iy,0)
            z1(nbnd)=cry(ix,iy,0)
            r2(nbnd)=crx(ix,iy,1)
            z2(nbnd)=cry(ix,iy,1)
            ixn = rightix(ix,iy)
            iyn = rightiy(ix,iy)
            call check_ok(nx, ny, ixn, iyn, idir)
            if (rightix(ixn,iyn).eq.nx+1) then
              idir = 1
              ixn = rightix(bottomix(ix,iy),bottomiy(ix,iy))
              iyn = rightiy(bottomix(ix,iy),bottomiy(ix,iy))
              if(iyn.eq.-1) then
! special case for 1d linear device
                 write(*,*) 'Special treatment: 1d linear device'
                 iyn = iy-1
              endif
              wall_end(nsegments) = nbnd
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd+1
              if(debug_init_wall.gt.1) then
                write(*,*)
     &           'init_wall(debug): ',
     &           'switching from idir=2 (right) to idir=1 (bottom)'
              endif
            else if (topiy(ixn,iyn).ne.ny+1) then
              idir = 3
              ixn = topix(rightix(ix,iy),rightiy(ix,iy))
              iyn = topiy(rightix(ix,iy),rightiy(ix,iy))
              wall_end(nsegments) = nbnd
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd+1
              if(debug_init_wall.gt.1) then
                write(*,*)
     &           'init_wall(debug): ',
     &           'switching from idir=2 (right) to idir=3 (top)'
              endif
            endif
          case(3)  ! top
            plate_area(nbnd) = gs(rightix(ix,iy),rightiy(ix,iy),0)
            r1(nbnd)=crx(ix,iy,1)
            z1(nbnd)=cry(ix,iy,1)
            r2(nbnd)=crx(ix,iy,3)
            z2(nbnd)=cry(ix,iy,3)
            ixn = topix(ix,iy)
            iyn = topiy(ix,iy)
            call check_ok(nx, ny, ixn, iyn, idir)
            if (topiy(ixn,iyn).eq.ny+1) then
              idir = 2
              ixn = topix(rightix(ix,iy),rightiy(ix,iy))
              iyn = topiy(rightix(ix,iy),rightiy(ix,iy))
              wall_end(nsegments) = nbnd
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd+1
              if(debug_init_wall.gt.1) then
                write(*,*)
     &           'init_wall(debug): ',
     &           'switching from idir=3 (top) to idir=2 (right)'
              endif
            else if (leftix(ixn,iyn).ne.-2) then
              idir = 0
              ixn = leftix(topix(ix,iy),topiy(ix,iy))
              iyn = leftiy(topix(ix,iy),topiy(ix,iy))
              wall_end(nsegments) = nbnd
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd+1
              if(debug_init_wall.gt.1) then
                write(*,*)
     &           'init_wall(debug): ',
     &           'switching from idir=3 (top) to idir=0 (left)'
              endif
            elseif (iy.eq.iyn) then
! special case for 1d linear device
              write(*,*) 'Special treatment: 1d linear device'
              idir = 2
              ixn = topix(rightix(ix,iy),rightiy(ix,iy))
              iyn = ny
              wall_end(nsegments) = nbnd
              nsegments = nsegments + 1
              wall_begin(nsegments) = nbnd+1
              if(debug_init_wall.gt.1) then
                write(*,*)
     &           'init_wall(debug): ',
     &           'switching from idir=3 (top) to idir=2 (right)'
              endif
            endif
          case default DIRVALUE
            call xerrab ('Invalid direction')
        end select DIRVALUE
        if(debug_init_wall.gt.0) then
          write(*,'(a,i3,1x,i3,1x,i1)')
     &     'init_wall(debug) continue: ixn, iyn, idir ',
     *     ixn, iyn, idir
        endif
      enddo
      nsegments = nsegments - 1

cxpb  We now identify the b2ah region numbers
      do ibnd = 1, nbnd
        ix = imapx_tmp(ibnd)
        iy = imapy_tmp(ibnd)
        if(iy.eq.ny) then
          ireg = cbirno
          do while (cbrbrk(ireg)*nx.lt.ix+0.5_R8 .and.
     &                     ireg.le.cbirno+cbnrno-1)
            ireg = ireg + 1
          enddo
          imapreg(ibnd) = ireg
        elseif (ix.eq.nx) then
          ireg = cbirea
          do while (cbrbrk(ireg)*ny.lt.iy+0.5_R8 .and.
     &                     ireg.le.cbirea+cbnrea-1)
            ireg = ireg + 1
          enddo
          imapreg(ibnd) = ireg
        elseif (iy.eq.-1) then
          ireg = cbirso
          do while (cbrbrk(ireg)*nx.lt.ix+0.5_R8 .and.
     &                     ireg.le.cbirso+cbnrso-1)
            ireg = ireg + 1
          enddo
          imapreg(ibnd) = ireg
        elseif (ix.eq.-1) then
          ireg = cbirwe
          do while (cbrbrk(ireg)*ny.lt.iy+0.5_R8 .and.
     &                     ireg.le.cbirwe+cbnrwe-1)
            ireg = ireg + 1
          enddo
          imapreg(ibnd) = ireg
        else
          imapreg(ibnd) = -1
        end if
      end do

      do is = nsegments+1, nwall_segments
        wall_begin(is) = nbnd+1
      enddo
#ifndef ALLOCATE_AND_NAMELIST
      call xertst (nbnd.le.2*(DEF_NXD+2)+(DEF_NCUT+2)*(DEF_NYD+2),
     & 'Increase size of B2 contribution to nwall in b2mod_wall')
#endif
cxpb Recopy geometrical arrays into namelist defaults
      imapx(1:nbnd) = imapx_tmp(1:nbnd)
      imapy(1:nbnd) = imapy_tmp(1:nbnd)
      xymap(-1:nx,-1:ny) = xymap_tmp(-1:nx,-1:ny)
      call subend()
      return

      entry get_wall_data(nx,ny,ns,BoRiS,wall_namelist_used)
      call subini('get_wall_data')
cxpb  Test namelist input
      if (ntargsp.gt.0) then
        write(hlp_frm,'(a,i2,a)') '(1x,a,',ntargsp,'(a2,1x))'
        write(*,hlp_frm) 'track_species : Bulk ', 
     .   track_species(1:ntargsp)
      else
        write(*,'(1x,a)') 'No tracked wall material species'
      end if
      if (wall_namelist_used) then
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
     .        imapx(xymap(ix,iy)).eq.ix .and.
     .        imapy(xymap(ix,iy)).eq.iy,
     .         'Inconsistent boundary mapping in namelist WALL!')
            endif
          end do
        end do
      end if
cxpb Fill up array data not given in namelist
      layer_composition = 0.0_R8
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
        layer_thickness=0.0_R8
        layer_depth=0.0_R8
        call subend()
        return
      endif
      call read_b2mod_wall(ntargsp)
      call ipgeti ('b2stbr_sput_phys_col', sput_phys_col)
c the TRIM sputtering files have 9 columns corresponding to
c 0 15 30 45 55 65 75 80 85 degree impact
      call xertst (0.le.sput_phys_col.and.sput_phys_col.le.9,
     1 'faulty b2stbr_sput_phys_col')
      call ipgeti ('b2stbr_sput_phys_model', sput_phys_model)
      call ipgeti ('b2stbr_refl_model', refl_model)
      if (sput_dst.ge.0) call read_sputter_data(ns,
     1 sput_phys_col,sput_phys_model,refl_model)

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
        layer_relconstituents(ibnd,1:layer_nconstituents(ibnd))=
     1   layer_nrelconstituents(ibnd,1:layer_nconstituents(ibnd))/sum
      end do

cxpb Compute target geometry data
      call ipgeti('b2stbr_output', output)
      call ipgetr('b2agmt_1d_width', width)
      call xertst(0.0_R8.lt.width, 'faulty input width')
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
cxpb  Compute initial layer thickness (deposited - eroded)
cxpb  and layer depth (sum of deposition),
cxpb  i.e. minimum depth over which material mixing occurs
        layer_thickness(ibnd)=0.0_R8
        layer_depth(ibnd)=0.0_R8
        if (coating_thickness(ibnd).gt.0.0_R8) then
cxpb  Figure out where to put the coating-bulk boundary
          if (hd_linlog.eq.1) then
            icoat = int(ndepth_nml*coating_thickness(ibnd)/
     .                       plate_thickness(ibnd) + onemeps)
          else
            icoat = int( (ndepth_nml-1)*
     .       log(coating_thickness(ibnd)*(hd_ratio-1.0_R8)/
     .           plate_thickness(ibnd) + 1.0_R8) / log(hd_ratio) +
     .           onemeps )+1
          end if
        else
          icoat = 1
        end if
        do i=1,icoat-1
          if (hd_linlog.eq.1) then
            target_depth(ibnd,i) = coating_thickness(ibnd) *
     .        float(i-1) / float(icoat - 1)
          else
            target_depth(ibnd,i) = coating_thickness(ibnd) *
     .       (hd_ratio ** (float(i-1) / float(icoat - 1))
     .                                                - 1.0_R8)/
     .       (hd_ratio - 1.0_R8)
          endif
        enddo
        target_depth(ibnd,icoat)=coating_thickness(ibnd)
        do i=icoat+1,ndepth_nml-1
          if (hd_linlog.eq.1) then
            target_depth(ibnd,i) = coating_thickness(ibnd) +
     .       (plate_thickness(ibnd)-coating_thickness(ibnd)) *
     .        float(i-1) / float(ndepth_nml - 1)
          else
            target_depth(ibnd,i) = coating_thickness(ibnd) +
     .       (plate_thickness(ibnd)-coating_thickness(ibnd)) *
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
          call add_layer_thickness(ns,ibnd)
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
        endif
        call compose_layer(ibnd,sput_frac)
        if (ndepth_nml.gt.1) then
          target_depth(ibnd,ndepth_nml) = plate_thickness(ibnd)
        else
          target_depth(ibnd,ndepth_nml) = 0.0_R8
        endif
      enddo

      error = .false.
      do is = 1, nsegments
        do ibnd = wall_begin(is), wall_end(is)
          if (ibnd.eq.wall_begin(is)) then
            ibndm = ibnd
            ikink = ndepth_nml+1
            ikink0 = ndepth_nml+1
          else
            ibndm = ibnd-1
          endif
          xa = r1(ibnd)
          xb = r2(ibnd)
          ya = z1(ibnd)
          yb = z2(ibnd)
          distance = sqrt((xb-xa)**2+(yb-ya)**2)
          target_width(ibnd,1) = distance
          targetcr(ibnd,1,1) = xa
          targetcz(ibnd,1,1) = ya
          targetcr(ibnd,1,2) = xb
          targetcz(ibnd,1,2) = yb
          if (ibnd.eq.wall_begin(is)) then
            normal_r = -(yb-ya)/target_width(ibnd,1)
            normal_z =  (xb-xa)/target_width(ibnd,1)
          endif
          nhalf_r = -(z2(ibnd)-z1(ibndm))
          nhalf_z =  (r2(ibnd)-r1(ibndm))
          nhalf = sqrt(nhalf_r**2 + nhalf_z**2)
          nhalf_r = nhalf_r/nhalf
          nhalf_z = nhalf_z/nhalf
          if (ibnd.eq.wall_begin(is)) then
            xaq = xa + normal_r*target_depth(ibnd,ndepth_nml)
            yaq = ya + normal_z*target_depth(ibnd,ndepth_nml)
          else
            xaq = targetcr(ibnd-1,min(ikink0,ndepth_nml),2)
            yaq = targetcz(ibnd-1,min(ikink0,ndepth_nml),2)
          endif
          xbq = xb + nhalf_r*target_depth(ibnd,ndepth_nml)
          ybq = yb + nhalf_z*target_depth(ibnd,ndepth_nml)
          if (abs(normal_r*nhalf_z-normal_z*nhalf_r).le.1.e-10_R8) then  ! parallel case
            xint = 0.0_R8
            yint = 0.0_R8
            dint = 1.0e30_R8
          else ! We are not interested with the actual intersection, but rather with the intersection
               ! to a line 0.1*AB closer
            xa2 = xa + 0.1_R8*(xb-xa)
            ya2 = ya + 0.1_R8*(yb-ya)
            if (nhalf_z.eq.0.0_R8) then  ! horizontal case
              yint = yb
              xint = xa2 + (normal_r/normal_z)*(yb-ya2)
            else
              yint = ((xb-xa2)*normal_z*nhalf_z+ya2*normal_r*nhalf_z
     .                                         -yb*nhalf_r*normal_z)/
     .               (normal_r*nhalf_z-normal_z*nhalf_r)
              xint = xb - (yb-yint)*nhalf_r/nhalf_z
            endif
            dint = sqrt((xint-xb)**2+(yint-yb)**2)
            if (((xint-xb)*(xbq-xb)+(yint-yb)*(ybq-yb)).lt.0.0_R8)
     >       dint = -dint
          end if
cxpb This does not properly take into account the case where the kink is receding away from the surface !!!
          if (dint.gt.0.0_R8 .and.
     &        dint.lt.target_depth(ibnd,min(ikink0,ndepth_nml))) then
            ikink = min(ikink0,ndepth_nml)
            do while (dint.lt.target_depth(ibnd,ikink))
              ikink = ikink - 1
            enddo
          else if (dint.lt.0.0_R8 .or.
     &             dint.gt.target_depth(ibnd,ndepth_nml)) then
            ikink = ndepth_nml+1
          else
            ikink = ikink0
            do while (dint.gt.target_depth(ibnd,ikink))
              ikink = ikink + 1
            enddo
          endif
          ik = ikink
          do i = 2, ndepth_nml
            if (ibnd.eq.wall_begin(is)) then
              xap = xa + normal_r*target_depth(ibnd,i)
              yap = ya + normal_z*target_depth(ibnd,i)
              xo = xa + normal_r*target_depth(ibnd,i-1)
              yo = ya + normal_z*target_depth(ibnd,i-1)
            else
              xap = targetcr(ibnd-1,i,2)
              yap = targetcz(ibnd-1,i,2)
              xo = targetcr(ibnd-1,i-1,2)
              yo = targetcz(ibnd-1,i-1,2)
            endif
            targetcr(ibnd,i,1) = xap
            targetcz(ibnd,i,1) = yap
cxpb Proceed carefully, checking at every step that we are not stepping backwards
cxpb First attempt: we assume the normal case
            if (i.le.ik) then
              xbp = xb + nhalf_r*target_depth(ibnd,i)
              ybp = yb + nhalf_z*target_depth(ibnd,i)
            else
              xbp = xap + 0.1_R8*(xb-xa)
              ybp = yap + 0.1_R8*(yb-ya)
            endif
            if (ibnd.ne.wall_begin(is)) then
cxpb We check we remain in the correct quadrant
              xbm = targetcr(ibnd,i-1,2)
              ybm = targetcz(ibnd,i-1,2)
cxpb We check: the segment OBp must be between OBm and OAp
cxpb and far enough away from the OAp line
              correct_quadrant =
     .         ((xbm-xo)*(ybp-yo)-(ybm-yo)*(xbp-xo))*
     .         ((xbp-xo)*(yap-yo)-(ybp-yo)*(xap-xo)).gt.0.0_R8
              distp = ((xbp-xo)*(yap-yo)+(ybp-yo)*(xo-xap))/
     .         sqrt((yap-yo)**2+(xo-xap)**2)
              if(.not.(correct_quadrant.and.distp/distance.ge.0.1_R8))
     .         then
                ik = min(ik,i)
                xbp = xo + (xbm-xo) + (xap-xo)
                ybp = yo + (ybm-yo) + (yap-yo)
              endif
cxpb We still need to check if we have not overstepped into a deeper target section
              do j = i+1, ndepth_nml
                if (overlap(xo,yo,xap,yap,xbp,ybp,xbm,ybm,
     .               targetcr(ibnd-1,j-1,1),targetcz(ibnd-1,j-1,1),
     .               targetcr(ibnd-1,j-1,2),targetcz(ibnd-1,j-1,2),
     .               targetcr(ibnd-1,j,2),targetcz(ibnd-1,j,2),
     .               targetcr(ibnd-1,j,1),targetcz(ibnd-1,j,1))) then
cxpb Do not know yet what to do here. Just a diagnostic in the meantime
                 error = .true.
                 write(6,'(a,i3,a1,i3,a,i3,a1,i3,a2)')
     .            'Target element (',ibnd,',',i,') overlapping with (',
     .                               ibnd-1,',',j,')!'
!                  xbp = xbp + 0.5_R8*
!     .             (targetcr(ibnd-1,j  ,2)-targetcr(ibnd-1,j  ,1)+
!     .              targetcr(ibnd-1,j-1,2)-targetcr(ibnd-1,j-1,1))
!                  ybp = ybp + 0.5_R8*
!     .             (targetcz(ibnd-1,j  ,2)-targetcz(ibnd-1,j  ,1)+
!     .              targetcz(ibnd-1,j-1,2)-targetcz(ibnd-1,j-1,1))
                endif
              enddo
            end if
            targetcr(ibnd,i,2) = xbp
            targetcz(ibnd,i,2) = ybp
            target_width(ibnd,i) = sqrt((xbp-xap)**2+(ybp-yap)**2)
          enddo
          ikink0 = ik
          normal_r = nhalf_r
          normal_z = nhalf_z
        end do
      end do
      if (error) then
        write(6,*) 'Found target elements overlap:'
        write(6,*)
     .   'Consider using thinner targets or straigthening them'
      endif
      if (ndepth_nml.gt.1) then
        do ibnd = 1, nbnd
          if (isymm.eq.0) then
            torl = width
          else if (isymm.eq.1 .or. isymm.eq.2) then
            torl = pi*(targetcr(ibnd,1,1)+targetcr(ibnd,1,2))
          else if (isymm.eq.3 .or. isymm.eq.4) then
            torl = pi*(targetcz(ibnd,1,1)+targetcz(ibnd,1,2))
          endif
          target_vol(ibnd,1) = target_width(ibnd,1)*torl*
     1     ((target_depth(ibnd,2)-target_depth(ibnd,1))/2.0_R8)
          do is = 2, ndepth_nml-1
            if (isymm.eq.0) then
              torl = width
            else if (isymm.eq.1 .or. isymm.eq.2) then
              torl = pi*(targetcr(ibnd,is,1)+targetcr(ibnd,is,2))
            else if (isymm.eq.3 .or. isymm.eq.4) then
              torl = pi*(targetcz(ibnd,is,1)+targetcz(ibnd,is,2))
            endif
            target_vol(ibnd,is) = target_width(ibnd,is)*torl*
     1       (target_depth(ibnd,is)-target_depth(ibnd,is-1))
          enddo
          if (isymm.eq.0) then
            torl = width
          else if (isymm.eq.1 .or. isymm.eq.2) then
            torl = pi*(targetcr(ibnd,ndepth_nml,1)+
     .               targetcr(ibnd,ndepth_nml,2))
          else if (isymm.eq.3 .or. isymm.eq.4) then
            torl = pi*(targetcz(ibnd,ndepth_nml,1)+
     .               targetcz(ibnd,ndepth_nml,2))
          endif
          target_vol(ibnd,ndepth_nml) = target_width(ibnd,ndepth_nml)*
     1     ((target_depth(ibnd,ndepth_nml)-
     .       target_depth(ibnd,ndepth_nml-1))
     2      /2.0_R8)*torl
        enddo
      else
        target_vol(1:nbnd,1) = target_width(1:nbnd,1)
      endif

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
     &                topiy(imapx(ibnd),imapy(ibnd)),1,1))/
     &            gs(topix(imapx(ibnd),imapy(ibnd)),
     &               topiy(imapx(ibnd),imapy(ibnd)),1)
              qin = qin + max(0.0_R8,
     &           -fhi_ext(topix(imapx(ibnd),imapy(ibnd)),
     &                topiy(imapx(ibnd),imapy(ibnd)),1,1))/
     &            gs(topix(imapx(ibnd),imapy(ibnd)),
     &               topiy(imapx(ibnd),imapy(ibnd)),1)
              do is=0,ns-1
                qin=qin+max(0.0_R8,
     &             -fna(topix(imapx(ibnd),imapy(ibnd)),
     &                  topiy(imapx(ibnd),imapy(ibnd)),1,1,is))
     &             /gs(topix(imapx(ibnd),imapy(ibnd)),
     &                 topiy(imapx(ibnd),imapy(ibnd)),1)*
     &             ((kinrgy(imapx(ibnd),imapy(ibnd),is)+
     &                   ti(imapx(ibnd),imapy(ibnd)))*(1.0_R8-BoRiS) +
     &              rpt(imapx(ibnd),imapy(ibnd),is)*ev)
              enddo
              do is=0,ns_ext-1
                qin=qin+max(0.0_R8,
     &             -fa_ext(topix(imapx(ibnd),imapy(ibnd)),
     &                     topiy(imapx(ibnd),imapy(ibnd)),1,1,is))
     &             /gs(topix(imapx(ibnd),imapy(ibnd)),
     &                 topiy(imapx(ibnd),imapy(ibnd)),1)*
     &             ((0.5_R8*am_ext(is)*mp*
     &               ua_ext(imapx(ibnd),imapy(ibnd),is)**2+
     &               ta_ext(imapx(ibnd),imapy(ibnd),is))*(1.0_R8-BoRiS)
     &              +pt_ext(imapx(ibnd),imapy(ibnd),is)*ev)
              enddo
              qel = max(0.0_R8,
     &           -(fhe(topix(imapx(ibnd),imapy(ibnd)),
     &                topiy(imapx(ibnd),imapy(ibnd)),1,1)+
     &            (1.0_R8-BoRiS)*fne(topix(imapx(ibnd),imapy(ibnd)),
     &                topiy(imapx(ibnd),imapy(ibnd)),1,1)*
     &                te(imapx(ibnd),imapy(ibnd))))/
     &            gs(topix(imapx(ibnd),imapy(ibnd)),
     &               topiy(imapx(ibnd),imapy(ibnd)),1)
            else if (topiy(imapx(ibnd),imapy(ibnd)).eq.ny+1) then
cxpb          North boundary
              qin = max(0.0_R8,fhi(imapx(ibnd),imapy(ibnd),1,1))
     &           /gs(imapx(ibnd),imapy(ibnd),1)
              qin = qin + max(0.0_R8,
     &           fhi_ext(imapx(ibnd),imapy(ibnd),1,1))
     &           /gs(imapx(ibnd),imapy(ibnd),1)
              do is=0,ns-1
                qin=qin+max(0.0_R8,fna(imapx(ibnd),imapy(ibnd),1,1,is))
     &             /gs(imapx(ibnd),imapy(ibnd),1)*
     &            ((kinrgy(imapx(ibnd),imapy(ibnd),is)+
     &              ti(imapx(ibnd),imapy(ibnd)))*(1.0_R8-BoRiS) +
     &              rpt(imapx(ibnd),imapy(ibnd),is)*ev)
              enddo
              do is=0,ns_ext-1
                qin=qin+max(0.0_R8,
     &             fa_ext(imapx(ibnd),imapy(ibnd),1,1,is))
     &             /gs(imapx(ibnd),imapy(ibnd),1)*
     &             ((0.5_R8*am_ext(is)*mp*
     &              ua_ext(imapx(ibnd),imapy(ibnd),is)**2+
     &              ta_ext(imapx(ibnd),imapy(ibnd),is))*(1.0_R8-BoRiS)
     &             +pt_ext(imapx(ibnd),imapy(ibnd),is)*ev)
              enddo
              qel = max(0.0_R8,fhe(imapx(ibnd),imapy(ibnd),1,1)+
     &          (1.0_R8-BoRiS)*fne(imapx(ibnd),imapy(ibnd),1,1)*
     &                          te(imapx(ibnd),imapy(ibnd)))
     &                         /gs(imapx(ibnd),imapy(ibnd),1)
            else if (leftix(imapx(ibnd),imapy(ibnd)).eq.-2) then
cxpb          West boundary
              qin = max(0.0_R8,-fhi(rightix(imapx(ibnd),imapy(ibnd)),
     &                          rightiy(imapx(ibnd),imapy(ibnd)),0,0))
     &           /gs(rightix(imapx(ibnd),imapy(ibnd)),
     &               rightiy(imapx(ibnd),imapy(ibnd)),0)
              qin = qin + max(0.0_R8,
     &                 -fhi_ext(rightix(imapx(ibnd),imapy(ibnd)),
     &                          rightiy(imapx(ibnd),imapy(ibnd)),0,0))
     &           /gs(rightix(imapx(ibnd),imapy(ibnd)),
     &               rightiy(imapx(ibnd),imapy(ibnd)),0)
              do is=0,ns-1
                qin=qin+max(0.0_R8,
     &             -fna(rightix(imapx(ibnd),imapy(ibnd)),
     &                  rightiy(imapx(ibnd),imapy(ibnd)),0,0,is))
     &              /gs(rightix(imapx(ibnd),imapy(ibnd)),
     &                  rightiy(imapx(ibnd),imapy(ibnd)),0)*
     &             ((kinrgy(imapx(ibnd),imapy(ibnd),is)+
     &                   ti(imapx(ibnd),imapy(ibnd)))*(1.0_R8-BoRiS) +
     &              rpt(imapx(ibnd),imapy(ibnd),is)*ev)
              enddo
              do is=0,ns_ext-1
                qin=qin+max(0.0_R8,
     &             -fa_ext(rightix(imapx(ibnd),imapy(ibnd)),
     &                     rightiy(imapx(ibnd),imapy(ibnd)),0,0,is))
     &              /gs(rightix(imapx(ibnd),imapy(ibnd)),
     &                  rightiy(imapx(ibnd),imapy(ibnd)),0)*
     &             ((0.5_R8*am_ext(is)*mp*
     &              ua_ext(imapx(ibnd),imapy(ibnd),is)**2+
     &              ta_ext(imapx(ibnd),imapy(ibnd),is))*(1.0_R8-BoRiS)
     &             +pt_ext(imapx(ibnd),imapy(ibnd),is)*ev)
              enddo
              qel = max(0.0_R8,
     &           -(fhe(rightix(imapx(ibnd),imapy(ibnd)),
     &                 rightiy(imapx(ibnd),imapy(ibnd)),0,0)+
     &             fne(rightix(imapx(ibnd),imapy(ibnd)),
     &                 rightiy(imapx(ibnd),imapy(ibnd)),0,0)*
     &              te(rightix(imapx(ibnd),imapy(ibnd)),
     &                 rightiy(imapx(ibnd),imapy(ibnd)))*
     &             (1.0_R8-BoRiS)))
     &           /gs(rightix(imapx(ibnd),imapy(ibnd)),
     &               rightiy(imapx(ibnd),imapy(ibnd)),0)
            else if (rightix(imapx(ibnd),imapy(ibnd)).eq.nx+1) then
cxpb          East boundary
              qin = max(0.0_R8,fhi(imapx(ibnd),imapy(ibnd),0,0))
     &           /gs(imapx(ibnd),imapy(ibnd),0)
              qin = qin + max(0.0_R8,
     &           fhi_ext(imapx(ibnd),imapy(ibnd),0,0))
     &           /gs(imapx(ibnd),imapy(ibnd),0)
              do is=0,ns-1
                qin = qin +
     &             max(0.0_R8,fna(imapx(ibnd),imapy(ibnd),0,0,is))
     &             /gs(imapx(ibnd),imapy(ibnd),0)*
     &            ((kinrgy(imapx(ibnd),imapy(ibnd),is)+
     &              ti(imapx(ibnd),imapy(ibnd)))*(1.0_R8-BoRiS) +
     &              rpt(imapx(ibnd),imapy(ibnd),is)*ev)
              enddo
              do is=0,ns_ext-1
                qin = qin +
     &             max(0.0_R8,fa_ext(imapx(ibnd),imapy(ibnd),0,0,is))
     &             /gs(imapx(ibnd),imapy(ibnd),0)*
     &            ((0.5_R8*am_ext(is)*mp*
     &              ua_ext(imapx(ibnd),imapy(ibnd),is)**2+
     &              ta_ext(imapx(ibnd),imapy(ibnd),is))*(1.0_R8-BoRiS)
     &             +pt_ext(imapx(ibnd),imapy(ibnd),is)*ev)
              enddo
              qel = max(0.0_R8,fhe(imapx(ibnd),imapy(ibnd),0,0)+
     &                         fne(imapx(ibnd),imapy(ibnd),0,0)*
     &                         te(imapx(ibnd),imapy(ibnd))*
     &          (1.0_R8-BoRiS))/gs(imapx(ibnd),imapy(ibnd),0)
            endif
cxpb Use 1-D equilibrium profile
            wrk(1:ndepth_nml) = target_temp(ibnd,1:ndepth_nml)
            call heatdiff1D(ibnd,ns,qin+qel,0.0_R8,wrk,output)
            target_temp(ibnd,1:ndepth_nml) = wrk(1:ndepth_nml)
          endif   ! plate_option
        enddo   !  ibnd
      endif   !  ndepth_nml
      call subend()
      return

* ..finalise time-dependent wall data
      entry wall_final(nx,ny,ns,wall_namelist_used)
      wall_namelist_used = wall_namelist_used .or. ntargsp.gt.0 .or.
     .  layer_thermal_contact.ne.1.0_R8 .or.
     .  maxval(target_temp(1:nbnd,1)).gt.plate_temp .or.
     .  maxval(thermal_evaporation(1:nbnd,0:ns-1)).gt.0.0_R8 .or.
     .  maxval(physical_sputtering(1:nbnd,0:ns-1)).gt.0.0_R8 .or.
     .  maxval(chemical_sputtering(1:nbnd,0:ns-1)).gt.0.0_R8 .or.
     .  maxval(RES_sputtering(1:nbnd,0:ns-1)).gt.0.0_R8 .or.
     .  maxval(backscattering(1:nbnd,0:ns-1)).gt.0.0_R8
      default_cooling = .true.
      do ibnd = 1, nbnd
        default_cooling = default_cooling .and.
     &   .not.inertial_cooling(ibnd)
      enddo
      if (.not.wall_namelist_used) return
      call write_b2mod_wall_namelist(99,nx,ny,ns,'    ')
      return

      contains

      subroutine check_ok(nx, ny, ix, iy, idir)
      implicit none
      integer, intent(in) :: nx, ny, ix, iy, idir
! check on position
      if(ix.eq.-2) then
         write(*,*) 'Stepped off the edge of the grid, ix, iy, idir = ',
     &        ix, iy, idir
         call xerrab('Invalid ix')
      endif
      if(ix.eq.nx+1) then
         write(*,*) 'Stepped off the edge of the grid, ix, iy, idir = ',
     &        ix, iy, idir
         call xerrab('Invalid ix')
      endif
      if(iy.eq.-2) then
         write(*,*) 'Stepped off the edge of the grid, ix, iy, idir = ',
     &        ix, iy, idir
         call xerrab('Invalid iy')
      endif
      if(iy.eq.ny+1) then
         write(*,*) 'Stepped off the edge of the grid, ix, iy, idir = ',
     &        ix, iy, idir
         call xerrab('Invalid iy')
      endif
      end subroutine check_ok

      end subroutine init_wall

      subroutine add_layer_thickness(ns,ibnd)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      use b2mod_b2cmpa
      use b2mod_sputter
      use b2mod_indirect
      implicit none
      integer ns, nz, ibnd
      integer i, is, itrack, icoat
      real (kind=R8) :: sum, sput_frac(0:ntrack)
      logical found
      character*24 filename
      integer, save :: ncall, hd_linlog, sput_frac_flag
      real (kind=R8), save :: hd_ratio,
     .  sput_mixed_alpha, sput_mixed_beta
      data ncall/0/, hd_linlog/1/, hd_ratio/100.0_R8/,
     .     sput_mixed_alpha/1.0_R8/, sput_mixed_beta/1.0_R8/,
     2     sput_frac_flag/0/

      real (kind=R8) :: mass_density
      real (kind=R8) :: onemeps
      parameter (onemeps = 1.0_R8 - 1.0e-5_R8)
      logical open_file
      intrinsic max
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
            if (coating_thickness(ibnd).gt.0.0_R8) then
              icoat = int(ndepth_nml*coating_thickness(ibnd)/
     .                   (plate_thickness(ibnd)-layer_thickness(ibnd)) +
     .                    onemeps)+2
            else
              icoat = 2
            endif
            do i = 2, icoat-1
              target_depth(ibnd,i+1) = target_depth(ibnd,2) +
     .         (coating_thickness(ibnd)-target_depth(ibnd,2)) *
     .         float(i-1) / float(icoat-1)
            enddo
            target_depth(ibnd,icoat)=coating_thickness(ibnd)
            do i = icoat+1, ndepth_nml-1
              target_depth(ibnd,i) = target_depth(ibnd,icoat) +
     .         (plate_thickness(ibnd)-target_depth(ibnd,icoat)) *
     .         float(i-icoat) / float(ndepth_nml - icoat)
            end do
          else
            target_depth(ibnd,2) = max(target_depth(ibnd,2),
     .       (plate_thickness(ibnd)-layer_thickness(ibnd)) *
     .         (hd_ratio ** (1.0_R8 / float(ndepth_nml - 1)) - 1.0_R8)/
     .         (hd_ratio - 1.0_R8))
            if (coating_thickness(ibnd).gt.0.0_R8) then
              icoat = int( (ndepth_nml-1)*
     .         log(coating_thickness(ibnd)*(hd_ratio-1.0_R8)/
     .            (plate_thickness(ibnd)-layer_thickness(ibnd))+1.0_R8)
     .            / log(hd_ratio) + onemeps )+2
            else
              icoat = 2
            endif
            do i = 2, icoat-1
              target_depth(ibnd,i) = target_depth(ibnd,2) +
     .         (coating_thickness(ibnd)-target_depth(ibnd,2)) *
     .         (hd_ratio ** (float(i-1) / float(icoat - 1))
     .                                                      - 1.0_R8)/
     .         (hd_ratio - 1.0_R8)
            end do
            target_depth(ibnd,icoat)=coating_thickness(ibnd)
            do i = icoat+1, ndepth_nml-1
              target_depth(ibnd,i) = target_depth(ibnd,icoat) +
     .         (plate_thickness(ibnd)-target_depth(ibnd,icoat)) *
     .         (hd_ratio ** (float(i-icoat) / float(ndepth_nml - icoat))
     .                                                      - 1.0_R8)/
     .         (hd_ratio - 1.0_R8)
            end do
          endif
        else
cxpb  Very deep erosion case, crosses coating boundary
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
            do i = 3, ndepth_nml-1
              target_depth(ibnd,i) = target_depth(ibnd,2) +
     .         (plate_thickness(ibnd)-target_depth(ibnd,2)) *
     .         (hd_ratio ** (float(i-2) / float(ndepth_nml - 2))
     .                                                       - 1.0_R8)/
     .         (hd_ratio - 1.0_R8)
            end do
          endif
          write(*,*)
     .     'Fully eroded coating. Read new bulk surface properties'
          filename='Surface_properties/'//trim(bulk_material_name(ibnd))
          if(open_file(99,filename,
     .       'data.local|modules/B2.5/Database')) then
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
cxpb Normalise surface constituents
          sum = 0.0_R8
          do i = 1, matsurf_nconstituents(ibnd)
            sum = sum + matsurf_nrelconstituents(ibnd,i)
          end do
          do i = 1, matsurf_nconstituents(ibnd)
            matsurf_relconstituents(ibnd,i) =
     .       matsurf_nrelconstituents(ibnd,i)/sum
          end do
        endif
      endif

cxpb We now identify the track index of the bulk or coating species
cxpb below the re-deposited layer.
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

cxpb  Code similar to b2stbr_phys.recycle in order to initialise the
cxpb  layer quantities
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
cxpb  Modify surface layer properties for backscattering
      layer_nconstituents(ibnd)=matsurf_nconstituents(ibnd)+ntargsp
cxpb  Contributions from the bulk
      layer_relconstituents(ibnd,1:matsurf_nconstituents(ibnd))=
     1  matsurf_relconstituents
     .   (ibnd,1:matsurf_nconstituents(ibnd))*
     2  sput_frac(0)
cxpb  Contributions from the deposited species
      do itrack=1,ntargsp
        layer_nzconstituents
     .     (ibnd,matsurf_nconstituents(ibnd)+itrack)=
     1     nint(zn(target_intcode(itrack)))
        layer_relconstituents
     .     (ibnd,matsurf_nconstituents(ibnd)+itrack)=
     .      sput_frac(itrack)
      enddo
      call compose_layer(ibnd,sput_frac)
cxpb  Update layer mass density
cxpb  with bulk and deposited contributions
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
cxpb  Compute layer composition
cxpb  This will need refinement including alloy formation
      nz = matsurf_nconstituents(ibnd)
cxpb  Contribution from the bulk
      layer_composition(ibnd,1:nz)=
     2  matsurf_relconstituents(ibnd,1:nz)*sput_frac(0)
      if (layer_nconstituents(ibnd).eq.nz) return
cxpb  Inventory of bulk or coating particles
      if (coating_thickness(ibnd).gt.0.0_R8 .and.
     &    target_depth(ibnd,1).lt.coating_thickness(ibnd)) then
cxpb  Coating case
        if (itrack_0(ibnd).ne.0) then
          particle_inventory(ibnd,itrack_0(ibnd))=
     1     (coating_thickness(ibnd)-layer_depth(ibnd))*
     2      plate_area(ibnd)*mass_density(ibnd,2)/
     3     (atomic_mass(matsurf_nzconstituents(ibnd,1))*mp)
        endif
      else
cxpb  Deep erosion or no coating case
        if (itrack_0(ibnd).ne.0) then
          particle_inventory(ibnd,itrack_0(ibnd))=
     1     (plate_thickness(ibnd)-layer_depth(ibnd))*
     2      plate_area(ibnd)*mass_density(ibnd,2)/
     3     (atomic_mass(matsurf_nzconstituents(ibnd,1))*mp)
        endif
      endif
cxpb  Additional contribution from deposited material
cxpb  Must check if Z of deposited material already exists in bulk
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
      real (kind=R8) :: T, layer_Kappa
      integer nz, itrack
      real (kind=R8) :: Tcels, Teff, a(5)
      real (kind=R8) :: bulk_value, coating_value, layer_value, frac
      external xerrab

      Tcels = T - 273.15_R8

cxpb First count the bulk contribution
      frac=layer_relconstituents(ibnd,1)/
     .   matsurf_relconstituents(ibnd,1)
      nz = matsurf_nconstituents(ibnd)
cxpb The bulk or coating surface parameters
      if (coating_thickness(ibnd).gt.0.0_R8 .and.
     &    target_depth(ibnd,1).lt.coating_thickness(ibnd)) then
        a(1) = coating_Kappa_A1(ibnd)
        a(2) = coating_Kappa_A2(ibnd)
        a(3) = coating_Kappa_A3(ibnd)
        a(4) = coating_Kappa_A4(ibnd)
        a(5) = coating_Kappa_A5(ibnd)
        if (coating_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (coating_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab ('Unrecognized temperature units (coat Kappa) !')
        endif
        coating_value = a(1) +
     .   (a(2) + (a(3) + (a(4) +  a(5) * Teff) * Teff) * Teff) * Teff
        layer_Kappa = coating_value*frac
      else
        a(1) = matbulk_Kappa_A1(ibnd)
        a(2) = matbulk_Kappa_A2(ibnd)
        a(3) = matbulk_Kappa_A3(ibnd)
        a(4) = matbulk_Kappa_A4(ibnd)
        a(5) = matbulk_Kappa_A5(ibnd)
        if (bulk_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (bulk_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab ('Unrecognized temperature units (bulk Kappa) !')
        endif
        bulk_value = a(1) +
     .   (a(2) + (a(3) + (a(4) +  a(5) * Teff) * Teff) * Teff) * Teff
        layer_Kappa = bulk_value*frac
      end if
      if (layer_nconstituents(ibnd).eq.nz) return
cxpb Loop over deposited species contributions
      do itrack = 1, ntargsp
        frac = layer_relconstituents(ibnd,nz+itrack)
        a(1) = species_Kappa_A1(itrack)
        a(2) = species_Kappa_A2(itrack)
        a(3) = species_Kappa_A3(itrack)
        a(4) = species_Kappa_A4(itrack)
        a(5) = species_Kappa_A5(itrack)
        if (species_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (species_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab
     .     ('Unrecognized temperature units (species Kappa) !')
        endif
        layer_value = a(1) +
     .   (a(2) + (a(3) + (a(4) +  a(5) * Teff) * Teff) * Teff) * Teff
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
      real (kind=R8) :: T, dlayer_KappabydT
      integer nz, itrack
      real (kind=R8) :: Tcels, Teff, a(2:5)
      real (kind=R8) :: bulk_value, coating_value, layer_value, frac

      Tcels = T - 273.15_R8

cxpb First count the bulk contribution
      frac=layer_relconstituents(ibnd,1)/
     .   matsurf_relconstituents(ibnd,1)
      nz = matsurf_nconstituents(ibnd)
cxpb The bulk or coating surface parameters
      if (coating_thickness(ibnd).gt.0.0_R8 .and.
     &    target_depth(ibnd,1).lt.coating_thickness(ibnd)) then
        a(2) = coating_Kappa_A2(ibnd)
        a(3) = coating_Kappa_A3(ibnd)
        a(4) = coating_Kappa_A4(ibnd)
        a(5) = coating_Kappa_A5(ibnd)
        if (coating_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (coating_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab
     .     ('Unrecognized temperature units (coat dKappadT) !')
        endif
        coating_value = a(2) +
     .   (2.0_R8*a(3) + (3.0_R8*a(4) +
     .                   4.0_R8*a(5) * Teff) * Teff) * Teff
        dlayer_KappabydT = coating_value*frac
      else
        a(2) = matbulk_Kappa_A2(ibnd)
        a(3) = matbulk_Kappa_A3(ibnd)
        a(4) = matbulk_Kappa_A4(ibnd)
        a(5) = matbulk_Kappa_A5(ibnd)
        if (bulk_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (bulk_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab
     .     ('Unrecognized temperature units (bulk dKappadT) !')
        endif
        bulk_value = a(2) +
     .   (2.0_R8*a(3) + (3.0_R8*a(4) +
     .                   4.0_R8*a(5) * Teff) * Teff) * Teff
        dlayer_KappabydT = bulk_value*frac
      end if
      if (layer_nconstituents(ibnd).eq.nz) return
cxpb Loop over deposited species contributions
      do itrack = 1, ntargsp
        frac = layer_relconstituents(ibnd,nz+itrack)
        a(2) = species_Kappa_A2(itrack)
        a(3) = species_Kappa_A3(itrack)
        a(4) = species_Kappa_A4(itrack)
        a(5) = species_Kappa_A5(itrack)
        if (species_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (species_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab
     .     ('Unrecognized temperature units (species dKappadT) !')
        endif
        layer_value = a(2) +
     .   (2.0_R8*a(3) + (3.0_R8*a(4) +
     .                   4.0_R8*a(5) * Teff) * Teff) * Teff
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
      real (kind=R8) :: T, layer_Cp
      integer nz, itrack, mf, cf, sf
      real (kind=R8) :: Tcels, Teff
      real (kind=R8) :: b(4)
      real (kind=R8) :: bulk_value, coating_value, layer_value, frac
      external xerrab

      Tcels = T - 273.15_R8

cxpb First count the bulk contribution
      frac=layer_relconstituents(ibnd,1)/
     .   matsurf_relconstituents(ibnd,1)
      nz = matsurf_nconstituents(ibnd)
cxpb The bulk or coating surface parameters
      if (coating_thickness(ibnd).gt.0.0_R8 .and.
     &    target_depth(ibnd,1).lt.coating_thickness(ibnd)) then
        b(1) = coating_Cp_B1(ibnd)
        b(2) = coating_Cp_B2(ibnd)
        b(3) = coating_Cp_B3(ibnd)
        b(4) = coating_Cp_B4(ibnd)
        cf = coating_format(ibnd)
        if (coating_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (coating_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab
     .     ('Unrecognized temperature units (coat Cp) !')
        endif
        if (cf.eq.1) then
          coating_value = b(1) +
     .     ( b(2) + ( b(3) + b(4) * Teff) * Teff ) * Teff
        else
          coating_value = b(1) + b(2) * Teff - b(3)/Teff**2
        end if
        layer_Cp = coating_value * frac
      else
        b(1) = matbulk_Cp_B1(ibnd)
        b(2) = matbulk_Cp_B2(ibnd)
        b(3) = matbulk_Cp_B3(ibnd)
        b(4) = matbulk_Cp_B4(ibnd)
        mf = matbulk_format(ibnd)
        if (bulk_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (bulk_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab
     .     ('Unrecognized temperature units (bulk Cp) !')
        endif
        if (mf.eq.1) then
          bulk_value = b(1) +
     .     ( b(2) + ( b(3) + b(4) * Teff) * Teff ) * Teff
        else
          bulk_value = b(1) + b(2) * Teff - b(3)/Teff**2
        end if
        layer_Cp = bulk_value * frac
      end if
cxpb Loop over deposited species contributions
      if (layer_nconstituents(ibnd).gt.nz) then
        do itrack = 1, ntargsp
          frac = layer_relconstituents(ibnd,nz+itrack)
          b(1) = species_Cp_B1(itrack)
          b(2) = species_Cp_B2(itrack)
          b(3) = species_Cp_B3(itrack)
          b(4) = species_Cp_B4(itrack)
          sf = species_cp_format(itrack)
          if (species_temp_unit(ibnd).eq.'K') then
            Teff = T
          else if (species_temp_unit(ibnd).eq.'C') then
            Teff = Tcels
          else
            call xerrab
     .       ('Unrecognized temperature units (species Cp) !')
          endif
          if (sf.eq.1) then
            layer_value = b(1) +
     .     ( b(2) + ( b(3) + b(4) * Teff) * Teff ) * Teff
          else
            layer_value = b(1) + b(2) * T - b(3)/T**2
          end if
          layer_Cp = layer_Cp + layer_value * frac
        end do
      end if

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

cxpb First count the bulk contribution
      frac=layer_relconstituents(ibnd,1)/
     .   matsurf_relconstituents(ibnd,1)
      nz = matsurf_nconstituents(ibnd)
cxpb The bulk or coating surface parameters
      layer_emissivity = matsurf_blackbody_emiss(ibnd)*frac
      if (layer_nconstituents(ibnd).eq.nz) return
cxpb Loop over deposited species contributions
      do itrack = 1, ntargsp
        frac = layer_relconstituents(ibnd,nz+itrack)
        layer_value = species_blackbody_emiss(itrack)
        layer_emissivity = layer_emissivity + layer_value*frac
      end do

      return
      end function layer_emissivity

      logical function overlap(x1,y1,x2,y2,x3,y3,x4,y4,
     .                         x5,y5,x6,y6,x7,y7,x8,y8)
      use b2mod_types
      implicit none
      real (kind=R8) :: x1,y1,x2,y2,x3,y3,x4,y4,
     .                  x5,y5,x6,y6,x7,y7,x8,y8
      logical inside_rect_strict, intersect
      external inside_rect_strict, intersect
cxpb We aim to find if the 1234 quadrangle has an overlap with the 5678 rectangle
cxpb Sharing a vertex or part of an edge is not considered an overlap here!

cxpb We first test to see if any of the corners of 1234 are inside 5678
      overlap = inside_rect_strict(x1,y1,x5,y5,x6,y6,x7,y7,x8,y8) .or.
     .          inside_rect_strict(x2,y2,x5,y5,x6,y6,x7,y7,x8,y8) .or.
     .          inside_rect_strict(x3,y3,x5,y5,x6,y6,x7,y7,x8,y8) .or.
     .          inside_rect_strict(x4,y4,x5,y5,x6,y6,x7,y7,x8,y8)
      if (overlap) return

cxpb We first test to see if any of the corners of 5678 are inside 1234
      overlap = inside_rect_strict(x5,y5,x1,y1,x2,y2,x3,y3,x4,y4) .or.
     .          inside_rect_strict(x6,y6,x1,y1,x2,y2,x3,y3,x4,y4) .or.
     .          inside_rect_strict(x7,y7,x1,y1,x2,y2,x3,y3,x4,y4) .or.
     .          inside_rect_strict(x8,y8,x1,y1,x2,y2,x3,y3,x4,y4)
      if (overlap) return

cxpb We must now test the intersection of side 12 with 5678
      overlap = intersect(x1,y1,x2,y2,x5,y5,x6,y6) .or.
     .          intersect(x1,y1,x2,y2,x6,y6,x7,y7) .or.
     .          intersect(x1,y1,x2,y2,x7,y7,x8,y8) .or.
     .          intersect(x1,y1,x2,y2,x8,y8,x5,y5)
      if (overlap) return

cxpb We must now test the intersection of side 23 with 5678
      overlap = intersect(x2,y2,x3,y3,x5,y5,x6,y6) .or.
     .          intersect(x2,y2,x3,y3,x6,y6,x7,y7) .or.
     .          intersect(x2,y2,x3,y3,x7,y7,x8,y8) .or.
     .          intersect(x2,y2,x3,y3,x8,y8,x5,y5)
      if (overlap) return

cxpb We must now test the intersection of side 34 with 5678
      overlap = intersect(x3,y3,x4,y4,x5,y5,x6,y6) .or.
     .          intersect(x3,y3,x4,y4,x6,y6,x7,y7) .or.
     .          intersect(x3,y3,x4,y4,x7,y7,x8,y8) .or.
     .          intersect(x3,y3,x4,y4,x8,y8,x5,y5)
      if (overlap) return

cxpb We must now test the intersection of side 41 with 5678
      overlap = intersect(x4,y4,x1,y1,x5,y5,x6,y6) .or.
     .          intersect(x4,y4,x1,y1,x6,y6,x7,y7) .or.
     .          intersect(x4,y4,x1,y1,x7,y7,x8,y8) .or.
     .          intersect(x4,y4,x1,y1,x8,y8,x5,y5)
      return
      end function overlap

      logical function intersect(x1,y1,x2,y2,x3,y3,x4,y4)
      use b2mod_types
      implicit none
      real (kind=R8) :: x1,y1,x2,y2,x3,y3,x4,y4
cxpb We aim to find if segment 12 intersects with segment 34

cxpb Check if vector 34 lies between 31 and 32 and vector 12 is between 13 and 14
      intersect = (((x4-x3)*(y2-y3)-(x2-x3)*(y4-y3))*
     .             ((x4-x3)*(y1-y3)-(x1-x3)*(y4-y3))).lt.0.0_R8.and.
     .            (((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))*
     .             ((x2-x1)*(y4-y1)-(x4-x1)*(y2-y1))).lt.0.0_R8
      return
      end function intersect

      logical function inside_rect_strict(x0,y0,x1,y1,x2,y2,x3,y3,x4,y4)
      use b2mod_types
      implicit none
      real (kind=R8) :: x0,y0,x1,y1,x2,y2,x3,y3,x4,y4
cxpb We aim to find if point 0 is strictly inside the rectangle 1234

      inside_rect_strict =
     . (((x0-x1)*(y1-y2)+(y0-y1)*(x2-x1)).gt.0.0_R8 .and.
     .  ((x0-x2)*(y2-y3)+(y0-y2)*(x3-x2)).gt.0.0_R8 .and.
     .  ((x0-x3)*(y3-y4)+(y0-y3)*(x4-x3)).gt.0.0_R8 .and.
     .  ((x0-x4)*(y4-y1)+(y0-y4)*(x1-x4)).gt.0.0_R8) .or.
     . (((x0-x1)*(y1-y2)+(y0-y1)*(x2-x1)).lt.0.0_R8 .and.
     .  ((x0-x2)*(y2-y3)+(y0-y2)*(x3-x2)).lt.0.0_R8 .and.
     .  ((x0-x3)*(y3-y4)+(y0-y3)*(x4-x3)).lt.0.0_R8 .and.
     .  ((x0-x4)*(y4-y1)+(y0-y4)*(x1-x4)).lt.0.0_R8)
      return
      end function inside_rect_strict
