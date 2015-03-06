      subroutine b2yt_wall_init(nncut1,nx1,ny1,nncut,nx,ny,nlim,nsts)
      use b2mod_types
      use b2mod_geo2
      use b2mod_wall
      use b2mod_wall2
      use b2mod_layer
      use b2mod_b2cmpb
#     include "DIMENSIONS.F"
      implicit none
      integer nncut1, nx1, ny1, nncut, nx, ny, nlim, nsts
      integer ix, iy, ixn, iyn, is, ibnd, ibnd1, idir, ireg, 
     .        nxtl1, nxtr1
      real (kind=R8) :: ub, lb, ub1, lb1
      external subini, subend, xertst, get_nxt

      call subini('b2yt_wall_init')
      call alloc_b2mod_layer(nx,ny,nncut,nlim,nsts)
      call alloc_b2mod_wall_init(nx,ny,nncut)
      call alloc_b2mod_wall2_init(nncut1,nx1,ny1,nlim,nsts)
      call init_wall(nx,ny)
cxpb Establish numbering of boundary wall cells
cxpb We must avoid corner cells and dead cells
      nbnd1 = 0
      wall_begin1 = 0
      wall_end1 = 0
      nsegments = 0

cxpb  We follow the same order as in init_wall
cxpb  Starting for HFS lower corner and going clockwise around the plasma
cxpb  keeping the plasma to our right (but we are moving in the guard cells!)
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
        call b2yt_nxt(nx1,nxtl1,nxtr1)
        cbrbrk(cbirso) = real(leftcut1(1))/real(nx1)
        cbrbrk(cbirso+1) = real(rightcut1(2))/real(nx1) 
        cbrbrk(cbirso+2) = real(nxtl1)/real(nx1) 
        cbrbrk(cbirso+3) = real(leftcut1(2))/real(nx1) 
        cbrbrk(cbirso+4) = real(rightcut1(2))/real(nx1) 
      endif
      cbrbrk(cbirso+cbnrso-1) = 2.00_R8
      if (nncut.eq.2) then
        cbrbrk(cbirno) = real(nxtl1)/real(nx1) 
      endif
      cbrbrk(cbirno+cbnrno-1) = 2.00_R8
      cbrbrk(cbirwe:cbirwe+cbnrwe-1) = 2.00_R8
      cbrbrk(cbirea:cbirea+cbnrea-1) = 2.00_R8

cxpb  Northern boundaries (excluding diverted island)
      iy = ny1
      ix = -1
      ixn = rightix1(ix,iy)
      iyn = rightiy1(ix,iy)
      idir = 2  ! 0 : left, 1: bottom, 2: right, 3: top
      nsegments = 1
      wall_begin1(1) = 1
      do while (xymap1(ixn,iyn).eq.0)
        nbnd1 = nbnd1 + 1
        imapx1(nbnd1) = ixn
        imapy1(nbnd1) = iyn
        xymap1(ixn,iyn)= nbnd1
        ix = ixn
        iy = iyn
        DIRVALUE: select case (idir)
          case(0)
            r1(nbnd1)=crx1(ix,iy,3)
            z1(nbnd1)=cry1(ix,iy,3)
            r2(nbnd1)=crx1(ix,iy,2)
            z2(nbnd1)=cry1(ix,iy,2)
            plate_area1(nbnd1) = gs1(topix1(ix,iy),topiy1(ix,iy),1)
            ixn = leftix1(ix,iy)
            iyn = leftiy1(ix,iy)
            if (leftix1(ixn,iyn).eq.-2) then 
              idir = 3
              ixn = leftix1(topix1(ix,iy),topiy1(ix,iy))
              iyn = leftiy1(topix1(ix,iy),topiy1(ix,iy))
              wall_end1(nsegments) = nbnd1
              nsegments = nsegments + 1
              wall_begin1(nsegments) = nbnd1+1
            endif 
          case(1)
            plate_area1(nbnd1) = gs1(ix,iy,0)
            r1(nbnd1)=crx1(ix,iy,2)
            z1(nbnd1)=cry1(ix,iy,2)
            r2(nbnd1)=crx1(ix,iy,0)
            z2(nbnd1)=cry1(ix,iy,0)
            ixn = bottomix1(ix,iy)
            iyn = bottomiy1(ix,iy)
            if (bottomiy1(ixn,iyn).eq.-2) then
              idir = 0
              ixn = bottomix1(leftix1(ix,iy),leftiy1(ix,iy))
              iyn = bottomiy1(leftix1(ix,iy),leftiy1(ix,iy))
              wall_end1(nsegments) = nbnd1
              nsegments = nsegments + 1
              wall_begin1(nsegments) = nbnd1+1
            endif 
          case(2)
            plate_area1(nbnd1) = gs1(ix,iy,1)
            r1(nbnd1)=crx1(ix,iy,0)
            z1(nbnd1)=cry1(ix,iy,0)
            r2(nbnd1)=crx1(ix,iy,1)
            z2(nbnd1)=cry1(ix,iy,1)
            ixn = rightix1(ix,iy)
            iyn = rightiy1(ix,iy)
            if (rightix1(ixn,iyn).eq.nx1+1) then
              idir = 1
              ixn = rightix1(bottomix1(ix,iy),bottomiy1(ix,iy))
              iyn = rightiy1(bottomix1(ix,iy),bottomiy1(ix,iy))
              wall_end1(nsegments) = nbnd1
              nsegments = nsegments + 1
              wall_begin1(nsegments) = nbnd1+1
            endif
          case(3)
            plate_area1(nbnd1) = gs1(rightix1(ix,iy),rightiy1(ix,iy),0)
            r1(nbnd1)=crx1(ix,iy,1)
            z1(nbnd1)=cry1(ix,iy,1)
            r2(nbnd1)=crx1(ix,iy,3)
            z2(nbnd1)=cry1(ix,iy,3)
            ixn = topix1(ix,iy)
            iyn = topiy1(ix,iy)
            if (topiy1(ixn,iyn).eq.ny1+1) then
              idir = 2
              ixn = topix1(rightix1(ix,iy),rightiy1(ix,iy))
              iyn = topiy1(rightix1(ix,iy),rightiy1(ix,iy))
              wall_end1(nsegments) = nbnd1
              nsegments = nsegments + 1
              wall_begin1(nsegments) = nbnd1+1
            endif
          case default DIRVALUE
            call xerrab ('Invalid direction')
        end select DIRVALUE
      enddo
      nsegments = nsegments - 1      

cxpb  We now identify the b2ah region numbers
      do ibnd = 1, nbnd1
        ix = imapx1(ibnd)
        iy = imapy1(ibnd)
        if(iy.eq.ny1) then
          ireg = cbirno
          do while (cbrbrk(ireg)*nx1.lt.ix+0.5_R8 .and. 
     &                     ireg.le.cbirno+cbnrno-1)
            ireg = ireg + 1
          enddo
          imapreg(ibnd) = ireg
        elseif (ix.eq.nx1) then
          ireg = cbirea
          do while (cbrbrk(ireg)*ny1.lt.iy+0.5_R8 .and. 
     &                     ireg.le.cbirea+cbnrea-1)
            ireg = ireg + 1
          enddo
          imapreg(ibnd) = ireg
        elseif (iy.eq.-1) then
          ireg = cbirso
          do while (cbrbrk(ireg)*nx1.lt.ix+0.5_R8 .and. 
     &                     ireg.le.cbirso+cbnrso-1)
            ireg = ireg + 1
          enddo
          imapreg(ibnd) = ireg
        elseif (ix.eq.-1) then
          ireg = cbirwe
          do while (cbrbrk(ireg)*ny1.lt.iy+0.5_R8 .and. 
     &                     ireg.le.cbirwe+cbnrwe-1)
            ireg = ireg + 1
          enddo
          imapreg(ibnd) = ireg
        else
          imapreg(ibnd) = -1
        endif
      enddo
      
      do is = nsegments+1, nwall_segments
        wall_begin1(is) = nbnd1+1
      enddo
#ifndef ALLOCATE_AND_NAMELIST
      call xertst (nbnd1.le.2*(DEF_NXD+2)+(DEF_NCUT+2)*(DEF_NYD+2),
     & 'Increase size of B2 contribution to nwall in b2mod_wall')
#endif
 
cxpb We map per segment
cxpb We interpolate according to the area along the segment
      do is= 1, nsegments
        acc_area(wall_begin(is)) = plate_area(wall_begin(is))
        do ibnd = wall_begin(is)+1, wall_end(is)
          acc_area(ibnd) = acc_area(ibnd-1) + plate_area(ibnd)
        enddo
        do ibnd = wall_begin(is), wall_end(is)-1
          acc_area(ibnd) = acc_area(ibnd)/acc_area(wall_end(is))
        end do
        acc_area(wall_end(is)) = 1.0_R8
        acc_area1(wall_begin1(is)) = plate_area1(wall_begin1(is))
        do ibnd = wall_begin1(is)+1, wall_end1(is)
          acc_area1(ibnd) = acc_area1(ibnd-1) + plate_area1(ibnd)
        enddo
        do ibnd = wall_begin1(is), wall_end1(is)-1
          acc_area1(ibnd) = acc_area1(ibnd)/acc_area1(wall_end1(is))
        end do
        acc_area1(wall_end1(is)) = 1.0_R8
cxpb The contr_area(ibnd,ibnd1) contains the fraction of the plate area from 
cxpb old wall element ibnd that is to be found as part of new wall element ibnd1
cxpb lb and ub and lower and upper interval bounds
        ibnd = wall_begin(is)
        lb = 0.0_R8
        ibnd1 = wall_begin1(is)
        lb1 = 0.0_R8
        do while (ibnd.le.wall_end(is) .and. ibnd1.le.wall_end1(is))
          ub = acc_area(ibnd)
          ub1 = acc_area1(ibnd1)
          if (lb.le.lb1 .and. ub.lt.ub1) then
            contr_area(ibnd,ibnd1) = (ub-lb1)/(ub-lb)
            ibnd = ibnd+1
            lb = ub
          else if (lb.le.lb1 .and. ub.gt.ub1) then
            contr_area(ibnd,ibnd1) = (ub1-lb1)/(ub-lb)
            ibnd1 = ibnd1 + 1
            lb1 = ub1
          else if (lb.le.lb1 .and. ub.eq.ub1) then
            contr_area(ibnd,ibnd1) = (ub-lb1)/(ub-lb)
            ibnd = ibnd + 1
            lb = ub
            ibnd1 = ibnd1 + 1
            lb1 = ub1
          else if (lb.gt.lb1 .and. ub.lt.ub1) then
            contr_area(ibnd,ibnd1) = 1.0_R8
            ibnd = ibnd + 1
            lb = ub
          else if (lb.gt.lb1 .and. ub.gt.ub1) then
            contr_area(ibnd,ibnd1) = (ub1-lb)/(ub-lb)
            ibnd1 = ibnd1 + 1
            lb1 = ub1
          else if (lb.gt.lb1 .and. ub.eq.ub1) then 
            contr_area(ibnd,ibnd1) = 1.0_R8
            ibnd = ibnd + 1
            lb = ub
            ibnd1 = ibnd1 + 1 
            lb1 = ub1
          end if
        end do
      end do  ! nsegments

      call subend()
      return
      end subroutine b2yt_wall_init

      subroutine b2yt_wall(nx1,ny1,nx,ny,ns,ns1,natm,natm1,nmol,nion,
     *   is_neutral1,iatm1,bundle_frac,zn1,l_chem_sput,
     *   itrack1,itrack2,itrack3)
      use b2mod_types
      use b2mod_wall
      use b2mod_wall2
      use b2mod_layer
      use b2mod_indirect
      use b2mod_b2cmpa
      use b2mod_sputter
#     include "DIMENSIONS.F"
      implicit none

      integer nx1, ny1, nx, ny, ns, ns1, natm, natm1, nmol, nion
      integer itrack1, itrack2, itrack3
      integer iatm1(natm)
      real (kind=R8) :: bundle_frac(0:ns1-1,0:ns-1), zn1(0:ns1-1)
      logical is_neutral1(0:ns1-1), l_chem_sput

      integer ix, iy, iz, is, ibnd, ibnd1, imax, i, i1, itrack
      integer ibegin, iend, iw, iw1
      integer ixa, ixb, iya, iyb, iza, ndepth1
      integer nsm, natmm, nbndm
      real (kind=R8) :: sum
      real (kind=R8), allocatable :: wall_data(:,:)
      real (kind=R8) :: data_tmp(nwall), xx(nwall+2), xx1(nwall+2) 
      real (kind=R8) :: track_data(nwall), track_source(nwall)
      real (kind=R8), allocatable :: yy(:), yy1(:)
      real (kind=R8), allocatable :: fa(:,:)
      real (kind=R8) :: fb(3,ndepth)
      logical wall_namelist_used, layer_model
      save layer_model
      data layer_model/.false./
      character*400 line
      character*17 numstr
      logical found

      intrinsic mod, max, min, sqrt
      external xertst, subini, subend
      external convert_bundle_r

      call subini('b2yt_wall')
      call alloc_b2mod_wall(ns,natm,nmol,nion)
      
      inquire(file='../b2.wall_save.parameters',
     1 exist=wall_namelist_used)
      if(wall_namelist_used) then
        call read_b2mod_wall_namelist(99,ns,'    ',wall_namelist_used)
        call xertst(ndepth_nml.le.ndepth .and. ndepth_nml.ge.1,
     .   'Invalid array size NDEPTH_NML!')
        ndepth1 = ndepth_nml
        call ipgeti('b2ytdr_ndepth1', ndepth1)
        call xertst(ndepth1.le.ndepth .and. ndepth1.ge.1,
     .   'Invalid array size NDEPTH1!')
        call alloc_b2mod_wall2(ns1,natm1,nmol,nion,ndepth1)
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
      else
        return
      endif

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
      do ibnd = 1, nbnd
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
      natmm = max(natm,natm1)
      nbndm = max(nbnd,nbnd1)
      allocate(wall_data(nwall,0:nsm-1+natmm+nmol+DEF_NION))
cxpb We map per segment
cxpb We interpolate according to the area along the segment
      allocate(yy(ndepth_nml),yy1(ndepth1))
cxpb For string-valued quantities, we choose to inherit from
cxpb the old element with the largest areal contribution
      do is = 1, nsegments
        do ibnd1 = wall_begin1(is), wall_end1(is)
          imax = wall_begin(is)
          do ibnd = wall_begin(is)+1, wall_end(is)
            if (contr_area(ibnd,ibnd1).gt.contr_area(imax,ibnd1)) 
     >       imax = ibnd
          end do
          surface_material_name1(ibnd1) = surface_material_name(imax)
          coating_material_name1(ibnd1) = coating_material_name(imax)
          bulk_material_name1(ibnd1) = bulk_material_name(imax)
          inertial_cooling1(ibnd1) = inertial_cooling(imax)
          layer_nconstituents1(ibnd1) = matsurf_nconstituents(imax)
          layer_nzconstituents1(ibnd1,1:layer_nconstituents1(ibnd1)) =
     1     matsurf_nzconstituents(imax,1:layer_nconstituents1(ibnd1))
          layer_nrelconstituents1(ibnd1,1:layer_nconstituents1(ibnd1)) =
     1     matsurf_nrelconstituents(imax,1:layer_nconstituents1(ibnd1))
cxpb  normalise the layer constituents relative concentrations
          sum = 0.0_R8
          do i = 1, layer_nconstituents1(ibnd1)
            sum = sum + layer_nrelconstituents1(ibnd1,i)
          enddo
          if (sum.gt.0.0_R8)
     1     layer_nrelconstituents1(ibnd1,1:layer_nconstituents1(ibnd1))=
     2     layer_nrelconstituents1(ibnd1,1:layer_nconstituents1(ibnd1))/
     3     sum
        end do
      end do  ! nsegments
cxpb Copy data arrays into new geometry
      call b2yt_wall_average(nwall,contr_area,plate_area,
     . plate_thickness,plate_thickness1)
      call b2yt_wall_average(nwall,contr_area,plate_area,
     . coating_thickness,coating_thickness1)
      call b2yt_wall_average(nwall,contr_area,plate_area,
     . plate_time_factor,plate_time_factor1)
      call b2yt_wall_average(nwall,contr_area,plate_area,
     . plate_time,plate_time1)
      call b2yt_wall_average(nwall,contr_area,plate_area,
     . backplate_temp,backplate_temp1)
      if (itrack1.gt.0) then
        data_tmp(1:nbnd) = deposition(1:nbnd,1)
        call b2yt_wall_sum(nwall,contr_area,data_tmp)
        deposition1(1:nbnd1,itrack1) = data_tmp(1:nbnd1)
        data_tmp(1:nbnd) = monolayer_deposition(1:nbnd,1)
        call b2yt_wall_sum(nwall,contr_area,data_tmp)
        monolayer_deposition1(1:nbnd1,itrack1) = data_tmp(1:nbnd1)
        data_tmp(1:nbnd) = erosion(1:nbnd,1)
        call b2yt_wall_sum(nwall,contr_area,data_tmp)
        erosion1(1:nbnd1,itrack1) = data_tmp(1:nbnd1)
        data_tmp(1:nbnd) = monolayer_erosion(1:nbnd,1)
        call b2yt_wall_sum(nwall,contr_area,data_tmp)
        monolayer_erosion1(1:nbnd1,itrack1) = data_tmp(1:nbnd1)
        do iw1 = 1, nbnd1
          layer_nconstituents1(iw1) = layer_nconstituents1(iw1)+1
          layer_nzconstituents1(iw1,layer_nconstituents1(iw1)) =
     .     zn1(sput_dst)
        end do
        do iw = 1, nbnd
          track_source(iw) = 
     .     layer_nrelconstituents(iw,matsurf_nconstituents(iw)+1)
        end do
        call b2yt_wall_average(nwall,contr_area,plate_area,
     .   track_source,track_data)
        do iw1 = 1, nbnd1
          layer_nrelconstituents1(iw1,layer_nconstituents1(iw1)) =
     .     track_data(iw1)
        end do
      end if
      if (itrack2.gt.0) then
        data_tmp(1:nbnd) = deposition(1:nbnd,2)
        call b2yt_wall_sum(nwall,contr_area,data_tmp)
        deposition1(1:nbnd1,itrack2) = data_tmp(1:nbnd1)
        data_tmp(1:nbnd) = monolayer_deposition(1:nbnd,2)
        call b2yt_wall_sum(nwall,contr_area,data_tmp)
        monolayer_deposition1(1:nbnd1,itrack2) = data_tmp(1:nbnd1)
        data_tmp(1:nbnd) = erosion(1:nbnd,2)
        call b2yt_wall_sum(nwall,contr_area,data_tmp)
        erosion1(1:nbnd1,itrack2) = data_tmp(1:nbnd1)
        data_tmp(1:nbnd) = monolayer_erosion(1:nbnd,2)
        call b2yt_wall_sum(nwall,contr_area,data_tmp)
        monolayer_erosion1(1:nbnd1,itrack2) = data_tmp(1:nbnd1)
        do iw1 = 1, nbnd1
          layer_nconstituents1(iw1) = layer_nconstituents1(iw1)+1
          if (itrack2.eq.2) then
            layer_nzconstituents1(iw1,layer_nconstituents1(iw1)) =
     .       zn1(sput_dst2)
          else if (itrack2.eq.1) then
            layer_nzconstituents1(iw1,layer_nconstituents1(iw1)) =
     .       zn1(sput_dst)
          end if
        end do
        do iw = 1, nbnd
          track_source(iw) = 
     .     layer_nrelconstituents(iw,matsurf_nconstituents(iw)+2)
        end do
        call b2yt_wall_average(nwall,contr_area,plate_area,
     .   track_source,track_data)
        do iw1 = 1, nbnd1
          layer_nrelconstituents1(iw1,layer_nconstituents1(iw1)) =
     .     track_data(iw1)
        end do
      end if
      if (itrack3.gt.0) then
        data_tmp(1:nbnd) = deposition(1:nbnd,3)
        call b2yt_wall_sum(nwall,contr_area,data_tmp)
        deposition1(1:nbnd1,itrack3) = data_tmp(1:nbnd1)
        data_tmp(1:nbnd) = monolayer_deposition(1:nbnd,3)
        call b2yt_wall_sum(nwall,contr_area,data_tmp)
        monolayer_deposition1(1:nbnd1,itrack3) = data_tmp(1:nbnd1)
        data_tmp(1:nbnd) = erosion(1:nbnd,3)
        call b2yt_wall_sum(nwall,contr_area,data_tmp)
        erosion1(1:nbnd1,itrack3) = data_tmp(1:nbnd1)
        data_tmp(1:nbnd) = monolayer_erosion(1:nbnd,3)
        call b2yt_wall_sum(nwall,contr_area,data_tmp)
        monolayer_erosion1(1:nbnd1,itrack3) = data_tmp(1:nbnd1)
        do iw1 = 1, nbnd1
          layer_nconstituents1(iw1) = layer_nconstituents1(iw1)+1
          if (itrack3.eq.3) then
            layer_nzconstituents1(iw1,layer_nconstituents1(iw1)) =
     .       zn1(sput_dst3)
          else if (itrack3.eq.2) then
            layer_nzconstituents1(iw1,layer_nconstituents1(iw1)) =
     .       zn1(sput_dst2)
          else if (itrack3.eq.1) then
            layer_nzconstituents1(iw1,layer_nconstituents1(iw1)) =
     .       zn1(sput_dst)
          endif
        end do
        do iw = 1, nbnd
          track_source(iw) = 
     .     layer_nrelconstituents(iw,matsurf_nconstituents(iw)+3)
        end do
        call b2yt_wall_average(nwall,contr_area,plate_area,
     .   track_source,track_data)
        do iw1 = 1, nbnd1
          layer_nrelconstituents1(iw1,layer_nconstituents1(iw1)) =
     .     track_data(iw1)
        end do
      endif
      do ibnd1 = 1, nbnd1
        sum = 0.0_R8
        do i = 1, layer_nconstituents1(ibnd1)
          sum = sum + layer_nrelconstituents1(ibnd1,i)
        end do
        layer_nrelconstituents1(ibnd1,1:layer_nconstituents1(ibnd1))=
     .   layer_nrelconstituents1(ibnd1,1:layer_nconstituents1(ibnd1))/
     .   sum
      end do
      wall_data(1:nbnd,0:ns-1+natm+nmol+nion) = 
     . physical_sputtering(1:nbnd,0:ns-1+natm+nmol+nion)
      call b2yt_wall_copy(ns,ns1,nsm,is_neutral1,bundle_frac,
     . natm,natm1,natmm,nmol,nion,iatm1,nwall,wall_data,
     . contr_area,plate_area)
      physical_sputtering1(1:nbnd1,0:ns1-1+natm1+nmol+nion) =
     . wall_data(1:nbnd1,0:ns1-1+natm1+nmol+nion)
      wall_data(1:nbnd,0:ns-1+natm+nmol+nion) = 
     . physical_sputtering_energy(1:nbnd,0:ns-1+natm+nmol+nion)
      call b2yt_wall_copy(ns,ns1,nsm,is_neutral1,bundle_frac,
     . natm,natm1,natmm,nmol,nion,iatm1,nwall,wall_data,
     . contr_area,plate_area)
      physical_sputtering_energy1(1:nbnd1,0:ns1-1+natm1+nmol+nion) =
     . wall_data(1:nbnd1,0:ns1-1+natm1+nmol+nion)
      if (l_chem_sput) then
        wall_data(1:nbnd,0:ns-1+natm+nmol+nion) = 
     .   chemical_sputtering(1:nbnd,0:ns-1+natm+nmol+nion)
        call b2yt_wall_copy(ns,ns1,nsm,is_neutral1,bundle_frac,
     .   natm,natm1,natmm,nmol,nion,iatm1,nwall,wall_data,
     .   contr_area,plate_area)
        chemical_sputtering1(1:nbnd1,0:ns1-1+natm1+nmol+nion) =
     .   wall_data(1:nbnd1,0:ns1-1+natm1+nmol+nion)
        wall_data(1:nbnd,0:ns-1+natm+nmol+nion) = 
     .   RES_sputtering(1:nbnd,0:ns-1+natm+nmol+nion)
        call b2yt_wall_copy(ns,ns1,nsm,is_neutral1,bundle_frac,
     .   natm,natm1,natmm,nmol,nion,iatm1,nwall,wall_data,
     .   contr_area,plate_area)
        RES_sputtering1(1:nbnd1,0:ns1-1+natm1+nmol+nion) =
     .   wall_data(1:nbnd1,0:ns1-1+natm1+nmol+nion)
      endif
      wall_data(1:nbnd,0:ns-1+natm+nmol+nion) = 
     . thermal_evaporation(1:nbnd,0:ns-1+natm+nmol+nion)
      call b2yt_wall_copy(ns,ns1,nsm,is_neutral1,bundle_frac,
     . natm,natm1,natmm,nmol,nion,iatm1,nwall,wall_data,
     . contr_area,plate_area)
      thermal_evaporation1(1:nbnd1,0:ns1-1+natm1+nmol+nion) =
     . wall_data(1:nbnd1,0:ns1-1+natm1+nmol+nion)
      wall_data(1:nbnd,0:ns-1+natm+nmol+nion) = 
     . backscattering(1:nbnd,0:ns-1+natm+nmol+nion)
      call b2yt_wall_copy(ns,ns1,nsm,is_neutral1,bundle_frac,
     . natm,natm1,natmm,nmol,nion,iatm1,nwall,wall_data,
     . contr_area,plate_area)
      backscattering1(1:nbnd1,0:ns1-1+natm1+nmol+nion) =
     . wall_data(1:nbnd1,0:ns1-1+natm1+nmol+nion)
      wall_data(1:nbnd,0:ns-1+natm+nmol+nion) = 
     . backscattering_energy(1:nbnd,0:ns-1+natm+nmol+nion)
      call b2yt_wall_copy(ns,ns1,nsm,is_neutral1,bundle_frac,
     . natm,natm1,natmm,nmol,nion,iatm1,nwall,wall_data,
     . contr_area,plate_area)
      backscattering_energy1(1:nbnd1,0:ns1-1+natm1+nmol+nion) =
     . wall_data(1:nbnd1,0:ns1-1+natm1+nmol+nion)
      do is = 1, nsegments
        do ibnd1 = wall_begin1(is), wall_end1(is)
cxpb We do a simple 2d interpolation of the target temperature,
cxpb but in computational space
          ibegin = wall_begin(is)
          found = .false.
          do while (.not.found.and.ibegin.le.wall_end(is))
            found = contr_area(ibegin,ibnd1).gt.0.0_R8
            if (.not.found) ibegin = ibegin+1
          end do
          found = ibegin.le.wall_end(is)
          iend = ibegin
          do while (found.and.iend.lt.wall_end(is))
            iend = iend+1
            found = contr_area(iend,ibnd1).gt.0.0_R8
          end do
          if (.not.found.and.iend.le.wall_end(is)) iend = iend-1
          if (ibegin.gt.wall_end(is).or.iend.gt.wall_end(is)) then
            write(6,*) 'ibegin, iend : ', ibegin, iend
            write(6,*) 'contr_area:', 
     .       (contr_area(iw,ibnd1),iw=wall_begin(is),wall_end(is))
            write(line,'(a,i5,a)') 'Wall element ',ibnd1,' not found!'
            call xerrab (trim(line))
          end if
          allocate(fa(1:iend-ibegin+3,1:ndepth_nml))
          fa = 0.0_R8
          fb = 0.0_R8
          if (ibegin.eq.wall_begin(is)) then
            xx(1) = 0.0_R8
            fa(1,1:ndepth_nml) = target_temp(1,1:ndepth_nml)
          else
            xx(1) = acc_area(ibegin-1)
            fa(1,1:ndepth_nml) = target_temp(ibegin-1,1:ndepth_nml)
          endif
          do iw = ibegin, iend
            xx(iw+2-ibegin) = (acc_area(iw)+xx(iw+1-ibegin))/2.0_R8
            fa(iw+2-ibegin,1:ndepth_nml) = target_temp(iw,1:ndepth_nml)
          end do
          if (iend.eq.wall_end(is)) then
            xx(iend+3-ibegin) = 1.0_R8
            fa(iend+3-ibegin,1:ndepth_nml) = 
     .       target_temp(iend,1:ndepth_nml)
          else
            xx(iend+3-ibegin) = acc_area(iend+1)
            fa(iend+3-ibegin,1:ndepth_nml) = 
     .       target_temp(iend+1,1:ndepth_nml)
          endif
          if (ibnd1.eq.wall_begin1(is)) then
            xx1(1) = 0.0_R8
          else
            xx1(1)= acc_area1(ibnd1-1)
          endif
          xx1(2) = (acc_area1(ibnd1)+xx1(1))/2.0_R8
          xx1(3) = acc_area1(ibnd1)
          yy(1) = 0.0_R8
          do iz = 2, ndepth_nml
            yy(iz) = (real(iz)-1.0_R8)/(real(ndepth_nml)-1.0_R8)
          end do
          yy1(1) = 0.0_R8
          do iz = 2, ndepth1
            yy1(iz) = (real(iz)-1.0_R8)/(real(ndepth1)-1.0_R8)
          end do
          call interp2d(iend-ibegin+3,ndepth_nml,3,ndepth1,
     .     xx(1:iend-ibegin+3),yy(1:ndepth_nml),
     .     xx1(1:3),yy1(1:ndepth1),fa,fb(1:3,1:ndepth1))
          target_temp1(ibnd1,1:ndepth1) = fb(2,1:ndepth1)
          deallocate(fa)
        end do          
      end do

      deallocate(yy,yy1)
      deallocate(wall_data)
      call subend()
      return

* ..finalise time-dependent wall data
      entry b2yt_wall_final(nx1,ny1,ns1)
      wall_namelist_used = ntargsp.gt.0 .or.
     .  layer_thermal_contact.ne.1.0_R8 .or. 
     .  maxval(target_temp1(1:nbnd1,1)).gt.plate_temp .or.
     .  maxval(thermal_evaporation1(1:nbnd1,0:ns1-1)).gt.0.0_R8 .or.
     .  maxval(physical_sputtering1(1:nbnd1,0:ns1-1)).gt.0.0_R8 .or.
     .  maxval(chemical_sputtering1(1:nbnd1,0:ns1-1)).gt.0.0_R8 .or.
     .  maxval(RES_sputtering1(1:nbnd1,0:ns1-1)).gt.0.0_R8 .or.
     .  maxval(backscattering1(1:nbnd1,0:ns1-1)).gt.0.0_R8
      if (.not.wall_namelist_used) return
      open(99,file='../b2.wall_save.parameters.2')
      write(6,*) 'Writing ../b2.wall_save.parameters.2'
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
	do ibnd1 = 1, nbnd1
        if (target_temp1(ibnd1,1).gt.plate_temp) 
     >     write(99,'(a12,i4,a4,1x,1es15.8,a1)') 
     .     'TARGET_TEMP(',ibnd1,',1)=',target_temp1(ibnd1,1),','
        enddo
      endif
      do is = 0, ns1-1
        if (maxval(thermal_evaporation1(1:nbnd1,is)).gt.0.0_R8) then
          do i = 1, (nbnd1-1)/nline+1
            if (i.eq.1) then
              write(line,'(a22,i3,a2)') 'THERMAL_EVAPORATION(1,',is,')='
            else
              line = ''
            endif
            do ibnd1 = (i-1)*nline+1, min(i*nline,nbnd1)
              write(numstr,'(1x,1es15.8,a1)') 
     .         thermal_evaporation1(ibnd1,is),','
              line = trim(line)//numstr
            enddo
            write(99,'(a)') trim(line)
          enddo
        endif
        if (maxval(physical_sputtering1(1:nbnd1,is)).gt.0.0_R8) then
          do i = 1, (nbnd1-1)/nline+1
            if (i.eq.1) then
              write(line,'(a22,i3,a2)') 'PHYSICAL_SPUTTERING(1,',is,')='
            else
              line = ''
            endif
            do ibnd1 = (i-1)*nline+1, min(i*nline,nbnd1)
              write(numstr,'(1x,1e15.8,a1)') 
     .         physical_sputtering1(ibnd1,is),','
              line = trim(line)//numstr
            enddo
            write(99,'(a)') trim(line)
          enddo
          do i = 1, (nbnd1-1)/nline+1
            if (i.eq.1) then
              write(line,'(a29,i3,a2)') 
     .         'PHYSICAL_SPUTTERING_ENERGY(1,',is,')='
            else
              line = ''
            endif
            do ibnd1 = (i-1)*nline+1, min(i*nline,nbnd1)
              write(numstr,'(1x,1es15.8,a1)') 
     .         physical_sputtering_energy1(ibnd1,is),','
              line = trim(line)//numstr
            enddo
            write(99,'(a)') trim(line)
          enddo
        endif
        if (maxval(chemical_sputtering1(1:nbnd1,is)).gt.0.0_R8) then
          do i = 1, (nbnd1-1)/nline+1
            if (i.eq.1) then
              write(line,'(a22,i3,a2)') 
     .         'CHEMICAL_SPUTTERING(1,',is,')='
            else
              line = ''
            endif
            do ibnd1 = (i-1)*nline+1, min(i*nline,nbnd1)
              write(numstr,'(1x,1es15.8,a1)') 
     .         chemical_sputtering1(ibnd1,is),','
              line = trim(line)//numstr
            enddo
            write(99,'(a)') trim(line)
          enddo
        endif
        if (maxval(RES_sputtering1(1:nbnd1,is)).gt.0.0_R8) then
          do i = 1, (nbnd1-1)/nline+1
            if (i.eq.1) then
              write(line,'(a17,i3,a2)') 
     .         'RES_SPUTTERING(1,',is,')='
            else
              line = ''
            endif
            do ibnd1 = (i-1)*nline+1, min(i*nline,nbnd1)
              write(numstr,'(1x,1es15.8,a1)') 
     .         RES_sputtering1(ibnd1,is),','
              line = trim(line)//numstr
            enddo
            write(99,'(a)') trim(line)
          enddo
        endif
        if (maxval(backscattering1(1:nbnd1,is)).gt.0.0_R8) then
          do i = 1, (nbnd1-1)/nline+1
            if (i.eq.1) then
              write(line,'(a17,i3,a2)') 
     .         'BACKSCATTERING(1,',is,')='
            else
              line = ''
            endif
            do ibnd1 = (i-1)*nline+1, min(i*nline,nbnd1)
              write(numstr,'(1x,1es15.8,a1)') 
     .         backscattering1(ibnd1,is),','
              line = trim(line)//numstr
            enddo
            write(99,'(a)') trim(line)
          enddo
          do i = 1, (nbnd1-1)/nline+1
            if (i.eq.1) then
              write(line,'(a24,i3,a2)') 
     .         'BACKSCATTERING_ENERGY(1,',is,')='
            else
              line = ''
            endif
            do ibnd1 = (i-1)*nline+1, min(i*nline,nbnd1)
              write(numstr,'(1x,1es15.8,a1)') 
     .         backscattering_energy1(ibnd1,is),','
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
        do ibnd1 = 1, nbnd1
          write(99,'(a6,i4,a2,i4,a1,1x,a6,i4,a2,i4,a1)') 
     .     'IMAPX(',ibnd1,')=',imapx1(ibnd1),',',
     .     'IMAPY(',ibnd1,')=',imapy1(ibnd1),','
          if (inertial_cooling1(ibnd1))
     .     write(99,'(a17,i4,a2,l1,a1)')
     .     'INERTIAL_COOLING(',ibnd1,')=',inertial_cooling1(ibnd1),','
          if (backplate_temp1(ibnd1).ne.plate_temp) 
     .     write(99,'(a15,i4,a2,1es15.8,a1)') 
     .     'BACKPLATE_TEMP(',ibnd1,')=',backplate_temp1(ibnd1),','
          if (plate_thickness1(ibnd1).ne.plate_thick)
     .     write(99,'(a16,i4,a2,1es15.8,a1)') 
     .     'PLATE_THICKNESS(',ibnd1,')=',plate_thickness1(ibnd1),','
          if (coating_thickness1(ibnd1).gt.0.0_R8)
     .     write(99,'(a18,i4,a2,1es15.8,a1)') 
     .     'COATING_THICKNESS(',ibnd1,')=',coating_thickness1(ibnd1),','
          if (plate_time_factor1(ibnd1).ne.1.0_R8)
     .     write(99,'(a18,i4,a2,1es15.8,a1)') 
     .     'PLATE_TIME_FACTOR(',ibnd1,')=',plate_time_factor1(ibnd1),','
          if (trim(surface_material_name1(ibnd1)).ne.'C') 
     .     write(99,'(a22,i4,a4,a6,a2)') 
     .     'SURFACE_MATERIAL_NAME(',ibnd1,')= ''',
     .      surface_material_name1(ibnd1),''','
          if (coating_thickness1(ibnd1).gt.0.0_R8) 
     .     write(99,'(a22,i4,a4,a6,a2)') 
     .     'COATING_MATERIAL_NAME(',ibnd1,')= ''',
     .      coating_material_name1(ibnd1),''','
          if (trim(bulk_material_name1(ibnd1)).ne.'C') 
     .     write(99,'(a19,i4,a4,a6,a2)') 
     .     'BULK_MATERIAL_NAME(',ibnd1,')= ''',
     .      bulk_material_name1(ibnd1),''','
          if (layer_model) then
            write(99,'(a20,i4,a2,i2,a1)')
     .       'LAYER_NCONSTITUENTS(',ibnd1,')=',
     .        layer_nconstituents1(ibnd1),','
            do i = 1, layer_nconstituents1(ibnd1)
              write(99,'(a21,i4,a1,i2,a4,i2,a1)')
     .         'LAYER_NZCONSTITUENTS(',ibnd1,',',i,')=  ',
     .          layer_nzconstituents1(ibnd1,i),','
              write(99,'(a23,i4,a1,i2,a2,1es15.8,a1)')
     .         'LAYER_NRELCONSTITUENTS(',ibnd1,',',i,')=',
     .          layer_nrelconstituents1(ibnd1,i),','
            end do
          end if
          do itrack=1,ntrack
            if (deposition1(ibnd1,itrack).ne.0.0_R8)
     .       write(99,'(a11,i4,a1,i1,a2,1es15.8,a1)')
     .        'DEPOSITION(',ibnd1,',',itrack,')=',
     1        deposition1(ibnd1,itrack),','
            if (monolayer_deposition1(ibnd1,itrack).ne.0.0_R8)
     .       write(99,'(a,i4,a,i1,a,1es15.8,a)') 
     .        'MONOLAYER_DEPOSITION(',ibnd1,',',itrack,')=',
     1         monolayer_deposition1(ibnd1,itrack),','
            if (erosion1(ibnd1,itrack).ne.0.0_R8)
     .       write(99,'(a8,i4,a1,i1,a2,1es15.8,a1)')
     .        'EROSION(',ibnd1,',',itrack,')=',
     1         erosion1(ibnd1,itrack),','
            if (monolayer_erosion1(ibnd1,itrack).ne.0.0_R8)
     .       write(99,'(a,i4,a,i1,a,1es15.8,a)') 
     .        'MONOLAYER_EROSION(',ibnd1,',',itrack,')=',
     1         monolayer_erosion1(ibnd1,itrack),','
          enddo
          if (plate_time1(ibnd1).ne.0.0_R8)
     .     write(99,'(a11,i4,a2,1es15.8,a1)')
     .      'PLATE_TIME(',ibnd1,')=',plate_time1(ibnd1),','
          if (plate_area1(ibnd1).ne.0.0_R8)
     .     write(99,'(a11,i4,a2,1es15.8,a1)')
     .      'PLATE_AREA(',ibnd1,')=',plate_area1(ibnd1),','
        enddo
        if (plate_model.ge.1.and.ndepth1.gt.1) then 
          do ibnd1 = 1, nbnd1
            do i = 1, (ndepth1-1)/nline+1
              if (i.eq.1) then
                write(line,'(a12,i4,a3,i4,a2)') 
     .           'TARGET_TEMP(',ibnd1,',1:',ndepth1,')='
              else
                line = ''
              endif
              do is = (i-1)*nline+1, min(i*nline,ndepth1)
                write(numstr,'(1x,1es15.8,a1)') 
     .           target_temp1(ibnd1,is),','
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

      subroutine b2yt_wall_copy(ns,ns1,nsm,is_neutral1,bundle_frac,
     .  natm,natm1,natmm,nmol,nion,iatm1,nwall,wall_data,
     .  contr_area,plate_area)
      use b2mod_types
      implicit none
      integer ns,ns1,nsm,natm,natm1,natmm,nmol,nion,nwall
      integer iatm1(natm)
      logical is_neutral1(0:ns1-1)
      real (kind=R8) bundle_frac(0:ns1-1,0:ns-1)
      real (kind=R8) contr_area(nwall,nwall)
      real (kind=R8) plate_area(nwall)
      real (kind=R8) wall_data(nwall,0:nsm-1+natmm+nmol+nion)
      integer iw, iw1, iatm
      real (kind=R8) wall_source(0:ns-1+natm+nmol+nion)
      real (kind=R8) wall_map(0:ns1-1+natm1+nmol+nion)
      real (kind=R8) wall_tmp(0:nsm-1)
      real (kind=R8) species_source(nwall)
      real (kind=R8) species_map(nwall)
      
cxpb We do the species mapping first
      do iw = 1, nwall
        wall_map = 0.0_R8
        wall_source(0:ns-1+natm+nmol+nion) = 
     .   wall_data(iw,0:ns-1+natm+nmol+nion)
        wall_tmp(0:ns-1) = wall_source(0:ns-1)
        call average_bundle
     .   (ns1,ns,nsm,is_neutral1,bundle_frac,wall_tmp)
        wall_map(0:ns1-1) = wall_tmp(0:ns1-1)
        do iatm = 1, natm
          if (iatm1(iatm).gt.0 .and. iatm1(iatm).le.natm1) then
            wall_map(ns1-1+iatm1(iatm)) = wall_source(ns-1+iatm)
          end if
        end do
        wall_map(ns1+natm1:ns1-1+natm1+nmol+nion) = 
     .   wall_source(ns+natm:ns-1+natm+nmol+nion)
        wall_data(iw,0:ns1-1+natm1+nmol+nion) =
     .   wall_map(0:ns1-1+natm1+nmol+nion)
      end do
cxpb Now for each new species, we do the space mapping
      do iatm = 0, ns1-1+natm1+nmol+nion
        species_map = 0.0_R8
        species_source(1:nwall) = wall_data(1:nwall,iatm)
        call b2yt_wall_average(nwall,contr_area,plate_area,
     .   species_source,species_map)
        wall_data(1:nwall,iatm) = species_map(1:nwall)
      end do

      return
      end subroutine b2yt_wall_copy

      subroutine b2yt_wall_average(nwall,contr_area,plate_area,
     . wall_source,wall_map)
      use b2mod_types
      implicit none
      integer nwall
      real (kind=R8) contr_area(nwall,nwall)
      real (kind=R8) plate_area(nwall)
      real (kind=R8) wall_source(nwall)
      real (kind=R8) wall_map(nwall)
      integer iw, iw1
      real (kind=R8) sum
      
      wall_map = 0.0_R8
      if (maxval(wall_source(1:nwall)).eq.minval(wall_source(1:nwall)))
     > then
        wall_map(1:nwall) = wall_source(1)
      else
        do iw1 = 1, nwall
          sum = 0.0_R8
          do iw = 1, nwall
            wall_map(iw1) = wall_map(iw1) + wall_source(iw)*
     .       contr_area(iw,iw1)*plate_area(iw)
            sum = sum + contr_area(iw,iw1)*plate_area(iw)
          end do
          if (sum.gt.0.0_R8) wall_map(iw1) = wall_map(iw1)/sum
        end do
      end if 

      return
      end subroutine b2yt_wall_average

      subroutine b2yt_wall_sum(nwall,contr_area,wall_data)
      use b2mod_types
      implicit none
      integer nwall
      real (kind=R8) contr_area(nwall,nwall)
      real (kind=R8) wall_data(nwall)
      real (kind=R8) wall_source(nwall)
      real (kind=R8) wall_map(nwall)
      integer iw, iw1
      
      wall_source(1:nwall) = wall_data(1:nwall)
      wall_map = 0.0_R8
      do iw1 = 1, nwall
        do iw = 1, nwall
          wall_map(iw1) = wall_map(iw1) + wall_source(iw)*
     .     contr_area(iw,iw1)
        end do
      end do
      wall_data(1:nwall) = wall_map(1:nwall)

      return
      end subroutine b2yt_wall_sum

      subroutine b2yt_wall_bundle_sum(ns,ns1,nsm,is_neutral1,
     .  bundle_frac,natm,natm1,natmm,nmol,nion,iatm1,
     .  nwall,wall_data,contr_area)
      use b2mod_types
      implicit none
      integer ns,ns1,nsm,natm,natm1,natmm,nmol,nion,nwall
      integer iatm1(natm)
      logical is_neutral1(0:ns1-1)
      real (kind=R8) bundle_frac(0:ns1-1,0:ns-1)
      real (kind=R8) contr_area(nwall,nwall)
      real (kind=R8) wall_data(nwall,0:nsm-1+natmm+nmol+nion)
      integer iw, iw1, iatm
      real (kind=R8) wall_source(0:ns-1+natm+nmol+nion)
      real (kind=R8) wall_map(0:ns1-1+natm1+nmol+nion)
      real (kind=R8) wall_tmp(0:nsm-1)
      real (kind=R8) species_source(nwall)

cxpb We do the species mapping first
      do iw = 1, nwall
        wall_map = 0.0_R8
        wall_source(0:ns-1+natm+nmol+nion) = 
     .   wall_data(iw,0:ns-1+natm+nmol+nion)
        wall_tmp(0:ns-1) = wall_source(0:ns-1)
        call sum_bundle
     .   (ns1,ns,nsm,is_neutral1,bundle_frac,wall_tmp)
        wall_map(0:ns1-1) = wall_tmp(0:ns1-1)
        do iatm = 1, natm
          if (iatm1(iatm).gt.0 .and. iatm1(iatm).le.natm1) then
            wall_map(ns1-1+iatm1(iatm)) = wall_source(ns-1+iatm)
          end if
        end do
        wall_map(ns1+natm1:ns1-1+natm1+nmol+nion) = 
     .   wall_source(ns+natm:ns-1+natm+nmol+nion)
        wall_data(iw,0:ns1-1+natm1+nmol+nion) =
     .   wall_map(0:ns1-1+natm1+nmol+nion)
      end do
cxpb Now for each new species, we do the space mapping
      do iatm = 0, ns1-1+natm1+nmol+nion
        species_source(1:nwall) = wall_data(1:nwall,iatm)
        call b2yt_wall_sum(nwall,contr_area,species_source)
        wall_data(1:nwall,iatm) = species_source(1:nwall)
      end do

      return
      end subroutine b2yt_wall_bundle_sum
      
