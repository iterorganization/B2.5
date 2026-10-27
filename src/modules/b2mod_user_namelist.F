      module b2mod_user_namelist
      use b2mod_types
      use b2mod_diag
      use b2mod_dimensions
      use b2us_geo
      use b2us_map
      use b2mod_switches
      use b2mod_dimensions
      implicit none
      integer npfrgrpd, nntrgrpd
      parameter (npfrgrpd=6, nntrgrpd=16)
      integer nhetrgtx
#ifdef ALLOCATE_AND_NAMELIST
      integer, allocatable :: lhetrgts(:),lpfrgrp(:),lntrgrp(:)
      integer, allocatable :: l_h_mol(:,:)
#else
      parameter (nhetrgtx=DEF_NLIM)
      integer lhetrgts(nhetrgtx),lpfrgrp(DEF_NLIM),lntrgrp(DEF_NLIM)
      integer l_h_mol(DEF_NMOL,3)
#endif
! VARIABLES RELATED TO DEFINITION OF MIDPLANES AND SEPARATRIX
! - [oi]mp: list of CVs belonging to the [OI]MP, defined through the segment rz[oi]mp
!           it should be built as an ordered list from core to SOL.
! - nromp: pre-defined maximum size for these lists
! - icsep[oi]mp: location of the first cell outside the separatrix in the [oi]mp list
! - ifsep[oi]mp: the face index of the separatrix intersection with the [oi]mp
! - n[oi]mp: number of CVs in the [oi]mp list
! - rz[o1]mp: R-Z coordinates of a segment whose intersection with the
!             plasma grid defines the OMP and IMP. rzomp(1,1:2) holds R-coordinates
!             of first and second point of the segment, rzomp(2,1:2) holds Z-coordinates
! - psi_[oi]mp: psi value of midplane cells
! - fs_ds_[oi]mp: distance along midplane of intersection with flux surface
! - ft_ds_[oi]mp: distance along midplane of flux tube cell center
! - ft_[oi]mp: index in list of flux tube cell intersection

      integer, parameter :: nromp=DEF_NYD*4
      integer, parameter :: ntarget_fclbl_max=32
      integer, parameter :: ntarget_face_group_list_max=4096
      integer, parameter :: target_fclbl_groups_str_len=4096
      integer, save :: omp(nromp), icsepomp, ifsepomp, nomp,
     ,                 imp(nromp), icsepimp, ifsepimp, nimp
      integer, save :: ntarget_fclbl,
     ,                 target_fclbl_list(ntarget_fclbl_max)
      integer, save :: ntarget_face_group
      integer, save :: target_face_group_sizes(ntarget_fclbl_max)
      integer, save :: target_face_group_list
      dimension target_face_group_list(ntarget_face_group_list_max)
      character*32, save :: target_fclbl_name_list(ntarget_fclbl_max)
      character*32, save :: target_fclbl_group_name_list
      dimension target_fclbl_group_name_list(ntarget_fclbl_max)
      character*32, save :: target_face_group_name_list
      dimension target_face_group_name_list(ntarget_fclbl_max)
      character*(target_fclbl_groups_str_len), save ::
     ,  target_fclbl_groups_str
      integer, allocatable :: ft_imp(:), ft_omp(:)
      real (kind=R8) :: rzomp(2,2), rzimp(2,2)
      real (kind=R8) :: psi_omp(nromp), psi_imp(nromp)
      real (kind=R8), save :: omp_pos(nromp), imp_pos(nromp)
      real (kind=R8), allocatable :: fs_ds_imp(:), fs_ds_omp(:)
      real (kind=R8), allocatable :: ft_ds_imp(:), ft_ds_omp(:)
      real (kind=R8) ::
     ,  fusion_power, Spmp_He_to_D, Spmp_nom
      integer :: lpfrb_i,lpfrb_o,lpfrt_i,lpfrt_o,
     ,  j_he_at,j_h_at(3),j_ne_at,lpfrs_pmp,
     ,  npfr_cvs,pfr_cvs(100),
     ,  npfrgrp,ipfrgrp(npfrgrpd),jpfrgrp(npfrgrpd),
     ,  nntrgrp,intrgrp(nntrgrpd),jntrgrp(nntrgrpd)
      real (kind=R8) :: Te_det_threshold, q95_albll,
     ,  R0_albll, B0_albll
      real (kind=R8) :: trgshp(DEF_NCUT*2)
      character*8 gpfrgrp(npfrgrpd), gntrgrp(nntrgrpd)
      character*80 user_filename
      logical scale_He, print_isat
      logical, save :: user_namelist_used = .true.
      logical, save :: written = .false.
      logical, save :: write_nml_user = .false.
      logical, save :: user_initialised = .false.

      namelist /user/
     ,  lhetrgts,lpfrb_i,lpfrb_o,lpfrt_i,lpfrt_o,
     ,  j_h_at,j_he_at,j_ne_at,
     ,  l_h_mol,fusion_power,Spmp_He_to_D,lpfrs_pmp,
     ,  npfr_cvs,pfr_cvs,
     ,  npfrgrp,lpfrgrp,ipfrgrp,jpfrgrp,gpfrgrp,
     ,  nntrgrp,lntrgrp,intrgrp,jntrgrp,gntrgrp,
     ,  Spmp_nom,scale_He,print_isat,write_nml_user,
     ,  filedata,user_filename,rzomp,rzimp,
     ,  ntarget_fclbl,target_fclbl_list,
     ,  target_fclbl_name_list,
     ,  target_fclbl_groups_str,
     ,  target_fclbl_group_name_list,
     ,  ntarget_face_group,
     ,  target_face_group_sizes,
     ,  target_face_group_list,
     ,  target_face_group_name_list,
     ,  Te_det_threshold,q95_albll,R0_albll,B0_albll,trgshp
c----------------------------------------------------------------------
c*** fusion_power : fusion power (including neutrons) [MW]
c*** Spmp_He_to_D : ratio of He to DT pumping speeds (typically, 0.8)
c*** Te_det_threshold : Temperature threshold at the target below which
c***                    a flux tube is considered detached (typically, 2.0 eV)
c*** q95_albll    : safety factor used to estimate the ballooning "alpha"
c***                parameter, if not assigned will be estimated as
c***                Lcon/(2*pi*R) at the core interface
c*** R0_albll     : major radius used to estimate the ballooning "alpha"
c***                parameter, if not assigned will be estimated as
c***                0.5*(Rmin-Rmax) at the core interface
c*** B0_albll     : major radius used to estimate the ballooning "alpha"
c***                parameter, if not assigned will be estimated as
c***                Bomp*R0/Romp at the core interface
c*** trgshp       : target shaping factor applied to plasma heat loads in
c***                peak heat flux estimate calculated by b2mod_usrtrc
c*** lhetrgts     : list of surface indices (EIRENE notation) which
c***                are used for calculation of helium enrichment
c*** j_h_at       : Hydrogen isotope index in EIRENE atoms list
c*** j_he_at      : Helium index in EIRENE atoms list
c*** j_ne_at      : Neon index in EIRENE atoms list
c*** l_h_mol      : Hydrogen isotope content in the molecules
c*** lpfrs_pmp    : Location of the pump.
c***                0 - no pump at all, 1 - lower PFR, 2 - lower outer
c*** lpfrb_[io]   : Obsolete. Use PFR_CVS instead.
c*** lpfrt_[io]   : Obsolete. Use PFR_CVS instead.
c*** npfr_cvs     : number of control volumes used for PFR neutral pressure
c***                calculation
c*** pfr_cvs      : list of those control volumes
c*** npfrgrpd     : max. number of surface groups for PFR flows
c*** npfrgrp      : actual number of surface groups for PFR flows
c*** lpfrgrp      : list of surface segments for PFR flows
c*** [ij]pfrgrp   : first and last positions in this list for each group
c*** gpfrgrp      : identifiers for the PFR surface groups
c*** nntrgrpd     : max. number of surface groups for neutral data
c*** nntrgrp      : actual number of surface groups for neutral data
c*** lntrgrp      : list of surface segments for neutral data
c*** [ij]ntrgrp   : first and last positions in this list for each group
c*** gntrgrp      : identifiers for the neutral data surface groups
c*** filedata     : flag to turn on/off individual tracing files
c***                (see b2mod_diag for their definition)
c----------------------------------------------------------------------

      contains

      subroutine alloc_b2mod_user(geo,m,ns,nlim,nmol,sw)
      implicit none
      integer, intent(in) :: ns, nlim, nmol
      type (geometry), intent (in) :: geo
      type(mapping), intent(inout) :: m
      type (switches), intent (inout) :: sw

      call xertst(1.le.ns, 'faulty argument ns')
      call xertst(0.le.nlim, 'faulty argument nlim')
      call xertst(0.le.nmol, 'faulty argument nmol')
      if (user_initialised) return
#ifdef ALLOCATE_AND_NAMELIST
      nhetrgtx=nlim
      allocate(lhetrgts(nhetrgtx),lpfrgrp(nlim),lntrgrp(nlim))
      allocate(l_h_mol(nmol,3))
#endif
      allocate(fs_ds_imp(m%nFs), fs_ds_omp(m%nFs))
      allocate(ft_ds_imp(m%nFt), ft_ds_omp(m%nFt))
      allocate(ft_imp(0:m%nFt), ft_omp(0:m%nFt))

c*** PFR specification - for He enrichment and neutral pressure
c***  Standard set-up:
      lpfrs_pmp=0            ! no pumping
      fusion_power=0.0_R8    ! no helium
      Spmp_He_to_D=1.0_R8    ! same pumping speed for He and D
      Spmp_nom=0.0_R8
      Te_det_threshold = 2.0_R8
      q95_albll = 0.0_R8
      B0_albll = 0.0_R8
      R0_albll = 0.0_R8
      trgshp = 1.0_R8
      j_h_at =0
      l_h_mol=0
      j_he_at=0
      j_ne_at=0
      lpfrb_i=0
      lpfrb_o=0
      lpfrt_i=0
      lpfrt_o=0
      npfr_cvs=0
      pfr_cvs=0
      npfrgrp=0
      ipfrgrp=0
      jpfrgrp=0
      nntrgrp=0
      intrgrp=0
      jntrgrp=0
      lhetrgts(1:nhetrgtx)=0
      lhetrgts(1:nsts)=-ltns(1:nsts)
      gpfrgrp=' '
      gntrgrp=' '
      scale_He=.false.
      print_isat=.false.
      rzomp = 0.0_R8
      rzimp = 0.0_R8
      ft_imp = 0
      ft_omp = 0
      fs_ds_imp = -999.0_R8
      fs_ds_omp = -999.0_R8
      ft_ds_imp = -999.0_R8
      ft_ds_omp = -999.0_R8
      ntarget_fclbl = 0
      target_fclbl_list = 0
      target_fclbl_name_list = ' '
      target_fclbl_groups_str = ' '
      target_fclbl_group_name_list = ' '
      ntarget_face_group = 0
      target_face_group_sizes = 0
      target_face_group_list = 0
      target_face_group_name_list = ' '
      user_filename='b2.user.parameters'
      call read_b2mod_user_namelist(ns,geo,m,sw)
      user_initialised = .true.

      return
      end subroutine alloc_b2mod_user

      subroutine dealloc_b2mod_user
      implicit none

      if (.not.user_initialised) return
#ifdef ALLOCATE_AND_NAMELIST
      deallocate(lhetrgts,lpfrgrp,lntrgrp)
      deallocate(l_h_mol)
#endif
      deallocate(fs_ds_omp,fs_ds_imp)
      deallocate(ft_ds_omp,ft_ds_imp)
      deallocate(ft_omp,ft_imp)
      user_initialised = .false.
      return
      end subroutine dealloc_b2mod_user

      subroutine read_b2mod_user_namelist(ns,geo,m,sw)
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
      integer, intent(in) :: ns
      type (geometry), intent (in) :: geo
      type(mapping), intent(inout) :: m
      type (switches), intent (inout) :: sw
      integer ic, iFc, iFs, iFs1, iFs2, iFt, iVx1, iVx2
      real (kind=R8) :: distint, dsep,
     . Xv1, Xv2, Yv1, Yv2, Xm1, Xm2, Ym1, Ym2, Xint, Yint

      logical :: file_ok, on_wall, wall_face(2), general_topology,
     &           reference_topology
      character*260 filename
      integer i
      logical done, ifail, found
      integer j, k, is, im
      logical old_style
      logical, save :: hydrogen, helium, neon
      external ipgeti, find_file, xertst

      call subini ('readuser')

*  -- USER --
      filename=user_filename
      call find_file(filename,file_ok)
      if(file_ok) then
        open(99,file=filename)
        write(*,*) 'Reading namelist user from ',trim(filename)
#ifdef USE_MPI
        read(99,user,err=999)
#else
        read(99,user)
#endif
        close(99)
        write(*,*) 'Read namelist user from ',trim(filename)
        user_namelist_used = .true.
      else
        write(*,*) 'No ',trim(filename)
        user_namelist_used = .false.
        filedata(ntrc_usr) = .false.
      endif
      call build_target_face_groups_from_fclbl_string(m)

! Explicit general, linear, and cylindrical topologies use one selected
! reference flux surface at both midplanes. Only GENERAL needs the new
! arbitrary-topology distance construction below.
      reference_topology = m%hasTopologicalData.and.
     &  (m%geometryID.eq.1.or.m%geometryID.eq.2.or.
     &   m%geometryID.eq.13)
      general_topology = m%hasTopologicalData.and.m%geometryID.eq.13

      nomp = 0
      omp = 0
! Define OMP list by finding CVs which are intersected by the rzomp segment;
! cv_intersections orders them along the segment (chord midpoints, running in
! increasing R). omp_pos holds these positions, oriented like the list; they
! are the abscissa of the OMP profiles (b2time.nc dsa, input profiles).
      call cv_intersections (geo, m, rzomp, omp, nomp, nromp, 1, ifail)
      call xertst(.not.ifail, 'Increase DEF_NYD in b2mod_dimensions')
      omp_pos = 0.0_R8
      if (nomp.gt.0) then
        call cv_segment_positions (geo, m, rzomp, omp, nomp,
     &                             omp_pos(1:nomp))
        if (rzomp(1,2).lt.rzomp(1,1)) omp_pos(1:nomp) = -omp_pos(1:nomp)
      end if
      do ic=2,nomp
        if (omp_pos(ic).lt.omp_pos(ic-1))
     &   call xerrab('OMP LIST NOT SORTED')
      end do
*   ..Definition of separatrix position in OMP list
      icsepomp = 0
      ic = 1
      found = .false.
      do while (ic.le.nomp .and..not.found)
        do i = m%cvFcP(omp(ic),1),
     &         m%cvFcP(omp(ic),1) + m%cvFcP(omp(ic),2) - 1
          iFc = m%cvFc(i)
          if ((m%fcFs(iFc).eq.m%iFssep.and.m%iFssep.ne.0).or.
     &        (.not.reference_topology.and.
     &         m%fcFs(iFc).eq.m%iFssep2.and.m%iFssep2.ne.0)) then
            icsepomp = ic + 1
            ifsepomp = iFc
            found = .true.
          end if
        end do
        ic = ic + 1
      end do
      if (.not.found.and.reference_topology) then
        call xerrab('read_b2mod_user_namelist: no OMP reference ' //
     &    'surface position')
      endif
! csc get psi for each OMP cell
      psi_omp = 0.0_R8
      do ic = 1, nomp
        psi_omp(ic) = geo%cvFpsi(omp(ic))
      end do
! xpb build the distance along the OMP for each flux surface
      if (general_topology) then
        call general_midplane_distances(m, geo, rzomp, omp, nomp,
     &    ifsepomp, .true., 'OMP', fs_ds_omp, ft_omp, ft_ds_omp)
      else
      Xm1 = rzomp(1,1)
      Xm2 = rzomp(1,2)
      Ym1 = rzomp(2,1)
      Ym2 = rzomp(2,2)
      do ic = 1, nomp
        do i = m%cvFcP(omp(ic),1),
     &         m%cvFcP(omp(ic),1) + m%cvFcP(omp(ic),2) - 1
          iFc = m%cvFc(i)
          iFs = m%fcFs(iFc)
          if (iFs.eq.0) cycle
          ! Does the flux-aligned face intersect the OMP?
          iVx1 = m%fcVx(iFc,1)
          iVx2 = m%fcVx(iFc,2)
          ! Find the intersection point of the two lines
          Xv1 = geo%vxX(iVx1)
          Xv2 = geo%vxX(iVx2)
          Yv1 = geo%vxY(iVx1)
          Yv2 = geo%vxY(iVx2)
          Yint =
     .     ((Xm1*Ym2)-(Xm2*Ym1)-
     .      (((Xv1*Yv2)-(Xv2*Yv1))/(Yv2-Yv1))*(Ym2-Ym1))/
     .     (((Xv1-Xv2)/(Yv2-Yv1))*(Ym2-Ym1)+(Xm1-Xm2))
          Xint = (((Xv1*YV2)-(Xv2*Yv1))-Yint*(Xv1-Xv2))/(Yv2-Yv1)
          ! Is the intersection point part of the cell face?
          distint = (Xint-Xv1)*(Xint-Xv2)+(Yint-Yv1)*(Yint-Yv2)
          if (distint.le.0.0_R8) then
            fs_ds_omp(iFs) = sqrt((Xint-Xm1)**2+(Yint-Ym1)**2)
          end if
        end do
      end do
      dsep = fs_ds_omp(m%fcFs(ifsepomp))
      do iFs = 1, m%nFs
        if (fs_ds_omp(iFs).eq.-999.0_R8) cycle
        if (Xm1.lt.Xm2) then
          fs_ds_omp(iFs) = fs_ds_omp(iFs) - dsep
        else
          fs_ds_omp(iFs) = dsep - fs_ds_omp(iFs)
        end if
      end do
      do ic = 1, nomp
        if (omp(ic).gt.m%nCi) cycle ! do not do guard cells
        iFs1 = 0
        iFs2 = 0
        iFt = m%cvFt(omp(ic))
        if (iFt.eq.0) cycle
        ft_omp(iFt) = ic
        wall_face = .false.
        do i = m%cvFcP(omp(ic),1),
     &         m%cvFcP(omp(ic),1) + m%cvFcP(omp(ic),2) - 1
          iFc = m%cvFc(i)
          iFs = m%fcFs(iFc)
          on_wall = m%fcLbl(iFc).ne.0
          if (on_wall) then
            if (iFs1.eq.0) wall_face(1) = .true.
            if (iFs1.ne.0 .and. iFs2.eq.0) wall_face(2) = .true.
          end if
          if (iFs1.eq.0 .and. iFs.ne.0) then
            if (wall_face(1)) then
              iFs2 = iFs
            else
              iFs1 = iFs
            end if
          else if (iFs2.eq.0 .and. iFs.ne.iFs1 .and. iFs.ne.0) then
            iFs2 = iFs
          end if
        end do
        if (iFs1.ne.0.and.iFs2.ne.0) then
          ft_ds_omp(iFt) = 0.5_R8*(fs_ds_omp(iFs1)+fs_ds_omp(iFs2))
        else if (wall_face(1)) then
          ft_ds_omp(iFt) = fs_ds_omp(iFs2)+0.5_R8*geo%cvHy(omp(ic))
        else
          ft_ds_omp(iFt) = fs_ds_omp(iFs1)+0.5_R8*geo%cvHy(omp(ic))
        end if
      end do
      endif
cxpb set default switches depending on midplane location
      if (sw%set_transport_iyref.eq.0) then
        sw%set_transport_iyref = icsepomp-1
      end if
      call xertst (1.le.sw%set_transport_iyref.and.
     .             sw%set_transport_iyref.le.nomp,
     .            'Invalid switch setting for set_transport_iyref')

cxpb convert old-style j_h_at and l_h_mol to deal with isotope mixes
      old_style = maxval(j_h_at).eq.0
      if (j_h_at(1).ne.0.and.j_h_at(2).eq.0.and.j_h_at(3).eq.0) then
        old_style = nint(am(eb2atcr(j_h_at(1)))).ne.1
      end if
c*** Set default values according to current species list
      hydrogen = .false.
      done = .false.
      is = 0
      do while (is.lt.ns)
        hydrogen = hydrogen .or.
     &   (is_neutral(is) .and. nint(zn(is)).eq.1)
        if (is_neutral(is).and.nint(zn(is)).eq.1) then
          if (old_style .and. nint(am(is)).ne.1 .and. .not.done) then
            if (j_h_at(1).ne.0) then
              i = nint(am(j_h_at(1)))
              j_h_at(i) = j_h_at(1)
              j_h_at(1) = 0
            else
              i = nint(am(is))
            end if
            do im = 1, nnmoli
              l_h_mol(im,i) = l_h_mol(im,1)
              l_h_mol(im,1) = 0
            end do
            done = .true.
          end if
          if (j_h_at(nint(am(is))).eq.0)
     >        j_h_at(nint(am(is))) = b2eatcr(is)
        end if
        is = is + 1
      end do
      helium = .false.
      is = 0
      do while (.not.helium .and. is.lt.ns)
        helium = is_neutral(is) .and. nint(zn(is)).eq.2
        if (helium .and. j_he_at.eq.0) j_he_at = b2eatcr(is)
        is = is + 1
      end do
      neon = .false.
      is = 0
      do while (.not.neon .and. is.lt.ns)
        neon = is_neutral(is) .and. nint(zn(is)).eq.10
        if (neon .and. j_ne_at.eq.0) j_ne_at = b2eatcr(is)
        is = is + 1
      end do

*   ..test user namelist input
      call xertst(0.0_R8.le.fusion_power,
     . 'faulty input fusion_power')
      call xertst(0.0_R8.le.Spmp_nom,
     . 'faulty input Spmp_nom')
      call xertst(0.0_R8.le.Spmp_He_to_D,
     . 'faulty input Spmp_He_to_D')
      call xertst(0.0_R8.lt.Te_det_threshold,
     . 'faulty input Te_det_threshold')
      call xertst(0.0_R8.le.q95_albll,
     . 'faulty input q95_albll')
      call xertst(0.0_R8.le.B0_albll,
     . 'faulty input B0_albll')
      call xertst(0.0_R8.le.R0_albll,
     . 'faulty input R0_albll')
      call xertst(1.0_R8.le.minval(trgshp),
     . 'faulty input trgshp')
      call xertst(0.le.npfr_cvs, 'faulty input npfr_cvs')
      if (npfr_cvs.gt.0)
     > call xertst(minval(pfr_cvs(1:npfr_cvs)).gt.0.and.
     &             maxval(pfr_cvs(1:npfr_cvs)).le.m%nCi,
     &            'faulty input pfr_cvs')
      call xertst (0.le.npfrgrp, 'faulty input npfrgrp')
      do i = 1, npfrgrp
        call xertst(1.le.ipfrgrp(i).and.ipfrgrp(i).le.jpfrgrp(i),
     .   'faulty inputs ipfrgrp, jpfrgrp')
      end do
      do i = 2, npfrgrp
        call xertst(jpfrgrp(i-1).lt.ipfrgrp(i),
     .   'faulty inputs ipfrgrp, jpfrgrp')
      end do
      call xertst (0.le.nntrgrp,
     . 'faulty input nntrgrp')
      do i = 1, nntrgrp
        call xertst(1.le.intrgrp(i).and.intrgrp(i).le.jntrgrp(i),
     .   'faulty inputs intrgrp, jntrgrp')
      end do
      do i = 2, nntrgrp
        call xertst(jntrgrp(i-1).lt.intrgrp(i),
     .   'faulty inputs intrgrp, jntrgrp')
      end do
      call xertst(.not.(hydrogen.and.maxval(j_h_at).eq.0),
     .   'faulty input j_h_at')
      do i = 1,3
        if (j_h_at(i).eq.0) cycle
        if (hydrogen) then
          call xertst(nint(zn(eb2atcr(j_h_at(i)))).eq.1,
     .     'faulty input j_h_at')
        end if
      end do
      do i = 1,2
        do j = i+1,3
          if (j_h_at(i).ne.0 .and. j_h_at(j).ne.0) then
            if (nint(am(eb2atcr(j_h_at(i)))).eq.
     .          nint(am(eb2atcr(j_h_at(j))))) then
              write(6,*)
     .         'Non-standard setting of hydrogen isotope species:'
              write(6,'(1x,a,i3)')
     .         'More than one hydrogen species defined with mass ',
     .          nint(am(eb2atcr(j_h_at(i))))
              write(6,'(1x,a,3i3)') 'j_h_at : ', (j_h_at(k),k=1,3)
              write(6,'(1x,a,3i3)') 'masses : ',
     .         (nint(am(eb2atcr(j_h_at(k)))),k=1,3)
            end if
          end if
        end do
      end do
      if (sw%b2stbr_neutrals_namelist.eq.1 .or. sw%use_eirene.ne.0) then
        do i = 1,3
          if (j_h_at(i).eq.0) cycle
          if (user_namelist_used) then
            if (maxval(l_h_mol).gt.0) then
              do is = 1, nnmoli
                call xertst(l_h_mol(is,i).eq.mlcmp(j_h_at(i),is),
     .           'faulty inputs l_h_mol, mlcmp')
              end do
            else
              do is = 1, nnmoli
                l_h_mol(is,i) = mlcmp(j_h_at(i),is)
              end do
            end if
          end if
        end do
      endif
      if (helium) call xertst(nint(zn(eb2atcr(j_he_at))).eq.2,
     . 'faulty input j_he_at')
      if (neon) call xertst(nint(zn(eb2atcr(j_ne_at))).eq.10,
     . 'faulty input j_ne_at')

*   ..Create list of IMP control volumes
! Same logic as OMP list, find CVs intersected by rzimp segment
      nimp = 0
      imp = 0
      call cv_intersections (geo, m, rzimp, imp, nimp, nromp, 1, ifail)
      call xertst(.not.ifail, 'Increase DEF_NYD in b2mod_dimensions')
      imp_pos = 0.0_R8
      if (nimp.gt.0) then
        call cv_segment_positions (geo, m, rzimp, imp, nimp,
     &                             imp_pos(1:nimp))
        if (rzimp(1,2).lt.rzimp(1,1)) imp_pos(1:nimp) = -imp_pos(1:nimp)
      end if
      do ic=2,nimp
        if (imp_pos(ic).lt.imp_pos(ic-1))
     &   call xerrab('IMP LIST NOT SORTED')
      end do
! Find position of separatrix in IMP list (same as for OMP)
      icsepimp = 0
      ic = 1
      found = .false.
      do while (ic.le.nimp .and..not.found)
        do i = m%cvFcP(imp(ic),1),
     &         m%cvFcP(imp(ic),1) + m%cvFcP(imp(ic),2) - 1
          iFc = m%cvFc(i)
          if ((m%fcFs(iFc).eq.m%iFssep.and.m%iFssep.ne.0).or.
     &        (.not.reference_topology.and.
     &         m%fcFs(iFc).eq.m%iFssep2.and.m%iFssep2.ne.0)) then
            icsepimp = ic + 1
            ifsepimp = iFc
            found = .true.
          end if
        end do
        ic = ic + 1
      end do
! csc get psi for each OMP cell
      psi_imp = 0.0_R8
      do ic = 1, nimp
        psi_imp(ic) = geo%cvFpsi(imp(ic))
      end do
! xpb build the distance along the IMP for each flux surface
      if (general_topology) then
        if (found) call general_midplane_distances(m, geo, rzimp,
     &    imp, nimp, ifsepimp, .false., 'IMP', fs_ds_imp,
     &    ft_imp, ft_ds_imp)
      else
      if (nimp.gt.0) then
        Xm1 = rzimp(1,1)
        Xm2 = rzimp(1,2)
        Ym1 = rzimp(2,1)
        Ym2 = rzimp(2,2)
      end if
      do ic = 1, nimp
        do i = m%cvFcP(imp(ic),1),
     &         m%cvFcP(imp(ic),1) + m%cvFcP(imp(ic),2) - 1
          iFc = m%cvFc(i)
          iFs = m%fcFs(iFc)
          if (iFs.eq.0) cycle
          ! Does the flux-aligned face intersect the IMP?
          iVx1 = m%fcVx(iFc,1)
          iVx2 = m%fcVx(iFc,2)
          ! Find the intersection point of the two lines
          Xv1 = geo%vxX(iVx1)
          Xv2 = geo%vxX(iVx2)
          Yv1 = geo%vxY(iVx1)
          Yv2 = geo%vxY(iVx2)
          Yint =
     .     ((Xm1*Ym2)-(Xm2*Ym1)-
     .      (((Xv1*Yv2)-(Xv2*Yv1))/(Yv2-Yv1))*(Ym2-Ym1))/
     .     (((Xv1-Xv2)/(Yv2-Yv1))*(Ym2-Ym1)+(Xm1-Xm2))
          Xint = (((Xv1*YV2)-(Xv2*Yv1))-Yint*(Xv1-Xv2))/(Yv2-Yv1)
          ! Is the intersection point part of the cell face?
          distint = (Xint-Xv1)*(Xint-Xv2)+(Yint-Yv1)*(Yint-Yv2)
          if (distint.le.0.0_R8) then
            fs_ds_imp(iFs) = sqrt((Xint-Xm1)**2+(Yint-Ym1)**2)
          end if
        end do
      end do
      if (ifsepimp.gt.0) dsep = fs_ds_imp(m%fcFs(ifsepimp))
      do iFs = 1, m%nFs
        if (fs_ds_imp(iFs).eq.-999.0_R8) cycle
        if (Xm1.gt.Xm2) then
          fs_ds_imp(iFs) = fs_ds_imp(iFs) - dsep
        else
          fs_ds_imp(iFs) = dsep - fs_ds_imp(iFs)
        end if
      end do
      do ic = 1, nimp
        if (imp(ic).gt.m%nCi) cycle ! do not do guard cells
        iFs1 = 0
        iFs2 = 0
        iFt = m%cvFt(imp(ic))
        if (iFt.eq.0) cycle
        ft_imp(iFt) = ic
        wall_face = .false.
        do i = m%cvFcP(imp(ic),1),
     &         m%cvFcP(imp(ic),1) + m%cvFcP(imp(ic),2) - 1
          iFc = m%cvFc(i)
          iFs = m%fcFs(iFc)
          on_wall = m%fcLbl(iFc).ne.0
          if (on_wall) then
            if (iFs1.eq.0) wall_face(1) = .true.
            if (iFs1.ne.0 .and. iFs2.eq.0) wall_face(2) = .true.
          end if
          if (iFs1.eq.0 .and. iFs.ne.0) then
            if (wall_face(1)) then
              iFs2 = iFs
            else
              iFs1 = iFs
            end if
          else if (iFs2.eq.0 .and. iFs.ne.iFs1 .and. iFs.ne.0) then
            iFs2 = iFs
          end if
        end do
        if (iFs1.ne.0.and.iFs2.ne.0) then
          ft_ds_imp(iFt) = 0.5_R8*(fs_ds_imp(iFs1)+fs_ds_imp(iFs2))
        else if (wall_face(1)) then
          ft_ds_imp(iFt) = fs_ds_imp(iFs2)+0.5_R8*geo%cvHy(imp(ic))
        else
          ft_ds_imp(iFt) = fs_ds_imp(iFs1)+0.5_R8*geo%cvHy(imp(ic))
        end if
      end do
      endif
      if (.not.found) then
        if (reference_topology) then
          write(*,'(a)') 'Warning: no IMP reference location'
        else
          write(*,'(a)') 'Warning: no IMP separatrix location'
        endif
        write(*,'(a)') '         check RZIMP in b2.user.parameters'
        write(*,'(a)') '         if IMP is needed.'
      endif
      if (m%isClassicalGrid.eq.1) call xertst(m%fsFcP(m%iFssep,2).gt.0,
     & 'No separatrix faces!')
      if (reference_topology) then
        write(*,*) 'REFERENCE FLUX SURFACE AT iFs=',m%iFssep
      else
        write(*,*) 'SEPARATRIX FLUX SURFACE AT iFs=',m%iFssep
      endif
      if (.not.reference_topology.and.m%iFssep2.gt.0) then
        write(*,*) 'SECOND HALF OF CDN SEPARATRIX FLUX SURFACE AT iFs=',
     &   m%iFssep2
      end if

      call subend ()
      return
#ifdef USE_MPI
  999 call xerrab ('Error while reading USER namelist!')
#endif
      end subroutine read_b2mod_user_namelist

**********************************************************************
!> Compute midplane flux-surface and flux-tube distances for an
!! explicit general topology. Flux-surface crossings come from the
!! cell faces, while the two bounding surfaces of a cell come from its
!! vertices. Only boundary faces of that same cell establish wall
!! contact.
      subroutine general_midplane_distances(m, geo, rz, midplane,
     &  nMidplane, referenceFace, outwardIncreasingR, midplaneName,
     &  fsDistance, ftIndex, ftDistance)
      implicit none

      type(mapping), intent(in) :: m
      type(geometry), intent(in) :: geo
      integer, intent(in) :: nMidplane, referenceFace
      integer, intent(in) :: midplane(nMidplane)
      logical, intent(in) :: outwardIncreasingR
      character*(*), intent(in) :: midplaneName
      real(kind=R8), intent(in) :: rz(2,2)
      real(kind=R8), intent(out) :: fsDistance(:), ftDistance(:)
      integer, intent(out) :: ftIndex(0:)

      integer :: ic, i, iCv, iFc, iFs, iFs1, iFs2, iFt, iVx
      real(kind=R8) :: distance, distanceTolerance, dsep,
     &  segmentLength
      logical :: intersects, collinear, onWall
      external xerrab

      fsDistance = -999.0_R8
      ftDistance = -999.0_R8
      ftIndex = 0

!   Validate the reference and the midplane segment before calculating
!   any intersections. A zero-length segment has no distance axis.
      if (referenceFace.lt.1.or.referenceFace.gt.m%nFc) then
        write(*,*) 'general_midplane_distances: ', midplaneName,
     &    ' reference face, nFc = ', referenceFace, m%nFc
        call xerrab('general_midplane_distances: reference face ' //
     &    'out of range')
      endif
      iFs = m%fcFs(referenceFace)
      if (iFs.ne.m%iFssep) then
        write(*,*) 'general_midplane_distances: ', midplaneName,
     &    ' reference face/surface/selected = ', referenceFace,
     &    iFs, m%iFssep
        call xerrab('general_midplane_distances: reference face ' //
     &    'does not lie on the selected surface')
      endif
      segmentLength = sqrt((rz(1,2)-rz(1,1))**2+
     &  (rz(2,2)-rz(2,1))**2)
      if (segmentLength.le.100.0_R8*epsilon(1.0_R8)*
     &    max(1.0_R8,maxval(abs(rz)))) then
        write(*,*) 'general_midplane_distances: degenerate ',
     &    midplaneName, ' segment ', rz
        call xerrab('general_midplane_distances: zero-length ' //
     &    'midplane segment')
      endif
      distanceTolerance = 1.0e-10_R8*max(1.0_R8,segmentLength)

!   Intersect every flux-aligned face of the midplane cells with the
!   user segment. The cross-product form supports vertical faces and
!   vertical midplane segments without a coordinate denominator.
      do ic = 1, nMidplane
        iCv = midplane(ic)
        do i = m%cvFcP(iCv,1),
     &         m%cvFcP(iCv,1) + m%cvFcP(iCv,2) - 1
          iFc = m%cvFc(i)
          iFs = m%fcFs(iFc)
          if (iFs.eq.0) cycle
          if (iFs.lt.1.or.iFs.gt.size(fsDistance)) then
            write(*,*) 'general_midplane_distances: ', midplaneName,
     &        ' face/surface/range = ', iFc, iFs, size(fsDistance)
            call xerrab('general_midplane_distances: face surface ' //
     &        'ID out of range')
          endif
          call midplane_face_intersection(geo, m, rz, iFc,
     &      intersects, collinear, distance)
          if (collinear) then
            write(*,*) 'general_midplane_distances: ', midplaneName,
     &        ' segment is collinear with face ', iFc
            call xerrab('general_midplane_distances: ambiguous ' //
     &        'collinear face intersection')
          endif
          if (.not.intersects) cycle
          if (fsDistance(iFs).eq.-999.0_R8) then
            fsDistance(iFs) = distance
          else if (abs(fsDistance(iFs)-distance).gt.
     &             distanceTolerance) then
            write(*,*) 'general_midplane_distances: ', midplaneName,
     &        ' surface has two crossings ', iFs,
     &        fsDistance(iFs), distance
            call xerrab('general_midplane_distances: ambiguous ' //
     &        'flux-surface intersection')
          endif
        enddo
      enddo

      if (fsDistance(m%iFssep).eq.-999.0_R8) then
        write(*,*) 'general_midplane_distances: ', midplaneName,
     &    ' selected surface did not cross segment ', m%iFssep
        call xerrab('general_midplane_distances: reference surface ' //
     &    'does not cross midplane segment')
      endif
      dsep = fsDistance(m%iFssep)
      do iFs = 1, m%nFs
        if (fsDistance(iFs).eq.-999.0_R8) cycle
        if ((outwardIncreasingR.and.rz(1,1).lt.rz(1,2)).or.
     &      (.not.outwardIncreasingR.and.rz(1,1).gt.rz(1,2))) then
          fsDistance(iFs) = fsDistance(iFs) - dsep
        else
          fsDistance(iFs) = dsep - fsDistance(iFs)
        endif
      enddo

!   Derive each intersected flux-tube center from the distinct surface
!   IDs on the current cell only. Wall contact is also a property of
!   the current cell; a boundary face belonging only to a neighboring
!   cell must not affect this result.
      do ic = 1, nMidplane
        iCv = midplane(ic)
        if (iCv.gt.m%nCi) cycle
        iFt = m%cvFt(iCv)
        if (iFt.eq.0.or.iFt.eq.m%nFt+1) cycle
        if (iFt.lt.1.or.iFt.gt.m%nFt) then
          write(*,*) 'general_midplane_distances: ', midplaneName,
     &      ' cell/tube/range = ', iCv, iFt, m%nFt
          call xerrab('general_midplane_distances: flux-tube ID ' //
     &      'out of range')
        endif

        iFs1 = 0
        iFs2 = 0
        do i = m%cvVxP(iCv,1),
     &         m%cvVxP(iCv,1) + m%cvVxP(iCv,2) - 1
          iVx = m%cvVx(i)
          iFs = m%vxFs(iVx)
          if (iFs.eq.0) cycle
          if (iFs.lt.1.or.iFs.gt.size(fsDistance)) then
            write(*,*) 'general_midplane_distances: ', midplaneName,
     &        ' cell/vertex/surface/range = ', iCv, iVx, iFs,
     &        size(fsDistance)
            call xerrab('general_midplane_distances: vertex surface ' //
     &        'ID out of range')
          endif
          if (iFs1.eq.0) then
            iFs1 = iFs
          else if (iFs.ne.iFs1.and.iFs2.eq.0) then
            iFs2 = iFs
          else if (iFs.ne.iFs1.and.iFs.ne.iFs2) then
            write(*,*) 'general_midplane_distances: ', midplaneName,
     &        ' cell has more than two surfaces ', iCv,
     &        iFs1, iFs2, iFs
            call xerrab('general_midplane_distances: more than two ' //
     &        'flux surfaces bound a midplane cell')
          endif
        enddo

        onWall = .false.
        do i = m%cvFcP(iCv,1),
     &         m%cvFcP(iCv,1) + m%cvFcP(iCv,2) - 1
          iFc = m%cvFc(i)
          onWall = onWall.or.m%fcLbl(iFc).ne.0
        enddo

        if (iFs1.ne.0.and.iFs2.ne.0) then
          if (fsDistance(iFs1).ne.-999.0_R8.and.
     &        fsDistance(iFs2).ne.-999.0_R8) then
            ftIndex(iFt) = ic
            ftDistance(iFt) = 0.5_R8*(fsDistance(iFs1)+
     &        fsDistance(iFs2))
          endif
        else if (iFs1.ne.0.and.onWall) then
          if (fsDistance(iFs1).ne.-999.0_R8) then
            ftIndex(iFt) = ic
            ftDistance(iFt) = fsDistance(iFs1)+0.5_R8*geo%cvHy(iCv)
          endif
        endif
      enddo

      return
      end subroutine general_midplane_distances

**********************************************************************
!> Return the distance from the first midplane endpoint to the unique
!! intersection with one face. Cross products avoid divisions by a
!! particular R or Z component.
      subroutine midplane_face_intersection(geo, m, rz, iFc,
     &  intersects, collinear, distance)
      implicit none

      type(mapping), intent(in) :: m
      type(geometry), intent(in) :: geo
      real(kind=R8), intent(in) :: rz(2,2)
      integer, intent(in) :: iFc
      logical, intent(out) :: intersects, collinear
      real(kind=R8), intent(out) :: distance

      integer :: iVx1, iVx2
      real(kind=R8) :: dmx, dmy, dvx, dvy, qx, qy,
     &  denom, denomTolerance, faceLength, offsetCross,
     &  offsetTolerance, parameterTolerance, segmentLength,
     &  tSegment, tFace
      external xerrab

      iVx1 = m%fcVx(iFc,1)
      iVx2 = m%fcVx(iFc,2)
      dmx = rz(1,2) - rz(1,1)
      dmy = rz(2,2) - rz(2,1)
      dvx = geo%vxX(iVx2) - geo%vxX(iVx1)
      dvy = geo%vxY(iVx2) - geo%vxY(iVx1)
      qx = geo%vxX(iVx1) - rz(1,1)
      qy = geo%vxY(iVx1) - rz(2,1)
      segmentLength = sqrt(dmx**2+dmy**2)
      faceLength = sqrt(dvx**2+dvy**2)
      if (faceLength.le.100.0_R8*epsilon(1.0_R8)*
     &    max(1.0_R8,abs(geo%vxX(iVx1)),
     &        abs(geo%vxY(iVx1)))) then
        write(*,*) 'midplane_face_intersection: degenerate face ',
     &    iFc, iVx1, iVx2
        call xerrab('midplane_face_intersection: zero-length face')
      endif

      denom = dmx*dvy - dmy*dvx
      denomTolerance = 100.0_R8*epsilon(1.0_R8)*
     &  segmentLength*faceLength
      offsetCross = qx*dmy - qy*dmx
      offsetTolerance = 100.0_R8*epsilon(1.0_R8)*
     &  segmentLength*max(1.0_R8,sqrt(qx**2+qy**2))
      intersects = .false.
      collinear = .false.
      distance = 0.0_R8
      if (abs(denom).le.denomTolerance) then
        collinear = abs(offsetCross).le.offsetTolerance
        return
      endif

      tSegment = (qx*dvy-qy*dvx)/denom
      tFace = offsetCross/denom
      parameterTolerance = 1000.0_R8*epsilon(1.0_R8)
      if (tSegment.ge.-parameterTolerance.and.
     &    tSegment.le.1.0_R8+parameterTolerance.and.
     &    tFace.ge.-parameterTolerance.and.
     &    tFace.le.1.0_R8+parameterTolerance) then
        intersects = .true.
        tSegment = max(0.0_R8,min(1.0_R8,tSegment))
        distance = tSegment*segmentLength
      endif

      return
      end subroutine midplane_face_intersection
      subroutine build_target_face_groups_from_fclbl_string(m)
      implicit none
      type(mapping), intent(in) :: m
      integer :: labels_flat(ntarget_face_group_list_max)
      integer :: sizes(ntarget_fclbl_max)
      integer :: starts(ntarget_fclbl_max)
      integer :: ngroups, nlabels
      integer :: g, ibc, iFc, lbl, nfaces, istart, iend, i
      logical :: duplicate

      if (len_trim(target_fclbl_groups_str).le.0) return
      if (.not.allocated(m%bcFc)) call xerrab(
     . 'target_fclbl_groups_str: boundary faces are not allocated')

      call parse_tfg_str(target_fclbl_groups_str,
     . ngroups,sizes,labels_flat,nlabels)

      if (ngroups.le.0) call xerrab(
     . 'target_fclbl_groups_str: no groups parsed')
      if (ngroups.gt.ntarget_fclbl_max) call xerrab(
     . 'target_fclbl_groups_str: too many groups')

      target_face_group_sizes = 0
      target_face_group_list = 0
      target_face_group_name_list = ' '

      ntarget_face_group = ngroups
      ntarget_fclbl = 0
      starts = 0
      istart = 1
      do g = 1, ngroups
        starts(g) = istart
        istart = istart + sizes(g)
      end do

      nfaces = 0
      do g = 1, ngroups
        istart = starts(g)
        iend = istart + sizes(g) - 1
        do ibc = 1, size(m%bcFc)
          iFc = m%bcFc(ibc)
          if (iFc.le.0) cycle
          lbl = m%fcLbl(iFc)
          if (lbl.le.0) cycle
          if (.not.int_in_list(lbl,labels_flat(istart:iend),
     .        sizes(g))) cycle
          duplicate = .false.
          if (target_face_group_sizes(g).gt.0) then
            duplicate = int_in_list(iFc,
     .        target_face_group_list(
     .          nfaces-target_face_group_sizes(g)+1:nfaces),
     .        target_face_group_sizes(g))
          end if
          if (duplicate) cycle
          if (nfaces+1.gt.ntarget_face_group_list_max) call xerrab(
     .     'target_fclbl_groups_str: too many grouped target faces')
          nfaces = nfaces + 1
          target_face_group_list(nfaces) = iFc
          target_face_group_sizes(g) = target_face_group_sizes(g) + 1
        end do
      end do

      do i = 1, ntarget_face_group
        if (len_trim(target_fclbl_group_name_list(i)).gt.0) then
          target_face_group_name_list(i) =
     .      target_fclbl_group_name_list(i)
        else
          target_face_group_name_list(i) = ' '
        end if
      end do

      write(*,*) 'Parsed target_fclbl_groups_str into',
     . ntarget_face_group,' groups'

      return
      end subroutine build_target_face_groups_from_fclbl_string

      subroutine parse_tfg_str(str,ngroups,sizes,labels,nlabels)
      implicit none
      character*(*) str
      integer, intent(out) :: ngroups, nlabels
      integer, intent(out) :: sizes(ntarget_fclbl_max)
      integer, intent(out) :: labels(ntarget_face_group_list_max)
      integer :: pos, slen, val, gsize
      logical :: close_group

      ngroups = 0
      nlabels = 0
      sizes = 0
      labels = 0
      slen = len_trim(str)
      pos = 1

      do while (pos.le.slen)
        call skip_group_separators(str,pos,slen)
        if (pos.gt.slen) exit
        if (str(pos:pos).ne.'[') call xerrab(
     .   'target_fclbl_groups_str: expected "["')
        ngroups = ngroups + 1
        if (ngroups.gt.ntarget_fclbl_max) call xerrab(
     .   'target_fclbl_groups_str: too many groups')
        pos = pos + 1
        gsize = 0
        close_group = .false.
        do while (.not.close_group)
          call skip_group_spaces(str,pos,slen)
          if (pos.gt.slen) call xerrab(
     .     'target_fclbl_groups_str: missing "]"')
          if (str(pos:pos).eq.']') then
            pos = pos + 1
            close_group = .true.
          else
            call read_group_integer(str,pos,slen,val)
            nlabels = nlabels + 1
            if (nlabels.gt.ntarget_face_group_list_max) call xerrab(
     .       'target_fclbl_groups_str: too many label entries')
            labels(nlabels) = val
            gsize = gsize + 1
            call skip_group_spaces(str,pos,slen)
            if (pos.le.slen .and. str(pos:pos).eq.',') pos = pos + 1
          end if
        end do
        if (gsize.le.0) call xerrab(
     .   'target_fclbl_groups_str: empty group is not allowed')
        sizes(ngroups) = gsize
      end do

      return
      end subroutine parse_tfg_str

      subroutine skip_group_spaces(str,pos,slen)
      implicit none
      character*(*) str
      integer, intent(inout) :: pos
      integer, intent(in) :: slen
      do while (pos.le.slen)
        if (str(pos:pos).eq.' ' .or. str(pos:pos).eq.char(9)) then
          pos = pos + 1
        else
          exit
        end if
      end do
      return
      end subroutine skip_group_spaces

      subroutine skip_group_separators(str,pos,slen)
      implicit none
      character*(*) str
      integer, intent(inout) :: pos
      integer, intent(in) :: slen
      do while (pos.le.slen)
        if (str(pos:pos).eq.' ' .or. str(pos:pos).eq.char(9) .or.
     .      str(pos:pos).eq.',') then
          pos = pos + 1
        else
          exit
        end if
      end do
      return
      end subroutine skip_group_separators

      subroutine read_group_integer(str,pos,slen,val)
      implicit none
      character*(*) str
      integer, intent(inout) :: pos
      integer, intent(in) :: slen
      integer, intent(out) :: val
      integer :: istart, iend, ios
      character*32 token

      call skip_group_spaces(str,pos,slen)
      if (pos.gt.slen) call xerrab(
     . 'target_fclbl_groups_str: expected integer value')
      istart = pos
      if (str(pos:pos).eq.'+' .or. str(pos:pos).eq.'-') pos = pos + 1
      do while (pos.le.slen)
        if (str(pos:pos).ge.'0' .and. str(pos:pos).le.'9') then
          pos = pos + 1
        else
          exit
        end if
      end do
      iend = pos - 1
      if (iend.lt.istart) call xerrab(
     . 'target_fclbl_groups_str: invalid integer token')
      token = ' '
      if (iend-istart+1.gt.len(token)) call xerrab(
     . 'target_fclbl_groups_str: integer token too long')
      token(1:iend-istart+1) = str(istart:iend)
      read(token,*,iostat=ios) val
      if (ios.ne.0) call xerrab(
     . 'target_fclbl_groups_str: failed integer parse')
      return
      end subroutine read_group_integer

      logical function int_in_list(val,list,n)
      implicit none
      integer, intent(in) :: val, n
      integer, intent(in) :: list(n)
      int_in_list = .false.
      if (n.le.0) return
      int_in_list = any(list(1:n).eq.val)
      return
      end function int_in_list

      subroutine write_b2mod_user_namelist
      implicit none

      if (written) return
      if (user_namelist_used.and.write_nml_user) write(*,user)
      written = .true.
      return
      end subroutine write_b2mod_user_namelist

      end module b2mod_user_namelist
