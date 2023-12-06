      module b2mod_user_namelist
      use b2mod_types
      use b2mod_diag
      use b2us_geo
      use b2us_map
      implicit none
#     include "DIMENSIONS.F"
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
! - icsep[oi]mp: location of the separatrix position in the [oi]mp list
! - ifsep[oi]mp: the face index of the separatrix intersection with the [oi]mp
! - n[oi]mp: number of CVs in the [oi]mp list
! - rz[o1]mp: R-Z coordinates of a segment whose intersection with the
!             plasma grid defines the OMP and IMP. rzomp(1,1:2) holds R-coordinates
!             of first and second point of the segment, rzomp(2,1:2) holds Z-coordinates
      integer, parameter :: nromp=DEF_NYD*4
      integer, save :: omp(nromp), icsepomp, ifsepomp, nomp,
     ,                 imp(nromp), icsepimp, ifsepimp, nimp
      real (kind=R8) :: rzomp(2,2), rzimp(2,2)
      real (kind=R8) ::
     ,  fusion_power, Spmp_He_to_D, Spmp_nom
      integer lpfrb_i,lpfrb_o,lpfrt_i,lpfrt_o,
     ,  j_he_at,j_h_at(3),j_ne_at,lpfrs_pmp,
     ,  npfrgrp,ipfrgrp(npfrgrpd),jpfrgrp(npfrgrpd),
     ,  nntrgrp,intrgrp(nntrgrpd),jntrgrp(nntrgrpd)
      character*8 gpfrgrp(npfrgrpd), gntrgrp(nntrgrpd)
      character*80 user_filename
      logical scale_He, print_isat
      logical, save :: user_namelist_used = .true.
      logical, save :: written = .false.
      logical, save :: write_nml_user = .false.
      logical, save :: user_initialised = .false.
      private :: user
      namelist /user/
     ,  lhetrgts,lpfrb_i,lpfrb_o,lpfrt_i,lpfrt_o,
     ,  j_h_at,j_he_at,j_ne_at,
     ,  l_h_mol,fusion_power,Spmp_He_to_D,lpfrs_pmp,
     ,  npfrgrp,lpfrgrp,ipfrgrp,jpfrgrp,gpfrgrp,
     ,  nntrgrp,lntrgrp,intrgrp,jntrgrp,gntrgrp,
     ,  Spmp_nom,scale_He,print_isat,write_nml_user,
     ,  filedata,user_filename,rzomp,rzimp
c----------------------------------------------------------------------
c*** fusion_power : fusion power (including neutrons) [MW]
c*** Spmp_He_to_D : ratio of He to DT pumping speeds (typically, 0.8)
c*** lhetrgts  :  list of surface indices (EIRENE notation) which
c***              are used for calculation of helium enrichment
c*** lpfrb_[io]:  B2 x-cell indices corresponding to the bottom of the
c***              bypass between the inner and outer divertor (counted
c***              from the bottom, if non-positive)
c*** lpfrt_[io]:  B2 x-cell indices corresponding to the top of the
c***              bypass between the inner and outer divertor (counted
c***              from the x-point, if non-positive)
c*** j_h_at    :  Hydrogen isotope index in EIRENE atoms list
c*** j_he_at   :  Helium index in EIRENE atoms list
c*** j_ne_at   :  Neon index in EIRENE atoms list
c*** l_h_mol   :  Hydrogen isotope content in the molecules
c*** lpfrs_pmp :  Location of the pump.
c***              0 no pump at all, 1 - lower PFR, 2 - lower outer
c*** npfrgrpd  :  max. number of surface groups for PFR flows
c*** npfrgrp   :  actual number of surface groups for PFR flows
c*** lpfrgrp   :  list of surface segments for PFR flows
c*** [ij]pfrgrp:  first and last positions in this list for each group
c*** gpfrgrp   :  identifiers for the PFR surface groups
c*** nntrgrpd  :  max. number of surface groups for neutral data
c*** nntrgrp   :  actual number of surface groups for neutral data
c*** lntrgrp   :  list of surface segments for neutral data
c*** [ij]ntrgrp:  first and last positions in this list for each group
c*** gntrgrp   :  identifiers for the neutral data surface groups
c*** filedata  :  flag to turn on/off individual tracing files
c***              (see b2mod_diag for their definition)
c----------------------------------------------------------------------

      contains

      subroutine alloc_b2mod_user(geo,m,ns,nlim,nmol)
      implicit none
      integer, intent(in) :: ns, nlim, nmol
      type (geometry), intent (in) :: geo
      type(mapping), intent(inout) :: m

      call xertst(1.le.ns, 'faulty argument ns')
      call xertst(0.le.nlim, 'faulty argument nlim')
      call xertst(0.le.nmol, 'faulty argument nmol')
      if (user_initialised) return
#ifdef ALLOCATE_AND_NAMELIST
      nhetrgtx=nlim
      allocate(lhetrgts(nhetrgtx),lpfrgrp(nlim),lntrgrp(nlim))
      allocate(l_h_mol(nmol,3))
#endif

c*** PFR specification - for He enrichment and neutral pressure
c***  Standard set-up:
      lpfrs_pmp=0            ! no pumping
      fusion_power=0.0_R8    ! no helium
      Spmp_He_to_D=1.0_R8    ! same pumping speed for He and D
      Spmp_nom=0.0_R8
      j_h_at =0
      l_h_mol=0
      j_he_at=0
      j_ne_at=0
      lpfrb_i=0
      lpfrb_o=0
      lpfrt_i=0
      lpfrt_o=0
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
      user_filename='b2.user.parameters'
      call read_b2mod_user_namelist(ns,geo,m)
      user_initialised = .true.

      return
      end subroutine alloc_b2mod_user

      subroutine dealloc_b2mod_user
      implicit none
#ifdef ALLOCATE_AND_NAMELIST
      deallocate(lhetrgts,lpfrgrp,lntrgrp)
      deallocate(l_h_mol)
#endif
      user_initialised = .false.
      return
      end subroutine dealloc_b2mod_user

      subroutine read_b2mod_user_namelist(ns,geo,m)
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
      integer, intent(in) :: ns
      type (geometry), intent (in) :: geo
      type(mapping), intent(inout) :: m
      integer ic, iFc

      logical :: file_ok
      character*260 filename
      integer i
      logical done, ifail, found
      integer, save :: use_eirene, neutrals_namelist
      integer j, k, is, im
      logical old_style
      logical, save :: hydrogen, helium, neon
      external ipgeti, find_file, xertst
#ifndef F2003
      integer newunit
      external newunit
#endif

      call subini ('readuser')
      call ipgeti ('b2mndr_eirene', use_eirene)
      call ipgeti ('b2stbr_neutrals_namelist', neutrals_namelist)

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

      nomp = 0
      omp = 0 !Need to add check for monotonicity in list
! Define OMP list by finding CVs which are intersected by the rzomp segment
      call cv_intersections (geo, m, rzomp, omp, nomp, nromp,
     & 1, ifail)
      call xertst(.not.ifail,'Increase nromp in b2.user.parameters')
! Here we need to guarantee ordering of CVs. With classical grids the conversion
! is done such that CVs are ordered from left to right and bottom to top wrt
! the structured computational domain. With fully unstructured grids this is not guaranteed.
! IDEA: project the cell center onto the segment and order CVs based on distance wrt one of the segment's extremes
! SAME FOR IMP
      do ic=2,nomp
        if (geo%cvX(omp(ic)).lt.geo%cvX(omp(ic-1)))
     &   call xerrab('OMP LIST NOT SORTED')
      end do
*   ..Definition of separatrix position in OMP list
! Separatrix is defined by checking when cvReg changes value inside the OMP list
      icsepomp = 0
      ic = 1
      found = .false.
      do while (ic.le.nomp .and..not.found)
        do i = m%cvFcP(omp(ic),1),
     &         m%cvFcP(omp(ic),1) + m%cvFcP(omp(ic),2) - 1
          iFc = m%cvFc(i)
          if (m%fcFs(iFc).eq.m%iFssep.and.m%iFssep.ne.0) then
            icsepomp = ic + 1
            ifsepomp = iFc
            found = .true.
          end if
        end do
        ic = ic + 1
      end do

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
#ifdef WG_TODO
      call xertst (abs(lpfrb_i).le.nx, 'faulty input lpfrb_i')
      call xertst (abs(lpfrb_o).le.nx, 'faulty input lpfrb_o')
      call xertst (abs(lpfrt_i).le.nx, 'faulty input lpfrt_i')
      call xertst (abs(lpfrt_o).le.nx, 'faulty input lpfrt_o')
#endif
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
      if (neutrals_namelist.eq.1 .or. use_eirene.ne.0) then
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
      call cv_intersections (geo, m, rzimp, imp, nimp, nromp,
     & 1, ifail)
      call xertst(.not.ifail,'Increase nromp in b2.user.parameters')
      do ic=2,nimp
        if (geo%cVx(imp(ic)).lt.geo%cVx(imp(ic-1)))
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
          if (m%fcFs(iFc).eq.m%iFssep.and.m%iFssep.ne.0) then
            icsepimp = ic + 1
            ifsepimp = iFc
            found = .true.
          end if
        end do
        ic = ic + 1
      end do
      if (.not.found) write(*,'(a/a/a)')
     &  'Warning: did not find separatrix location in IMP definition,',
     &  '         check specification of RZIMP in b2.user.parameters',
     &  '         if IMP is needed.'
      if (m%IsClassicalGrid.eq.1) call xertst(m%fsFcP(m%iFssep,2).gt.0,
     & 'No separatrix faces!')
      write(*,*) 'SEPARATRIX FLUX SURFACE AT iFs=',m%iFssep

      call subend ()
      return
#ifdef USE_MPI
  999 call xerrab ('Error while reading USER namelist!')
#endif
      end subroutine read_b2mod_user_namelist

      subroutine write_b2mod_user_namelist
      implicit none

      if (written) return
      if (user_namelist_used.and.write_nml_user) write(*,user)
      written = .true.
      return
      end subroutine write_b2mod_user_namelist

      end module b2mod_user_namelist

