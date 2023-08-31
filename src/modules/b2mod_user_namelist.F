      module b2mod_user_namelist
      use b2mod_types
      use b2mod_diag
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
      real (kind=R8) ::
     ,  fusion_power, Spmp_He_to_D, Spmp_nom
      integer lpfrb_i,lpfrb_o,lpfrt_i,lpfrt_o,
     ,  j_he_at,j_h_at(3),j_ne_at,lpfrs_pmp,
     ,  npfrgrp,ipfrgrp(npfrgrpd),jpfrgrp(npfrgrpd),
     ,  nntrgrp,intrgrp(nntrgrpd),jntrgrp(nntrgrpd)
      real (kind=R8) :: Te_det_threshold
      character*8 gpfrgrp(npfrgrpd), gntrgrp(nntrgrpd)
      character*80 user_filename
      logical scale_He, print_isat
      logical user_namelist_used, written
      logical, save :: write_nml_user = .true.
      logical, save :: user_initialised = .false.

      private :: user
      namelist /user/
     ,  lhetrgts,lpfrb_i,lpfrb_o,lpfrt_i,lpfrt_o,
     ,  j_h_at,j_he_at,j_ne_at,
     ,  l_h_mol,fusion_power,Spmp_He_to_D,lpfrs_pmp,
     ,  npfrgrp,lpfrgrp,ipfrgrp,jpfrgrp,gpfrgrp,
     ,  nntrgrp,lntrgrp,intrgrp,jntrgrp,gntrgrp,
     ,  Spmp_nom,scale_He,print_isat,write_nml_user,
     ,  filedata,user_filename,Te_det_threshold
c----------------------------------------------------------------------
c*** fusion_power : fusion power (including neutrons) [MW]
c*** Spmp_He_to_D : ratio of He to DT pumping speeds (typically, 0.8)
c*** Te_det_threshold : Temperature threshold at the target below which
c***                    a flux tube is considered detached (typically, 2.0 eV)
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

      subroutine alloc_b2mod_user(nx,ns,nlim,nmol,ntns)
      implicit none
      integer nx, ns, nlim, nmol, ntns
      external xertst

      call xertst(0.le.nx, 'faulty argument nx')
      call xertst(1.le.ns, 'faulty argument ns')
      call xertst(0.le.nlim, 'faulty argument nlim')
      call xertst(0.le.nmol, 'faulty argument nmol')
      call xertst(0.le.ntns, 'faulty argument ntns')
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
      Te_det_threshold = 2.0_R8
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
      user_filename='b2.user.parameters'
      call read_b2mod_user_namelist(nx,ns)
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

      subroutine read_b2mod_user_namelist(nx,ns)
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
      integer nx, ns

      logical :: file_ok
      character*260 filename
      integer i, j, k, is, im
      logical old_style, done
      integer, save :: use_eirene, neutrals_namelist
      logical, save :: hydrogen, helium, neon
      external ipgeti, find_file, xertst

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
      call xertst (abs(lpfrb_i).le.nx, 'faulty input lpfrb_i')
      call xertst (abs(lpfrb_o).le.nx, 'faulty input lpfrb_o')
      call xertst (abs(lpfrt_i).le.nx, 'faulty input lpfrt_i')
      call xertst (abs(lpfrt_o).le.nx, 'faulty input lpfrt_o')
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
      if (filedata(ntrc_uss)) call xertst(nlim.le.DEF_NLIM,
     . 'insufficient size of DEF_NLIM in b2mod_dimensions.F !')

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
