      module b2mod_user_namelist
      use b2mod_types
      use b2mod_diag
      use b2us_geo
      use b2us_map
      use b2mod_ad
     & , only : nncf, ffdata
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
      integer, parameter :: nromp=DEF_NYD*4
      integer, save :: omp(nromp), icsepomp, nomp, ifssep, imp(nromp),
     ,  icsepimp, nimp
      real (kind=R8) ::
     ,  fusion_power, Spmp_He_to_D, Spmp_nom, rzomp(2,2), rzimp(2,2)
      integer lpfrb_i,lpfrb_o,lpfrt_i,lpfrt_o,
     ,  j_he_at,j_h_at(3),j_ne_at,lpfrs_pmp,
     ,  npfrgrp,ipfrgrp(npfrgrpd),jpfrgrp(npfrgrpd),
     ,  nntrgrp,intrgrp(nntrgrpd),jntrgrp(nntrgrpd),
     ,  cfstart(nncf),cfend(nncf),cftype(nncf)
      character*8 gpfrgrp(npfrgrpd), gntrgrp(nntrgrpd)
      character*80 user_filename
      logical scale_He, print_isat
      logical, save :: user_namelist_used = .true.
      logical, save :: written = .false.
      logical, save :: write_nml_user = .false.
      logical, save :: user_initialised = .false.
      logical, save :: cfread(nncf), cfomp(nncf)
      integer, save :: ncfdata(nncf), ncf
      real (kind=R8), save :: cfdata(nncf,3,DEF_NLIM), cfweight(nncf)

      private :: user
      namelist /user/
     ,  lhetrgts,lpfrb_i,lpfrb_o,lpfrt_i,lpfrt_o,
     ,  j_h_at,j_he_at,j_ne_at,
     ,  l_h_mol,fusion_power,Spmp_He_to_D,lpfrs_pmp,
     ,  npfrgrp,lpfrgrp,ipfrgrp,jpfrgrp,gpfrgrp,
     ,  nntrgrp,lntrgrp,intrgrp,jntrgrp,gntrgrp,
     ,  Spmp_nom,scale_He,print_isat,write_nml_user,
     ,  filedata,user_filename,ncf,cfstart,cfend,
     ,  cfread, cftype, rzomp, rzimp, cfomp, cfweight
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

      subroutine alloc_b2mod_user(geo,m,ns,nlim,nmol,ntns)
      implicit none
      integer, intent(in) :: ns, nlim, nmol, ntns
      type (geometry), intent (in) :: geo
      type(mapping), intent(inout) :: m

      if (user_initialised) return
#ifdef ALLOCATE_AND_NAMELIST
      nhetrgtx=nlim
      allocate(lhetrgts(nhetrgtx),lpfrgrp(nlim),lntrgrp(nlim))
      allocate(l_h_mol(nmol,3))
#endif
      
c*** PFR specification - for He enrichment and neutral pressure
c***  Standard setup:
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
      user_filename='b2.user.parameters'
      ncf = 0
      cfstart = 0
      cfend = 0
      cfread = .false.
      cfomp = .false.
      cftype = 0
      cfdata = 0.0_R8
      ncfdata = 0
      rzomp = 0.0_R8
      rzimp = 0.0_R8
      cfweight = 1.0_R8
      call read_b2mod_user_namelist(ns,ntns,geo,m)
      user_initialised = .true.

      return
      end subroutine alloc_b2mod_user

      subroutine dealloc_b2mod_user
      implicit none
#ifdef ALLOCATE_AND_NAMELIST
      deallocate(lhetrgts,lpfrgrp,lntrgrp)
      deallocate(l_h_mol)
#endif
      if (allocated(ffdata)) deallocate(ffdata)
      user_initialised = .false.
      return
      end subroutine dealloc_b2mod_user

      subroutine read_b2mod_user_namelist(ns,ntns,geo,m)
      use b2mod_b2cmpa
      implicit none
      integer, intent(in) :: ns,ntns
      type (geometry), intent (in) :: geo
      type(mapping), intent(inout) :: m
      integer ic, noss, iss, indss,
     >        ncffc, incf, idb, icf,
     >        dx, ic1, ic2, ifc, ifcc, ifc1, ifc2,
     >        ii, ifc11, ifc22, jj

      logical :: file_ok
      character*260 filename, cffile
      character*1 str
      integer i, j, k, is, im, ndata
      logical old_style, done, ifail, found
      logical, save :: hydrogen, helium, neon
      external subini, subend

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

      if (ncf > 0) then
        write (*,*) ' Number of cost functions: ',ncf
        m%nCf = ncf
        m%mxnCf = count(m%fcLbl /= 0)
        call alloc_mapping_cf (m)
        m%cfFcP = 0

        open (newunit=idb,file='debug_user.out')        
        incf = 0
        do icf = 1, ncf
          if (cftype(icf).gt.5) call xerrab ('cftype>5 not coded ')
          if (cftype(1).eq.0 .and. ncf.eq.1) call xerrab(
     &     'cftype(1)=0 but no other cost functions are defined')
          if (cftype(icf).eq.0.and. icf.ne.1) call xerrab (
     &     'only cftype(1) can be zero')
          if ((cftype(icf).ge.1 .and. cftype(icf).le.3) .and.
     &         .not.cfread(icf)) call xerrab ('cftype=1,2,3 need to'//
     &         ' read experimental profile')
          if (cfread(icf)) then
            write (str,'(I1)') icf
            cffile = 'cf'//trim(str)//'.dat'
            call find_file(cffile,file_ok)
            if (.not. file_ok) call xerrab 
     &       ('Cannot find cost function file '//trim(cffile))
            open(unit=99,file=cffile,status='old',action='read')
            write(*,*) 'Reading cost function data from ',trim(cffile)
            ncfdata(icf) = 0
!csc File is of unknown length, so first count the number of lines
 5          read(99,*,end=15) str 
              ncfdata(icf) = ncfdata(icf) +1
            goto 5
 15         close(99)
            call xertst (ncfdata(icf).le.DEF_NLIM,
     &        'Cost function data limited to DEF_NLIM!')
            open(unit=99,file=cffile,status='old',action='read')
            do i=1,ncfdata(icf)
!csc The data in the file are r, z, data. Stored in cfdata according to
!    first dimension (# of cost functions), second dimension (r-z-data), third
!    dimension (number of records/data points)
              read(99,*) cfdata(icf,1:3,i) 
            enddo
            close(99)
          endif
          if (cfstart(icf) /= 0 .and. cfend(icf) /= 0 .and.
     &         cftype(icf).ne.0) then
            ncffc = 0
            m%cfFcP(icf,1) = incf +1 
            noss = cfend(icf) - cfstart(icf) + 1
            do iss = 1, noss
              indss = cfstart(icf)+iss-1
              write (*,*) ' Defining cost function ',icf
              write (*,*) ' surface index ',indss
! find all faces belonging to surface structure INDSS
              call find_faces (indss, incf, ncffc, 
     >             m%mxncf,m%cfFc, m%cfFcOr, m, idb)
            end do              ! iss	 
            m%cfFcP(icf,2) = ncffc
          else
            if (cftype(icf).eq.0 .and. cfstart(icf).eq.0 .or.
     &          cfend(icf).eq.0) call xerrab ('cfstart/cfend = 0 ')
            if (cftype(icf).ne.0 .and. .not.cfomp(icf)) 
     &        call xerrab ('cfstart/cfend = 0 ')
            m%cfFcP(icf,1) = 0
            m%cfFcP(icf,2) = 0
          endif
        end do
        close(idb)
      end if

      call cf_to_struct(m)

! now allocate vectors that will be used to interpolate data
      if (any(cfread)) then
        ndata = maxval(m%cfCvP(1:ncf,2),1)
        allocate(ffdata(ncf,ndata))
        ffdata = 0.0_R8
      endif

#ifdef WG_TODO
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
#ifdef B25_EIRENE
      do i = 1,3
        if (j_h_at(i).eq.0) cycle
        if (user_namelist_used) then
          if (maxval(l_h_mol).gt.0) then
            do is = 1, nnmoli
              call xertst(l_h_mol(is,i).eq.mlcmp(j_h_at(i),is),
     .         'faulty inputs l_h_mol, mlcmp')
            end do
          else
            do is = 1, nnmoli
              l_h_mol(is,i) = mlcmp(j_h_at(i),is)
            end do
          end if
        end if
      end do
#endif
      if (helium) call xertst(nint(zn(eb2atcr(j_he_at))).eq.2,
     . 'faulty input j_he_at')
      if (neon) call xertst(nint(zn(eb2atcr(j_ne_at))).eq.10,
     . 'faulty input j_ne_at')
#endif
      call xertst (ncf.le.nncf, 'increase nncf in b2mod_user_namelist')
*   ..Create list of OMP and IMP control volumes
      nomp = 0
      omp = 0
      call cv_intersections (m%nCv, geo, m, rzomp, omp, nomp, nromp,
     & 1, ifail)
      call xertst(.not.ifail,'Increase nromp in b2.user.parameters')
      do ic=2,nomp
        if (geo%cVx(omp(ic)).lt.geo%cVx(omp(ic-1)))
     & call xerrab('OMP LIST NOT SORTED')
      end do
      nimp = 0
      imp = 0
      call cv_intersections (m%nCv, geo, m, rzimp, imp, nimp, nromp,
     & 1, ifail)
      call xertst(.not.ifail,'Increase nromp in b2.user.parameters')
      do ic=2,nimp
        if (geo%cVx(imp(ic)).lt.geo%cVx(imp(ic-1)))
     & call xerrab('IMP LIST NOT SORTED')
      end do
      if (any(cfomp)) call xertst(nomp.gt.0,
     & 'Asking for OMP cost function but no OMP cells, check rzomp')
*   ..Redefine cost functions on midplane list
      do icf = 1, ncf
        if (cfomp(icf)) then
          write (*,*) ' Re-defining cost function ',icf,' on omp'
          dx = nomp - m%cfCvP(icf,2) !shift
          m%cfCvP(icf,2) = nomp
          if (icf + 1 .le. ncf) then
            ic1 = m%cfCvP(icf+1,1) ! starting CV in next CF
            ic2 = m%cfCvP(ncf,1) + m%cfCvP(ncf,2) - 1 ! last CV in last CF
*         ..Shift CV list
            call xertst(ic2+dx.le.m%mxnCf,
     &       ' Shifting of cfCv list not possible, increase mxnCf')
            m%cfCv(ic1+dx:ic2+dx) = m%cfCv(ic1:ic2)
*         ..Shift cfCvP(1,:) list
            do ic = icf + 1, ncf
              m%cfCvP(ic,1) = m%cfCvP(ic,1) + dx
            enddo
          endif
*       ..Insert OMP list
          m%cfCv(m%cfCvP(icf,1):m%cfCvP(icf,1)+m%cfCvP(icf,2)-1) =
     &      omp(1:nomp)
*       ..For now zero-out faces belonging to OMP (to be improved)
          m%cfFc(m%cfFcP(icf,1):m%cfFcP(icf,1)+m%cfFcP(icf,2)-1) = 0
          m%cfFcP(icf,2) = 0
        endif
      end do

*   ..Definition of separatrix position in OMP and IMP list (to be improved for DN)
      if (m%IsClassicalGrid) call xertst(m%nncut.le.2,
     & 'Separatrix position to be revised for DN')
      icsepomp = 0
      ic = 1
      found = .false.
      do while (ic.le.nomp .and..not.found)
        if (m%cvReg(omp(ic)).ne.m%cvReg(omp(ic+1))) then
          found = .true.
          icsepomp = ic+1
        endif
      ic = ic + 1 
      end do
      if (.not.found) call xerrab ('Did not find separatrix location. '
     &  // 'Check specification OMP in b2.user.paramters.')
      icsepimp = 0
      ic = 1
      found = .false.
      do while (ic.le.nimp .and..not.found)
        if (m%cvReg(imp(ic)).ne.m%cvReg(imp(ic+1))) then
          found = .true.
          icsepimp = ic+1
        endif
      ic = ic + 1 
      end do
      if (.not.found) write(*,*) ('Did not find separatrix location. '
     &  // 'Check specification IMP in b2.user.paramters.')
      ifssep = 0 !flux surface of separatrix
      ic1 = omp(icsepomp-1) !cv of flux tube just inside sep
      ic2 = omp(icsepomp) !cv of flux tube just outside sep
      ifc1 = m%cvFcP(ic1,1)
      ifc2 = ifc1 + m%cvFcP(ic1,2) - 1
      ifc11 = m%cvFcP(ic2,1)
      ifc22 = ifc11 + m%cvFcP(ic2,2) - 1
      ! find face belonging to both ic1 and ic2
      do ii = ifc1, ifc2
        ifc = m%cvFc(ii)
        do jj = ifc11, ifc22
          ifcc = m%cvFc(jj)
          if (ifc.eq.ifcc) exit
        enddo
        if (ifc.eq.ifcc) exit
      end do
      ! find this face in fsFc --> ii
      do ii = 1, size(m%fsFc,1)
        if (ifc.eq.m%fsFc(ii)) exit
      end do
      ! find ii in fsFcP(ifs,1):fsFcP(ifs,1)+fsFcP(ifs,2)
      do ifssep = 1, size(m%fsFcP,1)
      ifc1 = m%fsFcP(ifssep,1)
      ifc2 = ifc1 + m%fsFcP(ifssep,2) - 1
        if (ii.ge.ifc1 .and. ii.le.ifc2) exit
      enddo
      if (m%IsClassicalGrid) call xertst(m%fsFcP(ifssep,2).gt.0,
     & 'No separatrix faces!')
      write(*,*) 'SEPARATRIX FLUX SURFACE AT ifs=',ifssep

      call subend ()
      return
 999  call xerrab ('Error while reading USER namelist!')
      end subroutine read_b2mod_user_namelist

      subroutine write_b2mod_user_namelist(nx,ns)
      implicit none
      integer nx, ns

      if (written) return
      if (user_namelist_used.and.write_nml_user) write(*,user)
      written = .true.
      return
      end subroutine write_b2mod_user_namelist

***********************************************************************

      subroutine cf_to_struct (m)
      type(mapping), intent(inout) :: m
      integer :: ic, icf, ic1, incf, ifc

      m%nCf = ncf
      m%mxnCf = count(m%fcLbl /= 0)
      
      call alloc_mapping_cf(m)

      incf = 0 ! index in cfCv
      do ic = 1,m%nCf
        m%cfCvP(ic,1) = incf + 1
        m%cfCvP(ic,2) = m%cfFcP(ic,2)
        do icf = 1, m%cfFcP(ic,2)
          incf = incf + 1
          ifc = m%cfFc(m%cfFcP(ic,1)+icf-1)
          ic1 = min(m%fcCv(ifc,1), m%fcCv(ifc,2))
          if (ic1 <= m%nCi) then
            m%cfCv(incf) = ic1
          else
            call xerrab('cf_to_struct -- problem with boundary face') 
          end if
        end do
      end do


      end subroutine cf_to_struct


***********************************************************************

      end module b2mod_user_namelist
