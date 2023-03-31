      module b2mod_user_namelist
      use b2mod_types
      use b2mod_diag
      use b2us_geo
      use b2us_map
      use b2mod_ad
     & , only : nncf, b2rr, b2volOnCf, b2data, b2dataOnCf
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
! - [oi]mp: list of CVs beloning to the [OI]MP, defined through the segment rz[oi]mp
!           it should be built as an ordered list from core to SOL. 
! - nromp: pre-defined maximum size for these lists
! - icsep[oi]mp: location of the separatrix position in the [oi]omp list
! - n[oi]mp: number of CVs in the [oi]mp list
! - rz[o1]mp: R-Z coordinates of a segment whose intersection with the
!             plasma grid defines the OMP and IMP. rzomp(1,1:2) holds R-coordiantes 
!             of first and second point of the segment, rzomp(2,1:2) holds Z-coordinates
      integer, parameter :: nromp=DEF_NYD*4
      integer, save :: omp(nromp), icsepomp, nomp, ifssep, imp(nromp),
     ,  icsepimp, nimp
      real (kind=R8) :: rzomp(2,2), rzimp(2,2)
      real (kind=R8) ::
     ,  fusion_power, Spmp_He_to_D, Spmp_nom
      integer lpfrb_i,lpfrb_o,lpfrt_i,lpfrt_o,
     ,  j_he_at,j_h_at(3),j_ne_at,lpfrs_pmp,
     ,  npfrgrp,ipfrgrp(npfrgrpd),jpfrgrp(npfrgrpd),
     ,  nntrgrp,intrgrp(nntrgrpd),jntrgrp(nntrgrpd)
! VARIABLES RELATED TO COST FUNCTION DEFINITION
! - mxnCf: maximum number of cost functions allowed
! - ncf: actual number of cost functions
! - cftype: type of cost function, 0-9 possible options
! - cfdef: defines how the cost function domain is defined 1-4 possible options
! - cfstart/cfend: starting and ending element for the cost function
! - cf_reg: Listing for each cost function the corresponding domain CVs or FCs
! - cf_regP: Points for each cost function to the first index in the CF_REG list, and its number of items.
! - cfReg and cfRegP: same as above, but allocatable and used to pass CF data into mapping structure
! - cfweight: multiplier to the cost function
! - cfread: according to cftype, decide if reading experimental data from file for cost function
! - ncfdata: number of experimental data points read from file for each cost function
! - cfdata(nncf,3,DEF_NLIM): stores the experimental data points for each cost function
!   column 1 stores the curvilinear spatial coordinate (e.g. along OMP or target)
!   column 2 stores the experimental data value
!   column 3 stores the uncertianty
! - nvmx: maximum number of variables that can be used for optimization
! - nsigmx: maximum number of sigma values used in the Bayesian MAP cost function
! - nsigma: actual number of sigmas read and used for Bayesian MAP
! - prior_type: type of prior for each design parameter optimized in the MAP estimation
! - prior_par: parameters related to the type of prior (e.g. standard deviation for gaussian prior)
! - prior_range: range of validity of the prior, outside this interval the MAP cost function returns NaN
! - scale_sigma: tells of the sigma in the MAP cost function is intended as relative (true) or absolute (false)
! - mapToOMP: according to cftype, decide if cost function data is to be re-mapped to the OMP 
! - read_sigma: tells to read or not the experimental data uncertainty
! - sigma: stores the values of the standard deviation in  the MAP cost function
! - inf_opt: parameter that defines an 'infinity' bound for optimization variables (e.g. prior range)
      integer, parameter :: mxnCf=1000
      integer cfstart(nncf),cfend(nncf),cftype(nncf)
      logical, save :: cfread(nncf) = .false.
      integer, save :: cfdef(nncf)
      integer, save :: ncfdata(nncf), ncf, cf_reg(mxnCf),
     & cf_regp(nncf,2)
      real (kind=R8), save :: cfdata(nncf,3,DEF_NLIM), cfweight(nncf)
      integer nvmx, nsigmx
      parameter (nvmx = 50, nsigmx=10)
      integer, save :: prior_type(nvmx), nsigma !sc type of prior
      real (kind=r8), save :: prior_par(nvmx,2), prior_range(nvmx,2) !parameters for prior
      logical, save :: scale_sigma(nsigmx) = .true.
      logical, save :: mapToOMP(nncf) = .false.
      logical, save :: read_sigma(nncf) = .false.
      real (kind=r8), save :: sigma(nsigmx) = 0.0_R8
      double precision, parameter :: inf_opt = 1.0e30
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
     ,  filedata,user_filename,ncf,cfstart,cfend,
     ,  cftype, rzomp, rzimp, cfdef, cfweight,
     ,  prior_type, prior_par, prior_range, scale_sigma, nsigma,
     ,  sigma, cf_reg, cf_regP, mapToOMP, read_sigma
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
c FIXME add description of Cf quantities
c----------------------------------------------------------------------

      contains

      subroutine alloc_b2mod_user(geo,m,ns,nlim,nmol,ntns)
      implicit none
      integer, intent(in) :: ns, nlim, nmol, ntns
      type (geometry), intent (in) :: geo
      type(mapping), intent(inout) :: m

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
      cfdef = 0
      cftype = 0
      cfdata = 0.0_R8
      ncfdata = 0
      read_sigma = .false.
      rzomp = 0.0_R8
      rzimp = 0.0_R8
      cfweight = 1.0_R8
      nsigma = 0
      scale_sigma = .true.
      prior_type = -1
      prior_par = -1.0_R8
! sc  Set prior range to (-inf,+inf)
      prior_range(:,1) = -inf_opt*10.0
      prior_range(:,2) = inf_opt*10.0
      cf_reg = 0
      cf_regP = 0
      mapToOMP = .false.
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
      if (allocated(b2rr)) deallocate(b2rr)
      if (allocated(b2volOnCf)) deallocate(b2volOnCf)
      if (allocated(b2data)) deallocate(b2data)
      if (allocated(b2dataOnCf)) deallocate(b2dataOnCf)
      user_initialised = .false.
      return
      end subroutine dealloc_b2mod_user

      subroutine read_b2mod_user_namelist(ns,ntns,geo,m)
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
      integer, intent(in) :: ns,ntns
      type (geometry), intent (in) :: geo
      type(mapping), intent(inout) :: m
      integer ic, noss, iss, indss,
     >        ncffc, incf, idb, icf,
     >        ic1, ic2, icff, icv, ifc, ifcc, ifc1, ifc2, ifs,
     >        ii, ifc11, ifc22, jj

      logical :: file_ok
      character*260 filename, cffile
      character*1 str
      integer i, ndata
      integer, allocatable :: cfReg(:)
      logical done, ifail, found
      integer, save :: use_eirene, neutrals_namelist
#ifdef WG_TODO
      integer j, k, is, im
      logical old_style
      logical, save :: hydrogen, helium, neon
#endif
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
      call cv_intersections (m%nCv, geo, m, rzomp, omp, nomp, nromp,
     & 1, ifail)
      call xertst(.not.ifail,'Increase nromp in b2.user.parameters')
! Here we need to guarantee ordering of CVs. With classical grids the conversion
! is done such that CVs are ordered from left to right and bottom to top wrt
! the structured computational domain. With fully unstructured grids this is not guaranteed.
! IDEA: project the cell center onto the segment and order CVs based on distance wrt one of the segment's extremes
! SAME FOR IMP
      call xertst(m%isClassicalGrid.eq.1,
     & ' OMP list need to be revised for fully unstructured grids')
      do ic=2,nomp
        if (geo%cVx(omp(ic)).lt.geo%cVx(omp(ic-1)))
     & call xerrab('OMP LIST NOT SORTED')
      end do
*   ..Definition of separatrix position in OMP list (to be improved for DN)
! Separatrix is defined by checking when cvReg changes value inside the OMP list
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

      call xertst (ncf.le.nncf, 'increase nncf in b2mod_user_namelist')
      if (ncf > 0) then
        write (*,*) ' Number of cost functions: ',ncf
        m%nCf = ncf
        m%mxnCf = mxnCf
        call alloc_mapping_cf (m)
        m%cfRegP = 0
        allocate(cfReg(mxnCf))
        cfReg = 0
#ifdef F2003
        open (newunit=idb,file='debug_user.out')
#else
        idb = newunit()
        open (unit=idb,file='debug_user.out')
#endif
        incf = 0
        do icf = 1, ncf
          if (cftype(icf).gt.10) call xerrab ('cftype>10 not coded ')
          if (cftype(1).eq.0 .and. ncf.eq.1) call xerrab(
     &     'cftype(1)=0 but no other cost functions are defined')
          if (cftype(icf).eq.0.and. icf.ne.1) call xerrab (
     &     'only cftype(1) should be =0')
          if (cftype(icf).eq.6.and. icf.ne.1) call xerrab (
     &     'only cftype(1) should be =6')
          if ((cftype(icf).eq.0 .or. cftype(icf).eq.6) .and.
     &     ncf.lt.(cfend(icf)-cfstart(icf))+2) call xerrab (
     &     'cost functions specified for sum with cftype 0/6 > nCf')
          if ((cftype(icf).eq.0 .or. cftype(icf).eq.6) .and.
     &     cfstart(icf).lt.2 .and. cfend(icf).lt.2) call xerrab (
     &     'cost functions type 0/6 requires cfstart/cfend>1')
          if ((cftype(icf).ge.1 .and. cftype(icf).le.3) .or.
     &        (cftype(icf).ge.7 .and. cftype(icf).le.10))
     >         cfread(icf) = .true.
          if (cftype(icf).eq.5) m%cfOnCv(icf) = .false.
          if (cftype(icf).eq.5 .and. cfdef(icf).eq.3)
     &      call xerrab('cftype 5 must be defined using faces!')
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
!csc The data in the file are (r-rsep), data, sigma. Stored in cfdata according to
!    first dimension (# of cost functions), second dimension ((r-rsep)-data-sigma), third
!    dimension (number of records/data points)
              if (read_sigma(icf)) then
                read(99,*) cfdata(icf,1:3,i)
              else
                read(99,*) cfdata(icf,1:2,i)
              endif
            enddo
            close(99)
          endif
          if (cftype(icf).eq.0 .or. cftype(icf).eq.6) cycle !no need to specify domain here
          done  = .false.

! FIXME need to add flag for reversing lists?

          select case(cfdef(icf))

          case (1) !defined on OMP (only internal CVs)

            call xertst(m%cfOnCv(icf), 'This cost function type '//
     &       'cannot be evaluated on OMP')
            call xertst(nomp.gt.0, 'Asking for OMP cost function '//
     &       'but no OMP cells, check rzomp')
            write (*,*) ' defining cost function ',icf,' on omp'
            do icff = 1, icf-1
              if (cfdef(icff) .eq. 1) then
              m%cfRegP(icf,:) = m%cfRegP(icff,:)
              done = .true.
              endif
            end do
            if (done) cycle
            m%cfRegP(icf,1) = incf +1
            ndata = 0
            do iCv=1,nomp
              if (omp(iCv)>m%nCi) cycle !only use internal CVs
              m%cfReg(incf +1+ndata) = omp(iCv)
              ndata = ndata + 1
            enddo
            m%cfRegP(icf,2) = ndata
            incf = incf + ndata

          case (2) !defined using fclabel

            call xertst(cfstart(icf) /= 0 .and. cfend(icf) /= 0,
     &        'cfdef=2 requires cfstart/cfend>0')
! check if other CFs have same domain (must have same start/end, same type of definition (fcLabl)
! and also same type of domain (FcLabl can be used to identify either CVs or FCs)
            do icff = 1, icf-1
              if (cfstart(icf) .eq. cfstart(icff) .and.
     &            cfend(icf) .eq. cfend(icff) .and.
     &            cfdef(icff) .eq. 2 .and.
     &            (m%cfOnCv(icf) .and. m%cfOnCv(icff))) then
                m%cfRegP(icf,:) = m%cfRegP(icff,:)
                done = .true.
              endif
            end do
            if (done) cycle
! this domain has not yet been done
            m%cfRegP(icf,1) = incf +1
            ncffc = 0
            noss = cfend(icf) - cfstart(icf) + 1
            do iss = 1, noss
              indss = cfstart(icf)+iss-1
              write (*,*) ' Defining cost function ',icf
              write (*,*) ' surface label ',indss
! find all faces belonging to surface structure INDSS
              call find_faces (indss, incf, ncffc,
     >           m%mxncf,m%cfReg, m%cfFcOr, m, idb)
            end do              ! iss	
            m%cfRegP(icf,2) = ncffc
! if cost function is on CV e.g. target for CV-based CF --> substitute FCs with linked CVs
            if (m%cfOnCv(icf)) then
              do iFc = 1, m%cfRegP(icf,2)
                ifcc = m%cfReg(m%cfRegP(icf,1)+ifc-1)
                ic1 = min(m%fcCv(ifcc,1), m%fcCv(ifcc,2))
                ic2 = max(m%fcCv(ifcc,1), m%fcCv(ifcc,2))
                if (ic1 <= m%nCi) then
                  m%cfReg(m%cfRegP(icf,1)+ifc-1) = ic1
                elseif (ic2 <= m%nCi) then
                  m%cfReg(m%cfRegP(icf,1)+ifc-1) = ic2
                else
                  call xerrab('cfdef=2 -- problem with boundary face')
                end if
              end do
            endif

          case (3) !defined using cvlabel

            call xertst(cfstart(icf) /= 0 .and. cfend(icf) /= 0,
     &        'cfdef=3 requires cfstart/cfend>0')
! check if other CFs have same domain !FIXME check if also cfdef is the same
            do icff = 1, icf-1
              if (cfstart(icf) .eq. cfstart(icff) .and.
     &            cfend(icf) .eq. cfend(icff) .and.
     &            cfdef(icff) .eq. 3) then
                m%cfRegP(icf,:) = m%cfRegP(icff,:)
                done = .true.
              endif
            end do
            if (done) cycle
! this domain has not yet been done
            m%cfRegP(icf,1) = incf +1
            noss = cfend(icf) - cfstart(icf) + 1
            do iss = 1, noss
              indss = cfstart(icf)+iss-1
              write (*,*) ' Defining cost function ',icf
              write (*,*) ' cell label ',indss
              do icv = 1, m%nCv
                if (m%cvLbl(icv).eq.indss) then
                  m%cfReg(incf+1) = icv
                  incf = incf + 1
                endif
              end do
            end do
            ! FIXME no attempt to sort or check if they belong together somehow
            m%cfRegP(icf,2) = incf+1 - m%cfRegP(icf,1)

          case (4) !defined using cf_reg (No attempt whatsoever to check if they are sorted)

            call xertst(cf_RegP(icf,2).gt.0,'cf_RegP not defined')
            ! csc can we check somehow if this has already been filled?
            m%cfRegP(icf,1) = incf +1
            m%cfRegP(icf,2) = cf_RegP(icf,2)
            if (any(cf_Reg(cf_RegP(icf,1) : cf_RegP(icf,1) +
     &       cf_RegP(icf,2) - 1) .le. 0)) call xerrab('cf_Reg<=0')
            m%cfReg(incf +1:incf +1+m%cfRegP(icf,2)-1) =
     &        cf_Reg(cf_RegP(icf,1):cf_RegP(icf,1)+cf_RegP(icf,2)-1)
            incf = incf + cf_RegP(icf,2)

          case default

            call xerrab('cfDef undefined')

          end select

          call xertst(incf+1.le.m%mxncf,
     &      'Increase size of mxnCf in b2mod_user_namelist')
        end do !icf
        close(idb)
        if (allocated(cfReg)) deallocate(cfReg)
      end if

      if (any(cftype.eq.6)) then
        call xertst(nsigma.gt.0.and.nsigma.le.nsigmx,
     &   'b2mod_user_namelist: nsigma<=0 or nsigma>nsigmx')
! sc  prior_type refers to the prior type for the design variables:
!     0 - Uniform distribution.
!     1 - Uninformative, proper, Gaussian prior.
!     2 - Jeffrey's 1/sigma (only for sigma!!)
        call xertst(any(prior_type(1:nsigma).ge.0).and.
     &   any(prior_type(1:nsigma).le.2),
     &   'b2mod_user_namelist: prior_type<=0 or prior_type>=2')
! sc  prior_par refers to the prior parameters:
!     For prior_type = 1, prior_par(ii,1) = mu, prior_par(ii,2) = sigma
!     For prior_type = 2, not used?
        do icf = 1, nsigma
          if (prior_type(icf) .eq. 1)
     &      call xertst(prior_par(icf,2) .gt. 0.0_R8,
     &     'b2mod_user_namelist: sigma prior <=0')
          call xertst(prior_range(icf,1).lt.prior_range(icf,2),
     &     'b2mod_user_namelist: prior_range(.,1)>prior_range(.,2)')
        end do
      endif
! now allocate vector that will be used to interpolate solps onto data grid
      if (any(cfread)) then
        ndata = maxval(ncfdata(1:ncf),1)
        allocate(b2volOnCf(ncf,ndata)) !stores interpolated B2.5 volumes onto CF radial points
        allocate(b2dataOnCf(ndata)) !variable to store SOLPS data interpolated onto CF radial points
        ndata = maxval(m%cfRegP(1:ncf,2)) !max number of CVs in a cost function
        allocate(b2rr(ncf,ndata)) !store here radial distance coordiante of SOLPS data
        allocate(b2data(ndata)) !temporary variable to store SOLPS data for interpolation
        b2rr = 0.0_R8
        b2volOnCf = 0.0_R8
        b2data = 0.0
        b2dataOnCf = 0.0
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
#endif
*   ..Create list of IMP control volumes
! Same logic as OMP list, find CVs intersected by rzimp segment
      nimp = 0
      imp = 0
      call cv_intersections (m%nCv, geo, m, rzimp, imp, nimp, nromp,
     & 1, ifail)
      call xertst(.not.ifail,'Increase nromp in b2.user.parameters')
      do ic=2,nimp
        if (geo%cVx(imp(ic)).lt.geo%cVx(imp(ic-1)))
     & call xerrab('IMP LIST NOT SORTED')
      end do
! Find position of separatrix in IMP list (same as for OMP)
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
      if (.not.found) write(*,'(a/a/a)')
     &  'Warning: did not find separatrix location in IMP definition, ',
     &  '         check specification of RZIMP in b2.user.parameters',
     &  '         if IMP is needed.'
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
      end do
      ! find this face in fsFc --> ii
      do ii = 1, size(m%fsFc,1)
        if (ifc.eq.m%fsFc(ii)) exit
      end do
      ! find ii in fsFcP(ifs,1):fsFcP(ifs,1)+fsFcP(ifs,2)
      do ifs = 1, size(m%fsFcP,1)
        ifc1 = m%fsFcP(ifs,1)
        ifc2 = ifc1 + m%fsFcP(ifs,2) - 1
        if (ii.ge.ifc1 .and. ii.le.ifc2) exit
      enddo
      ifssep = ifs
      if (m%IsClassicalGrid.eq.1) call xertst(m%fsFcP(ifssep,2).gt.0,
     & 'No separatrix faces!')
      write(*,*) 'SEPARATRIX FLUX SURFACE AT ifs=',ifssep

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
