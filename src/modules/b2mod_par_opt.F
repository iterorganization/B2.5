      module b2mod_par_opt
      use b2mod_types
      use b2mod_user_namelist
     & , only : omp, nomp
      use b2mod_ad
     & , only : nncf, b2rr, b2volOnCf, b2data, b2dataOnCf
      use b2us_map
      implicit none
#     include "DIMENSIONS.F"
! VARIABLES RELATED TO COST FUNCTION DEFINITION
! - mxnCf: maximum number of cost functions allowed
! - ncf: actual number of cost functions
! - cftype: type of cost function, 0-9 possible options
! - cfdef: defines how the cost function domain is defined, 1-4 possible options
! - cfstart/cfend: starting and ending element for the cost function
! - cf_reg: listing for each cost function the corresponding domain CVs or FCs
! - cf_regP: points for each cost function to the first index in the CF_REG list, and its number of items.
! - cfReg and cfRegP: same as above, but allocatable and used to pass CF data into mapping structure
! - cfweight: multiplier to the cost function
! - cfread: according to cftype, decide if reading experimental data from file for cost function
! - ncfdata: number of experimental data points read from file for each cost function
! - cfdata(nncf,4,DEF_NLIM): stores the experimental data points for each cost function
!   column 1 stores the curvilinear spatial coordinate (e.g. along OMP or target)
!   column 2 stores the experimental data value
!   column 3 stores the standard deviation sigma
!   column 4 stores the mean
! - nvmx: maximum number of variables that can be used for optimization
! - nsigmx: maximum number of sigma and mean values used in the Bayesian MAP cost function
! - nsigma: actual number of sigmas read and used for Bayesian MAP
! - nmean: actual number of error means read and used for Bayesian MAP
! - prior_type: type of prior for each design parameter optimized in the MAP estimation
! - prior_par: parameters related to the type of prior (e.g. standard deviation for Gaussian prior)
! - prior_range: range of validity of the prior, outside this interval the MAP cost function returns NaN
! - scale_sigma: tells of the sigma in the MAP cost function is intended as relative (true) or absolute (false)
! - mapToOMP: according to cftype, decide if cost function data is to be re-mapped to the OMP 
! - read_sigma: tells to read or not the experimental data uncertainty
! - sigma: stores the values of the prediction error standard deviation in the MAP cost function
! - mean: stores the values of the prediction error mean in the MAP cost function
! - inf_opt: parameter that defines an 'infinity' bound for optimization variables (e.g. prior range)
      integer, parameter :: mxnCf=1000
      integer cfstart(nncf),cfend(nncf),cftype(nncf)
      logical, save :: cfread(nncf) = .false.
      integer, save :: cfdef(nncf)
      integer, save :: ncfdata(nncf), ncf, cf_reg(mxnCf),
     & cf_regp(nncf,2)
      real (kind=R8), save :: cfdata(nncf,4,DEF_NLIM), cfweight(nncf)
      integer nvmx, nsigmx
      parameter (nvmx = 50, nsigmx=10)
      integer, save :: prior_type(nvmx), nsigma !sc type of prior
      real (kind=r8), save :: prior_par(nvmx,2), prior_range(nvmx,2) !parameters for prior
      logical, save :: scale_sigma(nsigmx) = .true.
      logical, save :: mapToOMP(nncf) = .false.
      logical, save :: read_sigma(nncf) = .false.
      real (kind=r8), save :: sigma(nsigmx) = 0.0_R8
      real (kind=r8), parameter :: inf_opt = 1.0e30
      integer, save :: nmean = 0
      real (kind=r8), save :: mean(nsigmx) = 0.0_R8
! VARIABLES RELATED TO OPTIMIZATION
! - nnvar: number of optimization parameters (not including radially varying transport coefficients)
! - npar_opt: actual number of optimization parameters (including spatially varying transport coeff.)
! - par_opt_phys: array containing the values of the optimization parameters used to pass them to the
!                 optimization libraries through b2optim_tao/b2optim_ipopt, and to the routine for
!                 calculating the prior (if needed)
! - nsigma_opt stores only the number of sigmas that are actually optimized (which can be .le. nsigma).
!   similarly sigma_opt stores only such sigmas that are optimized. This differentiation is needed to
!   correctly pass the optimization variables to the optimization libraries
! - nmean_opt, mean_opt, do the same as for sigma but for the error mean
! - x0: initial guess of the optimization parameters (one for each actual parameter, so includes the
!        number of spatial points!)
! - xl/xu: lower/ipper bound of optimization parameters (one for each actual parameter, as x0)
! - nncon: number of nonlinear constraints (not used for now)
! - gl/gu: lowr/upper bound of such nonlinear constraints (not used for now)
! - nnjac: number of non-zeros in the jacobian of the nonlinear constraints (not used for now)
! - jcol/jrow: column/row indicex for each non-zero element of the jacobian
! - jj: value of the jacobian at jcol-jrow
      integer, save :: npar_opt = 0
      real (kind=R8), allocatable, save :: par_opt0(:), par_opt(:),
     &       par_opt_phys(:)
      logical, save :: flag_optim = .false.
      integer nnvar, nncon, nnjac
      real (kind=R8) :: x0(nvmx), xl(nvmx), xu(nvmx), gl(nvmx),
     &  gu(nvmx), jj(nvmx*nvmx),
     &  par_rescale(nvmx)
      integer :: jcol(nvmx*nvmx), jrow(nvmx*nvmx)
! sc some variables that allows to change optimization parameters, read from b2.optimization.parameters
      character*256, save :: limited_memory_update_type,
     &  hessian_approximation 
      real (kind=R8), save :: cpu_opt=0.0_R8, tol_opt=1.0e-7_R8
      integer, save :: maxiter=100, partype(nvmx), nsigma_opt = 0,
     &                 nmean_opt = 0
      integer, save :: paris(nvmx),parib(nvmx)
      logical, save :: sigma_opt(nsigmx), mean_opt(nsigmx)
      logical, save :: spatial_dep(nvmx)
      integer, save :: spatial_points(nvmx) = 0
      integer :: idd
      logical :: file_ok
      real (kind=R8) :: rangeold(nvmx,2), parold(nvmx,2)
      integer :: typeold(nvmx)
      private optimization
      namelist /optimization/ nnvar, nncon, nnjac, x0, xl, xu, gl, gu, 
     & jcol, jrow, jj, maxiter, cpu_opt, tol_opt, hessian_approximation,
     & limited_memory_update_type, partype, sigma_opt, spatial_dep,
     & spatial_points, par_rescale, paris, parib, ncf,cfstart,cfend,
     & cftype, cfdef, cfweight,prior_type, prior_par, prior_range, 
     & scale_sigma, nsigma, sigma, cf_reg, cf_regP, mapToOMP, 
     & read_sigma, nmean, mean, mean_opt

      contains

      subroutine read_b2mod_par_opt(ncon, nele_jac, ns, m)
      use b2mod_types
      implicit none
      integer, intent (out) :: ncon, nele_jac
      type(mapping), intent(inout) :: m
      integer, intent (in) :: ns
      ! csc local variables
      integer :: ii, isigma, ipp, i, ndata, iss, indss, icf, icff,
     &           noss, ncffc, incf, idb, ic1, ic2, icv, ifc, ifcc,
     &           ifc1, ifc2, jj, imean
      integer, allocatable :: cfReg(:)
      logical done
      character*1 str
      character*256 cffile, filename
#ifndef F2003
      integer newunit
      external newunit
#endif

      filename='b2.optimization.parameters'
      write (*,*) 'OPTIM: max number of readable cost functions :',nncf
      write (*,*) 'OPTIM: max number of readable parameters '//
     & 'and constraints:',nvmx
      write (*,*) 'OPTIM: infinity bound set to:',inf_opt
      ncf = 0
      cfstart = 0
      cfend = 0
      cfdef = 0
      cftype = -1
      cfdata = 0.0_R8
      ncfdata = 0
      read_sigma = .false.
      cfweight = 1.0_R8
      nsigma = 0
      nmean = 0
      scale_sigma = .true.
      prior_type = -1
      prior_par = -1.0_R8
! sc  Set prior range to (-inf,+inf)
      prior_range(:,1) = -inf_opt*10.0_R8
      prior_range(:,2) = inf_opt*10.0_R8
      cf_reg = 0
      cf_regP = 0
      mapToOMP = .false.
      x0 = inf_opt*10.0_R8
      xl = -inf_opt*10.0_R8
      xu = inf_opt*10.0_R8
      gl = -inf_opt*10.0_R8
      gu = inf_opt*10.0_R8
      par_rescale = 1.0_R8
      NELE_JAC = 0
      nnjac = 0
      jcol = 0
      jrow = 0
      jj = 0.0
      partype = -2
      sigma_opt = .true.
      mean_opt = .false.
      spatial_dep = .false.
      spatial_points = 0
      hessian_approximation = 'limited-memory'
      limited_memory_update_type = 'b2fgs'
      call find_file(filename,file_ok)
      if(file_ok) then
        open(99,file=filename)
        write(*,*) 'Reading namelist optimization from ',trim(filename)
#ifdef USE_MPI
        read(99,optimization)
#else
        read(99,optimization)
#endif
        close(99)
        write(*,*) 'Read namelist optimization from ',trim(filename)
      else
        call xerrab ('No '//trim(filename))
      endif
      write (*,*) 'OPTIM: read parameters'
      write(*,optimization)
! checks on cost function parameters
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
        open (newunit=idb,file='debug_optimization.out')
#else
        idb = newunit()
        open (unit=idb,file='debug_optimization.out')
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
! check if other CFs have same domain (must have same start/end, same type of definition (fcLbl)
! and also same type of domain (fcLbl can be used to identify either CVs or FCs)
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
! check if other CFs have same domain
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

          case (4) ! defined using cf_reg (No attempt whatsoever to check if they are sorted)

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

! now allocate vector that will be used to interpolate SOLPS onto data grid
      if (any(cfread)) then
        ndata = maxval(ncfdata(1:ncf),1)
        allocate(b2volOnCf(ncf,ndata)) !stores interpolated B2.5 volumes onto CF radial points
        allocate(b2dataOnCf(ndata)) !variable to store SOLPS data interpolated onto CF radial points
        ndata = maxval(m%cfRegP(1:ncf,2)) !max number of CVs in a cost function
        allocate(b2rr(ncf,ndata)) !store here radial distance coordinate of SOLPS data
        allocate(b2data(ndata)) !temporary variable to store SOLPS data for interpolation
        b2rr = 0.0_R8
        b2volOnCf = 0.0_R8
        b2data = 0.0
        b2dataOnCf = 0.0
      endif

! csc Tests on optimization parameters
      call xertst(nnvar.ge.0.and.nnvar.le.nvmx,
     & 'b2mod_par_opt: nnvar<0 or nnvar>nvmx')
      call xertst(nncon.le.nvmx,'b2mod_par_opt: increase size of nvmx')
      call xertst(nnjac.le.nvmx,'b2mod_par_opt: increase size of nvmx')
      call xertst(nncon.ge.0,'b2mod_par_opt: nncon<0')
      call xertst(nnjac.ge.0,'b2mod_par_opt: nnjac<0')
      call xertst(all(par_rescale.gt.0.D0),
     & 'b2mod_par_opt: par_rescale<=0')
      if (flag_optim) call xertst(nnvar.gt.0,
     & 'b2mod_par_opt: cannot optimize if nnvar=0!')
      if (nnvar .gt. 0) then
!     partype refers to the type of design variable:
!    -1  - mean, bias/mean of the prediction error for Bayesian inference. They must be the last parameters in the vector!
!     0  - sigma, standard deviation of the prediction error for Bayesian inference. They must be specified after the 'physical' parameters and before the means
!     1  - Dna, density driven particle diffusion coefficient
!     2  - Dpa, pressure driven particle diffusion coefficient
!     3  - X_i, ion heat diffusion coefficient
!     4  - X_e, electron heat diffusion coefficient
!     5  - vlax, x component of the anomalous velocity
!     6  - vlay, y component of the anomalous velocity FIXME not yet done
!     7  - vsa, anomalous viscosity
!     8  - sig, field-driven current radial conductivity
!     9  - alf, temperature-driven current radial conductivity
!     10 - enepar(parib(ipar), 1), with parib(ipar) identifies boundary
!     11 - enipar(parib(ipar), 1), with parib(ipar) identifies boundary
!     12 - conpar(paris(ipar), parib(ipar), 1), with paris(ipar) identifies species, parib(ipar) identifies boundary
!     13 - mompar(paris(ipar), parib(ipar), 1), with paris(ipar) identifies species, parib(ipar) identifies boundary
!     14 - potpar(parib(ipar), 1), with parib(ipar) identifies boundary
!     15 - enkpar(parib(ipar), 1), k boundary condition, parib(ipar) identifies boundary
!     16 - b2recyc(paris(ipar),parib(ipar)), with paris(ipar) identifies species, parib(ipar) identifies stratum
!     17 - b2tqna_ballooning
!     18 - b2tqna_ballooning_rescale
!     19 - keps_cd, k-model coefficient in Dna_ExB ~ Cd*sqrt(k) 
!     20 - keps_heat, multiplier Ce for X_e = Ce*Dna_ExB 
!     21 - keps_heat_i, multiplier Ci for X_i = Ci*Dna_ExB 
!     22 - keps_sig
!     23 - keps_alf
!     24 - keps_visc
!     25 - keps_dkt
!     26 - keps_dzt
!     27 - b2sikt_fac_diss
!     28 - b2sikt_fac_diss_core
!     29 - b2sikt_fac_sheath, multiplier Cl for Sk_diss = Cl*k*cs/Lconn (valid in SOL & PFR)
!     30 - b2sikt_fac_sheath_core, multiplier Cl for Sk_diss = Cl*k*cs/Lconn (valid in core)
!     31 - keps_shear
!     32 - b2tfhi_fsigkt
!     33 - b2sikt_fac_vis_RS
!     34 - b2tfhi_fflokt
!     35 - b2tfhi_fconkt
!     36 - b2tfhi_fflozt
!     37 - b2tfhi_fconzt
!     38 - b2tfhi_fkt_hie
!     39 - b2tfhe_vis_kt
       call xertst(any(partype(1:nnvar) .ge. -1) .and.
     &   any(partype(1:nnvar) .le. 39),
     &   'b2mod_par_opt: partype<-1 or partype>39')
       call xertst(any(spatial_points(1:nnvar) .ge. 0),
     &   'b2mod_par_opt: spatial_points<0')
       call xertst(sum(spatial_points(1:nnvar)) .le. nvmx,
     &   'b2mod_par_opt: sum(spatial_points)<=nvmx, increase nvmx')
       do ii=1,nnvar
        if (spatial_dep(ii) .and. 
     &      (partype(ii).lt.1 .or. partype(ii).gt.9)) !only parm_XXX can be spatially dependent
     >  call xerrab('b2mod_par_opt: spatial_points optimization '//
     &              'is available only for partype 1-9')
        if (spatial_dep(ii) .and. spatial_points(ii).le.1) then
          write(*,*) 'ipar, spatial_dep, spatial_points=',
     &      ii, spatial_dep(ii), spatial_points(ii)
          call xerrab(' b2mod_par_opt: seems you want to optimize '//
     &     'spatially dependent transport coefficients but you only '//
     &     'specify one!')
        endif
        if (partype(ii) .eq. 1 .or. partype(ii) .eq. 2 .or.
     &      partype(ii) .eq. 3 .or. partype(ii) .eq. 5 .or.
     &      partype(ii) .eq. 7 .or. partype(ii) .eq. 12 .or.
     &      partype(ii) .eq. 13 .or. partype(ii) .eq. 16) 
     >   call xertst((paris(ii) .ge. 0) .and. (paris(ii) .le. ns-1),
     &   'b2mod_par_opt: paris<0 or >ns-1')
        if (partype(ii) .ge. 10 .and. partype(ii) .le. 16)
     >   call xertst((parib(ii) .ge. 1) .and. (parib(ii) .le. m%nBc),
     &    'b2mod_par_opt: parib<0 or >nBc')
       enddo
      endif !nnvar>0

!     MAP estimation is turned on only when cftype(1)=6
      if (cftype(1) .eq. 6) then
!      First do some tests on means as they come last in vector partype
       call xertst(nmean.le.nsigmx,
     &   'b2mod_par_opt: nmean>nsigmx, please increase nsigmx')
       if (nmean .gt. 0 .and. any(mean_opt(1:nmean))) then
!       calculate first how many means are being optimized
        nmean_opt = 0
        do ii = 1, nmean
         if (mean_opt(ii)) nmean_opt = nmean_opt + 1
        end do
!       check how many means are specified as independent variables in partype
        ipp = 0
        do ii = 1, nnvar
          if (partype(ii).eq.-1) ipp = ipp + 1
        end do
        if (.not. (ipp.eq.nmean_opt)) then
          write(*,*) 'number of mean specified in partype = ',ipp
          write(*,*) 'number of mean optimized according'//
     &               ' to mean_opt = ',nmean_opt
          call xerrab("b2mod_par_opt: inconsistent number of mean "//
     &                "between partype and mean_opt")
        endif
!       check if means are at the end of the parameter vector, if not, issue error
        imean = nnvar - nmean_opt + 1 ! position of first mean in partype
        do ii = 1, nmean_opt
         call xertst(partype(ii+imean-1) .eq. -1,'b2mod_par_opt: '//
     &    'partype=-1 variables MUST be at the end of the vector')
        end do
       endif

!      Now do some tests on sigmas
       call xertst(nsigma.le.nsigmx,
     &   'b2mod_par_opt: nsigma>nsigmx, please increase nsigmx')
!      We expect one sigma for each CF used in the MAP
       ipp = abs(cfend(1)-cfstart(1))+1
       call xertst(nsigma.eq.ipp,
     &   'b2mod_par_opt: must specify a sigma for each CF used in MAP!')
       if (nsigma .gt. 0 .and. any(sigma_opt(1:nsigma))) then
!       calculate first how many sigmas are being optimized
        nsigma_opt = 0
        do ii = 1, nsigma
         if (sigma_opt(ii)) nsigma_opt = nsigma_opt + 1
        end do
!       check how many sigmas are specified as independent variables in partype
        ipp = 0
        do ii = 1, nnvar
          if (partype(ii).eq.0) ipp = ipp + 1
        end do
        if (.not. (ipp.eq.nsigma_opt)) then
          write(*,*) 'number of sigma specified in partype = ',ipp
          write(*,*) 'number of sigma optimized according'//
     &               ' to sigma_opt = ',nsigma_opt
          call xerrab("b2mod_par_opt: inconsistent number of sigma "//
     &                "between partype and sigma_opt")
        endif
!       check if sigmas are after the physical parameter vector and before the means, if not, issue error
        isigma = nnvar - nsigma_opt - nmean_opt + 1 ! position of first sigma in partype
        do ii = 1, nsigma_opt
         call xertst(partype(ii+isigma-1) .eq. 0,'b2mod_par_opt: '//
     &    'partype 0 variables MUST be after the physical parameters'//
     &    ' (PARTYPE>0) and before the error means (if any)')
        end do
       endif

!      Now do some tests on priors
!      prior_type refers to the prior type for the design variables:
!      0 - Uniform distribution.
!      1 - Uninformative, proper, Gaussian prior.
!      2 - Gamma distribution. (defined through mean and std)
!      3 - Jeffrey's 1/sigma (only for sigma!!)
       call xertst(any(prior_type(1:nnvar) .ge. 0).and.
     &  any(prior_type(1:nnvar) .le. 3),
     &  'b2mod_par_opt: prior_type<0 or prior_type>3')
!      prior_par refers to the prior parameters:
!      For prior_type = 1,2, prior_par(ii,1) = mu, prior_par(ii,2) = sigma
!      For prior_type = 3, not used?
       do ii = 1, nnvar
        if (prior_type(ii) .eq. 1) 
     &   call xertst(prior_par(ii,2) .gt. 0.0_R8,
     &    'b2mod_par_opt: sigma prior <=0')
        if (prior_type(ii) .eq. 2) then
         call xertst((prior_par(ii,1) .gt. 0.0_R8) .and. 
     &               (prior_par(ii,2) .gt. 0.0_R8),
     &    'b2mod_par_opt: mean and variance for gamma prior must be >0')
         if (max(prior_range(ii,1),xl(ii)).le.0.0_R8) then
           write(*,*) 'ii, prior_range(ii,1)',ii, prior_range(ii,1)
           write(*,*) 'ii, xl(ii)',ii, xl(ii)
           call xerrab('b2mod_par_opt: gamma distribution is '//
     &     'intended for parameters on domain (0 +inf)')
         endif
        endif
        call xertst(prior_range(ii,1).lt.prior_range(ii,2),
     &   'b2mod_par_opt: prior_range(.,1)>prior_range(.,2)')
       end do
      end if  !cftype(1).eq.6

      if (flag_optim) then
!     Now do some tests on constraints
!     For now no constraints are implemented apart from the trivial one in Ipopt
       call xertst(ncon.ge.0 .and. ncon.le.1,
     &  'b2mod_par_opt: ncon<0 or ncon>1 (not implemented yet)')
       if (nncon.gt.0) then
         write(*,*) 'WARNING! Specified ncon>0: one trivial '//
     &   'constraint is allowed for Ipopt. If different constraints '//
     &   'are required they must be implemented first in '//
     &   'b2optim_[ipopt-tao]'
       endif
       if (nnjac .gt. 0) then
        call xertst(any(jcol(1:nnjac) .gt. 0),'b2mod_par_opt: jcol<=0')
        call xertst(any(jrow(1:nnjac) .gt. 0),'b2mod_par_opt: jrow<=0')
        call xertst(maxval(jcol(1:nnjac)) .le. nncon,
     &  'b2mod_par_opt: jcol>ncon')
        call xertst(maxval(jrow(1:nnjac)) .le. nnvar,
     &  'b2mod_par_opt: jrow>nvar')
       endif
      endif

!     adjust some arrays when there are radially dependent coefficients
!     x0, xl, xu are supposed to have the correct dimension already, i.e. one for each spatial point
!     priors parameters are not
      npar_opt = 0
      idd = 1
      typeold = prior_type
      parold = prior_par
      rangeold = prior_range
      do ipp = 1, nnvar + (nsigma - nsigma_opt) + (nmean - nmean_opt)
       prior_type(idd) = typeold(ipp)
       prior_par(idd,1) = parold(ipp,1)
       prior_par(idd,2) = parold(ipp,2)
       prior_range(idd,1) = rangeold(ipp,1)
       prior_range(idd,2) = rangeold(ipp,2)
       if (spatial_dep(ipp)) then
        prior_type(idd:idd+spatial_points(ipp)-1) = typeold(ipp)
        prior_par(idd:idd+spatial_points(ipp)-1,1) = parold(ipp,1)
        prior_par(idd:idd+spatial_points(ipp)-1,2) = parold(ipp,2)
        prior_range(idd:idd+spatial_points(ipp)-1,1) = rangeold(ipp,1)
        prior_range(idd:idd+spatial_points(ipp)-1,2) = rangeold(ipp,2)
        npar_opt = npar_opt + spatial_points(ipp) - 1 
        idd = idd + spatial_points(ipp) - 1 !next index
       end if
       npar_opt = npar_opt + 1
       idd = idd + 1 !next index
      enddo
      npar_opt = npar_opt - (nsigma - nsigma_opt) - (nmean - nmean_opt)
      if (flag_optim .or. (cftype(1) .eq. 6)) then
        call xertst(any(x0(1:npar_opt) .lt. inf_opt*10.0_R8),
     &   'b2mod_par_opt: initial guess x0 MUST be specified for '//
     &   'all variables if optimizing or using a MAP cost function')
        do ipp = 1, npar_opt
         if ((x0(ipp).lt.xl(ipp)) .or. (x0(ipp).gt.xu(ipp))) then
          write (*,*) 'ipp, xl(ipp), x0(ipp), xu(ipp)', ipp, xl(ipp),
     &      x0(ipp), xu(ipp)
          call xerrab('b2mod_par_opt: x0 must be >= xl and <=xu')
         endif
        enddo
      endif
      isigma = npar_opt - nsigma_opt  - nmean_opt + 1 ! position of first sigma in x0
      do ii = 1, nsigma
!      only assign initial value from x0 if that sigma is being optimized or used in MAP
       if (sigma_opt(ii)) then
        call xertst(x0(isigma).gt.0.0_R8 .and. xl(isigma).gt.0.0_R8
     &              .and. xu(isigma).gt.0.0_R8, 'b2mod_par_opt: x0,'//
     &              'xl, and xu for sigmas MUST be >0')
        sigma(ii) = x0(isigma)
        isigma = isigma + 1 !next sigma in inital vector x0
       endif
      end do
      imean = npar_opt - nmean_opt + 1 ! position of first mean in x0
      do ii = 1, nmean
!      only assign initial value from x0 if that mean is being optimized or used in MAP
       if (mean_opt(ii)) then
        mean(ii) = x0(imean)
        imean = imean + 1 !next mean in inital vector x0
       endif
      end do
      ncon = nncon
      nele_jac = nnjac
      write (*,*) 'NPAR_OPT = ',npar_opt
      write (*,*) 'NCON = ',ncon
      write (*,*) 'NELE_JAC = ',NELE_JAC

      return
      end subroutine read_b2mod_par_opt

      subroutine dealloc_b2mod_par_opt
      implicit none
      if (allocated(b2rr)) deallocate(b2rr)
      if (allocated(b2volOnCf)) deallocate(b2volOnCf)
      if (allocated(b2data)) deallocate(b2data)
      if (allocated(b2dataOnCf)) deallocate(b2dataOnCf)
      end subroutine dealloc_b2mod_par_opt

      end module b2mod_par_opt
