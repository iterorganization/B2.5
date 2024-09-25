      module b2mod_par_opt
      use b2mod_types
      use b2mod_user_namelist
     & , only : omp, nomp
      use b2mod_ad
     & , only : nncf, b2rr, b2volOnCf, b2data, b2dataOnCf, b2psi
      use b2us_map
      use b2mod_dimensions
      use b2mod_transport_namelist
     & , only : flag_dna, flag_dpa, flag_hci, flag_hce, flag_vla,
     &          flag_vsa, flag_sig, flag_alf
      use b2mod_input_profile
     & , only : ndata
      implicit none
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
! - shift_cf_data: integer flag to shift certain cf data
! - shift_value: actual value of the shift [mm]
! - shift: working array used to pass shift values for CF and optimization
! - corr_model: type of correlation model used in the covariance of each CF used in MAP
! - corr_length: value of the correlation length for each covariance matrix [mm]
! - corr_cutoff: percentage caut-off of the covariance matrix. values below this wrt to the diagonal are set to zero
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
      integer, save :: shift_cf_data(nsigmx) = 0
      real (kind=r8), save :: shift_value(nsigmx) = 0.0_R8
      real (kind=r8), save :: shift(nsigmx) = 0.0_R8
      integer, save :: cf_to_shift(nsigmx) = 0
      integer, save :: nshift = 0, shift_prior_type(nsigmx)
      real (kind=r8), save :: shift_prior_par(nsigmx,2), 
     & shift_prior_range(nsigmx,2)
      integer, save :: corr_model(nncf) = 0
      real (kind=r8), save :: corr_length(nncf) = 0.0_R8
      integer, save :: ncorr = 0, corr_prior_type(nncf)
      real (kind=r8), save :: corr_prior_par(nncf,2), 
     & corr_prior_range(nncf,2), corr_cutoff = 0.01
! VARIABLES RELATED TO OPTIMIZATION
! - nnvar: number of optimization parameters (not including radially varying transport coefficients)
! - npar_opt: actual number of optimization parameters (including spatially varying transport coeff.)
! - par_opt_phys: array containing the values of the optimization parameters used to pass them to the
!                 optimization libraries through b2optim_tao, and to the routine for
!                 calculating the prior (if needed)
! - nsigma_opt stores only the number of sigmas that are actually optimized (which can be .le. nsigma).
! - sigma_opt: logical that tells the specific sigmas that are optimized. This differentiation is needed to
!   correctly pass the optimization variables to the optimization libraries
! - nmean_opt, mean_opt, do the same as for sigma but for the error mean
! - nshift_opt, shift_opt: same but for cost function shift
! - ncorr_opt, corr_opt: same but for cost function shift
! - x0: initial guess of the optimization parameters (one for each actual parameter, so includes the
!        number of spatial points!)
! - xl/xu: lower/ipper bound of optimization parameters (one for each actual parameter, as x0)
! - shift_L/U: lower/upper bound of the cf data shift [mm]
! - shiftopt: has similar meaning as shift_opt but 'reduced' to the actual number of unique shifts (e.g.
!             if some CF have the same shift value because they are the same diagnostic)
! - shift/shiftLL/shiftUU: has similar meaning as shift_value/shift_L/shift_U but 'reduced' to the actual 
!                          number of unique shifts
! - corr_L/U: lower/upper bound of the correlation parameter
      integer, save :: npar_opt = 0
      real (kind=R8), allocatable, save :: par_opt0(:), par_opt(:),
     &       par_opt_phys(:), xsave(:)
      logical, save :: flag_optim = .false., reset_gradient = .false.
      integer nnvar
      real (kind=R8) :: x0(nvmx), xl(nvmx), xu(nvmx), gl(nvmx),
     &  gu(nvmx), jj(nvmx*nvmx),
     &  par_rescale(nvmx)
      integer :: jcol(nvmx*nvmx), jrow(nvmx*nvmx)
! sc some variables that allows to change optimization parameters, read from b2.optimization.parameters
      character*256, save :: limited_memory_update_type,
     &  hessian_approximation 
      real (kind=R8), save :: cpu_opt=0.0_R8, tol_opt=1.0e-7_R8
      integer, save :: maxiter=100, partype(nvmx), nsigma_opt = 0,
     &                 nmean_opt = 0, nshift_opt = 0, ncorr_opt = 0
      integer, save :: paris(nvmx), parib(nvmx)
      logical, save :: sigma_opt(nsigmx), mean_opt(nsigmx),
     &                 shift_opt(nsigmx), shiftopt(nsigmx),
     &                 corr_opt(nncf), parallel_hf
      logical, save :: spatial_dep(nvmx)
      integer, save :: spatial_points(nvmx) = 0
      real (kind=R8), save :: shift_L(nsigmx) = -1000.0_R8, ![envisaged default max 1m shift]
     &                        shift_U(nsigmx) = 1000.0_R8,
     &                        shiftLL(nsigmx) = -1000.0_R8, !copy of the previous ones
     &                        shiftUU(nsigmx) = 1000.0_R8,
     &                        corr_L(nncf) = -inf_opt*10._R8,
     &                        corr_U(nncf) = inf_opt*10._R8,
     &                        corr_rescale(nncf) = 1.0_R8
      integer :: idd
      logical :: file_ok
      real (kind=R8) :: rangeold(nvmx,2), parold(nvmx,2)
      integer :: typeold(nvmx)
      private optimization
      namelist /optimization/ nnvar, x0, xl, xu, 
     & maxiter, cpu_opt, tol_opt, hessian_approximation,
     & limited_memory_update_type, partype, sigma_opt, spatial_dep,
     & spatial_points, par_rescale, paris, parib, ncf,cfstart,cfend,
     & cftype, cfdef, cfweight,prior_type, prior_par, prior_range, 
     & scale_sigma, nsigma, sigma, cf_reg, cf_regP, mapToOMP, 
     & read_sigma, nmean, mean, mean_opt, shift_cf_data, shift_value,
     & shift_opt, shift_L, shift_U, shift_prior_type, shift_prior_par,
     & shift_prior_range, corr_model, corr_length, corr_prior_type,
     & corr_prior_range, corr_prior_par, corr_opt, corr_L, corr_U,
     & corr_cutoff, corr_rescale, parallel_hf

      contains

      subroutine read_b2mod_par_opt(ns, m, sw)
      use b2mod_types
      use b2mod_switches
      implicit none
      type(mapping), intent(inout) :: m
      integer, intent (in) :: ns
      type(switches), intent (in)   :: sw
      ! csc local variables
      integer :: ii, isigma, ipp, i, numdata, iss, indss, icf, icff,
     &           noss, ncffc, incf, idb, ic1, ic2, icv, ifc, ifcc,
     &           jj, imean, curr_ind
      integer, allocatable :: cfReg(:), shiftcfdata(:)
      logical :: done, optimize
      character*3 ss
      character*256 cffile, filename
      real(kind=r8) :: dummy
#ifdef LEGACYCOMP
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
      shift_prior_type = -1
      corr_prior_type = -1
      prior_par = -1.0_R8
      shift_prior_par = -1.0_R8
      corr_prior_par = -1.0_R8
! sc  Set prior range to (-inf,+inf)
      prior_range(:,1) = -inf_opt*10.0_R8
      prior_range(:,2) = inf_opt*10.0_R8
      shift_prior_range(:,1) = -inf_opt*10.0_R8
      shift_prior_range(:,2) = inf_opt*10.0_R8
      corr_prior_range(:,1) = -inf_opt*10.0_R8
      corr_prior_range(:,2) = inf_opt*10.0_R8
      cf_reg = 0
      cf_regP = 0
      mapToOMP = .false.
      x0 = inf_opt*10.0_R8
      xl = -inf_opt*10.0_R8
      xu = inf_opt*10.0_R8
      par_rescale = 1.0_R8
      partype = -2
      sigma_opt = .true.
      mean_opt = .false.
      spatial_dep = .false.
      spatial_points = 0
      parallel_hf = .true.
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
#ifndef LEGACYCOMP
        open (newunit=idb,file='debug_optimization.out')
#else
        idb = newunit()
        open (unit=idb,file='debug_optimization.out')
#endif
        call xertst(any(shift_cf_data.ge.0),
     &   'b2mod_par_opt: shift_cf_data must be >=0')
        incf = 0
        do icf = 1, ncf
          if (cftype(icf).gt.12) call xerrab ('cftype>12 not coded ')
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
     &        (cftype(icf).ge.7 .and. cftype(icf).le.10) .or.
     &        (cftype(icf).ge.12 .and. cftype(icf).le.12))
     >         cfread(icf) = .true.
          if (cftype(icf).eq.5 .or. cftype(icf).eq.11 .or. 
     &        cftype(icf).eq.12) m%cfOnCv(icf) = .false. !heat flux CF needs faces
          if ((cftype(icf).eq.5 .or. cftype(icf).eq.11 .or. 
     &         cftype(icf).eq.11) .and. cfdef(icf).eq.3) call xerrab(
     &       'cftype 5, 11 and 12 must be defined using faces!')
          if (shift_cf_data(icf).gt.0 .and. .not.cfread(icf)) then
            write(*,*) 'icf, shift_cf_data(icf), cfread(icf),'//
     &        ' cftype(icf)'
            write(*,*) icf, shift_cf_data(icf), cfread(icf), cftype(icf)
            call xerrab('b2mod_par_opt: cannot shift a CF without data')
          endif
          if (shift_cf_data(icf).eq.0 .and. shift_opt(icf)) then
            write(*,*) 'icf, shift_cf_data(icf), shift_opt(icf)'
            write(*,*) icf, shift_cf_data(icf), shift_opt(icf)
            call xerrab('b2mod_par_opt: shift_opt=true requires '//
     &        'shift_cf_data>0')
          endif
          call xertst(corr_model(icf).ge.0 .and. corr_model(icf).le.1,
     &      'b2mod_par_opt: corr_model can only be 0 or 1')
          if (corr_model(icf).gt.0 .and. .not.cfread(icf)) then
            write(*,*) 'icf, corr_model(icf), cfread(icf),'//
     &        ' cftype(icf)'
            write(*,*) icf, corr_model(icf), cfread(icf), cftype(icf)
            call xerrab('b2mod_par_opt: cannot use a corr_model in CF'//
     &       ' without data')
          endif
          if (corr_model(icf).eq.0 .and. corr_opt(icf)) then
            write(*,*) 'icf, corr_model(icf), corr_opt(icf)'
            write(*,*) icf, corr_model(icf), corr_opt(icf)
            call xerrab('b2mod_par_opt: corr_opt=true requires '//
     &        'corr_model>0')
          endif
!csc reading files with CF data in it. From now on we assume data is stored as
!    1D radial profiles, and the abscissa is also stricly monotonically increasing. 
!    This facilitates shifting the CF and covariance matrix caclulation 
          if (cfread(icf)) then
            write (ss,'(I1)') icf
            if (icf.gt.9) write (ss,'(I2)') icf
            cffile = 'cf'//trim(ss)//'.dat'
            call find_file(cffile,file_ok)
            if (.not. file_ok) call xerrab
     &       ('Cannot find cost function file '//trim(cffile))
            open(unit=99,file=cffile,status='old',action='read')
            write(*,*) 'Reading cost function data from ',trim(cffile)
            ncfdata(icf) = 0
!csc File is of unknown length, so first count the number of lines
 5          read(99,*,end=15) ss
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
! Check if strictly monotonically increasing (sorted)
            dummy = cfdata(icf,1,1)
            do i=2,ncfdata(icf)
              if (.not.(cfdata(icf,1,i).gt.dummy)) then
                write(*,*) 'icf=',icf, ' ii=',i,' cfdata(icf,1,ii-1)=',
     &           dummy,' cfdata(icf,1,ii)=',cfdata(icf,1,i)
                call xerrab('b2mod_par_opt: abscissa in cfdata not '//
     &                      'stricly increasing!')
              endif
              dummy = cfdata(icf,1,i)
            enddo
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
            numdata = 0
            do iCv=1,nomp
              if (omp(iCv)>m%nCi) cycle !only use internal CVs
              m%cfReg(incf +1+numdata) = omp(iCv)
              numdata = numdata + 1
            enddo
            m%cfRegP(icf,2) = numdata
            incf = incf + numdata

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
            write(ss,"(I3)") icf
            call xertst(ncffc.gt.0, 
     &        'No faces found for cost function  = '//ss)
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
        !csc more testing and array preparation for shifting CF
        if (cftype(1).eq.6 .or. cftype(1).eq.0) call xertst(
     &    shift_cf_data(1).eq.0,'b2mod_par_opt: cftype 0 and 6 '//
     &    'cannot have a shift')
        nshift = 0
        allocate(shiftcfdata(nsigmx)) !working array
        shiftcfdata = shift_cf_data
        curr_ind = 0
        optimize = .false.
        shiftopt = .false.
        do icf=1, ncf
          if (shiftcfdata(icf).eq.0) cycle
          nshift = nshift + 1
          curr_ind = shiftcfdata(icf)
          shift(nshift) = shift_value(icf)
          cf_to_shift(icf) = nshift
          shiftLL(nshift) = shift_L(icf)
          shiftUU(nshift) = shift_U(icf)
          optimize = shift_opt(icf)
          if (optimize) shiftopt(nshift) = .true.
          ! search for all (next) CF that have same index
          do icff = icf+1, ncf
            if (shiftcfdata(icff).eq.curr_ind) then
              ! a cost function has the same shift
              cf_to_shift(icff) = nshift
              shiftcfdata(icff) = 0 ! set to zero to avoid duplication
              if (.not. (optimize .eqv. shift_opt(icff))) then
                write(*,*) 'ERROR: shift_opt inconsistent for '//
     &            'shift_cf_data=',curr_ind
                call xerrab('b2mod_par_opt: shift_opt inconsistent')
              endif
            endif
          end do
        end do !ncf
        if (allocated(shiftcfdata)) deallocate(shiftcfdata)
        ! csc some tests on corr_model
        ncorr = 0
        call xertst(corr_model(1).eq.0,'b2mod_par_opt: first cost '//
     &    'function cannot have a correlation model')
        do icf = 1, ncf
          if (corr_model(icf).gt.0) ncorr = ncorr+1
          if (corr_model(icf).eq.1)
     &      call xertst(corr_length(icf) .gt. 0.0_R8,
     &      'b2mod_par_opt: corr_length must be > 0.0 for corr_model=1')
        end do
      end if !if ncf>0

! now allocate vector that will be used to interpolate SOLPS onto data grid
      if (any(cfread)) then
        numdata = maxval(ncfdata(1:ncf),1)
        allocate(b2volOnCf(ncf,numdata)) !stores interpolated B2.5 volumes onto CF radial points
        allocate(b2dataOnCf(numdata)) !variable to store SOLPS data interpolated onto CF radial points
        numdata = maxval(m%cfRegP(1:ncf,2)) !max number of CVs in a cost function
        numdata = max(numdata,nomp) ! make sure this is the case when mapping to omp
        allocate(b2rr(ncf,numdata)) !store here radial distance coordinate of SOLPS data
        allocate(b2data(numdata)) !temporary variable to store SOLPS data for interpolation
        allocate(b2psi(ncf,numdata)) !store here psi of SOLPS data for interpolation
        b2rr = 0.0_R8
        b2psi = 0.0_R8
        b2volOnCf = 0.0_R8
        b2data = 0.0
        b2dataOnCf = 0.0
      endif

! csc Tests on optimization parameters
      call xertst(nnvar.ge.0.and.nnvar.le.nvmx,
     & 'b2mod_par_opt: nnvar<0 or nnvar>nvmx')
      call xertst(all(par_rescale.gt.0.0_R8),
     & 'b2mod_par_opt: par_rescale<=0')
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
! check that b2tqna_inputfile is used if any space-dependent transport sensitivity is needed
       if (any(spatial_dep(1:nnvar))) call xertst(
     &  sw%b2tqna_transport_inputfile.eq.1, 'Sensitivity of spatially'//
     &  ' dependent coefficients needs b2tqna_inputfile=1!')
       do ii=1,nnvar
! check that the corresponding transport coefficients are actually active before
! evaluating any sensitivity!!
        if ((partype(ii) .eq. 1) .and. .not.spatial_dep(ii)) 
     &    call xertst(flag_dna.eq.1, 
     &     'Sensitivity for parm_dna requires flag_dna=1!')
        if ((partype(ii) .eq. 2) .and. .not.spatial_dep(ii)) 
     &    call xertst(flag_dpa.eq.1, 
     &     'Sensitivity for parm_dpa requires flag_dpa=1!')
        if ((partype(ii) .eq. 3) .and. .not.spatial_dep(ii)) 
     &    call xertst(flag_hci.eq.1, 
     &     'Sensitivity for parm_hci requires flag_hci=1!')
        if ((partype(ii) .eq. 4) .and. .not.spatial_dep(ii)) 
     &    call xertst(flag_hce.eq.1, 
     &     'Sensitivity for parm_hce requires flag_hce=1!')
        if ((partype(ii) .eq. 6) .and. .not.spatial_dep(ii)) 
     &    call xertst(flag_vla.eq.1, 
     &     'Sensitivity for parm_vla requires flag_vla=1!')
        if ((partype(ii) .eq. 7) .and. .not.spatial_dep(ii)) 
     &    call xertst(flag_vsa.eq.1, 
     &     'Sensitivity for parm_vsa requires flag_vsa=1!')
        if ((partype(ii) .eq. 8) .and. .not.spatial_dep(ii)) 
     &    call xertst(flag_sig.eq.1, 
     &     'Sensitivity for parm_sig requires flag_sig=1!')
        if ((partype(ii) .eq. 9) .and. .not.spatial_dep(ii)) 
     &    call xertst(flag_alf.eq.1, 
     &     'Sensitivity for parm_alf requires flag_alf=1!')
! now some checks on spatially dependent coefficients
        if (spatial_dep(ii)) then
        !only parm_XXX can be spatially dependent
         call xertst(partype(ii).ge.1 .and. partype(ii).le.9, 
     &     'b2mod_par_opt: spatial_points optimization '//
     &     'is available only for partype 1-9')
         if (spatial_points(ii).le.1) then
          write(*,*) 'ipar, spatial_dep, spatial_points=',
     &      ii, spatial_dep(ii), spatial_points(ii)
          call xerrab(' b2mod_par_opt: seems you want to optimize '//
     &     'spatially dependent transport coefficients but you only '//
     &     'specify one!')
         endif
         if (spatial_points(ii).ne.ndata(1,partype(ii),paris(ii))) then
          write(*,*) 'ipar, spatial_dep, spatial_points, ndata=',
     &      ii, spatial_dep(ii), spatial_points(ii), 
     &       ndata(1,partype(ii),paris(ii))
          call xerrab(' b2mod_par_opt: spatial_points in '//
     &     'b2.optimization.parameters must be equal to ndata '//
     &     'in b2.transport.inputfile')
         endif
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

!      Now do some tests on priors for physical parameters (including sigma and mean for now)
       if (nnvar.gt.0) then
         write(*,*) 'Testing priors for physical parameters,sigma,mean'
         call test_priors(nnvar, prior_type(1:nnvar),
     &    prior_par(1:nnvar,1:2), prior_range(1:nnvar,1:2), xl(1:nnvar))
       endif
      end if  !cftype(1).eq.6

!     adjust some arrays when there are radially dependent coefficients
!     x0, xl, xu are supposed to have the correct dimension already, i.e. one for each spatial point
!     priors parameters are not
      npar_opt = 0
      idd = 1
      typeold = prior_type
      parold = prior_par
      rangeold = prior_range
      do ipp = 1, nnvar
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
!     similarly add the optimization part for the shifting of cfdata, not included yet
!     here X0 for shift is taken directly from the array shift, and same for lower/upper bounds
! FIXME generalize in a smilar way for mean and sigma, without need to specify partype etc. for them
      nshift_opt = 0
      do ii = 1, nshift
        if (shiftopt(ii)) then
          nshift_opt = nshift_opt+1
          call xertst(npar_opt.le.nvmx,
     &     'b2mod_par_opt: increase size of nvmx')
          npar_opt = npar_opt + 1
          x0(npar_opt) = shift(ii)
          xL(npar_opt) = shiftLL(ii)
          xU(npar_opt) = shiftUU(ii)
          par_rescale(npar_opt) = 1.0_R8 !fixed rescaling for shift, not expecting huge values
        endif
      enddo
      if (cftype(1).eq.6 .and. nshift_opt.gt.0) then
        ! test priors for shift. sigma and mean will be added here later
        write(*,*) 'Testing priors for shift_cf_data'
        call test_priors(nshift_opt, shift_prior_type(1:nshift_opt), 
     &   shift_prior_par(1:nshift_opt,1:2),
     &   shift_prior_range(1:nshift_opt,1:2), 
     &   xL(npar_opt-nshift_opt+1:npar_opt))
      endif
!     as above, add the contribution to sensitivity/optimization variables of the corr_model
      ncorr_opt = 0
      call xertst(.not.corr_opt(1),'b2mod_par_opt: corr_opt cannot '//
     & 'be .true. for icf=1')
      do ii=1, ncf
        if (corr_opt(ii)) then
          ncorr_opt = ncorr_opt+1
          call xertst(npar_opt.le.nvmx,
     &     'b2mod_par_opt: increase size of nvmx')
          npar_opt = npar_opt + 1
          x0(npar_opt) = corr_length(ii)
          if (.not. (corr_L(ii).gt.0.0_r8)) then
             write(*,*) 'icf=',icf,', corr_L(icf)=',corr_L(icf)
            call xerrab('b2mod_par_opt: corr_L must be >0.0')
          endif
          xL(npar_opt) = corr_L(ii)
          xU(npar_opt) = corr_U(ii)
          par_rescale(npar_opt) = corr_rescale(ii)
        endif
      enddo
      if (cftype(1).eq.6 .and. ncorr_opt.gt.0) then
        ! test priors for corr
        write(*,*) 'Testing priors for corr_length'
        call test_priors(ncorr_opt, corr_prior_type(1:ncorr_opt), 
     &   corr_prior_par(1:ncorr_opt,1:2),
     &   corr_prior_range(1:ncorr_opt,1:2), 
     &   xL(npar_opt-ncorr_opt+1:npar_opt))
      endif
      if (flag_optim) call xertst(npar_opt.gt.0,
     & 'b2mod_par_opt: cannot optimize if npar_opt=0!')
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
      isigma = npar_opt - nsigma_opt  - nmean_opt - nshift_opt - 
     &         ncorr_opt + 1 ! position of first sigma in x0
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
      imean = npar_opt - nmean_opt - nshift_opt - ncorr_opt + 1 ! position of first mean in x0
      do ii = 1, nmean
!      only assign initial value from x0 if that mean is being optimized or used in MAP
       if (mean_opt(ii)) then
        mean(ii) = x0(imean)
        imean = imean + 1 !next mean in inital vector x0
       endif
      end do
      write (*,*) 'NPAR_OPT = ',npar_opt

      return
      end subroutine read_b2mod_par_opt

      subroutine dealloc_b2mod_par_opt
      implicit none
      if (allocated(b2rr)) deallocate(b2rr)
      if (allocated(b2volOnCf)) deallocate(b2volOnCf)
      if (allocated(b2data)) deallocate(b2data)
      if (allocated(b2dataOnCf)) deallocate(b2dataOnCf)
      if (allocated(b2psi)) deallocate(b2psi)
      end subroutine dealloc_b2mod_par_opt

      subroutine test_priors(nn, pr_type, pr_par, pr_range, lowerB)
      implicit none
      integer, intent (IN) :: nn, pr_type(nn)
      real(kind=r8), intent(in) :: pr_par(nn,2), pr_range(nn,2),
     &  lowerB(nn)
      integer :: ii
!      prior_type refers to the prior type for the design variables:
!      0 - Uniform distribution.
!      1 - Uninformative, proper, Gaussian prior.
!      2 - Gamma distribution. (defined through mean and std)
!      3 - Jeffrey's 1/sigma (only for sigma!!)
       call xertst(any(pr_type(1:nn) .ge. 0).and.
     &  any(pr_type(1:nn) .le. 3),
     &  'b2mod_par_opt: prior_type<0 or prior_type>3')
!      prior_par refers to the prior parameters:
!      For prior_type = 1,2, prior_par(ii,1) = mu, prior_par(ii,2) = sigma
!      For prior_type = 3, not used?
       do ii = 1, nn
        if (pr_type(ii) .eq. 1) 
     &   call xertst(pr_par(ii,2) .gt. 0.0_R8,
     &    'b2mod_par_opt: sigma prior <=0')
        if (pr_type(ii) .eq. 2) then
         call xertst((pr_par(ii,1) .gt. 0.0_R8) .and. 
     &               (pr_par(ii,2) .gt. 0.0_R8),
     &    'b2mod_par_opt: mean and variance for gamma prior must be >0')
         if (max(pr_range(ii,1),lowerB(ii)).le.0.0_R8) then
           write(*,*) 'ii, prior_range(ii,1)',ii, pr_range(ii,1)
           write(*,*) 'ii, xl(ii)',ii, lowerB(ii)
           call xerrab('b2mod_par_opt: gamma distribution is '//
     &     'intended for parameters on domain (0 +inf)')
         endif
        endif
        call xertst(pr_range(ii,1).lt.pr_range(ii,2),
     &   'b2mod_par_opt: prior_range(.,1)>prior_range(.,2)')
       end do
      end subroutine test_priors

      end module b2mod_par_opt
