      program b2optim_ipopt
      use b2mod_types
      use b2mod_par_opt_diff
      use b2mod_main_diff
     , , only : b2mn_init_d, b2mn_step_d, b2mn_step, b2mn_fin_d
      use b2mod_ad_diff
     , , only : nncf
      use b2mod_user_namelist_diff
     , , only : nsigma, sigma
      use b2mod_transport_namelist_diff
      use b2mod_input_profile_diff
      use b2mod_boundary_namelist_diff
      use b2mod_switches_diff
      implicit none

!     include the Ipopt return codes
      include 'IpReturnCodes.inc'

!     Size of the problem (number of variables and equality constraints)
      integer :: ncon, nele_jac, nele_hess, IDX_STY
      parameter  (IDX_STY = 1 )
!sc N = number design variables
!   M number of constraints (it's 2 because the third one is the variable bound!
!   NELE_JAC # of nonzeros in Jacobian of Lagrangian
!   IDX_STY index style for matrix (1 for fortran normally)

!     Space for multipliers and constraints
      double precision, allocatable :: LAM(:), G(:)

!     Vector of variables
      double precision, allocatable :: X(:)

!     Vector of lower and upper bounds
      double precision, allocatable :: X_L(:), X_U(:), Z_L(:), Z_U(:)
      double precision, allocatable :: G_L(:), G_U(:)
!sc  X_L/X_U lower upper bound for each desiign variable
!    Z_L, Z_U initial value for the bound multiplier --> seems not needed to be initialized here
!    G_L,G_U lower/upper bounds for constraints

!     Private data for evaluation routines
!     This could be used to pass double precision and integer arrays untouched
!     to the evaluation subroutines EVAL_*
      double precision DAT(2)
      integer IDAT(1)

!     Place for storing the Ipopt Problem Handle
      integer*8 IPROBLEM
      integer*8 IPCREATE
      integer IERR
      integer IPSOLVE, IPADDSTROPTION
      integer IPADDNUMOPTION, IPADDINTOPTION
      integer IPOPENOUTPUTFILE
      double precision F
      integer ipar, icon, isigma, idir, i

      type(switches_diffv), save :: switchdiff

!     The following are the Fortran routines for computing the model
!     functions and their derivatives - their code can be found further
!     down in this file.
      external streql
      logical streql, hessian
!sc   these are the external subroutine names passed to IPOPT


!     ..START OF COMPUTATION

!     Allocate and initialize par_opt variables to be used in B2.5

      call b2mn_init_d
      call read_b2mod_par_opt(ncon, nele_jac)
      allocate(par_opt_phys(npar_opt))
!     Initialize derivatives for diffusion coefficients
!     FIXME initialization to 0.0_R8 should be done elsewhere?
!     FIXME if not differentiated wrt certain variables, then their derivative here will not exist and such the lines should be manually commented-out
!     FIXME remove from here for adjoint
      tdatad = 0.0_R8
      parm_dnad = 0.0_R8
      parm_dpad = 0.0_R8
      parm_hcid = 0.0_R8
      parm_hced = 0.0_R8
      parm_vlad = 0.0_R8
      parm_vsad = 0.0_R8
      parm_sigd = 0.0_R8
      parm_alfd = 0.0_R8
      enkpard = 0.0_R8
      switchdiff%b2sikt_fac_sheath = 0.0_R8
      switchdiff%b2sikt_fac_sheath_core = 0.0_R8
      switchdiff%keps_cd = 0.0_R8
      switchdiff%keps_heat = 0.0_R8
      switchdiff%keps_heat_i = 0.0_R8
      idir = 1 !this indicates the different directions of multidirectional derivative mode
      do ipar = 1, nnvar - nsigma_opt
        if (spatial_dep(ipar)) then
!         spatially dependent coefficient
          do i = 1, spatial_points(ipar)
            tdatad(idir,2,i,partype(ipar),1) = 1.0_R8
            idir = idir + 1
          end do
        else
!         non spatially dependent --> specified in b2tqna
          select case (partype(ipar))
          case (1) ! dna
            parm_dnad(idir,1) = 1.0_R8 ! FIXME to improve for multispecies
          case (2) ! dpa
            parm_dpad(idir,1) = 1.0_R8 ! FIXME to improve for multispecies
          case (3) ! hci
            parm_hcid(idir,1) = 1.0_R8 ! FIXME to improve for multispecies
          case (4) ! hce
            parm_hced(idir) = 1.0_R8
          case (5) ! vla
            parm_vlad(idir,1) = 1.0_R8
          case (7) ! vsa
            parm_vsad(idir,1) = 1.0_R8
          case (8) ! sig
            parm_sigd(idir) = 1.0_R8
          case (9) ! alf
            parm_alfd(idir) = 1.0_R8
          case (10) ! enkpar(1,1)
            enkpard(idir,1,1) = 1.0_R8
          case (11) ! b2sikt_fac_sheath
            switchdiff%b2sikt_fac_sheath(idir) = 1.0_R8
          case (12) ! b2sikt_fac_sheath_core
            switchdiff%b2sikt_fac_sheath_core(idir) = 1.0_R8
          case (13) ! keps_cd
            switchdiff%keps_cd(idir) = 1.0_R8
          case (14) ! keps_heat
            switchdiff%keps_heat(idir) = 1.0_R8
          case (15) ! keps_heat_i
            switchdiff%keps_heat_i(idir) = 1.0_R8
          case default
            write(*,*) partype(ipar)
            call xerrab ('partype out of bounds')
          end select
          idir = idir + 1
        endif
      end do
      if (nsigma_opt.gt.0) then
        sigmad = 0.0_R8
        isigma = npar_opt - nsigma_opt + 1
        do ipar = 1, nsigma
!         only if sigma is being optimized!
          if (sigma_opt(ipar)) then
            sigmad(isigma,ipar) = 1.0_R8
            isigma = isigma + 1
          endif
        end do
      endif
!     FIXME remove until here for adjoint
      par_opt_phys = 0.0_R8
      flag_optim  = .true. !csc this will tell in b2tqna to use par_opt_phys for parm_dna

      write (*,*) 'IPOPT, WARNING: setting only 1 hard-coded '//
     & 'trivial constraint'
      ncon = 1 ! FIXME to be updated whenever nonlinear or state constraints are needed
      nele_jac = 1

!     For the moment only steepest descent is available when specifiyng the Hessian
      hessian = .false.
      nele_hess = 0
      if (streql('exact',hessian_approximation)) then
        hessian = .true.
        write (*,*) 'IPOPT, WARNING: setting Hessian as '//
     &   'identity matrix for Steepest Descent method'
        nele_hess = npar_opt
      endif

      allocate(X(npar_opt),X_L(npar_opt),X_U(npar_opt))
      allocate(Z_L(npar_opt),Z_U(npar_opt))
      allocate(G(ncon),G_L(ncon),G_U(ncon),LAM(ncon))

      do ipar = 1, npar_opt
        X(ipar) = x0(ipar)
        X_L(ipar) = xl(ipar)
        X_U(ipar) = xu(ipar)
        if (X_L(ipar).lt.-inf_opt)
     &    write(*,*) 'IPOPT: warning, X_L(',ipar,') set to infty'
        if (X_U(ipar).gt.inf_opt)
     &    write(*,*) 'IPOPT: warning, X_U(',ipar,') set to infty'
      end do
      do icon = 1, ncon
        G_L(icon) = -inf_opt*10!gl(icon)
        G_U(icon) = 0.0!gu(icon)
        if (G_L(icon).lt.-inf_opt)
     &    write(*,*) 'IPOPT: warning, G_L(',icon,') set to infty'
        if (G_U(icon).gt.inf_opt)
     &    write(*,*) 'IPOPT: warning, G_U(',icon,') set to infty'
      end do
! FIXME Z_L/Z_U


!     First create a handle for the Ipopt problem (and read the options file)
!     If Hessian is provided or steepest descent is asked the additional argument needed
      if (hessian) then
        IPROBLEM = IPCREATE(npar_opt, X_L, X_U, ncon, G_L, G_U,
     &     NELE_JAC, NELE_HESS, IDX_STY, EV_F, EV_G, EV_GRAD_F,
     &     EV_JAC_G, EV_HESS)
      else
        IPROBLEM = IPCREATE(npar_opt, X_L, X_U, ncon, G_L, G_U,
     &      NELE_JAC, NELE_HESS, IDX_STY, EV_F, EV_G, EV_GRAD_F,
     &     EV_JAC_G)
      endif
      if (IPROBLEM.eq.0) then
         write(*,*) 'Error creating an Ipopt Problem handle.'
         stop
      endif

!     Open an output file
      IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 5)
      if (IERR.ne.0 ) then
         write(*,*) 'Error opening the Ipopt output file.'
         goto 9000
      endif

      IERR = IPADDSTROPTION(IPROBLEM, 'print_user_options', 'yes')
      if (IERR.ne.0 ) goto 9990

!     Set mu strategy
      IERR = IPADDSTROPTION(IPROBLEM, 'mu_strategy', 'adaptive')
      if (IERR.ne.0 ) goto 9990

!     Set scaling methods
!      IERR = IPADDSTROPTION(IPROBLEM, 'nlp_scaling_method',
!     & 'user-scaling')
!      if (IERR.ne.0 ) goto 9990

!sc IMPORTANT HERE these are needed for BFGS updates without hessian
      IERR = IPADDSTROPTION(IPROBLEM, 'hessian_approximation',
     & hessian_approximation)
      if (IERR.ne.0 ) goto 9990
      if (.not.hessian) then
        IERR = IPADDSTROPTION(IPROBLEM, 'limited_memory_update_type',
     &   limited_memory_update_type)
        if (IERR.ne.0 ) goto 9990
      endif

!     Set max number of iterations
      call xertst (maxiter.gt.0, 'faulty internal parameter maxiter')
      IERR = IPADDINTOPTION(IPROBLEM, 'max_iter', maxiter)
      if (IERR.ne.0 ) goto 9990

!      IERR = IPADDINTOPTION(IPROBLEM, 'print_level', 12)
!      if (IERR.ne.0 ) goto 9990

!     Set a minimum tolerance
      call xertst (tol_opt.gt.0.0,'faulty internal parameter tol_opt')
      IERR = IPADDNUMOPTION(IPROBLEM, 'tol', tol_opt)
      if (IERR.ne.0 ) goto 9990

!     Set a maximum cpu time for ipopt based on b2mndr_elapsed
      call xertst (cpu_opt.ge.0.0_R8, 
     & 'faulty internal parameter cpu_opt')
      if (cpu_opt.gt.0.0_R8) then
        IERR = IPADDNUMOPTION(IPROBLEM, 'max_cpu_time', cpu_opt)
        if (IERR.ne.0 ) goto 9990
      endif

!     Set infinity bounds
      IERR = IPADDNUMOPTION(IPROBLEM, 'nlp_upper_bound_inf', 
     & inf_opt*0.1)
      IERR = IPADDNUMOPTION(IPROBLEM, 'nlp_lower_bound_inf', 
     & -inf_opt*0.1)
      if (IERR.ne.0 ) goto 9990

!     Set a callback function to give you control once per iteration.
!     You can use it if you want to generate some output, or to stop
!     the optimization early.
      call IPSETCALLBACK(IPROBLEM, ITER_CB)

!     As a simple example, we pass the constants in the constraints to
!     the EVAL_C routine via the "private" DAT array.
      DAT(1) = 0.d0
      DAT(2) = 0.d0

!     Call optimization routine
      IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)

!     Output:
      if (IERR.eq.IP_SOLVE_SUCCEEDED) then
        write(*,*)
        write(*,*) 'Ipopt: the solution was found.'
        write(*,*)
      elseif (IERR.eq.IP_ITERATION_EXCEEDED) then
        write(*,*)
        write(*,*) 'Ipopt: reached maximum number of iterations.'
        write(*,*)
      elseif (IERR.eq.IP_CPUTIME_EXCEEDED) then
        write(*,*)
        write(*,*) 'Ipopt: reached maximum CPU time allowed.'
        write(*,*)
      else
        write(*,*)
        write(*,*) 'Ipopt: an error occoured.'
        write(*,*) 'The error code is ',IERR
        write(*,*)
      endif
      write(*,*)
      write(*,*) 'The final value of the objective function is ',F
      write(*,*)
      write(*,*) 'The optimal values of X are:'
      write(*,*)
      do ipar = 1, npar_opt
        write(*,*) 'X  (',ipar,') = ',X(ipar)
      enddo
      write(*,*)
      write(*,*) 'The multipliers for the lower bounds are:'
      write(*,*)
      do ipar = 1, npar_opt
        write(*,*) 'Z_L(',ipar,') = ',Z_L(ipar)
      enddo
      write(*,*)
      write(*,*) 'The multipliers for the upper bounds are:'
      write(*,*)
      do ipar = 1, npar_opt
        write(*,*) 'Z_U(',ipar,') = ',Z_U(ipar)
      enddo
      write(*,*)
      write(*,*) 'The multipliers for the equality constraints are:'
      write(*,*)
      do icon = 1, ncon
        write(*,*) 'LAM(',icon,') = ',LAM(icon)
      enddo
      write(*,*)
 9000 continue

!     Clean up
      call b2mn_fin_d
      call IPFREE(IPROBLEM)
      deallocate(X,X_L,X_U,Z_L,Z_U)
      deallocate(G,G_L,G_U,LAM)
      deallocate(par_opt_phys)
      stop 'b2optim'
 9990 continue
      write(*,*) 'Error setting an option'
      goto 9000

      contains

      subroutine EV_F(npar, XX, NEW_X, F, IDAT, DAT, IER)
      implicit none
      integer npar, NEW_X
      double precision F, XX(npar)
      real(kind=r8) j(nncf)
      double precision DAT(*)
      integer IDAT(*)
      integer IER, ipar, isigma
      character*3 str
      save j
!csc important here to check for NEW_X! so that it doesn't recalculate with same X
      if (new_x.eq.0) then
        do ipar = 1, npar - nsigma_opt
          par_opt_phys(ipar) = real(XX(ipar),kind=r8)
          write(str,"(I1)") ipar
          if (ipar.ge.10) write(str,"(I2)") ipar
          write(*,*) 'IPOPT: eval_f with x',trim(str),'= ',
     &     par_opt_phys(ipar)
        enddo
        isigma = npar - nsigma_opt + 1
        do ipar = 1, nsigma
          if (sigma_opt(ipar)) then
            sigma(ipar) = real(XX(isigma),kind=r8)
            write(str,"(I1)") isigma
            if (ipar.ge.10) write(str,"(I2)") isigma
            write(*,*) 'IPOPT: eval_f with x',trim(str),'= ',
     &       sigma(ipar)
            isigma = isigma + 1
          endif
        enddo
        call b2mn_step(j)
      endif
      F = DBLE(j(1))
      IER = 0
      return
      end subroutine EV_F

      subroutine EV_GRAD_F(npar, XX, NEW_X, GRAD, IDAT, DAT, IER)
      implicit none
      integer npar, NEW_X
      double precision GRAD(npar), XX(npar)
      real(kind=r8) j(nncf), jd(nncf)
      double precision DAT(*)
      integer IDAT(*)
      integer IER, ipar, isigma, idir, i
      character*3 str
      save jd
!csc important here to check for NEW_X! so that it doesn't recalculate with same X
! FIXME automatically choose with preprocessor data
      if (new_x.eq.0) then
        do ipar = 1, npar - nsigma_opt
          par_opt_phys(ipar) = real(XX(ipar),kind=r8)
          write(str,"(I1)") ipar
          if (ipar.ge.10) write(str,"(I2)") ipar
          write(*,*) 'IPOPT: eval_grad_f with x',trim(str),'= ',
     &     par_opt_phys(ipar)
        enddo
        isigma = npar - nsigma_opt + 1
        do ipar = 1, nsigma
          if (sigma_opt(ipar)) then
            sigma(ipar) = real(XX(isigma),kind=r8)
            write(str,"(I1)") isigma
            if (ipar.ge.10) write(str,"(I2)") isigma
            write(*,*) 'IPOPT: eval_grad_f with x',trim(str),'= ',
     &       sigma(ipar)
            isigma = isigma + 1
          endif
        enddo
        call b2mn_step_d(j,jd,switchdiff)
      endif
!     FIXME remove from here for adjoint
      do ipar = 1, npar
        grad(ipar) = DBLE(jd(1))
      end do
!     FIXME remove until here for adjoint
!     FIXME remove from here for tangent
      idir = 1 !this indicates the different directions of adjoint mode
      do ipar = 1, nnvar - nsigma_opt
        if (spatial_dep(ipar)) then
!         spatially dependent coefficient
          do i = 1, spatial_points(ipar)
            grad(idir) = DBLE(tdatab(2,i,partype(ipar),1))
            idir = idir + 1
          end do
        else
!         non spatially dependent --> specified in b2tqna
          select case (partype(ipar))
          case (1) ! dna
            grad(idir) = DBLE(parm_dnab(1)) ! FIXME to improve for multispecies
          case (2) ! dpa
            grad(idir) = DBLE(parm_dpab(1)) ! FIXME to improve for multispecies
          case (3) ! hci
            grad(idir) = DBLE(parm_hcib(1)) ! FIXME to improve for multispecies
          case (4) ! hce
            grad(idir) = DBLE(parm_hceb)
          case (5) ! vla
            grad(idir) = DBLE(parm_vlab(1))
          case (7) ! vsa
            grad(idir) = DBLE(parm_vsab(1))
          case (8) ! sig
            grad(idir) = DBLE(parm_sigb)
          case (9) ! alf
            grad(idir) = DBLE(parm_alfb)
          case (10) ! enkpar(1,1)
            grad(idir) = DBLE(enkparb(1,1))
          case (11) ! b2sikt_fac_sheath
            grad(idir) = DBLE(switchdiff%b2sikt_fac_sheath)
          case (12) ! b2sikt_fac_sheath_core
            grad(idir) = DBLE(switchdiff%b2sikt_fac_sheath_core)
          case (13) ! keps_cd
            grad(idir) = DBLE(switchdiff%keps_cd)
          case (14) ! keps_heat
            grad(idir) = DBLE(switchdiff%keps_heat)
          case (15) ! keps_heat_i
            grad(idir) = DBLE(switchdiff%keps_heat_i)
          case default
            write(*,*) partype(ipar)
            call xerrab ('partype out of bounds')
          end select
          idir = idir + 1
        endif
      end do
      if (nsigma_opt.gt.0) then
        do ipar = 1, nsigma
!         only if sigma is being optimized!
          if (sigma_opt(ipar)) then
            grad(idir) = DBLE(sigmab(ipar))
            idir = idir + 1
          endif
        end do
      endif
!     FIXME remove until here for tangent

      IER = 0
      return
      end subroutine EV_GRAD_F

      subroutine EV_G(npar, X, NEW_X, ncon, G, IDAT, DAT, IERR)
      implicit none
      integer npar, NEW_X, ncon
      double precision G(ncon), X(npar)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR, icon
!     s.t.  x1*x2*x3*x4                   >=  25
!           x1**2 + x2**2 + x3**2 + x4**2  =  40
!csc WARNING to be updated whenever nonlinear constraints are needed
      G(1) = -1.0

      IERR = 0
      return
      end subroutine EV_G

      subroutine EV_JAC_G(TASK, npar, X, NEW_X, ncon, NZ, ACON, AVAR, A,
     &    IDAT, DAT, IERR)
      integer TASK, npar, NEW_X, ncon, NZ
      double precision X(npar), A(NZ)
      integer ACON(NZ), AVAR(NZ), I
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
!     structure of Jacobian:
!
!     s.t.  x1*x2*x3*x4                   >=  25
!           x1**2 + x2**2 + x3**2 + x4**2  =  40
!csc  [x2*x3*x4   x1*x3*x4   x1*x2*x4   x1*x2*x3;
!         2*x1       2*x2       2*x3       2*x4]
      integer AVAR1(1), ACON1(1)
      data  AVAR1 /1/ !ROWS
      data  ACON1 /1/ !COLUMNS
      save  AVAR1, ACON1
!csc WARNING to be updated whenever nonlinear or state constraints are needed
      if( TASK.eq.0 ) then
        do I = 1, NZ
          AVAR(I) = 1!jrow(i)
          ACON(I) = 1!jcol(i)
        enddo
      else
        do I = 1, NZ
          A(I) = 0.0! FIXME jj(i) ![1,1]
        end do
      endif
      IERR = 0
      return
      end subroutine EV_JAC_G

      subroutine EV_HESS(TASK, npar, X, NEW_X, OBJFACT, ncon, LAM,
     &     NEW_LAM, NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
      implicit none
      integer TASK, npar, NEW_X, ncon, NEW_LAM, NNZH, i
      double precision X(npar), OBJFACT, LAM(ncon), HESS(NNZH)
      integer IRNH(NNZH), ICNH(NNZH)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
C
C     structure of Hessian:
C
C     L = F + lambda*h

c      integer IRNH1(10), ICNH1(10)
c      data  IRNH1 /1, 2, 2, 3, 3, 3, 4, 4, 4, 4/
c      data  ICNH1 /1, 1, 2, 1, 2, 3, 1, 2, 3, 4/
c      save  IRNH1, ICNH1
C     1  2  3  4
C   1 *
C   2    *
C   3       *
C   4          * 
      if( TASK.eq.0 ) then
         do i = 1, NNZH
            IRNH(i) = i
            ICNH(i) = i
         enddo
      else
         do i = 1, NNZH
            HESS(i) = 1.0d0
         enddo
      endif
      IERR = 0
      return
      end subroutine EV_HESS

      subroutine ITER_CB(ALG_MODE, ITER_COUNT, OBJVAL, INF_PR, INF_DU, 
     &    MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT, 
     &    DAT, ISTOP)
      implicit none
      integer ALG_MODE, ITER_COUNT, LS_TRIAL
      double precision OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
      double precision ALPHA_DU, ALPHA_PR
      double precision DAT(*)
      integer IDAT(*)
      integer ISTOP
      logical quitexist, quitexist_

      write(*,*) 'IPOPT: iteration ', ITER_COUNT
      write(*,*) 'IPOPT: objval ', OBJVAL

!     And set ISTOP to 1 if you want Ipopt to stop now.  Below is just a
!     simple example.
      inquire(file='.quit',exist=quitexist)
      inquire(file='_quit',exist=quitexist_)
      if(quitexist)  write(*,*) 'IPOPT: stopping because .quit seen'
      if(quitexist_) write(*,*) 'IPOPT: stopping because _quit seen'
      if (quitexist.or.quitexist_) ISTOP = 1

      return
      end subroutine ITER_CB

      end program b2optim_ipopt
