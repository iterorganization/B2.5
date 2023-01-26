      program b2optim_ipopt
      use b2mod_types
      use b2mod_par_opt_diff
      use b2mod_main_diff
     , , only : b2mn_init_b, b2mn_step_b, b2mn_step, b2mn_fin_b
      use b2mod_ad_diff
     , , only : nncf
      use b2us_data_diff
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
      integer ipar, icon

!     The following are the Fortran routines for computing the model
!     functions and their derivatives - their code can be found further
!     down in this file.
      external streql
      logical streql, hessian
!sc   these are the external subroutine names passed to IPOPT


!     ..START OF COMPUTATION

!     Allocate and initialize par_opt variables to be used in B2.5
      flag_optim  = .true.
      call b2mn_init_b(switch, geo, geob, mpg, mpgb, state, stateb,
     &  state_ext, state_extb)
      par_opt_phys = 0.0_R8
!     Initialize derivatives of estimated parameters
#ifdef TGT
      call xertst(npar_opt.le.nbdirsmax,
     & 'Increase size of nbdirsmax in diffsizes.F')
      call set_tgt_perturbation(switchd)
#endif

#ifdef ADJ
      par_opt_physb = 0.0_R8
#endif

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
        X(ipar) = x0(ipar)/par_rescale(ipar)
        X_L(ipar) = xl(ipar)/par_rescale(ipar)
        X_U(ipar) = xu(ipar)/par_rescale(ipar)
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
      call b2mn_fin_b(switch, geo, geob, mpg, mpgb, state, stateb,
     & state_ext, state_extb)
      call IPFREE(IPROBLEM)
      deallocate(X,X_L,X_U,Z_L,Z_U)
      deallocate(G,G_L,G_U,LAM)
      deallocate(par_opt_phys)
      deallocate(par_opt_physb)
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
          par_opt_phys(ipar) = real(XX(ipar),kind=r8)*par_rescale(ipar)
          write(str,"(I1)") ipar
          if (ipar.ge.10) write(str,"(I2)") ipar
          write(*,*) 'IPOPT: eval_f with x',trim(str),'= ',
     &     par_opt_phys(ipar)
        enddo
        isigma = npar - nsigma_opt + 1
        do ipar = 1, nsigma
          if (sigma_opt(ipar)) then
            sigma(ipar) = real(XX(isigma),kind=r8)*par_rescale(isigma)
            write(str,"(I1)") isigma
            if (ipar.ge.10) write(str,"(I2)") isigma
            write(*,*) 'IPOPT: eval_f with x',trim(str),'= ',
     &       sigma(ipar)
            isigma = isigma + 1
          endif
        enddo
        call b2mn_step(switch, geo, mpg, state, state_ext, j)
      endif
      F = DBLE(j(1))
      IER = 0
      return
      end subroutine EV_F

      subroutine EV_GRAD_F(npar, XX, NEW_X, GRAD, IDAT, DAT, IER)
      implicit none
      integer npar, NEW_X
      double precision GRAD(npar), XX(npar)
      real(kind=r8) j(nncf), jdiff(nncf), gradd(npar)
      double precision DAT(*)
      integer IDAT(*)
      integer IER, ipar, isigma
      character*3 str
      save jdiff
!csc important here to check for NEW_X! so that it doesn't recalculate with same X
! FIXME automatically choose with preprocessor data
      if (new_x.eq.0) then
        do ipar = 1, npar - nsigma_opt
          par_opt_phys(ipar) = real(XX(ipar),kind=r8)*par_rescale(ipar)
          write(str,"(I1)") ipar
          if (ipar.ge.10) write(str,"(I2)") ipar
          write(*,*) 'IPOPT: eval_grad_f with x',trim(str),'= ',
     &     par_opt_phys(ipar)
        enddo
        isigma = npar - nsigma_opt + 1
        do ipar = 1, nsigma
          if (sigma_opt(ipar)) then
            sigma(ipar) = real(XX(isigma),kind=r8)*par_rescale(isigma)
            write(str,"(I1)") isigma
            if (ipar.ge.10) write(str,"(I2)") isigma
            write(*,*) 'IPOPT: eval_grad_f with x',trim(str),'= ',
     &       sigma(ipar)
            isigma = isigma + 1
          endif
        enddo
! if forward, calculate the gradient using an 'effective' number of parameters which only includes the real physical parameters and not the sigmas
! because the gradient of the cost function wrt sigma is quite simple and only depends on the cost function. In this way we avoid iterating
! the forward problem over unnecessary directions
        call b2mn_step_b(switch, switchb, geo, geob, mpg, mpgb, state,
     &   stateb, state_ext, state_extb, j)
      endif
#ifdef TGT
      do ipar = 1, npar
        grad(ipar) = DBLE(jdiff(1))*par_rescale(ipar) ! rescale par to get order unity
        write (*,*) 'IPOPT GRAD:', grad(ipar)
      end do
#endif
#ifdef ADJ
      call set_adj_gradient(npar,gradd,switchb)
      do ipar = 1, npar
        grad(ipar) = DBLE(gradd(ipar))*par_rescale(ipar) ! rescale par to get order unity
        write (*,*) 'IPOPT GRAD:', grad(ipar)
      end do
#endif
      write (*,*) 'IPOPT GRADIENT NORM:', norm2(grad(1:npar))
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
      use b2us_data_diff
      use b2us_io_diff
     , , only : write_b2fstate
      use b2mod_version_diff
     , , only : newversion, cfverw
      use b2mod_b2cmpa_diff
      implicit none
      integer ALG_MODE, ITER_COUNT, LS_TRIAL
      double precision OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
      double precision ALPHA_DU, ALPHA_PR
      double precision DAT(*)
      integer IDAT(*)
      integer ISTOP, IDUM(0:2)
      logical quitexist, quitexist_
      character(22) :: opt_state_name
      character*120 label

      write(*,*) 'IPOPT: iteration ', ITER_COUNT
      write(*,*) 'IPOPT: objval ', OBJVAL

! Experimental: write intermediate state file?
      if (ITER_COUNT .gt. 0) then
        write(*,*) 'IPOPT: Saving intermediate optimization state'
        write (opt_state_name,'(a14,i4.4)') 'b2fstate.',ITER_COUNT
        call cfopen (99,trim(opt_state_name),'new','un*formatted')
        idum(0) = mpg%nCv
        idum(1) = mpg%nFc
        idum(2) = state%ns
        call cfverw (99, newversion)
        call cfwuin (99, 3, idum, 'nCv,nFc,ns')
        write (label,'(a46,i4)') 
     &    'b2optim_ipopt intermediate optimization state ',ITER_COUNT
        call cfwuch (99, 120, label, 'label')
        call b2wuzd_nodiff (99, newversion, state%ns, zamin, zamax, zn,
     &   am)
        call write_b2fstate (99, mpg%nCv, mpg%nFc, state%ns, state)
        close(99)
      endif

!     And set ISTOP to 1 if you want Ipopt to stop now.  Below is just a
!     simple example.
      inquire(file='.quit_opt',exist=quitexist)
      inquire(file='_quit_opt',exist=quitexist_)
      if(quitexist)  write(*,*) 'IPOPT: stopping because .quit seen'
      if(quitexist_) write(*,*) 'IPOPT: stopping because _quit seen'
      if (quitexist.or.quitexist_) ISTOP = 1

      return
      end subroutine ITER_CB

      end program b2optim_ipopt
