      program b2optim
C
      use b2mod_par_opt
      use b2mod_main_diff
     , , only : b2mn_init_d, b2mn_fin_d
      implicit none
C
C     include the Ipopt return codes
C
      include 'IpReturnCodes.inc'
C
C     Size of the problem (number of variables and equality constraints)
C

      integer :: nvar, ncon, nele_jac, NELE_HESS, IDX_STY
      parameter  (NELE_HESS = 0, IDX_STY = 1 )
csc N = number design variables
c   M number of constraints (it's 2 because the third one is the variable bound!
c   NELE_JAC # of nonzeros in Jacobian of Lagrangian
c   IDX_STY index style for matrix (1 for fortran normally)
C
C     Space for multipliers and constraints
C
      double precision, allocatable :: LAM(:), G(:)
C
C     Vector of variables
C
      double precision, allocatable :: X(:)
C
C     Vector of lower and upper bounds
C
      double precision, allocatable :: X_L(:), X_U(:), Z_L(:), Z_U(:)
      double precision, allocatable :: G_L(:), G_U(:)
Csc  X_L/X_U lower upper bound for each desiign variable
c    Z_L, Z_U initial value for the bound multiplier --> seems not needed to be initialized here
c    G_L,G_U lower/upper bounds for constraints
C
C     Private data for evaluation routines
C     This could be used to pass double precision and integer arrays untouched
C     to the evaluation subroutines EVAL_*
C
      double precision DAT(2)
      integer IDAT(1)
C
C     Place for storing the Ipopt Problem Handle
C
      integer*8 IPROBLEM
      integer*8 IPCREATE
C
      integer IERR
      integer IPSOLVE, IPADDSTROPTION
      integer IPADDNUMOPTION, IPADDINTOPTION
      integer IPOPENOUTPUTFILE
C
      double precision F
      integer ivar, icon
C
C     The following are the Fortran routines for computing the model
C     functions and their derivatives - their code can be found further
C     down in this file.
C
      external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, ITER_CB
!sc   these are the external subroutine names passed to IPOPT


!     ..START OF COMPUTATION

      call read_b2mod_par_opt(nvar, ncon, nele_jac)
      write (*,*) 'IPOPT, WARNING: setting only 1 hard-coded '//
     & 'trivial constraint'
      ncon = 1 ! FIXME to be updated whenever nonlinear constraints are needed
      nele_jac = 1

      allocate(X(nvar),X_L(nvar),X_U(nvar),Z_L(nvar),Z_U(nvar))
      allocate(G(ncon),G_L(ncon),G_U(ncon),LAM(ncon))

      do ivar = 1, nvar
        X(ivar) = x0(ivar)
        X_L(ivar) = xl(ivar)
        X_U(ivar) = xu(ivar)
        if (X_L(ivar).lt.-inf_opt)
     &    write(*,*) 'IPOPT: warning, X_L(',ivar,') set to infty'
        if (X_U(ivar).gt.inf_opt)
     &    write(*,*) 'IPOPT: warning, X_U(',ivar,') set to infty'
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


!     Allocate and initialize par_opt variables to be used in B2.5

      call b2mn_init_d
      npar_opt = nvar
      allocate(par_opt_phys(1:npar_opt))
      par_opt_phys = 0.0_R8
      flag_optim  = .true. !csc this will tell in b2tqna to use par_opt_phys for parm_dna

!     First create a handle for the Ipopt problem (and read the options file)
      IPROBLEM = IPCREATE(nvar, X_L, X_U, ncon, G_L, G_U, NELE_JAC, 
     &     NELE_HESS, IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G)
      if (IPROBLEM.eq.0) then
         write(*,*) 'Error creating an Ipopt Problem handle.'
         stop
      endif
C
C     Open an output file
C
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
c      IERR = IPADDSTROPTION(IPROBLEM, 'nlp_scaling_method',
c     & 'user-scaling')
c      if (IERR.ne.0 ) goto 9990

!csc IMPORTANT HERE these are needed for BFGS updates without hessian
      IERR = IPADDSTROPTION(IPROBLEM, 'hessian_approximation',
     & hessian_approximation)
      if (IERR.ne.0 ) goto 9990
      IERR = IPADDSTROPTION(IPROBLEM, 'limited_memory_update_type',
     & limited_memory_update_type)
      if (IERR.ne.0 ) goto 9990

!     Set max number of iterations
      call xertst (maxiter.gt.0,
     1  'faulty internal parameter maxiter')
      IERR = IPADDINTOPTION(IPROBLEM, 'max_iter', maxiter)
      if (IERR.ne.0 ) goto 9990

c      IERR = IPADDINTOPTION(IPROBLEM, 'print_level', 12)
c      if (IERR.ne.0 ) goto 9990

!     Set a minimum tolerance
      call xertst (tol_opt.gt.0.0,
     1  'faulty internal parameter tol_opt')
      IERR = IPADDNUMOPTION(IPROBLEM, 'tol', tol_opt)
      if (IERR.ne.0 ) goto 9990

!     Set a maximum cpu time for ipopt based on b2mndr_elapsed
      call xertst (cpu_opt.ge.0.0_R8,
     1  'faulty internal parameter cpu_opt')
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

CC     Set a callback function to give you control once per iteration.
CC     You can use it if you want to generate some output, or to stop
CC     the optimization early.
CC
      call IPSETCALLBACK(IPROBLEM, ITER_CB)

C
C     As a simple example, we pass the constants in the constraints to
C     the EVAL_C routine via the "private" DAT array.
C
      DAT(1) = 0.d0
      DAT(2) = 0.d0
C
C     Call optimization routine
C
      IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)
C
C     Output:
C
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
      do ivar = 1, nvar
        write(*,*) 'X  (',ivar,') = ',X(ivar)
      enddo
      write(*,*)
      write(*,*) 'The multipliers for the lower bounds are:'
      write(*,*)
      do ivar = 1, nvar
        write(*,*) 'Z_L(',ivar,') = ',Z_L(ivar)
      enddo
      write(*,*)
      write(*,*) 'The multipliers for the upper bounds are:'
      write(*,*)
      do ivar = 1, nvar
        write(*,*) 'Z_U(',ivar,') = ',Z_U(ivar)
      enddo
      write(*,*)
      write(*,*) 'The multipliers for the equality constraints are:'
      write(*,*)
      do icon = 1, ncon
        write(*,*) 'LAM(',icon,') = ',LAM(icon)
      enddo
      write(*,*)
C
 9000 continue
C
C     Clean up
C
      call b2mn_fin_d
      call IPFREE(IPROBLEM)
      deallocate(X,X_L,X_U,Z_L,Z_U)
      deallocate(G,G_L,G_U,LAM)
      deallocate(par_opt_phys)
      stop 'b2optim'
C
 9990 continue
      write(*,*) 'Error setting an option'
      goto 9000
      end program
C
C =============================================================================
C
C                    Computation of objective function
C
C =============================================================================
C
      subroutine EV_F(nvar, XX, NEW_X, F, IDAT, DAT, IER)
      use b2mod_types
      use b2mod_main_diff
     , , only : b2mn_step
      use b2mod_plasma_diff
      use b2mod_user_namelist_diff
     , , only : nncf
      use b2mod_par_opt
      implicit none
      integer nvar, NEW_X
      double precision F, XX(nvar)
      real(kind=r8) j(nncf)
      double precision DAT(*)
      integer IDAT(*)
      integer IER
      save j
!csc important here to check for NEW_X! so that it doesn't recalculate with same X
      if (new_x.eq.0) then
        par_opt_phys(1) = real(XX(1),kind=r8)
        write(*,*) 'IPOPT: eval_f with x= ', par_opt_phys(1)
        call b2mn_step(j)
      endif
      F = DBLE(j(1))
      IER = 0
      return
      end subroutine EV_F
C
C =============================================================================
C
C                Computation of gradient of objective function
C
C =============================================================================
C
      subroutine EV_GRAD_F(nvar, XX, NEW_X, GRAD, IDAT, DAT, IER)
      use b2mod_types
      use b2mod_main_diff
     , , only : b2mn_step_d
      use b2mod_user_namelist_diff
     , , only : nncf
      use b2mod_par_opt
      implicit none
      integer nvar, NEW_X
      double precision GRAD(nvar), XX(nvar)
      real(kind=r8) j(nncf), jd(nncf), jb(nncf)
      double precision DAT(*)
      integer IDAT(*)
      integer IER
      save jd
!csc important here to check for NEW_X! so that it doesn't recalculate with same X
! FIXME automatically choose with preprocessor data
      if (new_x.eq.0) then
        par_opt_phys(1) = real(XX(1),kind=r8)
        write(*,*) 'IPOPT: eval_grad_f with x= ', par_opt_phys(1)
        call b2mn_step_d(j,jd)
      endif
      grad(1) = DBLE(jd(1))

c      if (new_x.eq.0) then
c        call b2mn_init_b
c        ss = real(XX(1),kind=r8)
c        call b2mn_step_b(j, jb)
c        call b2mn_fin_b
c        grad(1) = DBLE(ssb)
c      endif
c      grad(1) = DBLE(ssb)

      IER = 0
      return
      end subroutine EV_GRAD_F

C
C =============================================================================
C
C                     Computation of equality constraints
C
C =============================================================================
C
      subroutine EV_G(nvar, X, NEW_X, ncon, G, IDAT, DAT, IERR)
      use b2mod_types
      implicit none
      integer nvar, NEW_X, ncon
      double precision G(ncon), X(nvar)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR, icon
C     s.t.  x1*x2*x3*x4                   >=  25
C           x1**2 + x2**2 + x3**2 + x4**2  =  40
!csc WARNING to be updated whenever nonlinear constraints are needed
      G(1) = -1.0

      IERR = 0
      return
      end subroutine EV_G
C
C =============================================================================
C
C                Computation of Jacobian of equality constraints
C
C =============================================================================
C
      subroutine EV_JAC_G(TASK, nvar, X, NEW_X, ncon, NZ, ACON, AVAR, A,
     1     IDAT, DAT, IERR)
      use b2mod_types
      use b2mod_par_opt
      integer TASK, nvar, NEW_X, ncon, NZ
      double precision X(nvar), A(NZ)
      integer ACON(NZ), AVAR(NZ), I
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
C
C     structure of Jacobian:
C
C     s.t.  x1*x2*x3*x4                   >=  25
C           x1**2 + x2**2 + x3**2 + x4**2  =  40
csc  [x2*x3*x4   x1*x3*x4   x1*x2*x4   x1*x2*x3;
c         2*x1       2*x2       2*x3       2*x4]
      integer AVAR1(1), ACON1(1)
      data  AVAR1 /1/ !ROWS
      data  ACON1 /1/ !COLUMNS
      save  AVAR1, ACON1
!csc WARNING to be updated whenever nonlinear constraints are needed
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

C
C =============================================================================
C
C                   Callback method called once per iteration
C
C =============================================================================
C
      subroutine ITER_CB(ALG_MODE, ITER_COUNT, OBJVAL, INF_PR, INF_DU,
     1     MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT,
     2     DAT, ISTOP)
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
C
C     And set ISTOP to 1 if you want Ipopt to stop now.  Below is just a
C     simple example.
C
      inquire(file='.quit',exist=quitexist)
      inquire(file='_quit',exist=quitexist_)
      if(quitexist)  write(*,*) 'IPOPT: stopping because .quit seen'
      if(quitexist_) write(*,*) 'IPOPT: stopping because _quit seen'
      if (quitexist.or.quitexist_) ISTOP = 1

      return
      end
