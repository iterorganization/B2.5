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
      integer     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
      parameter  (N = 1, M = 1, NELE_JAC = 1, NELE_HESS = 0)
      parameter  (IDX_STY = 1 )
csc N = number design variables
c   M number of constraints (it's 2 because the third one is the variable bound!
c   NELE_JAC # of nonzeros in Jacobian of Lagrangian
c   IDX_STY index style for matrix (1 for fortran normally)
C
C     Space for multipliers and constraints
C
      double precision LAM(M)
      double precision G(M)
C
C     Vector of variables
C
      double precision X(N)
C
C     Vector of lower and upper bounds
C
      double precision X_L(N), X_U(N), Z_L(N), Z_U(N)
      double precision G_L(M), G_U(M)
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
      integer i
      real(kind=R8), save :: b2opt_elapsed=0.0_R8
      double precision, save :: cpu_ipopt=0.0
      integer, save :: b2opt_maxiter=100
C
C     The following are the Fortran routines for computing the model
C     functions and their derivatives - their code can be found further
C     down in this file.
C
      external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, ITER_CB
csc   these are the external subroutine names passed to IPOPT
CC
CC     The next is an optional callback method.  It is called once per
CC     iteration.
CC
C      external ITER_CB
C
C     Set initial point and bounds:
C
      data X   / 1.0d0/
      data X_L / 0.5d0 /
      data X_U / 2.5d0 /
C
C     Set bounds for the constraints
C
      data G_L / 0.5d0 /
      data G_U / 2.5d0 /


C
C     Allocate and initialize par_opt variables to be used in B2.5
C
      call b2mn_init_d
      npar_opt = N
      allocate(par_opt_phys(1:npar_opt))
      par_opt_phys = 0.0_R8
      optimization  = .true. !csc this will tell in b2tqna to use par_opt_phys for parm_dna

C
C     First create a handle for the Ipopt problem (and read the options
C     file)
C
      IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,
     1     IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G)
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

!     Set mu strategy
      IERR = IPADDSTROPTION(IPROBLEM, 'mu_strategy', 'adaptive')
      if (IERR.ne.0 ) goto 9990
!csc IMPORTANT HERE these are needed for BFGS updates without hessian
      IERR = IPADDSTROPTION(IPROBLEM, 'hessian_approximation',
     & 'limited-memory')
      if (IERR.ne.0 ) goto 9990
      IERR = IPADDSTROPTION(IPROBLEM, 'limited_memory_update_type',
     & 'bfgs')
      if (IERR.ne.0 ) goto 9990

!     Set max number of iterations
      call ipgeti ('b2opt_maxiter', b2opt_maxiter)
      call xertst (b2opt_maxiter.gt.0,
     1  'faulty internal parameter b2opt_maxiter')
      IERR = IPADDINTOPTION(IPROBLEM, 'max_iter', b2opt_maxiter)
      if (IERR.ne.0 ) goto 9990

!     Set a minimum tolerance
      IERR = IPADDNUMOPTION(IPROBLEM, 'tol', 1.d-7)
      if (IERR.ne.0 ) goto 9990

!     Set a maximum cpu time for ipopt based on b2mndr_elapsed
      call ipgetr ('b2opt_elapsed', b2opt_elapsed)
      call xertst (b2opt_elapsed.ge.0.0_R8,
     1  'faulty internal parameter b2opt_elapsed')
      if (b2opt_elapsed.gt.0.0_R8) then
        cpu_ipopt = DBLE(b2opt_elapsed)
        IERR = IPADDNUMOPTION(IPROBLEM, 'max_cpu_time', cpu_ipopt)
        if (IERR.ne.0 ) goto 9990
      endif

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
      do i = 1, N
        write(*,*) 'X  (',i,') = ',X(i)
      enddo
      write(*,*)
      write(*,*) 'The multipliers for the lower bounds are:'
      write(*,*)
      do i = 1, N
        write(*,*) 'Z_L(',i,') = ',Z_L(i)
      enddo
      write(*,*)
      write(*,*) 'The multipliers for the upper bounds are:'
      write(*,*)
      do i = 1, N
        write(*,*) 'Z_U(',i,') = ',Z_U(i)
      enddo
      write(*,*)
      write(*,*) 'The multipliers for the equality constraints are:'
      write(*,*)
      do i = 1, M
        write(*,*) 'LAM(',i,') = ',LAM(i)
      enddo
      write(*,*)
C
 9000 continue
C
C     Clean up
C
      call b2mn_fin_d
      call IPFREE(IPROBLEM)
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
      subroutine EV_F(N, XX, NEW_X, F, IDAT, DAT, IER)
      use b2mod_types
      use b2mod_main_diff
     , , only : b2mn_step
      use b2mod_plasma_diff
      use b2mod_user_namelist_diff
     , , only : nncf
      use b2mod_par_opt
      implicit none
      integer N, NEW_X
      double precision F, XX(N)
      real(kind=r8) j(nncf)
      double precision DAT(*)
      integer IDAT(*)
      integer IER
      save j
!csc important here to check for NEW_X! so that it doesn't recalculate with same X
      if (new_x.eq.0) then
        par_opt_phys(1) = real(XX(1),kind=r8)
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
      subroutine EV_GRAD_F(N, XX, NEW_X, GRAD, IDAT, DAT, IER)
      use b2mod_types
      use b2mod_main_diff
     , , only : b2mn_step_d
      use b2mod_user_namelist_diff
     , , only : nncf
      use b2mod_par_opt
      implicit none
      integer N, NEW_X
      double precision GRAD(N), XX(N)
      real(kind=r8) j(nncf), jd(nncf), jb(nncf)
      double precision DAT(*)
      integer IDAT(*)
      integer IER
      save jd
!csc important here to check for NEW_X! so that it doesn't recalculate with same X
! FIXME automatically choose with preprocessor data
      if (new_x.eq.0) then
        par_opt_phys(1) = real(XX(1),kind=r8)
        call b2mn_step_d(j,jd)
        grad(1) = DBLE(jd(1))
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
      subroutine EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
      use b2mod_types
      implicit none
      integer N, NEW_X, M
      double precision G(M), X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
C     s.t.  x1*x2*x3*x4                   >=  25
C           x1**2 + x2**2 + x3**2 + x4**2  =  40
!csc important here to check for NEW_X! so that it doesn't recalculate with same X
      G(1) = X(1) - DAT(1)
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
      subroutine EV_JAC_G(TASK, N, X, NEW_X, M, NZ, ACON, AVAR, A,
     1     IDAT, DAT, IERR)
      use b2mod_types
      integer TASK, N, NEW_X, M, NZ
      double precision X(N), A(NZ)
      integer ACON(NZ), AVAR(NZ), I
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
!csc important here to check for NEW_X! so that it doesn't recalculate with same X
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
C
      if( TASK.eq.0 ) then
        do I = 1, NZ
          AVAR(I) = AVAR1(I)
          ACON(I) = ACON1(I)
        enddo
      else
        A(1) = 1.0 ![1,1]
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

      write(*,*) 'IPOPT: iteration ', ITER_COUNT
      write(*,*) 'IPOPT: objval ', OBJVAL
C
C     And set ISTOP to 1 if you want Ipopt to stop now.  Below is just a
C     simple example.
C
c      if (INF_PR.le.1D-04) ISTOP = 1

      return
      end
